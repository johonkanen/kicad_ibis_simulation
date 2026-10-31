"""Fast sweeps: write the SPICE netlist directly, convert each IBIS buffer only once, build the bit patterns in Python.

Why: kicad-cli's IBIS converter (KIBIS) costs ~15 s per driver at 150 bits and grows steeply with the bit count, and
a sweep re-converts the same buffers at every width / length point although only the lines change. Here:

1. Netlist: a `Netlist` collects R/L/C, TLINE, via and IBIS elements and writes the SPICE text itself (no schematic,
   no kicad-cli per point). The elements are the same as `sch.Sch` writes, so a committed baseline schematic made
   with `Sch` still matches what was simulated.
2. KIBIS cache: every IBIS element has a key (IBIS file content, component, pin, model, corner, bit rate). Missing keys
   are converted together, once, by kicad-cli on a small schematic (`KibisCache.resolve`); the subcircuits are stored
   in the cache folder and reused by every netlist that needs them.
3. Short conversion, own pattern: drivers are converted with a short PRBS (CONV_BITS). One rising and one falling edge
   segment of KIBIS' own switching coefficients Ku/Kd are then laid at every transition of the wanted pattern
   (`synth_kukd`): the KIBIS PRBS7 of any length, its complement (the N leg of a pair: KiCad 10's differential IBIS
   drivers come out with empty Ku/Kd), a shifted copy or a constant level (crosstalk aggressors). Built patterns are
   cached too.
4. Lines shorter than LUMP_TD (15 ps: neck-down pieces, stubs, short layer runs) become damped L/C pi sections of at
   most 2.5 ps (edges are 50+ ps). ngspice's lossless line (T) with a few-ps delay stops the run with "timestep too
   small" or makes it crawl: a CKE trunk of nine 3-13 ps lines stopped at 20 ns of 102 ns; lumped, it ran in 5 s.

Usage:
    n = Netlist('my point')
    n.ibis('U1', 'drv', 'models/fpga.ibs', 'FPGA', 'T62', 'out_40', ['dq_tx'], 'f0=2.8G n=150')
    n.tl('dq_tx', 'dq_rx', 40, 140e-12)
    n.ibis('U2', 'rx', 'models/dram.ibs', 'DRAM', 'C11', 'ODT_60', ['dq_rx'])
    n.directive('.tran 1p 55.6n')
    cache = KibisCache('work/kibis')
    cache.resolve(n.keys())                 # kicad-cli once for whatever is not cached yet
    n.write('work/point.cir', cache)        # die node of U2: xu2.die0, as with kicad-cli netlists

Validated against the kicad-cli path: same eyes to < 1 mV / 1 ps (see the board project's README).
"""
import hashlib, json, os, re, shutil, subprocess, sys, tempfile
import numpy as np

from . import spice, eye

CONV_BITS = 40          # PRBS bits KIBIS converts per driver: enough for 3-bit runs of both polarities
LUMP_TD = 15e-12        # lines shorter than this (s) -> damped L/C pi sections of <= LUMP_SEC each
LUMP_SEC = 1.0e-12     # section length
LUMP_RD = 4.0          # damping R = LUMP_RD * Z0 across each L: kills the ~160 GHz tank, >> |jwL| in the signal band


def _eng(v):
    return f'{v:.6g}'


# ----------------------------------------------------------------------------------------------------- pattern synth
def pwl_of(text, which):
    m = re.search(rf'^Vk{which} K{which.upper()} GND pwl \((.*?)\)', text, re.M | re.S)
    a = np.array(m.group(1).split(), float).reshape(-1, 2)
    t, idx = np.unique(a[::-1, 0], return_index=True)
    return t, a[::-1, 1][idx]


def source_bits(text, f0):
    """The bit pattern a KIBIS driver subcircuit was converted with, read back from its Ku (pull-up on = 1)."""
    tu, ku = pwl_of(text, 'u')
    n = int(round(tu[-1] * f0))
    return (np.interp((np.arange(n) + 0.98) / f0, tu, ku) > 0.5).astype(int)


def synth_kukd(text, f0, new_bits, dt=0.5e-12):
    """Rewrite a KIBIS driver's Ku/Kd for new_bits (bit k from k/f0). One rising and one falling segment, 3 bits
    long from the bit boundary, are cut from KIBIS' coefficients for the pattern it converted, and laid at every
    transition of new_bits; a segment is cut off at the next transition, as KIBIS does, and held after 3 bits."""
    ui = 1 / f0
    tu, ku = pwl_of(text, 'u'); td_, kd = pwl_of(text, 'd')
    src = source_bits(text, f0); n_src = len(src)
    seg_n = int(round(3 * ui / dt)); tt = np.arange(seg_n) * dt
    segs = {}
    for k in range(1, n_src - 3):
        e = 'r' if src[k - 1] == 0 and src[k] == 1 else 'f' if src[k - 1] == 1 and src[k] == 0 else None
        if e and e not in segs and src[k] == src[k + 1] == src[k + 2]:
            segs[e] = (np.interp(k * ui + tt, tu, ku), np.interp(k * ui + tt, td_, kd))
    if len(segs) < 2: raise RuntimeError('no 3-bit runs of both polarities in the converted pattern')
    n = len(new_bits); T = np.arange(int(round(n * ui / dt)) + 1) * dt
    U = np.empty_like(T); D = np.empty_like(T)
    lvl = {1: (segs['r'][0][-1], segs['r'][1][-1]), 0: (segs['f'][0][-1], segs['f'][1][-1])}
    U[:], D[:] = lvl[int(new_bits[0])]
    k = 0
    while k < n:
        j = k + 1
        while j < n and new_bits[j] == new_bits[k]: j += 1
        if k > 0:
            su, sd = segs['r' if new_bits[k] else 'f']
            i0, i1 = int(round(k * ui / dt)), min(int(round(j * ui / dt)), len(T))
            L = i1 - i0; m = min(L, seg_n)
            U[i0:i0 + m], D[i0:i0 + m] = su[:m], sd[:m]
            if L > m: U[i0 + m:i1], D[i0 + m:i1] = su[-1], sd[-1]
        k = j
    fmt = lambda t, v: ' '.join(f'{x:.6e} {y:.6e}' for x, y in zip(t, v))
    text = re.sub(r'^(Vku KU GND pwl \().*?\)', lambda mm: mm.group(1) + fmt(T, U) + ' )', text, flags=re.M | re.S)
    text = re.sub(r'^(Vkd KD GND pwl \().*?\)', lambda mm: mm.group(1) + fmt(T, D) + ' )', text, flags=re.M | re.S)
    return text


def kibis_pattern(text, f0, n):
    """The n-bit PRBS7 KIBIS would have sent for a PRBSDRIVER with n bits, continued from the converted one (the
    converted pattern is matched against eye.kibis_prbs7's rotations, so labels in eye.stimulus_bits stay valid)."""
    src = source_bits(text, f0)
    for s in range(-8, 9):
        if np.array_equal(eye.kibis_prbs7(len(src), s), src): return eye.kibis_prbs7(n, s)
    raise RuntimeError('converted driver pattern is not the KIBIS PRBS7')


# ---------------------------------------------------------------------------------------------------------- netlist
class Netlist:
    def __init__(self, title='fastnet'):
        self.title = title; self.lines = []; self.ib = []; self.directives = []; self.k = 0

    def _ref(self, p):
        self.k += 1; return f'{p}{self.k}'

    def r(self, a, b, v): self.lines.append(f'R{self._ref("")} {a} {b} {_eng(v)}')

    def l(self, a, b, v): self.lines.append(f'L{self._ref("")} {a} {b} {_eng(v)}')

    def c(self, a, b, v): self.lines.append(f'C{self._ref("")} {a} {b} {_eng(v)}')

    def tl(self, a, b, z0, td):
        """Lossless line; shorter than LUMP_TD: an L/C pi section with the same L = z0 td and C = td / z0, with z0
        across the L. Without it the pi section is an undamped tank near 1/(2 pi td) (100+ GHz): every edge rings it
        and ngspice crawls through the whole run at sub-ps steps. At signal frequencies (a few GHz) the L is a few
        ohms, so the resistor in parallel changes nothing there."""
        if td < LUMP_TD:
            n = max(1, int(np.ceil(td / LUMP_SEC))); d = td / n; nodes = [a] + [f'{a}_s{self.k}_{i}' for i in range(1, n)] + [b]
            for x, y in zip(nodes, nodes[1:]):
                self.c(x, '0', d / z0 / 2); self.l(x, y, z0 * d); self.r(x, y, LUMP_RD * z0); self.c(y, '0', d / z0 / 2)
        else:
            self.lines.append(f'T{self._ref("")} {a} 0 {b} 0 Z0={_eng(z0)} TD={_eng(td)}')

    def via(self, a, b, L, C):
        self.c(a, '0', C / 2); self.l(a, b, L); self.c(b, '0', C / 2)

    def coupled(self, nodes_in, nodes_out, Lm, Cm, length, rdc=0.0, dl=0.25e-3):
        """n coupled lines between nodes_in[i] and nodes_out[i] (length in m) as a lumped ladder of dl-long sections
        (0.25 mm, ~1.75 ps on a stripline: far shorter than any LPDDR4 edge): per section the self L with mutual K,
        C to ground (row sum of the Maxwell matrix) and C between lines, rdc (ohm/m) in series. Lm (H/m) and Cm (F/m,
        Maxwell: diagonal total, off-diagonal negative) from a 2D field solve. ngspice's CPL element stopped every
        LPDDR4 run with "timestep too small"; the ladder runs."""
        n = len(nodes_in); secs = max(1, int(round(length / dl))); d = length / secs
        prev = list(nodes_in)
        for k in range(secs):
            self.k += 1; tag = f'cp{self.k}'
            nxt = list(nodes_out) if k == secs - 1 else [f'{tag}_{i}' for i in range(n)]
            for i in range(n):
                mid = f'{tag}_r{i}'
                self.lines.append(f'R{tag}_{i} {prev[i]} {mid} {_eng(max(rdc * d, 1e-6))}')
                self.lines.append(f'L{tag}_{i} {mid} {nxt[i]} {_eng(Lm[i][i] * d)}')
                self.lines.append(f'C{tag}_g{i} {nxt[i]} 0 {_eng(sum(Cm[i]) * d)}')
            for i in range(n):
                for j in range(i + 1, n):
                    self.lines.append(f'K{tag}_{i}{j} L{tag}_{i} L{tag}_{j} {Lm[i][j] / (Lm[i][i] * Lm[j][j]) ** 0.5:.6f}')
                    self.lines.append(f'C{tag}_m{i}{j} {nxt[i]} {nxt[j]} {_eng(-Cm[i][j] * d)}')
            prev = nxt

    def directive(self, s): self.directives.append(s)

    def ibis(self, ref, kind, lib, comp, pin, model, nets, params='', stype=None, pattern=None):
        """kind 'drv' (PRBSDRIVER unless stype) or 'rx' (IBISDEVICE). nets: pin nets (one, or two for a diff
        device). pattern: None = the KIBIS PRBS7 of n bits (params n=), or ('inv', ref of another driver),
        ('shift', k), ('const', 0|1)."""
        st = stype or ('PRBSDRIVER' if kind == 'drv' else 'IBISDEVICE')
        self.ib.append(dict(ref=ref, kind=kind, lib=os.path.abspath(lib), comp=comp, pin=pin, model=model,
                            nets=list(nets), params=params, stype=st, pattern=pattern))

    def keys(self):
        return [KibisCache.key(e) for e in self.ib]

    def write(self, path, cache, tol=1e-3):
        """SPICE text with every IBIS element's subcircuit inline (drivers: their own pattern)."""
        out = [f'* {self.title} (kicad_ibis.fastnet)']
        subs, pats = {}, {}
        for e in self.ib:     # plain PRBS first, so ('inv', ref) can use the pattern of its reference
            if e['kind'] == 'drv' and not e['pattern']:
                f0, n = _f0_n(e['params']); text = cache.get(KibisCache.key(e))
                pats[e['ref']] = kibis_pattern(text, f0, n)
        for e in self.ib:
            k = KibisCache.key(e); text = cache.get(k)
            if e['kind'] == 'drv' and 'Vku KU' in text:
                f0, n = _f0_n(e['params']); p = e['pattern']
                bits = pats.get(e['ref']) if not p else 1 - pats[p[1]] if p[0] == 'inv' else \
                    np.roll(kibis_pattern(text, f0, n), p[1]) if p[0] == 'shift' else np.full(n, p[1])
                name = f'{e["ref"]}_s'
                body = cache.pattern(k, bits, f0, tol)
            else:
                name = f'K{k[:10]}'; body = None if name in subs else cache.clean(k)
            if body is not None:
                subs[name] = re.sub(r'^\.SUBCKT \S+', f'.SUBCKT {name}', body, count=1, flags=re.M | re.I)
            out.append(f'X{e["ref"]} 0 {" ".join(e["nets"])} {name}')
        txt = '\n'.join([out[0]] + list(subs.values()) + out[1:] + self.lines + self.directives + ['.end']) + '\n'
        open(path, 'w').write(txt)
        return path


def _f0_n(params):
    d = dict(re.findall(r'(\w+)=(\S+)', params))
    val = lambda s: float(re.match(r'[\d.eE+-]+', s).group()) * {'G': 1e9, 'M': 1e6, 'k': 1e3}.get(
        re.sub(r'[\d.eE+-]+', '', s)[:1], 1) if not s.endswith('Meg') else float(s[:-3]) * 1e6
    return val(d['f0']), int(d['n'])


# -------------------------------------------------------------------------------------------------------- the cache
class KibisCache:
    """KIBIS subcircuits by element key, in a folder (safe for parallel readers; writes are atomic renames)."""

    def __init__(self, root):
        self.root = os.path.abspath(root); os.makedirs(self.root, exist_ok=True)

    @staticmethod
    def key(e):
        h = hashlib.sha1(open(e['lib'], 'rb').read()).hexdigest()
        params = re.sub(r'\bn=\d+', f'n={CONV_BITS}', e['params']) if e['stype'] == 'PRBSDRIVER' else e['params']
        sig = json.dumps([h, e['comp'], e['pin'], e['model'], e['stype'], params.strip(), len(e['nets'])])
        return hashlib.sha1(sig.encode()).hexdigest()

    def path(self, k, ext='raw'): return os.path.join(self.root, f'{k}.{ext}')

    def get(self, k): return open(self.path(k)).read()

    def _put(self, p, text):
        fd, tmp = tempfile.mkstemp(dir=self.root); os.write(fd, text.encode()); os.close(fd); os.replace(tmp, p)

    def clean(self, k):
        p = self.path(k, 'clean')
        if not os.path.exists(p): self._put(p, spice.clean_cache(self.get(k)))
        return open(p).read()

    def pattern(self, k, bits, f0, tol=1e-3):
        h = hashlib.sha1(np.asarray(bits, np.int8).tobytes() + f'{f0:.6g}'.encode()).hexdigest()[:16]
        p = self.path(k, f'{h}.pat')
        if not os.path.exists(p): self._put(p, spice.clean_cache(synth_kukd(self.get(k), f0, bits), tol=tol))
        return open(p).read()

    def resolve(self, entries, batch=24, jobs=4):
        """Convert every element not cached yet: kicad-cli on small schematics (batch elements each, jobs at once,
        each with its own XDG cache folder: KIBIS writes to $XDG_CACHE_HOME/kicad/10.0/ibis)."""
        todo, seen = [], set()
        for e in entries:
            k = KibisCache.key(e)
            if k not in seen and not os.path.exists(self.path(k)): seen.add(k); todo.append((k, e))
        if not todo: return 0
        from concurrent.futures import ThreadPoolExecutor
        chunks = [todo[i:i + batch] for i in range(0, len(todo), batch)]
        with ThreadPoolExecutor(jobs) as ex:
            for msgs in ex.map(self._convert, chunks):
                if msgs: print('kicad-cli:', '; '.join(msgs), flush=True)
        missing = [k for k, _ in todo if not os.path.exists(self.path(k))]
        if missing: sys.exit(f'KIBIS conversion failed for {len(missing)} element(s)')
        return len(todo)

    def _convert(self, chunk):
        from .sch import Sch
        d = tempfile.mkdtemp(prefix='conv_', dir=self.root)
        try:
            s = Sch('.')
            s.text('.tran 1p 1n', 25.4, 15.24)
            refs = {}
            for k, e in chunk:
                s.row(k[:8])
                params = re.sub(r'\bn=\d+', f'n={CONV_BITS}', e['params']) if e['stype'] == 'PRBSDRIVER' else e['params']
                ref = s.ibis('drv' if e['kind'] == 'drv' else 'rx', (os.path.relpath(e['lib'], d), e['comp']),
                             e['pin'], e['model'],
                             f'n_{k[:8]}', params, e['stype'])
                refs[ref] = k
            sch = os.path.join(d, 'conv.kicad_sch'); s.write(sch)
            env = dict(os.environ, XDG_CACHE_HOME=os.path.join(d, 'xdg'))
            msgs = []
            for attempt in range(3):     # kicad-cli instances started at the same moment can fail silently
                r = subprocess.run(spice.KICAD_CLI + ['sch', 'export', 'netlist', '--format', 'spice', '-o', 'conv.cir',
                                                      'conv.kicad_sch'], cwd=d, capture_output=True, text=True, env=env)
                msgs = sorted(set(l for l in (r.stdout + r.stderr).splitlines() if l.strip()
                                  and 'Bitrate could be too high' not in l))
                if not r.returncode and not any('rror' in m or 'Invalid' in m for m in msgs): break
            cdir = os.path.join(d, 'xdg', 'kicad', '10.0', 'ibis')
            for ref, k in refs.items():
                f = os.path.join(cdir, f'{ref}.cache')
                if os.path.exists(f): self._put(self.path(k), open(f).read())
            return msgs
        finally:
            shutil.rmtree(d, ignore_errors=True)
