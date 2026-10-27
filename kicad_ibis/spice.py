"""KiCad schematic -> SPICE netlist (kicad-cli, which runs KiCad's IBIS converter KIBIS) -> ngspice (libngspice).

The KiCad 10 AppImage, extracted once (`kicad-10.x.AppImage --appimage-extract`), carries both kicad-cli and
libngspice with its code models; point KICAD_APPRUN at its AppRun. A system KiCad works too: set KICAD_CLI and
NGSPICE_LIB instead.
"""
import ctypes, glob, json, os, re, shutil, subprocess, sys
import numpy as np

APPRUN = os.environ.get('KICAD_APPRUN', '')
KICAD_CLI = [APPRUN, 'kicad-cli'] if APPRUN else [os.environ.get('KICAD_CLI', 'kicad-cli')]
LIBNG = os.environ.get('NGSPICE_LIB') or os.path.join(os.path.dirname(APPRUN), 'shared/lib/libngspice.so.0')
CODEMODELS = os.environ.get('NGSPICE_CM') or os.path.join(os.path.dirname(LIBNG), 'ngspice')
KSYM = os.environ.get('KICAD_SYMBOLS') or os.path.join(os.path.dirname(APPRUN), 'usr/share/kicad/symbols')
IBIS_CACHE = os.path.expanduser(os.environ.get('KIBIS_CACHE', '~/.cache/kicad/10.0/ibis'))


def check():
    """Exit with a hint when the KiCad/ngspice paths are not set up."""
    if not os.path.exists(LIBNG): sys.exit(f'libngspice not found ({LIBNG}): set KICAD_APPRUN to the extracted '
                                           'KiCad 10 AppImage AppRun, or NGSPICE_LIB')


# --------------------------------------------------------------------------------------------- KIBIS cache clean-up
def simplify(t, v, tol):
    """Douglas-Peucker on a PWL (iterative)."""
    keep = np.zeros(len(t), bool); keep[[0, -1]] = True
    stack = [(0, len(t) - 1)]
    while stack:
        a, b = stack.pop()
        if b - a < 2: continue
        tt, vv = t[a + 1:b], v[a + 1:b]
        lin = v[a] + (v[b] - v[a]) * (tt - t[a]) / (t[b] - t[a])
        i = int(np.argmax(np.abs(vv - lin)))
        if abs(vv[i] - lin[i]) > tol:
            m = a + 1 + i; keep[m] = True; stack += [(a, m), (m, b)]
    return t[keep], v[keep]


def clean_cache(text, dt=0.5e-12, tol=1e-3):
    """KIBIS writes the switching coefficients Ku/Kd as PWL sources with ~1e5 points per driver (0.8 ps steps, every
    IBIS waveform sample), which makes ngspice crawl, and with narrow solver spikes (to -0.9 / +1.3 for a few ps where
    the two-waveform system is ill-conditioned). Resample to 0.5 ps, median-filter over 2.5 ps, simplify to 1e-3."""
    def fix(m):
        a = np.array(m.group(2).split(), float).reshape(-1, 2)
        t, v = a[:, 0], a[:, 1]
        t, idx = np.unique(t[::-1], return_index=True); v = v[::-1][idx]       # last value wins at repeated times
        tu = np.arange(0, t[-1] + dt / 2, dt); vu = np.interp(tu, t, v)
        w = np.lib.stride_tricks.sliding_window_view(np.pad(vu, 2, mode='edge'), 5)
        ts, vs = simplify(tu, np.median(w, axis=1), tol)
        return m.group(1) + ' '.join(f'{x:.6e} {y:.6e}' for x, y in zip(ts, vs)) + ' )'
    return re.sub(r'^(Vk[ud] K[UD] GND pwl \()(.*?)\)', fix, text, flags=re.M)


def strip_package(path):
    """Zero the pin R/L/C in a copied KIBIS subcircuit (for checking a driver against its die-level V-t tables)."""
    s = open(path).read()
    s = re.sub(r'^(RPIN \S+ \S+) \S+', r'\1 1e-6', s, flags=re.M)
    s = re.sub(r'^(LPIN \S+ \S+) \S+', r'\1 1e-15', s, flags=re.M)
    s = re.sub(r'^(CPIN \S+ \S+) \S+', r'\1 1e-18', s, flags=re.M)
    open(path, 'w').write(s)


# ------------------------------------------------------------------------------------------------------------ export
def export(sch, raw_kukd=False):
    """sch (.kicad_sch) -> sch.cir next to it, with KIBIS' subcircuits copied into <name>_ibis/ and cleaned.

    kicad-cli runs in the schematic's folder, so Sim.Library paths resolve relative to it. KIBIS writes its
    subcircuits into one shared cache folder, so exports must not run in parallel (ngspice runs can). Returns the
    converter messages (only the last parser message survives in the headless CLI)."""
    d, name = os.path.split(os.path.abspath(sch)); tag = os.path.splitext(name)[0]
    for f in glob.glob(os.path.join(IBIS_CACHE, '*.cache')): os.remove(f)
    r = subprocess.run(KICAD_CLI + ['sch', 'export', 'netlist', '--format', 'spice', '-o', tag + '.cir', name],
                       cwd=d, capture_output=True, text=True)
    msgs = sorted(set(l for l in (r.stdout + r.stderr).splitlines() if l.strip()))
    if r.returncode or any('rror' in m or 'Invalid' in m for m in msgs):
        sys.exit(f'{tag}: kicad-cli failed:\n' + '\n'.join(msgs))
    idir = os.path.join(d, tag + '_ibis'); shutil.rmtree(idir, ignore_errors=True); os.makedirs(idir)
    cir = open(os.path.join(d, tag + '.cir')).read()
    for f in glob.glob(os.path.join(IBIS_CACHE, '*.cache')):
        dst = os.path.join(idir, os.path.basename(f))
        s = open(f).read()
        open(dst, 'w').write(s if raw_kukd else clean_cache(s)); cir = cir.replace(f, dst)
    cir = re.sub(r'(?<=\s)/(\w+)', r'\1', cir)            # KiCad net names '/x' -> 'x' (ngspice reads '/' as divide)
    open(os.path.join(d, tag + '.cir'), 'w').write(cir)
    return [m for m in msgs if m and 'Bitrate could be too high' not in m]


# ----------------------------------------------------------------------------------------------------------- ngspice
def ngspice(cir, vecs, out):
    """Run one netlist in a fresh process through libngspice (ctypes); write the vectors to out (ngspice wrdata,
    one time column, header with the vector names) and return (t, {vec: v}). Vectors are node names, e.g. 'dq_rx'
    or 'xu2.die0' (the die node inside a KIBIS subcircuit)."""
    check()
    code = f'''
import ctypes, glob, os
ng = ctypes.CDLL({LIBNG!r})
log = []
PF = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_void_p)
EX = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_int, ctypes.c_bool, ctypes.c_bool, ctypes.c_int, ctypes.c_void_p)
pf = PF(lambda s, i, u: log.append(s.decode(errors='replace')) or 0); ex = EX(lambda *a: 0)
ng.ngSpice_Init(pf, pf, ex, None, None, None, None)
c = lambda s: ng.ngSpice_Command(s.encode())
for cm in glob.glob(os.path.join({CODEMODELS!r}, '*.cm')): c('codemodel ' + cm)
c('set wr_singlescale'); c('set wr_vecnames'); c('source {os.path.abspath(cir)}'); c('run')
c('wrdata {os.path.abspath(out)} ' + {' '.join(f'v({v})' for v in vecs)!r})
open({os.path.abspath(out) + '.log'!r}, 'w').write('\\n'.join(log))
'''
    r = subprocess.run([sys.executable, '-I', '-c', code], capture_output=True, text=True)
    if r.returncode or not os.path.exists(out): sys.exit(f'ngspice failed on {cir}: {r.stderr[-2000:]}')
    return load(out, vecs)


def load(dat, vecs=None):
    """Read an ngspice wrdata file (set wr_singlescale, wr_vecnames) -> (t, {vec: v}). Names come from the header
    (v(x) -> x) unless vecs is given."""
    head = open(dat).readline().split()
    d = np.loadtxt(dat, skiprows=1)
    names = vecs or [re.sub(r'^v\((.*)\)$', r'\1', h) for h in head[1:]]
    return d[:, 0], {v: d[:, i + 1] for i, v in enumerate(names)}
