"""Make vendor IBIS files load in KiCad 10's IBIS parser (KIBIS), and read their V-t tables for validation.

    python3 -m kicad_ibis.ibisfix IN.ibs OUT.ibs          # pin header case, [R Series], package matrix -> pin RLC
    python3 -m kicad_ibis.ibisfix --corner min|max IN.ibs OUT.ibs   # min or max column as typ (corner_as_typ)

Problems found in KiCad 10.0.7 (all silent or with one terse message in the headless CLI):
  * [Pin] column names are matched case sensitively: 'R_Pin L_Pin C_Pin' (Winbond) is not read -> fix_pin_header
  * [R Series] models (and [Series Pin Mapping]) are rejected -> drop_series
  * [Define Package Model] matrices are ignored, so a file whose [Package] is all zero simulates without package
    -> package_to_pins writes the matrix diagonal (self R/L/C, no coupling) into the [Pin] columns
  * [Algorithmic Model] (AMI) blocks: analog simulation only; extract_models drops them
  * the corner keywords (vcc=MIN/MAX rpin=... on the IBIS element) are right only into the load the V-t tables were
    made with (R_fixture): into other loads the high level is off, e.g. an Efinix LPDDR4 driver +46 mV into 60 ohm
    and +129 mV into an open line at MAX -> corner_as_typ writes a file whose typ column is the min or max column;
    simulate that at typ (exact at every load). Check with dc_high (the level from the model's own V-I tables).
"""
import re, sys
import numpy as np


def fix_pin_header(t):
    return re.sub(r'^(\[Pin\].*)$', lambda m: re.sub(r'(?i)\b([RLC])_pin\b', r'\1_pin', m.group(1)), t, flags=re.M)


def drop_series(t):
    """Remove [Series Pin Mapping] and every [Model] whose Model_type is Series / Series_switch."""
    t = re.sub(r'\[Series Pin Mapping\].*?\n(?=\[)', '', t, flags=re.S)
    def keep(m): return '' if re.search(r'^Model_type\s+Series', m.group(0), re.M | re.I) else m.group(0)
    return re.sub(r'\[Model\].*?(?=\n\[Model\]|\n\[End\]|\n\[Define Package Model\])', keep, t, flags=re.S)


def extract_models(path, names, strip_ami=True):
    """[Model] blocks by name from a (large) IBIS file, read line by line -> {name: text}."""
    out, cur = {}, None
    with open(path, errors='replace') as f:
        for line in f:
            if line.startswith('[Model]'):
                n = line.split()[1]
                cur = n if n in names else None
                if cur: out[cur] = []
            elif line.startswith(('[End]', '[Define Package Model]', '[Model Selector]')):
                cur = None
            if cur: out[cur].append(line)
    missing = set(names) - set(out)
    if missing: raise KeyError(f'models not found in {path}: {missing}')
    return {n: re.sub(r'\[Algorithmic Model\].*?\[End Algorithmic Model\]\n', '', ''.join(v), flags=re.S)
            if strip_ami else ''.join(v) for n, v in out.items()}


def add_models(t, blocks):
    """Insert [Model] blocks (text) before the final [End]."""
    i = t.rfind('[End]')
    return t[:i] + ''.join(blocks) + t[i:]


def _matrix_diag(pkg, name):
    j = pkg.find(f'[{name} Matrix]')
    if j < 0: return {}
    body = pkg[j:]
    m = re.search(r'\n\[(?!Row\]|Bandwidth\])', body[5:])
    body = body[:m.start() + 5] if m else body
    sparse = 'sparse' in body.splitlines()[0].lower()      # [L Matrix] Sparse_matrix | Banded_matrix | Full_matrix
    d = {}
    for r in re.split(r'\[Row\]\s+', body)[1:]:
        lines = [l.split('|')[0].split() for l in r.strip().splitlines()]
        lines = [l for l in lines if l]
        pin = lines[0][0]
        if sparse:                                          # rows of 'pin value'
            for l in lines[1:]:
                if len(l) >= 2 and l[0] == pin: d[pin] = num(l[1])
        elif len(lines) > 1:                                # full / banded: the row starts at the diagonal
            d[pin] = num(lines[1][0])
    return d


def package_to_pins(t):
    """Write the [Define Package Model] diagonal R/L/C into the [Pin] R_pin/L_pin/C_pin columns, for pins that
    have only pin/signal/model columns. Pins missing from the matrix get NA (the [Package] typ values)."""
    i = t.find('[Define Package Model]')
    if i < 0: return t
    pkg = t[i:]
    L, C, R = (_matrix_diag(pkg, n) for n in ('Inductance', 'Capacitance', 'Resistance'))
    out, in_pin = [], False
    for line in t.splitlines(keepends=True):
        if line.startswith('[Pin]'):
            in_pin = True
            line = '[Pin] Signal_Name   Model_Name  R_pin  L_pin C_pin\n'
        elif line.startswith('['):
            in_pin = False
        elif in_pin and line.strip() and not line.startswith('|'):
            f = line.split()
            if len(f) == 3:
                p = f[0]
                rlc = (f'{R.get(p, 0):.4f}  {L[p] * 1e9:.4f}nH  {C[p] * 1e12:.4f}pF' if p in L and p in C
                       else 'NA  NA  NA')
                line = f'{f[0]:8s}{f[1]:10s}{f[2]:12s}{rlc}\n'
        out.append(line)
    return ''.join(out)


def fix(t):
    return package_to_pins(drop_series(fix_pin_header(t)))


# ------------------------------------------------------------------------------------------------ V-t tables
_UNIT = {'p': 1e-12, 'n': 1e-9, 'u': 1e-6, 'm': 1e-3, 'f': 1e-15}


def num(s):
    """IBIS number with optional scale suffix: '1.2n' -> 1.2e-9, '0.55V' -> 0.55."""
    return float(re.match(r'[-+.\deE]+', s).group()) * _UNIT.get(re.sub(r'[-+.\deE]+', '', s)[:1].lower(), 1)


def waveform(path, model, edge, vfix=0.0, col=1):
    """[Rising|Falling Waveform] of a model for the table with V_fixture == vfix; col 1/2/3 = typ/min/max.
    Returns (t, v) at the die, without package (as IBIS tables are)."""
    t = open(path, errors='replace').read()
    i = re.search(rf'^\[Model\]\s+{re.escape(model)}\s*$', t, re.M).start()
    j = t.find('[Model]', i + 8); blk = t[i:j if j > 0 else len(t)]
    for w in re.finditer(rf'\[{edge} Waveform\](.*?)(?=\n\[)', blk, re.S):
        body = w.group(1)
        vf = re.search(r'^V_fixture\s*=\s*(\S+)', body, re.M)
        if not vf or abs(num(vf.group(1).rstrip('Vv')) - vfix) > 1e-6: continue
        rows = [l.split() for l in body.splitlines() if l.strip() and not l.strip().startswith('|') and '=' not in l]
        return np.array([num(r[0].rstrip('sS')) for r in rows]), np.array([num(r[col].rstrip('Vv')) for r in rows])
    raise KeyError(f'{model} {edge} V_fixture={vfix}')


def t2080(t, v):
    """(start time, 20-80 % edge time) of a monotonic-ish edge."""
    a, b = v[0], v[-1]; lo, hi = a + 0.2 * (b - a), a + 0.8 * (b - a)
    x = (v - lo) / (hi - lo) if b > a else (lo - v) / (lo - hi)
    i1 = np.argmax(x >= 0); i2 = np.argmax(x >= 1)
    return t[i1], t[i2] - t[i1]



# ------------------------------------------------------------------------------------------------ corners
_ROW = re.compile(r'^(\s*(?:\[[^\]]+\]|[^\s\[|][^\s=]*)\s+)(\S+)(\s+)(\S+)(\s+)(\S+)(\s*(?:\|.*)?)$')


def corner_as_typ(t, corner):
    """The file text with the min or max column copied into the typ column, inside every [Model] and [Package]:
    each line of a key plus exactly three values (V-I and V-t rows, C_comp, [Voltage Range], [Pullup Reference] etc.,
    [Ramp], R_pkg/L_pkg/C_pkg; kept where the corner is NA) and V_fixture := V_fixture_<corner>. [Pin],
    [Model Selector] and [Diff Pin] are not touched. Simulate the result at typ instead of vcc=<CORNER>."""
    col = {'min': 3, 'max': 5}[corner]
    out, section = [], None
    for line in t.splitlines():
        kw = re.match(r'^\s*\[([^\]]+)\]', line)
        if kw:
            name = kw.group(1).strip().lower()
            if name in ('model', 'package'): section = name
            elif name in ('pin', 'model selector', 'diff pin', 'component', 'define package model', 'end'): section = None
        m = _ROW.match(line) if section else None
        if m and not re.match(r'^\s*\[(Model|Package|Pin)\]', line, re.I):
            g = list(m.groups())
            if g[col].upper() != 'NA': g[1] = g[col].ljust(len(g[1]))
            line = ''.join(g)
        out.append(line)
    t = '\n'.join(out) + '\n'

    def vfix(block):
        c = re.search(rf'^V_fixture_{corner}\s*=\s*(\S+)', block, re.M)
        return re.sub(r'^(V_fixture\s*=\s*)\S+', rf'\g<1>{c.group(1)}', block, flags=re.M) if c else block
    return re.sub(r'(\[(?:Rising|Falling) Waveform\].*?)(?=\n\[)', lambda m: vfix(m.group(1)), t, flags=re.S)


def _model_block(path, model):
    t = open(path, errors='replace').read()
    i = re.search(rf'^\[Model\]\s+{re.escape(model)}\s*$', t, re.M).start()
    j = t.find('[Model]', i + 8)
    return t[i:j if j > 0 else len(t)]


def r_fixture(path, model):
    """R_fixture (ohm) of a model's [Rising Waveform] tables: the load its V-t tables (and KiCad's corners) fit."""
    return num(re.search(r'^R_fixture\s*=\s*(\S+)', _model_block(path, model), re.M).group(1).rstrip('Ohm').rstrip('ohm'))


def dc_high(path, model, col, r_load):
    """DC high level (V) of a driver into r_load ohm to ground (1e9: open), solved from its own V-I tables ([Pullup],
    [POWER Clamp] relative to Vcc, [GND Clamp]) at col 1/2/3 = typ/min/max: the reference for any load."""
    blk = _model_block(path, model)

    def tab(name):
        m = re.search(rf'^\[{name}\](.*?)(?=^\[)', blk, re.S | re.M)
        if not m: return None
        rows = [l.split() for l in m.group(1).splitlines() if l.strip() and not l.strip().startswith('|')]
        a = np.array([[num(x.rstrip('VvAa')) if x.upper() != 'NA' else np.nan for x in r[:4]]
                      for r in rows if len(r) >= 4 and re.match(r'[-+.\d]', r[0])])
        for c in (2, 3): a[:, c] = np.where(np.isnan(a[:, c]), a[:, 1], a[:, c])
        return a
    ref = re.search(r'^\[Pullup Reference\]\s+(.*)$', blk, re.M) or re.search(r'^\[Voltage Range\]\s+(.*)$', blk, re.M)
    vcc = num(ref.group(1).split()[col - 1].rstrip('Vv'))
    v = np.linspace(-0.5, 2 * vcc, 20001); i_in = np.zeros_like(v)        # current into the pin (IBIS sign)
    for name, rel in (('Pullup', True), ('POWER Clamp', True), ('GND Clamp', False)):
        a = tab(name)
        if a is not None: i_in += np.interp(vcc - v if rel else v, a[:, 0], a[:, col])
    return float(v[np.argmin(abs(i_in + v / r_load))])

if __name__ == '__main__':
    a = sys.argv[1:]
    if len(a) == 4 and a[0] == '--corner' and a[1] in ('min', 'max'):
        open(a[3], 'w').write(corner_as_typ(open(a[2], errors='replace').read(), a[1])); print('wrote', a[3])
    elif len(a) == 2:
        open(a[1], 'w').write(fix(open(a[0], errors='replace').read())); print('wrote', a[1])
    else:
        sys.exit(__doc__)
