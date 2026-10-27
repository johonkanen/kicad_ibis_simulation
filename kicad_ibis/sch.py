"""Write a KiCad 10 simulation schematic from Python: IBIS driver/receiver, TLINE, R/L/C, text and .tran directives.

Every element is a normal KiCad symbol (Simulation_SPICE / Device libraries, copied into the file), connected by
net labels on the pin ends, so the file also opens, ERCs and simulates in Eeschema. Example:

    s = Sch('models')                     # Sim.Library paths are written relative to the schematic: models/x.ibs
    s.text('.tran 1p 100n', 25.4, 15.24)
    s.row('DQ write')
    s.ibis('drv', ('fpga.ibs', 'FPGA_COMPONENT'), 'T62', 'out_model_40', 'dq_tx', 'f0=1.6G n=150')
    s.tl('dq_tx', 'dq_rx', 50, 140e-12)
    s.ibis('rx', ('dram.ibs', 'DRAM'), 'C11', 'ODT_48', 'dq_rx')
    s.write('work/dq.kicad_sch')

KIBIS gives each IBIS symbol a subcircuit whose die node is 'x<ref>.die0' (ref lower case, e.g. xu2.die0); the
pin/ball node is the net name.
"""
import os, re, uuid

# KiCad symbol libraries (same rules as spice.py; no numpy needed here)
KSYM = os.environ.get('KICAD_SYMBOLS') or os.path.join(os.path.dirname(os.environ.get('KICAD_APPRUN', '')),
                                                       'usr/share/kicad/symbols')


def uid(): return str(uuid.uuid4())


def eng(v):
    """SPICE engineering notation: 1.6e9 -> 1.6G, 4e-13 -> 400f."""
    for e, s in ((1e9, 'G'), (1e6, 'Meg'), (1, ''), (1e-3, 'm'), (1e-6, 'u'), (1e-9, 'n'), (1e-12, 'p'), (1e-15, 'f')):
        if abs(v) >= e * 0.999: return f'{v / e:.4g}{s}'
    return f'{v:.4g}'


def _symfile(lib, name):
    return open(os.path.join(KSYM, f'{lib}.kicad_symdir', f'{name}.kicad_sym')).read()


def lib_symbol(lib, name):
    t = _symfile(lib, name)
    i = t.find(f'(symbol "{name}"')
    body = t[i:t.rstrip().rfind(')')].rstrip()
    return body.replace(f'(symbol "{name}"', f'(symbol "{lib}:{name}"', 1)


def pins_of(lib, name):
    return {m[2]: (float(m[0]), float(m[1])) for m in
            re.findall(r'\(pin \w+ \w+\s*\(at ([-\d.]+) ([-\d.]+) [-\d.]+\).*?\(number "(\w+)"', _symfile(lib, name),
                       re.S)}


class Sch:
    def __init__(self, models='models', project='ibis_sim', paper='A2'):
        self.models = models; self.project = project; self.paper = paper
        self.root = uid(); self.items = []; self.libs = {}; self.n = {}; self.x = 25.4; self.y = 25.4

    def prop(self, k, v, x, y, hide=True):
        return (f'\t\t(property "{k}" "{v}"\n\t\t\t(at {x:.2f} {y:.2f} 0)\n\t\t\t(effects\n\t\t\t\t(font\n'
                f'\t\t\t\t\t(size 1.27 1.27)\n\t\t\t\t)\n{"				(hide yes)" + chr(10) if hide else ""}\t\t\t)\n\t\t)\n')

    def sym(self, lib, name, pfx, nets, value='', props=None):
        """Place lib:name with reference prefix pfx; nets = {pin number: net}. Returns the reference."""
        key = f'{lib}:{name}'
        if key not in self.libs: self.libs[key] = (lib_symbol(lib, name), pins_of(lib, name))
        pins = self.libs[key][1]
        self.n[pfx] = self.n.get(pfx, 0) + 1; ref = f'{pfx}{self.n[pfx]}'
        x, y = self.x, self.y
        s = (f'\t(symbol\n\t\t(lib_id "{key}")\n\t\t(at {x:.2f} {y:.2f} 0)\n\t\t(unit 1)\n\t\t(exclude_from_sim no)\n'
             f'\t\t(in_bom yes)\n\t\t(on_board yes)\n\t\t(dnp no)\n\t\t(uuid "{uid()}")\n')
        s += self.prop('Reference', ref, x, y - 6.35, False) + self.prop('Value', value or name, x, y + 7.62, False)
        s += self.prop('Footprint', '', x, y) + self.prop('Datasheet', '', x, y)
        for k, v in (props or {}).items(): s += self.prop(k, v, x, y)
        for p in pins: s += f'\t\t(pin "{p}"\n\t\t\t(uuid "{uid()}")\n\t\t)\n'
        s += (f'\t\t(instances\n\t\t\t(project "{self.project}"\n\t\t\t\t(path "/{self.root}"\n'
              f'\t\t\t\t\t(reference "{ref}")\n\t\t\t\t\t(unit 1)\n\t\t\t\t)\n\t\t\t)\n\t\t)\n\t)\n')
        self.items.append(s)
        for p, net in nets.items():   # a local label on each pin end (no wires needed)
            px, py = pins[p]
            self.items.append(f'\t(label "{net}"\n\t\t(at {x + px:.2f} {y - py:.2f} 0)\n\t\t(effects\n\t\t\t(font\n'
                              f'\t\t\t\t(size 1.27 1.27)\n\t\t\t)\n\t\t\t(justify left bottom)\n\t\t)\n\t\t(uuid "{uid()}")\n\t)\n')
        self.x += 30.48
        return ref

    def text(self, s, x, y):
        """Sheet text; lines starting with '.' (.tran, .options, ...) are simulation directives."""
        s = s.replace('"', '\\"').replace('\n', '\\n')
        self.items.append(f'\t(text "{s}"\n\t\t(exclude_from_sim no)\n\t\t(at {x:.2f} {y:.2f} 0)\n\t\t(effects\n'
                          f'\t\t\t(font\n\t\t\t\t(size 1.27 1.27)\n\t\t\t)\n\t\t\t(justify left bottom)\n\t\t)\n'
                          f'\t\t(uuid "{uid()}")\n\t)\n')

    def row(self, title):
        """Start a new row of symbols (one channel per row) with a title."""
        self.y += 30.48; self.x = 25.4
        self.text(title, 25.4, self.y - 12.7)

    # elements
    def ibis(self, kind, lib, pin, model, net, params='', stype=None):
        """IBIS buffer. kind 'drv' (PRBSDRIVER unless stype, e.g. RECTDRIVER) or 'rx' (IBISDEVICE: a receiver or
        an ODT/terminator model). lib = (file under self.models, IBIS [Component] name). params e.g.
        'f0=1.6G n=150' (PRBS: f0 is the BIT RATE) plus corner overrides ' vcc=MAX rpin=MAX lpin=MAX cpin=MAX'
        (upper case: lower case silently gives typ)."""
        name = 'IBIS_DRIVER' if kind == 'drv' else 'IBIS_DEVICE'
        st = stype or ('PRBSDRIVER' if kind == 'drv' else 'IBISDEVICE')
        props = {'Sim.Device': 'IBIS', 'Sim.Type': st, 'Sim.Pins': '1=GND 2=IN/OUT',
                 'Sim.Library': f'{self.models}/{lib[0]}', 'Sim.Name': lib[1], 'Sim.Ibis.Pin': pin,
                 'Sim.Ibis.Model': model}
        if params: props['Sim.Params'] = params
        return self.sym('Simulation_SPICE', name, 'U', {'1': 'GND', '2': net}, f'{model}', props)

    def tl(self, a, b, z0, td):
        """Lossless line (SPICE T element), z0 ohm, delay td s."""
        return self.sym('Simulation_SPICE', 'TLINE', 'TL', {'1': a, '2': 'GND', '3': b, '4': 'GND'},
                        f'TLINE {z0:g}R {td * 1e12:.1f}ps',
                        {'Sim.Device': 'TLINE', 'Sim.Pins': '1=1+ 2=1- 3=2+ 4=2-',
                         'Sim.Params': f'z0={z0:g} td={eng(td)}'})

    def r(self, a, b, v): return self.sym('Device', 'R', 'R', {'1': a, '2': b}, eng(v))

    def c(self, a, b, v): return self.sym('Device', 'C', 'C', {'1': a, '2': b}, eng(v))

    def l(self, a, b, v): return self.sym('Device', 'L', 'L', {'1': a, '2': b}, eng(v))

    def via(self, a, b, L, C):
        """Via as a pi section: C/2, L, C/2."""
        self.c(a, 'GND', C / 2); self.l(a, b, L); self.c(b, 'GND', C / 2)

    def render(self):
        libs = ''.join('\t\t' + v[0].replace('\n', '\n\t') + '\n' for v in self.libs.values())
        return (f'(kicad_sch\n\t(version 20250114)\n\t(generator "eeschema")\n\t(generator_version "10.0")\n'
                f'\t(uuid "{self.root}")\n\t(paper "{self.paper}")\n\t(lib_symbols\n{libs}\t)\n' + ''.join(self.items) +
                '\t(sheet_instances\n\t\t(path "/"\n\t\t\t(page "1")\n\t\t)\n\t)\n\t(embedded_fonts no)\n)\n')

    def write(self, path):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        open(path, 'w').write(self.render())
