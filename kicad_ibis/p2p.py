"""Point-to-point IBIS link in one command: schematic -> kicad-cli/KIBIS -> ngspice -> eye.

    python3 -m kicad_ibis.p2p --drv models/fpga.ibs:FPGA_COMP:T62:out_40 --rx models/dram.ibs:DRAM:C11:ODT_48 \\
        --line 50:20:7.0 --via 0.5n:0.4p --rate 1600 --mask rect:140m:0.22 [--corner min] [--png eye.png]

--drv / --rx  FILE:COMPONENT:PIN:MODEL (IBIS [Component] name, pin name from its [Pin] table, model name)
--line        Z0:LENGTH_MM:PS_PER_MM, repeatable (segments in series, e.g. neck-down, route, neck-down)
--via         L:C, a pi via at each end (omit for none)
--rate        MT/s; the PRBS7 runs --bits bits (default 150) plus a 2 ns tail
--corner      typ | min | max: driver/receiver supply, package R/L/C (C_comp stays typ in KiCad 10)

Writes WORK/p2p.kicad_sch (opens in Eeschema), p2p.cir, p2p.dat, then shows the eye at the receiver pin and die.
"""
import argparse, os, sys
from . import sch, spice, eye
from .show_eye import parse_mask, si


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--drv', required=True); ap.add_argument('--rx', required=True)
    ap.add_argument('--line', action='append', required=True)
    ap.add_argument('--via'); ap.add_argument('--rate', type=float, required=True)
    ap.add_argument('--mask', type=parse_mask, default=('rect', 0.1, 0.25))
    ap.add_argument('--corner', default='typ', choices=('typ', 'min', 'max'))
    ap.add_argument('--bits', type=int, default=150); ap.add_argument('--work', default='work')
    ap.add_argument('--png'); ap.add_argument('--raw-kukd', action='store_true')
    a = ap.parse_args(argv)
    spice.check()
    ui = 1e-6 / a.rate
    work = os.path.abspath(a.work); os.makedirs(work, exist_ok=True)

    def part(s):
        f, comp, pin, model = s.split(':')
        return (os.path.relpath(os.path.abspath(f), work), comp), pin, model
    (dl, dpin, dmodel), (rl, rpin, rmodel) = part(a.drv), part(a.rx)
    c = a.corner.upper()
    corner = '' if a.corner == 'typ' else f' vcc={c} rpin={c} lpin={c} cpin={c}'

    s = sch.Sch(models='.', project='p2p')
    s.text(f'.tran 1p {a.bits * ui * 1e9 + 2:.2f}n\n.options method=gear maxord=2 reltol=1e-4', 25.4, 15.24)
    s.row(f'{dmodel} -> {" + ".join(a.line)} -> {rmodel}, {a.rate:g} MT/s, {a.corner}')
    # Sim.Library is resolved relative to the schematic folder (kicad-cli runs there)
    drv = s.ibis('drv', dl, dpin, dmodel, 'tx', f'f0={sch.eng(1 / ui)} n={a.bits}' + corner)
    net = 'tx'
    if a.via:
        L, C = (si(x) for x in a.via.split(':')); s.via(net, 'v1', L, C); net = 'v1'
    for i, seg in enumerate(a.line):
        z, l, d = (float(x) for x in seg.split(':'))
        nxt = f'n{i + 1}'; s.tl(net, nxt, z, l * d * 1e-12); net = nxt
    if a.via:
        s.via(net, 'rx', L, C)
    else:
        s.r(net, 'rx', 1e-3)
    rx = s.ibis('rx', rl, rpin, rmodel, 'rx', corner.strip())
    path = os.path.join(work, 'p2p.kicad_sch'); s.write(path)
    msgs = spice.export(path, a.raw_kukd)
    if msgs: print('kicad-cli:', '; '.join(msgs))
    dn, rn = f'x{drv.lower()}.die0', f'x{rx.lower()}.die0'
    t, w = spice.ngspice(os.path.join(work, 'p2p.cir'), [dn, 'rx', rn], os.path.join(work, 'p2p.dat'))
    from .show_eye import main as show
    return show([os.path.join(work, 'p2p.dat'), '--drv', dn, '--rx', 'rx', rn, '--rate', str(a.rate),
                 '--mask', f'{a.mask[0]}:{a.mask[1]}:{a.mask[2]}'] + (['--png', a.png] if a.png else []))


if __name__ == '__main__':
    sys.exit(main())
