"""Show eye diagrams from an ngspice wrdata file in an interactive matplotlib window (or save PNGs).

    python3 -m kicad_ibis.show_eye WAVE.dat --list
    python3 -m kicad_ibis.show_eye WAVE.dat --drv xu1.die0 --rx dq_rx xu2.die0 --rate 1600 --mask rect:140m:0.22
    python3 -m kicad_ibis.show_eye WAVE.dat ... --png eye.png          # save instead of show

WAVE.dat: written by kicad_ibis.spice.ngspice() or by ngspice itself ('set wr_singlescale', 'set wr_vecnames',
'wrdata WAVE.dat v(a) v(b) ...'). --drv is the PRBS driver's die node ('x<ref>.die0', e.g. xu1.die0), used to line
the bits up in time; each --rx node gets its own eye (several share one window). --mask shape:height:width_UI,
shape rect or diamond, height with an SI suffix (140m = 140 mV). --tail: simulated time after the last PRBS bit.
Traces are coloured by the bit sent (1 red, 0 blue): a red trace inside the blue band means a wrong label or a
closed eye, worth a look either way.
"""
import argparse, math, sys
from . import eye, spice

SI = {'f': 1e-15, 'p': 1e-12, 'n': 1e-9, 'u': 1e-6, 'm': 1e-3, 'k': 1e3, 'M': 1e6, 'G': 1e9}


def si(s):
    return float(s[:-1]) * SI[s[-1]] if s[-1] in SI else float(s)


def parse_mask(s):
    shape, h, w = s.split(':')
    if shape not in ('rect', 'diamond'): raise argparse.ArgumentTypeError('mask shape is rect or diamond')
    return shape, si(h), float(w)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dat')
    ap.add_argument('--list', action='store_true', help='list the vectors in the file')
    ap.add_argument('--drv'); ap.add_argument('--rx', nargs='+')
    ap.add_argument('--rate', type=float, help='MT/s (bit rate / 1e6)')
    ap.add_argument('--mask', type=parse_mask, default=('rect', 0.1, 0.25))
    ap.add_argument('--tail', type=si, default=2e-9); ap.add_argument('--skip', type=int, default=12)
    ap.add_argument('--png'); ap.add_argument('--title', default='')
    a = ap.parse_args(argv)
    t, w = spice.load(a.dat)
    if a.list or not (a.drv and a.rx and a.rate):
        print('vectors:', ' '.join(w)); return
    ui = 1e-6 / a.rate
    import matplotlib
    if a.png: matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    n = len(a.rx); cols = min(n, 2); rows = math.ceil(n / cols)
    fig, axs = plt.subplots(rows, cols, figsize=(8 * cols, 5 * rows), squeeze=False)
    for ax, node in zip(axs.ravel(), a.rx):
        ph, tr, lab, m = eye.analyse(t, w[a.drv], w[node], ui, a.mask, a.tail, a.skip)
        eye.draw(ax, ph, tr, lab, m, a.mask, ui, f'{a.title} {node}  ({a.rate:g} MT/s, delay {m["delay_ps"]:.0f} ps)')
        print(f"{node}: EH {m['eh_mV']:.0f} mV  EW {m['ew0_ps']:.0f} ps  margin {m['v_margin_mV']:+.0f} mV / "
              f"{m['t_margin_ps']:+.0f} ps  Vmax {m['vmax_mV']:.0f}  Vmin {m['vmin_mV']:.0f} mV")
    for ax in axs.ravel()[n:]: ax.set_visible(False)
    fig.tight_layout()
    if a.png: fig.savefig(a.png, dpi=110); print('wrote', a.png)
    else: plt.show()


if __name__ == '__main__':
    sys.exit(main())
