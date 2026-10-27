"""2D finite-difference field solve of a single trace: Z0 and delay (ps/mm) vs width. C only (lossless).

    python3 -m kicad_ibis.zsolve strip  --below 0.357 --above 0.1 --t 0.0152 --er-below 4.16 --er-above 4.6 -w 0.1 0.12
    python3 -m kicad_ibis.zsolve micro  --h 0.0994 --t 0.035 --er 4.1 --mask 0.03 --er-mask 3.5 -w 0.1 0.15

All dimensions in mm. Stripline: the trace sits on top of the 'below' dielectric (prepreg, er_below) and the
'above' dielectric (core, er_above) reaches the upper plane; the trace thickness is embedded in er_below.
Microstrip: trace on h of er over the plane, covered by a mask of thickness mask (also filling around the trace).
Checked against a commercial-style solve of a JLC 6-layer stack-up within ~2 %.
"""
import argparse
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

C0 = 299792458.0
EPS0 = 8.8541878128e-12


def solve(layers, w, t_cu, z_tr, dx=5e-6, dz=2.5e-6, span=3e-3):
    """layers: (thickness, er) from the lower plane (z=0) up to the upper plane (top boundary is a plane: model
    air as a tall er=1 layer). Trace: width w, thickness t_cu, bottom at z_tr. SI units. Returns (Z0, eps_eff)."""
    H = sum(t for t, _ in layers)
    nx, nz = int(span / dx) + 1, int(round(H / dz)) + 1
    x = (np.arange(nx) - (nx - 1) / 2) * dx; z = np.arange(nz) * dz
    zc = (z[:-1] + z[1:]) / 2
    er_c = np.empty(nz - 1); zt = 0
    for t, er in layers:
        er_c[(zc >= zt) & (zc < zt + t + 1e-12)] = er; zt += t
    cond = (np.abs(x)[:, None] <= w / 2 + 1e-12) & (z[None, :] >= z_tr - 1e-12) & (z[None, :] <= z_tr + t_cu + 1e-12)

    def cap(er_cells):
        E = np.repeat(er_cells[None, :], nx - 1, axis=0)
        Eh = np.zeros((nx - 1, nz)); Eh[:, 1:] += E / 2; Eh[:, :-1] += E / 2
        Eh[:, 0] *= 2; Eh[:, -1] *= 2
        Ev = np.zeros((nx, nz - 1)); Ev[1:, :] += E / 2; Ev[:-1, :] += E / 2
        Ev[0, :] *= 2; Ev[-1, :] *= 2
        gh = Eh * dz / dx; gv = Ev * dx / dz
        idx = np.arange(nx * nz).reshape(nx, nz)
        fixed = np.zeros((nx, nz), bool); fixed[:, 0] = True; fixed[:, -1] = True; fixed |= cond
        V = np.where(cond, 1.0, 0.0)
        rows, cols, vals = [], [], []
        diag = np.zeros((nx, nz))
        for (a, b, g) in ((idx[:-1, :], idx[1:, :], gh), (idx[:, :-1], idx[:, 1:], gv)):
            rows += [a.ravel(), b.ravel()]; cols += [b.ravel(), a.ravel()]; vals += [-g.ravel(), -g.ravel()]
            np.add.at(diag.ravel(), a.ravel(), g.ravel()); np.add.at(diag.ravel(), b.ravel(), g.ravel())
        A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(nx * nz,) * 2)
        A = A + sp.diags(diag.ravel())
        free = ~fixed.ravel()
        rhs = -A[free][:, ~free] @ V.ravel()[~free]
        Vf = V.ravel().copy(); Vf[free] = spla.spsolve(A[free][:, free].tocsc(), rhs)
        return (A @ Vf)[cond.ravel()].sum() * EPS0
    C = cap(er_c); Ca = cap(np.ones_like(er_c))
    return 1 / (C0 * np.sqrt(C * Ca)), C / Ca


def stripline(w, below, above, t, er_below, er_above):
    """Asymmetric stripline, mm in; returns (Z0 ohm, ps/mm)."""
    m = 1e-3
    z, ee = solve([(below * m, er_below), (t * m, er_below), (above * m, er_above)], w * m, t * m, below * m)
    return z, np.sqrt(ee) / C0 * 1e9


def microstrip(w, h, t, er, mask=0.0, er_mask=3.5, air=1.0):
    """Microstrip with optional solder mask, mm in; returns (Z0 ohm, ps/mm)."""
    m = 1e-3
    layers = [(h * m, er), (t * m, er_mask if mask else 1.0)] + ([(mask * m, er_mask)] if mask else []) + \
             [(air * m, 1.0)]
    z, ee = solve(layers, w * m, t * m, h * m)
    return z, np.sqrt(ee) / C0 * 1e9


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('kind', choices=('strip', 'micro'))
    ap.add_argument('-w', type=float, nargs='+', required=True, help='trace widths, mm')
    ap.add_argument('--below', type=float); ap.add_argument('--above', type=float)
    ap.add_argument('--er-below', type=float); ap.add_argument('--er-above', type=float)
    ap.add_argument('--h', type=float); ap.add_argument('--er', type=float)
    ap.add_argument('--t', type=float, default=0.035, help='copper thickness, mm')
    ap.add_argument('--mask', type=float, default=0.0); ap.add_argument('--er-mask', type=float, default=3.5)
    a = ap.parse_args()
    for w in a.w:
        z, d = (stripline(w, a.below, a.above, a.t, a.er_below, a.er_above) if a.kind == 'strip'
                else microstrip(w, a.h, a.t, a.er, a.mask, a.er_mask))
        print(f'w {w:.3f} mm  Z0 {z:5.1f} ohm  {d:.2f} ps/mm')
