"""Eye diagrams from a KIBIS PRBS simulation: bit labels, folding, Rx-mask margins, plots.

The waveform at the receiver is cut into 1-UI traces, each labelled with the bit that was sent. The labels come from
the KIBIS stimulus itself (a fixed PRBS7, see kibis_prbs7), lined up in time with the driver's own die waveform; the
receiver delay comes from a cross-correlation. Masks are (shape, height V, width in UI), shape 'rect' or 'diamond'.
"""
import numpy as np


def prbs7(n, taps=(7, 6), rot=0):
    """Fibonacci LFSR output, register initialised to all ones; n bits starting rot bits into the sequence."""
    r = [1] * 7; out = []
    for _ in range(127):
        b = r[taps[0] - 1] ^ r[taps[1] - 1]; out.append(r[-1]); r = [b] + r[:-1]
    q = np.array(out)
    return np.tile(q, (rot + n) // 127 + 2)[rot:rot + n]


def kibis_prbs7(n, shift=0):
    """The bits KiCad 10's PRBSDRIVER sends: x^7+x^6+1, 72 bits into prbs7(). Found by fitting every rotation and
    polarity of both PRBS7 polynomials to the driver waveforms of ~60 runs: always this one (100 % match on sane
    lines). Bit k of the returned array is the k-th full UI after the driver's first edge-grid point (see
    stimulus_bits); shift moves the alignment by whole bits."""
    return prbs7(n, rot=(72 + shift) % 127)


def edge_phase(t, vd, ui):
    """Phase (0..ui) of the bit grid from the driver's mid-level crossings, and the mid level."""
    lo, hi = np.percentile(vd, [2, 98]); mid = (lo + hi) / 2
    s = np.sign(vd - mid); i = np.nonzero(s[1:] != s[:-1])[0]
    tc = t[i] + (mid - vd[i]) * (t[i + 1] - t[i]) / (vd[i + 1] - vd[i])
    return np.angle(np.mean(np.exp(2j * np.pi * tc / ui))) / (2 * np.pi) * ui % ui, mid


def stimulus_bits(t, vd, ui, tail=2e-9, min_match=0.7):
    """Bit grid and bit values for a KIBIS PRBS driver, from its die waveform vd: edges at phi + k*ui, bit k sent
    in [phi + k*ui, phi + (k+1)*ui]. The values are kibis_prbs7(), aligned by the whole-bit shift that best matches
    the sampled driver levels. Do not use the sampled levels themselves: on a badly mismatched line the reflection
    back to the driver die flips isolated bits. tail = simulated time after the last stimulus bit (.tran length
    minus bits*UI). Raises if nothing matches (not a KIBIS PRBS7 stimulus, or a broken driver)."""
    phi, mid = edge_phase(t, vd, ui)
    n = int((t[-1] - tail - phi) / ui)
    raw = (np.interp(phi + (np.arange(n) + 0.5) * ui, t, vd) > mid).astype(int)
    score, bits = max(((np.mean(b == raw), b) for b in (kibis_prbs7(n, s) for s in range(-3, 4))),
                      key=lambda x: x[0])
    if score < min_match: raise ValueError(f'driver does not match the KIBIS PRBS7 ({score:.2f})')
    return phi, bits


def fold(t, v, phi, bits, ui, skip=12, dt=1e-12, max_delay=4e-9):
    """Cut the receiver waveform v into 1-UI traces aligned to the sent bits. The delay is the peak of the
    cross-correlation with the ideal bit stream (0..max_delay), refined by the whole-UI shift (+-2) that gives the
    largest eye opening (a slow, mismatched line broadens the peak). skip = start-up bits left out.
    Returns (phase axis, traces [n, samples], labels [n], delay)."""
    spb = int(round(ui / dt)); n = len(bits)
    tt = phi + np.arange(0, n * spb) * dt; vv = np.interp(tt, t, v)
    ideal = np.repeat(bits * 2 - 1.0, spb)
    s0 = skip * spb; L = int(max_delay / dt)
    a = vv[s0:] - vv[s0:].mean(); b = ideal[s0:]
    nf = 1 << int(np.ceil(np.log2(len(a) + len(b))))
    xc = np.fft.irfft(np.fft.rfft(a, nf) * np.conj(np.fft.rfft(b, nf)), nf)[:L]   # xc[lag] = sum a[i+lag] b[i]
    ph = np.arange(spb) * dt

    def cut(lag):
        traces, labels = [], []
        for k in range(skip, n):
            t0 = phi + lag + k * ui
            if t0 + ui > t[-1]: break
            traces.append(np.interp(t0 + ph, t, v)); labels.append(bits[k])
        return np.array(traces), np.array(labels)
    best = None
    for m in range(-2, 3):
        lag = int(np.argmax(xc)) * dt + m * ui
        if lag < 0: continue
        tr, lab = cut(lag); op = (tr[lab == 1].min(0) - tr[lab == 0].max(0)).max()
        if best is None or op > best[0]: best = (op, lag, tr, lab)
    return ph, best[2], best[3], best[1]


def metrics(ph, tr, lab, mask, ui):
    """Eye height/width and margins against the Rx mask. v_margin = opening inside the mask window minus the mask
    height (mask at its best position); t_margin = eye width at the mask height minus the mask width (best Vref).
    nan t_margin = no opening at the mask height."""
    shape, V, W = mask; dt = ph[1] - ph[0]; w = int(round(W * ui / dt))
    top, bot = tr[lab == 1].min(0), tr[lab == 0].max(0)
    op = top - bot
    eh = op.max()

    def widths(vrs, h):   # longest run of 'open' phases for each Vref
        ok = (top[None, :] >= vrs[:, None] + h / 2) & (bot[None, :] <= vrs[:, None] - h / 2)
        d = np.diff(np.pad(ok.astype(np.int8), ((0, 0), (1, 1))), axis=1)
        out = np.zeros(len(vrs))
        for i in range(len(vrs)):
            s, e = np.nonzero(d[i] == 1)[0], np.nonzero(d[i] == -1)[0]
            if len(s): out[i] = (e - s).max()
        return out * dt
    vr_scan = np.linspace(bot.min(), top.max(), 400)
    ew0 = widths(vr_scan, 0).max()
    wm = widths(vr_scan, V); ewm, vr_m = wm.max(), vr_scan[int(wm.argmax())]
    vm = -1e9
    for s in range(0, len(ph) - w):
        if shape == 'rect':
            o = top[s:s + w].min() - bot[s:s + w].max() - V
        else:
            x = np.abs(np.arange(w) - (w - 1) / 2) / ((w - 1) / 2)
            h = V / 2 * (1 - x)
            o = (top[s:s + w] - h).min() - (bot[s:s + w] + h).max()
        vm = max(vm, o)
    return dict(eh_mV=eh * 1e3, ew0_ps=ew0 * 1e12, vref_mV=vr_m * 1e3, ew_mask_ps=ewm * 1e12,
                t_margin_ps=(ewm - W * ui) * 1e12 if ewm > 0 else float('nan'), v_margin_mV=vm * 1e3,
                top=top, bot=bot)


def analyse(t, v_drv, v_rx, ui, mask, tail=2e-9, skip=12):
    """One receiver: labels, fold, metrics, plus the over/undershoot after start-up. Returns (ph, tr, lab, m)."""
    phi, bits = stimulus_bits(t, v_drv, ui, tail)
    ph, tr, lab, lag = fold(t, v_rx, phi, bits, ui, skip)
    m = metrics(ph, tr, lab, mask, ui)
    seg = v_rx[t > skip * ui]
    m.update(delay_ps=lag * 1e12, vmax_mV=seg.max() * 1e3, vmin_mV=seg.min() * 1e3)
    return ph, tr, lab, m


def draw(ax, ph, tr, lab, m, mask, ui, title='', colour_bits=True):
    """Two-UI eye on a matplotlib axis, with the mask at its widest point. 1 bits red, 0 bits blue."""
    import matplotlib.pyplot as plt
    x = np.concatenate([ph, ph + ui]) * 1e12
    for i in range(len(tr)):
        nxt = tr[i + 1] if i + 1 < len(tr) else tr[i]
        c = ('#d62728' if lab[i] else '#1f5fa8') if colour_bits else '#1f5fa8'
        ax.plot(x, np.concatenate([tr[i], nxt]) * 1e3, color=c, lw=0.5, alpha=0.3)
    shape, V, W = mask
    c = ph[int(np.argmax(m['top'] - m['bot']))] * 1e12; vr = m['vref_mV']; w = W * ui * 1e12; h = V * 1e3
    if shape == 'rect':
        ax.add_patch(plt.Rectangle((c - w / 2, vr - h / 2), w, h, fill=False, ec='k', lw=1.4))
    else:
        ax.add_patch(plt.Polygon([(c - w / 2, vr), (c, vr + h / 2), (c + w / 2, vr), (c, vr - h / 2)], fill=False,
                                 ec='k', lw=1.4))
    ax.set_xlim(0, 2 * ui * 1e12); ax.set_xlabel('time (ps), 2 UI'); ax.set_ylabel('mV'); ax.grid(alpha=0.3)
    ax.set_title(title, fontsize=9)
    txt = (f"EH {m['eh_mV']:.0f} mV  EW {m['ew0_ps']:.0f} ps  mask {shape} {V * 1e3:.0f} mV x {W:.2f} UI: "
           f"margin {m['v_margin_mV']:+.0f} mV / {m['t_margin_ps']:+.0f} ps")
    if 'vmax_mV' in m: txt += f"  Vmax {m['vmax_mV']:.0f}  Vmin {m['vmin_mV']:.0f} mV"
    ax.text(0.01, 0.02, txt, transform=ax.transAxes, fontsize=7.5)


def save(png, *a, **k):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=110)
    draw(ax, *a, colour_bits=False, **k); fig.tight_layout(); fig.savefig(png); plt.close(fig)
