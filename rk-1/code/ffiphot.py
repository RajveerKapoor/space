"""Simple aperture photometry on TESScut cutouts + single-event box search.
usage: ffiphot.py DIR TC_KNOWN DEPTH DUR_H out.png"""
import sys, glob, numpy as np, warnings
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from astropy.io import fits
from wotan import flatten
D = sys.argv[1]; tck = float(sys.argv[2]); dep = float(sys.argv[3]); dur = float(sys.argv[4]) / 24; out = sys.argv[5]
fs = sorted(glob.glob(f'{D}/s*/*.fits'), key=lambda f: int(f.split('tess-s')[1][:4]))
fig, axes = plt.subplots(len(fs), 1, figsize=(15, 2.6 * len(fs)))
allres = {}
for ax, f in zip(axes, fs):
    h = fits.open(f); d = h[1].data; sec = h[0].header['SECTOR']
    t = d['TIME']; cube = d['FLUX'].astype(float); q = d['QUALITY']
    ny0, nx0 = cube.shape[1:]; c0 = cube[:, ny0//2-3:ny0//2+4, nx0//2-3:nx0//2+4]
    ok = np.isfinite(t) & (q == 0) & np.all(np.isfinite(c0.reshape(len(t), -1)), axis=1)
    t, cube = t[ok], cube[ok]
    med = np.nanmedian(cube, 0); med = np.where(np.isfinite(med), med, np.nanmax(med))
    ny, nx = med.shape; cy, cx = ny // 2, nx // 2
    yy, xx = np.mgrid[0:ny, 0:nx]
    r = np.hypot(yy - cy, xx - cx)
    bkgmask = (med < np.percentile(med, 40)) & (r > 4)
    bkg = np.nanmedian(cube[:, bkgmask], axis=1)
    sub = cube - bkg[:, None, None]
    # aperture: brightest connected pixels near centre
    msub = np.nanmedian(sub, 0)
    ap = (r <= 2.0) & (msub > 0.25 * msub[cy, cx])
    lc = sub[:, ap].sum(1)
    lc = lc / np.nanmedian(lc)
    cad = np.nanmedian(np.diff(t)) * 1440
    fl = flatten(t, lc, method='biweight', window_length=1.0, break_tolerance=0.4)
    okf = np.isfinite(fl)
    dtt = np.diff(t); gi = np.where(dtt > 0.4)[0]
    edges = np.concatenate([[t[0]], t[gi], t[gi + 1], [t[-1]]])
    okf &= np.min(np.abs(t[:, None] - edges[None, :]), axis=1) > 0.25
    # box search with width = dur
    sig = 1.4826 * np.nanmedian(np.abs(fl[okf] - 1))
    best = []
    for tc in np.arange(t.min(), t.max(), 0.01):
        m = okf & (np.abs(t - tc) < dur / 2)
        if m.sum() < max(3, 0.6 * dur * 1440 / cad): continue
        dd = 1 - fl[m].mean(); s = dd / (sig / np.sqrt(m.sum()))
        best.append((s, tc, dd))
    best.sort(reverse=True)
    top = []
    for b in best:
        if all(abs(b[1] - x[1]) > 0.5 for x in top): top.append(b)
        if len(top) >= 4: break
    allres[sec] = dict(t=t, fl=fl, top=top, cad=cad, sig=sig, span=(t.min(), t.max()))
    print(f'S{sec}: cadence {cad:.1f} min, {len(t)} pts, {t.min():.2f}-{t.max():.2f}, sigma={sig*1e3:.2f} ppt, ap npix={ap.sum()}; '
          f'top events: ' + ', '.join(f'(snr {s:.1f} @ {tc:.2f} depth {dd*100:.2f}%)' for s, tc, dd in top))
    ax.plot(t, fl, '.', ms=1.5, color='k')
    ax.set_ylim(1 - 2.2 * dep, 1 + 1.2 * dep)
    for s, tc, dd in top[:2]:
        if s > 6: ax.axvline(tc, color='r', alpha=0.4)
    ax.axvline(tck, color='g', alpha=0.5) if t.min() < tck < t.max() else None
    ax.set_title(f'S{sec} ({cad:.0f}-min FFI cutout photometry)', fontsize=9)
plt.tight_layout(); plt.savefig(out, dpi=65)
import pickle; pickle.dump(allres, open(out.replace('.png', '.pkl'), 'wb'))
