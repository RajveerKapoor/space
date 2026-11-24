"""Pixel-level vetting: difference images & PSF-fit source location / true depth.
usage: pixvet.py TIC P T0 DUR_H  [neighborRA,neighborDec ...]"""
import sys, glob, numpy as np, warnings
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
from scipy.ndimage import median_filter
S = '../work'
tic = int(sys.argv[1]); P = float(sys.argv[2]); T0 = float(sys.argv[3]); DUR = float(sys.argv[4]) / 24
others = [tuple(map(float, a.split(','))) for a in sys.argv[5:]]
D = f'{S}/cand/{tic}'
tps = sorted(glob.glob(f'{D}/*_tp.fits'))
fig, axes = plt.subplots(3, len(tps), figsize=(4 * len(tps), 11), squeeze=False)

def gauss(params, yy, xx, n):
    # params: [x1,y1,A1, x2,y2,A2 ..., sigma, bkg]
    sig = params[-2]; out = np.full(xx.shape, params[-1])
    for k in range(n):
        x, y, A = params[3 * k:3 * k + 3]
        out += A * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * sig ** 2))
    return out

summary = []
for j, fn in enumerate(tps):
    h = fits.open(fn)
    d = h[1].data; t = d['TIME']; q = d['QUALITY']; cube = d['FLUX'].astype(float)
    sec = h[0].header['SECTOR']
    w = WCS(h[2].header)
    ok = np.isfinite(t) & (q == 0) & np.all(np.isfinite(cube.reshape(len(t), -1)), axis=1)
    t = t[ok]; cube = cube[ok]
    ra0, de0 = h[0].header['RA_OBJ'], h[0].header['DEC_OBJ']
    px_t = w.world_to_pixel_values(ra0, de0)
    px_o = [w.world_to_pixel_values(r, dd) for r, dd in others]
    ap = h[2].data
    ph = ((t - T0) / P + 0.5) % 1 - 0.5
    dt = ph * P
    epochs = np.round((t - T0) / P)
    intr = np.abs(dt) < DUR * 0.35
    oot = (np.abs(dt) > DUR * 0.75) & (np.abs(dt) < DUR * 0.75 + 0.35)
    if intr.sum() < 10:
        print(f'S{sec}: no in-transit data'); continue
    # local per-epoch differences to suppress drifts
    diffs, wts = [], []
    for e in np.unique(epochs[intr]):
        a = intr & (epochs == e); b = oot & (epochs == e)
        if a.sum() > 5 and b.sum() > 10:
            diffs.append(np.nanmean(cube[b], 0) - np.nanmean(cube[a], 0)); wts.append(a.sum())
    diff = np.average(diffs, axis=0, weights=wts)
    ootimg = np.nanmedian(cube[oot], 0)
    ny, nx = diff.shape; yy, xx = np.mgrid[0:ny, 0:nx]
    # 1) single-source fit to the difference image
    i0 = np.unravel_index(np.argmax(median_filter(diff, 2)), diff.shape)
    p0 = [i0[1], i0[0], diff.max(), 0.8, 0]
    r1 = least_squares(lambda p: (gauss(p, yy, xx, 1) - diff).ravel(), p0)
    xs, ys = r1.x[0], r1.x[1]
    # 2) two-star fit of OOT image at fixed WCS positions (target + first neighbor) + free sigma
    allpos = [px_t] + px_o
    def model_fixed(p, pos):
        sig, bkg = p[-2], p[-1]; out = np.full(xx.shape, bkg)
        for k, (x, y) in enumerate(pos):
            out += p[k] * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * sig ** 2))
        return out
    pp0 = [ootimg.max()] + [ootimg.max() / 10] * len(px_o) + [0.8, np.median(ootimg)]
    ro = least_squares(lambda p: (model_fixed(p, allpos) - ootimg).ravel(), pp0)
    sig_o = ro.x[-2]
    # 3) fit diff image with fixed positions & fixed sigma -> deficit amplitude per star
    def model_d(p):
        out = np.full(xx.shape, p[-1])
        for k, (x, y) in enumerate(allpos):
            out += p[k] * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * sig_o ** 2))
        return out
    rd = least_squares(lambda p: (model_d(p) - diff).ravel(), [0] * len(allpos) + [0])
    res = dict(sector=sec, src_x=xs, src_y=ys, tgt=px_t, nb=px_o,
               off_tgt=np.hypot(xs - px_t[0], ys - px_t[1]) * 21,
               off_nb=[np.hypot(xs - x, ys - y) * 21 for x, y in px_o],
               amp_oot=ro.x[:len(allpos)], amp_diff=rd.x[:len(allpos)], sig=sig_o)
    res['depth_each'] = res['amp_diff'] / res['amp_oot']
    summary.append(res)
    print(f"S{sec}: diff-source at ({xs:.2f},{ys:.2f}); target ({px_t[0]:.2f},{px_t[1]:.2f}) off={res['off_tgt']:.1f}\"; "
          f"neighbors off={[round(o,1) for o in res['off_nb']]}\"  PSFsig={sig_o:.2f}px")
    print(f"     OOT amps {np.round(ro.x[:len(allpos)],1)}  diff amps {np.round(rd.x[:len(allpos)],2)}  => depth per star {np.round(res['depth_each']*100,3)} %")
    for row, img, ttl in [(0, ootimg, f'S{sec} OOT'), (1, diff, 'OOT - IN (difference)'), (2, diff / np.sqrt(np.abs(ootimg) + 1), 'diff / sqrt(flux)')]:
        a = axes[row, j]; a.imshow(img, origin='lower', cmap='viridis')
        a.plot(*px_t, 'r*', ms=14, label='target'); [a.plot(x, y, 'wx', ms=10, mew=2) for x, y in px_o]
        if row == 1: a.plot(xs, ys, 'o', mfc='none', mec='magenta', ms=16, mew=2)
        a.contour(ap & 2 > 0, levels=[0.5], colors='orange', linewidths=1)
        a.set_title(ttl, fontsize=9)
plt.tight_layout(); plt.savefig(f'{D}/pixvet.png', dpi=70)
import pickle; pickle.dump(summary, open(f'{D}/pixvet.pkl', 'wb'))
