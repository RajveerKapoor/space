"""Automated transit-source localization.
usage: locate.py TIC DUR_H tc1 [tc2 ...]
Downloads TPFs for sectors containing the events, queries TIC neighbours, builds a
per-event difference image (local OOT - IN), fits a free 2-D Gaussian source, estimates
position errors by Monte-Carlo with per-pixel noise, and reports the nearest stars."""
import sys, os, glob, pickle, numpy as np, warnings
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from scipy.optimize import least_squares
sys.path.insert(0, '.')
from search import s3url, fetch
S = '../work'
tic = int(sys.argv[1]); DUR = float(sys.argv[2]) / 24; tcs = [float(x) for x in sys.argv[3:]]
D = f'{S}/cand/{tic}'; os.makedirs(D, exist_ok=True)
m, files = pickle.load(open(f'{S}/ticmap.pkl', 'rb'))
# neighbours
from astroquery.mast import Catalogs
import astropy.units as u
hdr_ra = hdr_dec = None
tp_paths = []
for s in sorted(m[tic]):
    fn = files[(tic, s)].replace('_lc.fits', '_tp.fits'); p = f'{D}/{fn}'
    tp_paths.append(p)
def ensure(p):
    if not os.path.exists(p):
        raw = fetch(s3url(os.path.basename(p).replace('_tp.fits', '_lc.fits')).replace('_lc.fits', '_tp.fits'))
        open(p, 'wb').write(raw)
neigh = None
def gauss(p, yy, xx):
    x, y, A, sg, b = p
    return b + A * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * sg ** 2))
rng = np.random.default_rng(1)
results = []
figs = []
for p in tp_paths:
    # check which events fall in this sector (use lc time range via filename sector)
    with_events = None
    ensure_needed = False
    # quick check with the LC file time span from ticmap sector -> need actual times: open TPF lazily
    ensure(p)
    h = fits.open(p)
    t = h[1].data['TIME']; q = h[1].data['QUALITY']
    tmin, tmax = np.nanmin(t), np.nanmax(t)
    evs = [tc for tc in tcs if tmin - 0.2 < tc < tmax + 0.2]
    if not evs:
        h.close(); os.remove(p); continue
    if neigh is None:
        ra0, de0 = h[0].header['RA_OBJ'], h[0].header['DEC_OBJ']
        r = Catalogs.query_region(f'{ra0} {de0}', radius=0.035 * u.deg, catalog='TIC')
        tm = h[0].header['TESSMAG']
        neigh = r[(r['Tmag'] < tm + 7.5)]
    w = WCS(h[2].header)
    cube = h[1].data['FLUX'].astype(float); err = h[1].data['FLUX_ERR'].astype(float)
    ok = np.isfinite(t) & (q == 0) & np.all(np.isfinite(cube.reshape(len(t), -1)), axis=1)
    t = t[ok]; cube = cube[ok]; err = err[ok]
    sec = h[0].header['SECTOR']
    for tc in evs:
        a = np.abs(t - tc) < DUR * 0.35
        b = (np.abs(t - tc) > DUR * 0.75) & (np.abs(t - tc) < DUR * 0.75 + 0.4)
        if a.sum() < 5 or b.sum() < 10: continue
        # linear baseline per pixel from OOT (both sides) evaluated at in-transit times
        X = np.vstack([np.ones(b.sum()), t[b] - tc]).T
        coef = np.linalg.lstsq(X, cube[b].reshape(b.sum(), -1), rcond=None)[0]
        base_in = (np.vstack([np.ones(a.sum()), t[a] - tc]).T @ coef).reshape((a.sum(),) + cube.shape[1:])
        diff = np.mean(base_in - cube[a], 0)
        # noise per pixel of the in-transit mean: from OOT residual scatter
        resid = cube[b].reshape(b.sum(), -1) - X @ coef
        sig = (np.std(resid, 0) / np.sqrt(a.sum()) * np.sqrt(1 + a.sum() / b.sum())).reshape(diff.shape)
        ny, nx = diff.shape; yy, xx = np.mgrid[0:ny, 0:nx]
        i0 = np.unravel_index(np.argmax(diff), diff.shape)
        fitf = lambda dimg: least_squares(lambda pp: ((gauss(pp, yy, xx) - dimg) / sig).ravel(),
                                          [i0[1], i0[0], max(dimg.max(), 1e-3), 0.8, 0],
                                          bounds=([-1, -1, 0, 0.4, -np.inf], [nx, ny, np.inf, 2.0, np.inf])).x
        best = fitf(diff)
        mc = np.array([fitf(diff + rng.normal(0, sig)) for _ in range(60)])
        ex, ey = np.std(mc[:, 0]), np.std(mc[:, 1])
        sra, sdec = w.pixel_to_world_values(best[0], best[1])
        rows = []
        for st in neigh:
            px, py = w.world_to_pixel_values(st['ra'], st['dec'])
            dpx = np.hypot(px - best[0], py - best[1]); nsig = np.hypot((px - best[0]) / max(ex, 1e-3), (py - best[1]) / max(ey, 1e-3))
            rows.append((dpx * 21, nsig, int(st['ID']), float(st['Tmag'])))
        rows.sort()
        snr_img = diff.sum() / np.sqrt((sig ** 2).sum())
        tgt = [r_ for r_ in rows if r_[2] == tic][0]
        print(f'S{sec} tc={tc:.3f}: source ({best[0]:.2f}±{ex:.2f}, {best[1]:.2f}±{ey:.2f}) px; diff-img SNR~{snr_img:.1f}; '
              f'target offset {tgt[0]:.1f}" ({tgt[1]:.1f} sigma)')
        print('    nearest stars:', [(f'{r_[0]:.1f}"', f'{r_[1]:.1f}sig', r_[2], round(r_[3], 2)) for r_ in rows[:4]])
        results.append(dict(sector=sec, tc=tc, x=best[0], y=best[1], ex=ex, ey=ey, tgt_off=tgt[0], tgt_sig=tgt[1], nearest=rows[:4]))
        fig, ax = plt.subplots(1, 2, figsize=(8, 4))
        ax[0].imshow(np.nanmedian(cube[b], 0), origin='lower'); ax[1].imshow(diff, origin='lower')
        for a_ in ax:
            for st in neigh:
                px, py = w.world_to_pixel_values(st['ra'], st['dec'])
                a_.plot(px, py, 'r*' if int(st['ID']) == tic else 'wx', ms=10)
            a_.set_xlim(-0.5, nx - 0.5); a_.set_ylim(-0.5, ny - 0.5)
        ax[1].errorbar(best[0], best[1], xerr=ex, yerr=ey, color='magenta', capsize=3)
        ax[0].set_title(f'TIC {tic} S{sec} OOT'); ax[1].set_title(f'diff tc={tc:.2f}')
        plt.tight_layout(); plt.savefig(f'{D}/locate_S{sec}_{tc:.1f}.png', dpi=70); plt.close()
    h.close()
pickle.dump(results, open(f'{D}/locate.pkl', 'wb'))
