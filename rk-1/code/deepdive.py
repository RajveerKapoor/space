"""Deep-dive on one TIC: fetch all 2-min LCs, detrend, BLS over combined data, fold, diagnostics.
usage: deepdive.py TIC [P_guess T0_guess]"""
import sys, io, os, pickle, numpy as np, warnings
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares
sys.path.insert(0, '.')
from search import s3url, fetch
from wotan import flatten
S = '../work'
tic = int(sys.argv[1])
D = f'{S}/cand/{tic}'; os.makedirs(D, exist_ok=True)
m, files = pickle.load(open(f'{S}/ticmap.pkl', 'rb'))
secs = sorted(m[tic])
T, F, E, SEC, C1, C2, SAP = [], [], [], [], [], [], []
for s in secs:
    fn = files[(tic, s)]; p = f'{D}/{fn}'
    if not os.path.exists(p):
        open(p, 'wb').write(fetch(s3url(fn)))
    with fits.open(p) as h:
        d = h[1].data
        t = d['TIME']; f = d['PDCSAP_FLUX']; e = d['PDCSAP_FLUX_ERR']; q = d['QUALITY']
        ok = np.isfinite(t) & np.isfinite(f) & (q == 0)
        med = np.nanmedian(f[ok])
        T.append(t[ok]); F.append(f[ok] / med); E.append(e[ok] / med); SEC.append(np.full(ok.sum(), s))
        C1.append(d['MOM_CENTR1'][ok] - np.nanmedian(d['MOM_CENTR1'][ok]))
        C2.append(d['MOM_CENTR2'][ok] - np.nanmedian(d['MOM_CENTR2'][ok]))
        SAP.append(d['SAP_FLUX'][ok] / np.nanmedian(d['SAP_FLUX'][ok]))
        print(f'S{s}: {ok.sum()} pts, {t[ok].min():.2f}-{t[ok].max():.2f}, CROWDSAP={h[1].header.get("CROWDSAP")}, FLFRCSAP={h[1].header.get("FLFRCSAP")}')
t = np.concatenate(T); f = np.concatenate(F); e = np.concatenate(E); sec = np.concatenate(SEC)
c1 = np.concatenate(C1); c2 = np.concatenate(C2); sap = np.concatenate(SAP)
fl = flatten(t, f, method='biweight', window_length=1.0, break_tolerance=0.4)
ok = np.isfinite(fl)
# clip 0.3 d around segment edges (gaps > 0.5 d)
dt = np.diff(t); gi = np.where(dt > 0.5)[0]
edges = np.concatenate([[t[0]], t[gi], t[gi + 1], [t[-1]]])
near_edge = np.min(np.abs(t[:, None] - edges[None, :]), axis=1) < 0.3
ok &= ~near_edge
print('clipped near-edge points:', near_edge.sum())
np.savez(f'{D}/lc.npz', t=t, f=f, e=e, sec=sec, fl=fl, c1=c1, c2=c2, sap=sap)
# bin to 10 min for BLS speed
tb = np.round(t / (10 / 1440)); u, inv = np.unique(tb, return_inverse=True)
bt = np.bincount(inv, t * ok) / np.maximum(np.bincount(inv, ok), 1)
bf = np.bincount(inv, np.where(ok, fl, 0)) / np.maximum(np.bincount(inv, ok), 1)
good = np.bincount(inv, ok) > 2
bt, bf = bt[good], bf[good]
sig = 1.4826*np.median(np.abs(bf-np.median(bf)))
bls = BoxLeastSquares(bt, bf, dy=np.full(len(bf), sig))
pmin, pmax = 1.0, min(60.0, (bt.max() - bt.min()) / 2)
if len(sys.argv) > 2:
    pg = float(sys.argv[2]); periods = np.linspace(pg * 0.995, pg * 1.005, 20000)
else:
    periods = np.exp(np.linspace(np.log(pmin), np.log(pmax), 300000))
durs = np.array([1.5, 2, 3, 4, 5, 6, 8]) / 24
r = bls.power(periods, durs, objective='likelihood')
i = np.argmax(r.power)
P, t0, dur, dep = r.period[i], r.transit_time[i], r.duration[i], r.depth[i]
print(f'BLS best: P={P:.5f} t0={t0:.4f} dur={dur*24:.2f}h depth={dep*1e6:.0f}ppm snr={r.depth_snr[i]:.1f} power={r.power[i]:.1f}')
order = np.argsort(r.power)[::-1]
shown = []
for j in order:
    if all(abs(r.period[j] / p - 1) > 0.01 for p in shown):
        shown.append(r.period[j]); print(f'   P={r.period[j]:.5f} pow={r.power[j]:.1f} dep={r.depth[j]*1e6:.0f}')
    if len(shown) >= 8: break
st = bls.compute_stats(P, dur, t0)
print('odd/even depth:', st['depth_odd'], st['depth_even'])
print('per-transit count:', st['per_transit_count'])
print('transit times:', np.round(st['transit_times'][st['per_transit_count'] > 0], 3))
pickle.dump(dict(P=P, t0=t0, dur=dur, dep=dep, periods=r.period, power=r.power), open(f'{D}/bls.pkl', 'wb'))
# plots
fig = plt.figure(figsize=(15, 11))
ax = fig.add_subplot(3, 1, 1); ax.plot(r.period, r.power, lw=0.5); ax.set_xscale('log'); ax.set_xlabel('P [d]'); ax.set_title(f'TIC {tic} BLS')
ph = ((t - t0) / P + 0.5) % 1 - 0.5
ax = fig.add_subplot(3, 2, 3); m_ = np.abs(ph * P) < 0.6
ax.plot(ph[m_ & ok] * P * 24, fl[m_ & ok], '.', ms=1, color='0.7')
bins = np.linspace(-0.6, 0.6, 80) * 24; idx = np.digitize(ph * P * 24, bins)
bm = [np.nanmedian(fl[(idx == k) & m_ & ok]) if ((idx == k) & m_ & ok).sum() > 3 else np.nan for k in range(1, len(bins))]
ax.plot(0.5 * (bins[1:] + bins[:-1]), bm, 'o', color='C3', ms=3); ax.set_xlabel('hours from mid-transit'); ax.set_title('phase-folded (2-min, detrended)')
ax = fig.add_subplot(3, 2, 4)
ph2 = ((t - t0) / P) % 1 - 0.5
m2 = np.abs(ph2 * P) < 0.6
ax.plot(ph2[m2 & ok] * P * 24, fl[m2 & ok], '.', ms=1, color='0.7'); ax.set_title('phase 0.5 (secondary eclipse?)'); ax.set_xlabel('hours from phase 0.5')
idx2 = np.digitize(ph2 * P * 24, bins)
bm2 = [np.nanmedian(fl[(idx2 == k) & m2 & ok]) if ((idx2 == k) & m2 & ok).sum() > 3 else np.nan for k in range(1, len(bins))]
ax.plot(0.5 * (bins[1:] + bins[:-1]), bm2, 'o', color='C3', ms=3)
tt = st['transit_times'][st['per_transit_count'] > 0]
n = len(tt)
for k, tc in enumerate(tt[:8]):
    ax = fig.add_subplot(3, 8, 17 + k)
    w = np.abs(t - tc) < 0.5
    ax.plot((t[w] - tc) * 24, fl[w], '.', ms=1, color='0.6')
    ax.set_title(f'{tc:.2f} S{sec[w][0] if w.any() else ""}', fontsize=8); ax.set_ylim(1 - 3 * dep, 1 + 1.5 * dep)
plt.tight_layout(); plt.savefig(f'{D}/deepdive.png', dpi=75)
