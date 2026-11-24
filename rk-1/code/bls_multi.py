"""Multi-sector BLS over the stored binned S97-108 light curves for stars with >=2 sectors."""
import sys, os, glob, pickle, numpy as np, time, warnings
warnings.filterwarnings('ignore')
from multiprocessing import Pool
sys.path.insert(0, '.')
from search import BIN
S = '../work'

def prep(t0, fb):
    from wotan import flatten
    tb = t0 + (np.arange(len(fb)) + 0.5) * BIN
    ok = np.isfinite(fb)
    tb, fb = tb[ok], fb[ok].astype(float)
    if len(tb) < 300: return None, None
    fl = flatten(tb, fb, method='biweight', window_length=1.0, break_tolerance=0.4)
    ok = np.isfinite(fl)
    tb, fl = tb[ok], fl[ok]
    dt = np.diff(tb); gi = np.where(dt > 0.5)[0]
    edges = np.concatenate([[tb[0]], tb[gi], tb[gi + 1], [tb[-1]]])
    keep = np.min(np.abs(tb[:, None] - edges[None, :]), axis=1) > 0.3
    tb, fl = tb[keep], fl[keep]
    # bin to 30 min
    k = np.floor(tb / (30 / 1440)).astype(np.int64); u, inv = np.unique(k, return_inverse=True)
    n = np.bincount(inv); bt = np.bincount(inv, tb) / n; bf = np.bincount(inv, fl) / n
    # clip upward outliers
    med = np.median(bf); mad = 1.4826 * np.median(np.abs(bf - med))
    g = bf < med + 4 * mad
    return bt[g], bf[g]

def run(item):
    tic, segs = item
    try:
        from astropy.timeseries import BoxLeastSquares
        T, F = [], []
        for t0, fb in segs:
            a, b = prep(t0, fb)
            if a is not None: T.append(a); F.append(b)
        if len(T) < 2: return None
        t = np.concatenate(T); f = np.concatenate(F)
        o = np.argsort(t); t, f = t[o], f[o]
        sig = 1.4826 * np.median(np.abs(f - np.median(f)))
        base = t.max() - t.min()
        pmax = min(120.0, base / 1.5)
        # frequency grid: oversample relative to shortest duration
        fmin, fmax = 1 / pmax, 1 / 8.0
        # Ofir-style grid: step ~ c f^(2/3) / baseline  (duration ~ P^(1/3))
        c = 0.018
        u = np.arange(fmin ** (1 / 3), fmax ** (1 / 3), c / (3 * base))
        freqs = u ** 3
        periods = 1 / freqs
        bls = BoxLeastSquares(t, f, dy=np.full(len(f), sig))
        r = bls.power(periods, np.array([2.0, 3.5, 6, 9]) / 24, objective='likelihood')
        pw = r.power
        # SDE on detrended power
        sde = (pw.max() - np.median(pw)) / np.std(pw)
        i = np.argmax(pw)
        P, tc, dur, dep, snr = r.period[i], r.transit_time[i], r.duration[i], r.depth[i], r.depth_snr[i]
        # number of transits with data
        ep = np.round((t - tc) / P); intr = np.abs(t - tc - ep * P) < dur / 2
        ntr = len(np.unique(ep[intr]))
        # per-transit depths consistency (max single-transit share of total SNR)
        dd = []
        for e in np.unique(ep[intr]):
            m = intr & (ep == e); dd.append((1 - f[m].mean()) * np.sqrt(m.sum()))
        dd = np.array(dd); frac = dd.max() / dd.sum() if dd.sum() > 0 else 1
        # odd-even
        odd = intr & (ep % 2 == 1); even = intr & (ep % 2 == 0)
        oe = np.nan
        if odd.sum() > 2 and even.sum() > 2:
            oe = abs(f[odd].mean() - f[even].mean()) / (sig * np.sqrt(1 / odd.sum() + 1 / even.sum()))
        # secondary at phase 0.5
        ph = ((t - tc) / P) % 1
        sec = np.abs(ph - 0.5) * P < dur / 2
        secd = (1 - f[sec].mean()) / (sig / np.sqrt(sec.sum())) if sec.sum() > 2 else np.nan
        return dict(tic=tic, P=P, tc=tc, dur_h=dur * 24, depth=dep, snr=snr, sde=sde, ntr=ntr, frac=frac,
                    oe=oe, sec_snr=secd, npts=len(t), base=base, sig30=sig)
    except Exception as ex:
        return dict(tic=tic, err=repr(ex)[:100])

if __name__ == '__main__':
    tics_ok = pickle.load(open(sys.argv[1], 'rb'))  # set of TICs to run
    nproc = int(sys.argv[2]); out = sys.argv[3]
    segs = {}
    for f in sorted(glob.glob(f'{S}/out/lcs_s*.pkl')):
        L = pickle.load(open(f, 'rb'))
        for tic, v in L.items():
            if tic in tics_ok: segs.setdefault(tic, []).append(v)
        del L
    items = [(k, v) for k, v in segs.items() if len(v) >= 2]
    import pandas as pd
    meta = pd.read_pickle(f'{S}/out/star_meta.pkl')
    pr = {k: (-(len(v) >= 3), float(meta.tmag.get(k, 20))) for k, v in items}
    items.sort(key=lambda kv: pr[kv[0]])
    print('stars to search:', len(items), flush=True)
    res = []; t0 = time.time()
    with Pool(nproc) as p:
        for k, r in enumerate(p.imap_unordered(run, items, chunksize=2)):
            if r: res.append(r)
            if k % 500 == 0:
                print(f'{k}/{len(items)} {time.time()-t0:.0f}s', flush=True)
                pickle.dump(res, open(out, 'wb'))
    pickle.dump(res, open(out, 'wb'))
    print('done', len(res), f'{time.time()-t0:.0f}s')
