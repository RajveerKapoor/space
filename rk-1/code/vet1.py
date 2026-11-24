"""Stage-1 automated vetting of single events using stored binned LCs."""
import sys, pickle, glob, numpy as np, pandas as pd, functools, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, '.')
from search import detrend, BIN
S = '../work'

ev = pd.read_pickle(f'{S}/out/all_events.pkl')
new = pd.read_pickle(f'{S}/out/new_events.pkl')

# ---- coincidence excess: events in other stars at same time (same sector)
ref = ev[(ev.snr > 7) & (ev.gapd > 0.1)]
coinc = {}
for sec, d in ref.groupby('sector'):
    t = np.sort(d.tc.values)
    coinc[sec] = t
def excess(r):
    t = coinc[r.sector]
    near = np.sum(np.abs(t - r.tc) < 0.1) - 1
    bg = (np.sum(np.abs(t - r.tc) < 3.0) - near - 1) / (5.8 / 0.2)
    return near, bg

cand = new[(new.snr >= 8) & (new.gapd >= 0.5) & (~new.tce_match) &
           ((new.rad < 3) | new.rad.isna()) & (new.crowd > 0.7)].copy()
cand[['n_coinc', 'bg_coinc']] = cand.apply(lambda r: pd.Series(excess(r)), axis=1)
cand = cand[cand.n_coinc <= np.maximum(3, 3 * cand.bg_coinc)]
# epochs repeated among candidate stars (same sector, within 0.06 d) -> systematics
rep = cand.apply(lambda r: ((cand.sector == r.sector) & (cand.tic != r.tic) & ((cand.tc - r.tc).abs() < 0.06)).sum(), axis=1)
cand = cand[rep < 2]
print('candidates after coincidence cut:', len(cand), cand.tic.nunique())

@functools.lru_cache(maxsize=2)
def load(sec):
    return pickle.load(open(f'{S}/out/lcs_s{sec:04d}.pkl', 'rb'))

out = []
for sec, d in cand.groupby('sector'):
    lcs = load(sec)
    for _, r in d.iterrows():
        t0, fb = lcs[r.tic]; fb = fb.astype(float)
        tb = t0 + (np.arange(len(fb)) + 0.5) * BIN
        fd = detrend(tb, fb)
        ok = np.isfinite(fd)
        sig = 1.4826 * np.nanmedian(np.abs(fd[ok] - np.nanmedian(fd[ok])))
        D = r.dur_h / 24
        # refine center & duration on detrended flux
        best = None
        for Dt in D * np.array([0.6, 0.8, 1.0, 1.25, 1.5]):
            for tcc in r.tc + np.linspace(-D / 2, D / 2, 13):
                inn = ok & (np.abs(tb - tcc) < Dt / 2)
                if inn.sum() < 3: continue
                dep = 1 - np.mean(fd[inn]); s = dep / (sig / np.sqrt(inn.sum()))
                if best is None or s > best[0]: best = (s, tcc, Dt, dep, inn.sum())
        if best is None: continue
        s, tcc, Dt, dep, nin = best
        # local raw-flux fit: quadratic baseline + box
        w = ok & np.isfinite(fb) & (np.abs(tb - tcc) < max(3 * Dt, 0.6))
        inn = np.abs(tb - tcc) < Dt / 2
        X = np.vstack([np.ones(w.sum()), tb[w] - tcc, (tb[w] - tcc) ** 2, inn[w].astype(float)]).T
        coef, *_ = np.linalg.lstsq(X, fb[w], rcond=None)
        res = fb[w] - X @ coef
        soot = np.std(res[~inn[w]]) if (~inn[w]).sum() > 5 else np.nan
        dep_raw = -coef[3] / coef[0]
        snr_raw = dep_raw / (soot / np.sqrt(max(inn[w].sum(), 1)))
        # step test: pre vs post windows of width Dt
        pre = ok & (tb > tcc - Dt / 2 - Dt) & (tb < tcc - Dt / 2)
        post = ok & (tb > tcc + Dt / 2) & (tb < tcc + Dt / 2 + Dt)
        step = (np.mean(fd[pre]) - np.mean(fd[post])) / dep if pre.sum() > 2 and post.sum() > 2 else np.nan
        # flare / upward excursion near event
        near = ok & (np.abs(tb - tcc) < 0.6)
        upmax = np.nanmax((fd[near] - 1) / sig) if near.sum() else np.nan
        # local scatter vs global
        loc = ok & (np.abs(tb - tcc) < 1.0) & ~inn
        lscat = np.std(fd[loc]) / sig if loc.sum() > 10 else np.nan
        out.append(dict(r) | dict(tc_ref=tcc, dur_ref_h=Dt * 24, dep_ref=dep, snr_ref=s,
                                  dep_raw=dep_raw, snr_raw=snr_raw, step=step, upmax=upmax,
                                  lscat=lscat, nin=nin))
V = pd.DataFrame(out)
V['rp_re'] = np.sqrt(np.clip(V.dep_ref, 0, None)) * V.rad * 109.1
V.to_pickle(f'{S}/out/vet1.pkl')
good = V[(V.snr_ref > 8) & (V.snr_raw > 7) & (V.step.abs() < 0.4) & (V.upmax < 5) & (V.lscat < 1.4)]
print('pass stage-1:', len(good), good.tic.nunique())
cols = ['tic', 'sector', 'tc_ref', 'dur_ref_h', 'dep_ref', 'snr_ref', 'snr_raw', 'step', 'upmax', 'lscat',
        'tmag', 'teff', 'logg', 'rad', 'rp_re', 'crowd', 'n_coinc', 'in_tce', 'n_star_events']
pd.set_option('display.width', 250)
print(good.sort_values('snr_ref', ascending=False)[cols].to_string())
good.to_pickle(f'{S}/out/vet1_good.pkl')
