"""Shape filter for BLS candidates: reject folds whose out-of-transit baseline is structured
(sinusoidal variability aliases), keep box-like transits."""
import sys, glob, pickle, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, '.')
from bls_multi import prep
S = '../work'
r = pd.DataFrame([x for x in pickle.load(open(f'{S}/out/bls_results_long.pkl', 'rb')) if 'P' in x])
meta = pd.read_pickle(f'{S}/out/star_meta.pkl')
r = r.join(meta[['tmag', 'rad', 'teff', 'logg', 'crowd', 'n']], on='tic')
c = r[(r.ntr >= 3) & (r.frac < 0.6) & (r.sde > 6.5) & (r.snr > 7.5) & (r.depth > 0) & (r.depth < 0.05)
      & ((r.oe < 3.5) | r.oe.isna()) & (r.sec_snr.abs() < 4)]
print('candidates to shape-test:', len(c))
want = set(c.tic)
segs = {}
for f in sorted(glob.glob(f'{S}/out/lcs_s*.pkl')):
    L = pickle.load(open(f, 'rb'))
    for t in want & set(L): segs.setdefault(t, []).append(L[t])
    del L
out = []
for _, x in c.iterrows():
    T, F = [], []
    for t0, fb in segs[x.tic]:
        a, b = prep(t0, fb)
        if a is not None: T.append(a); F.append(b)
    t = np.concatenate(T); f = np.concatenate(F)
    sig = 1.4826 * np.median(np.abs(f - np.median(f)))
    D = x.dur_h / 24
    dt = ((t - x.tc) / x.P + 0.5) % 1 * x.P - 0.5 * x.P
    intr = np.abs(dt) < D / 2
    dep = 1 - f[intr].mean(); dep_err = sig / np.sqrt(intr.sum())
    # OOT structure: bins of width D/2 from 0.75D to 4D on both sides
    edges = np.arange(0.75 * D, 4 * D + 1e-9, D / 2)
    vals = []
    for s in (-1, 1):
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (s * dt > lo) & (s * dt <= hi)
            if m.sum() > 3: vals.append(((f[m].mean() - 1) / (sig / np.sqrt(m.sum()))) ** 2)
    chi2r = np.mean(vals) if vals else np.nan
    # ingress shoulders: mean of the two D/2 bins adjacent to the transit
    sh = np.abs(dt) - D / 2
    shoulder = (1 - f[(sh > 0) & (sh < D / 2)].mean()) / dep if dep > 0 else np.nan
    out.append(dict(x) | dict(dep2=dep, dep2_snr=dep / dep_err, chi2r=chi2r, shoulder=shoulder))
o = pd.DataFrame(out)
o['rp'] = np.sqrt(o.dep2.clip(0)) * o.rad * 109.1
o.to_pickle(f'{S}/out/bls_shape.pkl')
g = o[(o.chi2r < 2.5) & (o.dep2_snr > 8) & (o.shoulder.abs() < 0.35)]
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
cols = ['tic', 'P', 'tc', 'dur_h', 'dep2', 'dep2_snr', 'sde', 'ntr', 'frac', 'oe', 'sec_snr', 'chi2r', 'shoulder', 'tmag', 'teff', 'rad', 'rp', 'crowd', 'n']
print('pass shape filter:', len(g))
print(g.sort_values('dep2_snr', ascending=False)[cols].round(3).to_string())
