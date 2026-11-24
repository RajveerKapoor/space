import pickle, pandas as pd, numpy as np, sys
S = '../work'
r = pd.DataFrame([x for x in pickle.load(open(f'{S}/out/bls_results_long.pkl', 'rb')) if 'P' in x])
print('results:', len(r))
meta = pd.read_pickle(f'{S}/out/star_meta.pkl')
r = r.join(meta[['tmag', 'rad', 'teff', 'logg', 'crowd', 'n']], on='tic')
r['rp'] = np.sqrt(r.depth.clip(0)) * r.rad * 109.1
# expected max duration for a circular orbit (b=0) given stellar density from logg/rad
M = 10 ** r.logg.fillna(4.4) * (r.rad.fillna(1) ** 2) / 10 ** 4.438
rho = M / r.rad.fillna(1) ** 3
r['tmax_h'] = 13 * (r.P / 365.25) ** (1 / 3) * rho ** (-1 / 3) * 1.15
# TCE ephemeris match
tce = pd.read_pickle(f'{S}/cat/tce_all.pkl'); tce = tce[tce.ticid.isin(set(r.tic))]
tg = {k: v for k, v in tce.groupby('ticid')}
def mt(x):
    d = tg.get(x.tic)
    if d is None: return False
    for _, y in d.iterrows():
        P, t0 = y.tce_period, y.tce_time0bt
        if not np.isfinite(P) or P <= 0: continue
        for mult in [1, 2, 0.5, 3, 1 / 3]:
            if abs(x.P / (P * mult) - 1) < 0.01: return True
        ph = ((x.tc - t0) / P + 0.5) % 1 - 0.5
        if abs(ph * P) < 0.15: return True
    return False
r['tce_match'] = r.apply(mt, axis=1)
g = r[(r.ntr >= 3) & (r.frac < 0.55) & (r.sde > 7) & (r.snr > 8) & (r.depth > 0) & (r.depth < 0.03)
      & ((r.oe < 3) | r.oe.isna()) & (r.sec_snr > -3) & (r.sec_snr < 3) & (r.dur_h <= r.tmax_h + 1)]
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
cols = ['tic', 'P', 'tc', 'dur_h', 'tmax_h', 'depth', 'snr', 'sde', 'ntr', 'frac', 'oe', 'sec_snr', 'tmag', 'teff', 'rad', 'rp', 'crowd', 'n', 'tce_match']
print('passing:', len(g), ' not TCE-matched:', (~g.tce_match).sum())
print(g[~g.tce_match].sort_values('sde', ascending=False)[cols].head(int(sys.argv[1]) if len(sys.argv) > 1 else 40).round(3).to_string())
g.to_pickle(f'{S}/out/bls_pass.pkl')
