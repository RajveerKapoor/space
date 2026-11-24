"""Aggregate single-event search results, reject systematics, cross-match catalogs."""
import pickle, glob, sys, numpy as np, pandas as pd, collections
S = '../work'
excl = pickle.load(open(f'{S}/cat/excl.pkl', 'rb'))
SNR_MIN = float(sys.argv[1]) if len(sys.argv) > 1 else 9.0

rows = []
for f in sorted(glob.glob(f'{S}/out/events_s*.pkl')):
    res = pickle.load(open(f, 'rb'))
    sec = int(f[-8:-4])
    # systematics map: event times shared by many stars
    allt = np.array([e[1] for m in res for e in m['events'] if e[0] > 7])
    bins = np.arange(allt.min() - 1, allt.max() + 1, 0.1) if len(allt) else np.array([0, 1])
    h, _ = np.histogram(allt, bins)
    thr = max(6, 4 * np.median(h[h > 0]) if (h > 0).any() else 6)
    bad = bins[:-1][h > thr]
    for m in res:
        for e in m['events']:
            snr, tc, wh, dep, gd = e
            nearbad = bool(len(bad)) and np.min(np.abs(bad + 0.05 - tc)) < 0.15
            rows.append(dict(tic=m['tic'], sector=sec, snr=snr, tc=tc, dur_h=wh, depth=dep, gapd=gd,
                             sysbad=nearbad, tmag=m['tmag'], teff=m['teff'], logg=m['logg'], rad=m['rad'],
                             crowd=m['crowd'], sig10=m['sig10'], nev=len(m['events']), ra=m['ra'], dec=m['dec']))
    print(f'sector {sec}: {len(res)} stars, {len(bad)} systematic epochs flagged (thr={thr:.0f})')

ev = pd.DataFrame(rows)
for k in ['toi', 'ctoi', 'ps', 'tce']:
    ev['in_' + k] = ev.tic.isin(excl[k])
ev.to_pickle(f'{S}/out/all_events.pkl')
good = ev[(ev.snr >= SNR_MIN) & (ev.gapd >= 0.25) & (~ev.sysbad) & (ev.depth < 0.05) & (ev.depth > 0)]
print('events passing basic cuts:', len(good), ' stars:', good.tic.nunique())
new = good[~(good.in_toi | good.in_ctoi | good.in_ps)]
print('  not in TOI/CTOI/confirmed:', len(new), ' stars:', new.tic.nunique())
# stars with very many events are variables/EBs -> flag
nper = new.groupby('tic').size()
new = new.assign(n_star_events=new.tic.map(nper))
new.to_pickle(f'{S}/out/new_events.pkl')
print(new.sort_values('snr', ascending=False).head(40).to_string())

# --- match events against TCE ephemerides (any sector/multisector run)
tce = pd.read_pickle(f'{S}/cat/tce_all.pkl')
tce = tce[tce.ticid.isin(set(new.tic))]
tg = {k: v for k, v in tce.groupby('ticid')}
def match_tce(r):
    d = tg.get(r.tic)
    if d is None: return False
    for _, x in d.iterrows():
        P, t0, dur = x.tce_period, x.tce_time0bt, x.tce_duration / 24
        if not np.isfinite(P) or P <= 0: continue
        ph = ((r.tc - t0) / P + 0.5) % 1 - 0.5
        if abs(ph * P) < max(dur, 0.1) or abs(abs(ph) - 0.5) * P < max(dur, 0.1):  # primary or secondary
            return True
    return False
new['tce_match'] = new.apply(match_tce, axis=1)
new.to_pickle(f'{S}/out/new_events.pkl')
un = new[~new.tce_match]
print('\n==== events NOT matching any TCE ephemeris:', len(un), 'stars', un.tic.nunique())
cols = ['tic','sector','snr','tc','dur_h','depth','gapd','tmag','teff','logg','rad','crowd','sig10','n_star_events','in_tce']
print(un.sort_values('snr', ascending=False)[cols].head(60).to_string())
