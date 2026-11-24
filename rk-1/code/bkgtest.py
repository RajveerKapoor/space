"""Background-event and centroid tests for a list of single events (needs raw LC files)."""
import sys, io, pickle, numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, '.')
from search import s3url, fetch
S = '../work'
m, files = pickle.load(open(f'{S}/ticmap.pkl', 'rb'))

def test(row):
    tic, sec, tc, dur = row
    try:
        from astropy.io import fits
        raw = fetch(s3url(files[(tic, sec)]))
        h = fits.open(io.BytesIO(raw)); d = h[1].data
        t = d['TIME']; q = d['QUALITY']
        cols = {k: np.array(d[k], float) for k in ['SAP_BKG', 'SAP_FLUX', 'PDCSAP_FLUX', 'MOM_CENTR1', 'MOM_CENTR2']}
        ok = np.isfinite(t) & (q == 0)
        for v in cols.values(): ok &= np.isfinite(v)
        t = t[ok]; c = {k: v[ok] for k, v in cols.items()}
        D = dur / 24
        a = np.abs(t - tc) < D * 0.35
        b = (np.abs(t - tc) > D * 0.75) & (np.abs(t - tc) < D * 0.75 + 0.4)
        if a.sum() < 4 or b.sum() < 10: return (tic, sec, tc, np.nan, np.nan, np.nan, np.nan)
        def lin(y):  # local linear baseline from OOT, evaluate in-transit mean deficit
            X = np.vstack([np.ones(b.sum()), t[b] - tc]).T
            co = np.linalg.lstsq(X, y[b], rcond=None)[0]
            pred = co[0] + co[1] * (t[a] - tc)
            res = y[b] - X @ co
            return np.mean(pred - y[a]), np.std(res) / np.sqrt(a.sum())
        dsap, esap = lin(c['SAP_FLUX'])
        dbk, ebk = lin(-c['SAP_BKG'])   # positive if background rose
        dc1, ec1 = lin(c['MOM_CENTR1']); dc2, ec2 = lin(c['MOM_CENTR2'])
        csig = np.hypot(dc1 / ec1, dc2 / ec2)
        return (tic, sec, tc, dsap / np.median(c['SAP_FLUX']), -dbk / dsap if dsap != 0 else np.nan, dsap / esap, csig)
    except Exception as ex:
        return (tic, sec, tc, np.nan, np.nan, np.nan, np.nan)

if __name__ == '__main__':
    df = pd.read_pickle(sys.argv[1]); out = sys.argv[2]
    rows = list(zip(df.tic.astype(int), df.sector.astype(int), df[sys.argv[3]].astype(float), df[sys.argv[4]].astype(float)))
    with Pool(24) as p:
        R = p.map(test, rows, chunksize=2)
    r = pd.DataFrame(R, columns=['tic', 'sector', 'tc', 'sap_depth', 'bkg_ratio', 'sap_snr', 'centroid_sig'])
    df = df.reset_index(drop=True)
    df[['sap_depth', 'bkg_ratio', 'sap_snr', 'centroid_sig']] = r[['sap_depth', 'bkg_ratio', 'sap_snr', 'centroid_sig']]
    df.to_pickle(out); print('done', len(df))
