"""Stream TESS SPOC 2-min light curves from the MAST S3 bucket, bin, detrend,
and run a single-event (box) transit search. Stores binned raw PDCSAP flux for
every target plus a per-target event summary."""
import io, os, sys, re, pickle, time, json, traceback
import numpy as np
import urllib.request
from multiprocessing import Pool

BIN = 10.0 / 1440.0            # 10-minute bins
WIDTHS_H = [1, 1.5, 2, 3, 4, 5, 6, 8, 10, 13]

def s3url(fn):
    m = re.match(r'tess\d+-s(\d{4})-(\d{16})-\d{4}-s_lc\.fits', fn)
    sec, t = m.group(1), m.group(2)
    return (f"https://stpubdata.s3.amazonaws.com/tess/public/tid/s{sec}/"
            f"{t[0:4]}/{t[4:8]}/{t[8:12]}/{t[12:16]}/{fn}")

def fetch(url, tries=4):
    for i in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=120).read()
        except Exception:
            time.sleep(2 * (i + 1))
    return None

def bin_lc(t, f):
    t0 = np.floor(t.min() * 10) / 10
    idx = ((t - t0) / BIN).astype(int)
    n = idx.max() + 1
    s = np.bincount(idx, weights=f, minlength=n)
    c = np.bincount(idx, minlength=n)
    with np.errstate(invalid='ignore', divide='ignore'):
        fb = s / c
    fb[c < 2] = np.nan
    return t0, fb.astype(np.float32)

def detrend(tb, fb):
    from wotan import flatten
    ok = np.isfinite(fb)
    out = np.full_like(fb, np.nan)
    if ok.sum() < 200:
        return out
    fl = flatten(tb[ok], fb[ok], method='biweight', window_length=1.0,
                 break_tolerance=0.4, edge_cutoff=0.0, return_trend=False)
    out[ok] = fl
    return out

def gap_distance(tb, ok):
    """distance (days) from each bin to the nearest missing stretch > 0.2 d or edge"""
    valid_t = tb[ok]
    if len(valid_t) < 2:
        return np.zeros_like(tb)
    d = np.diff(valid_t)
    gi = np.where(d > 0.2)[0]
    edges = np.concatenate([[valid_t[0]], valid_t[gi], valid_t[gi + 1], [valid_t[-1]]])
    return np.min(np.abs(tb[:, None] - edges[None, :]), axis=1)

def box_search(tb, fd):
    ok = np.isfinite(fd)
    x = np.where(ok, fd - 1.0, 0.0)
    w_ok = ok.astype(float)
    cs = np.concatenate([[0], np.cumsum(x)])
    cc = np.concatenate([[0], np.cumsum(w_ok)])
    mad = np.nanmedian(np.abs(fd[ok] - np.nanmedian(fd[ok])))
    sig_white = 1.4826 * mad
    gd = gap_distance(tb, ok)
    events = []
    for wh in WIDTHS_H:
        w = max(1, int(round(wh / 24 / BIN)))
        if w >= len(x):
            continue
        ssum = cs[w:] - cs[:-w]
        scnt = cc[w:] - cc[:-w]
        good = scnt >= 0.7 * w
        if good.sum() < 50:
            continue
        mean = np.where(good, ssum / np.maximum(scnt, 1), np.nan)
        m = mean[good]
        sig_red = 1.4826 * np.median(np.abs(m - np.median(m)))
        sig = max(sig_red, sig_white / np.sqrt(w))
        snr = -mean / sig
        snr[~good] = 0
        centers = tb[:len(snr)] + (w - 1) * BIN / 2
        # pick top distinct events for this width
        order = np.argsort(snr)[::-1]
        picked = []
        for i in order[:2000]:
            if snr[i] < 6:
                break
            if all(abs(centers[i] - p) > 0.5 for p in picked):
                picked.append(centers[i])
                ci = min(len(gd) - 1, i + w // 2)
                events.append((float(snr[i]), float(centers[i]), wh, float(-mean[i]),
                               float(gd[ci])))
            if len(picked) >= 5:
                break
    # merge across widths: keep best SNR per 0.5 d cluster
    events.sort(key=lambda e: -e[0])
    merged = []
    for e in events:
        if all(abs(e[1] - m[1]) > 0.5 for m in merged):
            merged.append(e)
    return merged[:8], float(sig_white)

def process(fn):
    try:
        raw = fetch(s3url(fn))
        if raw is None:
            return ('ERR', fn, 'download')
        from astropy.io import fits
        with fits.open(io.BytesIO(raw), memmap=False) as h:
            h0 = h[0].header; h1 = h[1].header; d = h[1].data
            t = np.array(d['TIME'], float); f = np.array(d['PDCSAP_FLUX'], float)
            q = np.array(d['QUALITY'])
            meta = dict(tic=int(h0.get('TICID', 0)), sector=int(h0.get('SECTOR', 0)),
                        tmag=h0.get('TESSMAG'), teff=h0.get('TEFF'), logg=h0.get('LOGG'),
                        rad=h0.get('RADIUS'), ra=h0.get('RA_OBJ'), dec=h0.get('DEC_OBJ'),
                        crowd=h1.get('CROWDSAP'), cdpp1=h1.get('CDPP1_0'), fn=fn)
        ok = np.isfinite(t) & np.isfinite(f) & (q == 0)
        if ok.sum() < 1000:
            return ('ERR', fn, 'few points')
        t = t[ok]; f = f[ok] / np.nanmedian(f[ok])
        t0, fb = bin_lc(t, f)
        tb = t0 + (np.arange(len(fb)) + 0.5) * BIN
        fd = detrend(tb, fb)
        # clip upward outliers
        ok2 = np.isfinite(fd)
        if ok2.sum() < 200:
            return ('ERR', fn, 'detrend')
        med = np.nanmedian(fd); mad = 1.4826 * np.nanmedian(np.abs(fd - med))
        fd[fd > med + 5 * mad] = np.nan
        events, sw = box_search(tb, fd)
        meta['sig10'] = sw
        meta['events'] = events
        return ('OK', meta, (t0, fb))
    except Exception as ex:
        return ('ERR', fn, repr(ex)[:200])

if __name__ == '__main__':
    sectors = [int(s) for s in sys.argv[1].split(',')]; outdir = sys.argv[2]; nproc = int(sys.argv[3])
    files = pickle.load(open(sys.argv[4], 'rb'))[1]
    os.makedirs(outdir, exist_ok=True)
    with Pool(nproc) as pool:
        for sector in sectors:
            if os.path.exists(f'{outdir}/events_s{sector:04d}.pkl'):
                continue
            fns = sorted(fn for (tic, sec), fn in files.items() if sec == sector)
            print(f'sector {sector}: {len(fns)} files', flush=True)
            res, lcs, errs = [], {}, []
            t_start = time.time()
            for k, r in enumerate(pool.imap_unordered(process, fns, chunksize=4)):
                if r[0] == 'OK':
                    res.append(r[1]); lcs[r[1]['tic']] = r[2]
                else:
                    errs.append(r[1:])
                if k % 2000 == 0:
                    print(f'  {k}/{len(fns)} {time.time()-t_start:.0f}s errs={len(errs)}', flush=True)
            pickle.dump(res, open(f'{outdir}/events_s{sector:04d}.pkl', 'wb'))
            pickle.dump(lcs, open(f'{outdir}/lcs_s{sector:04d}.pkl', 'wb'))
            json.dump(errs, open(f'{outdir}/errs_s{sector:04d}.json', 'w'))
            print(f'done sector {sector}: ok={len(res)} err={len(errs)} {time.time()-t_start:.0f}s', flush=True)
