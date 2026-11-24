"""Quick-look plot: full sector (raw binned + detrended) and zoom around event(s).
usage: quicklook.py out.png TIC:SECTOR:TC [TIC:SECTOR:TC ...]"""
import sys, pickle, numpy as np, functools
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, '.')
from search import detrend, BIN
S = '../work'

@functools.lru_cache(maxsize=4)
def load(sec):
    return pickle.load(open(f'{S}/out/lcs_s{sec:04d}.pkl', 'rb'))

out = sys.argv[1]; items = sys.argv[2:]
fig, axes = plt.subplots(len(items), 2, figsize=(15, 3.2 * len(items)), squeeze=False,
                         gridspec_kw=dict(width_ratios=[2.2, 1]))
for k, it in enumerate(items):
    tic, sec, tc = it.split(':'); tic = int(tic); sec = int(sec); tc = float(tc)
    t0, fb = load(sec)[tic]
    tb = t0 + (np.arange(len(fb)) + 0.5) * BIN
    fd = detrend(tb, fb.astype(float))
    a = axes[k, 0]
    a.plot(tb, fb, '.', ms=1.5, color='0.6')
    a.plot(tb, fd - 0.5 * (np.nanmax(fb) - np.nanmin(fb)) - 0.002, '.', ms=1.5, color='C0')
    a.axvline(tc, color='r', alpha=0.3)
    a.set_title(f'TIC {tic}  S{sec}  (grey raw PDCSAP 10-min, blue detrended offset)', fontsize=9)
    b = axes[k, 1]
    m = np.abs(tb - tc) < 1.5
    b.plot(tb[m] - tc, fb[m], '.', ms=3, color='0.6')
    b.plot(tb[m] - tc, fd[m] - 0.5 * (np.nanmax(fb[m]) - np.nanmin(fb[m])) - 0.001, '.', ms=3, color='C0')
    b.set_title(f'zoom tc={tc:.3f} BTJD', fontsize=9)
plt.tight_layout(); plt.savefig(out, dpi=80)
