"""MCMC fit of a single transit (circular orbit) with stellar-density prior -> period posterior."""
import sys, glob, numpy as np, warnings, pickle
warnings.filterwarnings('ignore')
import batman, emcee
from astropy.io import fits
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
S = '../work'
tic = int(sys.argv[1]); sec = int(sys.argv[2]); tc0 = float(sys.argv[3])
rho_mu, rho_sd = float(sys.argv[4]), float(sys.argv[5])      # solar units
u1m, u2m = float(sys.argv[6]), float(sys.argv[7])
Rstar = float(sys.argv[8])
fn = glob.glob(f'{S}/cand/{tic}/*s{sec:04d}*_lc.fits')[0]
d = fits.open(fn)[1].data
t = d['TIME']; f = d['PDCSAP_FLUX']; e = d['PDCSAP_FLUX_ERR']; q = d['QUALITY']
ok = np.isfinite(t) & np.isfinite(f) & (q == 0)
t, f, e = t[ok], f[ok], e[ok]
w = np.abs(t - tc0) < 0.6
t, f, e = t[w], f[w], e[w]
oot = np.abs(t - tc0) > 0.12
c = np.polyfit(t[oot] - tc0, f[oot], 2)
base = np.polyval(c, t - tc0)
y = f / base; ey = e / base
sig_oot = np.std(y[oot]); ey = np.full_like(y, sig_oot)
print(f'points {len(t)}, OOT scatter {sig_oot*1e3:.2f} ppt (2-min)')
params = batman.TransitParams()
params.t0 = tc0; params.per = 100.; params.rp = 0.15; params.a = 100.; params.inc = 89.9
params.ecc = 0.; params.w = 90.; params.limb_dark = 'quadratic'; params.u = [u1m, u2m]
mdl = batman.TransitModel(params, t, supersample_factor=3, exp_time=2 / 1440)
G = 6.674e-11; rho_sun = 1408.0
def model(th):
    t0, k, b, logP, rho, u1, u2 = th
    P = 10 ** logP
    aR = (rho * rho_sun * G * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
    inc = np.degrees(np.arccos(b / aR))
    params.t0 = t0; params.rp = k; params.per = P; params.a = aR; params.inc = inc; params.u = [u1, u2]
    return mdl.light_curve(params)
def lnp(th):
    t0, k, b, logP, rho, u1, u2 = th
    if not (abs(t0 - tc0) < 0.05 and 0.01 < k < 0.6 and 0 <= b < 1 + k and 0.3 < logP < 3.7 and 0.05 < rho < 500): return -np.inf
    if not (0 < u1 + u2 < 1 and u1 > 0): return -np.inf
    lp = -np.log(rho) - 0.5 * ((u1 - u1m) / 0.1) ** 2 - 0.5 * ((u2 - u2m) / 0.1) ** 2
    m = model(th)
    return lp - 0.5 * np.sum(((y - m) / ey) ** 2)
ndim, nw = 7, 40
p0 = np.array([tc0, 0.16, 0.5, 1.7, rho_mu, u1m, u2m])
pos = p0 + np.array([1e-3, 0.01, 0.2, 0.3, 0.1 * rho_mu, 0.02, 0.02]) * np.random.randn(nw, ndim)
pos[:, 2] = np.abs(pos[:, 2])
sam = emcee.EnsembleSampler(nw, ndim, lnp)
sam.run_mcmc(pos, 6000, progress=False)
ch = sam.get_chain(discard=2500, thin=5, flat=True)
lnl = sam.get_log_prob(discard=2500, thin=5, flat=True)
names = ['t0', 'k', 'b', 'logP', 'rho', 'u1', 'u2']
res = {}
for i, n in enumerate(names):
    lo, md, hi = np.percentile(ch[:, i], [16, 50, 84]); res[n] = (md, md - lo, hi - md)
    print(f'{n:5s} = {md:.5f} -{md-lo:.5f} +{hi-md:.5f}')
k = ch[:, 1]; P = 10 ** ch[:, 3]
rp = k * Rstar * 109.076
print('Rp/R* depth k^2 = %.4f' % np.median(k ** 2))
print('Rp [R_E] = %.2f (16-84%%: %.2f - %.2f);  Rp [R_J] = %.2f' % (np.median(rp), *np.percentile(rp, [16, 84]), np.median(rp) / 11.21))
print('P [d] median %.1f, 16-84%%: %.1f - %.1f, 2.5-97.5%%: %.1f - %.1f' % (np.median(P), *np.percentile(P, [16, 84]), *np.percentile(P, [2.5, 97.5])))
best = ch[np.argmax(lnl)]
mb = model(best)
# T14
th = best; Pb = 10 ** th[3]; aR = (th[4] * rho_sun * G * (Pb * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
T14 = Pb / np.pi * np.arcsin(np.sqrt((1 + th[1]) ** 2 - th[2] ** 2) / aR / np.sin(np.arccos(th[2] / aR))) * 24
print(f'best-fit: k={th[1]:.3f} b={th[2]:.2f} P={Pb:.1f} T14={T14:.2f} h  chi2={np.sum(((y-mb)/ey)**2):.1f} N={len(y)}')
pickle.dump(dict(chain=ch, names=names, t=t, y=y, best=best, res=res), open(f'{S}/cand/{tic}/fit_free.pkl', 'wb'))
fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
x = (t - best[0]) * 24
ax[0].plot(x, y, '.', ms=2, color='0.6', label='TESS 2-min PDCSAP')
nb = 10; mm = len(x) // nb * nb
ax[0].plot(x[:mm].reshape(-1, nb).mean(1), y[:mm].reshape(-1, nb).mean(1), 'o', ms=4, color='k', label='20-min bins')
tt = np.linspace(t.min(), t.max(), 3000); mdl2 = batman.TransitModel(params, tt); model(best)
ax[0].plot((tt - best[0]) * 24, mdl2.light_curve(params), 'r-', lw=2, label='best-fit model')
ax[0].legend(fontsize=8); ax[0].set_ylabel('relative flux'); ax[0].set_title(f'TIC {tic} single transit (S{sec})')
ax[1].plot(x, y - mb, '.', ms=2, color='0.5'); ax[1].axhline(0, color='r'); ax[1].set_xlabel('hours from mid-transit'); ax[1].set_ylabel('resid')
plt.tight_layout(); plt.savefig(f'{S}/cand/{tic}/fit_free.png', dpi=90)
fig, ax = plt.subplots(1, 3, figsize=(12, 3.2))
ax[0].hist(rp, 60, color='C0'); ax[0].set_xlabel('Rp [R_Earth]')
ax[1].hist(np.log10(P), 60, color='C1'); ax[1].set_xlabel('log10 P [d] (circular, rho* prior)')
ax[2].hist(ch[:, 2], 60, color='C2'); ax[2].set_xlabel('impact parameter b')
plt.tight_layout(); plt.savefig(f'{S}/cand/{tic}/fit_free_post.png', dpi=90)
