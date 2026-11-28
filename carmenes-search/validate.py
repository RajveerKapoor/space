"""Post-search catalog matches, activity checks and an external instrument test."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import f as f_dist

from survey import ROOT, OUT, nightly, fit, design

TARGET = "J02530+168"


def improvement(y, err, null, full):
    def solve(matrix):
        beta = np.linalg.lstsq(matrix/err[:, None], y/err, rcond=None)[0]
        return beta, np.sum(((y-matrix@beta)/err)**2)
    _, chi0 = solve(null)
    beta, chi1 = solve(full)
    n, k = full.shape
    extra = k-null.shape[1]
    return beta, chi0-chi1, f_dist.sf(((chi0-chi1)/extra)/(chi1/(n-k)), extra, n-k)


def main():
    rv = pd.read_csv(ROOT / "data/velocities.csv.gz")
    peaks = pd.read_csv(OUT / "blind-peaks.csv")
    stars = pd.read_csv(ROOT / "data/stars.csv.gz").set_index("carmenes_id")
    catalog = pd.read_csv(ROOT / "data/archive-cones.csv.gz")
    records = []
    for peak in peaks.query("0.2 <= flux_earth <= 1.6 and msini_earth < 10").itertuples():
        source = rv[rv.object == peak.target]
        star = stars.loc[peak.target]
        dra = (catalog.ra-star.ra)*np.cos(np.deg2rad(star.dec))
        neighbors = catalog[np.hypot(dra,catalog.dec-star.dec) < 0.1]
        matches = neighbors[np.abs(neighbors.pl_orbper/peak.period_d-1) < 0.015]
        rec = dict(target=peak.target, period_d=peak.period_d,
                   catalog_match=" | ".join(matches.pl_name))
        for indicator in ("halpha", "dlw", "ccf_bis"):
            t,y,e = nightly(source, indicator, "err_"+indicator, floor=0)
            t -= t.mean()
            rec[indicator+"_p"] = improvement(y,e,design(t,[]),design(t,[1/peak.period_d]))[2]
        records.append(rec)
    pd.DataFrame(records).to_csv(OUT/"temperate-vetting.csv", index=False, float_format="%.8g")

    group = rv[rv.object == TARGET]
    dates, y, err = nightly(group)
    epoch = dates.mean()
    t = dates-epoch
    freqs = 1/peaks[peaks.target == TARGET].period_d.to_numpy()
    coeff, res = fit(t,y,err,freqs)
    planet = design(t,freqs)[:,4:6] @ coeff[4:6]
    phase = (t*freqs[1]) % 1
    phase_data = pd.DataFrame(dict(bjd=dates,phase=phase,rv_minus_nuisance_ms=res+planet,error_ms=err))
    phase_data.to_csv(OUT/"teegarden-phase.csv",index=False,float_format="%.10f")

    # Freeze the frequency, phase and amplitude from DR1. Independent data use
    # instrument offsets, plus a free sinusoid at the inner planet's DR1 period.
    # Keep the red MAROON-X arm only, avoiding duplicate simultaneous spectra.
    other = pd.read_csv(ROOT/"data/teegarden-followup.csv.gz")
    other = other[~other.instrument.str.startswith(("CARM","MAROON-X_b"))].copy()
    nights = []
    for inst, rows in other.groupby("instrument"):
        rows = rows.assign(spec_flag=0)
        times, vals, errs = nightly(rows,"rv_ms","err_ms",floor=1.5)
        nights.extend(zip([inst]*len(times), times, vals, errs))
    oth = pd.DataFrame(nights,columns=["inst","bjd","y","err"])
    dt = oth.bjd.to_numpy()-epoch
    nuisance = np.column_stack([pd.get_dummies(oth.inst).to_numpy(float),
                                design(dt, freqs[:1])[:,2:]])
    pred = design(dt, freqs[1:2])[:,2:] @ coeff[4:6]
    YY, ee = oth.y.to_numpy(), oth.err.to_numpy()
    beta, delta, p = improvement(YY,ee,nuisance,np.column_stack([nuisance,pred]))
    # Amplitude fit relative to the DR1 waveform. Also test the actual frozen
    # amplitude, allowing exactly the same nuisance terms in both models.
    def chisq(values):
        nuisance_fit = np.linalg.lstsq(nuisance/ee[:,None],values/ee,rcond=None)[0]
        return np.sum(((values-nuisance@nuisance_fit)/ee)**2)
    print("Independent instrument nights",oth.inst.value_counts().to_dict())
    print("Template scale, improvement at fitted scale, nominal fixed-phase p, improvement at scale one",
          beta[-1],delta,p,chisq(YY)-chisq(YY-pred))
    print("DR1 nights, span, epoch",len(dates),np.ptp(t),epoch)
    fig,ax=plt.subplots(figsize=(8,4.5),layout="constrained")
    ax.errorbar(phase,res+planet,err,fmt=".",ms=3,color="#6994a0",ecolor="#cad7da",alpha=.7)
    centers=[];means=[];errors=[]
    for i in range(12):
        mask=(phase>=i/12)&(phase<(i+1)/12)
        weights=1/err[mask]**2
        centers.append(np.average(phase[mask],weights=weights))
        means.append(np.average((res+planet)[mask],weights=weights))
        errors.append(1/np.sqrt(weights.sum()))
    ax.errorbar(centers,means,errors,fmt="o",ms=5,color="#163d47",label="12 weighted phase bins")
    ph=np.linspace(0,1,400)
    ax.plot(ph,coeff[4]*np.cos(2*np.pi*ph)+coeff[5]*np.sin(2*np.pi*ph),
            color="#c6532c",lw=2,label="Circular fit")
    ax.set(xlim=(0,1),xlabel="Orbital phase",ylabel="Relative radial velocity (m/s)",
           title="Teegarden’s Star c  ·  233 observing nights folded at 11.405 days")
    ax.spines[["top","right"]].set_visible(False)
    ax.legend(frameon=False,loc="upper right")
    fig.savefig(OUT/"teegarden-phase.png",dpi=160)
    pd.DataFrame([dict(n_independent_nights=len(dt),scale=beta[-1],delta_chi2_fitted=delta,
                       nominal_p=p,delta_chi2_fixed=chisq(YY)-chisq(YY-pred))]).to_csv(
                       OUT/"independent-test.csv",index=False,float_format="%.10g")


if __name__ == "__main__":
    main()
