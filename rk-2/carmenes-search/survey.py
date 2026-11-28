"""Blind Doppler search. Only uses DR1 velocities, quality flags and stellar data."""
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"


def nightly(group, value="a_rv", error="err_a_rv", floor=1.5):
    g = group.loc[(group.spec_flag == 0) & group[value].notna() &
                  (group[error] > 0), ["bjd", value, error]].copy()
    # Check gross outliers against the target's own robust scatter.
    median = g[value].median()
    mad = 1.4826 * (g[value] - median).abs().median()
    g = g.loc[(g[value] - median).abs() <= 8 * np.hypot(mad, g[error])]
    w = 1 / g[error]**2
    weighted = pd.DataFrame(dict(w=w, t=g.bjd*w, y=g[value]*w))
    nights = weighted.groupby(np.floor(g.bjd)).sum()
    return ((nights.t/nights.w).to_numpy(), (nights.y/nights.w).to_numpy(),
            np.sqrt(1/nights.w + floor**2).to_numpy())


def design(t, freqs):
    cols = [np.ones_like(t), t / 365.25]
    for f in freqs:
        cols.extend([np.cos(2*np.pi*f*t), np.sin(2*np.pi*f*t)])
    return np.column_stack(cols)


def fit(t, y, err, freqs):
    matrix = design(t, freqs)
    coef = np.linalg.lstsq(matrix / err[:, None], y/err, rcond=None)[0]
    return coef, y - matrix @ coef


def refine(t, y, err, freqs, width):
    # Coordinate descent, with simultaneous linear amplitudes and a drift.
    freqs = list(freqs)
    for _ in range(3):
        for k, f in enumerate(freqs):
            def chi2(trial):
                trial_f = freqs.copy()
                trial_f[k] = trial
                return np.sum((fit(t, y, err, trial_f)[1]/err)**2)
            freqs[k] = minimize_scalar(chi2, method="bounded", bounds=(f-width, f+width),
                                       options={"xatol": 1e-10}).x
    return freqs


def run():
    stars = pd.read_csv(ROOT / "data/stars.csv.gz").set_index("carmenes_id")
    rv = pd.read_csv(ROOT / "data/velocities.csv.gz")
    candidates, ledger = [], []
    for target, obs in rv.groupby("object"):
        if target not in stars.index:
            continue
        star = stars.loc[target]
        if not star.survey.strip().startswith("GTO") or pd.notna(star.remarks):
            continue
        dates, y, err = nightly(obs)
        if len(dates) < 40 or np.ptp(dates) < 180:
            continue
        t = dates - dates.mean()
        maxperiod = min(200, np.ptp(t)/3)
        freqs, stages = [], []
        for stage in range(1, 5):
            residual = fit(t, y, err, freqs)[1]
            ls = LombScargle(t, residual, err)
            grid, power = ls.autopower(minimum_frequency=1/maxperiod, maximum_frequency=0.5,
                                       samples_per_peak=10)
            peak = np.argmax(power)
            fap = ls.false_alarm_probability(power[peak], method="baluev",
                                             minimum_frequency=grid[0], maximum_frequency=grid[-1])
            if fap >= 0.01:
                break
            freqs = refine(t, y, err, [*freqs, grid[peak]], 0.2/np.ptp(t))
            stages.append((stage, power[peak], fap))
        ledger.append(dict(target=target, nights=len(t), span_d=np.ptp(t), peaks=len(freqs)))
        coef, _ = fit(t, y, err, freqs)
        for i, (stage, power, fap) in enumerate(stages):
            period = 1/freqs[i]
            amp = np.hypot(*coef[2+2*i:4+2*i])
            axis = (star.mass*(period/365.25)**2)**(1/3)
            flux = star.radius**2*(star.teff/5772)**4/axis**2
            msini = amp/0.08946*(period/365.25)**(1/3)*star.mass**(2/3)
            candidates.append(dict(target=target, stage=stage, period_d=period, k_ms=amp,
                                   n=len(t), power=power, screen_fap=fap, axis_au=axis,
                                   flux_earth=flux, msini_earth=msini, rotation_d=star.p_rot))
    OUT.mkdir(exist_ok=True)
    pd.DataFrame(ledger).to_csv(OUT / "search-ledger.csv", index=False, float_format="%.9g")
    results = pd.DataFrame(candidates)
    results.to_csv(OUT / "blind-peaks.csv", index=False, float_format="%.10g")
    shortlist = results.query("0.2 <= flux_earth <= 1.6 and msini_earth < 10").sort_values("screen_fap")
    print(f"{len(ledger)} stars, {len(results)} peaks; exploratory temperate shortlist:")
    print(shortlist.to_string(index=False))


if __name__ == "__main__":
    run()
