"""Compare the blind recovery with published solutions and assess the irradiation scale."""
import argparse
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd

from survey import ROOT, OUT, nightly

TARGET = "J02530+168"
ARCHIVE = ROOT / "data/teegarden-c-archive.csv"


def refresh():
    fields = ("pl_name,default_flag,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_orbsmax,"
              "pl_orbsmaxerr1,pl_orbsmaxerr2,pl_msinie,pl_msinieerr1,pl_msinieerr2,pl_insol,"
              "pl_insolerr1,pl_insolerr2,st_teff,st_rad,st_mass,st_lum,disc_year,"
              "discoverymethod,tran_flag,rv_flag,pl_refname")
    query = f"select {fields} from ps where pl_name = 'Teegarden''s Star c'"
    url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?" + urlencode(dict(query=query, format="csv"))
    raw = urlopen(url, timeout=90).read()
    if not raw.startswith(b"pl_name,") or b"Teegarden" not in raw:
        raise ValueError("Unexpected NASA Archive response")
    ARCHIVE.write_bytes(raw)


def terms(time, period):
    phase = 2*np.pi*time/period
    return np.column_stack([np.cos(phase), np.sin(phase)])


def fit_at_period(data, period, companion_periods):
    """Use a separate offset per run and two free amplitudes per fixed frequency."""
    t = data.bjd.to_numpy()-2458000
    y = data.rv.to_numpy()
    err = data.error.to_numpy()
    null = np.column_stack([pd.get_dummies(data.instrument).to_numpy(float),
                            *(terms(t,p) for p in companion_periods)])
    def fit(matrix):
        beta = np.linalg.lstsq(matrix/err[:,None], y/err, rcond=None)[0]
        return beta, np.sum(((y-matrix@beta)/err)**2)
    _, null_chi2 = fit(null)
    beta, full_chi2 = fit(np.column_stack([null, terms(t,period)]))
    return np.hypot(*beta[-2:]), null_chi2-full_chi2


def main():
    published = pd.read_csv(ARCHIVE).query("default_flag == 1").iloc[0]
    peaks = pd.read_csv(OUT/"blind-peaks.csv")
    peaks = peaks[peaks.target == TARGET]
    blind = peaks.loc[(peaks.period_d-published.pl_orbper).abs().idxmin()]
    catalog = pd.read_csv(ROOT/"data/archive-cones.csv.gz").set_index("pl_name")
    companion_periods = [catalog.loc[f"Teegarden's Star {p}","pl_orbper"] for p in ("b","d")]
    discovery = pd.read_csv(ROOT/"data/velocities.csv.gz")
    cutoff = discovery.loc[discovery.object==TARGET,"bjd"].max()
    follow = pd.read_csv(ROOT/"data/teegarden-followup.csv.gz")
    # The two MAROON-X arms observed simultaneously. Keep the red arm.
    follow = follow[~follow.instrument.str.startswith("MAROON-X_b")]
    follow = follow[(follow.instrument!="CARM_VIS") | (follow.bjd>cutoff)]
    nights = []
    for inst, group in follow.groupby("instrument"):
        t, rv, err = nightly(group.assign(spec_flag=0),"rv_ms","err_ms",floor=1.5)
        nights.append(pd.DataFrame(dict(instrument=inst,bjd=t,rv=rv,error=err)))
    holdout = pd.concat(nights, ignore_index=True)
    rows = []
    for sample, data in [("other_instruments", holdout[holdout.instrument!="CARM_VIS"]),
                         ("other_plus_later_CARM", holdout)]:
        for source, period in [("blind",blind.period_d),("published",published.pl_orbper)]:
            amplitude, delta = fit_at_period(data,period,companion_periods)
            rows.append(dict(sample=sample,period_source=source,n_instrument_nights=len(data),
                             period_d=period, semiamplitude_ms=amplitude, delta_chi2=delta))
    checks = pd.DataFrame(rows)
    checks.to_csv(OUT/"teegarden-fixed-period-check.csv",index=False,float_format="%.9g")
    print(checks.to_string(index=False))
    # Kopparapu et al. 2014, one-Earth-mass runaway/maximum greenhouse boundaries.
    dt = published.st_teff-5780
    lum = 10**published.st_lum
    for name, coefficients in [
        ("inner", [1.107,1.332e-4,1.580e-8,-8.308e-12,-1.931e-15]),
        ("outer", [.356,6.171e-5,1.698e-9,-3.198e-12,-5.575e-16])
    ]:
        flux = np.polynomial.polynomial.polyval(dt,coefficients)
        print(f"{name} boundary {flux:.4f} Earth flux, {np.sqrt(lum/flux):.5f} AU")
    print(f"Published luminosity {lum:.7f} Suns; irradiation {lum/published.pl_orbsmax**2:.4f} Earth flux")
    print(f"Last date in blind RV archive {cutoff:.6f} BJD")


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--refresh",action="store_true")
    if parser.parse_args().refresh:
        refresh()
    main()
