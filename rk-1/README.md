# RK-1

**RK-1 (TIC 178994149.01)** is a single-transit giant-planet candidate orbiting a nearby M dwarf. It was found in public TESS data from Sector 98 and does not appear in any existing planet or candidate catalog.

Author: Rajveer Kapoor (independent researcher). This is a draft and has not been peer reviewed.

## The candidate

| Property | Value |
|---|---|
| Host star | TIC 178994149 = 2MASS J04390289−2919401 = Gaia DR3 4879180296770956416 |
| Host type | M3.5 V dwarf, 0.48 R☉, 116.1 pc |
| Mid-transit time | BJD_TDB 2460998.9041 ± 0.0020 (2025-11-19 09:41 UTC) |
| Depth | 26,500 (+2,200 / −1,800) ppm |
| Duration (T14) | 2.53 (+0.24 / −0.30) h |
| Radius (non-grazing) | 8.2 (+0.5 / −0.4) R⊕ = 0.73 R_J |
| Orbital period | > 49.6 d (one transit in seven TESS sectors, 2018–2026) |
| Signal-to-noise | ≈ 17 |

The signal is on the target star according to the difference image, shows no background or centroid change, and passes every vetting test applied (see the paper, Section 5). The transit is short for a long orbit, so either the orbit is eccentric (e ≳ 0.5) or the transit is grazing.

## What is still needed

This is a planet candidate, not a confirmed planet:

- **Period:** only one transit has been seen. A second transit from future TESS sectors or ground photometry would fix the period.
- **Mass:** radial velocities are needed to rule out a brown dwarf or low-mass star. A Saturn-mass planet at P = 100 d gives about 21 m/s, a Jupiter-mass planet about 71 m/s.
- **Host star:** the star is 1.1–1.4 mag brighter than the main sequence for its colour, which suggests youth or an unresolved companion. High-resolution imaging and a spectrum would test this.

## Repository contents

| Path | Contents |
|---|---|
| `paper/TIC178994149_paper.pdf` | Short formal paper (7 pages) |
| `paper/paper.html`, `paper/fig*.png` | Paper source and figures |
| `report/TIC178994149_discovery_report.pdf` | Illustrated discovery report |
| `report/index.html`, `report/fig_*.png` | Report source and figures |
| `code/` | Search, vetting and fitting pipeline (Python) |
| `data/tess2025312202959-s0098-…-s_lc.fits` | SPOC 2-min light curve containing the transit |
| `data/fit_free_density.pkl` | MCMC posterior, free stellar density (used in the paper) |
| `data/fit_density_prior.pkl` | MCMC posterior, TIC stellar-density prior |
| `data/period_window.pkl` | Periods allowed by TESS coverage |
| `data/ffi_photometry.pkl` | Aperture photometry from TESScut full-frame-image cutouts |

## How the search worked

1. `code/search.py` streams every SPOC 2-min light curve for TESS Sectors 97–108 (153,730 light curves of 84,380 stars) from the public MAST bucket on AWS, detrends each one, and runs a single-event box search (1–13 h).
2. `code/analyze.py` and `code/vet1.py` cross-match events against the TOI, CTOI, confirmed-planet and SPOC TCE catalogs and remove systematics.
3. `code/bkgtest.py`, `code/locate.py` and `code/pixvet.py` run the background, centroid and difference-image tests.
4. `code/bls_multi.py`, `code/bls_rank.py` and `code/bls_shape.py` run the multi-sector box-least-squares search over 11,496 stars.
5. `code/deepdive.py`, `code/ffiphot.py` and `code/fit_single*.py` do the per-candidate analysis and MCMC transit fit.

The scripts expect to be run from `code/`, with intermediate products in `../work/`. They need `numpy`, `scipy`, `astropy`, `astroquery`, `wotan`, `batman-package`, `emcee` and `matplotlib`. The catalogs (TOI, CTOI, NASA Exoplanet Archive, SPOC TCE lists) are not included; they are public on ExoFOP, the NASA Exoplanet Archive and MAST.

## Credits

This work uses data from the TESS mission (NASA), obtained from MAST at STScI, together with ExoFOP, the NASA Exoplanet Archive, and ESA's Gaia mission. The data search, analysis and drafting were carried out with the assistance of Claude, an AI model developed by Anthropic.
