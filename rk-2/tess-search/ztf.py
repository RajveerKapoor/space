"""Look for resolved observations of the TESS field during the predicted dips."""
import argparse
import io
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from astropy.utils import iers

from followup import EPOCH, PERIOD, DURATION
from photometry import OUT

DATA = Path(__file__).resolve().parent / "data" / "ztf-field.csv.gz"
API = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"
COLS = ["oid", "expid", "mjd", "mag", "magerr", "catflags", "filtercode"]
PALOMAR = EarthLocation.from_geodetic(-116.86*u.deg, 33.36*u.deg, 1700*u.m)


def fetch(stars):
    frames = []
    for star in stars.itertuples():
        params = {"POS": f"CIRCLE {star.RA_ICRS:.7f} {star.DE_ICRS:.7f} 0.00085",
                  "BAD_CATFLAGS_MASK": 32768, "FORMAT": "csv"}
        raw = urlopen(API + "?" + urlencode(params), timeout=150).read()
        rows = pd.read_csv(io.BytesIO(raw), usecols=COLS)
        frames.append(rows.assign(gaia_source=star.Source))
    DATA.parent.mkdir(exist_ok=True)
    pd.concat(frames, ignore_index=True).to_csv(DATA, index=False,
                                               compression={"method":"gzip","mtime":0})


def check():
    iers.conf.auto_download = False
    stars = pd.read_csv(OUT/"lead-gaia-field.csv", dtype={"Source":str}).set_index("Source")
    phot = pd.read_csv(DATA, dtype={"gaia_source":str})
    clean = phot.query("catflags == 0 and magerr > 0").dropna(subset=["mjd","mag","magerr"])
    clean = clean.sort_values("magerr").drop_duplicates(["gaia_source","expid","filtercode"])
    measurements = []
    for source, g in clean.groupby("gaia_source"):
        star = stars.loc[source]
        coord = SkyCoord(star.RA_ICRS,star.DE_ICRS,unit="deg")
        t = Time(g.mjd.to_numpy(),format="mjd",scale="utc",location=PALOMAR)
        bjd = (t.tdb + t.light_travel_time(coord)).jd
        cycle = np.rint((bjd-2457000-EPOCH)/PERIOD).astype(int)
        offset = bjd - (2457000+EPOCH+cycle*PERIOD)
        data = g.assign(bjd=bjd,cycle=cycle,offset=offset)
        for row in data[np.abs(offset) < DURATION/2].itertuples():
            outside = data[(data.oid==row.oid) & (np.abs(data.offset)>DURATION)]
            if len(outside) < 10:
                raise ValueError("Insufficient same-object, out-of-event baseline")
            baseline = outside.mag.median()
            measurements.append(dict(gaia_source=source, offset_arcsec=star.offset_arcsec,
                         clean_detections=len(data), expid=row.expid, cycle=row.cycle, bjd=row.bjd,
                         hours_from_center=row.offset*24, band=row.filtercode, mag=row.mag,
                         mag_err=row.magerr, baseline_n=len(outside), baseline_mag=baseline,
                         delta_mag=row.mag-baseline,
                         baseline_mad_mag=1.4826*(outside.mag-baseline).abs().median()))
    result = pd.DataFrame(measurements).sort_values("offset_arcsec")
    result.to_csv(OUT/"lead-ztf-in-window.csv", index=False,float_format="%.9f")
    print(result[["gaia_source","clean_detections","hours_from_center","delta_mag","mag_err",
                  "baseline_mad_mag"]].to_string(index=False))


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true", help="Replace the bundled IRSA photometry")
    if parser.parse_args().refresh:
        fetch(pd.read_csv(OUT/"lead-gaia-field.csv", dtype={"Source":str}))
    check()
