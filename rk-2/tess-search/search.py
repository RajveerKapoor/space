"""Rebuild the transit-event screen from public SPOC, ExoFOP and NASA tables."""
import io
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
MAST = "https://archive.stsci.edu/missions/tess/catalogs/tce/"
EXOFOP = "https://exofop.ipac.caltech.edu/tess/"
NASA = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+tic_id,pl_name,pl_orbper+from+pscomppars&format=csv"
FIELDS = ("ticid tce_plnt_num tce_time0bt tce_period tce_duration tce_depth tce_prad tce_insol "
          "tce_num_transits tce_model_snr tce_robstat tce_impact tce_dicco_msky "
          "tce_dicco_msky_err tce_bin_oedp_stat").split()


def read(url):
    return urlopen(url, timeout=90).read()


def main():
    filenames = sorted(f for f, end in re.findall(r'href="(tess\d+-s\d+-s(\d+)_dvr-tcestats\.csv)"',
                                                   read(MAST).decode()) if 81 <= int(end) <= 108)
    def load(name):
        table = pd.read_csv(io.BytesIO(read(MAST+name)), comment="#", usecols=FIELDS)
        return table.assign(source=name)
    with ThreadPoolExecutor(6) as pool:
        rows = pd.concat(pool.map(load, filenames), ignore_index=True)
    screen = rows.query("0.15 < tce_insol < 1.65 and 0.4 < tce_prad < 2.4 and "
                        "tce_model_snr >= 7 and tce_num_transits >= 2 and tce_impact < 1 and "
                        "tce_depth > 0 and tce_robstat >= 5").sort_values("tce_model_snr", ascending=False)
    OUT.mkdir(exist_ok=True)
    screen.to_csv(OUT/"flux-radius-screen.csv", index=False)
    excluded = set()
    for catalog in ["toi","ctoi"]:
        data = pd.read_csv(io.BytesIO(read(f"{EXOFOP}download_{catalog}.php?sort={catalog}&output=csv")),
                           usecols=["TIC ID"])
        excluded.update(data["TIC ID"].dropna().astype(int))
    planets = pd.read_csv(io.BytesIO(read(NASA)))
    excluded.update(planets.tic_id.str.extract(r"(\d+)")[0].dropna().astype(int))
    shortlist = screen.loc[~screen.ticid.isin(excluded)].query(
        "tce_num_transits >= 3 and tce_bin_oedp_stat < 9 and tce_robstat >= 7 and "
        "tce_dicco_msky_err > 0 and tce_dicco_msky/tce_dicco_msky_err < 3 and tce_dicco_msky < 15")
    shortlist.to_csv(OUT/"unmatched-shortlist.csv", index=False)
    print(len(filenames), "releases,", len(rows), "TCE rows,", len(screen), "screened in,", len(shortlist), "shortlisted")


if __name__ == "__main__":
    main()
