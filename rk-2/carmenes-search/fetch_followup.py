"""Post-detection reference tables. Run only after survey.py."""
import csv
import gzip
import io
from urllib.parse import urlencode
from urllib.request import urlopen

from fetch import ROOT


def download(url, params):
    return urlopen(url + "?" + urlencode(params), timeout=120).read()


def save(name, raw):
    (ROOT / "data" / name).write_bytes(gzip.compress(raw, mtime=0))


if __name__ == "__main__":
    # Published velocities from the primary and follow-up instruments.
    raw = download("https://vizier.cfa.harvard.edu/viz-bin/asu-tsv",
                   {"-source": "J/A+A/684/A117/tablea1", "-out.max": 10000})
    table = [row for row in csv.reader(io.StringIO(raw.decode()), delimiter="\t")
             if len(row) == 4 and row[1].strip().replace(".", "", 1).isdigit()]
    if len(table) < 300:
        raise ValueError("Incomplete VizieR result")
    stream = io.StringIO()
    writer = csv.writer(stream)
    writer.writerow(["instrument", "bjd", "rv_ms", "err_ms"])
    writer.writerows([[v.strip() for v in row] for row in table])
    save("teegarden-followup.csv.gz", stream.getvalue().encode())

    rows = list(csv.DictReader((ROOT / "results/blind-peaks.csv").open()))
    targets = {r["target"] for r in rows if 0.2 <= float(r["flux_earth"]) <= 1.6
               and float(r["msini_earth"]) < 10}
    with gzip.open(ROOT / "data/stars.csv.gz", "rt") as f:
        stars = [r for r in csv.DictReader(f) if r["carmenes_id"] in targets]
    # Wide positional cones allow for high proper motion and coordinate epochs.
    cones = [f"1=contains(point('ICRS',ra,dec),circle('ICRS',{r['ra']},{r['dec']},0.10))"
             for r in stars]
    fields = "pl_name,hostname,ra,dec,pl_orbper,pl_orbsmax,pl_bmasse,pl_bmassprov,sy_dist,disc_year"
    q = f"select {fields} from pscomppars where " + " or ".join(cones)
    raw = download("https://exoplanetarchive.ipac.caltech.edu/TAP/sync", {"query":q,"format":"csv"})
    if not raw.startswith(b"pl_name,"):
        raise ValueError("Invalid NASA response")
    save("archive-cones.csv.gz", raw)
