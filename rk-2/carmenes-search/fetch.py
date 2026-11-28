"""Download CARMENES DR1 observations and stellar parameters from the GAVO TAP."""
import csv
import gzip
import io
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
TAP = "https://dc.g-vo.org/tap/sync"
QUERIES = {
    "velocities": """SELECT object, bjd, a_rv, err_a_rv, spec_flag, halpha, err_halpha,
        dlw, err_dlw, ccf_bis, err_ccf_bis, berv, ccf_rv, err_ccf_rv
        FROM carmenes.rvs ORDER BY object, bjd""",
    "stars": """SELECT carmenes_id, name, survey, remarks, nobs, mass, err_mass,
        radius, err_radius, teff, err_teff, p_rot, ra, dec
        FROM carmenes.objects ORDER BY carmenes_id""",
}

def fetch(label, query):
    params = urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv",
                            MAXREC=100000, QUERY=" ".join(query.split())))
    data = urlopen(f"{TAP}?{params}", timeout=150).read()
    rows = list(csv.reader(io.StringIO(data.decode())))
    if not rows or rows[0][0] not in {"object", "carmenes_id"} or len(rows) >= 100001:
        raise ValueError("Invalid or truncated TAP response")
    path = ROOT / "data" / f"{label}.csv.gz"
    path.write_bytes(gzip.compress(data, mtime=0))
    print(label, len(rows)-1, flush=True)

if __name__ == "__main__":
    for label, query in QUERIES.items():
        fetch(label, query)
