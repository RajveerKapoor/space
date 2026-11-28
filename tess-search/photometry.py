"""Measure predicted dips in SPOC and full-frame light curves; keep extractions separate."""
import argparse
import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from tempfile import NamedTemporaryFile
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.io import fits

from search import OUT

# Full-frame follow-up after inspection of the short-cadence screen.
FOLLOWUP = {230091879, 125421413, 268319776, 198008979, 150392377}
S3 = "https://stpubdata.s3.amazonaws.com/tess/public/tid"
FP = re.compile(r"tess\d+-s\d+-\d{16}-\d{4}-s")


def observations(targets):
    req = dict(service="Mast.Caom.Filtered", format="json", pagesize=20000, page=1,
               params=dict(columns="obs_collection,obs_id,target_name,t_min,t_max,sequence_number,dataURL",
                           filters=[dict(paramName="target_name", values=list(map(str, targets))),
                                    dict(paramName="dataproduct_type", values=["timeseries"])]))
    url = "https://mast.stsci.edu/api/v0/invoke?" + urlencode(dict(request=json.dumps(req)))
    reply = json.load(urlopen(url, timeout=90))
    if reply["status"] != "COMPLETE" or len(reply["data"]) != reply["paging"]["rowsFiltered"]:
        raise RuntimeError("Incomplete MAST result")
    return reply["data"]


def product(item, tic):
    obs_id, uri = item["obs_id"], item.get("dataURL") or ""
    sector = item["sequence_number"]
    if item["obs_collection"] == "TESS" and FP.fullmatch(obs_id):
        digits = f"{tic:016d}"
        path = "/".join(digits[i:i+4] for i in range(0,16,4))
        return "SPOC", f"{S3}/s{sector:04}/{path}/{obs_id}_lc.fits"
    if tic in FOLLOWUP and uri.startswith(("mast:HLSP/qlp/","mast:HLSP/tess-spoc/")):
        path = uri.removeprefix("mast:HLSP/")
        if path.endswith(("_lc.fits","_llc.fits")):
            return path.split("/")[0], "https://archive.stsci.edu/hlsps/"+path
    return None


def epochs(time_min, time_max, epoch, period):
    first = int(np.ceil((time_min-epoch)/period))
    last = int(np.floor((time_max-epoch)/period))
    return epoch + period*np.arange(first, last+1)


def read_curve(url, cache):
    path = cache/url.rsplit("/",1)[1]
    if not path.exists():
        raw = urlopen(url, timeout=120).read()
        with NamedTemporaryFile(dir=cache, delete=False) as output:
            output.write(raw)
        Path(output.name).replace(path)
    with fits.open(path) as hdus:
        data = hdus[1].data
        col = next(c for c in ("PDCSAP_FLUX","KSPSAP_FLUX","DET_FLUX") if c in data.names)
        mask = (data["QUALITY"] == 0) & np.isfinite(data["TIME"]) & np.isfinite(data[col])
        t = data["TIME"][mask].astype(float)
        corrected = data[col][mask].astype(float)
        raw = data["SAP_FLUX"][mask].astype(float)
    return t, corrected/np.median(corrected)*1e6, raw/np.nanmedian(raw)*1e6


def depth(time, values, center, duration):
    dt = (time-center)/duration
    inside = np.abs(dt) < .38
    left = (dt > -2) & (dt < -.75)
    right = (dt < 2) & (dt > .75)
    count = .70*duration/np.median(np.diff(time))
    if inside.sum() < .75*count or min(left.sum(),right.sum()) < count:
        return np.nan
    side = left | right
    slope, intercept = np.polyfit(dt[side],values[side],1)
    return intercept-np.mean(values[inside]-slope*dt[inside])


def measure(task, cache):
    row, sector, pipeline, url, centers = task
    time, corrected, raw = read_curve(url, cache)
    dur, period, epoch = row.tce_duration/24, row.tce_period, row.tce_time0bt
    controls = np.arange(time.min()+2*dur,time.max()-2*dur,.8*dur)
    phases = (controls-epoch+period/2)%period-period/2
    blanks = np.array([depth(time,corrected,x,dur) for x in controls[np.abs(phases)>3*dur]])
    blanks = blanks[np.isfinite(blanks)]
    if len(blanks) < 12:
        return []
    sigma = 1.4826*np.median(np.abs(blanks-np.median(blanks)))
    if sigma <= 0:
        return []
    events = []
    for center in centers:
        dip = depth(time,corrected,center,dur)
        if np.isfinite(dip):
            events.append(dict(tic=row.ticid, period_d=period, btjd=center, sector=sector,
                               pipeline=pipeline, depth_ppm=dip, sap_depth_ppm=depth(time,raw,center,dur),
                               noise_ppm=sigma, local_snr=dip/sigma, file=url.rsplit("/",1)[1]))
    return events


def main(cache):
    cache.mkdir(parents=True,exist_ok=True)
    rows = pd.read_csv(OUT/"unmatched-shortlist.csv")
    index = observations(rows.ticid.unique())
    tasks = []
    for row in rows.itertuples():
        for item in index:
            if item["target_name"] != str(row.ticid):
                continue
            selected = product(item,row.ticid)
            if not selected:
                continue
            # MJD to BJD - 2457000, adequate for selecting data products only.
            times = epochs(item["t_min"]-56999.5,item["t_max"]-56999.5,row.tce_time0bt,row.tce_period)
            if len(times):
                tasks.append((row,item["sequence_number"],*selected,times))
    with ThreadPoolExecutor(6) as pool:
        events = [item for batch in pool.map(lambda task: measure(task,cache),tasks) for item in batch]
    result = pd.DataFrame(events).sort_values(["tic","period_d","btjd","pipeline"])
    result.to_csv(OUT/"event-photometry.csv",index=False,float_format="%.9f")
    print(len(result),"measurements of",len(result[["tic","period_d","btjd"]].drop_duplicates()),"distinct events")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache",type=Path,required=True)
    main(parser.parse_args().cache)
