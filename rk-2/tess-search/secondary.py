"""Remeasure the secondary epoch flagged in the discovery report using QLP flux."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from followup import EPOCH, PERIOD, DURATION, TIC
from photometry import OUT, observations, product, epochs, measure

OFFSET = 35.778  # Days after primary center in the SPOC DV summary.


def main(cache):
    peak = pd.Series(dict(ticid=TIC,tce_period=PERIOD,tce_time0bt=EPOCH+OFFSET,tce_duration=DURATION*24))
    events = []
    for item in observations([TIC]):
        selected = product(item,TIC)
        if not selected or selected[0] != "qlp":
            continue
        centers=epochs(item["t_min"]-56999.5,item["t_max"]-56999.5,peak.tce_time0bt,PERIOD)
        if len(centers):
            events.extend(measure((peak,item["sequence_number"],*selected,centers),cache))
    result=pd.DataFrame(events).sort_values("btjd")
    result.to_csv(OUT/"lead-secondary.csv",index=False,float_format="%.9f")
    w=1/result.noise_ppm**2
    print("Fixed secondary phase:",len(result),"events, depth and characteristic noise (ppm):",
          np.average(result.depth_ppm,weights=w),1/np.sqrt(w.sum()))
    print(result[["sector","depth_ppm","noise_ppm"]].to_string(index=False))


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--cache",type=Path,required=True)
    main(parser.parse_args().cache)
