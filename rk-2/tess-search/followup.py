"""Check coverage of the discovery fit and test its prediction in full-frame data."""
import argparse
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits

from photometry import depth, read_curve, OUT

TIC = 230091879
# Full precision from the original SPOC one-page Data Validation Summary.
EPOCH, PERIOD, DURATION = 1747.0873, 133.09719, 5.50008/24
DV = "tess2019199201929-s0014-s0086-0000000230091879-00971_dvt.fits"


def main(cache):
    path = cache/DV
    if not path.exists():
        query = urlencode({"uri":"mast:TESS/product/"+DV})
        path.write_bytes(urlopen("https://mast.stsci.edu/api/v0.1/Download/file?"+query, timeout=180).read())
    with fits.open(path) as hdus:
        d = hdus[1].data
        valid = d["TIME"][np.isfinite(d["LC_INIT"])]
    phot = pd.read_csv(OUT/"event-photometry.csv")
    target = phot[phot.tic == TIC].copy()
    rows = []
    for epoch in sorted(target.btjd.unique()):
        k = round((epoch-EPOCH)/PERIOD)
        center = EPOCH+k*PERIOD
        count = np.sum(np.abs(valid-center) < DURATION/2)
        rows.append(dict(cycle=k,btjd=center,sector=target[target.btjd==epoch].sector.iloc[0],
                         original_dv_finite_cadences=count))
    coverage = pd.DataFrame(rows)
    coverage.to_csv(OUT/"lead-coverage.csv",index=False)
    # QLP and TESS-SPOC reductions of one exposure share photons. Use QLP alone
    # for the test, and keep the other extraction as a consistency check.
    tests, controls, weights, profiles = [], [], [], []
    for row in coverage[coverage.original_dv_finite_cadences==0].itertuples():
        files = target[(target.sector==row.sector)&(target.pipeline=="qlp")].file
        filename = files.item()
        tstr = f"{TIC:016d}"
        url = ("https://archive.stsci.edu/hlsps/qlp/"+f"s{row.sector:04}/"+
               "/".join(tstr[i:i+4] for i in range(0,16,4))+"/"+filename)
        time,flux,raw = read_curve(url,cache)
        center = row.btjd
        trial = np.arange(time.min()+2*DURATION,time.max()-2*DURATION,.8*DURATION)
        trial = trial[np.abs(trial-center)>3*DURATION]
        offsets = np.array([depth(time,flux,c,DURATION) for c in trial])
        offsets = offsets[np.isfinite(offsets)]
        error = 1.4826*np.median(np.abs(offsets-np.median(offsets)))
        val = depth(time,flux,center,DURATION)
        tests.append(dict(sector=row.sector,cycle=row.cycle,btjd=center,qlp_depth_ppm=val,
                          raw_depth_ppm=depth(time,raw,center,DURATION),noise_ppm=error,
                          control_windows=len(offsets)))
        controls.append(offsets);weights.append(1/error**2)
        hours = (time-center)*24
        nearby = np.abs(hours)<2*DURATION*24
        side = nearby & (np.abs(hours)>.75*DURATION*24)
        slope,intercept = np.polyfit(hours[side],flux[side],1)
        profiles.append((row.sector,hours[nearby],(flux-(intercept+slope*hours))[nearby]))
    report = pd.DataFrame(tests)
    report.to_csv(OUT/"lead-extra-epochs.csv",index=False,float_format="%.9f")
    w = np.array(weights); combined = np.average(report.qlp_depth_ppm,weights=w)
    # Pair every control window in one sector with every control in the other.
    # Controls overlap. This rank is descriptive, not a survey false-alarm rate.
    if len(controls)!=2:
        raise RuntimeError("Expected the two coverage gaps in the archived discovery product")
    null=(w[0]*controls[0][:,None]+w[1]*controls[1][None,:])/w.sum()
    print(report.to_string(index=False))
    print("weighted depth, characteristic noise, control pairs, >= observed",
          combined,1/np.sqrt(w.sum()),null.size,(null>=combined).sum())

    fig,axes=plt.subplots(1,2,figsize=(10,4),sharey=True,layout="constrained")
    for ax,(sector,x,y) in zip(axes,profiles):
        half=DURATION*12
        ax.axvspan(-half,half,color="#ddaa76",alpha=.2,label="Predicted transit window")
        ax.axhline(0,color="#aaa",lw=.7)
        centers=np.arange(-10,10.01,.5)
        bins=[np.mean(y[np.abs(x-c)<.25]) if np.any(np.abs(x-c)<.25) else np.nan for c in centers]
        ax.plot(centers,np.array(bins)/1000,"o-",ms=3,lw=.8,color="#275e69")
        ax.set(xlabel="Hours from predicted center",title=f"TESS sector {sector}",xlim=(-10,10))
        ax.spines[["top","right"]].set_visible(False)
    axes[0].set_ylabel("Relative flux (parts per thousand)")
    fig.suptitle("TIC 230091879 · two predicted epochs outside the discovery fit")
    axes[1].legend(frameon=False,fontsize=8)
    fig.savefig(OUT/"lead-extra-epochs.png",dpi=165)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--cache",type=Path,required=True)
    main(parser.parse_args().cache)
