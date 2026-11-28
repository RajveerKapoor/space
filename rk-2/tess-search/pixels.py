"""Inspect independent full-frame cutouts at the two extra epochs in detector coordinates."""
import argparse
import io
import zipfile
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

from followup import DURATION
from photometry import OUT


def contrast(times, cube, center):
    t=(times-center)/DURATION
    inside=np.abs(t)<.38
    left=(t > -2)&(t < -.75)
    right=(t > .75)&(t < 2)
    count=.70*DURATION/np.median(np.diff(times))
    if inside.sum()<.75*count or min(left.sum(),right.sum())<count:
        return None
    side=left|right
    t0=t[side].mean()
    slope=(t[side]-t0) @ cube[side] / np.sum((t[side]-t0)**2)
    return cube[side].mean(axis=0)+(t[inside].mean()-t0)*slope-cube[inside].mean(axis=0)


def main(cache):
    gaia=pd.read_csv(OUT/"lead-gaia-field.csv",dtype={"Source":"string"})
    tests=pd.read_csv(OUT/"lead-extra-epochs.csv")
    fig,axes=plt.subplots(1,2,figsize=(8.5,4.1),layout="constrained")
    for ax,row in zip(axes,tests.itertuples()):
        path=cache/f"pixels-{row.sector}.fits"
        if not path.exists():
            url="https://mast.stsci.edu/tesscut/api/v0.1/astrocut?"+urlencode(
                dict(ra=247.856,dec=68.725,x=11,y=11,units="px",sector=row.sector))
            with zipfile.ZipFile(io.BytesIO(urlopen(url,timeout=150).read())) as archive:
                path.write_bytes(archive.read(next(f for f in archive.namelist() if f.endswith(".fits"))))
        with fits.open(path) as hdr:
            wcs=WCS(hdr[2].header)
            data=hdr[1].data
            mask=(data["QUALITY"]==0)&np.isfinite(data["TIME"])
            flux=data["FLUX"][mask].astype(float)
            ok=np.all(np.isfinite(flux),axis=(1,2))
            flux=flux[ok]
            times=data["TIME"][mask][ok]
        edge=np.ones(flux.shape[1:],bool)
        edge[1:-1,1:-1]=False
        flux-=np.median(flux[:,edge],axis=1)[:,None,None]
        cube=flux.reshape(len(times),-1)
        diff=contrast(times,cube,row.btjd)
        trials=np.arange(times.min()+2*DURATION,times.max()-2*DURATION,.8*DURATION)
        blanks=[contrast(times,cube,t) for t in trials if abs(t-row.btjd)>3*DURATION]
        blank=np.array([b for b in blanks if b is not None])
        med=np.median(blank,axis=0)
        sig=1.4826*np.median(np.abs(blank-med),axis=0)
        img=((diff-med)/sig).reshape(flux.shape[1:])
        im=ax.imshow(img,origin="lower",vmin=-3,vmax=3,cmap="RdBu_r",interpolation="nearest")
        x,y=wcs.world_to_pixel_values(gaia.RA_ICRS,gaia.DE_ICRS)
        ax.scatter(x[0],y[0],s=100,marker="+",color="#101010",label="Target")
        ax.scatter(x[1:],y[1:],s=48,facecolor="none",edgecolor="#272727",label="Gaia neighbors")
        ax.set(title=f"Sector {row.sector}",xlabel="Detector column",ylabel="Detector row")
        print(row.sector,"target pixel (column, row)",round(x[0],1),round(y[0],1),
              "highest contrast (row, column)",np.unravel_index(np.argmax(img),img.shape),
              "highest contrast",round(img.max(),1))
    axes[0].legend(loc="upper left",fontsize=8)
    fig.colorbar(im,ax=axes,label="Dimming / local pixel scatter",shrink=.8)
    fig.suptitle("TIC 230091879 · difference images at the two extra epochs")
    fig.savefig(OUT/"lead-pixels.png",dpi=145)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--cache",type=Path,required=True)
    main(parser.parse_args().cache)
