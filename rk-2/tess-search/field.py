"""Resolve the host, nearby Gaia sources and positional catalog matches after detection."""
import io
import json
from urllib.parse import urlencode

from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd

from search import OUT, EXOFOP, read

TIC = 230091879
VIZ = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"
MAST = "https://mast.stsci.edu/api/v0/invoke?"
TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?"
EB = "https://tessebs.villanova.edu/230091879"


def main():
    req = dict(service="Mast.Catalogs.Filtered.Tic",params=dict(columns="*",
               filters=[dict(paramName="ID",values=[str(TIC)])]),format="json",pagesize=5,page=1)
    star = json.loads(read(MAST+urlencode(dict(request=json.dumps(req)))))["data"][0]
    fields = ("ID TWOMASS GAIA ra dec RA_orig Dec_orig Tmag Teff e_Teff rad e_rad mass e_mass "
              "lum e_lum d e_d pmRA pmDEC").split()
    (OUT/"lead-star.json").write_text(json.dumps({k:star[k] for k in fields},indent=2)+"\n")
    center = SkyCoord(star["ra"],star["dec"],unit="deg")
    print("Sky region",center.get_constellation())
    for catalog in ("toi","ctoi"):
        frame=pd.read_csv(io.BytesIO(read(f"{EXOFOP}download_{catalog}.php?sort={catalog}&output=csv")),
                          usecols=["TIC ID","RA","Dec"]).dropna()
        units=(u.hourangle,u.deg) if catalog=="toi" else (u.deg,u.deg)
        coords=SkyCoord(frame.RA.to_numpy(),frame.Dec.to_numpy(),unit=units)
        print(catalog,"hosts within 2 arcmin:",frame.loc[center.separation(coords)<2*u.arcmin,"TIC ID"].tolist())
    query=("select pl_name from pscomppars where 1=contains(point('ICRS',ra,dec),"
           f"circle('ICRS',{star['ra']},{star['dec']},{2/60}))")
    hits=pd.read_csv(io.BytesIO(read(TAP+urlencode(dict(query=query,format="csv")))))
    print("Confirmed planets within 2 arcmin:",hits.pl_name.tolist())
    print("Live TESS EB catalog reports absent:", b"This TESS ID is not in the catalog" in read(EB))
    opts={"-source":"I/355/gaiadr3","-c":f"{star['RA_orig']} +{star['Dec_orig']}",
          "-c.rs":65,"-out":"Source,RA_ICRS,DE_ICRS,Gmag,RUWE,NSS,VarFlag","-sort":"Gmag","-out.max":200}
    raw=read(VIZ+urlencode(opts)).decode()
    lines=[s for s in raw.splitlines() if s.split("\t",1)[0].strip().isdigit()]
    if not lines:
        raise ValueError("No Gaia sources returned")
    sources=pd.read_csv(io.StringIO("\n".join(lines)),sep="\t",names=opts["-out"].split(","),
                        dtype={"Source":"string"},skipinitialspace=True)
    sky=SkyCoord(sources.RA_ICRS,sources.DE_ICRS,unit="deg")
    sources["offset_arcsec"]=sky[0].separation(sky).arcsec
    sources.to_csv(OUT/"lead-gaia-field.csv",index=False)
    var=read(VIZ+urlencode({"-source":"I/358/vclassre","Source":sources.Source.iloc[0]})).decode()
    print("Gaia variability class:",next(line.split("\t")[3].strip() for line in var.splitlines()
                                          if line.startswith(sources.Source.iloc[0])))


if __name__ == "__main__":
    main()
