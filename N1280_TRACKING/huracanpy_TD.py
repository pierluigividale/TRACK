import cartopy.crs as ccrs
import matplotlib.pyplot as plt

import huracanpy

BASE_DIR = "/Users/vidale/Python_scripts/TC-TRACKS"
#EXPERIMENTS = ["u-dd828"]
EXPERIMENTS = ["u-ch330"]
YEARS = range(1980, 1990)
domain_name = "NATL"

nat_tracks = []
for experiment in EXPERIMENTS:
    for year in YEARS:
        tag = f"{BASE_DIR}/{experiment}_{year}_{domain_name}"
        in_nc = f"{tag}_TCs.nc"
        print(in_nc)
        tracks = huracanpy.load(in_nc)
        nat_tracks.append(tracks)
nat_tracks
#nat_tracks.hrcn.plot_density()
plt.show()
