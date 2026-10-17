"""ERA5 reanalysis counterpart to huracanpy_TD.ipynb.

Computes and plots per-domain track density (full track + one panel per
WCSI nature category: TC/BC/Tr/MV/Vo/Ot), the same way huracanpy_TD.ipynb
does for the model experiments - but for the ERA5-tracked cyclones produced
by huracanpy_era5_batch_all.py, which writes its per-year, per-domain
{tag}_TCs.nc files using the exact same naming/location convention
(BASE_DIR/ERA5_{year}_{domain}_TCs.nc) the notebook's load_domain_data()
already expects. Because of that, this only needs to run huracanpy_TD.ipynb's
own plot_domain_density with experiment="ERA5" - see there for anything not
commented here.
"""

import os

import cartopy.crs as ccrs
import huracanpy
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator

from huracanpy_identify_batch_all import BASE_DIR, DOMAINS

EXPERIMENT = "ERA5"
YEARS = range(1980, 2014)  # matches huracanpy_era5_batch_all.py's default
GRID_SPACING = 10  # degrees - lon/lat gridline interval
GRID_LABEL_SPACING = 20  # degrees - interval of the lon/lat labels on the map border


def load_domain_data(experiment, years, domain_name):
    """Load and concatenate lon/lat/nature across years for one domain.

    Returns (None, None, None) if no _TCs.nc files exist yet for this
    domain/experiment (e.g. a domain with zero qualifying TCs, like MED).
    """
    lons, lats, nats = [], [], []
    for year in years:
        in_nc = f"{BASE_DIR}/{experiment}_{year}_{domain_name}_TCs.nc"
        if not os.path.exists(in_nc):
            continue
        tracks = huracanpy.load(in_nc)
        lons.append(tracks.lon)
        lats.append(tracks.lat)
        nats.append(tracks.nature)

    if not lons:
        return None, None, None

    return (
        xr.concat(lons, dim="record"),
        xr.concat(lats, dim="record"),
        xr.concat(nats, dim="record"),
    )


def plot_domain_density(domain_name, domain, experiment, years):
    lon_all, lat_all, nat_all = load_domain_data(experiment, years, domain_name)
    if lon_all is None:
        print(f"[{domain_name}] no _TCs.nc files found, skipping")
        return

    def subset_density(mask=None):
        lon_sub = lon_all if mask is None else lon_all.isel(record=mask)
        if lon_sub.sizes.get("record", 0) == 0:
            return None
        lat_sub = lat_all if mask is None else lat_all.isel(record=mask)
        return huracanpy.calc.density(lon_sub, lat_sub, method="histogram", bin_size=3, spherical=False)

    NATURE_PANELS = [
        ("Full track", subset_density()),
        ("TC", subset_density((nat_all == "TC").values)),
        ("BC", subset_density((nat_all == "BC").values)),
        ("Tr", subset_density((nat_all == "Tr").values)),
        ("MV", subset_density((nat_all == "MV").values)),
        ("Vo", subset_density((nat_all == "Vo").values)),
        ("Ot", subset_density((nat_all == "Ot").values)),
    ]

    vmax = max((float(D_panel.max()) for _, D_panel in NATURE_PANELS if D_panel is not None), default=0)
    if vmax <= 0:
        print(f"[{domain_name}] all-zero density, skipping plot")
        return
    vmin = vmax * 1e-2
    levels = np.logspace(np.log10(vmin), np.log10(vmax), 20)
    norm = LogNorm(vmin=vmin, vmax=vmax, clip=False)

    cmap = plt.get_cmap("terrain_r").copy()
    cmap.set_under("white")

    extent = domain["extent"]
    central_longitude = 0.5 * (extent[0] + extent[1])
    domain_crs = ccrs.PlateCarree(central_longitude=central_longitude)
    map_aspect = (extent[1] - extent[0]) / (extent[3] - extent[2])

    col_width = 4.2
    row_height = col_width / map_aspect
    fig_width = 3 * col_width + 2.3
    fig_height = 5 * row_height + 0.6

    fig = plt.figure(figsize=(fig_width, fig_height), layout="constrained")
    gs = GridSpec(3, 3, figure=fig, height_ratios=[3, 1, 1])
    ax_big = fig.add_subplot(gs[0, :], projection=domain_crs)
    axes_small = [
        fig.add_subplot(gs[1, 0], projection=domain_crs),
        fig.add_subplot(gs[1, 1], projection=domain_crs),
        fig.add_subplot(gs[1, 2], projection=domain_crs),
        fig.add_subplot(gs[2, 0], projection=domain_crs),
        fig.add_subplot(gs[2, 1], projection=domain_crs),
        fig.add_subplot(gs[2, 2], projection=domain_crs),
    ]

    mesh = None
    for ax, (title, D_panel) in zip([ax_big] + axes_small, NATURE_PANELS):
        ax.set_xlim(extent[0] - central_longitude, extent[1] - central_longitude)
        ax.set_ylim(extent[2], extent[3])
        ax.coastlines()
        ax.gridlines(
            xlocs=MultipleLocator(GRID_SPACING), ylocs=MultipleLocator(GRID_SPACING),
            color="0.3", linewidth=0.4, linestyle=":",
        )
        # Second, invisible gridliner just for the border labels - lets them be
        # sparser than the drawn lines so they stay legible on the small panels.
        ax.gridlines(
            xlocs=MultipleLocator(GRID_LABEL_SPACING), ylocs=MultipleLocator(GRID_LABEL_SPACING),
            draw_labels=["left", "bottom"], linewidth=0,
            xlabel_style={"size": 7}, ylabel_style={"size": 7},
        )
        if D_panel is None:
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center", va="center", fontsize=9, color="0.5")
        else:
            mesh = D_panel.plot(
                ax=ax, transform=ccrs.PlateCarree(), levels=levels, cmap=cmap, norm=norm, add_colorbar=False,
            )
        ax.set_title(title)

    fig.colorbar(mesh, ax=[ax_big] + axes_small, shrink=0.6, pad=0.02, label="TC Track Density (log scale)")

    fig.suptitle(f"{domain['title']} TC Track Density for {experiment} ({years.start}-{years.stop - 1})")
    fig.savefig(f"{BASE_DIR}/density_{experiment}_{years.start}-{years.stop - 1}_{domain_name}_MultiNature.png", dpi=300)
    plt.close(fig)


def main():
    for domain_name, domain in DOMAINS.items():
        plot_domain_density(domain_name, domain, EXPERIMENT, YEARS)


if __name__ == "__main__":
    main()
