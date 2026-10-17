"""Ensemble counterpart to huracanpy_TD.ipynb.

Combines several experiments' per-domain track densities into one
multi-panel plot per basin, laid out exactly like huracanpy_TD.ipynb's
per-experiment MultiNature plots (big "Full track" panel + one small panel
per WCSI nature category: TC/BC/Tr/MV/Vo/Ot). Each panel shows the ensemble
mean across experiments in colour (log scale, same style as the
per-experiment plots) and the ensemble spread across experiments (std dev,
smoothed black contour lines).

Densities are computed with an explicit lon_range=(0, 360)/lat_range=(-90, 90)
so every experiment lands on the identical histogram grid - required for the
per-gridpoint mean/std across experiments to be meaningful. huracanpy.load
already returns lon in 0-360 for all four experiments, so this matches what
huracanpy_TD.ipynb effectively used already (its lon_range=None default
resolves to the same (0, 360) grid whenever every experiment's lon.min() >= 0).

A nature category with zero storms for a given experiment still contributes
an all-zero density grid to that category's ensemble (a real "this member saw
none of these" value) - only an experiment missing _TCs.nc files entirely for
a domain is excluded from the ensemble.
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
from scipy.ndimage import gaussian_filter

from huracanpy_identify_batch_all import BASE_DIR, DOMAINS

# Which run's _TCs.nc files to read (a folder under TC-TRACKS), e.g.
# TC_RUN_DIR=TCIDENT_WC_RUN python huracanpy_TD_ensemble.py
BASE_DIR = os.path.join(BASE_DIR, os.environ.get("TC_RUN_DIR", "NEW_WCSI_RUN"))

EXPERIMENTS = ["u-ch330", "u-ck777", "u-ck778", "u-dd828","u-dc009"]
YEARS = range(1980, 2014)
BIN_SIZE = 3
LON_RANGE = (0, 360)
LAT_RANGE = (-90, 90)
GRID_SPACING = 10  # degrees - lon/lat gridline interval
GRID_LABEL_SPACING = 20  # degrees - interval of the lon/lat labels on the map border
SPREAD_SMOOTH_SIGMA = 1.5  # grid cells - see gaussian_filter call in plot_domain_ensemble

SKIP_DOMAINS = {"MED"}
# Small panels, in this order, for the natures that occur in the data (at most
# six): WCSI runs have TC/BC/Tr/MV/Vo/Ot, the tcident run only TC and NT
NATURES = ["TC", "BC", "Tr", "MV", "Vo", "Ot", "NT"]
NATURE_TITLES = {"NT": "Not TC (outside warm-core phase)"}


def load_experiment_data(experiment, years, domain_name):
    """Concatenate lon/lat/nature across years for one experiment/domain.

    Returns (None, None, None) if no _TCs.nc files exist yet for this
    experiment/domain (e.g. an experiment not yet run through the identify
    batch pipeline).
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


def density(lon, lat):
    return huracanpy.calc.density(
        lon, lat, method="histogram", bin_size=BIN_SIZE,
        lon_range=LON_RANGE, lat_range=LAT_RANGE, spherical=False,
    )


def experiment_panels(experiment, years, domain_name):
    """Full-track + per-nature density panels for one experiment.

    Returns None if the experiment has no _TCs.nc files for this domain.
    """
    lon_all, lat_all, nat_all = load_experiment_data(experiment, years, domain_name)
    if lon_all is None:
        return None

    panels = {"Full track": density(lon_all, lat_all)}
    for nature in NATURES:
        mask = (nat_all == nature).values
        if mask.any():
            panels[nature] = density(lon_all.isel(record=mask), lat_all.isel(record=mask))
    return panels


def plot_domain_ensemble(domain_name, domain, experiments, years):
    members = {}
    for experiment in experiments:
        panels = experiment_panels(experiment, years, domain_name)
        if panels is None:
            print(f"[{domain_name}] {experiment}: no _TCs.nc files found, excluding from ensemble")
            continue
        members[experiment] = panels

    if not members:
        print(f"[{domain_name}] no experiments have data, skipping")
        return

    n_members = len(members)
    panel_titles = ["Full track"] + [n for n in NATURES if any(n in panels for panels in members.values())]
    if len(panel_titles) > 7:
        raise ValueError(f"[{domain_name}] more natures than small panels: {panel_titles[1:]}")
    ens_mean, ens_spread = {}, {}
    for title in panel_titles:
        # a member with none of this nature contributes an all-zero grid
        stacked = xr.concat(
            [panels.get(title, xr.zeros_like(panels["Full track"])) for panels in members.values()],
            dim="experiment",
        )
        ens_mean[title] = stacked.mean(dim="experiment")
        ens_spread[title] = stacked.std(dim="experiment", ddof=1) if n_members > 1 else None

    vmax = max(float(D.max()) for D in ens_mean.values())
    if vmax <= 0:
        print(f"[{domain_name}] all-zero ensemble mean density, skipping plot")
        return
    vmin = vmax * 1e-2  # floor, as in huracanpy_TD.ipynb - avoids sparse-bin noise dominating the scale
    levels = np.logspace(np.log10(vmin), np.log10(vmax), 20)
    norm = LogNorm(vmin=vmin, vmax=vmax, clip=False)

    cmap = plt.get_cmap("terrain_r").copy()
    cmap.set_under("white")

    extent = domain["extent"]
    central_longitude = 0.5 * (extent[0] + extent[1])
    domain_crs = ccrs.PlateCarree(central_longitude=central_longitude)
    map_aspect = (extent[1] - extent[0]) / (extent[3] - extent[2])

    # Same layout math as huracanpy_TD.ipynb's plot_domain_density.
    col_width = 4.2
    row_height = col_width / map_aspect
    fig_width = 3 * col_width + 2.3
    fig_height = 5 * row_height + 0.6

    fig = plt.figure(figsize=(fig_width, fig_height), layout="constrained")
    gs = GridSpec(3, 3, figure=fig, height_ratios=[3, 1, 1])
    ax_big = fig.add_subplot(gs[0, :], projection=domain_crs)
    small_slots = [gs[1, 0], gs[1, 1], gs[1, 2], gs[2, 0], gs[2, 1], gs[2, 2]]
    axes_small = [fig.add_subplot(slot, projection=domain_crs) for slot in small_slots[:len(panel_titles) - 1]]

    mesh = None
    for ax, title in zip([ax_big] + axes_small, panel_titles):
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

        mesh = ens_mean[title].plot(
            ax=ax, transform=ccrs.PlateCarree(), levels=levels, cmap=cmap, norm=norm, add_colorbar=False,
        )

        spread = ens_spread[title]
        if spread is not None:
            # Gaussian-smooth the std field before contouring - the raw per-bin
            # std is noisy (each bin's spread depends on just a few tracks), which
            # makes the unsmoothed contour lines jagged. Only the contour input is
            # smoothed; the colour-filled ensemble mean above is left untouched.
            smoothed = gaussian_filter(spread.values, sigma=SPREAD_SMOOTH_SIGMA)
            smax = float(smoothed.max())
            if smax > 0:
                cs = ax.contour(
                    spread.lon, spread.lat, smoothed,
                    levels=np.linspace(smax * 0.2, smax, 4), colors="black", linewidths=0.8,
                    transform=ccrs.PlateCarree(),
                )
                ax.clabel(cs, inline=True, fontsize=6, fmt="%.1f")

        ax.set_title(NATURE_TITLES.get(title, title))  # after .plot() - xarray's plot() clears any pre-existing title

    fig.colorbar(
        mesh, ax=[ax_big] + axes_small, shrink=0.6, pad=0.02,
        label="Ensemble mean TC Track Density (log scale)",
    )

    if n_members <= 1:
        print(f"[{domain_name}] only one experiment has data ({next(iter(members))}); no spread contours")

    fig.suptitle(
        f"{domain['title']} TC Track Density Ensemble ({years.start}-{years.stop - 1})\n"
        f"Mean of {n_members} member(s): {', '.join(members.keys())} - black contours: ensemble std dev"
    )
    fig.savefig(
        f"{BASE_DIR}/density_ENSEMBLE_{years.start}-{years.stop - 1}_{domain_name}_MultiNature.png", dpi=300,
    )
    plt.show()


def main():
    for domain_name, domain in DOMAINS.items():
        if domain_name in SKIP_DOMAINS:
            continue
        plot_domain_ensemble(domain_name, domain, EXPERIMENTS, YEARS)


if __name__ == "__main__":
    main()
