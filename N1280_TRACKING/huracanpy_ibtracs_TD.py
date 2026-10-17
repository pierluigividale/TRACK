"""Observational counterpart to huracanpy_TD.ipynb, using IBTrACS.

Computes and plots per-domain track density (full track + one panel per nature
category), the same way huracanpy_TD.ipynb does for the model experiments -
except the per-nature panels are IBTrACS's own TS/SS/ET/DS/MX categories
(see huracanpy_ibtracs_batch_all.py) instead of the WCSI TC/BC/Tr/MV/Vo/Ot
categories, and everything is computed directly from the IBTrACS NetCDF rather
than from a batch script's saved per-year _TCs.nc files - IBTrACS is one file
covering every year, so there's nothing to loop-and-concatenate the way the
model version has to.
"""

import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator

import huracanpy

from huracanpy_ibtracs_batch_all import DOMAINS, IBTRACS_FILE, NATURE_MARKERS, YEARS, load_ibtracs

BIN_SIZE = 3  # degrees - matches huracanpy_TD.ipynb, for visual comparability
GRID_SPACING = 10  # degrees - lon/lat gridline interval
GRID_LABEL_SPACING = 20  # degrees - interval of the lon/lat labels on the map border
SKIP_DOMAINS = set()  # IBTrACS has ample storms in every domain, unlike the
                       # sparse model runs, so nothing is skipped by default
# Small-panel order (row-major over the 2x3 small slots), placed to sit in the
# same positions as the model plots' TC/BC/Tr/MV/Vo/Ot panels. Only loosely
# comparable: the IBTrACS nature codes and the phase-space (WCSI) categories
# are different classifications. TS vs TC and ET vs BC correspond reasonably
# well, DS is probably closest to the weak-vortex category (Vo); SS and MX
# have no real counterpart and just fill the Tr and Ot positions.
# None leaves a slot empty (the MV position).
SMALL_PANEL_NATURES = ["TS", "ET", "SS", None, "DS", "MX"]
assert set(SMALL_PANEL_NATURES) - {None} == set(NATURE_MARKERS)


def plot_domain_density(domain_name, domain, all_tracks, years):
    sel = all_tracks[
        (all_tracks.season >= years.start) & (all_tracks.season < years.stop)
        & (all_tracks.domain_basin == domain["basin"])
    ]
    if sel.empty:
        print(f"[{domain_name}] no tracks in basin {domain['basin']!r}, skipping")
        return

    def subset_density(mask=None):
        sub = sel if mask is None else sel[mask]
        if sub.empty:
            return None
        return huracanpy.calc.density(
            sub.lon, sub.lat, method="histogram", bin_size=BIN_SIZE, spherical=False
        )

    NATURE_PANELS = [("Full track", subset_density())] + [
        (code, subset_density((sel.nature == code).values))
        for code in SMALL_PANEL_NATURES if code is not None
    ]

    # Log-scaled shared colour range - see huracanpy_TD.ipynb for why (a linear
    # range washes out every nature category next to "Full track" and needs
    # re-tuning per basin; log scale adapts automatically).
    vmax = max((float(D.max()) for _, D in NATURE_PANELS if D is not None), default=0)
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

    # Same 3x3 grid as huracanpy_TD.ipynb (big "Full track" panel spanning row
    # 0, small panels below) but with only 5 nature categories instead of 6,
    # so one of the 6 small slots is left empty (see SMALL_PANEL_NATURES).
    col_width = 4.2
    row_height = col_width / map_aspect
    fig_width = 3 * col_width + 2.3
    fig_height = 5 * row_height + 0.6

    fig = plt.figure(figsize=(fig_width, fig_height), layout="constrained")
    gs = GridSpec(3, 3, figure=fig, height_ratios=[3, 1, 1])
    ax_big = fig.add_subplot(gs[0, :], projection=domain_crs)
    small_slots = [gs[1, 0], gs[1, 1], gs[1, 2], gs[2, 0], gs[2, 1], gs[2, 2]]
    axes_small = [
        fig.add_subplot(slot, projection=domain_crs)
        for slot, code in zip(small_slots, SMALL_PANEL_NATURES) if code is not None
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

    fig.colorbar(mesh, ax=[ax_big] + axes_small, shrink=0.6, pad=0.02, label="Storm Track Density (log scale)")

    fig.suptitle(f"{domain['title']} storm track density, IBTrACS ({years.start}-{years.stop - 1})")
    fig.savefig(f"density_IBTrACS_{years.start}-{years.stop - 1}_{domain_name}_MultiNature.png", dpi=300)
    plt.close(fig)


def main():
    print(f"Loading {IBTRACS_FILE} ...")
    all_tracks = load_ibtracs(IBTRACS_FILE)

    for domain_name, domain in DOMAINS.items():
        if domain_name in SKIP_DOMAINS:
            continue
        plot_domain_density(domain_name, domain, all_tracks, YEARS)


if __name__ == "__main__":
    main()
