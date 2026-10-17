"""Where do TCs form? Genesis maps and latitude distributions, per basin.

Companion to huracanpy_TC_interannual.py, written to look into why the model
ensemble's mean genesis latitude is 2-3 degrees closer to the equator than
IBTrACS in most basins. Uses exactly the same storms and genesis definition
as that script (it imports its settings and file-finding helpers):

  - storms: those reaching WIND_THRESHOLD_KT (34 kt) while tropical; every
    is_tc track for ERA5 unless REANALYSIS_USE_THRESHOLD is set
  - genesis: the first tropical point - first is_tc point for the models and
    ERA5 (WCSI structural criteria, no wind test), first nature == "TS" point
    for IBTrACS (usually the agency's first tropical-depression fix)

One figure per domain:

  top row     genesis density maps for IBTrACS, ERA5 and the model ensemble,
              on one colour scale, in genesis events per year per
              GRID_DEG x GRID_DEG box (models: averaged over member-years)
  bottom left distribution of genesis latitude (solid) and of the latitude
              where each storm first reaches 34 kt while tropical (dashed),
              as storms per year per degree of latitude
  bottom right
              model storms: genesis latitude vs latitude of first reaching
              34 kt (2-D histogram), IBTrACS storms as dots, 1:1 line. Points
              on the line deepened where they formed; points well above it
              drifted poleward between genesis and tropical-storm strength.

The reanalyses (REANALYSES: ERA5 and JRA3Q) are the key comparison: they go
through the same WCSI identification as the models, so a near-equatorial
genesis tail in the models but not in the reanalyses points at the model
itself rather than at the identification criteria. Their per-year files come
from huracanpy_era5_batch_all.py (run with --name JRA3Q for JRA3Q). Note the
JRA3Q input (JRA3Q_nolat-nwc-tcident_WCSI_nature.nc) holds only tracks that
already passed TRACK's own TC identification, whereas ERA5_all_nature.nc
holds every tracked vortex, so the two start from different candidate sets.
Reanalyses are left out of the dashed lines and the bottom-right panel
because their winds rarely reach 34 kt (see huracanpy_TC_interannual.py).

Latitudes in the bottom panels are absolute values (degrees poleward).
Also writes every storm's genesis and first-34-kt positions to a CSV.
"""

import os

import cartopy.crs as ccrs
import huracanpy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator

from huracanpy_TC_interannual import (
    BASE_DIR, COLOR_ENSEMBLE, COLOR_IBTRACS, DEFAULT_WIND_VARIABLE, DOMAINS, EXPERIMENTS,
    IBTRACS_DOMAINS, IBTRACS_FILE, IBTRACS_TC_NATURE, KT_TO_MS, REANALYSES, REANALYSIS_USE_THRESHOLD,
    SKIP_DOMAINS, WIND_THRESHOLD_KT, WIND_VARIABLES, YEARS, find_output, load_ibtracs,
)

GRID_DEG = 2.5  # genesis density box size (degrees)
LAT_BIN_DEG = 1  # latitude-distribution bin width (degrees)
STRENGTH_KT = 34  # "first reaches tropical-storm strength" latitude
GRID_SPACING = 10  # map gridlines, as in the track-density scripts
GRID_LABEL_SPACING = 20


# ---------------------------------------------------------------------------
# Per-storm genesis tables
# ---------------------------------------------------------------------------
def genesis_table(points, threshold_kt):
    """One row per counted track: genesis lat/lon (first tropical point) and
    lat/lon of the first tropical point with wind >= STRENGTH_KT (NaN if
    never). `points` has columns track_id, time, lat, lon, wind_kt, tropical.
    """
    tropical = points[points.tropical].sort_values(["track_id", "time"])
    if threshold_kt is None:
        counted = tropical.groupby("track_id").size() > 0
    else:
        counted = (tropical.wind_kt >= threshold_kt).groupby(tropical.track_id).any()
    first = tropical.groupby("track_id")[["lat", "lon"]].first()
    strong = tropical[tropical.wind_kt >= STRENGTH_KT].groupby("track_id")[["lat", "lon"]].first()
    table = first.add_prefix("genesis_").join(strong.add_prefix("ts_"))
    return table[counted.reindex(table.index, fill_value=False)].reset_index()


def model_genesis(experiment, domain_name, years, threshold_kt, wind_variable=DEFAULT_WIND_VARIABLE):
    """(storms, n_years): genesis table over all processed years of one
    experiment (or reanalysis), and the number of processed years (for
    per-year densities).
    """
    tables, n_years = [], 0
    for year in years:
        tag = f"{experiment}_{year}_{domain_name}"
        summary = find_output(f"{tag}_summary.parquet")
        if summary is None:
            continue
        n_years += 1
        tcs_nc = os.path.join(os.path.dirname(summary), f"{tag}_TCs.nc")  # see huracanpy_TC_interannual.model_storms
        if not os.path.exists(tcs_nc):
            continue
        tracks = huracanpy.load(tcs_nc)
        wind = tracks[wind_variable].where(tracks[wind_variable] < 1e10).values
        points = pd.DataFrame(dict(
            track_id=tracks.track_id.values, time=pd.to_datetime(tracks.time.values),
            lat=tracks.lat.values, lon=tracks.lon.values % 360,
            wind_kt=wind / KT_TO_MS, tropical=tracks.is_tc.values,
        ))
        tables.append(genesis_table(points, threshold_kt).assign(year=year))
    storms = pd.concat(tables, ignore_index=True) if tables else None
    return storms, n_years


def ibtracs_genesis(all_tracks, domain_name, years, threshold_kt):
    """Genesis table for one IBTrACS domain, years on the model's season
    convention (see huracanpy_TC_interannual.py). None if no such domain.
    """
    if domain_name not in IBTRACS_DOMAINS:
        return None
    sel = all_tracks[all_tracks.domain_basin == IBTRACS_DOMAINS[domain_name]["basin"]]
    points = pd.DataFrame(dict(
        track_id=sel.track_id.values, time=pd.DatetimeIndex(sel.time),
        lat=sel.lat.values, lon=sel.lon.values % 360,
        wind_kt=sel.wind.values, tropical=(sel.nature == IBTRACS_TC_NATURE).values,
    ))
    storms = genesis_table(points, threshold_kt)
    season = sel.groupby("track_id").season.first()
    year = season - 1 if DOMAINS[domain_name]["hemisphere"] == "S" else season
    storms["year"] = storms.track_id.map(year)
    return storms[storms.year.isin(years)].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def density_map(storms, n_years):
    """Genesis events per year per GRID_DEG box on a global 0-360 grid."""
    lon_edges = np.arange(0, 360 + GRID_DEG, GRID_DEG)
    lat_edges = np.arange(-90, 90 + GRID_DEG, GRID_DEG)
    counts, _, _ = np.histogram2d(storms.genesis_lon, storms.genesis_lat, bins=[lon_edges, lat_edges])
    return lon_edges, lat_edges, counts.T / n_years


def latitude_pdf(lat, n_years, edges):
    """Storms per year per degree of (absolute) latitude."""
    counts, _ = np.histogram(np.abs(lat.dropna()), bins=edges)
    return counts / n_years / LAT_BIN_DEG


def setup_map(ax, extent, central_longitude):
    ax.set_xlim(extent[0] - central_longitude, extent[1] - central_longitude)
    ax.set_ylim(extent[2], extent[3])
    ax.coastlines(linewidth=0.6)
    ax.gridlines(
        xlocs=MultipleLocator(GRID_SPACING), ylocs=MultipleLocator(GRID_SPACING),
        color="0.3", linewidth=0.4, linestyle=":",
    )
    ax.gridlines(
        xlocs=MultipleLocator(GRID_LABEL_SPACING), ylocs=MultipleLocator(GRID_LABEL_SPACING),
        draw_labels=["left", "bottom"], linewidth=0,
        xlabel_style={"size": 7}, ylabel_style={"size": 7},
    )
    ax.axhline(0, color="0.2", linewidth=0.8)  # equator (y is latitude in any PlateCarree)


def plot_domain(domain_name, domain, datasets, years, out_tag):
    """datasets: list of (label, colour, storms, n_years, is_reanalysis);
    storms may be None.
    """
    datasets = [d for d in datasets if d[2] is not None and len(d[2])]
    n_maps = len(datasets)
    extent = domain["extent"]
    central_longitude = 0.5 * (extent[0] + extent[1])
    domain_crs = ccrs.PlateCarree(central_longitude=central_longitude)
    map_aspect = (extent[1] - extent[0]) / (extent[3] - extent[2])

    map_width = 4.6
    fig = plt.figure(figsize=(n_maps * map_width + 1.2, map_width / map_aspect + 5.6), layout="constrained")
    # n_maps columns of maps on top; two equal panels below, each spanning half
    gs = GridSpec(2, 2 * n_maps, figure=fig, height_ratios=[map_width / map_aspect, 4.4])

    # --- genesis density maps, one shared colour scale ---------------------
    grids = [density_map(storms, n_years) for _, _, storms, n_years, _ in datasets]
    vmax = max(grid[2].max() for grid in grids)
    norm = LogNorm(vmin=vmax / 100, vmax=vmax)
    mesh = None
    map_axes = []
    for i, ((label, _, storms, n_years, _), (lon_edges, lat_edges, density)) in enumerate(zip(datasets, grids)):
        ax = fig.add_subplot(gs[0, 2 * i:2 * i + 2], projection=domain_crs)
        setup_map(ax, extent, central_longitude)
        mesh = ax.pcolormesh(
            lon_edges, lat_edges, np.ma.masked_equal(density, 0), cmap="Blues", norm=norm,
            transform=ccrs.PlateCarree(),
        )
        ax.set_title(f"{label}\n{len(storms) / n_years:.1f} storms/yr", fontsize=9)
        map_axes.append(ax)
    fig.colorbar(
        mesh, ax=map_axes, shrink=0.8, pad=0.01,
        label=f"Genesis per year per {GRID_DEG}°×{GRID_DEG}° box",
    )

    # --- latitude distributions ----------------------------------------------
    ax = fig.add_subplot(gs[1, :n_maps])
    edges = np.arange(0, 51, LAT_BIN_DEG)
    centres = 0.5 * (edges[:-1] + edges[1:])
    for label, colour, storms, n_years, is_reanalysis in datasets:
        ax.plot(
            centres, latitude_pdf(storms.genesis_lat, n_years, edges), color=colour, linewidth=2,
            label=f"{label}: genesis (mean {np.abs(storms.genesis_lat).mean():.1f}°)",
        )
        if not is_reanalysis:  # see module docstring
            ax.plot(
                centres, latitude_pdf(storms.ts_lat, n_years, edges), color=colour, linewidth=1.5,
                linestyle="--",
                label=f"{label}: first ≥ {STRENGTH_KT} kt (mean {np.abs(storms.ts_lat).mean():.1f}°)",
            )
    ax.set_xlim(0, 45)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(f"Latitude (°{domain['hemisphere']})")
    ax.set_ylabel("Storms per year per degree latitude")
    ax.grid(color="0.85", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=8, frameon=False)
    ax.set_title("Latitude of genesis and of first reaching tropical-storm strength", fontsize=9)

    # --- genesis vs first-34-kt latitude ---------------------------------------
    ax = fig.add_subplot(gs[1, n_maps:])
    lat_max = 45
    model = next((d for d in datasets if d[0].startswith("Ensemble")), None)
    if model is not None:
        storms = model[2].dropna(subset=["ts_lat"])
        h = ax.hist2d(
            np.abs(storms.genesis_lat), np.abs(storms.ts_lat), bins=np.arange(0, lat_max + 1, LAT_BIN_DEG),
            cmap="Blues", norm=LogNorm(), weights=np.full(len(storms), 1 / model[3]),
        )
        fig.colorbar(h[3], ax=ax, pad=0.01, label="Model storms per year per bin")
    obs = next((d for d in datasets if d[0] == "IBTrACS"), None)
    if obs is not None:
        storms = obs[2].dropna(subset=["ts_lat"])
        ax.scatter(
            np.abs(storms.genesis_lat), np.abs(storms.ts_lat), s=4, color=COLOR_IBTRACS, alpha=0.5,
            linewidths=0, label="IBTrACS storms",
        )
        ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.plot([0, lat_max], [0, lat_max], color="0.4", linewidth=0.8, linestyle=":")
    ax.set_xlim(0, lat_max)
    ax.set_ylim(0, lat_max)
    ax.set_xlabel(f"Genesis latitude (°{domain['hemisphere']})")
    ax.set_ylabel(f"Latitude of first ≥ {STRENGTH_KT} kt (°{domain['hemisphere']})")
    ax.set_title("Genesis vs first reaching tropical-storm strength", fontsize=9)

    threshold_text = "all TCs" if WIND_THRESHOLD_KT is None else f"storms reaching ≥ {WIND_THRESHOLD_KT} kt"
    fig.suptitle(
        f"{domain['title']}: TC genesis, {threshold_text}, {years.start}-{years.stop - 1}\n"
        f"Genesis = first is_tc point (models, reanalyses) / first TS point (IBTrACS)",
        fontsize=10,
    )
    suffix = "" if WIND_THRESHOLD_KT is None else f"_min{WIND_THRESHOLD_KT}kt"
    fig.savefig(
        f"{BASE_DIR}/TC_genesis_maps_{out_tag}{suffix}_{years.start}-{years.stop - 1}_{domain_name}.png", dpi=300,
    )
    plt.close(fig)


def main():
    print(f"Loading {IBTRACS_FILE} ...")
    all_tracks = load_ibtracs(IBTRACS_FILE)
    threshold = WIND_THRESHOLD_KT
    reanalysis_threshold = threshold if REANALYSIS_USE_THRESHOLD else None
    # only reanalyses with files in this run folder; the outputs are named after
    # them, so runs with other sets don't overwrite each other
    reanalysis_names = [
        name for name in REANALYSES
        if any(f.startswith(f"{name}_") and f.endswith("_summary.parquet") for f in os.listdir(BASE_DIR))
    ]
    out_tag = "-".join(reanalysis_names) or "no-reanalysis"

    all_storms = []
    for domain_name, domain in DOMAINS.items():
        if domain_name in SKIP_DOMAINS:
            continue
        obs = ibtracs_genesis(all_tracks, domain_name, YEARS, threshold)
        reanalyses = {
            name: model_genesis(name, domain_name, YEARS, reanalysis_threshold, WIND_VARIABLES[name])
            for name in reanalysis_names
        }

        members, member_years = [], 0
        for experiment in EXPERIMENTS:
            storms, n_years = model_genesis(experiment, domain_name, YEARS, threshold)
            if storms is not None:
                members.append(storms.assign(dataset=experiment))
                member_years += n_years
        ensemble = pd.concat(members, ignore_index=True) if members else None

        label_suffix = "" if reanalysis_threshold == threshold else " (all is_tc)"
        plot_domain(domain_name, domain, [
            ("IBTrACS", COLOR_IBTRACS, obs, len(YEARS), False),
            *[
                (name + label_suffix, REANALYSES[name], storms, n_years, True)
                for name, (storms, n_years) in reanalyses.items()
            ],
            (f"Ensemble ({len(members)} members)", COLOR_ENSEMBLE, ensemble, member_years, False),
        ], YEARS, out_tag)
        print(f"[{domain_name}] done")

        for storms, dataset in [(obs, "IBTrACS"), *[(st, name) for name, (st, _) in reanalyses.items()]]:
            if storms is not None:
                all_storms.append(storms.assign(dataset=dataset, domain=domain_name))
        if ensemble is not None:
            all_storms.append(ensemble.assign(domain=domain_name))

    if all_storms:
        suffix = "" if threshold is None else f"_min{threshold}kt"
        columns = ["domain", "dataset", "year", "track_id", "genesis_lat", "genesis_lon", "ts_lat", "ts_lon"]
        pd.concat(all_storms, ignore_index=True)[columns].to_csv(
            f"{BASE_DIR}/TC_genesis_points_{out_tag}{suffix}_{YEARS.start}-{YEARS.stop - 1}.csv", index=False,
        )


if __name__ == "__main__":
    main()
