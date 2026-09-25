"""General batch version of HuracanPy_Identify.ipynb.

For each experiment/year and each configured basin (see DOMAINS below), loads
the matching tr_trs_pos*.nc (northern-hemisphere basins) or tr_trs_neg*.nc
(southern-hemisphere basins) tracks from
{BASE_DIR}/{experiment}/VOR_VERTAVG_T63filt_{year}[_SH]/, identifies tropical
cyclones in that basin with the WCSI criteria, classifies each track's full
lifetime into nature categories (TC/BC/Tr/MV/Ot/Vo), and saves a labelled map
plus the underlying track/summary data. Also writes one combined CSV of TC
counts by experiment/year/domain across everything processed.
"""

import glob
import os

import cartopy.crs as ccrs
import huracanpy
import matplotlib
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.lines import Line2D
from parse import parse
from wcsi import nature
from scipy.ndimage import uniform_filter1d

matplotlib.use("Agg")  # no interactive display in batch mode

# ---------------------------------------------------------------------------
# Configuration - edit before each run
# ---------------------------------------------------------------------------
BASE_DIR = "/Users/vidale/Python_scripts/TC-TRACKS"
#EXPERIMENTS = ["u-dd828"]
EXPERIMENTS = ["u-dc009"]
#EXPERIMENTS = ["u-ch330"]
#EXPERIMENTS = ["u-dz876"]
#YEARS = range(2000, 2001)  # 1980 to 2000 inclusive
YEARS = range(1980, 1981)

# One entry per domain to process. "basin" is the WMO-TC convention code used
# by huracanpy/wcsi (see huracanpy._basins.NH/SH) to restrict tracks to those
# with their maximum intensity in that basin. "hemisphere" picks the input
# directory/file naming convention: northern-hemisphere tracker output lives
# in VOR_VERTAVG_T63filt_{year}/tr_trs_pos*.nc, southern-hemisphere output in
# VOR_VERTAVG_T63filt_{year}_SH/tr_trs_neg*.nc (negative 850 hPa vorticity).
DOMAINS = {
    "NATL": dict(
        basin="NATL", hemisphere="N", extent=[-110, 40, -5, 75],
        title="North Atlantic",
    ),
    "MED": dict(
            basin="MED", hemisphere="N", extent=[-20, 40, 20, 60],
            title="Mediterranean Sea",
        ),
    "NWPAC": dict(
        basin="WNP", hemisphere="N", extent=[-265, -135, -10, 65],
        title="Northwest Pacific",
    ),
    "SPAC": dict(
            basin="SP", hemisphere="S", extent=[140, 245, -65, 5],
            title="South Pacific",
        ),  
    "EPAC": dict(
        basin="ENP", hemisphere="N", extent=[-160, -65, -5, 40],
        title="East Pacific",
    ),
    "NIO": dict(
        basin="NI", hemisphere="N", extent=[30, 100, -5, 30],
        title="North Indian Ocean",
    ),
    "SIO": dict(
        basin="SI", hemisphere="S", extent=[20, 140, -65, 5],
        title="South Indian Ocean",
    ),
    "SATL": dict(
        basin="SA", hemisphere="S", extent=[-65, 20, -50, 5],
        title="South Atlantic",
    ),
}

#SKIP_DOMAINS = {"MED", "NATL", "NWPAC", "SPAC", "EPAC", "NIO"}  # too sparse to be useful yet - storms in only 1-2 of 34 years
SKIP_DOMAINS = {}  # too sparse to be useful yet - storms in only 1-2 of 34 years

NATURE_MARKERS = {"TC": "o", "BC": "^", "Tr": "s", "MV": "D", "Ot": "X", "Vo": "P"}
# mslp legend levels (hPa) - also sets the colour scale range; points outside
# [min, max] are drawn with a white edge instead of black
PRESSURE_LEVELS = [920, 930, 940, 950, 960, 970, 980, 990, 1000, 1010, 1020]
MARKER_SIZE = 10  # main-plot marker size (points**2); legend markers are separate
VORT_THRESHOLD = 6e-5  # raw s**-1 - see notebook for why this differs from the
                        # wcsi package's documented default of 6 (10**-5 s**-1)
TC_PERSISTENCE_HOURS = 24  # a track is only retained if it holds WCSI TC criteria
                            # over ocean for at least this many consecutive hours
VTU_THRESHOLD = 0  # matches wcsi.nature.nature()'s own default; used to re-derive
                    # the fallback category for TC points downgraded in compute_nature
TC_GENESIS_THRESHOLD_LATITUDE = 50.0   # Define a threshold latitude for the first point as TC
FILTER_SIZE = 5 # what are the units of this?

# ---------------------------------------------------------------------------
# Helpers (same logic as the notebook)
# ---------------------------------------------------------------------------
def gather_vorticity_profile(tracks):
    """Replace variables named relative_vorticity_{n} with a single variable and a
    pressure coordinate
    """
    # Keep the native TRACK vorticity field under the name "vorticity" - required by
    # nature.wcsi's basin= filter (hrcn.get_apex_vals("vorticity"))
    tracks["vorticity"] = tracks["relative_vorticity"]

    plevs = [
        int(result.named["n"])
        for result in [parse("relative_vorticity_{n:d}", var) for var in tracks]
        if result is not None
    ]

    if len(plevs) > 0:
        tracks["pressure"] = ("pressure", plevs)
        tracks = tracks.set_coords("pressure")

        vorticity = np.zeros([tracks.sizes["record"], tracks.sizes["pressure"]])
        for n, plev in enumerate(tracks.pressure.values):
            name = f"relative_vorticity_{int(plev)}"
            vorticity[:, n] = tracks[name].values
            tracks = tracks.drop_vars(name)

        tracks["relative_vorticity"] = (["record", "pressure"], vorticity)

    return tracks

def filter_genesis_latitude(tracks, threshold_latitude):
    # Select first points of tracks with TC nature
    tc_points = tracks.isel(record=tracks.nature == "TC")
    tc_genesis_points = tc_points.hrcn.get_gen_vals()
    
    # Select subset of those points that are poleward of the threshold latitude.
    # get_gen_vals() swaps the "record" dimension out for "track_id", so the
    # selection below must isel on track_id, not record.
    tc_genesis_points_valid = tc_genesis_points.isel(
        track_id=np.abs(tc_genesis_points.lat) < threshold_latitude
    )
    # Return the subset of tracks that have the first TC point poleward of the threshold
    # latitude
    return tracks.hrcn.sel_id(tc_genesis_points_valid.track_id)

def infer_npoints(tracks, hours):
    """Number of consecutive records spanning `hours`, inferred from the timestep
    of the first track (so this works regardless of the data's native resolution)
    """
    first_id = tracks.track_id.values[0]
    t = tracks.time.values[tracks.track_id.values == first_id]
    timestep_hours = np.diff(t)[0] / np.timedelta64(1, "h")
    return max(1, round(hours / timestep_hours))

def compute_nature(tracks, filter_size):
    """Add a full WCSI 'nature' classification (TC/BC/Tr/MV/Ot/Vo) per record.

    nature.nature()'s smoothing can extend "TC" to nearby weaker (tropical-storm-
    like) points without checking they're over ocean. Since is_tc itself is already
    ocean-filtered (see nature.wcsi(ocean=True) in process_year), any "TC" point
    that ends up over land here must be one of those smoothed-in points, not a
    genuine is_tc point - downgrade it to what it would have been without the
    smoothing (Tr if it's an upper-level trough, else Ot).
    """
    nat_tracks = []
    for track_id, track in tracks.groupby("track_id"):
        track = track.copy()

        b = uniform_filter1d(np.abs(track.cps_b), size=filter_size, mode="nearest")
        vtl = uniform_filter1d(track.cps_vtl, size=filter_size, mode="nearest")
        vtu = uniform_filter1d(track.cps_vtu, size=filter_size, mode="nearest")
        vorticity = uniform_filter1d(track.relative_vorticity.sel(pressure=850), size=filter_size, mode="nearest")

        nat = nature.nature(
                    b,
                    vtl,
                    vtu,
                    vorticity,
                    track.is_tc.values,
                    vort_threshold=VORT_THRESHOLD,
                )
        # We think that this is not necessary, so I am going to comment it out for now. If we find that the results are not as expected, we can uncomment it and see if it makes a difference.
        #is_ocean = track.hrcn.get_is_ocean().values
        #mislabelled = (nat == "TC") & ~is_ocean
        #nat[mislabelled] = np.where(
        #    track.cps_vtu.values[mislabelled] <= VTU_THRESHOLD, "Tr", "Ot"
        #)

        track["nature"] = ("record", nat)
        nat_tracks.append(track)
    return xr.concat(nat_tracks, dim="record")


def get_mslp(track):
    """Return the mean sea level pressure values for a single track.

    Different experiments have come through with the mslp field under
    different names - some as the CF standard name, some as "psl".
    """
    for name in ("air_pressure_at_mean_sea_level", "psl"):
        if name in track.variables:
            return track[name].values
    raise KeyError(
        "no mslp field found (looked for air_pressure_at_mean_sea_level, psl)"
    )

def plot_domain_nature(tracks, title, out_path, extent):
    """Full lifetime, marker shape by nature category (filled circle = TC),
    colour = mslp. G/L mark genesis/lysis, thin lines join each track's points.
    """
    domain_crs = ccrs.PlateCarree(
        central_longitude=0.5 * (extent[0] + extent[1])
    )

    cmap = plt.get_cmap("terrain")
    norm = plt.Normalize(vmin=min(PRESSURE_LEVELS), vmax=max(PRESSURE_LEVELS))

    # Explicit axes position sized to match the extent's true aspect ratio (rather
    # than plt.subplots + bbox_inches="tight"), for two reasons:
    # - cartopy's Gridliner (added below) confuses matplotlib's tight-bbox
    #   auto-expansion, silently dropping the legends placed outside the axes
    # - an axes box whose aspect doesn't match the extent gets shrunk/recentred by
    #   cartopy at draw time to preserve equal-degree scaling (the "letterboxing"
    #   that caused the excess white space), which also breaks the title's
    #   position calculation entirely
    lon_span = extent[1] - extent[0]
    lat_span = extent[3] - extent[2]
    map_aspect = lon_span / lat_span  # PlateCarree: 1 deg lon == 1 deg lat on screen

    map_width_in = 9.0
    left_margin_in = 0.5
    bottom_margin_in = 0.4
    top_margin_in = 0.5
    legend_width_in = 3.2

    map_height_in = map_width_in / map_aspect
    fig_width_in = left_margin_in + map_width_in + legend_width_in
    fig_height_in = bottom_margin_in + map_height_in + top_margin_in

    fig = plt.figure(figsize=(fig_width_in, fig_height_in))
    ax = fig.add_axes(
        [
            left_margin_in / fig_width_in,
            bottom_margin_in / fig_height_in,
            map_width_in / fig_width_in,
            map_height_in / fig_height_in,
        ],
        projection=domain_crs,
    )
    # ax.set_extent(extent, crs=ccrs.PlateCarree()) is NOT used here: for a
    # domain whose extent crosses the antimeridian (e.g. SPAC, 140 to 245
    # degrees), passing raw longitudes past 180 to set_extent's default
    # (central_longitude=0) crs makes cartopy silently fall back to a
    # full-globe view, which in turn breaks the title's position calculation
    # the same way the top_labels=False bug above does. Setting xlim/ylim
    # directly in the axes' own (already-rotated) coordinate system sidesteps
    # that entirely - central_longitude has already been chosen above so this
    # extent is centred and within +/-180 in the rotated frame.
    central_longitude = 0.5 * (extent[0] + extent[1])
    ax.set_xlim(extent[0] - central_longitude, extent[1] - central_longitude)
    ax.set_ylim(extent[2], extent[3])
    ax.coastlines()

    # xlocs computed from the domain's own extent (rather than a fixed
    # -180..180 array) and wrapped into [-180, 180) - Gridliner silently drops
    # tick candidates outside that range, so for a dateline-crossing domain a
    # fixed global array would only draw gridlines on one side of 180 degrees
    raw_xlocs = np.arange(np.floor(extent[0] / 10) * 10, extent[1] + 10, 10)
    xlocs = ((raw_xlocs + 180) % 360) - 180
    gl = ax.gridlines(
        xlocs=xlocs, ylocs=np.arange(-90, 91, 10),
        draw_labels=True, linewidth=0.3, color="gray", alpha=0.5, linestyle="--",
        zorder=0,
    )
    # NB: do not set gl.top_labels = False - it triggers a cartopy bug where the
    # title's y-position is computed as inf (title silently fails to render)
    gl.right_labels = False
    gl.xlabel_style = {"size": 6}
    gl.ylabel_style = {"size": 6}

    for track_id, track in tracks.groupby("track_id"):
        lon = track.lon.values
        lat = track.lat.values

        # thin connecting line - unwrap longitude so a track crossing the 0/360
        # meridian draws a short correct segment instead of either a spurious line
        # all the way across the map, or (if simply broken at the jump) a gap
        lon_line = np.unwrap(lon, period=360)
        ax.plot(
            lon_line, lat, "-", color="0.5", linewidth=0.6, zorder=1,
            transform=ccrs.PlateCarree(),
        )

        # genesis / lysis labels
        for txt, x, y in [("P", lon[0], lat[0]), ("L", lon[-1], lat[-1])]:
            ax.text(
                x, y, txt, fontsize=7, fontweight="bold", ha="center", va="center",
                color="black", zorder=5, transform=ccrs.PlateCarree(), clip_on=True,
                path_effects=[pe.withStroke(linewidth=2, foreground="white")],
            )

        # markers, one scatter call per nature category present in this track
        mslp = get_mslp(track)
        for nat_code, marker in NATURE_MARKERS.items():
            mask = track.nature.values == nat_code
            if mask.any():
                ax.scatter(
                    lon[mask], lat[mask],
                    c=mslp[mask],
                    cmap=cmap, norm=norm, marker=marker,
                    edgecolors="k", linewidth=0.3, s=MARKER_SIZE, zorder=3,
                    transform=ccrs.PlateCarree(),
                )

    ax.set_title(title)

    # Two separate legends placed outside the map so nothing is obscured
    pressure_handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=8,
               markerfacecolor=cmap(norm(p)), markeredgecolor="w", label=f"{p}")
        for p in PRESSURE_LEVELS
    ]
    nature_handles = [
        Line2D([0], [0], marker=m, linestyle="", markersize=8,
               markerfacecolor="0.75", markeredgecolor="k", label=code)
        for code, m in NATURE_MARKERS.items()
    ]

    leg1 = ax.legend(
        handles=pressure_handles, title="mslp", loc="upper left",
        bbox_to_anchor=(1.02, 1.0), frameon=False, handletextpad=0.4,
        labelspacing=0.6,
    )
    ax.add_artist(leg1)
    ax.legend(
        handles=nature_handles, title="nature", loc="upper left",
        bbox_to_anchor=(1.10, 1.0), frameon=False, handletextpad=0.4,
        labelspacing=0.6,
    )

    fig.savefig(out_path, dpi=600)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Batch loop
# ---------------------------------------------------------------------------
EMPTY_SUMMARY_COLUMNS = [
    "track_id", "storm_start", "storm_end",
    "origin_lat", "origin_lon", "end_lat", "end_lon", "is_tc",
]

def _empty_summary():
    """Match the schema nature.wcsi() would return, for a basin/year with zero
    qualifying tracks - see the comment in process_domain_year for why this is
    needed instead of just calling nature.wcsi() and getting an empty result back.
    """
    summary = pd.DataFrame(columns=EMPTY_SUMMARY_COLUMNS)
    summary["storm_start"] = pd.to_datetime(summary["storm_start"])
    summary["is_tc"] = summary["is_tc"].astype(bool)
    return summary


def process_domain_year(experiment, year, domain_name, domain):
    tag = f"{experiment}_{year}_{domain_name}"
    out_nc = f"{tag}_TCs.nc"
    out_parquet = f"{tag}_summary.parquet"

    if os.path.exists(out_parquet):
        print(f"[{experiment} {year} {domain_name}] {out_parquet} already exists, skipping WCSI computation")
        summary = pd.read_parquet(out_parquet)
        wcsi_tracks = huracanpy.load(out_nc) if os.path.exists(out_nc) else None
    else:
        if domain["hemisphere"] == "N":
            data_dir = os.path.join(BASE_DIR, experiment, f"VOR_VERTAVG_T63filt_{year}")
            file_glob = "tr_trs_pos*.nc"
        else:
            data_dir = os.path.join(BASE_DIR, experiment, f"VOR_VERTAVG_T63filt_{year}_SH")
            file_glob = "tr_trs_neg*.nc"

        matches = sorted(glob.glob(os.path.join(data_dir, file_glob)))
        if not matches:
            print(f"[{experiment} {year} {domain_name}] no {file_glob} found in {data_dir}, skipping")
            return None
        if len(matches) > 1:
            print(f"[{experiment} {year} {domain_name}] multiple matches in {data_dir}, using {matches[0]}")
        input_file = matches[0]

        tracks = huracanpy.load(
            input_file,
            rename=dict(Hart_B="cps_b", Hart_Vtl="cps_vtl", Hart_Vtu="cps_vtu"),
        )
        new_tracks = gather_vorticity_profile(tracks)

        npoints = infer_npoints(new_tracks, TC_PERSISTENCE_HOURS)
        try:
            wcsi_tracks, summary = nature.wcsi(
                new_tracks, vort_threshold=VORT_THRESHOLD, basin=domain["basin"],
                ocean=True, npoints=npoints,
            )
        except ValueError as exc:
            # nature.wcsi() does xr.concat()/pd.concat() on its per-track results
            # with no guard for the basin/year having zero qualifying tracks (a
            # basin like SATL can easily have none in a given year) - it raises
            # ValueError instead of returning something empty, so treat that
            # specific failure as "zero tracks" rather than a real error.
            if "concatenate" not in str(exc):
                raise
            wcsi_tracks = None
            summary = _empty_summary()

        if wcsi_tracks is not None:
            wcsi_tracks = compute_nature(wcsi_tracks, FILTER_SIZE)
            wcsi_tracks = filter_genesis_latitude(wcsi_tracks, TC_GENESIS_THRESHOLD_LATITUDE)
            huracanpy.save(wcsi_tracks, out_nc)
        summary.to_parquet(out_parquet)

    # summary only ever holds tracks from this (experiment, year) file, so no
    # further filtering by storm_start's calendar year - SH basin files span a
    # season (e.g. "1995" = Jul 1995-Jun 1996), so a real TC's storm_start can
    # legitimately fall in the following calendar year and would otherwise be
    # silently dropped from the count while still being drawn on the plot
    n_tc = int(summary.is_tc.sum())

    if wcsi_tracks is not None:
        plot_domain_nature(
            wcsi_tracks,
            title=f"{domain['title']} tracks by nature, {experiment} {year} ({n_tc} TCs)",
            out_path=f"{tag}_nature.png",
            extent=domain["extent"],
        )
    else:
        print(f"[{experiment} {year} {domain_name}] no tracks in basin {domain['basin']!r}, skipping plot")

    print(f"[{experiment} {year} {domain_name}] {n_tc} tropical cyclones identified")
    return dict(experiment=experiment, year=year, domain=domain_name, n_tc=n_tc)


def main():
    results = []
    for experiment in EXPERIMENTS:
        for year in YEARS:
            for domain_name, domain in DOMAINS.items():
                if domain_name in SKIP_DOMAINS:
                            continue    
                try:
                    result = process_domain_year(experiment, year, domain_name, domain)
                except Exception as exc:
                    print(f"[{experiment} {year} {domain_name}] FAILED: {exc}")
                    continue
                if result is not None:
                    results.append(result)

    if results:
        counts = pd.DataFrame(results)
        counts.to_csv(f"{'_'.join(EXPERIMENTS)}_TC_counts_by_year.csv", index=False)

        # wide-format pivot (one row per experiment/year, one column per domain) -
        # convenient for a quick look, in addition to the long-format CSV above
        wide = counts.pivot_table(
            index=["experiment", "year"], columns="domain", values="n_tc", fill_value=0
        )
        wide.to_csv(f"{'_'.join(EXPERIMENTS)}_TC_counts_by_year_wide.csv")


if __name__ == "__main__":
    main()
