"""Observational counterpart to huracanpy_identify_batch_all.py, using IBTrACS.

Produces the same style of per-domain, per-year "tracks by nature" maps as
huracanpy_identify_batch_all.py, but from IBTrACS observations instead of model
tracker output. IBTrACS does not carry the Hart cyclone-phase-space fields the
model script uses to run the WCSI classification, so there is no TC/BC/Tr/MV/Ot/Vo
split here. Instead each point is plotted using IBTrACS's own reported storm
nature: TS (tropical), SS (subtropical), ET (extratropical), DS (disturbance),
MX (mixed - agencies disagreed). "NR" (not reported) points are skipped, same as
any point whose nature isn't a plotted category.

By default this nature (and the matching lat/lon/pressure/wind) comes from
IBTrACS's consensus fields, i.e. whichever agency is officially responsible for
that basin. AGENCY_OVERRIDES below switches specific basins to a different
agency's own track instead - currently the South Indian Ocean ("SI"), sourced
from JTWC instead of BOM/La Reunion per the peer-reviewed literature. A
JTWC-sourced basin reports pressure less consistently than the consensus
fields do, so its mslp colour-coding will look sparser than other basins' -
that's a real gap in JTWC's own reporting, not a bug here.

Requires the classic IBTrACS "since1980" NetCDF (storm x date_time layout, with
the `nature` variable) - huracanpy.load(source="ibtracs") does NOT provide this,
its offline subsets only carry track_id/season/basin/time/lon/lat/wind/slp.
Download from:
https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/netcdf/IBTrACS.since1980.v04r01.nc
"""

import os

import cartopy.crs as ccrs
import matplotlib
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.lines import Line2D

matplotlib.use("Agg")  # no interactive display in batch mode

# ---------------------------------------------------------------------------
# Configuration - edit before each run
# ---------------------------------------------------------------------------
IBTRACS_FILE = "/Users/vidale/Python_scripts/IBTrACS.since1980.v04r01.nc"
YEARS = range(1980, 2014)  # matches the model script's default YEARS

# One entry per domain to process. "basin" is IBTrACS's own reported basin code
# (see the `basin` variable in the IBTrACS NetCDF: EP/NA/NI/SA/SI/SP/WP) - NOT
# the huracanpy/TRACK basin convention (NATL/WNP/ENP/...) used by
# huracanpy_identify_batch_all.py's DOMAINS. A track is assigned to a domain by
# the basin reported at its point of peak intensity (lowest wmo_pres, or highest
# wmo_wind if pressure was never reported), mirroring that script's "tracks with
# their maximum intensity in that basin" filtering.
DOMAINS = {
    "NATL": dict(basin="NA", extent=[-110, 40, -5, 75], title="North Atlantic"),
    "NWPAC": dict(basin="WP", extent=[-265, -135, -10, 65], title="Northwest Pacific"),
    "SPAC": dict(basin="SP", extent=[140, 245, -65, 5], title="South Pacific"),
    "EPAC": dict(basin="EP", extent=[-160, -65, -5, 40], title="East Pacific"),
    "NIO": dict(basin="NI", extent=[30, 100, -5, 30], title="North Indian Ocean"),
    "SIO": dict(basin="SI", extent=[20, 140, -65, 5], title="South Indian Ocean"),
    "SATL": dict(basin="SA", extent=[-65, 20, -50, 5], title="South Atlantic"),
    # MED is intentionally absent: IBTrACS's basin field has no separate
    # Mediterranean code, so medicanes can't be isolated the same way here.
}

SKIP_DOMAINS = set()  # IBTrACS has ample storms in every domain/year, unlike the
                       # sparse model runs, so nothing is skipped by default

# IBTrACS's generic lat/lon/pressure/wind/nature fields are a per-point
# consensus assembled from whichever agency is officially responsible for that
# basin (see the `wmo_agency` variable) - for the South Indian Ocean that's
# BOM/La Reunion. The peer-reviewed literature (e.g. Knapp et al.) considers
# JTWC's independent reanalysis (the `usa_*` fields, `usa_agency="jtwc_sh"`
# there) more reliable for that basin, so it's used here instead. Keyed by
# IBTrACS's own 2-letter basin code (DOMAINS' "basin" values above), not by
# domain_name, since the choice of reporting agency is inherently about the
# basin, and this is the natural place to add more overrides basin-by-basin.
#AGENCY_OVERRIDES = {
#    "SI": "jtwc",
#}
AGENCY_OVERRIDES = {
    "SI": "reunion",
}
DEFAULT_AGENCY = "wmo"

# usa_status (JTWC's own storm-status classification) uses a different
# vocabulary than IBTrACS's consensus `nature` field, so when a basin is
# switched to "jtwc" this maps it onto the same TS/SS/ET/DS/MX scheme used
# everywhere else. Statuses with no clear equivalent (inland, extrapolated
# position, unknown, dissipating, unreported) are left unmapped - those points
# simply won't be drawn, the same as an "NR" (not reported) point isn't.
USA_STATUS_TO_NATURE = {
    "TD": "TS", "TS": "TS", "TY": "TS", "ST": "TS", "TC": "TS", "HU": "TS", "HR": "TS",
    "SD": "SS", "SS": "SS",
    "EX": "ET", "PT": "ET",
    "DB": "DS", "WV": "DS", "LO": "DS", "MD": "DS",
}

NATURE_MARKERS = {"TS": "o", "SS": "^", "ET": "s", "DS": "D", "MX": "X"}
# mslp legend levels (hPa) - also sets the colour scale range; kept identical to
# huracanpy_identify_batch_all.py so obs and model maps are visually comparable
PRESSURE_LEVELS = [920, 930, 940, 950, 960, 970, 980, 990, 1000, 1010, 1020]
MARKER_SIZE = 10  # main-plot marker size (points**2); legend markers are separate


# ---------------------------------------------------------------------------
# IBTrACS loading
# ---------------------------------------------------------------------------
def load_ibtracs(path):
    """Load the classic IBTrACS NetCDF and flatten it to one row per track point.

    The raw file stores each storm as a fixed-length (storm, date_time) row
    padded with NaNs, doesn't decode its byte-string fields, and includes
    provisional/spur tracks - this returns a tidy per-record DataFrame with
    those handled: real (non-NaN) points only, "main" tracks only (provisional
    and spur-merge/split/other tracks dropped), and str instead of bytes.
    """
    ds = xr.open_dataset(path, decode_times=False)

    is_main = ds.track_type.values.astype(str) == "main"
    ds = ds.isel(storm=is_main)

    sid = ds.sid.values.astype(str)
    season = ds.season.values.astype(int)
    lat = ds.lat.values
    lon = ds.lon.values
    # time units: days since 1858-11-17 (see the file's `time` variable attrs)
    time_raw = ds.time.values
    time = (
        pd.Timestamp("1858-11-17")
        + pd.to_timedelta(time_raw.ravel(), unit="D")
    ).values.reshape(time_raw.shape)
    pres = ds.wmo_pres.values
    wind = ds.wmo_wind.values
    nature = ds.nature.values.astype(str)
    basin = ds.basin.values.astype(str)

    # JTWC's own reanalysis, for basins in AGENCY_OVERRIDES (see there for why).
    # units match the generic fields (mb, kts), so no conversion is needed.
    usa_lat = ds.usa_lat.values
    usa_lon = ds.usa_lon.values
    usa_pres = ds.usa_pres.values
    usa_wind = ds.usa_wind.values
    usa_nature = np.vectorize(lambda s: USA_STATUS_TO_NATURE.get(s, "NR"))(
        ds.usa_status.values.astype(str)
    )

    valid = ~np.isnan(lat)

    # Assign each track a single domain-basin at its point of peak intensity
    # (lowest pressure; falls back to highest wind, then to the first valid
    # point, for the storms IBTrACS never got an intensity estimate for). This
    # always uses the generic (consensus) fields, regardless of AGENCY_OVERRIDES
    # below - it's just a geographic label for which domain a storm belongs to,
    # not an intensity/nature judgement call.
    pres_for_argmin = np.where(valid, pres, np.inf)
    has_pres = np.isfinite(pres_for_argmin).any(axis=1)
    wind_for_argmax = np.where(valid, wind, -np.inf)
    has_wind = np.isfinite(wind_for_argmax).any(axis=1)
    first_valid = valid.argmax(axis=1)

    peak_idx = first_valid.copy()
    peak_idx[has_wind] = np.nanargmax(wind_for_argmax[has_wind], axis=1)
    peak_idx[has_pres] = np.nanargmin(pres_for_argmin[has_pres], axis=1)
    track_basin = basin[np.arange(len(peak_idx)), peak_idx]

    n_storms = len(sid)
    rows = []
    for i in range(n_storms):
        agency = AGENCY_OVERRIDES.get(track_basin[i], DEFAULT_AGENCY)
        if agency == "jtwc":
            i_lat, i_lon, i_pres, i_wind, i_nature = (
                usa_lat[i], usa_lon[i], usa_pres[i], usa_wind[i], usa_nature[i]
            )
        else:
            i_lat, i_lon, i_pres, i_wind, i_nature = (
                lat[i], lon[i], pres[i], wind[i], nature[i]
            )

        mask = ~np.isnan(i_lat)
        n = int(mask.sum())
        rows.append(pd.DataFrame(dict(
            track_id=np.repeat(sid[i], n),
            season=np.repeat(season[i], n),
            domain_basin=np.repeat(track_basin[i], n),
            time=time[i][mask],
            lat=i_lat[mask],
            lon=i_lon[mask],
            psl=i_pres[mask],
            wind=i_wind[mask],
            nature=i_nature[mask],
        )))
    return pd.concat(rows, ignore_index=True)


# ---------------------------------------------------------------------------
# Plotting (same approach as huracanpy_identify_batch_all.py's plot_domain_nature)
# ---------------------------------------------------------------------------
def plot_domain_nature(tracks, title, out_path, extent):
    """Full lifetime, marker shape by nature category, colour = mslp.

    P/L mark each track's first/last point, thin lines join each track's points.
    """
    domain_crs = ccrs.PlateCarree(
        central_longitude=0.5 * (extent[0] + extent[1])
    )

    cmap = plt.get_cmap("terrain")
    norm = plt.Normalize(vmin=min(PRESSURE_LEVELS), vmax=max(PRESSURE_LEVELS))

    lon_span = extent[1] - extent[0]
    lat_span = extent[3] - extent[2]
    map_aspect = lon_span / lat_span

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
    central_longitude = 0.5 * (extent[0] + extent[1])
    ax.set_xlim(extent[0] - central_longitude, extent[1] - central_longitude)
    ax.set_ylim(extent[2], extent[3])
    ax.coastlines()

    raw_xlocs = np.arange(np.floor(extent[0] / 10) * 10, extent[1] + 10, 10)
    xlocs = ((raw_xlocs + 180) % 360) - 180
    gl = ax.gridlines(
        xlocs=xlocs, ylocs=np.arange(-90, 91, 10),
        draw_labels=True, linewidth=0.3, color="gray", alpha=0.5, linestyle="--",
        zorder=0,
    )
    gl.right_labels = False
    gl.xlabel_style = {"size": 6}
    gl.ylabel_style = {"size": 6}

    for track_id, track in tracks.groupby("track_id"):
        track = track.sort_values("time")
        lon = track.lon.values
        lat = track.lat.values

        lon_line = np.unwrap(lon, period=360)
        ax.plot(
            lon_line, lat, "-", color="0.5", linewidth=0.6, zorder=1,
            transform=ccrs.PlateCarree(),
        )

        for txt, x, y in [("P", lon[0], lat[0]), ("L", lon[-1], lat[-1])]:
            ax.text(
                x, y, txt, fontsize=7, fontweight="bold", ha="center", va="center",
                color="black", zorder=5, transform=ccrs.PlateCarree(), clip_on=True,
                path_effects=[pe.withStroke(linewidth=2, foreground="white")],
            )

        mslp = track.psl.values
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
def summarize_tracks(tracks):
    rows = []
    for track_id, track in tracks.groupby("track_id"):
        track = track.sort_values("time")
        rows.append(dict(
            track_id=track_id,
            storm_start=track.time.iloc[0],
            storm_end=track.time.iloc[-1],
            origin_lat=track.lat.iloc[0],
            origin_lon=track.lon.iloc[0],
            end_lat=track.lat.iloc[-1],
            end_lon=track.lon.iloc[-1],
            min_pressure=track.psl.min(),
            max_wind=track.wind.max(),
            natures=",".join(sorted(set(track.nature) & set(NATURE_MARKERS))),
            is_ts=bool((track.nature == "TS").any()),
        ))
    return pd.DataFrame(rows)


def process_domain_year(year, domain_name, domain, all_tracks):
    tag = f"IBTrACS_{year}_{domain_name}"

    sel = all_tracks[
        (all_tracks.season == year) & (all_tracks.domain_basin == domain["basin"])
    ]
    n_storms = sel.track_id.nunique()

    if n_storms == 0:
        print(f"[IBTrACS {year} {domain_name}] no tracks in basin {domain['basin']!r}, skipping")
        return dict(year=year, domain=domain_name, n_storms=0, n_ts=0)

    summary = summarize_tracks(sel)
    summary.to_parquet(f"{tag}_summary.parquet")
    n_ts = int(summary.is_ts.sum())

    plot_domain_nature(
        sel,
        title=f"{domain['title']} tracks by nature, IBTrACS {year} ({n_ts} TS)",
        out_path=f"{tag}_nature.png",
        extent=domain["extent"],
    )

    print(f"[IBTrACS {year} {domain_name}] {n_storms} storms, {n_ts} reaching TS")
    return dict(year=year, domain=domain_name, n_storms=n_storms, n_ts=n_ts)


def main():
    print(f"Loading {IBTRACS_FILE} ...")
    all_tracks = load_ibtracs(IBTRACS_FILE)

    results = []
    for year in YEARS:
        for domain_name, domain in DOMAINS.items():
            if domain_name in SKIP_DOMAINS:
                continue
            try:
                result = process_domain_year(year, domain_name, domain, all_tracks)
            except Exception as exc:
                print(f"[IBTrACS {year} {domain_name}] FAILED: {exc}")
                continue
            if result is not None:
                results.append(result)

    if results:
        counts = pd.DataFrame(results)
        counts.to_csv("IBTrACS_TC_counts_by_year.csv", index=False)

        wide = counts.pivot_table(
            index="year", columns="domain", values="n_ts", fill_value=0
        )
        wide.to_csv("IBTrACS_TC_counts_by_year_wide.csv")


if __name__ == "__main__":
    main()
