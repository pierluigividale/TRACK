"""Method-6 WCSI tropical cyclone identification for TRACK output (self-contained).

Method 6 ("WCSI warm core only") is the identification used for the u-ch330 /
u-ck777 / u-ck778 / u-dc009 ensemble. For every experiment / year / hemisphere
this script:

  1. loads one TRACK output file (tr_trs_pos*.nc or tr_trs_neg*.nc, with Hart
     phase-space parameters, vorticity at several levels and the MSLP minimum
     added - the "hart.new" files);
  2. flags TC points with wcsi.wcsi() using the method-6 criteria below;
  3. classifies every point of each TC track (TC/BC/Tr/MV/Ot/Vo) with
     wcsi.nature();
  4. assigns each TC track to ONE ocean domain - the basin of its maximum
     vorticity during its TC phase, with the AUS and Central Pacific basins
     merged as IBTrACS counts them;
  5. writes, per domain, {experiment}_{year}_{domain}_TCs.nc (the TC tracks)
     and {experiment}_{year}_{domain}_summary.parquet (one row per track),
     plus optionally a map ({...}_nature.png) and counts CSVs for the run.

Those file names are what huracanpy_TD_ensemble.py, huracanpy_TC_interannual.py
and huracanpy_TC_genesis_maps.py read.

Method-6 criteria (fixed - see METHOD 6 block; change them only to test a
different method):
  - Hart CPS warm-core / symmetric criteria (wcsi defaults), no minimum
    850 hPa vorticity (VORT_THRESHOLD = 0);
  - warm-core vorticity test: 850 hPa minus upper-level relative vorticity
    > 6e-5 s**-1 (tcident's 6.0 in 1e-5 units);
  - the point is over ocean AND has a valid MSLP minimum (a cyclone, not just
    a vortex);
  - all of the above for 24 consecutive hours;
  - the track's first TC point equatorward of 50 degrees.

Replaces huracanpy_identify_batch_method6*.py, which reconfigured
huracanpy_identify_batch_wcsi_new.py (which itself borrowed from
huracanpy_identify_batch_all.py). Results are identical to that chain.

Requirements: huracanpy, wcsi from the current GitHub main branch
(pip install git+https://github.com/Huracan-project/wcsi - older releases
lack wcsi.wcsi(..., vort_warm_core_threshold=...)), numpy, pandas, xarray,
parse; cartopy + matplotlib only for --maps.

Usage:
  python huracanpy_identify_method6.py                      # settings below
  python huracanpy_identify_method6.py --experiments u-ck777 --years 2000-2001
  python huracanpy_identify_method6.py --out-dir /tmp/test --overwrite --maps

Already-processed experiment/year/domains (summary file present in OUT_DIR)
are SKIPPED unless --overwrite is given, so a re-run after changing settings
must use --overwrite or a new --out-dir.
"""

import argparse
import glob
import os
import re

import huracanpy
import numpy as np
import pandas as pd
import wcsi
from parse import parse

if not hasattr(wcsi, "wcsi") or "vort_warm_core_threshold" not in wcsi.wcsi.__code__.co_varnames:
    raise ImportError(
        f"The wcsi package at {wcsi.__file__} is too old for method 6 (needs "
        "wcsi.wcsi(..., vort_warm_core_threshold=...)). Install the current version: "
        "pip install --upgrade git+https://github.com/Huracan-project/wcsi  "
        "(or put a clone of it first on PYTHONPATH).")
import wcsi._nature  # noqa: E402

# ===========================================================================
# SETTINGS - edit for your data (all can also be given on the command line)
# ===========================================================================
DATA_DIR = "/Users/vidale/Python_scripts/TC-TRACKS"
OUT_DIR = "/Users/vidale/Python_scripts/TC-TRACKS/METHOD6_ENSEMBLE_RUN"
EXPERIMENTS = ["u-ch330", "u-ck777", "u-ck778", "u-dc009"]
YEARS = range(1980, 2015)  # last year exclusive

# Where each year's TRACK file is, relative to DATA_DIR, per hemisphere.
# {experiment} and {year} are filled in; glob wildcards are allowed. A missing
# file is reported and skipped. Set a hemisphere to None to skip it entirely.
INPUT_PATTERNS = {
    "N": "{experiment}/VOR_VERTAVG_T63filt_{year}/tr_trs_pos*.nc",
    "S": "{experiment}/VOR_VERTAVG_T63filt_{year}_SH/tr_trs_neg*.nc",
}

# Longitude of the MSLP minimum that TRACK stores next to the MSLP field
# (1e25, i.e. NaN once loaded, where it found no minimum). Its name depends on
# the order in which fields were added to the TRACK run (longitude_6 in our
# model files). "auto" picks the longitude_N variable immediately preceding
# the MSLP field and prints which one it used - check it once for new data.
MSLP_LON_VARIABLE = "auto"

# Upper level (hPa) for the warm-core test - the top relative_vorticity_{level}
# level in the files. Our files go up to 250 hPa (as tcident uses).
WARM_CORE_UPPER_LEVEL = 250

SKIP_DOMAINS = set()  # e.g. {"MED"}
MAKE_MAPS = False     # per-year-per-domain nature maps (slow; needs cartopy)

# ===========================================================================
# METHOD 6 - the identification criteria. Do not change for method-6 output.
# ===========================================================================
VORT_THRESHOLD = 0                # minimum 850 hPa vorticity, raw s**-1: none
VORT_WARM_CORE_THRESHOLD = 6e-5   # 850 minus upper-level vorticity, raw s**-1
COHERENT = False                  # vorticity defined at every level (wcsi default)
TC_PERSISTENCE_HOURS = 24         # criteria must hold this long, consecutively
TC_GENESIS_THRESHOLD_LATITUDE = 50.0  # first TC point must be equatorward of this
TRACK_FILL_VALUE = 1e25

# ===========================================================================
# Domains. "basin" is huracanpy's WMO-TC basin code.
# ===========================================================================
DOMAINS = {
    "NATL": dict(basin="NATL", hemisphere="N", extent=[-110, 40, -5, 75], title="North Atlantic"),
    "MED": dict(basin="MED", hemisphere="N", extent=[-20, 40, 20, 60], title="Mediterranean Sea"),
    "NWPAC": dict(basin="WNP", hemisphere="N", extent=[-265, -135, -10, 65], title="Northwest Pacific"),
    "SPAC": dict(basin="SP", hemisphere="S", extent=[130, 245, -65, 5], title="South Pacific"),
    "EPAC": dict(basin="ENP", hemisphere="N", extent=[-185, -65, -5, 40], title="East Pacific"),
    "NIO": dict(basin="NI", hemisphere="N", extent=[30, 100, -5, 30], title="North Indian Ocean"),
    "SIO": dict(basin="SI", hemisphere="S", extent=[20, 140, -65, 5], title="South Indian Ocean"),
    "SATL": dict(basin="SA", hemisphere="S", extent=[-65, 20, -50, 5], title="South Atlantic"),
}
# huracanpy basins with no domain of their own, merged as IBTrACS counts them:
# CP -> EPAC; AUS -> SIO west of AUS_SPLIT_LON (by longitude of peak), else SPAC.
AUS_SPLIT_LON = 135

SUMMARY_COLUMNS = ["track_id", "storm_start", "storm_end", "origin_lat", "origin_lon",
                   "end_lat", "end_lon", "is_tc"]
MSLP_NAMES = ("air_pressure_at_mean_sea_level", "psl")


# ---------------------------------------------------------------------------
# Workaround for a bug in wcsi._nature.smooth_excursions (empty sequences ->
# "need at least one array to concatenate"); partly fixed upstream. Same
# function with the empty cases guarded, otherwise unchanged, so results are
# identical either way. Delete once the installed wcsi is fully fixed.
# ---------------------------------------------------------------------------
def _smooth_excursions(track_id, criteria, *, min_count, starts=None, ends=None):
    if starts is None:
        starts = np.where(~criteria[:-1] & criteria[1:])[0] + 1
    if ends is None:
        ends = np.where(criteria[:-1] & ~criteria[1:])[0]

    if len(starts) == 0 and len(ends) == 0:
        if not criteria.any():
            return criteria.copy(), starts, ends
        return criteria.copy(), np.array([0]), np.array([len(criteria) - 1])
    if len(ends) and (len(starts) == 0 or starts[0] > ends[0]):
        starts = np.concat([[0], starts])
    if len(starts) > len(ends):
        ends = np.concat([ends, [len(criteria) - 1]])

    close = ((starts[1:] - ends[:-1]) < min_count) & (track_id[starts[1:]] == track_id[ends[:-1]])
    close_idx = np.where(close)[0] + 1

    joined = [np.arange(start, end + 1) for start, end in zip(starts[close_idx - 1], ends[close_idx])]
    idx = np.concat(joined) if joined else np.array([], dtype=int)

    starts = starts[np.concat([[True], ~close])]
    ends = ends[np.concat([~close, [True]])]

    new_criteria = criteria.copy()
    new_criteria[idx] = True
    return new_criteria, starts, ends


wcsi._nature.smooth_excursions = _smooth_excursions


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------
def find_input(experiment, year, hemisphere):
    pattern = INPUT_PATTERNS.get(hemisphere)
    if pattern is None:
        return None
    full = os.path.join(DATA_DIR, pattern.format(experiment=experiment, year=year))
    matches = sorted(glob.glob(full))
    if not matches:
        print(f"[{experiment} {year} {hemisphere}H] no input file matching {full} - skipping")
        return None
    if len(matches) > 1:
        print(f"[{experiment} {year} {hemisphere}H] {len(matches)} files match {full}, using {matches[0]}")
    return matches[0]


def load_tracks(path):
    """Load a TRACK file; gather relative_vorticity_{level} variables into one
    relative_vorticity(record, pressure) variable, keeping the tracked T63
    vorticity as "vorticity".
    """
    tracks = huracanpy.load(path, rename=dict(Hart_B="cps_b", Hart_Vtl="cps_vtl", Hart_Vtu="cps_vtu"))
    tracks["vorticity"] = tracks["relative_vorticity"]
    plevs = [int(r.named["n"]) for r in (parse("relative_vorticity_{n:d}", v) for v in tracks) if r is not None]
    if plevs:
        tracks["pressure"] = ("pressure", plevs)
        tracks = tracks.set_coords("pressure")
        vort = np.stack([tracks[f"relative_vorticity_{p}"].values for p in plevs], axis=1)
        tracks = tracks.drop_vars([f"relative_vorticity_{p}" for p in plevs])
        tracks["relative_vorticity"] = (["record", "pressure"], vort)
    return tracks


def mslp_lon_variable(tracks, path):
    if MSLP_LON_VARIABLE != "auto":
        return MSLP_LON_VARIABLE
    names = list(tracks.data_vars)
    for mslp in MSLP_NAMES:
        if mslp in names:
            i = names.index(mslp)
            if i >= 2 and re.fullmatch(r"longitude_\d+", names[i - 2]) and \
                    re.fullmatch(r"latitude_\d+", names[i - 1]):
                return names[i - 2]
    raise ValueError(
        f"{path}: could not work out which longitude_N is the MSLP minimum's "
        f"position (expected longitude_N, latitude_N, then one of {MSLP_NAMES}). "
        f"Set MSLP_LON_VARIABLE explicitly.")


def infer_npoints(tracks, hours):
    """Consecutive records spanning `hours`, from the first track's timestep."""
    first = tracks.track_id.values == tracks.track_id.values[0]
    step_h = np.diff(tracks.time.values[first])[0] / np.timedelta64(1, "h")
    return max(1, round(hours / step_h))


# ---------------------------------------------------------------------------
# Identification
# ---------------------------------------------------------------------------
def identify(tracks, path, label):
    """Return (tc_tracks with is_tc and nature added, basin per tc_tracks record),
    or (None, None) if no TC.
    """
    track_id = tracks.track_id.values
    lon_var = mslp_lon_variable(tracks, path)
    mslp_lon = tracks[lon_var].values
    no_mslp_min = ~np.isfinite(mslp_lon) | (mslp_lon >= TRACK_FILL_VALUE * 0.1)
    is_ocean = tracks.hrcn.get_is_ocean().values
    basin = tracks.hrcn.get_basin().values
    npoints = infer_npoints(tracks, TC_PERSISTENCE_HOURS)
    print(f"[{label}] {os.path.basename(path)}: {tracks.track_id.to_series().nunique()} tracks, "
          f"{100 * no_mslp_min.mean():.1f}% of points without an MSLP minimum ({lon_var}), "
          f"{npoints} points = {TC_PERSISTENCE_HOURS} h")
    if not no_mslp_min.any():
        print(f"[{label}] WARNING: every point has an MSLP minimum - is {lon_var} the right variable?")

    vort = tracks.relative_vorticity
    levels = list(vort.pressure.values)
    if 850 not in levels or WARM_CORE_UPPER_LEVEL not in levels:
        raise ValueError(f"{path}: need relative_vorticity at 850 and {WARM_CORE_UPPER_LEVEL} hPa, have {levels}")
    if WARM_CORE_UPPER_LEVEL != 200:
        # wcsi's warm-core test uses the level labelled 200 hPa
        vort = vort.assign_coords(pressure=[200 if p == WARM_CORE_UPPER_LEVEL else p for p in levels])

    is_tc = wcsi.wcsi(
        track_id, tracks.cps_b.values, tracks.cps_vtl.values, tracks.cps_vtu.values, vort,
        filter_all=is_ocean & ~no_mslp_min,
        npoints=npoints,
        vort_threshold=VORT_THRESHOLD,
        vort_warm_core_threshold=VORT_WARM_CORE_THRESHOLD,
        coherent=COHERENT,
    )
    tc_mask = np.isin(track_id, np.unique(track_id[is_tc]))
    if not tc_mask.any():
        return None, None
    tc_tracks = tracks.isel(record=tc_mask)
    tc_tracks["is_tc"] = ("record", is_tc[tc_mask])
    tc_tracks["nature"] = ("record", wcsi.nature(
        tc_tracks.track_id.values, tc_tracks.cps_b.values, tc_tracks.cps_vtl.values,
        tc_tracks.cps_vtu.values, no_mslp_min[tc_mask], is_tc[tc_mask],
    ))
    return tc_tracks, basin[tc_mask]


def gap_basin_domain(basin, apex_lon):
    if basin == "CP":
        return "EPAC"
    if basin == "AUS":
        return "SIO" if apex_lon % 360 < AUS_SPLIT_LON else "SPAC"
    return None


def assign_domains(tc_tracks, basin):
    """track_id -> domain: basin of the track's maximum |vorticity| during its
    TC phase (gap basins merged as IBTrACS does); NaN if in no domain.
    """
    points = pd.DataFrame(dict(
        track_id=tc_tracks.track_id.values, vorticity=np.abs(tc_tracks.vorticity.values),
        basin=basin, lon=tc_tracks.lon.values, is_tc=tc_tracks.is_tc.values,
    ))
    tc = points[points.is_tc]
    apex = tc.loc[tc.groupby("track_id").vorticity.idxmax()].set_index("track_id")
    domain_of = apex.basin.map({d["basin"]: name for name, d in DOMAINS.items()})
    gap = domain_of.isna()
    domain_of[gap] = [gap_basin_domain(b, x) for b, x in zip(apex.basin[gap], apex.lon[gap])]
    return domain_of


def summarize(tracks):
    points = pd.DataFrame(dict(
        track_id=tracks.track_id.values, time=pd.to_datetime(tracks.time.values),
        lat=tracks.lat.values, lon=tracks.lon.values, is_tc=tracks.is_tc.values,
    ))
    g = points.groupby("track_id", sort=False)
    return pd.DataFrame(dict(
        storm_start=g.time.first(), storm_end=g.time.last(),
        origin_lat=g.lat.first(), origin_lon=g.lon.first(),
        end_lat=g.lat.last(), end_lon=g.lon.last(), is_tc=g.is_tc.any(),
    )).reset_index()


def filter_genesis_latitude(tracks, threshold_latitude):
    """Keep tracks whose first TC-nature point is equatorward of the threshold."""
    gen = tracks.isel(record=tracks.nature == "TC").hrcn.get_gen_vals()
    keep = gen.isel(track_id=np.abs(gen.lat) < threshold_latitude)
    if keep.sizes["track_id"] == 0:
        return None
    return tracks.hrcn.sel_id(keep.track_id)


def process_file(experiment, year, hemisphere, path, domains, out_dir, make_maps):
    label = f"{experiment} {year} {hemisphere}H"
    tc_tracks, basin = identify(load_tracks(path), path, label)
    domain_of = assign_domains(tc_tracks, basin) if tc_tracks is not None else pd.Series(dtype=object)

    results = []
    for name, domain in domains.items():
        tag = f"{experiment}_{year}_{name}"
        ids = domain_of.index[domain_of == name].values
        dtracks = tc_tracks.hrcn.sel_id(ids) if len(ids) else None
        summary = summarize(dtracks) if dtracks is not None else pd.DataFrame(columns=SUMMARY_COLUMNS)
        n_tc = len(summary)
        tc_file = os.path.join(out_dir, f"{tag}_TCs.nc")
        if dtracks is not None:
            dtracks = filter_genesis_latitude(dtracks, TC_GENESIS_THRESHOLD_LATITUDE)
        if os.path.exists(tc_file):
            os.remove(tc_file)  # --overwrite: don't leave a stale file if now empty
        if dtracks is not None:
            huracanpy.save(dtracks, tc_file)
            if make_maps:
                plot_domain_nature(dtracks, f"{domain['title']} tracks by nature, {experiment} {year} ({n_tc} TCs)",
                                   os.path.join(out_dir, f"{tag}_nature.png"), domain["extent"])
        # summary last: its presence marks this experiment/year/domain as done
        summary.to_parquet(os.path.join(out_dir, f"{tag}_summary.parquet"))
        print(f"[{experiment} {year} {name}] {n_tc} tropical cyclones")
        results.append(dict(experiment=experiment, year=year, domain=name, n_tc=n_tc))
    n_out = int(domain_of.isna().sum())
    if n_out:
        print(f"[{label}] {n_out} TC tracks peak outside every domain (not saved)")
    return results


# ---------------------------------------------------------------------------
# Map (optional)
# ---------------------------------------------------------------------------
NATURE_MARKERS = {"TC": "o", "BC": "^", "Tr": "s", "MV": "D", "Ot": "X", "Vo": "P"}
PRESSURE_LEVELS = [920, 930, 940, 950, 960, 970, 980, 990, 1000, 1010, 1020]


def plot_domain_nature(tracks, title, out_path, extent):
    """Whole lifetime of each track: marker shape = nature, colour = MSLP (hPa),
    P/L = first/last point. Same figure as huracanpy_identify_batch_all.py.
    """
    import cartopy.crs as ccrs
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.patheffects as pe
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    central = 0.5 * (extent[0] + extent[1])
    cmap = plt.get_cmap("terrain")
    norm = plt.Normalize(vmin=min(PRESSURE_LEVELS), vmax=max(PRESSURE_LEVELS))

    # axes sized to the extent's aspect ratio, legends to the right (avoids
    # cartopy letterboxing and tight-bbox dropping the legends)
    map_w, left, bottom, top, legend_w = 9.0, 0.5, 0.4, 0.5, 3.2
    map_h = map_w * (extent[3] - extent[2]) / (extent[1] - extent[0])
    fig_w, fig_h = left + map_w + legend_w, bottom + map_h + top
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([left / fig_w, bottom / fig_h, map_w / fig_w, map_h / fig_h],
                      projection=ccrs.PlateCarree(central_longitude=central))
    # xlim/ylim in the rotated frame rather than set_extent (which breaks for
    # extents crossing the antimeridian, e.g. SPAC)
    ax.set_xlim(extent[0] - central, extent[1] - central)
    ax.set_ylim(extent[2], extent[3])
    ax.coastlines()
    xlocs = ((np.arange(np.floor(extent[0] / 10) * 10, extent[1] + 10, 10) + 180) % 360) - 180
    gl = ax.gridlines(xlocs=xlocs, ylocs=np.arange(-90, 91, 10), draw_labels=True,
                      linewidth=0.3, color="gray", alpha=0.5, linestyle="--", zorder=0)
    gl.right_labels = False  # (gl.top_labels = False breaks the title in cartopy)
    gl.xlabel_style = gl.ylabel_style = {"size": 6}

    pc = ccrs.PlateCarree()
    for _, track in tracks.groupby("track_id"):
        lon, lat = track.lon.values, track.lat.values
        ax.plot(np.unwrap(lon, period=360), lat, "-", color="0.5", linewidth=0.6, zorder=1, transform=pc)
        for txt, x, y in [("P", lon[0], lat[0]), ("L", lon[-1], lat[-1])]:
            ax.text(x, y, txt, fontsize=7, fontweight="bold", ha="center", va="center", zorder=5,
                    transform=pc, clip_on=True, path_effects=[pe.withStroke(linewidth=2, foreground="white")])
        mslp = next(track[n].values for n in MSLP_NAMES if n in track.variables)
        for code, marker in NATURE_MARKERS.items():
            m = track.nature.values == code
            if m.any():
                ax.scatter(lon[m], lat[m], c=mslp[m], cmap=cmap, norm=norm, marker=marker,
                           edgecolors="k", linewidth=0.3, s=10, zorder=3, transform=pc)
    ax.set_title(title)
    leg = ax.legend(handles=[Line2D([0], [0], marker="o", linestyle="", markersize=8,
                                    markerfacecolor=cmap(norm(p)), markeredgecolor="w", label=f"{p}")
                             for p in PRESSURE_LEVELS],
                    title="mslp", loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    ax.add_artist(leg)
    ax.legend(handles=[Line2D([0], [0], marker=m, linestyle="", markersize=8, markerfacecolor="0.75",
                              markeredgecolor="k", label=c) for c, m in NATURE_MARKERS.items()],
              title="nature", loc="upper left", bbox_to_anchor=(1.10, 1.0), frameon=False)
    fig.savefig(out_path, dpi=600)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def parse_years(text):
    """'1980-2014' (inclusive) or '1990,1995,2000'."""
    if "-" in text:
        a, b = text.split("-")
        return range(int(a), int(b) + 1)
    return [int(y) for y in text.split(",")]


def main():
    global DATA_DIR
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--experiments", nargs="+", default=EXPERIMENTS)
    p.add_argument("--years", type=parse_years, default=YEARS, help="e.g. 1980-2014 (inclusive) or 1990,2000")
    p.add_argument("--data-dir", default=DATA_DIR)
    p.add_argument("--out-dir", default=OUT_DIR)
    p.add_argument("--hemispheres", nargs="+", default=["N", "S"], choices=["N", "S"])
    p.add_argument("--overwrite", action="store_true", help="recompute even if outputs exist")
    p.add_argument("--maps", action="store_true", default=MAKE_MAPS, help="also write nature maps")
    args = p.parse_args()
    DATA_DIR = args.data_dir
    os.makedirs(args.out_dir, exist_ok=True)

    print(f"Method 6: vort_threshold={VORT_THRESHOLD}, warm core {VORT_WARM_CORE_THRESHOLD} "
          f"(850-{WARM_CORE_UPPER_LEVEL} hPa), closed MSLP + ocean, {TC_PERSISTENCE_HOURS} h, "
          f"genesis < {TC_GENESIS_THRESHOLD_LATITUDE} deg; output -> {args.out_dir}")

    results, skipped, failed = [], [], []
    for experiment in args.experiments:
        for year in args.years:
            for hemisphere in args.hemispheres:
                domains = {n: d for n, d in DOMAINS.items()
                           if d["hemisphere"] == hemisphere and n not in SKIP_DOMAINS}
                if not args.overwrite:
                    done = {n for n in domains
                            if os.path.exists(os.path.join(args.out_dir, f"{experiment}_{year}_{n}_summary.parquet"))}
                    if done == set(domains):
                        skipped.append(f"{experiment} {year} {hemisphere}H")
                        continue
                path = find_input(experiment, year, hemisphere)
                if path is None:
                    continue
                try:
                    results.extend(process_file(experiment, year, hemisphere, path, domains,
                                                args.out_dir, args.maps))
                except Exception as exc:
                    print(f"[{experiment} {year} {hemisphere}H] FAILED: {type(exc).__name__}: {exc}")
                    failed.append(f"{experiment} {year} {hemisphere}H")

    if results:
        counts = pd.DataFrame(results)
        stem = os.path.join(args.out_dir, f"{'_'.join(args.experiments)}_TC_counts_by_year")
        counts.to_csv(stem + ".csv", index=False)
        counts.pivot_table(index=["experiment", "year"], columns="domain", values="n_tc",
                           fill_value=0).to_csv(stem + "_wide.csv")
    print(f"\nProcessed {len({(r['experiment'], r['year']) for r in results})} experiment-years.")
    if skipped:
        print(f"SKIPPED {len(skipped)} experiment/year/hemispheres already in {args.out_dir} "
              f"(use --overwrite or a new --out-dir to recompute), e.g. {skipped[0]}")
    if failed:
        print(f"FAILED {len(failed)}: {', '.join(failed)}")


if __name__ == "__main__":
    main()
