"""Batch TC identification with the new wcsi API (closed-MSLP version).

Same job as huracanpy_identify_batch_all.py - per experiment/year/domain,
identify tropical cyclones in the TRACK output with the WCSI criteria,
classify every point's nature (TC/BC/Tr/MV/Ot/Vo), and save the TC tracks,
a per-track summary, a map and a counts CSV - but using the updated wcsi
package (wcsi.wcsi / wcsi.nature, following Leo Saffin's
wcsi_new_example.py) instead of the old wcsi.nature.wcsi:

  - No minimum-vorticity threshold by default (VORT_THRESHOLD = None).
    Instead a WCSI point must have a valid MSLP minimum, i.e. be a cyclone
    rather than just a vortex: TRACK stores the MSLP-minimum position as
    1e25 where it found none (MSLP_LON_VARIABLE), which xarray decodes to NaN
    when loading the model files - both are treated as "no minimum".
  - Extra criteria go in through filter_all: a WCSI point must be in the
    domain's basin (per point, huracanpy's WMO-TC convention), over ocean,
    and have a valid MSLP minimum - for npoints consecutive points.
  - "Vo" in the nature tag now means "no valid MSLP minimum".

Differences from huracanpy_identify_batch_all.py worth knowing:

  - Basin: the old script assigned each track to the basin of its maximum
    vorticity and then required an is_tc point in that basin. Here the
    basin is a per-point criterion, so a track that is a TC in two basins
    counts in both.
  - coherent (vorticity defined at every level) defaults to False in the new
    wcsi and in Leo's example; the old run used True. Set COHERENT to match
    the old run if you want to isolate the effect of the closed-MSLP change.
  - Each year/hemisphere file is loaded once and all its domains processed
    from it, rather than reloaded per domain.
  - Outputs go to OUT_DIR (a separate folder), not the working directory, so
    the existing NEW_WCSI_RUN results are untouched. The file names and the
    summary/_TCs.nc contents are the same as before (is_tc, nature, ...), so
    huracanpy_TC_interannual.py and huracanpy_TC_genesis_maps.py can read them
    by pointing their BASE_DIR at OUT_DIR.

Kept from the old pipeline: the same input file selection, the 24-hour
persistence (TC_PERSISTENCE_HOURS) and the genesis-latitude filter on the
saved TC tracks (TC_GENESIS_THRESHOLD_LATITUDE).

The script prints the fraction of points without a valid MSLP minimum for
every file it loads (typically about a quarter in the model hart.new files).
If that is ever 0%, the closed-MSLP criterion is doing nothing - and with no
vorticity threshold there would then be no intensity criterion at all.
"""

import glob
import os
import sys

# Use the updated wcsi from the local clone without replacing the installed
# (old) version, which huracanpy_identify_batch_all.py still needs. Set to None
# once the new version is installed.
WCSI_SOURCE = "/Users/vidale/Python_scripts/wcsi"
if WCSI_SOURCE is not None:
    sys.path.insert(0, WCSI_SOURCE)

import huracanpy  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import wcsi  # noqa: E402
import wcsi._nature  # noqa: E402

from huracanpy_identify_batch_all import (  # noqa: E402
    BASE_DIR,
    DOMAINS,
    TC_GENESIS_THRESHOLD_LATITUDE,
    TC_PERSISTENCE_HOURS,
    filter_genesis_latitude,
    gather_vorticity_profile,
    infer_npoints,
    plot_domain_nature,
)

# ---------------------------------------------------------------------------
# Configuration - edit before each run
# ---------------------------------------------------------------------------
OUT_DIR = os.path.join(BASE_DIR, "CLOSED_MSLP_RUN")
EXPERIMENTS = ["u-ch330", "u-ck777", "u-ck778", "u-dd828", "u-dc009"]
YEARS = range(1980, 2014)
SKIP_DOMAINS = set()

# None = no minimum 850 hPa vorticity (the point of the new version). If set,
# it is in raw s**-1 like the old run's 6e-5 - the tracker output is not in
# TRACK's 1e-5 units.
VORT_THRESHOLD = None
COHERENT = False  # new wcsi default; the old run effectively used True
# tcident-style warm-core test: 850 hPa minus upper-level vorticity must exceed
# this (raw s**-1, so 6e-5 = tcident's 6.0). None = off. wcsi computes it
# against 200 hPa, but these files have 250 hPa as their top level (tcident
# uses 850 - 250), so the WARM_CORE_UPPER_LEVEL level is passed to wcsi
# labelled as 200 for this test.
VORT_WARM_CORE_THRESHOLD = None
WARM_CORE_UPPER_LEVEL = 250
# Longitude of the MSLP minimum that TRACK added next to the MSLP field -
# 1e25 (NaN once loaded) where it found no minimum. In the model files the variable order is
# longitude_6, latitude_6, air_pressure_at_mean_sea_level, so it is _6
# (_7 belongs to wind_speed).
MSLP_LON_VARIABLE = "longitude_6"
TRACK_FILL_VALUE = 1e25
MAKE_MAPS = True


# ---------------------------------------------------------------------------
# Workaround for a bug in the new wcsi (reported to Leo Saffin; delete this
# block once fixed upstream): wcsi._nature.smooth_excursions does
# np.concat([...]) over the sequences to join, which raises "need at least one
# array to concatenate" when there are none - e.g. when no TC track has a short
# gap - and indexes starts[0]/ends[0] even when there are no sequences at all.
# Same function with those two cases guarded; otherwise unchanged.
# ---------------------------------------------------------------------------
def _smooth_excursions(track_id, criteria, *, min_count, starts=None, ends=None):
    if starts is None:
        starts = np.where(~criteria[:-1] & criteria[1:])[0] + 1
    if ends is None:
        ends = np.where(criteria[:-1] & ~criteria[1:])[0]

    if len(starts) == 0 and len(ends) == 0:
        if not criteria.any():
            return criteria.copy(), starts, ends
        # one sequence covering the whole array
        return criteria.copy(), np.array([0]), np.array([len(criteria) - 1])
    # first point is in a sequence: its start is missed
    if len(ends) and (len(starts) == 0 or starts[0] > ends[0]):
        starts = np.concat([[0], starts])
    # last point is in a sequence: its end is missed
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
# Identification
# ---------------------------------------------------------------------------
def find_input(experiment, year, hemisphere):
    """Same file selection as huracanpy_identify_batch_all.process_domain_year."""
    if hemisphere == "N":
        data_dir = os.path.join(BASE_DIR, experiment, f"VOR_VERTAVG_T63filt_{year}")
        file_glob = "tr_trs_pos*.nc"
    else:
        data_dir = os.path.join(BASE_DIR, experiment, f"VOR_VERTAVG_T63filt_{year}_SH")
        file_glob = "tr_trs_neg*.nc"
    matches = sorted(glob.glob(os.path.join(data_dir, file_glob)))
    if not matches:
        print(f"[{experiment} {year} {hemisphere}H] no {file_glob} found in {data_dir}, skipping")
        return None
    if len(matches) > 1:
        print(f"[{experiment} {year} {hemisphere}H] multiple matches in {data_dir}, using {matches[0]}")
    return matches[0]


def load_tracks(path):
    tracks = huracanpy.load(path, rename=dict(Hart_B="cps_b", Hart_Vtl="cps_vtl", Hart_Vtu="cps_vtu"))
    return gather_vorticity_profile(tracks)


def summarize(tracks, in_basin, is_tc):
    """Per-track summary, same columns as the old pipeline's summary.parquet,
    for every track with at least one point in the basin.
    """
    points = pd.DataFrame(dict(
        track_id=tracks.track_id.values, time=pd.to_datetime(tracks.time.values),
        lat=tracks.lat.values, lon=tracks.lon.values, in_basin=in_basin, is_tc=is_tc,
    ))
    g = points.groupby("track_id", sort=False)
    summary = pd.DataFrame(dict(
        storm_start=g.time.first(), storm_end=g.time.last(),
        origin_lat=g.lat.first(), origin_lon=g.lon.first(),
        end_lat=g.lat.last(), end_lon=g.lon.last(),
        in_basin=g.in_basin.any(), is_tc=g.is_tc.any(),
    ))
    return summary[summary.in_basin].drop(columns="in_basin").reset_index()


def process_file(experiment, year, path, domains):
    """Identify TCs for every domain of one year/hemisphere file."""
    tracks = load_tracks(path)
    track_id = tracks.track_id.values

    mslp_lon = tracks[MSLP_LON_VARIABLE].values
    no_mslp_min = ~np.isfinite(mslp_lon) | (mslp_lon >= TRACK_FILL_VALUE * 0.1)
    is_ocean = tracks.hrcn.get_is_ocean().values
    basin = tracks.hrcn.get_basin().values
    npoints = infer_npoints(tracks, TC_PERSISTENCE_HOURS)
    print(
        f"[{experiment} {year}] {os.path.basename(path)}: "
        f"{100 * no_mslp_min.mean():.1f}% of points without a valid MSLP minimum"
    )

    results = []
    for domain_name, domain in domains.items():
        try:
            results.append(process_domain(experiment, year, domain_name, domain, tracks, no_mslp_min, is_ocean, basin, npoints))
        except Exception as exc:
            print(f"[{experiment} {year} {domain_name}] FAILED: {exc}")
    return results


def process_domain(experiment, year, domain_name, domain, tracks, no_mslp_min, is_ocean, basin, npoints):
    """Identify TCs in one domain from one already-loaded file; returns its count row."""
    track_id = tracks.track_id.values
    tag = f"{experiment}_{year}_{domain_name}"
    in_basin = basin == domain["basin"]

    relative_vorticity = tracks.relative_vorticity
    if VORT_WARM_CORE_THRESHOLD is not None and WARM_CORE_UPPER_LEVEL != 200:
        # see VORT_WARM_CORE_THRESHOLD
        relative_vorticity = relative_vorticity.assign_coords(
            pressure=[200 if p == WARM_CORE_UPPER_LEVEL else p for p in relative_vorticity.pressure.values]
        )
    is_tc = wcsi.wcsi(
        track_id,
        tracks.cps_b.values,
        tracks.cps_vtl.values,
        tracks.cps_vtu.values,
        relative_vorticity,
        filter_all=in_basin & is_ocean & ~no_mslp_min,
        npoints=npoints,
        vort_threshold=VORT_THRESHOLD,
        vort_warm_core_threshold=VORT_WARM_CORE_THRESHOLD,
        coherent=COHERENT,
    )
    summary = summarize(tracks, in_basin, is_tc)
    summary.to_parquet(os.path.join(OUT_DIR, f"{tag}_summary.parquet"))
    n_tc = int(summary.is_tc.sum())

    # nature() classifies whole tracks, so run it on the TC tracks only -
    # a fresh subset per domain, keeping the full `tracks` for the next one
    tc_mask = np.isin(track_id, np.unique(track_id[is_tc]))
    if tc_mask.any():
        tc_tracks = tracks.isel(record=tc_mask)
        tc_tracks["is_tc"] = ("record", is_tc[tc_mask])
        tc_tracks["nature"] = ("record", wcsi.nature(
            tc_tracks.track_id.values,
            tc_tracks.cps_b.values,
            tc_tracks.cps_vtl.values,
            tc_tracks.cps_vtu.values,
            no_mslp_min[tc_mask],
            is_tc[tc_mask],
        ))
        tc_tracks = filter_genesis_latitude(tc_tracks, TC_GENESIS_THRESHOLD_LATITUDE)
        huracanpy.save(tc_tracks, os.path.join(OUT_DIR, f"{tag}_TCs.nc"))
        if MAKE_MAPS:
            plot_domain_nature(
                tc_tracks,
                title=f"{domain['title']} tracks by nature, {experiment} {year} ({n_tc} TCs)",
                out_path=os.path.join(OUT_DIR, f"{tag}_nature.png"),
                extent=domain["extent"],
            )

    print(f"[{experiment} {year} {domain_name}] {n_tc} tropical cyclones identified")
    return dict(experiment=experiment, year=year, domain=domain_name, n_tc=n_tc)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    results = []
    for experiment in EXPERIMENTS:
        for year in YEARS:
            for hemisphere in ("N", "S"):
                domains = {
                    name: domain for name, domain in DOMAINS.items()
                    if domain["hemisphere"] == hemisphere and name not in SKIP_DOMAINS
                }
                todo = {
                    name: domain for name, domain in domains.items()
                    if not os.path.exists(os.path.join(OUT_DIR, f"{experiment}_{year}_{name}_summary.parquet"))
                }
                if not todo:
                    print(f"[{experiment} {year} {hemisphere}H] already done, skipping")
                    continue
                path = find_input(experiment, year, hemisphere)
                if path is None:
                    continue
                try:
                    results.extend(process_file(experiment, year, path, todo))
                except Exception as exc:
                    print(f"[{experiment} {year} {hemisphere}H] FAILED: {exc}")

    if results:
        counts = pd.DataFrame(results)
        counts.to_csv(os.path.join(OUT_DIR, f"{'_'.join(EXPERIMENTS)}_TC_counts_by_year.csv"), index=False)
        counts.pivot_table(
            index=["experiment", "year"], columns="domain", values="n_tc", fill_value=0
        ).to_csv(os.path.join(OUT_DIR, f"{'_'.join(EXPERIMENTS)}_TC_counts_by_year_wide.csv"))


if __name__ == "__main__":
    main()
