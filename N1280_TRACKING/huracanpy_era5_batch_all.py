"""ERA5 reanalysis counterpart to huracanpy_identify_batch_all.py.

Produces the same style of per-domain, per-year "tracks by nature" maps as
huracanpy_identify_batch_all.py, but from ERA5-tracked cyclones instead of
model tracker output. Unlike the model script (and like the IBTrACS scripts),
there's no WCSI/CPS computation to run here: TC-TRACKS/ERA5_all_nature.nc was
already produced by wcsi's own `wcsi-nature` step (see wcsi/wcsi/scripts/
add_nature.py), so it already carries the exact same `nature` classification
(TC/BC/Tr/MV/Ot/Vo) and `is_tc` flag the model script computes for itself -
this just needs to load it, restrict to the same basins/domains, and plot it
the same way. Because of that, the plotting itself, DOMAINS, NATURE_MARKERS
etc. are all imported directly from huracanpy_identify_batch_all.py rather
than duplicated.

Unlike the model script's per-year-per-domain tracker files, ERA5_all_nature.nc
is one file covering the whole globe and the whole record (1940-2024) at once,
so basin assignment (huracanpy's own WMO-TC convention, via `.hrcn.add_basin()`
+ `.hrcn.get_apex_vals()`) and the per-track TC/basin selection that
wcsi.nature.wcsi() would otherwise do are done once for the full dataset up
front, rather than per year.

Outputs land in BASE_DIR (imported from huracanpy_identify_batch_all.py, i.e.
TC-TRACKS/) using the same {experiment}_{year}_{domain}_TCs.nc /
_summary.parquet naming as the model script (with experiment="ERA5"), so
huracanpy_TD.ipynb's density plots work on ERA5 too just by adding "ERA5" to
its EXPERIMENTS list - see huracanpy_era5_TD.py for a standalone equivalent.

Any other reanalysis already classified by wcsi-nature, with the same
variables as ERA5_all_nature.nc (e.g. JRA3Q), can be processed the same way
from the command line - outputs are then named {name}_{year}_{domain}_...:

    python huracanpy_era5_batch_all.py --name JRA3Q \
        --input TC-TRACKS/JRA3Q_nolat-nwc-tcident_WCSI_nature.nc \
        --out-dir TC-TRACKS/NEW_WCSI_RUN --no-maps
"""

import argparse

import os

import huracanpy
import numpy as np
import pandas as pd

from huracanpy_identify_batch_all import (
    BASE_DIR,
    DOMAINS,
    GAP_BASINS,
    TC_GENESIS_THRESHOLD_LATITUDE,
    filter_genesis_latitude,
    plot_domain_nature,
)

ERA5_FILE = os.path.join(BASE_DIR, "ERA5_all_nature.nc")
YEARS = range(1980, 2014)  # matches the model/IBTrACS scripts' default range;
                            # ERA5 itself covers 1940-2024 if you want to extend

SKIP_DOMAINS = set()  # ERA5 is dense in every domain/year, unlike the sparse
                       # model runs, so nothing is skipped by default


# ---------------------------------------------------------------------------
# ERA5 loading
# ---------------------------------------------------------------------------
def load_era5(path):
    """Load the pre-classified ERA5 tracks and tag every point with its basin."""
    tracks = huracanpy.load(path)
    tracks = tracks.rename(mslp="psl")  # matches get_mslp()'s expected names
    return tracks.hrcn.add_basin()  # WMO-TC convention: NATL/MED/WNP/SP/ENP/NI/SI/SA/...


def process_domain(tracks, apex, domain_name, domain):
    """All tracks whose peak-vorticity point is in this domain's basin, with a
    per-track summary and (for the ones that reach TC status while in that
    basin) the genesis-latitude-filtered set ready to plot - mirrors
    wcsi.nature.wcsi()'s own basin filter (see its `if basin is not None`
    block) plus huracanpy_identify_batch_all.py's filter_genesis_latitude,
    without re-deriving is_tc/nature, which ERA5_all_nature.nc already has.
    """
    basin = domain["basin"]
    candidate_ids = apex.track_id.values[apex.basin.values == basin]
    if len(candidate_ids) == 0:
        return None, None
    candidates = tracks.hrcn.sel_id(candidate_ids)

    rows = []
    tc_track_ids = []
    years = {}
    for track_id, track in candidates.groupby("track_id"):
        is_tc = bool((track.is_tc & (track.basin == basin)).values.any())
        t0 = pd.Timestamp(track.time.values[0])
        year = t0.year
        if domain["hemisphere"] == "S" and t0.month < 7:
            # matches the July-June season convention the model script's
            # separate "_SH" per-year files used
            year -= 1
        years[track_id] = year

        rows.append(dict(
            track_id=track_id, year=year,
            storm_start=t0, storm_end=pd.Timestamp(track.time.values[-1]),
            origin_lat=float(track.lat.values[0]), origin_lon=float(track.lon.values[0]),
            end_lat=float(track.lat.values[-1]), end_lon=float(track.lon.values[-1]),
            is_tc=is_tc,
        ))
        if is_tc:
            tc_track_ids.append(track_id)

    summary = pd.DataFrame(rows)

    if not tc_track_ids:
        return summary, None

    tc_tracks = candidates.hrcn.sel_id(tc_track_ids)
    tc_tracks = filter_genesis_latitude(tc_tracks, TC_GENESIS_THRESHOLD_LATITUDE)
    year_per_record = np.array([years[tid] for tid in tc_tracks.track_id.values])
    tc_tracks = tc_tracks.assign(year=("record", year_per_record))

    return summary, tc_tracks


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--name", default="ERA5", help="dataset name used in output file names")
    parser.add_argument("--input", default=ERA5_FILE, help="wcsi-nature classified track file")
    parser.add_argument("--out-dir", default=BASE_DIR, help="directory for the per-year outputs")
    parser.add_argument("--no-maps", action="store_true", help="skip the per-year tracks-by-nature maps")
    parser.add_argument(
        "--gap-basins", action="store_true",
        help="process only the AUS/CP basins no domain covers (GAP_BASINS), for merging later",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    name, out_dir = args.name, args.out_dir
    print(f"Loading {args.input} ...")
    tracks = load_era5(args.input)
    apex = tracks.hrcn.get_apex_vals("vorticity")

    results = []
    for domain_name, domain in (GAP_BASINS if args.gap_basins else DOMAINS).items():
        if domain_name in SKIP_DOMAINS:
            continue

        summary, tc_tracks = process_domain(tracks, apex, domain_name, domain)
        if summary is None:
            print(f"[{name} {domain_name}] no tracks in basin {domain['basin']!r}, skipping")
            continue

        for year in YEARS:
            tag = f"{name}_{year}_{domain_name}"
            year_summary = summary[summary.year == year].drop(columns="year")
            year_summary.to_parquet(os.path.join(out_dir, f"{tag}_summary.parquet"))
            n_tc = int(year_summary.is_tc.sum())

            year_tc_tracks = None
            if tc_tracks is not None:
                mask = tc_tracks.year.values == year
                if mask.any():
                    year_tc_tracks = tc_tracks.isel(record=mask)

            if year_tc_tracks is not None:
                huracanpy.save(year_tc_tracks, os.path.join(out_dir, f"{tag}_TCs.nc"))
                if not args.no_maps:
                    plot_domain_nature(
                        year_tc_tracks,
                        title=f"{domain['title']} tracks by nature, {name} {year} ({n_tc} TCs)",
                        out_path=os.path.join(out_dir, f"{tag}_nature.png"),
                        extent=domain["extent"],
                    )
            else:
                print(f"[{name} {year} {domain_name}] no TC tracks, skipping plot")

            print(f"[{name} {year} {domain_name}] {n_tc} tropical cyclones identified")
            results.append(dict(experiment=name, year=year, domain=domain_name, n_tc=n_tc))

    if results:
        counts = pd.DataFrame(results)
        suffix = "_gap_basins" if args.gap_basins else ""
        counts.to_csv(os.path.join(out_dir, f"{name}_TC_counts_by_year{suffix}.csv"), index=False)

        wide = counts.pivot_table(
            index=["experiment", "year"], columns="domain", values="n_tc", fill_value=0
        )
        wide.to_csv(os.path.join(out_dir, f"{name}_TC_counts_by_year{suffix}_wide.csv"))


if __name__ == "__main__":
    main()
