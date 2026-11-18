"""Maps and CSV tables for an existing method-6 output folder.

For every {experiment}_{year}_{domain}_TCs.nc written by
huracanpy_identify_method6.py, this script writes:
  {...}_nature.png  map of the TC tracks' whole lifetimes (marker = WCSI nature,
                    colour = MSLP), the same figure as `--maps`
  {...}_tracks.csv  one row per track point: track_id, time, lon, lat, the
                    TRACK T63 vorticity, relative vorticity at every level,
                    MSLP and 10 m wind speed (with the positions TRACK found
                    them at), Hart B / VTL / VTU, is_tc (WCSI TC criteria met
                    at that point) and nature (TC/BC/Tr/MV/Ot/Vo)

With --combine, it instead joins all years of each experiment and domain
(optionally limited by --years) and writes one set of files for the period:
  {experiment}_{y0}-{y1}_{domain}_allTCs.nc   the joined tracks (huracanpy NetCDF)
  {experiment}_{y0}-{y1}_{domain}_nature.png  all tracks of the period on one map
  {experiment}_{y0}-{y1}_{domain}_tracks.csv  as above, plus a "year" column
Track ids are renumbered year*1000000 + original id, so they stay unique
across years (the original TRACK id is the last 6 digits). y0-y1 are the
first and last years found; missing years in between are listed in the log.

It reads the NetCDF files only - nothing is recomputed - so it is quick to
re-run, and existing PNG/CSV files are skipped unless --overwrite.

Usage:
  python huracanpy_method6_maps_csv.py TC-TRACKS/DESTINE_METHOD6_RUN
  python huracanpy_method6_maps_csv.py OUT_DIR --experiments hist --no-maps
  python huracanpy_method6_maps_csv.py OUT_DIR --combine
  python huracanpy_method6_maps_csv.py OUT_DIR --combine --years 1990-1999 --experiments hist
"""

import argparse
import glob
import os
import re

import huracanpy
import numpy as np
import pandas as pd
import xarray as xr

import huracanpy_identify_method6 as m6

TC_FILE_RE = re.compile(r"^(?P<exp>.+)_(?P<year>\d{4})_(?P<domain>[A-Z]+)_TCs\.nc$")
FILL = 1e25


def tracks_table(tracks):
    """Per-point DataFrame with readable column names.

    Reads both _TCs.nc layouts: the current one (TRACK's own names - Hart_B /
    Hart_Vtl / Hart_Vtu, tracked T63 vorticity "relative_vorticity",
    relative_vorticity_{level}) and the one written before 2026-11-18 (cps_b /
    cps_vtl / cps_vtu, T63 vorticity "vorticity", profile merged into
    relative_vorticity(record, pressure)). The CSV columns are the same for both.
    """
    merged = "relative_vorticity" in tracks and "pressure" in tracks.relative_vorticity.dims
    t63 = tracks.vorticity if merged else tracks.relative_vorticity
    df = pd.DataFrame({"track_id": tracks.track_id.values,
                       "time": pd.to_datetime(tracks.time.values),
                       "lon": tracks.lon.values, "lat": tracks.lat.values,
                       "vorticity_T63": t63.values})
    if merged:
        for i, p in enumerate(tracks.pressure.values):
            df[f"relative_vorticity_{int(p)}"] = tracks.relative_vorticity.values[:, i]
    else:
        for v in tracks.data_vars:
            if re.fullmatch(r"relative_vorticity_\d+", v):
                df[v] = tracks[v].values
    # fields TRACK added with their own positions (longitude_N / latitude_N)
    names = list(tracks.variables)
    for field in ("air_pressure_at_mean_sea_level", "psl", "wind_speed", "wind_speed_10m"):
        if field not in tracks.variables:
            continue
        df[field] = tracks[field].values
        i = names.index(field)
        if i >= 2 and re.fullmatch(r"longitude_\d+", names[i - 2]):
            df[f"{field}_lon"] = tracks[names[i - 2]].values
            df[f"{field}_lat"] = tracks[names[i - 1]].values
    old = "cps_b" in tracks
    df["Hart_B"] = (tracks.cps_b if old else tracks.Hart_B).values
    df["Hart_VTL"] = (tracks.cps_vtl if old else tracks.Hart_Vtl).values
    df["Hart_VTU"] = (tracks.cps_vtu if old else tracks.Hart_Vtu).values
    df["is_tc"] = tracks.is_tc.values.astype(bool)
    df["nature"] = tracks.nature.values.astype(str)
    num = df.select_dtypes("number").columns
    df[num] = df[num].mask(df[num].abs() >= FILL * 0.1)  # TRACK's 1e25 "missing" -> empty
    return df


def combine_years(files_by_year):
    """Join {year: path} into one record-dimension dataset with unique track ids."""
    parts = []
    for year, f in sorted(files_by_year.items()):
        t = huracanpy.load(f)
        t = t.drop_vars([v for v in t.variables if "record" not in t[v].dims and v not in t.dims
                         and v != "pressure"])
        t["track_id"] = t.track_id.astype("int64") + year * 1_000_000
        t["year"] = ("record", np.full(t.sizes["record"], year, dtype="int32"))
        parts.append(t)
    return xr.concat(parts, dim="record", data_vars="minimal", coords="minimal", compat="override")


def run_combined(args):
    groups = {}  # (exp, domain) -> {year: path}
    for f in sorted(glob.glob(os.path.join(args.out_dir, "*_TCs.nc"))):
        mt = TC_FILE_RE.match(os.path.basename(f))
        if not mt or (args.experiments and mt["exp"] not in args.experiments):
            continue
        year = int(mt["year"])
        if args.years and year not in args.years:
            continue
        groups.setdefault((mt["exp"], mt["domain"]), {})[year] = f

    n_png = n_csv = n_nc = 0
    for (exp, dom), by_year in sorted(groups.items()):
        years = sorted(by_year)
        y0, y1 = years[0], years[-1]
        stem = os.path.join(args.out_dir, f"{exp}_{y0}-{y1}_{dom}")
        nc, png, csv = stem + "_allTCs.nc", stem + "_nature.png", stem + "_tracks.csv"
        want_nc = args.overwrite or not os.path.exists(nc)
        want_png = not args.no_maps and (args.overwrite or not os.path.exists(png))
        want_csv = not args.no_csv and (args.overwrite or not os.path.exists(csv))
        if not (want_nc or want_png or want_csv):
            continue
        tracks = combine_years(by_year)
        n_tc = len(np.unique(tracks.track_id.values))
        if want_nc:
            huracanpy.save(tracks, nc)
            n_nc += 1
        if want_csv:
            df = tracks_table(tracks)
            df.insert(1, "year", tracks.year.values)
            df.to_csv(csv, index=False, float_format="%.6g")
            n_csv += 1
        if want_png:
            domain = m6.DOMAINS[dom]
            m6.plot_domain_nature(
                tracks, f"{domain['title']} tracks by nature, {exp} {y0}-{y1} "
                        f"({n_tc} TCs in {len(years)} years)", png, domain["extent"])
            n_png += 1
        missing = sorted(set(range(y0, y1 + 1)) - set(years))
        if args.years:
            missing = [y for y in missing if y in args.years]
        # a year with a summary but no _TCs.nc was processed and had zero TCs
        zero = [y for y in missing if os.path.exists(os.path.join(args.out_dir, f"{exp}_{y}_{dom}_summary.parquet"))]
        absent = [y for y in missing if y not in zero]
        print(f"{os.path.basename(stem)}: {n_tc} TCs from {len(years)} years"
              + (f" (zero TCs in {zero})" if zero else "") + (f" (NOT PROCESSED: {absent})" if absent else "")
              + (" nc" if want_nc else "") + (" csv" if want_csv else "") + (" png" if want_png else ""),
              flush=True)
    print(f"wrote {n_nc} NetCDF, {n_csv} CSV and {n_png} PNG files in {args.out_dir}")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("out_dir")
    p.add_argument("--experiments", nargs="+")
    p.add_argument("--years", type=m6.parse_years, help="e.g. 1990-2014 (inclusive) or 1990,2000")
    p.add_argument("--combine", action="store_true",
                   help="one map/CSV/NetCDF per experiment and domain for all (selected) years")
    p.add_argument("--no-maps", action="store_true")
    p.add_argument("--no-csv", action="store_true")
    p.add_argument("--overwrite", action="store_true")
    args = p.parse_args()
    if args.combine:
        return run_combined(args)

    files = sorted(glob.glob(os.path.join(args.out_dir, "*_TCs.nc")))
    n_png = n_csv = 0
    for f in files:
        mt = TC_FILE_RE.match(os.path.basename(f))
        if not mt or (args.experiments and mt["exp"] not in args.experiments):
            continue
        if args.years and int(mt["year"]) not in args.years:
            continue
        stem = f[:-len("_TCs.nc")]
        png, csv = stem + "_nature.png", stem + "_tracks.csv"
        want_png = not args.no_maps and (args.overwrite or not os.path.exists(png))
        want_csv = not args.no_csv and (args.overwrite or not os.path.exists(csv))
        if not (want_png or want_csv):
            continue
        tracks = huracanpy.load(f)
        n_tc = len(np.unique(tracks.track_id.values))
        if want_csv:
            tracks_table(tracks).to_csv(csv, index=False, float_format="%.6g")
            n_csv += 1
        if want_png:
            domain = m6.DOMAINS[mt["domain"]]
            m6.plot_domain_nature(
                tracks, f"{domain['title']} tracks by nature, {mt['exp']} {mt['year']} ({n_tc} TCs)",
                png, domain["extent"])
            n_png += 1
        print(f"{os.path.basename(stem)}: {n_tc} TCs" + (" csv" if want_csv else "") + (" png" if want_png else ""),
              flush=True)
    print(f"wrote {n_csv} CSV and {n_png} PNG files in {args.out_dir}")


if __name__ == "__main__":
    main()
