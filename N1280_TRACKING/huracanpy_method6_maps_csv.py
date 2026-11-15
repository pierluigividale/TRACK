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

It reads the NetCDF files only - nothing is recomputed - so it is quick to
re-run, and existing PNG/CSV files are skipped unless --overwrite.

Usage:
  python huracanpy_method6_maps_csv.py TC-TRACKS/DESTINE_METHOD6_RUN
  python huracanpy_method6_maps_csv.py OUT_DIR --experiments hist --no-maps
"""

import argparse
import glob
import os
import re

import huracanpy
import numpy as np
import pandas as pd

import huracanpy_identify_method6 as m6

TC_FILE_RE = re.compile(r"^(?P<exp>.+)_(?P<year>\d{4})_(?P<domain>[A-Z]+)_TCs\.nc$")
FILL = 1e25


def tracks_table(tracks):
    """Per-point DataFrame with readable column names."""
    df = pd.DataFrame({"track_id": tracks.track_id.values,
                       "time": pd.to_datetime(tracks.time.values),
                       "lon": tracks.lon.values, "lat": tracks.lat.values,
                       "vorticity_T63": tracks.vorticity.values})
    if "relative_vorticity" in tracks and "pressure" in tracks.relative_vorticity.dims:
        for i, p in enumerate(tracks.pressure.values):
            df[f"relative_vorticity_{int(p)}"] = tracks.relative_vorticity.values[:, i]
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
    df["Hart_B"] = tracks.cps_b.values
    df["Hart_VTL"] = tracks.cps_vtl.values
    df["Hart_VTU"] = tracks.cps_vtu.values
    df["is_tc"] = tracks.is_tc.values.astype(bool)
    df["nature"] = tracks.nature.values.astype(str)
    num = df.select_dtypes("number").columns
    df[num] = df[num].mask(df[num].abs() >= FILL * 0.1)  # TRACK's 1e25 "missing" -> empty
    return df


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("out_dir")
    p.add_argument("--experiments", nargs="+")
    p.add_argument("--no-maps", action="store_true")
    p.add_argument("--no-csv", action="store_true")
    p.add_argument("--overwrite", action="store_true")
    args = p.parse_args()

    files = sorted(glob.glob(os.path.join(args.out_dir, "*_TCs.nc")))
    n_png = n_csv = 0
    for f in files:
        mt = TC_FILE_RE.match(os.path.basename(f))
        if not mt or (args.experiments and mt["exp"] not in args.experiments):
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
