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

With --combine, it instead joins tracks into one set of files, along time,
along the ensemble, or both:
  --combine years    per experiment and domain, all (selected) years:
                     {experiment}_{y0}-{y1}_{domain}_{allTCs.nc,nature.png,tracks.csv}
  --combine members  per year and domain, all members:
                     {ENSEMBLE}_{year}_{domain}_...
  --combine both     per domain, all members and all years:
                     {ENSEMBLE}_{y0}-{y1}_{domain}_...
(a bare --combine means "years", as before). allTCs.nc is the joined tracks
(huracanpy NetCDF), nature.png all of them on one map, tracks.csv as above
plus "year" (and "member") columns.

The members are --experiments, or the present-day ensemble
(huracanpy_identify_method6.EXPERIMENTS) when it is not given; ENSEMBLE is
--ensemble-name (default "ENSEMBLE"). The ensemble is not an xarray dimension
(tracks have different lengths and counts, so that would mean NaN padding and
huracanpy could not use it) but a "member" variable on the record dimension,
next to "year": tracks.groupby("member"), tracks.where(tracks.member == "u-dc009",
drop=True), etc.

Track ids are renumbered so they stay unique across what was joined:
  years:            year*1e6 + original id
  members / both:   member_index*1e10 + year*1e6 + original id
(member_index = position in the member list, from 0); the original TRACK id
is also kept as "track_id_orig". Members and both only use years that every
member has processed (a year with a summary but no _TCs.nc had zero TCs and
counts as processed), so a missing year of one member does not look like
fewer storms; dropped years are listed in the log. y0-y1 are the first and
last years used; missing years in between are listed in the log.

It reads the NetCDF files only - nothing is recomputed - so it is quick to
re-run, and existing PNG/CSV files are skipped unless --overwrite.

Usage:
  python huracanpy_method6_maps_csv.py TC-TRACKS/DESTINE_METHOD6_RUN
  python huracanpy_method6_maps_csv.py OUT_DIR --experiments hist --no-maps
  python huracanpy_method6_maps_csv.py OUT_DIR --combine
  python huracanpy_method6_maps_csv.py OUT_DIR --combine --years 1990-1999 --experiments hist
  python huracanpy_method6_maps_csv.py OUT_DIR --combine members --years 2005
  python huracanpy_method6_maps_csv.py OUT_DIR --combine both --experiments u-ch330 u-ck777 u-ck778 u-dc009 u-dd828 \\
      --ensemble-name ENS5 --output-dir OUT_DIR/combined
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
SUMMARY_OR_TC_RE = re.compile(r"^(?P<exp>.+)_(?P<year>\d{4})_(?P<domain>[A-Z]+)_(?P<kind>TCs\.nc|summary\.parquet)$")
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


def combine(parts):
    """Join [(member, member_index, year, path)] into one record-dimension
    dataset with unique track ids, plus per-point "year" (and "member" when
    member_index is not None) and "track_id_orig"."""
    out = []
    for member, idx, year, f in parts:
        t = huracanpy.load(f)
        t = t.drop_vars([v for v in t.variables if "record" not in t[v].dims and v not in t.dims
                         and v != "pressure"])
        orig = t.track_id.astype("int64")
        t["track_id_orig"] = orig
        t["track_id"] = orig + year * 1_000_000 + (idx * 10_000_000_000 if idx is not None else 0)
        t["year"] = ("record", np.full(t.sizes["record"], year, dtype="int32"))
        if idx is not None:
            t["member"] = ("record", np.full(t.sizes["record"], member, dtype=object))
        out.append(t)
    ds = xr.concat(out, dim="record", data_vars="minimal", coords="minimal", compat="override")
    if "member" in ds:
        ds["member"] = ds.member.astype(str)
    return ds


def year_ranges(years):
    """[1990, 1991, 1992, 1995] -> "1990-1992, 1995"."""
    out, years = [], sorted(years)
    for y in years:
        if out and y == out[-1][1] + 1:
            out[-1][1] = y
        else:
            out.append([y, y])
    return ", ".join(str(a) if a == b else f"{a}-{b}" for a, b in out)


def processed(out_dir, exp, year, dom):
    """True if huracanpy_identify_method6.py processed this year (a summary
    without a _TCs.nc means zero TCs)."""
    stem = os.path.join(out_dir, f"{exp}_{year}_{dom}")
    return os.path.exists(stem + "_TCs.nc") or os.path.exists(stem + "_summary.parquet")


def scan(out_dir, experiments, years):
    """{(exp, domain): {year: path}} for the _TCs.nc files, and
    {domain: years any of the experiments processed} (from summaries too)."""
    groups, domains = {}, {}
    for f in sorted(glob.glob(os.path.join(out_dir, "*_TCs.nc")) + glob.glob(os.path.join(out_dir, "*_summary.parquet"))):
        mt = SUMMARY_OR_TC_RE.match(os.path.basename(f))
        if not mt or (experiments and mt["exp"] not in experiments):
            continue
        year = int(mt["year"])
        if years and year not in years:
            continue
        domains.setdefault(mt["domain"], set()).add(year)
        if mt["kind"] == "TCs.nc":
            groups.setdefault((mt["exp"], mt["domain"]), {})[year] = f
    return groups, domains


def write_outputs(tracks, stem, title_fn, dom, args, want):
    """Write the NetCDF / CSV / PNG for one combined dataset; returns what was written."""
    want_nc, want_csv, want_png = want
    if want_nc:
        huracanpy.save(tracks, stem + "_allTCs.nc")
    if want_csv:
        df = tracks_table(tracks)
        df.insert(1, "year", tracks.year.values)
        if "member" in tracks:
            df.insert(1, "member", tracks.member.values)
        df.insert(len(df.columns), "track_id_orig", tracks.track_id_orig.values)
        df.to_csv(stem + "_tracks.csv", index=False, float_format="%.6g")
    if want_png:
        domain = m6.DOMAINS[dom]
        m6.plot_domain_nature(tracks, title_fn(domain["title"]), stem + "_nature.png", domain["extent"])
    return [k for k, w in zip(("nc", "csv", "png"), want) if w]


def wanted(stem, args):
    nc, png, csv = stem + "_allTCs.nc", stem + "_nature.png", stem + "_tracks.csv"
    return (args.overwrite or not os.path.exists(nc),
            not args.no_csv and (args.overwrite or not os.path.exists(csv)),
            not args.no_maps and (args.overwrite or not os.path.exists(png)))


def n_tracks(tracks):
    return len(np.unique(tracks.track_id.values))


def run_combined(args):
    dest = args.output_dir or args.out_dir
    os.makedirs(dest, exist_ok=True)
    counts = {"nc": 0, "csv": 0, "png": 0}
    if args.combine == "years":
        groups, _ = scan(args.out_dir, args.experiments, args.years)
        for (exp, dom), by_year in sorted(groups.items()):
            years = sorted(by_year)
            y0, y1 = years[0], years[-1]
            stem = os.path.join(dest, f"{exp}_{y0}-{y1}_{dom}")
            want = wanted(stem, args)
            if not any(want):
                continue
            tracks = combine([(exp, None, y, by_year[y]) for y in years])
            n_tc = n_tracks(tracks)
            done = write_outputs(tracks, stem, lambda t: f"{t} tracks by nature, {exp} {y0}-{y1} "
                                 f"({n_tc} TCs in {len(years)} years)", dom, args, want)
            for k in done:
                counts[k] += 1
            missing = sorted(set(range(y0, y1 + 1)) - set(years))
            if args.years:
                missing = [y for y in missing if y in args.years]
            zero = [y for y in missing if processed(args.out_dir, exp, y, dom)]
            absent = [y for y in missing if y not in zero]
            print(f"{os.path.basename(stem)}: {n_tc} TCs from {len(years)} years"
                  + (f" (zero TCs in {zero})" if zero else "") + (f" (NOT PROCESSED: {absent})" if absent else "")
                  + "".join(" " + k for k in done), flush=True)
    else:
        members = list(dict.fromkeys(args.experiments or m6.EXPERIMENTS))
        ens = args.ensemble_name
        groups, domains = scan(args.out_dir, members, args.years)
        print(f"{ens} = {', '.join(members)}")
        for dom in sorted(domains):
            seen = sorted(set(args.years) if args.years else domains[dom])
            common = [y for y in seen if all(processed(args.out_dir, m, y, dom) for m in members)]
            lacking = {m: [y for y in seen if not processed(args.out_dir, m, y, dom)] for m in members}
            if any(lacking.values()):
                print(f"[{dom}] years left out, not processed by: "
                      + "; ".join(f"{m} {year_ranges(ys)}" for m, ys in lacking.items() if ys), flush=True)
            if not common:
                continue
            batches = ([(y, [y]) for y in common] if args.combine == "members"
                       else [(f"{common[0]}-{common[-1]}", common)])
            for label, years in batches:
                stem = os.path.join(dest, f"{ens}_{label}_{dom}")
                want = wanted(stem, args)
                if not any(want):
                    continue
                parts = [(m, i, y, groups[(m, dom)][y]) for i, m in enumerate(members) for y in years
                         if y in groups.get((m, dom), {})]
                if not parts:
                    print(f"{os.path.basename(stem)}: zero TCs in every member, nothing written", flush=True)
                    continue
                tracks = combine(parts)
                n_tc = n_tracks(tracks)
                per = {m: len(np.unique(tracks.track_id.values[tracks.member.values == m])) for m in members}
                rate = n_tc / (len(members) * len(years))
                done = write_outputs(
                    tracks, stem,
                    lambda t: f"{t} tracks by nature, {ens} ({len(members)} members) {label}: {n_tc} TCs, "
                              + (f"{rate:.1f} per member" if len(years) == 1 else
                                 f"{len(members)} members x {len(years)} years = {rate:.1f} per member-year"),
                    dom, args, want)
                for k in done:
                    counts[k] += 1
                gaps = sorted(set(range(years[0], years[-1] + 1)) - set(years))
                print(f"{os.path.basename(stem)}: {n_tc} TCs ({', '.join(f'{m} {n}' for m, n in per.items())}),"
                      f" {rate:.1f} per member-year" + (f" (years missing: {gaps})" if gaps and len(years) > 1 else "")
                      + "".join(" " + k for k in done), flush=True)
    print(f"wrote {counts['nc']} NetCDF, {counts['csv']} CSV and {counts['png']} PNG files in {dest}")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("out_dir")
    p.add_argument("--experiments", nargs="+")
    p.add_argument("--years", type=m6.parse_years, help="e.g. 1990-2014 (inclusive) or 1990,2000")
    p.add_argument("--combine", nargs="?", const="years", choices=["years", "members", "both"],
                   help="join tracks along time (years; the default for a bare --combine), "
                        "the ensemble (members) or both")
    p.add_argument("--ensemble-name", default="ENSEMBLE",
                   help="file-name label for --combine members/both (default ENSEMBLE)")
    p.add_argument("--output-dir", help="where --combine writes (default: out_dir)")
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
