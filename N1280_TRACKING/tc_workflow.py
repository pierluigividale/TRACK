"""Audit and drive the TRACK-1.5.2 WCSI pipeline across multiple experiments.

Per (experiment, year, hemisphere) the pipeline is:
  convert   pp_to_nc_iris_WCSI.py + ll_to_gg_n1280.py, one SLURM job per month
            (SLURMscript_n1280_para_WCSI), building annual _geo/_mslp/
            _10m_wspeed/_uwinds_gg/_vwinds_gg .nc files in scrdir.
  track     track_n1280_para_WCSI.py (SLURMscript_n1280_track_para_WCSI):
            filtered vorticity, TRACK, Hart phase space, tcident, final
            NetCDF - submitted with --dependency=singleton so it waits for
            that year's convert jobs (same --job-name) to finish first.

`status` only reads on-disk markers (and squeue) - no side effects. `submit`
decides the minimal fix per incomplete (experiment, year, hemisphere) and
only acts with --yes; otherwise it just prints the sbatch calls it would run.

Usage:
  python tc_workflow.py status [--experiments ch330,ck777] [--years 1984-1990] [--hemisphere NH]
  python tc_workflow.py submit [--experiments ck777] [--years 1987] [--yes]
"""

import argparse
import glob
import os
import re
import subprocess

HOME = "/home/users/plvidale/TRACK-1.5.2/"
SCRIPTS_DIR = HOME + "scripts/"
TRUNC = "T63"

# Per-experiment source group and year range per hemisphere - verified
# directly on JASMIN (see plan). A hemisphere missing from `years` is not
# processed for that experiment at all. scrdir/permdir follow the existing
# submit scripts' convention and need no per-experiment override.
#
# NH range: capped at 2013, NOT 2014. All 5 experiments' raw archives have
# monthly files through Dec 2014, but the Dec 2014 file itself is truncated -
# the model runs actually stop on 2014-12-21 (confirmed identically for
# ck777 and ck778: annual geo.nc time axis ends 2014-12-21 00:00, and the
# raw apc archive is missing the normal third ~10-day December chunk). With
# a genuinely short final year, master_tracking()'s hardcoded -n=1,62,24
# chunking runs out of data mid-chunk and silently exits with no track
# output at all - every NH 2014 job fails identically for every experiment,
# wasting the full run before failing at the master step. Not fixable
# without either 2015 data (doesn't exist) or Kevin adjusting master's
# chunking to track a short final period - ask him before re-adding 2014.
# SH range: an SH "year" Y covers Jul(Y)-Jun(Y+1), so SH year=2014 would need
# 2015 data that doesn't exist - SH stays capped at 2013. ch330/ck777/ck778's
# raw data starts Mar 1979 (SH 1979 = Jul1979-Jun1980, fully present), but
# dd828/dc009's starts Dec 1979/Jan 1980 respectively - SH 1979 needs
# Jul-Nov 1979, which for those two is only on tape (slow to retrieve), so
# their SH range is left at 1980 for now.
#
# dz876/ec100 (coupled atmosphere-ocean, EERIE project, still running as of
# 2026-10-22) use a completely different raw-data convention - see
# pp_to_nc_iris_WCSI_EERIE.py and its dispatch in SLURMscript_n1280_para_WCSI
# (picked by run name, see EERIE_EXPERIMENTS below). Atmosphere grid is the
# same N1280 (1920x2560, checked directly), so TRACK's own config needs no
# change - only the conversion step differs.
#
# Their apparent raw-archive date ranges are NOT what's actually usable -
# most of dz876's directory listing (1950-1958) is an empty placeholder
# (owned by a different user, "lawrence", not the run's own "aosprey") left
# over from something else; the real run only starts Dec 1958. Verified by
# checking every month's file count directly, not just the directory
# listing, which is also how a real one-month gap (Dec 1959, both pc and p8
# streams) was caught - don't trust the raw directory listing for this pair
# without re-checking file-level completeness, same lesson as dd828/dc009.
#   dz876: real data Dec 1958 - Jul 1963, except Dec 1959 (missing) ->
#          full NH years 1960-1962, full SH seasons 1960-1962 (1959 excluded
#          by the Dec 1959 gap on both sides of the Jul-Jun boundary)
#   ec100: real data Dec 1949 - Jan 1951, no gaps in that span -> only NH
#          1950 is a full year; no SH season is complete yet
# Both are still running - re-check for newly-completed years before
# re-running `submit` next time, rather than assuming this range is final.
EXPERIMENTS = {
    "ch330": dict(source="n1280run", years={"NH": range(1980, 2014), "SH": range(1980, 2014)}),
    "ck777": dict(source="n1280run", years={"NH": range(1980, 2014), "SH": range(1980, 2014)}),
    "ck778": dict(source="n1280run", years={"NH": range(1980, 2014), "SH": range(1980, 2014)}),
    "dd828": dict(source="plvidale", years={"NH": range(1980, 2014), "SH": range(1980, 2014)}),
    "dc009": dict(source="torau", years={"NH": range(1980, 2014), "SH": range(1980, 2014)}),
    "dz876": dict(source="eerie", years={"NH": range(1960, 1963), "SH": range(1960, 1963)}),
    "ec100": dict(source="eerie", years={"NH": range(1950, 1951)}),
}

# Experiments whose conversion step goes through pp_to_nc_iris_WCSI_EERIE.py
# instead of pp_to_nc_iris_WCSI.py - SLURMscript_n1280_para_WCSI dispatches
# on this same name list, keep the two in sync.
EERIE_EXPERIMENTS = {run for run, spec in EXPERIMENTS.items() if spec["source"] == "eerie"}

# Raw .pp location per source group - read by pp_to_nc_iris_WCSI.py (kept
# here too so `status` can report on raw-data availability if ever needed).
SOURCE_DIRS = {
    "n1280run": lambda run, year: dict(
        apc="/gws/ssde/j25a/hrcm/n1280run/u-%s/%s/" % (run, "6hrly" if int(year) < 2000 else "apc"),
        ap7=None,  # levels are in the same apc/6hrly files for this group
    ),
    "plvidale": lambda run, year: dict(
        apc="/gws/ssde/j25a/hrcm/plvidale/u-%s/apc/" % run,
        ap7="/gws/ssde/j25a/hrcm/plvidale/u-%s/ap7/" % run,
    ),
    "torau": lambda run, year: dict(
        apc="/gws/ssde/j25a/hrcm/torau/u-%s/apc/" % run,
        ap7="/gws/ssde/j25a/hrcm/torau/u-%s/ap7/" % run,
    ),
    # Not a flat directory like the others - per-month subdirectories under
    # this base, named by a cycle_date offset one month ahead of the data
    # (see pp_to_nc_iris_WCSI_EERIE.py). Kept here as the base path only.
    "eerie": lambda run, year: dict(
        apc="/gws/ssde/j25a/hiresgw/aosprey/u-%s/" % run,
        ap7="/gws/ssde/j25a/hiresgw/aosprey/u-%s/" % run,
    ),
}

HEMI_SUFFIX = {"NH": "", "SH": "_SH"}  # matches track_n1280_para_WCSI.py's hemi variable
HEMI_ARG = {"NH": "_NH", "SH": "_SH"}  # matches the existing submit scripts' positional arg
TRACKNAME = {"NH": "tr_trs_pos", "SH": "tr_trs_neg"}


def scrdir(run):
    return "/work/scratch-pw4/plvidale/u-%s/" % run


def permdir(run):
    return "/gws/ssde/j25a/hrcm/N1280_analysis/u-%s/" % run


# ---------------------------------------------------------------------------
# Status
# ---------------------------------------------------------------------------
def convert_done(run, year, hemi):
    """Annual geo/mslp/10m_wspeed files exist, in scrdir (not yet cleaned up)
    or permdir (cleaned up after a fully successful year) - either means the
    convert stage finished for this (run, year, hemisphere).
    """
    suffix = HEMI_SUFFIX[hemi]
    for var in ("_geo", "_mslp", "_10m_wspeed"):
        name = "%s_%s%s%s.nc" % (run, year, var, suffix)
        if not (os.path.exists(scrdir(run) + name) or os.path.exists(permdir(run) + name)):
            return False
    return True


def vort_ready(run, year, hemi):
    """The T63-filtered, vertically-averaged vorticity TRACK needs - the
    input avg_filter_lev_vort() in track_n1280_para_WCSI.py produces.
    """
    suffix = HEMI_SUFFIX[hemi]
    name = "vor_%s_%s_avg_%s_filt%s.dat" % (run, year, TRUNC, suffix)
    return os.path.exists(permdir(run) + name)


def fully_done(run, year, hemi):
    """The final product: what huracanpy_identify_batch_all.py's tr_trs_pos*/
    tr_trs_neg*.nc glob actually reads.
    """
    trackdir = "%su-%s/VOR_VERTAVG_%sfilt_%s%s/" % (HOME, run, TRUNC, year, HEMI_SUFFIX[hemi])
    name = trackdir + TRACKNAME[hemi] + ".2day_addvor_addmslp_addwind.hart.new.nc"
    return os.path.exists(name)


# Job names: {run}_convert_{year} / {run}_track_{year} - readable in squeue
# (which experiment/stage/year at a glance) and, since run is now part of
# the name, collision-free between experiments processing the same year
# concurrently (the old n1280_TRACK_<year> scheme didn't encode run at all).
JOB_NAME_RE = re.compile(r"^(?P<run>[a-z0-9]+)_(?:convert|track)_(?P<year>\d+)$")


def running_jobs():
    """((run, year, hemi) precise matches, {(run, year)} imprecise fallback)
    for jobs currently in squeue.

    Precise: SLURMscript_n1280_track_para_WCSI (submitted directly, e.g. by
    submit_SLURMscript_n1280_track_para_WCSI, with no --job-name) echoes its
    args as the first line of {jobid}.out - "run year trunc scrdir permdir
    link_var_log add_names only_tcident hemi" - read that to attribute the
    job exactly, so submit() doesn't double-submit on top of it.

    Fallback: jobs matching JOB_NAME_RE that haven't started yet (no .out
    file: still waiting on their dependency) - hemisphere isn't in the name,
    so this is "don't touch this (run, year) in any hemisphere" rather than
    a fully precise match.
    """
    try:
        out = subprocess.run(
            ["squeue", "-u", "plvidale", "-h", "-o", "%i %j"],
            capture_output=True, text=True, check=True,
        ).stdout
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        print("warning: squeue check failed (%s); treating nothing as running" % exc)
        return set(), set()

    precise, year_fallback = set(), set()
    for line in out.splitlines():
        parts = line.split(None, 1)
        if not parts:
            continue
        jobid, name = parts[0], (parts[1] if len(parts) > 1 else "")
        out_path = os.path.join(SCRIPTS_DIR, "%s.out" % jobid)
        run = year = None
        hemi = "NH"
        if os.path.exists(out_path):
            try:
                with open(out_path) as f:
                    head = f.read(2000)
            except OSError:
                head = ""
            for tok in head.split():
                if tok in EXPERIMENTS and run is None:
                    run = tok
                elif tok.isdigit() and len(tok) == 4 and year is None:
                    year = tok
                elif tok == "_SH":
                    hemi = "SH"
        if run and year:
            precise.add((run, year, hemi))
        else:
            m = JOB_NAME_RE.match(name)
            if m:
                year_fallback.add((m.group("run"), m.group("year")))
    return precise, year_fallback


def audit(experiments, years_filter, hemispheres):
    """[(run, year, hemi, convert, vort, done, state)], state one of
    "ok" / "running" / "needs_convert" / "needs_track".
    """
    running_precise, running_year_fallback = running_jobs()
    rows = []
    for run in experiments:
        spec = EXPERIMENTS[run]
        for hemi, hemi_years in spec["years"].items():
            if hemi not in hemispheres:
                continue
            for year in hemi_years:
                if years_filter is not None and year not in years_filter:
                    continue
                conv = convert_done(run, year, hemi)
                vort = vort_ready(run, year, hemi)
                done = fully_done(run, year, hemi)
                is_running = (
                    (run, str(year), hemi) in running_precise
                    or (run, str(year)) in running_year_fallback
                )
                if done:
                    state = "ok"
                elif is_running:
                    state = "running"
                elif not conv:
                    state = "needs_convert"
                else:
                    state = "needs_track"
                rows.append((run, year, hemi, conv, vort, done, state))
    return rows


def print_status(rows):
    header = "%-8s %-6s %-4s %-9s %-11s %-6s %s" % (
        "run", "year", "hemi", "convert", "vort_ready", "done", "state",
    )
    print(header)
    print("-" * len(header))
    for run, year, hemi, conv, vort, done, state in rows:
        print(
            "%-8s %-6s %-4s %-9s %-11s %-6s %s"
            % (run, year, hemi, conv, vort, done, state)
        )
    n_total = len(rows)
    n_ok = sum(1 for r in rows if r[-1] == "ok")
    n_running = sum(1 for r in rows if r[-1] == "running")
    n_convert = sum(1 for r in rows if r[-1] == "needs_convert")
    n_track = sum(1 for r in rows if r[-1] == "needs_track")
    print(
        "\n%d total: %d ok, %d running, %d need convert(+track), %d need track only"
        % (n_total, n_ok, n_running, n_convert, n_track)
    )


# ---------------------------------------------------------------------------
# Submit
# ---------------------------------------------------------------------------
def months_for(hemi, year):
    """(year, month) pairs to convert for this (hemisphere, year) - matches
    submit_SLURMscript_n1280_para_WCSI's NH (Jan-Dec of year) / SH (Jul-Dec
    of year, Jan-Jun of year+1) split.
    """
    if hemi == "NH":
        return [(year, m) for m in range(1, 13)]
    return [(year, m) for m in range(7, 13)] + [(year + 1, m) for m in range(1, 7)]


SBATCH_JOBID_RE = re.compile(r"Submitted batch job (\d+)")


def sbatch(args, dry_run):
    """Runs sbatch (unless dry_run, which only prints it) and returns the new
    job's id, or None if dry_run or the id couldn't be parsed from sbatch's
    own output.
    """
    cmd = ["sbatch"] + args
    print(" ".join(cmd))
    if dry_run:
        return None
    result = subprocess.run(cmd, check=True, cwd=SCRIPTS_DIR, capture_output=True, text=True)
    print(result.stdout.strip())
    match = SBATCH_JOBID_RE.search(result.stdout)
    return match.group(1) if match else None


def submit_convert_and_track(run, year, hemi, dry_run):
    convert_name = "%s_convert_%s" % (run, year)
    # EERIE's raw monthly files are much larger (one pc-stream month is
    # ~100GB, several times the other experiments') - bump the conversion
    # job's memory above SLURMscript_n1280_para_WCSI's default 64GB as a
    # starting guess; adjust based on how the first real run actually does.
    convert_args = ["--mem=128000"] if run in EERIE_EXPERIMENTS else []
    job_ids = []
    for m_year, month in months_for(hemi, year):
        jid = sbatch(
            convert_args + ["--job-name=%s" % convert_name, "./SLURMscript_n1280_para_WCSI", run, str(m_year), str(month)],
            dry_run,
        )
        if jid:
            job_ids.append(jid)

    # The track job depends on exactly these convert jobs by id (afterany:
    # runs once they've all finished, regardless of their exit status - same
    # "wait, don't gate on success" behaviour as the old singleton scheme;
    # job ids rather than a shared name so convert/track can have distinct,
    # readable names without losing the dependency).
    if dry_run:
        print(
            "(track job below depends on the %d convert job(s) above via --dependency=afterany)"
            % len(job_ids or months_for(hemi, year))
        )
        dependency = None
    else:
        dependency = "--dependency=afterany:" + ":".join(job_ids) if job_ids else None
    submit_track_only(run, year, hemi, dry_run, jobname="%s_track_%s" % (run, year), dependency=dependency)


def submit_track_only(run, year, hemi, dry_run, jobname=None, dependency=None):
    jobname = jobname or ("%s_track_%s" % (run, year))
    args = []
    if dependency:
        args.append(dependency)
    args += [
        "--job-name=%s" % jobname, "./SLURMscript_n1280_track_para_WCSI",
        run, str(year), TRUNC, scrdir(run), permdir(run),
        "True", "addvor", "False", HEMI_ARG[hemi],
    ]
    sbatch(args, dry_run)


def submit(rows, dry_run):
    todo = [r for r in rows if r[-1] in ("needs_convert", "needs_track")]
    if not todo:
        print("Nothing to submit.")
        return
    for run, year, hemi, conv, vort, done, state in todo:
        if state == "needs_convert":
            submit_convert_and_track(run, year, hemi, dry_run)
        else:
            submit_track_only(run, year, hemi, dry_run)
    if dry_run:
        print("\n(dry run - pass --yes to actually submit these %d job(s))" % len(todo))
    else:
        print("\nSubmitted jobs for %d (experiment, year, hemisphere) combination(s)." % len(todo))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_years(spec):
    if spec is None:
        return None
    years = set()
    for part in spec.split(","):
        if "-" in part:
            lo, hi = part.split("-")
            years.update(range(int(lo), int(hi) + 1))
        else:
            years.add(int(part))
    return years


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["status", "submit"])
    parser.add_argument("--experiments", help="comma-separated, default all: " + ",".join(EXPERIMENTS))
    parser.add_argument("--years", help="e.g. 1984-1990,2000")
    parser.add_argument("--hemisphere", help="NH,SH (default both)")
    parser.add_argument("--yes", action="store_true", help="submit: actually run sbatch (default: dry run)")
    args = parser.parse_args()

    experiments = args.experiments.split(",") if args.experiments else list(EXPERIMENTS)
    for run in experiments:
        if run not in EXPERIMENTS:
            parser.error("unknown experiment %r (known: %s)" % (run, ", ".join(EXPERIMENTS)))
    years_filter = parse_years(args.years)
    hemispheres = args.hemisphere.split(",") if args.hemisphere else ["NH", "SH"]

    rows = audit(experiments, years_filter, hemispheres)
    print_status(rows)

    if args.command == "submit":
        print()
        submit(rows, dry_run=not args.yes)


if __name__ == "__main__":
    main()
