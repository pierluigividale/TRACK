"""Read-only status audit across the TC-tracking ensemble on JASMIN.

For each (experiment, year, hemisphere) in scope, reports whether the final
TRACK/Hart output exists, whether a SLURM job for it is currently running,
and - if neither - how far the per-month raw conversion got (the slow,
expensive stage) so you know whether a resubmit is cheap (just the TRACK/
Hart stage) or will have to redo raw pp conversion for some months too.

This does not submit or delete anything - it only reads the filesystem and
squeue. Run from ~/TRACK-1.5.2/scripts/ on jasmin-sci:

    python3 tc_status.py [experiment ...]

With no arguments, checks all experiments in EXPERIMENTS below.
"""
import glob
import os
import subprocess
import sys

HOME = os.path.expanduser("~/TRACK-1.5.2/")
SCRATCH = "/work/scratch-pw4/plvidale/"

# hemispheres + year ranges confirmed this session (NH years run Jan-Dec;
# SH years run Jul(year)-Jun(year+1), so SH's last usable year is one less
# than NH's when both share the same raw-data end point).
EXPERIMENTS = {
    "ch330": dict(hemispheres=["NH"], years=range(1980, 2015)),
    "ck777": dict(hemispheres=["NH"], years=range(1980, 2015)),
    "ck778": dict(hemispheres=["NH"], years=range(1980, 2015)),
    "dd828": dict(hemispheres=["NH", "SH"], years=range(1980, 2015)),
    "dc009": dict(hemispheres=["NH", "SH"], years=range(1980, 2015)),
}

MONTH_VARS = ["_uwinds_gg", "_vwinds_gg", "_mslp", "_geo", "_10m_wspeed"]


def running_jobs():
    try:
        out = subprocess.check_output(
            ["squeue", "-u", "plvidale", "-h", "-o", "%j"], text=True
        )
    except Exception:
        return set()
    return set(out.split())


def final_file(run, year, hemi):
    hemi_tag = "_SH" if hemi == "SH" else ""
    trackname = "tr_trs_neg" if hemi == "SH" else "tr_trs_pos"
    return os.path.join(
        HOME, "u-" + run, "VOR_VERTAVG_T63filt_%s%s" % (year, hemi_tag),
        trackname + ".2day_addvor_addmslp_addwind.hart.new.nc",
    )


def months_present(run, year, hemi):
    """(n_months_fully_converted, n_months_total) for this run+year+hemi's
    raw-conversion stage, counting a calendar month as done only when all 5
    MONTH_VARS exist for it. NH checks Jan-Dec(year); SH checks
    Jul(year)-Jun(year+1), matching link_variables()'s own windows.
    """
    scrdir = os.path.join(SCRATCH, "u-" + run) + "/"
    if hemi == "NH":
        months = [(year, m) for m in range(1, 13)]
    else:
        months = [(year, m) for m in range(7, 13)] + [(year + 1, m) for m in range(1, 7)]
    done = 0
    for y, m in months:
        tag = "%s_%04d%02d" % (run, y, m)
        if all(os.path.exists(scrdir + tag + v + ".nc") for v in MONTH_VARS):
            done += 1
    return done, len(months)


def main():
    experiments = sys.argv[1:] or list(EXPERIMENTS)
    running = running_jobs()

    for run in experiments:
        cfg = EXPERIMENTS[run]
        for hemi in cfg["hemispheres"]:
            missing_years = []
            running_years = []
            for year in cfg["years"]:
                job_name_guess = "%s_track_%s%s" % (run, year, "_SH" if hemi == "SH" else "")
                if os.path.exists(final_file(run, year, hemi)):
                    continue
                if job_name_guess in running or any(job_name_guess in j for j in running):
                    running_years.append(year)
                    continue
                missing_years.append(year)

            total = len(cfg["years"])
            n_missing = len(missing_years)
            n_running = len(running_years)
            n_done = total - n_missing - n_running
            print("[%s %s] %d/%d years done" % (run, hemi, n_done, total), end="")
            if n_running:
                print(", %d running" % n_running, end="")
            if n_missing:
                print(", %d missing" % n_missing)
                # for missing years, report conversion progress so you know
                # whether a resubmit is cheap (TRACK/Hart only) or needs
                # raw-conversion work redone for some months
                for year in missing_years:
                    done, tot = months_present(run, year, hemi)
                    flag = "" if done == tot else "  <-- conversion incomplete"
                    print("    %s %s %4d: convert %d/%d months%s" % (run, hemi, year, done, tot, flag))
            else:
                print()


if __name__ == "__main__":
    main()
