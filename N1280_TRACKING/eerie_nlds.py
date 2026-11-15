"""Shared EERIE raw-data / NLDS helpers for the TRACK WCSI workflow.

Used by:
  pp_to_nc_iris_WCSI_EERIE.py  (SLURM convert job, on a LOTUS compute node)
  tc_workflow.py stage/submit  (run interactively on a sci node)

Why staging is a separate step (found 2026-11-12): the NLDS server rejects
our token with HTTP 403 when called from LOTUS compute nodes, while the very
same ~/.nlds-token works from sci-vm-0x. So all NLDS calls (find/get/stat)
happen on a sci node via `tc_workflow.py stage`, which retrieves raw pp files
into STAGING_ROOT and marks each one verified with a ".complete" file. The
convert job only ever reads from the original location or from a verified
staged copy - it never talks to NLDS - and deletes the staged month after a
successful conversion.

The EERIE runs are ~100 years long and their raw pp files are removed from
disk on a rolling basis once archived (holdings ingested by the hiresgw
group, e.g. 5455 "hrcm__n1280o12__spinup__u-dz876" owned by lawrence). There
is deliberately no hardcoded holding id: files are located across ALL hiresgw
group holdings with `nlds find -A -p <path>`.
"""

import json
import os
import re
import shutil
import subprocess
import time

MONTH_ABBR = {
    1: "jan", 2: "feb", 3: "mar", 4: "apr", 5: "may", 6: "jun",
    7: "jul", 8: "aug", 9: "sep", 10: "oct", 11: "nov", 12: "dec",
}

RAW_ROOT = "/gws/ssde/j25a/hiresgw/aosprey/"
SCRATCH_ROOT = "/work/scratch-pw4/plvidale/"

NLDS_BIN = os.path.expanduser("~/.local/bin/nlds")
# Must be on a GWS (NLDS can't see /home). The per-month target dir must be
# group/world-writable: with plain drwxr-sr-x NLDS fails with "[Errno 13]
# Permission denied" (verified 2026-11-12; chmod 2777 fixed it). NLDS
# recreates the FULL original path underneath the target dir.
STAGING_ROOT = "/gws/ssde/j25a/hrcm/plvidale/nlds_staging/"
# Never stage if the hrcm GWS would be left with less than this free
# (gws_hrcm is ~96% of its 1PB quota). A month of p8+pc is ~450GB.
STAGING_MIN_FREE_BYTES = 5 * 1024**4
NLDS_TERMINAL_STATES = {"COMPLETE", "FAILED", "COMPLETE_WITH_ERRORS", "COMPLETE_WITH_WARNINGS"}
NLDS_ALL_STATES = NLDS_TERMINAL_STATES | {
    "INITIALISING", "ROUTING", "SPLITTING", "INDEXING", "CATALOG_PUTTING",
    "TRANSFER_PUTTING", "CATALOG_ROLLBACK", "CATALOG_GETTING", "ARCHIVE_GETTING",
    "TRANSFER_GETTING", "ARCHIVE_INIT", "CATALOG_ARCHIVE_AGGREGATING",
    "ARCHIVE_PUTTING", "CATALOG_ARCHIVE_UPDATING", "CATALOG_ARCHIVE_ROLLBACK", "SPLIT",
}


# ---------------------------------------------------------------------------
# Raw / output paths (no NLDS - safe on compute nodes)
# ---------------------------------------------------------------------------
def cycle_date(run, year, mon):
    """Each month's data is archived under a directory dated one month ahead
    of the data itself, e.g. 1950jan data sits under u-dz876/19500201T0000Z/."""
    year, mon = int(year), int(mon)
    if run == "dz876" and year == 1959 and mon == 12:
        # this one cycle bundled Nov+Dec 1959 output into the Nov cycle dir
        # instead of giving Dec its own - the normal +1-month dir (19600101)
        # has no atmosphere pp files at all; dz876a.p81959dec/pc1959dec etc.
        # actually sit in 19591201T0000Z alongside the 1959nov files.
        return "19591201T0000Z"
    dy, dm = year, mon + 1
    if dm == 13:
        dy, dm = dy + 1, 1
    return "%d%02d01T0000Z" % (dy, dm)


def raw_paths(run, year, mon):
    """{stream: original raw pp path}. p8 (~360GB): u/v winds, geo, mslp;
    pc (~100GB): 10m wind speed."""
    d = "%su-%s/%s/" % (RAW_ROOT, run, cycle_date(run, year, mon))
    tag = "%d%s" % (int(year), MONTH_ABBR[int(mon)])
    return {"p8": d + "%sa.p8%s" % (run, tag), "pc": d + "%sa.pc%s" % (run, tag)}


_CYCLE_DIR_RE = re.compile(r"^(\d{4})(\d{2})01T0000Z$")
# a raw file modified more recently than this may still be arriving from ARCHER2
RAW_SETTLE_SECONDS = 2 * 3600


def latest_raw_month(run, settle_seconds=RAW_SETTLE_SECONDS):
    """(year, month) of the newest month whose raw p8 AND pc files are both in
    RAW_ROOT/u-<run>/ (non-empty, not modified in the last settle_seconds, so
    a file still being transferred from ARCHER2 doesn't count) - i.e. how far
    the run has been delivered to JASMIN. None if there is none.

    Older months may since have been removed from disk (they are then on NLDS,
    which `stage` handles), but the newest output is always the last to go, so
    this is a reliable "latest available" marker for the year range.
    """
    top = "%su-%s/" % (RAW_ROOT, run)
    try:
        cycles = sorted((d for d in os.listdir(top) if _CYCLE_DIR_RE.match(d)), reverse=True)
    except OSError:
        return None
    now = time.time()
    for d in cycles:
        cy, cm = (int(x) for x in _CYCLE_DIR_RE.match(d).groups())
        year, mon = (cy, cm - 1) if cm > 1 else (cy - 1, 12)  # dir is dated one month ahead
        ok = True
        for path in raw_paths(run, year, mon).values():
            try:
                st = os.stat(path)
            except OSError:
                ok = False
                break
            if st.st_size == 0 or now - st.st_mtime < settle_seconds:
                ok = False
                break
        if ok:
            return year, mon
    return None


def dirout(run):
    return "%su-%s/" % (SCRATCH_ROOT, run)


def outputs(run, year, mon):
    """Monthly output file names written by pp_to_nc_iris_WCSI_EERIE.py (and
    then finalised by ll_to_gg_n1280.py)."""
    base = "%s%s_%d%02d_" % (dirout(run), run, int(year), int(mon))
    return {
        "u": base + "uwinds.nc", "u_gg": base + "uwinds_gg.nc",
        "v": base + "vwinds.nc", "v_gg": base + "vwinds_gg.nc",
        "geo": base + "geo.nc", "mslp": base + "mslp.nc", "wsp": base + "10m_wspeed.nc",
    }


def needed_vars(run, year, mon):
    """{var: True if still to be decoded from raw pp}. u/v count as done if
    either the regridded _gg file (ll_to_gg's final product) or the
    intermediate exists - intermediates are written atomically (temp name +
    rename), so an existing one is always complete."""
    o = outputs(run, year, mon)
    ex = os.path.exists
    return {
        "u": not (ex(o["u_gg"]) or ex(o["u"])),
        "v": not (ex(o["v_gg"]) or ex(o["v"])),
        "geo": not ex(o["geo"]),
        "mslp": not ex(o["mslp"]),
        "wsp": not ex(o["wsp"]),
    }


def needed_streams(run, year, mon):
    """{stream: original raw path} for the raw streams this month still
    needs - nothing if the month is already converted."""
    need = needed_vars(run, year, mon)
    paths = raw_paths(run, year, mon)
    out = {}
    if need["u"] or need["v"] or need["geo"] or need["mslp"]:
        out["p8"] = paths["p8"]
    if need["wsp"]:
        out["pc"] = paths["pc"]
    return out


def staging_month_dir(run, year, mon):
    return "%su-%s/%d%02d" % (STAGING_ROOT, run, int(year), int(mon))


def staging_target_dir(run, year, mon, original_path):
    """One NLDS target dir PER FILE, never shared between GETs. Found
    2026-11-12: with p8 and pc both targeting the same dir, the first GET to
    finish (pc) had NLDS hand the directories it created over to plvidale
    (drwxrwx--- plvidale:users); the still-running p8 GET - which writes as
    user `nlds` - then got "[Errno 13] Permission denied" renaming its
    .part.minio file to the final name, and its retry couldn't even enter
    the tree."""
    return staging_month_dir(run, year, mon) + "/" + os.path.basename(original_path)


def staged_path(run, year, mon, original_path):
    # NLDS recreates the full original path under the target dir
    return staging_target_dir(run, year, mon, original_path) + original_path


def is_staged(run, year, mon, original_path):
    """A staged copy only counts once `stage` has verified it (size checked
    against the NLDS catalogue) and written the .complete marker."""
    p = staged_path(run, year, mon, original_path)
    return os.path.exists(p) and os.path.exists(p + ".complete")


def locate(run, year, mon, original_path):
    """Readable path for a raw file: original location, else verified staged
    copy, else None."""
    if os.path.exists(original_path):
        return original_path
    if is_staged(run, year, mon, original_path):
        return staged_path(run, year, mon, original_path)
    return None


def remove_staged_month(run, year, mon, log=print):
    d = staging_month_dir(run, year, mon)
    if os.path.exists(d):
        log("removing staged NLDS copy %s" % d)
        shutil.rmtree(d, onerror=lambda f, p, e: log(
            "WARNING: could not remove %s: %s - delete by hand" % (p, e[1])))


def staging_usage_bytes(top=STAGING_ROOT):
    total = 0
    for root, _, files in os.walk(top):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return total


# ---------------------------------------------------------------------------
# NLDS client wrappers (sci nodes only - 403 from LOTUS compute nodes)
# ---------------------------------------------------------------------------
def _nlds(*args):
    return subprocess.run([NLDS_BIN, *args], capture_output=True, text=True)


def nlds_find(original_path):
    """(holding_id, size_bytes, storage_types) for original_path across all
    hiresgw group holdings, or None if it isn't archived."""
    r = _nlds("find", "-A", "-p", original_path, "-j")
    if r.returncode != 0 or not r.stdout.strip().startswith("{"):
        if "No matching files" in (r.stdout + r.stderr):
            return None
        raise RuntimeError("nlds find failed for %s:\n%s%s" % (original_path, r.stdout, r.stderr))
    for holding in json.loads(r.stdout)["data"]["holdings"].values():
        for tr in holding["transactions"].values():
            for f in tr["filelist"]:
                if f["original_path"] == original_path:
                    return (f["holding_id"], f["size"],
                            sorted({l["storage_type"] for l in f["locations"]}))
    return None


def job_label(run, year, mon, stream):
    # deterministic, so a re-run of `stage` recognises its own in-flight GETs
    return "stage_%s_%d%02d_%s" % (run, int(year), int(mon), stream)


def nlds_get(original_path, holding_id, target_dir, label):
    """Submit one GET (non-blocking) into a fresh per-file target dir (any
    leftover from an earlier failed attempt is cleared first - its subdirs
    may no longer be writable by NLDS). Returns the transaction id."""
    if os.path.exists(target_dir):
        shutil.rmtree(target_dir)
    os.makedirs(os.path.dirname(target_dir), exist_ok=True)
    os.makedirs(target_dir)
    os.chmod(target_dir, 0o2777)  # NLDS needs write access, see STAGING_ROOT
    r = _nlds("get", original_path, "-A", "-i", str(holding_id),
              "-r", target_dir + "/", "-b", label, "-j")
    if r.returncode != 0:
        raise RuntimeError("nlds get submission failed for %s:\n%s%s"
                           % (original_path, r.stdout, r.stderr))
    return json.loads(r.stdout)["transaction_id"]


_STATE_KV_RE = re.compile(r"state\s*:\s*(\S+)")


def latest_state(label):
    """State of the most recent GET with this job label, or None if there
    has never been one. Handles both of `nlds stat`'s output layouts
    (key/value for a single transaction, table for several)."""
    r = _nlds("stat", "-b", label, "-9", "-L", "1")
    out = r.stdout
    m = _STATE_KV_RE.search(out)
    if m:
        return m.group(1)
    for line in out.splitlines():
        if label in line.split():
            for tok in line.split():
                if tok in NLDS_ALL_STATES:
                    return tok
    return None


def errors_for(label):
    return _nlds("stat", "-b", label, "-9", "-L", "1", "-E").stdout
