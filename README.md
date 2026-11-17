# TRACK
A repository for my TRACK scripts
I have deposited a main script to run track, as well as a couple of service scripts, to prepare the data needed by TRACK. 
I have also deposited a number of namelists that are read by the main script and used to drive TRACK.

# The full workflow of Tropical Cyclone tracking and analysis

End-to-end workflow for identifying and analysing tropical cyclones in
high-resolution (N1280, ~10 km) atmosphere and coupled atmosphere-ocean model
output: TRACK-1.5.2 feature tracking on JASMIN, then TC identification and
ensemble analysis with [huracanpy](https://huracanpy.readthedocs.io/) on a
local machine.

This document is the map of the pipeline - what runs where, in what order,
and what a new simulation needs before it can be added. Script-level detail
(exact flags, file formats, known gotchas) lives in each script's own
docstring; this README tells you which script to open.

**Please note**: all scripts described here are contained in the N1280_TRACKING subdirectory on this repository. You can either copy them from PLV's home directory on JASMIN, under `/home/users/plvidale/TRACK-1.5.2/scripts`, or else download them from N1280_TRACKING on this repository.

## Quick start: TC identification on TRACK output you already have

If you already have TRACK output with the Hart parameters added (the
`tr_trs_pos*.hart.new.nc` / `tr_trs_neg*.hart.new.nc` files), you only need
two scripts, which work on a laptop or on JASMIN, for any TRACK output (our
N1280 runs, reanalyses, DestinE, ...):

```
pip install --upgrade git+https://github.com/Huracan-project/wcsi   # wcsi >= 2 required
pip install huracanpy parse pyarrow cartopy

# 1. method-6 TC identification -> {experiment}_{year}_{domain}_TCs.nc + _summary.parquet
python huracanpy_identify_method6.py --experiments <exp ...> --years 1990-2014 \
    --data-dir <folder holding the TRACK output> --out-dir <output folder> \
    --input-nh "<NH file pattern, with {experiment} and {year}>" \
    --input-sh "<SH file pattern; SH year Y = July Y - June Y+1>" \
    --warm-core-level <top vorticity level in the files, hPa>

# 2. a PNG map and a CSV of the tracks (Hart parameters, WCSI nature) per year and domain
python huracanpy_method6_maps_csv.py <output folder>

# 3. the same, but all years of each experiment and domain together
#    (one map, one CSV and one joined NetCDF file per experiment and domain)
python huracanpy_method6_maps_csv.py <output folder> --combine
python huracanpy_method6_maps_csv.py <output folder> --combine --years 1990-1999 --experiments hist
```

Worked example - Kevin Hodges' DestinE IFS-NEMO2 tracks on JASMIN (MSLP field
`psl`, vorticity up to 200 hPa):

```
python huracanpy_identify_method6.py --experiments hist cont --years 1990-2014 \
    --data-dir /gws/ssde/j25b/cmip6_track/DESTINe/IFS-NEMO2 --out-dir <output folder> \
    --warm-core-level 200 \
    --input-nh "{experiment}/TC/NH/IFS-NEMO_{experiment}_VOR_vertavg_{year}_T63/tr_trs_pos*.hart.new.nc" \
    --input-sh "{experiment}/TC/SH/IFS-NEMO_{experiment}_VOR_vertavg_jul-jun{year}*_T63/tr_trs_neg*.hart.new.nc"
```

Check the first line printed for each file: it names the variable used as the
MSLP-minimum position (detected automatically) and the percentage of points
with no MSLP minimum - typically 20-35%; 0% means the wrong variable. Years
already present in the output folder are skipped - use `--overwrite` or a new
`--out-dir` after changing anything. Do **not** use the older
`huracanpy_identify_batch_all.py` for identification (see Part 3).

## Two environments

The pipeline spans two places that don't share a filesystem:

- **JASMIN** (`TRACK-1.5.2/`, under `/home/users/<user>/`): converts raw UM
  (Unified Model) or EERIE-project pp output to NetCDF and runs TRACK-1.5.2
  itself, via SLURM batch jobs. This is where the compute-heavy tracking
  happens.
- **Local machine** (`TC-TRACKS/`, this repo): TC identification (which
  tracked vortices are tropical cyclones) and all plotting, with
  [huracanpy](https://huracanpy.readthedocs.io/).

Tracked output has to be `rsync`'d from JASMIN to the local machine's
`TC-TRACKS/` before any of the analysis scripts can see it.

## Part 1 - Tracking on JASMIN

### 1.1 Convert raw model output to NetCDF

`scripts/pp_to_nc_iris_WCSI.py <experiment> <year> <month>` (one SLURM job per
month) extracts the fields TRACK needs - u/v winds at several pressure
levels, mean sea level pressure, 10 m wind speed, and geopotential height at
three pressure levels - from the raw `.pp` files, using `iris`.

Every modelling centre/run lays out its raw archive differently, so the
script dispatches on a per-experiment `SOURCE_GROUP` table at the top of the
file: which directory the `.pp` files live in, and whether the level fields
(winds, geopotential) come from the same stream as everything else or a
separate one. **Adding a new standard (atmosphere-only) experiment means
adding one entry to this table** - look at the existing entries for the
pattern, and verify empirically (don't guess) which stream each field is
actually in: `iris.load_cube` raising `ConstraintMismatchError`, or a field
quietly coming from the wrong stream, are both easy mistakes here.

Coupled atmosphere-ocean (EERIE-project) experiments use a structurally
different raw-archive convention (`pp_to_nc_iris_WCSI_EERIE.py`): per-month
cycle directories offset one month from the data they contain, a third
stream, and - important - **different geopotential pressure levels** (see
1.3 below). `scripts/SLURMscript_n1280_para_WCSI` dispatches to this script
for the experiments in `tc_workflow.py`'s `EERIE_EXPERIMENTS` set.

`scripts/ll_to_gg_n1280.py` then resamples everything to a common
6-hourly cycle and regrids the winds to the Gaussian grid TRACK expects. If
10 m wind speed is natively on an offset cycle (e.g. 03/09/15/21 UTC instead
of 00/06/12/18) for a given experiment, add it to this script's
`WSPEED_SELHOUR` table rather than letting the resample step fail - this is
tolerable for TRACK's purposes (confirmed with Kevin Hodges), just needs the
right target hours so the step doesn't abort.

### 1.2 Run TRACK-1.5.2

`scripts/track_n1280_para_WCSI.py <experiment> <year> <truncation> <scrdir>
<permdir> ...` runs the actual tracker (vorticity filtering, `master`, the
Hart cyclone-phase-space calculation, `tcident`) for one experiment/year/
hemisphere. **Before tracking any experiment, both of these must exist
in `indat/` at the same time** (this is a standing setup, not created per
run):

```
indat/RUNDATIN.vor_global_hilat_NH_T63.in
indat/RUNDATIN.vor_global_hilat_SH_T63_A.in
```

`master` decides which sign of vorticity to track - positive (NH cyclonic)
or negative (SH cyclonic) - by checking indat/ for an `_A`-suffixed variant
of the base RUNDATIN name it's given.
Pass it the **base** name for both hemispheres (the SH call in
`track_n1280_para_WCSI.py` does *not* add `_A` itself); if it's passed the
`_A`-suffixed name directly, it probes for a nonexistent `<name>_A_A` file
and silently produces no track output at all for the SH. Confirmed the hard
way - don't change this without checking with Kevin first.

### 1.3 Non-standard pressure levels (coupled/EERIE experiments)

The Hart phase-space step's templates (`regZ_NH.in`, `regZ_SH.in`,
`hart_NH.in`, `hart_SH.in`) hardcode geopotential-height levels 900/600/300
hPa, which match the standard atmosphere-only experiments' extraction but
**not** every raw archive. If a new experiment's `.pp` files only offer
different levels (confirmed with Leo Saffin that the Hart parameters are
still valid computed at different levels, at some cost to cross-experiment
comparability), make `_EERIE`-suffixed copies of all four templates with the
levels substituted, and extend `EERIE_EXPERIMENTS` and the `eerie_suffix`
dispatch in `track_n1280_para_WCSI.py`'s `Hart()` function to use them for
that experiment. (This set is separate from, but should stay in sync with,
`tc_workflow.py`'s own `EERIE_EXPERIMENTS`, which controls the *conversion*
dispatch in 1.1.)

### 1.4 Orchestration

Running 1.1-1.2 by hand across many years/hemispheres/experiments doesn't
scale, and parts of the pipeline fail in ways SLURM itself won't flag (a
wrapper script can exit 0 even after the Python step inside it crashed).
`scripts/tc_workflow.py` is the audit-and-submit tool:

```
python tc_workflow.py status [--experiments ...] [--years ...] [--hemisphere NH|SH]
python tc_workflow.py stage  [--experiments ...] [--years ...] [--yes] [--wait]   # EERIE only, sci node
python tc_workflow.py submit [--experiments ...] [--years ...] [--hemisphere NH|SH] [--yes]
```

It submits one conversion job per month (`SLURMscript_n1280_para_WCSI`) and
a tracking job per year (`SLURMscript_n1280_track_para_WCSI`) that waits for
them.

**EERIE runs (dz876, ec100) - tape retrieval and automatic year ranges.**
These ~100-year coupled runs are still being delivered from ARCHER2, and their
raw pp files are removed from disk on a rolling basis once archived to
JASMIN's NLDS tape store. Two consequences:

- `stage` retrieves whatever raw files a year needs from NLDS into a staging
  area, verifies them and marks them complete; `submit` refuses an EERIE year
  until every raw file is on disk or staged. NLDS rejects requests from LOTUS
  compute nodes (HTTP 403), so `stage` must be run interactively on a sci
  node (`sci-vm-0x`); the conversion jobs never call NLDS themselves.
  NLDS helpers are in `eerie_nlds.py`.
- Their year range is not hand-set (`years="auto"` in `EXPERIMENTS`): it runs
  from a fixed start year up to the last complete NH year (needs December) /
  SH season (July-June, needs June) in the raw data delivered so far,
  judged from the newest month whose files are in
  `/gws/ssde/j25a/hiresgw/aosprey/u-<run>/`. `status` prints the latest month
  found; just re-run `status` and `submit` as ARCHER2 delivers more.

`status` is read-only and reports, per (experiment, year, hemisphere),
whether conversion, vorticity filtering and final tracking actually
succeeded (by checking for the real output files, not SLURM's own job
state). `submit` is a dry run by default - it prints the `sbatch` calls it
would make - and only actually submits with `--yes`. Adding a new
experiment means adding one entry to its `EXPERIMENTS` dict (source group,
year range per hemisphere) - everything else is driven from that.

Two things worth knowing before trusting a `status` "done" for a brand new
experiment: (1) **`sacct`/SLURM "COMPLETED" is not a reliable success
signal** for this pipeline - always check the real output files, which is
what `status` does; (2) a short/truncated final year in the raw archive
(the way we have defined TRACK's chunking in this workflow can't handle a partial year right now, but TRACK is fully capable of doing that, so this is an important modification to come to the script in the coming weeks) will fail silently rather
than erroring clearly - check the raw archive actually covers the full year
range before assuming a "missing" year is a pipeline bug.

## Part 2 - Bringing tracked output to the analysis machine

```
rsync -avP --exclude 'DJF*' \
    <jasmin-host>:/home/users/<user>/TRACK-1.5.2/u-<experiment>/VOR_* \
    TC-TRACKS/u-<experiment>/
```

(`--exclude 'DJF*'`, not `'DJF\*'` - inside single quotes the backslash is
literal and won't match anything.) The `DJF_MAX_*`/`DJF_MIN_*`
subdirectories are large intermediate output not needed for identification
or analysis.

Every later step expects tracked output at
`TC-TRACKS/u-<experiment>/VOR_VERTAVG_T<trunc>filt_<year>[_SH]/`, matching
the layout JASMIN itself uses - `rsync`ing `VOR_*` straight into
`TC-TRACKS/u-<experiment>/` reproduces that automatically.

## Part 3 - TC identification

TRACK finds vortices; a separate identification step decides which of them
are tropical cyclones. Several identification approaches were compared; the
one adopted (suggested by Kevin Hodges and Leo Saffin) is **method 6, "WCSI
warm core only"**:

- Hart cyclone-phase-space warm-core/symmetric criteria (Leo Saffin's
  [`wcsi`](https://github.com/Huracan-project/wcsi) package, version 2 or
  later), with **no** minimum 850 hPa vorticity threshold;
- a warm-core vorticity test: 850 hPa minus upper-level relative vorticity
  greater than 6e-5 s-1 (tcident's 6.0). The upper level is the top level in
  the TRACK files - 250 hPa for our N1280 runs, 200 hPa for the DestinE
  IFS-NEMO2 runs (`--warm-core-level`);
- the point must be over ocean and have a closed MSLP minimum (a cyclone,
  not just a vortex);
- all of the above for 24 consecutive hours, and the first TC point
  equatorward of 50 degrees.

Each TC track is assigned to one ocean basin (where its vorticity peaks
during the TC phase; Australian-region and Central Pacific storms are merged
into SIO/SPAC and EPAC as IBTrACS counts them).

**`huracanpy_identify_method6.py`** does all of this in one self-contained
script (no imports from the other scripts here). Settings are at the top of
the file and can all be overridden on the command line (`--help`; see the
Quick start above for an example). It writes every experiment's output into
one shared directory (`TC-TRACKS/METHOD6_ENSEMBLE_RUN/` by default):
`{experiment}_{year}_{domain}_TCs.nc` (the TC tracks, with `is_tc` and the
WCSI `nature` of every point) and `{experiment}_{year}_{domain}_summary.parquet`
(one row per track) - the layout Part 4's scripts expect. Already-processed
experiment/years are skipped unless `--overwrite` is given.

**`huracanpy_method6_maps_csv.py <output folder>`** then writes, for each
year and domain with at least one TC, a map (`_nature.png`: marker = WCSI
nature, colour = MSLP, P/L = first/last point) and a CSV of the track points
(`_tracks.csv`: position, vorticity at every level, MSLP and 10 m wind with the
positions TRACK found them at, Hart B/VTL/VTU, `is_tc`, `nature`; TRACK's 1e25
missing values are left empty). It only reads the NetCDF output, so it is
quick to re-run; existing files are skipped unless `--overwrite`. Options:

| Option | Effect |
|---|---|
| `--experiments hist cont` | only these experiments |
| `--years 1990-2014` or `--years 1990,2000` | only these years (inclusive range or list) |
| `--combine` | instead of one set of files per year, join all (selected) years of each experiment and domain: `{experiment}_{y0}-{y1}_{domain}_nature.png`, `_tracks.csv` (with an extra `year` column) and `_allTCs.nc` (the joined tracks). Useful when there are few TCs per year, e.g. DestinE. |
| `--no-maps`, `--no-csv` | skip the maps or the CSVs |
| `--overwrite` | redo existing files |

With `--combine`, track ids are renumbered year x 1000000 + the original TRACK
id, so they stay unique across years. The log lists years with zero TCs in that
domain (processed, but no `_TCs.nc` file because there was nothing to save)
separately from years that were never processed.

**NetCDF format.** All `_TCs.nc` and `_allTCs.nc` files are written by
huracanpy (`huracanpy.save`) in the standard CF-convention layout for
trajectories (a "contiguous ragged array", `featureType = "trajectory"`):
every track point along a `record` dimension (time, lon, lat, MSLP,
`relative_vorticity(record, pressure)`, `cps_b`/`cps_vtl`/`cps_vtu`, `is_tc`,
`nature`, ...), and per track along a `trajectory` dimension its `track_id`
and `rowSize` (number of points). `huracanpy.load(file)` reads them back with
`track_id` repeated on every point; other tools rebuild it by repeating each
`track_id` `rowSize` times. The CSVs already have it on every row.

`huracanpy_identify_batch_all.py` is the earlier identification script. It
is kept only because the Part 4, ERA5 and IBTrACS scripts import shared
definitions from it (basins, paths, helpers); running it needs the old wcsi
1.x and fails with wcsi 2 ("'function' object has no attribute 'wcsi'").

The ensemble/interannual/genesis scripts below also compare against ERA5
and JRA3Q reanalysis "pseudo-tracks" (`ERA5_<year>_<domain>_TCs.nc` etc.,
from `huracanpy_era5_batch_all.py`) and IBTrACS observations. The
reanalyses' TC classification is pre-computed once, independent of which
identification approach the *models* are run through (see
`huracanpy_era5_batch_all.py`'s docstring) - copy the existing
`ERA5_*`/`JRA3Q_*` files from another run directory (e.g. `NEW_WCSI_RUN`)
into the shared method-6 directory rather than regenerating them.

## Part 4 - Ensemble analysis

Three scripts, all driven by the same environment variables so they stay
pointed at the same run and the same subset of basins/years (plus the
reanalysis/observation maps below):

| Variable | Meaning | Default |
|---|---|---|
| `TC_RUN_DIR` | Which `TC-TRACKS/<dir>/` to read/write | `NEW_WCSI_RUN` |
| `TC_DOMAINS` | Comma-separated basin names to include, e.g. `SIO,SPAC,SATL` | all |
| `TC_YEARS` | `first-last`, e.g. `1980-1989` | 1980-2013 |

```
TC_RUN_DIR=METHOD6_ENSEMBLE_RUN python huracanpy_TD_ensemble.py        # track density
TC_RUN_DIR=METHOD6_ENSEMBLE_RUN python huracanpy_TC_genesis_maps.py    # genesis density + latitude
TC_RUN_DIR=METHOD6_ENSEMBLE_RUN python huracanpy_TC_interannual.py     # counts/ACE/genesis lat/LMI lat, interannual
```

- `huracanpy_TD_ensemble.py`: one multi-panel figure per basin (full track +
  per-nature-category density), ensemble mean in colour, ensemble spread
  (std dev across experiments) as smoothed black contours.
- `huracanpy_TC_genesis_maps.py`: genesis density maps (IBTrACS/ERA5/model
  ensemble on one colour scale) plus genesis-latitude and
  latitude-of-first-34kt distributions, per basin.
- `huracanpy_TC_interannual.py`: one time-series plot per basin per metric
  (storm count, ACE, genesis latitude, latitude of maximum intensity split
  by whether it happens during the TC phase) - ensemble mean line, ±1 std
  dev shaded band, IBTrACS/ERA5/JRA3Q for comparison, trend and correlation
  statistics in the legend and a `_stats.csv`.

**Track density for reanalyses and observations**, directly comparable with
the ensemble's `density_ENSEMBLE_<years>_<domain>_MultiNature.png`:

```
TC_RUN_DIR=METHOD6_ENSEMBLE_RUN python huracanpy_TD_reanalysis.py   # ERA5 and JRA3Q (or name them: ... ERA5)
python huracanpy_ibtracs_TD.py                                      # IBTrACS
```

- `huracanpy_TD_reanalysis.py` uses the ensemble script's own plotting
  (`plot_domain_ensemble(..., name="ERA5")`), so grid, colour scale, domains
  and the method-6 nature panels are the same; it writes
  `density_{ERA5,JRA3Q}_<years>_<domain>_MultiNature.png` into `TC_RUN_DIR`.
  There are no spread contours (one dataset).
- `huracanpy_ibtracs_TD.py` writes `density_IBTrACS_<years>_<domain>_MultiNature.png`
  into the current directory. WCSI natures can't be computed for observations,
  so its small panels show IBTrACS's own natures (TS, ET, SS, DS, MX); its
  years come from `YEARS` in `huracanpy_ibtracs_batch_all.py` (1980-2013).
- Each figure scales its colours to its own maximum - compare the colour bars,
  not just the colours.

All three read every experiment named in their own `EXPERIMENTS` list (kept
in sync across the three files) from `TC_RUN_DIR` - add a newly-identified
experiment there once its method-6 output exists in that shared directory.
An experiment with no files for a given domain is silently excluded from
that domain's ensemble rather than erroring, so a partially-tracked
experiment degrades gracefully rather than blocking the others.

## Adding a brand new simulation - checklist

1. Add a `SOURCE_GROUP` entry (1.1) describing where its raw `.pp` output
   lives and which stream(s) hold which fields - verify with `iris`
   directly, don't assume it matches an existing experiment.
2. If its geopotential levels aren't 900/600/300 hPa, add `_EERIE`-style
   regZ/hart templates (1.3) with the right levels.
3. Add it to `tc_workflow.py`'s `EXPERIMENTS` dict (1.4) and run `status`,
   then `submit` (dry-run first) to track it on JASMIN.
4. `rsync` the tracked output back (Part 2).
5. Sanity-check it alone (e.g. method-6 output with `--maps` for a few
   years) before trusting it in the ensemble.
6. Run `huracanpy_identify_method6.py --experiments u-<experiment>` to
   produce its method-6 output in the shared ensemble directory (optionally
   `huracanpy_method6_maps_csv.py` for maps and CSVs, per year or `--combine`d).
7. Add it to the `EXPERIMENTS` list in all three Part 4 scripts and
   re-run them.

## Requirements

- JASMIN account with access to TRACK-1.5.2, the relevant group workspace(s)
  for raw model output, and SLURM.
- On JASMIN: `module add jaspy` (provides `iris`, `cdo`); TRACK-1.5.2's own
  compiled binaries (`master`, `config.RUN`, `utils/bin/{hart,tcident,...}`).
- For identification and analysis (laptop or JASMIN): Python with
  `huracanpy`, `xarray`, `pandas`, `pyarrow`, `parse`, `scipy`, `cartopy`, and
  Leo Saffin's [`wcsi`](https://github.com/Huracan-project/wcsi) **version 2 or
  later**: `pip install --upgrade git+https://github.com/Huracan-project/wcsi`.
  `huracanpy_identify_method6.py` stops with an explicit message if the
  installed wcsi is too old.
- Default paths at the top of the scripts (`DATA_DIR`, `OUT_DIR`, `BASE_DIR`,
  `IBTRACS_FILE`) point at PLV's machine - change them or use the
  command-line options.
