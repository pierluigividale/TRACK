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
python tc_workflow.py submit [--experiments ...] [--years ...] [--hemisphere NH|SH] [--yes]
```

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
(TRACK's chunking can't handle a partial year) will fail silently rather
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
are tropical cyclones. Several identification approaches exist (see
`huracanpy_compare_4_methods.py`'s docstring for all six) - the one
currently recommended (suggested by Kevin Hodges and Leo Saffin) is
**"WCSI warm core only"**: a closed-MSLP WCSI criterion with no minimum
850 hPa vorticity threshold, selecting purely on the warm-core structure
test (850-250 hPa vorticity decaying with height). This is "method 6" in
`huracanpy_compare_4_methods.py` and the default settings in
`huracanpy_identify_batch_method6.py`.

Two tools, for two different purposes:

- **`huracanpy_compare_4_methods.py --exp=u-<experiment> [first_year
  last_year] [--sh] [--run=6]`** - runs one or more approaches on *one*
  experiment and compares them (storm counts, genesis-latitude
  distribution, overlap between approaches). Useful for sanity-checking a
  newly-tracked experiment, or for deciding which approach to use in the
  first place. Writes per-approach output under
  `TC-TRACKS/TEST_<experiment>_<NH|SH>H/<approach>/` - a layout for
  *comparing approaches on one experiment*, not for the ensemble scripts
  below.
- **`huracanpy_identify_batch_method6.py`** - bulk production, every
  experiment in its `EXPERIMENTS` list, every year, both hemispheres,
  method 6 only, writing every experiment's output into *one shared
  directory* (`TC-TRACKS/METHOD6_ENSEMBLE_RUN/` by default) - this is the
  layout Part 4's scripts expect: `{experiment}_{year}_{domain}_TCs.nc` /
  `_summary.parquet` files for every experiment side by side. Add a new
  experiment to `EXPERIMENTS` once it's tracked and `rsync`'d; the script
  skips anything it's already processed.

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
pointed at the same run and the same subset of basins/years:

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
5. Sanity-check it alone with `huracanpy_compare_4_methods.py --exp=...`
   (Part 3) before trusting it in the ensemble.
6. Add it to `huracanpy_identify_batch_method6.py`'s `EXPERIMENTS` and run
   it to produce method-6 output in the shared ensemble directory.
7. Add it to the `EXPERIMENTS` list in all three Part 4 scripts and
   re-run them.

## Requirements

- JASMIN account with access to TRACK-1.5.2, the relevant group workspace(s)
  for raw model output, and SLURM.
- On JASMIN: `module add jaspy` (provides `iris`, `cdo`); TRACK-1.5.2's own
  compiled binaries (`master`, `config.RUN`, `utils/bin/{hart,tcident,...}`).
- Locally: Python with `huracanpy`, `iris`, `cartopy`, `xarray`, `pandas`,
  `scipy`; Leo Saffin's [`wcsi`](https://wcsi.readthedocs.io/) package (used
  by `huracanpy_identify_batch_wcsi_new.py` for the closed-MSLP WCSI
  identification, referenced via a local path near the top of that file -
  point it at your own clone).
