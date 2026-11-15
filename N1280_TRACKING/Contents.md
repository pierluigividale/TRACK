# N1280_TRACKING - scripts to track and identify tropical cyclones

See the top-level README.md for how these fit together. Each script's own
docstring has the details.

## On JASMIN: conversion and tracking (TRACK-1.5.2)

| File | Purpose |
|---|---|
| `tc_workflow.py` | Audit / stage / submit tool for the whole JASMIN pipeline: `status`, `stage` (EERIE runs, NLDS tape retrieval, sci node only), `submit`. Experiments and year ranges are in its `EXPERIMENTS` table. |
| `eerie_nlds.py` | Shared helpers for the EERIE runs (raw-file paths, staging area, NLDS calls, latest month delivered from ARCHER2). |
| `pp_to_nc_iris_WCSI.py` | Converts one month of raw UM pp output to the NetCDF fields TRACK needs (standard atmosphere-only runs). |
| `pp_to_nc_iris_WCSI_EERIE.py` | Same for the coupled EERIE runs (dz876, ec100); reads raw files from disk or the NLDS staging area. |
| `ll_to_gg_n1280.py` | Resamples to a 6-hourly cycle and regrids winds to TRACK's Gaussian grid. |
| `track_n1280_para_WCSI.py` | Runs TRACK for one experiment/year/hemisphere: vorticity filtering, tracking, Hart phase space, tcident. |
| `SLURMscript_n1280_para_WCSI` | SLURM job for one month's conversion (submitted by `tc_workflow.py`). |
| `SLURMscript_n1280_track_para_WCSI` | SLURM job for one year's tracking (submitted by `tc_workflow.py`). |

## TC identification (laptop or JASMIN)

| File | Purpose |
|---|---|
| `huracanpy_identify_method6.py` | **The TC identification**: method 6 ("WCSI warm core only"), self-contained, command-line options for any TRACK output layout. |
| `huracanpy_method6_maps_csv.py` | Per year and domain: PNG map and CSV of the identified TC tracks (Hart parameters, WCSI nature). |
| `huracanpy_identify_batch_all.py` | Earlier identification script, kept as a library of shared definitions imported by the scripts below. Do not run it for identification (needs wcsi 1.x). |

## Analysis and comparison with observations / reanalyses

| File | Purpose |
|---|---|
| `huracanpy_TD_ensemble.py` | Ensemble track-density maps per basin. |
| `huracanpy_TC_genesis_maps.py` | Genesis density maps and genesis-latitude distributions. |
| `huracanpy_TC_interannual.py` | Interannual time series per basin (counts, ACE, genesis and LMI latitude). |
| `huracanpy_era5_batch_all.py`, `huracanpy_era5_TD.py` | TC identification and track density for ERA5 reanalysis tracks. |
| `huracanpy_ibtracs_batch_all.py`, `huracanpy_ibtracs_TD.py` | IBTrACS observations, same outputs. |
| `huracanpy_TD.py`, `huracanpy_TD.ipynb` | Single-experiment track density (notebook origin). |
