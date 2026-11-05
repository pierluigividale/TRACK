"""Bulk method-6 ("WCSI warm core only") identification, all ensemble members.

Reconfigures huracanpy_identify_batch_wcsi_new.py's process_file() with
method 6's settings (see huracanpy_compare_4_methods.py's run_wcsi() for the
6 WCSI warm core only approach: VORT_THRESHOLD=0, VORT_WARM_CORE_THRESHOLD=6e-5)
and runs it for every experiment/year/hemisphere into ONE shared output
directory, so huracanpy_TD_ensemble.py / huracanpy_TC_interannual.py /
huracanpy_TC_genesis_maps.py (which all expect every experiment's
{experiment}_{year}_{domain}_TCs.nc side by side in a single TC_RUN_DIR) see
a real multi-member ensemble - not huracanpy_compare_4_methods.py's own
per-experiment TEST_<EXP>_<NH|SH>H/6_WCSI_warm_core_only/ layout, which is
for comparing approaches on one experiment, not for ensemble consumption.

u-ch330 already has method-6 output (from huracanpy_compare_4_methods.py) -
copy its existing _TCs.nc/_summary.parquet into OUT_DIR once rather than
reprocessing it here; EXPERIMENTS below is only the other 4 members.

Usage: python huracanpy_identify_batch_method6.py
"""
import huracanpy_identify_batch_wcsi_new as w

w.OUT_DIR = "/Users/vidale/Python_scripts/TC-TRACKS/METHOD6_ENSEMBLE_RUN"
w.EXPERIMENTS = ["u-ck777", "u-ck778", "u-dd828", "u-dc009"]  # u-ch330 copied in separately, see docstring
w.YEARS = range(1980, 2014)
w.VORT_THRESHOLD = 0
w.VORT_WARM_CORE_THRESHOLD = 6e-5
w.MAKE_MAPS = False  # ensemble plots only - no per-year-per-domain nature maps

if __name__ == "__main__":
    w.main()
