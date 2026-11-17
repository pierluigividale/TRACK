"""Climatological track-density maps (MultiNature) for the reanalyses,
directly comparable with density_ENSEMBLE_*_MultiNature.png: same code
(huracanpy_TD_ensemble.plot_domain_ensemble), bins, colour scale, domains and
method-6 nature panels, applied to each reanalysis as a single dataset.

Reads <BASE_DIR>/<TC_RUN_DIR>/{ERA5,JRA3Q}_<year>_<domain>_TCs.nc (written by
huracanpy_identify_method6.py) and writes density_<name>_<y0>-<y1>_<domain>_MultiNature.png
in the same directory.

    TC_RUN_DIR=METHOD6_ENSEMBLE_RUN python huracanpy_TD_reanalysis.py [ERA5 JRA3Q]

For IBTrACS (observed natures, not WCSI) use huracanpy_ibtracs_TD.py.
"""
import os
import sys

os.environ.setdefault("TC_RUN_DIR", "METHOD6_ENSEMBLE_RUN")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import huracanpy_TD_ensemble as td

REANALYSES = sys.argv[1:] or ["ERA5", "JRA3Q"]

for name in REANALYSES:
    for domain_name, domain in td.DOMAINS.items():
        if domain_name in td.SKIP_DOMAINS:
            continue
        td.plot_domain_ensemble(domain_name, domain, [name], td.YEARS, name=name)
        plt.close("all")
        print(f"{name} {domain_name} done", flush=True)
