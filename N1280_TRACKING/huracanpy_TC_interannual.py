"""Interannual variability of tropical cyclone activity, per basin.

For each domain, one time-series plot per metric (see METRICS) for IBTrACS,
ERA5 and the model ensemble (ensemble mean as a line, +/- 1 std dev across
members as a shaded band):

  - counts:        number of TCs per year
  - ace:           accumulated cyclone energy per year
  - genesis_lat:   mean latitude of TC genesis
  - lmi_tc_lat:    mean LMI (latitude of maximum intensity), over the storms
                   whose maximum intensity happens while they are a TC
  - lmi_other_lat: the same, over the storms whose maximum intensity happens
                   as anything else (e.g. BC after extratropical transition)

The legend gives each series' linear trend (least squares, per decade, *
where p < 0.05) and the Pearson correlation of ERA5 and of the ensemble mean
with IBTrACS, both as is and after removing each series' own linear trend -
so the detrended r measures agreement on year-to-year variability alone, not
on a shared long-term drift. The same statistics are written to a _stats.csv
next to the plots, and the yearly values to a .csv.

Which storms count depends on WIND_THRESHOLD_KT:

  - 34 (default): named storms, i.e. the officially counted systems of each
    season - a track that reaches tropical-storm strength (>= 34 kt) while it
    is tropical. For IBTrACS that is the maximum reported wind over its
    nature == "TS" points; for the model experiments the maximum 10 m wind
    speed over the points flagged is_tc in the {experiment}_{year}_{domain}
    _TCs.nc files. The WCSI is_tc criteria themselves have no wind test, and
    IBTrACS's TS nature includes tropical depressions, so without this both
    sides count weaker systems than the published season totals do.
    The reanalyses (ERA5, JRA3Q) are NOT thresholded unless
    REANALYSIS_USE_THRESHOLD is set: their 10 m winds are far too weak for an absolute threshold to mean the same thing
    (see that setting), so for ERA5 every is_tc track counts.
  - None: every tropical system, no wind test - is_tc tracks for the models
    and ERA5 (for counts, taken from the per-year _summary.parquet files, the
    same count the batch scripts print in their "(N TCs)" map titles), and
    tracks with at least one nature == "TS" point for IBTrACS (the is_ts
    definition of huracanpy_ibtracs_batch_all.py).

Metric definitions, for the storms that count:

  - "Tropical" points are is_tc points for the models/ERA5 and nature == "TS"
    points for IBTrACS.
  - ACE: 1e-4 * sum of (max wind in kt)**2 over 6-hourly tropical points with
    max wind >= 34 kt, in 10**4 kt**2 - the standard definition, and the same
    34 kt cut for every dataset regardless of WIND_THRESHOLD_KT (it is part of
    what ACE is). IBTrACS is 3-hourly, so only its 00/06/12/18 UTC points are
    used; the model and ERA5 tracks are 6-hourly already. ERA5 ACE is low
    because its winds are (see above).
  - Genesis latitude: latitude of the first tropical point.
  - LMI: the latitude of the point of maximum wind (not minimum pressure)
    over the storm's whole track (not just its tropical part), split by the
    nature at that point - "TC" for the models
    and ERA5, "TS" for IBTrACS - versus anything else. Note this uses
    `nature`, not `is_tc`: in the model files every nature == "TC" point is
    is_tc but not the other way round. IBTrACS points with nature "NR" (not
    reported) are ignored, and IBTrACS reports wind for only a minority of
    its non-TS points, so its "other" series is thin.
  - Latitudes are absolute values (degrees poleward), so trends have the same
    sign convention in both hemispheres.
  - The model "wind" is the tracker's 10 m wind speed; IBTrACS's is each
    agency's own sustained wind (1-min for the US agencies, 10-min for most
    others), so ACE magnitudes are not strictly comparable across datasets.

IBTrACS is computed directly from its NetCDF rather than from
huracanpy_ibtracs_batch_all.py's CSV, because of the season shift below.

Southern-hemisphere seasons: the model/ERA5 files label a July-June season by
the year it STARTS ("1995" = Jul 1995-Jun 1996), while IBTrACS's `season`
labels it by the year it ENDS (1996). IBTrACS seasons are therefore shifted
back by one year for southern-hemisphere domains so all three series line up
on the same seasons - without this the correlations would be between
different seasons.
"""

import os

import huracanpy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator
from scipy import stats

from huracanpy_identify_batch_all import BASE_DIR, DOMAINS
from huracanpy_ibtracs_batch_all import DOMAINS as IBTRACS_DOMAINS
from huracanpy_ibtracs_batch_all import IBTRACS_FILE, load_ibtracs

# Which run's files to read and write (a folder under TC-TRACKS), e.g.
# TC_RUN_DIR=TCIDENT_WC_RUN python huracanpy_TC_interannual.py
BASE_DIR = os.path.join(BASE_DIR, os.environ.get("TC_RUN_DIR", "NEW_WCSI_RUN"))

EXPERIMENTS = ["u-ch330", "u-ck777", "u-ck778", "u-dc009"]  # keep in sync with huracanpy_TD_ensemble.py
# u-dd828 excluded: shorter radiation timestep experiment, confirmed outlier
# (r vs IBTrACS 0.24 vs 0.60-0.62 for the others in NATL counts) - excluded
# until the pressure-wind intensity check lands
# Reanalyses (name -> plot colour), each with {name}_{year}_{domain}_TCs.nc /
# _summary.parquet files from huracanpy_era5_batch_all.py (--name JRA3Q for
# JRA3Q). See huracanpy_TC_genesis_maps.py's docstring for how the two
# reanalyses' input track sets differ.
REANALYSES = {"ERA5": "#eb6834", "JRA3Q": "#1baf7a"}
YEARS = range(1980, 2015)  # NH now includes 2014; SH stays capped at 2013 (harmless skip)
# Optional override, e.g. TC_YEARS=1980-1989, for runs covering only part of the period
# (same convention as huracanpy_TD_ensemble.py).
# (per-year/per-domain lookups below already skip cleanly when a year has no data).
if os.environ.get("TC_YEARS"):
    _first, _last = (int(x) for x in os.environ["TC_YEARS"].split("-"))
    YEARS = range(_first, _last + 1)
SKIP_DOMAINS = {"MED"}  # no IBTrACS counterpart - see huracanpy_ibtracs_batch_all.py
# Optional: only these domains, e.g. TC_DOMAINS=SIO,SPAC,SATL for a Southern Hemisphere run folder
if os.environ.get("TC_DOMAINS"):
    SKIP_DOMAINS = set(SKIP_DOMAINS) | (set(DOMAINS) - set(os.environ["TC_DOMAINS"].split(",")))

# Minimum lifetime-maximum wind (knots, reached while tropical) for a track to
# be counted; 34 = named storms, 64 = hurricanes/typhoons, None = no wind test.
# See the module docstring.
WIND_THRESHOLD_KT = 34
# Reanalysis 10 m winds are much weaker than observed (only about half of
# ERA5's is_tc tracks even reach 34 kt in the East Pacific and North Indian
# Ocean), so by default the threshold is not applied to them.
REANALYSIS_USE_THRESHOLD = False
KT_TO_MS = 0.514444
WIND_VARIABLES = {name: "vmax10m" for name in REANALYSES}  # 10 m wind (m/s) in each _TCs.nc
DEFAULT_WIND_VARIABLE = "wind_speed"
ACE_THRESHOLD_KT = 34
TC_NATURE = "TC"  # models/ERA5 nature code for the LMI split
IBTRACS_TC_NATURE = "TS"

# Only the run folder itself. (Earlier this also fell back to the scripts'
# directory, where some experiments' older files still are - with a run folder
# that lacks a year on purpose, e.g. a merged domain whose parts aren't all
# processed, that would silently mix in results from another run.)
SUMMARY_DIRS = [BASE_DIR]

MIN_YEARS_FOR_CORRELATION = 10  # also the minimum for fitting a trend
TREND_SIGNIFICANCE = 0.05
MIN_MEMBER_YEARS = len(YEARS) // 2  # an experiment with fewer processed years than
                                    # this for a domain is left out of its ensemble

# One plot per metric. "latitude" metrics are means over the storms that have
# a value that year; the others are yearly totals (0 in a year with no storms).
METRICS = {
    "counts": dict(title="TC count", ylabel="Tropical cyclones per year", unit="", latitude=False),
    "ace": dict(title="ACE", ylabel="ACE (10$^4$ kt$^2$)", unit=" ×10⁴ kt²", latitude=False),
    "genesis_lat": dict(title="mean genesis latitude", ylabel="Mean genesis latitude", unit="°", latitude=True),
    "lmi_tc_lat": dict(
        title="mean LMI (latitude of maximum intensity), storms peaking as TC",
        ylabel="Mean LMI, peak as TC", unit="°", latitude=True,
    ),
    "lmi_other_lat": dict(
        title="mean LMI (latitude of maximum intensity), storms peaking as non-TC",
        ylabel="Mean LMI, peak as non-TC", unit="°", latitude=True,
    ),
}
LATITUDE_METRICS = [name for name, metric in METRICS.items() if metric["latitude"]]

COLOR_IBTRACS = "#1a1a19"
COLOR_REANALYSIS = REANALYSES["ERA5"]
COLOR_ENSEMBLE = "#2a78d6"


def find_output(filename):
    for summary_dir in SUMMARY_DIRS:
        path = os.path.join(summary_dir, filename)
        if os.path.exists(path):
            return path
    return None


# ---------------------------------------------------------------------------
# Per-storm tables: one row per track, with whether it counts and its
# ACE/genesis/LMI values - every metric is then a yearly aggregate of these.
# ---------------------------------------------------------------------------
def storm_table(points, threshold_kt):
    """Per-track metrics from a per-point DataFrame with columns track_id,
    time, lat, wind_kt, tropical (bool), is_tc_nature (bool), use_for_ace
    (bool), use_for_lmi (bool).
    """
    points = points.sort_values(["track_id", "time"])
    points = points.assign(abs_lat=points.lat.abs())
    tropical = points[points.tropical]

    reaches = tropical.wind_kt >= threshold_kt if threshold_kt is not None else pd.Series(True, tropical.index)
    counted = reaches.groupby(tropical.track_id).any()

    ace_points = points.use_for_ace & points.tropical & (points.wind_kt >= ACE_THRESHOLD_KT)
    ace = (1e-4 * points.wind_kt.where(ace_points) ** 2).groupby(points.track_id).sum()

    genesis_lat = tropical.groupby("track_id").abs_lat.first()

    lmi_candidates = points[points.use_for_lmi & points.wind_kt.notna()]
    lmi = lmi_candidates.loc[lmi_candidates.groupby("track_id").wind_kt.idxmax()].set_index("track_id")
    lmi_tc_lat = lmi.abs_lat.where(lmi.is_tc_nature)
    lmi_other_lat = lmi.abs_lat.where(~lmi.is_tc_nature)

    table = pd.DataFrame(dict(
        counted=counted, ace=ace, genesis_lat=genesis_lat, lmi_tc_lat=lmi_tc_lat, lmi_other_lat=lmi_other_lat,
    ))
    table["counted"] = table.counted.fillna(False).astype(bool)  # tracks with no tropical point
    return table


def model_storms(experiment, domain_name, years, threshold_kt):
    """(storms, summary_counts): per-track table across all years (with a
    `year` column) from the _TCs.nc files, and the yearly is_tc count from the
    _summary.parquet files. summary_counts is NaN for years with no summary
    (not run through the batch pipeline), as opposed to 0 for a processed
    year with no TCs; storms has no rows for either.
    """
    wind_variable = WIND_VARIABLES.get(experiment, DEFAULT_WIND_VARIABLE)
    summary_counts = pd.Series(np.nan, index=pd.Index(years, name="year"), dtype=float)
    tables = []
    for year in years:
        tag = f"{experiment}_{year}_{domain_name}"
        summary = find_output(f"{tag}_summary.parquet")
        if summary is None:
            continue
        summary_counts[year] = pd.read_parquet(summary, columns=["is_tc"]).is_tc.sum()

        # The batch scripts only write a _TCs.nc when the year has TC tracks.
        # Look for it only next to its own summary - searching every
        # SUMMARY_DIRS entry could pair this run's zero-TC year with another
        # run's _TCs.nc.
        tcs_nc = os.path.join(os.path.dirname(summary), f"{tag}_TCs.nc")
        if not os.path.exists(tcs_nc):
            continue
        tracks = huracanpy.load(tcs_nc)
        wind = tracks[wind_variable].where(tracks[wind_variable] < 1e10).values  # TRACK's 1e25 fill value
        points = pd.DataFrame(dict(
            track_id=tracks.track_id.values,
            time=pd.to_datetime(tracks.time.values),
            lat=tracks.lat.values,
            wind_kt=wind / KT_TO_MS,
            tropical=tracks.is_tc.values,
            is_tc_nature=tracks.nature.values == TC_NATURE,
            use_for_ace=True,  # already 6-hourly
            use_for_lmi=True,
        ))
        tables.append(storm_table(points, threshold_kt).assign(year=year))

    storms = pd.concat(tables) if tables else pd.DataFrame(columns=["year", "counted", *METRICS])
    return storms, summary_counts


def ibtracs_storms(all_tracks, domain_name, threshold_kt):
    """Per-track table for one IBTrACS domain, with `year` on the model's
    season-labelling convention (see module docstring). None if IBTrACS has
    no such domain.
    """
    if domain_name not in IBTRACS_DOMAINS:
        return None
    sel = all_tracks[all_tracks.domain_basin == IBTRACS_DOMAINS[domain_name]["basin"]]
    time = pd.DatetimeIndex(sel.time)
    points = pd.DataFrame(dict(
        track_id=sel.track_id.values,
        time=time,
        lat=sel.lat.values,
        wind_kt=sel.wind.values,
        tropical=(sel.nature == IBTRACS_TC_NATURE).values,
        is_tc_nature=(sel.nature == IBTRACS_TC_NATURE).values,
        use_for_ace=((time.hour % 6 == 0) & (time.minute == 0)),
        use_for_lmi=(sel.nature != "NR").values,
    ))
    storms = storm_table(points, threshold_kt)
    season = sel.groupby("track_id").season.first()
    storms["year"] = season - 1 if DOMAINS[domain_name]["hemisphere"] == "S" else season
    return storms


def yearly(storms, metric, years, processed_years=None):
    """Yearly series of one metric over the storms that count. Years outside
    processed_years (None = all) are NaN.
    """
    counted = storms.loc[storms.counted.astype(bool)]  # astype: an empty table has an untyped column
    if metric == "counts":
        series = counted.groupby("year").size()
    elif metric == "ace":
        series = counted.groupby("year").ace.sum()
    else:
        series = counted.groupby("year")[metric].mean()
    index = pd.Index(years, name="year")
    # a processed year with no storms is a real 0 for totals, but has no mean latitude
    series = series.reindex(index, fill_value=np.nan if METRICS[metric]["latitude"] else 0).astype(float)
    if processed_years is not None:
        series[~index.isin(processed_years)] = np.nan
    return series


def model_yearly(storms, summary_counts, metric, years, threshold_kt):
    if metric == "counts" and threshold_kt is None:
        # includes the tracks the batch scripts' genesis-latitude filter drops
        # from the _TCs.nc, so it matches their printed "(N TCs)" counts
        return summary_counts
    return yearly(storms, metric, years, processed_years=summary_counts.dropna().index)


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------
def linear_trend(series):
    """Least-squares fit over the years the series has: (slope per year,
    p-value), or None if too few years or a constant series.
    """
    valid = series.dropna()
    if len(valid) < MIN_YEARS_FOR_CORRELATION or valid.std() == 0:
        return None
    fit = stats.linregress(valid.index, valid.values)
    return fit.slope, fit.pvalue


def detrend(series):
    """Residuals from the series' own linear fit (NaN years stay NaN)."""
    valid = series.dropna()
    if len(valid) < MIN_YEARS_FOR_CORRELATION:
        return series * np.nan
    fit = stats.linregress(valid.index, valid.values)
    return series - (fit.intercept + fit.slope * series.index.values)


def correlation(series, reference):
    """Pearson r over the years both series have, or None if too few."""
    both = pd.concat([series, reference], axis=1).dropna()
    if len(both) < MIN_YEARS_FOR_CORRELATION or (both.std() == 0).any():
        return None
    return float(np.corrcoef(both.iloc[:, 0], both.iloc[:, 1])[0, 1])


def series_stats(series, reference):
    """Trend of `series`, and its correlation with `reference` before and
    after detrending both (None where not computable / no reference).
    """
    trend = linear_trend(series)
    r = r_detrended = None
    if reference is not None:
        r = correlation(series, reference)
        r_detrended = correlation(detrend(series), detrend(reference))
    return dict(
        mean=float(series.mean()) if series.notna().any() else None,
        trend_per_decade=None if trend is None else 10 * trend[0],
        trend_p=None if trend is None else trend[1],
        r=r, r_detrended=r_detrended,
    )


def legend_label(name, st, unit):
    """Series name, then a second line with its trend and correlations."""
    parts = []
    if st["trend_per_decade"] is not None:
        star = "*" if st["trend_p"] < TREND_SIGNIFICANCE else ""
        decimals = 2 if unit == "°" else 1
        parts.append(f"trend {st['trend_per_decade']:+.{decimals}f}{unit}/decade{star}")
    if st["r"] is not None:
        r_text = f"r = {st['r']:.2f}"
        if st["r_detrended"] is not None:
            r_text += f" ({st['r_detrended']:.2f} detrended)"
        parts.append(r_text)
    return name + ("\n" + ", ".join(parts) if parts else "")


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def output_stem(years, metric="counts"):
    """Common path prefix of the outputs - the wind threshold is part of the
    name so runs with different thresholds don't overwrite each other.
    """
    suffix = "" if WIND_THRESHOLD_KT is None else f"_min{WIND_THRESHOLD_KT}kt"
    return f"{BASE_DIR}/TC_{metric}_interannual{suffix}_{years.start}-{years.stop - 1}"


def plot_metric(domain_name, domain, metric, obs, reanalyses, members, years):
    """One metric's time-series plot for one domain; returns its stats rows.
    reanalyses: {name: yearly series}.
    """
    spec = METRICS[metric]
    threshold = WIND_THRESHOLD_KT
    reanalysis_threshold = threshold if REANALYSIS_USE_THRESHOLD else None
    wind_label = "" if threshold is None else f", ≥ {threshold} kt"
    reanalysis_suffix = ""
    if reanalysis_threshold != threshold and metric != "ace":  # ACE has its own fixed 34 kt cut
        reanalysis_suffix = " (no wind threshold)"

    ens_mean = members.mean(axis=1)
    ens_spread = members.std(axis=1, ddof=1)  # NaN where only one member has that year

    st_obs = series_stats(obs, None) if obs is not None else None
    st_reanalyses = {name: series_stats(series, obs) for name, series in reanalyses.items()}
    st_ensemble = series_stats(ens_mean, obs)

    fig, ax = plt.subplots(figsize=(10, 5.2), layout="constrained")

    if members.shape[1] > 1:
        ax.fill_between(
            ens_mean.index, ens_mean - ens_spread, ens_mean + ens_spread,
            color=COLOR_ENSEMBLE, alpha=0.2, linewidth=0, label="Ensemble mean ± 1 std dev",
        )
    if members.shape[1] > 0:
        ax.plot(
            ens_mean.index, ens_mean, color=COLOR_ENSEMBLE, linewidth=2,
            label=legend_label(f"Ensemble mean ({members.shape[1]} members)", st_ensemble, spec["unit"]),
        )
    for name, series in reanalyses.items():
        if series.notna().any():
            ax.plot(
                series.index, series, color=REANALYSES[name], linewidth=2,
                label=legend_label(name + reanalysis_suffix, st_reanalyses[name], spec["unit"]),
            )
    if obs is not None:
        ax.plot(
            obs.index, obs, color=COLOR_IBTRACS, linewidth=2, marker="o", markersize=3.5,
            label=legend_label("IBTrACS", st_obs, spec["unit"]),
        )

    ax.set_xlim(years.start, years.stop - 1)
    ax.set_xlabel("Season (year it starts)" if domain["hemisphere"] == "S" else "Year")
    if spec["latitude"]:
        ax.set_ylabel(f"{spec['ylabel']} (°{domain['hemisphere']})")
    else:
        ax.set_ylim(bottom=0)
        ax.set_ylabel(f"{spec['ylabel']}{wind_label if metric == 'counts' else ''}")
    if metric == "counts":
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(color="0.85", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    legend_title = f"Linear trend (* p < {TREND_SIGNIFICANCE})"
    if obs is not None:
        legend_title += "\nr: correlation with IBTrACS"
    ax.legend(
        title=legend_title, loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, alignment="left",
        labelspacing=1.0,
    )
    # suptitle rather than ax.set_title: centred on the whole figure, so the
    # longer metric titles aren't clipped by the legend column
    fig.suptitle(
        f"{domain['title']}: {spec['title']}{wind_label}, {years.start}-{years.stop - 1}\n"
        f"Ensemble: {', '.join(members.columns)}",
        fontsize=10,
    )

    fig.savefig(f"{output_stem(years, metric)}_{domain_name}.png", dpi=300)
    plt.close(fig)

    fmt = lambda x: "n/a" if x is None else f"{x:.2f}"
    print(
        f"[{domain_name} {metric}] r vs IBTrACS (detrended): "
        + "".join(f"{name} {fmt(st['r'])} ({fmt(st['r_detrended'])}), " for name, st in st_reanalyses.items())
        + f"ensemble mean {fmt(st_ensemble['r'])} ({fmt(st_ensemble['r_detrended'])})"
    )
    rows = [
        dict(domain=domain_name, metric=metric, dataset="ensemble_mean", **st_ensemble),
        *[dict(domain=domain_name, metric=metric, dataset=name, **st) for name, st in st_reanalyses.items()],
    ]
    if st_obs is not None:
        rows.insert(0, dict(domain=domain_name, metric=metric, dataset="IBTrACS", **st_obs))
    return rows, ens_mean, ens_spread


def process_domain(domain_name, domain, all_tracks, experiments, years):
    threshold = WIND_THRESHOLD_KT
    reanalysis_threshold = threshold if REANALYSIS_USE_THRESHOLD else None

    obs_storms = ibtracs_storms(all_tracks, domain_name, threshold)
    reanalysis_storms = {name: model_storms(name, domain_name, years, reanalysis_threshold) for name in REANALYSES}
    member_storms = {exp: model_storms(exp, domain_name, years, threshold) for exp in experiments}

    n_years = pd.Series({exp: summary.notna().sum() for exp, (_, summary) in member_storms.items()})
    for experiment in n_years.index[n_years < MIN_MEMBER_YEARS]:
        print(
            f"[{domain_name}] {experiment}: only {n_years[experiment]} of {len(years)} years have a "
            f"_summary.parquet, excluding from ensemble"
        )
        del member_storms[experiment]
    n_per_year = pd.DataFrame({exp: summary.notna() for exp, (_, summary) in member_storms.items()}).sum(axis=1)
    partial = n_per_year[n_per_year < len(member_storms)]
    if len(partial):
        print(
            f"[{domain_name}] ensemble mean uses fewer than {len(member_storms)} members in "
            f"{', '.join(f'{y} ({n})' for y, n in partial.items())}"
        )

    stats_rows, tables = [], []
    for metric in METRICS:
        obs = yearly(obs_storms, metric, years) if obs_storms is not None else None
        reanalyses = {
            name: model_yearly(*storms, metric, years, reanalysis_threshold)
            for name, storms in reanalysis_storms.items()
        }
        members = pd.DataFrame({
            exp: model_yearly(storms, summary, metric, years, threshold)
            for exp, (storms, summary) in member_storms.items()
        })
        rows, ens_mean, ens_spread = plot_metric(domain_name, domain, metric, obs, reanalyses, members, years)
        stats_rows.extend(rows)

        table = members.copy()
        for name, series in reversed(list(reanalyses.items())):
            table.insert(0, name, series)
        if obs is not None:
            table.insert(0, "IBTrACS", obs)
        table["ensemble_mean"] = ens_mean
        table["ensemble_std"] = ens_spread
        table.insert(0, "metric", metric)
        table.insert(0, "domain", domain_name)
        tables.append(table.reset_index())
    return tables, stats_rows


def main():
    print(f"Loading {IBTRACS_FILE} ...")
    all_tracks = load_ibtracs(IBTRACS_FILE)

    tables, stats_rows = [], []
    for domain_name, domain in DOMAINS.items():
        if domain_name in SKIP_DOMAINS:
            continue
        domain_tables, domain_stats = process_domain(domain_name, domain, all_tracks, EXPERIMENTS, YEARS)
        tables.extend(domain_tables)
        stats_rows.extend(domain_stats)

    if tables:
        stem = output_stem(YEARS, "metrics")
        pd.concat(tables, ignore_index=True).to_csv(f"{stem}.csv", index=False)
        pd.DataFrame(stats_rows).to_csv(f"{stem}_stats.csv", index=False, float_format="%.4g")


if __name__ == "__main__":
    main()
