"""Four TC identification approaches on u-ch330 1982 (NH), with the recomputed Hart parameters.

  1 tcident            tcident warm-core tracks as delivered (no latitude band)
  2 WCSI               WCSI on the full track set (closed MSLP + 850 hPa vorticity >= 6e-5)
  3 WCSI on tcident    the same WCSI criteria, applied to the tcident tracks only
  4 WCSI + warm core   approach 2 plus tcident's warm-core test (850 - 250 hPa vorticity > 6e-5)

Writes per-basin _TCs.nc files and nature maps to TC-TRACKS/TEST_1982/<approach>/,
and prints counts, overlaps and genesis statistics for storms reaching 34 kt
while tropical (is_tc for WCSI, the warm-core phase for tcident).
"""
import os, shutil, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, xarray as xr, huracanpy
import huracanpy_identify_batch_wcsi_new as w
import huracanpy_identify_batch_all as ib
import huracanpy_tcident_batch as tb

EXP, YEAR = "u-ch330", 1982
YD = f"{ib.BASE_DIR}/{EXP}/VOR_VERTAVG_T63filt_{YEAR}"
FULL = f"{YD}/tr_trs_pos.2day_addvor_addmslp_addwind.hart.new.nc"
TC_FULL = f"{YD}/tr_trs_pos.2day_addvor_addmslp_addwind.hart.tcident.new.nc"
TC_WC = f"{YD}/tr_trs_pos.2day_addvor_addmslp_addwind.hart.tcident.wc.new.nc"
OUT = f"{ib.BASE_DIR}/TEST_1982"
NH = {n: d for n, d in ib.DOMAINS.items() if d["hemisphere"] == "N"}
REPORT = ["NATL", "NWPAC", "EPAC", "NIO"]
KT = 0.514444
APPROACHES = ["1 tcident", "2 WCSI", "3 WCSI on tcident", "4 WCSI + warm core"]


def folder(label):
    path = os.path.join(OUT, label.replace(" ", "_").replace("+", "plus"))
    shutil.rmtree(path, ignore_errors=True)
    os.makedirs(path)
    return path


def run_wcsi(label, path, warm_core):
    w.OUT_DIR, w.MAKE_MAPS = folder(label), True
    w.VORT_THRESHOLD, w.VORT_WARM_CORE_THRESHOLD = 6e-5, (6e-5 if warm_core else None)
    w.process_file(EXP, YEAR, path, NH)
    return w.OUT_DIR


def run_tcident(label):
    out = folder(label)
    ti = tb.load_marked(TC_FULL, TC_WC)
    dom = tb.assign_domains(ti, NH)
    ib.NATURE_MARKERS.setdefault("NT", ".")  # draw tcident's non-warm-core points on the maps too
    for d, spec in NH.items():
        ids = dom.index[dom == d].values
        if not len(ids):
            continue
        tracks = ib.filter_genesis_latitude(ti.hrcn.sel_id(ids), ib.TC_GENESIS_THRESHOLD_LATITUDE)
        if not tracks.sizes["record"]:
            continue
        huracanpy.save(tracks, f"{out}/{EXP}_{YEAR}_{d}_TCs.nc")
        ib.plot_domain_nature(tracks, title=f"{spec['title']} tracks by nature, tcident {EXP} {YEAR}",
                              out_path=f"{out}/{EXP}_{YEAR}_{d}_nature.png", extent=spec["extent"])
    return out


# raw track index of every point of the full track set, to match storms across approaches
raw = huracanpy.load(FULL)
key = lambda t, la, lo: list(zip(np.asarray(t).astype("datetime64[m]").astype(np.int64), np.round(np.asarray(la, float) * 100).astype(int),
                                 np.round(np.asarray(lo, float) * 100).astype(int)))
key_to_raw = dict(zip(key(raw.time.values, raw.lat.values, raw.lon.values), raw.track_id.values))


def storms(out, domain):
    """{raw index: (genesis lat, max wind kt)} for storms reaching 34 kt while tropical."""
    p = f"{out}/{EXP}_{YEAR}_{domain}_TCs.nc"
    if not os.path.exists(p):
        return {}
    ds = huracanpy.load(p)
    lat, lon = ds.lat.values, ds.lon.values
    df = pd.DataFrame(dict(track=ds.track_id.values, raw=[key_to_raw.get(k, -1) for k in key(ds.time.values, lat, lon)],
                           lat=lat, w=ds.wind_speed.values / KT, tc=ds.is_tc.values.astype(bool)))
    out_ = {}
    for _, g in df.groupby("track", sort=False):
        tc = g[g.tc]
        if len(tc) and tc.w.max() >= 34:
            out_[int(g.raw.mode().iloc[0])] = (tc.lat.iloc[0], tc.w.max())
    return out_


import sys
if "--compare-only" in sys.argv:
    outs = {a: os.path.join(OUT, a.replace(" ", "_").replace("+", "plus")) for a in APPROACHES}
else:
  outs = {
    "1 tcident": run_tcident("1 tcident"),
    "2 WCSI": run_wcsi("2 WCSI", FULL, warm_core=False),
    "3 WCSI on tcident": run_wcsi("3 WCSI on tcident", TC_FULL, warm_core=False),
    "4 WCSI + warm core": run_wcsi("4 WCSI + warm core", FULL, warm_core=True),
  }
rows, overlap = [], []
for d in REPORT:
    sets = {a: storms(outs[a], d) for a in APPROACHES}
    for a, s in sets.items():
        g = np.abs([v[0] for v in s.values()])
        rows.append(dict(domain=d, approach=a, storms=len(s),
                         genesis_lat_mean=g.mean() if len(g) else np.nan,
                         pct_genesis_lt5=100 * (g < 5).mean() if len(g) else np.nan,
                         pct_genesis_gt30=100 * (g > 30).mean() if len(g) else np.nan,
                         median_max_wind=np.median([v[1] for v in s.values()]) if s else np.nan))
    for i, a in enumerate(APPROACHES):
        for b in APPROACHES[i + 1:]:
            A, B = set(sets[a]), set(sets[b])
            overlap.append(dict(domain=d, pair=f"{a.split()[0]} vs {b.split()[0]}", both=len(A & B), only_first=len(A - B), only_second=len(B - A)))
pd.set_option("display.width", 220)
r = pd.DataFrame(rows).round(1)
print(r.to_string(index=False))
o = pd.DataFrame(overlap)
print("\noverlaps (storms matched exactly on time + position):")
print(o.pivot_table(index="pair", columns="domain", values=["both", "only_first", "only_second"], aggfunc="first")
      .swaplevel(axis=1).sort_index(axis=1).reindex(columns=REPORT, level=0).to_string())
r.to_csv(f"{OUT}/summary.csv", index=False); o.to_csv(f"{OUT}/overlaps.csv", index=False)
