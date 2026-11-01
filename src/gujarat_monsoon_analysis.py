# Gujarat monsoon 2026: district rainfall, flood-risk and soil analysis
# Run from the repository root: python src/gujarat_monsoon_analysis.py



import os, glob, zipfile
def find_input(pattern):
    hits = glob.glob(pattern) + glob.glob(os.path.join("data", pattern))
    return hits[0] if hits else None

def unpack_zips():
    for z in glob.glob("*.zip") + glob.glob(os.path.join("data", "*.zip")):
        if zipfile.is_zipfile(z):
            with zipfile.ZipFile(z) as zf:
                zf.extractall("data" if os.path.isdir("data") else ".")

GRIB_PATH, KML_PATH = find_input("*.grib"), find_input("*.kml")
if GRIB_PATH is None:
    unpack_zips()                                   # e.g. data/data_grib.zip contains data.grib
    GRIB_PATH = find_input("*.grib")
if GRIB_PATH is None or KML_PATH is None:
    from google.colab import files                  # Colab: upload the data zip (or .grib) and the .kml
    up = files.upload()
    for name in up:
        if name.lower().endswith(".zip"):
            with zipfile.ZipFile(name) as zf:
                zf.extractall(".")
        elif name.lower().endswith(".kml"):
            KML_PATH = name
        elif name.lower().endswith(".grib"):
            GRIB_PATH = name
    GRIB_PATH = GRIB_PATH or find_input("*.grib")
    KML_PATH = KML_PATH or find_input("*.kml")
print("GRIB:", GRIB_PATH, "| KML:", KML_PATH)

import json, numpy as np, pandas as pd, eccodes, shutil
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt, matplotlib.colors as mcolors
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon as MplPolygon
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely.validation import make_valid
import shapely

OUT = "analysis_outputs"
SUB = {k: f"{OUT}/{k}" for k in ["tables", "timeseries", "district_maps", "district_grids", "heatmaps",
                                 "rankings", "grid_maps", "relationships", "insights"]}
for p in SUB.values():
    os.makedirs(p, exist_ok=True)
plt.rcParams.update({"figure.dpi": 100, "axes.grid": True, "grid.alpha": 0.3})

def read_grib(path):
    data, grid = {}, {}
    with open(path, "rb") as f:
        while True:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                break
            sn = eccodes.codes_get(gid, "shortName")
            vd = str(eccodes.codes_get(gid, "validityDate"))
            vt = eccodes.codes_get(gid, "validityTime")
            valid = pd.Timestamp(f"{vd[:4]}-{vd[4:6]}-{vd[6:]}") + pd.Timedelta(hours=vt // 100)
            ni, nj = eccodes.codes_get(gid, "Ni"), eccodes.codes_get(gid, "Nj")
            vals = eccodes.codes_get_array(gid, "values").reshape(nj, ni)
            if sn not in grid:
                grid[sn] = dict(name=eccodes.codes_get(gid, "name"), units=eccodes.codes_get(gid, "units"),
                                lat0=eccodes.codes_get(gid, "latitudeOfFirstGridPointInDegrees"),
                                lat1=eccodes.codes_get(gid, "latitudeOfLastGridPointInDegrees"),
                                lon0=eccodes.codes_get(gid, "longitudeOfFirstGridPointInDegrees"),
                                lon1=eccodes.codes_get(gid, "longitudeOfLastGridPointInDegrees"))
            data.setdefault(sn, {})[valid] = vals.astype(np.float32)
            eccodes.codes_release(gid)
    return data, grid

raw, grid_info = read_grib(GRIB_PATH)
ts = sorted(raw["2t"])
idx = pd.DatetimeIndex(ts)
g0 = grid_info["2t"]
nj, ni = raw["2t"][ts[0]].shape
lats = np.linspace(g0["lat0"], g0["lat1"], nj)
lons = np.linspace(g0["lon0"], g0["lon1"], ni)
print(f"Period {idx[0]} -> {idx[-1]} | {len(idx)} hours | grid {nj}x{ni}")
for sn, g in grid_info.items():
    print(f"  {sn:6s} {g['name']} [{g['units']}]")

def stack(var):
    return np.array([raw[var][t] for t in ts], dtype=np.float32)

# ---- physical units for every variable ----
VARS = {
    "t2m": stack("2t") - 273.15,            # °C
    "sp": stack("sp") / 100.0,              # hPa
    "tp": stack("tp") * 1000.0,             # mm per hour
    "cp": stack("cp") * 1000.0,             # mm per hour
    "crr": stack("crr") * 3600.0,           # mm/h rate
    "rain": (stack("ptype") > 0).astype(np.float32) * 100.0,   # % of hours with rain (per grid cell)
}
SLT = raw["slt"][ts[0]]                     # soil type (static)
print("Negative precipitation values:", int((VARS['tp'] < 0).sum() + (VARS['cp'] < 0).sum()))

KML_NS = "{http://www.opengis.net/kml/2.2}"
NAME_FIX = {"GANDH|NAGAR": "Gandhinagar", "G|R SOMNATH": "Gir Somnath", "DANGS": "The Dangs",
            "DEVBHUMI DWARKA": "Devbhumi Dwarka", "CHHOTA UDEPUR": "Chhota Udaipur"}

def ring_points(text):
    return [tuple(map(float, p.split(",")[:2])) for p in text.split()]

def load_kml_districts(path):
    root = ET.parse(path).getroot()
    parts = {}
    for pm in root.iter(KML_NS + "Placemark"):
        raw_name = (pm.findtext(KML_NS + "name") or "").strip()
        name = NAME_FIX.get(raw_name, raw_name.title())
        for poly in pm.iter(KML_NS + "Polygon"):
            outer = poly.find(f"{KML_NS}outerBoundaryIs/{KML_NS}LinearRing/{KML_NS}coordinates")
            if outer is None:
                continue
            holes = [ring_points(h.text) for h in
                     poly.findall(f"{KML_NS}innerBoundaryIs/{KML_NS}LinearRing/{KML_NS}coordinates")]
            p = Polygon(ring_points(outer.text), holes)
            parts.setdefault(name, []).append(p if p.is_valid else make_valid(p))
    return {n: unary_union(ps) for n, ps in parts.items()}

district_geoms = load_kml_districts(KML_PATH)
names = sorted(district_geoms)
polys = [district_geoms[n] for n in names]
gujarat = unary_union(polys)
print(f"{len(names)} districts loaded from KML")

LON, LAT = np.meshgrid(lons, lats)
flat_lon, flat_lat = LON.ravel(), LAT.ravel()
M = np.zeros((len(names), nj * ni), dtype=np.float32)
for i, g in enumerate(polys):
    M[i] = shapely.contains_xy(g, flat_lon, flat_lat)
# fallback: a district with no grid centre inside gets its nearest grid cell
for i, g in enumerate(polys):
    if M[i].sum() == 0:
        c = g.centroid
        M[i, np.argmin((flat_lon - c.x) ** 2 + (flat_lat - c.y) ** 2)] = 1
counts = M.sum(axis=1)
inside_cells = M.sum(axis=0) > 0                    # grid cells inside Gujarat
pd.DataFrame({"district": names, "grid_cells": counts.astype(int)}).to_csv(f"{SUB['tables']}/district_grid_cells.csv", index=False)
print("Districts with 3 or fewer grid cells:", [n for n, c in zip(names, counts) if c <= 3])

G = pd.DataFrame({k: v.reshape(len(idx), -1)[:, inside_cells].mean(axis=1) for k, v in VARS.items()}, index=idx)
G.columns = ["t2m", "sp", "tp", "cp", "crr", "rain"]

def district_hourly(arr):
    return pd.DataFrame((arr.reshape(len(idx), -1) @ M.T) / counts, index=idx, columns=names)

DH = {k: district_hourly(v) for k, v in VARS.items()}
T2_FLAT = VARS["t2m"].reshape(len(idx), -1)
DH["t2m_max"] = pd.DataFrame({n: T2_FLAT[:, M[i] > 0].max(axis=1) for i, n in enumerate(names)}, index=idx)
print("Hourly tables built:", list(G.columns), "| districts:", len(names))

RULES = {"daily": "D", "weekly": "W-MON", "monthly": "MS"}

def agg_gujarat(rule):
    r = G.resample(rule, label="left", closed="left")
    return pd.DataFrame({
        "t2m_mean_C": r.t2m.mean(), "t2m_max_C": r.t2m.max(), "t2m_min_C": r.t2m.min(),
        "sp_mean_hPa": r.sp.mean(), "precip_mm": r.tp.sum(), "convective_mm": r.cp.sum(),
        "convective_rate_mmh": r.crr.mean(), "rain_hours_pct": r.rain.mean(),
    })

GT = {k: agg_gujarat(v) for k, v in RULES.items()}
for k, df in GT.items():
    df.round(3).to_csv(f"{SUB['tables']}/gujarat_{k}_all_variables.csv")
    print(f"Gujarat {k}: {len(df)} rows")

def agg_district(df_h, rule, how):
    return getattr(df_h.resample(rule, label="left", closed="left"), how)()

DT = {}
for k, rule in RULES.items():
    DT[k] = {
        "precip": agg_district(DH["tp"], rule, "sum"),
        "conv": agg_district(DH["cp"], rule, "sum"),
        "t2m": agg_district(DH["t2m"], rule, "mean"),
        "t2m_max": agg_district(DH["t2m_max"], rule, "max"),
        "sp": agg_district(DH["sp"], rule, "mean"),
        "rain": agg_district(DH["rain"], rule, "mean"),
        "crr": agg_district(DH["crr"], rule, "mean"),
    }
    DT[k]["precip"].round(2).T.to_csv(f"{SUB['tables']}/district_{k}_precip_mm.csv")
    DT[k]["t2m"].round(2).T.to_csv(f"{SUB['tables']}/district_{k}_temp_C.csv")
print("District tables saved for daily, weekly and monthly")

def label_index(index, freq):
    return index.strftime("%d %b" if freq == "daily" else ("%d %b" if freq == "weekly" else "%b %Y"))

def plot_series(df, col, title, unit, fname, kind="line"):
    fig, ax = plt.subplots(figsize=(12, 4.5))
    if kind == "bar":
        ax.bar(df.index, df[col], width=0.8 if len(df) > 40 else (5 if len(df) > 8 else 20), color="tab:blue")
    else:
        ax.plot(df.index, df[col], lw=1.4, marker="o" if len(df) < 40 else None, ms=4)
    ax.set_title(title); ax.set_ylabel(unit)
    plt.tight_layout(); plt.savefig(fname, dpi=150); plt.close()

GUJ_PLOTS = [
    ("t2m_mean_C", "2 m temperature, Gujarat mean", "°C", "line"),
    ("t2m_max_C", "2 m temperature, Gujarat max", "°C", "line"),
    ("t2m_min_C", "2 m temperature, Gujarat min", "°C", "line"),
    ("sp_mean_hPa", "Surface pressure, Gujarat mean", "hPa", "line"),
    ("precip_mm", "Total precipitation, Gujarat", "mm", "bar"),
    ("convective_mm", "Convective precipitation, Gujarat", "mm", "bar"),
    ("convective_rate_mmh", "Convective rain rate, Gujarat mean", "mm/h", "line"),
    ("rain_hours_pct", "Share of hours with rain, Gujarat", "%", "line"),
]
for freq in RULES:
    for col, title, unit, kind in GUJ_PLOTS:
        plot_series(GT[freq], col, f"{title} ({freq})", unit,
                    f"{SUB['timeseries']}/gujarat_{freq}_{col}.png", kind)
# Combined stacked graph of all variables per frequency
for freq in RULES:
    df = GT[freq]
    fig, axes = plt.subplots(len(GUJ_PLOTS), 1, figsize=(13, 2.2 * len(GUJ_PLOTS)), sharex=True)
    for ax, (col, title, unit, kind) in zip(axes, GUJ_PLOTS):
        if kind == "bar":
            ax.bar(df.index, df[col], width=0.8 if freq == "daily" else (5 if freq == "weekly" else 20), color="tab:blue")
        else:
            ax.plot(df.index, df[col], lw=1.2)
        ax.set_ylabel(unit, fontsize=8); ax.set_title(title, fontsize=9, loc="left")
    fig.suptitle(f"All variables, Gujarat, {freq}", fontsize=13)
    plt.tight_layout(); plt.savefig(f"{SUB['timeseries']}/gujarat_{freq}_all_variables_stacked.png", dpi=150); plt.close()
print("Time-series PNGs saved")

LEVELS = {
    "precip_daily": [0, 1, 2, 5, 10, 20, 30, 50, 75, 100, 150, 250],
    "precip_weekly": [0, 10, 25, 50, 75, 100, 150, 200, 300, 400, 600, 800],
    "precip_monthly": [0, 20, 50, 100, 150, 200, 300, 400, 600, 800, 1000, 1500],
    "precip_period": [0, 100, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000],
    "temp": [16, 20, 22, 24, 26, 28, 30, 32, 34, 36, 40, 45],
    "pressure": [960, 970, 980, 985, 990, 995, 1000, 1005, 1010, 1015, 1020],
    "rain_pct": [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
}

def auto_levels(values, n=11):
    vmax = float(np.nanmax(values))
    return list(np.round(np.linspace(0, vmax * 1.001 + 1e-9, n), 4))

def draw_districts(ax, values, levels, cmap_name, annotate=False, fs=6):
    cmap = plt.get_cmap(cmap_name, len(levels) - 1)
    norm = mcolors.BoundaryNorm(levels, cmap.N)
    patches, cols = [], []
    for n, g in zip(names, polys):
        v = values.get(n, np.nan)
        colour = cmap(norm(v)) if np.isfinite(v) else (1, 1, 1, 1)
        for part in (g.geoms if hasattr(g, "geoms") else [g]):
            patches.append(MplPolygon(np.asarray(part.exterior.coords), closed=True))
            cols.append(colour)
    ax.add_collection(PatchCollection(patches, facecolor=cols, edgecolor="black", linewidth=0.4))
    if annotate:
        for n, g in zip(names, polys):
            v = values.get(n, np.nan)
            if np.isfinite(v):
                pt = g.representative_point()
                ax.text(pt.x, pt.y, f"{n}\n{v:.0f}", fontsize=fs, ha="center", va="center",
                        bbox=dict(fc="white", alpha=0.6, lw=0, pad=0.3))
    ax.set_xlim(68, 75); ax.set_ylim(20, 25.2); ax.set_aspect("equal")
    return cmap, norm

def district_map(values, title, levels, cmap_name, unit, fname):
    fig, ax = plt.subplots(figsize=(9, 8))
    cmap, norm = draw_districts(ax, values, levels, cmap_name, annotate=True)
    ax.set_title(title, fontsize=12); ax.set_xlabel("Longitude (°E)"); ax.set_ylabel("Latitude (°N)")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    fig.colorbar(sm, ax=ax, shrink=0.8, label=unit, ticks=levels, spacing="uniform")
    plt.savefig(fname, dpi=150, bbox_inches="tight"); plt.close()

def district_grid(table, levels, cmap_name, title, unit, fname, ncol=6, label_fmt="%d %b"):
    n = len(table)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 3.3, nrow * 3.4))
    axes = np.atleast_1d(axes).ravel()
    cmap = norm = None
    for ax, (t, row) in zip(axes, table.iterrows()):
        cmap, norm = draw_districts(ax, row.to_dict(), levels, cmap_name)
        ax.set_title(t.strftime(label_fmt), fontsize=9); ax.set_xticks([]); ax.set_yticks([])
    for ax in axes[n:]:
        ax.axis("off")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    fig.colorbar(sm, ax=axes.tolist(), shrink=0.5, label=unit, ticks=levels, spacing="uniform")
    fig.suptitle(title, fontsize=14)
    plt.savefig(fname, dpi=130, bbox_inches="tight"); plt.close()

month_names = {m: pd.Timestamp(2026, m, 1).strftime("%B %Y") for m in [7, 8, 9, 10]}
mp = DT["monthly"]["precip"]
for t in mp.index:
    district_map(mp.loc[t].to_dict(), f"Total precipitation by district: {t:%B %Y}",
                 LEVELS["precip_monthly"], "turbo", "mm", f"{SUB['district_maps']}/precip_monthly_{t:%m}.png")
    district_map(DT["monthly"]["t2m"].loc[t].to_dict(), f"Mean temperature by district: {t:%B %Y}",
                 LEVELS["temp"], "turbo", "°C", f"{SUB['district_maps']}/temp_monthly_{t:%m}.png")
    district_map(DT["monthly"]["sp"].loc[t].to_dict(), f"Mean surface pressure by district: {t:%B %Y}",
                 LEVELS["pressure"], "turbo", "hPa", f"{SUB['district_maps']}/pressure_monthly_{t:%m}.png")
    district_map(DT["monthly"]["rain"].loc[t].to_dict(), f"Share of hours with rain by district: {t:%B %Y}",
                 LEVELS["rain_pct"], "turbo", "%", f"{SUB['district_maps']}/rainhours_monthly_{t:%m}.png")

period_total = DH["tp"].sum()
district_map(period_total.to_dict(), "Total precipitation by district: full period",
             LEVELS["precip_period"], "turbo", "mm",
             f"{SUB['district_maps']}/precip_total_period.png")
period_conv = DH["cp"].sum()
district_map(period_conv.to_dict(), "Convective precipitation by district: full period",
             auto_levels(period_conv), "turbo", "mm", f"{SUB['district_maps']}/convective_total_period.png")
district_map(DH["sp"].mean().to_dict(), "Mean surface pressure by district: full period",
             LEVELS["pressure"], "turbo", "hPa", f"{SUB['district_maps']}/pressure_period.png")
district_map(DH["t2m"].mean().to_dict(), "Mean temperature by district: full period",
             LEVELS["temp"], "turbo", "°C", f"{SUB['district_maps']}/temp_period.png")
district_map(DH["t2m_max"].max().to_dict(), "Highest temperature by district: full period",
             LEVELS["temp"], "turbo", "°C", f"{SUB['district_maps']}/temp_max_period.png")
district_map(DH["rain"].mean().to_dict(), "Share of hours with rain by district: full period",
             LEVELS["rain_pct"], "turbo", "%", f"{SUB['district_maps']}/rainhours_period.png")
print("District maps saved")

dp = DT["daily"]["precip"]
for i, start in enumerate(range(0, len(dp), 24)):
    chunk = dp.iloc[start:start + 24]
    district_grid(chunk, LEVELS["precip_daily"], "turbo", f"Daily precipitation by district, part {i+1}",
                  "mm/day", f"{SUB['district_grids']}/precip_daily_part{i+1}.png")
wp = DT["weekly"]["precip"]
district_grid(wp, LEVELS["precip_weekly"], "turbo", "Weekly precipitation by district (week starting Monday)",
              "mm/week", f"{SUB['district_grids']}/precip_weekly_all.png", label_fmt="week %d %b")
district_grid(DT["weekly"]["t2m"], LEVELS["temp"], "turbo", "Weekly mean temperature by district",
              "°C", f"{SUB['district_grids']}/temp_weekly_all.png", label_fmt="week %d %b")
print("District grid PNGs saved")

def heatmap(table, title, fname, cmap="turbo", fmt="{:.0f}", unit=""):
    wide = table.shape[1] > 40
    fig, ax = plt.subplots(figsize=(min(45, max(7, 0.3 * table.shape[1] + 4)), max(7, 0.28 * table.shape[0])))
    im = ax.imshow(table.values, aspect="auto", cmap=cmap)
    ax.set_xticks(range(table.shape[1])); ax.set_xticklabels([c.strftime("%d %b") if hasattr(c, "strftime") else str(c) for c in table.columns], rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(table.shape[0])); ax.set_yticklabels(table.index, fontsize=7)
    if not wide:
        for i in range(table.shape[0]):
            for j in range(table.shape[1]):
                ax.text(j, i, fmt.format(table.values[i, j]), ha="center", va="center", fontsize=6, color="white")
    ax.set_title(title); ax.grid(False)
    plt.colorbar(im, ax=ax, shrink=0.6, label=unit)
    plt.tight_layout(); plt.savefig(fname, dpi=150); plt.close()

heatmap(DT["monthly"]["precip"].T.set_axis([c.strftime("%b") for c in DT["monthly"]["precip"].index], axis=1),
        "Monthly precipitation by district (mm)", f"{SUB['heatmaps']}/district_month_precip.png", unit="mm")
heatmap(DT["weekly"]["precip"].T, "Weekly precipitation by district (mm)",
        f"{SUB['heatmaps']}/district_week_precip.png", unit="mm")
heatmap(DT["daily"]["precip"].T, "Daily precipitation by district (mm)",
        f"{SUB['heatmaps']}/district_day_precip.png", unit="mm", fmt="{:.0f}")
heatmap(DT["monthly"]["t2m"].T.set_axis([c.strftime("%b") for c in DT["monthly"]["t2m"].index], axis=1),
        "Monthly mean temperature by district (°C)", f"{SUB['heatmaps']}/district_month_temp.png", cmap="coolwarm", fmt="{:.1f}", unit="°C")

summary = pd.DataFrame({
    "precip_total_mm": DH["tp"].sum(),
    "precip_daily_max_mm": DT["daily"]["precip"].max(),
    "rain_days_1mm": (DT["daily"]["precip"] >= 1).sum(),
    "heavy_days_20mm": (DT["daily"]["precip"] >= 20).sum(),
    "convective_total_mm": DH["cp"].sum(),
    "temp_mean_C": DH["t2m"].mean(),
    "temp_max_C": DH["t2m_max"].max(),
    "pressure_mean_hPa": DH["sp"].mean(),
    "rain_hours_pct": DH["rain"].mean(),
}).round(2)
summary["rain_per_rain_hour_mm"] = (summary.precip_total_mm / (summary.rain_hours_pct / 100 * len(idx))).round(3)
summary["convective_share_pct"] = (summary.convective_total_mm / summary.precip_total_mm * 100).round(1)
summary.index.name = "district"
summary.sort_values("precip_total_mm", ascending=False).to_csv(f"{SUB['tables']}/district_summary.csv")

for col, title, unit, asc in [
    ("precip_total_mm", "Total rainfall by district", "mm", True),
    ("temp_mean_C", "Mean temperature by district", "°C", True),
    ("rain_hours_pct", "Share of hours with rain", "%", True),
    ("convective_share_pct", "Convective share of rainfall", "%", True),
    ("rain_per_rain_hour_mm", "Rainfall intensity (mm per rain hour)", "mm", True),
    ("pressure_mean_hPa", "Mean surface pressure by district", "hPa", True),
]:
    s = summary[col].sort_values(ascending=asc)
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.barh(s.index, s.values, color=plt.cm.turbo(np.linspace(0.05, 0.95, len(s))))
    ax.set_xlabel(unit); ax.set_title(title); ax.grid(axis="x", alpha=0.3)
    plt.tight_layout(); plt.savefig(f"{SUB['rankings']}/rank_{col}.png", dpi=150); plt.close()
print("Heatmaps and rankings saved")

GRID_PARAMS = {
    "t2m": dict(title="2 m temperature (°C), monthly mean", agg="mean", unit="°C", levels=LEVELS["temp"], data=VARS["t2m"]),
    "sp": dict(title="Surface pressure (hPa), monthly mean", agg="mean", unit="hPa", levels=LEVELS["pressure"], data=VARS["sp"]),
    "tp": dict(title="Total precipitation (mm), monthly total", agg="sum", unit="mm", levels=LEVELS["precip_monthly"], data=VARS["tp"]),
    "cp": dict(title="Convective precipitation (mm), monthly total", agg="sum", unit="mm", levels=auto_levels(VARS["cp"].reshape(len(idx), -1).sum(0)), data=VARS["cp"]),
    "crr": dict(title="Convective rain rate (mm/h), monthly mean", agg="mean", unit="mm/h", levels=auto_levels(VARS["crr"].mean(0)), data=VARS["crr"]),
    "rain": dict(title="Share of hours with rain (%), monthly", agg="mean", unit="%", levels=LEVELS["rain_pct"], data=VARS["rain"]),
}
dlon, dlat = lons[1] - lons[0], lats[0] - lats[1]
edge_lon = np.append(lons - dlon / 2, lons[-1] + dlon / 2)
edge_lat = np.append(lats - dlat / 2, lats[-1] + dlat / 2)
gu_geoms = gujarat.geoms if hasattr(gujarat, "geoms") else [gujarat]
inside_grid = inside_cells.reshape(nj, ni)
MONTHS = [(7, "July 2026"), (8, "August 2026"), (9, "September 2026"), (10, "October 2026 (1-26 Oct)")]

def grid_field(arr, m, agg):
    sel = idx.month == m
    f = arr[sel].sum(axis=0) if agg == "sum" else arr[sel].mean(axis=0)
    return np.where(inside_grid, f, np.nan)

def grid_panel(ax, field, cmap, norm, title, label):
    ax.pcolormesh(edge_lon, edge_lat, field, cmap=cmap, norm=norm, shading="flat")
    for g in gu_geoms:
        x, y = g.exterior.xy
        ax.plot(x, y, color="black", lw=0.8)
    ax.set_xlim(68, 75); ax.set_ylim(20, 25); ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)
    ax.text(0.02, 0.02, label, transform=ax.transAxes, fontsize=9, fontweight="bold",
            bbox=dict(fc="white", alpha=0.85, lw=0))

for key, p in GRID_PARAMS.items():
    cmap = plt.get_cmap("turbo", len(p["levels"]) - 1).copy(); cmap.set_bad("white")
    norm = mcolors.BoundaryNorm(p["levels"], cmap.N)
    fields = {m: grid_field(p["data"], m, p["agg"]) for m, _ in MONTHS}
    fig, axes = plt.subplots(2, 2, figsize=(12, 11))
    for ax, (m, label) in zip(axes.ravel(), MONTHS):
        grid_panel(ax, fields[m], cmap, norm, f"{key}: {label}", f"mean {np.nanmean(fields[m]):.2f} {p['unit']}")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    fig.colorbar(sm, ax=axes, fraction=0.035, pad=0.03, label=p["title"], ticks=p["levels"], spacing="uniform")
    fig.suptitle(f"{p['title']}, Gujarat (same colour scale)", fontsize=13)
    plt.savefig(f"{SUB['grid_maps']}/{key}_4months.png", dpi=150, bbox_inches="tight"); plt.close()
    for m, label in MONTHS:
        fig, ax = plt.subplots(figsize=(7, 6.5))
        grid_panel(ax, fields[m], cmap, norm, f"{key}: {label}", f"mean {np.nanmean(fields[m]):.2f} {p['unit']}")
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
        fig.colorbar(sm, ax=ax, shrink=0.85, label=p["title"], ticks=p["levels"], spacing="uniform")
        plt.savefig(f"{SUB['grid_maps']}/{key}_{m:02d}.png", dpi=150, bbox_inches="tight"); plt.close()

# soil type (static)
cmap = plt.get_cmap("turbo", 5).copy(); cmap.set_bad("white")
norm = mcolors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], cmap.N)
fig, ax = plt.subplots(figsize=(7, 6.5))
grid_panel(ax, np.where(inside_grid, SLT, np.nan), cmap, norm, "slt: soil type (static)", "codes 0-4")
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
fig.colorbar(sm, ax=ax, shrink=0.85, ticks=[0, 1, 2, 3, 4], label="Soil type code (Code table 4.213)")
plt.savefig(f"{SUB['grid_maps']}/slt_static.png", dpi=150, bbox_inches="tight"); plt.close()
print("Grid maps saved")

D = GT["daily"]
corr = D[["t2m_mean_C", "sp_mean_hPa", "precip_mm", "convective_mm", "convective_rate_mmh", "rain_hours_pct"]].corr().round(2)
corr.to_csv(f"{SUB['tables']}/daily_correlation.csv")
fig, ax = plt.subplots(figsize=(8, 7))
im = ax.imshow(corr.values, cmap="RdBu_r", vmin=-1, vmax=1)
ax.set_xticks(range(len(corr))); ax.set_xticklabels(corr.columns, rotation=45, ha="right", fontsize=8)
ax.set_yticks(range(len(corr))); ax.set_yticklabels(corr.index, fontsize=8)
for i in range(len(corr)):
    for j in range(len(corr)):
        ax.text(j, i, f"{corr.values[i, j]:.2f}", ha="center", va="center", fontsize=8)
ax.set_title("Correlation between daily Gujarat variables"); ax.grid(False)
plt.colorbar(im, ax=ax, shrink=0.8)
plt.tight_layout(); plt.savefig(f"{SUB['relationships']}/daily_correlation_heatmap.png", dpi=150); plt.close()

def scatter(x, y, xl, yl, title, fname):
    fig, ax = plt.subplots(figsize=(7, 5.5))
    ax.scatter(x, y, s=18, alpha=0.7)
    m, b = np.polyfit(x, y, 1)
    xs = np.linspace(np.min(x), np.max(x), 50)
    ax.plot(xs, m * xs + b, color="red", lw=1.2)
    r = np.corrcoef(x, y)[0, 1]
    ax.set_title(f"{title} (r = {r:.2f})"); ax.set_xlabel(xl); ax.set_ylabel(yl)
    plt.tight_layout(); plt.savefig(fname, dpi=150); plt.close()

scatter(D.sp_mean_hPa, D.precip_mm, "Daily mean pressure (hPa)", "Daily precipitation (mm)",
        "Pressure vs rainfall", f"{SUB['relationships']}/scatter_pressure_rain.png")
scatter(D.t2m_mean_C, D.rain_hours_pct, "Daily mean temperature (°C)", "Share of hours with rain (%)",
        "Temperature vs rain hours", f"{SUB['relationships']}/scatter_temp_rainhours.png")
scatter(D.convective_rate_mmh, D.precip_mm, "Convective rain rate (mm/h)", "Daily precipitation (mm)",
        "Convective rate vs rainfall", f"{SUB['relationships']}/scatter_convective_rain.png")
scatter(summary.pressure_mean_hPa, summary.precip_total_mm, "Mean pressure by district (hPa)",
        "Total rainfall by district (mm)", "District pressure vs rainfall", f"{SUB['relationships']}/district_scatter_pressure_rain.png")
scatter(summary.temp_mean_C, summary.precip_total_mm, "Mean temperature by district (°C)",
        "Total rainfall by district (mm)", "District temperature vs rainfall", f"{SUB['relationships']}/district_scatter_temp_rain.png")

# Daily cumulative rain and dry/wet spell chart
wet = (D.precip_mm >= 1).astype(int)
fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
axes[0].plot(D.index, D.precip_mm.cumsum(), color="tab:blue")
axes[0].set_ylabel("Cumulative mm"); axes[0].set_title("Cumulative Gujarat rainfall")
axes[1].fill_between(D.index, 0, wet, step="mid", color="tab:green", alpha=0.6)
axes[1].set_yticks([0, 1]); axes[1].set_yticklabels(["dry (<1 mm)", "wet (>=1 mm)"]); axes[1].set_title("Wet and dry days")
plt.tight_layout(); plt.savefig(f"{SUB['relationships']}/cumulative_and_wet_dry_days.png", dpi=150); plt.close()
print("Relationship PNGs saved")

ins = []
ins.append(f"Total Gujarat rainfall over {idx[0]:%d %b} to {idx[-1]:%d %b %Y}: {G.tp.sum():.0f} mm (area mean).")
mon = GT["monthly"].precip_mm
ins.append("Monthly rainfall: " + ", ".join(f"{d:%b} {v:.0f} mm ({v / mon.sum() * 100:.0f}% of total)" for d, v in mon.items()) + ".")
top_day = D.precip_mm.idxmax()
ins.append(f"Wettest day in Gujarat was {top_day:%d %b} with {D.precip_mm.max():.1f} mm (area mean).")
wk = GT["weekly"].precip_mm
ins.append(f"Wettest week started {wk.idxmax():%d %b} with {wk.max():.0f} mm (area mean).")
ds_tot = DH["tp"].sum()
ins.append(f"Wettest district: {ds_tot.idxmax()} ({ds_tot.max():.0f} mm); driest: {ds_tot.idxmin()} ({ds_tot.min():.0f} mm). "
           f"Ratio wettest to driest: {ds_tot.max() / ds_tot.min():.1f}x.")
top5_share = ds_tot.nlargest(5).sum() / ds_tot.sum() * 100
ins.append(f"The five wettest districts received {top5_share:.0f}% of the district-summed rainfall.")
hot_t = G.t2m.idxmax()
ins.append(f"Highest temperature in any district: {DH['t2m_max'].max().max():.1f} °C in {DH['t2m_max'].max().idxmax()}; "
           f"hottest Gujarat-mean hour was {hot_t:%d %b %H:%M} UTC.")
tmean = DH["t2m"].mean()
ins.append(f"Warmest district on average: {tmean.idxmax()} ({tmean.max():.1f} °C); coolest: {tmean.idxmin()} ({tmean.min():.1f} °C).")
ins.append(f"Lowest Gujarat-mean pressure: {G.sp.min():.1f} hPa on {G.sp.idxmin():%d %b}; "
           f"pressure rose from {GT['monthly'].sp_mean_hPa.iloc[0]:.0f} hPa (Jul mean) to {GT['monthly'].sp_mean_hPa.iloc[-1]:.0f} hPa (Oct mean).")
# monsoon onset / withdrawal heuristic: 3-day rain total of at least 30 mm
r3 = D.precip_mm.rolling(3).sum()
wet_roll = r3[r3 >= 30]
if len(wet_roll):
    ins.append(f"Heavy 3-day spells (>=30 mm) occurred between {wet_roll.index[0]:%d %b} and {wet_roll.index[-1]:%d %b}; "
               f"the last one ended {wet_roll.index[-1]:%d %b}, after which rain was mostly light.")
# dry and wet spells
def longest_run(mask):
    best = cur = 0; end = None; start = None; cur_start = None
    for d, v in mask.items():
        if v:
            cur += 1
            if cur == 1: cur_start = d
            if cur > best: best, start, end = cur, cur_start, d
        else:
            cur = 0
    return best, start, end
b, s, e = longest_run(D.precip_mm < 1)
ins.append(f"Longest dry spell (all days below 1 mm): {b} days, {s:%d %b} to {e:%d %b}.")
b, s, e = longest_run(D.precip_mm >= 1)
ins.append(f"Longest wet spell (all days at least 1 mm): {b} days, {s:%d %b} to {e:%d %b}.")
heavy_days = (DT["daily"]["precip"] >= 20).sum().sort_values(ascending=False)
ins.append(f"Most heavy-rain days (>=20 mm): {heavy_days.index[0]} ({heavy_days.iloc[0]} days); "
           f"{(heavy_days == 0).sum()} districts had none.")
intensity = summary.rain_per_rain_hour_mm
ins.append(f"Most intense rainfall per rain hour: {intensity.idxmax()} ({intensity.max():.2f} mm); least: {intensity.idxmin()} ({intensity.min():.2f} mm).")
conv = summary.convective_share_pct
ins.append(f"Highest convective share of rainfall: {conv.idxmax()} ({conv.max():.0f}%); lowest: {conv.idxmin()} ({conv.min():.0f}%).")
hot_days = (DH["t2m_max"] >= 40).sum()
ins.append(f"Hours with a district maximum of at least 40 °C: most in {hot_days.idxmax()} ({hot_days.max()} hours); "
           f"{(hot_days > 0).sum()} of {len(names)} districts reached 40 °C at least once.")
ins.append(f"Daily correlation: pressure vs rainfall r = {corr.loc['sp_mean_hPa','precip_mm']:.2f}; "
           f"temperature vs rain-hour share r = {corr.loc['t2m_mean_C','rain_hours_pct']:.2f}; "
           f"convective rate vs rainfall r = {corr.loc['convective_rate_mmh','precip_mm']:.2f}.")
dist_r = np.corrcoef(summary.pressure_mean_hPa, summary.precip_total_mm)[0, 1]
ins.append(f"Across districts, mean pressure vs total rainfall r = {dist_r:.2f}.")
wk_cover = (DT["weekly"]["precip"] >= 100).sum(axis=1)
ins.append(f"Week with the most districts receiving at least 100 mm: week of {wk_cover.idxmax():%d %b} ({wk_cover.max()} districts).")
mon_dry = (DT["monthly"]["precip"] < 50).sum(axis=1)
ins.append("Districts under 50 mm by month: " + ", ".join(f"{d:%b} {v}" for d, v in mon_dry.items()) + ".")

with open(f"{SUB['insights']}/insights.md", "w") as f:
    f.write("# Findings from the data\n\n")
    for line in ins:
        f.write(f"- {line}\n")
for line in ins:
    print("-", line)

# Insight chart: rainfall intensity and convective share together
fig, ax = plt.subplots(figsize=(9, 8))
ax.scatter(summary.rain_hours_pct, summary.precip_total_mm, s=summary.convective_share_pct * 4 + 20,
           c=summary.temp_mean_C, cmap="turbo", alpha=0.85, edgecolor="black")
for n, row in summary.iterrows():
    ax.annotate(n, (row.rain_hours_pct, row.precip_total_mm), fontsize=7, xytext=(4, 3), textcoords="offset points")
ax.set_xlabel("Share of hours with rain (%)"); ax.set_ylabel("Total rainfall (mm)")
ax.set_title("Districts: how often it rains vs how much (bubble = convective share, colour = mean temp)")
plt.tight_layout(); plt.savefig(f"{SUB['insights']}/districts_frequency_vs_amount.png", dpi=150); plt.close()

with open(GRIB_PATH, "rb") as f:
    _gid = eccodes.codes_grib_new_from_file(f)
    prov = {k: eccodes.codes_get(_gid, k) for k in ["centreDescription", "generatingProcessIdentifier", "stepType", "gridType"]}
    eccodes.codes_release(_gid)
prov_text = (
    f"Producer: {prov['centreDescription']} (generating process {prov['generatingProcessIdentifier']})\n"
    f"Grid: {prov['gridType']}, {nj} x {ni} points, step type: {prov['stepType']}\n"
    f"Variables: {', '.join(sorted(raw))}\n"
    f"Period: {idx[0]} to {idx[-1]} ({len(idx)} hours)\n"
    "Note: the file has forecast-style step numbers. Confirm the exact product name from the download source.\n"
)
with open(f"{SUB['tables']}/data_provenance.txt", "w") as f:
    f.write(prov_text)
print(prov_text)

jul = idx.month == 7
tp_jul = VARS["tp"][jul]                                  # (hours, lat, lon)
jul_total = tp_jul.sum(axis=0).astype(float)
jul_total[~inside_grid] = np.nan
iy, ix = np.unravel_index(np.nanargmax(jul_total), jul_total.shape)
win = jul_total[max(iy - 1, 0):iy + 2, max(ix - 1, 0):ix + 2]
neigh_vals = np.delete(win.ravel(), np.argmax(win.ravel() == jul_total[iy, ix]))
neigh_vals = neigh_vals[np.isfinite(neigh_vals)]
cell_series = pd.Series(tp_jul[:, iy, ix], index=idx[jul])
neigh_series = pd.Series(np.nanmean(tp_jul[:, max(iy - 1, 0):iy + 2, max(ix - 1, 0):ix + 2].reshape(tp_jul.shape[0], -1), axis=1), index=idx[jul])
hotspot = {
    "lat": round(float(lats[iy]), 2), "lon": round(float(lons[ix]), 2),
    "july_total_mm": round(float(jul_total[iy, ix]), 1),
    "neighbour_mean_mm": round(float(neigh_vals.mean()), 1),
    "ratio_to_neighbours": round(float(jul_total[iy, ix] / neigh_vals.mean()), 2),
    "hours_with_rain": int((cell_series > 0.1).sum()),
    "share_from_top10_hours_pct": round(float(cell_series.nlargest(10).sum() / cell_series.sum() * 100), 1),
}
pd.Series(hotspot).to_csv(f"{SUB['tables']}/july_hotspot_check.csv", header=["value"])
print("July hotspot check:", hotspot)
print("Read: a ratio far above 1 with rain concentrated in a few hours suggests a single-cell spike. "
      "Spatially coherent rain gives a ratio close to 1 and rain spread over many hours.")
fig, ax = plt.subplots(figsize=(12, 4.5))
ax.plot(cell_series.index, cell_series.values, lw=0.8, label=f"Hotspot cell ({hotspot['lat']}N, {hotspot['lon']}E)")
ax.plot(neigh_series.index, neigh_series.values, lw=0.8, label="Neighbouring cells (mean)")
ax.set_title("July hourly rainfall: hotspot cell vs neighbours"); ax.set_ylabel("mm per hour"); ax.legend()
plt.tight_layout(); plt.savefig(f"{SUB['tables']}/july_hotspot_timeseries.png", dpi=150); plt.close()

max3 = DH["tp"].rolling(72, min_periods=72).sum().max()                     # largest 3-day total (mm)
heavy = (DT["daily"]["precip"] >= 20).sum()                                 # days with at least 20 mm
top5_share = DT["daily"]["precip"].apply(lambda s: s.nlargest(5).sum() / max(s.sum(), 1e-6) * 100)
risk = pd.DataFrame({
    "max_3day_mm": max3, "heavy_days_20mm": heavy,
    "rain_per_rain_hour_mm": summary.rain_per_rain_hour_mm,
    "top5_days_share_pct": top5_share,
}).round(2)
norm = (risk - risk.min()) / (risk.max() - risk.min())
risk["risk_index_0_to_1"] = norm.mean(axis=1).round(3)
risk = risk.sort_values("risk_index_0_to_1", ascending=False)
risk.index.name = "district"
risk.to_csv(f"{SUB['tables']}/district_flood_risk_index.csv")
print(risk.head(10).to_string())
print("Index = average of the four measures, each scaled 0 (lowest district) to 1 (highest district).")

district_map(risk["risk_index_0_to_1"].to_dict(), "Flood-risk index by district (0 = lowest, 1 = highest)",
             [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], "turbo", "index",
             f"{SUB['insights']}/flood_risk_index_map.png")
fig, ax = plt.subplots(figsize=(9, 9))
rs = risk.sort_values("risk_index_0_to_1")
ax.barh(rs.index, rs.risk_index_0_to_1, color=plt.cm.turbo(rs.risk_index_0_to_1.values))
ax.set_xlabel("Flood-risk index (0 to 1)"); ax.set_title("Flood-risk index by district"); ax.grid(axis="x", alpha=0.3)
plt.tight_layout(); plt.savefig(f"{SUB['insights']}/flood_risk_index_ranking.png", dpi=150); plt.close()

# Add a multi-year GRIB (same grid and variables) named baseline*.grib to compare 2026 with a normal year.
# Without it this cell only reports that no anomaly is computed.
import glob as _glob
base_files = _glob.glob("baseline*.grib")
if not base_files:
    print("No baseline file found. Add baseline*.grib (several past years, same grid) to enable anomaly maps.")
else:
    braw, _ = read_grib(base_files[0])
    bts = sorted(braw["tp"])
    bidx = pd.DatetimeIndex(bts)
    btp = np.array([braw["tp"][t] for t in bts], dtype=np.float32) * 1000.0
    rows = []
    for yr in sorted(set(bidx.year)):
        for m in [7, 8, 9, 10]:
            sel = (bidx.year == yr) & (bidx.month == m)
            if sel.any():
                rows.append((yr, m, btp[sel].sum(axis=0)))
    clim = {m: np.mean([r[2] for r in rows if r[1] == m], axis=0) for m in [7, 8, 9, 10]}
    anomaly_rows = {}
    for m, label in MONTHS:
        now = DT["monthly"]["precip"].loc[DT["monthly"]["precip"].index.month == m].iloc[0]
        base_district = (clim[m].ravel() @ M.T) / counts
        anomaly_rows[label] = pd.Series(now.values / np.maximum(base_district, 1e-6) * 100, index=names)
    anom = pd.DataFrame(anomaly_rows).round(1)
    anom.index.name = "district"
    anom.to_csv(f"{SUB['tables']}/district_monthly_anomaly_pct.csv")
    print("Monthly rainfall as % of baseline average (district):")
    print(anom.to_string())
    district_map(anom.iloc[:, 0].to_dict(), f"Rainfall as % of baseline: {MONTHS[0][1]}",
                 [0, 50, 75, 100, 125, 150, 200, 300, 500, 1000, 2000], "RdBu", "% of baseline",
                 f"{SUB['district_maps']}/anomaly_july_pct.png")

# Class names follow WMO code table 4.213 (Soil type), as published in the eccodes definitions.
SOIL_NAMES = {0: "Reserved / no class (likely water)", 1: "Sand", 2: "Loamy sand", 3: "Sandy loam", 4: "Silt loam"}
SOIL_COLOURS = {1: "#f4d03f", 2: "#e59866", 3: "#58d68d", 4: "#3498db"}
slt_grid = SLT.astype(int)
slt_in = np.where(inside_grid, slt_grid, -1)            # -1 = outside Gujarat

# cell-level variables inside Gujarat
cell_rain = VARS["tp"].sum(axis=0)
cell_rain_pct = VARS["rain"].mean(axis=0)
cell_conv_share = np.where(cell_rain > 0, VARS["cp"].sum(axis=0) / np.maximum(cell_rain, 1e-9) * 100, np.nan)
cell_temp = VARS["t2m"].mean(axis=0)
cell_pres = VARS["sp"].mean(axis=0)
cells = pd.DataFrame({
    "soil_code": slt_in[inside_grid],
    "rain_mm_total": cell_rain[inside_grid],
    "rain_hours_pct": cell_rain_pct[inside_grid],
    "convective_share_pct": cell_conv_share[inside_grid],
    "temp_mean_C": cell_temp[inside_grid],
    "pressure_mean_hPa": cell_pres[inside_grid],
})
cells["soil_name"] = cells.soil_code.map(SOIL_NAMES)
soil_class = cells.groupby("soil_code").agg(
    grid_cells=("soil_code", "size"),
    rain_mm_total_mean=("rain_mm_total", "mean"),
    rain_hours_pct_mean=("rain_hours_pct", "mean"),
    convective_share_pct_mean=("convective_share_pct", "mean"),
    temp_mean_C=("temp_mean_C", "mean"),
    pressure_mean_hPa=("pressure_mean_hPa", "mean"),
).round(2)
soil_class.insert(0, "name", [SOIL_NAMES[c] for c in soil_class.index])
soil_class.to_csv(f"{SUB['tables']}/soil_class_summary.csv")
print(soil_class.to_string())

# district composition and dominant class (excluding "reserved" cells)
comp = {}
dominant = {}
for i, n in enumerate(names):
    v = slt_grid.ravel()[M[i] > 0]
    comp[n] = {c: float(np.mean(v == c) * 100) for c in [0, 1, 2, 3, 4]}
    real = v[v > 0]
    dominant[n] = int(np.bincount(real).argmax()) if real.size else 0
comp = pd.DataFrame(comp).T.round(1)
comp.columns = [SOIL_NAMES[c] for c in comp.columns]
comp.index.name = "district"
comp["dominant_soil"] = [SOIL_NAMES[dominant[n]] for n in comp.index]
comp.to_csv(f"{SUB['tables']}/district_soil_composition.csv")

fig, ax = plt.subplots(figsize=(9, 9))
order = sorted(names, key=lambda n: dominant[n])
left = np.zeros(len(order))
for c in [1, 2, 3, 4, 0]:
    vals = np.array([comp.loc[n, SOIL_NAMES[c]] for n in order])
    ax.barh(order, vals, left=left, color=SOIL_COLOURS.get(c, "#bbbbbb"), label=SOIL_NAMES[c])
    left += vals
ax.set_xlabel("Share of grid cells in district (%)"); ax.set_title("Soil composition by district")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=8)
plt.tight_layout(); plt.savefig(f"{SUB['heatmaps']}/district_soil_composition.png", dpi=150); plt.close()

# categorical map: dominant soil class per district
fig, ax = plt.subplots(figsize=(9, 8))
for n, g in zip(names, polys):
    colour = SOIL_COLOURS.get(dominant[n], "#bbbbbb")
    for part in (g.geoms if hasattr(g, "geoms") else [g]):
        ax.add_patch(MplPolygon(np.asarray(part.exterior.coords), closed=True, fc=colour, ec="black", lw=0.4))
    pt = g.representative_point()
    ax.text(pt.x, pt.y, n, fontsize=6, ha="center", va="center")
ax.set_xlim(68, 75); ax.set_ylim(20, 25.2); ax.set_aspect("equal")
ax.set_title("Dominant soil type by district (WMO code table 4.213)")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(fc=SOIL_COLOURS[c], ec="black", label=SOIL_NAMES[c]) for c in [1, 2, 3, 4]],
          loc="lower left", fontsize=9)
plt.savefig(f"{SUB['district_maps']}/soil_dominant_map.png", dpi=150, bbox_inches="tight"); plt.close()

# ASSUMPTION (hydrology, not from this dataset): finer soils infiltrate less water and so produce more runoff.
RUNOFF_CLASS = {1: 1, 2: 1, 3: 2, 4: 3}     # 1 low, 2 moderate, 3 high
RUNOFF_NAME = {1: "low", 2: "moderate", 3: "high"}
suscept = {}
for i, n in enumerate(names):
    v = slt_grid.ravel()[M[i] > 0]
    v = v[v > 0]
    suscept[n] = float(np.mean([RUNOFF_CLASS[x] for x in v])) if v.size else np.nan
dist_suscept = pd.Series(suscept, name="runoff_susceptibility_1to3")

quad = risk[["risk_index_0_to_1"]].join(dist_suscept)
quad["dominant_soil"] = [SOIL_NAMES[dominant[n]] for n in quad.index]
quad.round(3).to_csv(f"{SUB['tables']}/district_soil_flood_quadrant.csv")

fig, ax = plt.subplots(figsize=(9, 7))
ax.scatter(quad.runoff_susceptibility_1to3, quad.risk_index_0_to_1, s=60, c="tab:blue", edgecolor="black")
for n, row in quad.iterrows():
    ax.annotate(n, (row.runoff_susceptibility_1to3, row.risk_index_0_to_1), fontsize=7, xytext=(4, 3), textcoords="offset points")
mx, my = quad.runoff_susceptibility_1to3.median(), quad.risk_index_0_to_1.median()
ax.axvline(mx, color="grey", ls="--", lw=0.8); ax.axhline(my, color="grey", ls="--", lw=0.8)
ax.set_xlabel("Runoff susceptibility from soil (1 = low, 3 = high; assumption-based)")
ax.set_ylabel("Flood-risk index (data-based, 0 to 1)")
ax.set_title("Soil runoff susceptibility vs rainfall-based flood risk")
ax.text(0.02, 0.97, "top-right: high rain risk AND high soil susceptibility", transform=ax.transAxes, va="top", fontsize=8)
plt.tight_layout(); plt.savefig(f"{SUB['insights']}/soil_runoff_vs_flood_risk.png", dpi=150); plt.close()

# rainfall and convection by soil class at cell level (confounded by location; see notes)
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
cls = cells[cells.soil_code > 0].groupby("soil_name")
for ax, col, title, unit in zip(axes, ["rain_mm_total", "rain_hours_pct", "convective_share_pct"],
                                ["Total rain per cell (Jul–Oct)", "Share of hours with rain", "Convective share of rain"],
                                ["mm", "%", "%"]):
    m = cls[col].mean(); sd = cls[col].std()
    ax.bar(m.index, m.values, yerr=sd.values, capsize=4, color=[SOIL_COLOURS[k] for k in [1, 2, 3, 4] if SOIL_NAMES[k] in m.index])
    ax.set_title(title); ax.set_ylabel(unit); ax.tick_params(axis="x", rotation=20)
plt.tight_layout(); plt.savefig(f"{SUB['relationships']}/soil_class_rain_comparison.png", dpi=150); plt.close()

soil_notes = []
soil_notes.append("Soil classes in this grid: " + ", ".join(f"{SOIL_NAMES[c].lower()} {soil_class.loc[c,'grid_cells']} cells"
                  for c in soil_class.index if c > 0) + f"; {int(soil_class.loc[0,'grid_cells']) if 0 in soil_class.index else 0} cells are reserved (mostly coastal).")
ds_sil = sorted([n for n in names if dominant[n] == 4], key=lambda n: -summary.loc[n, 'precip_total_mm'])
soil_notes.append(f"{len(ds_sil)} of {len(names)} districts are dominated by silt loam, including the wettest ones: "
                  + ", ".join(ds_sil[:5]) + ".")
kachchh_sand = comp.loc["Kachchh", "Sand"] + comp.loc["Kachchh", "Sandy loam"] if "Kachchh" in comp.index else np.nan
soil_notes.append(f"Kachchh has about {comp.loc['Kachchh', 'Sandy loam']:.0f}% sandy loam and {comp.loc['Kachchh', 'Sand']:.0f}% sand by grid cell, the most sandy district in the state.")
cls_rain = cells[cells.soil_code > 0].groupby("soil_name").rain_mm_total.mean().round(0)
soil_notes.append("Mean total rain per grid cell by soil class: " + "; ".join(f"{k} {v:.0f} mm" for k, v in cls_rain.items()) +
                  ". Soil and location are mixed together here, so these differences do not show that soil causes the rain pattern.")
soil_notes.append(f"Correlation between district runoff susceptibility and flood-risk index: "
                  f"{quad[['runoff_susceptibility_1to3','risk_index_0_to_1']].corr().iloc[0,1]:.2f} (assumption-based susceptibility).")
with open(f"{SUB['insights']}/soil_insights.md", "w") as f:
    f.write("# Soil type findings\n\n")
    for line in soil_notes:
        f.write(f"- {line}\n")
for line in soil_notes:
    print("-", line)

shutil.make_archive(OUT, "zip", OUT)
print(f"Saved {OUT}.zip  (PNG charts: {len(glob.glob(OUT + '/**/*.png', recursive=True))})")
# from google.colab import files
# files.download(f"{OUT}.zip")