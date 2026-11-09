"""Template for Fig. 4-5 PISCES-QUOTA model map figures.

This script is not run by `run_all_figures.py` because the NetCDF model output
files are not included in the uploaded experimental Excel workbook.

Expected workflow:
1. Place model output files in model_inputs/ or update the paths below.
2. Confirm the nitrogen fixation variable name and units.
3. Run this script to export PDF/SVG/TIFF/PNG at 600 dpi.

The plotting code keeps labels and annotations editable in PDF/SVG. For final
journal submission, check coastlines and map projections in Illustrator/Inkscape.
"""
from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

from ncc_style import mm_to_in, panel_label, save_figure, set_ncc_style

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
except Exception:  # pragma: no cover - optional dependency
    ccrs = None
    cfeature = None

# Update these paths and variable names after downloading the model data.
MODEL_FILES = {
    "control_historical": "model_inputs/control_historical_1996_2005.nc",
    "control_rcp85": "model_inputs/control_rcp85_2091_2100.nc",
    "short_term_rcp85": "model_inputs/short_term_adaptation_rcp85_2091_2100.nc",
    "long_term_rcp85": "model_inputs/long_term_adaptation_rcp85_2091_2100.nc",
}
NFIX_VAR = "nfix"       # Change to the actual variable name if needed.
LAT_NAME = "lat"
LON_NAME = "lon"
TIME_NAME = "time"
AREA_VAR = "area"       # Optional grid-cell area variable, if present.


def _guess_coord(ds: xr.Dataset, names: Sequence[str]) -> str:
    for name in names:
        if name in ds.coords or name in ds.variables:
            return name
    raise KeyError(f"None of these coordinates/variables were found: {names}")


def load_mean_nfix(path: str | Path, var: str = NFIX_VAR) -> xr.DataArray:
    """Load time-mean depth-integrated N2 fixation.

    If the variable has a depth dimension, this function sums over dimensions
    containing 'depth', 'lev' or 'z'. Adjust this if the uploaded model files
    already contain depth-integrated fields.
    """
    ds = xr.open_dataset(path)
    if var not in ds:
        raise KeyError(f"Variable {var!r} not found in {path}. Available: {list(ds.data_vars)}")
    da = ds[var]
    for dim in list(da.dims):
        if dim.lower() in {"depth", "deptht", "lev", "level", "z"}:
            da = da.sum(dim=dim, skipna=True)
    if TIME_NAME in da.dims:
        da = da.mean(dim=TIME_NAME, skipna=True)
    return da


def global_integral(da: xr.DataArray, area: xr.DataArray | None = None) -> float:
    """Return global integral. If area is unavailable, cosine-lat weights are used."""
    lat = _guess_coord(da.to_dataset(name="tmp"), [LAT_NAME, "latitude", "nav_lat", "y"])
    if area is not None:
        return float((da * area).sum(skipna=True))
    weights = np.cos(np.deg2rad(da[lat]))
    return float(da.weighted(weights).mean(skipna=True))


def _plot_map(ax, da: xr.DataArray, title: str, cmap: str, vmin=None, vmax=None):
    lon = _guess_coord(da.to_dataset(name="tmp"), [LON_NAME, "longitude", "nav_lon", "x"])
    lat = _guess_coord(da.to_dataset(name="tmp"), [LAT_NAME, "latitude", "nav_lat", "y"])
    if ccrs is not None:
        m = ax.pcolormesh(da[lon], da[lat], da, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
        ax.coastlines(linewidth=0.3)
        ax.add_feature(cfeature.LAND, facecolor="#D9D9D9", edgecolor="none", zorder=3)
        ax.set_global()
    else:
        m = ax.pcolormesh(da[lon], da[lat], da, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
    ax.set_title(title, pad=2)
    return m


def build_fig4(out_dir: str | Path = "outputs/figures", dpi: int = 600) -> list[Path]:
    set_ncc_style()
    ctrl_hist = load_mean_nfix(MODEL_FILES["control_historical"])
    ctrl_future = load_mean_nfix(MODEL_FILES["control_rcp85"])
    anomaly = ctrl_future - ctrl_hist

    projection = ccrs.Robinson() if ccrs is not None else None
    fig, axes = plt.subplots(1, 3, figsize=(mm_to_in(180), mm_to_in(65)), subplot_kw={"projection": projection} if projection else None)
    m0 = _plot_map(axes[0], ctrl_hist, "Historical control", cmap="viridis")
    m1 = _plot_map(axes[1], ctrl_future, "RCP8.5 control", cmap="viridis")
    vmax = float(np.nanmax(np.abs(anomaly.values)))
    m2 = _plot_map(axes[2], anomaly, "RCP8.5 - historical", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    for ax, label in zip(axes, ["a", "b", "c"]):
        panel_label(ax, label, x=0.02, y=1.02)
    fig.colorbar(m1, ax=axes[:2], shrink=0.65, pad=0.02, label="N₂ fixation")
    fig.colorbar(m2, ax=axes[2], shrink=0.65, pad=0.02, label="Anomaly")
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.08, top=0.90, wspace=0.10)
    return save_figure(fig, out_dir, "Fig4_model_control_maps", dpi=dpi)


def build_fig5(out_dir: str | Path = "outputs/figures", dpi: int = 600) -> list[Path]:
    set_ncc_style()
    ctrl_future = load_mean_nfix(MODEL_FILES["control_rcp85"])
    short_future = load_mean_nfix(MODEL_FILES["short_term_rcp85"])
    long_future = load_mean_nfix(MODEL_FILES["long_term_rcp85"])
    maps = [short_future - ctrl_future, long_future - ctrl_future, short_future - long_future]
    titles = ["Short-term - control", "Long-term - control", "Short-term - long-term"]
    vmax = max(float(np.nanmax(np.abs(m.values))) for m in maps)

    projection = ccrs.Robinson() if ccrs is not None else None
    fig, axes = plt.subplots(1, 3, figsize=(mm_to_in(180), mm_to_in(65)), subplot_kw={"projection": projection} if projection else None)
    last = None
    for ax, da, title, label in zip(axes, maps, titles, ["a", "b", "c"]):
        last = _plot_map(ax, da, title, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        panel_label(ax, label, x=0.02, y=1.02)
    fig.colorbar(last, ax=axes, shrink=0.65, pad=0.02, label="N₂ fixation anomaly")
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.08, top=0.90, wspace=0.10)
    return save_figure(fig, out_dir, "Fig5_model_adaptation_anomalies", dpi=dpi)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="outputs/figures")
    parser.add_argument("--dpi", default=600, type=int)
    parser.add_argument("--figure", choices=["4", "5", "both"], default="both")
    args = parser.parse_args()
    paths = []
    if args.figure in {"4", "both"}:
        paths += build_fig4(args.out, dpi=args.dpi)
    if args.figure in {"5", "both"}:
        paths += build_fig5(args.out, dpi=args.dpi)
    for p in paths:
        print(p)
