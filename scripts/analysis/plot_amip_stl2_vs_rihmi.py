#!/usr/bin/env python3
"""AMIP counterpart of the coupled 09C/10A/10B soil-temperature figure.

Adapted from a270270's plot_oifs_stl2_vs_rihmi.py so the two figures are directly
comparable: identical observations, identical QC=0 / >=150-month station screening,
identical nearest-land-cell matching, identical stl2 (7-28 cm, centre 17.5 cm) vs
observed 20 cm, identical metrics.  ONLY the experiment list changes -- every model
line here is AMIP, where the P series exists.  There is no coupled P run yet (that
would be 10C), and putting AMIP runs on coupled axes would compare two different
model configurations.

    N1  K1 base, snow scheme OFF   -- the baseline, AMIP analogue of 10A
    N2  ECE_SNOW_SCF=1 (tanh)      -- the broken winter, AMIP analogue of 10B
    P3  ECE_SNOW_SCF=3, SCALE=1    -- the observationally fitted curve
    P4  ECE_SNOW_SCF=3, SCALE=3    -- same curve, d_c x3 (uncalibrated bracket)

THREE DEVIATIONS FROM THE ORIGINAL, all forced and all stated:
  1. OBS_FILE is RIHMI-WDC_tpg.nc.  The original points at
     RIHMI-WDC_soil_temperature_v3_1963-2024.nc, which no longer exists after the
     2026-08-06 re-download; tpg.nc carries the same tsoil AND tsoil_qc, so the
     screening is unchanged.
  2. AMIP runs do not write lsm, so the land mask comes from the shared LSMF in
     runs.py -- the same remapped TCO95 grid every evaluator in this repo uses.
  3. Monthly files here are named atm_remapped_1m_<var>_1m_<y>-<y>.nc.

PERIOD CAVEAT, identical in kind to the coupled figure: the model is pre-industrial
and the observations are 1991-2020, so the whole model family reads cold by roughly
the industrial warming of Siberian soil.  That offset is common to every line, so
the RANKING and the SPREAD are the meaningful content, not the absolute bias.
"""
from __future__ import annotations

import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'

from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree

from runs import RT, LSMF

OBS_FILE = Path('/work/ab0246/a270092/obs/RIHMI-WDC/data/RIHMI-WDC_tpg.nc')
PLOTS = Path('/work/ab0246/a270092/postprocessing/'
             'investigation_awiesm3_high_lat_cold_bias_global_toa_positive/report/plots')
MONTHS = np.arange(1, 13)
MONTH_LABELS = np.array(list('JFMAMJJASOND'))


@dataclass(frozen=True)
class Experiment:
    key: str
    label: str
    first_year: int
    last_year: int
    run: str
    color: str

    @property
    def years_label(self) -> str:
        return f'{self.first_year}–{self.last_year}'


EXPERIMENTS = (
    Experiment('N1', 'N1 AMIP baseline (scheme off)', 1908, 1917,
               'amip_N1_snowdiag', '#5e3c99'),
    Experiment('N2', 'N2 AMIP + tanh depletion', 1908, 1917,
               'amip_N2_snowdiag_scf', '#1b9e77'),
    Experiment('P3', 'P3 AMIP + fitted curve', 1908, 1917,
               'amip_P3_scffit', '#e66101'),
    Experiment('P4', 'P4 AMIP + fitted curve, d_c x3', 1908, 1917,
               'amip_P4_scffit_x3', '#0571b0'),
)


def unit_xyz(lon, lat):
    lon_r = np.deg2rad((np.asarray(lon) + 180.0) % 360.0 - 180.0)
    lat_r = np.deg2rad(np.asarray(lat))
    return np.column_stack([np.cos(lat_r) * np.cos(lon_r),
                            np.cos(lat_r) * np.sin(lon_r), np.sin(lat_r)])


def read_observations() -> pd.DataFrame:
    """QC-screened monthly 1991-2020 station climatology at 20 cm."""
    ds = xr.open_dataset(OBS_FILE).sel(time=slice('1991-01-01', '2020-12-31'))
    temperature = ds['tsoil'].sel(depth=0.20)
    good = temperature.where(ds['tsoil_qc'].sel(depth=0.20) == 0)
    climatology = good.groupby('time.month').mean('time')
    count = good.groupby('time.month').count('time')
    selected = np.flatnonzero((count >= 150).all('month').values)
    values = climatology.isel(station=selected).transpose('station', 'month').values
    frame = pd.DataFrame(values, columns=MONTHS)
    frame.insert(0, 'lon', ds['lon'].isel(station=selected).values)
    frame.insert(1, 'lat', ds['lat'].isel(station=selected).values)
    frame.insert(2, 'station_id', ds['station_id'].isel(station=selected).values)
    ds.close()
    return frame


def station_grid_indices(station_lon, station_lat):
    """Nearest remapped grid cell with lsm >= 0.5 (shared mask; AMIP writes none)."""
    with xr.open_dataset(LSMF) as ds:
        lon = ds['lon'].values
        lat = ds['lat'].values
        mask = ds['lsm'].isel(time_counter=0).values
    grid_lon, grid_lat = np.meshgrid(lon, lat)
    land_y, land_x = np.where(mask >= 0.5)
    tree = cKDTree(unit_xyz(grid_lon[land_y, land_x], grid_lat[land_y, land_x]))
    chord, nearest = tree.query(unit_xyz(station_lon, station_lat), k=1)
    distance = np.rad2deg(2.0 * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0)))
    if np.max(distance) > 1.5:
        raise RuntimeError(f'Station-to-land distance reaches {np.max(distance):.2f} deg')
    return land_y[nearest], land_x[nearest], distance


def read_model_at_stations(experiment, iy, ix) -> np.ndarray:
    annual = []
    for year in range(experiment.first_year, experiment.last_year + 1):
        base = Path(RT) / experiment.run / 'outdata/oifs'
        path = base / f'atm_remapped_1m_stl2_1m_{year}-{year}.nc'
        if not path.exists():
            path = base / f'atm_remapped_1m_stl2_{year}-{year}.nc'
        if not path.exists():
            raise FileNotFoundError(path)
        with xr.open_dataset(path, decode_times=False) as ds:
            values = ds['stl2'].values[:, iy, ix] - 273.15
        if values.shape != (12, len(iy)):
            raise RuntimeError(f'Unexpected stl2 shape {values.shape}: {path}')
        annual.append(values)
    return np.mean(np.stack(annual, axis=0), axis=0).T


def metrics(bias) -> dict:
    return {'mean_bias_C': float(np.mean(bias)),
            'rmse_C': float(np.sqrt(np.mean(bias ** 2))),
            'mae_C': float(np.mean(np.abs(bias))),
            'DJF_bias_C': float(np.mean(bias[:, [11, 0, 1]])),
            'JJA_bias_C': float(np.mean(bias[:, [5, 6, 7]]))}


def plot_cycles(observations, models, results, station_count) -> Path:
    fig, axes = plt.subplots(2, 1, figsize=(12.7, 9.6), sharex=True)
    axes[0].plot(MONTHS, np.mean(observations, axis=0), color='black', lw=3,
                 marker='o', label='RIHMI-WDC 20 cm')
    for experiment in EXPERIMENTS:
        model = models[experiment.key]
        row = results.loc[experiment.key]
        axes[0].plot(MONTHS, np.mean(model, axis=0), color=experiment.color, lw=2.5,
                     marker='o',
                     label=(f'{experiment.label} '
                            f'(bias {row.mean_bias_C:+.2f} °C; '
                            f'RMSE {row.rmse_C:.2f} °C)'))
        axes[1].plot(MONTHS, np.mean(model - observations, axis=0),
                     color=experiment.color, lw=2.5, marker='o', label=experiment.label)
    axes[0].set_ylabel('Soil temperature [°C]')
    axes[0].legend(frameon=False, fontsize=9.5, ncol=2)
    axes[0].grid(alpha=0.25)
    axes[1].axhline(0, color='#555555', lw=1.2)
    axes[1].set_ylabel('OpenIFS − RIHMI-WDC [°C]')
    axes[1].set_xlabel('Month')
    axes[1].set_xticks(MONTHS, MONTH_LABELS)
    axes[1].grid(alpha=0.25)
    fig.suptitle('AMIP OpenIFS/HTESSEL Subsurface Soil Temperature versus RIHMI-WDC',
                 fontsize=17, fontweight='bold', y=0.985)
    fig.text(0.5, 0.94,
             f'OpenIFS stl2: 7–28 cm layer (centre 17.5 cm) · RIHMI-WDC: '
             f'observed 20 cm · equal-weight mean across {station_count} stations',
             ha='center', va='top', fontsize=10.5, color='#444444')
    fig.text(0.5, 0.91,
             'Observations: QC=0, 1991–2020 · AMIP models: final 10 years '
             '(1908–1917, pre-industrial) · climatological, not date-matched',
             ha='center', va='top', fontsize=9.6, color='#8a4b08', fontstyle='italic')
    fig.subplots_adjust(top=0.85, bottom=0.08, left=0.10, right=0.98, hspace=0.10)
    output = PLOTS / 'amip_stl2_vs_RIHMI_20cm_seasonal_cycle_N1_N2_P3_P4.png'
    fig.savefig(output, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return output


def main() -> None:
    obs = read_observations()
    observations = obs[MONTHS].to_numpy(dtype=float)
    iy, ix, distance = station_grid_indices(obs['lon'].to_numpy(), obs['lat'].to_numpy())
    print(f'{len(obs)} stations, max match distance {np.max(distance):.2f} deg', flush=True)
    models, rows = {}, []
    for experiment in EXPERIMENTS:
        print(f'Reading {experiment.label} ({experiment.years_label})', flush=True)
        model = read_model_at_stations(experiment, iy, ix)
        models[experiment.key] = model
        row = metrics(model - observations)
        row.update(experiment=experiment.key, model_years=experiment.years_label,
                   stations=len(obs))
        rows.append(row)
    results = pd.DataFrame(rows).set_index('experiment')
    results.to_csv(PLOTS / 'amip_stl2_vs_RIHMI_20cm_metrics_N1_N2_P3_P4.csv',
                   float_format='%.4f')
    out = plot_cycles(observations, models, results, len(obs))
    print(results.to_string(float_format=lambda v: f'{v:.3f}'), flush=True)
    print(f'Saved: {out}', flush=True)


if __name__ == '__main__':
    main()
