#!/usr/bin/env python3
"""
02_acquire_opendengue.py

Download and prepare multi-country dengue + weather data.

Steps:
  1. Download OpenDengue ZIP from figshare (v5)
  2. Extract and filter for Thailand, Brazil, Vietnam
  3. Download weather via meteostat for Bangkok, Sao Paulo, HCMC
  4. Save per-country CSVs: {country}_dengue_weekly.csv, {country}_raw_weather.csv

Output (in data/):
  thailand_dengue_weekly.csv, thailand_raw_weather.csv
  brazil_dengue_weekly.csv, brazil_raw_weather.csv
  vietnam_dengue_weekly.csv, vietnam_raw_weather.csv
"""

import os
import sys
import zipfile
import warnings
from datetime import datetime

import pandas as pd
import requests

warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
CODE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(CODE_DIR), 'data')
os.makedirs(DATA_DIR, exist_ok=True)

FIGSHARE_ZIP_URL = "https://figshare.com/ndownloader/articles/24259573/versions/5"
LOCAL_ZIP_CACHE = os.path.join(DATA_DIR, 'opendengue_figshare_v5.zip')

# ---------------------------------------------------------------------------
# Country configs
# ---------------------------------------------------------------------------
COUNTRIES = {
    'thailand': {
        'label': 'Thailand (Bangkok)',
        'iso3': 'THA',
        'names': ['Thailand', 'THA'],
        'adm1_filter': None,  # national level
        'weather_lat': 13.7563, 'weather_lon': 100.5018, 'weather_alt': 2,
        'year_range': (2012, 2022),
    },
    'brazil': {
        'label': 'Brazil (Sao Paulo)',
        'iso3': 'BRA',
        'names': ['Brazil', 'BRA', 'Brasil'],
        'adm1_filter': None,  # national level
        'weather_lat': -23.5505, 'weather_lon': -46.6333, 'weather_alt': 760,
        'year_range': (2012, 2022),
    },
    'vietnam': {
        'label': 'Vietnam (HCMC)',
        'iso3': 'VNM',
        'names': ['Vietnam', 'VNM', 'Viet Nam'],
        'adm1_filter': None,  # national level
        'weather_lat': 10.8231, 'weather_lon': 106.6297, 'weather_alt': 19,
        'year_range': (2012, 2022),
    },
}


# =============================================================================
# 1. DOWNLOAD OPENDENGUE ZIP
# =============================================================================

def load_opendengue():
    """Load OpenDengue data from cached ZIP, extracted CSV, or figshare download."""
    print("=" * 60)
    print("1. Loading OpenDengue data")
    print("=" * 60)

    # Try 1: Already-extracted CSV (e.g. Temporal_extract_V1_3.csv)
    import glob
    csv_candidates = glob.glob(os.path.join(DATA_DIR, '*emporal*extract*.csv'))
    csv_candidates += glob.glob(os.path.join(DATA_DIR, '*ational*extract*.csv'))
    if csv_candidates:
        best_csv = max(csv_candidates, key=os.path.getsize)
        print(f"  Found extracted CSV: {best_csv}")
        df = pd.read_csv(best_csv, low_memory=False)
        print(f"  Loaded: {len(df):,} rows x {len(df.columns)} cols")
        _print_country_summary(df)
        return df

    # Try 2: Cached ZIP
    if os.path.exists(LOCAL_ZIP_CACHE):
        size_mb = os.path.getsize(LOCAL_ZIP_CACHE) / 1e6
        print(f"  Using cached ZIP: {LOCAL_ZIP_CACHE} ({size_mb:.1f} MB)")
        return _load_from_zip(LOCAL_ZIP_CACHE)

    # Try 3: Download from figshare
    print(f"  Downloading from: {FIGSHARE_ZIP_URL}")
    try:
        resp = requests.get(FIGSHARE_ZIP_URL, timeout=300, stream=True)
        resp.raise_for_status()

        total = int(resp.headers.get('content-length', 0))
        downloaded = 0

        with open(LOCAL_ZIP_CACHE, 'wb') as f:
            for chunk in resp.iter_content(chunk_size=1024*1024):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total:
                        pct = downloaded / total * 100
                        print(f"  Downloaded {downloaded/1e6:.1f}/{total/1e6:.1f} MB ({pct:.0f}%)", end='\r')

        print(f"\n  Saved: {LOCAL_ZIP_CACHE}")
        return _load_from_zip(LOCAL_ZIP_CACHE)

    except Exception as e:
        raise RuntimeError(
            f"Could not load OpenDengue data: {e}\n"
            f"  Place the extracted CSV (e.g. Temporal_extract_V1_3.csv) in {DATA_DIR}/\n"
            f"  Or download the ZIP from https://figshare.com/articles/dataset/OpenDengue/24259573/5\n"
            f"  and save as {LOCAL_ZIP_CACHE}"
        )


def _load_from_zip(zip_path):
    """Extract and load the best CSV from an OpenDengue ZIP."""
    with zipfile.ZipFile(zip_path, 'r') as zf:
        csv_names = [n for n in zf.namelist() if n.lower().endswith('.csv')]
        if not csv_names:
            raise ValueError("No CSV files found in ZIP!")

        def score(name):
            lower = name.lower()
            s = 0
            if 'temporal' in lower: s += 20
            if 'national' in lower: s += 15
            if 'dengue' in lower: s += 10
            if any(x in lower for x in ('readme', 'metadata', 'codebook')): s -= 30
            s += zf.getinfo(name).file_size / 1e9
            return s

        best = sorted(csv_names, key=score, reverse=True)[0]
        print(f"  Reading from ZIP: {best}")

        with zf.open(best) as f:
            df = pd.read_csv(f, low_memory=False)

    print(f"  Loaded: {len(df):,} rows x {len(df.columns)} cols")
    _print_country_summary(df)
    return df


def _print_country_summary(df):
    """Print country column summary for diagnostics."""
    for col in df.columns:
        if any(x in col.lower() for x in ('country', 'adm0', 'iso')):
            vals = df[col].dropna().unique()
            print(f"  {col}: {len(vals)} unique values")
            sample = sorted(str(v) for v in vals)[:20]
            print(f"    Sample: {sample}")


# =============================================================================
# 3. FILTER COUNTRY DATA
# =============================================================================

def filter_country(df, country_key, cfg):
    """Filter OpenDengue dataframe for a specific country."""

    # Try to find country column
    country_col = None
    for col in df.columns:
        lower = col.lower()
        if 'adm0' in lower and ('name' in lower or 'es' not in lower):
            country_col = col
            break
    if country_col is None:
        for col in df.columns:
            if any(x in col.lower() for x in ('country', 'nation')):
                country_col = col
                break

    # Try ISO code column too
    iso_col = None
    for col in df.columns:
        lower = col.lower()
        if ('iso' in lower or 'adm0' in lower) and ('code' in lower or 'es' in lower or lower.endswith('_iso')):
            iso_col = col
            break
    # Fallback: adm0_es is often the ISO code in OpenDengue
    if iso_col is None:
        for col in df.columns:
            if col.lower() == 'adm0_es':
                iso_col = col
                break

    mask = pd.Series(False, index=df.index)

    if country_col:
        for name in cfg['names']:
            mask |= df[country_col].astype(str).str.lower() == name.lower()

    if iso_col:
        mask |= df[iso_col].astype(str).str.upper() == cfg['iso3']

    # Fallback: partial string match
    if mask.sum() == 0 and country_col:
        mask = df[country_col].astype(str).str.lower().str.contains(
            cfg['names'][0].lower(), na=False)

    df_country = df[mask].copy()
    print(f"    {cfg['label']}: {len(df_country):,} rows found")

    if len(df_country) == 0:
        print(f"    WARNING: No data found for {cfg['label']}")
        if country_col:
            vals = df[country_col].dropna().unique()
            print(f"    Available values in '{country_col}': {sorted(str(v) for v in vals)[:30]}")

    return df_country


def standardise_and_aggregate(df_country, cfg):
    """Standardise dates and aggregate to weekly."""

    # Find date column
    date_col = None
    for col in df_country.columns:
        lower = col.lower()
        if any(x in lower for x in ('start_date', 'calendar_start')):
            date_col = col
            break
    if date_col is None:
        for col in df_country.columns:
            if 'date' in col.lower():
                date_col = col
                break

    # Find case column
    case_col = None
    for col in df_country.columns:
        lower = col.lower()
        if any(x in lower for x in ('dengue_total', 'cases', 'case_count')):
            case_col = col
            break

    if date_col is None or case_col is None:
        print(f"    ERROR: Cannot identify date ({date_col}) or case ({case_col}) column")
        return None

    result = df_country[[date_col, case_col]].copy()
    result.columns = ['date', 'cases']
    result['date'] = pd.to_datetime(result['date'], errors='coerce')
    result['cases'] = pd.to_numeric(result['cases'], errors='coerce').fillna(0).astype(int)
    result = result.dropna(subset=['date'])

    # Filter year range
    result = result[(result['date'].dt.year >= cfg['year_range'][0]) &
                    (result['date'].dt.year <= cfg['year_range'][1])]

    if len(result) == 0:
        print(f"    No data in year range {cfg['year_range']}")
        return None

    result = result.sort_values('date')

    # Check temporal resolution and aggregate to weekly if needed
    if len(result) > 1:
        median_gap = result['date'].diff().median().days
        print(f"    Temporal resolution: median gap = {median_gap} days")

        if median_gap < 5:  # daily or sub-weekly
            print(f"    Aggregating to weekly...")
            result['week_start'] = result['date'].dt.to_period('W-SUN').dt.start_time
            result = result.groupby('week_start').agg({'cases': 'sum'}).reset_index()
            result = result.rename(columns={'week_start': 'date'})
        elif median_gap > 20:  # monthly — need to distribute
            print(f"    Data is monthly — distributing to weekly (uniform)")
            weekly_rows = []
            for _, row in result.iterrows():
                # Distribute monthly cases uniformly across ~4 weeks
                n_weeks = 4
                weekly_cases = int(row['cases'] / n_weeks)
                remainder = int(row['cases']) - weekly_cases * n_weeks
                for w in range(n_weeks):
                    week_date = row['date'] + pd.Timedelta(weeks=w)
                    c = weekly_cases + (1 if w < remainder else 0)
                    weekly_rows.append({'date': week_date, 'cases': c})
            result = pd.DataFrame(weekly_rows)

    print(f"    Final: {len(result)} weeks, {result['date'].min().date()} to {result['date'].max().date()}")
    print(f"    Cases: min={result['cases'].min()}, median={result['cases'].median():.0f}, max={result['cases'].max()}")

    return result


# =============================================================================
# 4. DOWNLOAD WEATHER DATA
# =============================================================================

def download_weather(country_key, cfg):
    """Download daily weather data using meteostat and aggregate to weekly (W-SUN)."""
    try:
        from meteostat import Point, Daily
    except ImportError:
        raise RuntimeError(
            "meteostat not installed. Install with: pip install meteostat"
        )

    try:
        location = Point(cfg['weather_lat'], cfg['weather_lon'], cfg['weather_alt'])
        start = datetime(cfg['year_range'][0], 1, 1)
        end = datetime(cfg['year_range'][1], 12, 31)

        wx = Daily(location, start, end).fetch().reset_index()
        print(f"    meteostat: {len(wx)} daily records retrieved")

        if len(wx) == 0:
            raise RuntimeError(
                f"No weather data returned for {cfg['label']}. "
                "Check coordinates or try a different meteostat station."
            )

        # Aggregate daily to weekly (W-SUN), matching 01_acquire_data.py
        wx['time'] = pd.to_datetime(wx['time'])
        wx = wx.set_index('time')

        weekly = wx.resample('W-SUN').agg({
            'tavg': 'mean',
            'prcp': 'sum'
        }).reset_index()

        weekly.columns = ['date', 'temp_mean', 'rainfall_total']

        # Fill missing values
        weekly['temp_mean'] = weekly['temp_mean'].interpolate().bfill().ffill()
        weekly['rainfall_total'] = weekly['rainfall_total'].fillna(0)

        print(f"    Aggregated to {len(weekly)} weekly records")
        return weekly

    except Exception as e:
        raise RuntimeError(
            f"Weather download failed for {cfg['label']}: {e}"
        )



# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("MULTI-COUNTRY DATA ACQUISITION")
    print("=" * 70)
    print(f"Data directory: {DATA_DIR}\n")

    # --- 1. Load OpenDengue data ---
    df_all = load_opendengue()

    # --- 3. Process each country ---
    print("\n" + "=" * 60)
    print("3. Extracting country-specific data")
    print("=" * 60)

    for country_key, cfg in COUNTRIES.items():
        print(f"\n  --- {cfg['label']} ---")

        # Filter
        df_country = filter_country(df_all, country_key, cfg)
        if len(df_country) == 0:
            continue

        # Standardise and aggregate to weekly
        df_weekly = standardise_and_aggregate(df_country, cfg)
        if df_weekly is None or len(df_weekly) < 52:
            print(f"    SKIP: insufficient data ({len(df_weekly) if df_weekly is not None else 0} weeks)")
            continue

        # Save dengue data
        dengue_file = os.path.join(DATA_DIR, f'{country_key}_dengue_weekly.csv')
        df_weekly.to_csv(dengue_file, index=False)
        print(f"    Saved: {dengue_file}")

    # --- 4. Download weather ---
    print("\n" + "=" * 60)
    print("4. Downloading weather data")
    print("=" * 60)

    for country_key, cfg in COUNTRIES.items():
        print(f"\n  --- {cfg['label']} ---")

        # Check if dengue data exists first
        dengue_file = os.path.join(DATA_DIR, f'{country_key}_dengue_weekly.csv')
        if not os.path.exists(dengue_file):
            print(f"    SKIP: no dengue data")
            continue

        weather = download_weather(country_key, cfg)

        weather_file = os.path.join(DATA_DIR, f'{country_key}_raw_weather.csv')
        weather.to_csv(weather_file, index=False)
        print(f"    Saved: {weather_file}")

    # --- Summary ---
    print("\n" + "=" * 70)
    print("ACQUISITION COMPLETE")
    print("=" * 70)
    for country_key in COUNTRIES:
        dengue_f = os.path.join(DATA_DIR, f'{country_key}_dengue_weekly.csv')
        weather_f = os.path.join(DATA_DIR, f'{country_key}_raw_weather.csv')
        d_ok = "✓" if os.path.exists(dengue_f) else "✗"
        w_ok = "✓" if os.path.exists(weather_f) else "✗"
        print(f"  {country_key}: dengue {d_ok}  weather {w_ok}")

    print("\nNext: Rscript code/51_prepare_multicountry_data.R")
