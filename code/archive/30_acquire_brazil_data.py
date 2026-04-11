#!/usr/bin/env python3
"""
30_acquire_brazil_data.py

Download and process data for dengue Rt estimation in Brazil:
1. Weekly dengue case counts + weather from InfoDengue API (São Paulo)
2. Meteorological data from Meteostat as backup

Output: data/brazil_raw_cases.csv and data/brazil_raw_weather.csv
"""

import os
import io
import time
import pandas as pd
import numpy as np
from datetime import datetime
import requests
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Create data directory if it doesn't exist
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')
os.makedirs(DATA_DIR, exist_ok=True)

# =============================================================================
# CONFIGURATION
# =============================================================================

# São Paulo municipality geocode (IBGE)
GEOCODE = 3550308
CITY_NAME = "São Paulo"

# Meteostat backup: São Paulo Congonhas Airport (WMO 83780)
METEOSTAT_WMO = "83780"
METEOSTAT_LAT = -23.6267
METEOSTAT_LON = -46.6544

# Study period
YEAR_START = 2012
YEAR_END = 2022

# API retry settings
MAX_RETRIES = 3
RETRY_DELAY = 5  # seconds


# =============================================================================
# 1. INFODENGUE API: DENGUE CASES + WEATHER
# =============================================================================

def download_infodengue_data():
    """
    Download weekly dengue data from InfoDengue alertcity API.

    The API returns weekly epidemiological data including case counts and
    climate variables (temperature, humidity) from weather stations.

    API docs: https://info.dengue.mat.br/api/
    """
    print("=" * 60)
    print(f"Downloading dengue data from InfoDengue ({CITY_NAME})...")
    print("=" * 60)

    # InfoDengue alertcity endpoint
    # We query year by year to avoid timeouts on large requests
    all_data = []

    for year in range(YEAR_START, YEAR_END + 1):
        url = "https://info.dengue.mat.br/api/alertcity"
        params = {
            "geocode": GEOCODE,
            "disease": "dengue",
            "format": "csv",
            "ew_start": 1,
            "ew_end": 52,
            "ey_start": year,
            "ey_end": year
        }

        success = False
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                print(f"  Fetching {year} (attempt {attempt}/{MAX_RETRIES})...")
                response = requests.get(url, params=params, timeout=60)
                response.raise_for_status()

                # Parse CSV response
                df_year = pd.read_csv(io.StringIO(response.text))

                if len(df_year) > 0:
                    all_data.append(df_year)
                    print(f"    Got {len(df_year)} weeks for {year}")
                    success = True
                    break
                else:
                    print(f"    Empty response for {year}")
                    success = True  # Empty is not an error
                    break

            except requests.exceptions.RequestException as e:
                print(f"    Error: {e}")
                if attempt < MAX_RETRIES:
                    print(f"    Retrying in {RETRY_DELAY} seconds...")
                    time.sleep(RETRY_DELAY)
                else:
                    print(f"    Failed after {MAX_RETRIES} attempts for {year}")

        if not success:
            print(f"  WARNING: Could not download data for {year}")

        # Brief pause between requests to be polite to the API
        time.sleep(1)

    if not all_data:
        print("  ERROR: No data downloaded from InfoDengue!")
        print("  Check API availability at: https://info.dengue.mat.br/api/")
        return None

    # Combine all years
    df = pd.concat(all_data, ignore_index=True)
    print(f"\n  Total records downloaded: {len(df)}")
    print(f"  Columns: {df.columns.tolist()}")

    return df


def process_infodengue_data(df):
    """
    Process raw InfoDengue data into weekly case counts and weather.

    InfoDengue columns of interest:
      - data_iniSE: epi week start date
      - casos: reported dengue cases
      - tempmin, tempmed, tempmax: temperature (°C)
      - umidmin, umidmed, umidmax: humidity (%)
    """
    print("\nProcessing InfoDengue data...")

    if df is None:
        print("  No data available.")
        return None, None

    # Parse date column
    df['date'] = pd.to_datetime(df['data_iniSE'])
    df = df.sort_values('date').reset_index(drop=True)

    # Drop duplicate weeks (can happen at year boundaries)
    df = df.drop_duplicates(subset=['date'], keep='first')

    # Filter to study period
    df = df[(df['date'].dt.year >= YEAR_START) & (df['date'].dt.year <= YEAR_END)]

    # --- Cases ---
    cases = df[['date', 'casos']].copy()
    cases.columns = ['date', 'cases']
    cases['cases'] = pd.to_numeric(cases['cases'], errors='coerce').fillna(0).astype(int)

    # --- Weather from InfoDengue ---
    # InfoDengue provides weather data from nearby stations
    weather_cols = {
        'date': 'date',
        'tempmin': 'temp_min',
        'tempmed': 'temp_mean',
        'tempmax': 'temp_max',
        'umidmin': 'humid_min',
        'umidmed': 'humid_mean',
        'umidmax': 'humid_max'
    }

    weather = pd.DataFrame()
    weather['date'] = df['date']

    for src_col, dst_col in weather_cols.items():
        if src_col == 'date':
            continue
        if src_col in df.columns:
            weather[dst_col] = pd.to_numeric(df[src_col], errors='coerce')
        else:
            print(f"  WARNING: Column '{src_col}' not found in InfoDengue data")

    # Interpolate missing weather values
    for col in weather.columns:
        if col != 'date':
            weather[col] = weather[col].interpolate(method='linear')

    print(f"  Cases: {len(cases)} weeks, range {cases['cases'].min()} - {cases['cases'].max()}")
    print(f"  Weather: {len(weather)} weeks")

    if 'temp_mean' in weather.columns:
        print(f"  Temperature range: {weather['temp_mean'].min():.1f} - {weather['temp_mean'].max():.1f}°C")
    if 'humid_mean' in weather.columns:
        print(f"  Humidity range: {weather['humid_mean'].min():.1f} - {weather['humid_mean'].max():.1f}%")

    return cases, weather


# =============================================================================
# 2. METEOSTAT WEATHER DATA (BACKUP)
# =============================================================================

def download_meteostat_weather():
    """
    Download daily weather data from Meteostat as backup.
    São Paulo Congonhas Airport (WMO: 83780)
    """
    print("\n" + "=" * 60)
    print("Downloading backup weather data from Meteostat...")
    print("=" * 60)

    try:
        from meteostat import Point, Daily

        sp = Point(METEOSTAT_LAT, METEOSTAT_LON, 802)  # elevation ~802m
        start = datetime(YEAR_START, 1, 1)
        end = datetime(YEAR_END, 12, 31)

        data = Daily(sp, start, end)
        df = data.fetch()

        if len(df) > 0:
            print(f"  Downloaded {len(df)} daily records from Meteostat")
            df = df.reset_index()
            return df
        else:
            print("  No data returned from Meteostat")
            return None

    except ImportError:
        print("  Meteostat not installed (pip install meteostat)")
        print("  Skipping backup weather download.")
        return None
    except Exception as e:
        print(f"  Error: {e}")
        return None


def process_meteostat_weather(df):
    """
    Aggregate daily Meteostat weather to weekly means/totals.
    """
    if df is None:
        return None

    print("\nProcessing Meteostat weather to weekly aggregates...")

    df['time'] = pd.to_datetime(df['time'])
    df = df.set_index('time')

    # Weekly aggregates: mean temperature, total precipitation
    agg_dict = {}
    if 'tavg' in df.columns:
        agg_dict['tavg'] = 'mean'
    if 'prcp' in df.columns:
        agg_dict['prcp'] = 'sum'

    if not agg_dict:
        print("  No usable columns in Meteostat data")
        return None

    weekly = df.resample('W-SUN').agg(agg_dict).reset_index()

    rename_map = {'time': 'date'}
    if 'tavg' in weekly.columns:
        rename_map['tavg'] = 'temp_mean'
    if 'prcp' in weekly.columns:
        rename_map['prcp'] = 'rainfall_total'

    weekly = weekly.rename(columns=rename_map)

    # Interpolate missing
    for col in weekly.columns:
        if col != 'date':
            weekly[col] = weekly[col].interpolate(method='linear')

    print(f"  Created {len(weekly)} weekly records")
    return weekly


# =============================================================================
# 3. COMBINE AND SAVE
# =============================================================================

def create_final_weather(infodengue_weather, meteostat_weather):
    """
    Create final weather dataset. Use InfoDengue data as primary source,
    with Meteostat as backup for missing values.

    For the Stan model we need:
      - temp_mean: weekly mean temperature (°C)
      - rainfall_total: weekly total rainfall (mm)

    InfoDengue provides temperature but not rainfall directly.
    We use humidity as an alternative climate covariate, or Meteostat
    rainfall if available.
    """
    print("\nCreating final weather dataset...")

    if infodengue_weather is None and meteostat_weather is None:
        print("  ERROR: No weather data available!")
        return None

    # Start with InfoDengue weather
    if infodengue_weather is not None:
        weather = infodengue_weather.copy()

        # If Meteostat has rainfall and InfoDengue doesn't, merge it in
        if meteostat_weather is not None and 'rainfall_total' in meteostat_weather.columns:
            if 'rainfall_total' not in weather.columns:
                weather = weather.merge(
                    meteostat_weather[['date', 'rainfall_total']],
                    on='date', how='left'
                )
                weather['rainfall_total'] = weather['rainfall_total'].interpolate(method='linear')
                print("  Added Meteostat rainfall to InfoDengue weather")

        # If we still don't have rainfall, use humidity as alternative
        if 'rainfall_total' not in weather.columns:
            if 'humid_mean' in weather.columns:
                print("  NOTE: No rainfall data available. Using humidity as climate covariate.")
                weather['rainfall_total'] = weather['humid_mean']
            else:
                print("  WARNING: No rainfall or humidity data. Setting rainfall to 0.")
                weather['rainfall_total'] = 0.0

    else:
        weather = meteostat_weather.copy()

    # Ensure required columns exist
    if 'temp_mean' not in weather.columns:
        print("  WARNING: No temperature data!")

    if 'rainfall_total' not in weather.columns:
        print("  WARNING: No rainfall data!")

    print(f"  Final weather: {len(weather)} weeks")
    return weather


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 70)
    print(f"DENGUE RT ESTIMATION - BRAZIL DATA ACQUISITION ({CITY_NAME})")
    print("=" * 70)

    # 1. Download from InfoDengue
    infodengue_raw = download_infodengue_data()
    cases, infodengue_weather = process_infodengue_data(infodengue_raw)

    # 2. Download backup weather from Meteostat
    meteostat_raw = download_meteostat_weather()
    meteostat_weather = process_meteostat_weather(meteostat_raw)

    # 3. Create final weather dataset
    weather = create_final_weather(infodengue_weather, meteostat_weather)

    # 4. Save outputs
    if cases is not None:
        cases.to_csv(os.path.join(DATA_DIR, 'brazil_raw_cases.csv'), index=False)
        print(f"\nSaved: data/brazil_raw_cases.csv ({len(cases)} rows)")
    else:
        print("\nERROR: No case data to save!")

    if weather is not None:
        weather.to_csv(os.path.join(DATA_DIR, 'brazil_raw_weather.csv'), index=False)
        print(f"Saved: data/brazil_raw_weather.csv ({len(weather)} rows)")
    else:
        print("WARNING: No weather data to save!")

    # 5. Summary statistics
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    if cases is not None:
        print(f"\n  City: {CITY_NAME} (geocode: {GEOCODE})")
        print(f"  Period: {cases['date'].min()} to {cases['date'].max()}")
        print(f"  Total weeks: {len(cases)}")
        print(f"\n  Case counts:")
        print(f"    Min:    {cases['cases'].min()}")
        print(f"    Median: {cases['cases'].median():.0f}")
        print(f"    Mean:   {cases['cases'].mean():.1f}")
        print(f"    Max:    {cases['cases'].max()}")
        print(f"    Total:  {cases['cases'].sum()}")

        # Annual breakdown
        cases['year'] = pd.to_datetime(cases['date']).dt.year
        annual = cases.groupby('year')['cases'].sum()
        print(f"\n  Annual totals:")
        for year, total in annual.items():
            print(f"    {year}: {total:,.0f}")

    if weather is not None:
        print(f"\n  Weather:")
        if 'temp_mean' in weather.columns:
            print(f"    Temperature: {weather['temp_mean'].mean():.1f}°C "
                  f"(range {weather['temp_mean'].min():.1f} - {weather['temp_mean'].max():.1f})")
        if 'rainfall_total' in weather.columns:
            print(f"    Rainfall/humidity: mean {weather['rainfall_total'].mean():.1f}, "
                  f"range {weather['rainfall_total'].min():.1f} - {weather['rainfall_total'].max():.1f}")

    print("\n" + "=" * 70)
    print("DATA ACQUISITION COMPLETE")
    print("=" * 70)
    print("\nNext step: Run 31_prepare_brazil_data.R to format data for Stan")
