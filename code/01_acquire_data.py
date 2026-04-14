#!/usr/bin/env python3
"""
01_acquire_data.py

Download and process data for dengue Rt estimation:
1. Weekly dengue case counts from data.gov.sg (MOH Weekly Infectious Diseases Bulletin)
2. Meteorological data from Meteostat (Singapore Changi Airport)

Output: data/raw_*.csv files ready for merging in 03_prepare_model_data.R
"""

import os
import pandas as pd
from datetime import datetime
import requests
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Create data directory if it doesn't exist
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')
os.makedirs(DATA_DIR, exist_ok=True)

# =============================================================================
# 1. DENGUE CASE DATA FROM DATA.GOV.SG
# =============================================================================

def download_dengue_data():
    """
    Download weekly dengue case counts from data.gov.sg
    Dataset: Weekly Infectious Diseases Bulletin

    If API fails, checks for manually downloaded file at:
    data/WeeklyInfectiousDiseaseBulletinCases.csv

    Manual download: https://data.gov.sg/datasets/d_ca168b2cb763640d72c4600a68f9909e/view
    """
    print("=" * 60)
    print("Downloading dengue case data from data.gov.sg...")
    print("=" * 60)

    # Check for manually downloaded file first
    manual_file = os.path.join(DATA_DIR, 'WeeklyInfectiousDiseaseBulletinCases.csv')
    if os.path.exists(manual_file):
        print(f"  Found manually downloaded file: {manual_file}")
        df = pd.read_csv(manual_file)
        print(f"  Loaded {len(df)} records from local file")
        return df

    # data.gov.sg API endpoint for the weekly infectious diseases bulletin
    # Dataset ID: d_ca168b2cb763640d72c4600a68f9909e
    url = "https://data.gov.sg/api/action/datastore_search"

    # The resource ID for the weekly infectious diseases bulletin
    # This may need to be updated if the API structure changes
    resource_id = "d_ca168b2cb763640d72c4600a68f9909e"

    all_records = []
    offset = 0
    limit = 1000

    while True:
        params = {
            "resource_id": resource_id,
            "limit": limit,
            "offset": offset
        }

        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            if data.get("success"):
                records = data.get("result", {}).get("records", [])
                if not records:
                    break
                all_records.extend(records)
                offset += limit
                print(f"  Downloaded {len(all_records)} records...")

                # Check if we've got all records
                total = data.get("result", {}).get("total", 0)
                if len(all_records) >= total:
                    break
            else:
                print(f"  API returned error: {data}")
                break

        except requests.exceptions.RequestException as e:
            print(f"  Error downloading data: {e}")
            print("  Attempting alternative download method...")
            break

    if all_records:
        df = pd.DataFrame(all_records)
        print(f"  Successfully downloaded {len(df)} records")
        return df

    # Alternative: Try direct CSV download
    print("  Trying direct CSV download...")
    csv_url = "https://data.gov.sg/dataset/d_ca168b2cb763640d72c4600a68f9909e/resource/d_ca168b2cb763640d72c4600a68f9909e/download/weekly-infectious-bulletin.csv"

    try:
        df = pd.read_csv(csv_url)
        print(f"  Successfully downloaded {len(df)} records via CSV")
        return df
    except Exception as e:
        print(f"  CSV download also failed: {e}")
        print("\n" + "!" * 60)
        print("MANUAL DOWNLOAD REQUIRED")
        print("!" * 60)
        print("The data.gov.sg API is rate-limiting requests.")
        print("Please manually download the data:")
        print("  1. Go to: https://data.gov.sg/datasets/d_ca168b2cb763640d72c4600a68f9909e/view")
        print("  2. Click 'Download' to get the CSV file")
        print("  3. Save it as: data/WeeklyInfectiousDiseaseBulletinCases.csv")
        print("  4. Re-run this script")
        print("!" * 60 + "\n")
        return None


def process_dengue_data(df):
    """
    Process raw dengue data into weekly case counts
    """
    print("\nProcessing dengue data...")

    if df is None:
        raise RuntimeError(
            "Dengue data download failed. Please manually download from:\n"
            "  https://data.gov.sg/datasets/d_ca168b2cb763640d72c4600a68f9909e/view\n"
            "Save as: data/WeeklyInfectiousDiseaseBulletinCases.csv"
        )

    # Print column names to understand structure
    print(f"  Columns: {df.columns.tolist()}")

    # Filter for dengue cases (dengue fever + dengue haemorrhagic fever)
    # Find the disease column
    disease_col = None
    for col in df.columns:
        if 'disease' in col.lower():
            disease_col = col
            break

    if disease_col:
        df_dengue = df[df[disease_col].str.lower().str.contains('dengue', na=False)]
    else:
        print("  Could not find disease column, using all data")
        df_dengue = df

    print(f"  Found {len(df_dengue)} dengue records")

    # Find the date/week column
    date_col = None
    for col in df.columns:
        if any(x in col.lower() for x in ['epi_week', 'week', 'date']):
            date_col = col
            break

    # Find the case count column
    case_col = None
    for col in df.columns:
        if any(x in col.lower() for x in ['no_of_cases', 'cases', 'count']):
            case_col = col
            break

    if date_col and case_col:
        # Aggregate by week
        df_dengue[case_col] = pd.to_numeric(df_dengue[case_col], errors='coerce')
        weekly = df_dengue.groupby(date_col)[case_col].sum().reset_index()
        weekly.columns = ['epi_week', 'cases']

        # Parse epi_week to get date (format: YYYY-Www or similar)
        weekly['date'] = weekly['epi_week'].apply(parse_epi_week)
        weekly = weekly.dropna(subset=['date'])
        weekly = weekly.sort_values('date')

        return weekly[['date', 'cases']]

    raise RuntimeError(
        "Could not find expected columns in dengue data. "
        f"Available columns: {df.columns.tolist()}"
    )


def parse_epi_week(epi_week):
    """
    Parse epidemiological week string to date
    Format examples: "2012-W01", "2012W01", "202201"
    """
    try:
        epi_week = str(epi_week)

        # Try various formats
        if '-W' in epi_week or '-w' in epi_week:
            # Format: YYYY-Www
            return datetime.strptime(epi_week + '-1', '%Y-W%W-%w')
        elif 'W' in epi_week.upper():
            # Format: YYYYWww
            parts = epi_week.upper().split('W')
            year = int(parts[0])
            week = int(parts[1])
            return datetime.strptime(f'{year}-W{week:02d}-1', '%Y-W%W-%w')
        elif len(epi_week) == 6:
            # Format: YYYYWW
            year = int(epi_week[:4])
            week = int(epi_week[4:])
            return datetime.strptime(f'{year}-W{week:02d}-1', '%Y-W%W-%w')
        else:
            return pd.to_datetime(epi_week)
    except:
        return None



# =============================================================================
# 2. METEOROLOGICAL DATA FROM METEOSTAT
# =============================================================================

def download_weather_data():
    """
    Download daily weather data from Meteostat for Singapore Changi Airport
    Station WMO: 48698
    """
    print("\n" + "=" * 60)
    print("Downloading weather data from Meteostat...")
    print("=" * 60)

    try:
        from meteostat import Point, Daily

        # Singapore Changi Airport coordinates
        singapore = Point(1.3521, 103.8198, 15)  # lat, lon, elevation

        # Date range
        start = datetime(2012, 1, 1)
        end = datetime(2022, 12, 31)

        # Get daily data
        df = Daily(singapore, start, end).fetch()

        if len(df) > 0:
            print(f"  Downloaded {len(df)} daily records")
            df = df.reset_index()
            return df
        else:
            print("  No data returned, trying alternative station...")

    except ImportError:
        print("  Meteostat not installed. Install with: pip install meteostat")
    except Exception as e:
        print(f"  Error downloading from Meteostat: {e}")

    raise RuntimeError(
        "Weather data download failed. Install meteostat (pip install meteostat) and retry."
    )



def process_weather_data(df):
    """
    Aggregate daily weather to weekly means/totals
    """
    print("\nProcessing weather data to weekly aggregates...")

    df['time'] = pd.to_datetime(df['time'])
    df = df.set_index('time')

    # Weekly aggregates: mean temperature, total rainfall
    weekly = df.resample('W-SUN').agg({
        'tavg': 'mean',   # Mean temperature
        'prcp': 'sum'     # Total precipitation
    }).reset_index()

    weekly.columns = ['date', 'temp_mean', 'rainfall_total']

    # Handle missing values
    weekly['temp_mean'] = weekly['temp_mean'].interpolate(method='linear')
    weekly['rainfall_total'] = weekly['rainfall_total'].fillna(0)

    print(f"  Created {len(weekly)} weekly records")

    return weekly


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("DENGUE RT ESTIMATION - DATA ACQUISITION")
    print("=" * 70)

    # 1. Download and process dengue data
    dengue_raw = download_dengue_data()
    dengue_df = process_dengue_data(dengue_raw)
    dengue_df.to_csv(os.path.join(DATA_DIR, 'raw_dengue_cases.csv'), index=False)
    print(f"\nSaved: data/raw_dengue_cases.csv ({len(dengue_df)} rows)")

    # 2. Download and process weather data
    weather_raw = download_weather_data()
    weather_df = process_weather_data(weather_raw)
    weather_df.to_csv(os.path.join(DATA_DIR, 'raw_weather.csv'), index=False)
    print(f"Saved: data/raw_weather.csv ({len(weather_df)} rows)")

    print("\n" + "=" * 70)
    print("DATA ACQUISITION COMPLETE")
    print("=" * 70)
    print("\nNext step: Run 03_prepare_model_data.R to merge and format data for Stan")
