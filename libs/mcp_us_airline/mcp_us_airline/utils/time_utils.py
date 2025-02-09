import pandas as pd
import pytz
from datetime import datetime
import timezonefinder

def convert_time_to_utc (time, airport_code, airport_dict):
    lat, lng = airport_dict.get(airport_code)
    tf = timezonefinder.TimezoneFinder()
    timezone_str = tf.certain_timezone_at(lat=lat, lng=lng)
    timezone = pytz.timezone(timezone_str)
    return timezone.normalize(timezone.localize(time)).astimezone(pytz.utc).timestamp()

def convert_to_datetime(row, time_col):
    """
    Converts separate year, month, day, and time columns into a single datetime object.
    :param row: DataFrame row
    :param time_col: Column name for the time (in HHMM format)
    :return: datetime object or None if data is missing
    """
    # Extract the components
    year = row['Year']
    month = row['Month']
    day = row['DayofMonth']
    time_val = row[time_col]
    
    # Handle missing values directly
    if pd.isnull([year, month, day, time_val]).any():
        return None
    
    # Ensure time is in HHMM format, then convert it to a full datetime string
    time_str = str(int(time_val)).zfill(4)  # Make sure time is in HHMM format
    time_str = f"{time_str[:2]}:{time_str[2:]}:00"  # Convert to HH:MM:SS

    # Combine year, month, day, and time into a single string and convert to datetime
    date_str = f"{int(year):04d}-{int(month):02d}-{int(day):02d} {time_str}"
    
    # Convert to datetime object
    return pd.to_datetime(date_str, format="%Y-%m-%d %H:%M:%S")