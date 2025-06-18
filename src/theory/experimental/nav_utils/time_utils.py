"""
GPS time conversion utilities.

This module provides functions for converting between different time formats
used in GPS systems.
"""

from datetime import datetime, timezone, timedelta


def generate_gps_prn_list(max_prn=32):
    """
    Generate a list of GPS PRN identifiers.
    
    Args:
        max_prn (int): Maximum PRN number to generate (default: 32)
        
    Returns:
        list: List of GPS PRN identifiers in the format 'Gxx'
    """
    # GPS: G01-G37 (33-37 are rarely used)
    # Galileo: E01-E36
    # GLONASS: R01-R24
    # BeiDou: C01-C37
    # QZSS: J01-J07
    # IRNSS: I01-I07
    # SBAS: S01-S39
    return [f"G{i:02d}" for i in range(1, max_prn + 1)]


def gps_to_datetime(gps_week, seconds_of_week):
    """
    Convert GPS week and seconds of week to datetime.
    
    Args:
        gps_week (int): GPS week number
        seconds_of_week (float): Seconds of the week
        
    Returns:
        datetime: Python datetime object in UTC
    """
    # https://github.com/GNSSpy-Project/gnsspy/blob/master/gnsspy/funcs/date.py
    return datetime(1980, 1, 6) + timedelta(weeks=gps_week, seconds=seconds_of_week)


def datetime_to_gps_seconds(obs_time: datetime):
    """
    Convert datetime object to GPS seconds of week.

    Args:
        obs_time: Python datetime object (assumes UTC if naive)

    Returns:
        tuple: (GPS week number, GPS seconds of week)
    """
    # For timezone-aware datetimes (convert to UTC first):
    if obs_time.tzinfo is not None:
        obs_time = obs_time.astimezone(timezone.utc).replace(tzinfo=None)

    # GPS epoch (January 6, 1980 00:00:00 UTC)
    gps_epoch = datetime(1980, 1, 6, 0, 0, 0)

    # Calculate time since GPS epoch
    delta = obs_time - gps_epoch

    # Total seconds since GPS epoch
    total_seconds = delta.total_seconds()

    # Calculate GPS week number
    gps_week = int(total_seconds // (7 * 86400))

    # Seconds into current GPS week
    seconds_of_week = total_seconds % (7 * 86400)

    return gps_week, seconds_of_week