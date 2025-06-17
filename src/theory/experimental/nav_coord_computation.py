import os
import sys
from datetime import datetime, timezone, timedelta

import georinex as gr
import numpy as np
import pandas as pd
import pymap3d as pm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.pyplot import tight_layout

# https://gage.upc.edu/en/learning-materials/library/gnss-format-descriptions
# https://server.gage.upc.edu/gLAB/HTML/GPS_Navigation_Rinex_v3.04.html

"""
These parameters are renewed periodically (typically every two hours for GPS)
and must not be used after the prescribed time (about four hours), because
the extrapolation error grows exponentially beyond this validity period.

GPS RINEX 3.04 Navigation Message Fields Explained

1. Satellite Identification & Metadata
    - sv: Satellite ID (e.g., "G01" for GPS PRN 1)
    - health: Satellite health status (0 = healthy, non-zero indicates problems)
    - GPSWeek: GPS week number (modulo 1024) when ephemeris was uploaded

2. Clock Parameters
    - SVclockBias (a₀): Clock bias in seconds (constant term)
    - SVclockDrift (a₁): Clock drift in sec/sec (linear term)
    - SVclockDriftRate (a₂): Clock drift rate in sec/sec² (quadratic term)
    - TGD: Total Group Delay (ionospheric correction term for L1/L2 signals) in seconds

3. Ephemeris Parameters (Keplerian Elements)
    - IODE: Issue of Data Ephemeris (changes when new data is uploaded)
    - IODC: Issue of Data Clock (should match IODE for valid ephemeris)
    - Toe: Time of Ephemeris (reference time for orbit parameters, in seconds of week)
    - TransTime: Transmission time of message (seconds of week)

    Orbit Shape Parameters
    - sqrtA: Square root of semi-major axis (√a, in m^(1/2))
    - Eccentricity (e): Orbit eccentricity (dimensionless)
    - M0: Mean anomaly at reference time (radians)
    - DeltaN (Δn): Mean motion correction from computed value (rad/s)

    Orbit Orientation Parameters
    - Omega0 (Ω₀): Longitude of ascending node at weekly epoch (radians)
    - OmegaDot (Ω̇): Rate of right ascension (rad/s)
    - Io (i₀): Inclination angle at reference time (radians)
    - IDOT (i̇): Rate of inclination angle (rad/s)
    - omega (ω): Argument of perigee (radians)

    Harmonic Correction Terms
    - Crs: Sine correction to orbit radius (m)
    - Crc: Cosine correction to orbit radius (m)
    - Cus: Sine correction to argument of latitude (rad)
    - Cuc: Cosine correction to argument of latitude (rad)
    - Cis: Sine correction to inclination (rad)
    - Cic: Cosine correction to inclination (rad)

4. Signal & Quality Parameters
    - CodesL2: Codes on L2 channel (0=invalid, 1=P-code, 2=C/A-code, etc.)
    - L2Pflag: L2 P-code data flag (0=unknown, 1=on, 2=off)
    - SVacc: User range accuracy estimate (URA index, converted to meters)
    - FitIntvl: Fit interval flag (0=4hrs, 1=longer than 4hrs)

5. Reserved Fields
    - spare0, spare1: Reserved for future use
"""


def generate_gps_prn_list(max_prn=32):
    # GPS: G01-G37 (33-37 are rarely used)
    # Galileo: E01-E36
    # GLONASS: R01-R24
    # BeiDou: C01-C37
    # QZSS: J01-J07
    # IRNSS: I01-I07
    # SBAS: S01-S39
    return [f"G{i:02d}" for i in range(1, max_prn + 1)]


def gps_to_datetime(gps_week, seconds_of_week):
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


def solve_kepler(M, e, tol=1e-10, max_iter=30):
    """Solve Kepler's equation M = E - e sin(E) using Newton-Raphson."""
    E = M  # Initial guess
    for _ in range(max_iter):
        delta = (E - e * np.sin(E) - M) / (1.0 - e * np.cos(E))
        E -= delta
        if abs(delta) < tol:
            break
    print('Newton-Raphson solver did not converge.')
    return E


def compute_satellite_position(eph, transmit_time):
    """
       Compute GPS satellite position in ECEF coordinates at given transmit time
       Source: [Section 3.3.1 Computation of GPS, Galileo and Beidou Coordinates,
       GNSS DATA PROCESSING - Volume I: Fundamentals and Algorithms, ESA]
       Args:
           eph: Series containing broadcast ephemeris parameters for one satellite
           transmit_time: seconds within the week GPS when signal was transmitted
       Returns:
           numpy array with [X, Y, Z] in meters
       """
    # Constants
    mu = 3.986005e14  # Earth's gravitational constant (m^3/s^2)
    omega_e = 7.2921151467e-5  # Earth rotation rate (rad/s)

    # Time from ephemeris reference epoch (seconds within the week)
    tk = transmit_time - eph['Toe']

    if isinstance(eph, pd.DataFrame):
        closest_time_arg = tk.abs().argmin()
        print(f'    Closest navigation message at {eph.index[closest_time_arg]}')
        eph = eph.iloc[closest_time_arg]
        tk = tk.iloc[closest_time_arg]

    assert np.abs(tk) <= eph['FitIntvl'] * 60 * 60, f'Observation time outside fit interval (4 hours): {np.abs(tk)}'

    if tk > 302400:
        tk = tk - 604800
    if tk < -302400:
        tk = tk + 604800

    # Corrected mean motion (rad/s)
    a = eph['sqrtA'] ** 2
    n0 = np.sqrt(mu) / np.sqrt(a ** 3)
    n = n0 + eph['DeltaN']

    # Mean anomaly (rad)
    Mk = eph['M0'] + n * tk

    # Solve Kepler's equation for eccentric anomaly Ek (rad)
    Ek = solve_kepler(Mk, eph['Eccentricity'], tol=1e-12, max_iter=30)

    # True anomaly (rad)
    nu_k = np.arctan2(np.sqrt(1 - eph['Eccentricity'] ** 2) * np.sin(Ek),
                      np.cos(Ek) - eph['Eccentricity'])

    # Argument of latitude (rad) from the argument of perigee
    phi_k = nu_k + eph['omega']

    # Second harmonic corrections
    du_k = eph['Cuc'] * np.cos(2 * phi_k) + eph['Cus'] * np.sin(2 * phi_k)  # Argument of latitude correction
    dr_k = eph['Crc'] * np.cos(2 * phi_k) + eph['Crs'] * np.sin(2 * phi_k)  # Radius correction
    di_k = eph['Cic'] * np.cos(2 * phi_k) + eph['Cis'] * np.sin(2 * phi_k)  # Inclination correction

    # Corrected arguments
    u_k = phi_k + du_k  # Argument of latitude
    r_k = eph['sqrtA'] ** 2 * (1 - eph['Eccentricity'] * np.cos(Ek)) + dr_k  # Radial distance
    i_k = eph['Io'] + di_k + eph['IDOT'] * tk  # Inclination

    # Position in orbital plane
    x_k_prime = r_k * np.cos(u_k)
    y_k_prime = r_k * np.sin(u_k)

    # Corrected longitude of ascending node (rad) using uses the right ascension at the beginning of the current week
    omega_k = eph['Omega0'] + (eph['OmegaDot'] - omega_e) * tk - omega_e * eph['Toe']

    # ECEF coordinates (m)
    x_k = x_k_prime * np.cos(omega_k) - y_k_prime * np.cos(i_k) * np.sin(omega_k)
    y_k = x_k_prime * np.sin(omega_k) + y_k_prime * np.cos(i_k) * np.cos(omega_k)
    z_k = y_k_prime * np.sin(i_k)

    return np.array([x_k, y_k, z_k])


def plot_satellite_worldmap(positions, obs_time=None, elevation_mask=0, receiver_lla=None):
    """
    Plot satellite positions on a world map with visibility indicators.

    Args:
        positions (dict): Dictionary of {PRN: ECEF_position} from compute_satellite_position()
        obs_time (datetime): Observation time (for title)
        elevation_mask (float): Minimum elevation angle to consider visible (degrees)
    """
    # Create figure with Plate Carrée projection
    fig = plt.figure(figsize=(15, 8), tight_layout=True)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    ax.set_global()  # Key addition for full world coverage

    # Add map features
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    ax.add_feature(cfeature.COASTLINE, edgecolor='black')
    # ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.gridlines(draw_labels=True, linestyle='--', alpha=0.7)

    # Convert ECEF to geodetic and mark positions
    for prn, ecef in positions.items():
        lat, lon, alt = pm.ecef2geodetic(ecef[0], ecef[1], ecef[2])
        ax.plot(lon, lat, 'ro', markersize=8, transform=ccrs.Geodetic())
        ax.text(lon + 2, lat, prn, transform=ccrs.Geodetic(), fontsize=10, weight='bold', color='red')

    # Add title and legend
    title = f"GPS Satellite Positions"
    if obs_time:
        title += f"\n{obs_time.strftime('%Y-%m-%d %H:%M:%S UTC')}"

    if receiver_lla is None:
        plt.title(title, fontsize=14)
        plt.legend(['Satellites'], loc='lower left')
        return fig, ax

    # Mark receiver position
    ax.plot(receiver_lla[1], receiver_lla[0], 'b*', markersize=12, transform=ccrs.Geodetic(), label='Receiver')

    # Calculate actual visibility
    visible = []
    for prn, ecef in positions.items():
        az, el, _ = pm.ecef2aer(ecef[0], ecef[1], ecef[2], *receiver_lla)
        if el > elevation_mask:
            visible.append(prn)

    # Highlight visible satellites
    for prn in visible:
        lat, lon, _ = pm.ecef2geodetic(*positions[prn])
        ax.plot(lon, lat, 'go', markersize=8, transform=ccrs.Geodetic())

    if elevation_mask > 0:
        title += f"\n(Elevation > {elevation_mask}°)"

    plt.title(title, fontsize=14)
    plt.legend(['All Satellites', 'Receiver', 'Visible'], loc='lower left')
    return fig, ax


# u-blox example
datadir = '/home/guillem/Code/rtklib-py/data/wine/static'
# nav_file = 'base_COM7___460800_250416_143415.nav'
nav_file = 'rover1_COM13___460800_250416_143419.nav'

# Load GPS RINEX file
nav_path = os.path.join(datadir, nav_file)
nav = gr.load(nav_path, use='G')

# Print navigation data
print(nav)

# Satellite and time of interest
prn_list = generate_gps_prn_list()
prn_active_list = [sv for sv in nav.sv.values if sv.startswith('G')]
print(f'Missing PRN: {list(set(prn_list) - set(prn_active_list))}')
obs_time = datetime(2025, 4, 16, 16, 35, 0)
obs_time_gps_week, obs_time_gps_seconds = datetime_to_gps_seconds(obs_time)

positions = {}

for prn in prn_active_list:
    print(f'Processing PRN {prn}')

    # Get broadcast ephemeris for specified PRN
    ephemeris = nav.sel(sv=prn).to_dataframe()
    ephemeris = ephemeris.dropna(thresh=len(ephemeris.columns) - 2)

    # Assert
    assert (ephemeris.health == 0).all(), 'Not all satellites are healthy'

    # Times
    initial_time = ephemeris.index[0]
    end_time = ephemeris.index[-1]
    print(f'    Initial time: {initial_time}')
    print(f'    End time: {end_time}')
    print(f'    Elapsed time: {end_time - initial_time}')
    print(f'    Number of navigation messages: {len(ephemeris)}')
    print(f'    Mean interval: {ephemeris.index.diff().mean()}')
    print(f'    Observation time: {obs_time}')

    # Compute position
    result = compute_satellite_position(ephemeris, obs_time_gps_seconds)
    print(f'    Position ECEF: {result}')
    print(f'    Position LLH: {pm.ecef2geodetic(*result)}')

    # Store
    positions[prn] = result

fig, ax = plot_satellite_worldmap(
    positions=positions,
    obs_time=obs_time,
    elevation_mask=30,
    receiver_lla=[41.1, 2.2, 0]
)
plt.show()
plt.close()

sys.exit()
