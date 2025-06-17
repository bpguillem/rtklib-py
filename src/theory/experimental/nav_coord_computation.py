import os
import sys
from datetime import datetime

import georinex as gr
import numpy as np
import pandas as pd

# https://gage.upc.edu/en/learning-materials/library/gnss-format-descriptions
# https://server.gage.upc.edu/gLAB/HTML/GPS_Navigation_Rinex_v3.04.html

"""
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


def compute_satellite_position(eph, transmit_time):
    """
       Compute GPS satellite position in ECEF coordinates at given transmit time
       Args:
           eph: Series containing broadcast ephemeris parameters for one satellite
           transmit_time: GPS datetime when signal was transmitted
       Returns:
           numpy array with [X, Y, Z] in meters
       """
    # Constants
    mu = 3.986005e14  # Earth's gravitational constant (m^3/s^2)
    omega_e = 7.2921151467e-5  # Earth rotation rate (rad/s)

    # Time from ephemeris reference epoch (seconds)
    tk = (transmit_time - eph['Toe']).total_seconds()

    # Corrected mean motion (rad/s)
    n0 = np.sqrt(mu) / (eph['sqrtA'] ** 3)
    n = n0 + eph['DeltaN']

    # Mean anomaly (rad)
    Mk = eph['M0'] + n * tk

    # Solve Kepler's equation for eccentric anomaly Ek (rad)
    Ek = Mk
    for _ in range(10):  # Simple iteration
        Ek_new = Mk + eph['Eccentricity'] * np.sin(Ek)
        if abs(Ek_new - Ek) < 1e-12:
            break
        Ek = Ek_new

    # True anomaly (rad)
    nu_k = np.arctan2(np.sqrt(1 - eph['Eccentricity'] ** 2) * np.sin(Ek),
                      np.cos(Ek) - eph['Eccentricity'])

    # Argument of latitude (rad)
    phi_k = nu_k + eph['Omega']

    # Second harmonic corrections
    du_k = eph['Cuc'] * np.cos(2 * phi_k) + eph['Cus'] * np.sin(2 * phi_k)  # Argument of latitude correction
    dr_k = eph['Crc'] * np.cos(2 * phi_k) + eph['Crs'] * np.sin(2 * phi_k)  # Radius correction
    di_k = eph['Cic'] * np.cos(2 * phi_k) + eph['Cis'] * np.sin(2 * phi_k)  # Inclination correction

    # Corrected arguments
    u_k = phi_k + du_k  # Argument of latitude
    r_k = eph['sqrtA'] ** 2 * (1 - eph['Eccentricity'] * np.cos(Ek)) + dr_k  # Radius
    i_k = eph['Io'] + di_k + eph['IDOT'] * tk  # Inclination

    # Position in orbital plane
    x_k_prime = r_k * np.cos(u_k)
    y_k_prime = r_k * np.sin(u_k)

    # Corrected longitude of ascending node (rad)
    omega_k = eph['Omega0'] + (eph['OmegaDot'] - omega_e) * tk - omega_e * eph['Toe'].timestamp()

    # ECEF coordinates (m)
    x_k = x_k_prime * np.cos(omega_k) - y_k_prime * np.cos(i_k) * np.sin(omega_k)
    y_k = x_k_prime * np.sin(omega_k) + y_k_prime * np.cos(i_k) * np.cos(omega_k)
    z_k = y_k_prime * np.sin(i_k)

    return np.array([x_k, y_k, z_k])


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
obs_time = datetime(2025, 4, 16, 16, 0, 0)

positions = {}

for prn in prn_active_list:
    print(f'Processing PRN {prn}')

    # Get ephemeris for specified PRN
    eph = nav.sel(sv=prn).to_dataframe()
    eph = eph.dropna(thresh=2)

    # Times
    initial_time = eph.index[0]
    end_time = eph.index[-1]
    print(f'Initial time: {initial_time}')
    print(f'End time: {end_time}')
    print(f'Duration: {end_time - initial_time}')
    print(f'Mean interval: {eph.index.diff().mean()}')

    # Convert Toe from time delta to datetime
    # Toe is the reference time for orbit calculations
    # TransTime helps verify ephemeris freshness
    eph['Toe'] = pd.to_datetime(eph['Toe'], unit='s')

    # Compute position
    result = compute_satellite_position(eph, obs_time)

    # Store
    positions[prn] = result

sys.exit()
