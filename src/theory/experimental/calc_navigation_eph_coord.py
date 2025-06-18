"""
GPS satellite position computation and visualization.

This script demonstrates how to compute GPS satellite positions from broadcast
ephemeris data and visualize them on 2D maps and 3D globes.
"""

import os
import sys
from datetime import datetime

import georinex as gr
import matplotlib.pyplot as plt
import pymap3d as pm

# Import functions from the nav_utils package
from src.theory.experimental.nav_utils import (
    generate_gps_prn_list,
    datetime_to_gps_seconds,
    compute_satellite_position,
    plot_satellite_worldmap_3d,
    plot_satellite_3d_plotly
)

# https://gage.upc.edu/en/learning-materials/library/gnss-format-descriptions
# https://server.gage.upc.edu/gLAB/HTML/GPS_Navigation_Rinex_v3.04.html

"""
Two different approaches are followed by the GPS/Galileo/Beidou and
Glonass satellites to account for satellite orbit perturbations. These ap-
proaches define what their messages contain.

In the case of the GPS, Galileo or Beidou satellites, the orbits are seen
as Keplerian in a first approximation, and the perturbations are treated as
temporal variations in the orbital elements.
Indeed, an extended set of 16 quasi-Keplerian parameters 
is broadcast to the user in the navigation message and regularly updated.
This expanded set consists of the six orbital elements (a(t), e(t), i(t), Ω(t),
ω(t), M (t)) and three rate parameters to account for the linear changes with
time (Ω, i, ∆n), three pairs of sinusoidal corrections (Cc , Cs ) (i.e. Cc cos(2φ),
Cs sin(2φ)), and the reference ephemeris epoch toe.

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

fig, ax = plot_satellite_worldmap_3d(positions=positions,
                                     obs_time=obs_time,
                                     elevation_mask=30,
                                     receiver_lla=[41.1, 2.2, 0])
plt.show()
plt.close()

plot_satellite_3d_plotly(positions, receiver_lla=[41.1, 2.2, 0])

sys.exit()
