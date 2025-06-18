"""
Satellite position computation utilities.

This module provides functions for computing GPS satellite positions
based on broadcast ephemeris data.
"""

import numpy as np
import pandas as pd
import pymap3d as pm


def solve_kepler(M, e, tol=1e-10, max_iter=30):
    """
    Solve Kepler's equation M = E - e sin(E) using Newton-Raphson.
    
    Args:
        M (float): Mean anomaly (radians)
        e (float): Eccentricity
        tol (float): Convergence tolerance
        max_iter (int): Maximum number of iterations
        
    Returns:
        float: Eccentric anomaly (radians)
    """
    E = M  # Initial guess
    for _ in range(max_iter):
        delta = (E - e * np.sin(E) - M) / (1.0 - e * np.cos(E))
        E -= delta
        if abs(delta) < tol:
            return E
    print('Newton-Raphson solver did not converge.')
    return E


def compute_satellite_position(eph, transmit_time):
    """
    Compute GPS satellite position in ECEF coordinates at given transmit time.
    
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