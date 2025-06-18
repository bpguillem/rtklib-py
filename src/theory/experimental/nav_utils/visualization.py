"""
Satellite visualization utilities.

This module provides functions for visualizing GPS satellite positions
on 2D maps and 3D globes.
"""

import numpy as np
import pymap3d as pm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.pyplot import tight_layout


def plot_satellite_worldmap(positions, obs_time=None, elevation_mask=0, receiver_lla=None):
    """
    Plot satellite positions on a world map with visibility indicators.

    Args:
        positions (dict): Dictionary of {PRN: ECEF_position} from compute_satellite_position()
        obs_time (datetime): Observation time (for title)
        elevation_mask (float): Minimum elevation angle to consider visible (degrees)
        receiver_lla (list): Receiver LLA
        
    Returns:
        tuple: (figure, axis) matplotlib objects
    """
    # Create figure with Plate Carrée projection
    fig = plt.figure(figsize=(15, 8), tight_layout=True)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    # Add map features
    ax.set_global()  # Full world coverage
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
        title += f"; (Elevation > {elevation_mask}°)"

    plt.title(title, fontsize=14)
    plt.legend(['All Satellites', 'Receiver', 'Visible'], loc='lower left')
    return fig, ax


def plot_satellite_worldmap_3d(positions, obs_time=None, elevation_mask=0, receiver_lla=None):
    """
    Plot satellite positions on a 3D globe with visibility indicators.

    Args:
        positions (dict): Dictionary of {PRN: ECEF_position} from compute_satellite_position()
        obs_time (datetime): Observation time (for title)
        elevation_mask (float): Minimum elevation angle to consider visible (degrees)
        receiver_lla (list): Receiver LLA [lat, lon, alt]
        
    Returns:
        tuple: (figure, axis) matplotlib objects
    """
    # Create figure with 3D axis
    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(111, projection='3d')

    # Define Earth radius
    earth_radius = 6371000  # meters

    # Create sphere coordinates
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    u, v = np.meshgrid(u, v)

    x = earth_radius * np.cos(u) * np.sin(v)
    y = earth_radius * np.sin(u) * np.sin(v)
    z = earth_radius * np.cos(v)

    # Plot Earth as a semi-transparent sphere
    ax.plot_surface(x, y, z, color='lightblue', alpha=0.6, edgecolor='gray')

    # Mark positions
    for prn, ecef in positions.items():
        ax.scatter(*ecef, color='red', s=50, label='Satellites' if prn == list(positions.keys())[0] else "")
        ax.text(*ecef, prn, color='red', fontsize=10, weight='bold')

    # Add title
    title = "GPS Satellite Positions (3D View)"
    if obs_time:
        title += f"\n{obs_time.strftime('%Y-%m-%d %H:%M:%S UTC')}"

    if receiver_lla is None:
        ax.set_title(title, fontsize=14)
        ax.legend()
        return fig, ax

    # Mark receiver position (convert LLA to ECEF)
    rx, ry, rz = pm.geodetic2ecef(receiver_lla[0], receiver_lla[1], receiver_lla[2])
    ax.scatter(rx, ry, rz, color='blue', s=100, marker='*', label='Receiver')

    # Calculate actual visibility
    visible = []
    for prn, ecef in positions.items():
        az, el, _ = pm.ecef2aer(ecef[0], ecef[1], ecef[2], *receiver_lla)
        if el > elevation_mask:
            visible.append(prn)
            # Draw line from receiver to visible satellite
            ax.plot([rx, ecef[0]], [ry, ecef[1]], [rz, ecef[2]],
                    'g--', alpha=0.3, label='Visible' if prn == visible[0] else "")

    # Highlight visible satellites
    for prn in visible:
        x, y, z = positions[prn]
        ax.scatter(x, y, z, color='green', s=50)

    if elevation_mask > 0:
        title += f"; (Elevation > {elevation_mask}°)"

    ax.set_title(title, fontsize=14)
    ax.legend(loc='upper left')

    # Equal aspect ratio
    ax.set_box_aspect([1, 1, 1])
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')

    # Set limits to match sat orbit altitude
    lim = 23000000
    ax.set_xlim([-lim, lim])
    ax.set_ylim([-lim, lim])
    ax.set_zlim([-lim, lim])

    ax.view_init(elev=20, azim=100)

    return fig, ax


def plot_satellite_3d_plotly(positions, receiver_lla=None):
    """
    Fully interactive 3D plot with Plotly.
    
    Args:
        positions (dict): Dictionary of {PRN: ECEF_position} from compute_satellite_position()
        receiver_lla (list): Receiver LLA [lat, lon, alt]
    """
    import plotly.graph_objects as go

    earth_radius = 6371  # km

    # Create sphere
    theta = np.linspace(0, 2 * np.pi, 100)
    phi = np.linspace(0, np.pi, 100)
    x = earth_radius * np.outer(np.cos(theta), np.sin(phi))
    y = earth_radius * np.outer(np.sin(theta), np.sin(phi))
    z = earth_radius * np.outer(np.ones(100), np.cos(phi))

    fig = go.Figure()

    # Add Earth
    fig.add_trace(go.Surface(x=x, y=y, z=z, colorscale='Blues', opacity=1.))

    # Add satellites
    for prn, ecef in positions.items():
        x, y, z = ecef[0] / 1000, ecef[1] / 1000, ecef[2] / 1000  # Convert to km
        fig.add_trace(go.Scatter3d(x=[x], y=[y], z=[z],
                                   mode='markers+text',
                                   marker=dict(size=4, color='red'),
                                   text=prn,
                                   textposition="top center",
                                   name=prn))

    # Add receiver
    if receiver_lla:
        rx, ry, rz = pm.geodetic2ecef(*receiver_lla)
        fig.add_trace(go.Scatter3d(x=[rx / 1000], y=[ry / 1000], z=[rz / 1000],
                                   mode='markers',
                                   marker=dict(size=4, color='blue', symbol='cross'),
                                   name='Receiver'))

    # Set axis limits and layout
    lim = 30000  # km
    fig.update_layout(
        scene=dict(
            xaxis=dict(range=[-lim, lim], title='X (km)'),
            yaxis=dict(range=[-lim, lim], title='Y (km)'),
            zaxis=dict(range=[-lim, lim], title='Z (km)'),
            aspectmode='cube'  # Maintain equal aspect ratio
        ),
        title="GPS Satellites (Interactive 3D)"
    )

    fig.show()