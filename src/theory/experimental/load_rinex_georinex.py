import os

import georinex as gr
import matplotlib.pyplot as plt
import pymap3d as pm


def process_obs_data(file_path, time_limits):
    """
    Process observation data from a RINEX file.

    Args:
        file_path (str): Path to the RINEX observation file
        time_limits (list): Time limits for data loading [start, end]

    Returns:
        tuple: (obs_data, dataframe) - The observation data and its DataFrame representation
    """
    print(f'Reading obs data from {os.path.basename(file_path)}...')
    obs_data = gr.load(file_path, tlim=time_limits)

    print(f'Position ECEF: {obs_data.position}')
    llh = pm.ecef2geodetic(*obs_data.position)
    print(f'Position LLH: {llh}')
    print(f'Position LLH: {obs_data.position_geodetic}')

    print(f'Tracking of satellites: {obs_data.sv.data}')

    # Convert to dataframe
    df = obs_data.to_dataframe()
    df = df.reset_index()

    return obs_data, df


def plot_satellite_cno(df, obs_columns):
    """
    Plot observables for each satellite.

    Args:
        df (DataFrame): DataFrame containing observation data
        obs_columns (list): List of column names to plot
    """
    print('Plotting...')

    unique_svs = df['sv'].unique()
    n_svs = len(unique_svs)

    # Plot each SV in its own subplot
    fig, axes = plt.subplots(n_svs, 1, figsize=(12, 3 * n_svs), tight_layout=True, sharex=True, sharey=True)
    axes = axes.flatten() if n_svs > 1 else [axes]

    for i, sv in enumerate(unique_svs):
        df_sv = df[df['sv'] == sv]
        for obs in obs_columns:
            axes[i].plot(df_sv['time'], df_sv[obs], ls='none', marker='.', label=obs)
        axes[i].set_title(f'Satellite {sv}')
        axes[i].set_ylabel('Observables')
        axes[i].legend(loc='upper left')
        axes[i].grid(True)

    plt.show()
    plt.close()


# u-blox example
datadir = '/home/guillem/Code/rtklib-py/data/wine/static'
nav_file = 'base_COM7___460800_250416_143415.nav'
rover_file = 'rover1_COM13___460800_250416_143419.obs'
base_file = 'base_COM7___460800_250416_143415.obs'

# Define time limits for data loading
time_limits = ['2025-04-16T14:34', '2025-04-17T14:34']

# Process rover data
rover_path = os.path.join(datadir, rover_file)
rover, rover_df = process_obs_data(rover_path, time_limits)

# List of observable columns (for reference)
observables_col = ['C1C', 'L1C', 'D1C', 'S1C', 'C2X', 'L2X', 'D2X', 'S2X']

# CNO columns to plot
obs_cno = ['S1C', 'S2X']

# Plot rover data
plot_satellite_cno(rover_df, obs_cno)

# Process base data
base_path = os.path.join(datadir, base_file)
base, base_df = process_obs_data(base_path, time_limits)

# Plot base data
plot_satellite_cno(base_df, obs_cno)
