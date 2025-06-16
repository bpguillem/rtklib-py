import matplotlib.pyplot as plt
import numpy as np

# GPS L1 signal parameters
GPS_L1_freq = 1575.42e6  # Hz
wavelength = 3e8 / GPS_L1_freq
d = 1.3  # meter separation between rovers

# Rover positions
rover_A = np.array([0, 0])
rover_B = np.array([d, 0])

# Incident angles in degrees
incident_angles_deg = [-10, -30, -60]
colors = plt.cm.viridis(np.linspace(0, 1, len(incident_angles_deg)))

# Plot setup
fig, ax = plt.subplots(figsize=(12, 10))

# Plot rovers
ax.plot(rover_A[0], rover_A[1], 'ro', markersize=10, label='Rover A')
ax.plot(rover_B[0], rover_B[1], 'go', markersize=10, label='Rover B')
ax.plot([rover_A[0], rover_B[0]], [rover_A[1], rover_B[1]], 'k:', alpha=0.3)

# Loop through each angle
for idx, angle_deg in enumerate(incident_angles_deg):
    theta = np.radians(angle_deg)
    k = np.array([np.sin(theta), -np.cos(theta)])

    # Draw wavefronts
    wavefront_spacing = wavelength
    num_wavefronts = 10
    for n in range(-2, num_wavefronts - 2):
        x0 = rover_A[0] + n * wavefront_spacing * k[0]
        y0 = rover_A[1] + n * wavefront_spacing * k[1]

        wavefront_length = 2
        wavefront_x = [x0 - wavefront_length * k[1], x0 + wavefront_length * k[1]]
        wavefront_y = [y0 + wavefront_length * k[0], y0 - wavefront_length * k[0]]

        ax.plot(wavefront_x, wavefront_y, '-', alpha=0.2, color=colors[idx])

    # Propagation direction arrow
    ax.arrow(rover_A[0] - 2 * k[0], rover_A[1] - 2 * k[1],
             2.5 * k[0], 2.5 * k[1],
             head_width=0.1, head_length=0.1, ls='--',
             fc='none', ec=colors[idx], alpha=0.8)

    # Wavefront through Rover A
    wavefront_at_A_x = [-2 * k[1], 2 * k[1]]
    wavefront_at_A_y = [2 * k[0], -2 * k[0]]
    # if idx == 0:
    #     ax.plot(wavefront_at_A_x, wavefront_at_A_y, '--', linewidth=1.5,
    #             label='Wavefront through Rover A', color='r', alpha=0.6)

    # Intersection of wavefront at A and direction from Rover B
    t_intersect = -np.dot(k, rover_B) / np.dot(k, k)
    intersection_point = rover_B + t_intersect * k
    additional_path = -t_intersect
    additional_wavelengths = additional_path / wavelength
    integer_cycles = int(additional_wavelengths)
    fractional_cycles = additional_wavelengths - integer_cycles

    # Line showing the additional path
    ax.plot([rover_B[0], intersection_point[0]],
            [rover_B[1], intersection_point[1]],
            '-', linewidth=2, color=colors[idx],
            label=f'θ = {angle_deg}°: {additional_path:.2f}m '
                  f'({integer_cycles}+{fractional_cycles:.2f}λ)')

    # Mark intersection
    ax.plot(intersection_point[0], intersection_point[1], 'x', color=colors[idx])

# Final touches
ax.set_xlim(-1, 2)
ax.set_ylim(-1.5, 1.5)
ax.set_aspect('equal')
ax.set_xlabel('X position (meters)')
ax.set_ylabel('Y position (meters)')
ax.set_title('Wavefronts and Phase Differences at Multiple Incident Angles')
ax.legend(loc='upper right', fontsize=9)
ax.grid(False)
plt.tight_layout()
plt.show()
