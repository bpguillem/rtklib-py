import matplotlib.pyplot as plt
import numpy as np

# GPS L1 signal parameters
GPS_L1_freq = 1575.42e6  # 1575.42 MHz
wavelength = 3e8 / GPS_L1_freq  # ~0.1903 meters
d = 1.3  # meter separation between rovers

# Example phase differencing

# Time parameters
t_sample = 0  # We'll look at a single time sample
phase_at_A = (2 * np.pi * GPS_L1_freq * t_sample - (2 * np.pi / wavelength) * 0) % (2 * np.pi)
phase_at_B = (2 * np.pi * GPS_L1_freq * t_sample - (2 * np.pi / wavelength) * d) % (2 * np.pi)

# Calculate phase difference
total_phase_diff = (phase_at_B - phase_at_A) % (2 * np.pi)
integer_cycles = int(d / wavelength)
fractional_phase = (d / wavelength - integer_cycles) * 2 * np.pi

# Print results in cycles + radians
print(f"GPS L1 Signal Parameters:")
print(f"Frequency: {GPS_L1_freq / 1e6} MHz")
print(f"Wavelength: {wavelength:.4f} meters")
print(f"\nRover Measurements (1 meter separation):")
print(f"Rover A phase: {phase_at_A:.4f} rad ({phase_at_A / (2 * np.pi):.4f} cycles)")
print(f"Rover B phase: {phase_at_B:.4f} rad ({phase_at_B / (2 * np.pi):.4f} cycles)")
print(f"\nPhase Difference Breakdown:")
print(f"Total measured phase difference: {total_phase_diff:.4f} rad")
print(f"Integer cycles difference: {integer_cycles} full wavelengths")
print(f"Fractional phase difference: {fractional_phase:.4f} rad ({fractional_phase / (2 * np.pi):.4f} cycles)")
print(f"Calculated from geometry: {d / wavelength:.4f} total wavelength units")

# Verification
calculated_total = integer_cycles * 2 * np.pi + fractional_phase
print(f"\nVerification:")
print(f"Integer cycles (2π) + fractional: {integer_cycles}*2π + {fractional_phase:.4f} = {calculated_total:.4f} rad")
print(f"Should match: {2 * np.pi * d / wavelength:.4f} rad")

# 2D PLOT
# Parameters
theta = np.radians(30)  # Incidence angle (30 degrees)

# Rover positions
rover_A = np.array([0, 0])
rover_B = np.array([d, 0])

# Wave propagation direction vector
k = np.array([np.sin(theta), -np.cos(theta)])  # Coming from above at angle theta

# Create figure
fig, ax = plt.subplots(figsize=(10, 8))

# Plot rovers
ax.plot(rover_A[0], rover_A[1], 'ro', markersize=10, label='Rover A')
ax.plot(rover_B[0], rover_B[1], 'go', markersize=10, label='Rover B')

# Draw wavefronts (perpendicular to k)
wavefront_spacing = wavelength
num_wavefronts = 10

for n in range(-2, num_wavefronts - 2):
    # Calculate wavefront position (lines perpendicular to k)
    x0 = rover_A[0] + n * wavefront_spacing * k[0]
    y0 = rover_A[1] + n * wavefront_spacing * k[1]

    # Draw wavefront (perpendicular to propagation direction)
    wavefront_length = 2
    wavefront_x = [x0 - wavefront_length * k[1], x0 + wavefront_length * k[1]]
    wavefront_y = [y0 + wavefront_length * k[0], y0 - wavefront_length * k[0]]

    if n == 0:
        ax.plot(wavefront_x, wavefront_y, 'b-', alpha=0.5, label='Wavefronts')
    else:
        ax.plot(wavefront_x, wavefront_y, 'b-', alpha=0.5)

# Draw propagation direction arrow
ax.arrow(rover_A[0] - 2 * k[0], rover_A[1] - 2 * k[1],
         3 * k[0], 3 * k[1],
         head_width=0.1, head_length=0.1, ls='--', fc='none', ec='b', alpha=0.5, label='Propagation direction')

# Find perpendicular intersection with rover A
# The wavefront passing through rover A is when n=0
wavefront_at_A_x = [-2 * k[1], 2 * k[1]]
wavefront_at_A_y = [2 * k[0], -2 * k[0]]

# Draw perpendicular line to wavefront through rover A
ax.plot(wavefront_at_A_x, wavefront_at_A_y, 'r--', linewidth=2, label='Perpendicular through Rover A')

# Perpendicular wavefront through Rover B
# Find the wavefront that intersects Rover B
n_B = (rover_B[0] * k[0] + rover_B[1] * k[1]) / wavefront_spacing
x0_B = n_B * wavefront_spacing * k[0]
y0_B = n_B * wavefront_spacing * k[1]

# Draw wavefront through Rover B
wavefront_B_x = [x0_B - 2 * k[1], x0_B + 2 * k[1]]
wavefront_B_y = [y0_B + 2 * k[0], y0_B - 2 * k[0]]
ax.plot(wavefront_B_x, wavefront_B_y, 'g--', linewidth=2, label='Perpendicular through Rover B')

# Add dashed line between rovers
ax.plot([rover_A[0], rover_B[0]], [rover_A[1], rover_B[1]], 'k:', alpha=0.5)

# # Annotate the angle
# # angle_arc = np.linspace(-theta, 0, 30)
# # ax.plot(0.5 * np.cos(angle_arc), 0.5 * np.sin(angle_arc), 'k-')
# ax.text(0.3, -0.2, f'θ = {np.degrees(theta):.0f}°', fontsize=12)

# Find intersection of wave direction through Rover B with wavefront through A
# Parametric line: r = rover_B + t*k
# We need to find t where this intersects the wavefront through A (k·r = 0)
t_intersect = -np.dot(k, rover_B) / np.dot(k, k)
intersection_point = rover_B + t_intersect * k

# Draw the additional path (in wave propagation direction)
additional_path = -t_intersect  # Since t_intersect is negative
additional_wavelengths = additional_path / wavelength

# Calculate phase differencing
integer_cycles = int(additional_path / wavelength)
fractional_cycles = additional_wavelengths - integer_cycles

ax.plot([rover_B[0], intersection_point[0]],
        [rover_B[1], intersection_point[1]],
        'k-', linewidth=2,
        label=f'Additional path: {additional_path:.3f}m ({integer_cycles} + {fractional_cycles:.3f} cycles)')

# Mark the intersection point
ax.plot(intersection_point[0], intersection_point[1], 'kx', markersize=10)

# Annotate the phase difference
ax.text((rover_B[0] + intersection_point[0]) / 2,
        (rover_B[1] + intersection_point[1]) / 2, 'Δϕ')
# f'Δϕ = {additional_wavelengths * 360:.1f}°\n= {additional_wavelengths:.3f} cycles',
# color='k', ha='center', va='top')

# Formatting
ax.set_xlim(-1, 2)
ax.set_ylim(-1.5, 1.5)
ax.set_aspect('equal')
ax.set_xlabel('X position (meters)')
ax.set_ylabel('Y position (meters)')
ax.set_title(f'GPS Wavefronts at θ = {np.degrees(theta):.0f}°')  # \n'
# f'Additional path: {additional_path:.3f}m ({additional_wavelengths:.3f} wavelengths)')
ax.grid(False)
ax.legend(loc='upper right')

plt.tight_layout()
plt.show()
