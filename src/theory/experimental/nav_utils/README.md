# Navigation Utilities Package

This package contains reusable functions for GPS data processing and visualization.

## Structure

The package is organized into the following modules:

- `time_utils.py`: Functions for GPS time conversion
- `satellite_position.py`: Functions for computing satellite positions from ephemeris data
- `visualization.py`: Functions for visualizing satellite positions on 2D maps and 3D globes

## Usage

You can import individual functions:

```python
from src.theory.experimental.nav_utils import generate_gps_prn_list, compute_satellite_position
```

Or import all functions:

```python
from src.theory.experimental.nav_utils import (
    generate_gps_prn_list,
    gps_to_datetime,
    datetime_to_gps_seconds,
    solve_kepler,
    compute_satellite_position,
    plot_satellite_worldmap,
    plot_satellite_worldmap_3d,
    plot_satellite_3d_plotly
)
```

## Example

See `nav_coord_computation_refactored.py` for a complete example of how to use these functions.

## Functions

### Time Utilities

- `generate_gps_prn_list(max_prn=32)`: Generate a list of GPS PRN identifiers
- `gps_to_datetime(gps_week, seconds_of_week)`: Convert GPS week and seconds to datetime
- `datetime_to_gps_seconds(obs_time)`: Convert datetime to GPS week and seconds

### Satellite Position

- `solve_kepler(M, e, tol=1e-10, max_iter=30)`: Solve Kepler's equation using Newton-Raphson
- `compute_satellite_position(eph, transmit_time)`: Compute GPS satellite position in ECEF coordinates

### Visualization

- `plot_satellite_worldmap(positions, obs_time=None, elevation_mask=0, receiver_lla=None)`: Plot satellites on a 2D world map
- `plot_satellite_worldmap_3d(positions, obs_time=None, elevation_mask=0, receiver_lla=None)`: Plot satellites on a 3D globe
- `plot_satellite_3d_plotly(positions, receiver_lla=None)`: Create an interactive 3D plot with Plotly