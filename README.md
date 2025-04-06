# S-Curve Motion Profile Calculator

This Python script generates S-curve motion profiles for smooth motion control applications. It calculates and visualizes position, velocity, acceleration, and jerk profiles based on given motion parameters.

## Features

- Calculates complete S-curve motion profiles
- Handles both trapezoidal and triangular velocity profiles
- Visualizes all motion parameters (Position, Velocity, Acceleration, Jerk)
- Includes phase transition markers and maximum value indicators
- Robust error handling for invalid inputs

## Requirements

- Python 3.x
- NumPy
- Matplotlib
- SciPy

## Usage

```python
from calculator import calculate_scurve_times, plot_scurve_profile

# Define motion parameters
v_max = 200.0    # Maximum velocity (units/s)
a_max = 1000.0   # Maximum acceleration (units/s²)
d_max = 1000.0   # Maximum deceleration (units/s²)
jerk = 10000.0   # Maximum jerk (units/s³)
distance = 3000.0 # Total distance to move (units)

# Calculate motion profile
results = calculate_scurve_times(v_max, a_max, d_max, jerk, distance)

# Plot the results
plot_scurve_profile(results)
```

## Function Parameters

### calculate_scurve_times()
- `V_max`: Maximum velocity
- `A_max`: Maximum acceleration
- `D_max`: Maximum deceleration
- `Jerk`: Maximum jerk
- `Distance`: Total distance to travel

Returns a dictionary containing:
- Time intervals for each motion phase
- Peak velocity achieved
- Whether maximum velocity was reached
- Total motion time
- Input parameters used

### plot_scurve_profile()
- `results`: Dictionary returned by calculate_scurve_times()
- `num_points`: Number of points to plot (default: 500)

## Examples

The script includes three example scenarios:
1. Long move reaching maximum velocity (trapezoidal profile)
2. Short move with triangular velocity profile
3. High acceleration/deceleration case

## Motion Profile Phases

An S-curve motion consists of up to seven phases:
1. Jerk up to max acceleration (tj1)
2. Constant acceleration (ta)
3. Jerk down to zero acceleration (tj2)
4. Constant velocity (tv)
5. Jerk down to max deceleration (tj3)
6. Constant deceleration (td)
7. Jerk up to zero deceleration (tj4)

## Output Visualization

The script generates four subplots showing:
1. Position profile
2. Velocity profile (with V_peak indicator)
3. Acceleration profile (with A_max/D_max indicators)
4. Jerk profile

Vertical lines indicate phase transitions for easy analysis.
