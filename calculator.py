import math
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

# --- Previous calculate_scurve_times function ---
def calculate_scurve_times(V_max, A_max, D_max, Jerk, Distance):
    """
    Calculates the elapsed times for different phases of an S-curve motion profile.
    (Code from the previous response - keeping it concise here)
    ... [Same function code as before] ...
    """
    # --- Input Validation ---
    if not all(val > 0 for val in [V_max, A_max, D_max, Jerk, Distance]):
        return {"error": "Inputs (V_max, A_max, D_max, Jerk, Distance) must all be positive."}
    if A_max**2 / Jerk == 0 or D_max**2 / Jerk == 0:
         # Avoid potential division by zero later if A_max or D_max is effectively zero relative to Jerk
         # Although the initial check covers A/D > 0, extremely small values could cause issues.
         # Let's refine the conditions where V_max comparisons happen
         pass # Keep original logic for now

    # --- Pre-calculate Jerk Times ---
    tj1 = A_max / Jerk
    tj2 = tj1
    tj3 = D_max / Jerk
    tj4 = tj3

    # --- Step 1: Calculate minimum time and distance for ideal trapezoid (reaching V_max) ---
    Ta_ideal = V_max / A_max + A_max / Jerk
    Td_ideal = V_max / D_max + D_max / Jerk
    Dist_a_ideal = V_max / 2 * Ta_ideal # Approximation, better if adjusted for tri vs trap accel
    Dist_d_ideal = V_max / 2 * Td_ideal # Approximation, better if adjusted for tri vs trap decel

    # Refined calculation for Dist_a_ideal based on accel profile shape IF V_max were reached
    if V_max <= A_max**2 / Jerk: # Triangular accel needed to reach V_max
        _Ta_tri = 2 * math.sqrt(V_max / Jerk)
        Dist_a_ideal = Jerk * (_Ta_tri/2)**3 / 3 + V_max/2 * (_Ta_tri/2) # Complex integral result? simpler: V_max * _Ta_tri * 1/2 ? Needs check.
        # Let's stick to Dist = AvgVel * Time for simplicity in MIN DIST calc.
        # Avg Vel during ramp = V_max / 2. So Dist_a = V_max/2 * Ta
        # Ta_ideal was already calculated considering the shape. Let's use the general formula.
        Ta_ideal = V_max / A_max + A_max / Jerk # Use the general form
        Dist_a_ideal = V_max / 2 * Ta_ideal # Recalculate using general time

    else: # Trapezoidal accel needed to reach V_max
        Ta_ideal = V_max / A_max + A_max / Jerk
        Dist_a_ideal = V_max / 2 * Ta_ideal

    # Refined calculation for Dist_d_ideal based on decel profile shape IF V_max were reached
    if V_max <= D_max**2 / Jerk: # Triangular decel needed from V_max
        _Td_tri = 2 * math.sqrt(V_max / Jerk)
        Td_ideal = V_max / D_max + D_max / Jerk # Use the general form
        Dist_d_ideal = V_max / 2 * Td_ideal # Recalculate using general time
    else: # Trapezoidal decel needed from V_max
        Td_ideal = V_max / D_max + D_max / Jerk
        Dist_d_ideal = V_max / 2 * Td_ideal

    Dist_min = Dist_a_ideal + Dist_d_ideal

    # --- Step 2: Determine Profile Type and Calculate Actual Times ---
    results = {
        "error": None,
        "inputs": {"V_max": V_max, "A_max": A_max, "D_max": D_max, "Jerk": Jerk, "Distance": Distance},
        "reached_max_v": None,
        "v_peak": None,
        "Ta": 0.0, "Tv": 0.0, "Td": 0.0, "T_total": 0.0,
        "tj1": 0.0, "ta": 0.0, "tj2": 0.0,
        "tj3": 0.0, "td": 0.0, "tj4": 0.0,
    }

    if Distance >= Dist_min and V_max > 0: # Added V_max > 0 check
        results["reached_max_v"] = True
        results["v_peak"] = V_max
        Ta = Ta_ideal
        Td = Td_ideal
        Dist_v = Distance - Dist_min
        Tv = Dist_v / V_max
        T_total = Ta + Tv + Td
        results["Ta"] = Ta
        results["Tv"] = Tv
        results["Td"] = Td
        results["T_total"] = T_total

    else: # Case 2: Triangular Velocity Profile or Zero Velocity target edge case
        results["reached_max_v"] = False
        results["Tv"] = 0.0
        if Distance <= 1e-9: # Handle zero or near-zero distance
             results["v_peak"] = 0.0
             # All times are zero, handled by initial dict values
             # T_total is already 0.0
        else:
            # Solve quadratic equation for V_peak
            a_quad = (A_max + D_max) / (2 * A_max * D_max) if (A_max > 0 and D_max > 0) else 0
            b_quad = (A_max + D_max) / (2 * Jerk) if Jerk > 0 else 0
            c_quad = -Distance

            if a_quad == 0 or b_quad == 0: # Avoid division by zero if A/D/J were problematic (should be caught earlier)
                 results["error"] = "Calculation error: Invalid coefficients for V_peak."
                 return results

            discriminant = b_quad**2 - 4 * a_quad * c_quad

            if discriminant < -1e-9: # Allow for small numerical errors
                results["error"] = f"Calculation error: Negative discriminant ({discriminant}). Motion might be impossible with these parameters."
                return results
            discriminant = max(0, discriminant) # Ensure non-negative

            v_peak = (-b_quad + math.sqrt(discriminant)) / (2 * a_quad)
            results["v_peak"] = v_peak

            # Recalculate Ta and Td based on the actual V_peak
            if v_peak <= A_max**2 / Jerk:
                 Ta = 2 * math.sqrt(v_peak / Jerk) if v_peak > 0 else 0
            else:
                 Ta = v_peak / A_max + A_max / Jerk

            if v_peak <= D_max**2 / Jerk:
                 Td = 2 * math.sqrt(v_peak / Jerk) if v_peak > 0 else 0
            else:
                 Td = v_peak / D_max + D_max / Jerk

            T_total = Ta + Td
            results["Ta"] = Ta
            results["Td"] = Td
            results["T_total"] = T_total


    # --- Step 3: Calculate Sub-Phase Times (Jerk/Constant Accel/Decel) ---
    Ta_final = results["Ta"]
    Td_final = results["Td"]
    v_peak_final = results["v_peak"] if results["v_peak"] is not None else 0.0

    # Accel sub-phases
    if Ta_final > 1e-9:
        # Check accel profile shape at actual v_peak
        if v_peak_final < A_max**2 / Jerk: # Triangular acceleration
             tj1_actual = math.sqrt(v_peak_final / Jerk) if v_peak_final > 0 else 0
             ta_actual = 0.0
             tj2_actual = tj1_actual
             # Ensure Ta consistency: Ta_final should be tj1+tj2 = 2*sqrt(v_peak/Jerk)
        else: # Trapezoidal acceleration
             tj1_actual = A_max / Jerk
             tj2_actual = A_max / Jerk
             ta_actual = Ta_final - tj1_actual - tj2_actual

        results["tj1"] = tj1_actual
        results["ta"] = max(0.0, ta_actual) # Clamp negative due to float errors
        results["tj2"] = tj2_actual
        # Simple correction if ta is near zero
        if abs(results["ta"]) < 1e-9 : results["tj2"] = results["tj1"]

    # Decel sub-phases
    if Td_final > 1e-9:
         # Check decel profile shape at actual v_peak
        if v_peak_final < D_max**2 / Jerk: # Triangular deceleration
            tj3_actual = math.sqrt(v_peak_final / Jerk) if v_peak_final > 0 else 0
            td_actual = 0.0
            tj4_actual = tj3_actual
        else: # Trapezoidal deceleration
            tj3_actual = D_max / Jerk
            tj4_actual = D_max / Jerk
            td_actual = Td_final - tj3_actual - tj4_actual

        results["tj3"] = tj3_actual
        results["td"] = max(0.0, td_actual) # Clamp negative due to float errors
        results["tj4"] = tj4_actual
        # Simple correction if td is near zero
        if abs(results["td"]) < 1e-9 : results["tj4"] = results["tj3"]

    return results

# --- NEW Plotting Function ---
def plot_scurve_profile(results, num_points=500):
    """
    Plots the Jerk, Acceleration, Velocity, and Position profiles.

    Args:
        results (dict): The dictionary returned by calculate_scurve_times.
        num_points (int): Number of points to generate for the plot curves.
    """
    if results.get("error"):
        print(f"Cannot plot profile due to calculation error: {results['error']}")
        return
    if results['T_total'] <= 1e-9:
        print("Cannot plot profile: Total time is zero or near zero.")
        return

    # Extract parameters from results
    T_total = results['T_total']
    tj1, ta, tj2 = results['tj1'], results['ta'], results['tj2']
    tv = results['Tv']
    tj3, td, tj4 = results['tj3'], results['td'], results['tj4']
    J = results['inputs']['Jerk']
    A_max = results['inputs']['A_max']
    D_max = results['inputs']['D_max'] # Note: D_max is positive value

    # Define phase end times
    t1 = tj1
    t2 = t1 + ta
    t3 = t2 + tj2  # End of Accel phase (Ta)
    t4 = t3 + tv   # End of Const Vel phase
    t5 = t4 + tj3
    t6 = t5 + td
    t7 = t6 + tj4  # End of motion (T_total)

    # Ensure T_total consistency (correct for potential float inaccuracies)
    t7 = T_total

    # Create time vector
    t = np.linspace(0, T_total, num_points)
    jerk = np.zeros_like(t)

    # Build Jerk Profile based on phases
    # Phase 1: +J
    jerk[(t >= 0) & (t < t1)] = J
    # Phase 2: 0 (Constant Accel) - Jerk is 0
    # Phase 3: -J
    jerk[(t >= t2) & (t < t3)] = -J
    # Phase 4: 0 (Constant Velocity) - Jerk is 0
    # Phase 5: -J (Start Decel)
    jerk[(t >= t4) & (t < t5)] = -J
    # Phase 6: 0 (Constant Decel) - Jerk is 0
    # Phase 7: +J (End Decel)
    jerk[(t >= t6) & (t <= t7)] = J # Use <= for the last segment

    # Integrate Jerk to get Acceleration
    accel = cumulative_trapezoid(jerk, t, initial=0)

    # Integrate Acceleration to get Velocity
    velocity = cumulative_trapezoid(accel, t, initial=0)

    # Integrate Velocity to get Position
    position = cumulative_trapezoid(velocity, t, initial=0)

    # --- Plotting ---
    fig, axs = plt.subplots(4, 1, figsize=(10, 12), sharex=True) # 4 rows, 1 column

    # Plot Position
    axs[0].plot(t, position, label='Position', color='blue')
    axs[0].set_ylabel('Position (units)')
    axs[0].grid(True)
    axs[0].set_title(f'S-Curve Motion Profile (Total Time: {T_total:.3f} s)')
    axs[0].legend()

    # Plot Velocity
    axs[1].plot(t, velocity, label='Velocity', color='green')
    axs[1].set_ylabel('Velocity (units/s)')
    axs[1].grid(True)
    # Add V_peak line
    if results['v_peak'] is not None:
      axs[1].axhline(results['v_peak'], color='r', linestyle='--', lw=1, label=f'V_peak ({results["v_peak"]:.2f})')
    axs[1].legend()


    # Plot Acceleration
    axs[2].plot(t, accel, label='Acceleration', color='red')
    axs[2].set_ylabel('Acceleration (units/s^2)')
    axs[2].grid(True)
     # Add A_max/-D_max lines if relevant
    if results['ta'] > 1e-9: axs[2].axhline( A_max, color='orange', linestyle='--', lw=1, label=f'A_max ({A_max:.2f})')
    if results['td'] > 1e-9: axs[2].axhline(-D_max, color='purple', linestyle='--', lw=1, label=f'-D_max ({-D_max:.2f})')
    axs[2].legend()

    # Plot Jerk
    axs[3].plot(t, jerk, label='Jerk', color='purple')
    axs[3].set_xlabel('Time (s)')
    axs[3].set_ylabel('Jerk (units/s^3)')
    axs[3].grid(True)
    axs[3].legend()

    # Add vertical lines for phase transitions for clarity
    phase_times = [t1, t2, t3, t4, t5, t6, t7]
    labels = ['tj1', 'ta', 'tj2', 'tv', 'tj3', 'td', 'tj4']
    current_label_time = 0
    for i, time_pt in enumerate(phase_times):
        if time_pt > current_label_time + 1e-6 and time_pt < T_total - 1e-6 : # Only draw if phase exists and not at start/end
             for ax in axs:
                ax.axvline(time_pt, color='gray', linestyle=':', lw=1)
             # Add phase label slightly after the line
             axs[0].text(time_pt + 0.01 * T_total, axs[0].get_ylim()[1]*0.9, labels[i], color='gray', rotation=90, verticalalignment='top')
        current_label_time = time_pt


    plt.tight_layout() # Adjust subplot params for a tight layout
    plt.show()


# --- Example Usage (with plotting) ---
# Scenario 1: Long move, should reach V_max (Trapezoidal Velocity)
v = 200.0  # units/s
a = 1000.0   # units/s^2
d = 1000.0   # units/s^2
j = 10000.0  # units/s^3
dist = 3000.0 # units

print("--- Scenario 1 ---")
motion_times1 = calculate_scurve_times(v, a, d, j, dist)
print(json.dumps(motion_times1, indent=4))
plot_scurve_profile(motion_times1) # Call the plotting function
print("-" * 20)

# Scenario 2: Short move, should NOT reach V_max (Triangular Velocity)
v = 200.0  # units/s
a = 1000.0   # units/s^2
d = 1000.0   # units/s^2
j = 10000.0  # units/s^3
dist = 57.0 # units

print("--- Scenario 2 ---")
motion_times2 = calculate_scurve_times(v, a, d, j, dist)
print(json.dumps(motion_times2, indent=4))
plot_scurve_profile(motion_times2)
print("-" * 20)

# Scenario 3: Very high acceleration/deceleration relative to V_max
v = 200.0  # units/s
a = 1000.0   # units/s^2
d = 1000.0   # units/s^2
j = 10000.0  # units/s^3
dist = 16.0 # units

print("--- Scenario 3 ---")
motion_times3 = calculate_scurve_times(v, a, d, j, dist)
print(json.dumps(motion_times3, indent=4))
plot_scurve_profile(motion_times3)
print("-" * 20)
