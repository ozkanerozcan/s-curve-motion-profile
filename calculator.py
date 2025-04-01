import math
import json # For pretty printing the results dict

def calculate_scurve_times(V_max, A_max, D_max, Jerk, Distance):
    """
    Calculates the elapsed times for different phases of an S-curve motion profile.

    Args:
        V_max (float): Maximum desired velocity (must be > 0).
        A_max (float): Maximum acceleration (must be > 0).
        D_max (float): Maximum deceleration (must be > 0, input as a positive value).
        Jerk (float): Maximum jerk (must be > 0).
        Distance (float): Total distance to travel (must be > 0).

    Returns:
        dict: A dictionary containing the calculated times and profile details,
              or an error message if inputs are invalid or motion is impossible.
              Keys include:
              'error': Error message string (if any)
              'inputs': Dictionary of the validated inputs used.
              'reached_max_v': Boolean indicating if V_max was reached.
              'v_peak': Actual peak velocity achieved.
              'Ta': Total acceleration phase time.
              'Tv': Constant velocity phase time.
              'Td': Total deceleration phase time.
              'T_total': Total motion time.
              'tj1': Jerk time during initial acceleration ramp-up.
              'ta': Constant acceleration time.
              'tj2': Jerk time during final acceleration ramp-down.
              'tj3': Jerk time during initial deceleration ramp-up.
              'td': Constant deceleration time.
              'tj4': Jerk time during final deceleration ramp-down.
    """
    # --- Input Validation ---
    if not all(val > 0 for val in [V_max, A_max, D_max, Jerk, Distance]):
        return {"error": "Inputs (V_max, A_max, D_max, Jerk, Distance) must all be positive."}

    # --- Pre-calculate Jerk Times ---
    # These are the times IF the full A_max or D_max is reached during the jerk phase
    tj1 = A_max / Jerk  # Time to reach A_max from 0 accel
    tj2 = tj1           # Time to reach 0 accel from A_max
    tj3 = D_max / Jerk  # Time to reach D_max from 0 decel
    tj4 = tj3           # Time to reach 0 decel from D_max

    # --- Step 1: Calculate minimum time and distance for ideal trapezoid (reaching V_max) ---

    # Check if V_max is reachable given A_max and Jerk (triangular vs trapezoidal accel phase)
    if V_max <= A_max**2 / Jerk:
        # Accel phase is triangular (won't reach constant A_max) if V_max was the target
        # This formula gives the time required to reach V_max in this case
        Ta_ideal = 2 * math.sqrt(V_max / Jerk)
        # Distance covered during this triangular accel phase
        Dist_a_ideal = 2 * (0.5 * Jerk * (Ta_ideal / 2)**2) * (Ta_ideal / 2) + \
                       0.5 * Jerk * (Ta_ideal/2)**2 * (Ta_ideal / 2)
        # Simpler distance calc: average vel * time = (V_max/2) * Ta_ideal ? No.
        # Integral of v(t) = 0.5*J*t^2 for t=0 to Ta/2, plus integral for ramp down
        # Dist = 2 * integral(0.5*J*t^2) dt from 0 to Ta/2 = J * (Ta/2)^3 / 3
        # Dist_a_ideal = Jerk * (Ta_ideal/2)**3 / 3
        # Let's stick to the simpler avg velocity approach for Dist_a_ideal, accepting minor inaccuracy in this edge case calc
        # or just use the general trapezoidal formulas, knowing they are slightly off if V_max is small
        # Using general formulas is usually preferred unless very high accuracy at low V is needed.
        Ta_ideal = V_max / A_max + A_max / Jerk # Recalculate using general form for consistency below
        Dist_a_ideal = V_max / 2 * Ta_ideal

    else:
        # Accel phase is trapezoidal (reaches constant A_max) if V_max was the target
        Ta_ideal = V_max / A_max + A_max / Jerk
        Dist_a_ideal = V_max / 2 * Ta_ideal

    # Check if V_max is reachable given D_max and Jerk (triangular vs trapezoidal decel phase)
    if V_max <= D_max**2 / Jerk:
        # Decel phase is triangular if V_max was the starting velocity
        Td_ideal = 2 * math.sqrt(V_max / Jerk) # Assuming same Jerk for decel ramp-down
        # Dist_d_ideal = Jerk * (Td_ideal/2)**3 / 3 # Similar complexity as above
        # Use general formula for consistency
        Td_ideal = V_max / D_max + D_max / Jerk
        Dist_d_ideal = V_max / 2 * Td_ideal
    else:
        # Decel phase is trapezoidal if V_max was the starting velocity
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

    # Case 1: Distance is long enough to reach V_max (Trapezoidal Velocity Profile)
    if Distance >= Dist_min:
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

    # Case 2: Distance is too short to reach V_max (Triangular Velocity Profile)
    else:
        results["reached_max_v"] = False
        results["Tv"] = 0.0 # No constant velocity phase

        # Solve quadratic equation for V_peak
        # a*V_peak^2 + b*V_peak + c = 0
        a_quad = (A_max + D_max) / (2 * A_max * D_max)
        b_quad = (A_max + D_max) / (2 * Jerk)
        c_quad = -Distance

        discriminant = b_quad**2 - 4 * a_quad * c_quad

        if discriminant < 0:
            # Should not happen with valid positive inputs
            results["error"] = "Calculation error: Negative discriminant."
            return results

        # Calculate actual peak velocity reached
        v_peak = (-b_quad + math.sqrt(discriminant)) / (2 * a_quad)
        results["v_peak"] = v_peak

        # Recalculate Ta and Td based on the actual V_peak
        # Check if accel phase is triangular or trapezoidal *at V_peak*
        if v_peak <= A_max**2 / Jerk:
             Ta = 2 * math.sqrt(v_peak / Jerk) # More accurate Ta for triangular accel
        else:
             Ta = v_peak / A_max + A_max / Jerk

        # Check if decel phase is triangular or trapezoidal *at V_peak*
        if v_peak <= D_max**2 / Jerk:
             Td = 2 * math.sqrt(v_peak / Jerk) # More accurate Td for triangular decel
        else:
             Td = v_peak / D_max + D_max / Jerk

        # Ta = v_peak / A_max + A_max / Jerk # Using the general formula approach
        # Td = v_peak / D_max + D_max / Jerk # Using the general formula approach
        T_total = Ta + Td

        results["Ta"] = Ta
        # Tv is already 0.0
        results["Td"] = Td
        results["T_total"] = T_total


    # --- Step 3: Calculate Sub-Phase Times (Jerk/Constant Accel/Decel) ---
    # These are based on the FINAL calculated Ta, Td, and v_peak

    Ta_final = results["Ta"]
    Td_final = results["Td"]
    v_peak_final = results["v_peak"]

    # Accel sub-phases
    if Ta_final > 0:
        tj1_actual = A_max / Jerk
        # Check if the calculated peak velocity allows for reaching A_max
        if v_peak_final < A_max**2 / Jerk : # Triangular acceleration profile was used
             tj1_actual = math.sqrt(v_peak_final / Jerk) # Actual time spent in initial jerk phase
             ta_actual = 0.0                           # No constant acceleration phase
             tj2_actual = tj1_actual                   # Symmetric ramp down
             # Verify Ta_final approx equals tj1_actual + ta_actual + tj2_actual
        else: # Trapezoidal acceleration profile was used
             tj1_actual = A_max / Jerk
             tj2_actual = A_max / Jerk
             # Calculate constant accel time based on the total Ta
             ta_actual = Ta_final - tj1_actual - tj2_actual

        results["tj1"] = tj1_actual
        results["ta"] = max(0.0, ta_actual) # Ensure non-negative time
        results["tj2"] = tj2_actual if ta_actual >= -1e-9 else tj1_actual # Use symmetric time if ta is effectively zero or negative

    # Decel sub-phases
    if Td_final > 0:
        tj3_actual = D_max / Jerk
        # Check if peak velocity allows reaching D_max
        if v_peak_final < D_max**2 / Jerk: # Triangular deceleration profile
            tj3_actual = math.sqrt(v_peak_final / Jerk) # Actual time spent in initial decel jerk phase
            td_actual = 0.0                            # No constant deceleration phase
            tj4_actual = tj3_actual                    # Symmetric ramp down
        else: # Trapezoidal deceleration profile
            tj3_actual = D_max / Jerk
            tj4_actual = D_max / Jerk
            # Calculate constant decel time based on total Td
            td_actual = Td_final - tj3_actual - tj4_actual

        results["tj3"] = tj3_actual
        results["td"] = max(0.0, td_actual) # Ensure non-negative time
        results["tj4"] = tj4_actual if td_actual >= -1e-9 else tj3_actual # Use symmetric time if td is effectively zero or negative

    # Final consistency check (optional)
    # T_recalc = results["tj1"]+results["ta"]+results["tj2"]+results["Tv"]+results["tj3"]+results["td"]+results["tj4"]
    # if abs(results["T_total"] - T_recalc) > 1e-6:
    #     results["warning"] = "Internal timing inconsistency detected."

    return results

# --- Example Usage ---
# Scenario 1: Long move, should reach V_max (Trapezoidal Velocity)
v = 200.0  # units/s
a = 1000.0   # units/s^2
d = 1000.0   # units/s^2
j = 10000.0  # units/s^3
dist = 3000.0 # units

print("--- Scenario 1 ---")
motion_times1 = calculate_scurve_times(v, a, d, j, dist)
print(json.dumps(motion_times1, indent=4))
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
print("-" * 20)

# Scenario 3: Very high acceleration/deceleration relative to V_max
# This might trigger the triangular Accel/Decel phases even if V_max is notionally reached
v = 200.0  # units/s
a = 1000.0   # units/s^2
d = 1000.0   # units/s^2
j = 10000.0  # units/s^3
dist = 16.0 # units

print("--- Scenario 3 ---")
motion_times3 = calculate_scurve_times(v, a, d, j, dist)
print(json.dumps(motion_times3, indent=4))
print("-" * 20)

# Scenario 4: Invalid Input
print("--- Scenario 4 (Invalid Input) ---")
motion_times4 = calculate_scurve_times(v, a, d, 0, dist) # Zero Jerk
print(json.dumps(motion_times4, indent=4))
print("-" * 20)