import numpy as np


# --------------------------------------------------------------------------------------
# Physical parameters
# --------------------------------------------------------------------------------------

V_up = 5.0e-6
V_down = 5.0e-6

L = 0.05
D = 0.0254
area = np.pi * D**2 / 4.0

porosity = 0.10
k = 1.0e-17
mu = 1.8e-5
b = 1.0e5

P_up_initial = 1.2e6
P_down_initial = 1.0e6


# --------------------------------------------------------------------------------------
# Fixed spatial grid
# --------------------------------------------------------------------------------------

N = 40

dx = L / N

V_up_eff = V_up + porosity * area * dx / 2.0
V_down_eff = V_down + porosity * area * dx / 2.0


# --------------------------------------------------------------------------------------
# Base stable timestep
# --------------------------------------------------------------------------------------

safety_factor = 0.8

dt_limit = (
    porosity * mu * dx**2
    / (2.0 * k * (P_up_initial + b))
)

dt_base = safety_factor * dt_limit


# --------------------------------------------------------------------------------------
# Common comparison time
# --------------------------------------------------------------------------------------

target_time = 500.0

base_num_steps = int(round(target_time / dt_base))

comparison_time = base_num_steps * dt_base


# --------------------------------------------------------------------------------------
# Time refinement levels
# --------------------------------------------------------------------------------------

dt_factors = [1.0, 0.5, 0.25, 0.125, 0.0625]

print(f"Fixed mesh: N = {N}")
print(f"dx = {dx:.6e} m")
print(f"Comparison time = {comparison_time:.6f} s")

print()

# --------------------------------------------------------------------------------------
# Run one timestep refinement level
# --------------------------------------------------------------------------------------

def run_simulation(dt_factor):

    dt = dt_base * dt_factor
    num_steps = int(base_num_steps / dt_factor)

    # Initial pressure
    P = np.full(N + 1, P_down_initial)
    P[0] = P_up_initial

    # Transformed pressure
    phi = (P + b)**2

    for n in range(num_steps):

        # Interior coefficients
        r = (
            k * (P[1:-1] + b) * dt
            / (porosity * mu * dx**2)
        )

        # Boundary coefficients
        r_up = (
            k * area * (P[0] + b) * dt
            / (mu * dx * V_up_eff)
        )

        r_down = (
            k * area * (P[-1] + b) * dt
            / (mu * dx * V_down_eff)
        )

        # Stability checks
        if np.any(r < 0.0) or np.any(r > 0.5):
            raise ValueError(
                f"Interior stability limit violated at step {n}."
            )

        if not 0.0 <= r_up <= 1.0:
            raise ValueError(
                f"Invalid upstream coefficient at step {n}."
            )

        if not 0.0 <= r_down <= 1.0:
            raise ValueError(
                f"Invalid downstream coefficient at step {n}."
            )

        # Explicit update
        phi_new = phi.copy()

        phi_new[1:-1] = (
            phi[1:-1]
            + r * (
                phi[2:]
                - 2.0 * phi[1:-1]
                + phi[:-2]
            )
        )

        phi_new[0] = (
            phi[0]
            + r_up * (phi[1] - phi[0])
        )

        phi_new[-1] = (
            phi[-1]
            + r_down * (phi[-2] - phi[-1])
        )

        # State validity
        if not np.all(np.isfinite(phi_new)):
            raise ValueError(
                f"Non-finite phi at step {n}."
            )

        if np.any(phi_new <= 0.0):
            raise ValueError(
                f"Non-positive phi at step {n}."
            )

        # Recover and accept pressure
        P = np.sqrt(phi_new) - b
        phi = phi_new

    return {
        "dt_factor": dt_factor,
        "dt": dt,
        "num_steps": num_steps,
        "final_time": num_steps * dt,
        "P_up": P[0],
        "P_down": P[-1],
        "delta_P": P[0] - P[-1],
    }

# --------------------------------------------------------------------------------------
# Time-refinement study
# --------------------------------------------------------------------------------------

dt_factors = [1.0, 0.5, 0.25, 0.125, 0.0625]

results = []

for dt_factor in dt_factors:

    print(f"Running dt factor = {dt_factor:.4f}...")

    result = run_simulation(dt_factor)

    results.append(result)

print()
print(
    f"{'factor':>8} "
    f"{'dt [s]':>12} "
    f"{'steps':>10} "
    f"{'Pu [MPa]':>12} "
    f"{'Pd [MPa]':>12} "
    f"{'dP [kPa]':>12}"
)

print("-" * 72)

for result in results:

    print(
        f"{result['dt_factor']:8.4f} "
        f"{result['dt']:12.6e} "
        f"{result['num_steps']:10d} "
        f"{result['P_up'] / 1e6:12.6f} "
        f"{result['P_down'] / 1e6:12.6f} "
        f"{result['delta_P'] / 1e3:12.6f}"
    )

# --------------------------------------------------------------------------------------
# Observed temporal convergence order
# --------------------------------------------------------------------------------------
#
# If the temporal error behaves as:
#
#     E(dt) ~ C * dt**p
#
# then successive solution differences can be used even when the exact
# transient solution is unknown.
#
# For timestep refinement by a factor of two:
#
#     D1 = |Q(dt)   - Q(dt/2)|
#     D2 = |Q(dt/2) - Q(dt/4)|
#
# and, in the asymptotic regime:
#
#     p = log(D1 / D2) / log(2)
#
# Forward Euler is first order in time, so p should approach 1 as dt decreases.
# --------------------------------------------------------------------------------------

def self_convergence_order(Q1, Q2, Q3):
    D1 = abs(Q1 - Q2)
    D2 = abs(Q2 - Q3)

    return np.log(D1 / D2) / np.log(2.0)


print()
print("Observed temporal convergence order")
print()

print(
    f"{'Levels':>20} "
    f"{'p Pu':>12} "
    f"{'p Pd':>12} "
    f"{'p dP':>12}"
)

print("-" * 60)

for i in range(len(results) - 2):

    coarse = results[i]
    medium = results[i + 1]
    fine = results[i + 2]

    p_up = self_convergence_order(
        coarse["P_up"],
        medium["P_up"],
        fine["P_up"],
    )

    p_down = self_convergence_order(
        coarse["P_down"],
        medium["P_down"],
        fine["P_down"],
    )

    p_delta = self_convergence_order(
        coarse["delta_P"],
        medium["delta_P"],
        fine["delta_P"],
    )

    label = (
        f"{coarse['dt_factor']:.3f} -> "
        f"{medium['dt_factor']:.3f} -> "
        f"{fine['dt_factor']:.3f}"
    )

    print(
        f"{label:>20} "
        f"{p_up:12.4f} "
        f"{p_down:12.4f} "
        f"{p_delta:12.4f}"
    )