import numpy as np


# --------------------------------------------------------------------------------------
# Physical parameters
# --------------------------------------------------------------------------------------

V_up = 5.0e-6       # Upstream tank volume [m^3]
V_down = 5.0e-6     # Downstream tank volume [m^3]

L = 0.05             # Core length [m]
D = 0.0254           # Core diameter [m]

area = np.pi * D**2 / 4.0

porosity = 0.10      # Connected porosity [-]
k = 1.0e-17          # Intrinsic permeability [m^2]
mu = 1.8e-5          # Gas viscosity [Pa s]
b = 1.0e5            # Klinkenberg coefficient [Pa]

P_up_initial = 1.2e6
P_down_initial = 1.0e6

safety_factor = 0.8
t_end = 5000.0


# --------------------------------------------------------------------------------------
# Run one mesh
# --------------------------------------------------------------------------------------

def run_simulation(N):

    # Spatial grid
    dx = L / N
    x = np.linspace(0.0, L, N + 1)

    # Effective boundary storage volumes
    V_up_eff = V_up + porosity * area * dx / 2.0
    V_down_eff = V_down + porosity * area * dx / 2.0

    # Initial pressure distribution
    P = np.full(N + 1, P_down_initial)
    P[0] = P_up_initial

    # Pressure transformation
    phi = (P + b)**2

    # Stable explicit timestep
    P_max = np.max(P)

    dt_limit = (
        porosity * mu * dx**2
        / (2.0 * k * (P_max + b))
    )

    dt = safety_factor * dt_limit

    # Number of timesteps
    num_steps = int(np.ceil(t_end / dt))

    # Initial discrete inventory
    inventory_initial = (
        V_up * P[0]
        + V_down * P[-1]
        + porosity * area * dx * (
            0.5 * P[0]
            + np.sum(P[1:-1])
            + 0.5 * P[-1]
        )
    )

    inventory_current = inventory_initial
    max_inventory_error = 0.0

    # ----------------------------------------------------------------------------------
    # Explicit time integration
    # ----------------------------------------------------------------------------------

    for n in range(num_steps):

        # Pressure-dependent interior coefficients
        r = (
            k * (P[1:-1] + b) * dt
            / (porosity * mu * dx**2)
        )

        # Pressure-dependent boundary coefficients
        r_up = (
            k * area * (P[0] + b) * dt
            / (mu * dx * V_up_eff)
        )

        r_down = (
            k * area * (P[-1] + b) * dt
            / (mu * dx * V_down_eff)
        )

        # Stability checks
        if not np.all(np.isfinite(r)):
            raise ValueError(
                f"N = {N}: non-finite interior coefficient at step {n}."
            )

        if np.any(r < 0.0) or np.any(r > 0.5):
            raise ValueError(
                f"N = {N}: interior stability limit violated at step {n}."
            )

        if not np.isfinite(r_up) or not 0.0 <= r_up <= 1.0:
            raise ValueError(
                f"N = {N}: invalid upstream boundary coefficient."
            )

        if not np.isfinite(r_down) or not 0.0 <= r_down <= 1.0:
            raise ValueError(
                f"N = {N}: invalid downstream boundary coefficient."
            )

        # New transformed-pressure state
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
                f"N = {N}: non-finite phi at step {n}."
            )

        if np.any(phi_new <= 0.0):
            raise ValueError(
                f"N = {N}: non-positive phi at step {n}."
            )

        # Recover pressure
        P_new = np.sqrt(phi_new) - b

        # Accept timestep
        phi = phi_new
        P = P_new

        # Current discrete inventory
        inventory_current = (
            V_up * P[0]
            + V_down * P[-1]
            + porosity * area * dx * (
                0.5 * P[0]
                + np.sum(P[1:-1])
                + 0.5 * P[-1]
            )
        )

        inventory_error = (
            inventory_current - inventory_initial
        ) / inventory_initial

        max_inventory_error = max(
            max_inventory_error,
            abs(inventory_error)
        )

    # ----------------------------------------------------------------------------------
    # Equilibrium diagnostics
    # ----------------------------------------------------------------------------------

    P_eq_analytical = (
        V_up * P_up_initial
        + (V_down + porosity * area * L) * P_down_initial
    ) / (
        V_up
        + V_down
        + porosity * area * L
    )

    P_eq_discrete = (
        inventory_initial
        / (V_up + V_down + porosity * area * L)
    )

    P_eq_numerical = 0.5 * (P[0] + P[-1])

    final_inventory_error = (
        inventory_current - inventory_initial
    ) / inventory_initial



    return {
        "N": N,
        "dx": dx,
        "dt": dt,
        "num_steps": num_steps,
        "P_eq_analytical": P_eq_analytical,
        "P_eq_discrete": P_eq_discrete,
        "P_eq_numerical": P_eq_numerical,
        "error_analytical": P_eq_numerical - P_eq_analytical,
        "error_discrete": P_eq_numerical - P_eq_discrete,
        "final_inventory_error": final_inventory_error,
        "max_inventory_error": max_inventory_error,
    }


# --------------------------------------------------------------------------------------
# Grid-refinement study
# --------------------------------------------------------------------------------------

mesh_sizes = [10, 20, 40, 80]

results = []

for N in mesh_sizes:

    print(f"Running N = {N}...")

    result = run_simulation(N)

    results.append(result)


# --------------------------------------------------------------------------------------
# Print refinement table
# --------------------------------------------------------------------------------------

print()
print(
    f"{'N':>5} "
    f"{'dx [m]':>12} "
    f"{'dt [s]':>12} "
    f"{'steps':>10} "
    f"{'Peq num [MPa]':>15} "
    f"{'Err anal [kPa]':>15} "
    f"{'Err disc [kPa]':>15} "
    f"{'Max inv err':>15}"
)

print("-" * 105)

for result in results:

    print(
        f"{result['N']:5d} "
        f"{result['dx']:12.6e} "
        f"{result['dt']:12.6e} "
        f"{result['num_steps']:10d} "
        f"{result['P_eq_numerical'] / 1e6:15.6f} "
        f"{result['error_analytical'] / 1e3:15.6f} "
        f"{result['error_discrete'] / 1e3:15.6f} "
        f"{result['max_inventory_error']:15.6e}"
    )



# --------------------------------------------------------------------------------------
# Observed convergence order
# --------------------------------------------------------------------------------------
#
# For an error that behaves approximately as:
#
#     E ~ C * dx**p
#
# the observed order p is estimated from two successive meshes as:
#
#     p = log(E_coarse / E_fine) / log(dx_coarse / dx_fine)
#
# Interpretation:
#     p ~ 1  -> first-order convergence: halving dx approximately halves the error
#     p ~ 2  -> second-order convergence: halving dx approximately quarters the error
#
# In this coupled tank-core FDM, the equilibrium and inventory diagnostics
# converge at approximately first order, even though the interior centred
# second-difference stencil itself is second order in space. The overall
# observed order is influenced by the boundary treatment, initial
# discontinuity, storage representation, and time discretization.
# --------------------------------------------------------------------------------------

def observed_order(error_coarse, error_fine, dx_coarse, dx_fine):
     #Estimate convergence order p assuming error scales as E ~ C * dx**p.
    return (
        np.log(abs(error_coarse) / abs(error_fine))
        / np.log(dx_coarse / dx_fine)
    )


print()
print("Observed convergence order")
print()

print(
    f"{'Refinement':>12} "
    f"{'p analytical':>15} "
    f"{'p discrete':>15} "
    f"{'p inventory':>15}"
)

print("-" * 62)

for i in range(len(results) - 1):

    coarse = results[i]
    fine = results[i + 1]

    p_analytical = observed_order(
        coarse["error_analytical"],
        fine["error_analytical"],
        coarse["dx"],
        fine["dx"],
    )

    p_discrete = observed_order(
        coarse["error_discrete"],
        fine["error_discrete"],
        coarse["dx"],
        fine["dx"],
    )

    p_inventory = observed_order(
        coarse["max_inventory_error"],
        fine["max_inventory_error"],
        coarse["dx"],
        fine["dx"],
    )

    refinement_label = (
        f"{coarse['N']} -> {fine['N']}"
    )

    print(
        f"{refinement_label:>12} "
        f"{p_analytical:15.4f} "
        f"{p_discrete:15.4f} "
        f"{p_inventory:15.4f}"
    )