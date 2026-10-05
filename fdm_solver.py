import numpy as np


def run_fdm(
    *,
    V_up,
    V_down,
    L,
    D,
    porosity,
    k,
    mu,
    b,
    P_up_initial,
    P_down_initial,
    N=40,
    t_end=5000.0,
    safety_factor=0.8,
):
    """
    Run the explicit storage-consistent FDM pulse-decay model.

    Returns a dictionary containing the time history, tank pressures,
    differential pressure, final pressure profile, timestep information,
    and inventory diagnostics.
    """

    area = np.pi * D**2 / 4.0

    dx = L / N
    x = np.linspace(0.0, L, N + 1)

    V_up_eff = V_up + porosity * area * dx / 2.0
    V_down_eff = V_down + porosity * area * dx / 2.0

    # Initial pressure field
    P = np.full(N + 1, P_down_initial, dtype=float)
    P[0] = P_up_initial

    phi = (P + b)**2


    # --------------------------------------------------------------------------
    # Discrete pressure-volume inventory
    # --------------------------------------------------------------------------

    def calculate_inventory(P_state):
        """Return the discrete pressure-volume inventory [Pa m^3]."""

        return (
            V_up * P_state[0]
            + V_down * P_state[-1]
            + porosity * area * dx * (
                0.5 * P_state[0]
                + np.sum(P_state[1:-1])
                + 0.5 * P_state[-1]
            )
        )

    inventory_initial = calculate_inventory(P)

    # --------------------------------------------------------------------------
    # Zero-permeability limiting case
    # --------------------------------------------------------------------------

    if k == 0.0:

        time = np.array([0.0, t_end])

        P_up_history = np.array([
            P_up_initial,
            P_up_initial,
        ])

        P_down_history = np.array([
            P_down_initial,
            P_down_initial,
        ])

        delta_P_history = P_up_history - P_down_history

        inventory_history = np.array([
            inventory_initial,
            inventory_initial,
        ])

        inventory_error_history = np.zeros(2)

        return {
            "x": x,
            "time": time,
            "P_up": P_up_history,
            "P_down": P_down_history,
            "delta_P": delta_P_history,
            "P_final": P.copy(),
            "dx": dx,
            "dt": None,
            "num_steps": 0,
            "inventory": inventory_history,
            "inventory_error": inventory_error_history,
            "max_inventory_error": 0.0,
        }

    # --------------------------------------------------------------------------
    # Explicit timestep
    # --------------------------------------------------------------------------

    P_max = np.max(P)

    dt_limit = (
        porosity * mu * dx**2
        / (2.0 * k * (P_max + b))
    )

    dt = safety_factor * dt_limit

    num_steps = int(np.ceil(t_end / dt))


    # --------------------------------------------------------------------------
    # Solution history
    # --------------------------------------------------------------------------

    time = np.zeros(num_steps + 1)
    P_up_history = np.zeros(num_steps + 1)
    P_down_history = np.zeros(num_steps + 1)
    inventory_history = np.zeros(num_steps + 1)

    P_up_history[0] = P[0]
    P_down_history[0] = P[-1]
    inventory_history[0] = inventory_initial

    # --------------------------------------------------------------------------
    # Explicit time integration
    # --------------------------------------------------------------------------

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


        # ----------------------------------------------------------------------
        # Stability checks
        # ----------------------------------------------------------------------

        if not np.all(np.isfinite(r)):
            raise ValueError(
                f"Non-finite interior FDM coefficient at step {n}."
            )

        if np.any(r < 0.0) or np.any(r > 0.5):
            raise ValueError(
                f"Interior stability limit violated at step {n}: "
                f"min(r) = {np.min(r):.6f}, "
                f"max(r) = {np.max(r):.6f}"
            )

        if not np.isfinite(r_up) or not 0.0 <= r_up <= 1.0:
            raise ValueError(
                f"Invalid upstream boundary coefficient at step {n}: "
                f"r_up = {r_up:.6f}"
            )

        if not np.isfinite(r_down) or not 0.0 <= r_down <= 1.0:
            raise ValueError(
                f"Invalid downstream boundary coefficient at step {n}: "
                f"r_down = {r_down:.6f}"
            )


        # ----------------------------------------------------------------------
        # Explicit update
        # ----------------------------------------------------------------------

        phi_new = phi.copy()

        # Interior nodes
        phi_new[1:-1] = (
            phi[1:-1]
            + r * (
                phi[2:]
                - 2.0 * phi[1:-1]
                + phi[:-2]
            )
        )

        # Upstream boundary
        phi_new[0] = (
            phi[0]
            + r_up * (phi[1] - phi[0])
        )

        # Downstream boundary
        phi_new[-1] = (
            phi[-1]
            + r_down * (phi[-2] - phi[-1])
        )


        # ----------------------------------------------------------------------
        # New-state validity checks
        # ----------------------------------------------------------------------

        if not np.all(np.isfinite(phi_new)):
            raise ValueError(
                f"Non-finite transformed pressure at step {n}."
            )

        if np.any(phi_new <= 0.0):
            raise ValueError(
                f"Non-positive transformed pressure at step {n}."
            )


        # ----------------------------------------------------------------------
        # Recover and accept pressure
        # ----------------------------------------------------------------------

        P_new = np.sqrt(phi_new) - b

        phi = phi_new
        P = P_new


        # ----------------------------------------------------------------------
        # Store results
        # ----------------------------------------------------------------------

        time[n + 1] = (n + 1) * dt
        P_up_history[n + 1] = P[0]
        P_down_history[n + 1] = P[-1]
        inventory_history[n + 1] = calculate_inventory(P)

    # --------------------------------------------------------------------------
    # Derived histories and diagnostics
    # --------------------------------------------------------------------------

    delta_P_history = P_up_history - P_down_history

    inventory_error_history = (
        inventory_history - inventory_initial
    ) / inventory_initial

    max_inventory_error = np.max(
        np.abs(inventory_error_history)
    )


    return {
    "x": x,
    "time": time,
    "P_up": P_up_history,
    "P_down": P_down_history,
    "delta_P": delta_P_history,
    "P_final": P.copy(),
    "dx": dx,
    "dt": dt,
    "num_steps": num_steps,
    "inventory": inventory_history,
    "inventory_error": inventory_error_history,
    "max_inventory_error": max_inventory_error,
}