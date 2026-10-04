import numpy as np
import matplotlib.pyplot as plt

# Tank volumes
V_up = 5.0e-6 # Upstream tank volume [m^3]
V_down = 5.0e-6 # Downstream tank volume [m^3]


# --------------------------------------------------------------------------------------
# Core geometry and gas properties
# --------------------------------------------------------------------------------------

L = 0.05 # Core Length [m]
D = 0.0254 # Core diameter [m]

area = np.pi * D**2 / 4

porosity = 0.10 # Connected porosity [-]
k = 1.0e-17 # Intrinsic permeability [m^2]
mu = 1.8e-5 # Gas viscosity [Pa * s]

# --------------------------------------------------------------------------------------
# Spatial grid
# --------------------------------------------------------------------------------------

N = 40
dx = L / N

x = np.linspace(0.0, L, N + 1)

# --------------------------------------------------------------------------------------
# Effective boundary storage volumes
# --------------------------------------------------------------------------------------

V_up_eff = V_up + porosity * area * dx / 2.0
V_down_eff = V_down + porosity * area * dx / 2.0


# --------------------------------------------------------------------------------------
# Initial pressures
# --------------------------------------------------------------------------------------

P_up_initial = 1.2e6 # Upstream pressure [Pa]
P_down_initial = 1.0e6 # Downstream pressure [Pa]

# --------------------------------------------------------------------------------------
# Initial pressure distribution
# --------------------------------------------------------------------------------------

P = np.full(N + 1, P_down_initial)

P[0] = P_up_initial


# --------------------------------------------------------------------------------------
# Pressure transformation
# phi = (P + b)^2
# --------------------------------------------------------------------------------------
b = 1.0e5 # Klinkenberg coefficient [Pa]

phi = (P + b)**2

P_recovered = np.sqrt(phi) - b

# --------------------------------------------------------------------------------------
# Explicit time-step stability limit
# --------------------------------------------------------------------------------------

safety_factor = 0.8

P_max = np.max(P)

dt_limit = (
    porosity * mu * dx**2
    / (2.0 * k * (P_max + b))
)

dt = safety_factor * dt_limit

# --------------------------------------------------------------------------------------
# Simulation time
# --------------------------------------------------------------------------------------

t_end = 5000.0  # Total simulation time [s]

num_steps = int(np.ceil(t_end / dt))

# --------------------------------------------------------------------------------------
# Solution history
# --------------------------------------------------------------------------------------

time = np.zeros(num_steps + 1)
P_up_history = np.zeros(num_steps + 1)
P_down_history = np.zeros(num_steps + 1)

P_up_history[0] = P[0]
P_down_history[0] = P[-1]

# Verify the corresponding maximum diffusion coefficient
r_max_bound = (
    k * (P_max + b) * dt
    / (porosity * mu * dx**2)
)

# --------------------------------------------------------------------------------------
# Explicit time integration
# --------------------------------------------------------------------------------------

for n in range(num_steps):

    # Recalculate pressure-dependent interior coefficients
    r = (
        k * (P[1:-1] + b) * dt
        / (porosity * mu * dx**2)
    )

    # Recalculate pressure-dependent boundary coefficients
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
        raise ValueError("Non-finite interior FDM coefficient detected.")

    if np.any(r < 0.0) or np.any(r > 0.5):
        raise ValueError(
            f"Interior stability limit violated at step {n}: "
            f"min(r) = {np.min(r):.6f}, max(r) = {np.max(r):.6f}"
        )

    if not np.isfinite(r_up) or not 0.0 <= r_up <= 1.0:
        raise ValueError(
            f"Upstream boundary coefficient invalid at step {n}: "
            f"r_up = {r_up:.6f}"
        )

    if not np.isfinite(r_down) or not 0.0 <= r_down <= 1.0:
        raise ValueError(
            f"Downstream boundary coefficient invalid at step {n}: "
            f"r_down = {r_down:.6f}"
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

    # Validate new state
    if not np.all(np.isfinite(phi_new)):
        raise ValueError(f"Non-finite phi detected at step {n}.")

    if np.any(phi_new <= 0.0):
        raise ValueError(f"Non-positive phi detected at step {n}.")

    # Recover pressure
    P_new = np.sqrt(phi_new) - b

    # Accept the timestep
    phi = phi_new
    P = P_new

    # Store results
    time[n + 1] = (n + 1) * dt
    P_up_history[n + 1] = P[0]
    P_down_history[n + 1] = P[-1]

# --------------------------------------------------------------------------------------
# Transient solution checks
# --------------------------------------------------------------------------------------

delta_P_history = P_up_history - P_down_history

if np.any(np.diff(P_up_history) > 0.0):
    raise ValueError("Upstream pressure increased during pulse decay.")

if np.any(np.diff(P_down_history) < 0.0):
    raise ValueError("Downstream pressure decreased during pulse decay.")

if np.any(np.diff(delta_P_history) > 0.0):
    raise ValueError("Differential pressure increased during pulse decay.")


# --------------------------------------------------------------------------------------
# Diagnostics
# --------------------------------------------------------------------------------------

print()
print(f"Number of time steps: {num_steps}")
print(f"Final time:           {time[-1]:.3f} s")

print()
print(f"Initial Pu: {P_up_history[0] / 1e6:.6f} MPa")
print(f"Final Pu:   {P_up_history[-1] / 1e6:.6f} MPa")

print()
print(f"Initial Pd: {P_down_history[0] / 1e6:.6f} MPa")
print(f"Final Pd:   {P_down_history[-1] / 1e6:.6f} MPa")

print()
print(f"Final delta P: {delta_P_history[-1] / 1e3:.3f} kPa")

print()
print("Transient monotonicity checks passed.")

# --------------------------------------------------------------------------------------
# Plot pressure histories
# --------------------------------------------------------------------------------------

plt.figure()

plt.plot(
    time,
    P_up_history / 1e6,
    label="Upstream pressure"
)

plt.plot(
    time,
    P_down_history / 1e6,
    label="Downstream pressure"
)

plt.xlabel("Time [s]")
plt.ylabel("Pressure [MPa]")
plt.title("Tank pressure histories")
plt.grid(True)
plt.legend()
plt.tight_layout()


# --------------------------------------------------------------------------------------
# Plot differential pressure
# --------------------------------------------------------------------------------------

plt.figure()

plt.plot(
    time,
    delta_P_history / 1e3
)

plt.xlabel("Time [s]")
plt.ylabel("Differential pressure [kPa]")
plt.title("Pulse-decay differential pressure")
plt.grid(True)
plt.tight_layout()

plt.show()