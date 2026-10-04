import numpy as np

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
dx = L/N

x = np.linspace(0.0, L, N + 1)

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

# Verify the corresponding maximum diffusion coefficient
r_max_bound = (
    k * (P_max + b) *dt
    / (porosity * mu * dx**2)
)

# --------------------------------------------------------------------------------------
# Interior FDM coefficients
# --------------------------------------------------------------------------------------

r = (
    k * (P[1:-1] + b) *dt
    / (porosity * mu * dx**2)
)

# --------------------------------------------------------------------------------------
# First complete explicit FDM step
# Interior nodes + storage-consistent tank boundaries
# --------------------------------------------------------------------------------------

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

# --------------------------------------------------------------------------------------
# Storage-consistent boundary updates
# --------------------------------------------------------------------------------------

V_up_eff = V_up + porosity * area * dx / 2.0
V_down_eff = V_down + porosity * area * dx / 2.0

r_up = (
    k * area * (P[0] + b) * dt
    / (mu * dx * V_up_eff)
)

r_down = (
    k * area * (P[-1] + b) * dt
    / (mu * dx * V_up_eff)
)

P_new = np.sqrt(phi_new) - b

# --------------------------------------------------------------------------------------
# Storage-consistent boundary updates
# --------------------------------------------------------------------------------------

phi_new[0] = (
    phi[0]
    + r_up * (phi[1] - phi[0])
)

phi_new[-1] = (
    phi[-1]
    + r_down * (phi[-2] - phi[-1])
)

P_new = np.sqrt(phi_new) - b

# --------------------------------------------------------------------------------------
# Diagnostics
# --------------------------------------------------------------------------------------

print()
print(f"Effective upstream storage: {V_up_eff:.6e} m^3")
print(f"Effective downstream storage: {V_down_eff:.6e} m^3")

print()
print(f"Upstream boundary r: {r_up:.6f}")
print(f"Downstream boundary r: {r_down:.6f}")

print()
print("Pressure after one complete FDM step:")
print(P_new)

print(f"Upstream pressure:   {P[0] / 1e6:.6f} -> {P_new[0] / 1e6:.6f} MPa")
print(f"Node 1 pressure:     {P[1] / 1e6:.6f} -> {P_new[1] / 1e6:.6f} MPa")
print(f"Downstream pressure: {P[-1] / 1e6:.6f} -> {P_new[-1] / 1e6:.6f} MPa")