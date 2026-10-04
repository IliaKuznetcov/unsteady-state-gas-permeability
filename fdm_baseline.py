import numpy as np

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
# Basic verification
# --------------------------------------------------------------------------------------

print(f"Core length: {L:.6f} m")
print(f"Core diameter: {D:.6f} m")
print(f"Core area: {area:.6f} m^2")

print()

print(f"Number of intervals: {N}")
print(f"Number of nodes: {len(x)}")
print(f"Grid spacing dx: {dx:.06e} m")

print()

print("Initial pressure array:")
print(P)

print()
print(f"Upstream pressure: {P[0] / 1.e6:.3f} MPa")
print(f"Downstream pressure: {P[-1] / 1e6:.3f} MPa")

# --------------------------------------------------------------------------------------
# Verify transformation
# --------------------------------------------------------------------------------------

transformation_error = np.max(np.abs(P_recovered - P))

print()
print("Initial transformed pressure array:")
print(phi)

print()
print(f"Maximum transformation error: {transformation_error:.6e} Pa")


# Verify the corresponding maximum diffusion coefficient

print()
print(f"Maximum pressure: {P_max / 1.e6:.3f} MPa")
print(f"Stability limit dt: {dt_limit:.6e} s")
print(f"Selected dt: {dt:.6e} s")
print(f"Global r bound: {r_max_bound:.6f}")