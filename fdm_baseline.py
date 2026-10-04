import numpy as np

# --------------------------------------------------------------------------------------
# Core geometry
# --------------------------------------------------------------------------------------

L = 0.05 # Core Length [m]
D = 0.0254 # Core diameter [m]

area = np.pi * D**2 / 4

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
