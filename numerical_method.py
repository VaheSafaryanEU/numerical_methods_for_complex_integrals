import numpy as np

def contour_integral_trapezoidal(f, z_points):
    n = len(z_points)
    integral = 0.0 + 0.0j
    for k in range(n - 1):
        dz = z_points[k + 1] - z_points[k]
        integral += (f(z_points[k]) + f(z_points[k + 1])) * dz / 2
    dz = z_points[0] - z_points[-1]
    integral += (f(z_points[-1]) + f(z_points[0])) * dz / 2
    return integral

def contour_integral_simpson(f, z_points):
    n = len(z_points)
    if n % 2 == 0:
        raise ValueError("Number of points must be odd for Simpson's rule.")

    integral = 0.0 + 0.0j
    for k in range(0, n - 2, 2):
        dz = z_points[k + 2] - z_points[k]
        integral += (f(z_points[k]) + 4 * f(z_points[k + 1]) + f(z_points[k + 2])) * dz / 6
    dz = z_points[0] - z_points[-2]
    integral += (f(z_points[-2]) + 4 * f(z_points[-1]) + f(z_points[0])) * dz / 6
    return integral

def f(z):
    return np.exp(z) / z

num_points = 1000
t = np.linspace(0, 2*np.pi, num_points)
radius = 1.0
z_points = radius * np.exp(1j * t)

integral_trap = contour_integral_trapezoidal(f, z_points)
print(f"Trapezoidal rule: {integral_trap}")

if num_points % 2 == 0:
    num_points += 1
    t = np.linspace(0, 2*np.pi, num_points)
    z_points = radius * np.exp(1j * t)

integral_simp = contour_integral_simpson(f, z_points)
print(f"Simpson's rule: {integral_simp}")
