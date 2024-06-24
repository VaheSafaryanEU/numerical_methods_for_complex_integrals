import numpy as np
import scipy.integrate as integrate

def contour_integral_trapezoidal(f, z_points):
    """
    Evaluate the contour integral of f along a given contour using the trapezoidal rule.
    
    Parameters:
    f : callable
        The integrand function f(z).
    z_points : array-like
        A sequence of complex points defining the contour.

    Returns:
    complex
        The approximate value of the contour integral.
    """
    n = len(z_points)
    integral = 0.0 + 0.0j
    for k in range(n - 1):
        dz = z_points[k + 1] - z_points[k]
        integral += (f(z_points[k]) + f(z_points[k + 1])) * dz / 2
    return integral

def contour_integral_simpson(f, z_points):
    """
    Evaluate the contour integral of f along a given contour using Simpson's rule.
    
    Parameters:
    f : callable
        The integrand function f(z).
    z_points : array-like
        A sequence of complex points defining the contour.

    Returns:
    complex
        The approximate value of the contour integral.
    """
    n = len(z_points)
    if n % 2 == 0:
        raise ValueError("Number of points must be odd for Simpson's rule.")
    
    integral = 0.0 + 0.0j
    for k in range(0, n - 2, 2):
        dz = z_points[k + 2] - z_points[k]
        integral += (f(z_points[k]) + 4 * f(z_points[k + 1]) + f(z_points[k + 2])) * dz / 6
    return integral

# Example usage
def f(z):
    return np.exp(z) / z

# Define the contour as a sequence of complex points (e.g., a circle around the origin)
num_points = 1000
t = np.linspace(0, 2*np.pi, num_points)
radius = 1.0
z_points = radius * np.exp(1j * t)

# Evaluate the contour integral using the trapezoidal rule
integral_trap = contour_integral_trapezoidal(f, z_points)
print(f"Trapezoidal rule: {integral_trap}")

# Evaluate the contour integral using Simpson's rule
# Ensure the number of points is odd for Simpson's rule
if num_points % 2 == 0:
    num_points += 1
    t = np.linspace(0, 2*np.pi, num_points)
    z_points = radius * np.exp(1j * t)

integral_simp = contour_integral_simpson(f, z_points)
print(f"Simpson's rule: {integral_simp}")

