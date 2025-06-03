#ifndef FRINZFITTING_HPP
#define FRINZFITTING_HPP

#include <vector>
#include <cmath> // For std::abs
#include <string>

// Structure to hold the results of quadratic fitting
struct QuadraticFitResult {
    double peak_x = 0.0; // x-coordinate of the peak of the quadratic
    bool success = false;
    double a = 0.0, b = 0.0, c = 0.0; // Coefficients of y = ax^2 + bx + c
    std::string message;
};

// Fits a quadratic function y = ax^2 + bx + c to three points (x_coords[i], y_values[i])
// x_coords and y_values must both have a size of 3.
// Returns the x-coordinate of the peak (-b / 2a).
QuadraticFitResult fit_quadratic_to_3_points(
    const std::vector<double>& x_coords,
    const std::vector<double>& y_values
);

// Fits a quadratic function y = ax^2 + bx + c to N points (x_coords[i], y_values[i])
// using the method of least squares.
// x_coords and y_values must have the same size, and size must be >= 3.
// Returns the x-coordinate of the peak (-b / 2a) if 'a' is negative.
QuadraticFitResult fit_quadratic_least_squares(
    const std::vector<double>& x_coords,
    const std::vector<double>& y_values
);

#endif // FRINZFITTING_HPP
