#include "frinZfitting.hpp"
#include <stdexcept> // For std::runtime_error if needed, though returning struct

QuadraticFitResult fit_quadratic_to_3_points(
    const std::vector<double>& x_coords,
    const std::vector<double>& y_values) {

    QuadraticFitResult result;
    const double epsilon = 1e-9; // 浮動小数点比較のための許容誤差

    if (x_coords.size() != 3 || y_values.size() != 3) {
        result.success = false;
        result.message = "Input vectors must have exactly 3 points.";
        return result;
    }

    const double x1 = x_coords[0], y1 = y_values[0];
    const double x2 = x_coords[1], y2 = y_values[1];
    const double x3 = x_coords[2], y3 = y_values[2];

    // 係数計算時の分母となる値
    const double den12 = x1 - x2;
    const double den23 = x2 - x3;
    const double den13 = x1 - x3;

    if (std::abs(den12) < epsilon || std::abs(den23) < epsilon || std::abs(den13) < epsilon) {
        result.success = false;
        result.message = "X-coordinates are identical or too close for a unique quadratic fit.";
        return result;
    }

    // Solve for a, b, c in y = ax^2 + bx + c
    // Using Cramer's rule or substitution. For 3 points:
    // y1 = a*x1^2 + b*x1 + c
    // y2 = a*x2^2 + b*x2 + c
    // y3 = a*x3^2 + b*x3 + c
    // This can be solved systematically. A common formula for 'a' is:
    // a = ( (y1-y2)/(x1-x2) - (y2-y3)/(x2-x3) ) / (x1-x3)
    // b = (y1-y2)/(x1-x2) - a*(x1+x2)
    // c = y1 - a*x1^2 - b*x1

    // 計算の明確化と再計算を避けるための中間項
    const double term1 = (y1 - y2) / den12;
    const double term2 = (y2 - y3) / den23;


    result.a = (term1 - term2) / den13;
    result.b = term1 - result.a * (x1 + x2);
    result.c = y1 - result.a * x1 * x1 - result.b * x1;

    if (std::abs(result.a) < epsilon) { 
        result.success = false;
        result.message = "Coefficient 'a' (" + std::to_string(result.a) + ") is near zero; quadratic is degenerate. Peak cannot be determined from -b/2a.";
        return result;
    }

    // リクエストは「二次関数の最大値を返す」なので、
    // a > 0 (下に凸) の場合は、定義域全体での最大値は存在しない（頂点は最小値）。
    if (result.a > 0) {
        result.success = false;
        result.message = "係数 'a' が正です。二次関数は下に凸であり、この関数が求めるべき最大値を与える頂点は存在しません (頂点は最小値を示します)。";
        // result.peak_x = -result.b / (2.0 * result.a); // 最小値のx座標は計算可能だが、ここでは返さない
        return result;
    }


    result.peak_x = -result.b / (2.0 * result.a);
    result.success = true;
    result.message = "Quadratic fit successful.";
    return result;
}

QuadraticFitResult fit_quadratic_least_squares(
    const std::vector<double>& x_coords,
    const std::vector<double>& y_values) {

    QuadraticFitResult result;
    const double epsilon = 1e-9;
    const size_t N = x_coords.size();

    if (N < 3 || N != y_values.size()) {
        result.success = false;
        result.message = "Input vectors must have the same size, and at least 3 points for least squares fit.";
        return result;
    }

    double s0 = static_cast<double>(N);
    double s1 = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0;
    double t0 = 0.0, t1 = 0.0, t2 = 0.0;

    for (size_t i = 0; i < N; ++i) {
        double x = x_coords[i];
        double y = y_values[i];
        double x_sq = x * x;
        s1 += x;
        s2 += x_sq;
        s3 += x_sq * x;
        s4 += x_sq * x_sq;
        t0 += y;
        t1 += x * y;
        t2 += x_sq * y;
    }

    // Solve the system for c, b, a:
    // S0*c + S1*b + S2*a = T0
    // S1*c + S2*b + S3*a = T1
    // S2*c + S3*b + S4*a = T2

    // Denominator D
    double D = s0 * (s2 * s4 - s3 * s3) -
               s1 * (s1 * s4 - s2 * s3) +
               s2 * (s1 * s3 - s2 * s2);

    if (std::abs(D) < epsilon) {
        result.success = false;
        result.message = "Denominator D (" + std::to_string(D) + ") is near zero; matrix is singular or ill-conditioned for " + std::to_string(N) + "-point fit.";
        return result;
    }

    // Numerators for c, b, a
    double Dc_num = t0 * (s2 * s4 - s3 * s3) - s1 * (t1 * s4 - t2 * s3) + s2 * (t1 * s3 - t2 * s2);
    double Db_num = s0 * (t1 * s4 - t2 * s3) - t0 * (s1 * s4 - s2 * s3) + s2 * (s1 * t2 - s2 * t1);
    double Da_num = s0 * (s2 * t2 - s3 * t1) - s1 * (s1 * t2 - s2 * t1) + t0 * (s1 * s3 - s2 * s2);

    result.c = Dc_num / D;
    result.b = Db_num / D;
    result.a = Da_num / D;

    if (std::abs(result.a) < epsilon) {
        result.success = false;
        result.message = "Coefficient 'a' is near zero; quadratic is degenerate.";
        return result;
    }

    if (result.a > 0) {
        result.success = false;
        result.message = "Coefficient 'a' is positive. Quadratic is convex (opens upwards), seeking maximum.";
        return result;
    }

    result.peak_x = -result.b / (2.0 * result.a);
    result.success = true;
    result.message = std::to_string(N) + "-point least squares quadratic fit successful.";
    return result;
}
