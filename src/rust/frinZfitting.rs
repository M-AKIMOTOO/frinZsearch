use std::f64;

use crate::frinZerror::FrinZError;

// C++版の QuadraticFitResult 構造体に対応
#[derive(Debug, Default)]
pub struct QuadraticFitResult {
    pub peak_x: f64, // 2次関数の頂点のx座標
    pub success: bool,
    pub a: f64, // y = ax^2 + bx + c の係数 a
    pub b: f64, // y = ax^2 + bx + c の係数 b
    pub c: f64, // y = ax^2 + bx + c の係数 c
    pub message: String,
}

/// 最小二乗法を用いてN点に2次関数 y = ax^2 + bx + c をフィッティングする
/// x_coords と y_values は同じサイズで、かつ3点以上必要。
/// 頂点のx座標 (-b / 2a) を返す。'a' が負の場合のみ成功とする。
pub fn fit_quadratic_least_squares(
    x_coords: &[f64],
    y_values: &[f64],
) -> Result<QuadraticFitResult, FrinZError> {

    let n = x_coords.len();
    let epsilon = 1e-9;

    if n < 3 || n != y_values.len() {
        return Err(FrinZError::Fit("入力ベクトルは同じサイズで、かつ最小二乗法には3点以上必要です。".to_string()));
    }

    let s0 = n as f64;
    let mut s1 = 0.0;
    let mut s2 = 0.0;
    let mut s3 = 0.0;
    let mut s4 = 0.0;
    let mut t0 = 0.0;
    let mut t1 = 0.0;
    let mut t2 = 0.0;

    for i in 0..n {
        let x = x_coords[i];
        let y = y_values[i];
        let x_sq = x * x;
        s1 += x;
        s2 += x_sq;
        s3 += x_sq * x;
        s4 += x_sq * x_sq;
        t0 += y;
        t1 += x * y;
        t2 += x_sq * y;
    }

    // 連立方程式を解くための行列式
    // S0*c + S1*b + S2*a = T0
    // S1*c + S2*b + S3*a = T1
    // S2*c + S3*b + S4*a = T2

    // Denominator D
    let d = s0 * (s2 * s4 - s3 * s3) -
            s1 * (s1 * s4 - s2 * s3) +
            s2 * (s1 * s3 - s2 * s2);

    if d.abs() < epsilon {
        return Err(FrinZError::Fit(format!("分母 D ({}) がほぼゼロです。行列が特異または条件が悪いです。", d)));
    }

    // Numerators for c, b, a (using Cramer's rule implicitly)
    let dc_num = t0 * (s2 * s4 - s3 * s3) - s1 * (t1 * s4 - t2 * s3) + s2 * (t1 * s3 - t2 * s2);
    let db_num = s0 * (t1 * s4 - t2 * s3) - t0 * (s1 * s4 - s2 * s3) + s2 * (s1 * t2 - s2 * t1);
    let da_num = s0 * (s2 * t2 - s3 * t1) - s1 * (s1 * t2 - s2 * t1) + t0 * (s1 * s3 - s2 * s2);

    let a = da_num / d;
    let b = db_num / d;
    let c = dc_num / d;

    if a.abs() < epsilon {
        return Err(FrinZError::Fit("係数 'a' がほぼゼロです。2次関数が退化しています。".to_string()));
    }

    if a > 0.0 {
        return Err(FrinZError::Fit("係数 'a' が正です。2次関数は下に凸であり、最大値は存在しません。".to_string()));
    }

    let peak_x = -b / (2.0 * a);

    Ok(QuadraticFitResult {
        peak_x,
        success: true,
        a,
        b,
        c,
        message: format!("{}-point least squares quadratic fit successful.", n),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fit_quadratic_least_squares_basic() {
        // y = -x^2 + 2x + 3 の頂点は x = 1, y = 4
        let x_coords = vec![0.0, 1.0, 2.0];
        let y_values = vec![3.0, 4.0, 3.0];

        let result = fit_quadratic_least_squares(&x_coords, &y_values);
        assert!(result.is_ok());
        let fit_result = result.unwrap();

        assert!(fit_result.success);
        assert!((fit_result.peak_x - 1.0).abs() < 1e-9);
        assert!((fit_result.a - (-1.0)).abs() < 1e-9);
        assert!((fit_result.b - 2.0).abs() < 1e-9);
        assert!((fit_result.c - 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_fit_quadratic_least_squares_more_points() {
        // y = -2x^2 + 8x + 1 の頂点は x = 2, y = 9
        let x_coords = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let y_values = vec![1.0, 7.0, 9.0, 7.0, 1.0];

        let result = fit_quadratic_least_squares(&x_coords, &y_values);
        assert!(result.is_ok());
        let fit_result = result.unwrap();

        assert!(fit_result.success);
        assert!((fit_result.peak_x - 2.0).abs() < 1e-9);
        assert!((fit_result.a - (-2.0)).abs() < 1e-9);
        assert!((fit_result.b - 8.0).abs() < 1e-9);
        assert!((fit_result.c - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_fit_quadratic_least_squares_positive_a() {
        // y = x^2 の頂点は x = 0, y = 0 (下に凸)
        let x_coords = vec![-1.0, 0.0, 1.0];
        let y_values = vec![1.0, 0.0, 1.0];

        let result = fit_quadratic_least_squares(&x_coords, &y_values);
        assert!(result.is_err());
        if let Err(FrinZError::Fit(msg)) = result {
            assert!(msg.contains("係数 'a' が正です"));
        } else {
            panic!("予期しないエラータイプ: {:?}", result);
        }
    }

    #[test]
    fn test_fit_quadratic_least_squares_degenerate() {
        // 直線の場合 (a=0)
        let x_coords = vec![0.0, 1.0, 2.0];
        let y_values = vec![0.0, 1.0, 2.0];

        let result = fit_quadratic_least_squares(&x_coords, &y_values);
        assert!(result.is_err());
        if let Err(FrinZError::Fit(msg)) = result {
            assert!(msg.contains("係数 'a' がほぼゼロです"));
        } else {
            panic!("予期しないエラータイプ: {:?}", result);
        }
    }

    #[test]
    fn test_fit_quadratic_least_squares_insufficient_points() {
        let x_coords = vec![0.0, 1.0];
        let y_values = vec![0.0, 1.0];

        let result = fit_quadratic_least_squares(&x_coords, &y_values);
        assert!(result.is_err());
        if let Err(FrinZError::Fit(msg)) = result {
            assert!(msg.contains("3点以上必要"));
        } else {
            panic!("予期しないエラータイプ: {:?}", result);
        }
    }
}
