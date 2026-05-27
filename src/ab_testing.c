#include "ab_testing.h"
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ===================================================================
 * Normal CDF — Abramowitz & Stegun approximation 26.2.17
 * =================================================================== */

double ab_normal_cdf(double x) {
    if (x < -8.0) return 0.0;
    if (x > 8.0) return 1.0;

    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    int sign = x >= 0 ? 1 : -1;
    double t = 1.0 / (1.0 + p * fabs(x));
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x / 2.0);
    return 0.5 * (1.0 + sign * y);
}

/* ===================================================================
 * Gamma function — Lanczos approximation
 * =================================================================== */

double ab_gamma_func(double x) {
    if (x < 0.5) {
        return M_PI / (sin(M_PI * x) * ab_gamma_func(1.0 - x));
    }
    x -= 1.0;
    int g = 7;
    double coefs[] = {
        0.99999999999980993,   676.5203681218851,  -1259.1392167224028,
        771.32342877765313,  -176.61502916214059,   12.507343278686905,
        -0.13857109526572012,  9.9843695780195716e-6, 1.5056327351493116e-7,
    };
    double s = coefs[0];
    for (int i = 1; i < g + 2; i++) {
        s += coefs[i] / (x + i);
    }
    double t_val = x + g + 0.5;
    return sqrt(2.0 * M_PI) * pow(t_val, x + 0.5) * exp(-t_val) * s;
}

/* ===================================================================
 * Lower regularised incomplete gamma P(s, x) — series expansion
 * =================================================================== */

double ab_lower_incomplete_gamma(double s, double x) {
    if (x < 0.0) return 0.0;
    if (x == 0.0) return 0.0;

    double total = 1.0 / s;
    double term = 1.0 / s;
    for (int n = 1; n < 200; n++) {
        term *= x / (s + n);
        total += term;
    }
    return total * pow(x, s) * exp(-x) / ab_gamma_func(s);
}

/* ===================================================================
 * Regularised incomplete beta I_x(a, b) — continued fraction (Lentz)
 * =================================================================== */

double ab_incomplete_beta(double a, double b, double x) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;

    double lbeta = lgamma(a) + lgamma(b) - lgamma(a + b);
    double front = exp(a * log(x) + b * log(1.0 - x) - lbeta) / a;

    double f = 1.0, c = 1.0, d;
    double num = (a + 1.0) * x;
    d = num != 0.0 ? (1.0 - num / (a + 1.0)) : 1.0;
    if (fabs(d) < 1e-30) d = 1e-30;
    d = 1.0 / d;
    f = d;

    for (int m = 1; m <= 100; m++) {
        /* Even step */
        num = m * (b - m) * x / ((a + 2.0 * m - 1.0) * (a + 2.0 * m));
        d = 1.0 + num * d;
        if (fabs(d) < 1e-30) d = 1e-30;
        c = 1.0 + num / c;
        if (fabs(c) < 1e-30) c = 1e-30;
        d = 1.0 / d;
        f *= c * d;

        /* Odd step */
        num = -(a + m) * (a + b + m) * x / ((a + 2.0 * m) * (a + 2.0 * m + 1.0));
        d = 1.0 + num * d;
        if (fabs(d) < 1e-30) d = 1e-30;
        c = 1.0 + num / c;
        if (fabs(c) < 1e-30) c = 1e-30;
        d = 1.0 / d;
        double delta = c * d;
        f *= delta;

        if (fabs(delta - 1.0) < 1e-10) break;
    }

    return front * f;
}

/* ===================================================================
 * Critical z-value via binary search
 * =================================================================== */

double ab_z_for_confidence(double confidence) {
    double alpha_half = (1.0 - confidence) / 2.0;
    double target = 1.0 - alpha_half;
    double lo = 0.0, hi = 8.0;
    for (int i = 0; i < 100; i++) {
        double mid = (lo + hi) / 2.0;
        if (ab_normal_cdf(mid) < target) lo = mid;
        else hi = mid;
    }
    return (lo + hi) / 2.0;
}

/* ===================================================================
 * Approximate critical t-value — Cornish-Fisher expansion
 * =================================================================== */

static double t_critical(double df, double confidence) {
    double z = ab_z_for_confidence(confidence);
    double g1 = (z * z * z + z) / (4.0 * df);
    double g2 = (5.0 * pow(z, 5) + 16.0 * z * z * z + 3.0 * z) / (96.0 * df * df);
    return z + g1 + g2;
}

/* ===================================================================
 * Two-tailed p-value for t-distribution
 * =================================================================== */

static double two_tailed_p(double t, double df) {
    if (df > 30.0) {
        return 2.0 * (1.0 - ab_normal_cdf(fabs(t)));
    }
    double x = df / (df + t * t);
    return ab_incomplete_beta(df / 2.0, 0.5, x);
}

/* ===================================================================
 * Chi-squared test (two-proportion, Yates-corrected)
 * =================================================================== */

ab_chi2_result_t ab_chi_squared_test(int ctrl_successes, int ctrl_total,
                                      int treat_successes, int treat_total,
                                      double alpha) {
    ab_chi2_result_t r = {0, 0, 1, 0};
    if (ctrl_total <= 0 || treat_total <= 0) return r;

    int a = ctrl_successes;
    int b = ctrl_total - ctrl_successes;
    int c = treat_successes;
    int d = treat_total - treat_successes;
    long long n = a + b + c + d;
    if (n == 0) return r;

    long long ad_minus_bc = (long long)a * d - (long long)b * c;
    double diff = fabs((double)ad_minus_bc) - n / 2.0;
    if (diff < 0.0) diff = 0.0;

    double denom = (double)(a + b) * (c + d) * (a + c) * (b + d);
    r.chi2 = n * diff * diff / denom;
    r.degrees_of_freedom = 1;
    r.p_value = 1.0 - ab_lower_incomplete_gamma(0.5, r.chi2 / 2.0);
    r.significant = r.p_value < alpha ? 1 : 0;
    return r;
}

/* ===================================================================
 * Welch's t-test
 * =================================================================== */

ab_ttest_result_t ab_welch_t_test(const double *control, int n_ctrl,
                                   const double *treatment, int n_treat,
                                   double alpha) {
    ab_ttest_result_t r = {0, 1.0, 0.0, 0};

    if (n_ctrl < 2 || n_treat < 2) return r;

    double m1 = 0.0, m2 = 0.0;
    for (int i = 0; i < n_ctrl; i++) m1 += control[i];
    for (int i = 0; i < n_treat; i++) m2 += treatment[i];
    m1 /= n_ctrl;
    m2 /= n_treat;

    double v1 = 0.0, v2 = 0.0;
    for (int i = 0; i < n_ctrl; i++) { double d = control[i] - m1; v1 += d * d; }
    for (int i = 0; i < n_treat; i++) { double d = treatment[i] - m2; v2 += d * d; }
    v1 /= (n_ctrl - 1);
    v2 /= (n_treat - 1);

    double se = sqrt(v1 / n_ctrl + v2 / n_treat);
    if (se == 0.0) {
        r.t_statistic = fabs(m2 - m1) > 0.0 ? 1.0 / 0.0 : 0.0;
        r.p_value = fabs(m2 - m1) > 0.0 ? 0.0 : 1.0;
        r.degrees_of_freedom = n_ctrl + n_treat - 2;
        r.significant = fabs(m2 - m1) > 0.0 ? 1 : 0;
        return r;
    }

    r.t_statistic = (m2 - m1) / se;

    double vn1 = v1 / n_ctrl, vn2 = v2 / n_treat;
    double num = (vn1 + vn2) * (vn1 + vn2);
    double denom = vn1 * vn1 / (n_ctrl - 1) + vn2 * vn2 / (n_treat - 1);
    r.degrees_of_freedom = denom > 0.0 ? num / denom : (n_ctrl + n_treat - 2);

    r.p_value = two_tailed_p(r.t_statistic, r.degrees_of_freedom);
    r.significant = r.p_value < alpha ? 1 : 0;
    return r;
}

/* ===================================================================
 * Wilson-score confidence interval for a proportion
 * =================================================================== */

ab_ci_t ab_proportion_ci(int successes, int total, double confidence) {
    ab_ci_t ci = {0, 0, confidence};
    if (total <= 0) return ci;
    double p_hat = (double)successes / total;
    double z = ab_z_for_confidence(confidence);
    double n = total;
    double denom = 1.0 + z * z / n;
    double centre = (p_hat + z * z / (2.0 * n)) / denom;
    double spread = z * sqrt(p_hat * (1.0 - p_hat) / n + z * z / (4.0 * n * n)) / denom;
    ci.lower = centre - spread;
    ci.upper = centre + spread;
    if (ci.lower < 0.0) ci.lower = 0.0;
    if (ci.upper > 1.0) ci.upper = 1.0;
    return ci;
}

/* ===================================================================
 * Confidence interval for a mean (t-distribution)
 * =================================================================== */

ab_ci_t ab_mean_ci(const double *values, int n, double confidence) {
    ab_ci_t ci = {0, 0, confidence};
    if (n < 2) return ci;
    double m = 0.0;
    for (int i = 0; i < n; i++) m += values[i];
    m /= n;
    double var = 0.0;
    for (int i = 0; i < n; i++) { double d = values[i] - m; var += d * d; }
    var /= (n - 1);
    double se = sqrt(var / n);
    double tc = t_critical(n - 1, confidence);
    ci.lower = m - tc * se;
    ci.upper = m + tc * se;
    return ci;
}
