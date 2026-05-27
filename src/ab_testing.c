/*
 * ab_testing — statistical tests from scratch (no external deps)
 *
 * C99 port of ab-testing Python library.
 * Implements chi-squared test, Welch's t-test, confidence intervals,
 * with all special functions (gamma, erf, incomplete beta) from scratch.
 */

#include "ab_testing.h"
#include <math.h>
#include <float.h>

/* ─── Standard normal CDF (Abramowitz & Stegun 26.2.17) ────── */

double ab_normal_cdf(double x)
{
    if (x < -8.0) return 0.0;
    if (x >  8.0) return 1.0;

    static const double a1 =  0.254829592;
    static const double a2 = -0.284496736;
    static const double a3 =  1.421413741;
    static const double a4 = -1.453152027;
    static const double a5 =  1.061405429;
    static const double p  =  0.3275911;

    int sign = (x >= 0) ? 1 : -1;
    double t = 1.0 / (1.0 + p * fabs(x));
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t
               * exp(-x * x / 2.0);
    return 0.5 * (1.0 + sign * y);
}

/* ─── Gamma function (Lanczos approximation) ───────────────── */

double ab_gamma_func(double x)
{
    if (x < 0.5)
        return 3.14159265358979323846 / (sin(3.14159265358979323846 * x) * ab_gamma_func(1.0 - x));

    x -= 1.0;
    static const double coefs[] = {
         0.99999999999980993,   676.5203681218851,   -1259.1392167224028,
         771.32342877765313,  -176.61502916214059,    12.507343278686905,
        -0.13857109526572012, 9.9843695780195716e-6,  1.5056327351493116e-7,
    };

    double s = coefs[0];
    for (int i = 1; i <= 8; i++)
        s += coefs[i] / (x + (double)i);

    double t = x + 7.5;
    return sqrt(2.0 * 3.14159265358979323846) * pow(t, x + 0.5) * exp(-t) * s;
}

/* ─── Lower regularised incomplete gamma P(s, x) ──────────── */

double ab_lower_incomplete_gamma(double s, double x)
{
    if (x < 0.0) return 0.0;
    if (x == 0.0) return 0.0;

    /* Series expansion */
    double total = 1.0 / s;
    double term  = 1.0 / s;
    for (int n = 1; n < 200; n++) {
        term  *= x / (s + (double)n);
        total += term;
    }
    return total * pow(x, s) * exp(-x) / ab_gamma_func(s);
}

/* ─── Regularised incomplete beta I_x(a, b) ────────────────── */

double ab_incomplete_beta(double a, double b, double x)
{
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;

    double lbeta = lgamma(a) + lgamma(b) - lgamma(a + b);
    double front = exp(a * log(x) + b * log(1.0 - x) - lbeta) / a;

    /* Lentz's continued fraction */
    double f = 1.0, c = 1.0, d;
    d = 1.0 - (a + 1.0) * x / (a + 1.0);
    if (fabs(d) < 1e-30) d = 1e-30;
    d = 1.0 / d;
    f = d;

    for (int m = 1; m <= 100; m++) {
        /* Even step */
        double num = (double)m * (b - (double)m) * x
                   / ((a + 2.0 * m - 1.0) * (a + 2.0 * m));
        d = 1.0 + num * d; if (fabs(d) < 1e-30) d = 1e-30;
        c = 1.0 + num / c; if (fabs(c) < 1e-30) c = 1e-30;
        d = 1.0 / d;
        f *= c * d;

        /* Odd step */
        num = -(a + (double)m) * (a + b + (double)m) * x
            / ((a + 2.0 * m) * (a + 2.0 * m + 1.0));
        d = 1.0 + num * d; if (fabs(d) < 1e-30) d = 1e-30;
        c = 1.0 + num / c; if (fabs(c) < 1e-30) c = 1e-30;
        d = 1.0 / d;
        double delta = c * d;
        f *= delta;

        if (fabs(delta - 1.0) < 1e-10) break;
    }

    return front * f;
}

/* ─── Z critical value via binary search ───────────────────── */

double ab_z_for_confidence(double confidence)
{
    double alpha_half = (1.0 - confidence) / 2.0;
    double target = 1.0 - alpha_half;
    double lo = 0.0, hi = 8.0;
    for (int i = 0; i < 100; i++) {
        double mid = (lo + hi) * 0.5;
        if (ab_normal_cdf(mid) < target)
            lo = mid;
        else
            hi = mid;
    }
    return (lo + hi) * 0.5;
}

/* ─── Approximate critical t-value (Cornish-Fisher) ────────── */

double ab_t_critical(double df, double confidence)
{
    double z = ab_z_for_confidence(confidence);
    double g1 = (z * z * z + z) / (4.0 * df);
    double g2 = (5.0 * z * z * z * z * z + 16.0 * z * z * z + 3.0 * z)
              / (96.0 * df * df);
    return z + g1 + g2;
}

/* ─── Chi-squared test (Yates-corrected) ───────────────────── */

AbChiSquaredResult ab_chi_squared_test(
    int ctrl_succ, int ctrl_tot,
    int treat_succ, int treat_tot,
    double alpha)
{
    AbChiSquaredResult r = {0, 0, 1, 0};

    if (ctrl_tot <= 0 || treat_tot <= 0) return r;

    double a = (double)ctrl_succ;
    double b = (double)(ctrl_tot - ctrl_succ);
    double c = (double)treat_succ;
    double d = (double)(treat_tot - treat_succ);
    double n = a + b + c + d;
    if (n == 0.0) return r;

    double numer = fabs(a * d - b * c) - n / 2.0;
    r.chi2 = n * numer * numer / ((a + b) * (c + d) * (a + c) * (b + d));
    r.degrees_of_freedom = 1;
    r.p_value = 1.0 - ab_lower_incomplete_gamma(0.5, r.chi2 / 2.0);
    r.significant = (r.p_value < alpha) ? 1 : 0;

    return r;
}

/* ─── Welch's two-sample t-test ────────────────────────────── */

AbTTestResult ab_t_test(
    const double *ctrl, size_t n1,
    const double *treat, size_t n2,
    double alpha)
{
    AbTTestResult r = {0, 0, 0, 0};

    if (n1 < 2 || n2 < 2) return r;

    /* Means */
    double m1 = 0, m2 = 0;
    for (size_t i = 0; i < n1; i++) m1 += ctrl[i];
    for (size_t i = 0; i < n2; i++) m2 += treat[i];
    m1 /= (double)n1;
    m2 /= (double)n2;

    /* Variances */
    double v1 = 0, v2 = 0;
    for (size_t i = 0; i < n1; i++) v1 += (ctrl[i] - m1) * (ctrl[i] - m1);
    for (size_t i = 0; i < n2; i++) v2 += (treat[i] - m2) * (treat[i] - m2);
    v1 /= (double)(n1 - 1);
    v2 /= (double)(n2 - 1);

    double se = sqrt(v1 / (double)n1 + v2 / (double)n2);

    if (se == 0.0) {
        r.t_statistic = (fabs(m2 - m1) > 0) ? HUGE_VAL : 0.0;
        r.p_value = (fabs(m2 - m1) > 0) ? 0.0 : 1.0;
        r.degrees_of_freedom = (double)(n1 + n2 - 2);
        r.significant = (fabs(m2 - m1) > 0) ? 1 : 0;
        return r;
    }

    r.t_statistic = (m2 - m1) / se;

    /* Welch–Satterthwaite df */
    double v1n = v1 / (double)n1, v2n = v2 / (double)n2;
    double num = (v1n + v2n) * (v1n + v2n);
    double denom = v1n * v1n / (double)(n1 - 1) + v2n * v2n / (double)(n2 - 1);
    r.degrees_of_freedom = (denom > 0) ? num / denom : (double)(n1 + n2 - 2);

    /* Two-tailed p-value */
    if (r.degrees_of_freedom > 30.0) {
        r.p_value = 2.0 * (1.0 - ab_normal_cdf(fabs(r.t_statistic)));
    } else {
        double x = r.degrees_of_freedom
                 / (r.degrees_of_freedom + r.t_statistic * r.t_statistic);
        r.p_value = ab_incomplete_beta(r.degrees_of_freedom / 2.0, 0.5, x);
    }

    r.significant = (r.p_value < alpha) ? 1 : 0;
    return r;
}

/* ─── Confidence intervals ─────────────────────────────────── */

AbConfidenceInterval ab_proportion_ci(int successes, int total, double confidence)
{
    AbConfidenceInterval r = {0, 0, confidence};
    if (total <= 0) return r;

    double p_hat = (double)successes / (double)total;
    double z = ab_z_for_confidence(confidence);
    double n = (double)total;
    double denom = 1.0 + z * z / n;
    double centre = (p_hat + z * z / (2.0 * n)) / denom;
    double spread = z * sqrt(p_hat * (1.0 - p_hat) / n + z * z / (4.0 * n * n)) / denom;

    r.lower = centre - spread;
    r.upper = centre + spread;
    if (r.lower < 0.0) r.lower = 0.0;
    if (r.upper > 1.0) r.upper = 1.0;
    return r;
}

AbConfidenceInterval ab_mean_ci(const double *values, size_t n, double confidence)
{
    AbConfidenceInterval r = {0, 0, confidence};
    if (n < 2) return r;

    double m = 0;
    for (size_t i = 0; i < n; i++) m += values[i];
    m /= (double)n;

    double var = 0;
    for (size_t i = 0; i < n; i++) var += (values[i] - m) * (values[i] - m);
    var /= (double)(n - 1);

    double se = sqrt(var / (double)n);
    double tc = ab_t_critical((double)(n - 1), confidence);
    r.lower = m - tc * se;
    r.upper = m + tc * se;
    return r;
}
