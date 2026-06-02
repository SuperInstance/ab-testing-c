#define _POSIX_C_SOURCE 199309L
#include "ab_testing.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ── Utility: sample mean ──────────────────────────────────────────── */

double sample_mean(const double *data, int n) {
    if (n <= 0 || !data) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += data[i];
    return sum / n;
}

/* ── Utility: sample variance ──────────────────────────────────────── */

double sample_variance(const double *data, int n) {
    if (n <= 1 || !data) return 0.0;
    double mean = sample_mean(data, n);
    double ss = 0.0;
    for (int i = 0; i < n; i++) {
        double d = data[i] - mean;
        ss += d * d;
    }
    return ss / (n - 1);
}

/* ── Utility: sample std dev (Bessel-corrected) ────────────────────── */

double sample_std(const double *data, int n) {
    return sqrt(sample_variance(data, n));
}

/* ── Utility: sample median ────────────────────────────────────────── */

static int dbl_cmp(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

double sample_median(const double *data, int n) {
    if (n <= 0 || !data) return 0.0;
    double *sorted = (double *)malloc((size_t)n * sizeof(double));
    memcpy(sorted, data, (size_t)n * sizeof(double));
    qsort(sorted, (size_t)n, sizeof(double), dbl_cmp);
    double med = (n % 2 == 0)
        ? (sorted[n/2 - 1] + sorted[n/2]) / 2.0
        : sorted[n/2];
    free(sorted);
    return med;
}

/* ── Utility: min/max ──────────────────────────────────────────────── */

double sample_min(const double *data, int n) {
    if (n <= 0 || !data) return 0.0;
    double m = data[0];
    for (int i = 1; i < n; i++) if (data[i] < m) m = data[i];
    return m;
}

double sample_max(const double *data, int n) {
    if (n <= 0 || !data) return 0.0;
    double m = data[0];
    for (int i = 1; i < n; i++) if (data[i] > m) m = data[i];
    return m;
}

/* ── Utility: standard error ───────────────────────────────────────── */

double standard_error(const double *data, int n) {
    if (n <= 0 || !data) return 0.0;
    return sample_std(data, n) / sqrt((double)n);
}

/* ── Gamma function (Lanczos approximation, g=7) ──────────────────── */

double gamma_func(double x) {
    static const double c[] = {
        0.99999999999980993, 676.5203681218851, -1259.1392167224028,
        771.32342877765313, -176.61502916214059, 12.507343278686905,
        -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7
    };
    if (x < 0.5) {
        return 3.14159265358979323846 / (sin(3.14159265358979323846 * x) * gamma_func(1.0 - x));
    }
    x -= 1.0;
    double ag = c[0];
    for (int i = 1; i < 9; i++) {
        ag += c[i] / (x + (double)i);
    }
    double t = x + 7.5;
    return sqrt(2.0 * 3.14159265358979323846) * pow(t, x + 0.5) * exp(-t) * ag;
}

/* ── Lower incomplete gamma function (series expansion) ─────────────── */

double incomplete_gamma(double s, double x) {
    if (x < 0.0) return 0.0;
    double sum = 1.0 / s;
    double term = 1.0 / s;
    for (int n = 1; n < 200; n++) {
        term *= x / (s + n);
        sum += term;
        if (fabs(term) < fabs(sum) * 1e-12) break;
    }
    return pow(x, s) * exp(-x) * sum;
}

/* ── erf approximation (Abramowitz & Stegun 7.1.26) ─────────────────── */

double erf_approx(double x) {
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    int sign = (x >= 0) ? 1 : -1;
    x = fabs(x);
    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);
    return sign * y;
}

/* ── Normal CDF ──────────────────────────────────────────────────────── */

double normal_cdf(double x) {
    return 0.5 * (1.0 + erf_approx(x / sqrt(2.0)));
}

/* ── Normal PDF ──────────────────────────────────────────────────────── */

double normal_pdf(double x) {
    return exp(-x * x / 2.0) / sqrt(2.0 * 3.14159265358979323846);
}

/* ── t-distribution CDF (approximation) ─────────────────────────────── */

double t_cdf(double t_val, double df) {
    if (df > 30) {
        return normal_cdf(t_val);
    }
    double z = t_val * (1.0 - 1.0 / (4.0 * df)) / sqrt(1.0 + t_val * t_val / (2.0 * df));
    return normal_cdf(z);
}

/* ── Cohen's d ──────────────────────────────────────────────────────── */

double cohens_d(const double *s1, int n1, const double *s2, int n2) {
    if (n1 <= 0 || n2 <= 0 || !s1 || !s2) return 0.0;
    double m1 = sample_mean(s1, n1);
    double m2 = sample_mean(s2, n2);
    double v1 = sample_variance(s1, n1);
    double v2 = sample_variance(s2, n2);
    /* Pooled std */
    double pooled = sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (double)(n1 + n2 - 2));
    if (pooled <= 0.0) return 0.0;
    return (m1 - m2) / pooled;
}

/* ── Chi-squared test (2x2 contingency table) ──────────────────────── */

TestResult chi_squared_test(int a, int b, int n1, int n2) {
    TestResult result = {0, 0, 0, 0};
    int total = n1 + n2;
    if (total == 0 || n1 == 0 || n2 == 0) return result;

    double p1 = (double)a / n1;
    double p2 = (double)b / n2;
    double p_pool = (double)(a + b) / total;

    double obs_diff = p1 - p2;
    double se = sqrt(p_pool * (1.0 - p_pool) * (1.0 / n1 + 1.0 / n2));
    if (se == 0.0) {
        result.statistic = 0.0;
        result.p_value = 1.0;
        result.df = 1.0;
        result.significant = 0;
        return result;
    }

    result.statistic = (obs_diff * obs_diff) / (se * se);
    result.df = 1.0;

    double gamma_val = gamma_func(result.df / 2.0);
    if (gamma_val == 0.0) return result;
    double lower_gamma = incomplete_gamma(result.df / 2.0, result.statistic / 2.0);
    double cdf = lower_gamma / gamma_val;
    result.p_value = 1.0 - cdf;
    if (result.p_value < 0.0) result.p_value = 0.0;
    if (result.p_value > 1.0) result.p_value = 1.0;
    result.significant = (result.p_value < 0.05) ? 1 : 0;

    return result;
}

/* ── Welch's t-test ──────────────────────────────────────────────────── */

TestResult welch_t_test(const double *s1, int n1, const double *s2, int n2) {
    TestResult result = {0, 0, 0, 0};
    if (n1 < 2 || n2 < 2 || !s1 || !s2) return result;

    double m1 = sample_mean(s1, n1);
    double m2 = sample_mean(s2, n2);
    double v1 = sample_variance(s1, n1);
    double v2 = sample_variance(s2, n2);

    double se = sqrt(v1 / n1 + v2 / n2);
    if (se == 0.0) {
        result.statistic = 0.0;
        result.p_value = 1.0;
        result.df = 1.0;
        result.significant = 0;
        return result;
    }

    result.statistic = (m1 - m2) / se;

    double num = (v1 / n1 + v2 / n2) * (v1 / n1 + v2 / n2);
    double den = (v1 / n1) * (v1 / n1) / (n1 - 1) + (v2 / n2) * (v2 / n2) / (n2 - 1);
    result.df = (den > 0.0) ? num / den : 1.0;

    double p_one_tail = 1.0 - t_cdf(fabs(result.statistic), result.df);
    result.p_value = 2.0 * p_one_tail;
    if (result.p_value < 0.0) result.p_value = 0.0;
    if (result.p_value > 1.0) result.p_value = 1.0;
    result.significant = (result.p_value < 0.05) ? 1 : 0;

    return result;
}

/* ── Confidence Intervals ───────────────────────────────────────────── */

ConfidenceInterval proportion_ci(int successes, int n, double confidence) {
    ConfidenceInterval ci = {0, 0, 0, 0};
    if (n <= 0) return ci;
    double p = (double)successes / n;
    double alpha = 1.0 - confidence;
    double z = 1.96;
    if (fabs(alpha - 0.10) < 0.01) z = 1.645;
    if (fabs(alpha - 0.01) < 0.01) z = 2.576;
    double se = sqrt(p * (1.0 - p) / n);
    ci.mean = p;
    ci.lower = p - z * se;
    ci.upper = p + z * se;
    ci.confidence = confidence;
    if (ci.lower < 0.0) ci.lower = 0.0;
    if (ci.upper > 1.0) ci.upper = 1.0;
    return ci;
}

ConfidenceInterval mean_ci(const double *sample, int n, double confidence) {
    ConfidenceInterval ci = {0, 0, 0, 0};
    if (n < 2 || !sample) return ci;
    double m = sample_mean(sample, n);
    double s = sample_std(sample, n);
    double se_val = s / sqrt((double)n);
    double t_crit = 1.96;
    if (n < 30) t_crit = 2.0 + 2.0 / n;
    ci.mean = m;
    ci.lower = m - t_crit * se_val;
    ci.upper = m + t_crit * se_val;
    ci.confidence = confidence;
    return ci;
}

/* ── Benchmarks ─────────────────────────────────────────────────────── */

double benchmark_welch_t(int n) {
    double *a = (double *)malloc((size_t)n * sizeof(double));
    double *b = (double *)malloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) {
        a[i] = (double)i * 0.1;
        b[i] = (double)i * 0.1 + 1.0;
    }
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (int i = 0; i < 1000; i++) {
        welch_t_test(a, n, b, n);
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    free(b);
    free(a);
    return (double)(end.tv_sec - start.tv_sec) +
           (double)(end.tv_nsec - start.tv_nsec) / 1e9;
}

double benchmark_chi_squared(int n) {
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (int i = 0; i < n; i++) {
        chi_squared_test(100 + i % 50, 120 + i % 30, 1000, 1000);
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    return (double)(end.tv_sec - start.tv_sec) +
           (double)(end.tv_nsec - start.tv_nsec) / 1e9;
}
