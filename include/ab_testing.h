#ifndef AB_TESTING_H
#define AB_TESTING_H

#include <stddef.h>
#include <math.h>

/* ── Statistical Results ───────────────────────────────────────────── */

typedef struct {
    double statistic;    /* chi-squared or t statistic */
    double p_value;
    double df;           /* degrees of freedom */
    int significant;     /* 1 if p < alpha */
} TestResult;

typedef struct {
    double lower;
    double upper;
    double mean;
    double confidence;
} ConfidenceInterval;

/* ── Chi-Squared Test ──────────────────────────────────────────────── */

/* Chi-squared goodness-of-fit / independence test on 2x2 contingency table.
 * a, b = counts for group 1 and 2 (successes)
 * n1, n2 = total counts for group 1 and 2
 */
TestResult chi_squared_test(int a, int b, int n1, int n2);

/* ── Welch's t-Test ────────────────────────────────────────────────── */

/* Welch's t-test for two independent samples with unequal variances. */
TestResult welch_t_test(
    const double *sample1, int n1,
    const double *sample2, int n2
);

/* ── Confidence Intervals ──────────────────────────────────────────── */

/* Proportion confidence interval (normal approximation) */
ConfidenceInterval proportion_ci(int successes, int n, double confidence);

/* Mean confidence interval (t-distribution approximation) */
ConfidenceInterval mean_ci(const double *sample, int n, double confidence);

/* ── Utility Functions ─────────────────────────────────────────────── */

/* Sample mean */
double sample_mean(const double *data, int n);

/* Sample standard deviation (using Bessel's correction) */
double sample_std(const double *data, int n);

/* Gamma function (Stirling/Lanczos approximation) */
double gamma_func(double x);

/* Incomplete gamma function for chi-squared p-value */
double incomplete_gamma(double s, double x);

/* erf approximation for t-test p-value */
double erf_approx(double x);

/* Normal CDF using erf */
double normal_cdf(double x);

/* t-distribution CDF (approximation) */
double t_cdf(double t, double df);

#endif /* AB_TESTING_H */
