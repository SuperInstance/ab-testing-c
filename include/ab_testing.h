#ifndef AB_TESTING_H
#define AB_TESTING_H

#include <stddef.h>
#include <math.h>

/* ── Statistical Results ───────────────────────────────────────────── */

/** Result from a hypothesis test (chi-squared or Welch's t-test) */
typedef struct {
    double statistic;    /**< Chi-squared or t statistic */
    double p_value;      /**< Two-tailed p-value */
    double df;           /**< Degrees of freedom */
    int significant;     /**< 1 if p < 0.05 */
} TestResult;

/** Confidence interval for a statistic */
typedef struct {
    double lower;        /**< Lower bound */
    double upper;        /**< Upper bound */
    double mean;         /**< Point estimate */
    double confidence;   /**< Confidence level (e.g. 0.95) */
} ConfidenceInterval;

/* ── Chi-Squared Test ──────────────────────────────────────────────── */

/**
 * Chi-squared goodness-of-fit test on a 2x2 contingency table.
 *
 * @param a   Successes in group 1
 * @param b   Successes in group 2
 * @param n1  Total observations in group 1 (must be > 0)
 * @param n2  Total observations in group 2 (must be > 0)
 * @return    Test result with statistic, p-value, significance
 */
TestResult chi_squared_test(int a, int b, int n1, int n2);

/* ── Welch's t-Test ────────────────────────────────────────────────── */

/**
 * Welch's t-test for two independent samples with unequal variances.
 * Handles the case of unequal sample sizes and variances.
 *
 * @param sample1  First sample array (must not be NULL, n1 >= 2)
 * @param n1       Size of first sample
 * @param sample2  Second sample array (must not be NULL, n2 >= 2)
 * @param n2       Size of second sample
 * @return         Test result with statistic, p-value, df, significance
 */
TestResult welch_t_test(
    const double *sample1, int n1,
    const double *sample2, int n2
);

/* ── Confidence Intervals ──────────────────────────────────────────── */

/**
 * Proportion confidence interval using normal approximation.
 * @param successes   Number of successes
 * @param n           Total observations (must be > 0)
 * @param confidence  Confidence level (e.g. 0.95)
 * @return            Confidence interval
 */
ConfidenceInterval proportion_ci(int successes, int n, double confidence);

/**
 * Mean confidence interval using t-distribution approximation.
 * @param sample      Data array (must not be NULL, n >= 2)
 * @param n           Sample size
 * @param confidence  Confidence level (e.g. 0.95)
 * @return            Confidence interval
 */
ConfidenceInterval mean_ci(const double *sample, int n, double confidence);

/* ── Utility Functions ─────────────────────────────────────────────── */

/**
 * Compute the arithmetic mean of a sample.
 * @param data  Data array
 * @param n     Number of elements (returns 0.0 if n <= 0 or data is NULL)
 * @return      Sample mean
 */
double sample_mean(const double *data, int n);

/**
 * Compute the sample standard deviation with Bessel's correction.
 * @param data  Data array
 * @param n     Number of elements (returns 0.0 if n <= 1 or data is NULL)
 * @return      Sample standard deviation
 */
double sample_std(const double *data, int n);

/**
 * Compute the sample variance with Bessel's correction.
 * @param data  Data array
 * @param n     Number of elements
 * @return      Sample variance
 */
double sample_variance(const double *data, int n);

/**
 * Compute the median of a sample.
 * @param data  Data array (does not modify input)
 * @param n     Number of elements
 * @return      Median value (0.0 if n <= 0)
 */
double sample_median(const double *data, int n);

/**
 * Gamma function using the Lanczos approximation (g=7).
 * Accurate to ~15 significant digits for positive reals.
 * @param x  Input value
 * @return   Γ(x)
 */
double gamma_func(double x);

/**
 * Lower incomplete gamma function (series expansion).
 * Used for computing chi-squared p-values.
 * @param s  Shape parameter
 * @param x  Upper limit of integration
 * @return   γ(s, x)
 */
double incomplete_gamma(double s, double x);

/**
 * Error function approximation (Abramowitz & Stegun 7.1.26).
 * Maximum error: 1.5×10⁻⁷
 * @param x  Input value
 * @return   erf(x)
 */
double erf_approx(double x);

/**
 * Standard normal CDF using erf.
 * @param x  Input value
 * @return   Φ(x) = 0.5 * (1 + erf(x/√2))
 */
double normal_cdf(double x);

/**
 * Normal PDF (probability density function).
 * @param x  Input value
 * @return   φ(x) = (1/√2π) * exp(-x²/2)
 */
double normal_pdf(double x);

/**
 * t-distribution CDF (approximation).
 * Uses normal approximation for df > 30.
 * @param t_val  t statistic
 * @param df     Degrees of freedom
 * @return       P(T ≤ t)
 */
double t_cdf(double t_val, double df);

/**
 * Compute Cohen's d effect size between two samples.
 * @param s1    First sample
 * @param n1    Size of first sample
 * @param s2    Second sample
 * @param n2    Size of second sample
 * @return      Cohen's d (standardized mean difference)
 */
double cohens_d(const double *s1, int n1, const double *s2, int n2);

/**
 * Compute the standard error of the mean.
 * @param data  Sample data
 * @param n     Sample size
 * @return      SEM = std / √n
 */
double standard_error(const double *data, int n);

/**
 * Compute the min of a sample.
 * @param data  Sample data
 * @param n     Sample size
 * @return      Minimum value
 */
double sample_min(const double *data, int n);

/**
 * Compute the max of a sample.
 * @param data  Sample data
 * @param n     Sample size
 * @return      Maximum value
 */
double sample_max(const double *data, int n);

/**
 * Run a benchmark timing Welch's t-test.
 * @param n  Size of each sample
 * @return   Elapsed time in seconds
 */
double benchmark_welch_t(int n);

/**
 * Run a benchmark timing chi-squared test.
 * @param n  Number of iterations
 * @return   Elapsed time in seconds
 */
double benchmark_chi_squared(int n);

#endif /* AB_TESTING_H */
