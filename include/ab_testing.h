#ifndef AB_TESTING_H
#define AB_TESTING_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ─── Chi-squared test (two-proportion, Yates-corrected) ───── */

typedef struct {
    double chi2;
    double p_value;
    int    degrees_of_freedom;
    int    significant;  /* 1 if p < alpha, 0 otherwise */
} AbChiSquaredResult;

AbChiSquaredResult ab_chi_squared_test(
    int control_successes, int control_total,
    int treatment_successes, int treatment_total,
    double alpha);

/* ─── Welch's t-test ───────────────────────────────────────── */

typedef struct {
    double t_statistic;
    double p_value;
    double degrees_of_freedom;
    int    significant;
} AbTTestResult;

AbTTestResult ab_t_test(
    const double *control, size_t n_control,
    const double *treatment, size_t n_treatment,
    double alpha);

/* ─── Confidence intervals ─────────────────────────────────── */

typedef struct {
    double lower;
    double upper;
    double confidence;
} AbConfidenceInterval;

/* Wilson-score CI for a proportion */
AbConfidenceInterval ab_proportion_ci(
    int successes, int total, double confidence);

/* CI for a mean using t-distribution */
AbConfidenceInterval ab_mean_ci(
    const double *values, size_t n, double confidence);

/* ─── Internal math utilities (exposed for testing) ────────── */

double ab_normal_cdf(double x);
double ab_gamma_func(double x);
double ab_lower_incomplete_gamma(double s, double x);
double ab_incomplete_beta(double a, double b, double x);
double ab_z_for_confidence(double confidence);
double ab_t_critical(double df, double confidence);

#ifdef __cplusplus
}
#endif

#endif /* AB_TESTING_H */
