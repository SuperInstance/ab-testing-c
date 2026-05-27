#ifndef AB_TESTING_H
#define AB_TESTING_H

#include <stddef.h>

/* ===================================================================
 * Statistical tests — no external dependencies
 * =================================================================== */

/* --- Chi-squared test (two-proportion, Yates-corrected) --- */

typedef struct {
    double chi2;
    double p_value;
    int degrees_of_freedom;
    int significant;  /* 1 if p < alpha */
} ab_chi2_result_t;

ab_chi2_result_t ab_chi_squared_test(int ctrl_successes, int ctrl_total,
                                      int treat_successes, int treat_total,
                                      double alpha);

/* --- Welch's t-test --- */

typedef struct {
    double t_statistic;
    double p_value;
    double degrees_of_freedom;
    int significant;
} ab_ttest_result_t;

ab_ttest_result_t ab_welch_t_test(const double *control, int n_ctrl,
                                   const double *treatment, int n_treat,
                                   double alpha);

/* --- Confidence intervals --- */

typedef struct {
    double lower;
    double upper;
    double confidence;
} ab_ci_t;

/* Wilson-score CI for a proportion */
ab_ci_t ab_proportion_ci(int successes, int total, double confidence);

/* CI for a mean using t-distribution */
ab_ci_t ab_mean_ci(const double *values, int n, double confidence);

/* --- Internal helpers (exposed for testing) --- */

double ab_normal_cdf(double x);
double ab_gamma_func(double x);
double ab_lower_incomplete_gamma(double s, double x);
double ab_incomplete_beta(double a, double b, double x);
double ab_z_for_confidence(double confidence);

#endif /* AB_TESTING_H */
