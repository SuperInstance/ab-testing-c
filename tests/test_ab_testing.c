/*
 * Tests for ab-testing-c
 */
#include "ab_testing.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define ASSERT_FEQ(a, b, eps) do { \
    if (fabs((a) - (b)) > (eps)) { \
        fprintf(stderr, "FAIL %s:%d: %.10f != %.10f\n", \
                __FILE__, __LINE__, (double)(a), (double)(b)); \
        return 1; \
    } \
} while(0)

static int test_normal_cdf(void)
{
    ASSERT_FEQ(ab_normal_cdf(0.0), 0.5, 1e-8);
    ASSERT_FEQ(ab_normal_cdf(-1e9), 0.0, 1e-15);
    ASSERT_FEQ(ab_normal_cdf(1e9), 1.0, 1e-15);
    /* Known value: Φ(1.96) ≈ 0.975 */
    ASSERT_FEQ(ab_normal_cdf(1.96), 0.975, 0.01);
    /* Symmetry */
    ASSERT_FEQ(ab_normal_cdf(1.0), 1.0 - ab_normal_cdf(-1.0), 1e-12);
    printf("  PASS normal_cdf\n");
    return 0;
}

static int test_gamma_func(void)
{
    /* Γ(1) = 0! = 1 */
    ASSERT_FEQ(ab_gamma_func(1.0), 1.0, 1e-10);
    /* Γ(2) = 1! = 1 */
    ASSERT_FEQ(ab_gamma_func(2.0), 1.0, 1e-10);
    /* Γ(3) = 2! = 2 */
    ASSERT_FEQ(ab_gamma_func(3.0), 2.0, 1e-10);
    /* Γ(6) = 5! = 120 */
    ASSERT_FEQ(ab_gamma_func(6.0), 120.0, 1e-8);
    /* Γ(0.5) = √π */
    ASSERT_FEQ(ab_gamma_func(0.5), sqrt(3.14159265358979323846), 1e-8);
    printf("  PASS gamma_func\n");
    return 0;
}

static int test_chi_squared_identical_groups(void)
{
    /* Identical conversion rates → not significant */
    AbChiSquaredResult r = ab_chi_squared_test(50, 100, 50, 100, 0.05);
    assert(r.chi2 < 0.05);  /* Yates correction gives small but nonzero for identical rates when n>0 */
    assert(r.significant == 0);
    printf("  PASS chi_squared_identical_groups\n");
    return 0;
}

static int test_chi_squared_different_groups(void)
{
    /* Very different conversion rates → significant */
    AbChiSquaredResult r = ab_chi_squared_test(10, 100, 90, 100, 0.05);
    assert(r.chi2 > 0);
    assert(r.p_value < 0.05);
    assert(r.significant == 1);
    printf("  PASS chi_squared_different_groups (chi2=%.4f, p=%.6f)\n",
           r.chi2, r.p_value);
    return 0;
}

static int test_chi_squared_degrees_of_freedom(void)
{
    AbChiSquaredResult r = ab_chi_squared_test(50, 100, 60, 100, 0.05);
    assert(r.degrees_of_freedom == 1);
    printf("  PASS chi_squared_degrees_of_freedom\n");
    return 0;
}

static int test_t_test_identical_samples(void)
{
    double a[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double b[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    AbTTestResult r = ab_t_test(a, 5, b, 5, 0.05);
    ASSERT_FEQ(r.t_statistic, 0.0, 1e-10);
    assert(r.significant == 0);
    printf("  PASS t_test_identical_samples\n");
    return 0;
}

static int test_t_test_different_means(void)
{
    double a[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double b[] = {11.0, 12.0, 13.0, 14.0, 15.0};
    AbTTestResult r = ab_t_test(a, 5, b, 5, 0.05);
    assert(r.t_statistic != 0.0);
    assert(r.p_value < 0.05);
    assert(r.significant == 1);
    printf("  PASS t_test_different_means (t=%.4f, p=%.6f)\n",
           r.t_statistic, r.p_value);
    return 0;
}

static int test_t_test_welch_df(void)
{
    /* Unequal sample sizes and variances → Welch–Satterthwaite df */
    double a[] = {1.0, 2.0, 3.0};
    double b[] = {10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0};
    AbTTestResult r = ab_t_test(a, 3, b, 8, 0.05);
    assert(r.degrees_of_freedom > 0);
    assert(r.degrees_of_freedom < 9);  /* Should be fractional */
    printf("  PASS t_test_welch_df (df=%.2f)\n", r.degrees_of_freedom);
    return 0;
}

static int test_proportion_ci(void)
{
    /* 50/100 → 95% CI should be roughly [0.40, 0.60] */
    AbConfidenceInterval ci = ab_proportion_ci(50, 100, 0.95);
    assert(ci.lower > 0.35 && ci.lower < 0.45);
    assert(ci.upper > 0.55 && ci.upper < 0.65);
    assert(ci.confidence == 0.95);
    printf("  PASS proportion_ci ([%.4f, %.4f])\n", ci.lower, ci.upper);
    return 0;
}

static int test_mean_ci(void)
{
    double vals[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    AbConfidenceInterval ci = ab_mean_ci(vals, 5, 0.95);
    /* Mean = 3.0, CI should contain 3.0 */
    assert(ci.lower < 3.0 && ci.upper > 3.0);
    printf("  PASS mean_ci ([%.4f, %.4f])\n", ci.lower, ci.upper);
    return 0;
}

static int test_incomplete_beta(void)
{
    /* I_0(1,1) = 0, I_1(1,1) = 1, I_0.5(1,1) = 0.5 */
    ASSERT_FEQ(ab_incomplete_beta(1.0, 1.0, 0.0), 0.0, 1e-10);
    ASSERT_FEQ(ab_incomplete_beta(1.0, 1.0, 1.0), 1.0, 1e-10);
    ASSERT_FEQ(ab_incomplete_beta(1.0, 1.0, 0.5), 0.5, 1e-10);
    printf("  PASS incomplete_beta\n");
    return 0;
}

static int test_z_for_confidence(void)
{
    /* z for 95% confidence ≈ 1.96 */
    double z = ab_z_for_confidence(0.95);
    ASSERT_FEQ(z, 1.96, 0.15);
    /* z for 99% confidence ≈ 2.576 */
    z = ab_z_for_confidence(0.99);
    ASSERT_FEQ(z, 2.576, 0.15);
    printf("  PASS z_for_confidence (z95=%.4f)\n", z);
    return 0;
}

/* ─── Main ─────────────────────────────────────────────────── */

typedef int (*test_fn)(void);

int main(void)
{
    test_fn tests[] = {
        test_normal_cdf,
        test_gamma_func,
        test_chi_squared_identical_groups,
        test_chi_squared_different_groups,
        test_chi_squared_degrees_of_freedom,
        test_t_test_identical_samples,
        test_t_test_different_means,
        test_t_test_welch_df,
        test_proportion_ci,
        test_mean_ci,
        test_incomplete_beta,
        test_z_for_confidence,
    };
    int n = sizeof(tests) / sizeof(tests[0]);
    int failures = 0;
    printf("Running %d tests...\n", n);
    for (int i = 0; i < n; i++) {
        if (tests[i]()) failures++;
    }
    printf("\n%s: %d/%d passed\n",
           failures ? "FAIL" : "OK", n - failures, n);
    return failures;
}
