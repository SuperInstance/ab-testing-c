#include "ab_testing.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>

#define ASSERT_FEQ(a, b) assert(fabs((a) - (b)) < 1e-2)

static void test_normal_cdf(void) {
    ASSERT_FEQ(ab_normal_cdf(0.0), 0.5);
    assert(ab_normal_cdf(1.96) > 0.97);
    assert(ab_normal_cdf(-1.96) < 0.03);
    assert(ab_normal_cdf(10.0) > 0.9999);
    assert(ab_normal_cdf(-10.0) < 0.0001);
    printf("  normal_cdf: PASS\n");
}

static void test_gamma_func(void) {
    /* Γ(1) = 1, Γ(2) = 1, Γ(3) = 2, Γ(4) = 6 */
    ASSERT_FEQ(ab_gamma_func(1.0), 1.0);
    ASSERT_FEQ(ab_gamma_func(2.0), 1.0);
    ASSERT_FEQ(ab_gamma_func(3.0), 2.0);
    ASSERT_FEQ(ab_gamma_func(4.0), 6.0);
    ASSERT_FEQ(ab_gamma_func(5.0), 24.0);
    printf("  gamma_func: PASS\n");
}

static void test_chi_squared(void) {
    /* Clear significant difference: 50/100 vs 70/100 */
    ab_chi2_result_t r = ab_chi_squared_test(50, 100, 70, 100, 0.05);
    assert(r.chi2 > 0.0);
    assert(r.degrees_of_freedom == 1);
    assert(r.significant == 1);
    printf("  chi_squared_significant: PASS (chi2=%.4f, p=%.6f)\n", r.chi2, r.p_value);

    /* No difference: 50/100 vs 50/100 */
    ab_chi2_result_t r2 = ab_chi_squared_test(50, 100, 50, 100, 0.05);
    assert(r2.significant == 0);
    assert(r2.chi2 == 0.0);
    printf("  chi_squared_nosig: PASS\n");
}

static void test_welch_t(void) {
    /* Clearly different groups */
    double ctrl[] = {10.0, 11.0, 10.5, 10.2, 9.8};
    double treat[] = {15.0, 16.0, 14.5, 15.5, 15.2};
    ab_ttest_result_t r = ab_welch_t_test(ctrl, 5, treat, 5, 0.05);
    assert(r.significant == 1);
    assert(r.t_statistic > 5.0);
    assert(r.p_value < 0.01);
    printf("  welch_t_significant: PASS (t=%.4f, p=%.6f)\n", r.t_statistic, r.p_value);

    /* Same group */
    ab_ttest_result_t r2 = ab_welch_t_test(ctrl, 5, ctrl, 5, 0.05);
    assert(r2.significant == 0);
    ASSERT_FEQ(r2.t_statistic, 0.0);
    printf("  welch_t_nosig: PASS\n");
}

static void test_proportion_ci(void) {
    ab_ci_t ci = ab_proportion_ci(50, 100, 0.95);
    assert(ci.lower < 0.5);
    assert(ci.upper > 0.5);
    assert(ci.upper - ci.lower < 0.3);
    printf("  proportion_ci: PASS ([%.4f, %.4f])\n", ci.lower, ci.upper);
}

static void test_mean_ci(void) {
    double vals[] = {10.0, 11.0, 12.0, 10.5, 11.5};
    ab_ci_t ci = ab_mean_ci(vals, 5, 0.95);
    assert(ci.lower < 11.0);
    assert(ci.upper > 11.0);
    printf("  mean_ci: PASS ([%.4f, %.4f])\n", ci.lower, ci.upper);
}

static void test_z_for_confidence(void) {
    double z = ab_z_for_confidence(0.95);
    assert(z > 1.8 && z < 2.0);
    double z99 = ab_z_for_confidence(0.99);
    assert(z99 > 2.3 && z99 < 2.8);
    printf("  z_for_confidence: PASS\n");
}

int main(void) {
    printf("=== AB Testing C Tests ===\n");
    test_normal_cdf();
    test_gamma_func();
    test_chi_squared();
    test_welch_t();
    test_proportion_ci();
    test_mean_ci();
    test_z_for_confidence();
    printf("All tests passed!\n");
    return 0;
}
