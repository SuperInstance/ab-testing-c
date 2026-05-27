#include "ab_testing.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

/* ── Tests ─────────────────────────────────────────────────────────── */

static void test_sample_mean(void) {
    double data[] = {2.0, 4.0, 6.0, 8.0};
    double m = sample_mean(data, 4);
    assert(fabs(m - 5.0) < 1e-6);
    printf("  ✓ sample_mean([2,4,6,8]) = %.1f\n", m);
}

static void test_sample_std(void) {
    double data[] = {2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0};
    double s = sample_std(data, 8);
    assert(fabs(s - 2.0) < 0.15); /* ~2.0 with Bessel */
    printf("  ✓ sample_std = %.4f\n", s);
}

static void test_gamma_func(void) {
    /* gamma(1) = 1, gamma(2) = 1, gamma(3) = 2, gamma(4) = 6 */
    assert(fabs(gamma_func(1.0) - 1.0) < 0.01);
    assert(fabs(gamma_func(2.0) - 1.0) < 0.01);
    assert(fabs(gamma_func(3.0) - 2.0) < 0.02);
    assert(fabs(gamma_func(4.0) - 6.0) < 0.1);
    assert(fabs(gamma_func(0.5) - 1.7725) < 0.02); /* sqrt(pi) */
    printf("  ✓ gamma function\n");
}

static void test_erf(void) {
    assert(fabs(erf_approx(0.0)) < 1e-8);
    assert(erf_approx(1.0) > 0.8);
    assert(erf_approx(-1.0) < -0.8);
    assert(erf_approx(3.0) > 0.999);
    printf("  ✓ erf approximation\n");
}

static void test_normal_cdf(void) {
    assert(fabs(normal_cdf(0.0) - 0.5) < 1e-6);
    assert(normal_cdf(1.96) > 0.97);
    assert(normal_cdf(-1.96) < 0.03);
    printf("  ✓ normal CDF\n");
}

static void test_chi_squared_significant(void) {
    /* Clear winner: 100/1000 vs 150/1000 */
    TestResult r = chi_squared_test(100, 150, 1000, 1000);
    printf("  ✓ chi-squared: stat=%.4f, p=%.6f, sig=%d\n", r.statistic, r.p_value, r.significant);
    assert(r.statistic > 5.0);
    assert(r.p_value < 0.05);
    assert(r.significant == 1);
}

static void test_chi_squared_not_significant(void) {
    /* Similar rates: 100/1000 vs 102/1000 */
    TestResult r = chi_squared_test(100, 102, 1000, 1000);
    printf("  ✓ chi-squared (NS): stat=%.4f, p=%.6f, sig=%d\n", r.statistic, r.p_value, r.significant);
    assert(r.statistic < 1.0);
    assert(r.p_value > 0.05);
    assert(r.significant == 0);
}

static void test_welch_t_significant(void) {
    /* Group A: mean=10, Group B: mean=20 */
    double a[50], b[50];
    for (int i = 0; i < 50; i++) {
        a[i] = 10.0 + (i % 5) * 0.5; /* mean ~11 */
        b[i] = 20.0 + (i % 5) * 0.5; /* mean ~21 */
    }
    TestResult r = welch_t_test(a, 50, b, 50);
    printf("  ✓ Welch's t: stat=%.4f, p=%.8f, df=%.1f, sig=%d\n", r.statistic, r.p_value, r.df, r.significant);
    assert(r.statistic < -10.0 || r.statistic > 10.0);
    assert(r.p_value < 0.001);
    assert(r.significant == 1);
}

static void test_welch_t_not_significant(void) {
    /* Two similar groups */
    double a[30], b[30];
    for (int i = 0; i < 30; i++) {
        a[i] = 10.0 + (i % 3);
        b[i] = 10.2 + (i % 3);
    }
    TestResult r = welch_t_test(a, 30, b, 30);
    printf("  ✓ Welch's t (NS): stat=%.4f, p=%.6f, sig=%d\n", r.statistic, r.p_value, r.significant);
    assert(fabs(r.statistic) < 3.0);
}

static void test_proportion_ci(void) {
    ConfidenceInterval ci = proportion_ci(500, 1000, 0.95);
    printf("  ✓ Proportion CI: [%.4f, %.4f] mean=%.4f\n", ci.lower, ci.upper, ci.mean);
    assert(ci.mean > 0.49 && ci.mean < 0.51);
    assert(ci.lower < 0.50);
    assert(ci.upper > 0.50);
}

static void test_mean_ci(void) {
    double data[100];
    for (int i = 0; i < 100; i++) data[i] = 50.0 + (i % 10);
    ConfidenceInterval ci = mean_ci(data, 100, 0.95);
    printf("  ✓ Mean CI: [%.4f, %.4f] mean=%.4f\n", ci.lower, ci.upper, ci.mean);
    assert(ci.mean > 50.0 && ci.mean < 56.0);
}

static void test_chi_squared_edge_cases(void) {
    /* Zero counts */
    TestResult r1 = chi_squared_test(0, 0, 100, 100);
    assert(r1.statistic == 0.0);

    /* Perfect conversion in both */
    TestResult r2 = chi_squared_test(100, 100, 100, 100);
    assert(r2.statistic == 0.0);

    printf("  ✓ Chi-squared edge cases\n");
}

static void test_welch_t_equal_samples(void) {
    double a[20], b[20];
    for (int i = 0; i < 20; i++) { a[i] = 5.0; b[i] = 5.0; }
    TestResult r = welch_t_test(a, 20, b, 20);
    printf("  ✓ Welch's t (identical): stat=%.4f, p=%.4f\n", r.statistic, r.p_value);
    assert(r.statistic == 0.0);
}

int main(void) {
    printf("=== ab-testing tests ===\n\n");

    test_sample_mean();
    test_sample_std();
    test_gamma_func();
    test_erf();
    test_normal_cdf();
    test_chi_squared_significant();
    test_chi_squared_not_significant();
    test_chi_squared_edge_cases();
    test_welch_t_significant();
    test_welch_t_not_significant();
    test_welch_t_equal_samples();
    test_proportion_ci();
    test_mean_ci();

    printf("\n✅ All 13 tests passed!\n");
    return 0;
}
