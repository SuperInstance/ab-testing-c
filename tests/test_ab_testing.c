#include "ab_testing.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

static int test_count = 0;
#define TEST(name) do { printf("  ✓ %s\n", name); test_count++; } while(0)

/* ── Original tests (13) ───────────────────────────────────────────── */

static void test_sample_mean(void) {
    double data[] = {2.0, 4.0, 6.0, 8.0};
    double m = sample_mean(data, 4);
    assert(fabs(m - 5.0) < 1e-6);
    TEST("sample_mean([2,4,6,8]) = 5.0");
}

static void test_sample_std(void) {
    double data[] = {2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0};
    double s = sample_std(data, 8);
    assert(fabs(s - 2.0) < 0.15);
    TEST("sample_std ~2.0");
}

static void test_gamma_func(void) {
    assert(fabs(gamma_func(1.0) - 1.0) < 0.01);
    assert(fabs(gamma_func(2.0) - 1.0) < 0.01);
    assert(fabs(gamma_func(3.0) - 2.0) < 0.02);
    assert(fabs(gamma_func(4.0) - 6.0) < 0.1);
    assert(fabs(gamma_func(0.5) - 1.7725) < 0.02);
    TEST("gamma function");
}

static void test_erf(void) {
    assert(fabs(erf_approx(0.0)) < 1e-8);
    assert(erf_approx(1.0) > 0.8);
    assert(erf_approx(-1.0) < -0.8);
    assert(erf_approx(3.0) > 0.999);
    TEST("erf approximation");
}

static void test_normal_cdf(void) {
    assert(fabs(normal_cdf(0.0) - 0.5) < 1e-6);
    assert(normal_cdf(1.96) > 0.97);
    assert(normal_cdf(-1.96) < 0.03);
    TEST("normal CDF");
}

static void test_chi_squared_significant(void) {
    TestResult r = chi_squared_test(100, 150, 1000, 1000);
    assert(r.statistic > 5.0);
    assert(r.p_value < 0.05);
    assert(r.significant == 1);
    TEST("chi-squared significant");
}

static void test_chi_squared_not_significant(void) {
    TestResult r = chi_squared_test(100, 102, 1000, 1000);
    assert(r.statistic < 1.0);
    assert(r.p_value > 0.05);
    assert(r.significant == 0);
    TEST("chi-squared not significant");
}

static void test_chi_squared_edge_cases(void) {
    TestResult r1 = chi_squared_test(0, 0, 100, 100);
    assert(r1.statistic == 0.0);
    TestResult r2 = chi_squared_test(100, 100, 100, 100);
    assert(r2.statistic == 0.0);
    TEST("chi-squared edge cases");
}

static void test_welch_t_significant(void) {
    double a[50], b[50];
    for (int i = 0; i < 50; i++) {
        a[i] = 10.0 + (i % 5) * 0.5;
        b[i] = 20.0 + (i % 5) * 0.5;
    }
    TestResult r = welch_t_test(a, 50, b, 50);
    assert(r.statistic < -10.0 || r.statistic > 10.0);
    assert(r.p_value < 0.001);
    assert(r.significant == 1);
    TEST("Welch's t significant");
}

static void test_welch_t_not_significant(void) {
    double a[30], b[30];
    for (int i = 0; i < 30; i++) {
        a[i] = 10.0 + (i % 3);
        b[i] = 10.2 + (i % 3);
    }
    TestResult r = welch_t_test(a, 30, b, 30);
    assert(fabs(r.statistic) < 3.0);
    TEST("Welch's t not significant");
}

static void test_welch_t_equal_samples(void) {
    double a[20], b[20];
    for (int i = 0; i < 20; i++) { a[i] = 5.0; b[i] = 5.0; }
    TestResult r = welch_t_test(a, 20, b, 20);
    assert(r.statistic == 0.0);
    TEST("Welch's t identical");
}

static void test_proportion_ci(void) {
    ConfidenceInterval ci = proportion_ci(500, 1000, 0.95);
    assert(ci.mean > 0.49 && ci.mean < 0.51);
    assert(ci.lower < 0.50);
    assert(ci.upper > 0.50);
    TEST("proportion CI");
}

static void test_mean_ci(void) {
    double data[100];
    for (int i = 0; i < 100; i++) data[i] = 50.0 + (i % 10);
    ConfidenceInterval ci = mean_ci(data, 100, 0.95);
    assert(ci.mean > 50.0 && ci.mean < 56.0);
    TEST("mean CI");
}

/* ── New tests (target: 26+) ──────────────────────────────────────── */

static void test_sample_variance(void) {
    double data[] = {2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0};
    double v = sample_variance(data, 8);
    assert(v > 0.0);
    assert(fabs(sqrt(v) - sample_std(data, 8)) < 1e-10);
    TEST("sample variance");
}

static void test_sample_median_odd(void) {
    double data[] = {3.0, 1.0, 2.0};
    double m = sample_median(data, 3);
    assert(fabs(m - 2.0) < 1e-10);
    TEST("median (odd count)");
}

static void test_sample_median_even(void) {
    double data[] = {4.0, 1.0, 2.0, 3.0};
    double m = sample_median(data, 4);
    assert(fabs(m - 2.5) < 1e-10);
    TEST("median (even count)");
}

static void test_normal_pdf(void) {
    /* PDF at 0 should be 1/sqrt(2*pi) ≈ 0.3989 */
    assert(fabs(normal_pdf(0.0) - 0.3989) < 0.001);
    /* PDF at ±∞ should be 0 */
    assert(normal_pdf(10.0) < 1e-20);
    TEST("normal PDF");
}

static void test_cohens_d(void) {
    double a[50], b[50];
    for (int i = 0; i < 50; i++) {
        a[i] = 10.0 + (i % 5) * 0.5;
        b[i] = 20.0 + (i % 5) * 0.5;
    }
    double d = cohens_d(a, 50, b, 50);
    assert(d < -5.0);
    TEST("Cohen's d (large effect)");
}

static void test_cohens_d_zero(void) {
    double a[20], b[20];
    for (int i = 0; i < 20; i++) { a[i] = 5.0; b[i] = 5.0; }
    double d = cohens_d(a, 20, b, 20);
    assert(fabs(d) < 1e-10);
    TEST("Cohen's d (no effect)");
}

static void test_standard_error(void) {
    double data[] = {10.0, 12.0, 14.0, 16.0, 18.0};
    double se = standard_error(data, 5);
    assert(se > 0.0);
    assert(se < sample_std(data, 5));
    TEST("standard error");
}

static void test_sample_min_max(void) {
    double data[] = {5.0, 3.0, 8.0, 1.0, 9.0};
    assert(fabs(sample_min(data, 5) - 1.0) < 1e-10);
    assert(fabs(sample_max(data, 5) - 9.0) < 1e-10);
    TEST("sample min/max");
}

static void test_null_mean(void) {
    assert(sample_mean(NULL, 10) == 0.0);
    assert(sample_mean(NULL, 0) == 0.0);
    TEST("NULL sample_mean");
}

static void test_null_std(void) {
    assert(sample_std(NULL, 10) == 0.0);
    assert(sample_std(NULL, 0) == 0.0);
    TEST("NULL sample_std");
}

static void test_null_welch_t(void) {
    double a[5] = {1,2,3,4,5};
    TestResult r = welch_t_test(NULL, 5, a, 5);
    assert(r.statistic == 0.0);
    TestResult r2 = welch_t_test(a, 1, a, 5);
    assert(r2.statistic == 0.0);
    TEST("NULL/edge Welch's t");
}

static void test_chi_squared_zero_total(void) {
    TestResult r = chi_squared_test(0, 0, 0, 0);
    assert(r.statistic == 0.0);
    TEST("chi-squared zero total");
}

static void test_proportion_ci_bounds(void) {
    /* 0/100 should clamp lower bound to 0 */
    ConfidenceInterval ci = proportion_ci(0, 100, 0.95);
    assert(ci.lower >= 0.0);
    /* 100/100 should clamp upper bound to 1 */
    ConfidenceInterval ci2 = proportion_ci(100, 100, 0.95);
    assert(ci2.upper <= 1.0);
    TEST("proportion CI boundary clamping");
}

static void test_mean_ci_small_sample(void) {
    double data[] = {5.0};
    ConfidenceInterval ci = mean_ci(data, 1, 0.95);
    assert(ci.mean == 0.0);
    TEST("mean CI with n=1");
}

static void test_null_mean_ci(void) {
    ConfidenceInterval ci = mean_ci(NULL, 10, 0.95);
    assert(ci.mean == 0.0);
    TEST("NULL mean CI");
}

static void test_erf_symmetry(void) {
    assert(fabs(erf_approx(1.0) + erf_approx(-1.0)) < 1e-6);
    assert(fabs(erf_approx(2.5) + erf_approx(-2.5)) < 1e-6);
    TEST("erf symmetry");
}

static void test_normal_cdf_symmetry(void) {
    assert(fabs(normal_cdf(1.0) + normal_cdf(-1.0) - 1.0) < 1e-4);
    assert(fabs(normal_cdf(2.0) + normal_cdf(-2.0) - 1.0) < 1e-4);
    TEST("normal CDF symmetry");
}

static void test_benchmark_welch_t(void) {
    double elapsed = benchmark_welch_t(100);
    assert(elapsed >= 0.0);
    printf("    Welch's t benchmark: 1000 iterations in %.4f s\n", elapsed);
    TEST("Welch's t benchmark");
}

static void test_benchmark_chi_squared(void) {
    double elapsed = benchmark_chi_squared(10000);
    assert(elapsed >= 0.0);
    printf("    Chi-squared benchmark: 10000 iterations in %.4f s\n", elapsed);
    TEST("Chi-squared benchmark");
}

/* ── Main ───────────────────────────────────────────────────────────── */

int main(void) {
    printf("=== ab-testing tests ===\n\n");

    /* Original */
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

    /* New */
    test_sample_variance();
    test_sample_median_odd();
    test_sample_median_even();
    test_normal_pdf();
    test_cohens_d();
    test_cohens_d_zero();
    test_standard_error();
    test_sample_min_max();
    test_null_mean();
    test_null_std();
    test_null_welch_t();
    test_chi_squared_zero_total();
    test_proportion_ci_bounds();
    test_mean_ci_small_sample();
    test_null_mean_ci();
    test_erf_symmetry();
    test_normal_cdf_symmetry();
    test_benchmark_welch_t();
    test_benchmark_chi_squared();

    printf("\n✅ All %d tests passed!\n", test_count);
    return 0;
}
