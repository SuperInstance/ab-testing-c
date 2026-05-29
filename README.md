# ab-testing-c

Statistical A/B testing in pure C — chi-squared tests, Welch's t-test, and confidence intervals. All math from scratch, zero external dependencies beyond the C standard library.

## What This Gives You

- **Chi-squared test** for comparing two proportions (with Yates correction)
- **Welch's t-test** for comparing two independent samples with unequal variances
- **Confidence intervals** for proportions (normal approximation) and means (t-distribution)
- **Utility functions**: sample mean, standard deviation, gamma function, erf, normal CDF — all from scratch
- Runs anywhere C99 compiles — embedded systems, edge computing, kernels

## Quick Start

```c
#include "ab_testing.h"

/* Chi-squared: compare conversion rates */
TestResult r = chi_squared_test(150, 80, 1000, 1000);
if (r.significant) {
    printf("Significant difference (p=%.4f)\n", r.p_value);
}

/* Welch's t-test: compare continuous metrics */
double control[] = {1.0, 2.0, 3.0, 4.0, 5.0};
double treatment[] = {10.0, 11.0, 12.0, 13.0, 14.0};
TestResult t = welch_t_test(control, 5, treatment, 5);

/* Confidence interval for a proportion */
ConfidenceInterval ci = proportion_ci(500, 1000, 0.95);
printf("Rate: %.1f%% [%.1f%%, %.1f%%]\n", ci.mean * 100, ci.lower * 100, ci.upper * 100);
```

## API Reference

### Test Functions

| Function | Parameters | Returns |
|----------|-----------|---------|
| `chi_squared_test(a, b, n1, n2)` | Successes in each group, total counts | `TestResult` |
| `welch_t_test(s1, n1, s2, n2)` | Sample arrays and sizes | `TestResult` |

### Confidence Intervals

| Function | Parameters | Returns |
|----------|-----------|---------|
| `proportion_ci(successes, n, confidence)` | Success count, total, confidence level | `ConfidenceInterval` |
| `mean_ci(sample, n, confidence)` | Sample array, size, confidence level | `ConfidenceInterval` |

### Result Types

```c
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
```

### Utility Functions

| Function | Description |
|----------|-------------|
| `sample_mean(data, n)` | Arithmetic mean |
| `sample_std(data, n)` | Bessel-corrected standard deviation |
| `gamma_func(x)` | Lanczos approximation (g=7) |
| `erf_func(x)` | Abramowitz-Stegun approximation |
| `normal_cdf(x)` | Standard normal CDF |

## How It Fits

C port of [ab-testing-rs](https://github.com/SuperInstance/ab-testing-rs). Part of the [SuperInstance OpenConstruct](https://github.com/SuperInstance/OpenConstruct) ecosystem. Use this when you need statistical testing in environments where you can't pull in scipy or a Rust runtime.

## Testing

**31 assertions** across tests for chi-squared, Welch's t-test, proportion CI, mean CI, and utility functions.

## Installation

```bash
# Build and test
make test

# Just build
make

# Clean
make clean
```

Requires a C99 compiler and `math.h`. No other dependencies.
