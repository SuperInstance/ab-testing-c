# ab-testing-c

C port of [ab-testing](https://github.com/SuperInstance/ab-testing) — statistical A/B testing with from-scratch math.

## What's included

- **Chi-squared test** for 2×2 contingency tables (with Yates correction)
- **Welch's t-test** for unequal-variance independent samples
- **Confidence intervals** for proportions and means
- **Utility functions**: gamma, erf, normal CDF — all from scratch, no external deps

## Build & Test

```bash
make test    # build and run tests
make clean
```

## API

```c
#include "ab_testing.h"

// Chi-squared: compare conversion rates
TestResult r = chi_squared_test(convs_a, convs_b, n_a, n_b);
// r.statistic, r.p_value, r.significant

// Welch's t-test: compare continuous metrics
TestResult r = welch_t_test(sample1, n1, sample2, n2);

// Confidence intervals
ConfidenceInterval ci = proportion_ci(500, 1000, 0.95);
```

## Why C?

For embedded experiment analysis, edge computing, and systems where you can't pull in scipy. All math implemented from scratch — Lanczos gamma, Abramowitz-Stegun erf, normal/t CDF approximations.

## License

MIT
