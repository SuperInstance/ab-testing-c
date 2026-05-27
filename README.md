# ab-testing-c

C99 port of the [ab-testing](https://github.com/SuperInstance/ab-testing) Python library — statistical tests for A/B experiments, implemented from scratch with no external dependencies.

Designed for embedded systems and bare-metal environments where you need statistical analysis without pulling in a math library.

## What It Does

- **Chi-squared test** — two-proportion comparison with Yates correction
- **Welch's t-test** — two-sample comparison with unequal variances
- **Wilson-score confidence interval** — for proportions
- **Mean confidence interval** — using t-distribution with Cornish-Fisher approximation
- All special functions (gamma, incomplete beta, normal CDF) implemented from scratch

## C API

```c
#include <ab_testing.h>

/* Chi-squared test: 50/100 vs 70/100 */
AbChiSquaredResult r = ab_chi_squared_test(50, 100, 70, 100, 0.05);
printf("chi2=%.4f p=%.6f significant=%d\n", r.chi2, r.p_value, r.significant);

/* Welch's t-test */
double control[]  = {1.2, 1.5, 1.3, 1.4, 1.6};
double treatment[] = {1.8, 2.1, 1.9, 2.0, 2.2};
AbTTestResult t = ab_t_test(control, 5, treatment, 5, 0.05);

/* Confidence interval for a proportion */
AbConfidenceInterval ci = ab_proportion_ci(45, 100, 0.95);
printf("95%% CI: [%.4f, %.4f]\n", ci.lower, ci.upper);
```

## Build

```sh
make        # builds libab_testing.a
make test   # builds and runs tests
make clean
```

No dependencies beyond C99, `libc`, and `libm`.

## Port Notes

Ported from Python with these changes:
- Replaced Python's `@dataclass` with C structs
- Implemented all special functions (Lanczos gamma, Abramowitz & Stegun normal CDF, Lentz's continued fraction for incomplete beta)
- All arithmetic is `double` precision
- Welch–Satterthwaite degrees of freedom computed exactly as in the Python version
- Cornish-Fisher expansion for critical t-values

## License

MIT
