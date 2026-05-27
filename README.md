# ab-testing-c

C port of [ab-testing](https://github.com/SuperInstance/ab-testing) — chi-squared test, Welch's t-test, confidence intervals, and statistical significance. No external dependencies.

## What It Does

- **Chi-squared test**: Two-proportion test with Yates correction, p-value via incomplete gamma function
- **Welch's t-test**: Unequal-variance two-sample t-test with Welch–Satterthwaite df
- **Confidence intervals**: Wilson-score CI for proportions, t-distribution CI for means
- **Math primitives**: Normal CDF (Abramowitz & Stegun), gamma function (Lanczos), incomplete beta (Lentz continued fraction)

## Building

```bash
make lib    # static library
make test   # build and run tests
```

## License

MIT
