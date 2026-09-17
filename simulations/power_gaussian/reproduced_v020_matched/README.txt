Current-source reproduction of Section 5.1 (Gaussian case).

Package source: svpChange 0.2.0.  Design: n = 1000, patterns none/up/updown/rand1,
jumps 0.1 to 2.0 by 0.1, 100 replicates, seed = 123, tolerance = 2 or 3
observations (round(0.0025*n)).  PELT uses the Manual 2 log(n) penalty.

To preserve the historical parameterization, SVP BIC uses subtests = "right",
SVP BIC calibrated uses subtests = "right", and the multiscale extension uses
subtests = "both" with gamma = 1.8 log(n).  The original results.csv is not
overwritten.
