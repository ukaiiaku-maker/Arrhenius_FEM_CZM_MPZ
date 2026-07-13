#!/usr/bin/env python3
"""Numerical check of the minimum-Gumbel moments implied by Exp(1) hazard action."""
import math
import numpy as np

rng = np.random.default_rng(12345)
mu = 20.0
beta = 1.7
H = rng.exponential(1.0, size=500_000)
K = mu + beta * np.log(H)
mean_expected = mu - 0.5772156649015329 * beta
std_expected = math.pi / math.sqrt(6.0) * beta
assert abs(K.mean() - mean_expected) < 0.015
assert abs(K.std(ddof=1) - std_expected) < 0.015
print("PASS: Exp(1) hazard action maps to the expected minimum-type Gumbel moments")
