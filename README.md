This file provides the source code to reproduce the simulation results and real data analyses demonstrated in the following work:

 

1. local submit code

* [A Two-Step Projection-Based Goodness-of-Fit Test for Ultra-High Dimensional Sparse Regressions](https://arxiv.org/abs/2412.10721):
This paper proposes a two-step strategy that first constructs multiple test statistics based on projected predictors from distinct projections to mitigate the dimensionality problem. The second step employs p-value combination methods (e.g., the minimum p-value and the Fisher combination) to form the final test statistics. This projection-based approach significantly mitigates the dimensionality problem, enabling the tests to detect local alternatives converging to the null at the rate as if the predictor were univariate.

2. global submit code

* [Asymptotic Distribution-Free Tests for Ultra-high Dimensional Parametric Regressions via Projected Empirical Processes and p-value Combination](https://arxiv.org/abs/2601.00541v1):
This paper constructs a Cramer-von Mises type test based on a martingale-transformed, projected residual-marked empirical process, which renders the test asymptotically distribution-free. Furthermore, it proposes a novel hybrid test that aggregates empirical process-based and local smoothing tests using Cauchy combination, which is powerful against both low-frequency and high-frequency alternatives.
