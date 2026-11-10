

# bayesammi 1.0.0

## Changes visible to users

- Every function keeps its name, its arguments and the names of its
  results, so existing scripts run unchanged. Results print, summarise
  and plot (`print()`, `summary()`, `plot()`), and any component
  converts to a table with
  `data.table::as.data.table(x, component = "...")`.
- `bayes_ammi()`: new arguments `.k` (number of bilinear terms, default
  2 as in Crossa et al. 2011; earlier versions fitted all min(G, E) -
  1), `.burnin`, `.thin`, `.chains`, `.priors` (`bam_priors()`),
  `.level`, `.seed` (the random number state of the session is restored)
  and `.align` (sign of every term); `.nIter` now has a default (1000).
  New components `summary` (posterior mean, SD, Monte Carlo error, HPD
  interval, effective sample size, potential scale reduction factor of
  every parameter), `scores`, `draws` (`coda::mcmc.list`), `diagnostics`
  (Gelman-Rubin, Raftery-Lewis, Heidelberger-Welch), `ls`, `settings`.
- `tau1` now holds the draws of the error variance sigma^2 (it held a
  constant precision).
- `biplots()`: HPD regions of the biplot coordinates lambda_k^(1/2)
  alpha_ik of the paper, at one or more levels (`pb = c(0.05, 0.01)`);
  new arguments `.which`, `.grid`; `biplot_data` gains `axis1_lower`,
  `axis1_upper`, `axis1_significant` and `group`.
- New `bam_priors()`, `bam_interpret()`; `plot()` of a `bayes_ammi()`
  result draws traces, histograms, the biplot of posterior means and the
  main effects with HPD intervals.
- Genotypes and environments keep their labels in all tables.

## Beyond the article

- Unbalanced trials: missing plots and whole genotype-environment cells
  are handled by data augmentation inside the sampler (earlier versions
  required balanced data).
- `bayes_ammi()`: new arguments `.model` (`"AMMI"` or `"GGE"`, the sites
  regression model of Perez-Elizalde et al. 2012), `.het` (an error
  variance for every environment; the environment scores are then drawn
  from their Bingham-von Mises-Fisher conditional) and `.waic`; new
  component `fit` (DIC, WAIC and the contribution of every plot to
  WAIC).
- Decisions: `bam_predict()` (posterior of every cell and genotype mean,
  including empty cells), `bam_winners()` (probability that each
  genotype wins in each environment, groups of environments),
  `bam_superiority()` (probabilities of superior performance and
  stability, Dias et al. 2022), `bam_ranks()` and `bam_stability()`
  (Bayesian AMMI stability value, interaction variance, yield stability
  index).
- Model choice: `bam_compare()` (DIC, WAIC and the standard error of
  WAIC differences) and `bam_ppc()` (posterior predictive checks).
- `as_draws_df()`, `as_draws()` and `as.mcmc.list()` methods, for the
  **posterior**, **bayesplot** and **coda** packages.
- `bam_app()` and `bam_app_object()`: a built-in **shiny** application
  with tabs for the interpretation, the plot, the table and the printed
  output, each with copy buttons, and the R code of the fit.
- All plots are **ggplot2** with `scales::pretty_breaks()`. `bam_save()`
  writes lossless WebP; the HTML vignettes and the application use
  lossless WebP images (about half the size of PNG).
- Bayesian stability analysis: `bayes_stability()` and
  `bam_stab_priors()` fit the Finlay-Wilkinson (1963),
  Eberhart-Russell (1966) and Shukla (1972) models by a partially
  collapsed Gibbs sampler (compiled, with an R reference giving the same
  draws), with the environmental index estimated jointly, unbalanced
  data, slopes, deviation and stability variances, the safety-first
  index of Eskridge (1990), DIC and WAIC. The decision tools,
  `bam_compare()` and `bam_ppc()` accept these fits.
- `bam_indices()`: posterior distributions of Wricke’s ecovalence,
  Shukla’s stability variance, Lin and Binns’ superiority,
  Annicchiarico’s reliability and Eskridge’s safety-first index, from
  any fit. `bam_ammi_indices()`: ASV, SIPC, EV, Za, WAAS and WAASY from
  AMMI or GGE fits. Each index has HPD intervals, probabilities of being
  best and expected ranks.
- Citations: every printed result and interpretation ends with the
  sources of its method; `bam_cite()` gives them as text or BibTeX;
  `citation("bayesammi")` lists the package and every method;
  `inst/REFERENCES.bib` has all references. The application shows the
  sources under each interpretation.
- The application gains a Stability menu (regression models, stability
  indices, AMMI indices).
- Seven new vignettes: unbalanced trials and environment-specific
  variances, the GGE model, decision tools, the Bayesian stability
  models, posterior stability indices, the application, and the
  references of the methods.

## Corrections

- `bayes_ammi()`: the singular values were computed once before the loop
  and never sampled (`lambdas1` was constant); they are now drawn every
  iteration from their truncated normal conditional.
- `bayes_ammi()`: the residuals for sigma^2 compared every cell with the
  environment of a data row instead of its own environment (with `Maiz`,
  always site 1); the within-cell variances were summed in place of the
  sums of squares; the covariance matrices of the main effects were
  square roots of variances; the environment effects used B E in place
  of B G.
- `bayes_ammi()`: the scores were drawn on the unconstrained Stiefel
  manifold and then orthogonalised; they are now von Mises-Fisher draws
  on the constrained spheres (Crossa et al. 2011, eq. 8-9).
- `bayes_ammi()`: the prior of sigma^2 was improper (0 degrees of
  freedom); it is now proper.
- `biplots()`: the contours were of the singular vectors with a fixed
  bandwidth of 0.2, and “stable” was decided by the bounding box of a
  contour.
- `Maiz` help page: it described another data set; the data have 720
  rows (9 genotypes, 20 sites, 4 replicates).

## Internal

- The sampler is written in C++ (**Rcpp**, **RcppArmadillo**) and draws
  its random numbers from R; a plain R version is kept as a reference,
  and both give the same draws from the same seed.
- The tidyverse, `rlang`, `rstiefel`, `tmvtnorm`, `mvtnorm` and `ks` are
  no longer needed; `lme4` is suggested (only `ge_model()` uses it;
  `ge_ammi()` computes the error variance in closed form). The
  diagnostics use `coda`.
- Tests against Tables 1 and B1 of Crossa et al. (2011), simulated data,
  the building blocks of the sampler and the results of version 0.3.0
  for the tables.

# bayesammi 0.3.0

- Added the `biplots()` function.

# bayesammi 0.2.0

- Updated the DESCRIPTION file to meet the requirements of CRAN.

# bayesammi 0.1.0

- First version of bayesammi.
