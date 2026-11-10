# Crossa et al. (2011), Tables 1 and B1: posterior means of the maize data. Monte Carlo estimates on both sides, so
# the comparison allows a quarter of the posterior standard deviation. The second scores are compared in size only:
# their published means depend on how often the chain switched sign (see the vignette on validation).
test_that("Tables 1 and B1 of Crossa et al. (2011) are reproduced", {
  skip_on_cran()
  pub <- utils::read.csv(system.file("extdata", "crossa2011.csv", package = "bayesammi"))
  fit <- bayes_ammi(.data = Maiz, .y = y, .gen = entry, .env = site, .rep = rep, .nIter = 6000, .k = 2,
                    .burnin = 1000, .chains = 2, .seed = 2011)
  s <- fit$summary
  # genotype and environment labels are 1..9 and 1..20, so the parameter names match the table
  m <- merge(pub, s[, list(parameter, est = mean, est_sd = sd)], by = "parameter")
  expect_equal(nrow(m), nrow(pub))
  first <- !grepl(",2\\]$", m$parameter)
  dev <- abs(m$est - m$mean) / m$sd
  expect_true(all(dev[first] < 0.25), label = paste(m$parameter[first][dev[first] >= 0.25], collapse = " "))
  expect_true(all(abs(m$est_sd - m$sd)[first] / m$sd[first] < 0.15))
  # second scores: posterior SD of the same order
  expect_true(all(m$est_sd[!first] / m$sd[!first] > 0.5 & m$est_sd[!first] / m$sd[!first] < 2))
})
