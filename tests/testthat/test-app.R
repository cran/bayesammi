test_that("the application fits a model and fills every page", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  app <- bam_app_object()
  expect_s3_class(app, "shiny.appobj")
  shiny::testServer(app, {
    session$setInputs(source = "maiz", col_y = "y", col_gen = "entry", col_env = "site", col_rep = "rep", hb = TRUE,
                      model = "AMMI", k = 2, het = FALSE, iter = 300, burnin = 50, thin = 1, chains = 2, seed = 1,
                      level = 0.95, lambda_var = 1e8, nu = 1, fit_plot = "trace", pb = c("0.05"), which = "both",
                      stable = TRUE, bp_thin = 1, s = 0.2, sp_plot = "marginal", pr_plot = "cells")
    session$setInputs(run = 1)
    expect_s3_class(fit(), "bam_bayes_ammi")
    expect_s3_class(bp(), "bam_biplots")
    expect_match(code(), "bayes_ammi\\(")
    expect_match(code(), ".het    = FALSE", fixed = TRUE)
    session$setInputs(cmp_k = c("1", "2"), cmp_het = FALSE, compare = 1)
    expect_s3_class(cmp(), "bam_compare")
    expect_equal(nrow(cmp()$table), 2L)
    session$setInputs(st_model = "ER", st_alpha = 0.05, st_run = 1)
    expect_s3_class(stab(), "bam_bayes_stability")
    session$setInputs(ind_src = "stab", ind_prob = 0.25, ind_alpha = 0.05, wy = 50)
    expect_match(.bam_app_code_stab(stab_args(), "maiz"), "bayes_stability\\(")
  })
})

test_that("plots for the application are WebP when magick can write it", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
  f <- bayesammi:::.bam_plot_file(p, 400, 300)
  expect_true(file.exists(f$src))
  if (requireNamespace("magick", quietly = TRUE) && isTRUE(magick::magick_config()$webp)) {
    expect_equal(f$contentType, "image/webp")
    w <- tempfile(fileext = ".webp")
    bam_save(plot = p, filename = w, width = 4, height = 3, dpi = 100)
    expect_equal(magick::image_info(magick::image_read(w))$format, "WEBP")
  }
})
