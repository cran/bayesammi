## ----setup, include = FALSE---------------------------------------------------
# HTML figures are drawn as PNG and kept as WebP: lossless (the same pixels in about half the space of PNG), or
# lossy at quality 90 when that is less than half the size, as for traces and densities of thousands of draws. The
# PDF figures stay vector graphics. Any failure of magick (not installed, no WebP, or a broken installation whose
# system libraries are missing) keeps the PNG
webp <- function(path) {
  if (knitr::is_latex_output() || !grepl(pattern = "[.]png$", x = path)) return(path)
  tryCatch(
      {
        if (!requireNamespace("magick", quietly = TRUE) || !isTRUE(magick::magick_config()$webp)) return(path)
        img <- magick::image_read(path)
        out <- sub(pattern = "[.]png$", replacement = ".webp", x = path)
        tmp <- tempfile(fileext = ".webp")
        magick::image_write(
            image   = img
          , path    = out
          , format  = "webp"
          , quality = 100
          , defines = c(`webp:lossless` = "true")
          )
        magick::image_write(image = img, path = tmp, format = "webp", quality = 90)
        if (file.size(tmp) < file.size(out) / 2) file.copy(from = tmp, to = out, overwrite = TRUE)
        unlink(c(path, tmp))
        out
      }
    , error = function(e) path
    )
}
knitr::opts_chunk$set(collapse = FALSE, comment = "", fig.width = 7, fig.height = 5,
                      dev = if (knitr::is_latex_output()) "cairo_pdf" else "ragg_png", dpi = 96,
                      fig.process = webp)
if (requireNamespace("bayesammi", quietly = TRUE)) library(bayesammi) else
  pkgload::load_all(export_all = FALSE, helpers = FALSE, quiet = TRUE)


## ----fit----------------------------------------------------------------------
fit <-
  bayes_ammi(
      .data   = Maiz
    , .y      = y
    , .gen    = entry
    , .env    = site
    , .rep    = rep
    , .nIter  = 4000
    , .k      = 2
    , .burnin = 1000
    , .thin   = 1
    , .chains = 3
    , .priors = bam_priors(
          mu_mean     = 0
        , mu_var      = 1e+8
        , tau_mean    = 0
        , tau_var     = 1e+8
        , delta_mean  = 0
        , delta_var   = 1e+8
        , lambda_mean = 0
        , lambda_var  = 1e+8
        , nu          = 1
        , s2          = NULL
        )
    , .level  = 0.95
    , .seed   = 5
    , .align  = TRUE
    , .model  = "AMMI"
    , .het    = FALSE
    , .waic   = TRUE
    )


## ----trace, fig.height = 6----------------------------------------------------
plot(x = fit, which = "trace", pars = c("mu", "sigma2", "lambda[1]", "lambda[2]", "alpha[8,1]"))


## ----gelman-------------------------------------------------------------------
fit$diagnostics$gelman
summary(object = fit$summary$rhat)


## ----raftery------------------------------------------------------------------
fit$diagnostics$raftery[[1]]


## ----heidel-------------------------------------------------------------------
fit$diagnostics$heidel[[1]]


## ----ess----------------------------------------------------------------------
fit$summary[block %in% c("mu", "sigma2", "lambda"), list(parameter, mean, sd, mc_se, ess, rhat)]


## ----coda---------------------------------------------------------------------
coda::autocorr.diag(mcmc.obj = fit$draws[, c("lambda[1]", "lambda[2]")], lags = c(1, 5, 10))


## ----posterior----------------------------------------------------------------
dr <- posterior::as_draws_df(x = fit)
posterior::summarise_draws(
    .x = posterior::subset_draws(x = dr, variable = c("mu", "sigma2", "lambda[1]", "lambda[2]"))
  , "mean", "sd", "rhat", "ess_bulk", "ess_tail"
  )


## ----bayesplot, fig.height = 3.4----------------------------------------------
bayesplot::mcmc_rank_overlay(x = fit$draws, pars = c("lambda[1]", "lambda[2]"))


## ----int----------------------------------------------------------------------
bam_interpret(
    x             = fit
  , higher_better = TRUE
  , top           = 5
  )$findings[1:2]

