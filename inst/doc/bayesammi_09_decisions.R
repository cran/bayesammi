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
knitr::opts_chunk$set(collapse = FALSE, comment = "", fig.width = 7, fig.height = 4.8,
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
    , .chains = 2
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
    , .seed   = 9
    , .align  = TRUE
    , .model  = "AMMI"
    , .het    = FALSE
    , .waic   = TRUE
    )


## ----sup----------------------------------------------------------------------
sp <-
  bam_superiority(
      fit           = fit
    , .s            = 0.2
    , higher_better = TRUE
    )
sp


## ----sup-plots, fig.height = 4.4----------------------------------------------
plot(x = sp, which = "marginal")
plot(x = sp, which = "pairwise")


## ----sup-cond, fig.height = 4.8-----------------------------------------------
plot(x = sp, which = "conditional")


## ----ranks--------------------------------------------------------------------
rk <-
  bam_ranks(
      fit           = fit
    , higher_better = TRUE
    , .env          = NULL
    )
rk


## ----ranks-plot, fig.height = 4.6---------------------------------------------
plot(x = rk)


## ----ranks-env----------------------------------------------------------------
bam_ranks(
    fit           = fit
  , higher_better = TRUE
  , .env          = "12"
  )$summary[1:3]


## ----stab---------------------------------------------------------------------
st <-
  bam_stability(
      fit           = fit
    , higher_better = TRUE
    )
st


## ----stab-plot, fig.height = 5------------------------------------------------
plot(x = st)


## ----pred---------------------------------------------------------------------
pr <-
  bam_predict(
      fit    = fit
    , .level = 0.95
    )
pr$genotypes


## ----pred-plot, fig.height = 4.4----------------------------------------------
plot(x = pr, which = "genotypes")


## ----win----------------------------------------------------------------------
bam_winners(
    fit           = fit
  , higher_better = TRUE
  )$mega


## ----compare------------------------------------------------------------------
fits <-
  lapply(
      X   = 1:3
    , FUN = function(K) {
        bayes_ammi(
            .data   = Maiz
          , .y      = y
          , .gen    = entry
          , .env    = site
          , .rep    = rep
          , .nIter  = 2000
          , .k      = K
          , .burnin = 500
          , .thin   = 1
          , .chains = 1
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
          , .seed   = 9
          , .align  = TRUE
          , .model  = "AMMI"
          , .het    = TRUE
          , .waic   = TRUE
          )
      }
    )
names(fits) <- paste0("K", 1:3)
cp <- do.call(what = bam_compare, args = fits)
cp


## ----compare-plot, fig.height = 3.6-------------------------------------------
plot(x = cp)


## ----kk, include = FALSE------------------------------------------------------
t <- cp$table


## ----ppc----------------------------------------------------------------------
pp <-
  bam_ppc(
      fit    = fit
    , .draws = 500
    , .seed  = 2
    )
pp


## ----ppc-plot, fig.height = 5.5-----------------------------------------------
plot(x = pp)


## ----bayesplot, fig.height = 3.6----------------------------------------------
dr <- posterior::as_draws_df(x = fit)
posterior::summarise_draws(.x = posterior::subset_draws(x = dr, variable = c("lambda[1]", "lambda[2]", "sigma2")))
bayesplot::mcmc_areas(x = fit$draws, pars = c("tau[4]", "tau[5]", "tau[6]"), prob = 0.95)


## ----int----------------------------------------------------------------------
bam_interpret(
    x             = sp
  , higher_better = TRUE
  , top           = 5
  )
bam_interpret(
    x             = st
  , higher_better = TRUE
  , top           = 5
  )

