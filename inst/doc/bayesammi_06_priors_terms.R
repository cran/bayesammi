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
knitr::opts_chunk$set(collapse = FALSE, comment = "", fig.width = 7, fig.height = 4.4,
                      dev = if (knitr::is_latex_output()) "cairo_pdf" else "ragg_png", dpi = 96,
                      fig.process = webp)
if (requireNamespace("bayesammi", quietly = TRUE)) library(bayesammi) else
  pkgload::load_all(export_all = FALSE, helpers = FALSE, quiet = TRUE)


## ----vague--------------------------------------------------------------------
f0 <-
  bayes_ammi(
      .data   = Maiz
    , .y      = y
    , .gen    = entry
    , .env    = site
    , .rep    = rep
    , .nIter  = 3000
    , .k      = 2
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
    , .seed   = 6
    , .align  = TRUE
    , .model  = "AMMI"
    , .het    = FALSE
    , .waic   = TRUE
    )
f0$summary[block %in% c("lambda", "sigma2"), list(parameter, mean, sd, lower, upper)]


## ----informative--------------------------------------------------------------
f1 <-
  bayes_ammi(
      .data   = Maiz
    , .y      = y
    , .gen    = entry
    , .env    = site
    , .rep    = rep
    , .nIter  = 3000
    , .k      = 2
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
        , lambda_mean = 5000
        , lambda_var  = 250^2
        , nu          = 50
        , s2          = 8e5
        )
    , .level  = 0.95
    , .seed   = 6
    , .align  = TRUE
    , .model  = "AMMI"
    , .het    = FALSE
    , .waic   = TRUE
    )
pm <- rbind(vague       = f0$summary[block %in% c("lambda", "sigma2"), mean],
            informative = f1$summary[block %in% c("lambda", "sigma2"), mean])
colnames(pm) <- f0$summary[block %in% c("lambda", "sigma2"), parameter]
round(pm)


## ----k------------------------------------------------------------------------
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
          , .nIter  = 3000
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
          , .seed   = 6
          , .align  = TRUE
          , .model  = "AMMI"
          , .het    = FALSE
          , .waic   = TRUE
          )
      }
    )
data.table::rbindlist(l = lapply(X = fits, FUN = function(f) {
  f$summary[block == "lambda", list(K = f$settings$k, parameter, mean, sd, lower, upper)]
}))


## ----klow, include = FALSE----------------------------------------------------
lw <- vapply(X = fits, FUN = function(f) { l <- f$summary[block == "lambda"]
l$lower[nrow(l)] }, FUN.VALUE = numeric(length = 1))


## ----s2-----------------------------------------------------------------------
vapply(
    X         = fits
  , FUN       = function(f) f$summary[parameter == "sigma2", mean]
  , FUN.VALUE = numeric(length = 1)
  )


## ----waic---------------------------------------------------------------------
names(fits) <- paste0("K", 1:3)
cp <- do.call(what = bam_compare, args = fits)
cp


## ----waicw, include = FALSE---------------------------------------------------
best <- cp$table$model[1]

