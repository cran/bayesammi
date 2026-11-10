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
knitr::opts_chunk$set(collapse = FALSE, comment = "", fig.width = 7, fig.height = 4.6,
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
    , .nIter  = 30000
    , .k      = 2
    , .burnin = 5000
    , .thin   = 2
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
    , .seed   = 2011
    , .align  = TRUE
    , .model  = "AMMI"
    , .het    = FALSE
    , .waic   = TRUE
    )
fit$settings$seconds


## ----table1-------------------------------------------------------------------
pub <- utils::read.csv(file = system.file("extdata", "crossa2011.csv", package = "bayesammi"))
s <- fit$summary
t1 <- merge(x = pub[pub$table == "1", ], y = s[, list(parameter, mean, sd, lower, upper)],
            by = "parameter", suffixes = c(".paper", ".here"), sort = FALSE)
t1$dev_sd <- round(x = (t1$mean.here - t1$mean.paper) / t1$sd.paper, digits = 3)
num <- c("mean.paper", "mean.here", "sd.paper", "sd.here", "lower.paper", "lower.here", "upper.paper", "upper.here")
t1[num] <- lapply(X = t1[num], FUN = signif, digits = 5)
t1[, c("parameter", "mean.paper", "mean.here", "sd.paper", "sd.here", "dev_sd")]
t1[, c("parameter", "lower.paper", "lower.here", "upper.paper", "upper.here")]


## ----s2-----------------------------------------------------------------------
fit_s2 <-
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
        , s2          = 1e8
        )
    , .level  = 0.95
    , .seed   = 1
    , .align  = TRUE
    , .model  = "AMMI"
    , .het    = FALSE
    , .waic   = TRUE
    )
c(default = s[parameter == "sigma2", mean], s2_1e8 = fit_s2$summary[parameter == "sigma2", mean])


## ----b1-----------------------------------------------------------------------
b1 <-
  merge(
      x        = pub[pub$table == "B1" & grepl(pattern = ",1\\]$", x = pub$parameter), ]
    , y        = s[, list(parameter, mean, sd, lower, upper)]
    , by       = "parameter"
    , suffixes = c(".paper", ".here")
    , sort     = FALSE
    )
b1$dev_sd <- round(x = (b1$mean.here - b1$mean.paper) / b1$sd.paper, digits = 3)
b1[c("mean.here", "sd.here")] <- lapply(X = b1[c("mean.here", "sd.here")], FUN = round, digits = 4)
b1[, c("parameter", "mean.paper", "mean.here", "sd.paper", "sd.here", "dev_sd")]


## ----b2-----------------------------------------------------------------------
b2 <-
  merge(
      x        = pub[pub$table == "B1" & grepl(pattern = ",2\\]$", x = pub$parameter), ]
    , y        = s[, list(parameter, mean, sd)]
    , by       = "parameter"
    , suffixes = c(".paper", ".here")
    , sort     = FALSE
    )
b2[c("mean.here", "sd.here")] <- lapply(X = b2[c("mean.here", "sd.here")], FUN = round, digits = 4)
b2[, c("parameter", "mean.paper", "mean.here", "sd.paper", "sd.here")]


## ----noalign------------------------------------------------------------------
fit_free <-
  bayes_ammi(
      .data   = Maiz
    , .y      = y
    , .gen    = entry
    , .env    = site
    , .rep    = rep
    , .nIter  = 30000
    , .k      = 2
    , .burnin = 5000
    , .thin   = 2
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
    , .seed   = 2011
    , .align  = FALSE
    , .model  = "AMMI"
    , .het    = FALSE
    , .waic   = TRUE
    )
d <- fit_free$draws[[1]][, "alpha[8,2]"]
c(share_positive = mean(d > 0), mean = mean(d), paper = pub$mean[pub$parameter == "alpha[8,2]"])


## ----share, include = FALSE---------------------------------------------------
am <- s[parameter == "alpha[8,2]", mean]
pp <- pub$mean[pub$parameter == "alpha[8,2]"]


## ----fig1, fig.height = 5-----------------------------------------------------
plot(x = fit, which = "trace", pars = c("lambda[1]", "lambda[2]"))
plot(x = fit, which = "density", pars = c("lambda[1]", "lambda[2]"))


## ----fig2, fig.height = 6-----------------------------------------------------
plot(x = fit, which = "biplot")


## ----fig34, fig.width = 9, fig.height = 5.5-----------------------------------
bp <-
  biplots(
      model         = fit
    , burnin        = 0
    , thin          = 0.25
    , pb            = c(0.05, 0.01)
    , plot_stable   = TRUE
    , plot_unstable = TRUE
    , ncolors       = 5
    , .which        = "both"
    , .grid         = 100
    )
bp$biplot
bp$biplot_data[axis1_significant == TRUE, list(byplot, name, x, axis1_lower, axis1_upper, group)]

