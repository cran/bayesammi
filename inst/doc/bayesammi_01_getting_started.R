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
knitr::opts_chunk$set(collapse = FALSE, comment = "", fig.width = 7, fig.height = 4.2,
                      dev = if (knitr::is_latex_output()) "cairo_pdf" else "ragg_png", dpi = 96,
                      fig.process = webp)
if (requireNamespace("bayesammi", quietly = TRUE)) library(bayesammi) else
  pkgload::load_all(export_all = FALSE, helpers = FALSE, quiet = TRUE)


## ----data---------------------------------------------------------------------
data(Maiz)
str(Maiz)
m <-
  ge_mean(
      .data = Maiz
    , .y    = y
    , .gen  = entry
    , .env  = site
    )
round(m$ge_means[, 1:6])
m$grand_mean


## ----ls-----------------------------------------------------------------------
ls_fit <-
  ge_ammi(
      .data = Maiz
    , .y    = y
    , .gen  = entry
    , .env  = site
    , .rep  = rep
    )
round(ls_fit$lambdas)
ls_fit$sigma2


## ----fit----------------------------------------------------------------------
fit <-
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
    , .seed   = 1
    , .align  = TRUE
    , .model  = "AMMI"
    , .het    = FALSE
    , .waic   = TRUE
    )
fit


## ----scores-------------------------------------------------------------------
head(x = fit$summary[block == "alpha"], n = 9)


## ----trace, fig.height = 5.5--------------------------------------------------
plot(x = fit, which = "trace")


## ----biplot, fig.height = 6---------------------------------------------------
plot(x = fit, which = "biplot")


## ----hpd, fig.width = 9, fig.height = 5.5-------------------------------------
bp <-
  biplots(
      model         = fit
    , burnin        = 0
    , thin          = 1
    , pb            = c(0.05, 0.01)
    , plot_stable   = TRUE
    , plot_unstable = TRUE
    , ncolors       = 5
    , .which        = "both"
    , .grid         = 100
    )
bp$biplot
bp$biplot_data[byplot == "Genotype"]


## ----int----------------------------------------------------------------------
bam_interpret(
    x             = fit
  , higher_better = TRUE
  , top           = 5
  )

