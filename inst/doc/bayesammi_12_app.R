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


## ----launch, eval = FALSE-----------------------------------------------------
# bam_app(
#     launch.browser = TRUE
#   )


## ----webp---------------------------------------------------------------------
isTRUE(tryCatch(requireNamespace("magick", quietly = TRUE) && isTRUE(magick::magick_config()$webp),
                error = function(e) FALSE))


## ----shot-model, echo = FALSE, out.width = "100%"-----------------------------
knitr::include_graphics(path = if (knitr::is_latex_output()) "figures/app-model.jpg" else "figures/app-model.webp")


## ----shot-stab, echo = FALSE, out.width = "100%"------------------------------
knitr::include_graphics(path = if (knitr::is_latex_output()) "figures/app-stability.jpg" else "figures/app-stability.webp")


## ----code---------------------------------------------------------------------
settings <-
  list(
      .y      = "y"
    , .gen    = "entry"
    , .env    = "site"
    , .rep    = "rep"
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
    , .model  = "AMMI"
    , .het    = TRUE
    )
cat(bayesammi:::.bam_app_code(a = settings, source = "maiz"))


## ----deploy, eval = FALSE-----------------------------------------------------
# # app.R on a Shiny server
# bayesammi::bam_app_object()

