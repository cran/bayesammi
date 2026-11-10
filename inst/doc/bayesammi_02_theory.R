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


## ----priors-------------------------------------------------------------------
bam_priors(
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


## ----contrasts----------------------------------------------------------------
K <-
  matrix_k(
      n = 4
    )
round(x = K %*% t(K), digits = 12)
round(x = K %*% rep(1, 4), digits = 12)


## ----vmf----------------------------------------------------------------------
set.seed(seed = 1)
X <- t(replicate(n = 2000, expr = bayesammi:::.bam_rvmf(v = c(10, 0, 0))))
c(mean_first = mean(X[, 1]),
  theory     = besselI(x = 10, nu = 1.5) / besselI(x = 10, nu = 0.5))

