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
knitr::opts_chunk$set(collapse = FALSE, comment = "", fig.width = 7, fig.height = 6,
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
    , .seed   = 11
    , .align  = TRUE
    , .model  = "AMMI"
    , .het    = FALSE
    , .waic   = TRUE
    )


## ----ind----------------------------------------------------------------------
ind <-
  bam_indices(
      fit           = fit
    , higher_better = TRUE
    , .prob         = 0.25
    , .alpha        = 0.05
    )
ind


## ----ind-plot, fig.height = 6-------------------------------------------------
plot(x = ind)


## ----agree, include = FALSE---------------------------------------------------
m <- ind$means
rc <- round(x = cor(x = m$wricke, y = m$lin_binns, method = "spearman"), digits = 2)
neg <- ind$indices[index == "shukla" & lower < 0, .N]


## ----ai-----------------------------------------------------------------------
ai <-
  bam_ammi_indices(
      fit           = fit
    , higher_better = TRUE
    , .weight_y     = 50
    )
ai


## ----ai-plot, fig.height = 6--------------------------------------------------
plot(x = ai)


## ----waasy--------------------------------------------------------------------
data.table::rbindlist(
    l = lapply(
        X   = c(30, 50, 70)
      , FUN = function(wy) {
          a <-
            bam_ammi_indices(
                fit           = fit
              , higher_better = TRUE
              , .weight_y     = wy
              )
          a$indices[index == "WAASY"][order(-mean)][1:3, list(weight = wy, Gen, mean, prob_best)]
        }
      )
  )


## ----from-stab----------------------------------------------------------------
fs <-
  bayes_stability(
      .data   = Maiz
    , .y      = y
    , .gen    = entry
    , .env    = site
    , .rep    = rep
    , .model  = "ER"
    , .nIter  = 2000
    , .burnin = 500
    , .thin   = 1
    , .chains = 1
    , .priors = bam_stab_priors(
          mu_var = NULL
        , g_var  = NULL
        , h_var  = NULL
        , b_var  = 100
        , nu     = 1
        , s2     = NULL
        , nu_d   = 1
        , s2_d   = NULL
        )
    , .level  = 0.95
    , .alpha  = 0.05
    , .seed   = 11
    , .waic   = TRUE
    )
bam_indices(
    fit           = fs
  , higher_better = TRUE
  , .prob         = 0.25
  , .alpha        = 0.05
  )$means


## ----int----------------------------------------------------------------------
bam_interpret(
    x             = ind
  , higher_better = TRUE
  , top           = 5
  )

