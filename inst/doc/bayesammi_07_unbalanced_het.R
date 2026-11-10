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


## ----unbal--------------------------------------------------------------------
set.seed(seed = 21)
mz   <- Maiz[-sample(x = nrow(Maiz), size = 60), ]
gone <- (mz$entry == 2 & mz$site == 5) | (mz$entry == 7 & mz$site == 11) | (mz$entry == 9 & mz$site == 18)
mz   <- mz[!gone, ]
table(plots = table(factor(mz$entry):factor(mz$site)))


## ----unbal-fit----------------------------------------------------------------
fu <-
  bayes_ammi(
      .data   = mz
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
    , .seed   = 3
    , .align  = TRUE
    , .model  = "AMMI"
    , .het    = FALSE
    , .waic   = TRUE
    )
fu


## ----unbal-pred---------------------------------------------------------------
pu <-
  bam_predict(
      fit    = fu
    , .level = 0.95
    )
pu$cells[n == 0]
Maiz_cells <- data.table::as.data.table(x = Maiz)[, list(full = mean(y)), by = list(
    Gen = factor(entry)
  , Env = factor(site)
  )]
merge(x = pu$cells[n == 0], y = Maiz_cells, by = c("Gen", "Env"))[, list(
    Gen
  , Env
  , mean
  , lower
  , upper
  , full
  )]


## ----cover, include = FALSE---------------------------------------------------
mm <- merge(x = pu$cells[n == 0], y = Maiz_cells, by = c("Gen", "Env"))
inside <- sum(mm$lower <= mm$full & mm$full <= mm$upper)
miss1  <- mm[!(lower <= full & full <= upper)]
sdr    <- mean(pu$cells[n == 0, sd]) / mean(pu$cells[n == 4, sd])


## ----unbal-plot, fig.height = 5-----------------------------------------------
plot(x = pu, which = "cells")


## ----het----------------------------------------------------------------------
fh <-
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
    , .seed   = 4
    , .align  = TRUE
    , .model  = "AMMI"
    , .het    = TRUE
    , .waic   = TRUE
    )
fh$summary[block == "sigma2", list(parameter, mean, lower, upper, rhat)]


## ----het-range, include = FALSE-----------------------------------------------
v <- fh$summary[block == "sigma2", mean]
rr <- max(v) / min(v)


## ----het-plot, fig.height = 4.6-----------------------------------------------
ggplot2::ggplot(
    data    = fh$summary[block == "sigma2"]
  , mapping = ggplot2::aes(x = reorder(label, mean), y = mean, ymin = lower, ymax = upper)
  ) +
  ggplot2::geom_pointrange(colour = "#2c7fb8") +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 6), labels = scales::label_comma()) +
  ggplot2::labs(x = "Environment", y = "Error variance (posterior mean, 95% HPD)") +
  ggplot2::theme_bw()


## ----het-common---------------------------------------------------------------
fc <-
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
    , .seed   = 4
    , .align  = TRUE
    , .model  = "AMMI"
    , .het    = FALSE
    , .waic   = TRUE
    )
cp <-
  bam_compare(
      common          = fc
    , environment_var = fh
    )
cp


## ----dw, include = FALSE------------------------------------------------------
dw  <- cp$table$dWAIC[2]
se  <- cp$table$dWAIC_se[2]
l2  <- cp$table$lambda_last_lower


## ----het-ppc------------------------------------------------------------------
pp <-
  bam_ppc(
      fit    = fc
    , .draws = 500
    , .seed  = 1
    )
pp


## ----ppcp, include = FALSE----------------------------------------------------
pc <- pp$checks
pv <- pc[statistic == "largest / smallest environment variance", p_value]

