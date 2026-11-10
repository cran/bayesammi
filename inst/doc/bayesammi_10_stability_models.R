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


## ----priors-------------------------------------------------------------------
bam_stab_priors(
    mu_var = NULL
  , g_var  = NULL
  , h_var  = NULL
  , b_var  = 100
  , nu     = 1
  , s2     = NULL
  , nu_d   = 1
  , s2_d   = NULL
  )


## ----fit----------------------------------------------------------------------
fs <-
  bayes_stability(
      .data   = Maiz
    , .y      = y
    , .gen    = entry
    , .env    = site
    , .rep    = rep
    , .model  = "ER"
    , .nIter  = 3000
    , .burnin = 500
    , .thin   = 1
    , .chains = 2
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
    , .seed   = 1
    , .waic   = TRUE
    )
fs


## ----facts, include = FALSE---------------------------------------------------
s <- fs$stability
up <- s[slope_lower > 1, Gen]
dn <- s[slope_upper < 1, Gen]
best <- s[which.max(prob_best_safety)]
lst <- function(x) if (length(x) == 0) "none" else if (length(x) == 1) x else paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])


## ----slopes, fig.height = 5.2-------------------------------------------------
plot(x = fs, which = "slopes")


## ----lines, fig.height = 5.2--------------------------------------------------
plot(x = fs, which = "lines")


## ----deviations, fig.height = 4.4---------------------------------------------
plot(x = fs, which = "deviations")


## ----ls-----------------------------------------------------------------------
data.table::data.table(
    Gen        = s$Gen
  , slope_ls   = round(x = fs$ls$slope, digits = 3)
  , slope      = round(x = s$slope, digits = 3)
  , dev_var_ls = round(fs$ls$dev_var)
  , dev_var    = round(s$dev_var)
  )


## ----lsr, include = FALSE-----------------------------------------------------
rs <- round(x = cor(x = fs$ls$slope, y = s$slope), digits = 3)
rd <- round(x = cor(x = fs$ls$dev_var, y = s$dev_var, method = "spearman"), digits = 2)


## ----sf-----------------------------------------------------------------------
s[order(-safety_first), list(
    Gen
  , mean
  , safety_first
  , safety_first_lower
  , safety_first_upper
  , prob_best_safety
  )]


## ----models-------------------------------------------------------------------
ff <-
  bayes_stability(
      .data   = Maiz
    , .y      = y
    , .gen    = entry
    , .env    = site
    , .rep    = rep
    , .model  = "FW"
    , .nIter  = 3000
    , .burnin = 500
    , .thin   = 1
    , .chains = 2
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
    , .seed   = 1
    , .waic   = TRUE
    )
fh <-
  bayes_stability(
      .data   = Maiz
    , .y      = y
    , .gen    = entry
    , .env    = site
    , .rep    = rep
    , .model  = "Shukla"
    , .nIter  = 3000
    , .burnin = 500
    , .thin   = 1
    , .chains = 2
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
    , .seed   = 1
    , .waic   = TRUE
    )
cp <-
  bam_compare(
      ER     = fs
    , FW     = ff
    , Shukla = fh
    )
cp


## ----cpf, include = FALSE-----------------------------------------------------
t <- cp$table


## ----decide-------------------------------------------------------------------
bam_superiority(
    fit           = fs
  , .s            = 0.2
  , higher_better = TRUE
  )$marginal


## ----int----------------------------------------------------------------------
bam_interpret(
    x             = fs
  , higher_better = TRUE
  , top           = 5
  )


## ----cite---------------------------------------------------------------------
bam_cite(
    x     = fs
  , style = c("text", "bibtex")[1]
  )

