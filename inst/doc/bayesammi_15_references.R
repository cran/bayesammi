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


## ----cit, eval = FALSE--------------------------------------------------------
# citation(package = "bayesammi")


## ----cite-one-----------------------------------------------------------------
g <-
  ge_ammi(
      .data = Maiz
    , .y    = y
    , .gen  = entry
    , .env  = site
    , .rep  = rep
    )
bam_cite(
    x     = g
  , style = c("text", "bibtex")[2]
  )


## ----bib----------------------------------------------------------------------
basename(system.file("REFERENCES.bib", package = "bayesammi"))


## ----list, echo = FALSE, results = "asis"-------------------------------------
refs <- bayesammi:::.bam_refs_list()
groups <-
  list(
      `Bayesian AMMI and GGE models` = c("Crossa2011", "PerezElizalde2012", "Viele2000", "Hoff2009", "Yan2000")
    , `Bayesian stability models`    = c("FinlayWilkinson1963", "EberhartRussell1966", "Shukla1972", "LianCampos2016",
                                      "Oliveira2018", "EdwardsJannink2006", "Eskridge1990")
    , `Stability indices`            = c("Wricke1962", "LinBinns1988", "Annicchiarico1992", "Purchase2000", "Sneller1997",
                              "Zobel1994", "Zali2012", "Olivoto2019", "Farshadfar2008")
    , `Decisions`                    = c("Dias2022")
    , `Model choice and checking`    = c("Spiegelhalter2002", "Watanabe2010", "Vehtari2017", "Gelman1996", "Gabry2019")
    )
for (gname in names(groups)) {
  cat("## ", gname, "\n\n", sep = "")
  cat(paste0("* ", vapply(refs[groups[[gname]]], `[[`, "", "text")), sep = "\n")
  cat("\n")
}

