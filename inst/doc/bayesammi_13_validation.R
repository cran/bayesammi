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


## ----tables-------------------------------------------------------------------
ge_eff(
    .data = Maiz
  , .y    = y
  , .gen  = entry
  , .env  = site
  )$ge_svd$d[1:3]


## ----vmf----------------------------------------------------------------------
set.seed(seed = 7)
p <- 18
vapply(
    X         = c(1, 10, 100)
  , FUN       = function(kap) {
      X <- t(replicate(n = 3000, expr = bayesammi:::.bam_rvmf(v = c(kap, rep(0, p - 1)))))
      c(
          kappa     = kap
        , simulated = mean(X[, 1])
        , theory    = besselI(x = kap, nu = p / 2, expon.scaled = TRUE) / besselI(
              x            = kap
            , nu           = p / 2 - 1
            , expon.scaled = TRUE
            )
        )
    }
  , FUN.VALUE = numeric(length = 3)
  )


## ----same---------------------------------------------------------------------
ns <- asNamespace(ns = "bayesammi")
cd <-
  ns$.bam_prep(
      "check"
    , Maiz
    , quote(y)
    , quote(entry)
    , quote(site)
    , quote(rep)
    , environment()
    , allow_missing = TRUE
    )
ls <- ns$.bam_ls(Y = cd$Y, K = 2)
pr <-
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
pr$s2 <- (cd$ssw + 4 * ls$rss) / 600
st <- ns$.bam_start(ls, pr$s2, 1)
set.seed(seed = 5)
t_r <- system.time(r <- ns$.bam_chain(cd$Y, 4, cd$ssw, 2, pr, st, list(alpha = ls$alpha), 500, 1:500))
set.seed(seed = 5)
t_c <- system.time(cc <- ns$.bam_run_chain(cd, 2, FALSE, FALSE, pr, st, ls$alpha, 500, rep(1L, 500), TRUE, FALSE))
c(scalars = max(abs(cc$scalar - r$scalar)), alpha = max(abs(cc$alpha - r$alpha)), gamma = max(abs(cc$gamma - r$gamma)))
c(R = t_r[["elapsed"]], compiled = t_c[["elapsed"]])


## ----same-stab----------------------------------------------------------------
vapply(
    X         = 0:2
  , FUN       = function(code) {
      ls <- ns$.bam_stab_ls(cd, code)
      pr <-
        ns$.bam_stab_fill_priors(
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
          , cd
          , ls
          )
      st <- ns$.bam_stab_start(ls, 1L, code)
      set.seed(seed = 9)
      r  <- ns$.bam_stab_chain(cd$S, cd$Q, cd$N * 1, code, pr, st, 200, rep(1L, 200))
      set.seed(seed = 9)
      cc <-
        ns$.bam_stab_cpp(
            cd$S
          , cd$Q
          , cd$N * 1
          , code
          , pr
          , st
          , 200L
          , rep(1L, 200)
          , cd$data$y
          , cd$data$cell - 1L
          , FALSE
          )
      max(abs(cc$mu - r$mu), abs(cc$g - r$g), abs(cc$b - r$b), abs(cc$d - r$d))
    }
  , FUN.VALUE = numeric(length = 1)
  )


## ----sim----------------------------------------------------------------------
set.seed(seed = 8)
G <- 8
E <- 12
B <- 3
A  <- qr.Q(qr(cbind(1, matrix(data = rnorm(G * 2), nrow = G))))[, 2:3]
Gm <- qr.Q(qr(cbind(1, matrix(data = rnorm(E * 2), nrow = E))))[, 2:3]
tau <- as.vector(scale(x = rnorm(n = G, mean = 0, sd = 6), scale = FALSE))
del <- as.vector(scale(x = rnorm(n = E, mean = 0, sd = 10), scale = FALSE))
M <- 100 + outer(X = tau, Y = del, FUN = "+") + A %*% (c(80, 40) * t(Gm))
sim <- expand.grid(rep = 1:B, gen = 1:G, env = 1:E)
sim$y <- M[cbind(sim$gen, sim$env)] + rnorm(n = nrow(sim), mean = 0, sd = 5)
fs <-
  bayes_ammi(
      .data   = sim
    , .y      = y
    , .gen    = gen
    , .env    = env
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
    , .seed   = 9
    , .align  = TRUE
    , .model  = "AMMI"
    , .het    = FALSE
    , .waic   = TRUE
    )
s <- fs$summary
c(tau_covered   = mean(s[block == "tau", lower] <= tau & tau <= s[block == "tau", upper]),
  delta_covered = mean(s[block == "delta", lower] <= del & del <= s[block == "delta", upper]))
s[block %in% c("sigma2", "lambda"), list(parameter, truth = c(25, 80, 40), mean, lower, upper)]


## ----ref----------------------------------------------------------------------
bam_reference(
    source   = "Crossa"
  , evaluate = TRUE
  , level    = 1
  )[, list(quantity, value, computed, agrees)]

