#' Biplots with highest posterior density regions
#'
#' @description
#' \code{biplots()} draws the Bayesian AMMI biplot of the first two bilinear terms with a highest posterior density
#' (HPD) region around every genotype and environment (Crossa et al. 2011, Fig. 3 and 4). A region that contains the
#' origin marks a genotype or environment whose interaction is not different from zero (stable); regions that
#' overlap mark genotypes or environments that cannot be told apart, a homogeneous group.
#'
#' @details
#' The points are the biplot coordinates of the paper,
#' \deqn{\left(\lambda_1^{1/2}\alpha_{i1}, \lambda_2^{1/2}\alpha_{i2}\right), \qquad
#'   \left(\lambda_1^{1/2}\gamma_{j1}, \lambda_2^{1/2}\gamma_{j2}\right),}{(lambda_1^(1/2) alpha_i1, lambda_2^(1/2)
#'   alpha_i2), (lambda_1^(1/2) gamma_j1, lambda_2^(1/2) gamma_j2),}
#' computed draw by draw. For each genotype and environment a bivariate kernel density \eqn{\hat{f}}{f} is estimated
#' from the draws (normal kernels, normal reference bandwidths); the \eqn{100(1 - p)}{100(1 - p)}\% HPD region is
#' \eqn{\{\mathbf{x} : \hat{f}(\mathbf{x}) \ge c_p\}}{{x : f(x) >= c_p}} with \eqn{c_p}{c_p} the
#' \eqn{p}{p}-quantile of \eqn{\hat{f}}{f} at the draws, so that the region holds \eqn{100(1 - p)}{100(1 - p)}\% of
#' them. The region contains the origin when \eqn{\hat{f}(\mathbf{0}) \ge c_p}{f(0) >= c_p} (\code{stable}). The
#' criterion of the paper, the marginal HPD interval of the first coordinate excluding zero, is given as
#' \code{axis1_significant}.
#'
#' \strong{Changes from the earlier versions.} The regions were drawn for the singular vectors without the factor
#' \eqn{\lambda_k^{1/2}}{lambda_k^(1/2)}, with a fixed bandwidth of 0.2 whatever the scale; a genotype counted as
#' stable when the bounding box of its contour covered the origin.
#'
#' @param model the result of \code{\link{bayes_ammi}()} (with \code{.k >= 2}).
#' @param burnin proportion of the stored draws of every chain discarded at its start (default 0.3).
#' @param thin proportion of the remaining draws used (default 0.2: every fifth draw).
#' @param pb one or more probabilities \eqn{p}{p}: the regions hold \eqn{100(1 - p)}{100(1 - p)}\% of the posterior
#'   (default 0.05, the 95\% region; \code{c(0.05, 0.01)} draws the 95\% and 99\% regions of the paper). The
#'   classification uses the first value.
#' @param plot_stable,plot_unstable draw the regions that contain the origin, and those that do not (default both).
#' @param ncolors number of colours used for the regions (default 5).
#' @param .which \code{"both"} (default), \code{"genotype"} or \code{"environment"}.
#' @param .grid number of grid points per axis of the density estimate (default 100).
#' @return An object of class \code{c("bam_biplots", "bam_result", "list")}, with the components of the earlier
#'   version:
#' \describe{
#'   \item{\code{contorns}}{\code{data.table} of the region boundaries: \code{byplot} (Genotype, Environment),
#'     \code{ID} (number), \code{name}, \code{prob} (the region, for instance "95\%"), \code{piece} (a region can have
#'     several pieces), \code{x}, \code{y}, \code{stable}.}
#'   \item{\code{biplot_data}}{\code{data.table}, one row per genotype and environment: \code{byplot}, \code{ID},
#'     \code{name}, \code{x}, \code{y} (posterior means of the coordinates), \code{stable} (region contains the
#'     origin), \code{axis1_lower}, \code{axis1_upper} (marginal HPD interval of the first coordinate),
#'     \code{axis1_significant} (it excludes zero), \code{group} (sign of the first coordinate when significant:
#'     positive, negative, or not significant).}
#'   \item{\code{biplot}}{the \code{ggplot}.}
#' }
#'   Methods: \code{print()}, \code{plot()}, \code{\link{bam_interpret}()}.
#' @references
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#'
#' Hyndman, R. J. (1996). Computing and graphing highest density regions. \emph{The American Statistician},
#' \strong{50}(2), 120--126.
#' @author Julian Garcia Abadillo Velasco (\email{garciaabadillo.j@@ufl.edu}), Diego Jarquin
#'   (\email{diego.jarquin@@gmail.com}), Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bayes_ammi}}, \code{\link{bam_interpret}}
#' @examples
#' data(Maiz)
#' fm1 <-
#'   bayes_ammi(
#'       .data   = Maiz
#'     , .y      = y
#'     , .gen    = entry
#'     , .env    = site
#'     , .rep    = rep
#'     , .nIter  = 1000
#'     , .k      = 2
#'     , .burnin = 200
#'     , .thin   = 1
#'     , .chains = 1
#'     , .priors = bam_priors(
#'           mu_mean     = 0
#'         , mu_var      = 1e+8
#'         , tau_mean    = 0
#'         , tau_var     = 1e+8
#'         , delta_mean  = 0
#'         , delta_var   = 1e+8
#'         , lambda_mean = 0
#'         , lambda_var  = 1e+8
#'         , nu          = 1
#'         , s2          = NULL
#'         )
#'     , .level  = 0.95
#'     , .seed   = 1
#'     , .align  = TRUE
#'     , .model  = c("AMMI", "GGE")
#'     , .het    = FALSE
#'     , .waic   = TRUE
#'     )
#' bp <-
#'   biplots(
#'       model         = fm1
#'     , burnin        = 0
#'     , thin          = 1
#'     , pb            = c(0.05, 0.01)
#'     , plot_stable   = TRUE
#'     , plot_unstable = TRUE
#'     , ncolors       = 5
#'     , .which        = "both"
#'     , .grid         = 100
#'     )
#' bp$biplot
#' bp$biplot_data
#' @export
biplots <- function(model, burnin = 0.3, thin = 0.2, pb = 0.05, plot_stable = TRUE, plot_unstable = TRUE,
                    ncolors = 5, .which = c("both", "genotype", "environment"), .grid = 100L) {
  UseMethod("biplots")
}

#' @rdname biplots
#' @export
biplots.default <- function(model, burnin = 0.3, thin = 0.2, pb = 0.05, plot_stable = TRUE, plot_unstable = TRUE,
                            ncolors = 5, .which = c("both", "genotype", "environment"), .grid = 100L) {
  fn <- "biplots"
  .bam_check_class(model, "bam_bayes_ammi", "model", fn)
  if (model$settings$k < 2L) .bam_stop(fn, "the biplot needs at least two bilinear terms (bayes_ammi(.k >= 2)).")
  burnin  <- .bam_check_number(burnin, "burnin", fn, min = 0, max = 0.9)
  thin    <- .bam_check_number(thin, "thin", fn, min = 1e-6, max = 1)
  pb      <- sort(.bam_check_number(pb, "pb", fn, min = 1e-4, max = 0.5, scalar = FALSE), decreasing = TRUE)
  plot_stable   <- .bam_check_flag(plot_stable, "plot_stable", fn)
  plot_unstable <- .bam_check_flag(plot_unstable, "plot_unstable", fn)
  ncolors <- .bam_check_number(ncolors, "ncolors", fn, min = 1, integer = TRUE)
  .which  <- .bam_check_choice(.which[1L], c("both", "genotype", "environment"), ".which", fn)
  .grid   <- .bam_check_number(.grid, ".grid", fn, min = 25, integer = TRUE)
  m   <- attr(model, "meta")
  Gl  <- m$gen_labels
  El  <- m$env_labels
  K   <- model$settings$k
  # the draws used: burn-in and thinning within every chain
  sel <- lapply(model$draws, function(ch) {
    n <- nrow(ch)
    i <- seq.int(floor(burnin * n) + 1L, n, by = max(1L, round(1 / thin)))
    unclass(ch)[i, , drop = FALSE]
  })
  D <- do.call(rbind, sel)
  if (nrow(D) < 50L) .bam_stop(fn, "only ", nrow(D), " draws are left after burn-in and thinning; at least 50 are ",
                               "needed (run bayes_ammi() longer or use larger `thin`).")
  l1 <- sqrt(D[, "lambda[1]"])
  l2 <- sqrt(D[, "lambda[2]"])
  coords <- function(prefix, labs, type) {
    lapply(seq_along(labs), function(i) {
      list(type = type, ID = i, name = labs[i],
           x = D[, sprintf("%s[%s,1]", prefix, labs[i])] * l1, y = D[, sprintf("%s[%s,2]", prefix, labs[i])] * l2)
    })
  }
  sets <- c(if (.which != "environment") coords("alpha", Gl, "Genotype"),
            if (.which != "genotype") coords("gamma", El, "Environment"))
  res <- lapply(sets, .bam_hpd_region, pb = pb, grid = .grid, level1 = 1 - pb[1L])
  cont <- data.table::rbindlist(lapply(res, `[[`, "contour"))
  bd   <- data.table::rbindlist(lapply(res, `[[`, "summary"))
  cont[, stable := bd$stable[match(paste(cont$byplot, cont$ID), paste(bd$byplot, bd$ID))]]
  out <- list(contorns = cont, biplot_data = bd)
  out$biplot <- .bam_plot_hpd(cont, bd, pb, plot_stable, plot_unstable, ncolors)
  .bam_result(out, "biplots",
              label = sprintf("HPD regions of the AMMI biplot (%s), %d draws", paste0(100 * (1 - pb), "%",
                                                                                      collapse = ", "), nrow(D)),
              call = match.call(), meta = list(pb = pb, draws = nrow(D), level = 1 - pb[1L], k = K))
}

# the HPD regions of one point: contours at every level, the summary row
.bam_hpd_region <- function(s, pb, grid, level1) {
  x <- s$x
  y <- s$y
  bw <- function(v) {
    h <- MASS::bandwidth.nrd(v)
    if (!is.finite(h) || h <= 0) h <- 4 * 1.06 * max(stats::sd(v), 1e-8) * length(v)^(-1 / 5)
    h
  }
  h  <- c(bw(x), bw(y))
  lx <- range(x, 0) + c(-1, 1) * 1.5 * h[1L]
  ly <- range(y, 0) + c(-1, 1) * 1.5 * h[2L]
  kd <- MASS::kde2d(x, y, h = h, n = grid, lims = c(lx, ly))
  fx <- .bam_interp2(kd, x, y)
  f0 <- .bam_interp2(kd, 0, 0)
  cp <- stats::quantile(fx, probs = pb, names = FALSE, type = 1L)
  cl <- lapply(seq_along(pb), function(j) {
    lines <- grDevices::contourLines(kd$x, kd$y, kd$z, levels = cp[j])
    if (!length(lines)) return(NULL)
    data.table::rbindlist(lapply(seq_along(lines), function(p) {
      data.table::data.table(prob = sprintf("%g%%", 100 * (1 - pb[j])), piece = p, x = lines[[p]]$x,
                             y = lines[[p]]$y)
    }))
  })
  cont <- data.table::rbindlist(cl)
  if (nrow(cont)) cont <- cbind(data.table::data.table(byplot = s$type, ID = s$ID, name = s$name), cont)
  iv <- .bam_hpd(matrix(x), level1)
  sig <- iv[1L] > 0 || iv[2L] < 0
  list(contour = cont,
       summary = data.table::data.table(byplot = s$type, ID = s$ID, name = s$name, x = mean(x), y = mean(y),
                                        stable = f0 >= cp[1L], axis1_lower = iv[1L], axis1_upper = iv[2L],
                                        axis1_significant = sig,
                                        group = if (!sig) "not significant" else if (iv[1L] > 0) "positive" else
                                          "negative"))
}

# bilinear interpolation of a kde2d grid at points (0 outside the grid)
.bam_interp2 <- function(kd, x, y) {
  nx <- length(kd$x)
  ny <- length(kd$y)
  out <- numeric(length(x))
  ok  <- x >= kd$x[1L] & x <= kd$x[nx] & y >= kd$y[1L] & y <= kd$y[ny]
  if (!any(ok)) return(out)
  i <- pmin(findInterval(x[ok], kd$x), nx - 1L)
  j <- pmin(findInterval(y[ok], kd$y), ny - 1L)
  tx <- (x[ok] - kd$x[i]) / (kd$x[i + 1L] - kd$x[i])
  ty <- (y[ok] - kd$y[j]) / (kd$y[j + 1L] - kd$y[j])
  z  <- kd$z
  out[ok] <- (1 - tx) * (1 - ty) * z[cbind(i, j)] + tx * (1 - ty) * z[cbind(i + 1L, j)] +
    (1 - tx) * ty * z[cbind(i, j + 1L)] + tx * ty * z[cbind(i + 1L, j + 1L)]
  out
}

.bam_plot_hpd <- function(cont, bd, pb, plot_stable, plot_unstable, ncolors) {
  keep <- c(if (plot_stable) TRUE, if (plot_unstable) FALSE)
  cd <- cont[stable %in% keep]
  pd <- bd[stable %in% keep]
  cd[, colour := factor((ID - 1L) %% ncolors)]
  pd[, colour := factor((ID - 1L) %% ncolors)]
  cd[, grp := paste(byplot, ID, prob, piece)]
  cd[, byplot := factor(byplot, levels = c("Genotype", "Environment"))]
  pd[, byplot := factor(byplot, levels = c("Genotype", "Environment"))]
  lim <- max(abs(c(cd$x, cd$y, pd$x, pd$y, 0)), na.rm = TRUE) * 1.03
  pal <- rep_len(.bam_palette, ncolors)
  levs <- sprintf("%g%%", 100 * (1 - pb))
  cd[, prob := factor(prob, levels = levs)]
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey55", linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey55", linetype = "dashed", linewidth = 0.3)
  if (nrow(cd)) {
    p <- p + ggplot2::geom_polygon(data = cd, mapping = ggplot2::aes(x = x, y = y, group = grp, fill = colour,
                                                                     alpha = prob),
                                   colour = "grey25", linewidth = 0.2)
  }
  p + ggplot2::geom_point(data = pd, mapping = ggplot2::aes(x = x, y = y, colour = colour), size = 1.8) +
    ggrepel::geom_text_repel(data = pd, mapping = ggplot2::aes(x = x, y = y, label = name, colour = colour),
                             size = 3, family = "serif", max.overlaps = 100L, show.legend = FALSE) +
    ggplot2::facet_wrap(ggplot2::vars(byplot), ncol = 2L) +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::scale_alpha_manual(values = stats::setNames(seq(0.25, 0.08, length.out = length(levs)), levs),
                                name = "HPD region") +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggplot2::labs(x = expression(lambda[1]^{1/2} ~ "times first score"),
                  y = expression(lambda[2]^{1/2} ~ "times second score"),
                  title = "Bayesian AMMI biplot with highest posterior density regions") +
    bam_theme()
}

.bam_print_biplots <- function(x, digits = 3L, ...) {
  cat(attr(x, "label"), "\n\n")
  print(x$biplot_data, digits = digits, nrows = 500L)
  invisible(x)
}

.bam_plot_biplots <- function(x, ...) x$biplot

utils::globalVariables(c("stable", "ID", "byplot", "colour", "grp", "prob", "piece", "name", "x", "y"))
