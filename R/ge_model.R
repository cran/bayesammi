#' Linear mixed model of a multi-environment trial
#'
#' @description
#' \code{ge_model()} fits the two-way model with interaction and random replicates within environments by restricted
#' maximum likelihood (\code{lme4::lmer()}). Its residual variance is the error variance of the least-squares AMMI
#' fit (\code{\link{ge_ammi}()}, which computes it in closed form).
#'
#' @details
#' \deqn{y_{ijl} = \mu + \delta_j + \tau_i + (\tau\delta)_{ij} + b_{l(j)} + e_{ijl}, \qquad b_{l(j)} \sim N(0,
#'   \sigma^2_b), \quad e_{ijl} \sim N(0, \sigma^2).}{y_ijl = mu + delta_j + tau_i + (tau delta)_ij + b_l(j) + e_ijl,
#'   b ~ N(0, sigma_b^2), e ~ N(0, sigma^2).}
#' The package \pkg{lme4} is needed (it is suggested, not required, since version 1.0.0).
#'
#' @param .data a data frame with one row per plot.
#' @param .y,.gen,.env,.rep the columns of the response, genotypes, environments and replicates: bare names, strings
#'   or vectors.
#' @return The \code{lmerMod} object of \code{lme4::lmer(Y ~ Env + Gen + Env:Gen + (1 | Env:Rep))}.
#' @references
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#'
#' Bates, D., Maechler, M., Bolker, B. and Walker, S. (2015). Fitting linear mixed-effects models using lme4.
#' \emph{Journal of Statistical Software}, \strong{67}(1), 1--48. \doi{10.18637/jss.v067.i01}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{ge_ammi}}
#' @examples
#' data(Maiz)
#' if (requireNamespace("lme4", quietly = TRUE)) {
#'   fm <-
#'     ge_model(
#'         .data = Maiz
#'       , .y    = y
#'       , .gen  = entry
#'       , .env  = site
#'       , .rep  = rep
#'       )
#'   stats::sigma(fm)^2
#' }
#' @export
ge_model <- function(.data, .y, .gen, .env, .rep) UseMethod("ge_model")

#' @rdname ge_model
#' @export
ge_model.default <- function(.data, .y, .gen, .env, .rep) {
  fn <- "ge_model"
  if (!requireNamespace("lme4", quietly = TRUE)) .bam_stop(fn, "the package lme4 is needed: install.packages(\"lme4\").")
  cd <- .bam_prep(fn, .data, substitute(.y), substitute(.gen), substitute(.env), substitute(.rep), parent.frame())
  d  <- data.frame(Y = cd$data$y, Env = cd$data$Env, Gen = cd$data$Gen, Rep = cd$data$Rep)
  lme4::lmer(Y ~ Env + Gen + Env:Gen + (1 | Env:Rep), data = d)
}
