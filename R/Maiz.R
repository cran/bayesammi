#' Maize multi-environment trial
#'
#' @description
#' Grain yield of 9 maize (\emph{Zea mays} L.) genotypes in 20 sites, a randomized complete block design with 4
#' replicates in every site: the example of Crossa et al. (2011), also used by Burgueno et al. (2008).
#'
#' @format A \code{data.frame} with 720 rows and 7 columns:
#' \describe{
#'   \item{\code{site}}{site (environment), 1 to 20.}
#'   \item{\code{entry}}{genotype, 1 to 9.}
#'   \item{\code{fold10}, \code{fold5}, \code{fold3}}{fold numbers of 10-, 5- and 3-fold cross-validation (by
#'     genotype).}
#'   \item{\code{rep}}{replicate within the site, 1 to 4.}
#'   \item{\code{y}}{grain yield (kg/ha).}
#' }
#' @details
#' The grand mean is 4858 kg/ha; the least-squares singular values of the interaction are 5923 and 3070 for the first
#' two terms (Crossa et al. 2011).
#'
#' \strong{Change from the earlier versions.} The help page described another data set (1320 rows; columns Gen,
#' Institute, Rep, Block, Env, Yield).
#' @source International Maize and Wheat Improvement Center (CIMMYT).
#' @references
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#'
#' Burgueno, J., Crossa, J., Cornelius, P. L. and Yang, R.-C. (2008). Using factor analytic models for joining
#' environments and genotypes without crossover genotype x environment interaction. \emph{Crop Science},
#' \strong{48}(4), 1291--1305.
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com}), Jose Crossa (\email{j.crossa@@cgiar.org})
#' @examples
#' data(Maiz)
#' str(Maiz)
#' ge_mean(
#'     .data = Maiz
#'   , .y    = y
#'   , .gen  = entry
#'   , .env  = site
#'   )$grand_mean
"Maiz"
