#' @name    ge_ammi
#' @aliases ge_ammi
#' @title   AMMI of Genotype by Environment Interaction Model
#' @description Performs Additive Main Effects and Multiplication Interaction Analysis of Genotype by Environment Interaction Model
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#' @param .rep   Replication Factor
#'
#' @return Genotype by Environment Interaction Model
#'
#' @author
#' \enumerate{
#'     \item Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'    }
#'
#' @references
#'  Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J.M., Viele, K., Liu, G., and Cornelius, P.L. (2011)
#'  Bayesian Estimation of the Additive Main Effects and Multiplicative Interaction Model
#'  \emph{Crop Science}, \strong{51}, 1458–1469.
#'  (\href{https://dl.sciencesocieties.org/publications/cs/abstracts/51/4/1458?access=0&view=pdf}{doi: 10.2135/cropsci2010.06.0343})
#'
#' @import lme4
#' @import rlang
#' @import tidyr
#' @importFrom  dplyr select group_by count
#' @importFrom magrittr %>%
#' @importFrom stats sigma
#'
#' @export
#'
#' @examples
#'
#' data(Maiz)
#' fm1 <-
#'    ge_ammi(
#'       .data  = Maiz
#'      , .y    = y
#'      , .gen  = entry
#'      , .env  = site
#'      , .rep  = rep
#'      )
#'
#'
#'
#'

ge_ammi <- function(.data, .y, .gen, .env, .rep) {
  UseMethod("ge_ammi")
}


#' @export
#' @rdname ge_ammi

ge_ammi.default <-
  function(.data, .y, .gen, .env, .rep){

    .y    <- deparse(substitute(.y))
    .gen  <- deparse(substitute(.gen))
    .env  <- deparse(substitute(.env))
    .rep  <- deparse(substitute(.rep))

    df1 <- tibble::as_tibble(data.frame(
        Env = factor(.data[[.env]])
      , Gen = factor(.data[[.gen]])
      , Rep = factor(.data[[.rep]])
      , Y   = .data[[.y]]
    ))

    g   <- length(levels(df1$Gen))
    e   <- length(levels(df1$Env))
    Rep <- length(levels(df1$Rep))
    k   <- min(g, e) - 1

    ge_means <-
      ge_mean(
         .data = df1
        , .y    = Y
        , .gen  = Gen
        , .env  = Env
        )

    ge_effects <-
      ge_eff(
          .data = df1
        , .y    = Y
        , .gen  = Gen
        , .env  = Env
      )


    g_effects <-
      g_eff(
        .data = df1
        , .y    = Y
        , .gen  = Gen
        , .env  = Env
      )

    e_effects <-
      e_eff(
        .data = df1
        , .y    = Y
        , .gen  = Gen
        , .env  = Env
      )

    ge_model <-
      ge_model(
         .data = df1
        , .y    = Y
        , .gen  = Gen
        , .env  = Env
        , .rep  = Rep
      )

    return(list(
      g        = g
    , e        = e
    , Rep      = Rep
    , k        = k
    , ge_means = ge_means
    , mu       = ge_means$grand_mean
    , sigma2   = sigma(ge_model)^2
    , tau      = 1/sigma(ge_model)^2
    , tao      = g_effects
    , delta    = e_effects
    , lambdas  = ge_effects$ge_svd$d[1:k]
    , alphas   = ge_effects$ge_svd$u[,1:k]
    , gammas   = ge_effects$ge_svd$v[,1:k]
      ))
}
