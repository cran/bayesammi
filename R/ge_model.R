#' @name    ge_model
#' @aliases ge_model
#' @title   Genotype by Environment Interaction Model
#' @description Calcuates Genotype by Environment Interaction Model
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
#'  (\href{https://acsess.onlinelibrary.wiley.com/doi/abs/10.2135/cropsci2010.06.0343}{doi: 10.2135/cropsci2010.06.0343})
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
#'    ge_model(
#'       .data  = Maiz
#'      , .y    = y
#'      , .gen  = entry
#'      , .env  = site
#'      , .rep  = rep
#'      )
#'
#'

ge_model <- function(.data, .y, .gen, .env, .rep) {
  UseMethod("ge_model")
}


#' @export
#' @rdname ge_model

ge_model.default <-
  function(.data, .y, .gen, .env, .rep){

    .y    <- deparse(substitute(.y))
    .gen  <- deparse(substitute(.gen))
    .env  <- deparse(substitute(.env))
    .rep  <- deparse(substitute(.rep))

    df1 <- data.frame(
            Y = .data[[.y]]
          , Env = factor(.data[[.env]])
          , Gen = factor(.data[[.gen]])
          , Rep = factor(.data[[.rep]])
          )

     ge_fm <-
       lme4::lmer(Y ~ Env + Gen + Env:Gen + (1|Env:Rep), data = df1)
      return(ge_fm)
}
