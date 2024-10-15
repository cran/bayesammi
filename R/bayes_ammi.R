#' @name    bayes_ammi
#' @aliases bayes_ammi
#' @title   Bayesian Estimation of the Additive Main Effects and Multiplicative Interaction Model
#' @description Performs Bayesian Estimation of the Additive Main Effects and Multiplicative Interaction Model
#'
#' @param .data  data.frame
#' @param .y     Response Variable
#' @param .gen   Genotypes Factor
#' @param .env   Environment Factor
#' @param .rep   Replication Factor
#' @param .nIter Number of Iterations
#'
#' @return Genotype by Environment Interaction Model
#'
#' @author
#' \enumerate{
#'     \item Muhammad Yaseen (\email{myaseen208@gmail.com})
#'     \item Jose Crossa (\email{j.crossa@cgiar.org})
#'     \item Sergio Perez-Elizalde (\email{sergiop@colpos.mx})
#'     \item Diego Jarquin (\email{diego.jarquin@gmail.com})
#'     \item Jose Miguel Cotes
#'     \item Kert Viele
#'     \item Genzhou Liu
#'     \item Paul L. Cornelius
#'    }
#'
#' @references
#'  Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J.M., Viele, K., Liu, G., and Cornelius, P.L. (2011)
#'  Bayesian Estimation of the Additive Main Effects and Multiplicative Interaction Model
#'  \emph{Crop Science}, \strong{51}, 1458–1469.
#'  (\href{https://acsess.onlinelibrary.wiley.com/doi/abs/10.2135/cropsci2010.06.0343}{doi: 10.2135/cropsci2010.06.0343})
#'
#'
#' @import ggplot2
#' @import lme4
#' @import mvtnorm
#' @import rlang
#' @import rstiefel
#' @import scales
#' @import tidyr
#' @import tibble
#' @import tmvtnorm
#' @importFrom  dplyr select group_by count
#' @importFrom magrittr %>%
#' @importFrom MASS Null mvrnorm
#' @importFrom stats rexp rgamma rnorm runif
#'
#' @export
#'
#' @examples
#'
#' data(Maiz)
#'
#' fm1 <-
#'  bayes_ammi(
#'      .data  = Maiz
#'    , .y     = y
#'    , .gen   = entry
#'    , .env   = site
#'    , .rep   = rep
#'    , .nIter = 20
#'   )
#' names(fm1)
#' fm1$mu1
#' fm1$tau1
#' fm1$tao1
#' fm1$delta1
#' fm1$lambdas1
#' fm1$alphas1
#' fm1$gammas1
#'
#' library(ggplot2)
#'
#' Plot1Mu <-
#'   ggplot(data = fm1$mu1, mapping = aes(x = 1:nrow(fm1$mu1), y = mu)) +
#'   geom_line(color = "blue") +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = expression(mu), x = "Iterations") +
#'   theme_bw()
#' print(Plot1Mu)
#'
#' Plot2Mu <-
#'   ggplot(data = fm1$mu1, mapping = aes(mu)) +
#'   geom_histogram() +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = "Frequency", x = expression(mu)) +
#'   theme_bw()
#' print(Plot2Mu)
#'
#'
#' Plot1Sigma2 <-
#'   ggplot(data = fm1$tau1, mapping = aes(x = 1:nrow(fm1$tau1), y = tau)) +
#'   geom_line(color = "blue") +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = expression(sigma^2), x = "Iterations") +
#'   theme_bw()
#' print(Plot1Sigma2)
#'
#'
#' Plot2Sigma2 <-
#'   ggplot(data = fm1$tau1, mapping = aes(tau)) +
#'   geom_histogram() +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = "Frequency", x = expression(sigma^2)) +
#'   theme_bw()
#' print(Plot2Sigma2)
#'
#'
#' # Plot of Alphas
#' Plot1Alpha1 <-
#'   ggplot(data = fm1$tao1, mapping = aes(x = 1:nrow(fm1$tao1), y = tao1)) +
#'   geom_line(color = "blue") +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = expression(alpha[1]), x = "Iterations") +
#'   theme_bw()
#' print(Plot1Alpha1)
#'
#' Plot2Alpha1 <-
#'   ggplot(data = fm1$tao1, mapping = aes(tao1)) +
#'   geom_histogram() +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = "Frequency", x = expression(alpha[1])) +
#'   theme_bw()
#' print(Plot2Alpha1)
#'
#' Plot1Alpha2 <-
#'   ggplot(data = fm1$tao1, mapping = aes(x = 1:nrow(fm1$tao1), y = tao2)) +
#'   geom_line(color = "blue") +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = expression(alpha[2]), x = "Iterations") +
#'   theme_bw()
#' print(Plot1Alpha2)
#'
#' Plot2Alpha2 <-
#'   ggplot(data = fm1$tao1, mapping = aes(tao2)) +
#'   geom_histogram() +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = "Frequency", x = expression(alpha[2])) +
#'   theme_bw()
#' print(Plot2Alpha2)
#'
#' # Plot of Betas
#' Plot1Beta1 <-
#'   ggplot(data = fm1$delta1, mapping = aes(x = 1:nrow(fm1$delta1), y = delta1)) +
#'   geom_line(color = "blue") +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = expression(beta[1]), x = "Iterations") +
#'   theme_bw()
#' print(Plot1Beta1)
#'
#' Plot2Beta1 <-
#'   ggplot(data = fm1$delta1, mapping = aes(delta1)) +
#'   geom_histogram() +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = "Frequency", x = expression(beta[1])) +
#'   theme_bw()
#' print(Plot2Beta1)
#'
#'
#' Plot1Beta2 <-
#'   ggplot(data = fm1$delta1, mapping = aes(x = 1:nrow(fm1$delta1), y = delta2)) +
#'   geom_line(color = "blue") +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = expression(beta[2]), x = "Iterations") +
#'   theme_bw()
#' print(Plot1Beta2)
#'
#' Plot2Beta2 <-
#'   ggplot(data = fm1$delta1, mapping = aes(delta2)) +
#'   geom_histogram() +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = "Frequency", x = expression(beta[2])) +
#'   theme_bw()
#' print(Plot2Beta2)
#'
#'
#' Plot1Beta3 <-
#'   ggplot(data = fm1$delta1, mapping = aes(x = 1:nrow(fm1$delta1), y = delta3)) +
#'   geom_line(color = "blue") +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = expression(beta[3]), x = "Iterations") +
#'   theme_bw()
#' print(Plot1Beta3)
#'
#' Plot2Beta3 <-
#'   ggplot(data = fm1$delta1, mapping = aes(delta3)) +
#'   geom_histogram() +
#'   scale_x_continuous(labels = scales::comma) +
#'   scale_y_continuous(labels = scales::comma) +
#'   labs(y = "Frequency", x = expression(beta[3])) +
#'   theme_bw()
#' print(Plot2Beta3)
#'
#'
#' BiplotAMMI <-
#'   ggplot(data = fm1$alphas0, mapping = aes(x = alphas1, y = alphas2)) +
#'   geom_point() +
#'   geom_hline(yintercept = 0) +
#'   geom_vline(xintercept = 0) +
#'   geom_text(aes(label = 1:nrow(fm1$alphas0)),
#'             vjust = "inward", hjust = "inward") +
#'   geom_point(data = fm1$gammas0, mapping = aes(x = gammas1, y = gammas2)) +
#'   geom_segment(data = fm1$gammas0,
#'                aes(x = 0, y = 0, xend = gammas1, yend = gammas2),
#'                arrow = arrow(length = unit(0.2, "cm"))
#'                , alpha = 0.75, color = "red") +
#'   geom_text(data = fm1$gammas0,
#'             aes(x = gammas1, y = gammas2,
#'                 label = paste0("E", 1:nrow(fm1$gammas0))),
#'             vjust = "inward", hjust = "inward") +
#'   scale_x_continuous(
#'     limits = c(-max(abs(c(range(fm1$alphas0[, 1:2], fm1$gammas0[, 1:2]))))
#'                , max(abs(c(range(fm1$alphas0[, 1:2], fm1$gammas0[, 1:2])))))) +
#'   scale_y_continuous(
#'     limits = c(-max(abs(c(range(fm1$alphas0[, 1:2], fm1$gammas0[, 1:2]))))
#'                , max(abs(c(range(fm1$alphas0[, 1:2], fm1$gammas0[, 1:2])))))) +
#'   labs(title = "MCO Method", x = expression(PC[1]), y = expression(PC[2])) +
#'   theme_bw() +
#'   theme(plot.title = element_text(hjust = 0.5))
#'
#' print(BiplotAMMI)
#'
#'
#' BiplotBayesAMMI <-
#'   ggplot(data = fm1$alphas1, mapping = aes(x = alphas1, y = alphas2)) +
#'   geom_point() +
#'   geom_hline(yintercept = 0) +
#'   geom_vline(xintercept = 0) +
#'   geom_text(aes(label = 1:nrow(fm1$alphas1)),
#'             vjust = "inward", hjust = "inward") +
#'   geom_point(data = fm1$gammas1, mapping = aes(x = gammas1, y = gammas2)) +
#'   geom_segment(data = fm1$gammas1,
#'                aes(x = 0, y = 0, xend = gammas1, yend = gammas2),
#'                arrow = arrow(length = unit(0.2, "cm"))
#'                , alpha = 0.75, color = "red") +
#'   geom_text(data = fm1$gammas1,
#'             aes(x = gammas1, y = gammas2,
#'                 label = paste0("E", 1:nrow(fm1$gammas1))),
#'             vjust = "inward", hjust = "inward") +
#'   scale_x_continuous(
#'     limits = c(-max(abs(c(range(fm1$alphas1[, 1:2], fm1$gammas1[, 1:2]))))
#'                , max(abs(c(range(fm1$alphas1[, 1:2], fm1$gammas1[, 1:2])))))) +
#'   scale_y_continuous(
#'     limits = c(-max(abs(c(range(fm1$alphas1[, 1:2], fm1$gammas1[, 1:2]))))
#'                , max(abs(c(range(fm1$alphas1[, 1:2], fm1$gammas1[, 1:2])))))) +
#'   labs(title = "Bayesian Method", x = expression(PC[1]), y = expression(PC[2])) +
#'   theme_bw() +
#'   theme(plot.title = element_text(hjust = 0.5))
#'
#' print(BiplotBayesAMMI)
#'
#'
# if(getRversion() >= "2.15.1"){
#   utils::globalVariables(
#     c(
#       "alphas0"
#       , "c0"
#       , "delta0"
#       , "gammas0"
#       , "lambdas0"
#       , "ge_means0"
#       , "mu0"
#       , "n0"
#       , "tao0"
#       , "tau0"
#     )
#   )
# }
#


bayes_ammi <- function(.data, .y, .gen, .env, .rep, .nIter){
  UseMethod("bayes_ammi")
}


#' @export
#' @rdname bayes_ammi


bayes_ammi.default <- function(.data, .y, .gen, .env, .rep, .nIter){

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


fm1   <-
    ge_ammi(
          .data = df1
        , .y    = Y
        , .gen  = Gen
        , .env  = Env
        , .rep  = Rep
      )

  Envs         <- df1$Env
  Rep          <- fm1$Rep
  alphas       <- fm1$alphas
  gammas       <- fm1$gammas
  ge_means     <- fm1$ge_means$ge_means
  sigma2       <- fm1$sigma2
  tau          <- fm1$tau
  tao          <- fm1$tao
  lambdas      <- fm1$lambdas
  delta        <- fm1$delta
  mu           <- fm1$ge_means$grand_mean
  g            <- fm1$g
  e            <- fm1$e
  k            <- min(c(g, e)) - 1
  #

  ge_variances <- ge_var(.data = df1, .y = Y, .gen = Gen, .env = Env)$ge_variances
  ssy                 <- sum(ge_variances, na.rm = TRUE)

  alphas0             <- alphas
  gammas0             <- gammas
  lambdas0            <- lambdas
  grand_mean          <- mu

  kg <- t(matrix_k(g))
  ke <- t(matrix_k(e))

  # Values for hyper-parameters
  sigma2_u    <- 1*10^15
  alfa_sig_ap <- 0
  beta_sig_ap <- 0
  mu_tao_ap   <- rep(0, g)
  sig_tao_ap  <- 1*10^15
  mu_del_ap   <- rep(0, e)
  sig_del_ap  <- 1*10^15
  mu_lamb_ap  <- 0
  sig_lamb_ap <- 1*10^15

  g_means  <- matrix(rowMeans(ge_means), ncol = 1)
  tao      <- matrix(g_means - mean(ge_means), ncol = 1)
  e_means  <- colMeans(ge_means)
  delta    <- matrix(e_means - mean(ge_means), ncol = 1)

  phi          <- 0
  uni          <- 1
  med_lamb     <- (Rep*tau*(sum((alphas %*% t(gammas))*ge_means)) + Rep*tau*lambdas0)/(Rep*tau + Rep*tau)
  u_neg        <- - med_lamb/sqrt((Rep*tau + Rep*tau)^(-1))
  u_pos        <- (lambdas[-length(lambdas)] - med_lamb[-1])/sqrt((Rep*tau + Rep*tau)^(-1))
  alpha_star   <- (u_neg + sqrt(u_neg^2 + 4))/2
  z            <- c(runif(n = 1, min = u_neg, max = u_pos), rexp(n = (k-1), rate = alpha_star) + u_neg[-1])
  lambdas_m    <- sqrt((Rep*tau + Rep*tau)^(-1))*z + med_lamb

  suma_parci3    <- function(lambdas){
    suma_parcial <- alphas %*% diag(lambdas) %*% t(gammas)
    return(suma_parcial)
  }

  # suma_parci3(lambdas)

  matriz_desv_model <- function(nada){
    iter  <- 1
    df1se2 <- NULL
    while(iter <= g) {
      df1s <- df1[df1$Gen == iter, ]
      iter2 <- 1
      df1se1 <- NULL
      while(iter2 <= e) {
        df1se  <- ((mean(df1s[df1s$Env == Envs[iter2], ]$Y) - mu - tao[iter] - delta[iter2] - suma_parci3(lambdas)[iter,iter2])^2)
        df1se1 <- cbind(df1se1, df1se)
        iter2  <- iter2+1
      }
      df1se2 <- rbind(df1se2, df1se1)
      iter   <- iter + 1
    }
    return(df1se2)
  }

  # matriz_desv_model(1)

  rmf.matrix.gibbs4 <- function(M, X) {
    rscol <- max(2, min(round(log(dim(M)[1])), dim(M)[2]))
    sM    <- svd(M)
    H     <- sM$u %*% diag(sM$d)
    Y     <- X %*% sM$v
    m     <- dim(H)[1]
    R     <- dim(H)[2]

   for (iter in 1:round(R/rscol)) {
      r <- 1:R
      N <- rstiefel::NullC(Y[, -r])
      y <- rstiefel::rmf.matrix(t(N) %*% H[, r])
      Y[, r] <- N %*% y
   }

    entrada  <- Y %*% t(sM$v)
    entrada1 <- cbind(rep(1, m), entrada)
    salida2  <- orthnorm(entrada1, basis = TRUE, norm = TRUE)
    salida3  <- salida2[ ,2:(k+1)]
    return(salida3)
  }

  mu1      <- NULL
  tau1     <- NULL
  tao1     <- NULL
  delta1   <- NULL
  lambdas1 <- NULL
  alphas1  <- NULL
  gammas1  <- NULL


    times <- 1
    for(time in 1:.nIter) {
    mu    <- rnorm(
               n    = 1
             , mean = (Rep*g*e*sigma2_u*grand_mean)/(Rep*g*e*sigma2_u + sigma2)
             , sd   = sqrt((sigma2*sigma2_u)/(Rep*g*e*sigma2_u+sigma2))
             )

    tao   <- kg %*%
                MASS::mvrnorm(
                 n      = 1
               , mu     = (Rep*e*sig_tao_ap*t(kg)%*%(g_means) + sigma2*t(kg)%*%mu_tao_ap)/(Rep*e*sig_tao_ap+sigma2)
               , Sigma  = sqrt((sig_tao_ap*sigma2)/(Rep*e*sig_tao_ap+sigma2)*diag(rep(1, g-1)))
               )

    delta <- ke %*%
                MASS::mvrnorm(
                  n     = 1
                , mu    = (Rep*e*sig_del_ap*t(ke)%*%(e_means) + sigma2*t(ke)%*%mu_del_ap)/(Rep*e*sig_del_ap+sigma2)
                , Sigma = sqrt((sig_del_ap*sigma2)/(Rep*e*sig_del_ap+sigma2)*diag(rep(1, e-1)))
                )

    mean_lambI <- diag(t(gammas) %*% t(ge_means) %*% alphas)
    lambdas    <- lambdas_m
    alphas     <- rmf.matrix.gibbs4(M = (Rep * ge_means %*% gammas %*% diag(lambdas))/sigma2,    X = alphas)[ ,1:k]
    gammas     <- rmf.matrix.gibbs4(M = (Rep * t(ge_means) %*% alphas %*% diag(lambdas))/sigma2, X = gammas)[ ,1:k]
    sigma2     <- 1/rgamma(1, Rep*g*e/2 + alfa_sig_ap, .5*(Rep*sum(matriz_desv_model(1))+ssy)+beta_sig_ap)


    for(i in 1:ncol(alphas)) {
      maxal0 <- which.max(abs(alphas0[ ,i]))
      maxga0 <- which.max(abs(gammas0[ ,i]))
      sgn1 <- (alphas0[ ,i] < 0) != (alphas[ ,i] < 0)
      sgn2 <- (gammas0[ ,i] < 0) != (gammas[ ,i] < 0)
      sgn3 <- (alphas0[maxal0,i]<0)!=(alphas[maxal0,i]<0) && (abs(alphas[maxal0,i])>19.5/20)
      sgn4 <- (gammas0[maxga0,i]<0)!=(gammas[maxga0,i]<0) && (abs(alphas[maxal0,i])>19.5/20)
      if(any(all(sgn1), all(sgn2), sgn3, sgn4)){
        alphas[ ,i] <- -1*alphas[ ,i]
        gammas[ ,i] <- -1*gammas[ ,i]
      }
    }


    alphas2  <- matrix(alphas, nrow = 1)
    gammas2  <- matrix(gammas, nrow = 1)

    mu1      <-  rbind(mu1, mu)
    tau1     <-  rbind(tau1, tau)
    tao1     <-  cbind(tao1, tao)
    delta1   <-  cbind(delta1, delta)
    lambdas1 <-  cbind(lambdas1, lambdas)
    alphas1  <-  matrix(alphas2, ncol = k, byrow = TRUE)
    gammas1  <-  matrix(gammas2, ncol = k, byrow = TRUE)

     times   <- times + 1
     }

    mu1      <-  tibble::as_tibble(mu1)
    tau1     <-  tibble::as_tibble(tau1)
    tao1     <-  tibble::as_tibble(t(tao1))
    delta1   <-  tibble::as_tibble(t(delta1))
    lambdas1 <-  tibble::as_tibble(t(lambdas1))
    alphas1  <-  tibble::as_tibble(alphas1)
    gammas1  <-  tibble::as_tibble(gammas1)

    colnames(mu1)      <- c("mu")
    colnames(tau1)     <- c("tau")
    colnames(tao1)     <- paste0("tao", 1:g)
    colnames(delta1)   <- paste0("delta", 1:e)
    colnames(lambdas1) <- paste0("lambdas", 1:k)
    colnames(alphas1)  <- paste0("alphas", 1:k)
    colnames(gammas1)  <- paste0("gammas", 1:k)

    alphas0  <-  tibble::as_tibble(alphas0)
    gammas0  <-  tibble::as_tibble(gammas0)

    colnames(alphas0)  <- paste0("alphas", 1:k)
    colnames(gammas0)  <- paste0("gammas", 1:k)

    return(list(
        mu1         = mu1
      , tau1        = tau1
      , tao1        = tao1
      , delta1      = delta1
      , lambdas1    = lambdas1
      , alphas1     = alphas1
      , gammas1     = gammas1
      , alphas0     = alphas0
      , gammas0     = gammas0
    ))
}
