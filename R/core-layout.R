#-------------------------------------------------------------------------------------------------------------------
# The layout of the draws: one row per parameter (name, block, label, axis), the same for the summary table, the
# coda draws and every function that reads them, whatever the model (AMMI or GGE) and the variances (one or one per
# environment).
#-------------------------------------------------------------------------------------------------------------------
.bam_layout <- function(G, E, K, gge = FALSE, het = FALSE) {
  kk  <- seq_len(K)
  sig <- if (het) sprintf("sigma2[%s]", E) else "sigma2"
  tau <- if (gge) character(0) else sprintf("tau[%s]", G)
  data.table::data.table(
      parameter = c("mu", sig, tau, sprintf("delta[%s]", E), sprintf("lambda[%d]", kk),
                    sprintf("alpha[%s,%d]", rep(G, K), rep(kk, each = length(G))),
                    sprintf("gamma[%s,%d]", rep(E, K), rep(kk, each = length(E))))
    , block = c("mu", rep("sigma2", length(sig)), rep("tau", length(tau)), rep("delta", length(E)), rep("lambda", K),
                rep("alpha", length(G) * K), rep("gamma", length(E) * K))
    , label = c("", if (het) E else "", if (gge) character(0) else G, E, rep("", K), rep(G, K), rep(E, K))
    , axis  = c(rep(NA_integer_, 1L + length(sig) + length(tau) + length(E)), kk, rep(kk, each = length(G)),
                rep(kk, each = length(E)))
    )
}

# all draws of a fit, chains one after the other (matrix with the parameter names as columns)
.bam_all <- function(fit) do.call(rbind, lapply(fit$draws, unclass))

# the draws of one block: a matrix (draws x labels) or, for the scores, an array (draws x labels x axes)
.bam_block <- function(all, lay, block) {
  i <- which(lay$block == block)
  M <- all[, lay$parameter[i], drop = FALSE]
  if (block %in% c("alpha", "gamma")) {
    K <- max(lay$axis[i])
    return(array(M, c(nrow(M), length(i) / K, K)))
  }
  M
}

# the draws of the cell means eta_ij of environment j (draws x genotypes)
.bam_eta_env <- function(all, lay, j, G) {
  E   <- sum(lay$block == "delta")
  K   <- sum(lay$block == "lambda")
  mu  <- all[, "mu"]
  tau <- if (any(lay$block == "tau")) .bam_block(all, lay, "tau") else matrix(0, nrow(all), G)
  del <- .bam_block(all, lay, "delta")[, j]
  lam <- .bam_block(all, lay, "lambda")
  A   <- .bam_block(all, lay, "alpha")
  Gm  <- .bam_block(all, lay, "gamma")
  out <- mu + tau + del
  for (k in seq_len(K)) out <- out + A[, , k] * (lam[, k] * Gm[, j, k])
  out
}
