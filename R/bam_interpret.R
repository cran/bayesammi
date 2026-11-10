#' Interpretation of a Bayesian AMMI analysis
#'
#' @description
#' \code{bam_interpret()} turns a result of \code{\link{bayes_ammi}()}, \code{\link{biplots}()} or
#' \code{\link{ge_ammi}()} into a plain-language interpretation: whether the chains have converged, which genotypes
#' and environments differ from the average (main effects whose HPD interval excludes zero), how strong the bilinear
#' terms are, which genotypes and environments interact significantly and in which direction, which are stable, and
#' what to check before acting on the results. Every statement is computed from the numbers of the result.
#'
#' @details
#' The rules:
#' \itemize{
#'   \item convergence: potential scale reduction factor below 1.1 (two or more chains), effective sample size of
#'     at least 400, Raftery-Lewis dependence factor below 5 and the Heidelberger-Welch tests passed (Crossa et
#'     al. 2011);
#'   \item main effects: a genotype or environment differs from the average when the HPD interval of its effect
#'     excludes zero;
#'   \item bilinear terms: a term whose HPD interval of \eqn{\lambda_k}{lambda_k} reaches down to near zero (lower
#'     limit below a tenth of the posterior mean) is weakly supported;
#'   \item interaction: a genotype (environment) interacts significantly on the first term when the HPD interval of
#'     its first score excludes zero; genotypes and environments with significant scores of the same sign interact
#'     positively (the genotype does better there than its main effect predicts), of opposite sign negatively;
#'   \item stability (\code{\link{biplots}()}): a genotype is stable when its HPD region contains the origin.
#' }
#'
#' @param x a result of \code{\link{bayes_ammi}()}, \code{\link{biplots}()} or \code{\link{ge_ammi}()}.
#' @param higher_better \code{TRUE} (default) when larger values of the trait are better.
#' @param top number of genotypes or environments named in lists (default 5).
#' @param ... unused.
#' @return An object of class \code{"bam_interpretation"}: a list with \code{analysis}, \code{title},
#'   \code{headline}, \code{findings}, \code{recommendations} and \code{cautions} (sentences), \code{facts} (the
#'   numbers and lists behind the sentences), \code{guide} (how the method is read) and \code{sources} (the full
#'   references of the method; \code{\link{bam_cite}()} gives them as BibTeX). \code{print()} shows the text;
#'   \code{format(x, style = "markdown")} returns it as Markdown.
#' @references
#' Crossa, J., Perez-Elizalde, S., Jarquin, D., Cotes, J. M., Viele, K., Liu, G. and Cornelius, P. L. (2011).
#' Bayesian estimation of the additive main effects and multiplicative interaction model. \emph{Crop Science},
#' \strong{51}(4), 1458--1469. \doi{10.2135/cropsci2010.06.0343}
#' @author Muhammad Yaseen (\email{myaseen208@@gmail.com})
#' @seealso \code{\link{bayes_ammi}}, \code{\link{biplots}}
#' @examples
#' data(Maiz)
#' fm1 <-
#'   bayes_ammi(
#'       .data   = Maiz
#'     , .y      = y
#'     , .gen    = entry
#'     , .env    = site
#'     , .rep    = rep
#'     , .nIter  = 500
#'     , .k      = 2
#'     , .burnin = 100
#'     , .thin   = 1
#'     , .chains = 2
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
#' txt <-
#'   bam_interpret(
#'       x             = fm1
#'     , higher_better = TRUE
#'     , top           = 5
#'     )
#' txt
#' cat(format(x = txt, style = "markdown"))
#' @export
bam_interpret <- function(x, higher_better = TRUE, top = 5L, ...) {
  fn  <- "bam_interpret"
  higher_better <- .bam_check_flag(higher_better, "higher_better", fn)
  top <- .bam_check_number(top, "top", fn, min = 1, integer = TRUE)
  if (!inherits(x, "bam_result")) .bam_stop(fn, "`x` must be a result of a function of the bayesammi package.")
  cls <- .bam_cls(x)
  f <- get0(paste0(".bam_int_", cls), envir = asNamespace("bayesammi"), inherits = FALSE)
  if (is.null(f)) .bam_stop(fn, "no interpretation is defined for ", cls, " results.")
  ctx <- list(hb = higher_better, top = top, good = if (higher_better) "above" else "below",
              bad = if (higher_better) "below" else "above")
  out <- f(x, ctx)
  out$analysis <- cls
  out$guide    <- .bam_guides[[cls]]
  out$sources  <- unname(vapply(.bam_refs_list()[unique(.bam_sources_of(x))], `[[`, "", "text"))
  structure(out[c("analysis", "title", "headline", "findings", "recommendations", "cautions", "facts", "guide",
                  "sources")],
            class = "bam_interpretation")
}

#' @rdname bam_interpret
#' @export
print.bam_interpretation <- function(x, ...) {
  cat(format(x, style = "text"), sep = "\n")
  invisible(x)
}

#' @rdname bam_interpret
#' @param style \code{"text"} (default) or \code{"markdown"}.
#' @param width line width of the text (default \code{getOption("width")}).
#' @export
format.bam_interpretation <- function(x, style = c("text", "markdown"), width = getOption("width", 80L), ...) {
  style <- match.arg(style)
  wrap <- function(s, prefix) {
    if (!length(s)) return(character(0))
    unlist(lapply(s, function(z) {
      w <- strwrap(z, width = max(40L, width - 2L) - nchar(prefix))
      c(paste0(prefix, w[1L]), if (length(w) > 1L) paste0(strrep(" ", nchar(prefix)), w[-1L]))
    }))
  }
  sec <- function(h, s) if (length(s)) c("", if (style == "markdown") paste("###", h) else paste0(h, ":"),
                                          wrap(s, "- ")) else character(0)
  c(if (style == "markdown") paste("##", x$title) else x$title,
    "", wrap(x$headline, ""),
    sec("Findings", x$findings), sec("Recommendations", x$recommendations), sec("Cautions", x$cautions),
    sec("Sources", x$sources))
}

.bam_list <- function(x, max = Inf) {
  x <- as.character(x)
  if (!length(x)) return("none")
  more <- length(x) - max
  if (more > 0) x <- c(x[seq_len(max)], sprintf("%d more", more))
  if (length(x) == 1L) return(x)
  paste0(paste(x[-length(x)], collapse = ", "), " and ", x[length(x)])
}
.bam_n <- function(v, d = 4L) {
  vapply(v, function(z) format(signif(z, d), big.mark = ",", trim = TRUE, scientific = FALSE), "")
}

# convergence facts of a bayes_ammi result
.bam_conv_facts <- function(x) {
  s  <- x$summary
  sc <- s[block %in% c("mu", "sigma2", "tau", "delta", "lambda", "slope", "sigma2_d")]
  dg <- x$diagnostics
  rl <- if (!is.null(dg$raftery)) {
    unlist(lapply(dg$raftery, function(r) if (is.matrix(r$resmatrix)) r$resmatrix[, "I"] else NA_real_))
  }
  hw <- if (!is.null(dg$heidel)) {
    unlist(lapply(dg$heidel, function(h) c(h[, "stest"], h[, "htest"])))
  }
  list(max_rhat = if (all(is.na(s$rhat))) NA_real_ else max(s$rhat, na.rm = TRUE),
       min_ess = min(s$ess), min_ess_scalar = min(sc$ess), max_dependence = if (any(is.finite(rl))) max(rl, na.rm = TRUE) else NA_real_,
       heidel_failed = if (length(hw)) sum(hw == 0, na.rm = TRUE) else NA_integer_,
       heidel_tests = if (length(hw)) sum(!is.na(hw)) else NA_integer_, chains = x$settings$chains,
       draws = x$settings$kept * x$settings$chains)
}

.bam_int_bayes_ammi <- function(x, ctx) {
  s   <- x$summary
  K   <- x$settings$k
  lev <- 100 * x$settings$level
  cf  <- .bam_conv_facts(x)
  # Heidelberger-Welch: each test has level 0.05, so a few of the many tests fail by chance; more than 10% is a
  # sign of non-stationary chains
  hw_bad <- !is.na(cf$heidel_failed) && cf$heidel_tests > 0 && cf$heidel_failed / cf$heidel_tests > 0.1
  ok  <- (is.na(cf$max_rhat) || cf$max_rhat < 1.1) && cf$min_ess_scalar >= 400 &&
    (is.na(cf$max_dependence) || cf$max_dependence < 5) && !hw_bad
  lam <- s[block == "lambda"]
  weak <- lam$axis[lam$lower < 0.1 * lam$mean]
  sig_pos <- function(b) s[block == b & lower > 0, label]
  sig_neg <- function(b) s[block == b & upper < 0, label]
  gen_good <- if (ctx$hb) sig_pos("tau") else sig_neg("tau")
  gen_bad  <- if (ctx$hb) sig_neg("tau") else sig_pos("tau")
  env_hi   <- sig_pos("delta")
  env_lo   <- sig_neg("delta")
  a1 <- s[block == "alpha" & axis == 1L]
  g1 <- s[block == "gamma" & axis == 1L]
  ga_pos <- a1[lower > 0, label]
  ga_neg <- a1[upper < 0, label]
  ge_pos <- g1[lower > 0, label]
  ge_neg <- g1[upper < 0, label]
  stable_good <- setdiff(gen_good, c(ga_pos, ga_neg))
  mu  <- s[block == "mu"]
  s2  <- s[block == "sigma2"]
  cv  <- 100 * sqrt(s2$mean) / mu$mean
  f <- c(cf, list(converged = ok, lambda = stats::setNames(lam$mean, paste0("lambda", lam$axis)),
                  lambda_lower = lam$lower, weak_terms = weak, genotypes_above = gen_good,
                  genotypes_below = gen_bad, environments_high = env_hi, environments_low = env_lo,
                  axis1_genotypes_positive = ga_pos, axis1_genotypes_negative = ga_neg,
                  axis1_environments_positive = ge_pos, axis1_environments_negative = ge_neg,
                  stable_and_good = stable_good, cv = cv))
  conv <- if (ok) {
    sprintf("The chains have converged by the criteria of Crossa et al. (2011)%s; the smallest effective sample size of the main parameters is %s.",
            if (is.na(cf$max_rhat)) " that apply to one chain" else sprintf(" (largest potential scale reduction factor %.3f)", cf$max_rhat),
            .bam_n(cf$min_ess_scalar, 3L))
  } else {
    sprintf("Convergence is not established: %s.", paste(c(
        if (!is.na(cf$max_rhat) && cf$max_rhat >= 1.1) sprintf("potential scale reduction factor up to %.2f", cf$max_rhat)
      , if (cf$min_ess_scalar < 400) sprintf("effective sample size as small as %s", .bam_n(cf$min_ess_scalar, 3L))
      , if (!is.na(cf$max_dependence) && cf$max_dependence >= 5) sprintf("Raftery-Lewis dependence factor up to %.1f", cf$max_dependence)
      , if (hw_bad) sprintf("%d of %d Heidelberger-Welch tests failed", cf$heidel_failed, cf$heidel_tests)
      ), collapse = "; "))
  }
  inter <- c(
      sprintf("On the first bilinear term, genotypes with significantly positive scores: %s; negative: %s.",
              .bam_list(ga_pos, ctx$top), .bam_list(ga_neg, ctx$top))
    , sprintf("Environments with significantly positive first scores: %s; negative: %s.", .bam_list(ge_pos, ctx$top),
              .bam_list(ge_neg, ctx$top))
    , if (length(c(ga_pos, ga_neg)) && length(c(ge_pos, ge_neg)))
      "Genotypes and environments with significant first scores of the same sign interact positively: the genotype yields more there than its main effect predicts; of opposite sign, less.")
  list(title = "Bayesian AMMI analysis",
       headline = sprintf("%d bilinear term(s): lambda_1 = %s%s; %s", K, .bam_n(lam$mean[1L]),
                          if (K > 1L) sprintf(", lambda_2 = %s", .bam_n(lam$mean[2L])) else "",
                          if (ok) "the chains have converged." else "convergence is not yet established."),
       findings = c(
           conv
         , if (!is.na(cf$heidel_failed)) sprintf("%d of %d Heidelberger-Welch stationarity and half-width tests failed%s.",
                                                 cf$heidel_failed, cf$heidel_tests,
                                                 if (hw_bad) "" else " (about 5% fail by chance at level 0.05)")
         , sprintf("Overall mean %s (%g%% HPD %s to %s); error variance %s, a coefficient of variation of %.1f%%.",
                   .bam_n(mu$mean), lev, .bam_n(mu$lower), .bam_n(mu$upper), .bam_n(s2$mean), cv)
         , sprintf("Genotypes significantly %s average (HPD of the effect excludes zero): %s; %s average: %s.", ctx$good,
                   .bam_list(gen_good, ctx$top), ctx$bad, .bam_list(gen_bad, ctx$top))
         , sprintf("Environments with significantly high means: %s; low: %s.", .bam_list(env_hi, ctx$top),
                   .bam_list(env_lo, ctx$top))
         , sprintf("Singular values with %g%% HPD intervals: %s.", lev,
                   paste(sprintf("lambda_%d = %s (%s to %s)", lam$axis, .bam_n(lam$mean), .bam_n(lam$lower),
                                 .bam_n(lam$upper)), collapse = "; "))
         , inter),
       recommendations = c(
           if (length(stable_good)) sprintf("%s are %s average and do not interact significantly on the first term: candidates for wide recommendation.",
                                            .bam_list(stable_good, ctx$top), ctx$good)
         , if (length(c(ga_pos, ga_neg))) "Recommend the genotypes with significant first scores for the environments whose scores have the same sign (biplots() shows the regions)."
         , "Draw biplots() to see which genotypes and environments form homogeneous groups (overlapping HPD regions)."),
       cautions = c(
           if (!ok) "Run longer chains (and two or more chains) before relying on the estimates."
         , if (cf$chains < 2L) "With one chain the Gelman-Rubin diagnostic is not available; use .chains = 2 or more."
         , if (length(weak)) sprintf("The HPD interval of lambda_%s reaches down to near zero: the term is weakly supported and its scores are uncertain.",
                                     paste(weak, collapse = ", lambda_"))
         , if (!isTRUE(x$settings$align)) "Signs were not aligned (.align = FALSE): posterior means of scores of weak terms shrink towards zero."),
       facts = f)
}

.bam_int_biplots <- function(x, ctx) {
  b  <- x$biplot_data
  lv <- 100 * attr(x, "meta")$level
  gen <- b[byplot == "Genotype"]
  env <- b[byplot == "Environment"]
  f <- list(stable_genotypes = gen[stable == TRUE, name], unstable_genotypes = gen[stable == FALSE, name],
            stable_environments = env[stable == TRUE, name], unstable_environments = env[stable == FALSE, name],
            positive_genotypes = gen[group == "positive", name], negative_genotypes = gen[group == "negative", name],
            positive_environments = env[group == "positive", name], negative_environments = env[group == "negative", name])
  list(title = "Biplot with highest posterior density regions",
       headline = sprintf("%d of %d genotypes and %d of %d environments have a %g%% HPD region that excludes the origin (significant interaction).",
                          length(f$unstable_genotypes), nrow(gen), length(f$unstable_environments), nrow(env), lv),
       findings = c(
           if (nrow(gen)) sprintf("Genotypes whose region contains the origin (interaction not different from zero, stable): %s.",
                                  .bam_list(f$stable_genotypes, ctx$top))
         , if (nrow(gen)) sprintf("Groups on the first term (marginal HPD interval excludes zero): positive %s; negative %s.",
                                  .bam_list(f$positive_genotypes, ctx$top), .bam_list(f$negative_genotypes, ctx$top))
         , if (nrow(env)) sprintf("Environments on the first term: positive %s; negative %s; the others do not differ from zero.",
                                  .bam_list(f$positive_environments, ctx$top), .bam_list(f$negative_environments, ctx$top))
         , "Overlapping regions mark genotypes (environments) whose interaction parameters do not differ: a homogeneous group."),
       recommendations = c(
           if (length(f$positive_genotypes) && length(f$positive_environments))
             sprintf("Genotypes %s are adapted to environments %s (both positive on the first term).",
                     .bam_list(f$positive_genotypes, ctx$top), .bam_list(f$positive_environments, ctx$top))
         , if (length(f$negative_genotypes) && length(f$negative_environments))
             sprintf("Genotypes %s are adapted to environments %s (both negative on the first term).",
                     .bam_list(f$negative_genotypes, ctx$top), .bam_list(f$negative_environments, ctx$top))),
       cautions = c("The regions are estimated from the draws by a kernel density; with few draws their outline is rough.",
                    "The second term is usually estimated with much more uncertainty than the first: long regions along the second axis."),
       facts = f)
}

.bam_int_ge_ammi <- function(x, ctx) {
  pct <- 100 * x$lambdas^2 / sum(x$lambdas^2)
  f <- list(lambdas = x$lambdas, percent = pct, sigma2 = x$sigma2)
  list(title = "Least-squares AMMI",
       headline = sprintf("The first two of %d interaction axes hold %.1f%% and %.1f%% of the interaction sum of squares.",
                          x$k, pct[1L], pct[2L]),
       findings = sprintf("Singular values: %s; error variance %s.", paste(.bam_n(utils::head(x$lambdas, 4L)), collapse = ", "),
                          .bam_n(x$sigma2)),
       recommendations = "Use bayes_ammi() for the posterior distributions and biplots() for inference on the biplot.",
       cautions = "Least-squares singular values are biased upwards (Crossa et al. 2011); the Bayesian estimates correct this.",
       facts = f)
}

.bam_int_predict <- function(x, ctx) {
  g  <- x$genotypes[order(if (ctx$hb) -mean else mean)]
  ms <- x$cells[n == 0L]
  f <- list(best = utils::head(as.character(g$Gen), ctx$top), missing = nrow(ms),
            widest = as.character(x$cells$Gen[which.max(x$cells$upper - x$cells$lower)]))
  list(title = "Predicted cell means",
       headline = sprintf("Highest posterior genotype means: %s; %d cell(s) were not observed and are predicted.",
                          .bam_list(f$best), f$missing),
       findings = c(sprintf("Genotype means (posterior mean, HPD): %s.",
                            paste(sprintf("%s %s (%s to %s)", g$Gen, .bam_n(g$mean), .bam_n(g$lower), .bam_n(g$upper)),
                                  collapse = "; ")),
                    if (nrow(ms)) sprintf("Predicted missing cells: %s.", paste(sprintf("%s in %s: %s (%s to %s)", ms$Gen,
                      ms$Env, .bam_n(ms$mean), .bam_n(ms$lower), .bam_n(ms$upper)), collapse = "; "))),
       recommendations = "Predictions of missing cells borrow the interaction pattern of the genotype and environment; confirm important ones in the field.",
       cautions = "Predictions of unobserved cells are wider than those of observed ones; do not read them as measurements.",
       facts = f)
}

.bam_int_winners <- function(x, ctx) {
  w <- x$winners
  sure <- w[Probability >= 0.5]
  f <- list(mega = x$mega, n_groups = nrow(x$mega), clear = sure$Env, unclear = w[Probability < 0.5, Env])
  list(title = "Which genotype wins where",
       headline = sprintf("%d genotype(s) are the most probable winners; %d of %d environments have a winner with probability of at least 0.5.",
                          nrow(x$mega), nrow(sure), nrow(w)),
       findings = c(sprintf("%s wins most probably in %s (mean probability %.2f).", x$mega$Winner, x$mega$Environments,
                            x$mega$mean_probability),
                    if (length(f$unclear)) sprintf("No clear winner (largest probability below 0.5) in: %s.", .bam_list(f$unclear, ctx$top))),
       recommendations = c(if (nrow(x$mega) > 1L) "Several winners: the environments may form mega-environments; confirm the grouping over years."
                           else "One genotype wins everywhere: one recommendation for the region is supported."),
       cautions = "A probability of winning near 1/G means the genotypes cannot be told apart in that environment.",
       facts = f)
}

# "a 0.89, b 0.63, ..." sorted by the probability; the genotypes below 0.005 are gathered at the end
.bam_prob_list <- function(g, p) {
  o <- order(-p)
  g <- g[o]
  p <- p[o]
  big <- p >= 0.005
  out <- paste(sprintf("%s %.2f", g[big], p[big]), collapse = ", ")
  if (any(!big)) out <- paste0(out, if (any(big)) "; " else "", "under 0.005: ", paste(g[!big], collapse = ", "))
  out
}

.bam_int_superiority <- function(x, ctx) {
  mg <- x$marginal
  ns <- attr(x, "meta")$selected
  f <- list(top_joint = utils::head(mg$Gen, ctx$top), perf = mg[order(-prob_performance)]$Gen[1L],
            stab = mg[order(-prob_stability)]$Gen[1L])
  list(title = "Probabilities of superior performance and stability",
       headline = sprintf("Selecting %d genotype(s): the highest joint probability of superior performance and stability is that of %s (%.2f).",
                          ns, mg$Gen[1L], mg$prob_joint[1L]),
       findings = c(sprintf("Probability of superior performance, largest first: %s.", .bam_prob_list(mg$Gen, mg$prob_performance)),
                    sprintf("Probability of superior stability, largest first: %s.", .bam_prob_list(mg$Gen, mg$prob_stability)),
                    sprintf("Most likely superior in performance: %s; in stability: %s.", f$perf, f$stab)),
       recommendations = {
         cand <- utils::head(mg[prob_joint >= 0.05]$Gen, 3L)
         if (length(cand)) sprintf("For a wide recommendation consider %s (joint probability of at least 0.05); the probabilities by environment (component conditional) show where each genotype is superior.",
                                   .bam_list(cand))
         else sprintf("No genotype is likely to be both superior and stable; choose on performance (%s) or stability (%s) according to the target.",
                      f$perf, f$stab)
       },
       cautions = "The joint probability multiplies the two probabilities as in Dias et al. (2022), which treats them as independent.",
       facts = f)
}

.bam_int_ranks <- function(x, ctx) {
  s <- x$summary
  f <- list(order = s$Gen, prob_first = stats::setNames(s$prob_first, s$Gen))
  list(title = "Posterior ranks",
       headline = sprintf("%s has the best expected rank (%.1f; probability %.2f of being first).", s$Gen[1L],
                          s$mean_rank[1L], s$prob_first[1L]),
       findings = sprintf("Expected ranks with HPD intervals: %s.",
                          paste(sprintf("%s %.1f (%g to %g)", s$Gen, s$mean_rank, s$lower, s$upper), collapse = "; ")),
       recommendations = "Genotypes whose rank intervals overlap cannot be ranked with confidence; select among them on other traits.",
       cautions = "Ranks depend on the set of genotypes in the trial.",
       facts = f)
}

.bam_int_stability <- function(x, ctx) {
  s <- x$stability
  f <- list(best_ysi = utils::head(s$Gen, ctx$top), most_stable = s$Gen[which.max(s$prob_most_stable)])
  list(title = "Bayesian stability",
       headline = sprintf("Best yield stability index (expected ranks of mean and ASV): %s; most probably the most stable: %s (%.2f).",
                          .bam_list(f$best_ysi), f$most_stable, max(s$prob_most_stable)),
       findings = sprintf("ASV with HPD intervals: %s.", paste(sprintf("%s %s (%s to %s)", s$Gen, .bam_n(s$ASV, 3L),
                                                                      .bam_n(s$ASV_lower, 3L), .bam_n(s$ASV_upper, 3L)), collapse = "; ")),
       recommendations = sprintf("Consider %s: they combine a good mean with a small AMMI stability value.", .bam_list(utils::head(f$best_ysi, 3L))),
       cautions = "The ASV uses the first two terms; with one term it is the absolute first score.",
       facts = f)
}

.bam_int_compare <- function(x, ctx) {
  t <- x$table
  sed   <- if (all(is.na(t$dWAIC_se))) t$se else t$dWAIC_se
  close <- t$model[t$dWAIC <= 2 * sed]
  f <- list(best = t$model[1L], close = close)
  list(title = "Model comparison",
       headline = sprintf("%s has the smallest WAIC%s.", t$model[1L],
                          if (length(close) > 1L) sprintf("; %s within two standard errors of the difference", .bam_list(setdiff(close, t$model[1L]))) else ""),
       findings = sprintf("%s: WAIC %s (%s), DIC %s, lower HPD limit of the last singular value %s.",
                          t$model, .bam_n(t$WAIC, 5L),
                          ifelse(t$dWAIC == 0, "the best", sprintf("%s above the best, standard error %s",
                                                                   .bam_n(t$dWAIC, 3L), .bam_n(sed, 3L))),
                          .bam_n(t$DIC, 5L), .bam_n(t$lambda_last_lower, 3L)),
       recommendations = "Among models within two standard errors of the best, prefer the simplest; a last singular value whose HPD interval reaches zero adds little.",
       cautions = "WAIC and DIC compare predictive accuracy for new plots of the same genotypes and environments.",
       facts = f)
}

.bam_int_ppc <- function(x, ctx) {
  c <- x$checks
  bad <- c[p_value < 0.025 | p_value > 0.975, statistic]
  f <- list(failed = bad)
  list(title = "Posterior predictive checks",
       headline = if (length(bad)) sprintf("The model does not reproduce: %s.", .bam_list(bad)) else
         "Replicated trials reproduce every checked feature of the observed trial.",
       findings = sprintf("%s: observed %s, replicated mean %s, p = %.3f.", c$statistic, .bam_n(c$observed), .bam_n(c$mean_rep), c$p_value),
       recommendations = c(if ("largest / smallest environment variance" %in% bad) "Fit environment-specific variances (.het = TRUE).",
                           if ("first singular value" %in% bad) "Check the number of multiplicative terms (bam_compare())."),
       cautions = "A p-value near 0 or 1 marks a misfit; values in between do not prove the model right.",
       facts = f)
}

.bam_int_bayes_stability <- function(x, ctx) {
  s <- x$stability
  st <- x$settings
  cf <- .bam_conv_facts(x)
  ok <- (is.na(cf$max_rhat) || cf$max_rhat < 1.1) && cf$min_ess_scalar >= 400
  conv <- list(converged = ok, max_rhat = cf$max_rhat,
               sentence = sprintf("%s (largest potential scale reduction factor %s; smallest effective sample size of the main parameters %s).",
                                  if (ok) "The chains have converged" else "The chains may not have converged",
                                  if (is.na(cf$max_rhat)) "not available with one chain" else format(round(cf$max_rhat, 3L)),
                                  format(round(cf$min_ess_scalar))))
  slopes <- "slope" %in% names(s)
  sf <- if (ctx$hb) s[order(-safety_first)] else s[order(safety_first)]
  f <- list(model = st$model, converged = conv$converged, max_rhat = conv$max_rhat,
            responsive = if (slopes) s[slope_lower > 1, Gen] else character(0),
            unresponsive = if (slopes) s[slope_upper < 1, Gen] else character(0),
            smallest_variance = s[which.min(dev_var), Gen], safety_first = utils::head(sf$Gen, ctx$top),
            best_safety = s[which.max(prob_best_safety), Gen])
  vname <- if (st$model == "Shukla") "stability variance" else "deviation variance"
  list(title = sprintf("Bayesian %s stability model",
                       c(ER = "Eberhart-Russell", FW = "Finlay-Wilkinson", Shukla = "Shukla")[[st$model]]),
       headline = sprintf("Best by the safety-first index (alpha = %g): %s; most probably the best: %s (%.2f).",
                          st$alpha, .bam_list(f$safety_first), f$best_safety, max(s$prob_best_safety)),
       findings = c(conv$sentence,
                    if (slopes) sprintf("Slopes whose HPD interval lies above 1 (respond more than average to good environments): %s; below 1 (respond less, relatively better in poor environments): %s.",
                                        .bam_list(f$responsive), .bam_list(f$unresponsive)),
                    sprintf("Smallest %s: %s (%s); largest: %s (%s).", vname, f$smallest_variance,
                            .bam_n(min(s$dev_var), 3L), s[which.max(dev_var), Gen], .bam_n(max(s$dev_var), 3L)),
                    if (slopes) sprintf("Share of each genotype's variation explained by the regression (R2): %s.",
                                        paste(sprintf("%s %.2f", s$Gen, s$R2), collapse = ", "))),
       recommendations = c(
         if (slopes) "For a wide recommendation choose genotypes with a high mean, a slope near 1 and a small deviation variance; for favourable environments a slope above 1 with a high mean."
         else "Prefer genotypes with a high mean and a small stability variance.",
         sprintf("The safety-first index combines the two: %s lead.", .bam_list(utils::head(f$safety_first, 3L)))),
       cautions = c(if (!isTRUE(f$converged)) "The chains have not converged by the criteria used; run longer chains before using the estimates.",
                    "The environmental index is estimated from the same trial; with few environments the slopes are uncertain (see the HPD intervals)."),
       facts = f)
}

.bam_int_indices_common <- function(x, ctx, title) {
  ii <- x$indices
  m <- x$means
  best <- ii[, .SD[which.max(prob_best)], by = index]
  f <- list(overall = utils::head(m$Gen, ctx$top), best = stats::setNames(as.character(best$Gen), best$index),
            agree = length(unique(best$Gen)) == 1L)
  list(title = title,
       headline = sprintf("Best on average over the indices (mean expected rank): %s.", .bam_list(f$overall)),
       findings = c(sprintf("Most probably the best on %s: %s (probability %.2f).", best$index, best$Gen,
                            best$prob_best),
                    if (f$agree) "All indices point to the same genotype." else
                      "The indices disagree on the best genotype: they weigh mean and stability differently."),
       recommendations = "Prefer genotypes that rank well on several indices and whose HPD intervals separate from the others; overlapping intervals mean the trial cannot tell them apart.",
       cautions = "Stability indices describe the environments of this trial; confirm over years before recommending.",
       facts = f)
}
.bam_int_indices <- function(x, ctx) .bam_int_indices_common(x, ctx, "Posterior stability indices")
.bam_int_ammi_indices <- function(x, ctx) .bam_int_indices_common(x, ctx, "Posterior AMMI-based stability indices")

.bam_guides <- list(
    bayes_ammi = "Bayesian AMMI: every parameter has a posterior distribution, summarised by its mean, standard deviation and highest posterior density (HPD) interval. Main effects whose HPD interval excludes zero differ from the average. Singular values measure the strength of each bilinear term. Genotype and environment scores with HPD intervals excluding zero interact significantly; same signs mean positive interaction. Check convergence first: potential scale reduction factor near 1, large effective sample sizes, Raftery-Lewis dependence factor below 5, Heidelberger-Welch tests passed.",
    biplots = "Biplot of lambda^(1/2) times the scores of the first two terms with HPD regions. A region containing the origin: no significant interaction (stable). Regions that do not overlap: genotypes (environments) that differ in interaction. The marginal HPD interval of the first coordinate decides significance on the first term, as in Crossa et al. (2011).",
    predict = "Posterior predictions of every cell mean, observed or not, with HPD intervals; genotype means over the environments.",
    winners = "Posterior probability that each genotype is best in each environment; the most probable winners group the environments.",
    superiority = "Dias et al. (2022): probability that a genotype is among the selected best (performance) and among the most stable (smallest interaction variance), and their product.",
    ranks = "Posterior distribution of the rank of every genotype: expected rank, HPD interval and probability of being first.",
    stability = "Posterior of the AMMI stability value and of the interaction variance of every genotype; the yield stability index adds the expected ranks of mean and ASV.",
    compare = "WAIC and DIC: smaller predicts new data better; a difference within two of its standard errors is not decisive.",
    ppc = "Posterior predictive p-values: the share of replicated trials with a statistic at least as large as the observed one; values near 0 or 1 mark a misfit.",
    bayes_stability = "Bayesian joint regression (Finlay-Wilkinson, Eberhart-Russell) or heterogeneous interaction variance (Shukla) model: the slope on the environmental index measures adaptability (1 = average response), the deviation (or stability) variance measures stability, and the safety-first index of Eskridge (1990) combines mean and stability. HPD intervals and probabilities replace significance tests.",
    indices = "Classical stability indices (Wricke, Shukla, Lin-Binns, Annicchiarico, Eskridge) computed on every posterior draw of the cell means: posterior means, HPD intervals, probability of being the best and expected rank.",
    ammi_indices = "AMMI-based stability indices (ASV, SIPC, EV, Za, WAAS) and the mean-stability index WAASY computed on every posterior draw of the singular values and scores.",
    ge_ammi = "Least-squares AMMI: singular value decomposition of the interaction effects. The share of each axis in the interaction sum of squares measures its importance; the estimates start the Bayesian sampler."
  )

utils::globalVariables(c("axis", "group"))
