// Gibbs sampler of the Bayesian joint regression (Finlay-Wilkinson, Eberhart-Russell) and heterogeneous interaction
// variance (Shukla) models:
//
//   y_ijl = mu + g_i + (1 + b_i) h_j + d_ij + e_ijl,   e_ijl ~ N(0, sigma^2),   d_ij ~ N(0, sigma^2_d,i)
//
// with sum(g) = sum(h) = sum(b) = 0 (orthonormal contrasts). model 0 = Eberhart-Russell (a deviation variance for
// every genotype), 1 = Finlay-Wilkinson (one deviation variance), 2 = Shukla (no slopes, b = 0; sigma^2_d,i is the
// stability variance of genotype i).
//
// Data enter as the sums S, sums of squares Q and counts N of the observed plots of every cell, so missing plots and
// empty cells need no imputation. The linear parameters are drawn with the deviations d integrated out (cell mean
// variance sigma^2_d,i + sigma^2 / n_ij), then the deviations given them, then the variances: a partially collapsed
// Gibbs sampler that mixes well when the deviations are large. Random numbers come from R in a fixed order, the same
// as the R reference .bam_stab_chain(), so set.seed() reproduces a run and the two give the same draws.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

static arma::mat contrasts_s(arma::uword n) {
  arma::mat m(n - 1, n, arma::fill::zeros);
  for (arma::uword r = 0; r < n - 1; ++r) {
    double nr = n - (r + 1.0);
    double c  = std::sqrt(nr * (nr + 1.0));
    m(r, r) = nr / c;
    for (arma::uword j = r + 1; j < n; ++j) m(r, j) = -1.0 / c;
  }
  return m;
}

// draw x ~ N(P^{-1} K s, P^{-1}) with P = K diag(w) K' + I / v, returned as K' x (a sum-zero vector)
static arma::vec draw_contrast(const arma::mat& K, const arma::vec& w, const arma::vec& s, double v) {
  const arma::uword m = K.n_rows;
  arma::mat P = K * arma::diagmat(w) * K.t();
  P.diag() += 1.0 / v;
  arma::vec rhs = K * s;
  arma::mat R = arma::chol(P);
  arma::vec mean = arma::solve(arma::trimatu(R), arma::solve(arma::trimatl(R.t()), rhs));
  arma::vec z(m);
  for (arma::uword k = 0; k < m; ++k) z(k) = norm_rand();
  arma::vec x = mean + arma::solve(arma::trimatu(R), z);
  return K.t() * x;
}

// [[Rcpp::export(name = ".bam_stab_cpp")]]
List bam_stab_cpp(const arma::mat& S, const arma::mat& Q, const arma::mat& N, int model, List pr, List start,
                  int n_total, IntegerVector keep, const arma::vec& y, const IntegerVector& cell, bool waic) {
  RNGScope scope;
  const arma::uword G = S.n_rows, E = S.n_cols, nobs = y.n_elem;
  const double mu_var = pr["mu_var"], g_var = pr["g_var"], h_var = pr["h_var"], b_var = pr["b_var"];
  const double nu = pr["nu"], s2 = pr["s2"], nu_d = pr["nu_d"], s2_d = pr["s2_d"];
  const arma::mat KG = contrasts_s(G), KE = contrasts_s(E);
  const bool slopes = model != 2, common = model == 1;
  // observed cell means and the within-cell sum of squares
  arma::mat Yb(G, E, arma::fill::zeros);
  double ssw = 0.0, ntot = 0.0;
  for (arma::uword j = 0; j < E; ++j) for (arma::uword i = 0; i < G; ++i) if (N(i, j) > 0) {
    Yb(i, j) = S(i, j) / N(i, j);
    ssw += Q(i, j) - S(i, j) * S(i, j) / N(i, j);
    ntot += N(i, j);
  }
  // state
  double mu = as<double>(start["mu"]), sig2 = as<double>(start["sigma2"]);
  arma::vec g = as<arma::vec>(start["g"]), h = as<arma::vec>(start["h"]), b = as<arma::vec>(start["b"]);
  arma::vec sd2 = as<arma::vec>(start["sd2"]);
  if (!slopes) b.zeros();
  arma::mat D(G, E, arma::fill::zeros);
  // storage
  int n_keep = 0;
  for (int it = 0; it < n_total; ++it) n_keep += keep[it];
  const arma::uword nsd = common ? 1 : G;
  arma::mat out_g(n_keep, G), out_h(n_keep, E), out_b(n_keep, G), out_sd(n_keep, nsd), out_d(n_keep, G * E);
  arma::vec out_mu(n_keep), out_s2(n_keep), dev(n_keep);
  arma::vec lse(nobs), lmean(nobs, arma::fill::zeros), lm2(nobs, arma::fill::zeros);
  lse.fill(-arma::datum::inf);
  int j_keep = 0;
  for (int it = 0; it < n_total; ++it) {
    // weights of the cell means with the deviations integrated out
    arma::mat W(G, E, arma::fill::zeros);
    for (arma::uword j = 0; j < E; ++j) for (arma::uword i = 0; i < G; ++i) if (N(i, j) > 0)
      W(i, j) = 1.0 / (sd2(common ? 0 : i) + sig2 / N(i, j));
    arma::vec c = 1.0 + b;
    // overall mean
    {
      arma::mat r = Yb - arma::repmat(g, 1, E) - c * h.t();
      double prec = arma::accu(W) + 1.0 / mu_var;
      mu = arma::accu(W % r) / prec + norm_rand() / std::sqrt(prec);
    }
    // genotype effects
    {
      arma::mat r = Yb - mu - c * h.t();
      g = draw_contrast(KG, arma::sum(W, 1), arma::sum(W % r, 1), g_var);
    }
    // environment effects (the environmental index)
    {
      arma::mat r = Yb - mu - arma::repmat(g, 1, E);
      arma::mat Wc = W.each_col() % (c % c);
      arma::mat WR = W % r;
      arma::mat Rc = WR.each_col() % c;
      h = draw_contrast(KE, arma::sum(Wc, 0).t(), arma::sum(Rc, 0).t(), h_var);
    }
    // slopes
    if (slopes) {
      arma::mat r = Yb - mu - arma::repmat(g, 1, E) - arma::repmat(h.t(), G, 1);
      arma::mat Wh = W.each_row() % (h % h).t();
      arma::mat WR = W % r;
      arma::mat Rh = WR.each_row() % h.t();
      b = draw_contrast(KG, arma::sum(Wh, 1), arma::sum(Rh, 1), b_var);
      c = 1.0 + b;
    }
    arma::mat eta = mu + arma::repmat(g, 1, E) + c * h.t();
    // deviations
    for (arma::uword j = 0; j < E; ++j) for (arma::uword i = 0; i < G; ++i) {
      double v = sd2(common ? 0 : i);
      if (N(i, j) > 0) {
        double a = N(i, j) / sig2, prec = a + 1.0 / v;
        D(i, j) = a * (Yb(i, j) - eta(i, j)) / prec + norm_rand() / std::sqrt(prec);
      } else {
        D(i, j) = norm_rand() * std::sqrt(v);
      }
    }
    arma::mat m = eta + D;
    // error variance
    {
      double ss = ssw;
      for (arma::uword j = 0; j < E; ++j) for (arma::uword i = 0; i < G; ++i) if (N(i, j) > 0)
        ss += N(i, j) * (Yb(i, j) - m(i, j)) * (Yb(i, j) - m(i, j));
      sig2 = (nu * s2 + ss) / R::rchisq(nu + ntot);
    }
    // deviation variances
    if (common) {
      sd2(0) = (nu_d * s2_d + arma::accu(D % D)) / R::rchisq(nu_d + G * E);
    } else {
      for (arma::uword i = 0; i < G; ++i)
        sd2(i) = (nu_d * s2_d + arma::accu(D.row(i) % D.row(i))) / R::rchisq(nu_d + E);
    }
    if (keep[it]) {
      out_mu(j_keep) = mu;
      out_s2(j_keep) = sig2;
      out_g.row(j_keep) = g.t();
      out_h.row(j_keep) = h.t();
      out_b.row(j_keep) = b.t();
      out_sd.row(j_keep) = sd2.t();
      out_d.row(j_keep) = arma::vectorise(D).t();
      double Dv = 0.0;
      for (arma::uword j = 0; j < E; ++j) for (arma::uword i = 0; i < G; ++i) if (N(i, j) > 0)
        Dv += N(i, j) * std::log(2.0 * M_PI * sig2) +
          (Q(i, j) - 2.0 * m(i, j) * S(i, j) + N(i, j) * m(i, j) * m(i, j)) / sig2;
      dev(j_keep) = Dv;
      if (waic) {
        const double cnt = j_keep + 1.0;
        for (arma::uword o = 0; o < nobs; ++o) {
          double r = y(o) - m(cell[o]);
          double l = -0.5 * std::log(2.0 * M_PI * sig2) - 0.5 * r * r / sig2;
          double mx = std::max(lse(o), l);
          lse(o) = mx + std::log(std::exp(lse(o) - mx) + std::exp(l - mx));
          double dl = l - lmean(o);
          lmean(o) += dl / cnt;
          lm2(o) += dl * (l - lmean(o));
        }
      }
      ++j_keep;
    }
  }
  return List::create(Named("mu") = out_mu, Named("sigma2") = out_s2, Named("g") = out_g, Named("h") = out_h,
                      Named("b") = out_b, Named("sd2") = out_sd, Named("d") = out_d, Named("deviance") = dev,
                      Named("lse") = lse, Named("lmean") = lmean, Named("lm2") = lm2);
}
