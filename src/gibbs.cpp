// Gibbs sampler of the Bayesian AMMI and GGE (SREG) models (Crossa et al. 2011; Perez-Elizalde et al. 2012).
//
// Data enter as the sums S, sums of squares Q and counts N of the observed plots of every genotype-environment cell.
// Cells with fewer than B plots (missing plots or whole missing cells) are completed by data augmentation: at every
// iteration the missing plots are drawn from N(eta_ij, sigma_j^2), after which the conditional distributions of the
// balanced design apply exactly.
//
// The random numbers come from R (norm_rand, unif_rand, rbeta, rchisq, exp_rand) so that set.seed() reproduces a run;
// in the balanced design with one error variance the order of the draws is that of the R reference sampler
// (.bam_chain in R/core-gibbs.R), so the two give identical draws.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// one draw from N(m, s^2) truncated to (lo, hi)
static double rtnorm1(double m, double s, double lo, double hi) {
  double a = (lo - m) / s, b = (hi - m) / s;
  if (a > 0) {
    double pa = R::pnorm(a, 0.0, 1.0, 0, 0), pb = R::pnorm(b, 0.0, 1.0, 0, 0);
    if (pa - pb <= 0) return lo + s * exp_rand() / a;
    double u = pb + (pa - pb) * unif_rand();
    return m + s * R::qnorm(u, 0.0, 1.0, 0, 0);
  }
  double pa = R::pnorm(a, 0.0, 1.0, 1, 0), pb = R::pnorm(b, 0.0, 1.0, 1, 0);
  if (pb - pa <= 0) return hi - s * exp_rand() / (-b);
  double u = pa + (pb - pa) * unif_rand();
  return m + s * R::qnorm(u, 0.0, 1.0, 1, 0);
}

// one draw from the von Mises-Fisher distribution on the unit sphere of R^p with parameter v (Wood 1994)
static arma::vec rvmf(const arma::vec& v) {
  const arma::uword p = v.n_elem;
  const double kap = arma::norm(v);
  arma::vec x(p);
  if (p == 1) {
    x(0) = (unif_rand() < R::plogis(2.0 * v(0), 0.0, 1.0, 1, 0)) ? 1.0 : -1.0;
    return x;
  }
  if (kap < 1e-12) {
    for (arma::uword i = 0; i < p; ++i) x(i) = norm_rand();
    return x / arma::norm(x);
  }
  arma::vec mu = v / kap;
  const double m1 = p - 1.0;
  const double b  = m1 / (2.0 * kap + std::sqrt(4.0 * kap * kap + m1 * m1));
  const double x0 = (1.0 - b) / (1.0 + b);
  const double cc = kap * x0 + m1 * std::log(1.0 - x0 * x0);
  double w;
  while (true) {
    double z = R::rbeta(m1 / 2.0, m1 / 2.0);
    w = (1.0 - (1.0 + b) * z) / (1.0 - (1.0 - b) * z);
    if (kap * w + m1 * std::log(1.0 - x0 * w) - cc >= std::log(unif_rand())) break;
  }
  arma::vec u(p);
  for (arma::uword i = 0; i < p; ++i) u(i) = norm_rand();
  u -= arma::dot(u, mu) * mu;
  u /= arma::norm(u);
  return w * mu + std::sqrt(std::max(0.0, 1.0 - w * w)) * u;
}

// One Gibbs sweep for x on the unit sphere of R^p with density proportional to exp(c'x + x'Ax) (the vector
// Bingham-von Mises-Fisher distribution; Hoff 2009): in the eigenbasis of A, y = V'x, every coordinate y_i is drawn
// given the direction u of the others, from the density exp(l_i t^2 + q (1 - t^2) + d_i t + s sqrt(1 - t^2))
// (1 - t^2)^((p - 3) / 2) on (-1, 1), by slice sampling with shrinkage (Neal 2003).
static arma::vec rbmf_sweep(const arma::vec& x, const arma::vec& c, const arma::mat& A) {
  const arma::uword p = x.n_elem;
  if (p == 1) {
    arma::vec out(1);
    out(0) = (unif_rand() < R::plogis(2.0 * c(0), 0.0, 1.0, 1, 0)) ? 1.0 : -1.0;
    return out;
  }
  arma::vec l;
  arma::mat V;
  arma::eig_sym(l, V, A);
  arma::vec y = V.t() * x, d = V.t() * c;
  const double jac = (p - 3.0) / 2.0;
  for (arma::uword i = 0; i < p; ++i) {
    double r2 = 1.0 - y(i) * y(i);
    if (r2 < 1e-14) continue;
    double r = std::sqrt(r2), q = 0.0, sd = 0.0;
    for (arma::uword j = 0; j < p; ++j) if (j != i) {
      double u = y(j) / r;
      q  += l(j) * u * u;
      sd += d(j) * u;
    }
    auto f = [&](double t) {
      double o = 1.0 - t * t;
      return l(i) * t * t + q * o + d(i) * t + sd * std::sqrt(o) + jac * std::log(o);
    };
    double t0 = y(i), lev = f(t0) - exp_rand(), lo = -1.0, hi = 1.0, t1 = t0;
    for (int tries = 0; tries < 200; ++tries) {
      t1 = lo + (hi - lo) * unif_rand();
      if (f(t1) > lev) break;
      if (t1 < t0) lo = t1; else hi = t1;
      t1 = t0;
    }
    double scale = std::sqrt(std::max(0.0, 1.0 - t1 * t1)) / r;
    for (arma::uword j = 0; j < p; ++j) if (j != i) y(j) *= scale;
    y(i) = t1;
  }
  y /= arma::norm(y);
  return V * y;
}

// orthonormal basis of the orthogonal complement of the columns of A: Gram-Schmidt (two passes) of A, then of the
// unit vectors in order (the same algorithm as .bam_null() in R)
static arma::mat null_basis(const arma::mat& A) {
  const arma::uword p = A.n_rows;
  arma::mat Q(p, 0);
  auto add = [&](arma::vec v, double tol) {
    for (int pass = 0; pass < 2; ++pass) if (Q.n_cols > 0) v -= Q * (Q.t() * v);
    double nv = arma::norm(v);
    if (nv > tol) { Q.insert_cols(Q.n_cols, v / nv); return true; }
    return false;
  };
  for (arma::uword c = 0; c < A.n_cols; ++c) add(A.col(c), 1e-10);
  const arma::uword r = Q.n_cols;
  for (arma::uword i = 0; i < p && Q.n_cols < p; ++i) {
    arma::vec e(p, arma::fill::zeros);
    e(i) = 1.0;
    add(e, 1e-6);
  }
  return Q.cols(r, p - 1);
}

// the columns of M other than k, with a column of ones in front when `one` is true
static arma::mat others(const arma::mat& M, arma::uword k, bool one) {
  arma::mat out(M.n_rows, 0);
  if (one) out.insert_cols(0, arma::ones<arma::vec>(M.n_rows));
  for (arma::uword j = 0; j < M.n_cols; ++j) if (j != k) out.insert_cols(out.n_cols, M.col(j));
  return out;
}

// orthonormal contrasts K_n ((n - 1) x n), as matrix_k() in R
static arma::mat contrasts(arma::uword n) {
  arma::mat m(n - 1, n, arma::fill::zeros);
  for (arma::uword r = 0; r < n - 1; ++r) {
    double nr = n - (r + 1.0);
    double c  = std::sqrt(nr * (nr + 1.0));
    m(r, r) = nr / c;
    for (arma::uword j = r + 1; j < n; ++j) m(r, j) = -1.0 / c;
  }
  return m;
}

// one draw from the von Mises distribution with mean direction m and concentration kap (Best and Fisher 1979)
static double rvonmises(double m, double kap) {
  if (kap < 1e-8) return M_PI * (2.0 * unif_rand() - 1.0);
  const double t = 1.0 + std::sqrt(1.0 + 4.0 * kap * kap);
  const double rho = (t - std::sqrt(2.0 * t)) / (2.0 * kap);
  const double r = (1.0 + rho * rho) / (2.0 * rho);
  double f;
  while (true) {
    double u1 = unif_rand(), u2 = unif_rand();
    double z = std::cos(M_PI * u1);
    f = (1.0 + r * z) / (r + z);
    double c = kap * (r - f);
    if (c * (2.0 - c) - u2 > 0 || std::log(c / u2) + 1.0 - c >= 0) break;
  }
  double u3 = unif_rand();
  return (u3 > 0.5 ? 1.0 : -1.0) * std::acos(f) + m;
}

// Two columns at a time (Hoff 2009): when the columns fill the space allowed by the constraints, each column is fixed
// by the others up to its sign and a column-wise sampler cannot move. Given the other columns, every pair (x_k, x_l)
// is an orthonormal pair of the 2-dimensional complement N, x_k = N (cos t, sin t)', x_l = e N (-sin t, cos t)',
// e = +1 or -1; with vMF parameters c_k, c_l its density is exp(k_e cos(t - phi_e)): e is drawn with probability
// proportional to I_0(k_e), then t from the von Mises distribution.
static void pair_sweep(arma::mat& X, const arma::mat& C, bool one) {
  const arma::uword K = X.n_cols;
  for (arma::uword k = 0; k + 1 < K; ++k) for (arma::uword l = k + 1; l < K; ++l) {
    arma::mat O;
    if (one) O = arma::ones(X.n_rows, 1);
    for (arma::uword m = 0; m < K; ++m) if (m != k && m != l) O = arma::join_rows(O, X.col(m));
    arma::mat N = null_basis(O);
    if (N.n_cols != 2) continue;
    arma::vec a = N.t() * C.col(k), b = N.t() * C.col(l);
    double ap = a(0) + b(1), bp = a(1) - b(0);
    double am = a(0) - b(1), bm = a(1) + b(0);
    double kp = std::sqrt(ap * ap + bp * bp), km = std::sqrt(am * am + bm * bm);
    double lp = std::log(R::bessel_i(kp, 0.0, 2.0)) + kp, lm = std::log(R::bessel_i(km, 0.0, 2.0)) + km;
    double pp = 1.0 / (1.0 + std::exp(lm - lp));
    double e = (unif_rand() < pp) ? 1.0 : -1.0;
    double th = (e > 0) ? rvonmises(std::atan2(bp, ap), kp) : rvonmises(std::atan2(bm, am), km);
    X.col(k)     = N.col(0) * std::cos(th) + N.col(1) * std::sin(th);
    X.col(l)     = e * (-N.col(0) * std::sin(th) + N.col(1) * std::cos(th));
  }
}

// [[Rcpp::export(name = ".bam_gibbs_cpp")]]
List bam_gibbs_cpp(const arma::mat& S, const arma::mat& Q, const arma::mat& N, int B, int K, bool gge, bool het,
                   List pr, List start, const arma::mat& ref_alpha, bool align, int n_total,
                   const IntegerVector& keep, const arma::vec& y, const arma::uvec& cell, bool waic) {
  const arma::uword G = S.n_rows, E = S.n_cols;
  const double n = (double) B * G * E;
  const bool balanced = arma::all(arma::vectorise(N) == (double) B);
  const arma::mat KG = contrasts(G), KE = contrasts(E);
  // priors
  const double mu0 = pr["mu_mean"], vmu = pr["mu_var"], vtau = pr["tau_var"], vdel = pr["delta_var"];
  const double l0 = pr["lambda_mean"], vl = pr["lambda_var"], nu = pr["nu"];
  const arma::vec mtau = KG * as<arma::vec>(pr["tau_mean"]);
  const arma::vec mdel = KE * as<arma::vec>(pr["delta_mean"]);
  const arma::vec s2p  = as<arma::vec>(pr["s2"]);
  // state
  arma::vec sig2 = as<arma::vec>(start["sigma2"]);
  arma::vec lam  = as<arma::vec>(start["lambda"]);
  arma::mat A    = as<arma::mat>(start["alpha"]);
  arma::mat Gm   = as<arma::mat>(start["gamma"]);
  double mu = 0.0;
  arma::vec tau(G, arma::fill::zeros), del(E, arma::fill::zeros);
  if (start.containsElementNamed("mu")) {
    mu  = as<double>(start["mu"]);
    tau = as<arma::vec>(start["tau"]);
    del = as<arma::vec>(start["delta"]);
  }
  // completed data
  arma::mat Y = S / (double) B;
  arma::vec SSW(E, arma::fill::zeros);
  for (arma::uword j = 0; j < E; ++j) for (arma::uword i = 0; i < G; ++i)
    if (N(i, j) > 0) SSW(j) += Q(i, j) - S(i, j) * S(i, j) / N(i, j);
  // storage
  const int nsig = het ? E : 1;
  const int ntau = gge ? 0 : G;
  const int nk = sum(keep);
  arma::mat sc(nk, 1 + nsig + ntau + E + K);
  arma::cube aa(nk, G, K), gg(nk, E, K);
  arma::vec dev(nk);
  const arma::uword nobs = y.n_elem;
  arma::vec lse(nobs), lmean(nobs, arma::fill::zeros), lm2(nobs, arma::fill::zeros);
  lse.fill(-arma::datum::inf);
  int j_keep = 0;
  for (int it = 0; it < n_total; ++it) {
    // data augmentation of the missing plots
    if (!balanced) {
      arma::mat eta = mu + arma::repmat(tau, 1, E) + arma::repmat(del.t(), G, 1) + A * arma::diagmat(lam) * Gm.t();
      SSW.zeros();
      for (arma::uword j = 0; j < E; ++j) {
        double sd = std::sqrt(het ? sig2(j) : sig2(0));
        for (arma::uword i = 0; i < G; ++i) {
          double tot = S(i, j), sq = Q(i, j);
          for (int m = (int) N(i, j); m < B; ++m) {
            double x = eta(i, j) + sd * norm_rand();
            tot += x;
            sq  += x * x;
          }
          Y(i, j) = tot / B;
          SSW(j) += sq - tot * tot / B;
        }
      }
    }
    if (!het) {
      // ---------------------------------------------------------------------------------------------------------------
      // one error variance: the order of the R reference sampler
      // ---------------------------------------------------------------------------------------------------------------
      const double s2 = sig2(0);
      const double c_mu = s2 / vmu;
      mu = (n * arma::mean(arma::vectorise(Y)) + c_mu * mu0) / (n + c_mu) + std::sqrt(s2 / (n + c_mu)) * norm_rand();
      if (!gge) {
        const double c_t = s2 / vtau;
        arma::vec yr = KG * arma::mean(Y, 1), t(G - 1);
        for (arma::uword i = 0; i < G - 1; ++i)
          t(i) = (B * E * yr(i) + c_t * mtau(i)) / (B * E + c_t) + std::sqrt(s2 / (B * E + c_t)) * norm_rand();
        tau = KG.t() * t;
      }
      const double c_d = s2 / vdel;
      arma::vec yc = KE * arma::mean(Y, 0).t(), d(E - 1);
      for (arma::uword j = 0; j < E - 1; ++j)
        d(j) = (B * G * yc(j) + c_d * mdel(j)) / (B * G + c_d) + std::sqrt(s2 / (B * G + c_d)) * norm_rand();
      del = KE.t() * d;
      const double c_l = s2 / vl;
      for (int k = 0; k < K; ++k) {
        arma::vec yg = Y * Gm.col(k);
        double lo = (k < K - 1) ? lam(k + 1) : 0.0;
        double hi = (k > 0) ? lam(k - 1) : R_PosInf;
        lam(k) = rtnorm1((B * arma::dot(A.col(k), yg) + c_l * l0) / (B + c_l), std::sqrt(s2 / (B + c_l)), lo, hi);
        arma::mat H = null_basis(others(A, k, true));
        A.col(k) = H * rvmf(H.t() * yg * (B * lam(k) / s2));
        H = null_basis(others(Gm, k, !gge));
        Gm.col(k) = H * rvmf(H.t() * (Y.t() * A.col(k)) * (B * lam(k) / s2));
        if (align && arma::dot(A.col(k), ref_alpha.col(k)) < 0) {
          A.col(k)  = -A.col(k);
          Gm.col(k) = -Gm.col(k);
        }
      }
      // pair moves for score matrices whose columns fill the constrained space (one free dimension per column)
      if (K > 1) {
        if (G - 1 == (arma::uword) K) {
          arma::mat C = (Y * Gm * arma::diagmat(lam * (double) B)) / s2;
          pair_sweep(A, C, true);
        }
        if (E - (gge ? 0 : 1) == (arma::uword) K) {
          arma::mat C = (Y.t() * A * arma::diagmat(lam * (double) B)) / s2;
          pair_sweep(Gm, C, !gge);
        }
        if (align) for (int k = 0; k < K; ++k) if (arma::dot(A.col(k), ref_alpha.col(k)) < 0) {
          A.col(k)  = -A.col(k);
          Gm.col(k) = -Gm.col(k);
        }
      }
      arma::mat res = Y - mu - arma::repmat(tau, 1, E) - arma::repmat(del.t(), G, 1) - A * arma::diagmat(lam) * Gm.t();
      sig2(0) = (nu * s2p(0) + arma::accu(SSW) + B * arma::accu(arma::square(res))) / R::rchisq(nu + n);
    } else {
      // ---------------------------------------------------------------------------------------------------------------
      // an error variance for every environment: weights w_j = 1 / sigma_j^2
      // ---------------------------------------------------------------------------------------------------------------
      arma::vec w = 1.0 / sig2;
      const double W = arma::accu(w);
      arma::mat Int = A * arma::diagmat(lam) * Gm.t();
      arma::mat Z = Y - Int;
      if (!gge) {
        // genotype effects: data precision B W for every genotype (independent of mu and delta, as tau sums to zero)
        arma::vec zt = (Z * w) / W, t(G - 1);
        arma::vec kz = KG * zt;
        const double prec = B * W + 1.0 / vtau;
        for (arma::uword i = 0; i < G - 1; ++i) t(i) = (B * W * kz(i) + mtau(i) / vtau) / prec + norm_rand() / std::sqrt(prec);
        tau = KG.t() * t;
      }
      // overall mean and environment effects jointly: eta_j = mu + delta_j, data precision B G w_j
      arma::vec cm = arma::mean(Z.each_col() - tau, 0).t();
      arma::mat M(E, E);
      M.col(0).ones();
      M.cols(1, E - 1) = KE.t();
      arma::vec P = B * G * w;
      arma::mat prec = M.t() * arma::diagmat(P) * M;
      prec(0, 0) += 1.0 / vmu;
      for (arma::uword j = 1; j < E; ++j) prec(j, j) += 1.0 / vdel;
      arma::vec pm(E);
      pm(0) = mu0 / vmu;
      pm.subvec(1, E - 1) = mdel / vdel;
      arma::vec rhs = M.t() * (P % cm) + pm;
      arma::mat U = arma::chol(prec);
      arma::vec mean = arma::solve(arma::trimatu(U), arma::solve(arma::trimatl(U.t()), rhs));
      arma::vec z(E);
      for (arma::uword j = 0; j < E; ++j) z(j) = norm_rand();
      arma::vec phi = mean + arma::solve(arma::trimatu(U), z);
      mu  = phi(0);
      del = KE.t() * phi.subvec(1, E - 1);
      arma::mat F = mu + arma::repmat(tau, 1, E) + arma::repmat(del.t(), G, 1);
      for (int k = 0; k < K; ++k) {
        arma::mat Rk = Y - F - A * arma::diagmat(lam) * Gm.t() + lam(k) * A.col(k) * Gm.col(k).t();
        arma::vec gw = Gm.col(k) % w;
        double pl = B * arma::dot(Gm.col(k), gw) + 1.0 / vl;
        double ll = B * arma::dot(A.col(k), Rk * gw) + l0 / vl;
        double lo = (k < K - 1) ? lam(k + 1) : 0.0;
        double hi = (k > 0) ? lam(k - 1) : R_PosInf;
        lam(k) = rtnorm1(ll / pl, 1.0 / std::sqrt(pl), lo, hi);
        // genotype scores: von Mises-Fisher (the quadratic term is constant on the sphere)
        arma::mat H = null_basis(others(A, k, true));
        A.col(k) = H * rvmf(H.t() * (Rk * gw) * (B * lam(k)));
        // environment scores: Bingham-von Mises-Fisher (the quadratic term differs between environments)
        H = null_basis(others(Gm, k, !gge));
        arma::vec lin = (Rk.t() * A.col(k)) % w * (B * lam(k));
        arma::mat Aq = H.t() * arma::diagmat(-0.5 * B * lam(k) * lam(k) * w) * H;
        arma::vec xg = H.t() * Gm.col(k);
        xg /= arma::norm(xg);
        Gm.col(k) = H * rbmf_sweep(xg, H.t() * lin, Aq);
        if (align && arma::dot(A.col(k), ref_alpha.col(k)) < 0) {
          A.col(k)  = -A.col(k);
          Gm.col(k) = -Gm.col(k);
        }
      }
      arma::mat res = Y - F - A * arma::diagmat(lam) * Gm.t();
      for (arma::uword j = 0; j < E; ++j) {
        double sj = (s2p.n_elem == E) ? s2p(j) : s2p(0);
        sig2(j) = (nu * sj + SSW(j) + B * arma::accu(arma::square(res.col(j)))) / R::rchisq(nu + B * G);
      }
    }
    if (keep[it]) {
      arma::rowvec row(1 + nsig + ntau + E + K);
      int c = 0;
      row(c++) = mu;
      for (int j = 0; j < nsig; ++j) row(c++) = sig2(j);
      for (int i = 0; i < ntau; ++i) row(c++) = tau(i);
      for (arma::uword j = 0; j < E; ++j) row(c++) = del(j);
      for (int k = 0; k < K; ++k) row(c++) = lam(k);
      sc.row(j_keep) = row;
      for (int k = 0; k < K; ++k) {
        aa.slice(k).row(j_keep) = A.col(k).t();
        gg.slice(k).row(j_keep) = Gm.col(k).t();
      }
      // deviance of the observed plots and the pointwise log-likelihood (WAIC)
      arma::mat eta = mu + arma::repmat(tau, 1, E) + arma::repmat(del.t(), G, 1) + A * arma::diagmat(lam) * Gm.t();
      double D = 0.0;
      for (arma::uword j = 0; j < E; ++j) {
        double s = het ? sig2(j) : sig2(0);
        for (arma::uword i = 0; i < G; ++i) {
          double nij = N(i, j);
          if (nij == 0) continue;
          D += nij * std::log(2.0 * M_PI * s) + (Q(i, j) - 2.0 * eta(i, j) * S(i, j) + nij * eta(i, j) * eta(i, j)) / s;
        }
      }
      dev(j_keep) = D;
      if (waic) {
        const double cnt = j_keep + 1.0;
        for (arma::uword o = 0; o < nobs; ++o) {
          arma::uword ci = cell(o), jj = ci / G;
          double s = het ? sig2(jj) : sig2(0);
          double r = y(o) - eta(ci);
          double l = -0.5 * std::log(2.0 * M_PI * s) - 0.5 * r * r / s;
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
  return List::create(Named("scalar") = sc, Named("alpha") = aa, Named("gamma") = gg, Named("deviance") = dev,
                      Named("lse") = lse, Named("lmean") = lmean, Named("lm2") = lm2,
                      Named("accept") = NA_REAL);
}
