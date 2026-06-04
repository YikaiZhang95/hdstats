// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;
using namespace std;

// --------------------------- small utilities ---------------------------

// Branch-free clamp: compiles to SIMD min/max instructions
inline double clip_huber(double r, double h) {
  return r > h ? h : (r < -h ? -h : r);
}

// return a with the sign of b (Fortran SIGN(a,b))
inline double sign_copy(double a, double b) {
  return (b >= 0.0) ? std::abs(a) : -std::abs(a);
}

// elementwise max with 0 for matrix pf (column-major traversal)
static void clamp_nonneg(NumericMatrix pf) {
  const int p = pf.nrow(), q = pf.ncol();
  for (int l = 0; l < q; ++l)
    for (int j = 0; j < p; ++j)
      if (pf(j, l) < 0.0) pf(j, l) = 0.0;
}

// elementwise max with 0 for vector pf2
static void clamp_nonneg_vec(NumericVector v) {
  for (int i = 0; i < v.size(); ++i)
    if (v[i] < 0.0) v[i] = 0.0;
}

// gradient helper used for strong rule / KKT checks
static NumericVector huber_drv_internal(const NumericMatrix& X,
                                        const NumericVector& r,
                                        const double hval) {
  const int n = X.nrow(), p = X.ncol();
  const double ninv = 1.0 / static_cast<double>(n);
  NumericVector vl(p);
  std::vector<double> dl(n);
  for (int i = 0; i < n; ++i) dl[i] = clip_huber(r[i], hval);
  for (int j = 0; j < p; ++j) {
    double s = 0.0;
    const double* colj = &X(0, j);
    for (int i = 0; i < n; ++i) s += colj[i] * dl[i];
    vl[j] = -s * ninv;
  }
  return vl;
}

// --------------------------- chkvars ---------------------------
static IntegerVector chkvars_cpp(const NumericMatrix& X) {
  const int n = X.nrow(), p = X.ncol();
  IntegerVector ju(p);
  for (int j = 0; j < p; ++j) {
    const double* col = &X(0, j);
    double mu = 0.0;
    for (int i = 0; i < n; ++i) mu += col[i];
    mu /= static_cast<double>(n);
    double ss = 0.0;
    for (int i = 0; i < n; ++i) {
      const double z = col[i] - mu;
      ss += z * z;
    }
    ju[j] = (ss > 0.0) ? 1 : 0;
  }
  return ju;
}

// --------------------------- Standard ---------------------------
static void standardize_cpp(NumericMatrix& X,
                            const IntegerVector& ju,
                            const int isd,
                            NumericVector& xmean,
                            NumericVector& xnorm,
                            NumericVector& maj) {
  const int n = X.nrow(), p = X.ncol();
  xmean = NumericVector(p);
  xnorm = NumericVector(p);
  maj   = NumericVector(p);

  for (int j = 0; j < p; ++j) {
    if (ju[j] == 0) {
      xmean[j] = 0.0; xnorm[j] = 1.0; maj[j] = 0.0;
      continue;
    }
    double mu = 0.0;
    for (int i = 0; i < n; ++i) mu += X(i, j);
    mu /= static_cast<double>(n);
    xmean[j] = mu;
    for (int i = 0; i < n; ++i) X(i, j) -= mu;
    double ss = 0.0;
    for (int i = 0; i < n; ++i) ss += X(i, j) * X(i, j);
    const double msq = ss / static_cast<double>(n);
    const double sd  = (msq > 0.0) ? std::sqrt(msq) : 0.0;
    xnorm[j] = (sd > 0.0) ? sd : 1.0;
    if (isd == 1 && sd > 0.0) {
      const double invsd = 1.0 / sd;
      for (int i = 0; i < n; ++i) X(i, j) *= invsd;
      maj[j] = 1.0;
    } else {
      maj[j] = msq;
    }
  }
}

// --------------------------- CD update helpers ---------------------------

// Soft-threshold update for one variable.
// Uses pre-maintained dl[] for the gradient; updates r[] and dl[] in-place.
// No loop-carried dependency → both passes are SIMD-vectorizable.
// Returns d^2 (0 if unchanged).
inline double update_coord(const double* xk, int n,
                           double& bk,
                           double* r_ptr, double* dl_ptr,
                           double majk, double hval,
                           double al_pfk, double lam2_pf2k, double ninv) {
  // Pass 1: dot product using pre-clipped dl (vectorizable reduction)
  double u = 0.0;
  for (int i = 0; i < n; ++i) u -= dl_ptr[i] * xk[i];
  u = majk * bk - u * ninv;

  const double v = std::abs(u) - al_pfk;
  const double new_bk = (v > 0.0) ? sign_copy(v, u) / (majk + lam2_pf2k) : 0.0;
  const double d = new_bk - bk;
  bk = new_bk;
  if (d == 0.0) return 0.0;

  // Pass 2: SAXPY on r and clip-update dl — NO loop-carried dependency → vectorizable
  for (int i = 0; i < n; ++i) {
    const double ri = r_ptr[i] - xk[i] * d;
    r_ptr[i] = ri;
    dl_ptr[i] = clip_huber(ri, hval);
  }
  return d * d;
}

// Intercept update.
// Recomputes sum_dl from scratch (vectorizable reduction), then updates r and dl.
// Returns d^2 (0 if no change).
inline double update_intercept(double& b0, double* r_ptr, double* dl_ptr, int n,
                                double hval, double ninv, double mval) {
  // Vectorizable reduction over pre-maintained dl
  double sum_dl = 0.0;
  for (int i = 0; i < n; ++i) sum_dl += dl_ptr[i];

  const double d = sum_dl * ninv / mval;
  if (d == 0.0) return 0.0;
  b0 += d;

  // SAXPY + clip-update — vectorizable
  for (int i = 0; i < n; ++i) {
    const double ri = r_ptr[i] - d;
    r_ptr[i] = ri;
    dl_ptr[i] = clip_huber(ri, hval);
  }
  return d * d;
}

// --------------------------- huber_path core ---------------------------

static List huber_path_core(double alpha,
                            double lam2_in,
                            double hval,
                            NumericVector maj,
                            double mval,
                            const NumericMatrix& X,
                            const NumericVector& y,
                            const IntegerVector& ju,
                            int pfncol,
                            const NumericMatrix& pf,
                            const NumericVector& pf2,
                            int dfmax,
                            int pmax,
                            int nlam,
                            double flmin_in,
                            const NumericVector& ulam,
                            double eps,
                            int maxit,
                            int istrong) {

  const double BIG   = 9.9e30;
  const double MFL   = 1.0e-6;
  const int    MNLAM = 6;

  const int n = X.nrow();
  const int p = X.ncol();
  const double ninv = 1.0 / static_cast<double>(n);

  NumericVector r = clone(y);
  double* r_ptr = r.begin();
  NumericVector ga(p), vl(p);
  NumericVector b0_out(nlam);
  NumericMatrix beta_out(pmax, nlam);
  IntegerVector m(pmax);
  IntegerVector mm(p);
  IntegerVector nbeta_out(nlam);
  NumericVector alam(nlam);
  IntegerVector jxx(p);

  double al = 0.0, al0 = 0.0;
  double alf = 0.01;
  int npass = 0, nalam = 0, jerr = 0;

  const double* Xbase = X.begin();
  auto colptr = [&](int j)->const double* { return Xbase + (size_t)j*n; };

  for (int j = 0; j < p; ++j) maj[j] *= mval;

  std::vector<double> dl(n);
  double* dl_ptr = dl.data();
  for (int i = 0; i < n; ++i) dl_ptr[i] = clip_huber(r_ptr[i], hval);

  if (istrong == 1) {
    std::fill(jxx.begin(), jxx.end(), 0);
  } else {
    std::fill(jxx.begin(), jxx.end(), 1);
  }

  int mnl = std::min(MNLAM, nlam);
  double flmin = flmin_in;
  if (flmin < 1.0) {
    flmin = std::max(MFL, flmin);
    if (nlam > 1) alf = std::pow(flmin, 1.0 / (static_cast<double>(nlam) - 1.0));
  }

  std::vector<double> b(p, 0.0), oldb_all(p, 0.0);
  // Precomputed per-variable penalty terms (updated once per lambda)
  std::vector<double> al_pf(p), lam2_pf2(p);
  double b0 = 0.0, oldb0 = 0.0;
  int ni = 0;

  for (int l = 0; l < nlam; ++l) {
    const int pfl = (pfncol == 1 ? 0 : l);

    al0 = al;
    if (flmin >= 1.0) {
      al = ulam[l];
    } else {
      if (l > 1) {
        al = al * alf;
      } else if (l == 0) {
        al = BIG;
      } else {  // l == 1
        al0 = 0.0;
        vl  = huber_drv_internal(X, r, hval);
        for (int j = 0; j < p; ++j) ga[j] = std::abs(vl[j]);
        double maxv = 0.0;
        for (int j = 0; j < p; ++j) {
          if (ju[j] == 0) continue;
          const double pf_j = pf(j, pfl);
          if (pf_j > 0.0) maxv = std::max(maxv, ga[j] / pf_j);
        }
        al = maxv * alf;
      }
    }

    double lam2 = lam2_in;
    if (alpha != -1.0) lam2 = al * (1.0 - alpha) * 0.5 / alpha;

    // Precompute per-variable penalty terms once per lambda — avoids
    // repeated pf(k,pfl) matrix access and al/lam2 multiplications in hot loops
    for (int k = 0; k < p; ++k) {
      al_pf[k]    = al * pf(k, pfl);
      lam2_pf2[k] = pf2[k] * lam2;
    }

    const double tlam = 2.0 * al - al0;
    for (int j = 0; j < p; ++j) {
      if (jxx[j] == 1) continue;
      if (ga[j] > pf(j, pfl) * tlam) jxx[j] = 1;
    }

    // ---------------- outer loop ----------------
    for (;;) {
      oldb0 = b0;
      for (int t = 0; t < ni; ++t) oldb_all[m[t] - 1] = b[m[t] - 1];

      // ------------- middle loop: all strong-set variables -------------
      for (;;) {
        npass += 1;
        double dif = 0.0;

        for (int k = 0; k < p; ++k) {
          if (jxx[k] == 0 || ju[k] == 0) continue;
          const double sq_d = update_coord(colptr(k), n, b[k],
                                           r_ptr, dl_ptr,
                                           maj[k], hval,
                                           al_pf[k], lam2_pf2[k], ninv);
          dif = std::max(dif, sq_d);
          if (sq_d != 0.0 && mm[k] == 0) {
            ni += 1;
            if (ni > pmax) { jerr = -10000 - (l + 1); goto done_path; }
            mm[k] = ni;
            m[ni - 1] = k + 1;
          }
        }

        dif = std::max(dif, update_intercept(b0, r_ptr, dl_ptr, n, hval, ninv, mval));

        if (dif < eps) break;
        if (npass > maxit) { jerr = -1; goto done_path; }

        // ------------- inner loop: active set only -------------
        for (;;) {
          npass += 1;
          double dif = 0.0;

          for (int t = 0; t < ni; ++t) {
            const int k = m[t] - 1;
            dif = std::max(dif, update_coord(colptr(k), n, b[k],
                                             r_ptr, dl_ptr,
                                             maj[k], hval,
                                             al_pf[k], lam2_pf2[k], ninv));
          }

          dif = std::max(dif, update_intercept(b0, r_ptr, dl_ptr, n, hval, ninv, mval));

          if (dif < eps) break;
          if (npass > maxit) { jerr = -1; goto done_path; }
        } // inner
      } // middle

      // ----- outer convergence check -----
      {
        int jx = 0;
        if ((b0 - oldb0) * (b0 - oldb0) >= eps) jx = 1;
        if (!jx) {
          for (int t = 0; t < ni; ++t) {
            const int k = m[t] - 1;
            const double diff = b[k] - oldb_all[k];
            if (diff * diff >= eps) { jx = 1; break; }
          }
        }
        if (jx) continue;

        // KKT check for discarded variables
        vl = huber_drv_internal(X, r, hval);
        for (int j = 0; j < p; ++j) ga[j] = std::abs(vl[j]);
        for (int j = 0; j < p; ++j) {
          if (jxx[j] == 1) continue;
          if (ga[j] > al_pf[j]) { jxx[j] = 1; jx = 1; }
        }
        if (jx) continue;
      }

      break; // outer converged
    } // outer

    // ---- save results for this lambda ----
    if (ni > pmax) { jerr = -10000 - (l + 1); goto done_path; }
    for (int t = 0; t < ni; ++t) beta_out(t, l) = b[m[t] - 1];
    nbeta_out[l] = ni;
    b0_out[l]    = b0;
    alam[l]      = al;
    nalam        = l + 1;
    if (l >= (mnl - 1) && flmin < 1.0) {
      int me = 0;
      for (int t = 0; t < ni; ++t) if (beta_out(t, l) != 0.0) ++me;
      if (me > dfmax) break;
    }
  } // for each lambda

done_path:
  return List::create(
    _["nalam"] = nalam,
    _["b0"]    = b0_out,
    _["beta"]  = beta_out,
    _["ibeta"] = clone(m),
    _["nbeta"] = nbeta_out,
    _["alam"]  = alam,
    _["npass"] = npass,
    _["jerr"]  = jerr
  );
}

// --------------------------- Exported: huber_path_cpp ---------------------------

// [[Rcpp::export]]
List huber_path_cpp(double alpha,
                    double lam2,
                    double hval,
                    NumericVector maj,
                    double mval,
                    NumericMatrix X,
                    NumericVector y,
                    IntegerVector ju,
                    int pfncol,
                    NumericMatrix pf,
                    NumericVector pf2,
                    int dfmax,
                    int pmax,
                    int nlam,
                    double flmin,
                    NumericVector ulam,
                    double eps,
                    int maxit,
                    int istrong) {
  const int n = X.nrow(), p = X.ncol();
  if (y.size() != n) stop("Length of y must equal nrow(X).");
  if (maj.size() != p) stop("Length of maj must equal ncol(X).");
  if (pf.nrow() != p) stop("pf must have nrow = ncol(X).");
  if (!(pfncol == 1 || pf.ncol() == nlam))
    stop("pfncol must be 1 or equal to nlam and match pf.ncol().");
  if (pf2.size() != p) stop("pf2 must have length ncol(X).");
  if (ulam.size() != nlam && flmin >= 1.0)
    stop("When flmin >= 1, length(ulam) must equal nlam.");

  // clone maj to avoid mutating the caller's vector (core scales it by mval)
  return huber_path_core(alpha, lam2, hval, clone(maj), mval,
                         X, y, ju, pfncol, pf, pf2,
                         dfmax, pmax, nlam, flmin, ulam,
                         eps, maxit, istrong);
}

// --------------------------- Exported: huber_cd_cpp ---------------------------

// [[Rcpp::export]]
List huber_cd_cpp(double alpha,
                  double lam2,
                  double hval,
                  int nobs,
                  int nvars,
                  NumericMatrix X,
                  NumericVector y,
                  IntegerVector jd,
                  int pfncol,
                  NumericMatrix pf,
                  NumericVector pf2,
                  int dfmax,
                  int pmax,
                  int nlam,
                  double flmin,
                  NumericVector ulam,
                  double eps,
                  int isd,
                  int maxit,
                  int istrong) {

  if (X.nrow() != nobs || X.ncol() != nvars) stop("X must be nobs x nvars.");
  if (y.size() != nobs) stop("y must have length nobs.");
  if (pf.nrow() != nvars) stop("pf must have nrow=nvars.");
  if (!(pfncol == 1 || pf.ncol() == nlam))
    stop("pfncol must be 1 or equal to nlam and match pf.ncol().");
  if (pf2.size() != nvars) stop("pf2 must have length nvars.");

  IntegerVector ju = chkvars_cpp(X);

  if (jd.size() > 0 && jd[0] > 0) {
    const int deln = jd[0];
    if (jd.size() < deln + 1) stop("jd length inconsistent with jd[1].");
    for (int k = 0; k < deln; ++k) {
      const int idx1 = jd[k + 1];
      if (idx1 >= 1 && idx1 <= nvars) ju[idx1 - 1] = 0;
    }
  }

  {
    int anyju = 0;
    for (int j = 0; j < nvars; ++j) if (ju[j] != 0) { anyju = 1; break; }
    if (!anyju)
      return List::create(_["nalam"]=0, _["b0"]=NumericVector(0), _["beta"]=NumericMatrix(pmax, 0),
                          _["ibeta"]=IntegerVector(pmax), _["nbeta"]=IntegerVector(0),
                          _["alam"]=NumericVector(0), _["npass"]=0, _["jerr"]=7777);
  }

  {
    double maxpf = R_NegInf, maxpf2 = R_NegInf;
    for (int j = 0; j < pf.nrow(); ++j)
      for (int c = 0; c < pf.ncol(); ++c)
        if (pf(j, c) > maxpf) maxpf = pf(j, c);
    for (int j = 0; j < pf2.size(); ++j)
      if (pf2[j] > maxpf2) maxpf2 = pf2[j];
    if (!(maxpf > 0.0) || !(maxpf2 > 0.0))
      return List::create(_["nalam"]=0, _["b0"]=NumericVector(0), _["beta"]=NumericMatrix(pmax, 0),
                          _["ibeta"]=IntegerVector(pmax), _["nbeta"]=IntegerVector(0),
                          _["alam"]=NumericVector(0), _["npass"]=0, _["jerr"]=10000);
  }

  clamp_nonneg(pf);
  clamp_nonneg_vec(pf2);

  NumericMatrix Xstd = clone(X);
  NumericVector xmean, xnorm, maj;
  standardize_cpp(Xstd, ju, isd, xmean, xnorm, maj);

  List res = huber_path_core(alpha, lam2, hval, clone(maj), 1.0,
                             Xstd, y, ju, pfncol, pf, pf2,
                             dfmax, pmax, nlam, flmin, ulam,
                             eps, maxit, istrong);

  int jerr = as<int>(res["jerr"]);
  if (jerr > 0) return res;

  int nalam = as<int>(res["nalam"]);
  IntegerVector ibeta = res["ibeta"];
  IntegerVector nbeta = res["nbeta"];
  NumericMatrix beta  = res["beta"];
  NumericVector b0    = res["b0"];

  if (nalam > 0) {
    for (int l = 0; l < nalam; ++l) {
      int nk = nbeta[l];
      if (isd == 1) {
        for (int j = 0; j < nk; ++j) {
          const int col = ibeta[j] - 1;
          if (xnorm[col] != 0.0) beta(j, l) /= xnorm[col];
        }
      }
      double adj = 0.0;
      for (int j = 0; j < nk; ++j) adj += beta(j, l) * xmean[ibeta[j] - 1];
      b0[l] -= adj;
    }
  }

  NumericMatrix beta_full(nvars, nalam);
  std::fill(beta_full.begin(), beta_full.end(), 0.0);
  for (int l = 0; l < nalam; ++l) {
    const int nk = nbeta[l];
    for (int j = 0; j < nk; ++j) beta_full(ibeta[j] - 1, l) = beta(j, l);
  }

  res["beta_full"] = beta_full;
  res["beta"] = beta;
  res["b0"]   = b0;
  return res;
}

// --------------------------- Exported: huber_drv_cpp ---------------------------

// [[Rcpp::export]]
NumericVector huber_drv_cpp(const NumericMatrix& X,
                            const NumericVector& r,
                            const double hval) {
  const int n = X.nrow();
  if (r.size() != n) stop("Length of r must match nrow(X).");
  return huber_drv_internal(X, r, hval);
}
