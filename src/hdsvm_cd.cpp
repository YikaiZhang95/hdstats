// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace Rcpp;

// =====================================================================
// hdsvm: penalized SVM via uniform-kernel convolution-smoothed hinge.
// Rcpp translation of hdsvm_cd.f90 / utilities.f90.
// =====================================================================

inline double sign_copy(double a, double b) {
  return (b >= 0.0) ? std::abs(a) : -std::abs(a);
}

// Smoothed-hinge loss derivative wrt margin r:
//   dl = -1                       if r <= 1 - h
//   dl = 0.5/h * (r - (1+h))      if r <  1 + h
//   dl = 0                        otherwise
inline double svm_dl(double r, double hinv, double onemh, double oneph) {
  if (r <= onemh) return -1.0;
  if (r <  oneph) return 0.5 * hinv * (r - oneph);
  return 0.0;
}

// --------------------------- chkvars ---------------------------
static IntegerVector chkvars_cpp(const NumericMatrix& X) {
  const int n = X.nrow(), p = X.ncol();
  IntegerVector ju(p);
  for (int j = 0; j < p; ++j) {
    const double* col = &X(0, j);
    const double t = col[0];
    int nz = 0;
    for (int i = 1; i < n; ++i) if (col[i] != t) { nz = 1; break; }
    ju[j] = nz;
  }
  return ju;
}

// --------------------------- standardize ---------------------------
static void standardize_cpp(NumericMatrix& X, const IntegerVector& ju, int isd,
                            NumericVector& xmean, NumericVector& xnorm,
                            NumericVector& maj) {
  const int n = X.nrow(), p = X.ncol();
  xmean = NumericVector(p); xnorm = NumericVector(p); maj = NumericVector(p);
  for (int j = 0; j < p; ++j) {
    xnorm[j] = 1.0;
    if (ju[j] == 0) continue;
    double mu = 0.0;
    for (int i = 0; i < n; ++i) mu += X(i, j);
    mu /= n;
    xmean[j] = mu;
    for (int i = 0; i < n; ++i) X(i, j) -= mu;
    double ss = 0.0;
    for (int i = 0; i < n; ++i) ss += X(i, j) * X(i, j);
    double msq = ss / n;
    maj[j] = msq;
    if (isd == 1) {
      double sd = std::sqrt(msq);
      xnorm[j] = sd;
      if (sd > 0.0) { double inv = 1.0 / sd; for (int i = 0; i < n; ++i) X(i, j) *= inv; }
      maj[j] = 1.0;
    }
  }
}

// --------------------------- svm_drv ---------------------------
// vl_j = (1/n) * sum_i (dl_i * y_i) * x_ij
static void svm_drv(const NumericMatrix& X, const double* y, const double* r,
                    double* vl, double hinv, double onemh, double oneph) {
  const int n = X.nrow(), p = X.ncol();
  const double ninv = 1.0 / n;
  std::vector<double> dly(n);
  for (int i = 0; i < n; ++i) dly[i] = svm_dl(r[i], hinv, onemh, oneph) * y[i];
  for (int j = 0; j < p; ++j) {
    const double* xj = &X(0, j);
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += dly[i] * xj[i];
    vl[j] = s * ninv;
  }
}

// --------------------------- objective ---------------------------
// xi_i = max(0, 1 - y_i*(ka_i + intcpt));  obj = lam2/2*bb + mean(xi) + lam1*bt
static double objfun_enetsvm(double intcpt, double bb, double bt,
                             const double* ka, const double* y,
                             double lam1, double lam2, int n) {
  double s = 0.0;
  for (int i = 0; i < n; ++i) {
    double xi = 1.0 - y[i] * (ka[i] + intcpt);
    if (xi > 0.0) s += xi;
  }
  return 0.5 * lam2 * bb + s / n + lam1 * bt;
}

// Brent's method (ported from R's fmin / Fortran opt_int_enetsvm)
static double opt_int_enetsvm(double lmin, double lmax, int n, double bt,
                              const double* ka, double bb, const double* y,
                              double lam1, double lam2) {
  const double gold = (3.0 - std::sqrt(5.0)) * 0.5;
  double eps = std::numeric_limits<double>::epsilon();
  double tol = std::pow(eps, 0.25);
  eps = std::sqrt(eps);
  double tol3 = tol / 3.0;

  double a = lmin, b = lmax;
  double v = a + gold * (b - a), w = v, x = v;
  double d = 0.0, e = 0.0;
  double fx = objfun_enetsvm(x, bb, bt, ka, y, lam1, lam2, n);
  double fv = fx, fw = fx;

  for (;;) {
    double xm = (a + b) * 0.5;
    double tol1 = eps * std::abs(x) + tol3;
    double t2 = tol1 * 2.0;
    if (std::abs(x - xm) <= t2 - (b - a) * 0.5) break;
    double p = 0.0, q = 0.0, r = 0.0;
    if (std::abs(e) > tol1) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = (q - r) * 2.0;
      if (q > 0.0) p = -p; else q = -q;
      r = e; e = d;
    }
    if (std::abs(p) >= std::abs(q * 0.5 * r) || p <= q * (a - x) || p >= q * (b - x)) {
      e = (x < xm) ? (b - x) : (a - x);
      d = gold * e;
    } else {
      d = p / q;
      double u = x + d;
      if (u - a < t2 || b - u < t2) { d = tol1; if (x >= xm) d = -d; }
    }
    double u;
    if (std::abs(d) >= tol1) u = x + d;
    else if (d > 0.0) u = x + tol1;
    else u = x - tol1;
    double fu = objfun_enetsvm(u, bb, bt, ka, y, lam1, lam2, n);
    if (fu <= fx) {
      if (u < x) b = x; else a = x;
      v = w; w = x; x = u; fv = fw; fw = fx; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) { v = w; fv = fw; w = u; fw = fu; }
      else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
    }
  }
  return x;
}

// --------------------------- path core ---------------------------
static List hdsvm_path(double alpha, double lam2_in, double hval_in,
                       NumericVector maj0, const NumericMatrix& X,
                       const NumericVector& y, const IntegerVector& ju,
                       int pfncol, const NumericMatrix& pf, const NumericVector& pf2,
                       int dfmax, int pmax, int nlam, double flmin_in,
                       const NumericVector& ulam, double eps, int maxit,
                       double sigma, int is_exact) {

  const double BIG = 9.9e30, MFL = 1.0e-6;
  const int MNLAM = 6, hval_len = 4, mproj = 2;
  const double KKTeps = 1e-3, hvaleps = 1e-6;

  const int n = X.nrow(), p = X.ncol();
  const double ninv = 1.0 / n;
  const double projeps = n * eps;
  const double* yp = y.begin();

  std::vector<double> r(n), dl(n);
  for (int i = 0; i < n; ++i) r[i] = yp[i];

  std::vector<double> b(p + 1, 0.0), oldbeta(p + 1, 0.0), maj(p, 0.0);
  std::vector<double> ga(p, 0.0), vl(p, 0.0), sx(p, 0.0), ka(n), theta(n, 0.0);
  std::vector<int> m(pmax, 0), mm(p, 0), eset(p, 0), sset(n, 0), ss(n, 0);

  NumericVector b0_out(nlam), alam(nlam);
  NumericMatrix beta_out(pmax, nlam);
  IntegerVector nbeta_out(nlam), npass(nlam);
  int ni = 0, nalam = 0, jerr = 0;
  double al = 0.0, al0 = 0.0, alf = 0.01;

  int mnl = std::min(MNLAM, nlam);
  double flmin = flmin_in;
  if (flmin < 1.0) { flmin = std::max(MFL, flmin); if (nlam > 1) alf = std::pow(flmin, 1.0 / (nlam - 1.0)); }

  const double delta = hval_in;

  // ---- initial intercept-only fit + lambda_max gradient ----
  for (int i = 0; i < n; ++i) ka[i] = 0.0;
  double intercept = opt_int_enetsvm(-1e2, 1e2, n, 0.0, ka.data(), 0.0, yp, 0.0, 0.0);
  {
    std::vector<double> rtmp(n);
    for (int i = 0; i < n; ++i) rtmp[i] = yp[i] * intercept;
    svm_drv(X, yp, rtmp.data(), vl.data(), 1.0/1e-9, 1.0 - 1e-9, 1.0 + 1e-9);
  }
  for (int j = 0; j < p; ++j) ga[j] = std::abs(vl[j]);

  int pfl = 0;
  if (pfncol == 1) {
    al0 = 0.0;
    for (int j = 0; j < p; ++j)
      if (ju[j] != 0 && pf(j, 0) > 0.0) al0 = std::max(al0, ga[j] / pf(j, 0));
  }

  // ===================== lambda loop =====================
  for (int l = 0; l < nlam; ++l) {
    if (pfncol != 1) {
      pfl = l; al0 = 0.0;
      for (int j = 0; j < p; ++j)
        if (ju[j] != 0 && pf(j, pfl) > 0.0) al0 = std::max(al0, ga[j] / pf(j, pfl));
    }
    double hval = delta;
    int hval_id = 0;

    if (flmin >= 1.0) al = ulam[l];
    else if (l > 1) al = al * alf;
    else if (l == 0) al = BIG;
    else al = al0 * alf;

    if (al > al0) {
      for (int j = 1; j <= p; ++j) b[j] = 0.0;
      b[0] = intercept;
      for (int i = 0; i < n; ++i) r[i] = yp[i] * b[0];
      for (int t = 0; t < ni; ++t) beta_out(t, l) = 0.0;
      nbeta_out[l] = 0; b0_out[l] = b[0]; alam[l] = al; nalam = l + 1;
      continue;
    }

    double lam2 = lam2_in;
    if (alpha != -1.0) lam2 = al * (1.0 - alpha) * 0.5 / alpha;

    // ---------------- hval annealing loop ----------------
    for (;;) {
      hval_id += 1;
      double hinv = 1.0 / hval, mval = hinv * 0.5;
      for (int j = 0; j < p; ++j) maj[j] = maj0[j] * mval;
      double onemh = 1.0 - hval, oneph = 1.0 + hval;

      oldbeta[0] = b[0];
      for (int t = 0; t < ni; ++t) oldbeta[m[t]] = b[m[t]];

      // ---------- middle loop ----------
      for (;;) {
        npass[l] += 1;
        double dif = 0.0;

        for (int k = 1; k <= p; ++k) {
          if (ju[k - 1] == 0) continue;
          const double* xk = &X(0, k - 1);
          double oldb = b[k], u = 0.0;
          for (int i = 0; i < n; ++i) u += svm_dl(r[i], hinv, onemh, oneph) * yp[i] * xk[i];
          u = maj[k - 1] * b[k] - u * ninv;
          double v = std::abs(u) - al * pf(k - 1, pfl);
          b[k] = (v > 0.0) ? sign_copy(v, u) / (pf2[k - 1] * lam2 + maj[k - 1]) : 0.0;
          double d = b[k] - oldb;
          if (d != 0.0) {
            dif = std::max(dif, d * d);
            for (int i = 0; i < n; ++i) r[i] += yp[i] * xk[i] * d;
            if (mm[k - 1] == 0) {
              ni += 1;
              if (ni > pmax) { jerr = -10000 - (l + 1); goto done; }
              mm[k - 1] = ni; m[ni - 1] = k;
            }
          }
        }

        // intercept
        {
          double d = 0.0;
          for (int i = 0; i < n; ++i) d += svm_dl(r[i], hinv, onemh, oneph) * yp[i];
          d = -d * ninv / mval;
          if (d != 0.0) { b[0] += d; for (int i = 0; i < n; ++i) r[i] += yp[i] * d; dif = std::max(dif, d * d); }
        }

        if (dif < eps) break;
        { long sp = 0; for (int q = 0; q < nlam; ++q) sp += npass[q]; if (sp > maxit) { jerr = -1; goto done; } }

        // ---------- inner loop ----------
        for (;;) {
          npass[l] += 1;
          double difi = 0.0;
          for (int t = 0; t < ni; ++t) {
            int k = m[t];
            const double* xk = &X(0, k - 1);
            double oldb = b[k], u = 0.0;
            for (int i = 0; i < n; ++i) u += svm_dl(r[i], hinv, onemh, oneph) * yp[i] * xk[i];
            u = maj[k - 1] * b[k] - u * ninv;
            double v = std::abs(u) - al * pf(k - 1, pfl);
            b[k] = (v > 0.0) ? sign_copy(v, u) / (pf2[k - 1] * lam2 + maj[k - 1]) : 0.0;
            double d = b[k] - oldb;
            if (d != 0.0) { difi = std::max(difi, d * d); for (int i = 0; i < n; ++i) r[i] += yp[i] * xk[i] * d; }
          }
          {
            double d = 0.0;
            for (int i = 0; i < n; ++i) d += svm_dl(r[i], hinv, onemh, oneph) * yp[i];
            d = -d * ninv / mval;
            if (d != 0.0) { b[0] += d; for (int i = 0; i < n; ++i) r[i] += yp[i] * d; difi = std::max(difi, d * d); }
          }
          if (difi < eps) break;
          { long sp = 0; for (int q = 0; q < nlam; ++q) sp += npass[q]; if (sp > maxit) { jerr = -1; goto done; } }
        }
      } // middle

      // ---- exact intercept optimization via Brent ----
      double dif;
      {
        double ab = 0.0, bb = 0.0;
        for (int j = 1; j <= p; ++j) { ab += std::abs(b[j]); bb += b[j] * b[j]; }
        for (int i = 0; i < n; ++i) ka[i] = yp[i] * r[i] - b[0];
        double obj0 = objfun_enetsvm(b[0], bb, ab, ka.data(), yp, al, lam2, n);
        double b00 = opt_int_enetsvm(-1e2, 1e2, n, ab, ka.data(), bb, yp, al, lam2);
        double obj1 = objfun_enetsvm(b00, bb, ab, ka.data(), yp, al, lam2, n);
        if (obj1 < obj0) { double d0 = b00 - b[0]; for (int i = 0; i < n; ++i) r[i] += yp[i] * d0; b[0] = b00; }
      }
      if (ni > pmax) { jerr = -10000 - (l + 1); goto done; }

      // ---- KKT check ----
      dif = (ni > 0) ? 0.0 : 1e300;
      for (int t = 0; t < ni; ++t) { double q = std::abs(oldbeta[m[t]] - b[m[t]]); dif = std::max(dif, q); }
      svm_drv(X, yp, r.data(), vl.data(), 1.0/1e-12, 1.0 - 1e-12, 1.0 + 1e-12);
      double dif_kkt = 0.0;
      for (int j = 1; j <= p; ++j) {
        double d_kkt = 0.0;
        if (std::abs(b[j]) > 0.0) {
          double u = vl[j - 1] + lam2 * pf2[j - 1] * b[j];
          double v = std::abs(u) - al * pf(j - 1, pfl);
          d_kkt = sign_copy(v, u);
        } else {
          double v = std::abs(vl[j - 1]) - al * pf(j - 1, pfl);
          if (v > 0.0) d_kkt = v;
        }
        if (d_kkt != 0.0) dif_kkt = std::max(dif_kkt, d_kkt * d_kkt);
      }
      {
        double uo = std::max(al, 1.0);
        if (dif_kkt * dif_kkt / uo / uo < KKTeps && dif * dif < hvaleps) break;
      }

      if (is_exact == 0) break;

      // ---------------- exact projection (ADMM) ----------------
      if (hval < 0.125 * 0.125) {
        for (int j = 0; j < p; ++j) eset[j] = 0;
        int si = 0;
        for (int i = 0; i < n; ++i) { sset[i] = 0; ss[i] = 0; theta[i] = 0.0; }
        for (int mp = 0; mp < mproj; ++mp) {
          for (int i = 0; i < n; ++i)
            if (std::abs(1.0 - r[i]) < hval && ss[i] == 0) { sset[si++] = i; ss[i] = 1; }
          for (int k = 1; k <= p; ++k) if (std::abs(b[k]) < hval) eset[k - 1] = 1;

          for (;;) { // proj middle
            npass[l] += 1;
            double difp = 0.0;
            for (int k = 1; k <= p; ++k) {
              if (ju[k - 1] == 0) continue;
              const double* xk = &X(0, k - 1);
              double oldb = b[k];
              if (eset[k - 1] == 1) {
                b[k] = 0.0;
              } else {
                double u = 0.0;
                for (int i = 0; i < n; ++i) u += svm_dl(r[i], hinv, onemh, oneph) * yp[i] * xk[i];
                double mb = 0.0;
                if (si > 0) {
                  double s2 = 0.0, mbacc = 0.0;
                  for (int t = 0; t < si; ++t) {
                    int idx = sset[t];
                    s2 += xk[idx] * xk[idx];
                    mbacc += xk[idx] * (theta[t] + sigma * (1.0 - r[idx])) / yp[idx];
                  }
                  sx[k - 1] = s2;
                  mb = sigma * s2 * b[k] + mbacc;
                }
                u = mb + maj[k - 1] * b[k] - u * ninv;
                double v = std::abs(u) - al * pf(k - 1, pfl);
                if (v > 0.0)
                  b[k] = sign_copy(v, u) / (pf2[k - 1] * lam2 + maj[k - 1] + (si > 0 ? sigma * sx[k - 1] : 0.0));
                else b[k] = 0.0;
              }
              double d = b[k] - oldb;
              if (d != 0.0) {
                difp = std::max(difp, d * d);
                for (int i = 0; i < n; ++i) r[i] += yp[i] * xk[i] * d;
                if (mm[k - 1] == 0) { ni += 1; if (ni > pmax) { jerr = -10000 - (l + 1); goto done; } mm[k - 1] = ni; m[ni - 1] = k; }
              }
            }
            // intercept (proj middle)
            {
              double d = 0.0;
              for (int i = 0; i < n; ++i) d += svm_dl(r[i], hinv, onemh, oneph) * yp[i];
              d = -d * ninv;
              if (si > 0) { double acc = 0.0; for (int t = 0; t < si; ++t) acc += (theta[t] + sigma * (1.0 - r[sset[t]])) / yp[sset[t]]; d = (d + acc) / (sigma * si + mval); }
              else d = d / mval;
              if (d != 0.0) { b[0] += d; for (int i = 0; i < n; ++i) r[i] += yp[i] * d; difp = std::max(difp, d * d); }
            }
            if (si > 0) for (int t = 0; t < si; ++t) theta[t] += sigma * (1.0 - r[sset[t]]);
            if (difp < eps) break;
            { long sp = 0; for (int q = 0; q < nlam; ++q) sp += npass[q]; if (sp > maxit) { jerr = -1; goto done; } }

            for (;;) { // proj inner
              npass[l] += 1;
              double difpi = 0.0;
              for (int t = 0; t < ni; ++t) {
                int k = m[t];
                const double* xk = &X(0, k - 1);
                double oldb = b[k];
                if (eset[k - 1] == 1) { b[k] = 0.0; }
                else {
                  double u = 0.0;
                  for (int i = 0; i < n; ++i) u += svm_dl(r[i], hinv, onemh, oneph) * yp[i] * xk[i];
                  double mb = 0.0;
                  if (si > 0) {
                    double s2 = 0.0, mbacc = 0.0;
                    for (int tt = 0; tt < si; ++tt) { int idx = sset[tt]; s2 += xk[idx] * xk[idx]; mbacc += xk[idx] * (theta[tt] + sigma * (1.0 - r[idx])) / yp[idx]; }
                    sx[k - 1] = s2; mb = sigma * s2 * b[k] + mbacc;
                  }
                  u = mb + maj[k - 1] * b[k] - u * ninv;
                  double v = std::abs(u) - al * pf(k - 1, pfl);
                  if (v > 0.0) b[k] = sign_copy(v, u) / (pf2[k - 1] * lam2 + maj[k - 1] + (si > 0 ? sigma * sx[k - 1] : 0.0));
                  else b[k] = 0.0;
                }
                double d = b[k] - oldb;
                if (d != 0.0) { difpi = std::max(difpi, d * d); for (int i = 0; i < n; ++i) r[i] += yp[i] * xk[i] * d; }
              }
              {
                double d = 0.0;
                for (int i = 0; i < n; ++i) d += svm_dl(r[i], hinv, onemh, oneph) * yp[i];
                d = -d * ninv;
                if (si > 0) { double acc = 0.0; for (int t = 0; t < si; ++t) acc += (theta[t] + sigma * (1.0 - r[sset[t]])) / yp[sset[t]]; d = (d + acc) / (sigma * si + mval); }
                else d = d / mval;
                if (d != 0.0) { b[0] += d; for (int i = 0; i < n; ++i) r[i] += yp[i] * d; difpi = std::max(difpi, d * d); }
              }
              if (si > 0) for (int t = 0; t < si; ++t) theta[t] += sigma * (1.0 - r[sset[t]]);
              if (difpi < projeps) break;
              { long sp = 0; for (int q = 0; q < nlam; ++q) sp += npass[q]; if (sp > maxit) { jerr = -1; goto done; } }
            }
          }
        }
      }

      if (hval_id == hval_len) break;
      hval *= 0.125;
    } // hval loop

    // ---- save ----
    if (ni > pmax) { jerr = -10000 - (l + 1); goto done; }
    for (int t = 0; t < ni; ++t) beta_out(t, l) = b[m[t]];
    nbeta_out[l] = ni; b0_out[l] = b[0]; alam[l] = al; nalam = l + 1;
    if (l + 1 < mnl || flmin >= 1.0) continue;
    int me = 0; for (int t = 0; t < ni; ++t) if (beta_out(t, l) != 0.0) ++me;
    if (me > dfmax) break;
  }

done:
  IntegerVector ibeta(pmax);
  for (int t = 0; t < pmax; ++t) ibeta[t] = m[t];
  return List::create(_["nalam"] = nalam, _["b0"] = b0_out, _["beta"] = beta_out,
                      _["ibeta"] = ibeta, _["nbeta"] = nbeta_out, _["alam"] = alam,
                      _["npass"] = npass, _["jerr"] = jerr);
}

// --------------------------- exported wrapper ---------------------------
// [[Rcpp::export]]
List hdsvm_cd_cpp(double alpha, double lam2, double hval, int nobs, int nvars,
                  NumericMatrix X, NumericVector y, IntegerVector jd, int pfncol,
                  NumericMatrix pf, NumericVector pf2, int dfmax, int pmax, int nlam,
                  double flmin, NumericVector ulam, double eps, int isd, int maxit,
                  double sigma, int is_exact) {
  if (X.nrow() != nobs || X.ncol() != nvars) stop("X must be nobs x nvars.");
  if (y.size() != nobs) stop("y must have length nobs.");
  if (pf.nrow() != nvars) stop("pf must have nrow=nvars.");
  if (!(pfncol == 1 || pf.ncol() == nlam)) stop("pfncol must be 1 or nlam.");
  if (pf2.size() != nvars) stop("pf2 must have length nvars.");

  IntegerVector ju = chkvars_cpp(X);
  if (jd.size() > 0 && jd[0] > 0) {
    const int deln = jd[0];
    for (int k = 0; k < deln; ++k) { int idx = jd[k + 1]; if (idx >= 1 && idx <= nvars) ju[idx - 1] = 0; }
  }
  { int any = 0; for (int j = 0; j < nvars; ++j) if (ju[j]) { any = 1; break; }
    if (!any) return List::create(_["nalam"]=0,_["jerr"]=7777); }

  for (int c = 0; c < pf.ncol(); ++c) for (int j = 0; j < nvars; ++j) if (pf(j,c) < 0.0) pf(j,c) = 0.0;
  for (int j = 0; j < nvars; ++j) if (pf2[j] < 0.0) pf2[j] = 0.0;

  NumericMatrix Xstd = clone(X);
  NumericVector xmean, xnorm, maj;
  standardize_cpp(Xstd, ju, isd, xmean, xnorm, maj);

  List res = hdsvm_path(alpha, lam2, hval, maj, Xstd, y, ju, pfncol, pf, pf2,
                        dfmax, pmax, nlam, flmin, ulam, eps, maxit, sigma, is_exact);

  int jerr = as<int>(res["jerr"]);
  if (jerr > 0) return res;

  int nalam = as<int>(res["nalam"]);
  IntegerVector ibeta = res["ibeta"]; IntegerVector nbeta = res["nbeta"];
  NumericMatrix beta = res["beta"]; NumericVector b0 = res["b0"];
  for (int l = 0; l < nalam; ++l) {
    int nk = nbeta[l];
    if (isd == 1) for (int j = 0; j < nk; ++j) { int col = ibeta[j] - 1; if (xnorm[col] != 0.0) beta(j, l) /= xnorm[col]; }
    double adj = 0.0; for (int j = 0; j < nk; ++j) adj += beta(j, l) * xmean[ibeta[j] - 1];
    b0[l] -= adj;
  }
  NumericMatrix beta_full(nvars, nalam);
  for (int l = 0; l < nalam; ++l) { int nk = nbeta[l]; for (int j = 0; j < nk; ++j) beta_full(ibeta[j] - 1, l) = beta(j, l); }
  res["beta_full"] = beta_full; res["beta"] = beta; res["b0"] = b0;
  return res;
}
