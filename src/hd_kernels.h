// Inner loops shared by the coordinate-descent kernels.
//
// The reductions use several independent accumulators so that the latency of
// the floating-point add chain does not bound the throughput, and the clamps
// are written in the form that compilers turn into min/max instructions.  The
// summation order is fixed, so results are reproducible across platforms and
// no compiler-specific flags or pragmas are needed.
#ifndef HDSTATS_KERNELS_H
#define HDSTATS_KERNELS_H

#include <cmath>

namespace hd {

// sum_i a_i * b_i
inline double dot(const double* a, const double* b, int n) {
  double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0, s5 = 0.0, s6 = 0.0, s7 = 0.0;
  int i = 0;
  for (; i + 7 < n; i += 8) {
    s0 += a[i] * b[i];         s1 += a[i + 1] * b[i + 1];
    s2 += a[i + 2] * b[i + 2]; s3 += a[i + 3] * b[i + 3];
    s4 += a[i + 4] * b[i + 4]; s5 += a[i + 5] * b[i + 5];
    s6 += a[i + 6] * b[i + 6]; s7 += a[i + 7] * b[i + 7];
  }
  for (; i < n; ++i) s0 += a[i] * b[i];
  return ((s0 + s1) + (s2 + s3)) + ((s4 + s5) + (s6 + s7));
}

// sum_i a_i
inline double sum(const double* a, int n) {
  double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
  int i = 0;
  for (; i + 3 < n; i += 4) { s0 += a[i]; s1 += a[i + 1]; s2 += a[i + 2]; s3 += a[i + 3]; }
  for (; i < n; ++i) s0 += a[i];
  return (s0 + s1) + (s2 + s3);
}

// clamp r to [-h, h]
inline double clip(double r, double h) {
  double t = (r < h) ? r : h;
  return (t > -h) ? t : -h;
}

// clamp t to [lo, hi]
inline double clamp(double t, double lo, double hi) {
  t = (t < hi) ? t : hi;
  return (t > lo) ? t : lo;
}

// ---------------------------------------------------------------- Huber ----
// The Huber kernel keeps dl_i = clip(r_i, h) = psi(r_i) in sync with the
// residuals r_i = y_i - x_i'b - b0.

// r_i -= x_i * d ; dl_i = clip(r_i, h)
inline void huber_axpy(const double* x, double* r, double* dl, int n, double d, double h) {
  int i = 0;
  for (; i + 3 < n; i += 4) {
    const double r0 = r[i] - x[i] * d,         r1 = r[i + 1] - x[i + 1] * d;
    const double r2 = r[i + 2] - x[i + 2] * d, r3 = r[i + 3] - x[i + 3] * d;
    r[i] = r0; r[i + 1] = r1; r[i + 2] = r2; r[i + 3] = r3;
    dl[i] = clip(r0, h); dl[i + 1] = clip(r1, h); dl[i + 2] = clip(r2, h); dl[i + 3] = clip(r3, h);
  }
  for (; i < n; ++i) { const double ri = r[i] - x[i] * d; r[i] = ri; dl[i] = clip(ri, h); }
}

// r_i -= d ; dl_i = clip(r_i, h)   (intercept move)
inline void huber_shift(double* r, double* dl, int n, double d, double h) {
  int i = 0;
  for (; i + 3 < n; i += 4) {
    const double r0 = r[i] - d, r1 = r[i + 1] - d, r2 = r[i + 2] - d, r3 = r[i + 3] - d;
    r[i] = r0; r[i + 1] = r1; r[i + 2] = r2; r[i + 3] = r3;
    dl[i] = clip(r0, h); dl[i + 1] = clip(r1, h); dl[i + 2] = clip(r2, h); dl[i + 3] = clip(r3, h);
  }
  for (; i < n; ++i) { const double ri = r[i] - d; r[i] = ri; dl[i] = clip(ri, h); }
}

// ------------------------------------------------------------------ SVM ----
// The SVM kernel works with margins r_i = y_i (x_i'b + b0) and keeps
// dly_i = dl(r_i) * y_i, where dl is the derivative of the smoothed hinge loss:
//   dl = -1 if r <= 1 - h ; dl = (r - (1 + h)) / (2h) if 1 - h < r < 1 + h ; dl = 0 otherwise.
inline double svm_dl(double r, double hinv, double oneph) {
  return clamp(0.5 * hinv * (r - oneph), -1.0, 0.0);
}

inline void svm_refresh(const double* r, const double* y, double* dly, int n, double hinv, double oneph) {
  for (int i = 0; i < n; ++i) dly[i] = svm_dl(r[i], hinv, oneph) * y[i];
}

// r_i += y_i x_i d ; dly_i = dl(r_i) y_i
inline void svm_axpy(const double* x, const double* y, double* r, double* dly, int n,
                     double d, double hinv, double oneph) {
  int i = 0;
  for (; i + 3 < n; i += 4) {
    const double r0 = r[i] + y[i] * x[i] * d,             r1 = r[i + 1] + y[i + 1] * x[i + 1] * d;
    const double r2 = r[i + 2] + y[i + 2] * x[i + 2] * d, r3 = r[i + 3] + y[i + 3] * x[i + 3] * d;
    r[i] = r0; r[i + 1] = r1; r[i + 2] = r2; r[i + 3] = r3;
    dly[i] = svm_dl(r0, hinv, oneph) * y[i];         dly[i + 1] = svm_dl(r1, hinv, oneph) * y[i + 1];
    dly[i + 2] = svm_dl(r2, hinv, oneph) * y[i + 2]; dly[i + 3] = svm_dl(r3, hinv, oneph) * y[i + 3];
  }
  for (; i < n; ++i) { const double ri = r[i] + y[i] * x[i] * d; r[i] = ri; dly[i] = svm_dl(ri, hinv, oneph) * y[i]; }
}

// r_i += y_i d ; dly_i = dl(r_i) y_i   (intercept move)
inline void svm_shift(const double* y, double* r, double* dly, int n, double d, double hinv, double oneph) {
  int i = 0;
  for (; i + 3 < n; i += 4) {
    const double r0 = r[i] + y[i] * d,         r1 = r[i + 1] + y[i + 1] * d;
    const double r2 = r[i + 2] + y[i + 2] * d, r3 = r[i + 3] + y[i + 3] * d;
    r[i] = r0; r[i + 1] = r1; r[i + 2] = r2; r[i + 3] = r3;
    dly[i] = svm_dl(r0, hinv, oneph) * y[i];         dly[i + 1] = svm_dl(r1, hinv, oneph) * y[i + 1];
    dly[i + 2] = svm_dl(r2, hinv, oneph) * y[i + 2]; dly[i + 3] = svm_dl(r3, hinv, oneph) * y[i + 3];
  }
  for (; i < n; ++i) { const double ri = r[i] + y[i] * d; r[i] = ri; dly[i] = svm_dl(ri, hinv, oneph) * y[i]; }
}

// ------------------------------------------------------- quantile (QR) ----
// The quantile kernel keeps dl_i = dl(r_i) for the residuals r_i, where dl is
// the derivative of the smoothed check loss:
//   dl = 1 - tau if r < -h ; dl = -r/(2h) - tau + 1/2 if -h <= r <= h ; dl = -tau if r > h.
inline double qr_dl(double r, double hinv, double tau, double lo, double hi) {
  return clamp(-r * 0.5 * hinv - tau + 0.5, lo, hi);
}

inline void qr_refresh(const double* r, double* dl, int n, double hinv, double tau, double lo, double hi) {
  for (int i = 0; i < n; ++i) dl[i] = qr_dl(r[i], hinv, tau, lo, hi);
}

// r_i -= x_i d ; dl_i = dl(r_i)
inline void qr_axpy(const double* x, double* r, double* dl, int n, double d,
                    double hinv, double tau, double lo, double hi) {
  int i = 0;
  for (; i + 3 < n; i += 4) {
    const double r0 = r[i] - x[i] * d,         r1 = r[i + 1] - x[i + 1] * d;
    const double r2 = r[i + 2] - x[i + 2] * d, r3 = r[i + 3] - x[i + 3] * d;
    r[i] = r0; r[i + 1] = r1; r[i + 2] = r2; r[i + 3] = r3;
    dl[i] = qr_dl(r0, hinv, tau, lo, hi);     dl[i + 1] = qr_dl(r1, hinv, tau, lo, hi);
    dl[i + 2] = qr_dl(r2, hinv, tau, lo, hi); dl[i + 3] = qr_dl(r3, hinv, tau, lo, hi);
  }
  for (; i < n; ++i) { const double ri = r[i] - x[i] * d; r[i] = ri; dl[i] = qr_dl(ri, hinv, tau, lo, hi); }
}

// r_i -= d ; dl_i = dl(r_i)   (intercept move)
inline void qr_shift(double* r, double* dl, int n, double d, double hinv, double tau, double lo, double hi) {
  int i = 0;
  for (; i + 3 < n; i += 4) {
    const double r0 = r[i] - d, r1 = r[i + 1] - d, r2 = r[i + 2] - d, r3 = r[i + 3] - d;
    r[i] = r0; r[i + 1] = r1; r[i + 2] = r2; r[i + 3] = r3;
    dl[i] = qr_dl(r0, hinv, tau, lo, hi);     dl[i + 1] = qr_dl(r1, hinv, tau, lo, hi);
    dl[i + 2] = qr_dl(r2, hinv, tau, lo, hi); dl[i + 3] = qr_dl(r3, hinv, tau, lo, hi);
  }
  for (; i < n; ++i) { const double ri = r[i] - d; r[i] = ri; dl[i] = qr_dl(ri, hinv, tau, lo, hi); }
}

} // namespace hd

#endif
