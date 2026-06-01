#pragma once
#include <math.h>

namespace chem {

// Trig-free multiple-angle recurrence shared by the dihedral_fast_floored kernel
// (generated from energy-kernels.lisp) and the dihedralKernelCompare.cc
// agreement test.
//
// Given the un-normalized torsion components t1 (~ |c1||c2| sin phi) and
// t2 (~ |c1||c2| cos phi), compute cos(n*phi) and sin(n*phi) with NO
// transcendentals: normalize to cos(phi)/sin(phi) (one sqrt, floored so it never
// divides by zero at a linear torsion), then rotate n times by (cos phi, sin phi).
// n is the small integer periodicity (1..6).  This is the runtime-n analog of the
// old SinCos<N>::sinNPhiCosNPhi template.  cos(phi)=t2/|.|, sin(phi)=t1/|.| matches
// phi=atan2(t1,t2), so it agrees with the kernel's Blondel-Karplus dphi/dq.
inline void cos_sin_n_phi(int n, double t1, double t2, double &cn, double &sn) {
  const double FLOOR = 1e-12; // linear-torsion guard (matches energyComponent.h VERYSMALL)
  double inv = 1.0 / sqrt(t1 * t1 + t2 * t2 + FLOOR);
  double c = t2 * inv, s = t1 * inv; // cos(phi), sin(phi)
  double rc = 1.0, rs = 0.0;         // (cos 0, sin 0)
  for (int k = 0; k < n; ++k) {
    double nc = rc * c - rs * s;
    double ns = rs * c + rc * s;
    rc = nc;
    rs = ns;
  }
  cn = rc;
  sn = rs;
}
inline double cos_n_phi(double n, double t1, double t2) {
  double cn, sn;
  cos_sin_n_phi((int)n, t1, t2, cn, sn);
  return cn;
}
inline double sin_n_phi(double n, double t1, double t2) {
  double cn, sn;
  cos_sin_n_phi((int)n, t1, t2, cn, sn);
  return sn;
}

} // namespace chem
