// File: mm_angle.c
//
// Drop-in replacement for angle.c. Exposes the same
// template <typename HESSIAN> struct Angle interface (PositionSize,
// energy/gradient/hessian and *_fd siblings, taking const angle_term&),
// but the *analytical* methods drive the legacy Mathematica-derived
// _Angle_termCode.cc via the ANGLE_* macro family. The *_fd methods are
// central-difference references that call the local energy(), so the FD
// reference automatically matches the macro-based energy law here.
//
// Use as an A/B test against the new symbolic kernel: swap the include in
// energyAngle.cc from
//   #include "cando/chem/energyKernels/angle.c"
// to
//   #include "cando/chem/energyKernels/mm_angle.c"
// rebuild, and re-run mmin with only the angle component enabled.

#include "angle.h"

template <typename HESSIAN>
struct Angle {
  static constexpr size_t PositionSize = 9;
  static std::string description() { return "macrokernel-angle"; }

  double energy(const angle_term& term, double* position, double* energy_accumulate) {
    double energy_added = 0.0;
    // _Angle_termCode.cc expects locals named kt, t0, I1, I2, I3 in scope.
    double kt = term.kt;
    double t0 = term.t0;
    int    I1 = term.i3x1;
    int    I2 = term.i3x2;
    int    I3 = term.i3x3;
    #define USE_EXPLICIT_DECLARES 1
    #define DECLARE_FLOAT(xx) double xx;
    #define ArcCos(xx) acos(xx)
    #define ANGLE_SET_PARAMETER(xx)
    #define ANGLE_APPLY_ATOM_MASK(ii1,ii2,ii3)
    #define ANGLE_SET_POSITION(xx,ii,of) xx = position[ii+of];
    #define ANGLE_ENERGY_ACCUMULATE(eee) { energy_added += (eee); *energy_accumulate += (eee); }
    #define ANGLE_FORCE_ACCUMULATE(ii,oo,ff)
    #define ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(ii1,oo1,ii2,oo2,vv)
    #define ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(ii1,oo1,ii2,oo2,vv)
    #define ANGLE_DEBUG_INTERACTIONS(ii1,ii2,ii3)
    int calcForce = 0;
    int calcDiagonalHessian = 0;
    int calcOffDiagonalHessian = 0;
    (void)calcForce; (void)calcDiagonalHessian; (void)calcOffDiagonalHessian;
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-variable"
    #pragma clang diagnostic ignored "-Wunused-label"
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic ignored "-Wunused-label"
    // _Angle_termDeclares.cc does NOT declare x1..z3; declare them here.
    // The y1 local also shadows libm's Bessel y1(double).
    double x1, y1, z1, x2, y2, z2, x3, y3, z3;
    int IllegalAngle;
    #include <cando/chem/energy_functions/_Angle_termDeclares.cc>
    #include <cando/chem/energy_functions/_Angle_termCode.cc>
    #pragma GCC diagnostic pop
    #pragma clang diagnostic pop
    #undef USE_EXPLICIT_DECLARES
    #undef DECLARE_FLOAT
    #undef ArcCos
    #undef ANGLE_SET_PARAMETER
    #undef ANGLE_APPLY_ATOM_MASK
    #undef ANGLE_SET_POSITION
    #undef ANGLE_ENERGY_ACCUMULATE
    #undef ANGLE_FORCE_ACCUMULATE
    #undef ANGLE_DIAGONAL_HESSIAN_ACCUMULATE
    #undef ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE
    #undef ANGLE_DEBUG_INTERACTIONS
    return energy_added;
  }

  double gradient(const angle_term& term, double* position, double* energy_accumulate, double* force) {
    double energy_added = 0.0;
    double kt = term.kt;
    double t0 = term.t0;
    int    I1 = term.i3x1;
    int    I2 = term.i3x2;
    int    I3 = term.i3x3;
    #define USE_EXPLICIT_DECLARES 1
    #define DECLARE_FLOAT(xx) double xx;
    #define ArcCos(xx) acos(xx)
    #define ANGLE_SET_PARAMETER(xx)
    #define ANGLE_APPLY_ATOM_MASK(ii1,ii2,ii3)
    #define ANGLE_SET_POSITION(xx,ii,of) xx = position[ii+of];
    #define ANGLE_ENERGY_ACCUMULATE(eee) { energy_added += (eee); *energy_accumulate += (eee); }
    #define ANGLE_FORCE_ACCUMULATE(ii,oo,ff) force[ii+oo] += (ff);
    #define ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(ii1,oo1,ii2,oo2,vv)
    #define ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(ii1,oo1,ii2,oo2,vv)
    #define ANGLE_DEBUG_INTERACTIONS(ii1,ii2,ii3)
    #define ANGLE_CALC_FORCE 1
    int calcForce = 1;
    int calcDiagonalHessian = 0;
    int calcOffDiagonalHessian = 0;
    (void)calcDiagonalHessian; (void)calcOffDiagonalHessian;
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-variable"
    #pragma clang diagnostic ignored "-Wunused-label"
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic ignored "-Wunused-label"
    double x1, y1, z1, x2, y2, z2, x3, y3, z3;
    int IllegalAngle;
    #include <cando/chem/energy_functions/_Angle_termDeclares.cc>
    #include <cando/chem/energy_functions/_Angle_termCode.cc>
    #pragma GCC diagnostic pop
    #pragma clang diagnostic pop
    #undef USE_EXPLICIT_DECLARES
    #undef DECLARE_FLOAT
    #undef ArcCos
    #undef ANGLE_SET_PARAMETER
    #undef ANGLE_APPLY_ATOM_MASK
    #undef ANGLE_SET_POSITION
    #undef ANGLE_ENERGY_ACCUMULATE
    #undef ANGLE_FORCE_ACCUMULATE
    #undef ANGLE_DIAGONAL_HESSIAN_ACCUMULATE
    #undef ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE
    #undef ANGLE_DEBUG_INTERACTIONS
    #undef ANGLE_CALC_FORCE
    return energy_added;
  }

  double hessian(const angle_term& term, double* position, double* energy_accumulate, double* force,
                 HESSIAN hessian, double* dvec, double* hdvec) {
    double energy_added = 0.0;
    double kt = term.kt;
    double t0 = term.t0;
    int    I1 = term.i3x1;
    int    I2 = term.i3x2;
    int    I3 = term.i3x3;
    #define USE_EXPLICIT_DECLARES 1
    #define DECLARE_FLOAT(xx) double xx;
    #define ArcCos(xx) acos(xx)
    #define ANGLE_SET_PARAMETER(xx)
    #define ANGLE_APPLY_ATOM_MASK(ii1,ii2,ii3)
    #define ANGLE_SET_POSITION(xx,ii,of) xx = position[ii+of];
    #define ANGLE_ENERGY_ACCUMULATE(eee) { energy_added += (eee); *energy_accumulate += (eee); }
    #define ANGLE_FORCE_ACCUMULATE(ii,oo,ff) force[ii+oo] += (ff);
    #define ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(ii1,oo1,ii2,oo2,vv) \
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, ii1, oo1, ii2, oo2, (vv) )
    #define ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(ii1,oo1,ii2,oo2,vv) \
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, ii1, oo1, ii2, oo2, (vv) )
    #define ANGLE_DEBUG_INTERACTIONS(ii1,ii2,ii3)
    #define ANGLE_CALC_FORCE 1
    #define ANGLE_CALC_DIAGONAL_HESSIAN 1
    #define ANGLE_CALC_OFF_DIAGONAL_HESSIAN 1
    int calcForce = 1;
    int calcDiagonalHessian = 1;
    int calcOffDiagonalHessian = 1;
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-variable"
    #pragma clang diagnostic ignored "-Wunused-label"
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic ignored "-Wunused-label"
    double x1, y1, z1, x2, y2, z2, x3, y3, z3;
    int IllegalAngle;
    #include <cando/chem/energy_functions/_Angle_termDeclares.cc>
    #include <cando/chem/energy_functions/_Angle_termCode.cc>
    #pragma GCC diagnostic pop
    #pragma clang diagnostic pop
    #undef USE_EXPLICIT_DECLARES
    #undef DECLARE_FLOAT
    #undef ArcCos
    #undef ANGLE_SET_PARAMETER
    #undef ANGLE_APPLY_ATOM_MASK
    #undef ANGLE_SET_POSITION
    #undef ANGLE_ENERGY_ACCUMULATE
    #undef ANGLE_FORCE_ACCUMULATE
    #undef ANGLE_DIAGONAL_HESSIAN_ACCUMULATE
    #undef ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE
    #undef ANGLE_DEBUG_INTERACTIONS
    #undef ANGLE_CALC_FORCE
    #undef ANGLE_CALC_DIAGONAL_HESSIAN
    #undef ANGLE_CALC_OFF_DIAGONAL_HESSIAN
    return energy_added;
  }

  // ---------------------------------------------------------------------
  // Finite-difference reference implementations
  // ---------------------------------------------------------------------

  void energy_fd(const angle_term& term, double* position, double* energy_accumulate) {
    energy(term, position, energy_accumulate);
  }

  void gradient_fd(const angle_term& term, double* position, double* energy_accumulate, double* force) {
    constexpr size_t PS = PositionSize;
    const double h = 1.0e-5;
    const double inv2h = 1.0/(2.0*h);
    double e0 = 0.0;
    energy(term, position, &e0);
    if (energy_accumulate) { *energy_accumulate += e0; }
    auto fd_dE_d = [&](size_t idx) {
      double saved = position[idx];
      double e_plus = 0.0, e_minus = 0.0;
      position[idx] = saved + h;
      energy(term, position, &e_plus);
      position[idx] = saved - h;
      energy(term, position, &e_minus);
      position[idx] = saved;
      return (e_plus - e_minus) * inv2h;
    };
    size_t atom[PS] = { (size_t)term.i3x1, (size_t)term.i3x1, (size_t)term.i3x1,
                        (size_t)term.i3x2, (size_t)term.i3x2, (size_t)term.i3x2,
                        (size_t)term.i3x3, (size_t)term.i3x3, (size_t)term.i3x3 };
    size_t off [PS] = { 0, 1, 2, 0, 1, 2, 0, 1, 2 };
    for (size_t r = 0; r < PS; ++r) {
      KernelGradientAcc(atom[r], off[r], fd_dE_d(atom[r] + off[r]));
    }
  }

  void hessian_fd(const angle_term& term, double* position, double* energy_accumulate, double* force,
                  HESSIAN hessian, double* dvec, double* hdvec) {
    constexpr size_t PS = PositionSize;
    const double h = 1.0e-5;
    const double inv2h = 1.0/(2.0*h);
    const double invh2 = 1.0/(h*h);
    double e0 = 0.0;
    energy(term, position, &e0);
    if (energy_accumulate) { *energy_accumulate += e0; }
    auto fd_dE_d = [&](size_t idx) {
      double saved = position[idx];
      double e_plus = 0.0, e_minus = 0.0;
      position[idx] = saved + h;
      energy(term, position, &e_plus);
      position[idx] = saved - h;
      energy(term, position, &e_minus);
      position[idx] = saved;
      return (e_plus - e_minus) * inv2h;
    };
    auto fd_d2E_dii = [&](size_t idx) {
      double saved = position[idx];
      double e_plus = 0.0, e_minus = 0.0;
      position[idx] = saved + h;
      energy(term, position, &e_plus);
      position[idx] = saved - h;
      energy(term, position, &e_minus);
      position[idx] = saved;
      return (e_plus + e_minus - 2.0*e0) * invh2;
    };
    auto fd_d2E_dij = [&](size_t idx_i, size_t idx_j) {
      double saved_i = position[idx_i];
      double saved_j = position[idx_j];
      double e_pp=0.0, e_pm=0.0, e_mp=0.0, e_mm=0.0;
      position[idx_i] = saved_i + h; position[idx_j] = saved_j + h;
      energy(term, position, &e_pp);
      position[idx_j] = saved_j - h;
      energy(term, position, &e_pm);
      position[idx_i] = saved_i - h; position[idx_j] = saved_j + h;
      energy(term, position, &e_mp);
      position[idx_j] = saved_j - h;
      energy(term, position, &e_mm);
      position[idx_i] = saved_i; position[idx_j] = saved_j;
      return (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    };
    size_t atom[PS] = { (size_t)term.i3x1, (size_t)term.i3x1, (size_t)term.i3x1,
                        (size_t)term.i3x2, (size_t)term.i3x2, (size_t)term.i3x2,
                        (size_t)term.i3x3, (size_t)term.i3x3, (size_t)term.i3x3 };
    size_t off [PS] = { 0, 1, 2, 0, 1, 2, 0, 1, 2 };
    for (size_t r = 0; r < PS; ++r) {
      KernelGradientAcc(atom[r], off[r], fd_dE_d(atom[r] + off[r]));
    }
    for (size_t r = 0; r < PS; ++r) {
      double hval = fd_d2E_dii(atom[r] + off[r]);
      KernelHessDiagAcc( PS, hessian, dvec, hdvec, atom[r], off[r], atom[r], off[r], hval );
    }
    for (size_t r = 1; r < PS; ++r) {
      for (size_t c = 0; c < r; ++c) {
        double hval = fd_d2E_dij(atom[r] + off[r], atom[c] + off[c]);
        KernelHessOffDiagAcc( PS, hessian, dvec, hdvec, atom[r], off[r], atom[c], off[c], hval );
      }
    }
  }

};
