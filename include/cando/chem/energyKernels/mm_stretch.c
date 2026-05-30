// File: mm_stretch.c
//
// Drop-in replacement for stretch.c. Exposes the same
// template <typename HESSIAN> struct Stretch interface (PositionSize,
// energy/gradient/hessian and *_fd siblings, taking const TermStretch&),
// but the *analytical* methods drive the legacy Mathematica-derived
// _stretch_termCode.cc via the STRETCH_* macro family. The *_fd methods
// are central-difference references that call the local energy(), so the
// FD reference automatically matches the macro-based energy law here.
//
// Use as an A/B test against the new symbolic kernel: swap the include in
// energyStretch.cc from
//   #include "cando/chem/energyKernels/stretch.c"
// to
//   #include "cando/chem/energyKernels/mm_stretch.c"
// rebuild, and re-run mmin. If chem:run-test-calls now reports zero
// mismatches on the jostled initial state, the bug is in the symbolic
// derivation in stretch.c, not in the FD machinery or the kernel-to-matrix
// scatter.

#include "stretch.h"

template <typename HESSIAN>
struct Stretch {
  static constexpr size_t PositionSize = 6;
  static std::string description() { return "macrokernel-stretch"; }

  double energy(const TermStretch& term, double* position, double* energy_accumulate) {
    double energy_added = 0.0;
    // _stretch_termCode.cc expects locals named kb, r0, I1, I2 in scope.
    // Pull them out of term so the legacy macros work unchanged.
    double kb = term.kb;
    double r0 = term.r0;
    int    I1 = term.i3x1;
    int    I2 = term.i3x2;
    #define USE_EXPLICIT_DECLARES 1
    #define DECLARE_FLOAT(xx) double xx;
    #define STRETCH_SET_PARAMETER(xx)
    #define STRETCH_APPLY_ATOM_MASK(ii1,ii2)
    #define STRETCH_SET_POSITION(xx,ii,of) xx = position[ii+of];
    #define STRETCH_ENERGY_ACCUMULATE(eee) { energy_added += (eee); *energy_accumulate += (eee); }
    #define STRETCH_FORCE_ACCUMULATE(ii,oo,ff)
    #define STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(ii1,oo1,ii2,oo2,vv)
    #define STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(ii1,oo1,ii2,oo2,vv)
    #define STRETCH_DEBUG_INTERACTIONS(ii1,ii2)
    // STRETCH_CALC_FORCE / *_HESSIAN intentionally undefined: the macro
    // file #ifdef-guards skip the corresponding blocks at compile time.
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
    // _stretch_termDeclares.cc does NOT declare x1..z2; declare them here.
    // The y1 local also shadows libm's Bessel y1(double), which is the
    // actual cause of the "non-object type ... not assignable" error.
    double x1, y1, z1, x2, y2, z2;
    #include <cando/chem/energy_functions/_Stretch_termDeclares.cc>
    #include <cando/chem/energy_functions/_Stretch_termCode.cc>
    #pragma GCC diagnostic pop
    #pragma clang diagnostic pop
    #undef USE_EXPLICIT_DECLARES
    #undef DECLARE_FLOAT
    #undef STRETCH_SET_PARAMETER
    #undef STRETCH_APPLY_ATOM_MASK
    #undef STRETCH_SET_POSITION
    #undef STRETCH_ENERGY_ACCUMULATE
    #undef STRETCH_FORCE_ACCUMULATE
    #undef STRETCH_DIAGONAL_HESSIAN_ACCUMULATE
    #undef STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE
    #undef STRETCH_DEBUG_INTERACTIONS
    return energy_added;
  }

  double gradient(const TermStretch& term, double* position, double* energy_accumulate, double* force) {
    double energy_added = 0.0;
    double kb = term.kb;
    double r0 = term.r0;
    int    I1 = term.i3x1;
    int    I2 = term.i3x2;
    #define USE_EXPLICIT_DECLARES 1
    #define DECLARE_FLOAT(xx) double xx;
    #define STRETCH_SET_PARAMETER(xx)
    #define STRETCH_APPLY_ATOM_MASK(ii1,ii2)
    #define STRETCH_SET_POSITION(xx,ii,of) xx = position[ii+of];
    #define STRETCH_ENERGY_ACCUMULATE(eee) { energy_added += (eee); *energy_accumulate += (eee); }
    // Legacy convention: STRETCH_FORCE_ACCUMULATE receives a FORCE value
    // (fx1 = -gx1 etc.), so we add it directly to force[].
    #define STRETCH_FORCE_ACCUMULATE(ii,oo,ff) force[ii+oo] += (ff);
    #define STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(ii1,oo1,ii2,oo2,vv)
    #define STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(ii1,oo1,ii2,oo2,vv)
    #define STRETCH_DEBUG_INTERACTIONS(ii1,ii2)
    #define STRETCH_CALC_FORCE 1
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
    double x1, y1, z1, x2, y2, z2;
    #include <cando/chem/energy_functions/_Stretch_termDeclares.cc>
    #include <cando/chem/energy_functions/_Stretch_termCode.cc>
    #pragma GCC diagnostic pop
    #pragma clang diagnostic pop
    #undef USE_EXPLICIT_DECLARES
    #undef DECLARE_FLOAT
    #undef STRETCH_SET_PARAMETER
    #undef STRETCH_APPLY_ATOM_MASK
    #undef STRETCH_SET_POSITION
    #undef STRETCH_ENERGY_ACCUMULATE
    #undef STRETCH_FORCE_ACCUMULATE
    #undef STRETCH_DIAGONAL_HESSIAN_ACCUMULATE
    #undef STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE
    #undef STRETCH_DEBUG_INTERACTIONS
    #undef STRETCH_CALC_FORCE
    return energy_added;
  }

  double hessian(const TermStretch& term, double* position, double* energy_accumulate, double* force,
                 HESSIAN hessian, double* dvec, double* hdvec) {
    double energy_added = 0.0;
    double kb = term.kb;
    double r0 = term.r0;
    int    I1 = term.i3x1;
    int    I2 = term.i3x2;
    #define USE_EXPLICIT_DECLARES 1
    #define DECLARE_FLOAT(xx) double xx;
    #define STRETCH_SET_PARAMETER(xx)
    #define STRETCH_APPLY_ATOM_MASK(ii1,ii2)
    #define STRETCH_SET_POSITION(xx,ii,of) xx = position[ii+of];
    #define STRETCH_ENERGY_ACCUMULATE(eee) { energy_added += (eee); *energy_accumulate += (eee); }
    #define STRETCH_FORCE_ACCUMULATE(ii,oo,ff) force[ii+oo] += (ff);
    // For HESSIAN=double* this writes hessian[(r,c)] (and the symmetric off-
    // diagonal entry) AND updates hdvec via the dvec dot product. For
    // HESSIAN=NoHessian the templated overload only updates hdvec.
    #define STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(ii1,oo1,ii2,oo2,vv) \
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, ii1, oo1, ii2, oo2, (vv) )
    #define STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(ii1,oo1,ii2,oo2,vv) \
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, ii1, oo1, ii2, oo2, (vv) )
    #define STRETCH_DEBUG_INTERACTIONS(ii1,ii2)
    #define STRETCH_CALC_FORCE 1
    #define STRETCH_CALC_DIAGONAL_HESSIAN 1
    #define STRETCH_CALC_OFF_DIAGONAL_HESSIAN 1
    int calcForce = 1;
    int calcDiagonalHessian = 1;
    int calcOffDiagonalHessian = 1;
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-variable"
    #pragma clang diagnostic ignored "-Wunused-label"
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic ignored "-Wunused-label"
    double x1, y1, z1, x2, y2, z2;
    #include <cando/chem/energy_functions/_Stretch_termDeclares.cc>
    #include <cando/chem/energy_functions/_Stretch_termCode.cc>
    #pragma GCC diagnostic pop
    #pragma clang diagnostic pop
    #undef USE_EXPLICIT_DECLARES
    #undef DECLARE_FLOAT
    #undef STRETCH_SET_PARAMETER
    #undef STRETCH_APPLY_ATOM_MASK
    #undef STRETCH_SET_POSITION
    #undef STRETCH_ENERGY_ACCUMULATE
    #undef STRETCH_FORCE_ACCUMULATE
    #undef STRETCH_DIAGONAL_HESSIAN_ACCUMULATE
    #undef STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE
    #undef STRETCH_DEBUG_INTERACTIONS
    #undef STRETCH_CALC_FORCE
    #undef STRETCH_CALC_DIAGONAL_HESSIAN
    #undef STRETCH_CALC_OFF_DIAGONAL_HESSIAN
    return energy_added;
  }

  // ---------------------------------------------------------------------
  // Finite-difference reference implementations
  //
  // Identical layout to stretch.c::*_fd: they perturb the live position
  // array by +/- h, call energy() at each perturbed point, and combine the
  // values into the FD gradient / second derivative. Because they call the
  // local energy() above, the FD reference is automatically consistent
  // with the macro-based energy law in this file.
  // ---------------------------------------------------------------------

  void energy_fd(const TermStretch& term, double* position, double* energy_accumulate) {
    energy(term, position, energy_accumulate);
  }

  void gradient_fd(const TermStretch& term, double* position, double* energy_accumulate, double* force) {
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
    KernelGradientAcc(term.i3x1, 0, fd_dE_d(term.i3x1 + 0));
    KernelGradientAcc(term.i3x1, 1, fd_dE_d(term.i3x1 + 1));
    KernelGradientAcc(term.i3x1, 2, fd_dE_d(term.i3x1 + 2));
    KernelGradientAcc(term.i3x2, 0, fd_dE_d(term.i3x2 + 0));
    KernelGradientAcc(term.i3x2, 1, fd_dE_d(term.i3x2 + 1));
    KernelGradientAcc(term.i3x2, 2, fd_dE_d(term.i3x2 + 2));
    (void)PS;
  }

  void hessian_fd(const TermStretch& term, double* position, double* energy_accumulate, double* force,
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
    // Gradient via central differences.
    KernelGradientAcc(term.i3x1, 0, fd_dE_d(term.i3x1 + 0));
    KernelGradientAcc(term.i3x1, 1, fd_dE_d(term.i3x1 + 1));
    KernelGradientAcc(term.i3x1, 2, fd_dE_d(term.i3x1 + 2));
    KernelGradientAcc(term.i3x2, 0, fd_dE_d(term.i3x2 + 0));
    KernelGradientAcc(term.i3x2, 1, fd_dE_d(term.i3x2 + 1));
    KernelGradientAcc(term.i3x2, 2, fd_dE_d(term.i3x2 + 2));
    // Diagonal of the Hessian.
    KernelHessDiagAcc( PS, hessian, dvec, hdvec, term.i3x1, 0, term.i3x1, 0, fd_d2E_dii(term.i3x1 + 0));
    KernelHessDiagAcc( PS, hessian, dvec, hdvec, term.i3x1, 1, term.i3x1, 1, fd_d2E_dii(term.i3x1 + 1));
    KernelHessDiagAcc( PS, hessian, dvec, hdvec, term.i3x1, 2, term.i3x1, 2, fd_d2E_dii(term.i3x1 + 2));
    KernelHessDiagAcc( PS, hessian, dvec, hdvec, term.i3x2, 0, term.i3x2, 0, fd_d2E_dii(term.i3x2 + 0));
    KernelHessDiagAcc( PS, hessian, dvec, hdvec, term.i3x2, 1, term.i3x2, 1, fd_d2E_dii(term.i3x2 + 1));
    KernelHessDiagAcc( PS, hessian, dvec, hdvec, term.i3x2, 2, term.i3x2, 2, fd_d2E_dii(term.i3x2 + 2));
    // Off-diagonal: same set of (row,col) pairs that stretch.c::hessian_fd writes.
    KernelHessOffDiagAcc( PS, hessian, dvec, hdvec, term.i3x1, 1, term.i3x1, 0, fd_d2E_dij(term.i3x1+1, term.i3x1+0));
    KernelHessOffDiagAcc( PS, hessian, dvec, hdvec, term.i3x1, 2, term.i3x1, 0, fd_d2E_dij(term.i3x1+2, term.i3x1+0));
    KernelHessOffDiagAcc( PS, hessian, dvec, hdvec, term.i3x1, 2, term.i3x1, 1, fd_d2E_dij(term.i3x1+2, term.i3x1+1));
    KernelHessOffDiagAcc( PS, hessian, dvec, hdvec, term.i3x2, 0, term.i3x1, 0, fd_d2E_dij(term.i3x2+0, term.i3x1+0));
    KernelHessOffDiagAcc( PS, hessian, dvec, hdvec, term.i3x2, 0, term.i3x1, 1, fd_d2E_dij(term.i3x2+0, term.i3x1+1));
    KernelHessOffDiagAcc( PS, hessian, dvec, hdvec, term.i3x2, 0, term.i3x1, 2, fd_d2E_dij(term.i3x2+0, term.i3x1+2));
    KernelHessOffDiagAcc( PS, hessian, dvec, hdvec, term.i3x2, 1, term.i3x1, 0, fd_d2E_dij(term.i3x2+1, term.i3x1+0));
    KernelHessOffDiagAcc( PS, hessian, dvec, hdvec, term.i3x2, 1, term.i3x1, 1, fd_d2E_dij(term.i3x2+1, term.i3x1+1));
    KernelHessOffDiagAcc( PS, hessian, dvec, hdvec, term.i3x2, 1, term.i3x1, 2, fd_d2E_dij(term.i3x2+1, term.i3x1+2));
    KernelHessOffDiagAcc( PS, hessian, dvec, hdvec, term.i3x2, 1, term.i3x2, 0, fd_d2E_dij(term.i3x2+1, term.i3x2+0));
    KernelHessOffDiagAcc( PS, hessian, dvec, hdvec, term.i3x2, 2, term.i3x1, 0, fd_d2E_dij(term.i3x2+2, term.i3x1+0));
    KernelHessOffDiagAcc( PS, hessian, dvec, hdvec, term.i3x2, 2, term.i3x1, 1, fd_d2E_dij(term.i3x2+2, term.i3x1+1));
    KernelHessOffDiagAcc( PS, hessian, dvec, hdvec, term.i3x2, 2, term.i3x1, 2, fd_d2E_dij(term.i3x2+2, term.i3x1+2));
    KernelHessOffDiagAcc( PS, hessian, dvec, hdvec, term.i3x2, 2, term.i3x2, 0, fd_d2E_dij(term.i3x2+2, term.i3x2+0));
    KernelHessOffDiagAcc( PS, hessian, dvec, hdvec, term.i3x2, 2, term.i3x2, 1, fd_d2E_dij(term.i3x2+2, term.i3x2+1));
  }

};
