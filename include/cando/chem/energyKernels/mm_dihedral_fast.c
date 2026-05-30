// File: mm_dihedral_fast.c
//
// Drop-in replacement for dihedral_fast.c. Exposes the same
// template <typename HESSIAN> struct Dihedral interface (PositionSize,
// energy/gradient/hessian and *_fd siblings, taking const dihedral_term&),
// but the *analytical* methods drive the legacy Mathematica-derived
// _Dihedral_termCode.cc via the DIHEDRAL_* macro family. The *_fd methods
// are central-difference references that call the local energy(), so the
// FD reference automatically matches the macro-based energy law here.
//
// Use as an A/B test against the new symbolic kernel: swap the include in
// energyDihedral.cc from
//   #include "cando/chem/energyKernels/dihedral_fast.c"
// to
//   #include "cando/chem/energyKernels/mm_dihedral_fast.c"
// rebuild, and re-run mmin with only the dihedral component enabled.
//
// Parameter-name remap: the new dihedral_term carries (v, n, sinphase,
// cosphase); the legacy macro file uses (V, DN, IN, sinPhase, cosPhase).
// We pull values out of the term struct into locals with the legacy names.

#include "dihedral_fast.h"

template <typename HESSIAN>
struct Dihedral {
  static constexpr size_t PositionSize = 12;
  static std::string description() { return "macrokernel-dihedral"; }

  double energy(const dihedral_term& term, double* position, double* energy_accumulate) {
    double energy_added = 0.0;
    // _Dihedral_termCode.cc expects: V, DN (double mult.), IN (int mult.),
    // sinPhase, cosPhase, I1..I4 in scope. Pull from term struct.
    double V        = term.v;
    double DN       = term.n;
    int    IN       = (int)term.n;
    double sinPhase = term.sinphase;
    double cosPhase = term.cosphase;
    int    I1       = term.i3x1;
    int    I2       = term.i3x2;
    int    I3       = term.i3x3;
    int    I4       = term.i3x4;
    #define USE_EXPLICIT_DECLARES 1
    #define DECLARE_FLOAT(xx) double xx;
    #define ZERO_SMALL_LEN(xxx,yyy)
    #define VEC_CONST(xxx) (xxx)
    #define MAX(xxx,yyy) fmax(xxx,yyy)
    #define MIN(xxx,yyy) fmin(xxx,yyy)
    #define DO_sinNPhiCosNPhi(II,SN,CN,SP,CP) sinNPhiCosNPhi(II,SN,CN,SP,CP)
    #define DIHEDRAL_SET_PARAMETER(xx)
    #define DIHEDRAL_APPLY_ATOM_MASK(ii1,ii2,ii3,ii4)
    #define DIHEDRAL_SET_POSITION(xx,ii,of) xx = position[ii+of];
    #define DIHEDRAL_ENERGY_ACCUMULATE(eee) { energy_added += (eee); *energy_accumulate += (eee); }
    #define DIHEDRAL_FORCE_ACCUMULATE(ii,oo,ff)
    #define DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(ii1,oo1,ii2,oo2,vv)
    #define DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(ii1,oo1,ii2,oo2,vv)
    #define DIHEDRAL_DEBUG_INTERACTIONS(ii1,ii2,ii3,ii4)
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
    // _Dihedral_termDeclares.cc does NOT declare x1..z4; declare them here.
    // y1 also shadows libm's Bessel y1(double).
    double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
    double SinNPhi, CosNPhi;
    // EraseLinearDihedral scales V in tzz1100 = EraseLinearDihedral * V.
    // Production sets it per-term; for the standard dihedral energy law
    // (no linear-term erasure) use 1.0.
    double EraseLinearDihedral = 1.0;
    #include <cando/chem/energy_functions/_Dihedral_termDeclares.cc>
    #include <cando/chem/energy_functions/_Dihedral_termCode.cc>
    #pragma GCC diagnostic pop
    #pragma clang diagnostic pop
    #undef USE_EXPLICIT_DECLARES
    #undef DECLARE_FLOAT
    #undef ZERO_SMALL_LEN
    #undef VEC_CONST
    #undef MAX
    #undef MIN
    #undef DO_sinNPhiCosNPhi
    #undef DIHEDRAL_SET_PARAMETER
    #undef DIHEDRAL_APPLY_ATOM_MASK
    #undef DIHEDRAL_SET_POSITION
    #undef DIHEDRAL_ENERGY_ACCUMULATE
    #undef DIHEDRAL_FORCE_ACCUMULATE
    #undef DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
    #undef DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
    #undef DIHEDRAL_DEBUG_INTERACTIONS
    return energy_added;
  }

  double gradient(const dihedral_term& term, double* position, double* energy_accumulate, double* force) {
    double energy_added = 0.0;
    double V        = term.v;
    double DN       = term.n;
    int    IN       = (int)term.n;
    double sinPhase = term.sinphase;
    double cosPhase = term.cosphase;
    int    I1       = term.i3x1;
    int    I2       = term.i3x2;
    int    I3       = term.i3x3;
    int    I4       = term.i3x4;
    #define USE_EXPLICIT_DECLARES 1
    #define DECLARE_FLOAT(xx) double xx;
    #define ZERO_SMALL_LEN(xxx,yyy)
    #define VEC_CONST(xxx) (xxx)
    #define MAX(xxx,yyy) fmax(xxx,yyy)
    #define MIN(xxx,yyy) fmin(xxx,yyy)
    #define DO_sinNPhiCosNPhi(II,SN,CN,SP,CP) sinNPhiCosNPhi(II,SN,CN,SP,CP)
    #define DIHEDRAL_SET_PARAMETER(xx)
    #define DIHEDRAL_APPLY_ATOM_MASK(ii1,ii2,ii3,ii4)
    #define DIHEDRAL_SET_POSITION(xx,ii,of) xx = position[ii+of];
    #define DIHEDRAL_ENERGY_ACCUMULATE(eee) { energy_added += (eee); *energy_accumulate += (eee); }
    #define DIHEDRAL_FORCE_ACCUMULATE(ii,oo,ff) force[ii+oo] += (ff);
    #define DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(ii1,oo1,ii2,oo2,vv)
    #define DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(ii1,oo1,ii2,oo2,vv)
    #define DIHEDRAL_DEBUG_INTERACTIONS(ii1,ii2,ii3,ii4)
    #define DIHEDRAL_CALC_FORCE 1
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
    double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
    double SinNPhi, CosNPhi;
    double EraseLinearDihedral = 1.0;
    #include <cando/chem/energy_functions/_Dihedral_termDeclares.cc>
    #include <cando/chem/energy_functions/_Dihedral_termCode.cc>
    #pragma GCC diagnostic pop
    #pragma clang diagnostic pop
    #undef USE_EXPLICIT_DECLARES
    #undef DECLARE_FLOAT
    #undef ZERO_SMALL_LEN
    #undef VEC_CONST
    #undef MAX
    #undef MIN
    #undef DO_sinNPhiCosNPhi
    #undef DIHEDRAL_SET_PARAMETER
    #undef DIHEDRAL_APPLY_ATOM_MASK
    #undef DIHEDRAL_SET_POSITION
    #undef DIHEDRAL_ENERGY_ACCUMULATE
    #undef DIHEDRAL_FORCE_ACCUMULATE
    #undef DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
    #undef DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
    #undef DIHEDRAL_DEBUG_INTERACTIONS
    #undef DIHEDRAL_CALC_FORCE
    return energy_added;
  }

  double hessian(const dihedral_term& term, double* position, double* energy_accumulate, double* force,
                 HESSIAN hessian, double* dvec, double* hdvec) {
    double energy_added = 0.0;
    double V        = term.v;
    double DN       = term.n;
    int    IN       = (int)term.n;
    double sinPhase = term.sinphase;
    double cosPhase = term.cosphase;
    int    I1       = term.i3x1;
    int    I2       = term.i3x2;
    int    I3       = term.i3x3;
    int    I4       = term.i3x4;
    #define USE_EXPLICIT_DECLARES 1
    #define DECLARE_FLOAT(xx) double xx;
    #define ZERO_SMALL_LEN(xxx,yyy)
    #define VEC_CONST(xxx) (xxx)
    #define MAX(xxx,yyy) fmax(xxx,yyy)
    #define MIN(xxx,yyy) fmin(xxx,yyy)
    #define DO_sinNPhiCosNPhi(II,SN,CN,SP,CP) sinNPhiCosNPhi(II,SN,CN,SP,CP)
    #define DIHEDRAL_SET_PARAMETER(xx)
    #define DIHEDRAL_APPLY_ATOM_MASK(ii1,ii2,ii3,ii4)
    #define DIHEDRAL_SET_POSITION(xx,ii,of) xx = position[ii+of];
    #define DIHEDRAL_ENERGY_ACCUMULATE(eee) { energy_added += (eee); *energy_accumulate += (eee); }
    #define DIHEDRAL_FORCE_ACCUMULATE(ii,oo,ff) force[ii+oo] += (ff);
    #define DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(ii1,oo1,ii2,oo2,vv) \
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, ii1, oo1, ii2, oo2, (vv) )
    #define DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(ii1,oo1,ii2,oo2,vv) \
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, ii1, oo1, ii2, oo2, (vv) )
    #define DIHEDRAL_DEBUG_INTERACTIONS(ii1,ii2,ii3,ii4)
    #define DIHEDRAL_CALC_FORCE 1
    #define DIHEDRAL_CALC_DIAGONAL_HESSIAN 1
    #define DIHEDRAL_CALC_OFF_DIAGONAL_HESSIAN 1
    int calcForce = 1;
    int calcDiagonalHessian = 1;
    int calcOffDiagonalHessian = 1;
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-variable"
    #pragma clang diagnostic ignored "-Wunused-label"
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic ignored "-Wunused-label"
    double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
    double SinNPhi, CosNPhi;
    double EraseLinearDihedral = 1.0;
    #include <cando/chem/energy_functions/_Dihedral_termDeclares.cc>
    #include <cando/chem/energy_functions/_Dihedral_termCode.cc>
    #pragma GCC diagnostic pop
    #pragma clang diagnostic pop
    #undef USE_EXPLICIT_DECLARES
    #undef DECLARE_FLOAT
    #undef ZERO_SMALL_LEN
    #undef VEC_CONST
    #undef MAX
    #undef MIN
    #undef DO_sinNPhiCosNPhi
    #undef DIHEDRAL_SET_PARAMETER
    #undef DIHEDRAL_APPLY_ATOM_MASK
    #undef DIHEDRAL_SET_POSITION
    #undef DIHEDRAL_ENERGY_ACCUMULATE
    #undef DIHEDRAL_FORCE_ACCUMULATE
    #undef DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
    #undef DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
    #undef DIHEDRAL_DEBUG_INTERACTIONS
    #undef DIHEDRAL_CALC_FORCE
    #undef DIHEDRAL_CALC_DIAGONAL_HESSIAN
    #undef DIHEDRAL_CALC_OFF_DIAGONAL_HESSIAN
    return energy_added;
  }

  // ---------------------------------------------------------------------
  // Finite-difference reference implementations
  // ---------------------------------------------------------------------

  void energy_fd(const dihedral_term& term, double* position, double* energy_accumulate) {
    energy(term, position, energy_accumulate);
  }

  void gradient_fd(const dihedral_term& term, double* position, double* energy_accumulate, double* force) {
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
                        (size_t)term.i3x3, (size_t)term.i3x3, (size_t)term.i3x3,
                        (size_t)term.i3x4, (size_t)term.i3x4, (size_t)term.i3x4 };
    size_t off [PS] = { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2 };
    for (size_t r = 0; r < PS; ++r) {
      KernelGradientAcc(atom[r], off[r], fd_dE_d(atom[r] + off[r]));
    }
  }

  void hessian_fd(const dihedral_term& term, double* position, double* energy_accumulate, double* force,
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
                        (size_t)term.i3x3, (size_t)term.i3x3, (size_t)term.i3x3,
                        (size_t)term.i3x4, (size_t)term.i3x4, (size_t)term.i3x4 };
    size_t off [PS] = { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2 };
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
