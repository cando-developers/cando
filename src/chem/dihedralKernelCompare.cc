/*
    File: dihedralKernelCompare.cc

    Test-only translation unit that compiles BOTH dihedral kernels
    (dihedral_fast and dihedral_fast_floored) into one binary and compares their
    energy and gradient term-by-term.  This lets a single build validate that the
    trig-free Blondel-Karplus kernel agrees with the AD kernel -- no recompile to
    swap, which the production A/B selector (CANDO_DIHEDRAL_KERNEL in
    energyDihedral.cc) would otherwise require.

    Both generated kernel headers define an identical `struct dihedral_term`, so
    each kernel is #included in its own namespace to avoid a duplicate definition.
    The per-namespace `term_type` alias resolves to whichever struct that kernel's
    functor actually uses, which keeps this robust to CANDO_DIHEDRAL_KERNEL.
*/

#include <cando/main/extension.h>        // DOUBLE (vecreal), NVECTOR_*_TYPE
#include <cando/chem/energyDihedral.h>   // EnergyDihedral_O, KernelGradientAcc (via energyComponent.h)
#include <cando/chem/nVector.h>
#include <cando/chem/energyKernels/dihedralRecurrence.h> // cos_n_phi / sin_n_phi (floored kernel)
#include <clasp/core/wrappers.h>
#include <cmath> // std::isfinite

namespace chem {

namespace dih_fast {
#include "cando/chem/energyKernels/dihedral_fast.c"
using term_type = dihedral_term;
} // namespace dih_fast

namespace dih_floored {
#include "cando/chem/energyKernels/dihedral_fast_floored.c"
using term_type = dihedral_term;
} // namespace dih_floored

};

namespace chem {
// Compare dihedral_fast vs dihedral_fast_floored over every term of COMP at the
// coordinate vector COORDS.  Returns two values: the maximum absolute energy
// difference and the maximum absolute force-component difference across all
// dihedral terms.  Both kernels live in this one binary, so no rebuild/swap is
// needed to A/B them.
CL_DEFUN core::T_mv chem__compare_dihedral_kernels(chem::EnergyDihedral_sp comp, chem::NVector_sp coords) {
  dih_fast::Dihedral_Fast<double *> ka;
  dih_floored::Dihedral_Fast_Floored<double *> kb;
  double pos[12], fa[12], fb[12];
  double maxdE = 0.0;
  double maxdF = 0.0;
  for (auto di = comp->_Terms.begin(); di != comp->_Terms.end(); ++di) {
    size_t base[4] = {(size_t)di->term.i3x1, (size_t)di->term.i3x2, (size_t)di->term.i3x3, (size_t)di->term.i3x4};
    for (int a = 0; a < 4; ++a)
      for (int c = 0; c < 3; ++c)
        pos[a * 3 + c] = coords[base[a] + c];
    for (int i = 0; i < 12; ++i) {
      fa[i] = 0.0;
      fb[i] = 0.0;
    }
    dih_fast::term_type ta(di->term.v, di->term.n, di->term.sinphase, di->term.cosphase, 0, 3, 6, 9);
    dih_floored::term_type tb(di->term.v, di->term.n, di->term.sinphase, di->term.cosphase, 0, 3, 6, 9);
    double ea = 0.0, eb = 0.0;
    ka.gradient(ta, pos, &ea, fa);
    kb.gradient(tb, pos, &eb, fb);
    double de = fabs(ea - eb);
    if (de > maxdE)
      maxdE = de;
    for (int i = 0; i < 12; ++i) {
      double df = fabs(fa[i] - fb[i]);
      if (df > maxdF)
        maxdF = df;
    }
  }
  return Values(core::clasp_make_double_float(maxdE), core::clasp_make_double_float(maxdF));
}

// Evaluate both kernels on a single (near-)linear dihedral and return the max
// absolute force component of each: (max|force_fast|, max|force_floored|).
//
// Atoms 1-2-3 lie on a line offset only by COLINEAR_EPS in y, so the plane
// normal c1 = (p2-p1)x(p3-p2) -> 0.  With COLINEAR_EPS = 0 they are exactly
// collinear: dihedral_fast's 1/(t1^2+t2^2) factor becomes 1/0 -> non-finite,
// while dihedral_fast_floored stays finite (the floored denominators, and the
// c1=0 factor that multiplies them, give a clean bounded force).  Demonstrates
// the linear-dihedral robustness the floored kernel was built for.
//
// Call inside ext:with-float-traps-masked -- the fast kernel divides by zero.
CL_DEFUN core::T_mv chem__dihedral_kernel_near_linear(double colinear_eps) {
  dih_fast::Dihedral_Fast<double *> ka;
  dih_floored::Dihedral_Fast_Floored<double *> kb;
  double pos[12] = {0.0, 0.0,          0.0,   // atom 1
                    1.0, 0.0,          0.0,   // atom 2
                    2.0, colinear_eps, 0.0,   // atom 3 (collinear with 1,2 when eps=0)
                    3.0, 0.0,          1.0};  // atom 4 (out of the 1-2-3 line)
  // V=10, n=2, phase=0  ->  sinPhase=0, cosPhase=1
  dih_fast::term_type ta(10.0, 2.0, 0.0, 1.0, 0, 3, 6, 9);
  dih_floored::term_type tb(10.0, 2.0, 0.0, 1.0, 0, 3, 6, 9);
  double fa[12], fb[12];
  for (int i = 0; i < 12; ++i) {
    fa[i] = 0.0;
    fb[i] = 0.0;
  }
  double ea = 0.0, eb = 0.0;
  ka.gradient(ta, pos, &ea, fa);
  kb.gradient(tb, pos, &eb, fb);
  double maxfa = 0.0, maxfb = 0.0;
  for (int i = 0; i < 12; ++i) {
    // propagate non-finite (inf/nan) rather than letting a NaN comparison drop it
    if (!std::isfinite(fa[i]))
      maxfa = INFINITY;
    else if (std::fabs(fa[i]) > maxfa)
      maxfa = std::fabs(fa[i]);
    if (!std::isfinite(fb[i]))
      maxfb = INFINITY;
    else if (std::fabs(fb[i]) > maxfb)
      maxfb = std::fabs(fb[i]);
  }
  return Values(core::clasp_make_double_float(maxfa), core::clasp_make_double_float(maxfb));
}

} // namespace chem
