
#define I1 i3x1
#define I2 i3x2
#define STRETCH_SET_PARAMETER(xx)
#define STRETCH_APPLY_ATOM_MASK(ii1,ii2)
#define STRETCH_SET_POSITION(xx,ii,oo) xx = position[ii+oo];
#define USE_EXPLICIT_DECLARES 1
#define power2(xx) (xx*xx)
#define mysqrt(xx) sqrt(xx)
#define reciprocal(xxx) (1.0/xxx)
#define DECLARE_FLOAT(xx) double xx;
#define STRETCH_ENERGY_ACCUMULATE(eee) *energy_accumulate = eee;
#define STRETCH_FORCE_ACCUMULATE(ii,oo,ff) force[ii+oo] = (-(ff));
#define STRETCH_DEBUG_INTERACTIONS(ii1,ii2)
#define STRETCH_CALC_FORCE 1
#define STRETCH_CALC_DIAGONAL_HESSIAN 1
#define STRETCH_CALC_OFF_DIAGONAL_HESSIAN 1
#define	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE 	KernelDiagHessAcc
#define	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE KernelOffDiagHessAcc
int calcForce = (force!=NULL);
int calcDiagonalHessian = (hdvec!=NULL);
int calcOffDiagonalHessian = (hdvec!=NULL);
  double x1;
  double y1;
  double z1;
  double x2;
  double y2;
  double z2;
#include "../energy_functions/_Stretch_termDeclares.cc"
#include "../energy_functions/_Stretch_termCode.cc"

#undef STRETCH_SET_PARAMETER
#undef STRETCH_APPLY_ATOM_MASK
#undef STRETCH_SET_POSITION
#undef DECLARE_FLOAT
#undef USE_EXPLICIT_DECLARES
#undef power2
#undef mysqrt
#undef reciprocal
#undef STRETCH_ENERGY_ACCUMULATE
#undef STRETCH_FORCE_ACCUMULATE
#undef STRETCH_DEBUG_INTERACTIONS
#undef STRETCH_CALC_FORCE
#undef STRETCH_CALC_DIAGONAL_HESSIAN
#undef STRETCH_CALC_OFF_DIAGONAL_HESSIAN
#undef	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE
#undef	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#undef I1
#undef I2
