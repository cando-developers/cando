#define I1 i3x1
#define I2 i3x2
#define I3 i3x3
#define ANGLE_SET_PARAMETER(xx)
#define ANGLE_APPLY_ATOM_MASK(ii1,ii2,ii3)
#define ANGLE_SET_POSITION(xx,ii,oo) xx = position[ii+oo];
#define USE_EXPLICIT_DECLARES 1
#define power2(xx) (xx*xx)
#define mysqrt(xx) sqrt(xx)
#define reciprocal(xxx) (1.0/xxx)
#define ArcCos(xxx) acos(xxx)
#define DECLARE_FLOAT(xx) double xx;
#define ANGLE_ENERGY_ACCUMULATE(eee) *energy_accumulate = eee;
#define ANGLE_FORCE_ACCUMULATE(ii,oo,ff) force[ii+oo] = (-(ff));
#define ANGLE_DEBUG_INTERACTIONS(ii1,ii2,ii3)
#define ANGLE_CALC_FORCE 1
#define ANGLE_CALC_DIAGONAL_HESSIAN 1
#define ANGLE_CALC_OFF_DIAGONAL_HESSIAN 1
#define	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE 	KernelDiagHessAcc
#define	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE KernelOffDiagHessAcc
int calcForce = (force!=NULL);
int calcDiagonalHessian = (hdvec!=NULL);
int calcOffDiagonalHessian = (hdvec!=NULL);
  double x1;
  double y1;
  double z1;
  double x2;
  double y2;
  double z2;
  double x3;
  double y3;
  double z3;
int IllegalAngle;
#include "../energy_functions/_Angle_termDeclares.cc"
#include "../energy_functions/_Angle_termCode.cc"

#undef ANGLE_SET_PARAMETER
#undef ANGLE_APPLY_ATOM_MASK
#undef ANGLE_SET_POSITION
#undef DECLARE_FLOAT
#undef USE_EXPLICIT_DECLARES
#undef power2
#undef mysqrt
#undef reciprocal
#undef ANGLE_ENERGY_ACCUMULATE
#undef ANGLE_FORCE_ACCUMULATE
#undef ANGLE_DEBUG_INTERACTIONS
#undef ANGLE_CALC_FORCE
#undef ANGLE_CALC_DIAGONAL_HESSIAN
#undef ANGLE_CALC_OFF_DIAGONAL_HESSIAN
#undef ANGLE_DIAGONAL_HESSIAN_ACCUMULATE
#undef ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#undef I1
#undef I2
#undef I3
