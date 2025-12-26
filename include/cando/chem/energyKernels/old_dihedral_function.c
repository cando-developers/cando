#define IN n
#define DN n
#define sinPhase sinphase
#define cosPhase cosphase
#define V v
#define I1 i3x1
#define I2 i3x2
#define I3 i3x3
#define I4 i3x4
#define DIHEDRAL_SET_PARAMETER(xx)
#define DIHEDRAL_APPLY_ATOM_MASK(ii1,ii2,ii3,ii4)
#define DIHEDRAL_SET_POSITION(xx,ii,oo) xx = position[ii+oo];
#define USE_EXPLICIT_DECLARES 1
#define power2(xx) (xx*xx)
#define mysqrt(xx) sqrt(xx)
#define reciprocal(xxx) (1.0/xxx)
#define DECLARE_FLOAT(xx) double xx;
#define DIHEDRAL_ENERGY_ACCUMULATE(eee) *energy_accumulate = eee;
#define DIHEDRAL_FORCE_ACCUMULATE(ii,oo,ff) force[ii+oo] = (-(ff));
#define DIHEDRAL_DEBUG_INTERACTIONS(ii1,ii2,ii3,ii4)
#define DIHEDRAL_CALC_FORCE 1
#define DIHEDRAL_CALC_DIAGONAL_HESSIAN 1
#define DIHEDRAL_CALC_OFF_DIAGONAL_HESSIAN 1
#define	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE 	KernelDiagHessAcc
#define	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE KernelOffDiagHessAcc
#define ZERO_SMALL_LEN(xxx,yyy)
#define VEC_CONST(xxx) (xxx)
#define MAX(xxx,yyy) fmax(xxx,yyy)
#define MIN(xxx,yyy) fmin(xxx,yyy)
#define DO_sinNPhiCosNPhi(IN,SinNPhi,CosNPhi,SinPhi,CosPhi) sinNPhiCosNPhi(IN,&SinNPhi,&CosNPhi,SinPhi,CosPhi)

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
  double x4;
  double y4;
  double z4;
double SinNPhi;
double CosNPhi;
size_t EraseLinearDihedral;

#include "../energy_functions/_Dihedral_termDeclares.cc"
#include "../energy_functions/_Dihedral_termCode.cc"

#undef DIHEDRAL_SET_PARAMETER
#undef DIHEDRAL_APPLY_ATOM_MASK
#undef DIHEDRAL_SET_POSITION
#undef DECLARE_FLOAT
#undef USE_EXPLICIT_DECLARES
#undef power2
#undef mysqrt
#undef reciprocal
#undef DIHEDRAL_ENERGY_ACCUMULATE
#undef DIHEDRAL_FORCE_ACCUMULATE
#undef DIHEDRAL_DEBUG_INTERACTIONS
#undef DIHEDRAL_CALC_FORCE
#undef DIHEDRAL_CALC_DIAGONAL_HESSIAN
#undef DIHEDRAL_CALC_OFF_DIAGONAL_HESSIAN
#undef DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE
#undef DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#undef I1
#undef I2
#undef I3
#undef I4
#undef ZERO_SMALL_LEN

#undef VEC_CONST
#undef MAX
#undef MIN

#undef DO_sinNPhiCosNPhi
#undef IN
#undef DN
#undef sinPhase
#undef cosPhase
#undef V
