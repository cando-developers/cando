{Name -> "Oozp", AdditionalCDeclares -> "", Input -> {x1, y1, z1, kb, za}, 
 Output -> {Energy, fx1, fy1, fz1, dhx1x1, dhy1y1, dhz1z1, ohx1y1, ohx1z1, 
   ohy1z1, OozpDeviation}, EnergyFunction -> EnergyFunction, 
 DerivativeVariables -> {x1, y1, z1}, 
 Rules -> {CCode["OOZP_SET_PARAMETER(kb);"], 
   CCode["OOZP_SET_PARAMETER(za);"], CCode["OOZP_SET_PARAMETER(I1);"], 
   CCode["OOZP_SET_PARAMETER(I2);"], CCode["OOZP_SET_POSITION(x1,I1,0);"], 
   CCode["OOZP_SET_POSITION(y1,I1,1);"], 
   CCode["OOZP_SET_POSITION(z1,I1,2);"], -za -> tx5, tx5 + z1 -> tx6, 
   power2[tx6] -> tx7, kb*tx7 -> Energy, 
   CCode["OOZP_ENERGY_ACCUMULATE(Energy);"], 
   CCode["#ifdef OOZP_CALC_FORCE //["], CCode["if ( calcForce ) {"], 
   0 -> gx1, -gx1 -> fx1, CCode["OOZP_FORCE_ACCUMULATE(I1, 0, fx1 );"], 
   gx1 -> gy1, -gy1 -> fy1, CCode["OOZP_FORCE_ACCUMULATE(I1, 1, fy1 );"], 
   tx6 -> tx8, 2.*kb*tx8 -> gz1, -gz1 -> fz1, 
   CCode["OOZP_FORCE_ACCUMULATE(I1, 2, fz1 );"], 
   CCode["#ifdef OOZP_CALC_DIAGONAL_HESSIAN //["], 
   CCode["if ( calcDiagonalHessian ) {"], gy1 -> dhx1x1, 
   CCode["OOZP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);"], 
   dhx1x1 -> dhy1y1, 
   CCode["OOZP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);"], 
   2.*kb -> dhz1z1, 
   CCode["OOZP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);"], 
   CCode["#ifdef OOZP_CALC_OFF_DIAGONAL_HESSIAN //["], 
   CCode["if ( calcOffDiagonalHessian ) {"], dhy1y1 -> ohx1y1, 
   CCode["OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);"], 
   ohx1y1 -> ohx1z1, 
   CCode["OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);"], 
   ohx1z1 -> ohy1z1, 
   CCode["OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);"], 
   CCode["} /*if calcOffDiagonalHessian */ "], 
   CCode["#endif /* OOZP_CALC_OFF_DIAGONAL_HESSIAN ]*/"], 
   CCode["} /*calcDiagonalHessian */"], 
   CCode["#endif /* OOZP_CALC_DIAGONAL_HESSIAN ]*/"], 
   CCode["} /*calcForce */"], CCode["#endif /* OOZP_CALC_FORCE ]*/"]}, 
 HessianStructure -> {{5, 8, 9}, {8, 6, 10}, {9, 10, 7}}, 
 CCode -> "NOT USED", MathCode -> "Block[{dhx1x1,dhy1y1,dhz1z1,Energy,fx1,fy1\
,fz1,gx1,gy1,gz1,ohx1y1,ohx1z1,ohy1z1,OozpDeviation,tx5,tx6,tx7,tx8,xxxDummy}\
,\n\t(*OOZP_SET_PARAMETER(kb);*)\n\t(*OOZP_SET_PARAMETER(za);*)\n\t(*OOZP_SET\
_PARAMETER(I1);*)\n\t(*OOZP_SET_PARAMETER(I2);*)\n\t(*OOZP_SET_POSITION(x1,I1\
,0);*)\n\t(*OOZP_SET_POSITION(y1,I1,1);*)\n\t(*OOZP_SET_POSITION(z1,I1,2);*)\
\n\ttx5=-za; (*Rule 8*)\n\ttx6=tx5 + z1; (*Rule 9*)\n\ttx7= mathPower2[tx6]; \
(*Rule 10*)\n\tEnergy=kb tx7; (*Rule \
11*)\n\t(*OOZP_ENERGY_ACCUMULATE(Energy);*)\n\t(*#ifdef OOZP_CALC_FORCE \
//[*)\n\t(*if ( calcForce ) {*)\n\tgx1=0; (*Rule 15*)\n\tfx1=-gx1; (*Rule \
16*)\n\t(*OOZP_FORCE_ACCUMULATE(I1, 0, fx1 );*)\n\tgy1=gx1; (*Rule \
18*)\n\tfy1=-gy1; (*Rule 19*)\n\t(*OOZP_FORCE_ACCUMULATE(I1, 1, fy1 \
);*)\n\ttx8=tx6; (*Rule 21*)\n\tgz1=2. kb tx8; (*Rule 22*)\n\tfz1=-gz1; \
(*Rule 23*)\n\t(*OOZP_FORCE_ACCUMULATE(I1, 2, fz1 );*)\n\t(*#ifdef \
OOZP_CALC_DIAGONAL_HESSIAN //[*)\n\t(*if ( calcDiagonalHessian ) \
{*)\n\tdhx1x1=gy1; (*Rule 27*)\n\t(*OOZP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, \
I1, 0, dhx1x1);*)\n\tdhy1y1=dhx1x1; (*Rule \
29*)\n\t(*OOZP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, \
dhy1y1);*)\n\tdhz1z1=2. kb; (*Rule \
31*)\n\t(*OOZP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, \
dhz1z1);*)\n\t(*#ifdef OOZP_CALC_OFF_DIAGONAL_HESSIAN //[*)\n\t(*if ( \
calcOffDiagonalHessian ) {*)\n\tohx1y1=dhy1y1; (*Rule \
35*)\n\t(*OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, \
ohx1y1);*)\n\tohx1z1=ohx1y1; (*Rule \
37*)\n\t(*OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, \
ohx1z1);*)\n\tohy1z1=ohx1z1; (*Rule \
39*)\n\t(*OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, \
ohy1z1);*)\n\t(*} /*if calcOffDiagonalHessian */ *)\n\t(*#endif /* \
OOZP_CALC_OFF_DIAGONAL_HESSIAN ]*/*)\n\t(*} /*calcDiagonalHessian \
*/*)\n\t(*#endif /* OOZP_CALC_DIAGONAL_HESSIAN ]*/*)\n\t(*} /*calcForce \
*/*)\n\t(*#endif /* OOZP_CALC_FORCE ]*/*)\n{Energy, fx1, fy1, fz1, dhx1x1, \
dhy1y1, dhz1z1, ohx1y1, ohx1z1, ohy1z1, OozpDeviation}]\n"}
