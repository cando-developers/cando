{Name -> "Erep", AdditionalCDeclares -> "", 
 Input -> {x1, y1, z1, x2, y2, z2, crep}, 
 Output -> {Erep, fx1, fy1, fz1, fx2, fy2, fz2, dhx1x1, dhy1y1, dhz1z1, 
   dhx2x2, dhy2y2, dhz2z2, ohx1y1, ohx1z1, ohx1x2, ohx1y2, ohx1z2, ohy1z1, 
   ohy1x2, ohy1y2, ohy1z2, ohz1x2, ohz1y2, ohz1z2, ohx2y2, ohx2z2, ohy2z2, 
   ErepDistance}, EnergyFunction -> 
  crep*Log[Sqrt[(x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2]], 
 DerivativeVariables -> {x1, y1, z1, x2, y2, z2}, 
 Rules -> {CCode["EREP_SET_PARAMETER(crep);"], 
   CCode["EREP_SET_PARAMETER(I1);"], CCode["EREP_SET_PARAMETER(I2);"], 
   CCode["EREP_APPLY_ATOM_MASK(I1,I2);"], 
   CCode["EREP_SET_POSITION(x1,I1,0);"], 
   CCode["EREP_SET_POSITION(y1,I1,1);"], 
   CCode["EREP_SET_POSITION(z1,I1,2);"], 
   CCode["EREP_SET_POSITION(x2,I2,0);"], 
   CCode["EREP_SET_POSITION(y2,I2,1);"], 
   CCode["EREP_SET_POSITION(z2,I2,2);"], -x2 -> tx1, -y2 -> tx2, -z2 -> tx3, 
   tx1 + x1 -> tx4, tx2 + y1 -> tx5, tx3 + z1 -> tx6, power2[tx4] -> tx7, 
   power2[tx5] -> tx8, power2[tx6] -> tx9, tx7 + tx8 + tx9 -> tx10, 
   mysqrt[tx10] -> ErepDistance, CCode["MAYBE_BAIL(ErepDistance);"], 
   Log[ErepDistance] -> tx11, crep*tx11 -> Erep, 
   CCode["EREP_ENERGY_ACCUMULATE(Erep);"], 
   CCode["#ifdef EREP_CALC_FORCE //["], CCode["if ( calcForce ) {"], 
   reciprocal[tx10] -> tx12, crep*tx12*tx4 -> gx1, -gx1 -> fx1, 
   CCode["EREP_FORCE_ACCUMULATE(I1, 0, fx1 );"], crep*tx12*tx5 -> gy1, 
   -gy1 -> fy1, CCode["EREP_FORCE_ACCUMULATE(I1, 1, fy1 );"], 
   crep*tx12*tx6 -> gz1, -gz1 -> fz1, 
   CCode["EREP_FORCE_ACCUMULATE(I1, 2, fz1 );"], fx1 -> gx2, -gx2 -> fx2, 
   CCode["EREP_FORCE_ACCUMULATE(I2, 0, fx2 );"], fy1 -> gy2, -gy2 -> fy2, 
   CCode["EREP_FORCE_ACCUMULATE(I2, 1, fy2 );"], fz1 -> gz2, -gz2 -> fz2, 
   CCode["EREP_FORCE_ACCUMULATE(I2, 2, fz2 );"], 
   CCode["#ifdef EREP_CALC_DIAGONAL_HESSIAN //["], 
   CCode["if ( calcDiagonalHessian ) {"], power2[tx12] -> tx13, 
   crep*tx12 -> tx14, -2.*crep*tx13*tx7 -> tx15, tx14 + tx15 -> dhx1x1, 
   CCode["EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);"], 
   -2.*crep*tx13*tx8 -> tx16, tx14 + tx16 -> dhy1y1, 
   CCode["EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);"], 
   -2.*crep*tx13*tx9 -> tx17, tx14 + tx17 -> dhz1z1, 
   CCode["EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);"], 
   dhx1x1 -> dhx2x2, 
   CCode["EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);"], 
   dhy1y1 -> dhy2y2, 
   CCode["EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);"], 
   dhz1z1 -> dhz2z2, 
   CCode["EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);"], 
   CCode["#ifdef EREP_CALC_OFF_DIAGONAL_HESSIAN //["], 
   CCode["if ( calcOffDiagonalHessian ) {"], -2.*crep*tx13*tx4*tx5 -> ohx1y1, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);"], 
   -2.*crep*tx13*tx4*tx6 -> ohx1z1, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);"], 
   -tx14 -> tx18, 2.*crep*tx13*tx7 -> tx19, tx18 + tx19 -> ohx1x2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);"], 
   2.*crep*tx13*tx4*tx5 -> ohx1y2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);"], 
   2.*crep*tx13*tx4*tx6 -> ohx1z2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);"], 
   -2.*crep*tx13*tx5*tx6 -> ohy1z1, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);"], 
   ohx1y2 -> ohy1x2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);"], 
   2.*crep*tx13*tx8 -> tx20, tx18 + tx20 -> ohy1y2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);"], 
   2.*crep*tx13*tx5*tx6 -> ohy1z2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);"], 
   ohx1z2 -> ohz1x2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);"], 
   ohy1z2 -> ohz1y2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);"], 
   2.*crep*tx13*tx9 -> tx21, tx18 + tx21 -> ohz1z2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);"], 
   ohx1y1 -> ohx2y2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);"], 
   ohx1z1 -> ohx2z2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);"], 
   ohy1z1 -> ohy2z2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);"], 
   CCode["} /*if calcOffDiagonalHessian */ "], 
   CCode["#endif /* EREP_CALC_OFF_DIAGONAL_HESSIAN ]*/"], 
   CCode["} /*calcDiagonalHessian */"], 
   CCode["#endif /* EREP_CALC_DIAGONAL_HESSIAN ]*/"], 
   CCode["} /*calcForce */"], CCode["#endif /* EREP_CALC_FORCE ]*/"], 
   CCode["SKIP_term:  (void)0;"]}, HessianStructure -> 
  {{8, 14, 15, 16, 17, 18}, {14, 9, 19, 20, 21, 22}, 
   {15, 19, 10, 23, 24, 25}, {16, 20, 23, 11, 26, 27}, 
   {17, 21, 24, 26, 12, 28}, {18, 22, 25, 27, 28, 13}}, CCode -> "NOT USED", 
 MathCode -> "Block[{dhx1x1,dhx2x2,dhy1y1,dhy2y2,dhz1z1,dhz2z2,Erep,ErepDista\
nce,fx1,fx2,fy1,fy2,fz1,fz2,gx1,gx2,gy1,gy2,gz1,gz2,ohx1x2,ohx1y1,ohx1y2,ohx1\
z1,ohx1z2,ohx2y2,ohx2z2,ohy1x2,ohy1y2,ohy1z1,ohy1z2,ohy2z2,ohz1x2,ohz1y2,ohz1\
z2,tx1,tx10,tx11,tx12,tx13,tx14,tx15,tx16,tx17,tx18,tx19,tx2,tx20,tx21,tx3,tx\
4,tx5,tx6,tx7,tx8,tx9,xxxDummy},\n\t(*EREP_SET_PARAMETER(crep);*)\n\t(*EREP_S\
ET_PARAMETER(I1);*)\n\t(*EREP_SET_PARAMETER(I2);*)\n\t(*EREP_APPLY_ATOM_MASK(\
I1,I2);*)\n\t(*EREP_SET_POSITION(x1,I1,0);*)\n\t(*EREP_SET_POSITION(y1,I1,1);\
*)\n\t(*EREP_SET_POSITION(z1,I1,2);*)\n\t(*EREP_SET_POSITION(x2,I2,0);*)\n\t(\
*EREP_SET_POSITION(y2,I2,1);*)\n\t(*EREP_SET_POSITION(z2,I2,2);*)\n\ttx1=-x2; \
(*Rule 11*)\n\ttx2=-y2; (*Rule 12*)\n\ttx3=-z2; (*Rule 13*)\n\ttx4=tx1 + x1; \
(*Rule 14*)\n\ttx5=tx2 + y1; (*Rule 15*)\n\ttx6=tx3 + z1; (*Rule 16*)\n\ttx7= \
mathPower2[tx4]; (*Rule 17*)\n\ttx8= mathPower2[tx5]; (*Rule 18*)\n\ttx9= \
mathPower2[tx6]; (*Rule 19*)\n\ttx10=tx7 + tx8 + tx9; (*Rule \
20*)\n\tErepDistance= mathSqrt[tx10]; (*Rule \
21*)\n\t(*MAYBE_BAIL(ErepDistance);*)\n\ttx11=Log[ErepDistance]; (*Rule \
23*)\n\tErep=crep tx11; (*Rule \
24*)\n\t(*EREP_ENERGY_ACCUMULATE(Erep);*)\n\t(*#ifdef EREP_CALC_FORCE \
//[*)\n\t(*if ( calcForce ) {*)\n\ttx12= mathReciprocal[tx10]; (*Rule \
28*)\n\tgx1=crep tx12 tx4; (*Rule 29*)\n\tfx1=-gx1; (*Rule \
30*)\n\t(*EREP_FORCE_ACCUMULATE(I1, 0, fx1 );*)\n\tgy1=crep tx12 tx5; (*Rule \
32*)\n\tfy1=-gy1; (*Rule 33*)\n\t(*EREP_FORCE_ACCUMULATE(I1, 1, fy1 \
);*)\n\tgz1=crep tx12 tx6; (*Rule 35*)\n\tfz1=-gz1; (*Rule \
36*)\n\t(*EREP_FORCE_ACCUMULATE(I1, 2, fz1 );*)\n\tgx2=fx1; (*Rule \
38*)\n\tfx2=-gx2; (*Rule 39*)\n\t(*EREP_FORCE_ACCUMULATE(I2, 0, fx2 \
);*)\n\tgy2=fy1; (*Rule 41*)\n\tfy2=-gy2; (*Rule \
42*)\n\t(*EREP_FORCE_ACCUMULATE(I2, 1, fy2 );*)\n\tgz2=fz1; (*Rule \
44*)\n\tfz2=-gz2; (*Rule 45*)\n\t(*EREP_FORCE_ACCUMULATE(I2, 2, fz2 \
);*)\n\t(*#ifdef EREP_CALC_DIAGONAL_HESSIAN //[*)\n\t(*if ( \
calcDiagonalHessian ) {*)\n\ttx13= mathPower2[tx12]; (*Rule 49*)\n\ttx14=crep \
tx12; (*Rule 50*)\n\ttx15=-2. crep tx13 tx7; (*Rule 51*)\n\tdhx1x1=tx14 + \
tx15; (*Rule 52*)\n\t(*EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, \
dhx1x1);*)\n\ttx16=-2. crep tx13 tx8; (*Rule 54*)\n\tdhy1y1=tx14 + tx16; \
(*Rule 55*)\n\t(*EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, \
dhy1y1);*)\n\ttx17=-2. crep tx13 tx9; (*Rule 57*)\n\tdhz1z1=tx14 + tx17; \
(*Rule 58*)\n\t(*EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, \
dhz1z1);*)\n\tdhx2x2=dhx1x1; (*Rule \
60*)\n\t(*EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, \
dhx2x2);*)\n\tdhy2y2=dhy1y1; (*Rule \
62*)\n\t(*EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, \
dhy2y2);*)\n\tdhz2z2=dhz1z1; (*Rule \
64*)\n\t(*EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, \
dhz2z2);*)\n\t(*#ifdef EREP_CALC_OFF_DIAGONAL_HESSIAN //[*)\n\t(*if ( \
calcOffDiagonalHessian ) {*)\n\tohx1y1=-2. crep tx13 tx4 tx5; (*Rule \
68*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, \
ohx1y1);*)\n\tohx1z1=-2. crep tx13 tx4 tx6; (*Rule \
70*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, \
ohx1z1);*)\n\ttx18=-tx14; (*Rule 72*)\n\ttx19=2. crep tx13 tx7; (*Rule \
73*)\n\tohx1x2=tx18 + tx19; (*Rule \
74*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, \
ohx1x2);*)\n\tohx1y2=2. crep tx13 tx4 tx5; (*Rule \
76*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, \
ohx1y2);*)\n\tohx1z2=2. crep tx13 tx4 tx6; (*Rule \
78*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, \
ohx1z2);*)\n\tohy1z1=-2. crep tx13 tx5 tx6; (*Rule \
80*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, \
ohy1z1);*)\n\tohy1x2=ohx1y2; (*Rule \
82*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, \
ohy1x2);*)\n\ttx20=2. crep tx13 tx8; (*Rule 84*)\n\tohy1y2=tx18 + tx20; \
(*Rule 85*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, \
ohy1y2);*)\n\tohy1z2=2. crep tx13 tx5 tx6; (*Rule \
87*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, \
ohy1z2);*)\n\tohz1x2=ohx1z2; (*Rule \
89*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, \
ohz1x2);*)\n\tohz1y2=ohy1z2; (*Rule \
91*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, \
ohz1y2);*)\n\ttx21=2. crep tx13 tx9; (*Rule 93*)\n\tohz1z2=tx18 + tx21; \
(*Rule 94*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, \
ohz1z2);*)\n\tohx2y2=ohx1y1; (*Rule \
96*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, \
ohx2y2);*)\n\tohx2z2=ohx1z1; (*Rule \
98*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, \
ohx2z2);*)\n\tohy2z2=ohy1z1; (*Rule \
100*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, \
ohy2z2);*)\n\t(*} /*if calcOffDiagonalHessian */ *)\n\t(*#endif /* \
EREP_CALC_OFF_DIAGONAL_HESSIAN ]*/*)\n\t(*} /*calcDiagonalHessian \
*/*)\n\t(*#endif /* EREP_CALC_DIAGONAL_HESSIAN ]*/*)\n\t(*} /*calcForce \
*/*)\n\t(*#endif /* EREP_CALC_FORCE ]*/*)\n\t(*SKIP_term:  (void)0;*)\n{Erep, \
fx1, fy1, fz1, fx2, fy2, fz2, dhx1x1, dhy1y1, dhz1z1, dhx2x2, dhy2y2, dhz2z2, \
ohx1y1, ohx1z1, ohx1x2, ohx1y2, ohx1z2, ohy1z1, ohy1x2, ohy1y2, ohy1z2, \
ohz1x2, ohz1y2, ohz1z2, ohx2y2, ohx2z2, ohy2z2, ErepDistance}]\n"}
