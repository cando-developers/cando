{Name -> "Erep", AdditionalCDeclares -> "", 
 Input -> {x1, y1, z1, x2, y2, z2, crep}, 
 Output -> {Erep, fx1, fy1, fz1, fx2, fy2, fz2, dhx1x1, dhy1y1, dhz1z1, 
   dhx2x2, dhy2y2, dhz2z2, ohx1y1, ohx1z1, ohx1x2, ohx1y2, ohx1z2, ohy1z1, 
   ohy1x2, ohy1y2, ohy1z2, ohz1x2, ohz1y2, ohz1z2, ohx2y2, ohx2z2, ohy2z2, 
   ErepDistance}, EnergyFunction -> 
  crep*Log[Sqrt[(x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2]], 
 DerivativeVariables -> {x1, y1, z1, x2, y2, z2}, 
 Rules -> {CCode["EREP_SET_PARAMETER(I1);"], 
   CCode["EREP_SET_PARAMETER(I2);"], CCode["EREP_SET_PARAMETER(crep);"], 
   CCode["EREP_SET_POSITION(x1,I1,0);"], 
   CCode["EREP_SET_POSITION(y1,I1,1);"], 
   CCode["EREP_SET_POSITION(z1,I1,2);"], 
   CCode["EREP_SET_POSITION(x2,I2,0);"], 
   CCode["EREP_SET_POSITION(y2,I2,1);"], 
   CCode["EREP_SET_POSITION(z2,I2,2);"], -x2 -> tx22, -y2 -> tx23, 
   -z2 -> tx24, tx22 + x1 -> tx25, tx23 + y1 -> tx26, tx24 + z1 -> tx27, 
   power2[tx25] -> tx28, power2[tx26] -> tx29, power2[tx27] -> tx30, 
   tx28 + tx29 + tx30 -> tx31, mysqrt[tx31] -> ErepDistance, 
   CCode["MAYBE_BAIL(ErepDistance);"], Log[ErepDistance] -> tx32, 
   crep*tx32 -> Erep, CCode["EREP_ENERGY_ACCUMULATE(Erep);"], 
   CCode["#ifdef EREP_CALC_FORCE //["], CCode["if ( calcForce ) {"], 
   reciprocal[tx31] -> tx33, crep*tx33 -> tzz46, tx25*tzz46 -> gx1, 
   -gx1 -> fx1, CCode["EREP_FORCE_ACCUMULATE(I1, 0, fx1 );"], 
   tx26*tzz46 -> gy1, -gy1 -> fy1, 
   CCode["EREP_FORCE_ACCUMULATE(I1, 1, fy1 );"], tx27*tzz46 -> gz1, 
   -gz1 -> fz1, CCode["EREP_FORCE_ACCUMULATE(I1, 2, fz1 );"], fx1 -> gx2, 
   -gx2 -> fx2, CCode["EREP_FORCE_ACCUMULATE(I2, 0, fx2 );"], fy1 -> gy2, 
   -gy2 -> fy2, CCode["EREP_FORCE_ACCUMULATE(I2, 1, fy2 );"], fz1 -> gz2, 
   -gz2 -> fz2, CCode["EREP_FORCE_ACCUMULATE(I2, 2, fz2 );"], 
   CCode["#ifdef EREP_CALC_DIAGONAL_HESSIAN //["], 
   CCode["if ( calcDiagonalHessian ) {"], power2[tx33] -> tx34, 
   tzz46 -> tx35, crep*tx34 -> tzz43, -2.*tzz43 -> tzz45, tx28*tzz45 -> tx36, 
   tx35 + tx36 -> dhx1x1, 
   CCode["EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);"], 
   tx29*tzz45 -> tx37, tx35 + tx37 -> dhy1y1, 
   CCode["EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);"], 
   tx30*tzz45 -> tx38, tx35 + tx38 -> dhz1z1, 
   CCode["EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);"], 
   dhx1x1 -> dhx2x2, 
   CCode["EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);"], 
   dhy1y1 -> dhy2y2, 
   CCode["EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);"], 
   dhz1z1 -> dhz2z2, 
   CCode["EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);"], 
   CCode["#ifdef EREP_CALC_OFF_DIAGONAL_HESSIAN //["], 
   CCode["if ( calcOffDiagonalHessian ) {"], tx26*tzz45 -> tzz48, 
   tx25*tzz48 -> ohx1y1, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);"], 
   tx25*tx27*tzz45 -> ohx1z1, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);"], 
   2.*tzz43 -> tzz44, tx28*tzz44 -> tx39, -tx35 -> tx40, 
   tx39 + tx40 -> ohx1x2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);"], 
   tx25*tx26*tzz44 -> ohx1y2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);"], 
   tx27*tzz44 -> tzz47, tx25*tzz47 -> ohx1z2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);"], 
   tx27*tzz48 -> ohy1z1, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);"], 
   ohx1y2 -> ohy1x2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);"], 
   tx29*tzz44 -> tx41, tx40 + tx41 -> ohy1y2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);"], 
   tx26*tzz47 -> ohy1z2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);"], 
   ohx1z2 -> ohz1x2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);"], 
   ohy1z2 -> ohz1y2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);"], 
   tx30*tzz44 -> tx42, tx40 + tx42 -> ohz1z2, 
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
   CCode["} /*calcForce */"], CCode["#endif /* EREP_CALC_FORCE ]*/"]}, 
 HessianStructure -> {{8, 14, 15, 16, 17, 18}, {14, 9, 19, 20, 21, 22}, 
   {15, 19, 10, 23, 24, 25}, {16, 20, 23, 11, 26, 27}, 
   {17, 21, 24, 26, 12, 28}, {18, 22, 25, 27, 28, 13}}, CCode -> "NOT USED", 
 MathCode -> "Block[{dhx1x1,dhx2x2,dhy1y1,dhy2y2,dhz1z1,dhz2z2,Erep,ErepDista\
nce,fx1,fx2,fy1,fy2,fz1,fz2,gx1,gx2,gy1,gy2,gz1,gz2,ohx1x2,ohx1y1,ohx1y2,ohx1\
z1,ohx1z2,ohx2y2,ohx2z2,ohy1x2,ohy1y2,ohy1z1,ohy1z2,ohy2z2,ohz1x2,ohz1y2,ohz1\
z2,tx22,tx23,tx24,tx25,tx26,tx27,tx28,tx29,tx30,tx31,tx32,tx33,tx34,tx35,tx36\
,tx37,tx38,tx39,tx40,tx41,tx42,tzz43,tzz44,tzz45,tzz46,tzz47,tzz48,xxxDummy},\
\n\t(*EREP_SET_PARAMETER(I1);*)\n\t(*EREP_SET_PARAMETER(I2);*)\n\t(*EREP_SET_\
PARAMETER(crep);*)\n\t(*EREP_SET_POSITION(x1,I1,0);*)\n\t(*EREP_SET_POSITION(\
y1,I1,1);*)\n\t(*EREP_SET_POSITION(z1,I1,2);*)\n\t(*EREP_SET_POSITION(x2,I2,0\
);*)\n\t(*EREP_SET_POSITION(y2,I2,1);*)\n\t(*EREP_SET_POSITION(z2,I2,2);*)\n\
\ttx22=-x2; (*Rule 10*)\n\ttx23=-y2; (*Rule 11*)\n\ttx24=-z2; (*Rule \
12*)\n\ttx25=tx22 + x1; (*Rule 13*)\n\ttx26=tx23 + y1; (*Rule \
14*)\n\ttx27=tx24 + z1; (*Rule 15*)\n\ttx28= mathPower2[tx25]; (*Rule \
16*)\n\ttx29= mathPower2[tx26]; (*Rule 17*)\n\ttx30= mathPower2[tx27]; (*Rule \
18*)\n\ttx31=tx28 + tx29 + tx30; (*Rule 19*)\n\tErepDistance= mathSqrt[tx31]; \
(*Rule 20*)\n\t(*MAYBE_BAIL(ErepDistance);*)\n\ttx32=Log[ErepDistance]; \
(*Rule 22*)\n\tErep=crep tx32; (*Rule \
23*)\n\t(*EREP_ENERGY_ACCUMULATE(Erep);*)\n\t(*#ifdef EREP_CALC_FORCE \
//[*)\n\t(*if ( calcForce ) {*)\n\ttx33= mathReciprocal[tx31]; (*Rule \
27*)\n\ttzz46=crep tx33; (*Rule 28*)\n\tgx1=tx25 tzz46; (*Rule \
29*)\n\tfx1=-gx1; (*Rule 30*)\n\t(*EREP_FORCE_ACCUMULATE(I1, 0, fx1 \
);*)\n\tgy1=tx26 tzz46; (*Rule 32*)\n\tfy1=-gy1; (*Rule \
33*)\n\t(*EREP_FORCE_ACCUMULATE(I1, 1, fy1 );*)\n\tgz1=tx27 tzz46; (*Rule \
35*)\n\tfz1=-gz1; (*Rule 36*)\n\t(*EREP_FORCE_ACCUMULATE(I1, 2, fz1 \
);*)\n\tgx2=fx1; (*Rule 38*)\n\tfx2=-gx2; (*Rule \
39*)\n\t(*EREP_FORCE_ACCUMULATE(I2, 0, fx2 );*)\n\tgy2=fy1; (*Rule \
41*)\n\tfy2=-gy2; (*Rule 42*)\n\t(*EREP_FORCE_ACCUMULATE(I2, 1, fy2 \
);*)\n\tgz2=fz1; (*Rule 44*)\n\tfz2=-gz2; (*Rule \
45*)\n\t(*EREP_FORCE_ACCUMULATE(I2, 2, fz2 );*)\n\t(*#ifdef \
EREP_CALC_DIAGONAL_HESSIAN //[*)\n\t(*if ( calcDiagonalHessian ) {*)\n\ttx34= \
mathPower2[tx33]; (*Rule 49*)\n\ttx35=tzz46; (*Rule 50*)\n\ttzz43=crep tx34; \
(*Rule 51*)\n\ttzz45=-2. tzz43; (*Rule 52*)\n\ttx36=tx28 tzz45; (*Rule \
53*)\n\tdhx1x1=tx35 + tx36; (*Rule \
54*)\n\t(*EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, \
dhx1x1);*)\n\ttx37=tx29 tzz45; (*Rule 56*)\n\tdhy1y1=tx35 + tx37; (*Rule \
57*)\n\t(*EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, \
dhy1y1);*)\n\ttx38=tx30 tzz45; (*Rule 59*)\n\tdhz1z1=tx35 + tx38; (*Rule \
60*)\n\t(*EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, \
dhz1z1);*)\n\tdhx2x2=dhx1x1; (*Rule \
62*)\n\t(*EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, \
dhx2x2);*)\n\tdhy2y2=dhy1y1; (*Rule \
64*)\n\t(*EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, \
dhy2y2);*)\n\tdhz2z2=dhz1z1; (*Rule \
66*)\n\t(*EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, \
dhz2z2);*)\n\t(*#ifdef EREP_CALC_OFF_DIAGONAL_HESSIAN //[*)\n\t(*if ( \
calcOffDiagonalHessian ) {*)\n\ttzz48=tx26 tzz45; (*Rule 70*)\n\tohx1y1=tx25 \
tzz48; (*Rule 71*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, \
ohx1y1);*)\n\tohx1z1=tx25 tx27 tzz45; (*Rule \
73*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, \
ohx1z1);*)\n\ttzz44=2. tzz43; (*Rule 75*)\n\ttx39=tx28 tzz44; (*Rule \
76*)\n\ttx40=-tx35; (*Rule 77*)\n\tohx1x2=tx39 + tx40; (*Rule \
78*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, \
ohx1x2);*)\n\tohx1y2=tx25 tx26 tzz44; (*Rule \
80*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, \
ohx1y2);*)\n\ttzz47=tx27 tzz44; (*Rule 82*)\n\tohx1z2=tx25 tzz47; (*Rule \
83*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, \
ohx1z2);*)\n\tohy1z1=tx27 tzz48; (*Rule \
85*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, \
ohy1z1);*)\n\tohy1x2=ohx1y2; (*Rule \
87*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, \
ohy1x2);*)\n\ttx41=tx29 tzz44; (*Rule 89*)\n\tohy1y2=tx40 + tx41; (*Rule \
90*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, \
ohy1y2);*)\n\tohy1z2=tx26 tzz47; (*Rule \
92*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, \
ohy1z2);*)\n\tohz1x2=ohx1z2; (*Rule \
94*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, \
ohz1x2);*)\n\tohz1y2=ohy1z2; (*Rule \
96*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, \
ohz1y2);*)\n\ttx42=tx30 tzz44; (*Rule 98*)\n\tohz1z2=tx40 + tx42; (*Rule \
99*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, \
ohz1z2);*)\n\tohx2y2=ohx1y1; (*Rule \
101*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, \
ohx2y2);*)\n\tohx2z2=ohx1z1; (*Rule \
103*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, \
ohx2z2);*)\n\tohy2z2=ohy1z1; (*Rule \
105*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, \
ohy2z2);*)\n\t(*} /*if calcOffDiagonalHessian */ *)\n\t(*#endif /* \
EREP_CALC_OFF_DIAGONAL_HESSIAN ]*/*)\n\t(*} /*calcDiagonalHessian \
*/*)\n\t(*#endif /* EREP_CALC_DIAGONAL_HESSIAN ]*/*)\n\t(*} /*calcForce \
*/*)\n\t(*#endif /* EREP_CALC_FORCE ]*/*)\n{Erep, fx1, fy1, fz1, fx2, fy2, \
fz2, dhx1x1, dhy1y1, dhz1z1, dhx2x2, dhy2y2, dhz2z2, ohx1y1, ohx1z1, ohx1x2, \
ohx1y2, ohx1z2, ohy1z1, ohy1x2, ohy1y2, ohy1z2, ohz1x2, ohz1y2, ohz1z2, \
ohx2y2, ohx2z2, ohy2z2, ErepDistance}]\n"}
