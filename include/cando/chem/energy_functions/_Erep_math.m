{Name -> "Erep", AdditionalCDeclares -> "", 
 Input -> {x1, y1, z1, x2, y2, z2, crep}, 
 Output -> {Erep, fx1, fy1, fz1, fx2, fy2, fz2, dhx1x1, dhy1y1, dhz1z1, 
   dhx2x2, dhy2y2, dhz2z2, ohx1y1, ohx1z1, ohx1x2, ohx1y2, ohx1z2, ohy1z1, 
   ohy1x2, ohy1y2, ohy1z2, ohz1x2, ohz1y2, ohz1z2, ohx2y2, ohx2z2, ohy2z2, 
   ErepDistance}, EnergyFunction -> 
  crep/Sqrt[(x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2], 
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
   CCode["MAYBE_BAIL(ErepDistance);"], reciprocal[ErepDistance] -> tx32, 
   crep*tx32 -> Erep, CCode["EREP_ENERGY_ACCUMULATE(Erep);"], 
   CCode["#ifdef EREP_CALC_FORCE //["], CCode["if ( calcForce ) {"], 
   tx31 -> tx43, tx32 -> tx44, reciprocal[tx43] -> tx45, tx44*tx45 -> tx33, 
   crep*tx33 -> tzz49, -tzz49 -> tzz52, tx25*tzz52 -> gx1, -gx1 -> fx1, 
   CCode["EREP_FORCE_ACCUMULATE(I1, 0, fx1 );"], tx26*tzz52 -> gy1, 
   -gy1 -> fy1, CCode["EREP_FORCE_ACCUMULATE(I1, 1, fy1 );"], 
   tx27*tzz52 -> gz1, -gz1 -> fz1, 
   CCode["EREP_FORCE_ACCUMULATE(I1, 2, fz1 );"], tx25*tzz49 -> gx2, 
   -gx2 -> fx2, CCode["EREP_FORCE_ACCUMULATE(I2, 0, fx2 );"], 
   tx26*tzz49 -> gy2, -gy2 -> fy2, 
   CCode["EREP_FORCE_ACCUMULATE(I2, 1, fy2 );"], tx27*tzz49 -> gz2, 
   -gz2 -> fz2, CCode["EREP_FORCE_ACCUMULATE(I2, 2, fz2 );"], 
   CCode["#ifdef EREP_CALC_DIAGONAL_HESSIAN //["], 
   CCode["if ( calcDiagonalHessian ) {"], power2[tx43] -> tx46, 
   reciprocal[tx46] -> tx47, tx44*tx47 -> tx34, tzz52 -> tx35, 
   crep*tx34 -> tzz48, 3*tzz48 -> tzz50, tx28*tzz50 -> tx36, 
   tx35 + tx36 -> dhx1x1, 
   CCode["EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);"], 
   tx29*tzz50 -> tx37, tx35 + tx37 -> dhy1y1, 
   CCode["EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);"], 
   tx30*tzz50 -> tx38, tx35 + tx38 -> dhz1z1, 
   CCode["EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);"], 
   dhx1x1 -> dhx2x2, 
   CCode["EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);"], 
   dhy1y1 -> dhy2y2, 
   CCode["EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);"], 
   dhz1z1 -> dhz2z2, 
   CCode["EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);"], 
   CCode["#ifdef EREP_CALC_OFF_DIAGONAL_HESSIAN //["], 
   CCode["if ( calcOffDiagonalHessian ) {"], tx25*tx26*tzz50 -> ohx1y1, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);"], 
   tx27*tzz50 -> tzz53, tx25*tzz53 -> ohx1z1, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);"], 
   tzz49 -> tx39, -3*tzz48 -> tzz51, tx28*tzz51 -> tx40, 
   tx39 + tx40 -> ohx1x2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);"], 
   tx26*tzz51 -> tzz54, tx25*tzz54 -> ohx1y2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);"], 
   tx25*tx27*tzz51 -> ohx1z2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);"], 
   tx26*tzz53 -> ohy1z1, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);"], 
   ohx1y2 -> ohy1x2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);"], 
   tx29*tzz51 -> tx41, tx39 + tx41 -> ohy1y2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);"], 
   tx27*tzz54 -> ohy1z2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);"], 
   ohx1z2 -> ohz1x2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);"], 
   ohy1z2 -> ohz1y2, 
   CCode["EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);"], 
   tx30*tzz51 -> tx42, tx39 + tx42 -> ohz1z2, 
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
,tx37,tx38,tx39,tx40,tx41,tx42,tx43,tx44,tx45,tx46,tx47,tzz48,tzz49,tzz50,tzz\
51,tzz52,tzz53,tzz54,xxxDummy},\n\t(*EREP_SET_PARAMETER(I1);*)\n\t(*EREP_SET_\
PARAMETER(I2);*)\n\t(*EREP_SET_PARAMETER(crep);*)\n\t(*EREP_SET_POSITION(x1,I\
1,0);*)\n\t(*EREP_SET_POSITION(y1,I1,1);*)\n\t(*EREP_SET_POSITION(z1,I1,2);*)\
\n\t(*EREP_SET_POSITION(x2,I2,0);*)\n\t(*EREP_SET_POSITION(y2,I2,1);*)\n\t(*E\
REP_SET_POSITION(z2,I2,2);*)\n\ttx22=-x2; (*Rule 10*)\n\ttx23=-y2; (*Rule \
11*)\n\ttx24=-z2; (*Rule 12*)\n\ttx25=tx22 + x1; (*Rule 13*)\n\ttx26=tx23 + \
y1; (*Rule 14*)\n\ttx27=tx24 + z1; (*Rule 15*)\n\ttx28= mathPower2[tx25]; \
(*Rule 16*)\n\ttx29= mathPower2[tx26]; (*Rule 17*)\n\ttx30= mathPower2[tx27]; \
(*Rule 18*)\n\ttx31=tx28 + tx29 + tx30; (*Rule 19*)\n\tErepDistance= \
mathSqrt[tx31]; (*Rule 20*)\n\t(*MAYBE_BAIL(ErepDistance);*)\n\ttx32= \
mathReciprocal[ErepDistance]; (*Rule 22*)\n\tErep=crep tx32; (*Rule \
23*)\n\t(*EREP_ENERGY_ACCUMULATE(Erep);*)\n\t(*#ifdef EREP_CALC_FORCE \
//[*)\n\t(*if ( calcForce ) {*)\n\ttx43=tx31; (*Rule 27*)\n\ttx44=tx32; \
(*Rule 28*)\n\ttx45= mathReciprocal[tx43]; (*Rule 29*)\n\ttx33=tx44 tx45; \
(*Rule 30*)\n\ttzz49=crep tx33; (*Rule 31*)\n\ttzz52=-tzz49; (*Rule \
32*)\n\tgx1=tx25 tzz52; (*Rule 33*)\n\tfx1=-gx1; (*Rule \
34*)\n\t(*EREP_FORCE_ACCUMULATE(I1, 0, fx1 );*)\n\tgy1=tx26 tzz52; (*Rule \
36*)\n\tfy1=-gy1; (*Rule 37*)\n\t(*EREP_FORCE_ACCUMULATE(I1, 1, fy1 \
);*)\n\tgz1=tx27 tzz52; (*Rule 39*)\n\tfz1=-gz1; (*Rule \
40*)\n\t(*EREP_FORCE_ACCUMULATE(I1, 2, fz1 );*)\n\tgx2=tx25 tzz49; (*Rule \
42*)\n\tfx2=-gx2; (*Rule 43*)\n\t(*EREP_FORCE_ACCUMULATE(I2, 0, fx2 \
);*)\n\tgy2=tx26 tzz49; (*Rule 45*)\n\tfy2=-gy2; (*Rule \
46*)\n\t(*EREP_FORCE_ACCUMULATE(I2, 1, fy2 );*)\n\tgz2=tx27 tzz49; (*Rule \
48*)\n\tfz2=-gz2; (*Rule 49*)\n\t(*EREP_FORCE_ACCUMULATE(I2, 2, fz2 \
);*)\n\t(*#ifdef EREP_CALC_DIAGONAL_HESSIAN //[*)\n\t(*if ( \
calcDiagonalHessian ) {*)\n\ttx46= mathPower2[tx43]; (*Rule 53*)\n\ttx47= \
mathReciprocal[tx46]; (*Rule 54*)\n\ttx34=tx44 tx47; (*Rule \
55*)\n\ttx35=tzz52; (*Rule 56*)\n\ttzz48=crep tx34; (*Rule 57*)\n\ttzz50=3 \
tzz48; (*Rule 58*)\n\ttx36=tx28 tzz50; (*Rule 59*)\n\tdhx1x1=tx35 + tx36; \
(*Rule 60*)\n\t(*EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, \
dhx1x1);*)\n\ttx37=tx29 tzz50; (*Rule 62*)\n\tdhy1y1=tx35 + tx37; (*Rule \
63*)\n\t(*EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, \
dhy1y1);*)\n\ttx38=tx30 tzz50; (*Rule 65*)\n\tdhz1z1=tx35 + tx38; (*Rule \
66*)\n\t(*EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, \
dhz1z1);*)\n\tdhx2x2=dhx1x1; (*Rule \
68*)\n\t(*EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, \
dhx2x2);*)\n\tdhy2y2=dhy1y1; (*Rule \
70*)\n\t(*EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, \
dhy2y2);*)\n\tdhz2z2=dhz1z1; (*Rule \
72*)\n\t(*EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, \
dhz2z2);*)\n\t(*#ifdef EREP_CALC_OFF_DIAGONAL_HESSIAN //[*)\n\t(*if ( \
calcOffDiagonalHessian ) {*)\n\tohx1y1=tx25 tx26 tzz50; (*Rule \
76*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, \
ohx1y1);*)\n\ttzz53=tx27 tzz50; (*Rule 78*)\n\tohx1z1=tx25 tzz53; (*Rule \
79*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, \
ohx1z1);*)\n\ttx39=tzz49; (*Rule 81*)\n\ttzz51=-3 tzz48; (*Rule \
82*)\n\ttx40=tx28 tzz51; (*Rule 83*)\n\tohx1x2=tx39 + tx40; (*Rule \
84*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, \
ohx1x2);*)\n\ttzz54=tx26 tzz51; (*Rule 86*)\n\tohx1y2=tx25 tzz54; (*Rule \
87*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, \
ohx1y2);*)\n\tohx1z2=tx25 tx27 tzz51; (*Rule \
89*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, \
ohx1z2);*)\n\tohy1z1=tx26 tzz53; (*Rule \
91*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, \
ohy1z1);*)\n\tohy1x2=ohx1y2; (*Rule \
93*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, \
ohy1x2);*)\n\ttx41=tx29 tzz51; (*Rule 95*)\n\tohy1y2=tx39 + tx41; (*Rule \
96*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, \
ohy1y2);*)\n\tohy1z2=tx27 tzz54; (*Rule \
98*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, \
ohy1z2);*)\n\tohz1x2=ohx1z2; (*Rule \
100*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, \
ohz1x2);*)\n\tohz1y2=ohy1z2; (*Rule \
102*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, \
ohz1y2);*)\n\ttx42=tx30 tzz51; (*Rule 104*)\n\tohz1z2=tx39 + tx42; (*Rule \
105*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, \
ohz1z2);*)\n\tohx2y2=ohx1y1; (*Rule \
107*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, \
ohx2y2);*)\n\tohx2z2=ohx1z1; (*Rule \
109*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, \
ohx2z2);*)\n\tohy2z2=ohy1z1; (*Rule \
111*)\n\t(*EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, \
ohy2z2);*)\n\t(*} /*if calcOffDiagonalHessian */ *)\n\t(*#endif /* \
EREP_CALC_OFF_DIAGONAL_HESSIAN ]*/*)\n\t(*} /*calcDiagonalHessian \
*/*)\n\t(*#endif /* EREP_CALC_DIAGONAL_HESSIAN ]*/*)\n\t(*} /*calcForce \
*/*)\n\t(*#endif /* EREP_CALC_FORCE ]*/*)\n{Erep, fx1, fy1, fz1, fx2, fy2, \
fz2, dhx1x1, dhy1y1, dhz1z1, dhx2x2, dhy2y2, dhz2z2, ohx1y1, ohx1z1, ohx1x2, \
ohx1y2, ohx1z2, ohy1z1, ohy1x2, ohy1y2, ohy1z2, ohz1x2, ohz1y2, ohz1z2, \
ohx2y2, ohx2z2, ohy2z2, ErepDistance}]\n"}
