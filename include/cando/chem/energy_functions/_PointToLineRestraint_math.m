{Name -> "PointToLineRestraint", AdditionalCDeclares -> "", 
 Input -> {x1, y1, z1, xa, ya, za, xb, yb, zb, ka, ra}, 
 Output -> {Energy, fx1, fy1, fz1, dhx1x1, dhy1y1, dhz1z1, ohx1y1, ohx1z1, 
   ohy1z1, PointToLineDeviation}, EnergyFunction -> EnergyFunction, 
 DerivativeVariables -> {x1, y1, z1}, 
 Rules -> {CCode["POINT_TO_LINE_RESTRAINT_SET_PARAMETER(ka);"], 
   CCode["POINT_TO_LINE_RESTRAINT_SET_PARAMETER(ra);"], 
   CCode["POINT_TO_LINE_RESTRAINT_SET_PARAMETER(xa);"], 
   CCode["POINT_TO_LINE_RESTRAINT_SET_PARAMETER(ya);"], 
   CCode["POINT_TO_LINE_RESTRAINT_SET_PARAMETER(za);"], 
   CCode["POINT_TO_LINE_RESTRAINT_SET_PARAMETER(xb);"], 
   CCode["POINT_TO_LINE_RESTRAINT_SET_PARAMETER(yb);"], 
   CCode["POINT_TO_LINE_RESTRAINT_SET_PARAMETER(zb);"], 
   CCode["POINT_TO_LINE_RESTRAINT_SET_PARAMETER(I1);"], 
   CCode["POINT_TO_LINE_RESTRAINT_SET_POSITION(x1,I1,0);"], 
   CCode["POINT_TO_LINE_RESTRAINT_SET_POSITION(y1,I1,1);"], 
   CCode["POINT_TO_LINE_RESTRAINT_SET_POSITION(z1,I1,2);"], -xa -> tx100, 
   tx100*y1 -> tx101, xb*y1 -> tx102, -ya -> tx103, x1*ya -> tx104, 
   tx103*xb -> tx105, -(x1*yb) -> tx106, xa*yb -> tx107, xa*z1 -> tx108, 
   -(xb*z1) -> tx109, tx103*z1 -> tx110, yb*z1 -> tx111, -za -> tx112, 
   tx112*x1 -> tx113, xb*za -> tx114, y1*za -> tx115, tx112*yb -> tx116, 
   x1*zb -> tx117, tx100*zb -> tx118, -(y1*zb) -> tx119, ya*zb -> tx120, 
   tx113 -> tx121, tx105 -> tx122, tx101 -> tx123, tx116 -> tx124, 
   tx110 -> tx125, tx118 -> tx126, tx102 + tx104 + tx106 + tx107 + tx122 + 
     tx123 -> tx127, tx111 + tx115 + tx119 + tx120 + tx124 + tx125 -> tx128, 
   tx108 + tx109 + tx114 + tx117 + tx121 + tx126 -> tx129, 
   tx100 + xb -> tx130, tx103 + yb -> tx131, tx112 + zb -> tx132, 
   power2[tx127] -> tx133, power2[tx128] -> tx134, power2[tx129] -> tx135, 
   power2[tx130] -> tx136, power2[tx131] -> tx137, power2[tx132] -> tx138, 
   tx133 + tx134 + tx135 -> tx139, tx136 + tx137 + tx138 -> tx140, 
   mysqrt[tx139] -> tx141, mysqrt[tx140] -> tx199, 
   reciprocal[tx199] -> tx142, -ra -> tx143, tx141*tx142 -> tx144, 
   tx143 + tx144 -> tx145, power2[tx145] -> tx146, ka*tx146 -> Energy, 
   CCode["POINT_TO_LINE_RESTRAINT_ENERGY_ACCUMULATE(Energy);"], 
   CCode["#ifdef POINT_TO_LINE_RESTRAINT_CALC_FORCE //["], 
   CCode["if ( calcForce ) {"], -yb -> tx147, tx147 + ya -> tx148, 
   2.*tx129 -> tzz216, tx132*tzz216 -> tx149, 2.*tx127 -> tzz218, 
   tx148*tzz218 -> tx150, reciprocal[tx141] -> tx151, tx149 + tx150 -> tx152, 
   tx142*tx145 -> tzz203, ka*tzz203 -> tzz204, tx151*tzz204 -> tzz211, 
   tx152*tzz211 -> gx1, -gx1 -> fx1, 
   CCode["POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 0, fx1 );"], 
   -zb -> tx153, tx153 + za -> tx154, tx130*tzz218 -> tx155, 
   2.*tx128 -> tzz217, tx154*tzz217 -> tx156, tx155 + tx156 -> tx157, 
   tx157*tzz211 -> gy1, -gy1 -> fy1, 
   CCode["POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 1, fy1 );"], 
   -xb -> tx158, tx158 + xa -> tx159, tx131*tzz217 -> tx160, 
   tx159*tzz216 -> tx161, tx160 + tx161 -> tx162, tx162*tzz211 -> gz1, 
   -gz1 -> fz1, 
   CCode["POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 2, fz1 );"], 
   CCode["#ifdef POINT_TO_LINE_RESTRAINT_CALC_DIAGONAL_HESSIAN //["], 
   CCode["if ( calcDiagonalHessian ) {"], power2[tx148] -> tx163, 
   2.*tx138 -> tx164, 2.*tx163 -> tx165, tx139 -> tx200, tx151 -> tx201, 
   reciprocal[tx200] -> tx202, tx201*tx202 -> tx166, tx202 -> tx167, 
   reciprocal[tx140] -> tx168, power2[tx152] -> tx169, 
   tx164 + tx165 -> tx170, tx166*tzz204 -> tzz207, -0.5*tzz207 -> tzz210, 
   tx169*tzz210 -> tx171, tx167*tx168 -> tzz206, ka*tzz206 -> tzz208, 
   0.5*tzz208 -> tzz209, tx169*tzz209 -> tx172, tx201*tzz204 -> tzz205, 
   tx170*tzz205 -> tx173, tx171 + tx172 + tx173 -> dhx1x1, 
   CCode["POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, \
dhx1x1);"], power2[tx154] -> tx174, 2.*tx136 -> tx175, 2.*tx174 -> tx176, 
   power2[tx157] -> tx177, tx175 + tx176 -> tx178, tx177*tzz210 -> tx179, 
   tx177*tzz209 -> tx180, tx178*tzz205 -> tx181, 
   tx179 + tx180 + tx181 -> dhy1y1, CCode["POINT_TO_LINE_RESTRAINT_DIAGONAL_H\
ESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);"], power2[tx159] -> tx182, 
   2.*tx137 -> tx183, 2.*tx182 -> tx184, power2[tx162] -> tx185, 
   tx183 + tx184 -> tx186, tx185*tzz210 -> tx187, tx185*tzz209 -> tx188, 
   tx186*tzz205 -> tx189, tx187 + tx188 + tx189 -> dhz1z1, 
   CCode["POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, \
dhz1z1);"], 
   CCode["#ifdef POINT_TO_LINE_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN //["], 
   CCode["if ( calcOffDiagonalHessian ) {"], 2.*tzz205 -> tzz212, 
   tx130*tx148*tzz212 -> tx190, tx152*tx157 -> tzz215, 
   tzz210*tzz215 -> tx191, tzz209*tzz215 -> tx192, 
   tx190 + tx191 + tx192 -> ohx1y1, CCode["POINT_TO_LINE_RESTRAINT_OFF_DIAGON\
AL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);"], tx132*tx159*tzz212 -> tx193, 
   tx162*tzz210 -> tzz213, tx152*tzz213 -> tx194, tx162*tzz209 -> tzz214, 
   tx152*tzz214 -> tx195, tx193 + tx194 + tx195 -> ohx1z1, 
   CCode["POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, \
2, ohx1z1);"], tx131*tx154*tzz212 -> tx196, tx157*tzz213 -> tx197, 
   tx157*tzz214 -> tx198, tx196 + tx197 + tx198 -> ohy1z1, 
   CCode["POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, \
2, ohy1z1);"], CCode["} /*if calcOffDiagonalHessian */ "], 
   CCode["#endif /* POINT_TO_LINE_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN ]*/"], 
   CCode["} /*calcDiagonalHessian */"], 
   CCode["#endif /* POINT_TO_LINE_RESTRAINT_CALC_DIAGONAL_HESSIAN ]*/"], 
   CCode["} /*calcForce */"], 
   CCode["#endif /* POINT_TO_LINE_RESTRAINT_CALC_FORCE ]*/"]}, 
 HessianStructure -> {{5, 8, 9, 0, 0, 0}, {8, 6, 10, 0, 0, 0}, 
   {9, 10, 7, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, 
   {0, 0, 0, 0, 0, 0}}, CCode -> "NOT USED", MathCode -> "Block[{dhx1x1,dhy1y\
1,dhz1z1,Energy,fx1,fy1,fz1,gx1,gy1,gz1,ohx1y1,ohx1z1,ohy1z1,PointToLineDevia\
tion,tx100,tx101,tx102,tx103,tx104,tx105,tx106,tx107,tx108,tx109,tx110,tx111,\
tx112,tx113,tx114,tx115,tx116,tx117,tx118,tx119,tx120,tx121,tx122,tx123,tx124\
,tx125,tx126,tx127,tx128,tx129,tx130,tx131,tx132,tx133,tx134,tx135,tx136,tx13\
7,tx138,tx139,tx140,tx141,tx142,tx143,tx144,tx145,tx146,tx147,tx148,tx149,tx1\
50,tx151,tx152,tx153,tx154,tx155,tx156,tx157,tx158,tx159,tx160,tx161,tx162,tx\
163,tx164,tx165,tx166,tx167,tx168,tx169,tx170,tx171,tx172,tx173,tx174,tx175,t\
x176,tx177,tx178,tx179,tx180,tx181,tx182,tx183,tx184,tx185,tx186,tx187,tx188,\
tx189,tx190,tx191,tx192,tx193,tx194,tx195,tx196,tx197,tx198,tx199,tx200,tx201\
,tx202,tzz203,tzz204,tzz205,tzz206,tzz207,tzz208,tzz209,tzz210,tzz211,tzz212,\
tzz213,tzz214,tzz215,tzz216,tzz217,tzz218,xxxDummy},\n\t(*POINT_TO_LINE_RESTR\
AINT_SET_PARAMETER(ka);*)\n\t(*POINT_TO_LINE_RESTRAINT_SET_PARAMETER(ra);*)\n\
\t(*POINT_TO_LINE_RESTRAINT_SET_PARAMETER(xa);*)\n\t(*POINT_TO_LINE_RESTRAINT\
_SET_PARAMETER(ya);*)\n\t(*POINT_TO_LINE_RESTRAINT_SET_PARAMETER(za);*)\n\t(*\
POINT_TO_LINE_RESTRAINT_SET_PARAMETER(xb);*)\n\t(*POINT_TO_LINE_RESTRAINT_SET\
_PARAMETER(yb);*)\n\t(*POINT_TO_LINE_RESTRAINT_SET_PARAMETER(zb);*)\n\t(*POIN\
T_TO_LINE_RESTRAINT_SET_PARAMETER(I1);*)\n\t(*POINT_TO_LINE_RESTRAINT_SET_POS\
ITION(x1,I1,0);*)\n\t(*POINT_TO_LINE_RESTRAINT_SET_POSITION(y1,I1,1);*)\n\t(*\
POINT_TO_LINE_RESTRAINT_SET_POSITION(z1,I1,2);*)\n\ttx100=-xa; (*Rule \
13*)\n\ttx101=tx100 y1; (*Rule 14*)\n\ttx102=xb y1; (*Rule 15*)\n\ttx103=-ya; \
(*Rule 16*)\n\ttx104=x1 ya; (*Rule 17*)\n\ttx105=tx103 xb; (*Rule \
18*)\n\ttx106=-(x1 yb); (*Rule 19*)\n\ttx107=xa yb; (*Rule 20*)\n\ttx108=xa \
z1; (*Rule 21*)\n\ttx109=-(xb z1); (*Rule 22*)\n\ttx110=tx103 z1; (*Rule \
23*)\n\ttx111=yb z1; (*Rule 24*)\n\ttx112=-za; (*Rule 25*)\n\ttx113=tx112 x1; \
(*Rule 26*)\n\ttx114=xb za; (*Rule 27*)\n\ttx115=y1 za; (*Rule \
28*)\n\ttx116=tx112 yb; (*Rule 29*)\n\ttx117=x1 zb; (*Rule \
30*)\n\ttx118=tx100 zb; (*Rule 31*)\n\ttx119=-(y1 zb); (*Rule \
32*)\n\ttx120=ya zb; (*Rule 33*)\n\ttx121=tx113; (*Rule 34*)\n\ttx122=tx105; \
(*Rule 35*)\n\ttx123=tx101; (*Rule 36*)\n\ttx124=tx116; (*Rule \
37*)\n\ttx125=tx110; (*Rule 38*)\n\ttx126=tx118; (*Rule 39*)\n\ttx127=tx102 + \
tx104 + tx106 + tx107 + tx122 + tx123; (*Rule 40*)\n\ttx128=tx111 + tx115 + \
tx119 + tx120 + tx124 + tx125; (*Rule 41*)\n\ttx129=tx108 + tx109 + tx114 + \
tx117 + tx121 + tx126; (*Rule 42*)\n\ttx130=tx100 + xb; (*Rule \
43*)\n\ttx131=tx103 + yb; (*Rule 44*)\n\ttx132=tx112 + zb; (*Rule \
45*)\n\ttx133= mathPower2[tx127]; (*Rule 46*)\n\ttx134= mathPower2[tx128]; \
(*Rule 47*)\n\ttx135= mathPower2[tx129]; (*Rule 48*)\n\ttx136= \
mathPower2[tx130]; (*Rule 49*)\n\ttx137= mathPower2[tx131]; (*Rule \
50*)\n\ttx138= mathPower2[tx132]; (*Rule 51*)\n\ttx139=tx133 + tx134 + tx135; \
(*Rule 52*)\n\ttx140=tx136 + tx137 + tx138; (*Rule 53*)\n\ttx141= \
mathSqrt[tx139]; (*Rule 54*)\n\ttx199= mathSqrt[tx140]; (*Rule 55*)\n\ttx142= \
mathReciprocal[tx199]; (*Rule 56*)\n\ttx143=-ra; (*Rule 57*)\n\ttx144=tx141 \
tx142; (*Rule 58*)\n\ttx145=tx143 + tx144; (*Rule 59*)\n\ttx146= \
mathPower2[tx145]; (*Rule 60*)\n\tEnergy=ka tx146; (*Rule \
61*)\n\t(*POINT_TO_LINE_RESTRAINT_ENERGY_ACCUMULATE(Energy);*)\n\t(*#ifdef \
POINT_TO_LINE_RESTRAINT_CALC_FORCE //[*)\n\t(*if ( calcForce ) \
{*)\n\ttx147=-yb; (*Rule 65*)\n\ttx148=tx147 + ya; (*Rule 66*)\n\ttzz216=2. \
tx129; (*Rule 67*)\n\ttx149=tx132 tzz216; (*Rule 68*)\n\ttzz218=2. tx127; \
(*Rule 69*)\n\ttx150=tx148 tzz218; (*Rule 70*)\n\ttx151= \
mathReciprocal[tx141]; (*Rule 71*)\n\ttx152=tx149 + tx150; (*Rule \
72*)\n\ttzz203=tx142 tx145; (*Rule 73*)\n\ttzz204=ka tzz203; (*Rule \
74*)\n\ttzz211=tx151 tzz204; (*Rule 75*)\n\tgx1=tx152 tzz211; (*Rule \
76*)\n\tfx1=-gx1; (*Rule \
77*)\n\t(*POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 0, fx1 \
);*)\n\ttx153=-zb; (*Rule 79*)\n\ttx154=tx153 + za; (*Rule \
80*)\n\ttx155=tx130 tzz218; (*Rule 81*)\n\ttzz217=2. tx128; (*Rule \
82*)\n\ttx156=tx154 tzz217; (*Rule 83*)\n\ttx157=tx155 + tx156; (*Rule \
84*)\n\tgy1=tx157 tzz211; (*Rule 85*)\n\tfy1=-gy1; (*Rule \
86*)\n\t(*POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 1, fy1 \
);*)\n\ttx158=-xb; (*Rule 88*)\n\ttx159=tx158 + xa; (*Rule \
89*)\n\ttx160=tx131 tzz217; (*Rule 90*)\n\ttx161=tx159 tzz216; (*Rule \
91*)\n\ttx162=tx160 + tx161; (*Rule 92*)\n\tgz1=tx162 tzz211; (*Rule \
93*)\n\tfz1=-gz1; (*Rule \
94*)\n\t(*POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 2, fz1 \
);*)\n\t(*#ifdef POINT_TO_LINE_RESTRAINT_CALC_DIAGONAL_HESSIAN //[*)\n\t(*if \
( calcDiagonalHessian ) {*)\n\ttx163= mathPower2[tx148]; (*Rule \
98*)\n\ttx164=2. tx138; (*Rule 99*)\n\ttx165=2. tx163; (*Rule \
100*)\n\ttx200=tx139; (*Rule 101*)\n\ttx201=tx151; (*Rule 102*)\n\ttx202= \
mathReciprocal[tx200]; (*Rule 103*)\n\ttx166=tx201 tx202; (*Rule \
104*)\n\ttx167=tx202; (*Rule 105*)\n\ttx168= mathReciprocal[tx140]; (*Rule \
106*)\n\ttx169= mathPower2[tx152]; (*Rule 107*)\n\ttx170=tx164 + tx165; \
(*Rule 108*)\n\ttzz207=tx166 tzz204; (*Rule 109*)\n\ttzz210=-0.5 tzz207; \
(*Rule 110*)\n\ttx171=tx169 tzz210; (*Rule 111*)\n\ttzz206=tx167 tx168; \
(*Rule 112*)\n\ttzz208=ka tzz206; (*Rule 113*)\n\ttzz209=0.5 tzz208; (*Rule \
114*)\n\ttx172=tx169 tzz209; (*Rule 115*)\n\ttzz205=tx201 tzz204; (*Rule \
116*)\n\ttx173=tx170 tzz205; (*Rule 117*)\n\tdhx1x1=tx171 + tx172 + tx173; \
(*Rule 118*)\n\t(*POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, \
I1, 0, dhx1x1);*)\n\ttx174= mathPower2[tx154]; (*Rule 120*)\n\ttx175=2. \
tx136; (*Rule 121*)\n\ttx176=2. tx174; (*Rule 122*)\n\ttx177= \
mathPower2[tx157]; (*Rule 123*)\n\ttx178=tx175 + tx176; (*Rule \
124*)\n\ttx179=tx177 tzz210; (*Rule 125*)\n\ttx180=tx177 tzz209; (*Rule \
126*)\n\ttx181=tx178 tzz205; (*Rule 127*)\n\tdhy1y1=tx179 + tx180 + tx181; \
(*Rule 128*)\n\t(*POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, \
I1, 1, dhy1y1);*)\n\ttx182= mathPower2[tx159]; (*Rule 130*)\n\ttx183=2. \
tx137; (*Rule 131*)\n\ttx184=2. tx182; (*Rule 132*)\n\ttx185= \
mathPower2[tx162]; (*Rule 133*)\n\ttx186=tx183 + tx184; (*Rule \
134*)\n\ttx187=tx185 tzz210; (*Rule 135*)\n\ttx188=tx185 tzz209; (*Rule \
136*)\n\ttx189=tx186 tzz205; (*Rule 137*)\n\tdhz1z1=tx187 + tx188 + tx189; \
(*Rule 138*)\n\t(*POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, \
I1, 2, dhz1z1);*)\n\t(*#ifdef \
POINT_TO_LINE_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN //[*)\n\t(*if ( \
calcOffDiagonalHessian ) {*)\n\ttzz212=2. tzz205; (*Rule 142*)\n\ttx190=tx130 \
tx148 tzz212; (*Rule 143*)\n\ttzz215=tx152 tx157; (*Rule \
144*)\n\ttx191=tzz210 tzz215; (*Rule 145*)\n\ttx192=tzz209 tzz215; (*Rule \
146*)\n\tohx1y1=tx190 + tx191 + tx192; (*Rule \
147*)\n\t(*POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, \
1, ohx1y1);*)\n\ttx193=tx132 tx159 tzz212; (*Rule 149*)\n\ttzz213=tx162 \
tzz210; (*Rule 150*)\n\ttx194=tx152 tzz213; (*Rule 151*)\n\ttzz214=tx162 \
tzz209; (*Rule 152*)\n\ttx195=tx152 tzz214; (*Rule 153*)\n\tohx1z1=tx193 + \
tx194 + tx195; (*Rule \
154*)\n\t(*POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, \
2, ohx1z1);*)\n\ttx196=tx131 tx154 tzz212; (*Rule 156*)\n\ttx197=tx157 \
tzz213; (*Rule 157*)\n\ttx198=tx157 tzz214; (*Rule 158*)\n\tohy1z1=tx196 + \
tx197 + tx198; (*Rule \
159*)\n\t(*POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, \
2, ohy1z1);*)\n\t(*} /*if calcOffDiagonalHessian */ *)\n\t(*#endif /* \
POINT_TO_LINE_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN ]*/*)\n\t(*} \
/*calcDiagonalHessian */*)\n\t(*#endif /* \
POINT_TO_LINE_RESTRAINT_CALC_DIAGONAL_HESSIAN ]*/*)\n\t(*} /*calcForce \
*/*)\n\t(*#endif /* POINT_TO_LINE_RESTRAINT_CALC_FORCE ]*/*)\n{Energy, fx1, \
fy1, fz1, dhx1x1, dhy1y1, dhz1z1, ohx1y1, ohx1z1, ohy1z1, \
PointToLineDeviation}]\n"}
