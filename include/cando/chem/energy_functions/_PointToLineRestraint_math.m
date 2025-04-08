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
   CCode["POINT_TO_LINE_RESTRAINT_APPLY_ATOM_MASK(I1);"], 
   CCode["POINT_TO_LINE_RESTRAINT_SET_POSITION(x1,I1,0);"], 
   CCode["POINT_TO_LINE_RESTRAINT_SET_POSITION(y1,I1,1);"], 
   CCode["POINT_TO_LINE_RESTRAINT_SET_POSITION(z1,I1,2);"], -xa -> tx1, 
   -(xa*y1) -> tx2, xb*y1 -> tx3, -ya -> tx4, x1*ya -> tx5, -(xb*ya) -> tx6, 
   -(x1*yb) -> tx7, xa*yb -> tx8, xa*z1 -> tx9, -z1 -> tzz132, 
   tzz132*xb -> tx10, tzz132*ya -> tx11, yb*z1 -> tx12, -za -> tzz131, 
   tzz131 -> tx13, tx13*x1 -> tx14, xb*za -> tx15, y1*za -> tx16, 
   tx13*yb -> tx17, x1*zb -> tx18, -zb -> tzz130, tzz130*xa -> tx19, 
   tzz130*y1 -> tx20, ya*zb -> tx21, tx14 -> tx22, tx4*xb -> tx23, 
   tx1*y1 -> tx24, tx17 -> tx25, tx4*z1 -> tx26, tx1*zb -> tx27, 
   tx21 + tx25 -> tzz137, tx20 + tzz137 -> tzz138, tx16 + tzz138 -> tzz140, 
   tx12 + tzz140 -> tzz142, tx26 + tzz142 -> tx28, tx7 + tx8 -> tzz133, 
   tx5 + tzz133 -> tzz134, tx3 + tzz134 -> tzz135, 
   tx23 + tx24 + tzz135 -> tx29, tx22 + tx9 -> tzz136, 
   tx18 + tzz136 -> tzz139, tx15 + tzz139 -> tzz141, tx10 + tzz141 -> tzz143, 
   tx27 + tzz143 -> tx30, tx1 + xb -> tx31, tx4 + yb -> tx32, 
   tx13 + zb -> tx33, power2[tx28] -> tx34, power2[tx29] -> tx35, 
   power2[tx30] -> tx36, power2[tx31] -> tx37, power2[tx32] -> tx38, 
   power2[tx33] -> tx39, tx34 + tx35 + tx36 -> tx40, 
   tx37 + tx38 + tx39 -> tx41, mysqrt[tx40] -> tx42, mysqrt[tx41] -> tx110, 
   reciprocal[tx110] -> tx43, -ra -> tx44, tx42*tx43 -> tx45, 
   tx44 + tx45 -> tx46, power2[tx46] -> tx47, ka*tx47 -> Energy, 
   CCode["POINT_TO_LINE_RESTRAINT_ENERGY_ACCUMULATE(Energy);"], 
   CCode["#ifdef POINT_TO_LINE_RESTRAINT_CALC_FORCE //["], 
   CCode["if ( calcForce ) {"], tx11 + tzz142 -> tx48, 
   tx2 + tx6 + tzz135 -> tx49, tx19 + tzz143 -> tx50, power2[tx48] -> tx51, 
   power2[tx49] -> tx52, power2[tx50] -> tx53, tx51 + tx52 + tx53 -> tx54, 
   -yb -> tx55, mysqrt[tx54] -> tx56, tx55 + ya -> tx57, 2.*tx50 -> tzz127, 
   tx33*tzz127 -> tx58, tx43*tx56 -> tx59, 2.*tx49 -> tzz128, 
   tx57*tzz128 -> tx60, reciprocal[tx56] -> tx61, tx44 + tx59 -> tx62, 
   tx58 + tx60 -> tx63, tx43*tx62 -> tzz114, ka*tzz114 -> tzz115, 
   tx61*tzz115 -> tzz122, tx63*tzz122 -> gx1, -gx1 -> fx1, 
   CCode["POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 0, fx1 );"], 
   tzz130 -> tx64, tx64 + za -> tx65, tx31*tzz128 -> tx66, 2.*tx48 -> tzz129, 
   tx65*tzz129 -> tx67, tx66 + tx67 -> tx68, tx68*tzz122 -> gy1, -gy1 -> fy1, 
   CCode["POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 1, fy1 );"], 
   -xb -> tx69, tx69 + xa -> tx70, tx32*tzz129 -> tx71, tx70*tzz127 -> tx72, 
   tx71 + tx72 -> tx73, tx73*tzz122 -> gz1, -gz1 -> fz1, 
   CCode["POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 2, fz1 );"], 
   CCode["#ifdef POINT_TO_LINE_RESTRAINT_CALC_DIAGONAL_HESSIAN //["], 
   CCode["if ( calcDiagonalHessian ) {"], power2[tx57] -> tx74, 
   2.*tx39 -> tx75, 2.*tx74 -> tx76, reciprocal[tx41] -> tx77, tx54 -> tx111, 
   reciprocal[tx111] -> tx112, tx61 -> tx113, tx112*tx113 -> tx78, 
   tx112 -> tx79, power2[tx63] -> tx80, tx75 + tx76 -> tx81, 
   tx78*tzz115 -> tzz116, -0.5*tzz116 -> tzz121, tx80*tzz121 -> tx82, 
   tx77*tx79 -> tzz117, ka*tzz117 -> tzz119, 0.5*tzz119 -> tzz120, 
   tx80*tzz120 -> tx83, tx113*tzz115 -> tzz118, tx81*tzz118 -> tx84, 
   tx82 + tx83 + tx84 -> dhx1x1, CCode["POINT_TO_LINE_RESTRAINT_DIAGONAL_HESS\
IAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);"], power2[tx65] -> tx85, 
   2.*tx37 -> tx86, 2.*tx85 -> tx87, power2[tx68] -> tx88, 
   tx86 + tx87 -> tx89, tx88*tzz121 -> tx90, tx88*tzz120 -> tx91, 
   tx89*tzz118 -> tx92, tx90 + tx91 + tx92 -> dhy1y1, 
   CCode["POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, \
dhy1y1);"], power2[tx70] -> tx93, 2.*tx38 -> tx94, 2.*tx93 -> tx95, 
   power2[tx73] -> tx96, tx94 + tx95 -> tx97, tx96*tzz121 -> tx98, 
   tx96*tzz120 -> tx99, tx97*tzz118 -> tx100, tx100 + tx98 + tx99 -> dhz1z1, 
   CCode["POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, \
dhz1z1);"], 
   CCode["#ifdef POINT_TO_LINE_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN //["], 
   CCode["if ( calcOffDiagonalHessian ) {"], 2.*tzz118 -> tzz123, 
   tx31*tx57*tzz123 -> tx101, tx63*tx68 -> tzz126, tzz121*tzz126 -> tx102, 
   tzz120*tzz126 -> tx103, tx101 + tx102 + tx103 -> ohx1y1, 
   CCode["POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, \
1, ohx1y1);"], tx33*tx70*tzz123 -> tx104, tx73*tzz121 -> tzz124, 
   tx63*tzz124 -> tx105, tx73*tzz120 -> tzz125, tx63*tzz125 -> tx106, 
   tx104 + tx105 + tx106 -> ohx1z1, CCode["POINT_TO_LINE_RESTRAINT_OFF_DIAGON\
AL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);"], tx32*tx65*tzz123 -> tx107, 
   tx68*tzz124 -> tx108, tx68*tzz125 -> tx109, tx107 + tx108 + tx109 -> 
    ohy1z1, CCode["POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1\
, 1, I1, 2, ohy1z1);"], CCode["} /*if calcOffDiagonalHessian */ "], 
   CCode["#endif /* POINT_TO_LINE_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN ]*/"], 
   CCode["} /*calcDiagonalHessian */"], 
   CCode["#endif /* POINT_TO_LINE_RESTRAINT_CALC_DIAGONAL_HESSIAN ]*/"], 
   CCode["} /*calcForce */"], 
   CCode["#endif /* POINT_TO_LINE_RESTRAINT_CALC_FORCE ]*/"], 
   CCode["POINT_TO_LINE_RESTRAINT_DEBUG_INTERACTIONS(I1);"], 
   CCode["SKIP_term: (void)0;"]}, HessianStructure -> 
  {{5, 8, 9, 0, 0, 0}, {8, 6, 10, 0, 0, 0}, {9, 10, 7, 0, 0, 0}, 
   {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}}, 
 CCode -> "NOT USED", MathCode -> "Block[{dhx1x1,dhy1y1,dhz1z1,Energy,fx1,fy1\
,fz1,gx1,gy1,gz1,ohx1y1,ohx1z1,ohy1z1,PointToLineDeviation,tx1,tx10,tx100,tx1\
01,tx102,tx103,tx104,tx105,tx106,tx107,tx108,tx109,tx11,tx110,tx111,tx112,tx1\
13,tx12,tx13,tx14,tx15,tx16,tx17,tx18,tx19,tx2,tx20,tx21,tx22,tx23,tx24,tx25,\
tx26,tx27,tx28,tx29,tx3,tx30,tx31,tx32,tx33,tx34,tx35,tx36,tx37,tx38,tx39,tx4\
,tx40,tx41,tx42,tx43,tx44,tx45,tx46,tx47,tx48,tx49,tx5,tx50,tx51,tx52,tx53,tx\
54,tx55,tx56,tx57,tx58,tx59,tx6,tx60,tx61,tx62,tx63,tx64,tx65,tx66,tx67,tx68,\
tx69,tx7,tx70,tx71,tx72,tx73,tx74,tx75,tx76,tx77,tx78,tx79,tx8,tx80,tx81,tx82\
,tx83,tx84,tx85,tx86,tx87,tx88,tx89,tx9,tx90,tx91,tx92,tx93,tx94,tx95,tx96,tx\
97,tx98,tx99,tzz114,tzz115,tzz116,tzz117,tzz118,tzz119,tzz120,tzz121,tzz122,t\
zz123,tzz124,tzz125,tzz126,tzz127,tzz128,tzz129,tzz130,tzz131,tzz132,tzz133,t\
zz134,tzz135,tzz136,tzz137,tzz138,tzz139,tzz140,tzz141,tzz142,tzz143,xxxDummy\
},\n\t(*POINT_TO_LINE_RESTRAINT_SET_PARAMETER(ka);*)\n\t(*POINT_TO_LINE_RESTR\
AINT_SET_PARAMETER(ra);*)\n\t(*POINT_TO_LINE_RESTRAINT_SET_PARAMETER(xa);*)\n\
\t(*POINT_TO_LINE_RESTRAINT_SET_PARAMETER(ya);*)\n\t(*POINT_TO_LINE_RESTRAINT\
_SET_PARAMETER(za);*)\n\t(*POINT_TO_LINE_RESTRAINT_SET_PARAMETER(xb);*)\n\t(*\
POINT_TO_LINE_RESTRAINT_SET_PARAMETER(yb);*)\n\t(*POINT_TO_LINE_RESTRAINT_SET\
_PARAMETER(zb);*)\n\t(*POINT_TO_LINE_RESTRAINT_SET_PARAMETER(I1);*)\n\t(*POIN\
T_TO_LINE_RESTRAINT_APPLY_ATOM_MASK(I1);*)\n\t(*POINT_TO_LINE_RESTRAINT_SET_P\
OSITION(x1,I1,0);*)\n\t(*POINT_TO_LINE_RESTRAINT_SET_POSITION(y1,I1,1);*)\n\t\
(*POINT_TO_LINE_RESTRAINT_SET_POSITION(z1,I1,2);*)\n\ttx1=-xa; (*Rule \
14*)\n\ttx2=-(xa y1); (*Rule 15*)\n\ttx3=xb y1; (*Rule 16*)\n\ttx4=-ya; \
(*Rule 17*)\n\ttx5=x1 ya; (*Rule 18*)\n\ttx6=-(xb ya); (*Rule \
19*)\n\ttx7=-(x1 yb); (*Rule 20*)\n\ttx8=xa yb; (*Rule 21*)\n\ttx9=xa z1; \
(*Rule 22*)\n\ttzz132=-z1; (*Rule 23*)\n\ttx10=tzz132 xb; (*Rule \
24*)\n\ttx11=tzz132 ya; (*Rule 25*)\n\ttx12=yb z1; (*Rule 26*)\n\ttzz131=-za; \
(*Rule 27*)\n\ttx13=tzz131; (*Rule 28*)\n\ttx14=tx13 x1; (*Rule \
29*)\n\ttx15=xb za; (*Rule 30*)\n\ttx16=y1 za; (*Rule 31*)\n\ttx17=tx13 yb; \
(*Rule 32*)\n\ttx18=x1 zb; (*Rule 33*)\n\ttzz130=-zb; (*Rule \
34*)\n\ttx19=tzz130 xa; (*Rule 35*)\n\ttx20=tzz130 y1; (*Rule 36*)\n\ttx21=ya \
zb; (*Rule 37*)\n\ttx22=tx14; (*Rule 38*)\n\ttx23=tx4 xb; (*Rule \
39*)\n\ttx24=tx1 y1; (*Rule 40*)\n\ttx25=tx17; (*Rule 41*)\n\ttx26=tx4 z1; \
(*Rule 42*)\n\ttx27=tx1 zb; (*Rule 43*)\n\ttzz137=tx21 + tx25; (*Rule \
44*)\n\ttzz138=tx20 + tzz137; (*Rule 45*)\n\ttzz140=tx16 + tzz138; (*Rule \
46*)\n\ttzz142=tx12 + tzz140; (*Rule 47*)\n\ttx28=tx26 + tzz142; (*Rule \
48*)\n\ttzz133=tx7 + tx8; (*Rule 49*)\n\ttzz134=tx5 + tzz133; (*Rule \
50*)\n\ttzz135=tx3 + tzz134; (*Rule 51*)\n\ttx29=tx23 + tx24 + tzz135; (*Rule \
52*)\n\ttzz136=tx22 + tx9; (*Rule 53*)\n\ttzz139=tx18 + tzz136; (*Rule \
54*)\n\ttzz141=tx15 + tzz139; (*Rule 55*)\n\ttzz143=tx10 + tzz141; (*Rule \
56*)\n\ttx30=tx27 + tzz143; (*Rule 57*)\n\ttx31=tx1 + xb; (*Rule \
58*)\n\ttx32=tx4 + yb; (*Rule 59*)\n\ttx33=tx13 + zb; (*Rule 60*)\n\ttx34= \
mathPower2[tx28]; (*Rule 61*)\n\ttx35= mathPower2[tx29]; (*Rule 62*)\n\ttx36= \
mathPower2[tx30]; (*Rule 63*)\n\ttx37= mathPower2[tx31]; (*Rule 64*)\n\ttx38= \
mathPower2[tx32]; (*Rule 65*)\n\ttx39= mathPower2[tx33]; (*Rule \
66*)\n\ttx40=tx34 + tx35 + tx36; (*Rule 67*)\n\ttx41=tx37 + tx38 + tx39; \
(*Rule 68*)\n\ttx42= mathSqrt[tx40]; (*Rule 69*)\n\ttx110= mathSqrt[tx41]; \
(*Rule 70*)\n\ttx43= mathReciprocal[tx110]; (*Rule 71*)\n\ttx44=-ra; (*Rule \
72*)\n\ttx45=tx42 tx43; (*Rule 73*)\n\ttx46=tx44 + tx45; (*Rule 74*)\n\ttx47= \
mathPower2[tx46]; (*Rule 75*)\n\tEnergy=ka tx47; (*Rule \
76*)\n\t(*POINT_TO_LINE_RESTRAINT_ENERGY_ACCUMULATE(Energy);*)\n\t(*#ifdef \
POINT_TO_LINE_RESTRAINT_CALC_FORCE //[*)\n\t(*if ( calcForce ) \
{*)\n\ttx48=tx11 + tzz142; (*Rule 80*)\n\ttx49=tx2 + tx6 + tzz135; (*Rule \
81*)\n\ttx50=tx19 + tzz143; (*Rule 82*)\n\ttx51= mathPower2[tx48]; (*Rule \
83*)\n\ttx52= mathPower2[tx49]; (*Rule 84*)\n\ttx53= mathPower2[tx50]; (*Rule \
85*)\n\ttx54=tx51 + tx52 + tx53; (*Rule 86*)\n\ttx55=-yb; (*Rule \
87*)\n\ttx56= mathSqrt[tx54]; (*Rule 88*)\n\ttx57=tx55 + ya; (*Rule \
89*)\n\ttzz127=2. tx50; (*Rule 90*)\n\ttx58=tx33 tzz127; (*Rule \
91*)\n\ttx59=tx43 tx56; (*Rule 92*)\n\ttzz128=2. tx49; (*Rule \
93*)\n\ttx60=tx57 tzz128; (*Rule 94*)\n\ttx61= mathReciprocal[tx56]; (*Rule \
95*)\n\ttx62=tx44 + tx59; (*Rule 96*)\n\ttx63=tx58 + tx60; (*Rule \
97*)\n\ttzz114=tx43 tx62; (*Rule 98*)\n\ttzz115=ka tzz114; (*Rule \
99*)\n\ttzz122=tx61 tzz115; (*Rule 100*)\n\tgx1=tx63 tzz122; (*Rule \
101*)\n\tfx1=-gx1; (*Rule \
102*)\n\t(*POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 0, fx1 \
);*)\n\ttx64=tzz130; (*Rule 104*)\n\ttx65=tx64 + za; (*Rule \
105*)\n\ttx66=tx31 tzz128; (*Rule 106*)\n\ttzz129=2. tx48; (*Rule \
107*)\n\ttx67=tx65 tzz129; (*Rule 108*)\n\ttx68=tx66 + tx67; (*Rule \
109*)\n\tgy1=tx68 tzz122; (*Rule 110*)\n\tfy1=-gy1; (*Rule \
111*)\n\t(*POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 1, fy1 \
);*)\n\ttx69=-xb; (*Rule 113*)\n\ttx70=tx69 + xa; (*Rule 114*)\n\ttx71=tx32 \
tzz129; (*Rule 115*)\n\ttx72=tx70 tzz127; (*Rule 116*)\n\ttx73=tx71 + tx72; \
(*Rule 117*)\n\tgz1=tx73 tzz122; (*Rule 118*)\n\tfz1=-gz1; (*Rule \
119*)\n\t(*POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 2, fz1 \
);*)\n\t(*#ifdef POINT_TO_LINE_RESTRAINT_CALC_DIAGONAL_HESSIAN //[*)\n\t(*if \
( calcDiagonalHessian ) {*)\n\ttx74= mathPower2[tx57]; (*Rule \
123*)\n\ttx75=2. tx39; (*Rule 124*)\n\ttx76=2. tx74; (*Rule 125*)\n\ttx77= \
mathReciprocal[tx41]; (*Rule 126*)\n\ttx111=tx54; (*Rule 127*)\n\ttx112= \
mathReciprocal[tx111]; (*Rule 128*)\n\ttx113=tx61; (*Rule 129*)\n\ttx78=tx112 \
tx113; (*Rule 130*)\n\ttx79=tx112; (*Rule 131*)\n\ttx80= mathPower2[tx63]; \
(*Rule 132*)\n\ttx81=tx75 + tx76; (*Rule 133*)\n\ttzz116=tx78 tzz115; (*Rule \
134*)\n\ttzz121=-0.5 tzz116; (*Rule 135*)\n\ttx82=tx80 tzz121; (*Rule \
136*)\n\ttzz117=tx77 tx79; (*Rule 137*)\n\ttzz119=ka tzz117; (*Rule \
138*)\n\ttzz120=0.5 tzz119; (*Rule 139*)\n\ttx83=tx80 tzz120; (*Rule \
140*)\n\ttzz118=tx113 tzz115; (*Rule 141*)\n\ttx84=tx81 tzz118; (*Rule \
142*)\n\tdhx1x1=tx82 + tx83 + tx84; (*Rule \
143*)\n\t(*POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, \
dhx1x1);*)\n\ttx85= mathPower2[tx65]; (*Rule 145*)\n\ttx86=2. tx37; (*Rule \
146*)\n\ttx87=2. tx85; (*Rule 147*)\n\ttx88= mathPower2[tx68]; (*Rule \
148*)\n\ttx89=tx86 + tx87; (*Rule 149*)\n\ttx90=tx88 tzz121; (*Rule \
150*)\n\ttx91=tx88 tzz120; (*Rule 151*)\n\ttx92=tx89 tzz118; (*Rule \
152*)\n\tdhy1y1=tx90 + tx91 + tx92; (*Rule \
153*)\n\t(*POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, \
dhy1y1);*)\n\ttx93= mathPower2[tx70]; (*Rule 155*)\n\ttx94=2. tx38; (*Rule \
156*)\n\ttx95=2. tx93; (*Rule 157*)\n\ttx96= mathPower2[tx73]; (*Rule \
158*)\n\ttx97=tx94 + tx95; (*Rule 159*)\n\ttx98=tx96 tzz121; (*Rule \
160*)\n\ttx99=tx96 tzz120; (*Rule 161*)\n\ttx100=tx97 tzz118; (*Rule \
162*)\n\tdhz1z1=tx100 + tx98 + tx99; (*Rule \
163*)\n\t(*POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, \
dhz1z1);*)\n\t(*#ifdef POINT_TO_LINE_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN \
//[*)\n\t(*if ( calcOffDiagonalHessian ) {*)\n\ttzz123=2. tzz118; (*Rule \
167*)\n\ttx101=tx31 tx57 tzz123; (*Rule 168*)\n\ttzz126=tx63 tx68; (*Rule \
169*)\n\ttx102=tzz121 tzz126; (*Rule 170*)\n\ttx103=tzz120 tzz126; (*Rule \
171*)\n\tohx1y1=tx101 + tx102 + tx103; (*Rule \
172*)\n\t(*POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, \
1, ohx1y1);*)\n\ttx104=tx33 tx70 tzz123; (*Rule 174*)\n\ttzz124=tx73 tzz121; \
(*Rule 175*)\n\ttx105=tx63 tzz124; (*Rule 176*)\n\ttzz125=tx73 tzz120; (*Rule \
177*)\n\ttx106=tx63 tzz125; (*Rule 178*)\n\tohx1z1=tx104 + tx105 + tx106; \
(*Rule 179*)\n\t(*POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, \
0, I1, 2, ohx1z1);*)\n\ttx107=tx32 tx65 tzz123; (*Rule 181*)\n\ttx108=tx68 \
tzz124; (*Rule 182*)\n\ttx109=tx68 tzz125; (*Rule 183*)\n\tohy1z1=tx107 + \
tx108 + tx109; (*Rule \
184*)\n\t(*POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, \
2, ohy1z1);*)\n\t(*} /*if calcOffDiagonalHessian */ *)\n\t(*#endif /* \
POINT_TO_LINE_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN ]*/*)\n\t(*} \
/*calcDiagonalHessian */*)\n\t(*#endif /* \
POINT_TO_LINE_RESTRAINT_CALC_DIAGONAL_HESSIAN ]*/*)\n\t(*} /*calcForce \
*/*)\n\t(*#endif /* POINT_TO_LINE_RESTRAINT_CALC_FORCE \
]*/*)\n\t(*POINT_TO_LINE_RESTRAINT_DEBUG_INTERACTIONS(I1);*)\n\t(*SKIP_term: \
(void)0;*)\n{Energy, fx1, fy1, fz1, dhx1x1, dhy1y1, dhz1z1, ohx1y1, ohx1z1, \
ohy1z1, PointToLineDeviation}]\n"}
