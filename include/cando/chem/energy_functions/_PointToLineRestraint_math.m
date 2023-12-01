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
   tx1*y1 -> tx2, xb*y1 -> tx3, -ya -> tx4, x1*ya -> tx5, tx4*xb -> tx6, 
   -(x1*yb) -> tx7, xa*yb -> tx8, xa*z1 -> tx9, -(xb*z1) -> tx10, 
   tx4*z1 -> tx11, yb*z1 -> tx12, -za -> tx13, tx13*x1 -> tx14, 
   xb*za -> tx15, y1*za -> tx16, tx13*yb -> tx17, x1*zb -> tx18, 
   tx1*zb -> tx19, -(y1*zb) -> tx20, ya*zb -> tx21, tx14 -> tx22, 
   tx6 -> tx23, tx2 -> tx24, tx17 -> tx25, tx11 -> tx26, tx19 -> tx27, 
   tx12 + tx16 + tx20 + tx21 + tx25 + tx26 -> tx28, 
   tx23 + tx24 + tx3 + tx5 + tx7 + tx8 -> tx29, 
   tx10 + tx15 + tx18 + tx22 + tx27 + tx9 -> tx30, tx1 + xb -> tx31, 
   tx4 + yb -> tx32, tx13 + zb -> tx33, power2[tx28] -> tx34, 
   power2[tx29] -> tx35, power2[tx30] -> tx36, power2[tx31] -> tx37, 
   power2[tx32] -> tx38, power2[tx33] -> tx39, tx34 + tx35 + tx36 -> tx40, 
   tx37 + tx38 + tx39 -> tx41, mysqrt[tx40] -> tx42, mysqrt[tx41] -> tx100, 
   reciprocal[tx100] -> tx43, -ra -> tx44, tx42*tx43 -> tx45, 
   tx44 + tx45 -> tx46, power2[tx46] -> tx47, ka*tx47 -> Energy, 
   CCode["POINT_TO_LINE_RESTRAINT_ENERGY_ACCUMULATE(Energy);"], 
   CCode["#ifdef POINT_TO_LINE_RESTRAINT_CALC_FORCE //["], 
   CCode["if ( calcForce ) {"], -yb -> tx48, tx48 + ya -> tx49, 
   2.*tx30*tx33 -> tx50, 2.*tx29*tx49 -> tx51, reciprocal[tx42] -> tx52, 
   tx50 + tx51 -> tx53, ka*tx43*tx46*tx52*tx53 -> gx1, -gx1 -> fx1, 
   CCode["POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 0, fx1 );"], 
   -zb -> tx54, tx54 + za -> tx55, 2.*tx29*tx31 -> tx56, 
   2.*tx28*tx55 -> tx57, tx56 + tx57 -> tx58, ka*tx43*tx46*tx52*tx58 -> gy1, 
   -gy1 -> fy1, 
   CCode["POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 1, fy1 );"], 
   -xb -> tx59, tx59 + xa -> tx60, 2.*tx28*tx32 -> tx61, 
   2.*tx30*tx60 -> tx62, tx61 + tx62 -> tx63, ka*tx43*tx46*tx52*tx63 -> gz1, 
   -gz1 -> fz1, 
   CCode["POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 2, fz1 );"], 
   CCode["#ifdef POINT_TO_LINE_RESTRAINT_CALC_DIAGONAL_HESSIAN //["], 
   CCode["if ( calcDiagonalHessian ) {"], power2[tx49] -> tx64, 
   2.*tx39 -> tx65, 2.*tx64 -> tx66, tx40 -> tx101, 
   reciprocal[tx101] -> tx102, tx52 -> tx103, tx102*tx103 -> tx67, 
   tx102 -> tx68, reciprocal[tx41] -> tx69, power2[tx53] -> tx70, 
   tx65 + tx66 -> tx71, -0.5*ka*tx43*tx46*tx67*tx70 -> tx72, 
   0.5*ka*tx68*tx69*tx70 -> tx73, ka*tx103*tx43*tx46*tx71 -> tx74, 
   tx72 + tx73 + tx74 -> dhx1x1, CCode["POINT_TO_LINE_RESTRAINT_DIAGONAL_HESS\
IAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);"], power2[tx55] -> tx75, 
   2.*tx37 -> tx76, 2.*tx75 -> tx77, power2[tx58] -> tx78, 
   tx76 + tx77 -> tx79, -0.5*ka*tx43*tx46*tx67*tx78 -> tx80, 
   0.5*ka*tx68*tx69*tx78 -> tx81, ka*tx103*tx43*tx46*tx79 -> tx82, 
   tx80 + tx81 + tx82 -> dhy1y1, CCode["POINT_TO_LINE_RESTRAINT_DIAGONAL_HESS\
IAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);"], power2[tx60] -> tx83, 
   2.*tx38 -> tx84, 2.*tx83 -> tx85, power2[tx63] -> tx86, 
   tx84 + tx85 -> tx87, -0.5*ka*tx43*tx46*tx67*tx86 -> tx88, 
   0.5*ka*tx68*tx69*tx86 -> tx89, ka*tx103*tx43*tx46*tx87 -> tx90, 
   tx88 + tx89 + tx90 -> dhz1z1, CCode["POINT_TO_LINE_RESTRAINT_DIAGONAL_HESS\
IAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);"], 
   CCode["#ifdef POINT_TO_LINE_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN //["], 
   CCode["if ( calcOffDiagonalHessian ) {"], 
   2.*ka*tx103*tx31*tx43*tx46*tx49 -> tx91, 
   -0.5*ka*tx43*tx46*tx53*tx58*tx67 -> tx92, 0.5*ka*tx53*tx58*tx68*tx69 -> 
    tx93, tx91 + tx92 + tx93 -> ohx1y1, CCode["POINT_TO_LINE_RESTRAINT_OFF_DI\
AGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);"], 
   2.*ka*tx103*tx33*tx43*tx46*tx60 -> tx94, 
   -0.5*ka*tx43*tx46*tx53*tx63*tx67 -> tx95, 0.5*ka*tx53*tx63*tx68*tx69 -> 
    tx96, tx94 + tx95 + tx96 -> ohx1z1, CCode["POINT_TO_LINE_RESTRAINT_OFF_DI\
AGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);"], 
   2.*ka*tx103*tx32*tx43*tx46*tx55 -> tx97, 
   -0.5*ka*tx43*tx46*tx58*tx63*tx67 -> tx98, 0.5*ka*tx58*tx63*tx68*tx69 -> 
    tx99, tx97 + tx98 + tx99 -> ohy1z1, CCode["POINT_TO_LINE_RESTRAINT_OFF_DI\
AGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);"], 
   CCode["} /*if calcOffDiagonalHessian */ "], 
   CCode["#endif /* POINT_TO_LINE_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN ]*/"], 
   CCode["} /*calcDiagonalHessian */"], 
   CCode["#endif /* POINT_TO_LINE_RESTRAINT_CALC_DIAGONAL_HESSIAN ]*/"], 
   CCode["} /*calcForce */"], 
   CCode["#endif /* POINT_TO_LINE_RESTRAINT_CALC_FORCE ]*/"], 
   CCode["SKIP_term:  (void)0;"]}, HessianStructure -> 
  {{5, 8, 9, 0, 0, 0}, {8, 6, 10, 0, 0, 0}, {9, 10, 7, 0, 0, 0}, 
   {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}}, 
 CCode -> "NOT USED", MathCode -> "Block[{dhx1x1,dhy1y1,dhz1z1,Energy,fx1,fy1\
,fz1,gx1,gy1,gz1,ohx1y1,ohx1z1,ohy1z1,PointToLineDeviation,tx1,tx10,tx100,tx1\
01,tx102,tx103,tx11,tx12,tx13,tx14,tx15,tx16,tx17,tx18,tx19,tx2,tx20,tx21,tx2\
2,tx23,tx24,tx25,tx26,tx27,tx28,tx29,tx3,tx30,tx31,tx32,tx33,tx34,tx35,tx36,t\
x37,tx38,tx39,tx4,tx40,tx41,tx42,tx43,tx44,tx45,tx46,tx47,tx48,tx49,tx5,tx50,\
tx51,tx52,tx53,tx54,tx55,tx56,tx57,tx58,tx59,tx6,tx60,tx61,tx62,tx63,tx64,tx6\
5,tx66,tx67,tx68,tx69,tx7,tx70,tx71,tx72,tx73,tx74,tx75,tx76,tx77,tx78,tx79,t\
x8,tx80,tx81,tx82,tx83,tx84,tx85,tx86,tx87,tx88,tx89,tx9,tx90,tx91,tx92,tx93,\
tx94,tx95,tx96,tx97,tx98,tx99,xxxDummy},\n\t(*POINT_TO_LINE_RESTRAINT_SET_PAR\
AMETER(ka);*)\n\t(*POINT_TO_LINE_RESTRAINT_SET_PARAMETER(ra);*)\n\t(*POINT_TO\
_LINE_RESTRAINT_SET_PARAMETER(xa);*)\n\t(*POINT_TO_LINE_RESTRAINT_SET_PARAMET\
ER(ya);*)\n\t(*POINT_TO_LINE_RESTRAINT_SET_PARAMETER(za);*)\n\t(*POINT_TO_LIN\
E_RESTRAINT_SET_PARAMETER(xb);*)\n\t(*POINT_TO_LINE_RESTRAINT_SET_PARAMETER(y\
b);*)\n\t(*POINT_TO_LINE_RESTRAINT_SET_PARAMETER(zb);*)\n\t(*POINT_TO_LINE_RE\
STRAINT_SET_PARAMETER(I1);*)\n\t(*POINT_TO_LINE_RESTRAINT_APPLY_ATOM_MASK(I1)\
;*)\n\t(*POINT_TO_LINE_RESTRAINT_SET_POSITION(x1,I1,0);*)\n\t(*POINT_TO_LINE_\
RESTRAINT_SET_POSITION(y1,I1,1);*)\n\t(*POINT_TO_LINE_RESTRAINT_SET_POSITION(\
z1,I1,2);*)\n\ttx1=-xa; (*Rule 14*)\n\ttx2=tx1 y1; (*Rule 15*)\n\ttx3=xb y1; \
(*Rule 16*)\n\ttx4=-ya; (*Rule 17*)\n\ttx5=x1 ya; (*Rule 18*)\n\ttx6=tx4 xb; \
(*Rule 19*)\n\ttx7=-(x1 yb); (*Rule 20*)\n\ttx8=xa yb; (*Rule 21*)\n\ttx9=xa \
z1; (*Rule 22*)\n\ttx10=-(xb z1); (*Rule 23*)\n\ttx11=tx4 z1; (*Rule \
24*)\n\ttx12=yb z1; (*Rule 25*)\n\ttx13=-za; (*Rule 26*)\n\ttx14=tx13 x1; \
(*Rule 27*)\n\ttx15=xb za; (*Rule 28*)\n\ttx16=y1 za; (*Rule \
29*)\n\ttx17=tx13 yb; (*Rule 30*)\n\ttx18=x1 zb; (*Rule 31*)\n\ttx19=tx1 zb; \
(*Rule 32*)\n\ttx20=-(y1 zb); (*Rule 33*)\n\ttx21=ya zb; (*Rule \
34*)\n\ttx22=tx14; (*Rule 35*)\n\ttx23=tx6; (*Rule 36*)\n\ttx24=tx2; (*Rule \
37*)\n\ttx25=tx17; (*Rule 38*)\n\ttx26=tx11; (*Rule 39*)\n\ttx27=tx19; (*Rule \
40*)\n\ttx28=tx12 + tx16 + tx20 + tx21 + tx25 + tx26; (*Rule \
41*)\n\ttx29=tx23 + tx24 + tx3 + tx5 + tx7 + tx8; (*Rule 42*)\n\ttx30=tx10 + \
tx15 + tx18 + tx22 + tx27 + tx9; (*Rule 43*)\n\ttx31=tx1 + xb; (*Rule \
44*)\n\ttx32=tx4 + yb; (*Rule 45*)\n\ttx33=tx13 + zb; (*Rule 46*)\n\ttx34= \
mathPower2[tx28]; (*Rule 47*)\n\ttx35= mathPower2[tx29]; (*Rule 48*)\n\ttx36= \
mathPower2[tx30]; (*Rule 49*)\n\ttx37= mathPower2[tx31]; (*Rule 50*)\n\ttx38= \
mathPower2[tx32]; (*Rule 51*)\n\ttx39= mathPower2[tx33]; (*Rule \
52*)\n\ttx40=tx34 + tx35 + tx36; (*Rule 53*)\n\ttx41=tx37 + tx38 + tx39; \
(*Rule 54*)\n\ttx42= mathSqrt[tx40]; (*Rule 55*)\n\ttx100= mathSqrt[tx41]; \
(*Rule 56*)\n\ttx43= mathReciprocal[tx100]; (*Rule 57*)\n\ttx44=-ra; (*Rule \
58*)\n\ttx45=tx42 tx43; (*Rule 59*)\n\ttx46=tx44 + tx45; (*Rule 60*)\n\ttx47= \
mathPower2[tx46]; (*Rule 61*)\n\tEnergy=ka tx47; (*Rule \
62*)\n\t(*POINT_TO_LINE_RESTRAINT_ENERGY_ACCUMULATE(Energy);*)\n\t(*#ifdef \
POINT_TO_LINE_RESTRAINT_CALC_FORCE //[*)\n\t(*if ( calcForce ) \
{*)\n\ttx48=-yb; (*Rule 66*)\n\ttx49=tx48 + ya; (*Rule 67*)\n\ttx50=2. tx30 \
tx33; (*Rule 68*)\n\ttx51=2. tx29 tx49; (*Rule 69*)\n\ttx52= \
mathReciprocal[tx42]; (*Rule 70*)\n\ttx53=tx50 + tx51; (*Rule 71*)\n\tgx1=ka \
tx43 tx46 tx52 tx53; (*Rule 72*)\n\tfx1=-gx1; (*Rule \
73*)\n\t(*POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 0, fx1 \
);*)\n\ttx54=-zb; (*Rule 75*)\n\ttx55=tx54 + za; (*Rule 76*)\n\ttx56=2. tx29 \
tx31; (*Rule 77*)\n\ttx57=2. tx28 tx55; (*Rule 78*)\n\ttx58=tx56 + tx57; \
(*Rule 79*)\n\tgy1=ka tx43 tx46 tx52 tx58; (*Rule 80*)\n\tfy1=-gy1; (*Rule \
81*)\n\t(*POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 1, fy1 \
);*)\n\ttx59=-xb; (*Rule 83*)\n\ttx60=tx59 + xa; (*Rule 84*)\n\ttx61=2. tx28 \
tx32; (*Rule 85*)\n\ttx62=2. tx30 tx60; (*Rule 86*)\n\ttx63=tx61 + tx62; \
(*Rule 87*)\n\tgz1=ka tx43 tx46 tx52 tx63; (*Rule 88*)\n\tfz1=-gz1; (*Rule \
89*)\n\t(*POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 2, fz1 \
);*)\n\t(*#ifdef POINT_TO_LINE_RESTRAINT_CALC_DIAGONAL_HESSIAN //[*)\n\t(*if \
( calcDiagonalHessian ) {*)\n\ttx64= mathPower2[tx49]; (*Rule 93*)\n\ttx65=2. \
tx39; (*Rule 94*)\n\ttx66=2. tx64; (*Rule 95*)\n\ttx101=tx40; (*Rule \
96*)\n\ttx102= mathReciprocal[tx101]; (*Rule 97*)\n\ttx103=tx52; (*Rule \
98*)\n\ttx67=tx102 tx103; (*Rule 99*)\n\ttx68=tx102; (*Rule 100*)\n\ttx69= \
mathReciprocal[tx41]; (*Rule 101*)\n\ttx70= mathPower2[tx53]; (*Rule \
102*)\n\ttx71=tx65 + tx66; (*Rule 103*)\n\ttx72=-0.5 ka tx43 tx46 tx67 tx70; \
(*Rule 104*)\n\ttx73=0.5 ka tx68 tx69 tx70; (*Rule 105*)\n\ttx74=ka tx103 \
tx43 tx46 tx71; (*Rule 106*)\n\tdhx1x1=tx72 + tx73 + tx74; (*Rule \
107*)\n\t(*POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, \
dhx1x1);*)\n\ttx75= mathPower2[tx55]; (*Rule 109*)\n\ttx76=2. tx37; (*Rule \
110*)\n\ttx77=2. tx75; (*Rule 111*)\n\ttx78= mathPower2[tx58]; (*Rule \
112*)\n\ttx79=tx76 + tx77; (*Rule 113*)\n\ttx80=-0.5 ka tx43 tx46 tx67 tx78; \
(*Rule 114*)\n\ttx81=0.5 ka tx68 tx69 tx78; (*Rule 115*)\n\ttx82=ka tx103 \
tx43 tx46 tx79; (*Rule 116*)\n\tdhy1y1=tx80 + tx81 + tx82; (*Rule \
117*)\n\t(*POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, \
dhy1y1);*)\n\ttx83= mathPower2[tx60]; (*Rule 119*)\n\ttx84=2. tx38; (*Rule \
120*)\n\ttx85=2. tx83; (*Rule 121*)\n\ttx86= mathPower2[tx63]; (*Rule \
122*)\n\ttx87=tx84 + tx85; (*Rule 123*)\n\ttx88=-0.5 ka tx43 tx46 tx67 tx86; \
(*Rule 124*)\n\ttx89=0.5 ka tx68 tx69 tx86; (*Rule 125*)\n\ttx90=ka tx103 \
tx43 tx46 tx87; (*Rule 126*)\n\tdhz1z1=tx88 + tx89 + tx90; (*Rule \
127*)\n\t(*POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, \
dhz1z1);*)\n\t(*#ifdef POINT_TO_LINE_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN \
//[*)\n\t(*if ( calcOffDiagonalHessian ) {*)\n\ttx91=2. ka tx103 tx31 tx43 \
tx46 tx49; (*Rule 131*)\n\ttx92=-0.5 ka tx43 tx46 tx53 tx58 tx67; (*Rule \
132*)\n\ttx93=0.5 ka tx53 tx58 tx68 tx69; (*Rule 133*)\n\tohx1y1=tx91 + tx92 \
+ tx93; (*Rule \
134*)\n\t(*POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, \
1, ohx1y1);*)\n\ttx94=2. ka tx103 tx33 tx43 tx46 tx60; (*Rule \
136*)\n\ttx95=-0.5 ka tx43 tx46 tx53 tx63 tx67; (*Rule 137*)\n\ttx96=0.5 ka \
tx53 tx63 tx68 tx69; (*Rule 138*)\n\tohx1z1=tx94 + tx95 + tx96; (*Rule \
139*)\n\t(*POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, \
2, ohx1z1);*)\n\ttx97=2. ka tx103 tx32 tx43 tx46 tx55; (*Rule \
141*)\n\ttx98=-0.5 ka tx43 tx46 tx58 tx63 tx67; (*Rule 142*)\n\ttx99=0.5 ka \
tx58 tx63 tx68 tx69; (*Rule 143*)\n\tohy1z1=tx97 + tx98 + tx99; (*Rule \
144*)\n\t(*POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, \
2, ohy1z1);*)\n\t(*} /*if calcOffDiagonalHessian */ *)\n\t(*#endif /* \
POINT_TO_LINE_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN ]*/*)\n\t(*} \
/*calcDiagonalHessian */*)\n\t(*#endif /* \
POINT_TO_LINE_RESTRAINT_CALC_DIAGONAL_HESSIAN ]*/*)\n\t(*} /*calcForce \
*/*)\n\t(*#endif /* POINT_TO_LINE_RESTRAINT_CALC_FORCE ]*/*)\n\t(*SKIP_term:  \
(void)0;*)\n{Energy, fx1, fy1, fz1, dhx1x1, dhy1y1, dhz1z1, ohx1y1, ohx1z1, \
ohy1z1, PointToLineDeviation}]\n"}
