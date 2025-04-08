// C-code
	POINT_TO_LINE_RESTRAINT_SET_PARAMETER(ka);
	POINT_TO_LINE_RESTRAINT_SET_PARAMETER(ra);
	POINT_TO_LINE_RESTRAINT_SET_PARAMETER(xa);
	POINT_TO_LINE_RESTRAINT_SET_PARAMETER(ya);
	POINT_TO_LINE_RESTRAINT_SET_PARAMETER(za);
	POINT_TO_LINE_RESTRAINT_SET_PARAMETER(xb);
	POINT_TO_LINE_RESTRAINT_SET_PARAMETER(yb);
	POINT_TO_LINE_RESTRAINT_SET_PARAMETER(zb);
	POINT_TO_LINE_RESTRAINT_SET_PARAMETER(I1);
	POINT_TO_LINE_RESTRAINT_APPLY_ATOM_MASK(I1);
	POINT_TO_LINE_RESTRAINT_SET_POSITION(x1,I1,0);
	POINT_TO_LINE_RESTRAINT_SET_POSITION(y1,I1,1);
	POINT_TO_LINE_RESTRAINT_SET_POSITION(z1,I1,2);
	 tx1 = -xa; 		/* rule 14 */
	 tx2 = -(xa*y1); 		/* rule 15 */
	 tx3 = xb*y1; 		/* rule 16 */
	 tx4 = -ya; 		/* rule 17 */
	 tx5 = x1*ya; 		/* rule 18 */
	 tx6 = -(xb*ya); 		/* rule 19 */
	 tx7 = -(x1*yb); 		/* rule 20 */
	 tx8 = xa*yb; 		/* rule 21 */
	 tx9 = xa*z1; 		/* rule 22 */
	 tzz132 = -z1; 		/* rule 23 */
	 tx10 = tzz132*xb; 		/* rule 24 */
	 tx11 = tzz132*ya; 		/* rule 25 */
	 tx12 = yb*z1; 		/* rule 26 */
	 tzz131 = -za; 		/* rule 27 */
	 tx13 = tzz131; 		/* rule 28 */
	 tx14 = tx13*x1; 		/* rule 29 */
	 tx15 = xb*za; 		/* rule 30 */
	 tx16 = y1*za; 		/* rule 31 */
	 tx17 = tx13*yb; 		/* rule 32 */
	 tx18 = x1*zb; 		/* rule 33 */
	 tzz130 = -zb; 		/* rule 34 */
	 tx19 = tzz130*xa; 		/* rule 35 */
	 tx20 = tzz130*y1; 		/* rule 36 */
	 tx21 = ya*zb; 		/* rule 37 */
	 tx22 = tx14; 		/* rule 38 */
	 tx23 = tx4*xb; 		/* rule 39 */
	 tx24 = tx1*y1; 		/* rule 40 */
	 tx25 = tx17; 		/* rule 41 */
	 tx26 = tx4*z1; 		/* rule 42 */
	 tx27 = tx1*zb; 		/* rule 43 */
	 tzz137 = tx21 + tx25; 		/* rule 44 */
	 tzz138 = tx20 + tzz137; 		/* rule 45 */
	 tzz140 = tx16 + tzz138; 		/* rule 46 */
	 tzz142 = tx12 + tzz140; 		/* rule 47 */
	 tx28 = tx26 + tzz142; 		/* rule 48 */
	 tzz133 = tx7 + tx8; 		/* rule 49 */
	 tzz134 = tx5 + tzz133; 		/* rule 50 */
	 tzz135 = tx3 + tzz134; 		/* rule 51 */
	 tx29 = tx23 + tx24 + tzz135; 		/* rule 52 */
	 tzz136 = tx22 + tx9; 		/* rule 53 */
	 tzz139 = tx18 + tzz136; 		/* rule 54 */
	 tzz141 = tx15 + tzz139; 		/* rule 55 */
	 tzz143 = tx10 + tzz141; 		/* rule 56 */
	 tx30 = tx27 + tzz143; 		/* rule 57 */
	 tx31 = tx1 + xb; 		/* rule 58 */
	 tx32 = tx4 + yb; 		/* rule 59 */
	 tx33 = tx13 + zb; 		/* rule 60 */
	 tx34 = power2(tx28); 		/* rule 61 */
	 tx35 = power2(tx29); 		/* rule 62 */
	 tx36 = power2(tx30); 		/* rule 63 */
	 tx37 = power2(tx31); 		/* rule 64 */
	 tx38 = power2(tx32); 		/* rule 65 */
	 tx39 = power2(tx33); 		/* rule 66 */
	 tx40 = tx34 + tx35 + tx36; 		/* rule 67 */
	 tx41 = tx37 + tx38 + tx39; 		/* rule 68 */
	 tx42 = mysqrt(tx40); 		/* rule 69 */
	 tx110 = mysqrt(tx41); 		/* rule 70 */
	 tx43 = reciprocal(tx110); 		/* rule 71 */
	 tx44 = -ra; 		/* rule 72 */
	 tx45 = tx42*tx43; 		/* rule 73 */
	 tx46 = tx44 + tx45; 		/* rule 74 */
	 tx47 = power2(tx46); 		/* rule 75 */
	 Energy = ka*tx47; 		/* rule 76 */
	POINT_TO_LINE_RESTRAINT_ENERGY_ACCUMULATE(Energy);
	#ifdef POINT_TO_LINE_RESTRAINT_CALC_FORCE //[
	if ( calcForce ) {
	 tx48 = tx11 + tzz142; 		/* rule 80 */
	 tx49 = tx2 + tx6 + tzz135; 		/* rule 81 */
	 tx50 = tx19 + tzz143; 		/* rule 82 */
	 tx51 = power2(tx48); 		/* rule 83 */
	 tx52 = power2(tx49); 		/* rule 84 */
	 tx53 = power2(tx50); 		/* rule 85 */
	 tx54 = tx51 + tx52 + tx53; 		/* rule 86 */
	 tx55 = -yb; 		/* rule 87 */
	 tx56 = mysqrt(tx54); 		/* rule 88 */
	 tx57 = tx55 + ya; 		/* rule 89 */
	 tzz127 = 2.*tx50; 		/* rule 90 */
	 tx58 = tx33*tzz127; 		/* rule 91 */
	 tx59 = tx43*tx56; 		/* rule 92 */
	 tzz128 = 2.*tx49; 		/* rule 93 */
	 tx60 = tx57*tzz128; 		/* rule 94 */
	 tx61 = reciprocal(tx56); 		/* rule 95 */
	 tx62 = tx44 + tx59; 		/* rule 96 */
	 tx63 = tx58 + tx60; 		/* rule 97 */
	 tzz114 = tx43*tx62; 		/* rule 98 */
	 tzz115 = ka*tzz114; 		/* rule 99 */
	 tzz122 = tx61*tzz115; 		/* rule 100 */
	 gx1 = tx63*tzz122; 		/* rule 101 */
	 fx1 = -gx1; 		/* rule 102 */
	POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx64 = tzz130; 		/* rule 104 */
	 tx65 = tx64 + za; 		/* rule 105 */
	 tx66 = tx31*tzz128; 		/* rule 106 */
	 tzz129 = 2.*tx48; 		/* rule 107 */
	 tx67 = tx65*tzz129; 		/* rule 108 */
	 tx68 = tx66 + tx67; 		/* rule 109 */
	 gy1 = tx68*tzz122; 		/* rule 110 */
	 fy1 = -gy1; 		/* rule 111 */
	POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx69 = -xb; 		/* rule 113 */
	 tx70 = tx69 + xa; 		/* rule 114 */
	 tx71 = tx32*tzz129; 		/* rule 115 */
	 tx72 = tx70*tzz127; 		/* rule 116 */
	 tx73 = tx71 + tx72; 		/* rule 117 */
	 gz1 = tx73*tzz122; 		/* rule 118 */
	 fz1 = -gz1; 		/* rule 119 */
	POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 2, fz1 );
	#ifdef POINT_TO_LINE_RESTRAINT_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx74 = power2(tx57); 		/* rule 123 */
	 tx75 = 2.*tx39; 		/* rule 124 */
	 tx76 = 2.*tx74; 		/* rule 125 */
	 tx77 = reciprocal(tx41); 		/* rule 126 */
	 tx111 = tx54; 		/* rule 127 */
	 tx112 = reciprocal(tx111); 		/* rule 128 */
	 tx113 = tx61; 		/* rule 129 */
	 tx78 = tx112*tx113; 		/* rule 130 */
	 tx79 = tx112; 		/* rule 131 */
	 tx80 = power2(tx63); 		/* rule 132 */
	 tx81 = tx75 + tx76; 		/* rule 133 */
	 tzz116 = tx78*tzz115; 		/* rule 134 */
	 tzz121 = -0.5*tzz116; 		/* rule 135 */
	 tx82 = tx80*tzz121; 		/* rule 136 */
	 tzz117 = tx77*tx79; 		/* rule 137 */
	 tzz119 = ka*tzz117; 		/* rule 138 */
	 tzz120 = 0.5*tzz119; 		/* rule 139 */
	 tx83 = tx80*tzz120; 		/* rule 140 */
	 tzz118 = tx113*tzz115; 		/* rule 141 */
	 tx84 = tx81*tzz118; 		/* rule 142 */
	 dhx1x1 = tx82 + tx83 + tx84; 		/* rule 143 */
	POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx85 = power2(tx65); 		/* rule 145 */
	 tx86 = 2.*tx37; 		/* rule 146 */
	 tx87 = 2.*tx85; 		/* rule 147 */
	 tx88 = power2(tx68); 		/* rule 148 */
	 tx89 = tx86 + tx87; 		/* rule 149 */
	 tx90 = tx88*tzz121; 		/* rule 150 */
	 tx91 = tx88*tzz120; 		/* rule 151 */
	 tx92 = tx89*tzz118; 		/* rule 152 */
	 dhy1y1 = tx90 + tx91 + tx92; 		/* rule 153 */
	POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx93 = power2(tx70); 		/* rule 155 */
	 tx94 = 2.*tx38; 		/* rule 156 */
	 tx95 = 2.*tx93; 		/* rule 157 */
	 tx96 = power2(tx73); 		/* rule 158 */
	 tx97 = tx94 + tx95; 		/* rule 159 */
	 tx98 = tx96*tzz121; 		/* rule 160 */
	 tx99 = tx96*tzz120; 		/* rule 161 */
	 tx100 = tx97*tzz118; 		/* rule 162 */
	 dhz1z1 = tx100 + tx98 + tx99; 		/* rule 163 */
	POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	#ifdef POINT_TO_LINE_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tzz123 = 2.*tzz118; 		/* rule 167 */
	 tx101 = tx31*tx57*tzz123; 		/* rule 168 */
	 tzz126 = tx63*tx68; 		/* rule 169 */
	 tx102 = tzz121*tzz126; 		/* rule 170 */
	 tx103 = tzz120*tzz126; 		/* rule 171 */
	 ohx1y1 = tx101 + tx102 + tx103; 		/* rule 172 */
	POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tx104 = tx33*tx70*tzz123; 		/* rule 174 */
	 tzz124 = tx73*tzz121; 		/* rule 175 */
	 tx105 = tx63*tzz124; 		/* rule 176 */
	 tzz125 = tx73*tzz120; 		/* rule 177 */
	 tx106 = tx63*tzz125; 		/* rule 178 */
	 ohx1z1 = tx104 + tx105 + tx106; 		/* rule 179 */
	POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx107 = tx32*tx65*tzz123; 		/* rule 181 */
	 tx108 = tx68*tzz124; 		/* rule 182 */
	 tx109 = tx68*tzz125; 		/* rule 183 */
	 ohy1z1 = tx107 + tx108 + tx109; 		/* rule 184 */
	POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	} /*if calcOffDiagonalHessian */ 
	#endif /* POINT_TO_LINE_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* POINT_TO_LINE_RESTRAINT_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* POINT_TO_LINE_RESTRAINT_CALC_FORCE ]*/
	POINT_TO_LINE_RESTRAINT_DEBUG_INTERACTIONS(I1);
	SKIP_term: (void)0;
