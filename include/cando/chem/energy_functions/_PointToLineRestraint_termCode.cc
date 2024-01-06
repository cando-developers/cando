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
	 tx2 = tx1*y1; 		/* rule 15 */
	 tx3 = xb*y1; 		/* rule 16 */
	 tx4 = -ya; 		/* rule 17 */
	 tx5 = x1*ya; 		/* rule 18 */
	 tx6 = tx4*xb; 		/* rule 19 */
	 tx7 = -(x1*yb); 		/* rule 20 */
	 tx8 = xa*yb; 		/* rule 21 */
	 tx9 = xa*z1; 		/* rule 22 */
	 tx10 = -(xb*z1); 		/* rule 23 */
	 tx11 = tx4*z1; 		/* rule 24 */
	 tx12 = yb*z1; 		/* rule 25 */
	 tx13 = -za; 		/* rule 26 */
	 tx14 = tx13*x1; 		/* rule 27 */
	 tx15 = xb*za; 		/* rule 28 */
	 tx16 = y1*za; 		/* rule 29 */
	 tx17 = tx13*yb; 		/* rule 30 */
	 tx18 = x1*zb; 		/* rule 31 */
	 tx19 = tx1*zb; 		/* rule 32 */
	 tx20 = -(y1*zb); 		/* rule 33 */
	 tx21 = ya*zb; 		/* rule 34 */
	 tx22 = tx14; 		/* rule 35 */
	 tx23 = tx6; 		/* rule 36 */
	 tx24 = tx2; 		/* rule 37 */
	 tx25 = tx17; 		/* rule 38 */
	 tx26 = tx11; 		/* rule 39 */
	 tx27 = tx19; 		/* rule 40 */
	 tx28 = tx12 + tx16 + tx20 + tx21 + tx25 + tx26; 		/* rule 41 */
	 tx29 = tx23 + tx24 + tx3 + tx5 + tx7 + tx8; 		/* rule 42 */
	 tx30 = tx10 + tx15 + tx18 + tx22 + tx27 + tx9; 		/* rule 43 */
	 tx31 = tx1 + xb; 		/* rule 44 */
	 tx32 = tx4 + yb; 		/* rule 45 */
	 tx33 = tx13 + zb; 		/* rule 46 */
	 tx34 = power2(tx28); 		/* rule 47 */
	 tx35 = power2(tx29); 		/* rule 48 */
	 tx36 = power2(tx30); 		/* rule 49 */
	 tx37 = power2(tx31); 		/* rule 50 */
	 tx38 = power2(tx32); 		/* rule 51 */
	 tx39 = power2(tx33); 		/* rule 52 */
	 tx40 = tx34 + tx35 + tx36; 		/* rule 53 */
	 tx41 = tx37 + tx38 + tx39; 		/* rule 54 */
	 tx42 = mysqrt(tx40); 		/* rule 55 */
	 tx100 = mysqrt(tx41); 		/* rule 56 */
	 tx43 = reciprocal(tx100); 		/* rule 57 */
	 tx44 = -ra; 		/* rule 58 */
	 tx45 = tx42*tx43; 		/* rule 59 */
	 tx46 = tx44 + tx45; 		/* rule 60 */
	 tx47 = power2(tx46); 		/* rule 61 */
	 Energy = ka*tx47; 		/* rule 62 */
	POINT_TO_LINE_RESTRAINT_ENERGY_ACCUMULATE(Energy);
	#ifdef POINT_TO_LINE_RESTRAINT_CALC_FORCE //[
	if ( calcForce ) {
	 tx48 = -yb; 		/* rule 66 */
	 tx49 = tx48 + ya; 		/* rule 67 */
	 tx50 = 2.*tx30*tx33; 		/* rule 68 */
	 tx51 = 2.*tx29*tx49; 		/* rule 69 */
	 tx52 = reciprocal(tx42); 		/* rule 70 */
	 tx53 = tx50 + tx51; 		/* rule 71 */
	 gx1 = ka*tx43*tx46*tx52*tx53; 		/* rule 72 */
	 fx1 = -gx1; 		/* rule 73 */
	POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx54 = -zb; 		/* rule 75 */
	 tx55 = tx54 + za; 		/* rule 76 */
	 tx56 = 2.*tx29*tx31; 		/* rule 77 */
	 tx57 = 2.*tx28*tx55; 		/* rule 78 */
	 tx58 = tx56 + tx57; 		/* rule 79 */
	 gy1 = ka*tx43*tx46*tx52*tx58; 		/* rule 80 */
	 fy1 = -gy1; 		/* rule 81 */
	POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx59 = -xb; 		/* rule 83 */
	 tx60 = tx59 + xa; 		/* rule 84 */
	 tx61 = 2.*tx28*tx32; 		/* rule 85 */
	 tx62 = 2.*tx30*tx60; 		/* rule 86 */
	 tx63 = tx61 + tx62; 		/* rule 87 */
	 gz1 = ka*tx43*tx46*tx52*tx63; 		/* rule 88 */
	 fz1 = -gz1; 		/* rule 89 */
	POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 2, fz1 );
	#ifdef POINT_TO_LINE_RESTRAINT_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx64 = power2(tx49); 		/* rule 93 */
	 tx65 = 2.*tx39; 		/* rule 94 */
	 tx66 = 2.*tx64; 		/* rule 95 */
	 tx101 = tx40; 		/* rule 96 */
	 tx102 = reciprocal(tx101); 		/* rule 97 */
	 tx103 = tx52; 		/* rule 98 */
	 tx67 = tx102*tx103; 		/* rule 99 */
	 tx68 = tx102; 		/* rule 100 */
	 tx69 = reciprocal(tx41); 		/* rule 101 */
	 tx70 = power2(tx53); 		/* rule 102 */
	 tx71 = tx65 + tx66; 		/* rule 103 */
	 tx72 = -0.5*ka*tx43*tx46*tx67*tx70; 		/* rule 104 */
	 tx73 = 0.5*ka*tx68*tx69*tx70; 		/* rule 105 */
	 tx74 = ka*tx103*tx43*tx46*tx71; 		/* rule 106 */
	 dhx1x1 = tx72 + tx73 + tx74; 		/* rule 107 */
	POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx75 = power2(tx55); 		/* rule 109 */
	 tx76 = 2.*tx37; 		/* rule 110 */
	 tx77 = 2.*tx75; 		/* rule 111 */
	 tx78 = power2(tx58); 		/* rule 112 */
	 tx79 = tx76 + tx77; 		/* rule 113 */
	 tx80 = -0.5*ka*tx43*tx46*tx67*tx78; 		/* rule 114 */
	 tx81 = 0.5*ka*tx68*tx69*tx78; 		/* rule 115 */
	 tx82 = ka*tx103*tx43*tx46*tx79; 		/* rule 116 */
	 dhy1y1 = tx80 + tx81 + tx82; 		/* rule 117 */
	POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx83 = power2(tx60); 		/* rule 119 */
	 tx84 = 2.*tx38; 		/* rule 120 */
	 tx85 = 2.*tx83; 		/* rule 121 */
	 tx86 = power2(tx63); 		/* rule 122 */
	 tx87 = tx84 + tx85; 		/* rule 123 */
	 tx88 = -0.5*ka*tx43*tx46*tx67*tx86; 		/* rule 124 */
	 tx89 = 0.5*ka*tx68*tx69*tx86; 		/* rule 125 */
	 tx90 = ka*tx103*tx43*tx46*tx87; 		/* rule 126 */
	 dhz1z1 = tx88 + tx89 + tx90; 		/* rule 127 */
	POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	#ifdef POINT_TO_LINE_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tx91 = 2.*ka*tx103*tx31*tx43*tx46*tx49; 		/* rule 131 */
	 tx92 = -0.5*ka*tx43*tx46*tx53*tx58*tx67; 		/* rule 132 */
	 tx93 = 0.5*ka*tx53*tx58*tx68*tx69; 		/* rule 133 */
	 ohx1y1 = tx91 + tx92 + tx93; 		/* rule 134 */
	POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tx94 = 2.*ka*tx103*tx33*tx43*tx46*tx60; 		/* rule 136 */
	 tx95 = -0.5*ka*tx43*tx46*tx53*tx63*tx67; 		/* rule 137 */
	 tx96 = 0.5*ka*tx53*tx63*tx68*tx69; 		/* rule 138 */
	 ohx1z1 = tx94 + tx95 + tx96; 		/* rule 139 */
	POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx97 = 2.*ka*tx103*tx32*tx43*tx46*tx55; 		/* rule 141 */
	 tx98 = -0.5*ka*tx43*tx46*tx58*tx63*tx67; 		/* rule 142 */
	 tx99 = 0.5*ka*tx58*tx63*tx68*tx69; 		/* rule 143 */
	 ohy1z1 = tx97 + tx98 + tx99; 		/* rule 144 */
	POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	} /*if calcOffDiagonalHessian */ 
	#endif /* POINT_TO_LINE_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* POINT_TO_LINE_RESTRAINT_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* POINT_TO_LINE_RESTRAINT_CALC_FORCE ]*/
	POINT_TO_LINE_RESTRAINT_DEBUG_INTERACTIONS(I1);
	SKIP_term: (void)0;
