// C-code
	STRETCH_SET_PARAMETER(kb);
	STRETCH_SET_PARAMETER(r0);
	STRETCH_SET_PARAMETER(I1);
	STRETCH_SET_PARAMETER(I2);
	STRETCH_APPLY_ATOM_MASK(I1,I2);
	STRETCH_SET_POSITION(x1,I1,0);
	STRETCH_SET_POSITION(y1,I1,1);
	STRETCH_SET_POSITION(z1,I1,2);
	STRETCH_SET_POSITION(x2,I2,0);
	STRETCH_SET_POSITION(y2,I2,1);
	STRETCH_SET_POSITION(z2,I2,2);
	 tx1 = -x1; 		/* rule 12 */
	 tx2 = -y1; 		/* rule 13 */
	 tx3 = -z1; 		/* rule 14 */
	 tx4 = tx1 + x2; 		/* rule 15 */
	 tx5 = tx2 + y2; 		/* rule 16 */
	 tx6 = tx3 + z2; 		/* rule 17 */
	 tx7 = power2(tx4); 		/* rule 18 */
	 tx8 = power2(tx5); 		/* rule 19 */
	 tx9 = power2(tx6); 		/* rule 20 */
	 tx10 = tx7 + tx8 + tx9; 		/* rule 21 */
	 tx11 = -r0; 		/* rule 22 */
	 tx12 = mysqrt(tx10); 		/* rule 23 */
	 StretchDeviation = tx11 + tx12; 		/* rule 24 */
	 tx13 = power2(StretchDeviation); 		/* rule 25 */
	 Energy = kb*tx13; 		/* rule 26 */
	STRETCH_ENERGY_ACCUMULATE(Energy);
	#ifdef STRETCH_CALC_FORCE //[
	if ( calcForce ) {
	 tx14 = StretchDeviation; 		/* rule 30 */
	 tx15 = reciprocal(tx12); 		/* rule 31 */
	 tzz47 = kb*tx14; 		/* rule 32 */
	 tzz50 = tx15*tzz47; 		/* rule 33 */
	 tzz59 = -2.*tzz50; 		/* rule 34 */
	 gx1 = tx4*tzz59; 		/* rule 35 */
	 fx1 = -gx1; 		/* rule 36 */
	STRETCH_FORCE_ACCUMULATE(I1, 0, fx1 );
	 gy1 = tx5*tzz59; 		/* rule 38 */
	 fy1 = -gy1; 		/* rule 39 */
	STRETCH_FORCE_ACCUMULATE(I1, 1, fy1 );
	 gz1 = tx6*tzz59; 		/* rule 41 */
	 fz1 = -gz1; 		/* rule 42 */
	STRETCH_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tzz58 = 2.*tzz50; 		/* rule 44 */
	 gx2 = tx4*tzz58; 		/* rule 45 */
	 fx2 = -gx2; 		/* rule 46 */
	STRETCH_FORCE_ACCUMULATE(I2, 0, fx2 );
	 gy2 = tx5*tzz58; 		/* rule 48 */
	 fy2 = -gy2; 		/* rule 49 */
	STRETCH_FORCE_ACCUMULATE(I2, 1, fy2 );
	 gz2 = tx6*tzz58; 		/* rule 51 */
	 fz2 = -gz2; 		/* rule 52 */
	STRETCH_FORCE_ACCUMULATE(I2, 2, fz2 );
	#ifdef STRETCH_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx44 = tx10; 		/* rule 56 */
	 tx45 = tx15; 		/* rule 57 */
	 tx46 = reciprocal(tx44); 		/* rule 58 */
	 tx16 = tx45*tx46; 		/* rule 59 */
	 tx17 = tx46; 		/* rule 60 */
	 tx18 = 2.*tx45*tzz47; 		/* rule 61 */
	 tzz48 = tx16*tzz47; 		/* rule 62 */
	 tzz54 = -2.*tzz48; 		/* rule 63 */
	 tx19 = tx7*tzz54; 		/* rule 64 */
	 tzz49 = kb*tx17; 		/* rule 65 */
	 tzz51 = 2.*tzz49; 		/* rule 66 */
	 tx20 = tx7*tzz51; 		/* rule 67 */
	 dhx1x1 = tx18 + tx19 + tx20; 		/* rule 68 */
	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx21 = tx8*tzz54; 		/* rule 70 */
	 tx22 = tx8*tzz51; 		/* rule 71 */
	 dhy1y1 = tx18 + tx21 + tx22; 		/* rule 72 */
	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx23 = tx9*tzz54; 		/* rule 74 */
	 tx24 = tx9*tzz51; 		/* rule 75 */
	 dhz1z1 = tx18 + tx23 + tx24; 		/* rule 76 */
	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 dhx2x2 = dhx1x1; 		/* rule 78 */
	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 dhy2y2 = dhy1y1; 		/* rule 80 */
	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 dhz2z2 = dhz1z1; 		/* rule 82 */
	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	#ifdef STRETCH_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tzz57 = tx4*tx5; 		/* rule 86 */
	 tx25 = tzz54*tzz57; 		/* rule 87 */
	 tx26 = tzz51*tzz57; 		/* rule 88 */
	 ohx1y1 = tx25 + tx26; 		/* rule 89 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tzz56 = tx4*tx6; 		/* rule 91 */
	 tx27 = tzz54*tzz56; 		/* rule 92 */
	 tx28 = tzz51*tzz56; 		/* rule 93 */
	 ohx1z1 = tx27 + tx28; 		/* rule 94 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx29 = -2.*tx45*tzz47; 		/* rule 96 */
	 tzz52 = 2.*tzz48; 		/* rule 97 */
	 tx30 = tx7*tzz52; 		/* rule 98 */
	 tzz53 = -2.*tzz49; 		/* rule 99 */
	 tx31 = tx7*tzz53; 		/* rule 100 */
	 ohx1x2 = tx29 + tx30 + tx31; 		/* rule 101 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx32 = tzz52*tzz57; 		/* rule 103 */
	 tx33 = tzz53*tzz57; 		/* rule 104 */
	 ohx1y2 = tx32 + tx33; 		/* rule 105 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx34 = tzz52*tzz56; 		/* rule 107 */
	 tx35 = tzz53*tzz56; 		/* rule 108 */
	 ohx1z2 = tx34 + tx35; 		/* rule 109 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tzz55 = tx5*tx6; 		/* rule 111 */
	 tx36 = tzz54*tzz55; 		/* rule 112 */
	 tx37 = tzz51*tzz55; 		/* rule 113 */
	 ohy1z1 = tx36 + tx37; 		/* rule 114 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 ohy1x2 = ohx1y2; 		/* rule 116 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx38 = tx8*tzz52; 		/* rule 118 */
	 tx39 = tx8*tzz53; 		/* rule 119 */
	 ohy1y2 = tx29 + tx38 + tx39; 		/* rule 120 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx40 = tzz52*tzz55; 		/* rule 122 */
	 tx41 = tzz53*tzz55; 		/* rule 123 */
	 ohy1z2 = tx40 + tx41; 		/* rule 124 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 ohz1x2 = ohx1z2; 		/* rule 126 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 ohz1y2 = ohy1z2; 		/* rule 128 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx42 = tx9*tzz52; 		/* rule 130 */
	 tx43 = tx9*tzz53; 		/* rule 131 */
	 ohz1z2 = tx29 + tx42 + tx43; 		/* rule 132 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 ohx2y2 = ohx1y1; 		/* rule 134 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 ohx2z2 = ohx1z1; 		/* rule 136 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 ohy2z2 = ohy1z1; 		/* rule 138 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	} /*if calcOffDiagonalHessian */ 
	#endif /* STRETCH_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* STRETCH_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* STRETCH_CALC_FORCE ]*/
	STRETCH_DEBUG_INTERACTIONS(I1,I2);
	SKIP_term: (void)0;
