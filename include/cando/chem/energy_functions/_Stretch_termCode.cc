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
	 gx1 = -2.*kb*tx14*tx15*tx4; 		/* rule 32 */
	 fx1 = -gx1; 		/* rule 33 */
	STRETCH_FORCE_ACCUMULATE(I1, 0, fx1 );
	 gy1 = -2.*kb*tx14*tx15*tx5; 		/* rule 35 */
	 fy1 = -gy1; 		/* rule 36 */
	STRETCH_FORCE_ACCUMULATE(I1, 1, fy1 );
	 gz1 = -2.*kb*tx14*tx15*tx6; 		/* rule 38 */
	 fz1 = -gz1; 		/* rule 39 */
	STRETCH_FORCE_ACCUMULATE(I1, 2, fz1 );
	 gx2 = 2.*kb*tx14*tx15*tx4; 		/* rule 41 */
	 fx2 = -gx2; 		/* rule 42 */
	STRETCH_FORCE_ACCUMULATE(I2, 0, fx2 );
	 gy2 = 2.*kb*tx14*tx15*tx5; 		/* rule 44 */
	 fy2 = -gy2; 		/* rule 45 */
	STRETCH_FORCE_ACCUMULATE(I2, 1, fy2 );
	 gz2 = 2.*kb*tx14*tx15*tx6; 		/* rule 47 */
	 fz2 = -gz2; 		/* rule 48 */
	STRETCH_FORCE_ACCUMULATE(I2, 2, fz2 );
	#ifdef STRETCH_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx44 = tx10; 		/* rule 52 */
	 tx45 = tx15; 		/* rule 53 */
	 tx46 = reciprocal(tx44); 		/* rule 54 */
	 tx16 = tx45*tx46; 		/* rule 55 */
	 tx17 = tx46; 		/* rule 56 */
	 tx18 = 2.*kb*tx14*tx45; 		/* rule 57 */
	 tx19 = -2.*kb*tx14*tx16*tx7; 		/* rule 58 */
	 tx20 = 2.*kb*tx17*tx7; 		/* rule 59 */
	 dhx1x1 = tx18 + tx19 + tx20; 		/* rule 60 */
	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx21 = -2.*kb*tx14*tx16*tx8; 		/* rule 62 */
	 tx22 = 2.*kb*tx17*tx8; 		/* rule 63 */
	 dhy1y1 = tx18 + tx21 + tx22; 		/* rule 64 */
	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx23 = -2.*kb*tx14*tx16*tx9; 		/* rule 66 */
	 tx24 = 2.*kb*tx17*tx9; 		/* rule 67 */
	 dhz1z1 = tx18 + tx23 + tx24; 		/* rule 68 */
	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 dhx2x2 = dhx1x1; 		/* rule 70 */
	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 dhy2y2 = dhy1y1; 		/* rule 72 */
	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 dhz2z2 = dhz1z1; 		/* rule 74 */
	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	#ifdef STRETCH_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tx25 = -2.*kb*tx14*tx16*tx4*tx5; 		/* rule 78 */
	 tx26 = 2.*kb*tx17*tx4*tx5; 		/* rule 79 */
	 ohx1y1 = tx25 + tx26; 		/* rule 80 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tx27 = -2.*kb*tx14*tx16*tx4*tx6; 		/* rule 82 */
	 tx28 = 2.*kb*tx17*tx4*tx6; 		/* rule 83 */
	 ohx1z1 = tx27 + tx28; 		/* rule 84 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx29 = -2.*kb*tx14*tx45; 		/* rule 86 */
	 tx30 = 2.*kb*tx14*tx16*tx7; 		/* rule 87 */
	 tx31 = -2.*kb*tx17*tx7; 		/* rule 88 */
	 ohx1x2 = tx29 + tx30 + tx31; 		/* rule 89 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx32 = 2.*kb*tx14*tx16*tx4*tx5; 		/* rule 91 */
	 tx33 = -2.*kb*tx17*tx4*tx5; 		/* rule 92 */
	 ohx1y2 = tx32 + tx33; 		/* rule 93 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx34 = 2.*kb*tx14*tx16*tx4*tx6; 		/* rule 95 */
	 tx35 = -2.*kb*tx17*tx4*tx6; 		/* rule 96 */
	 ohx1z2 = tx34 + tx35; 		/* rule 97 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tx36 = -2.*kb*tx14*tx16*tx5*tx6; 		/* rule 99 */
	 tx37 = 2.*kb*tx17*tx5*tx6; 		/* rule 100 */
	 ohy1z1 = tx36 + tx37; 		/* rule 101 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 ohy1x2 = ohx1y2; 		/* rule 103 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx38 = 2.*kb*tx14*tx16*tx8; 		/* rule 105 */
	 tx39 = -2.*kb*tx17*tx8; 		/* rule 106 */
	 ohy1y2 = tx29 + tx38 + tx39; 		/* rule 107 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx40 = 2.*kb*tx14*tx16*tx5*tx6; 		/* rule 109 */
	 tx41 = -2.*kb*tx17*tx5*tx6; 		/* rule 110 */
	 ohy1z2 = tx40 + tx41; 		/* rule 111 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 ohz1x2 = ohx1z2; 		/* rule 113 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 ohz1y2 = ohy1z2; 		/* rule 115 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx42 = 2.*kb*tx14*tx16*tx9; 		/* rule 117 */
	 tx43 = -2.*kb*tx17*tx9; 		/* rule 118 */
	 ohz1z2 = tx29 + tx42 + tx43; 		/* rule 119 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 ohx2y2 = ohx1y1; 		/* rule 121 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 ohx2z2 = ohx1z1; 		/* rule 123 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 ohy2z2 = ohy1z1; 		/* rule 125 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	} /*if calcOffDiagonalHessian */ 
	#endif /* STRETCH_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* STRETCH_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* STRETCH_CALC_FORCE ]*/
	SKIP_term:  (void)0;
