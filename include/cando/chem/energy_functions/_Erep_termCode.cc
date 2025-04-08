// C-code
	EREP_SET_PARAMETER(crep);
	EREP_SET_PARAMETER(I1);
	EREP_SET_PARAMETER(I2);
	EREP_APPLY_ATOM_MASK(I1,I2);
	EREP_SET_POSITION(x1,I1,0);
	EREP_SET_POSITION(y1,I1,1);
	EREP_SET_POSITION(z1,I1,2);
	EREP_SET_POSITION(x2,I2,0);
	EREP_SET_POSITION(y2,I2,1);
	EREP_SET_POSITION(z2,I2,2);
	 tx1 = -x2; 		/* rule 11 */
	 tx2 = -y2; 		/* rule 12 */
	 tx3 = -z2; 		/* rule 13 */
	 tx4 = tx1 + x1; 		/* rule 14 */
	 tx5 = tx2 + y1; 		/* rule 15 */
	 tx6 = tx3 + z1; 		/* rule 16 */
	 tx7 = power2(tx4); 		/* rule 17 */
	 tx8 = power2(tx5); 		/* rule 18 */
	 tx9 = power2(tx6); 		/* rule 19 */
	 tx10 = tx7 + tx8 + tx9; 		/* rule 20 */
	 ErepDistance = mysqrt(tx10); 		/* rule 21 */
	MAYBE_BAIL(ErepDistance);
	 tx11 = Log(ErepDistance); 		/* rule 23 */
	 Erep = crep*tx11; 		/* rule 24 */
	EREP_ENERGY_ACCUMULATE(Erep);
	#ifdef EREP_CALC_FORCE //[
	if ( calcForce ) {
	 tx12 = reciprocal(tx10); 		/* rule 28 */
	 tzz23 = crep*tx12; 		/* rule 29 */
	 gx1 = tx4*tzz23; 		/* rule 30 */
	 fx1 = -gx1; 		/* rule 31 */
	EREP_FORCE_ACCUMULATE(I1, 0, fx1 );
	 gy1 = tx5*tzz23; 		/* rule 33 */
	 fy1 = -gy1; 		/* rule 34 */
	EREP_FORCE_ACCUMULATE(I1, 1, fy1 );
	 gz1 = tx6*tzz23; 		/* rule 36 */
	 fz1 = -gz1; 		/* rule 37 */
	EREP_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tzz26 = -tzz23; 		/* rule 39 */
	 gx2 = tx4*tzz26; 		/* rule 40 */
	 fx2 = -gx2; 		/* rule 41 */
	EREP_FORCE_ACCUMULATE(I2, 0, fx2 );
	 gy2 = tx5*tzz26; 		/* rule 43 */
	 fy2 = -gy2; 		/* rule 44 */
	EREP_FORCE_ACCUMULATE(I2, 1, fy2 );
	 gz2 = tx6*tzz26; 		/* rule 46 */
	 fz2 = -gz2; 		/* rule 47 */
	EREP_FORCE_ACCUMULATE(I2, 2, fz2 );
	#ifdef EREP_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx13 = power2(tx12); 		/* rule 51 */
	 tx14 = tzz23; 		/* rule 52 */
	 tzz22 = crep*tx13; 		/* rule 53 */
	 tzz25 = -2.*tzz22; 		/* rule 54 */
	 tx15 = tx7*tzz25; 		/* rule 55 */
	 dhx1x1 = tx14 + tx15; 		/* rule 56 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx16 = tx8*tzz25; 		/* rule 58 */
	 dhy1y1 = tx14 + tx16; 		/* rule 59 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx17 = tx9*tzz25; 		/* rule 61 */
	 dhz1z1 = tx14 + tx17; 		/* rule 62 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 dhx2x2 = dhx1x1; 		/* rule 64 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 dhy2y2 = dhy1y1; 		/* rule 66 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 dhz2z2 = dhz1z1; 		/* rule 68 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	#ifdef EREP_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tzz29 = tx4*tx5; 		/* rule 72 */
	 ohx1y1 = tzz25*tzz29; 		/* rule 73 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tzz27 = tx6*tzz25; 		/* rule 75 */
	 ohx1z1 = tx4*tzz27; 		/* rule 76 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx18 = tzz26; 		/* rule 78 */
	 tzz24 = 2.*tzz22; 		/* rule 79 */
	 tx19 = tx7*tzz24; 		/* rule 80 */
	 ohx1x2 = tx18 + tx19; 		/* rule 81 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 ohx1y2 = tzz24*tzz29; 		/* rule 83 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tzz28 = tx6*tzz24; 		/* rule 85 */
	 ohx1z2 = tx4*tzz28; 		/* rule 86 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 ohy1z1 = tx5*tzz27; 		/* rule 88 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 ohy1x2 = ohx1y2; 		/* rule 90 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx20 = tx8*tzz24; 		/* rule 92 */
	 ohy1y2 = tx18 + tx20; 		/* rule 93 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 ohy1z2 = tx5*tzz28; 		/* rule 95 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 ohz1x2 = ohx1z2; 		/* rule 97 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 ohz1y2 = ohy1z2; 		/* rule 99 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx21 = tx9*tzz24; 		/* rule 101 */
	 ohz1z2 = tx18 + tx21; 		/* rule 102 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 ohx2y2 = ohx1y1; 		/* rule 104 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 ohx2z2 = ohx1z1; 		/* rule 106 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 ohy2z2 = ohy1z1; 		/* rule 108 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	} /*if calcOffDiagonalHessian */ 
	#endif /* EREP_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* EREP_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* EREP_CALC_FORCE ]*/
	EREP_DEBUG_INTERACTIONS(I1,I2);
	SKIP_term: (void)0;
