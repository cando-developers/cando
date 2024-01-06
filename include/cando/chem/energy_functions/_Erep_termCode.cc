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
	 gx1 = crep*tx12*tx4; 		/* rule 29 */
	 fx1 = -gx1; 		/* rule 30 */
	EREP_FORCE_ACCUMULATE(I1, 0, fx1 );
	 gy1 = crep*tx12*tx5; 		/* rule 32 */
	 fy1 = -gy1; 		/* rule 33 */
	EREP_FORCE_ACCUMULATE(I1, 1, fy1 );
	 gz1 = crep*tx12*tx6; 		/* rule 35 */
	 fz1 = -gz1; 		/* rule 36 */
	EREP_FORCE_ACCUMULATE(I1, 2, fz1 );
	 gx2 = fx1; 		/* rule 38 */
	 fx2 = -gx2; 		/* rule 39 */
	EREP_FORCE_ACCUMULATE(I2, 0, fx2 );
	 gy2 = fy1; 		/* rule 41 */
	 fy2 = -gy2; 		/* rule 42 */
	EREP_FORCE_ACCUMULATE(I2, 1, fy2 );
	 gz2 = fz1; 		/* rule 44 */
	 fz2 = -gz2; 		/* rule 45 */
	EREP_FORCE_ACCUMULATE(I2, 2, fz2 );
	#ifdef EREP_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx13 = power2(tx12); 		/* rule 49 */
	 tx14 = crep*tx12; 		/* rule 50 */
	 tx15 = -2.*crep*tx13*tx7; 		/* rule 51 */
	 dhx1x1 = tx14 + tx15; 		/* rule 52 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx16 = -2.*crep*tx13*tx8; 		/* rule 54 */
	 dhy1y1 = tx14 + tx16; 		/* rule 55 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx17 = -2.*crep*tx13*tx9; 		/* rule 57 */
	 dhz1z1 = tx14 + tx17; 		/* rule 58 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 dhx2x2 = dhx1x1; 		/* rule 60 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 dhy2y2 = dhy1y1; 		/* rule 62 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 dhz2z2 = dhz1z1; 		/* rule 64 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	#ifdef EREP_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 ohx1y1 = -2.*crep*tx13*tx4*tx5; 		/* rule 68 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 ohx1z1 = -2.*crep*tx13*tx4*tx6; 		/* rule 70 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx18 = -tx14; 		/* rule 72 */
	 tx19 = 2.*crep*tx13*tx7; 		/* rule 73 */
	 ohx1x2 = tx18 + tx19; 		/* rule 74 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 ohx1y2 = 2.*crep*tx13*tx4*tx5; 		/* rule 76 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 ohx1z2 = 2.*crep*tx13*tx4*tx6; 		/* rule 78 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 ohy1z1 = -2.*crep*tx13*tx5*tx6; 		/* rule 80 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 ohy1x2 = ohx1y2; 		/* rule 82 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx20 = 2.*crep*tx13*tx8; 		/* rule 84 */
	 ohy1y2 = tx18 + tx20; 		/* rule 85 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 ohy1z2 = 2.*crep*tx13*tx5*tx6; 		/* rule 87 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 ohz1x2 = ohx1z2; 		/* rule 89 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 ohz1y2 = ohy1z2; 		/* rule 91 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx21 = 2.*crep*tx13*tx9; 		/* rule 93 */
	 ohz1z2 = tx18 + tx21; 		/* rule 94 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 ohx2y2 = ohx1y1; 		/* rule 96 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 ohx2z2 = ohx1z1; 		/* rule 98 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 ohy2z2 = ohy1z1; 		/* rule 100 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	} /*if calcOffDiagonalHessian */ 
	#endif /* EREP_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* EREP_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* EREP_CALC_FORCE ]*/
	EREP_DEBUG_INTERACTIONS(I1,I2);
	SKIP_term: (void)0;
