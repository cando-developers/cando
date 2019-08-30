// C-code
	EREP_SET_PARAMETER(I1);
	EREP_SET_PARAMETER(I2);
	EREP_SET_PARAMETER(crep);
	EREP_SET_POSITION(x1,I1,0);
	EREP_SET_POSITION(y1,I1,1);
	EREP_SET_POSITION(z1,I1,2);
	EREP_SET_POSITION(x2,I2,0);
	EREP_SET_POSITION(y2,I2,1);
	EREP_SET_POSITION(z2,I2,2);
	 tx22 = -x2; 		/* rule 10 */
	 tx23 = -y2; 		/* rule 11 */
	 tx24 = -z2; 		/* rule 12 */
	 tx25 = tx22 + x1; 		/* rule 13 */
	 tx26 = tx23 + y1; 		/* rule 14 */
	 tx27 = tx24 + z1; 		/* rule 15 */
	 tx28 = power2(tx25); 		/* rule 16 */
	 tx29 = power2(tx26); 		/* rule 17 */
	 tx30 = power2(tx27); 		/* rule 18 */
	 tx31 = tx28 + tx29 + tx30; 		/* rule 19 */
	 ErepDistance = mysqrt(tx31); 		/* rule 20 */
	MAYBE_BAIL(ErepDistance);
	 tx32 = reciprocal(ErepDistance); 		/* rule 22 */
	 Erep = crep*tx32; 		/* rule 23 */
	EREP_ENERGY_ACCUMULATE(Erep);
	#ifdef EREP_CALC_FORCE //[
	if ( calcForce ) {
	 tx43 = tx31; 		/* rule 27 */
	 tx44 = tx32; 		/* rule 28 */
	 tx45 = reciprocal(tx43); 		/* rule 29 */
	 tx33 = tx44*tx45; 		/* rule 30 */
	 tzz49 = crep*tx33; 		/* rule 31 */
	 tzz52 = -tzz49; 		/* rule 32 */
	 gx1 = tx25*tzz52; 		/* rule 33 */
	 fx1 = -gx1; 		/* rule 34 */
	EREP_FORCE_ACCUMULATE(I1, 0, fx1 );
	 gy1 = tx26*tzz52; 		/* rule 36 */
	 fy1 = -gy1; 		/* rule 37 */
	EREP_FORCE_ACCUMULATE(I1, 1, fy1 );
	 gz1 = tx27*tzz52; 		/* rule 39 */
	 fz1 = -gz1; 		/* rule 40 */
	EREP_FORCE_ACCUMULATE(I1, 2, fz1 );
	 gx2 = tx25*tzz49; 		/* rule 42 */
	 fx2 = -gx2; 		/* rule 43 */
	EREP_FORCE_ACCUMULATE(I2, 0, fx2 );
	 gy2 = tx26*tzz49; 		/* rule 45 */
	 fy2 = -gy2; 		/* rule 46 */
	EREP_FORCE_ACCUMULATE(I2, 1, fy2 );
	 gz2 = tx27*tzz49; 		/* rule 48 */
	 fz2 = -gz2; 		/* rule 49 */
	EREP_FORCE_ACCUMULATE(I2, 2, fz2 );
	#ifdef EREP_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx46 = power2(tx43); 		/* rule 53 */
	 tx47 = reciprocal(tx46); 		/* rule 54 */
	 tx34 = tx44*tx47; 		/* rule 55 */
	 tx35 = tzz52; 		/* rule 56 */
	 tzz48 = crep*tx34; 		/* rule 57 */
	 tzz50 = 3*tzz48; 		/* rule 58 */
	 tx36 = tx28*tzz50; 		/* rule 59 */
	 dhx1x1 = tx35 + tx36; 		/* rule 60 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx37 = tx29*tzz50; 		/* rule 62 */
	 dhy1y1 = tx35 + tx37; 		/* rule 63 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx38 = tx30*tzz50; 		/* rule 65 */
	 dhz1z1 = tx35 + tx38; 		/* rule 66 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 dhx2x2 = dhx1x1; 		/* rule 68 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 dhy2y2 = dhy1y1; 		/* rule 70 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 dhz2z2 = dhz1z1; 		/* rule 72 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	#ifdef EREP_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 ohx1y1 = tx25*tx26*tzz50; 		/* rule 76 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tzz53 = tx27*tzz50; 		/* rule 78 */
	 ohx1z1 = tx25*tzz53; 		/* rule 79 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx39 = tzz49; 		/* rule 81 */
	 tzz51 = -3*tzz48; 		/* rule 82 */
	 tx40 = tx28*tzz51; 		/* rule 83 */
	 ohx1x2 = tx39 + tx40; 		/* rule 84 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tzz54 = tx26*tzz51; 		/* rule 86 */
	 ohx1y2 = tx25*tzz54; 		/* rule 87 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 ohx1z2 = tx25*tx27*tzz51; 		/* rule 89 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 ohy1z1 = tx26*tzz53; 		/* rule 91 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 ohy1x2 = ohx1y2; 		/* rule 93 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx41 = tx29*tzz51; 		/* rule 95 */
	 ohy1y2 = tx39 + tx41; 		/* rule 96 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 ohy1z2 = tx27*tzz54; 		/* rule 98 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 ohz1x2 = ohx1z2; 		/* rule 100 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 ohz1y2 = ohy1z2; 		/* rule 102 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx42 = tx30*tzz51; 		/* rule 104 */
	 ohz1z2 = tx39 + tx42; 		/* rule 105 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 ohx2y2 = ohx1y1; 		/* rule 107 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 ohx2z2 = ohx1z1; 		/* rule 109 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 ohy2z2 = ohy1z1; 		/* rule 111 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	} /*if calcOffDiagonalHessian */ 
	#endif /* EREP_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* EREP_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* EREP_CALC_FORCE ]*/
