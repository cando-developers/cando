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
	 tx32 = Log(ErepDistance); 		/* rule 22 */
	 Erep = crep*tx32; 		/* rule 23 */
	EREP_ENERGY_ACCUMULATE(Erep);
	#ifdef EREP_CALC_FORCE //[
	if ( calcForce ) {
	 tx33 = reciprocal(tx31); 		/* rule 27 */
	 tzz46 = crep*tx33; 		/* rule 28 */
	 gx1 = tx25*tzz46; 		/* rule 29 */
	 fx1 = -gx1; 		/* rule 30 */
	EREP_FORCE_ACCUMULATE(I1, 0, fx1 );
	 gy1 = tx26*tzz46; 		/* rule 32 */
	 fy1 = -gy1; 		/* rule 33 */
	EREP_FORCE_ACCUMULATE(I1, 1, fy1 );
	 gz1 = tx27*tzz46; 		/* rule 35 */
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
	 tx34 = power2(tx33); 		/* rule 49 */
	 tx35 = tzz46; 		/* rule 50 */
	 tzz43 = crep*tx34; 		/* rule 51 */
	 tzz45 = -2.*tzz43; 		/* rule 52 */
	 tx36 = tx28*tzz45; 		/* rule 53 */
	 dhx1x1 = tx35 + tx36; 		/* rule 54 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx37 = tx29*tzz45; 		/* rule 56 */
	 dhy1y1 = tx35 + tx37; 		/* rule 57 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx38 = tx30*tzz45; 		/* rule 59 */
	 dhz1z1 = tx35 + tx38; 		/* rule 60 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 dhx2x2 = dhx1x1; 		/* rule 62 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 dhy2y2 = dhy1y1; 		/* rule 64 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 dhz2z2 = dhz1z1; 		/* rule 66 */
	EREP_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	#ifdef EREP_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tzz48 = tx26*tzz45; 		/* rule 70 */
	 ohx1y1 = tx25*tzz48; 		/* rule 71 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 ohx1z1 = tx25*tx27*tzz45; 		/* rule 73 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tzz44 = 2.*tzz43; 		/* rule 75 */
	 tx39 = tx28*tzz44; 		/* rule 76 */
	 tx40 = -tx35; 		/* rule 77 */
	 ohx1x2 = tx39 + tx40; 		/* rule 78 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 ohx1y2 = tx25*tx26*tzz44; 		/* rule 80 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tzz47 = tx27*tzz44; 		/* rule 82 */
	 ohx1z2 = tx25*tzz47; 		/* rule 83 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 ohy1z1 = tx27*tzz48; 		/* rule 85 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 ohy1x2 = ohx1y2; 		/* rule 87 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx41 = tx29*tzz44; 		/* rule 89 */
	 ohy1y2 = tx40 + tx41; 		/* rule 90 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 ohy1z2 = tx26*tzz47; 		/* rule 92 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 ohz1x2 = ohx1z2; 		/* rule 94 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 ohz1y2 = ohy1z2; 		/* rule 96 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx42 = tx30*tzz44; 		/* rule 98 */
	 ohz1z2 = tx40 + tx42; 		/* rule 99 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 ohx2y2 = ohx1y1; 		/* rule 101 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 ohx2z2 = ohx1z1; 		/* rule 103 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 ohy2z2 = ohy1z1; 		/* rule 105 */
	EREP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	} /*if calcOffDiagonalHessian */ 
	#endif /* EREP_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* EREP_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* EREP_CALC_FORCE ]*/
