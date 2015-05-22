// C-code
	ANCHOR_RESTRAINT_SET_PARAMETER(ka);
	ANCHOR_RESTRAINT_SET_PARAMETER(xa);
	ANCHOR_RESTRAINT_SET_PARAMETER(ya);
	ANCHOR_RESTRAINT_SET_PARAMETER(za);
	ANCHOR_RESTRAINT_SET_PARAMETER(I1);
	ANCHOR_RESTRAINT_SET_POSITION(x1,I1,0);
	ANCHOR_RESTRAINT_SET_POSITION(y1,I1,1);
	ANCHOR_RESTRAINT_SET_POSITION(z1,I1,2);
	 tx10 = -xa; 		/* rule 9 */
	 tx11 = -ya; 		/* rule 10 */
	 tx12 = -za; 		/* rule 11 */
	 tx13 = tx10 + x1; 		/* rule 12 */
	 tx14 = tx11 + y1; 		/* rule 13 */
	 tx15 = tx12 + z1; 		/* rule 14 */
	 tx16 = power2(tx13); 		/* rule 15 */
	 tx17 = power2(tx14); 		/* rule 16 */
	 tx18 = power2(tx15); 		/* rule 17 */
	 AnchorDeviation = tx16 + tx17 + tx18; 		/* rule 18 */
	 Energy = AnchorDeviation*ka; 		/* rule 19 */
	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE(Energy);
	#ifdef ANCHOR_RESTRAINT_CALC_FORCE //[
	if ( calcForce ) {
	 gx1 = 2.*ka*tx13; 		/* rule 23 */
	 fx1 = -gx1; 		/* rule 24 */
	ANCHOR_RESTRAINT_FORCE_ACCUMULATE(I1, 0, fx1 );
	 gy1 = 2.*ka*tx14; 		/* rule 26 */
	 fy1 = -gy1; 		/* rule 27 */
	ANCHOR_RESTRAINT_FORCE_ACCUMULATE(I1, 1, fy1 );
	 gz1 = 2.*ka*tx15; 		/* rule 29 */
	 fz1 = -gz1; 		/* rule 30 */
	ANCHOR_RESTRAINT_FORCE_ACCUMULATE(I1, 2, fz1 );
	#ifdef ANCHOR_RESTRAINT_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 dhx1x1 = 2.*ka; 		/* rule 34 */
	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 dhy1y1 = dhx1x1; 		/* rule 36 */
	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 dhz1z1 = dhy1y1; 		/* rule 38 */
	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	#ifdef ANCHOR_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 ohx1y1 = 0; 		/* rule 42 */
	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 ohx1z1 = ohx1y1; 		/* rule 44 */
	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 ohy1z1 = ohx1z1; 		/* rule 46 */
	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	} /*if calcOffDiagonalHessian */ 
	#endif /* ANCHOR_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* ANCHOR_RESTRAINT_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* ANCHOR_RESTRAINT_CALC_FORCE ]*/
