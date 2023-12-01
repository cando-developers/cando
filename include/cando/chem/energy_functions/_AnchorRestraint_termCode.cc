// C-code
	ANCHOR_RESTRAINT_SET_PARAMETER(ka);
	ANCHOR_RESTRAINT_SET_PARAMETER(xa);
	ANCHOR_RESTRAINT_SET_PARAMETER(ya);
	ANCHOR_RESTRAINT_SET_PARAMETER(za);
	ANCHOR_RESTRAINT_SET_PARAMETER(I1);
	ANCHOR_RESTRAINT_APPLY_ATOM_MASK(I1);
	ANCHOR_RESTRAINT_SET_POSITION(x1,I1,0);
	ANCHOR_RESTRAINT_SET_POSITION(y1,I1,1);
	ANCHOR_RESTRAINT_SET_POSITION(z1,I1,2);
	 tx1 = -xa; 		/* rule 10 */
	 tx2 = -ya; 		/* rule 11 */
	 tx3 = -za; 		/* rule 12 */
	 tx4 = tx1 + x1; 		/* rule 13 */
	 tx5 = tx2 + y1; 		/* rule 14 */
	 tx6 = tx3 + z1; 		/* rule 15 */
	 tx7 = power2(tx4); 		/* rule 16 */
	 tx8 = power2(tx5); 		/* rule 17 */
	 tx9 = power2(tx6); 		/* rule 18 */
	 AnchorDeviation = tx7 + tx8 + tx9; 		/* rule 19 */
	 Energy = AnchorDeviation*ka; 		/* rule 20 */
	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE(Energy);
	#ifdef ANCHOR_RESTRAINT_CALC_FORCE //[
	if ( calcForce ) {
	 gx1 = 2.*ka*tx4; 		/* rule 24 */
	 fx1 = -gx1; 		/* rule 25 */
	ANCHOR_RESTRAINT_FORCE_ACCUMULATE(I1, 0, fx1 );
	 gy1 = 2.*ka*tx5; 		/* rule 27 */
	 fy1 = -gy1; 		/* rule 28 */
	ANCHOR_RESTRAINT_FORCE_ACCUMULATE(I1, 1, fy1 );
	 gz1 = 2.*ka*tx6; 		/* rule 30 */
	 fz1 = -gz1; 		/* rule 31 */
	ANCHOR_RESTRAINT_FORCE_ACCUMULATE(I1, 2, fz1 );
	#ifdef ANCHOR_RESTRAINT_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 dhx1x1 = 2.*ka; 		/* rule 35 */
	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 dhy1y1 = dhx1x1; 		/* rule 37 */
	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 dhz1z1 = dhy1y1; 		/* rule 39 */
	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	#ifdef ANCHOR_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 ohx1y1 = 0; 		/* rule 43 */
	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 ohx1z1 = ohx1y1; 		/* rule 45 */
	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 ohy1z1 = ohx1z1; 		/* rule 47 */
	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	} /*if calcOffDiagonalHessian */ 
	#endif /* ANCHOR_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* ANCHOR_RESTRAINT_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* ANCHOR_RESTRAINT_CALC_FORCE ]*/
	SKIP_term:  (void)0;
