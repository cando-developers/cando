// C-code
	OOZP_SET_PARAMETER(kb);
	OOZP_SET_PARAMETER(za);
	OOZP_SET_PARAMETER(I1);
	OOZP_SET_POSITION(x1,I1,0);
	OOZP_SET_POSITION(y1,I1,1);
	OOZP_SET_POSITION(z1,I1,2);
	 tx5 = -za; 		/* rule 8 */
	 tx6 = tx5 + z1; 		/* rule 9 */
	 tx7 = power2(tx6); 		/* rule 10 */
	 Energy = kb*tx7; 		/* rule 11 */
	OOZP_ENERGY_ACCUMULATE(Energy);
	#ifdef OOZP_CALC_FORCE //[
	if ( calcForce ) {
	 gx1 = 0; 		/* rule 15 */
	 fx1 = -gx1; 		/* rule 16 */
	OOZP_FORCE_ACCUMULATE(I1, 0, fx1 );
	 gy1 = gx1; 		/* rule 18 */
	 fy1 = -gy1; 		/* rule 19 */
	OOZP_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx8 = tx6; 		/* rule 21 */
	 gz1 = 2.*kb*tx8; 		/* rule 22 */
	 fz1 = -gz1; 		/* rule 23 */
	OOZP_FORCE_ACCUMULATE(I1, 2, fz1 );
	#ifdef OOZP_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 dhx1x1 = gy1; 		/* rule 27 */
	OOZP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 dhy1y1 = dhx1x1; 		/* rule 29 */
	OOZP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 dhz1z1 = 2.*kb; 		/* rule 31 */
	OOZP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	#ifdef OOZP_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 ohx1y1 = dhy1y1; 		/* rule 35 */
	OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 ohx1z1 = ohx1y1; 		/* rule 37 */
	OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 ohy1z1 = ohx1z1; 		/* rule 39 */
	OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	} /*if calcOffDiagonalHessian */ 
	#endif /* OOZP_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* OOZP_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* OOZP_CALC_FORCE ]*/
