// C-code
	OOZP_SET_PARAMETER(kb);
	OOZP_SET_PARAMETER(za);
	OOZP_SET_PARAMETER(I1);
	OOZP_SET_PARAMETER(I2);
	OOZP_APPLY_ATOM_MASK(I1,I2);
	OOZP_SET_POSITION(x1,I1,0);
	OOZP_SET_POSITION(y1,I1,1);
	OOZP_SET_POSITION(z1,I1,2);
	 tx1 = -za; 		/* rule 9 */
	 tx2 = tx1 + z1; 		/* rule 10 */
	 tx3 = power2(tx2); 		/* rule 11 */
	 Energy = kb*tx3; 		/* rule 12 */
	OOZP_ENERGY_ACCUMULATE(Energy);
	#ifdef OOZP_CALC_FORCE //[
	if ( calcForce ) {
	 gx1 = 0; 		/* rule 16 */
	 fx1 = -gx1; 		/* rule 17 */
	OOZP_FORCE_ACCUMULATE(I1, 0, fx1 );
	 gy1 = gx1; 		/* rule 19 */
	 fy1 = -gy1; 		/* rule 20 */
	OOZP_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx4 = tx2; 		/* rule 22 */
	 gz1 = 2.*kb*tx4; 		/* rule 23 */
	 fz1 = -gz1; 		/* rule 24 */
	OOZP_FORCE_ACCUMULATE(I1, 2, fz1 );
	#ifdef OOZP_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 dhx1x1 = gy1; 		/* rule 28 */
	OOZP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 dhy1y1 = dhx1x1; 		/* rule 30 */
	OOZP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 dhz1z1 = 2.*kb; 		/* rule 32 */
	OOZP_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	#ifdef OOZP_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 ohx1y1 = dhy1y1; 		/* rule 36 */
	OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 ohx1z1 = ohx1y1; 		/* rule 38 */
	OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 ohy1z1 = ohx1z1; 		/* rule 40 */
	OOZP_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	} /*if calcOffDiagonalHessian */ 
	#endif /* OOZP_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* OOZP_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* OOZP_CALC_FORCE ]*/
	OOZP_DEBUG_INTERACTIONS(I1,I2);
	SKIP_term: (void)0;
