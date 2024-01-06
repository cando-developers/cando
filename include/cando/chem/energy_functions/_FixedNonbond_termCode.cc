// C-code
	// FNONBOND_SET_PARAMETER(xf);
	// FNONBOND_SET_PARAMETER(yf);
	// FNONBOND_SET_PARAMETER(zf);
	// FNONBOND_SET_PARAMETER(dQ1Q2);
	// FNONBOND_SET_PARAMETER(dA);
	// FNONBOND_SET_PARAMETER(dC);
	// FNONBOND_SET_PARAMETER(I1);
	FNONBOND_SET_POSITION(x1,I1,0);
	FNONBOND_SET_POSITION(y1,I1,1);
	FNONBOND_SET_POSITION(z1,I1,2);
	 tx1 = -xf; 		/* rule 11 */
	 tx2 = -yf; 		/* rule 12 */
	 tx3 = -zf; 		/* rule 13 */
	 tx4 = tx1 + x1; 		/* rule 14 */
	 tx5 = tx2 + y1; 		/* rule 15 */
	 tx6 = tx3 + z1; 		/* rule 16 */
	 tx7 = power2(tx4); 		/* rule 17 */
	 tx8 = power2(tx5); 		/* rule 18 */
	 tx9 = power2(tx6); 		/* rule 19 */
	 tx10 = tx7 + tx8 + tx9; 		/* rule 20 */
	 NonbondDistance = mysqrt(tx10); 		/* rule 21 */
	 tx63 = reciprocal(tx10); 		/* rule 22 */
	 tx61 = power2(tx63); 		/* rule 23 */
	 tx62 = tx63; 		/* rule 24 */
	 tx52 = tx61*tx62; 		/* rule 25 */
	 tx11 = power2(tx52); 		/* rule 26 */
	 tx53 = power2(tx62); 		/* rule 27 */
	 tx54 = tx62; 		/* rule 28 */
	 tx12 = tx53*tx54; 		/* rule 29 */
	 tx13 = dA*tx11; 		/* rule 30 */
	 tx14 = -(dC*tx12); 		/* rule 31 */
	 Efvdw = tx13 + tx14; 		/* rule 32 */
	FNONBOND_EFVDW_ENERGY_ACCUMULATE(Efvdw);
	 tx15 = reciprocal(NonbondDistance); 		/* rule 34 */
	 Efeel = dQ1Q2*tx15; 		/* rule 35 */
	FNONBOND_EFEEL_ENERGY_ACCUMULATE(Efeel);
	 Energy = Efeel + Efvdw; 		/* rule 37 */
	FNONBOND_ENERGY_ACCUMULATE(Energy);
	#ifdef FNONBOND_CALC_FORCE //[
	if ( calcForce ) {
	 tx55 = power2(tx61); 		/* rule 41 */
	 tx16 = tx52*tx55; 		/* rule 42 */
	 tx17 = power2(tx53); 		/* rule 43 */
	 tx56 = tx10; 		/* rule 44 */
	 tx57 = tx15; 		/* rule 45 */
	 tx58 = reciprocal(tx56); 		/* rule 46 */
	 tx18 = tx57*tx58; 		/* rule 47 */
	 tx19 = -12.*dA*tx16*tx4; 		/* rule 48 */
	 tx20 = 6.*dC*tx17*tx4; 		/* rule 49 */
	 tx21 = -(dQ1Q2*tx18*tx4); 		/* rule 50 */
	 gx1 = tx19 + tx20 + tx21; 		/* rule 51 */
	 fx1 = -gx1; 		/* rule 52 */
	FNONBOND_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx22 = -12.*dA*tx16*tx5; 		/* rule 54 */
	 tx23 = 6.*dC*tx17*tx5; 		/* rule 55 */
	 tx24 = -(dQ1Q2*tx18*tx5); 		/* rule 56 */
	 gy1 = tx22 + tx23 + tx24; 		/* rule 57 */
	 fy1 = -gy1; 		/* rule 58 */
	FNONBOND_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx25 = -12.*dA*tx16*tx6; 		/* rule 60 */
	 tx26 = 6.*dC*tx17*tx6; 		/* rule 61 */
	 tx27 = -(dQ1Q2*tx18*tx6); 		/* rule 62 */
	 gz1 = tx25 + tx26 + tx27; 		/* rule 63 */
	 fz1 = -gz1; 		/* rule 64 */
	FNONBOND_FORCE_ACCUMULATE(I1, 2, fz1 );
	#ifdef FNONBOND_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx28 = power2(tx55); 		/* rule 68 */
	 tx29 = tx52*tx53; 		/* rule 69 */
	 tx59 = power2(tx56); 		/* rule 70 */
	 tx60 = reciprocal(tx59); 		/* rule 71 */
	 tx30 = tx57*tx60; 		/* rule 72 */
	 tx31 = -12.*dA*tx16; 		/* rule 73 */
	 tx32 = 6.*dC*tx17; 		/* rule 74 */
	 tx33 = -(dQ1Q2*tx18); 		/* rule 75 */
	 tx34 = 168*dA*tx28*tx7; 		/* rule 76 */
	 tx35 = -48*dC*tx29*tx7; 		/* rule 77 */
	 tx36 = 3*dQ1Q2*tx30*tx7; 		/* rule 78 */
	 dhx1x1 = tx31 + tx32 + tx33 + tx34 + tx35 + tx36; 		/* rule 79 */
	FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx37 = 168*dA*tx28*tx8; 		/* rule 81 */
	 tx38 = -48*dC*tx29*tx8; 		/* rule 82 */
	 tx39 = 3*dQ1Q2*tx30*tx8; 		/* rule 83 */
	 dhy1y1 = tx31 + tx32 + tx33 + tx37 + tx38 + tx39; 		/* rule 84 */
	FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx40 = 168*dA*tx28*tx9; 		/* rule 86 */
	 tx41 = -48*dC*tx29*tx9; 		/* rule 87 */
	 tx42 = 3*dQ1Q2*tx30*tx9; 		/* rule 88 */
	 dhz1z1 = tx31 + tx32 + tx33 + tx40 + tx41 + tx42; 		/* rule 89 */
	FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	#ifdef FNONBOND_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tx43 = 168*dA*tx28*tx4*tx5; 		/* rule 93 */
	 tx44 = -48*dC*tx29*tx4*tx5; 		/* rule 94 */
	 tx45 = 3*dQ1Q2*tx30*tx4*tx5; 		/* rule 95 */
	 ohx1y1 = tx43 + tx44 + tx45; 		/* rule 96 */
	FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tx46 = 168*dA*tx28*tx4*tx6; 		/* rule 98 */
	 tx47 = -48*dC*tx29*tx4*tx6; 		/* rule 99 */
	 tx48 = 3*dQ1Q2*tx30*tx4*tx6; 		/* rule 100 */
	 ohx1z1 = tx46 + tx47 + tx48; 		/* rule 101 */
	FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx49 = 168*dA*tx28*tx5*tx6; 		/* rule 103 */
	 tx50 = -48*dC*tx29*tx5*tx6; 		/* rule 104 */
	 tx51 = 3*dQ1Q2*tx30*tx5*tx6; 		/* rule 105 */
	 ohy1z1 = tx49 + tx50 + tx51; 		/* rule 106 */
	FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	} /*if calcOffDiagonalHessian */ 
	#endif /* FNONBOND_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* FNONBOND_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* FNONBOND_CALC_FORCE ]*/
