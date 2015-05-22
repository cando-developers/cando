// C-code
	// FNONBOND_SET_PARAMETER(xf);
	// FNONBOND_SET_PARAMETER(yf);
	// FNONBOND_SET_PARAMETER(zf);
	// FNONBOND_SET_PARAMETER(I1);
	// FNONBOND_SET_PARAMETER(dQ1Q2);
	// FNONBOND_SET_PARAMETER(dA);
	// FNONBOND_SET_PARAMETER(dC);
	FNONBOND_SET_POSITION(x1,I1,0);
	FNONBOND_SET_POSITION(y1,I1,1);
	FNONBOND_SET_POSITION(z1,I1,2);
	 tx52 = -xf; 		/* rule 11 */
	 tx53 = -yf; 		/* rule 12 */
	 tx54 = -zf; 		/* rule 13 */
	 tx55 = tx52 + x1; 		/* rule 14 */
	 tx56 = tx53 + y1; 		/* rule 15 */
	 tx57 = tx54 + z1; 		/* rule 16 */
	 tx58 = power2(tx55); 		/* rule 17 */
	 tx59 = power2(tx56); 		/* rule 18 */
	 tx60 = power2(tx57); 		/* rule 19 */
	 tx61 = tx58 + tx59 + tx60; 		/* rule 20 */
	 NonbondDistance = mysqrt(tx61); 		/* rule 21 */
	 tx114 = reciprocal(tx61); 		/* rule 22 */
	 tx112 = power2(tx114); 		/* rule 23 */
	 tx113 = tx114; 		/* rule 24 */
	 tx103 = tx112*tx113; 		/* rule 25 */
	 tx62 = power2(tx103); 		/* rule 26 */
	 tx104 = power2(tx113); 		/* rule 27 */
	 tx105 = tx113; 		/* rule 28 */
	 tx63 = tx104*tx105; 		/* rule 29 */
	 tx64 = dA*tx62; 		/* rule 30 */
	 tx65 = -(dC*tx63); 		/* rule 31 */
	 Efvdw = tx64 + tx65; 		/* rule 32 */
	FNONBOND_EFVDW_ENERGY_ACCUMULATE(Efvdw);
	 tx66 = reciprocal(NonbondDistance); 		/* rule 34 */
	 Efeel = dQ1Q2*tx66; 		/* rule 35 */
	FNONBOND_EFEEL_ENERGY_ACCUMULATE(Efeel);
	 Energy = Efeel + Efvdw; 		/* rule 37 */
	FNONBOND_ENERGY_ACCUMULATE(Energy);
	#ifdef FNONBOND_CALC_FORCE //[
	if ( calcForce ) {
	 tx106 = power2(tx112); 		/* rule 41 */
	 tx67 = tx103*tx106; 		/* rule 42 */
	 tx68 = power2(tx104); 		/* rule 43 */
	 tx107 = tx61; 		/* rule 44 */
	 tx108 = tx66; 		/* rule 45 */
	 tx109 = reciprocal(tx107); 		/* rule 46 */
	 tx69 = tx108*tx109; 		/* rule 47 */
	 tx70 = -12.*dA*tx55*tx67; 		/* rule 48 */
	 tx71 = 6.*dC*tx55*tx68; 		/* rule 49 */
	 tx72 = -(dQ1Q2*tx55*tx69); 		/* rule 50 */
	 gx1 = tx70 + tx71 + tx72; 		/* rule 51 */
	 fx1 = -gx1; 		/* rule 52 */
	FNONBOND_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx73 = -12.*dA*tx56*tx67; 		/* rule 54 */
	 tx74 = 6.*dC*tx56*tx68; 		/* rule 55 */
	 tx75 = -(dQ1Q2*tx56*tx69); 		/* rule 56 */
	 gy1 = tx73 + tx74 + tx75; 		/* rule 57 */
	 fy1 = -gy1; 		/* rule 58 */
	FNONBOND_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx76 = -12.*dA*tx57*tx67; 		/* rule 60 */
	 tx77 = 6.*dC*tx57*tx68; 		/* rule 61 */
	 tx78 = -(dQ1Q2*tx57*tx69); 		/* rule 62 */
	 gz1 = tx76 + tx77 + tx78; 		/* rule 63 */
	 fz1 = -gz1; 		/* rule 64 */
	FNONBOND_FORCE_ACCUMULATE(I1, 2, fz1 );
	#ifdef FNONBOND_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx79 = power2(tx106); 		/* rule 68 */
	 tx80 = tx103*tx104; 		/* rule 69 */
	 tx110 = power2(tx107); 		/* rule 70 */
	 tx111 = reciprocal(tx110); 		/* rule 71 */
	 tx81 = tx108*tx111; 		/* rule 72 */
	 tx82 = -12.*dA*tx67; 		/* rule 73 */
	 tx83 = 6.*dC*tx68; 		/* rule 74 */
	 tx84 = -(dQ1Q2*tx69); 		/* rule 75 */
	 tx85 = 168*dA*tx58*tx79; 		/* rule 76 */
	 tx86 = -48*dC*tx58*tx80; 		/* rule 77 */
	 tx87 = 3*dQ1Q2*tx58*tx81; 		/* rule 78 */
	 dhx1x1 = tx82 + tx83 + tx84 + tx85 + tx86 + tx87; 		/* rule 79 */
	FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx88 = 168*dA*tx59*tx79; 		/* rule 81 */
	 tx89 = -48*dC*tx59*tx80; 		/* rule 82 */
	 tx90 = 3*dQ1Q2*tx59*tx81; 		/* rule 83 */
	 dhy1y1 = tx82 + tx83 + tx84 + tx88 + tx89 + tx90; 		/* rule 84 */
	FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx91 = 168*dA*tx60*tx79; 		/* rule 86 */
	 tx92 = -48*dC*tx60*tx80; 		/* rule 87 */
	 tx93 = 3*dQ1Q2*tx60*tx81; 		/* rule 88 */
	 dhz1z1 = tx82 + tx83 + tx84 + tx91 + tx92 + tx93; 		/* rule 89 */
	FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	#ifdef FNONBOND_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tx94 = 168*dA*tx55*tx56*tx79; 		/* rule 93 */
	 tx95 = -48*dC*tx55*tx56*tx80; 		/* rule 94 */
	 tx96 = 3*dQ1Q2*tx55*tx56*tx81; 		/* rule 95 */
	 ohx1y1 = tx94 + tx95 + tx96; 		/* rule 96 */
	FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tx97 = 168*dA*tx55*tx57*tx79; 		/* rule 98 */
	 tx98 = -48*dC*tx55*tx57*tx80; 		/* rule 99 */
	 tx99 = 3*dQ1Q2*tx55*tx57*tx81; 		/* rule 100 */
	 ohx1z1 = tx97 + tx98 + tx99; 		/* rule 101 */
	FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx100 = 168*dA*tx56*tx57*tx79; 		/* rule 103 */
	 tx101 = -48*dC*tx56*tx57*tx80; 		/* rule 104 */
	 tx102 = 3*dQ1Q2*tx56*tx57*tx81; 		/* rule 105 */
	 ohy1z1 = tx100 + tx101 + tx102; 		/* rule 106 */
	FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	} /*if calcOffDiagonalHessian */ 
	#endif /* FNONBOND_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* FNONBOND_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* FNONBOND_CALC_FORCE ]*/
