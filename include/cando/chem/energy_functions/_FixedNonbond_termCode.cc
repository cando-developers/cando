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
	 tzz72 = dA*tx16; 		/* rule 48 */
	 tzz78 = -12.*tzz72; 		/* rule 49 */
	 tx19 = tx4*tzz78; 		/* rule 50 */
	 tzz71 = dC*tx17; 		/* rule 51 */
	 tzz76 = 6.*tzz71; 		/* rule 52 */
	 tx20 = tx4*tzz76; 		/* rule 53 */
	 tzz70 = dQ1Q2*tx18; 		/* rule 54 */
	 tzz77 = -tzz70; 		/* rule 55 */
	 tx21 = tx4*tzz77; 		/* rule 56 */
	 gx1 = tx19 + tx20 + tx21; 		/* rule 57 */
	 fx1 = -gx1; 		/* rule 58 */
	FNONBOND_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx22 = tx5*tzz78; 		/* rule 60 */
	 tx23 = tx5*tzz76; 		/* rule 61 */
	 tx24 = tx5*tzz77; 		/* rule 62 */
	 gy1 = tx22 + tx23 + tx24; 		/* rule 63 */
	 fy1 = -gy1; 		/* rule 64 */
	FNONBOND_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx25 = tx6*tzz78; 		/* rule 66 */
	 tx26 = tx6*tzz76; 		/* rule 67 */
	 tx27 = tx6*tzz77; 		/* rule 68 */
	 gz1 = tx25 + tx26 + tx27; 		/* rule 69 */
	 fz1 = -gz1; 		/* rule 70 */
	FNONBOND_FORCE_ACCUMULATE(I1, 2, fz1 );
	#ifdef FNONBOND_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx28 = power2(tx55); 		/* rule 74 */
	 tx29 = tx52*tx53; 		/* rule 75 */
	 tx59 = power2(tx56); 		/* rule 76 */
	 tx60 = reciprocal(tx59); 		/* rule 77 */
	 tx30 = tx57*tx60; 		/* rule 78 */
	 tx31 = tzz78; 		/* rule 79 */
	 tx32 = tzz76; 		/* rule 80 */
	 tx33 = tzz77; 		/* rule 81 */
	 tzz66 = dA*tx28; 		/* rule 82 */
	 tzz67 = 168*tzz66; 		/* rule 83 */
	 tx34 = tx7*tzz67; 		/* rule 84 */
	 tzz65 = dC*tx29; 		/* rule 85 */
	 tzz69 = -48*tzz65; 		/* rule 86 */
	 tx35 = tx7*tzz69; 		/* rule 87 */
	 tzz64 = dQ1Q2*tx30; 		/* rule 88 */
	 tzz68 = 3*tzz64; 		/* rule 89 */
	 tx36 = tx7*tzz68; 		/* rule 90 */
	 tzz79 = tx32 + tx33; 		/* rule 91 */
	 tzz80 = tx31 + tzz79; 		/* rule 92 */
	 dhx1x1 = tx34 + tx35 + tx36 + tzz80; 		/* rule 93 */
	FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx37 = tx8*tzz67; 		/* rule 95 */
	 tx38 = tx8*tzz69; 		/* rule 96 */
	 tx39 = tx8*tzz68; 		/* rule 97 */
	 dhy1y1 = tx37 + tx38 + tx39 + tzz80; 		/* rule 98 */
	FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx40 = tx9*tzz67; 		/* rule 100 */
	 tx41 = tx9*tzz69; 		/* rule 101 */
	 tx42 = tx9*tzz68; 		/* rule 102 */
	 dhz1z1 = tx40 + tx41 + tx42 + tzz80; 		/* rule 103 */
	FNONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	#ifdef FNONBOND_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tzz75 = tx4*tx5; 		/* rule 107 */
	 tx43 = tzz67*tzz75; 		/* rule 108 */
	 tx44 = tzz69*tzz75; 		/* rule 109 */
	 tx45 = tzz68*tzz75; 		/* rule 110 */
	 ohx1y1 = tx43 + tx44 + tx45; 		/* rule 111 */
	FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tzz74 = tx4*tx6; 		/* rule 113 */
	 tx46 = tzz67*tzz74; 		/* rule 114 */
	 tx47 = tzz69*tzz74; 		/* rule 115 */
	 tx48 = tzz68*tzz74; 		/* rule 116 */
	 ohx1z1 = tx46 + tx47 + tx48; 		/* rule 117 */
	FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tzz73 = tx5*tx6; 		/* rule 119 */
	 tx49 = tzz67*tzz73; 		/* rule 120 */
	 tx50 = tzz69*tzz73; 		/* rule 121 */
	 tx51 = tzz68*tzz73; 		/* rule 122 */
	 ohy1z1 = tx49 + tx50 + tx51; 		/* rule 123 */
	FNONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	} /*if calcOffDiagonalHessian */ 
	#endif /* FNONBOND_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* FNONBOND_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* FNONBOND_CALC_FORCE ]*/
