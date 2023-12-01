// C-code
	NONBOND_SET_PARAMETER(dQ1Q2);
	NONBOND_SET_PARAMETER(dA);
	NONBOND_SET_PARAMETER(dC);
	NONBOND_SET_PARAMETER(I1);
	NONBOND_SET_PARAMETER(I2);
	NONBOND_APPLY_ATOM_MASK(I1,I2);
	NONBOND_SET_POSITION(x1,I1,0);
	NONBOND_SET_POSITION(y1,I1,1);
	NONBOND_SET_POSITION(z1,I1,2);
	NONBOND_SET_POSITION(x2,I2,0);
	NONBOND_SET_POSITION(y2,I2,1);
	NONBOND_SET_POSITION(z2,I2,2);
	 tx1 = -x1; 		/* rule 13 */
	 DeltaX = tx1 + x2; 		/* rule 14 */
	 tx2 = -y1; 		/* rule 15 */
	 DeltaY = tx2 + y2; 		/* rule 16 */
	 tx3 = -z1; 		/* rule 17 */
	 DeltaZ = tx3 + z2; 		/* rule 18 */
	 tx4 = PBX(DeltaX); 		/* rule 19 */
	 tx5 = PBY(DeltaY); 		/* rule 20 */
	 tx6 = PBZ(DeltaZ); 		/* rule 21 */
	 tx7 = power2(tx4); 		/* rule 22 */
	 tx8 = power2(tx5); 		/* rule 23 */
	 tx9 = power2(tx6); 		/* rule 24 */
	 DistanceSquared = tx7 + tx8 + tx9; 		/* rule 25 */
	BAIL_OUT_IF_CUTOFF(DistanceSquared)
	 tx93 = reciprocal(DistanceSquared); 		/* rule 27 */
	 tx91 = power2(tx93); 		/* rule 28 */
	 tx92 = tx93; 		/* rule 29 */
	 tx81 = tx91*tx92; 		/* rule 30 */
	 tx10 = power2(tx81); 		/* rule 31 */
	 tx82 = power2(tx92); 		/* rule 32 */
	 tx83 = tx92; 		/* rule 33 */
	 tx11 = tx82*tx83; 		/* rule 34 */
	 tx12 = dA*tx10; 		/* rule 35 */
	 tx13 = -(dC*tx11); 		/* rule 36 */
	 Evdw = tx12 + tx13; 		/* rule 37 */
	NONBOND_EVDW_ENERGY_ACCUMULATE(Evdw);
	 tx84 = mysqrt(DistanceSquared); 		/* rule 39 */
	 tx14 = reciprocal(tx84); 		/* rule 40 */
	 Eeel = dQ1Q2*tx14; 		/* rule 41 */
	NONBOND_EEEL_ENERGY_ACCUMULATE(Eeel);
	 Energy = Eeel + Evdw; 		/* rule 43 */
	NONBOND_ENERGY_ACCUMULATE(Energy);
	#ifdef NONBOND_CALC_FORCE //[
	if ( calcForce ) {
	 tx85 = power2(tx91); 		/* rule 47 */
	 tx15 = tx81*tx85; 		/* rule 48 */
	 tx16 = power2(tx82); 		/* rule 49 */
	 tx86 = DistanceSquared; 		/* rule 50 */
	 tx87 = tx14; 		/* rule 51 */
	 tx88 = reciprocal(tx86); 		/* rule 52 */
	 tx17 = tx87*tx88; 		/* rule 53 */
	 tx18 = 12.*dA*tx15*tx4; 		/* rule 54 */
	 tx19 = -6.*dC*tx16*tx4; 		/* rule 55 */
	 tx20 = dQ1Q2*tx17*tx4; 		/* rule 56 */
	 gx1 = tx18 + tx19 + tx20; 		/* rule 57 */
	 fx1 = -gx1; 		/* rule 58 */
	NONBOND_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx21 = 12.*dA*tx15*tx5; 		/* rule 60 */
	 tx22 = -6.*dC*tx16*tx5; 		/* rule 61 */
	 tx23 = dQ1Q2*tx17*tx5; 		/* rule 62 */
	 gy1 = tx21 + tx22 + tx23; 		/* rule 63 */
	 fy1 = -gy1; 		/* rule 64 */
	NONBOND_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx24 = 12.*dA*tx15*tx6; 		/* rule 66 */
	 tx25 = -6.*dC*tx16*tx6; 		/* rule 67 */
	 tx26 = dQ1Q2*tx17*tx6; 		/* rule 68 */
	 gz1 = tx24 + tx25 + tx26; 		/* rule 69 */
	 fz1 = -gz1; 		/* rule 70 */
	NONBOND_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tx27 = -tx20; 		/* rule 72 */
	 tx28 = -12.*dA*tx15*tx4; 		/* rule 73 */
	 tx29 = 6.*dC*tx16*tx4; 		/* rule 74 */
	 gx2 = tx27 + tx28 + tx29; 		/* rule 75 */
	 fx2 = -gx2; 		/* rule 76 */
	NONBOND_FORCE_ACCUMULATE(I2, 0, fx2 );
	 tx30 = -tx23; 		/* rule 78 */
	 tx31 = -12.*dA*tx15*tx5; 		/* rule 79 */
	 tx32 = 6.*dC*tx16*tx5; 		/* rule 80 */
	 gy2 = tx30 + tx31 + tx32; 		/* rule 81 */
	 fy2 = -gy2; 		/* rule 82 */
	NONBOND_FORCE_ACCUMULATE(I2, 1, fy2 );
	 tx33 = -tx26; 		/* rule 84 */
	 tx34 = -12.*dA*tx15*tx6; 		/* rule 85 */
	 tx35 = 6.*dC*tx16*tx6; 		/* rule 86 */
	 gz2 = tx33 + tx34 + tx35; 		/* rule 87 */
	 fz2 = -gz2; 		/* rule 88 */
	NONBOND_FORCE_ACCUMULATE(I2, 2, fz2 );
	#ifdef NONBOND_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx36 = power2(tx85); 		/* rule 92 */
	 tx37 = tx81*tx82; 		/* rule 93 */
	 tx89 = power2(tx86); 		/* rule 94 */
	 tx90 = reciprocal(tx89); 		/* rule 95 */
	 tx38 = tx87*tx90; 		/* rule 96 */
	 tx39 = -12.*dA*tx15; 		/* rule 97 */
	 tx40 = 6.*dC*tx16; 		/* rule 98 */
	 tx41 = -(dQ1Q2*tx17); 		/* rule 99 */
	 tx42 = 168*dA*tx36*tx7; 		/* rule 100 */
	 tx43 = -48*dC*tx37*tx7; 		/* rule 101 */
	 tx44 = 3*dQ1Q2*tx38*tx7; 		/* rule 102 */
	 dhx1x1 = tx39 + tx40 + tx41 + tx42 + tx43 + tx44; 		/* rule 103 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx45 = 168*dA*tx36*tx8; 		/* rule 105 */
	 tx46 = -48*dC*tx37*tx8; 		/* rule 106 */
	 tx47 = 3*dQ1Q2*tx38*tx8; 		/* rule 107 */
	 dhy1y1 = tx39 + tx40 + tx41 + tx45 + tx46 + tx47; 		/* rule 108 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx48 = 168*dA*tx36*tx9; 		/* rule 110 */
	 tx49 = -48*dC*tx37*tx9; 		/* rule 111 */
	 tx50 = 3*dQ1Q2*tx38*tx9; 		/* rule 112 */
	 dhz1z1 = tx39 + tx40 + tx41 + tx48 + tx49 + tx50; 		/* rule 113 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 dhx2x2 = dhx1x1; 		/* rule 115 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 dhy2y2 = dhy1y1; 		/* rule 117 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 dhz2z2 = dhz1z1; 		/* rule 119 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	#ifdef NONBOND_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tx51 = 168*dA*tx36*tx4*tx5; 		/* rule 123 */
	 tx52 = -48*dC*tx37*tx4*tx5; 		/* rule 124 */
	 tx53 = 3*dQ1Q2*tx38*tx4*tx5; 		/* rule 125 */
	 ohx1y1 = tx51 + tx52 + tx53; 		/* rule 126 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tx54 = 168*dA*tx36*tx4*tx6; 		/* rule 128 */
	 tx55 = -48*dC*tx37*tx4*tx6; 		/* rule 129 */
	 tx56 = 3*dQ1Q2*tx38*tx4*tx6; 		/* rule 130 */
	 ohx1z1 = tx54 + tx55 + tx56; 		/* rule 131 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx57 = 12.*dA*tx15; 		/* rule 133 */
	 tx58 = -6.*dC*tx16; 		/* rule 134 */
	 tx59 = dQ1Q2*tx17; 		/* rule 135 */
	 tx60 = -168*dA*tx36*tx7; 		/* rule 136 */
	 tx61 = 48*dC*tx37*tx7; 		/* rule 137 */
	 tx62 = -3*dQ1Q2*tx38*tx7; 		/* rule 138 */
	 ohx1x2 = tx57 + tx58 + tx59 + tx60 + tx61 + tx62; 		/* rule 139 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx63 = -168*dA*tx36*tx4*tx5; 		/* rule 141 */
	 tx64 = 48*dC*tx37*tx4*tx5; 		/* rule 142 */
	 tx65 = -3*dQ1Q2*tx38*tx4*tx5; 		/* rule 143 */
	 ohx1y2 = tx63 + tx64 + tx65; 		/* rule 144 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx66 = -168*dA*tx36*tx4*tx6; 		/* rule 146 */
	 tx67 = 48*dC*tx37*tx4*tx6; 		/* rule 147 */
	 tx68 = -3*dQ1Q2*tx38*tx4*tx6; 		/* rule 148 */
	 ohx1z2 = tx66 + tx67 + tx68; 		/* rule 149 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tx69 = 168*dA*tx36*tx5*tx6; 		/* rule 151 */
	 tx70 = -48*dC*tx37*tx5*tx6; 		/* rule 152 */
	 tx71 = 3*dQ1Q2*tx38*tx5*tx6; 		/* rule 153 */
	 ohy1z1 = tx69 + tx70 + tx71; 		/* rule 154 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 ohy1x2 = ohx1y2; 		/* rule 156 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx72 = -168*dA*tx36*tx8; 		/* rule 158 */
	 tx73 = 48*dC*tx37*tx8; 		/* rule 159 */
	 tx74 = -3*dQ1Q2*tx38*tx8; 		/* rule 160 */
	 ohy1y2 = tx57 + tx58 + tx59 + tx72 + tx73 + tx74; 		/* rule 161 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx75 = -168*dA*tx36*tx5*tx6; 		/* rule 163 */
	 tx76 = 48*dC*tx37*tx5*tx6; 		/* rule 164 */
	 tx77 = -3*dQ1Q2*tx38*tx5*tx6; 		/* rule 165 */
	 ohy1z2 = tx75 + tx76 + tx77; 		/* rule 166 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 ohz1x2 = ohx1z2; 		/* rule 168 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 ohz1y2 = ohy1z2; 		/* rule 170 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx78 = -168*dA*tx36*tx9; 		/* rule 172 */
	 tx79 = 48*dC*tx37*tx9; 		/* rule 173 */
	 tx80 = -3*dQ1Q2*tx38*tx9; 		/* rule 174 */
	 ohz1z2 = tx57 + tx58 + tx59 + tx78 + tx79 + tx80; 		/* rule 175 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 ohx2y2 = ohx1y1; 		/* rule 177 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 ohx2z2 = ohx1z1; 		/* rule 179 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 ohy2z2 = ohy1z1; 		/* rule 181 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	} /*if calcOffDiagonalHessian */ 
	#endif /* NONBOND_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* NONBOND_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* NONBOND_CALC_FORCE ]*/
	SKIP_term:  (void)0;
