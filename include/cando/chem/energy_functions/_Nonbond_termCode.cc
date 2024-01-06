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
	 tzz98 = dA*tx15; 		/* rule 54 */
	 tzz109 = 12.*tzz98; 		/* rule 55 */
	 tx18 = tx4*tzz109; 		/* rule 56 */
	 tzz97 = dC*tx16; 		/* rule 57 */
	 tzz111 = -6.*tzz97; 		/* rule 58 */
	 tx19 = tx4*tzz111; 		/* rule 59 */
	 tzz108 = dQ1Q2*tx17; 		/* rule 60 */
	 tx20 = tx4*tzz108; 		/* rule 61 */
	 gx1 = tx18 + tx19 + tx20; 		/* rule 62 */
	 fx1 = -gx1; 		/* rule 63 */
	NONBOND_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx21 = tx5*tzz109; 		/* rule 65 */
	 tx22 = tx5*tzz111; 		/* rule 66 */
	 tx23 = tx5*tzz108; 		/* rule 67 */
	 gy1 = tx21 + tx22 + tx23; 		/* rule 68 */
	 fy1 = -gy1; 		/* rule 69 */
	NONBOND_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx24 = tx6*tzz109; 		/* rule 71 */
	 tx25 = tx6*tzz111; 		/* rule 72 */
	 tx26 = tx6*tzz108; 		/* rule 73 */
	 gz1 = tx24 + tx25 + tx26; 		/* rule 74 */
	 fz1 = -gz1; 		/* rule 75 */
	NONBOND_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tx27 = -tx20; 		/* rule 77 */
	 tzz112 = -12.*tzz98; 		/* rule 78 */
	 tx28 = tx4*tzz112; 		/* rule 79 */
	 tzz110 = 6.*tzz97; 		/* rule 80 */
	 tx29 = tx4*tzz110; 		/* rule 81 */
	 gx2 = tx27 + tx28 + tx29; 		/* rule 82 */
	 fx2 = -gx2; 		/* rule 83 */
	NONBOND_FORCE_ACCUMULATE(I2, 0, fx2 );
	 tx30 = -tx23; 		/* rule 85 */
	 tx31 = tx5*tzz112; 		/* rule 86 */
	 tx32 = tx5*tzz110; 		/* rule 87 */
	 gy2 = tx30 + tx31 + tx32; 		/* rule 88 */
	 fy2 = -gy2; 		/* rule 89 */
	NONBOND_FORCE_ACCUMULATE(I2, 1, fy2 );
	 tx33 = -tx26; 		/* rule 91 */
	 tx34 = tx6*tzz112; 		/* rule 92 */
	 tx35 = tx6*tzz110; 		/* rule 93 */
	 gz2 = tx33 + tx34 + tx35; 		/* rule 94 */
	 fz2 = -gz2; 		/* rule 95 */
	NONBOND_FORCE_ACCUMULATE(I2, 2, fz2 );
	#ifdef NONBOND_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx36 = power2(tx85); 		/* rule 99 */
	 tx37 = tx81*tx82; 		/* rule 100 */
	 tx89 = power2(tx86); 		/* rule 101 */
	 tx90 = reciprocal(tx89); 		/* rule 102 */
	 tx38 = tx87*tx90; 		/* rule 103 */
	 tx39 = tzz112; 		/* rule 104 */
	 tx40 = tzz110; 		/* rule 105 */
	 tx41 = -tzz108; 		/* rule 106 */
	 tzz96 = dA*tx36; 		/* rule 107 */
	 tzz102 = 168*tzz96; 		/* rule 108 */
	 tx42 = tx7*tzz102; 		/* rule 109 */
	 tzz95 = dC*tx37; 		/* rule 110 */
	 tzz106 = -48*tzz95; 		/* rule 111 */
	 tx43 = tx7*tzz106; 		/* rule 112 */
	 tzz94 = dQ1Q2*tx38; 		/* rule 113 */
	 tzz104 = 3*tzz94; 		/* rule 114 */
	 tx44 = tx7*tzz104; 		/* rule 115 */
	 tzz115 = tx40 + tx41; 		/* rule 116 */
	 tzz116 = tx39 + tzz115; 		/* rule 117 */
	 dhx1x1 = tx42 + tx43 + tx44 + tzz116; 		/* rule 118 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx45 = tx8*tzz102; 		/* rule 120 */
	 tx46 = tx8*tzz106; 		/* rule 121 */
	 tx47 = tx8*tzz104; 		/* rule 122 */
	 dhy1y1 = tx45 + tx46 + tx47 + tzz116; 		/* rule 123 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx48 = tx9*tzz102; 		/* rule 125 */
	 tx49 = tx9*tzz106; 		/* rule 126 */
	 tx50 = tx9*tzz104; 		/* rule 127 */
	 dhz1z1 = tx48 + tx49 + tx50 + tzz116; 		/* rule 128 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 dhx2x2 = dhx1x1; 		/* rule 130 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 dhy2y2 = dhy1y1; 		/* rule 132 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 dhz2z2 = dhz1z1; 		/* rule 134 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	#ifdef NONBOND_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tzz101 = tx4*tx5; 		/* rule 138 */
	 tx51 = tzz101*tzz102; 		/* rule 139 */
	 tx52 = tzz101*tzz106; 		/* rule 140 */
	 tx53 = tzz101*tzz104; 		/* rule 141 */
	 ohx1y1 = tx51 + tx52 + tx53; 		/* rule 142 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tzz100 = tx4*tx6; 		/* rule 144 */
	 tx54 = tzz100*tzz102; 		/* rule 145 */
	 tx55 = tzz100*tzz106; 		/* rule 146 */
	 tx56 = tzz100*tzz104; 		/* rule 147 */
	 ohx1z1 = tx54 + tx55 + tx56; 		/* rule 148 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx57 = tzz109; 		/* rule 150 */
	 tx58 = tzz111; 		/* rule 151 */
	 tx59 = tzz108; 		/* rule 152 */
	 tzz107 = -168*tzz96; 		/* rule 153 */
	 tx60 = tx7*tzz107; 		/* rule 154 */
	 tzz103 = 48*tzz95; 		/* rule 155 */
	 tx61 = tx7*tzz103; 		/* rule 156 */
	 tzz105 = -3*tzz94; 		/* rule 157 */
	 tx62 = tx7*tzz105; 		/* rule 158 */
	 tzz113 = tx58 + tx59; 		/* rule 159 */
	 tzz114 = tx57 + tzz113; 		/* rule 160 */
	 ohx1x2 = tx60 + tx61 + tx62 + tzz114; 		/* rule 161 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx63 = tzz101*tzz107; 		/* rule 163 */
	 tx64 = tzz101*tzz103; 		/* rule 164 */
	 tx65 = tzz101*tzz105; 		/* rule 165 */
	 ohx1y2 = tx63 + tx64 + tx65; 		/* rule 166 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx66 = tzz100*tzz107; 		/* rule 168 */
	 tx67 = tzz100*tzz103; 		/* rule 169 */
	 tx68 = tzz100*tzz105; 		/* rule 170 */
	 ohx1z2 = tx66 + tx67 + tx68; 		/* rule 171 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tzz99 = tx5*tx6; 		/* rule 173 */
	 tx69 = tzz102*tzz99; 		/* rule 174 */
	 tx70 = tzz106*tzz99; 		/* rule 175 */
	 tx71 = tzz104*tzz99; 		/* rule 176 */
	 ohy1z1 = tx69 + tx70 + tx71; 		/* rule 177 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 ohy1x2 = ohx1y2; 		/* rule 179 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx72 = tx8*tzz107; 		/* rule 181 */
	 tx73 = tx8*tzz103; 		/* rule 182 */
	 tx74 = tx8*tzz105; 		/* rule 183 */
	 ohy1y2 = tx72 + tx73 + tx74 + tzz114; 		/* rule 184 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx75 = tzz107*tzz99; 		/* rule 186 */
	 tx76 = tzz103*tzz99; 		/* rule 187 */
	 tx77 = tzz105*tzz99; 		/* rule 188 */
	 ohy1z2 = tx75 + tx76 + tx77; 		/* rule 189 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 ohz1x2 = ohx1z2; 		/* rule 191 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 ohz1y2 = ohy1z2; 		/* rule 193 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx78 = tx9*tzz107; 		/* rule 195 */
	 tx79 = tx9*tzz103; 		/* rule 196 */
	 tx80 = tx9*tzz105; 		/* rule 197 */
	 ohz1z2 = tx78 + tx79 + tx80 + tzz114; 		/* rule 198 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 ohx2y2 = ohx1y1; 		/* rule 200 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 ohx2z2 = ohx1z1; 		/* rule 202 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 ohy2z2 = ohy1z1; 		/* rule 204 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	} /*if calcOffDiagonalHessian */ 
	#endif /* NONBOND_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* NONBOND_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* NONBOND_CALC_FORCE ]*/
	NONBOND_DEBUG_INTERACTIONS(I1,I2);
	SKIP_term: (void)0;
