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
	 tx131 = reciprocal(DistanceSquared); 		/* rule 27 */
	 tx129 = power2(tx131); 		/* rule 28 */
	 tx130 = tx131; 		/* rule 29 */
	 tx119 = tx129*tx130; 		/* rule 30 */
	 tx10 = power2(tx119); 		/* rule 31 */
	 tx120 = power2(tx130); 		/* rule 32 */
	 tx121 = tx130; 		/* rule 33 */
	 tx11 = tx120*tx121; 		/* rule 34 */
	 tx12 = dA*tx10; 		/* rule 35 */
	 tx13 = -(dC*tx11); 		/* rule 36 */
	 tx14 = tx12 + tx13; 		/* rule 37 */
	 Evdw = tx14*vdwScale; 		/* rule 38 */
	NONBOND_EVDW_ENERGY_ACCUMULATE(Evdw);
	 tx15 = reciprocal(DIELECTRIC); 		/* rule 40 */
	 tx122 = mysqrt(DistanceSquared); 		/* rule 41 */
	 tx16 = reciprocal(tx122); 		/* rule 42 */
	 tzz132 = eelScale*tx15; 		/* rule 43 */
	 tzz133 = dQ1Q2*tzz132; 		/* rule 44 */
	 Eeel = tx16*tzz133; 		/* rule 45 */
	NONBOND_EEEL_ENERGY_ACCUMULATE(Eeel);
	 Energy = Eeel + Evdw; 		/* rule 47 */
	NONBOND_ENERGY_ACCUMULATE(Energy);
	#ifdef NONBOND_CALC_FORCE //[
	if ( calcForce ) {
	 tx123 = power2(tx129); 		/* rule 51 */
	 tx17 = tx119*tx123; 		/* rule 52 */
	 tx18 = power2(tx120); 		/* rule 53 */
	 tzz138 = dA*tx17; 		/* rule 54 */
	 tzz149 = 12.*tzz138; 		/* rule 55 */
	 tx19 = tx4*tzz149; 		/* rule 56 */
	 tzz137 = dC*tx18; 		/* rule 57 */
	 tzz151 = -6.*tzz137; 		/* rule 58 */
	 tx20 = tx4*tzz151; 		/* rule 59 */
	 tx124 = DistanceSquared; 		/* rule 60 */
	 tx125 = tx16; 		/* rule 61 */
	 tx126 = reciprocal(tx124); 		/* rule 62 */
	 tx21 = tx125*tx126; 		/* rule 63 */
	 tx22 = tx19 + tx20; 		/* rule 64 */
	 tzz148 = tx21*tzz133; 		/* rule 65 */
	 tx23 = tx4*tzz148; 		/* rule 66 */
	 tx24 = tx22*vdwScale; 		/* rule 67 */
	 gx1 = tx23 + tx24; 		/* rule 68 */
	 fx1 = -gx1; 		/* rule 69 */
	NONBOND_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx25 = tx5*tzz149; 		/* rule 71 */
	 tx26 = tx5*tzz151; 		/* rule 72 */
	 tx27 = tx25 + tx26; 		/* rule 73 */
	 tx28 = tx5*tzz148; 		/* rule 74 */
	 tx29 = tx27*vdwScale; 		/* rule 75 */
	 gy1 = tx28 + tx29; 		/* rule 76 */
	 fy1 = -gy1; 		/* rule 77 */
	NONBOND_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx30 = tx6*tzz149; 		/* rule 79 */
	 tx31 = tx6*tzz151; 		/* rule 80 */
	 tx32 = tx30 + tx31; 		/* rule 81 */
	 tx33 = tx6*tzz148; 		/* rule 82 */
	 tx34 = tx32*vdwScale; 		/* rule 83 */
	 gz1 = tx33 + tx34; 		/* rule 84 */
	 fz1 = -gz1; 		/* rule 85 */
	NONBOND_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tzz152 = -12.*tzz138; 		/* rule 87 */
	 tx35 = tx4*tzz152; 		/* rule 88 */
	 tzz150 = 6.*tzz137; 		/* rule 89 */
	 tx36 = tx4*tzz150; 		/* rule 90 */
	 tx37 = tx35 + tx36; 		/* rule 91 */
	 tx38 = -tx23; 		/* rule 92 */
	 tx39 = tx37*vdwScale; 		/* rule 93 */
	 gx2 = tx38 + tx39; 		/* rule 94 */
	 fx2 = -gx2; 		/* rule 95 */
	NONBOND_FORCE_ACCUMULATE(I2, 0, fx2 );
	 tx40 = tx5*tzz152; 		/* rule 97 */
	 tx41 = tx5*tzz150; 		/* rule 98 */
	 tx42 = tx40 + tx41; 		/* rule 99 */
	 tx43 = -tx28; 		/* rule 100 */
	 tx44 = tx42*vdwScale; 		/* rule 101 */
	 gy2 = tx43 + tx44; 		/* rule 102 */
	 fy2 = -gy2; 		/* rule 103 */
	NONBOND_FORCE_ACCUMULATE(I2, 1, fy2 );
	 tx45 = tx6*tzz152; 		/* rule 105 */
	 tx46 = tx6*tzz150; 		/* rule 106 */
	 tx47 = tx45 + tx46; 		/* rule 107 */
	 tx48 = -tx33; 		/* rule 108 */
	 tx49 = tx47*vdwScale; 		/* rule 109 */
	 gz2 = tx48 + tx49; 		/* rule 110 */
	 fz2 = -gz2; 		/* rule 111 */
	NONBOND_FORCE_ACCUMULATE(I2, 2, fz2 );
	#ifdef NONBOND_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx50 = power2(tx123); 		/* rule 115 */
	 tx51 = tx119*tx120; 		/* rule 116 */
	 tx52 = tzz152; 		/* rule 117 */
	 tx53 = tzz150; 		/* rule 118 */
	 tzz136 = dA*tx50; 		/* rule 119 */
	 tzz142 = 168*tzz136; 		/* rule 120 */
	 tx54 = tx7*tzz142; 		/* rule 121 */
	 tzz135 = dC*tx51; 		/* rule 122 */
	 tzz146 = -48*tzz135; 		/* rule 123 */
	 tx55 = tx7*tzz146; 		/* rule 124 */
	 tx127 = power2(tx124); 		/* rule 125 */
	 tx128 = reciprocal(tx127); 		/* rule 126 */
	 tx56 = tx125*tx128; 		/* rule 127 */
	 tzz154 = tx52 + tx53; 		/* rule 128 */
	 tx57 = tx54 + tx55 + tzz154; 		/* rule 129 */
	 tx58 = -tzz148; 		/* rule 130 */
	 tzz134 = tx56*tzz133; 		/* rule 131 */
	 tzz144 = 3*tzz134; 		/* rule 132 */
	 tx59 = tx7*tzz144; 		/* rule 133 */
	 tx60 = tx57*vdwScale; 		/* rule 134 */
	 dhx1x1 = tx58 + tx59 + tx60; 		/* rule 135 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx61 = tx8*tzz142; 		/* rule 137 */
	 tx62 = tx8*tzz146; 		/* rule 138 */
	 tx63 = tx61 + tx62 + tzz154; 		/* rule 139 */
	 tx64 = tx8*tzz144; 		/* rule 140 */
	 tx65 = tx63*vdwScale; 		/* rule 141 */
	 dhy1y1 = tx58 + tx64 + tx65; 		/* rule 142 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx66 = tx9*tzz142; 		/* rule 144 */
	 tx67 = tx9*tzz146; 		/* rule 145 */
	 tx68 = tx66 + tx67 + tzz154; 		/* rule 146 */
	 tx69 = tx9*tzz144; 		/* rule 147 */
	 tx70 = tx68*vdwScale; 		/* rule 148 */
	 dhz1z1 = tx58 + tx69 + tx70; 		/* rule 149 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 dhx2x2 = dhx1x1; 		/* rule 151 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 dhy2y2 = dhy1y1; 		/* rule 153 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 dhz2z2 = dhz1z1; 		/* rule 155 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	#ifdef NONBOND_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tzz141 = tx4*tx5; 		/* rule 159 */
	 tx71 = tzz141*tzz142; 		/* rule 160 */
	 tx72 = tzz141*tzz146; 		/* rule 161 */
	 tx73 = tx71 + tx72; 		/* rule 162 */
	 tx74 = tzz141*tzz144; 		/* rule 163 */
	 tx75 = tx73*vdwScale; 		/* rule 164 */
	 ohx1y1 = tx74 + tx75; 		/* rule 165 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tzz140 = tx4*tx6; 		/* rule 167 */
	 tx76 = tzz140*tzz142; 		/* rule 168 */
	 tx77 = tzz140*tzz146; 		/* rule 169 */
	 tx78 = tx76 + tx77; 		/* rule 170 */
	 tx79 = tzz140*tzz144; 		/* rule 171 */
	 tx80 = tx78*vdwScale; 		/* rule 172 */
	 ohx1z1 = tx79 + tx80; 		/* rule 173 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx81 = tzz149; 		/* rule 175 */
	 tx82 = tzz151; 		/* rule 176 */
	 tzz147 = -168*tzz136; 		/* rule 177 */
	 tx83 = tx7*tzz147; 		/* rule 178 */
	 tzz143 = 48*tzz135; 		/* rule 179 */
	 tx84 = tx7*tzz143; 		/* rule 180 */
	 tzz153 = tx81 + tx82; 		/* rule 181 */
	 tx85 = tx83 + tx84 + tzz153; 		/* rule 182 */
	 tx86 = tzz148; 		/* rule 183 */
	 tzz145 = -3*tzz134; 		/* rule 184 */
	 tx87 = tx7*tzz145; 		/* rule 185 */
	 tx88 = tx85*vdwScale; 		/* rule 186 */
	 ohx1x2 = tx86 + tx87 + tx88; 		/* rule 187 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx89 = tzz141*tzz147; 		/* rule 189 */
	 tx90 = tzz141*tzz143; 		/* rule 190 */
	 tx91 = tx89 + tx90; 		/* rule 191 */
	 tx92 = tzz141*tzz145; 		/* rule 192 */
	 tx93 = tx91*vdwScale; 		/* rule 193 */
	 ohx1y2 = tx92 + tx93; 		/* rule 194 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx94 = tzz140*tzz147; 		/* rule 196 */
	 tx95 = tzz140*tzz143; 		/* rule 197 */
	 tx96 = tx94 + tx95; 		/* rule 198 */
	 tx97 = tzz140*tzz145; 		/* rule 199 */
	 tx98 = tx96*vdwScale; 		/* rule 200 */
	 ohx1z2 = tx97 + tx98; 		/* rule 201 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tzz139 = tx5*tx6; 		/* rule 203 */
	 tx99 = tzz139*tzz142; 		/* rule 204 */
	 tx100 = tzz139*tzz146; 		/* rule 205 */
	 tx101 = tx100 + tx99; 		/* rule 206 */
	 tx102 = tzz139*tzz144; 		/* rule 207 */
	 tx103 = tx101*vdwScale; 		/* rule 208 */
	 ohy1z1 = tx102 + tx103; 		/* rule 209 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 ohy1x2 = ohx1y2; 		/* rule 211 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx104 = tx8*tzz147; 		/* rule 213 */
	 tx105 = tx8*tzz143; 		/* rule 214 */
	 tx106 = tx104 + tx105 + tzz153; 		/* rule 215 */
	 tx107 = tx8*tzz145; 		/* rule 216 */
	 tx108 = tx106*vdwScale; 		/* rule 217 */
	 ohy1y2 = tx107 + tx108 + tx86; 		/* rule 218 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx109 = tzz139*tzz147; 		/* rule 220 */
	 tx110 = tzz139*tzz143; 		/* rule 221 */
	 tx111 = tx109 + tx110; 		/* rule 222 */
	 tx112 = tzz139*tzz145; 		/* rule 223 */
	 tx113 = tx111*vdwScale; 		/* rule 224 */
	 ohy1z2 = tx112 + tx113; 		/* rule 225 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 ohz1x2 = ohx1z2; 		/* rule 227 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 ohz1y2 = ohy1z2; 		/* rule 229 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx114 = tx9*tzz147; 		/* rule 231 */
	 tx115 = tx9*tzz143; 		/* rule 232 */
	 tx116 = tx114 + tx115 + tzz153; 		/* rule 233 */
	 tx117 = tx9*tzz145; 		/* rule 234 */
	 tx118 = tx116*vdwScale; 		/* rule 235 */
	 ohz1z2 = tx117 + tx118 + tx86; 		/* rule 236 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 ohx2y2 = ohx1y1; 		/* rule 238 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 ohx2z2 = ohx1z1; 		/* rule 240 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 ohy2z2 = ohy1z1; 		/* rule 242 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	} /*if calcOffDiagonalHessian */ 
	#endif /* NONBOND_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* NONBOND_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* NONBOND_CALC_FORCE ]*/
	NONBOND_DEBUG_INTERACTIONS(I1,I2);
	SKIP_term: (void)0;
