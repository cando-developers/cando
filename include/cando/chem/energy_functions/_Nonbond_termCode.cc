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
	 tzz152 = -6.*tzz137; 		/* rule 58 */
	 tx20 = tx4*tzz152; 		/* rule 59 */
	 tx124 = DistanceSquared; 		/* rule 60 */
	 tx125 = tx16; 		/* rule 61 */
	 tx126 = reciprocal(tx124); 		/* rule 62 */
	 tx21 = tx125*tx126; 		/* rule 63 */
	 tx22 = tx19 + tx20; 		/* rule 64 */
	 tzz139 = tx21*tzz133; 		/* rule 65 */
	 tx23 = tx4*tzz139; 		/* rule 66 */
	 tx24 = tx22*vdwScale; 		/* rule 67 */
	 gx1 = tx23 + tx24; 		/* rule 68 */
	 fx1 = -gx1; 		/* rule 69 */
	NONBOND_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx25 = tx5*tzz149; 		/* rule 71 */
	 tx26 = tx5*tzz152; 		/* rule 72 */
	 tx27 = tx25 + tx26; 		/* rule 73 */
	 tx28 = tx5*tzz139; 		/* rule 74 */
	 tx29 = tx27*vdwScale; 		/* rule 75 */
	 gy1 = tx28 + tx29; 		/* rule 76 */
	 fy1 = -gy1; 		/* rule 77 */
	NONBOND_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx30 = tx6*tzz149; 		/* rule 79 */
	 tx31 = tx6*tzz152; 		/* rule 80 */
	 tx32 = tx30 + tx31; 		/* rule 81 */
	 tx33 = tx6*tzz139; 		/* rule 82 */
	 tx34 = tx32*vdwScale; 		/* rule 83 */
	 gz1 = tx33 + tx34; 		/* rule 84 */
	 fz1 = -gz1; 		/* rule 85 */
	NONBOND_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tzz153 = -12.*tzz138; 		/* rule 87 */
	 tx35 = tx4*tzz153; 		/* rule 88 */
	 tzz150 = 6.*tzz137; 		/* rule 89 */
	 tx36 = tx4*tzz150; 		/* rule 90 */
	 tx37 = tx35 + tx36; 		/* rule 91 */
	 tzz151 = -tzz139; 		/* rule 92 */
	 tx38 = tx4*tzz151; 		/* rule 93 */
	 tx39 = tx37*vdwScale; 		/* rule 94 */
	 gx2 = tx38 + tx39; 		/* rule 95 */
	 fx2 = -gx2; 		/* rule 96 */
	NONBOND_FORCE_ACCUMULATE(I2, 0, fx2 );
	 tx40 = tx5*tzz153; 		/* rule 98 */
	 tx41 = tx5*tzz150; 		/* rule 99 */
	 tx42 = tx40 + tx41; 		/* rule 100 */
	 tx43 = tx5*tzz151; 		/* rule 101 */
	 tx44 = tx42*vdwScale; 		/* rule 102 */
	 gy2 = tx43 + tx44; 		/* rule 103 */
	 fy2 = -gy2; 		/* rule 104 */
	NONBOND_FORCE_ACCUMULATE(I2, 1, fy2 );
	 tx45 = tx6*tzz153; 		/* rule 106 */
	 tx46 = tx6*tzz150; 		/* rule 107 */
	 tx47 = tx45 + tx46; 		/* rule 108 */
	 tx48 = tx6*tzz151; 		/* rule 109 */
	 tx49 = tx47*vdwScale; 		/* rule 110 */
	 gz2 = tx48 + tx49; 		/* rule 111 */
	 fz2 = -gz2; 		/* rule 112 */
	NONBOND_FORCE_ACCUMULATE(I2, 2, fz2 );
	#ifdef NONBOND_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx50 = power2(tx123); 		/* rule 116 */
	 tx51 = tx119*tx120; 		/* rule 117 */
	 tx52 = tzz153; 		/* rule 118 */
	 tx53 = tzz150; 		/* rule 119 */
	 tzz136 = dA*tx50; 		/* rule 120 */
	 tzz143 = 168*tzz136; 		/* rule 121 */
	 tx54 = tx7*tzz143; 		/* rule 122 */
	 tzz135 = dC*tx51; 		/* rule 123 */
	 tzz147 = -48*tzz135; 		/* rule 124 */
	 tx55 = tx7*tzz147; 		/* rule 125 */
	 tx127 = power2(tx124); 		/* rule 126 */
	 tx128 = reciprocal(tx127); 		/* rule 127 */
	 tx56 = tx125*tx128; 		/* rule 128 */
	 tzz155 = tx52 + tx53; 		/* rule 129 */
	 tx57 = tx54 + tx55 + tzz155; 		/* rule 130 */
	 tx58 = tzz151; 		/* rule 131 */
	 tzz134 = tx56*tzz133; 		/* rule 132 */
	 tzz145 = 3*tzz134; 		/* rule 133 */
	 tx59 = tx7*tzz145; 		/* rule 134 */
	 tx60 = tx57*vdwScale; 		/* rule 135 */
	 dhx1x1 = tx58 + tx59 + tx60; 		/* rule 136 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx61 = tx8*tzz143; 		/* rule 138 */
	 tx62 = tx8*tzz147; 		/* rule 139 */
	 tx63 = tx61 + tx62 + tzz155; 		/* rule 140 */
	 tx64 = tx8*tzz145; 		/* rule 141 */
	 tx65 = tx63*vdwScale; 		/* rule 142 */
	 dhy1y1 = tx58 + tx64 + tx65; 		/* rule 143 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx66 = tx9*tzz143; 		/* rule 145 */
	 tx67 = tx9*tzz147; 		/* rule 146 */
	 tx68 = tx66 + tx67 + tzz155; 		/* rule 147 */
	 tx69 = tx9*tzz145; 		/* rule 148 */
	 tx70 = tx68*vdwScale; 		/* rule 149 */
	 dhz1z1 = tx58 + tx69 + tx70; 		/* rule 150 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 dhx2x2 = dhx1x1; 		/* rule 152 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 dhy2y2 = dhy1y1; 		/* rule 154 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 dhz2z2 = dhz1z1; 		/* rule 156 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	#ifdef NONBOND_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tzz142 = tx4*tx5; 		/* rule 160 */
	 tx71 = tzz142*tzz143; 		/* rule 161 */
	 tx72 = tzz142*tzz147; 		/* rule 162 */
	 tx73 = tx71 + tx72; 		/* rule 163 */
	 tx74 = tzz142*tzz145; 		/* rule 164 */
	 tx75 = tx73*vdwScale; 		/* rule 165 */
	 ohx1y1 = tx74 + tx75; 		/* rule 166 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tzz141 = tx4*tx6; 		/* rule 168 */
	 tx76 = tzz141*tzz143; 		/* rule 169 */
	 tx77 = tzz141*tzz147; 		/* rule 170 */
	 tx78 = tx76 + tx77; 		/* rule 171 */
	 tx79 = tzz141*tzz145; 		/* rule 172 */
	 tx80 = tx78*vdwScale; 		/* rule 173 */
	 ohx1z1 = tx79 + tx80; 		/* rule 174 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx81 = tzz149; 		/* rule 176 */
	 tx82 = tzz152; 		/* rule 177 */
	 tzz148 = -168*tzz136; 		/* rule 178 */
	 tx83 = tx7*tzz148; 		/* rule 179 */
	 tzz144 = 48*tzz135; 		/* rule 180 */
	 tx84 = tx7*tzz144; 		/* rule 181 */
	 tzz154 = tx81 + tx82; 		/* rule 182 */
	 tx85 = tx83 + tx84 + tzz154; 		/* rule 183 */
	 tx86 = tzz139; 		/* rule 184 */
	 tzz146 = -3*tzz134; 		/* rule 185 */
	 tx87 = tx7*tzz146; 		/* rule 186 */
	 tx88 = tx85*vdwScale; 		/* rule 187 */
	 ohx1x2 = tx86 + tx87 + tx88; 		/* rule 188 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx89 = tzz142*tzz148; 		/* rule 190 */
	 tx90 = tzz142*tzz144; 		/* rule 191 */
	 tx91 = tx89 + tx90; 		/* rule 192 */
	 tx92 = tzz142*tzz146; 		/* rule 193 */
	 tx93 = tx91*vdwScale; 		/* rule 194 */
	 ohx1y2 = tx92 + tx93; 		/* rule 195 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx94 = tzz141*tzz148; 		/* rule 197 */
	 tx95 = tzz141*tzz144; 		/* rule 198 */
	 tx96 = tx94 + tx95; 		/* rule 199 */
	 tx97 = tzz141*tzz146; 		/* rule 200 */
	 tx98 = tx96*vdwScale; 		/* rule 201 */
	 ohx1z2 = tx97 + tx98; 		/* rule 202 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tzz140 = tx5*tx6; 		/* rule 204 */
	 tx99 = tzz140*tzz143; 		/* rule 205 */
	 tx100 = tzz140*tzz147; 		/* rule 206 */
	 tx101 = tx100 + tx99; 		/* rule 207 */
	 tx102 = tzz140*tzz145; 		/* rule 208 */
	 tx103 = tx101*vdwScale; 		/* rule 209 */
	 ohy1z1 = tx102 + tx103; 		/* rule 210 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 ohy1x2 = ohx1y2; 		/* rule 212 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx104 = tx8*tzz148; 		/* rule 214 */
	 tx105 = tx8*tzz144; 		/* rule 215 */
	 tx106 = tx104 + tx105 + tzz154; 		/* rule 216 */
	 tx107 = tx8*tzz146; 		/* rule 217 */
	 tx108 = tx106*vdwScale; 		/* rule 218 */
	 ohy1y2 = tx107 + tx108 + tx86; 		/* rule 219 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx109 = tzz140*tzz148; 		/* rule 221 */
	 tx110 = tzz140*tzz144; 		/* rule 222 */
	 tx111 = tx109 + tx110; 		/* rule 223 */
	 tx112 = tzz140*tzz146; 		/* rule 224 */
	 tx113 = tx111*vdwScale; 		/* rule 225 */
	 ohy1z2 = tx112 + tx113; 		/* rule 226 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 ohz1x2 = ohx1z2; 		/* rule 228 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 ohz1y2 = ohy1z2; 		/* rule 230 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx114 = tx9*tzz148; 		/* rule 232 */
	 tx115 = tx9*tzz144; 		/* rule 233 */
	 tx116 = tx114 + tx115 + tzz154; 		/* rule 234 */
	 tx117 = tx9*tzz146; 		/* rule 235 */
	 tx118 = tx116*vdwScale; 		/* rule 236 */
	 ohz1z2 = tx117 + tx118 + tx86; 		/* rule 237 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 ohx2y2 = ohx1y1; 		/* rule 239 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 ohx2z2 = ohx1z1; 		/* rule 241 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 ohy2z2 = ohy1z1; 		/* rule 243 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	} /*if calcOffDiagonalHessian */ 
	#endif /* NONBOND_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* NONBOND_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* NONBOND_CALC_FORCE ]*/
	NONBOND_DEBUG_INTERACTIONS(I1,I2);
	SKIP_term: (void)0;
