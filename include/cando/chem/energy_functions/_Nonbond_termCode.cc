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
	 R = mysqrt(DistanceSquared); 		/* rule 26 */
	BAIL_OUT_IF_CUTOFF(DistanceSquared)
	 tx130 = reciprocal(DistanceSquared); 		/* rule 28 */
	 tx128 = power2(tx130); 		/* rule 29 */
	 tx129 = tx130; 		/* rule 30 */
	 tx119 = tx128*tx129; 		/* rule 31 */
	 tx10 = power2(tx119); 		/* rule 32 */
	 tx120 = power2(tx129); 		/* rule 33 */
	 tx121 = tx129; 		/* rule 34 */
	 tx11 = tx120*tx121; 		/* rule 35 */
	 tx12 = dA*tx10; 		/* rule 36 */
	 tx13 = -(dC*tx11); 		/* rule 37 */
	 tx14 = tx12 + tx13; 		/* rule 38 */
	 Evdw = tx14*vdwScale; 		/* rule 39 */
	NONBOND_EVDW_ENERGY_ACCUMULATE(Evdw);
	 tx15 = reciprocal(DIELECTRIC); 		/* rule 41 */
	 tx16 = reciprocal(R); 		/* rule 42 */
	 Eeel = dQ1Q2*eelScale*tx15*tx16; 		/* rule 43 */
	NONBOND_EEEL_ENERGY_ACCUMULATE(Eeel);
	 Energy = Eeel + Evdw; 		/* rule 45 */
	NONBOND_ENERGY_ACCUMULATE(Energy);
	#ifdef NONBOND_CALC_FORCE //[
	if ( calcForce ) {
	 tx122 = power2(tx128); 		/* rule 49 */
	 tx17 = tx119*tx122; 		/* rule 50 */
	 tx18 = power2(tx120); 		/* rule 51 */
	 tx19 = 12.*dA*tx17*tx4; 		/* rule 52 */
	 tx20 = -6.*dC*tx18*tx4; 		/* rule 53 */
	 tx123 = DistanceSquared; 		/* rule 54 */
	 tx124 = tx16; 		/* rule 55 */
	 tx125 = reciprocal(tx123); 		/* rule 56 */
	 tx21 = tx124*tx125; 		/* rule 57 */
	 tx22 = tx19 + tx20; 		/* rule 58 */
	 tx23 = dQ1Q2*eelScale*tx15*tx21*tx4; 		/* rule 59 */
	 tx24 = tx22*vdwScale; 		/* rule 60 */
	 gx1 = tx23 + tx24; 		/* rule 61 */
	 fx1 = -gx1; 		/* rule 62 */
	NONBOND_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx25 = 12.*dA*tx17*tx5; 		/* rule 64 */
	 tx26 = -6.*dC*tx18*tx5; 		/* rule 65 */
	 tx27 = tx25 + tx26; 		/* rule 66 */
	 tx28 = dQ1Q2*eelScale*tx15*tx21*tx5; 		/* rule 67 */
	 tx29 = tx27*vdwScale; 		/* rule 68 */
	 gy1 = tx28 + tx29; 		/* rule 69 */
	 fy1 = -gy1; 		/* rule 70 */
	NONBOND_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx30 = 12.*dA*tx17*tx6; 		/* rule 72 */
	 tx31 = -6.*dC*tx18*tx6; 		/* rule 73 */
	 tx32 = tx30 + tx31; 		/* rule 74 */
	 tx33 = dQ1Q2*eelScale*tx15*tx21*tx6; 		/* rule 75 */
	 tx34 = tx32*vdwScale; 		/* rule 76 */
	 gz1 = tx33 + tx34; 		/* rule 77 */
	 fz1 = -gz1; 		/* rule 78 */
	NONBOND_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tx35 = -12.*dA*tx17*tx4; 		/* rule 80 */
	 tx36 = 6.*dC*tx18*tx4; 		/* rule 81 */
	 tx37 = tx35 + tx36; 		/* rule 82 */
	 tx38 = -tx23; 		/* rule 83 */
	 tx39 = tx37*vdwScale; 		/* rule 84 */
	 gx2 = tx38 + tx39; 		/* rule 85 */
	 fx2 = -gx2; 		/* rule 86 */
	NONBOND_FORCE_ACCUMULATE(I2, 0, fx2 );
	 tx40 = -12.*dA*tx17*tx5; 		/* rule 88 */
	 tx41 = 6.*dC*tx18*tx5; 		/* rule 89 */
	 tx42 = tx40 + tx41; 		/* rule 90 */
	 tx43 = -tx28; 		/* rule 91 */
	 tx44 = tx42*vdwScale; 		/* rule 92 */
	 gy2 = tx43 + tx44; 		/* rule 93 */
	 fy2 = -gy2; 		/* rule 94 */
	NONBOND_FORCE_ACCUMULATE(I2, 1, fy2 );
	 tx45 = -12.*dA*tx17*tx6; 		/* rule 96 */
	 tx46 = 6.*dC*tx18*tx6; 		/* rule 97 */
	 tx47 = tx45 + tx46; 		/* rule 98 */
	 tx48 = -tx33; 		/* rule 99 */
	 tx49 = tx47*vdwScale; 		/* rule 100 */
	 gz2 = tx48 + tx49; 		/* rule 101 */
	 fz2 = -gz2; 		/* rule 102 */
	NONBOND_FORCE_ACCUMULATE(I2, 2, fz2 );
	#ifdef NONBOND_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx50 = power2(tx122); 		/* rule 106 */
	 tx51 = tx119*tx120; 		/* rule 107 */
	 tx52 = -12.*dA*tx17; 		/* rule 108 */
	 tx53 = 6.*dC*tx18; 		/* rule 109 */
	 tx54 = 168*dA*tx50*tx7; 		/* rule 110 */
	 tx55 = -48*dC*tx51*tx7; 		/* rule 111 */
	 tx126 = power2(tx123); 		/* rule 112 */
	 tx127 = reciprocal(tx126); 		/* rule 113 */
	 tx56 = tx124*tx127; 		/* rule 114 */
	 tx57 = tx52 + tx53 + tx54 + tx55; 		/* rule 115 */
	 tx58 = -(dQ1Q2*eelScale*tx15*tx21); 		/* rule 116 */
	 tx59 = 3*dQ1Q2*eelScale*tx15*tx56*tx7; 		/* rule 117 */
	 tx60 = tx57*vdwScale; 		/* rule 118 */
	 dhx1x1 = tx58 + tx59 + tx60; 		/* rule 119 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx61 = 168*dA*tx50*tx8; 		/* rule 121 */
	 tx62 = -48*dC*tx51*tx8; 		/* rule 122 */
	 tx63 = tx52 + tx53 + tx61 + tx62; 		/* rule 123 */
	 tx64 = 3*dQ1Q2*eelScale*tx15*tx56*tx8; 		/* rule 124 */
	 tx65 = tx63*vdwScale; 		/* rule 125 */
	 dhy1y1 = tx58 + tx64 + tx65; 		/* rule 126 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx66 = 168*dA*tx50*tx9; 		/* rule 128 */
	 tx67 = -48*dC*tx51*tx9; 		/* rule 129 */
	 tx68 = tx52 + tx53 + tx66 + tx67; 		/* rule 130 */
	 tx69 = 3*dQ1Q2*eelScale*tx15*tx56*tx9; 		/* rule 131 */
	 tx70 = tx68*vdwScale; 		/* rule 132 */
	 dhz1z1 = tx58 + tx69 + tx70; 		/* rule 133 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 dhx2x2 = dhx1x1; 		/* rule 135 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 dhy2y2 = dhy1y1; 		/* rule 137 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 dhz2z2 = dhz1z1; 		/* rule 139 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	#ifdef NONBOND_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tx71 = 168*dA*tx4*tx5*tx50; 		/* rule 143 */
	 tx72 = -48*dC*tx4*tx5*tx51; 		/* rule 144 */
	 tx73 = tx71 + tx72; 		/* rule 145 */
	 tx74 = 3*dQ1Q2*eelScale*tx15*tx4*tx5*tx56; 		/* rule 146 */
	 tx75 = tx73*vdwScale; 		/* rule 147 */
	 ohx1y1 = tx74 + tx75; 		/* rule 148 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tx76 = 168*dA*tx4*tx50*tx6; 		/* rule 150 */
	 tx77 = -48*dC*tx4*tx51*tx6; 		/* rule 151 */
	 tx78 = tx76 + tx77; 		/* rule 152 */
	 tx79 = 3*dQ1Q2*eelScale*tx15*tx4*tx56*tx6; 		/* rule 153 */
	 tx80 = tx78*vdwScale; 		/* rule 154 */
	 ohx1z1 = tx79 + tx80; 		/* rule 155 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx81 = 12.*dA*tx17; 		/* rule 157 */
	 tx82 = -6.*dC*tx18; 		/* rule 158 */
	 tx83 = -168*dA*tx50*tx7; 		/* rule 159 */
	 tx84 = 48*dC*tx51*tx7; 		/* rule 160 */
	 tx85 = tx81 + tx82 + tx83 + tx84; 		/* rule 161 */
	 tx86 = dQ1Q2*eelScale*tx15*tx21; 		/* rule 162 */
	 tx87 = -3*dQ1Q2*eelScale*tx15*tx56*tx7; 		/* rule 163 */
	 tx88 = tx85*vdwScale; 		/* rule 164 */
	 ohx1x2 = tx86 + tx87 + tx88; 		/* rule 165 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx89 = -168*dA*tx4*tx5*tx50; 		/* rule 167 */
	 tx90 = 48*dC*tx4*tx5*tx51; 		/* rule 168 */
	 tx91 = tx89 + tx90; 		/* rule 169 */
	 tx92 = -3*dQ1Q2*eelScale*tx15*tx4*tx5*tx56; 		/* rule 170 */
	 tx93 = tx91*vdwScale; 		/* rule 171 */
	 ohx1y2 = tx92 + tx93; 		/* rule 172 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx94 = -168*dA*tx4*tx50*tx6; 		/* rule 174 */
	 tx95 = 48*dC*tx4*tx51*tx6; 		/* rule 175 */
	 tx96 = tx94 + tx95; 		/* rule 176 */
	 tx97 = -3*dQ1Q2*eelScale*tx15*tx4*tx56*tx6; 		/* rule 177 */
	 tx98 = tx96*vdwScale; 		/* rule 178 */
	 ohx1z2 = tx97 + tx98; 		/* rule 179 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tx99 = 168*dA*tx5*tx50*tx6; 		/* rule 181 */
	 tx100 = -48*dC*tx5*tx51*tx6; 		/* rule 182 */
	 tx101 = tx100 + tx99; 		/* rule 183 */
	 tx102 = 3*dQ1Q2*eelScale*tx15*tx5*tx56*tx6; 		/* rule 184 */
	 tx103 = tx101*vdwScale; 		/* rule 185 */
	 ohy1z1 = tx102 + tx103; 		/* rule 186 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 ohy1x2 = ohx1y2; 		/* rule 188 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx104 = -168*dA*tx50*tx8; 		/* rule 190 */
	 tx105 = 48*dC*tx51*tx8; 		/* rule 191 */
	 tx106 = tx104 + tx105 + tx81 + tx82; 		/* rule 192 */
	 tx107 = -3*dQ1Q2*eelScale*tx15*tx56*tx8; 		/* rule 193 */
	 tx108 = tx106*vdwScale; 		/* rule 194 */
	 ohy1y2 = tx107 + tx108 + tx86; 		/* rule 195 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx109 = -168*dA*tx5*tx50*tx6; 		/* rule 197 */
	 tx110 = 48*dC*tx5*tx51*tx6; 		/* rule 198 */
	 tx111 = tx109 + tx110; 		/* rule 199 */
	 tx112 = -3*dQ1Q2*eelScale*tx15*tx5*tx56*tx6; 		/* rule 200 */
	 tx113 = tx111*vdwScale; 		/* rule 201 */
	 ohy1z2 = tx112 + tx113; 		/* rule 202 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 ohz1x2 = ohx1z2; 		/* rule 204 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 ohz1y2 = ohy1z2; 		/* rule 206 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx114 = -168*dA*tx50*tx9; 		/* rule 208 */
	 tx115 = 48*dC*tx51*tx9; 		/* rule 209 */
	 tx116 = tx114 + tx115 + tx81 + tx82; 		/* rule 210 */
	 tx117 = -3*dQ1Q2*eelScale*tx15*tx56*tx9; 		/* rule 211 */
	 tx118 = tx116*vdwScale; 		/* rule 212 */
	 ohz1z2 = tx117 + tx118 + tx86; 		/* rule 213 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 ohx2y2 = ohx1y1; 		/* rule 215 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 ohx2z2 = ohx1z1; 		/* rule 217 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 ohy2z2 = ohy1z1; 		/* rule 219 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	} /*if calcOffDiagonalHessian */ 
	#endif /* NONBOND_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* NONBOND_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* NONBOND_CALC_FORCE ]*/
	NONBOND_DEBUG_INTERACTIONS(I1,I2);
	SKIP_term: (void)0;
