// C-code
	NONBOND_SET_PARAMETER(I1);
	NONBOND_SET_PARAMETER(I2);
	NONBOND_SET_PARAMETER(dQ1Q2);
	NONBOND_SET_PARAMETER(dA);
	NONBOND_SET_PARAMETER(dC);
	NONBOND_SET_POSITION(x1,I1,0);
	NONBOND_SET_POSITION(y1,I1,1);
	NONBOND_SET_POSITION(z1,I1,2);
	NONBOND_SET_POSITION(x2,I2,0);
	NONBOND_SET_POSITION(y2,I2,1);
	NONBOND_SET_POSITION(z2,I2,2);
	 tx81 = -x1; 		/* rule 12 */
	 DeltaX = tx81 + x2; 		/* rule 13 */
	 tx82 = -y1; 		/* rule 14 */
	 DeltaY = tx82 + y2; 		/* rule 15 */
	 tx83 = -z1; 		/* rule 16 */
	 DeltaZ = tx83 + z2; 		/* rule 17 */
	 tx84 = PBX(DeltaX); 		/* rule 18 */
	 tx85 = PBY(DeltaY); 		/* rule 19 */
	 tx86 = PBZ(DeltaZ); 		/* rule 20 */
	 tx87 = power2(tx84); 		/* rule 21 */
	 tx88 = power2(tx85); 		/* rule 22 */
	 tx89 = power2(tx86); 		/* rule 23 */
	 DistanceSquared = tx87 + tx88 + tx89; 		/* rule 24 */
	BAIL_OUT_IF_CUTOFF(DistanceSquared)
	 tx192 = reciprocal(DistanceSquared); 		/* rule 26 */
	 tx171 = power2(tx192); 		/* rule 27 */
	 tx172 = tx192; 		/* rule 28 */
	 tx161 = tx171*tx172; 		/* rule 29 */
	 tx90 = power2(tx161); 		/* rule 30 */
	 tx162 = power2(tx172); 		/* rule 31 */
	 tx163 = tx172; 		/* rule 32 */
	 tx91 = tx162*tx163; 		/* rule 33 */
	 tx92 = dA*tx90; 		/* rule 34 */
	 tx93 = -(dC*tx91); 		/* rule 35 */
	 Evdw = tx92 + tx93; 		/* rule 36 */
	NONBOND_EVDW_ENERGY_ACCUMULATE(Evdw);
	 tx164 = mysqrt(DistanceSquared); 		/* rule 38 */
	 tx94 = reciprocal(tx164); 		/* rule 39 */
	 Eeel = dQ1Q2*tx94; 		/* rule 40 */
	NONBOND_EEEL_ENERGY_ACCUMULATE(Eeel);
	 Energy = Eeel + Evdw; 		/* rule 42 */
	NONBOND_ENERGY_ACCUMULATE(Energy);
	#ifdef NONBOND_CALC_FORCE //[
	if ( calcForce ) {
	 tx165 = power2(tx171); 		/* rule 46 */
	 tx95 = tx161*tx165; 		/* rule 47 */
	 tx96 = power2(tx162); 		/* rule 48 */
	 tx166 = DistanceSquared; 		/* rule 49 */
	 tx167 = tx94; 		/* rule 50 */
	 tx168 = reciprocal(tx166); 		/* rule 51 */
	 tx97 = tx167*tx168; 		/* rule 52 */
	 tzz177 = dA*tx95; 		/* rule 53 */
	 tzz188 = 12.*tzz177; 		/* rule 54 */
	 tx98 = tx84*tzz188; 		/* rule 55 */
	 tzz176 = dC*tx96; 		/* rule 56 */
	 tzz190 = -6.*tzz176; 		/* rule 57 */
	 tx99 = tx84*tzz190; 		/* rule 58 */
	 tzz187 = dQ1Q2*tx97; 		/* rule 59 */
	 tx100 = tx84*tzz187; 		/* rule 60 */
	 gx1 = tx100 + tx98 + tx99; 		/* rule 61 */
	 fx1 = -gx1; 		/* rule 62 */
	NONBOND_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx101 = tx85*tzz188; 		/* rule 64 */
	 tx102 = tx85*tzz190; 		/* rule 65 */
	 tx103 = tx85*tzz187; 		/* rule 66 */
	 gy1 = tx101 + tx102 + tx103; 		/* rule 67 */
	 fy1 = -gy1; 		/* rule 68 */
	NONBOND_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx104 = tx86*tzz188; 		/* rule 70 */
	 tx105 = tx86*tzz190; 		/* rule 71 */
	 tx106 = tx86*tzz187; 		/* rule 72 */
	 gz1 = tx104 + tx105 + tx106; 		/* rule 73 */
	 fz1 = -gz1; 		/* rule 74 */
	NONBOND_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tx107 = -tx100; 		/* rule 76 */
	 tzz191 = -12.*tzz177; 		/* rule 77 */
	 tx108 = tx84*tzz191; 		/* rule 78 */
	 tzz189 = 6.*tzz176; 		/* rule 79 */
	 tx109 = tx84*tzz189; 		/* rule 80 */
	 gx2 = tx107 + tx108 + tx109; 		/* rule 81 */
	 fx2 = -gx2; 		/* rule 82 */
	NONBOND_FORCE_ACCUMULATE(I2, 0, fx2 );
	 tx110 = -tx103; 		/* rule 84 */
	 tx111 = tx85*tzz191; 		/* rule 85 */
	 tx112 = tx85*tzz189; 		/* rule 86 */
	 gy2 = tx110 + tx111 + tx112; 		/* rule 87 */
	 fy2 = -gy2; 		/* rule 88 */
	NONBOND_FORCE_ACCUMULATE(I2, 1, fy2 );
	 tx113 = -tx106; 		/* rule 90 */
	 tx114 = tx86*tzz191; 		/* rule 91 */
	 tx115 = tx86*tzz189; 		/* rule 92 */
	 gz2 = tx113 + tx114 + tx115; 		/* rule 93 */
	 fz2 = -gz2; 		/* rule 94 */
	NONBOND_FORCE_ACCUMULATE(I2, 2, fz2 );
	#ifdef NONBOND_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx116 = power2(tx165); 		/* rule 98 */
	 tx117 = tx161*tx162; 		/* rule 99 */
	 tx169 = power2(tx166); 		/* rule 100 */
	 tx170 = reciprocal(tx169); 		/* rule 101 */
	 tx118 = tx167*tx170; 		/* rule 102 */
	 tzz175 = dA*tx116; 		/* rule 103 */
	 tzz181 = 168*tzz175; 		/* rule 104 */
	 tx119 = tx87*tzz181; 		/* rule 105 */
	 tzz174 = dC*tx117; 		/* rule 106 */
	 tzz185 = -48*tzz174; 		/* rule 107 */
	 tx120 = tx87*tzz185; 		/* rule 108 */
	 tzz173 = dQ1Q2*tx118; 		/* rule 109 */
	 tzz183 = 3*tzz173; 		/* rule 110 */
	 tx121 = tx87*tzz183; 		/* rule 111 */
	 tx122 = tzz191; 		/* rule 112 */
	 tx123 = tzz189; 		/* rule 113 */
	 tx124 = -tzz187; 		/* rule 114 */
	 tzz195 = tx123 + tx124; 		/* rule 115 */
	 tzz196 = tx122 + tzz195; 		/* rule 116 */
	 dhx1x1 = tx119 + tx120 + tx121 + tzz196; 		/* rule 117 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx125 = tx88*tzz181; 		/* rule 119 */
	 tx126 = tx88*tzz185; 		/* rule 120 */
	 tx127 = tx88*tzz183; 		/* rule 121 */
	 dhy1y1 = tx125 + tx126 + tx127 + tzz196; 		/* rule 122 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx128 = tx89*tzz181; 		/* rule 124 */
	 tx129 = tx89*tzz185; 		/* rule 125 */
	 tx130 = tx89*tzz183; 		/* rule 126 */
	 dhz1z1 = tx128 + tx129 + tx130 + tzz196; 		/* rule 127 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 dhx2x2 = dhx1x1; 		/* rule 129 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 dhy2y2 = dhy1y1; 		/* rule 131 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 dhz2z2 = dhz1z1; 		/* rule 133 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	#ifdef NONBOND_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tzz180 = tx84*tx85; 		/* rule 137 */
	 tx131 = tzz180*tzz181; 		/* rule 138 */
	 tx132 = tzz180*tzz185; 		/* rule 139 */
	 tx133 = tzz180*tzz183; 		/* rule 140 */
	 ohx1y1 = tx131 + tx132 + tx133; 		/* rule 141 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tzz179 = tx84*tx86; 		/* rule 143 */
	 tx134 = tzz179*tzz181; 		/* rule 144 */
	 tx135 = tzz179*tzz185; 		/* rule 145 */
	 tx136 = tzz179*tzz183; 		/* rule 146 */
	 ohx1z1 = tx134 + tx135 + tx136; 		/* rule 147 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tzz186 = -168*tzz175; 		/* rule 149 */
	 tx137 = tx87*tzz186; 		/* rule 150 */
	 tzz182 = 48*tzz174; 		/* rule 151 */
	 tx138 = tx87*tzz182; 		/* rule 152 */
	 tzz184 = -3*tzz173; 		/* rule 153 */
	 tx139 = tx87*tzz184; 		/* rule 154 */
	 tx140 = tzz188; 		/* rule 155 */
	 tx141 = tzz190; 		/* rule 156 */
	 tx142 = tzz187; 		/* rule 157 */
	 tzz193 = tx141 + tx142; 		/* rule 158 */
	 tzz194 = tx140 + tzz193; 		/* rule 159 */
	 ohx1x2 = tx137 + tx138 + tx139 + tzz194; 		/* rule 160 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx143 = tzz180*tzz186; 		/* rule 162 */
	 tx144 = tzz180*tzz182; 		/* rule 163 */
	 tx145 = tzz180*tzz184; 		/* rule 164 */
	 ohx1y2 = tx143 + tx144 + tx145; 		/* rule 165 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx146 = tzz179*tzz186; 		/* rule 167 */
	 tx147 = tzz179*tzz182; 		/* rule 168 */
	 tx148 = tzz179*tzz184; 		/* rule 169 */
	 ohx1z2 = tx146 + tx147 + tx148; 		/* rule 170 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tzz178 = tx85*tx86; 		/* rule 172 */
	 tx149 = tzz178*tzz181; 		/* rule 173 */
	 tx150 = tzz178*tzz185; 		/* rule 174 */
	 tx151 = tzz178*tzz183; 		/* rule 175 */
	 ohy1z1 = tx149 + tx150 + tx151; 		/* rule 176 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 ohy1x2 = ohx1y2; 		/* rule 178 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx152 = tx88*tzz186; 		/* rule 180 */
	 tx153 = tx88*tzz182; 		/* rule 181 */
	 tx154 = tx88*tzz184; 		/* rule 182 */
	 ohy1y2 = tx152 + tx153 + tx154 + tzz194; 		/* rule 183 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx155 = tzz178*tzz186; 		/* rule 185 */
	 tx156 = tzz178*tzz182; 		/* rule 186 */
	 tx157 = tzz178*tzz184; 		/* rule 187 */
	 ohy1z2 = tx155 + tx156 + tx157; 		/* rule 188 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 ohz1x2 = ohx1z2; 		/* rule 190 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 ohz1y2 = ohy1z2; 		/* rule 192 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx158 = tx89*tzz186; 		/* rule 194 */
	 tx159 = tx89*tzz182; 		/* rule 195 */
	 tx160 = tx89*tzz184; 		/* rule 196 */
	 ohz1z2 = tx158 + tx159 + tx160 + tzz194; 		/* rule 197 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 ohx2y2 = ohx1y1; 		/* rule 199 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 ohx2z2 = ohx1z1; 		/* rule 201 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 ohy2z2 = ohy1z1; 		/* rule 203 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	} /*if calcOffDiagonalHessian */ 
	#endif /* NONBOND_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* NONBOND_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* NONBOND_CALC_FORCE ]*/
