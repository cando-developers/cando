// C-code
	POINT_TO_LINE_RESTRAINT_SET_PARAMETER(ka);
	POINT_TO_LINE_RESTRAINT_SET_PARAMETER(ra);
	POINT_TO_LINE_RESTRAINT_SET_PARAMETER(xa);
	POINT_TO_LINE_RESTRAINT_SET_PARAMETER(ya);
	POINT_TO_LINE_RESTRAINT_SET_PARAMETER(za);
	POINT_TO_LINE_RESTRAINT_SET_PARAMETER(xb);
	POINT_TO_LINE_RESTRAINT_SET_PARAMETER(yb);
	POINT_TO_LINE_RESTRAINT_SET_PARAMETER(zb);
	POINT_TO_LINE_RESTRAINT_SET_PARAMETER(I1);
	POINT_TO_LINE_RESTRAINT_SET_POSITION(x1,I1,0);
	POINT_TO_LINE_RESTRAINT_SET_POSITION(y1,I1,1);
	POINT_TO_LINE_RESTRAINT_SET_POSITION(z1,I1,2);
	 tx100 = -xa; 		/* rule 13 */
	 tx101 = tx100*y1; 		/* rule 14 */
	 tx102 = xb*y1; 		/* rule 15 */
	 tx103 = -ya; 		/* rule 16 */
	 tx104 = x1*ya; 		/* rule 17 */
	 tx105 = tx103*xb; 		/* rule 18 */
	 tx106 = -(x1*yb); 		/* rule 19 */
	 tx107 = xa*yb; 		/* rule 20 */
	 tx108 = xa*z1; 		/* rule 21 */
	 tx109 = -(xb*z1); 		/* rule 22 */
	 tx110 = tx103*z1; 		/* rule 23 */
	 tx111 = yb*z1; 		/* rule 24 */
	 tx112 = -za; 		/* rule 25 */
	 tx113 = tx112*x1; 		/* rule 26 */
	 tx114 = xb*za; 		/* rule 27 */
	 tx115 = y1*za; 		/* rule 28 */
	 tx116 = tx112*yb; 		/* rule 29 */
	 tx117 = x1*zb; 		/* rule 30 */
	 tx118 = tx100*zb; 		/* rule 31 */
	 tx119 = -(y1*zb); 		/* rule 32 */
	 tx120 = ya*zb; 		/* rule 33 */
	 tx121 = tx113; 		/* rule 34 */
	 tx122 = tx105; 		/* rule 35 */
	 tx123 = tx101; 		/* rule 36 */
	 tx124 = tx116; 		/* rule 37 */
	 tx125 = tx110; 		/* rule 38 */
	 tx126 = tx118; 		/* rule 39 */
	 tx127 = tx102 + tx104 + tx106 + tx107 + tx122 + tx123; 		/* rule 40 */
	 tx128 = tx111 + tx115 + tx119 + tx120 + tx124 + tx125; 		/* rule 41 */
	 tx129 = tx108 + tx109 + tx114 + tx117 + tx121 + tx126; 		/* rule 42 */
	 tx130 = tx100 + xb; 		/* rule 43 */
	 tx131 = tx103 + yb; 		/* rule 44 */
	 tx132 = tx112 + zb; 		/* rule 45 */
	 tx133 = power2(tx127); 		/* rule 46 */
	 tx134 = power2(tx128); 		/* rule 47 */
	 tx135 = power2(tx129); 		/* rule 48 */
	 tx136 = power2(tx130); 		/* rule 49 */
	 tx137 = power2(tx131); 		/* rule 50 */
	 tx138 = power2(tx132); 		/* rule 51 */
	 tx139 = tx133 + tx134 + tx135; 		/* rule 52 */
	 tx140 = tx136 + tx137 + tx138; 		/* rule 53 */
	 tx141 = mysqrt(tx139); 		/* rule 54 */
	 tx199 = mysqrt(tx140); 		/* rule 55 */
	 tx142 = reciprocal(tx199); 		/* rule 56 */
	 tx143 = -ra; 		/* rule 57 */
	 tx144 = tx141*tx142; 		/* rule 58 */
	 tx145 = tx143 + tx144; 		/* rule 59 */
	 tx146 = power2(tx145); 		/* rule 60 */
	 Energy = ka*tx146; 		/* rule 61 */
	POINT_TO_LINE_RESTRAINT_ENERGY_ACCUMULATE(Energy);
	#ifdef POINT_TO_LINE_RESTRAINT_CALC_FORCE //[
	if ( calcForce ) {
	 tx147 = -yb; 		/* rule 65 */
	 tx148 = tx147 + ya; 		/* rule 66 */
	 tzz216 = 2.*tx129; 		/* rule 67 */
	 tx149 = tx132*tzz216; 		/* rule 68 */
	 tzz218 = 2.*tx127; 		/* rule 69 */
	 tx150 = tx148*tzz218; 		/* rule 70 */
	 tx151 = reciprocal(tx141); 		/* rule 71 */
	 tx152 = tx149 + tx150; 		/* rule 72 */
	 tzz203 = tx142*tx145; 		/* rule 73 */
	 tzz204 = ka*tzz203; 		/* rule 74 */
	 tzz211 = tx151*tzz204; 		/* rule 75 */
	 gx1 = tx152*tzz211; 		/* rule 76 */
	 fx1 = -gx1; 		/* rule 77 */
	POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx153 = -zb; 		/* rule 79 */
	 tx154 = tx153 + za; 		/* rule 80 */
	 tx155 = tx130*tzz218; 		/* rule 81 */
	 tzz217 = 2.*tx128; 		/* rule 82 */
	 tx156 = tx154*tzz217; 		/* rule 83 */
	 tx157 = tx155 + tx156; 		/* rule 84 */
	 gy1 = tx157*tzz211; 		/* rule 85 */
	 fy1 = -gy1; 		/* rule 86 */
	POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx158 = -xb; 		/* rule 88 */
	 tx159 = tx158 + xa; 		/* rule 89 */
	 tx160 = tx131*tzz217; 		/* rule 90 */
	 tx161 = tx159*tzz216; 		/* rule 91 */
	 tx162 = tx160 + tx161; 		/* rule 92 */
	 gz1 = tx162*tzz211; 		/* rule 93 */
	 fz1 = -gz1; 		/* rule 94 */
	POINT_TO_LINE_RESTRAINT_FORCE_ACCUMULATE(I1, 2, fz1 );
	#ifdef POINT_TO_LINE_RESTRAINT_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx163 = power2(tx148); 		/* rule 98 */
	 tx164 = 2.*tx138; 		/* rule 99 */
	 tx165 = 2.*tx163; 		/* rule 100 */
	 tx200 = tx139; 		/* rule 101 */
	 tx201 = tx151; 		/* rule 102 */
	 tx202 = reciprocal(tx200); 		/* rule 103 */
	 tx166 = tx201*tx202; 		/* rule 104 */
	 tx167 = tx202; 		/* rule 105 */
	 tx168 = reciprocal(tx140); 		/* rule 106 */
	 tx169 = power2(tx152); 		/* rule 107 */
	 tx170 = tx164 + tx165; 		/* rule 108 */
	 tzz207 = tx166*tzz204; 		/* rule 109 */
	 tzz210 = -0.5*tzz207; 		/* rule 110 */
	 tx171 = tx169*tzz210; 		/* rule 111 */
	 tzz206 = tx167*tx168; 		/* rule 112 */
	 tzz208 = ka*tzz206; 		/* rule 113 */
	 tzz209 = 0.5*tzz208; 		/* rule 114 */
	 tx172 = tx169*tzz209; 		/* rule 115 */
	 tzz205 = tx201*tzz204; 		/* rule 116 */
	 tx173 = tx170*tzz205; 		/* rule 117 */
	 dhx1x1 = tx171 + tx172 + tx173; 		/* rule 118 */
	POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx174 = power2(tx154); 		/* rule 120 */
	 tx175 = 2.*tx136; 		/* rule 121 */
	 tx176 = 2.*tx174; 		/* rule 122 */
	 tx177 = power2(tx157); 		/* rule 123 */
	 tx178 = tx175 + tx176; 		/* rule 124 */
	 tx179 = tx177*tzz210; 		/* rule 125 */
	 tx180 = tx177*tzz209; 		/* rule 126 */
	 tx181 = tx178*tzz205; 		/* rule 127 */
	 dhy1y1 = tx179 + tx180 + tx181; 		/* rule 128 */
	POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx182 = power2(tx159); 		/* rule 130 */
	 tx183 = 2.*tx137; 		/* rule 131 */
	 tx184 = 2.*tx182; 		/* rule 132 */
	 tx185 = power2(tx162); 		/* rule 133 */
	 tx186 = tx183 + tx184; 		/* rule 134 */
	 tx187 = tx185*tzz210; 		/* rule 135 */
	 tx188 = tx185*tzz209; 		/* rule 136 */
	 tx189 = tx186*tzz205; 		/* rule 137 */
	 dhz1z1 = tx187 + tx188 + tx189; 		/* rule 138 */
	POINT_TO_LINE_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	#ifdef POINT_TO_LINE_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tzz212 = 2.*tzz205; 		/* rule 142 */
	 tx190 = tx130*tx148*tzz212; 		/* rule 143 */
	 tzz215 = tx152*tx157; 		/* rule 144 */
	 tx191 = tzz210*tzz215; 		/* rule 145 */
	 tx192 = tzz209*tzz215; 		/* rule 146 */
	 ohx1y1 = tx190 + tx191 + tx192; 		/* rule 147 */
	POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tx193 = tx132*tx159*tzz212; 		/* rule 149 */
	 tzz213 = tx162*tzz210; 		/* rule 150 */
	 tx194 = tx152*tzz213; 		/* rule 151 */
	 tzz214 = tx162*tzz209; 		/* rule 152 */
	 tx195 = tx152*tzz214; 		/* rule 153 */
	 ohx1z1 = tx193 + tx194 + tx195; 		/* rule 154 */
	POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx196 = tx131*tx154*tzz212; 		/* rule 156 */
	 tx197 = tx157*tzz213; 		/* rule 157 */
	 tx198 = tx157*tzz214; 		/* rule 158 */
	 ohy1z1 = tx196 + tx197 + tx198; 		/* rule 159 */
	POINT_TO_LINE_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	} /*if calcOffDiagonalHessian */ 
	#endif /* POINT_TO_LINE_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* POINT_TO_LINE_RESTRAINT_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* POINT_TO_LINE_RESTRAINT_CALC_FORCE ]*/
