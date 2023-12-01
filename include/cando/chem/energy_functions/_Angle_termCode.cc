// C-code
	ANGLE_SET_PARAMETER(kt);
	ANGLE_SET_PARAMETER(t0);
	ANGLE_SET_PARAMETER(I1);
	ANGLE_SET_PARAMETER(I2);
	ANGLE_SET_PARAMETER(I3);
	ANGLE_APPLY_ATOM_MASK(I1,I2,I3);
	ANGLE_SET_POSITION(x1,I1,0);
	ANGLE_SET_POSITION(y1,I1,1);
	ANGLE_SET_POSITION(z1,I1,2);
	ANGLE_SET_POSITION(x2,I2,0);
	ANGLE_SET_POSITION(y2,I2,1);
	ANGLE_SET_POSITION(z2,I2,2);
	ANGLE_SET_POSITION(x3,I3,0);
	ANGLE_SET_POSITION(y3,I3,1);
	ANGLE_SET_POSITION(z3,I3,2);
	 tx1 = -x2; 		/* rule 16 */
	 tx2 = -y2; 		/* rule 17 */
	 tx3 = -z2; 		/* rule 18 */
	 tx4 = tx1 + x1; 		/* rule 19 */
	 tx5 = tx1 + x3; 		/* rule 20 */
	 tx6 = tx2 + y1; 		/* rule 21 */
	 tx7 = tx2 + y3; 		/* rule 22 */
	 tx8 = tx3 + z1; 		/* rule 23 */
	 tx9 = tx3 + z3; 		/* rule 24 */
	 tx10 = power2(tx4); 		/* rule 25 */
	 tx11 = power2(tx5); 		/* rule 26 */
	 tx12 = power2(tx6); 		/* rule 27 */
	 tx13 = power2(tx7); 		/* rule 28 */
	 tx14 = power2(tx8); 		/* rule 29 */
	 tx15 = power2(tx9); 		/* rule 30 */
	 tx16 = tx10 + tx12 + tx14; 		/* rule 31 */
	 tx17 = tx11 + tx13 + tx15; 		/* rule 32 */
	 tx18 = tx4*tx5; 		/* rule 33 */
	 tx19 = tx6*tx7; 		/* rule 34 */
	 tx20 = tx8*tx9; 		/* rule 35 */
	 tx450 = mysqrt(tx16); 		/* rule 36 */
	 tx21 = reciprocal(tx450); 		/* rule 37 */
	 tx451 = mysqrt(tx17); 		/* rule 38 */
	 tx22 = reciprocal(tx451); 		/* rule 39 */
	 tx23 = tx18 + tx19 + tx20; 		/* rule 40 */
	 DotNormAbNormCb = tx21*tx22*tx23; 		/* rule 41 */
	 tx24 = -t0; 		/* rule 42 */
	 tx25 = ArcCos(DotNormAbNormCb); 		/* rule 43 */
	 AngleDeviation = tx24 + tx25; 		/* rule 44 */
	IllegalAngle=false;
	if(fabs(DotNormAbNormCb)>(1.0-VERYSMALL)) IllegalAngle=true;
	 tx26 = power2(AngleDeviation); 		/* rule 47 */
	 Energy = kt*tx26; 		/* rule 48 */
	ANGLE_ENERGY_ACCUMULATE(Energy);
	#ifdef ANGLE_CALC_FORCE //[
	if ( calcForce ) {
	 tx27 = reciprocal(tx16); 		/* rule 52 */
	 tx28 = reciprocal(tx17); 		/* rule 53 */
	 tx29 = power2(tx23); 		/* rule 54 */
	 tx452 = tx16; 		/* rule 55 */
	 tx453 = tx21; 		/* rule 56 */
	 tx454 = reciprocal(tx452); 		/* rule 57 */
	 tx30 = tx453*tx454; 		/* rule 58 */
	 tx31 = -(tx27*tx28*tx29); 		/* rule 59 */
	 tx32 = 1 + tx31; 		/* rule 60 */
	 tx33 = -(tx22*tx23*tx30*tx4); 		/* rule 61 */
	 tx34 = tx22*tx453*tx5; 		/* rule 62 */
	 tx455 = mysqrt(tx32); 		/* rule 63 */
	 tx35 = reciprocal(tx455); 		/* rule 64 */
	 tx36 = tx33 + tx34; 		/* rule 65 */
	 gx1 = -2.*AngleDeviation*kt*tx35*tx36; 		/* rule 66 */
	 fx1 = -gx1; 		/* rule 67 */
	ANGLE_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx37 = -(tx22*tx23*tx30*tx6); 		/* rule 69 */
	 tx38 = tx22*tx453*tx7; 		/* rule 70 */
	 tx39 = tx37 + tx38; 		/* rule 71 */
	 gy1 = -2.*AngleDeviation*kt*tx35*tx39; 		/* rule 72 */
	 fy1 = -gy1; 		/* rule 73 */
	ANGLE_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx40 = -(tx22*tx23*tx30*tx8); 		/* rule 75 */
	 tx41 = tx22*tx453*tx9; 		/* rule 76 */
	 tx42 = tx40 + tx41; 		/* rule 77 */
	 gz1 = -2.*AngleDeviation*kt*tx35*tx42; 		/* rule 78 */
	 fz1 = -gz1; 		/* rule 79 */
	ANGLE_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tx43 = -x1; 		/* rule 81 */
	 tx44 = 2.*x2; 		/* rule 82 */
	 tx45 = -x3; 		/* rule 83 */
	 tx456 = tx17; 		/* rule 84 */
	 tx457 = tx22; 		/* rule 85 */
	 tx458 = reciprocal(tx456); 		/* rule 86 */
	 tx46 = tx457*tx458; 		/* rule 87 */
	 tx47 = tx43 + tx44 + tx45; 		/* rule 88 */
	 tx48 = tx23*tx30*tx4*tx457; 		/* rule 89 */
	 tx49 = tx453*tx457*tx47; 		/* rule 90 */
	 tx50 = tx23*tx453*tx46*tx5; 		/* rule 91 */
	 tx51 = tx48 + tx49 + tx50; 		/* rule 92 */
	 gx2 = -2.*AngleDeviation*kt*tx35*tx51; 		/* rule 93 */
	 fx2 = -gx2; 		/* rule 94 */
	ANGLE_FORCE_ACCUMULATE(I2, 0, fx2 );
	 tx52 = -y1; 		/* rule 96 */
	 tx53 = 2.*y2; 		/* rule 97 */
	 tx54 = -y3; 		/* rule 98 */
	 tx55 = tx52 + tx53 + tx54; 		/* rule 99 */
	 tx56 = tx453*tx457*tx55; 		/* rule 100 */
	 tx57 = tx23*tx30*tx457*tx6; 		/* rule 101 */
	 tx58 = tx23*tx453*tx46*tx7; 		/* rule 102 */
	 tx59 = tx56 + tx57 + tx58; 		/* rule 103 */
	 gy2 = -2.*AngleDeviation*kt*tx35*tx59; 		/* rule 104 */
	 fy2 = -gy2; 		/* rule 105 */
	ANGLE_FORCE_ACCUMULATE(I2, 1, fy2 );
	 tx60 = -z1; 		/* rule 107 */
	 tx61 = 2.*z2; 		/* rule 108 */
	 tx62 = -z3; 		/* rule 109 */
	 tx63 = tx60 + tx61 + tx62; 		/* rule 110 */
	 tx64 = tx453*tx457*tx63; 		/* rule 111 */
	 tx65 = tx23*tx30*tx457*tx8; 		/* rule 112 */
	 tx66 = tx23*tx453*tx46*tx9; 		/* rule 113 */
	 tx67 = tx64 + tx65 + tx66; 		/* rule 114 */
	 gz2 = -2.*AngleDeviation*kt*tx35*tx67; 		/* rule 115 */
	 fz2 = -gz2; 		/* rule 116 */
	ANGLE_FORCE_ACCUMULATE(I2, 2, fz2 );
	 tx68 = tx4*tx453*tx457; 		/* rule 118 */
	 tx69 = -tx50; 		/* rule 119 */
	 tx70 = tx68 + tx69; 		/* rule 120 */
	 gx3 = -2.*AngleDeviation*kt*tx35*tx70; 		/* rule 121 */
	 fx3 = -gx3; 		/* rule 122 */
	ANGLE_FORCE_ACCUMULATE(I3, 0, fx3 );
	 tx71 = -tx58; 		/* rule 124 */
	 tx72 = tx453*tx457*tx6; 		/* rule 125 */
	 tx73 = tx71 + tx72; 		/* rule 126 */
	 gy3 = -2.*AngleDeviation*kt*tx35*tx73; 		/* rule 127 */
	 fy3 = -gy3; 		/* rule 128 */
	ANGLE_FORCE_ACCUMULATE(I3, 1, fy3 );
	 tx74 = -tx66; 		/* rule 130 */
	 tx75 = tx453*tx457*tx8; 		/* rule 131 */
	 tx76 = tx74 + tx75; 		/* rule 132 */
	 gz3 = -2.*AngleDeviation*kt*tx35*tx76; 		/* rule 133 */
	 fz3 = -gz3; 		/* rule 134 */
	ANGLE_FORCE_ACCUMULATE(I3, 2, fz3 );
	#ifdef ANGLE_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx459 = power2(tx452); 		/* rule 138 */
	 tx460 = reciprocal(tx459); 		/* rule 139 */
	 tx77 = tx453*tx460; 		/* rule 140 */
	 tx78 = power2(tx27); 		/* rule 141 */
	 tx79 = -2.*tx18*tx30*tx457; 		/* rule 142 */
	 tx80 = -(tx23*tx30*tx457); 		/* rule 143 */
	 tx81 = -2.*tx23*tx27*tx28*tx5; 		/* rule 144 */
	 tx82 = 3*tx10*tx23*tx457*tx77; 		/* rule 145 */
	 tx83 = 2.*tx28*tx29*tx4*tx78; 		/* rule 146 */
	 tx461 = tx32; 		/* rule 147 */
	 tx462 = tx35; 		/* rule 148 */
	 tx463 = reciprocal(tx461); 		/* rule 149 */
	 tx84 = tx462*tx463; 		/* rule 150 */
	 tx85 = tx463; 		/* rule 151 */
	 tx86 = power2(tx36); 		/* rule 152 */
	 tx87 = tx79 + tx80 + tx82; 		/* rule 153 */
	 tx88 = tx81 + tx83; 		/* rule 154 */
	 tx89 = 2.*kt*tx85*tx86; 		/* rule 155 */
	 tx90 = -2.*AngleDeviation*kt*tx462*tx87; 		/* rule 156 */
	 tx91 = AngleDeviation*kt*tx36*tx84*tx88; 		/* rule 157 */
	 dhx1x1 = tx89 + tx90 + tx91; 		/* rule 158 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx92 = -2.*tx19*tx30*tx457; 		/* rule 160 */
	 tx93 = -2.*tx23*tx27*tx28*tx7; 		/* rule 161 */
	 tx94 = 3*tx12*tx23*tx457*tx77; 		/* rule 162 */
	 tx95 = 2.*tx28*tx29*tx6*tx78; 		/* rule 163 */
	 tx96 = power2(tx39); 		/* rule 164 */
	 tx97 = tx80 + tx92 + tx94; 		/* rule 165 */
	 tx98 = tx93 + tx95; 		/* rule 166 */
	 tx99 = 2.*kt*tx85*tx96; 		/* rule 167 */
	 tx100 = -2.*AngleDeviation*kt*tx462*tx97; 		/* rule 168 */
	 tx101 = AngleDeviation*kt*tx39*tx84*tx98; 		/* rule 169 */
	 dhy1y1 = tx100 + tx101 + tx99; 		/* rule 170 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx102 = -2.*tx20*tx30*tx457; 		/* rule 172 */
	 tx103 = 3*tx14*tx23*tx457*tx77; 		/* rule 173 */
	 tx104 = 2.*tx28*tx29*tx78*tx8; 		/* rule 174 */
	 tx105 = -2.*tx23*tx27*tx28*tx9; 		/* rule 175 */
	 tx106 = tx104 + tx105; 		/* rule 176 */
	 tx107 = power2(tx42); 		/* rule 177 */
	 tx108 = tx102 + tx103 + tx80; 		/* rule 178 */
	 tx109 = -2.*AngleDeviation*kt*tx108*tx462; 		/* rule 179 */
	 tx110 = AngleDeviation*kt*tx106*tx42*tx84; 		/* rule 180 */
	 tx111 = 2.*kt*tx107*tx85; 		/* rule 181 */
	 dhz1z1 = tx109 + tx110 + tx111; 		/* rule 182 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 tx464 = power2(tx456); 		/* rule 184 */
	 tx465 = reciprocal(tx464); 		/* rule 185 */
	 tx112 = tx457*tx465; 		/* rule 186 */
	 tx113 = power2(tx28); 		/* rule 187 */
	 tx114 = 2.*tx453*tx457; 		/* rule 188 */
	 tx115 = 3*tx11*tx112*tx23*tx453; 		/* rule 189 */
	 tx116 = -(tx23*tx453*tx46); 		/* rule 190 */
	 tx117 = 2.*tx18*tx23*tx30*tx46; 		/* rule 191 */
	 tx118 = -2.*tx23*tx27*tx28*tx47; 		/* rule 192 */
	 tx119 = 2.*tx30*tx4*tx457*tx47; 		/* rule 193 */
	 tx120 = -2.*tx113*tx27*tx29*tx5; 		/* rule 194 */
	 tx121 = 2.*tx453*tx46*tx47*tx5; 		/* rule 195 */
	 tx122 = -2.*tx28*tx29*tx4*tx78; 		/* rule 196 */
	 tx123 = tx118 + tx120 + tx122; 		/* rule 197 */
	 tx124 = power2(tx51); 		/* rule 198 */
	 tx125 = tx114 + tx115 + tx116 + tx117 + tx119 + tx121 + tx80 + tx82; 		/* rule 199 */
	 tx126 = -2.*AngleDeviation*kt*tx125*tx462; 		/* rule 200 */
	 tx127 = AngleDeviation*kt*tx123*tx51*tx84; 		/* rule 201 */
	 tx128 = 2.*kt*tx124*tx85; 		/* rule 202 */
	 dhx2x2 = tx126 + tx127 + tx128; 		/* rule 203 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 tx129 = 3*tx112*tx13*tx23*tx453; 		/* rule 205 */
	 tx130 = 2.*tx19*tx23*tx30*tx46; 		/* rule 206 */
	 tx131 = -2.*tx23*tx27*tx28*tx55; 		/* rule 207 */
	 tx132 = 2.*tx30*tx457*tx55*tx6; 		/* rule 208 */
	 tx133 = -2.*tx113*tx27*tx29*tx7; 		/* rule 209 */
	 tx134 = 2.*tx453*tx46*tx55*tx7; 		/* rule 210 */
	 tx135 = -2.*tx28*tx29*tx6*tx78; 		/* rule 211 */
	 tx136 = tx131 + tx133 + tx135; 		/* rule 212 */
	 tx137 = power2(tx59); 		/* rule 213 */
	 tx138 = tx114 + tx116 + tx129 + tx130 + tx132 + tx134 + tx80 + tx94; 		/* rule 214 */
	 tx139 = -2.*AngleDeviation*kt*tx138*tx462; 		/* rule 215 */
	 tx140 = AngleDeviation*kt*tx136*tx59*tx84; 		/* rule 216 */
	 tx141 = 2.*kt*tx137*tx85; 		/* rule 217 */
	 dhy2y2 = tx139 + tx140 + tx141; 		/* rule 218 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 tx142 = 3*tx112*tx15*tx23*tx453; 		/* rule 220 */
	 tx143 = 2.*tx20*tx23*tx30*tx46; 		/* rule 221 */
	 tx144 = -2.*tx23*tx27*tx28*tx63; 		/* rule 222 */
	 tx145 = 2.*tx30*tx457*tx63*tx8; 		/* rule 223 */
	 tx146 = -2.*tx28*tx29*tx78*tx8; 		/* rule 224 */
	 tx147 = -2.*tx113*tx27*tx29*tx9; 		/* rule 225 */
	 tx148 = 2.*tx453*tx46*tx63*tx9; 		/* rule 226 */
	 tx149 = tx144 + tx146 + tx147; 		/* rule 227 */
	 tx150 = power2(tx67); 		/* rule 228 */
	 tx151 = tx103 + tx114 + tx116 + tx142 + tx143 + tx145 + tx148 + tx80; 		/* rule 229 */
	 tx152 = -2.*AngleDeviation*kt*tx151*tx462; 		/* rule 230 */
	 tx153 = AngleDeviation*kt*tx149*tx67*tx84; 		/* rule 231 */
	 tx154 = 2.*kt*tx150*tx85; 		/* rule 232 */
	 dhz2z2 = tx152 + tx153 + tx154; 		/* rule 233 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	 tx155 = -2.*tx23*tx27*tx28*tx4; 		/* rule 235 */
	 tx156 = -2.*tx18*tx453*tx46; 		/* rule 236 */
	 tx157 = 2.*tx113*tx27*tx29*tx5; 		/* rule 237 */
	 tx158 = tx115 + tx116 + tx156; 		/* rule 238 */
	 tx159 = tx155 + tx157; 		/* rule 239 */
	 tx160 = power2(tx70); 		/* rule 240 */
	 tx161 = -2.*AngleDeviation*kt*tx158*tx462; 		/* rule 241 */
	 tx162 = AngleDeviation*kt*tx159*tx70*tx84; 		/* rule 242 */
	 tx163 = 2.*kt*tx160*tx85; 		/* rule 243 */
	 dhx3x3 = tx161 + tx162 + tx163; 		/* rule 244 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 0, dhx3x3);
	 tx164 = -2.*tx19*tx453*tx46; 		/* rule 246 */
	 tx165 = -2.*tx23*tx27*tx28*tx6; 		/* rule 247 */
	 tx166 = 2.*tx113*tx27*tx29*tx7; 		/* rule 248 */
	 tx167 = tx116 + tx129 + tx164; 		/* rule 249 */
	 tx168 = tx165 + tx166; 		/* rule 250 */
	 tx169 = power2(tx73); 		/* rule 251 */
	 tx170 = -2.*AngleDeviation*kt*tx167*tx462; 		/* rule 252 */
	 tx171 = AngleDeviation*kt*tx168*tx73*tx84; 		/* rule 253 */
	 tx172 = 2.*kt*tx169*tx85; 		/* rule 254 */
	 dhy3y3 = tx170 + tx171 + tx172; 		/* rule 255 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 1, dhy3y3);
	 tx173 = -2.*tx20*tx453*tx46; 		/* rule 257 */
	 tx174 = -2.*tx23*tx27*tx28*tx8; 		/* rule 258 */
	 tx175 = 2.*tx113*tx27*tx29*tx9; 		/* rule 259 */
	 tx176 = tx116 + tx142 + tx173; 		/* rule 260 */
	 tx177 = tx174 + tx175; 		/* rule 261 */
	 tx178 = power2(tx76); 		/* rule 262 */
	 tx179 = -2.*AngleDeviation*kt*tx176*tx462; 		/* rule 263 */
	 tx180 = AngleDeviation*kt*tx177*tx76*tx84; 		/* rule 264 */
	 tx181 = 2.*kt*tx178*tx85; 		/* rule 265 */
	 dhz3z3 = tx179 + tx180 + tx181; 		/* rule 266 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I3, 2, dhz3z3);
	#ifdef ANGLE_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tx182 = -(tx30*tx457*tx5*tx6); 		/* rule 270 */
	 tx183 = -(tx30*tx4*tx457*tx7); 		/* rule 271 */
	 tx184 = 3*tx23*tx4*tx457*tx6*tx77; 		/* rule 272 */
	 tx185 = tx182 + tx183 + tx184; 		/* rule 273 */
	 tx186 = -2.*AngleDeviation*kt*tx185*tx462; 		/* rule 274 */
	 tx187 = 2.*kt*tx36*tx39*tx85; 		/* rule 275 */
	 tx188 = AngleDeviation*kt*tx36*tx84*tx98; 		/* rule 276 */
	 ohx1y1 = tx186 + tx187 + tx188; 		/* rule 277 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tx189 = -(tx30*tx457*tx5*tx8); 		/* rule 279 */
	 tx190 = 3*tx23*tx4*tx457*tx77*tx8; 		/* rule 280 */
	 tx191 = -(tx30*tx4*tx457*tx9); 		/* rule 281 */
	 tx192 = tx189 + tx190 + tx191; 		/* rule 282 */
	 tx193 = -2.*AngleDeviation*kt*tx192*tx462; 		/* rule 283 */
	 tx194 = AngleDeviation*kt*tx106*tx36*tx84; 		/* rule 284 */
	 tx195 = 2.*kt*tx36*tx42*tx85; 		/* rule 285 */
	 ohx1z1 = tx193 + tx194 + tx195; 		/* rule 286 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx196 = -(tx453*tx457); 		/* rule 288 */
	 tx197 = tx18*tx30*tx457; 		/* rule 289 */
	 tx198 = tx23*tx30*tx457; 		/* rule 290 */
	 tx199 = tx11*tx453*tx46; 		/* rule 291 */
	 tx200 = -(tx18*tx23*tx30*tx46); 		/* rule 292 */
	 tx201 = -(tx30*tx4*tx457*tx47); 		/* rule 293 */
	 tx202 = -3*tx10*tx23*tx457*tx77; 		/* rule 294 */
	 tx203 = tx196 + tx197 + tx198 + tx199 + tx200 + tx201 + tx202; 		/* rule 295 */
	 tx204 = -2.*AngleDeviation*kt*tx203*tx462; 		/* rule 296 */
	 tx205 = AngleDeviation*kt*tx123*tx36*tx84; 		/* rule 297 */
	 tx206 = 2.*kt*tx36*tx51*tx85; 		/* rule 298 */
	 ohx1x2 = tx204 + tx205 + tx206; 		/* rule 299 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx207 = -(tx30*tx4*tx457*tx55); 		/* rule 301 */
	 tx208 = tx30*tx457*tx5*tx6; 		/* rule 302 */
	 tx209 = -(tx23*tx30*tx4*tx46*tx7); 		/* rule 303 */
	 tx210 = tx453*tx46*tx5*tx7; 		/* rule 304 */
	 tx211 = -3*tx23*tx4*tx457*tx6*tx77; 		/* rule 305 */
	 tx212 = tx207 + tx208 + tx209 + tx210 + tx211; 		/* rule 306 */
	 tx213 = -2.*AngleDeviation*kt*tx212*tx462; 		/* rule 307 */
	 tx214 = AngleDeviation*kt*tx136*tx36*tx84; 		/* rule 308 */
	 tx215 = 2.*kt*tx36*tx59*tx85; 		/* rule 309 */
	 ohx1y2 = tx213 + tx214 + tx215; 		/* rule 310 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx216 = -(tx30*tx4*tx457*tx63); 		/* rule 312 */
	 tx217 = tx30*tx457*tx5*tx8; 		/* rule 313 */
	 tx218 = -3*tx23*tx4*tx457*tx77*tx8; 		/* rule 314 */
	 tx219 = -(tx23*tx30*tx4*tx46*tx9); 		/* rule 315 */
	 tx220 = tx453*tx46*tx5*tx9; 		/* rule 316 */
	 tx221 = tx216 + tx217 + tx218 + tx219 + tx220; 		/* rule 317 */
	 tx222 = -2.*AngleDeviation*kt*tx221*tx462; 		/* rule 318 */
	 tx223 = AngleDeviation*kt*tx149*tx36*tx84; 		/* rule 319 */
	 tx224 = 2.*kt*tx36*tx67*tx85; 		/* rule 320 */
	 ohx1z2 = tx222 + tx223 + tx224; 		/* rule 321 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tx225 = -tx199; 		/* rule 323 */
	 tx226 = tx453*tx457; 		/* rule 324 */
	 tx227 = -(tx10*tx30*tx457); 		/* rule 325 */
	 tx228 = tx18*tx23*tx30*tx46; 		/* rule 326 */
	 tx229 = tx225 + tx226 + tx227 + tx228; 		/* rule 327 */
	 tx230 = -2.*AngleDeviation*kt*tx229*tx462; 		/* rule 328 */
	 tx231 = AngleDeviation*kt*tx159*tx36*tx84; 		/* rule 329 */
	 tx232 = 2.*kt*tx36*tx70*tx85; 		/* rule 330 */
	 ohx1x3 = tx230 + tx231 + tx232; 		/* rule 331 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 0, ohx1x3);
	 tx233 = -tx210; 		/* rule 333 */
	 tx234 = -(tx30*tx4*tx457*tx6); 		/* rule 334 */
	 tx235 = tx23*tx30*tx4*tx46*tx7; 		/* rule 335 */
	 tx236 = tx233 + tx234 + tx235; 		/* rule 336 */
	 tx237 = -2.*AngleDeviation*kt*tx236*tx462; 		/* rule 337 */
	 tx238 = AngleDeviation*kt*tx168*tx36*tx84; 		/* rule 338 */
	 tx239 = 2.*kt*tx36*tx73*tx85; 		/* rule 339 */
	 ohx1y3 = tx237 + tx238 + tx239; 		/* rule 340 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 1, ohx1y3);
	 tx240 = -tx220; 		/* rule 342 */
	 tx241 = -(tx30*tx4*tx457*tx8); 		/* rule 343 */
	 tx242 = tx23*tx30*tx4*tx46*tx9; 		/* rule 344 */
	 tx243 = tx240 + tx241 + tx242; 		/* rule 345 */
	 tx244 = -2.*AngleDeviation*kt*tx243*tx462; 		/* rule 346 */
	 tx245 = AngleDeviation*kt*tx177*tx36*tx84; 		/* rule 347 */
	 tx246 = 2.*kt*tx36*tx76*tx85; 		/* rule 348 */
	 ohx1z3 = tx244 + tx245 + tx246; 		/* rule 349 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 2, ohx1z3);
	 tx247 = -(tx30*tx457*tx7*tx8); 		/* rule 351 */
	 tx248 = 3*tx23*tx457*tx6*tx77*tx8; 		/* rule 352 */
	 tx249 = -(tx30*tx457*tx6*tx9); 		/* rule 353 */
	 tx250 = tx247 + tx248 + tx249; 		/* rule 354 */
	 tx251 = -2.*AngleDeviation*kt*tx250*tx462; 		/* rule 355 */
	 tx252 = AngleDeviation*kt*tx106*tx39*tx84; 		/* rule 356 */
	 tx253 = 2.*kt*tx39*tx42*tx85; 		/* rule 357 */
	 ohy1z1 = tx251 + tx252 + tx253; 		/* rule 358 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 tx254 = -(tx30*tx457*tx47*tx6); 		/* rule 360 */
	 tx255 = -(tx23*tx30*tx46*tx5*tx6); 		/* rule 361 */
	 tx256 = tx30*tx4*tx457*tx7; 		/* rule 362 */
	 tx257 = tx210 + tx211 + tx254 + tx255 + tx256; 		/* rule 363 */
	 tx258 = -2.*AngleDeviation*kt*tx257*tx462; 		/* rule 364 */
	 tx259 = AngleDeviation*kt*tx123*tx39*tx84; 		/* rule 365 */
	 tx260 = 2.*kt*tx39*tx51*tx85; 		/* rule 366 */
	 ohy1x2 = tx258 + tx259 + tx260; 		/* rule 367 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx261 = tx19*tx30*tx457; 		/* rule 369 */
	 tx262 = tx13*tx453*tx46; 		/* rule 370 */
	 tx263 = -(tx19*tx23*tx30*tx46); 		/* rule 371 */
	 tx264 = -(tx30*tx457*tx55*tx6); 		/* rule 372 */
	 tx265 = -3*tx12*tx23*tx457*tx77; 		/* rule 373 */
	 tx266 = tx196 + tx198 + tx261 + tx262 + tx263 + tx264 + tx265; 		/* rule 374 */
	 tx267 = -2.*AngleDeviation*kt*tx266*tx462; 		/* rule 375 */
	 tx268 = AngleDeviation*kt*tx136*tx39*tx84; 		/* rule 376 */
	 tx269 = 2.*kt*tx39*tx59*tx85; 		/* rule 377 */
	 ohy1y2 = tx267 + tx268 + tx269; 		/* rule 378 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx270 = -(tx30*tx457*tx6*tx63); 		/* rule 380 */
	 tx271 = tx30*tx457*tx7*tx8; 		/* rule 381 */
	 tx272 = -3*tx23*tx457*tx6*tx77*tx8; 		/* rule 382 */
	 tx273 = -(tx23*tx30*tx46*tx6*tx9); 		/* rule 383 */
	 tx274 = tx453*tx46*tx7*tx9; 		/* rule 384 */
	 tx275 = tx270 + tx271 + tx272 + tx273 + tx274; 		/* rule 385 */
	 tx276 = -2.*AngleDeviation*kt*tx275*tx462; 		/* rule 386 */
	 tx277 = AngleDeviation*kt*tx149*tx39*tx84; 		/* rule 387 */
	 tx278 = 2.*kt*tx39*tx67*tx85; 		/* rule 388 */
	 ohy1z2 = tx276 + tx277 + tx278; 		/* rule 389 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 tx279 = tx23*tx30*tx46*tx5*tx6; 		/* rule 391 */
	 tx280 = tx233 + tx234 + tx279; 		/* rule 392 */
	 tx281 = -2.*AngleDeviation*kt*tx280*tx462; 		/* rule 393 */
	 tx282 = AngleDeviation*kt*tx159*tx39*tx84; 		/* rule 394 */
	 tx283 = 2.*kt*tx39*tx70*tx85; 		/* rule 395 */
	 ohy1x3 = tx281 + tx282 + tx283; 		/* rule 396 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 0, ohy1x3);
	 tx284 = -tx262; 		/* rule 398 */
	 tx285 = -(tx12*tx30*tx457); 		/* rule 399 */
	 tx286 = tx19*tx23*tx30*tx46; 		/* rule 400 */
	 tx287 = tx226 + tx284 + tx285 + tx286; 		/* rule 401 */
	 tx288 = -2.*AngleDeviation*kt*tx287*tx462; 		/* rule 402 */
	 tx289 = AngleDeviation*kt*tx168*tx39*tx84; 		/* rule 403 */
	 tx290 = 2.*kt*tx39*tx73*tx85; 		/* rule 404 */
	 ohy1y3 = tx288 + tx289 + tx290; 		/* rule 405 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 1, ohy1y3);
	 tx291 = -tx274; 		/* rule 407 */
	 tx292 = -(tx30*tx457*tx6*tx8); 		/* rule 408 */
	 tx293 = tx23*tx30*tx46*tx6*tx9; 		/* rule 409 */
	 tx294 = tx291 + tx292 + tx293; 		/* rule 410 */
	 tx295 = -2.*AngleDeviation*kt*tx294*tx462; 		/* rule 411 */
	 tx296 = AngleDeviation*kt*tx177*tx39*tx84; 		/* rule 412 */
	 tx297 = 2.*kt*tx39*tx76*tx85; 		/* rule 413 */
	 ohy1z3 = tx295 + tx296 + tx297; 		/* rule 414 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 2, ohy1z3);
	 tx298 = -(tx30*tx457*tx47*tx8); 		/* rule 416 */
	 tx299 = -(tx23*tx30*tx46*tx5*tx8); 		/* rule 417 */
	 tx300 = tx30*tx4*tx457*tx9; 		/* rule 418 */
	 tx301 = tx218 + tx220 + tx298 + tx299 + tx300; 		/* rule 419 */
	 tx302 = -2.*AngleDeviation*kt*tx301*tx462; 		/* rule 420 */
	 tx303 = AngleDeviation*kt*tx123*tx42*tx84; 		/* rule 421 */
	 tx304 = 2.*kt*tx42*tx51*tx85; 		/* rule 422 */
	 ohz1x2 = tx302 + tx303 + tx304; 		/* rule 423 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 tx305 = -(tx30*tx457*tx55*tx8); 		/* rule 425 */
	 tx306 = -(tx23*tx30*tx46*tx7*tx8); 		/* rule 426 */
	 tx307 = tx30*tx457*tx6*tx9; 		/* rule 427 */
	 tx308 = tx272 + tx274 + tx305 + tx306 + tx307; 		/* rule 428 */
	 tx309 = -2.*AngleDeviation*kt*tx308*tx462; 		/* rule 429 */
	 tx310 = AngleDeviation*kt*tx136*tx42*tx84; 		/* rule 430 */
	 tx311 = 2.*kt*tx42*tx59*tx85; 		/* rule 431 */
	 ohz1y2 = tx309 + tx310 + tx311; 		/* rule 432 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx312 = tx20*tx30*tx457; 		/* rule 434 */
	 tx313 = tx15*tx453*tx46; 		/* rule 435 */
	 tx314 = -(tx20*tx23*tx30*tx46); 		/* rule 436 */
	 tx315 = -3*tx14*tx23*tx457*tx77; 		/* rule 437 */
	 tx316 = -(tx30*tx457*tx63*tx8); 		/* rule 438 */
	 tx317 = tx196 + tx198 + tx312 + tx313 + tx314 + tx315 + tx316; 		/* rule 439 */
	 tx318 = -2.*AngleDeviation*kt*tx317*tx462; 		/* rule 440 */
	 tx319 = AngleDeviation*kt*tx149*tx42*tx84; 		/* rule 441 */
	 tx320 = 2.*kt*tx42*tx67*tx85; 		/* rule 442 */
	 ohz1z2 = tx318 + tx319 + tx320; 		/* rule 443 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 tx321 = tx23*tx30*tx46*tx5*tx8; 		/* rule 445 */
	 tx322 = tx240 + tx241 + tx321; 		/* rule 446 */
	 tx323 = -2.*AngleDeviation*kt*tx322*tx462; 		/* rule 447 */
	 tx324 = AngleDeviation*kt*tx159*tx42*tx84; 		/* rule 448 */
	 tx325 = 2.*kt*tx42*tx70*tx85; 		/* rule 449 */
	 ohz1x3 = tx323 + tx324 + tx325; 		/* rule 450 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 0, ohz1x3);
	 tx326 = tx23*tx30*tx46*tx7*tx8; 		/* rule 452 */
	 tx327 = tx291 + tx292 + tx326; 		/* rule 453 */
	 tx328 = -2.*AngleDeviation*kt*tx327*tx462; 		/* rule 454 */
	 tx329 = AngleDeviation*kt*tx168*tx42*tx84; 		/* rule 455 */
	 tx330 = 2.*kt*tx42*tx73*tx85; 		/* rule 456 */
	 ohz1y3 = tx328 + tx329 + tx330; 		/* rule 457 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 1, ohz1y3);
	 tx331 = -(tx14*tx30*tx457); 		/* rule 459 */
	 tx332 = -tx313; 		/* rule 460 */
	 tx333 = tx20*tx23*tx30*tx46; 		/* rule 461 */
	 tx334 = tx226 + tx331 + tx332 + tx333; 		/* rule 462 */
	 tx335 = -2.*AngleDeviation*kt*tx334*tx462; 		/* rule 463 */
	 tx336 = AngleDeviation*kt*tx177*tx42*tx84; 		/* rule 464 */
	 tx337 = 2.*kt*tx42*tx76*tx85; 		/* rule 465 */
	 ohz1z3 = tx335 + tx336 + tx337; 		/* rule 466 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 2, ohz1z3);
	 tx338 = tx30*tx4*tx457*tx55; 		/* rule 468 */
	 tx339 = tx453*tx46*tx5*tx55; 		/* rule 469 */
	 tx340 = tx30*tx457*tx47*tx6; 		/* rule 470 */
	 tx341 = tx453*tx46*tx47*tx7; 		/* rule 471 */
	 tx342 = 3*tx112*tx23*tx453*tx5*tx7; 		/* rule 472 */
	 tx343 = tx184 + tx235 + tx279 + tx338 + tx339 + tx340 + tx341 + tx342; 		/* rule 473 */
	 tx344 = -2.*AngleDeviation*kt*tx343*tx462; 		/* rule 474 */
	 tx345 = AngleDeviation*kt*tx136*tx51*tx84; 		/* rule 475 */
	 tx346 = 2.*kt*tx51*tx59*tx85; 		/* rule 476 */
	 ohx2y2 = tx344 + tx345 + tx346; 		/* rule 477 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 tx347 = tx30*tx4*tx457*tx63; 		/* rule 479 */
	 tx348 = tx453*tx46*tx5*tx63; 		/* rule 480 */
	 tx349 = tx30*tx457*tx47*tx8; 		/* rule 481 */
	 tx350 = tx453*tx46*tx47*tx9; 		/* rule 482 */
	 tx351 = 3*tx112*tx23*tx453*tx5*tx9; 		/* rule 483 */
	 tx352 = tx190 + tx242 + tx321 + tx347 + tx348 + tx349 + tx350 + tx351; 		/* rule 484 */
	 tx353 = -2.*AngleDeviation*kt*tx352*tx462; 		/* rule 485 */
	 tx354 = AngleDeviation*kt*tx149*tx51*tx84; 		/* rule 486 */
	 tx355 = 2.*kt*tx51*tx67*tx85; 		/* rule 487 */
	 ohx2z2 = tx353 + tx354 + tx355; 		/* rule 488 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 tx356 = -3*tx11*tx112*tx23*tx453; 		/* rule 490 */
	 tx357 = tx10*tx30*tx457; 		/* rule 491 */
	 tx358 = tx18*tx453*tx46; 		/* rule 492 */
	 tx359 = tx23*tx453*tx46; 		/* rule 493 */
	 tx360 = -(tx453*tx46*tx47*tx5); 		/* rule 494 */
	 tx361 = tx196 + tx200 + tx356 + tx357 + tx358 + tx359 + tx360; 		/* rule 495 */
	 tx362 = -2.*AngleDeviation*kt*tx361*tx462; 		/* rule 496 */
	 tx363 = AngleDeviation*kt*tx159*tx51*tx84; 		/* rule 497 */
	 tx364 = 2.*kt*tx51*tx70*tx85; 		/* rule 498 */
	 ohx2x3 = tx362 + tx363 + tx364; 		/* rule 499 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 0, ohx2x3);
	 tx365 = -tx341; 		/* rule 501 */
	 tx366 = tx30*tx4*tx457*tx6; 		/* rule 502 */
	 tx367 = tx453*tx46*tx5*tx6; 		/* rule 503 */
	 tx368 = -3*tx112*tx23*tx453*tx5*tx7; 		/* rule 504 */
	 tx369 = tx209 + tx365 + tx366 + tx367 + tx368; 		/* rule 505 */
	 tx370 = -2.*AngleDeviation*kt*tx369*tx462; 		/* rule 506 */
	 tx371 = AngleDeviation*kt*tx168*tx51*tx84; 		/* rule 507 */
	 tx372 = 2.*kt*tx51*tx73*tx85; 		/* rule 508 */
	 ohx2y3 = tx370 + tx371 + tx372; 		/* rule 509 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 1, ohx2y3);
	 tx373 = -tx350; 		/* rule 511 */
	 tx374 = tx30*tx4*tx457*tx8; 		/* rule 512 */
	 tx375 = tx453*tx46*tx5*tx8; 		/* rule 513 */
	 tx376 = -3*tx112*tx23*tx453*tx5*tx9; 		/* rule 514 */
	 tx377 = tx219 + tx373 + tx374 + tx375 + tx376; 		/* rule 515 */
	 tx378 = -2.*AngleDeviation*kt*tx377*tx462; 		/* rule 516 */
	 tx379 = AngleDeviation*kt*tx177*tx51*tx84; 		/* rule 517 */
	 tx380 = 2.*kt*tx51*tx76*tx85; 		/* rule 518 */
	 ohx2z3 = tx378 + tx379 + tx380; 		/* rule 519 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 2, ohx2z3);
	 tx381 = tx30*tx457*tx6*tx63; 		/* rule 521 */
	 tx382 = tx453*tx46*tx63*tx7; 		/* rule 522 */
	 tx383 = tx30*tx457*tx55*tx8; 		/* rule 523 */
	 tx384 = tx453*tx46*tx55*tx9; 		/* rule 524 */
	 tx385 = 3*tx112*tx23*tx453*tx7*tx9; 		/* rule 525 */
	 tx386 = tx248 + tx293 + tx326 + tx381 + tx382 + tx383 + tx384 + tx385; 		/* rule 526 */
	 tx387 = -2.*AngleDeviation*kt*tx386*tx462; 		/* rule 527 */
	 tx388 = AngleDeviation*kt*tx149*tx59*tx84; 		/* rule 528 */
	 tx389 = 2.*kt*tx59*tx67*tx85; 		/* rule 529 */
	 ohy2z2 = tx387 + tx388 + tx389; 		/* rule 530 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	 tx390 = -tx339; 		/* rule 532 */
	 tx391 = tx4*tx453*tx46*tx7; 		/* rule 533 */
	 tx392 = tx255 + tx366 + tx368 + tx390 + tx391; 		/* rule 534 */
	 tx393 = -2.*AngleDeviation*kt*tx392*tx462; 		/* rule 535 */
	 tx394 = AngleDeviation*kt*tx159*tx59*tx84; 		/* rule 536 */
	 tx395 = 2.*kt*tx59*tx70*tx85; 		/* rule 537 */
	 ohy2x3 = tx393 + tx394 + tx395; 		/* rule 538 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 0, ohy2x3);
	 tx396 = -3*tx112*tx13*tx23*tx453; 		/* rule 540 */
	 tx397 = tx12*tx30*tx457; 		/* rule 541 */
	 tx398 = tx19*tx453*tx46; 		/* rule 542 */
	 tx399 = -(tx453*tx46*tx55*tx7); 		/* rule 543 */
	 tx400 = tx196 + tx263 + tx359 + tx396 + tx397 + tx398 + tx399; 		/* rule 544 */
	 tx401 = -2.*AngleDeviation*kt*tx400*tx462; 		/* rule 545 */
	 tx402 = AngleDeviation*kt*tx168*tx59*tx84; 		/* rule 546 */
	 tx403 = 2.*kt*tx59*tx73*tx85; 		/* rule 547 */
	 ohy2y3 = tx401 + tx402 + tx403; 		/* rule 548 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 1, ohy2y3);
	 tx404 = -tx384; 		/* rule 550 */
	 tx405 = tx30*tx457*tx6*tx8; 		/* rule 551 */
	 tx406 = tx453*tx46*tx7*tx8; 		/* rule 552 */
	 tx407 = -3*tx112*tx23*tx453*tx7*tx9; 		/* rule 553 */
	 tx408 = tx273 + tx404 + tx405 + tx406 + tx407; 		/* rule 554 */
	 tx409 = -2.*AngleDeviation*kt*tx408*tx462; 		/* rule 555 */
	 tx410 = AngleDeviation*kt*tx177*tx59*tx84; 		/* rule 556 */
	 tx411 = 2.*kt*tx59*tx76*tx85; 		/* rule 557 */
	 ohy2z3 = tx409 + tx410 + tx411; 		/* rule 558 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 2, ohy2z3);
	 tx412 = -tx348; 		/* rule 560 */
	 tx413 = tx4*tx453*tx46*tx9; 		/* rule 561 */
	 tx414 = tx299 + tx374 + tx376 + tx412 + tx413; 		/* rule 562 */
	 tx415 = -2.*AngleDeviation*kt*tx414*tx462; 		/* rule 563 */
	 tx416 = AngleDeviation*kt*tx159*tx67*tx84; 		/* rule 564 */
	 tx417 = 2.*kt*tx67*tx70*tx85; 		/* rule 565 */
	 ohz2x3 = tx415 + tx416 + tx417; 		/* rule 566 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 0, ohz2x3);
	 tx418 = -tx382; 		/* rule 568 */
	 tx419 = tx453*tx46*tx6*tx9; 		/* rule 569 */
	 tx420 = tx306 + tx405 + tx407 + tx418 + tx419; 		/* rule 570 */
	 tx421 = -2.*AngleDeviation*kt*tx420*tx462; 		/* rule 571 */
	 tx422 = AngleDeviation*kt*tx168*tx67*tx84; 		/* rule 572 */
	 tx423 = 2.*kt*tx67*tx73*tx85; 		/* rule 573 */
	 ohz2y3 = tx421 + tx422 + tx423; 		/* rule 574 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 1, ohz2y3);
	 tx424 = -3*tx112*tx15*tx23*tx453; 		/* rule 576 */
	 tx425 = tx14*tx30*tx457; 		/* rule 577 */
	 tx426 = tx20*tx453*tx46; 		/* rule 578 */
	 tx427 = -(tx453*tx46*tx63*tx9); 		/* rule 579 */
	 tx428 = tx196 + tx314 + tx359 + tx424 + tx425 + tx426 + tx427; 		/* rule 580 */
	 tx429 = -2.*AngleDeviation*kt*tx428*tx462; 		/* rule 581 */
	 tx430 = AngleDeviation*kt*tx177*tx67*tx84; 		/* rule 582 */
	 tx431 = 2.*kt*tx67*tx76*tx85; 		/* rule 583 */
	 ohz2z3 = tx429 + tx430 + tx431; 		/* rule 584 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 2, ohz2z3);
	 tx432 = -tx367; 		/* rule 586 */
	 tx433 = -tx391; 		/* rule 587 */
	 tx434 = tx342 + tx432 + tx433; 		/* rule 588 */
	 tx435 = -2.*AngleDeviation*kt*tx434*tx462; 		/* rule 589 */
	 tx436 = AngleDeviation*kt*tx168*tx70*tx84; 		/* rule 590 */
	 tx437 = 2.*kt*tx70*tx73*tx85; 		/* rule 591 */
	 ohx3y3 = tx435 + tx436 + tx437; 		/* rule 592 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 1, ohx3y3);
	 tx438 = -tx375; 		/* rule 594 */
	 tx439 = -tx413; 		/* rule 595 */
	 tx440 = tx351 + tx438 + tx439; 		/* rule 596 */
	 tx441 = -2.*AngleDeviation*kt*tx440*tx462; 		/* rule 597 */
	 tx442 = AngleDeviation*kt*tx177*tx70*tx84; 		/* rule 598 */
	 tx443 = 2.*kt*tx70*tx76*tx85; 		/* rule 599 */
	 ohx3z3 = tx441 + tx442 + tx443; 		/* rule 600 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 2, ohx3z3);
	 tx444 = -tx406; 		/* rule 602 */
	 tx445 = -tx419; 		/* rule 603 */
	 tx446 = tx385 + tx444 + tx445; 		/* rule 604 */
	 tx447 = -2.*AngleDeviation*kt*tx446*tx462; 		/* rule 605 */
	 tx448 = AngleDeviation*kt*tx177*tx73*tx84; 		/* rule 606 */
	 tx449 = 2.*kt*tx73*tx76*tx85; 		/* rule 607 */
	 ohy3z3 = tx447 + tx448 + tx449; 		/* rule 608 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 2, ohy3z3);
	} /*if calcOffDiagonalHessian */ 
	#endif /* ANGLE_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* ANGLE_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* ANGLE_CALC_FORCE ]*/
	SKIP_term:  (void)0;
