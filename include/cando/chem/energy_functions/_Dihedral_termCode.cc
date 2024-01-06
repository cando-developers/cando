// C-code
	DIHEDRAL_SET_PARAMETER(sinPhase);
	DIHEDRAL_SET_PARAMETER(cosPhase);
	DIHEDRAL_SET_PARAMETER(V);
	DIHEDRAL_SET_PARAMETER(DN);
	DIHEDRAL_SET_PARAMETER(IN);
	DIHEDRAL_SET_PARAMETER(I1);
	DIHEDRAL_SET_PARAMETER(I2);
	DIHEDRAL_SET_PARAMETER(I3);
	DIHEDRAL_SET_PARAMETER(I4);
	DIHEDRAL_APPLY_ATOM_MASK(I1,I2,I3,I4);
	DIHEDRAL_SET_POSITION(x1,I1,0);
	DIHEDRAL_SET_POSITION(y1,I1,1);
	DIHEDRAL_SET_POSITION(z1,I1,2);
	DIHEDRAL_SET_POSITION(x2,I2,0);
	DIHEDRAL_SET_POSITION(y2,I2,1);
	DIHEDRAL_SET_POSITION(z2,I2,2);
	DIHEDRAL_SET_POSITION(x3,I3,0);
	DIHEDRAL_SET_POSITION(y3,I3,1);
	DIHEDRAL_SET_POSITION(z3,I3,2);
	DIHEDRAL_SET_POSITION(x4,I4,0);
	DIHEDRAL_SET_POSITION(y4,I4,1);
	DIHEDRAL_SET_POSITION(z4,I4,2);
	 tx1 = -(x2*y1); 		/* rule 23 */
	 tx2 = x3*y1; 		/* rule 24 */
	 tx3 = x1*y2; 		/* rule 25 */
	 tx4 = -(x3*y2); 		/* rule 26 */
	 tx5 = -(x1*y3); 		/* rule 27 */
	 tx6 = x2*y3; 		/* rule 28 */
	 tx7 = x2*z1; 		/* rule 29 */
	 tx8 = -(x3*z1); 		/* rule 30 */
	 tx9 = -(y2*z1); 		/* rule 31 */
	 tx10 = y3*z1; 		/* rule 32 */
	 tx11 = -(x1*z2); 		/* rule 33 */
	 tx12 = x3*z2; 		/* rule 34 */
	 tx13 = y1*z2; 		/* rule 35 */
	 tx14 = -(y3*z2); 		/* rule 36 */
	 tx15 = x1*z3; 		/* rule 37 */
	 tx16 = -(x2*z3); 		/* rule 38 */
	 tx17 = -(y1*z3); 		/* rule 39 */
	 tx18 = y2*z3; 		/* rule 40 */
	 tx19 = tx1 + tx2 + tx3 + tx4 + tx5 + tx6; 		/* rule 41 */
	 tx20 = tx11 + tx12 + tx15 + tx16 + tx7 + tx8; 		/* rule 42 */
	 tx21 = tx10 + tx13 + tx14 + tx17 + tx18 + tx9; 		/* rule 43 */
	 tx22 = power2(tx19); 		/* rule 44 */
	 tx23 = power2(tx20); 		/* rule 45 */
	 tx24 = power2(tx21); 		/* rule 46 */
	 tx25 = tx22 + tx23 + tx24; 		/* rule 47 */
	 LenA = mysqrt(tx25); 		/* rule 48 */
	 tx26 = x4*y2; 		/* rule 49 */
	 tx27 = -(x4*y3); 		/* rule 50 */
	 tx28 = -(x2*y4); 		/* rule 51 */
	 tx29 = x3*y4; 		/* rule 52 */
	 tx30 = -(x4*z2); 		/* rule 53 */
	 tx31 = y4*z2; 		/* rule 54 */
	 tx32 = x4*z3; 		/* rule 55 */
	 tx33 = -(y4*z3); 		/* rule 56 */
	 tx34 = x2*z4; 		/* rule 57 */
	 tx35 = -(x3*z4); 		/* rule 58 */
	 tx36 = -(y2*z4); 		/* rule 59 */
	 tx37 = y3*z4; 		/* rule 60 */
	 tx38 = tx12 + tx16 + tx30 + tx32 + tx34 + tx35; 		/* rule 61 */
	 tx39 = tx14 + tx18 + tx31 + tx33 + tx36 + tx37; 		/* rule 62 */
	 tx40 = tx26 + tx27 + tx28 + tx29 + tx4 + tx6; 		/* rule 63 */
	 tx41 = power2(tx38); 		/* rule 64 */
	 tx42 = power2(tx39); 		/* rule 65 */
	 tx43 = power2(tx40); 		/* rule 66 */
	 tx44 = tx41 + tx42 + tx43; 		/* rule 67 */
	 LenB = mysqrt(tx44); 		/* rule 68 */
	 ReciprocalLenA = reciprocal(LenA); 		/* rule 69 */
	 ReciprocalLenB = reciprocal(LenB); 		/* rule 70 */
	ZERO_SMALL_LEN(LenA,ReciprocalLenA);
	ZERO_SMALL_LEN(LenB,ReciprocalLenB);
	 RecLenARecLenB = ReciprocalLenA*ReciprocalLenB; 		/* rule 73 */
EraseLinearDihedral = VEC_CONST(1.0);
	EraseLinearDihedral = (RecLenARecLenB==0.0) ? 0.0 : EraseLinearDihedral;
	 tx45 = tx20*tx38; 		/* rule 76 */
	 tx46 = tx21*tx39; 		/* rule 77 */
	 tx47 = tx19*tx40; 		/* rule 78 */
	 tx48 = tx45 + tx46 + tx47; 		/* rule 79 */
	 CosPhi = RecLenARecLenB*tx48; 		/* rule 80 */
	 tx49 = -x3; 		/* rule 81 */
	 tx50 = -y3; 		/* rule 82 */
	 tx51 = -z3; 		/* rule 83 */
	 tx52 = tx49 + x2; 		/* rule 84 */
	 tx53 = tx50 + y2; 		/* rule 85 */
	 tx54 = tx51 + z2; 		/* rule 86 */
	 tx55 = power2(tx52); 		/* rule 87 */
	 tx56 = power2(tx53); 		/* rule 88 */
	 tx57 = power2(tx54); 		/* rule 89 */
	 tx58 = tx49 + x4; 		/* rule 90 */
	 tx59 = tx50 + y4; 		/* rule 91 */
	 tx60 = tx51 + z4; 		/* rule 92 */
	 tx61 = tx55 + tx56 + tx57; 		/* rule 93 */
	 tx62 = tx21*tx58; 		/* rule 94 */
	 tx63 = tx20*tx59; 		/* rule 95 */
	 tx64 = tx19*tx60; 		/* rule 96 */
	 tx65 = mysqrt(tx61); 		/* rule 97 */
	 tx66 = tx62 + tx63 + tx64; 		/* rule 98 */
	 SinPhi = RecLenARecLenB*tx65*tx66; 		/* rule 99 */
	CosPhi=MAX(-1.0,MIN(1.0,CosPhi));
/*CosNPhi = mathCosNPhi[IN,SinPhi,CosPhi];*/
/*SinNPhi = mathSinNPhi[IN,SinPhi,CosPhi];*/
	DO_sinNPhiCosNPhi(IN, SinNPhi, CosNPhi, SinPhi, CosPhi);
	 tx67 = CosNPhi*cosPhase; 		/* rule 104 */
	 tx68 = SinNPhi*sinPhase; 		/* rule 105 */
	 DihedralDeviation = 1. + tx67 + tx68; 		/* rule 106 */
	 Energy = DihedralDeviation*EraseLinearDihedral*V; 		/* rule 107 */
	DIHEDRAL_ENERGY_ACCUMULATE(Energy);
	#ifdef DIHEDRAL_CALC_FORCE //[
	if (calcForce ) {
	 tx69 = cosPhase*SinNPhi; 		/* rule 111 */
	 tx70 = -(CosNPhi*sinPhase); 		/* rule 112 */
	 tx71 = tx69 + tx70; 		/* rule 113 */
	 DeDPhi = -(DN*EraseLinearDihedral*tx71*V); 		/* rule 114 */
	 tx72 = reciprocal(tx25); 		/* rule 115 */
	 gx1 = -(DeDPhi*tx21*tx65*tx72); 		/* rule 116 */
	 fx1 = -gx1; 		/* rule 117 */
	DIHEDRAL_FORCE_ACCUMULATE(I1, 0, fx1 );
	 gy1 = -(DeDPhi*tx20*tx65*tx72); 		/* rule 119 */
	 fy1 = -gy1; 		/* rule 120 */
	DIHEDRAL_FORCE_ACCUMULATE(I1, 1, fy1 );
	 gz1 = -(DeDPhi*tx19*tx65*tx72); 		/* rule 122 */
	 fz1 = -gz1; 		/* rule 123 */
	DIHEDRAL_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tx73 = -x2; 		/* rule 125 */
	 tx74 = -y2; 		/* rule 126 */
	 tx75 = -z2; 		/* rule 127 */
	 tx76 = tx73 + x1; 		/* rule 128 */
	 tx77 = tx74 + y1; 		/* rule 129 */
	 tx78 = tx75 + z1; 		/* rule 130 */
	 tx79 = tx52*tx58; 		/* rule 131 */
	 tx80 = tx53*tx59; 		/* rule 132 */
	 tx81 = tx54*tx60; 		/* rule 133 */
	 tx82 = tx52*tx76; 		/* rule 134 */
	 tx83 = tx53*tx77; 		/* rule 135 */
	 tx84 = tx54*tx78; 		/* rule 136 */
	 tx85 = reciprocal(tx44); 		/* rule 137 */
	 tx86 = reciprocal(tx65); 		/* rule 138 */
	 tx87 = tx79 + tx80 + tx81; 		/* rule 139 */
	 tx88 = tx82 + tx83 + tx84; 		/* rule 140 */
	 tx89 = tx21*tx65*tx72; 		/* rule 141 */
	 tx90 = -(tx39*tx85*tx86*tx87); 		/* rule 142 */
	 tx91 = tx21*tx72*tx86*tx88; 		/* rule 143 */
	 tx92 = tx89 + tx90 + tx91; 		/* rule 144 */
	 gx2 = DeDPhi*tx92; 		/* rule 145 */
	 fx2 = -gx2; 		/* rule 146 */
	DIHEDRAL_FORCE_ACCUMULATE(I2, 0, fx2 );
	 tx93 = tx20*tx65*tx72; 		/* rule 148 */
	 tx94 = -(tx38*tx85*tx86*tx87); 		/* rule 149 */
	 tx95 = tx20*tx72*tx86*tx88; 		/* rule 150 */
	 tx96 = tx93 + tx94 + tx95; 		/* rule 151 */
	 gy2 = DeDPhi*tx96; 		/* rule 152 */
	 fy2 = -gy2; 		/* rule 153 */
	DIHEDRAL_FORCE_ACCUMULATE(I2, 1, fy2 );
	 tx97 = tx19*tx65*tx72; 		/* rule 155 */
	 tx98 = -(tx40*tx85*tx86*tx87); 		/* rule 156 */
	 tx99 = tx19*tx72*tx86*tx88; 		/* rule 157 */
	 tx100 = tx97 + tx98 + tx99; 		/* rule 158 */
	 gz2 = DeDPhi*tx100; 		/* rule 159 */
	 fz2 = -gz2; 		/* rule 160 */
	DIHEDRAL_FORCE_ACCUMULATE(I2, 2, fz2 );
	 tx101 = -(tx39*tx65*tx85); 		/* rule 162 */
	 tx102 = tx39*tx85*tx86*tx87; 		/* rule 163 */
	 tx103 = -tx91; 		/* rule 164 */
	 tx104 = tx101 + tx102 + tx103; 		/* rule 165 */
	 gx3 = DeDPhi*tx104; 		/* rule 166 */
	 fx3 = -gx3; 		/* rule 167 */
	DIHEDRAL_FORCE_ACCUMULATE(I3, 0, fx3 );
	 tx105 = -(tx38*tx65*tx85); 		/* rule 169 */
	 tx106 = tx38*tx85*tx86*tx87; 		/* rule 170 */
	 tx107 = -tx95; 		/* rule 171 */
	 tx108 = tx105 + tx106 + tx107; 		/* rule 172 */
	 gy3 = DeDPhi*tx108; 		/* rule 173 */
	 fy3 = -gy3; 		/* rule 174 */
	DIHEDRAL_FORCE_ACCUMULATE(I3, 1, fy3 );
	 tx109 = -(tx40*tx65*tx85); 		/* rule 176 */
	 tx110 = tx40*tx85*tx86*tx87; 		/* rule 177 */
	 tx111 = -tx99; 		/* rule 178 */
	 tx112 = tx109 + tx110 + tx111; 		/* rule 179 */
	 gz3 = DeDPhi*tx112; 		/* rule 180 */
	 fz3 = -gz3; 		/* rule 181 */
	DIHEDRAL_FORCE_ACCUMULATE(I3, 2, fz3 );
	 gx4 = DeDPhi*tx39*tx65*tx85; 		/* rule 183 */
	 fx4 = -gx4; 		/* rule 184 */
	DIHEDRAL_FORCE_ACCUMULATE(I4, 0, fx4 );
	 gy4 = DeDPhi*tx38*tx65*tx85; 		/* rule 186 */
	 fy4 = -gy4; 		/* rule 187 */
	DIHEDRAL_FORCE_ACCUMULATE(I4, 1, fy4 );
	 gz4 = DeDPhi*tx40*tx65*tx85; 		/* rule 189 */
	 fz4 = -gz4; 		/* rule 190 */
	DIHEDRAL_FORCE_ACCUMULATE(I4, 2, fz4 );
	#ifdef DIHEDRAL_CALC_DIAGONAL_HESSIAN //[
	if (calcDiagonalHessian) {
	 tx113 = power2(y2); 		/* rule 194 */
	 tx114 = power2(y3); 		/* rule 195 */
	 tx115 = power2(z2); 		/* rule 196 */
	 tx116 = power2(z3); 		/* rule 197 */
	 tx117 = tx113*tx49; 		/* rule 198 */
	 tx118 = tx115*tx49; 		/* rule 199 */
	 tx119 = tx2*tx50; 		/* rule 200 */
	 tx120 = tx114*tx73; 		/* rule 201 */
	 tx121 = tx116*tx73; 		/* rule 202 */
	 tx122 = tx7*tx75; 		/* rule 203 */
	 tx123 = tx113*x1; 		/* rule 204 */
	 tx124 = tx114*x1; 		/* rule 205 */
	 tx125 = tx115*x1; 		/* rule 206 */
	 tx126 = tx116*x1; 		/* rule 207 */
	 tx127 = tx6*y1; 		/* rule 208 */
	 tx128 = tx1*y2; 		/* rule 209 */
	 tx129 = tx2*y2; 		/* rule 210 */
	 tx130 = tx6*y2; 		/* rule 211 */
	 tx131 = -2.*tx3*y3; 		/* rule 212 */
	 tx132 = x3*y2*y3; 		/* rule 213 */
	 tx133 = tx12*z1; 		/* rule 214 */
	 tx134 = -2.*tx15*z2; 		/* rule 215 */
	 tx135 = tx12*z3; 		/* rule 216 */
	 tx136 = tx7*z3; 		/* rule 217 */
	 tx137 = tx8*z3; 		/* rule 218 */
	 tx138 = x2*z2*z3; 		/* rule 219 */
	 tx139 = power2(DN); 		/* rule 220 */
	 tx140 = tx117 + tx118 + tx119 + tx120 + tx121 + tx122 + tx123 + tx124 + tx125 + tx126 + tx127 + tx128 + tx129 + tx130 + tx131 + tx132 + tx133 + tx134 + tx135 + tx136 + tx137 + tx138; 		/* rule 221 */
	 tx141 = power2(tx72); 		/* rule 222 */
	 tx142 = tx67 + tx68; 		/* rule 223 */
	 tx143 = 2.*DeDPhi*tx140*tx141*tx21*tx65; 		/* rule 224 */
	 tx144 = -(tx139*tx141*tx142*tx24*tx61*V); 		/* rule 225 */
	 dhx1x1 = tx143 + tx144; 		/* rule 226 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx145 = power2(x2); 		/* rule 228 */
	 tx146 = power2(x3); 		/* rule 229 */
	 tx147 = tx115*tx50; 		/* rule 230 */
	 tx148 = tx145*tx50; 		/* rule 231 */
	 tx149 = tx10*tx51; 		/* rule 232 */
	 tx150 = tx3*tx73; 		/* rule 233 */
	 tx151 = tx116*tx74; 		/* rule 234 */
	 tx152 = tx146*tx74; 		/* rule 235 */
	 tx153 = tx6*x1; 		/* rule 236 */
	 tx154 = -2.*tx2*x2; 		/* rule 237 */
	 tx155 = tx3*x3; 		/* rule 238 */
	 tx156 = tx5*x3; 		/* rule 239 */
	 tx157 = tx6*x3; 		/* rule 240 */
	 tx158 = tx115*y1; 		/* rule 241 */
	 tx159 = tx116*y1; 		/* rule 242 */
	 tx160 = tx145*y1; 		/* rule 243 */
	 tx161 = tx146*y1; 		/* rule 244 */
	 tx162 = x2*x3*y2; 		/* rule 245 */
	 tx163 = tx18*z1; 		/* rule 246 */
	 tx164 = tx10*z2; 		/* rule 247 */
	 tx165 = tx18*z2; 		/* rule 248 */
	 tx166 = tx9*z2; 		/* rule 249 */
	 tx167 = -2.*tx13*z3; 		/* rule 250 */
	 tx168 = y3*z2*z3; 		/* rule 251 */
	 tx169 = tx147 + tx148 + tx149 + tx150 + tx151 + tx152 + tx153 + tx154 + tx155 + tx156 + tx157 + tx158 + tx159 + tx160 + tx161 + tx162 + tx163 + tx164 + tx165 + tx166 + tx167 + tx168; 		/* rule 252 */
	 tx170 = 2.*DeDPhi*tx141*tx169*tx20*tx65; 		/* rule 253 */
	 tx171 = -(tx139*tx141*tx142*tx23*tx61*V); 		/* rule 254 */
	 dhy1y1 = tx170 + tx171; 		/* rule 255 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx172 = tx15*tx49; 		/* rule 257 */
	 tx173 = tx113*tx51; 		/* rule 258 */
	 tx174 = tx145*tx51; 		/* rule 259 */
	 tx175 = tx13*tx74; 		/* rule 260 */
	 tx176 = tx114*tx75; 		/* rule 261 */
	 tx177 = tx146*tx75; 		/* rule 262 */
	 tx178 = tx12*x1; 		/* rule 263 */
	 tx179 = tx11*x2; 		/* rule 264 */
	 tx180 = tx12*x2; 		/* rule 265 */
	 tx181 = tx15*x2; 		/* rule 266 */
	 tx182 = -2.*tx7*x3; 		/* rule 267 */
	 tx183 = tx18*y1; 		/* rule 268 */
	 tx184 = -2.*tx10*y2; 		/* rule 269 */
	 tx185 = tx13*y3; 		/* rule 270 */
	 tx186 = tx17*y3; 		/* rule 271 */
	 tx187 = tx18*y3; 		/* rule 272 */
	 tx188 = tx113*z1; 		/* rule 273 */
	 tx189 = tx114*z1; 		/* rule 274 */
	 tx190 = tx145*z1; 		/* rule 275 */
	 tx191 = tx146*z1; 		/* rule 276 */
	 tx192 = y2*y3*z2; 		/* rule 277 */
	 tx193 = x2*x3*z3; 		/* rule 278 */
	 tx194 = tx172 + tx173 + tx174 + tx175 + tx176 + tx177 + tx178 + tx179 + tx180 + tx181 + tx182 + tx183 + tx184 + tx185 + tx186 + tx187 + tx188 + tx189 + tx190 + tx191 + tx192 + tx193; 		/* rule 279 */
	 tx195 = 2.*DeDPhi*tx141*tx19*tx194*tx65; 		/* rule 280 */
	 tx196 = -(tx139*tx141*tx142*tx22*tx61*V); 		/* rule 281 */
	 dhz1z1 = tx195 + tx196; 		/* rule 282 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 tx197 = power2(y1); 		/* rule 284 */
	 tx198 = power2(z1); 		/* rule 285 */
	 tx199 = -tx123; 		/* rule 286 */
	 tx200 = -tx124; 		/* rule 287 */
	 tx201 = -tx125; 		/* rule 288 */
	 tx202 = -tx126; 		/* rule 289 */
	 tx203 = 2.*tx129; 		/* rule 290 */
	 tx204 = 2.*tx133; 		/* rule 291 */
	 tx205 = tx197*tx49; 		/* rule 292 */
	 tx206 = tx198*tx49; 		/* rule 293 */
	 tx207 = tx3*tx50; 		/* rule 294 */
	 tx208 = tx12*tx51; 		/* rule 295 */
	 tx209 = tx51*tx7; 		/* rule 296 */
	 tx210 = tx2*tx74; 		/* rule 297 */
	 tx211 = tx6*tx74; 		/* rule 298 */
	 tx212 = tx114*x2; 		/* rule 299 */
	 tx213 = tx116*x2; 		/* rule 300 */
	 tx214 = tx197*x2; 		/* rule 301 */
	 tx215 = tx198*x2; 		/* rule 302 */
	 tx216 = tx113*x3; 		/* rule 303 */
	 tx217 = tx115*x3; 		/* rule 304 */
	 tx218 = -(tx3*y1); 		/* rule 305 */
	 tx219 = x2*y1*y2; 		/* rule 306 */
	 tx220 = tx1*y3; 		/* rule 307 */
	 tx221 = tx2*y3; 		/* rule 308 */
	 tx222 = 2.*tx3*y3; 		/* rule 309 */
	 tx223 = tx4*y3; 		/* rule 310 */
	 tx224 = x1*y1*y3; 		/* rule 311 */
	 tx225 = tx11*z1; 		/* rule 312 */
	 tx226 = tx15*z1; 		/* rule 313 */
	 tx227 = 2.*tx15*z2; 		/* rule 314 */
	 tx228 = tx16*z2; 		/* rule 315 */
	 tx229 = tx7*z2; 		/* rule 316 */
	 tx230 = tx8*z2; 		/* rule 317 */
	 tx231 = tx11*z3; 		/* rule 318 */
	 tx232 = x3*z1*z3; 		/* rule 319 */
	 tx233 = tx117 + tx118 + tx122 + tx123 + tx125 + tx128 + tx130 + tx138 + tx203 + tx204 + tx205 + tx206 + tx207 + tx209 + tx214 + tx215 + tx218 + tx220 + tx224 + tx225 + tx226 + tx231; 		/* rule 320 */
	 tx234 = tx199 + tx200 + tx201 + tx202 + tx208 + tx209 + tx210 + tx211 + tx212 + tx213 + tx216 + tx217 + tx219 + tx220 + tx221 + tx222 + tx223 + tx227 + tx228 + tx229 + tx230 + tx232; 		/* rule 321 */
	 tx235 = power2(y4); 		/* rule 322 */
	 tx236 = power2(z4); 		/* rule 323 */
	 tx237 = -2.*tx129; 		/* rule 324 */
	 tx238 = -2.*tx133; 		/* rule 325 */
	 tx239 = -tx226; 		/* rule 326 */
	 tx240 = tx26*tx50; 		/* rule 327 */
	 tx241 = tx29*tx50; 		/* rule 328 */
	 tx242 = tx21*tx233*tx61; 		/* rule 329 */
	 tx243 = tx197*tx73; 		/* rule 330 */
	 tx244 = tx198*tx73; 		/* rule 331 */
	 tx245 = tx235*tx73; 		/* rule 332 */
	 tx246 = tx236*tx73; 		/* rule 333 */
	 tx247 = tx34*tx75; 		/* rule 334 */
	 tx248 = tx21*tx234*tx88; 		/* rule 335 */
	 tx249 = tx197*x3; 		/* rule 336 */
	 tx250 = tx198*x3; 		/* rule 337 */
	 tx251 = tx235*x3; 		/* rule 338 */
	 tx252 = tx236*x3; 		/* rule 339 */
	 tx253 = tx113*x4; 		/* rule 340 */
	 tx254 = tx114*x4; 		/* rule 341 */
	 tx255 = tx115*x4; 		/* rule 342 */
	 tx256 = tx116*x4; 		/* rule 343 */
	 tx257 = tx3*y1; 		/* rule 344 */
	 tx258 = tx5*y1; 		/* rule 345 */
	 tx259 = tx28*y2; 		/* rule 346 */
	 tx260 = tx29*y2; 		/* rule 347 */
	 tx261 = -2.*tx26*y3; 		/* rule 348 */
	 tx262 = tx3*y3; 		/* rule 349 */
	 tx263 = tx26*y4; 		/* rule 350 */
	 tx264 = tx27*y4; 		/* rule 351 */
	 tx265 = tx4*y4; 		/* rule 352 */
	 tx266 = tx6*y4; 		/* rule 353 */
	 tx267 = 2.*tx266; 		/* rule 354 */
	 tx268 = tx15*z2; 		/* rule 355 */
	 tx269 = -2.*tx32*z2; 		/* rule 356 */
	 tx270 = x1*z1*z2; 		/* rule 357 */
	 tx271 = tx30*z3; 		/* rule 358 */
	 tx272 = tx34*z3; 		/* rule 359 */
	 tx273 = 2.*tx272; 		/* rule 360 */
	 tx274 = tx35*z3; 		/* rule 361 */
	 tx275 = -(tx12*z4); 		/* rule 362 */
	 tx276 = tx12*z4; 		/* rule 363 */
	 tx277 = -(tx32*z4); 		/* rule 364 */
	 tx278 = x4*z2*z4; 		/* rule 365 */
	 tx279 = tx267; 		/* rule 366 */
	 tx280 = tx273; 		/* rule 367 */
	 tx281 = tx242 + tx248; 		/* rule 368 */
	 tx282 = tx127 + tx136 + tx199 + tx201 + tx211 + tx216 + tx217 + tx219 + tx228 + tx229 + tx237 + tx238 + tx239 + tx243 + tx244 + tx249 + tx250 + tx257 + tx258 + tx262 + tx268 + tx270; 		/* rule 369 */
	 tx283 = tx117 + tx118 + tx120 + tx121 + tx130 + tx132 + tx135 + tx138 + tx241 + tx247 + tx253 + tx254 + tx255 + tx256 + tx259 + tx260 + tx261 + tx266 + tx269 + tx272 + tx274 + tx276; 		/* rule 370 */
	 tx284 = tx120 + tx121 + tx132 + tx135 + tx240 + tx241 + tx245 + tx246 + tx251 + tx252 + tx254 + tx256 + tx263 + tx264 + tx265 + tx271 + tx274 + tx275 + tx277 + tx278 + tx279 + tx280; 		/* rule 371 */
	 tx285 = power2(tx85); 		/* rule 372 */
	 tx907 = tx61; 		/* rule 373 */
	 tx908 = tx86; 		/* rule 374 */
	 tx909 = reciprocal(tx907); 		/* rule 375 */
	 tx286 = tx908*tx909; 		/* rule 376 */
	 tx287 = 2.*tx140*tx141*tx21*tx65; 		/* rule 377 */
	 tx288 = tx140*tx21*tx286*tx72; 		/* rule 378 */
	 tx289 = -(tx283*tx286*tx39*tx85); 		/* rule 379 */
	 tx290 = -2.*tx141*tx281*tx908; 		/* rule 380 */
	 tx291 = -2.*tx284*tx285*tx39*tx87*tx908; 		/* rule 381 */
	 tx292 = 2.*tx141*tx21*tx282*tx88*tx908; 		/* rule 382 */
	 tx293 = tx287 + tx288 + tx289 + tx290 + tx291 + tx292; 		/* rule 383 */
	 tx294 = power2(tx92); 		/* rule 384 */
	 tx295 = DeDPhi*tx293; 		/* rule 385 */
	 tx296 = -(tx139*tx142*tx294*V); 		/* rule 386 */
	 dhx2x2 = tx295 + tx296; 		/* rule 387 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 tx297 = power2(x1); 		/* rule 389 */
	 tx298 = 2.*tx153; 		/* rule 390 */
	 tx299 = -tx158; 		/* rule 391 */
	 tx300 = -tx159; 		/* rule 392 */
	 tx301 = -tx160; 		/* rule 393 */
	 tx302 = -tx161; 		/* rule 394 */
	 tx303 = 2.*tx164; 		/* rule 395 */
	 tx304 = tx3*tx49; 		/* rule 396 */
	 tx305 = tx198*tx50; 		/* rule 397 */
	 tx306 = tx297*tx50; 		/* rule 398 */
	 tx307 = tx13*tx51; 		/* rule 399 */
	 tx308 = tx49*tx6; 		/* rule 400 */
	 tx309 = tx10*tx75; 		/* rule 401 */
	 tx310 = tx18*tx75; 		/* rule 402 */
	 tx311 = tx1*x1; 		/* rule 403 */
	 tx312 = tx2*x1; 		/* rule 404 */
	 tx313 = 2.*tx2*x2; 		/* rule 405 */
	 tx314 = tx3*x2; 		/* rule 406 */
	 tx315 = tx4*x2; 		/* rule 407 */
	 tx316 = tx5*x2; 		/* rule 408 */
	 tx317 = tx1*x3; 		/* rule 409 */
	 tx318 = tx116*y2; 		/* rule 410 */
	 tx319 = tx146*y2; 		/* rule 411 */
	 tx320 = tx198*y2; 		/* rule 412 */
	 tx321 = tx297*y2; 		/* rule 413 */
	 tx322 = tx115*y3; 		/* rule 414 */
	 tx323 = tx145*y3; 		/* rule 415 */
	 tx324 = x1*x3*y3; 		/* rule 416 */
	 tx325 = -(tx13*z1); 		/* rule 417 */
	 tx326 = y2*z1*z2; 		/* rule 418 */
	 tx327 = tx10*z3; 		/* rule 419 */
	 tx328 = 2.*tx13*z3; 		/* rule 420 */
	 tx329 = tx14*z3; 		/* rule 421 */
	 tx330 = tx9*z3; 		/* rule 422 */
	 tx331 = y1*z1*z3; 		/* rule 423 */
	 tx332 = tx299 + tx300 + tx301 + tx302 + tx304 + tx308 + tx309 + tx310 + tx313 + tx314 + tx315 + tx316 + tx318 + tx319 + tx322 + tx323 + tx324 + tx326 + tx327 + tx328 + tx329 + tx330; 		/* rule 424 */
	 tx333 = tx147 + tx148 + tx150 + tx158 + tx160 + tx162 + tx165 + tx166 + tx298 + tx303 + tx304 + tx305 + tx306 + tx307 + tx311 + tx312 + tx317 + tx320 + tx321 + tx325 + tx330 + tx331; 		/* rule 425 */
	 tx334 = power2(x4); 		/* rule 426 */
	 tx335 = -2.*tx153; 		/* rule 427 */
	 tx336 = -2.*tx164; 		/* rule 428 */
	 tx337 = -tx312; 		/* rule 429 */
	 tx338 = tx31*tx51; 		/* rule 430 */
	 tx339 = tx37*tx51; 		/* rule 431 */
	 tx340 = tx20*tx333*tx907; 		/* rule 432 */
	 tx341 = tx26*tx73; 		/* rule 433 */
	 tx342 = tx198*tx74; 		/* rule 434 */
	 tx343 = tx236*tx74; 		/* rule 435 */
	 tx344 = tx297*tx74; 		/* rule 436 */
	 tx345 = tx334*tx74; 		/* rule 437 */
	 tx346 = tx20*tx332*tx88; 		/* rule 438 */
	 tx347 = tx2*x2; 		/* rule 439 */
	 tx348 = -2.*tx29*x2; 		/* rule 440 */
	 tx349 = tx26*x3; 		/* rule 441 */
	 tx350 = 2.*tx349; 		/* rule 442 */
	 tx351 = tx27*x3; 		/* rule 443 */
	 tx352 = tx28*x3; 		/* rule 444 */
	 tx353 = -(tx29*x4); 		/* rule 445 */
	 tx354 = -(tx6*x4); 		/* rule 446 */
	 tx355 = tx6*x4; 		/* rule 447 */
	 tx356 = x1*x2*y1; 		/* rule 448 */
	 tx357 = tx198*y3; 		/* rule 449 */
	 tx358 = tx236*y3; 		/* rule 450 */
	 tx359 = tx297*y3; 		/* rule 451 */
	 tx360 = tx334*y3; 		/* rule 452 */
	 tx361 = tx115*y4; 		/* rule 453 */
	 tx362 = tx116*y4; 		/* rule 454 */
	 tx363 = tx145*y4; 		/* rule 455 */
	 tx364 = tx146*y4; 		/* rule 456 */
	 tx365 = x2*x4*y4; 		/* rule 457 */
	 tx366 = tx13*z1; 		/* rule 458 */
	 tx367 = tx17*z1; 		/* rule 459 */
	 tx368 = tx36*z2; 		/* rule 460 */
	 tx369 = tx37*z2; 		/* rule 461 */
	 tx370 = tx13*z3; 		/* rule 462 */
	 tx371 = -2.*tx31*z3; 		/* rule 463 */
	 tx372 = tx14*z4; 		/* rule 464 */
	 tx373 = tx18*z4; 		/* rule 465 */
	 tx374 = 2.*tx373; 		/* rule 466 */
	 tx375 = tx31*z4; 		/* rule 467 */
	 tx376 = tx33*z4; 		/* rule 468 */
	 tx377 = tx350; 		/* rule 469 */
	 tx378 = tx374; 		/* rule 470 */
	 tx379 = tx340 + tx346; 		/* rule 471 */
	 tx380 = tx155 + tx163 + tx299 + tx301 + tx310 + tx314 + tx315 + tx322 + tx323 + tx326 + tx335 + tx336 + tx337 + tx342 + tx344 + tx347 + tx356 + tx357 + tx359 + tx366 + tx367 + tx370; 		/* rule 472 */
	 tx381 = tx147 + tx148 + tx151 + tx152 + tx157 + tx162 + tx165 + tx168 + tx339 + tx341 + tx348 + tx349 + tx351 + tx355 + tx361 + tx362 + tx363 + tx364 + tx368 + tx369 + tx371 + tx373; 		/* rule 473 */
	 tx382 = tx151 + tx152 + tx157 + tx168 + tx338 + tx339 + tx343 + tx345 + tx351 + tx352 + tx353 + tx354 + tx358 + tx360 + tx362 + tx364 + tx365 + tx372 + tx375 + tx376 + tx377 + tx378; 		/* rule 474 */
	 tx383 = 2.*tx141*tx169*tx20*tx65; 		/* rule 475 */
	 tx384 = tx169*tx20*tx286*tx72; 		/* rule 476 */
	 tx385 = -(tx286*tx38*tx381*tx85); 		/* rule 477 */
	 tx386 = -2.*tx141*tx379*tx908; 		/* rule 478 */
	 tx387 = -2.*tx285*tx38*tx382*tx87*tx908; 		/* rule 479 */
	 tx388 = 2.*tx141*tx20*tx380*tx88*tx908; 		/* rule 480 */
	 tx389 = tx383 + tx384 + tx385 + tx386 + tx387 + tx388; 		/* rule 481 */
	 tx390 = power2(tx96); 		/* rule 482 */
	 tx391 = DeDPhi*tx389; 		/* rule 483 */
	 tx392 = -(tx139*tx142*tx390*V); 		/* rule 484 */
	 dhy2y2 = tx391 + tx392; 		/* rule 485 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 tx393 = 2.*tx181; 		/* rule 487 */
	 tx394 = 2.*tx183; 		/* rule 488 */
	 tx395 = -tx188; 		/* rule 489 */
	 tx396 = -tx189; 		/* rule 490 */
	 tx397 = -tx190; 		/* rule 491 */
	 tx398 = -tx191; 		/* rule 492 */
	 tx399 = tx13*tx50; 		/* rule 493 */
	 tx400 = tx18*tx50; 		/* rule 494 */
	 tx401 = tx197*tx51; 		/* rule 495 */
	 tx402 = tx297*tx51; 		/* rule 496 */
	 tx403 = tx49*tx7; 		/* rule 497 */
	 tx404 = tx12*tx73; 		/* rule 498 */
	 tx405 = tx15*tx73; 		/* rule 499 */
	 tx406 = -(tx7*x1); 		/* rule 500 */
	 tx407 = tx11*x3; 		/* rule 501 */
	 tx408 = tx15*x3; 		/* rule 502 */
	 tx409 = tx16*x3; 		/* rule 503 */
	 tx410 = 2.*tx7*x3; 		/* rule 504 */
	 tx411 = tx10*y1; 		/* rule 505 */
	 tx412 = tx9*y1; 		/* rule 506 */
	 tx413 = 2.*tx10*y2; 		/* rule 507 */
	 tx414 = tx13*y2; 		/* rule 508 */
	 tx415 = tx14*y2; 		/* rule 509 */
	 tx416 = tx17*y2; 		/* rule 510 */
	 tx417 = tx9*y3; 		/* rule 511 */
	 tx418 = x1*x3*z1; 		/* rule 512 */
	 tx419 = tx114*z2; 		/* rule 513 */
	 tx420 = tx146*z2; 		/* rule 514 */
	 tx421 = tx197*z2; 		/* rule 515 */
	 tx422 = tx297*z2; 		/* rule 516 */
	 tx423 = x1*x2*z2; 		/* rule 517 */
	 tx424 = tx113*z3; 		/* rule 518 */
	 tx425 = tx145*z3; 		/* rule 519 */
	 tx426 = y1*y3*z3; 		/* rule 520 */
	 tx427 = tx173 + tx174 + tx175 + tx179 + tx180 + tx188 + tx190 + tx192 + tx393 + tx394 + tx399 + tx401 + tx402 + tx403 + tx406 + tx407 + tx411 + tx412 + tx417 + tx418 + tx421 + tx422; 		/* rule 521 */
	 tx428 = tx395 + tx396 + tx397 + tx398 + tx399 + tx400 + tx404 + tx405 + tx407 + tx408 + tx409 + tx410 + tx413 + tx414 + tx415 + tx416 + tx419 + tx420 + tx423 + tx424 + tx425 + tx426; 		/* rule 522 */
	 tx429 = -2.*tx181; 		/* rule 523 */
	 tx430 = -2.*tx183; 		/* rule 524 */
	 tx431 = -tx411; 		/* rule 525 */
	 tx432 = tx32*tx49; 		/* rule 526 */
	 tx433 = tx34*tx49; 		/* rule 527 */
	 tx434 = tx19*tx427*tx907; 		/* rule 528 */
	 tx435 = tx31*tx74; 		/* rule 529 */
	 tx436 = tx197*tx75; 		/* rule 530 */
	 tx437 = tx235*tx75; 		/* rule 531 */
	 tx438 = tx297*tx75; 		/* rule 532 */
	 tx439 = tx334*tx75; 		/* rule 533 */
	 tx440 = tx19*tx428*tx88; 		/* rule 534 */
	 tx441 = tx7*x1; 		/* rule 535 */
	 tx442 = tx8*x1; 		/* rule 536 */
	 tx443 = tx30*x2; 		/* rule 537 */
	 tx444 = tx32*x2; 		/* rule 538 */
	 tx445 = -2.*tx34*x3; 		/* rule 539 */
	 tx446 = tx7*x3; 		/* rule 540 */
	 tx447 = tx12*x4; 		/* rule 541 */
	 tx448 = 2.*tx447; 		/* rule 542 */
	 tx449 = tx16*x4; 		/* rule 543 */
	 tx450 = tx34*x4; 		/* rule 544 */
	 tx451 = tx35*x4; 		/* rule 545 */
	 tx452 = tx10*y2; 		/* rule 546 */
	 tx453 = -2.*tx37*y2; 		/* rule 547 */
	 tx454 = tx31*y3; 		/* rule 548 */
	 tx455 = 2.*tx454; 		/* rule 549 */
	 tx456 = tx33*y3; 		/* rule 550 */
	 tx457 = tx36*y3; 		/* rule 551 */
	 tx458 = -(tx18*y4); 		/* rule 552 */
	 tx459 = tx18*y4; 		/* rule 553 */
	 tx460 = -(tx37*y4); 		/* rule 554 */
	 tx461 = y1*y2*z1; 		/* rule 555 */
	 tx462 = tx197*z3; 		/* rule 556 */
	 tx463 = tx235*z3; 		/* rule 557 */
	 tx464 = tx297*z3; 		/* rule 558 */
	 tx465 = tx334*z3; 		/* rule 559 */
	 tx466 = tx113*z4; 		/* rule 560 */
	 tx467 = tx114*z4; 		/* rule 561 */
	 tx468 = tx145*z4; 		/* rule 562 */
	 tx469 = tx146*z4; 		/* rule 563 */
	 tx470 = y2*y4*z4; 		/* rule 564 */
	 tx471 = tx448; 		/* rule 565 */
	 tx472 = tx455; 		/* rule 566 */
	 tx473 = tx434 + tx440; 		/* rule 567 */
	 tx474 = tx178 + tx185 + tx395 + tx397 + tx404 + tx414 + tx415 + tx423 + tx424 + tx425 + tx429 + tx430 + tx431 + tx436 + tx438 + tx441 + tx442 + tx446 + tx452 + tx461 + tx462 + tx464; 		/* rule 568 */
	 tx475 = tx173 + tx174 + tx176 + tx177 + tx180 + tx187 + tx192 + tx193 + tx432 + tx435 + tx443 + tx444 + tx445 + tx447 + tx453 + tx454 + tx456 + tx459 + tx466 + tx467 + tx468 + tx469; 		/* rule 569 */
	 tx476 = tx176 + tx177 + tx187 + tx193 + tx432 + tx433 + tx437 + tx439 + tx449 + tx450 + tx451 + tx456 + tx457 + tx458 + tx460 + tx463 + tx465 + tx467 + tx469 + tx470 + tx471 + tx472; 		/* rule 570 */
	 tx477 = 2.*tx141*tx19*tx194*tx65; 		/* rule 571 */
	 tx478 = tx19*tx194*tx286*tx72; 		/* rule 572 */
	 tx479 = -(tx286*tx40*tx475*tx85); 		/* rule 573 */
	 tx480 = -2.*tx141*tx473*tx908; 		/* rule 574 */
	 tx481 = -2.*tx285*tx40*tx476*tx87*tx908; 		/* rule 575 */
	 tx482 = 2.*tx141*tx19*tx474*tx88*tx908; 		/* rule 576 */
	 tx483 = power2(tx100); 		/* rule 577 */
	 tx484 = tx477 + tx478 + tx479 + tx480 + tx481 + tx482; 		/* rule 578 */
	 tx485 = DeDPhi*tx484; 		/* rule 579 */
	 tx486 = -(tx139*tx142*tx483*V); 		/* rule 580 */
	 dhz2z2 = tx485 + tx486; 		/* rule 581 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	 tx487 = -tx253; 		/* rule 583 */
	 tx488 = -tx254; 		/* rule 584 */
	 tx489 = -tx255; 		/* rule 585 */
	 tx490 = -tx256; 		/* rule 586 */
	 tx491 = -tx263; 		/* rule 587 */
	 tx492 = -2.*tx266; 		/* rule 588 */
	 tx493 = -tx266; 		/* rule 589 */
	 tx494 = -2.*tx272; 		/* rule 590 */
	 tx495 = tx235*tx49; 		/* rule 591 */
	 tx496 = tx236*tx49; 		/* rule 592 */
	 tx497 = tx235*x2; 		/* rule 593 */
	 tx498 = tx236*x2; 		/* rule 594 */
	 tx499 = tx26*y3; 		/* rule 595 */
	 tx500 = 2.*tx499; 		/* rule 596 */
	 tx501 = tx29*y3; 		/* rule 597 */
	 tx502 = x2*y2*y4; 		/* rule 598 */
	 tx503 = x4*y3*y4; 		/* rule 599 */
	 tx504 = tx32*z2; 		/* rule 600 */
	 tx505 = 2.*tx504; 		/* rule 601 */
	 tx506 = tx34*z2; 		/* rule 602 */
	 tx507 = tx16*z4; 		/* rule 603 */
	 tx508 = tx30*z4; 		/* rule 604 */
	 tx509 = tx32*z4; 		/* rule 605 */
	 tx510 = x3*z3*z4; 		/* rule 606 */
	 tx511 = tx500; 		/* rule 607 */
	 tx512 = tx505; 		/* rule 608 */
	 tx513 = tx208 + tx212 + tx213 + tx223 + tx260 + tx276 + tx488 + tx490 + tx491 + tx492 + tx494 + tx495 + tx496 + tx497 + tx498 + tx499 + tx501 + tx503 + tx504 + tx508 + tx509 + tx510; 		/* rule 609 */
	 tx514 = tx208 + tx211 + tx212 + tx213 + tx216 + tx217 + tx223 + tx228 + tx265 + tx275 + tx487 + tx488 + tx489 + tx490 + tx493 + tx501 + tx502 + tx506 + tx507 + tx510 + tx511 + tx512; 		/* rule 610 */
	 tx515 = tx39*tx513*tx907; 		/* rule 611 */
	 tx516 = tx39*tx514*tx87; 		/* rule 612 */
	 tx517 = tx515 + tx516; 		/* rule 613 */
	 tx518 = -2.*tx283*tx285*tx39*tx65; 		/* rule 614 */
	 tx519 = -2.*tx285*tx517*tx908; 		/* rule 615 */
	 tx520 = power2(tx104); 		/* rule 616 */
	 tx521 = tx288 + tx289 + tx291 + tx292 + tx518 + tx519; 		/* rule 617 */
	 tx522 = DeDPhi*tx521; 		/* rule 618 */
	 tx523 = -(tx139*tx142*tx520*V); 		/* rule 619 */
	 dhx3x3 = tx522 + tx523; 		/* rule 620 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 0, dhx3x3);
	 tx524 = -2.*tx349; 		/* rule 622 */
	 tx525 = -tx361; 		/* rule 623 */
	 tx526 = -tx362; 		/* rule 624 */
	 tx527 = -tx363; 		/* rule 625 */
	 tx528 = -tx364; 		/* rule 626 */
	 tx529 = -2.*tx373; 		/* rule 627 */
	 tx530 = -tx373; 		/* rule 628 */
	 tx531 = -tx375; 		/* rule 629 */
	 tx532 = tx236*tx50; 		/* rule 630 */
	 tx533 = tx334*tx50; 		/* rule 631 */
	 tx534 = tx26*x2; 		/* rule 632 */
	 tx535 = tx29*x2; 		/* rule 633 */
	 tx536 = 2.*tx535; 		/* rule 634 */
	 tx537 = tx28*x4; 		/* rule 635 */
	 tx538 = tx29*x4; 		/* rule 636 */
	 tx539 = tx4*x4; 		/* rule 637 */
	 tx540 = tx236*y2; 		/* rule 638 */
	 tx541 = tx334*y2; 		/* rule 639 */
	 tx542 = x3*x4*y3; 		/* rule 640 */
	 tx543 = tx31*z3; 		/* rule 641 */
	 tx544 = 2.*tx543; 		/* rule 642 */
	 tx545 = tx37*z3; 		/* rule 643 */
	 tx546 = y2*z2*z4; 		/* rule 644 */
	 tx547 = y4*z3*z4; 		/* rule 645 */
	 tx548 = tx536; 		/* rule 646 */
	 tx549 = tx544; 		/* rule 647 */
	 tx550 = tx308 + tx318 + tx319 + tx329 + tx355 + tx369 + tx524 + tx526 + tx528 + tx529 + tx531 + tx532 + tx533 + tx535 + tx537 + tx538 + tx540 + tx541 + tx542 + tx543 + tx545 + tx547; 		/* rule 648 */
	 tx551 = tx308 + tx310 + tx315 + tx318 + tx319 + tx322 + tx323 + tx329 + tx354 + tx372 + tx525 + tx526 + tx527 + tx528 + tx530 + tx534 + tx539 + tx542 + tx545 + tx546 + tx548 + tx549; 		/* rule 649 */
	 tx552 = tx38*tx550*tx907; 		/* rule 650 */
	 tx553 = tx38*tx551*tx87; 		/* rule 651 */
	 tx554 = tx552 + tx553; 		/* rule 652 */
	 tx555 = -2.*tx285*tx38*tx381*tx65; 		/* rule 653 */
	 tx556 = -2.*tx285*tx554*tx908; 		/* rule 654 */
	 tx557 = power2(tx108); 		/* rule 655 */
	 tx558 = tx384 + tx385 + tx387 + tx388 + tx555 + tx556; 		/* rule 656 */
	 tx559 = DeDPhi*tx558; 		/* rule 657 */
	 tx560 = -(tx139*tx142*tx557*V); 		/* rule 658 */
	 dhy3y3 = tx559 + tx560; 		/* rule 659 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 1, dhy3y3);
	 tx561 = -2.*tx447; 		/* rule 661 */
	 tx562 = -tx447; 		/* rule 662 */
	 tx563 = -tx450; 		/* rule 663 */
	 tx564 = -2.*tx454; 		/* rule 664 */
	 tx565 = -tx466; 		/* rule 665 */
	 tx566 = -tx467; 		/* rule 666 */
	 tx567 = -tx468; 		/* rule 667 */
	 tx568 = -tx469; 		/* rule 668 */
	 tx569 = tx235*tx51; 		/* rule 669 */
	 tx570 = tx334*tx51; 		/* rule 670 */
	 tx571 = tx32*x3; 		/* rule 671 */
	 tx572 = tx34*x3; 		/* rule 672 */
	 tx573 = 2.*tx572; 		/* rule 673 */
	 tx574 = tx31*y2; 		/* rule 674 */
	 tx575 = tx37*y2; 		/* rule 675 */
	 tx576 = 2.*tx575; 		/* rule 676 */
	 tx577 = tx14*y4; 		/* rule 677 */
	 tx578 = tx36*y4; 		/* rule 678 */
	 tx579 = tx37*y4; 		/* rule 679 */
	 tx580 = tx235*z2; 		/* rule 680 */
	 tx581 = tx334*z2; 		/* rule 681 */
	 tx582 = x2*x4*z2; 		/* rule 682 */
	 tx583 = y3*y4*z3; 		/* rule 683 */
	 tx584 = x3*x4*z4; 		/* rule 684 */
	 tx585 = tx573; 		/* rule 685 */
	 tx586 = tx576; 		/* rule 686 */
	 tx587 = tx400 + tx409 + tx419 + tx420 + tx444 + tx459 + tx561 + tx563 + tx564 + tx566 + tx568 + tx569 + tx570 + tx571 + tx572 + tx575 + tx578 + tx579 + tx580 + tx581 + tx583 + tx584; 		/* rule 687 */
	 tx588 = tx400 + tx404 + tx409 + tx415 + tx419 + tx420 + tx424 + tx425 + tx449 + tx458 + tx562 + tx565 + tx566 + tx567 + tx568 + tx571 + tx574 + tx577 + tx582 + tx583 + tx585 + tx586; 		/* rule 688 */
	 tx589 = tx40*tx587*tx907; 		/* rule 689 */
	 tx590 = tx40*tx588*tx87; 		/* rule 690 */
	 tx591 = tx589 + tx590; 		/* rule 691 */
	 tx592 = -2.*tx285*tx40*tx475*tx65; 		/* rule 692 */
	 tx593 = -2.*tx285*tx591*tx908; 		/* rule 693 */
	 tx594 = power2(tx112); 		/* rule 694 */
	 tx595 = tx478 + tx479 + tx481 + tx482 + tx592 + tx593; 		/* rule 695 */
	 tx596 = DeDPhi*tx595; 		/* rule 696 */
	 tx597 = -(tx139*tx142*tx594*V); 		/* rule 697 */
	 dhz3z3 = tx596 + tx597; 		/* rule 698 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I3, 2, dhz3z3);
	 tx598 = DeDPhi*tx518; 		/* rule 700 */
	 tx599 = -(tx139*tx142*tx285*tx42*tx907*V); 		/* rule 701 */
	 dhx4x4 = tx598 + tx599; 		/* rule 702 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 0, dhx4x4);
	 tx600 = DeDPhi*tx555; 		/* rule 704 */
	 tx601 = -(tx139*tx142*tx285*tx41*tx907*V); 		/* rule 705 */
	 dhy4y4 = tx600 + tx601; 		/* rule 706 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I4, 1, I4, 1, dhy4y4);
	 tx602 = DeDPhi*tx592; 		/* rule 708 */
	 tx603 = -(tx139*tx142*tx285*tx43*tx907*V); 		/* rule 709 */
	 dhz4z4 = tx602 + tx603; 		/* rule 710 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I4, 2, I4, 2, dhz4z4);
	#ifdef DIHEDRAL_CALC_OFF_DIAGONAL_HESSIAN //[
	if (calcOffDiagonalHessian) { 
	 tx604 = tx140*tx20; 		/* rule 714 */
	 tx605 = tx169*tx21; 		/* rule 715 */
	 tx606 = tx604 + tx605; 		/* rule 716 */
	 tx607 = DeDPhi*tx141*tx606*tx65; 		/* rule 717 */
	 tx608 = -(tx139*tx141*tx142*tx20*tx21*tx907*V); 		/* rule 718 */
	 ohx1y1 = tx607 + tx608; 		/* rule 719 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tx609 = tx140*tx19; 		/* rule 721 */
	 tx610 = tx194*tx21; 		/* rule 722 */
	 tx611 = tx609 + tx610; 		/* rule 723 */
	 tx612 = DeDPhi*tx141*tx611*tx65; 		/* rule 724 */
	 tx613 = -(tx139*tx141*tx142*tx19*tx21*tx907*V); 		/* rule 725 */
	 ohx1z1 = tx612 + tx613; 		/* rule 726 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx614 = -2.*tx140*tx141*tx21*tx65; 		/* rule 728 */
	 tx615 = tx141*tx281*tx908; 		/* rule 729 */
	 tx616 = tx614 + tx615; 		/* rule 730 */
	 tx617 = DeDPhi*tx616; 		/* rule 731 */
	 tx618 = tx139*tx142*tx89*tx92*V; 		/* rule 732 */
	 ohx1x2 = tx617 + tx618; 		/* rule 733 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx619 = tx20*tx233*tx907; 		/* rule 735 */
	 tx620 = tx21*tx332*tx88; 		/* rule 736 */
	 tx621 = tx619 + tx620; 		/* rule 737 */
	 tx622 = -(tx141*tx606*tx65); 		/* rule 738 */
	 tx623 = tx141*tx621*tx908; 		/* rule 739 */
	 tx624 = tx622 + tx623; 		/* rule 740 */
	 tx625 = DeDPhi*tx624; 		/* rule 741 */
	 tx626 = tx139*tx142*tx89*tx96*V; 		/* rule 742 */
	 ohx1y2 = tx625 + tx626; 		/* rule 743 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx627 = tx19*tx233*tx907; 		/* rule 745 */
	 tx628 = tx21*tx428*tx88; 		/* rule 746 */
	 tx629 = tx627 + tx628; 		/* rule 747 */
	 tx630 = -(tx141*tx611*tx65); 		/* rule 748 */
	 tx631 = tx141*tx629*tx908; 		/* rule 749 */
	 tx632 = tx630 + tx631; 		/* rule 750 */
	 tx633 = DeDPhi*tx632; 		/* rule 751 */
	 tx634 = tx100*tx139*tx142*tx89*V; 		/* rule 752 */
	 ohx1z2 = tx633 + tx634; 		/* rule 753 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tx635 = -(DeDPhi*tx615); 		/* rule 755 */
	 tx636 = tx104*tx139*tx142*tx89*V; 		/* rule 756 */
	 ohx1x3 = tx635 + tx636; 		/* rule 757 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 0, ohx1x3);
	 tx637 = -(DeDPhi*tx623); 		/* rule 759 */
	 tx638 = tx108*tx139*tx142*tx89*V; 		/* rule 760 */
	 ohx1y3 = tx637 + tx638; 		/* rule 761 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 1, ohx1y3);
	 tx639 = -(DeDPhi*tx631); 		/* rule 763 */
	 tx640 = tx112*tx139*tx142*tx89*V; 		/* rule 764 */
	 ohx1z3 = tx639 + tx640; 		/* rule 765 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 2, ohx1z3);
	 ohx1x4 = tx139*tx142*tx46*tx72*tx85*tx907*V; 		/* rule 767 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 0, ohx1x4);
	 ohx1y4 = tx139*tx142*tx21*tx38*tx72*tx85*tx907*V; 		/* rule 769 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 1, ohx1y4);
	 ohx1z4 = tx139*tx142*tx21*tx40*tx72*tx85*tx907*V; 		/* rule 771 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 2, ohx1z4);
	 tx641 = tx169*tx19; 		/* rule 773 */
	 tx642 = tx194*tx20; 		/* rule 774 */
	 tx643 = tx641 + tx642; 		/* rule 775 */
	 tx644 = DeDPhi*tx141*tx643*tx65; 		/* rule 776 */
	 tx645 = -(tx139*tx141*tx142*tx19*tx20*tx907*V); 		/* rule 777 */
	 ohy1z1 = tx644 + tx645; 		/* rule 778 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 tx646 = tx21*tx333*tx907; 		/* rule 780 */
	 tx647 = tx20*tx234*tx88; 		/* rule 781 */
	 tx648 = tx646 + tx647; 		/* rule 782 */
	 tx649 = tx141*tx648*tx908; 		/* rule 783 */
	 tx650 = tx622 + tx649; 		/* rule 784 */
	 tx651 = DeDPhi*tx650; 		/* rule 785 */
	 tx652 = tx139*tx142*tx92*tx93*V; 		/* rule 786 */
	 ohy1x2 = tx651 + tx652; 		/* rule 787 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx653 = -2.*tx141*tx169*tx20*tx65; 		/* rule 789 */
	 tx654 = tx141*tx379*tx908; 		/* rule 790 */
	 tx655 = tx653 + tx654; 		/* rule 791 */
	 tx656 = DeDPhi*tx655; 		/* rule 792 */
	 tx657 = tx139*tx142*tx93*tx96*V; 		/* rule 793 */
	 ohy1y2 = tx656 + tx657; 		/* rule 794 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx658 = tx19*tx333*tx907; 		/* rule 796 */
	 tx659 = tx20*tx428*tx88; 		/* rule 797 */
	 tx660 = tx658 + tx659; 		/* rule 798 */
	 tx661 = -(tx141*tx643*tx65); 		/* rule 799 */
	 tx662 = tx141*tx660*tx908; 		/* rule 800 */
	 tx663 = tx661 + tx662; 		/* rule 801 */
	 tx664 = DeDPhi*tx663; 		/* rule 802 */
	 tx665 = tx100*tx139*tx142*tx93*V; 		/* rule 803 */
	 ohy1z2 = tx664 + tx665; 		/* rule 804 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 tx666 = -(DeDPhi*tx649); 		/* rule 806 */
	 tx667 = tx104*tx139*tx142*tx93*V; 		/* rule 807 */
	 ohy1x3 = tx666 + tx667; 		/* rule 808 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 0, ohy1x3);
	 tx668 = -(DeDPhi*tx654); 		/* rule 810 */
	 tx669 = tx108*tx139*tx142*tx93*V; 		/* rule 811 */
	 ohy1y3 = tx668 + tx669; 		/* rule 812 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 1, ohy1y3);
	 tx670 = -(DeDPhi*tx662); 		/* rule 814 */
	 tx671 = tx112*tx139*tx142*tx93*V; 		/* rule 815 */
	 ohy1z3 = tx670 + tx671; 		/* rule 816 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 2, ohy1z3);
	 ohy1x4 = tx139*tx142*tx20*tx39*tx72*tx85*tx907*V; 		/* rule 818 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 0, ohy1x4);
	 ohy1y4 = tx139*tx142*tx45*tx72*tx85*tx907*V; 		/* rule 820 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 1, ohy1y4);
	 ohy1z4 = tx139*tx142*tx20*tx40*tx72*tx85*tx907*V; 		/* rule 822 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 2, ohy1z4);
	 tx672 = tx21*tx427*tx907; 		/* rule 824 */
	 tx673 = tx19*tx234*tx88; 		/* rule 825 */
	 tx674 = tx672 + tx673; 		/* rule 826 */
	 tx675 = tx141*tx674*tx908; 		/* rule 827 */
	 tx676 = tx630 + tx675; 		/* rule 828 */
	 tx677 = DeDPhi*tx676; 		/* rule 829 */
	 tx678 = tx139*tx142*tx92*tx97*V; 		/* rule 830 */
	 ohz1x2 = tx677 + tx678; 		/* rule 831 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 tx679 = tx20*tx427*tx907; 		/* rule 833 */
	 tx680 = tx19*tx332*tx88; 		/* rule 834 */
	 tx681 = tx679 + tx680; 		/* rule 835 */
	 tx682 = tx141*tx681*tx908; 		/* rule 836 */
	 tx683 = tx661 + tx682; 		/* rule 837 */
	 tx684 = DeDPhi*tx683; 		/* rule 838 */
	 tx685 = tx139*tx142*tx96*tx97*V; 		/* rule 839 */
	 ohz1y2 = tx684 + tx685; 		/* rule 840 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx686 = -2.*tx141*tx19*tx194*tx65; 		/* rule 842 */
	 tx687 = tx141*tx473*tx908; 		/* rule 843 */
	 tx688 = tx686 + tx687; 		/* rule 844 */
	 tx689 = DeDPhi*tx688; 		/* rule 845 */
	 tx690 = tx100*tx139*tx142*tx97*V; 		/* rule 846 */
	 ohz1z2 = tx689 + tx690; 		/* rule 847 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 tx691 = -(DeDPhi*tx675); 		/* rule 849 */
	 tx692 = tx104*tx139*tx142*tx97*V; 		/* rule 850 */
	 ohz1x3 = tx691 + tx692; 		/* rule 851 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 0, ohz1x3);
	 tx693 = -(DeDPhi*tx682); 		/* rule 853 */
	 tx694 = tx108*tx139*tx142*tx97*V; 		/* rule 854 */
	 ohz1y3 = tx693 + tx694; 		/* rule 855 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 1, ohz1y3);
	 tx695 = -(DeDPhi*tx687); 		/* rule 857 */
	 tx696 = tx112*tx139*tx142*tx97*V; 		/* rule 858 */
	 ohz1z3 = tx695 + tx696; 		/* rule 859 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 2, ohz1z3);
	 ohz1x4 = tx139*tx142*tx19*tx39*tx72*tx85*tx907*V; 		/* rule 861 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 0, ohz1x4);
	 ohz1y4 = tx139*tx142*tx19*tx38*tx72*tx85*tx907*V; 		/* rule 863 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 1, ohz1y4);
	 ohz1z4 = tx139*tx142*tx47*tx72*tx85*tx907*V; 		/* rule 865 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 2, ohz1z4);
	 tx697 = tx20*tx282; 		/* rule 867 */
	 tx698 = tx283*tx38; 		/* rule 868 */
	 tx699 = tx284*tx38; 		/* rule 869 */
	 tx700 = tx21*tx380; 		/* rule 870 */
	 tx701 = tx381*tx39; 		/* rule 871 */
	 tx702 = tx382*tx39; 		/* rule 872 */
	 tx703 = tx697 + tx700; 		/* rule 873 */
	 tx704 = tx698 + tx701; 		/* rule 874 */
	 tx705 = tx699 + tx702; 		/* rule 875 */
	 tx706 = -tx623; 		/* rule 876 */
	 tx707 = -tx649; 		/* rule 877 */
	 tx708 = tx141*tx606*tx65; 		/* rule 878 */
	 tx709 = 0.5*tx286*tx606*tx72; 		/* rule 879 */
	 tx710 = -0.5*tx286*tx704*tx85; 		/* rule 880 */
	 tx711 = -(tx285*tx705*tx87*tx908); 		/* rule 881 */
	 tx712 = tx141*tx703*tx88*tx908; 		/* rule 882 */
	 tx713 = tx706 + tx707 + tx708 + tx709 + tx710 + tx711 + tx712; 		/* rule 883 */
	 tx714 = DeDPhi*tx713; 		/* rule 884 */
	 tx715 = -(tx139*tx142*tx92*tx96*V); 		/* rule 885 */
	 ohx2y2 = tx714 + tx715; 		/* rule 886 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 tx716 = tx19*tx282; 		/* rule 888 */
	 tx717 = tx283*tx40; 		/* rule 889 */
	 tx718 = tx284*tx40; 		/* rule 890 */
	 tx719 = tx21*tx474; 		/* rule 891 */
	 tx720 = tx39*tx475; 		/* rule 892 */
	 tx721 = tx39*tx476; 		/* rule 893 */
	 tx722 = tx716 + tx719; 		/* rule 894 */
	 tx723 = tx717 + tx720; 		/* rule 895 */
	 tx724 = tx718 + tx721; 		/* rule 896 */
	 tx725 = -tx631; 		/* rule 897 */
	 tx726 = tx141*tx611*tx65; 		/* rule 898 */
	 tx727 = -tx675; 		/* rule 899 */
	 tx728 = 0.5*tx286*tx611*tx72; 		/* rule 900 */
	 tx729 = -0.5*tx286*tx723*tx85; 		/* rule 901 */
	 tx730 = -(tx285*tx724*tx87*tx908); 		/* rule 902 */
	 tx731 = tx141*tx722*tx88*tx908; 		/* rule 903 */
	 tx732 = tx725 + tx726 + tx727 + tx728 + tx729 + tx730 + tx731; 		/* rule 904 */
	 tx733 = DeDPhi*tx732; 		/* rule 905 */
	 tx734 = -(tx100*tx139*tx142*tx92*V); 		/* rule 906 */
	 ohx2z2 = tx733 + tx734; 		/* rule 907 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 tx735 = -tx288; 		/* rule 909 */
	 tx736 = tx283*tx286*tx39*tx85; 		/* rule 910 */
	 tx737 = tx285*tx517*tx908; 		/* rule 911 */
	 tx738 = 2.*tx284*tx285*tx39*tx87*tx908; 		/* rule 912 */
	 tx739 = -2.*tx141*tx21*tx282*tx88*tx908; 		/* rule 913 */
	 tx740 = tx615 + tx735 + tx736 + tx737 + tx738 + tx739; 		/* rule 914 */
	 tx741 = DeDPhi*tx740; 		/* rule 915 */
	 tx742 = -(tx104*tx139*tx142*tx92*V); 		/* rule 916 */
	 ohx2x3 = tx741 + tx742; 		/* rule 917 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 0, ohx2x3);
	 tx743 = tx39*tx550*tx907; 		/* rule 919 */
	 tx744 = tx38*tx514*tx87; 		/* rule 920 */
	 tx745 = tx743 + tx744; 		/* rule 921 */
	 tx746 = -tx712; 		/* rule 922 */
	 tx747 = -0.5*tx286*tx606*tx72; 		/* rule 923 */
	 tx748 = 0.5*tx286*tx704*tx85; 		/* rule 924 */
	 tx749 = tx285*tx745*tx908; 		/* rule 925 */
	 tx750 = tx285*tx705*tx87*tx908; 		/* rule 926 */
	 tx751 = tx623 + tx746 + tx747 + tx748 + tx749 + tx750; 		/* rule 927 */
	 tx752 = DeDPhi*tx751; 		/* rule 928 */
	 tx753 = -(tx108*tx139*tx142*tx92*V); 		/* rule 929 */
	 ohx2y3 = tx752 + tx753; 		/* rule 930 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 1, ohx2y3);
	 tx754 = tx39*tx587*tx907; 		/* rule 932 */
	 tx755 = tx40*tx514*tx87; 		/* rule 933 */
	 tx756 = tx754 + tx755; 		/* rule 934 */
	 tx757 = -0.5*tx286*tx611*tx72; 		/* rule 935 */
	 tx758 = -tx731; 		/* rule 936 */
	 tx759 = 0.5*tx286*tx723*tx85; 		/* rule 937 */
	 tx760 = tx285*tx756*tx908; 		/* rule 938 */
	 tx761 = tx285*tx724*tx87*tx908; 		/* rule 939 */
	 tx762 = tx631 + tx757 + tx758 + tx759 + tx760 + tx761; 		/* rule 940 */
	 tx763 = DeDPhi*tx762; 		/* rule 941 */
	 tx764 = -(tx112*tx139*tx142*tx92*V); 		/* rule 942 */
	 ohx2z3 = tx763 + tx764; 		/* rule 943 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 2, ohx2z3);
	 tx765 = -(DeDPhi*tx737); 		/* rule 945 */
	 tx766 = tx101*tx139*tx142*tx92*V; 		/* rule 946 */
	 ohx2x4 = tx765 + tx766; 		/* rule 947 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 0, ohx2x4);
	 tx767 = -(DeDPhi*tx749); 		/* rule 949 */
	 tx768 = tx105*tx139*tx142*tx92*V; 		/* rule 950 */
	 ohx2y4 = tx767 + tx768; 		/* rule 951 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 1, ohx2y4);
	 tx769 = -(DeDPhi*tx760); 		/* rule 953 */
	 tx770 = tx109*tx139*tx142*tx92*V; 		/* rule 954 */
	 ohx2z4 = tx769 + tx770; 		/* rule 955 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 2, ohx2z4);
	 tx771 = tx19*tx380; 		/* rule 957 */
	 tx772 = tx381*tx40; 		/* rule 958 */
	 tx773 = tx382*tx40; 		/* rule 959 */
	 tx774 = tx20*tx474; 		/* rule 960 */
	 tx775 = tx38*tx475; 		/* rule 961 */
	 tx776 = tx38*tx476; 		/* rule 962 */
	 tx777 = tx771 + tx774; 		/* rule 963 */
	 tx778 = tx772 + tx775; 		/* rule 964 */
	 tx779 = tx773 + tx776; 		/* rule 965 */
	 tx780 = tx141*tx643*tx65; 		/* rule 966 */
	 tx781 = -tx662; 		/* rule 967 */
	 tx782 = -tx682; 		/* rule 968 */
	 tx783 = 0.5*tx286*tx643*tx72; 		/* rule 969 */
	 tx784 = -0.5*tx286*tx778*tx85; 		/* rule 970 */
	 tx785 = -(tx285*tx779*tx87*tx908); 		/* rule 971 */
	 tx786 = tx141*tx777*tx88*tx908; 		/* rule 972 */
	 tx787 = tx780 + tx781 + tx782 + tx783 + tx784 + tx785 + tx786; 		/* rule 973 */
	 tx788 = DeDPhi*tx787; 		/* rule 974 */
	 tx789 = -(tx100*tx139*tx142*tx96*V); 		/* rule 975 */
	 ohy2z2 = tx788 + tx789; 		/* rule 976 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	 tx790 = tx38*tx513*tx907; 		/* rule 978 */
	 tx791 = tx39*tx551*tx87; 		/* rule 979 */
	 tx792 = tx790 + tx791; 		/* rule 980 */
	 tx793 = tx285*tx792*tx908; 		/* rule 981 */
	 tx794 = tx649 + tx746 + tx747 + tx748 + tx750 + tx793; 		/* rule 982 */
	 tx795 = DeDPhi*tx794; 		/* rule 983 */
	 tx796 = -(tx104*tx139*tx142*tx96*V); 		/* rule 984 */
	 ohy2x3 = tx795 + tx796; 		/* rule 985 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 0, ohy2x3);
	 tx797 = -tx384; 		/* rule 987 */
	 tx798 = tx286*tx38*tx381*tx85; 		/* rule 988 */
	 tx799 = tx285*tx554*tx908; 		/* rule 989 */
	 tx800 = 2.*tx285*tx38*tx382*tx87*tx908; 		/* rule 990 */
	 tx801 = -2.*tx141*tx20*tx380*tx88*tx908; 		/* rule 991 */
	 tx802 = tx654 + tx797 + tx798 + tx799 + tx800 + tx801; 		/* rule 992 */
	 tx803 = DeDPhi*tx802; 		/* rule 993 */
	 tx804 = -(tx108*tx139*tx142*tx96*V); 		/* rule 994 */
	 ohy2y3 = tx803 + tx804; 		/* rule 995 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 1, ohy2y3);
	 tx805 = tx38*tx587*tx907; 		/* rule 997 */
	 tx806 = tx40*tx551*tx87; 		/* rule 998 */
	 tx807 = tx805 + tx806; 		/* rule 999 */
	 tx808 = -0.5*tx286*tx643*tx72; 		/* rule 1000 */
	 tx809 = -tx786; 		/* rule 1001 */
	 tx810 = 0.5*tx286*tx778*tx85; 		/* rule 1002 */
	 tx811 = tx285*tx807*tx908; 		/* rule 1003 */
	 tx812 = tx285*tx779*tx87*tx908; 		/* rule 1004 */
	 tx813 = tx662 + tx808 + tx809 + tx810 + tx811 + tx812; 		/* rule 1005 */
	 tx814 = DeDPhi*tx813; 		/* rule 1006 */
	 tx815 = -(tx112*tx139*tx142*tx96*V); 		/* rule 1007 */
	 ohy2z3 = tx814 + tx815; 		/* rule 1008 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 2, ohy2z3);
	 tx816 = -(DeDPhi*tx793); 		/* rule 1010 */
	 tx817 = tx101*tx139*tx142*tx96*V; 		/* rule 1011 */
	 ohy2x4 = tx816 + tx817; 		/* rule 1012 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 0, ohy2x4);
	 tx818 = -(DeDPhi*tx799); 		/* rule 1014 */
	 tx819 = tx105*tx139*tx142*tx96*V; 		/* rule 1015 */
	 ohy2y4 = tx818 + tx819; 		/* rule 1016 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 1, ohy2y4);
	 tx820 = -(DeDPhi*tx811); 		/* rule 1018 */
	 tx821 = tx109*tx139*tx142*tx96*V; 		/* rule 1019 */
	 ohy2z4 = tx820 + tx821; 		/* rule 1020 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 2, ohy2z4);
	 tx822 = tx40*tx513*tx907; 		/* rule 1022 */
	 tx823 = tx39*tx588*tx87; 		/* rule 1023 */
	 tx824 = tx822 + tx823; 		/* rule 1024 */
	 tx825 = tx285*tx824*tx908; 		/* rule 1025 */
	 tx826 = tx675 + tx757 + tx758 + tx759 + tx761 + tx825; 		/* rule 1026 */
	 tx827 = DeDPhi*tx826; 		/* rule 1027 */
	 tx828 = -(tx100*tx104*tx139*tx142*V); 		/* rule 1028 */
	 ohz2x3 = tx827 + tx828; 		/* rule 1029 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 0, ohz2x3);
	 tx829 = tx40*tx550*tx907; 		/* rule 1031 */
	 tx830 = tx38*tx588*tx87; 		/* rule 1032 */
	 tx831 = tx829 + tx830; 		/* rule 1033 */
	 tx832 = tx285*tx831*tx908; 		/* rule 1034 */
	 tx833 = tx682 + tx808 + tx809 + tx810 + tx812 + tx832; 		/* rule 1035 */
	 tx834 = DeDPhi*tx833; 		/* rule 1036 */
	 tx835 = -(tx100*tx108*tx139*tx142*V); 		/* rule 1037 */
	 ohz2y3 = tx834 + tx835; 		/* rule 1038 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 1, ohz2y3);
	 tx836 = -tx478; 		/* rule 1040 */
	 tx837 = tx286*tx40*tx475*tx85; 		/* rule 1041 */
	 tx838 = tx285*tx591*tx908; 		/* rule 1042 */
	 tx839 = 2.*tx285*tx40*tx476*tx87*tx908; 		/* rule 1043 */
	 tx840 = -2.*tx141*tx19*tx474*tx88*tx908; 		/* rule 1044 */
	 tx841 = tx687 + tx836 + tx837 + tx838 + tx839 + tx840; 		/* rule 1045 */
	 tx842 = DeDPhi*tx841; 		/* rule 1046 */
	 tx843 = -(tx100*tx112*tx139*tx142*V); 		/* rule 1047 */
	 ohz2z3 = tx842 + tx843; 		/* rule 1048 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 2, ohz2z3);
	 tx844 = -(DeDPhi*tx825); 		/* rule 1050 */
	 tx845 = tx100*tx101*tx139*tx142*V; 		/* rule 1051 */
	 ohz2x4 = tx844 + tx845; 		/* rule 1052 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 0, ohz2x4);
	 tx846 = -(DeDPhi*tx832); 		/* rule 1054 */
	 tx847 = tx100*tx105*tx139*tx142*V; 		/* rule 1055 */
	 ohz2y4 = tx846 + tx847; 		/* rule 1056 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 1, ohz2y4);
	 tx848 = -(DeDPhi*tx838); 		/* rule 1058 */
	 tx849 = tx100*tx109*tx139*tx142*V; 		/* rule 1059 */
	 ohz2z4 = tx848 + tx849; 		/* rule 1060 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 2, ohz2z4);
	 tx850 = -(tx285*tx65*tx704); 		/* rule 1062 */
	 tx851 = -tx749; 		/* rule 1063 */
	 tx852 = -tx793; 		/* rule 1064 */
	 tx853 = tx709 + tx710 + tx711 + tx712 + tx850 + tx851 + tx852; 		/* rule 1065 */
	 tx854 = DeDPhi*tx853; 		/* rule 1066 */
	 tx855 = -(tx104*tx108*tx139*tx142*V); 		/* rule 1067 */
	 ohx3y3 = tx854 + tx855; 		/* rule 1068 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 1, ohx3y3);
	 tx856 = -(tx285*tx65*tx723); 		/* rule 1070 */
	 tx857 = -tx760; 		/* rule 1071 */
	 tx858 = -tx825; 		/* rule 1072 */
	 tx859 = tx728 + tx729 + tx730 + tx731 + tx856 + tx857 + tx858; 		/* rule 1073 */
	 tx860 = DeDPhi*tx859; 		/* rule 1074 */
	 tx861 = -(tx104*tx112*tx139*tx142*V); 		/* rule 1075 */
	 ohx3z3 = tx860 + tx861; 		/* rule 1076 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 2, ohx3z3);
	 tx862 = 2.*tx283*tx285*tx39*tx65; 		/* rule 1078 */
	 tx863 = tx737 + tx862; 		/* rule 1079 */
	 tx864 = DeDPhi*tx863; 		/* rule 1080 */
	 tx865 = tx101*tx104*tx139*tx142*V; 		/* rule 1081 */
	 ohx3x4 = tx864 + tx865; 		/* rule 1082 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 0, ohx3x4);
	 tx866 = tx285*tx65*tx704; 		/* rule 1084 */
	 tx867 = tx749 + tx866; 		/* rule 1085 */
	 tx868 = DeDPhi*tx867; 		/* rule 1086 */
	 tx869 = tx104*tx105*tx139*tx142*V; 		/* rule 1087 */
	 ohx3y4 = tx868 + tx869; 		/* rule 1088 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 1, ohx3y4);
	 tx870 = tx285*tx65*tx723; 		/* rule 1090 */
	 tx871 = tx760 + tx870; 		/* rule 1091 */
	 tx872 = DeDPhi*tx871; 		/* rule 1092 */
	 tx873 = tx104*tx109*tx139*tx142*V; 		/* rule 1093 */
	 ohx3z4 = tx872 + tx873; 		/* rule 1094 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 2, ohx3z4);
	 tx874 = -(tx285*tx65*tx778); 		/* rule 1096 */
	 tx875 = -tx811; 		/* rule 1097 */
	 tx876 = -tx832; 		/* rule 1098 */
	 tx877 = tx783 + tx784 + tx785 + tx786 + tx874 + tx875 + tx876; 		/* rule 1099 */
	 tx878 = DeDPhi*tx877; 		/* rule 1100 */
	 tx879 = -(tx108*tx112*tx139*tx142*V); 		/* rule 1101 */
	 ohy3z3 = tx878 + tx879; 		/* rule 1102 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 2, ohy3z3);
	 tx880 = tx793 + tx866; 		/* rule 1104 */
	 tx881 = DeDPhi*tx880; 		/* rule 1105 */
	 tx882 = tx101*tx108*tx139*tx142*V; 		/* rule 1106 */
	 ohy3x4 = tx881 + tx882; 		/* rule 1107 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 0, ohy3x4);
	 tx883 = 2.*tx285*tx38*tx381*tx65; 		/* rule 1109 */
	 tx884 = tx799 + tx883; 		/* rule 1110 */
	 tx885 = DeDPhi*tx884; 		/* rule 1111 */
	 tx886 = tx105*tx108*tx139*tx142*V; 		/* rule 1112 */
	 ohy3y4 = tx885 + tx886; 		/* rule 1113 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 1, ohy3y4);
	 tx887 = tx285*tx65*tx778; 		/* rule 1115 */
	 tx888 = tx811 + tx887; 		/* rule 1116 */
	 tx889 = DeDPhi*tx888; 		/* rule 1117 */
	 tx890 = tx108*tx109*tx139*tx142*V; 		/* rule 1118 */
	 ohy3z4 = tx889 + tx890; 		/* rule 1119 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 2, ohy3z4);
	 tx891 = tx825 + tx870; 		/* rule 1121 */
	 tx892 = DeDPhi*tx891; 		/* rule 1122 */
	 tx893 = tx101*tx112*tx139*tx142*V; 		/* rule 1123 */
	 ohz3x4 = tx892 + tx893; 		/* rule 1124 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 0, ohz3x4);
	 tx894 = tx832 + tx887; 		/* rule 1126 */
	 tx895 = DeDPhi*tx894; 		/* rule 1127 */
	 tx896 = tx105*tx112*tx139*tx142*V; 		/* rule 1128 */
	 ohz3y4 = tx895 + tx896; 		/* rule 1129 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 1, ohz3y4);
	 tx897 = 2.*tx285*tx40*tx475*tx65; 		/* rule 1131 */
	 tx898 = tx838 + tx897; 		/* rule 1132 */
	 tx899 = DeDPhi*tx898; 		/* rule 1133 */
	 tx900 = tx109*tx112*tx139*tx142*V; 		/* rule 1134 */
	 ohz3z4 = tx899 + tx900; 		/* rule 1135 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 2, ohz3z4);
	 tx901 = DeDPhi*tx850; 		/* rule 1137 */
	 tx902 = -(tx139*tx142*tx285*tx38*tx39*tx907*V); 		/* rule 1138 */
	 ohx4y4 = tx901 + tx902; 		/* rule 1139 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 1, ohx4y4);
	 tx903 = DeDPhi*tx856; 		/* rule 1141 */
	 tx904 = -(tx139*tx142*tx285*tx39*tx40*tx907*V); 		/* rule 1142 */
	 ohx4z4 = tx903 + tx904; 		/* rule 1143 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 2, ohx4z4);
	 tx905 = DeDPhi*tx874; 		/* rule 1145 */
	 tx906 = -(tx139*tx142*tx285*tx38*tx40*tx907*V); 		/* rule 1146 */
	 ohy4z4 = tx905 + tx906; 		/* rule 1147 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 1, I4, 2, ohy4z4);
	} /*calcOffDiagonalHessian*/
	#endif // DIHEDRAL_CALC_OFF_DIAGONAL_HESSIAN ]
	} /*calcDiagonalHessian*/
	#endif // DIHEDRAL_CALC_DIAGONAL_HESSIAN ]
	} /*calcForce*/
	#endif // DIHEDRAL_CALC_FORCE ]
	DIHEDRAL_DEBUG_INTERACTIONS(I1,I2,I3,I4);
	SKIP_term: (void)0;
