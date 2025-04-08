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
	 tzz502 = tx27*tx29; 		/* rule 59 */
	 tx31 = -(tx28*tzz502); 		/* rule 60 */
	 tx32 = 1 + tx31; 		/* rule 61 */
	 tzz475 = tx23*tx30; 		/* rule 62 */
	 tzz528 = tx22*tzz475; 		/* rule 63 */
	 tzz532 = -tzz528; 		/* rule 64 */
	 tx33 = tx4*tzz532; 		/* rule 65 */
	 tzz529 = tx22*tx453; 		/* rule 66 */
	 tx34 = tx5*tzz529; 		/* rule 67 */
	 tx455 = mysqrt(tx32); 		/* rule 68 */
	 tx35 = reciprocal(tx455); 		/* rule 69 */
	 tx36 = tx33 + tx34; 		/* rule 70 */
	 tzz466 = AngleDeviation*kt; 		/* rule 71 */
	 tzz468 = -2.*tzz466; 		/* rule 72 */
	 tzz488 = tx35*tzz468; 		/* rule 73 */
	 gx1 = tx36*tzz488; 		/* rule 74 */
	 fx1 = -gx1; 		/* rule 75 */
	ANGLE_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx37 = tx6*tzz532; 		/* rule 77 */
	 tx38 = tx7*tzz529; 		/* rule 78 */
	 tx39 = tx37 + tx38; 		/* rule 79 */
	 gy1 = tx39*tzz488; 		/* rule 80 */
	 fy1 = -gy1; 		/* rule 81 */
	ANGLE_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx40 = tx8*tzz532; 		/* rule 83 */
	 tx41 = tx9*tzz529; 		/* rule 84 */
	 tx42 = tx40 + tx41; 		/* rule 85 */
	 gz1 = tx42*tzz488; 		/* rule 86 */
	 fz1 = -gz1; 		/* rule 87 */
	ANGLE_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tx43 = -x1; 		/* rule 89 */
	 tx44 = 2.*x2; 		/* rule 90 */
	 tx45 = -x3; 		/* rule 91 */
	 tx456 = tx17; 		/* rule 92 */
	 tx457 = tx22; 		/* rule 93 */
	 tx458 = reciprocal(tx456); 		/* rule 94 */
	 tx46 = tx457*tx458; 		/* rule 95 */
	 tx47 = tx43 + tx44 + tx45; 		/* rule 96 */
	 tzz469 = tx30*tx457; 		/* rule 97 */
	 tzz496 = tx4*tzz469; 		/* rule 98 */
	 tx48 = tx23*tzz496; 		/* rule 99 */
	 tzz495 = tx453*tx457; 		/* rule 100 */
	 tx49 = tx47*tzz495; 		/* rule 101 */
	 tzz467 = tx453*tx46; 		/* rule 102 */
	 tzz494 = tx5*tzz467; 		/* rule 103 */
	 tx50 = tx23*tzz494; 		/* rule 104 */
	 tx51 = tx48 + tx49 + tx50; 		/* rule 105 */
	 gx2 = tx51*tzz488; 		/* rule 106 */
	 fx2 = -gx2; 		/* rule 107 */
	ANGLE_FORCE_ACCUMULATE(I2, 0, fx2 );
	 tx52 = -y1; 		/* rule 109 */
	 tx53 = 2.*y2; 		/* rule 110 */
	 tx54 = -y3; 		/* rule 111 */
	 tx55 = tx52 + tx53 + tx54; 		/* rule 112 */
	 tx56 = tx55*tzz495; 		/* rule 113 */
	 tzz486 = tx6*tzz469; 		/* rule 114 */
	 tx57 = tx23*tzz486; 		/* rule 115 */
	 tzz485 = tx7*tzz467; 		/* rule 116 */
	 tx58 = tx23*tzz485; 		/* rule 117 */
	 tx59 = tx56 + tx57 + tx58; 		/* rule 118 */
	 gy2 = tx59*tzz488; 		/* rule 119 */
	 fy2 = -gy2; 		/* rule 120 */
	ANGLE_FORCE_ACCUMULATE(I2, 1, fy2 );
	 tx60 = -z1; 		/* rule 122 */
	 tx61 = 2.*z2; 		/* rule 123 */
	 tx62 = -z3; 		/* rule 124 */
	 tx63 = tx60 + tx61 + tx62; 		/* rule 125 */
	 tx64 = tx63*tzz495; 		/* rule 126 */
	 tzz484 = tx8*tzz469; 		/* rule 127 */
	 tx65 = tx23*tzz484; 		/* rule 128 */
	 tzz483 = tx9*tzz467; 		/* rule 129 */
	 tx66 = tx23*tzz483; 		/* rule 130 */
	 tx67 = tx64 + tx65 + tx66; 		/* rule 131 */
	 gz2 = tx67*tzz488; 		/* rule 132 */
	 fz2 = -gz2; 		/* rule 133 */
	ANGLE_FORCE_ACCUMULATE(I2, 2, fz2 );
	 tx68 = tx4*tzz495; 		/* rule 135 */
	 tzz474 = -tzz467; 		/* rule 136 */
	 tzz506 = tx5*tzz474; 		/* rule 137 */
	 tx69 = tx23*tzz506; 		/* rule 138 */
	 tx70 = tx68 + tx69; 		/* rule 139 */
	 gx3 = tx70*tzz488; 		/* rule 140 */
	 fx3 = -gx3; 		/* rule 141 */
	ANGLE_FORCE_ACCUMULATE(I3, 0, fx3 );
	 tx71 = tx6*tzz495; 		/* rule 143 */
	 tzz500 = tx7*tzz474; 		/* rule 144 */
	 tx72 = tx23*tzz500; 		/* rule 145 */
	 tx73 = tx71 + tx72; 		/* rule 146 */
	 gy3 = tx73*tzz488; 		/* rule 147 */
	 fy3 = -gy3; 		/* rule 148 */
	ANGLE_FORCE_ACCUMULATE(I3, 1, fy3 );
	 tx74 = tx8*tzz495; 		/* rule 150 */
	 tzz492 = tx9*tzz474; 		/* rule 151 */
	 tx75 = tx23*tzz492; 		/* rule 152 */
	 tx76 = tx74 + tx75; 		/* rule 153 */
	 gz3 = tx76*tzz488; 		/* rule 154 */
	 fz3 = -gz3; 		/* rule 155 */
	ANGLE_FORCE_ACCUMULATE(I3, 2, fz3 );
	#ifdef ANGLE_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx459 = power2(tx452); 		/* rule 159 */
	 tx460 = reciprocal(tx459); 		/* rule 160 */
	 tx77 = tx453*tx460; 		/* rule 161 */
	 tx78 = power2(tx27); 		/* rule 162 */
	 tzz476 = -tzz469; 		/* rule 163 */
	 tx79 = tx23*tzz476; 		/* rule 164 */
	 tzz482 = -2.*tx28; 		/* rule 165 */
	 tzz489 = tx27*tzz482; 		/* rule 166 */
	 tzz490 = tx23*tzz489; 		/* rule 167 */
	 tx80 = tx5*tzz490; 		/* rule 168 */
	 tx81 = -2.*tx5*tzz496; 		/* rule 169 */
	 tzz478 = tx457*tx77; 		/* rule 170 */
	 tzz479 = tx23*tzz478; 		/* rule 171 */
	 tzz510 = 3*tzz479; 		/* rule 172 */
	 tx82 = tx10*tzz510; 		/* rule 173 */
	 tzz507 = tx29*tx78; 		/* rule 174 */
	 tzz527 = tx28*tzz507; 		/* rule 175 */
	 tzz530 = 2.*tzz527; 		/* rule 176 */
	 tx83 = tx4*tzz530; 		/* rule 177 */
	 tx461 = tx32; 		/* rule 178 */
	 tx462 = tx35; 		/* rule 179 */
	 tx463 = reciprocal(tx461); 		/* rule 180 */
	 tx84 = tx462*tx463; 		/* rule 181 */
	 tx85 = tx463; 		/* rule 182 */
	 tx86 = power2(tx36); 		/* rule 183 */
	 tx87 = tx79 + tx81 + tx82; 		/* rule 184 */
	 tx88 = tx80 + tx83; 		/* rule 185 */
	 tzz472 = kt*tx85; 		/* rule 186 */
	 tzz473 = 2.*tzz472; 		/* rule 187 */
	 tx89 = tx86*tzz473; 		/* rule 188 */
	 tzz471 = tx462*tzz468; 		/* rule 189 */
	 tx90 = tx87*tzz471; 		/* rule 190 */
	 tzz470 = tx84*tzz466; 		/* rule 191 */
	 tzz487 = tx36*tzz470; 		/* rule 192 */
	 tx91 = tx88*tzz487; 		/* rule 193 */
	 dhx1x1 = tx89 + tx90 + tx91; 		/* rule 194 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx92 = tx7*tzz490; 		/* rule 196 */
	 tx93 = -2.*tx7*tzz486; 		/* rule 197 */
	 tx94 = tx12*tzz510; 		/* rule 198 */
	 tx95 = tx6*tzz530; 		/* rule 199 */
	 tx96 = power2(tx39); 		/* rule 200 */
	 tx97 = tx79 + tx93 + tx94; 		/* rule 201 */
	 tx98 = tx92 + tx95; 		/* rule 202 */
	 tx99 = tx96*tzz473; 		/* rule 203 */
	 tx100 = tx97*tzz471; 		/* rule 204 */
	 tzz497 = tx39*tzz470; 		/* rule 205 */
	 tx101 = tx98*tzz497; 		/* rule 206 */
	 dhy1y1 = tx100 + tx101 + tx99; 		/* rule 207 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx102 = tx14*tzz510; 		/* rule 209 */
	 tx103 = tx8*tzz530; 		/* rule 210 */
	 tx104 = tx9*tzz490; 		/* rule 211 */
	 tx105 = -2.*tx9*tzz484; 		/* rule 212 */
	 tx106 = tx103 + tx104; 		/* rule 213 */
	 tx107 = power2(tx42); 		/* rule 214 */
	 tx108 = tx102 + tx105 + tx79; 		/* rule 215 */
	 tx109 = tx108*tzz471; 		/* rule 216 */
	 tzz501 = tx42*tzz470; 		/* rule 217 */
	 tx110 = tx106*tzz501; 		/* rule 218 */
	 tx111 = tx107*tzz473; 		/* rule 219 */
	 dhz1z1 = tx109 + tx110 + tx111; 		/* rule 220 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 tx464 = power2(tx456); 		/* rule 222 */
	 tx465 = reciprocal(tx464); 		/* rule 223 */
	 tx112 = tx457*tx465; 		/* rule 224 */
	 tx113 = power2(tx28); 		/* rule 225 */
	 tx114 = 2.*tzz495; 		/* rule 226 */
	 tzz480 = tx23*tx453; 		/* rule 227 */
	 tzz481 = tx112*tzz480; 		/* rule 228 */
	 tzz509 = 3*tzz481; 		/* rule 229 */
	 tx115 = tx11*tzz509; 		/* rule 230 */
	 tx116 = tx23*tzz474; 		/* rule 231 */
	 tx117 = tx47*tzz490; 		/* rule 232 */
	 tzz544 = 2.*tx47; 		/* rule 233 */
	 tx118 = tzz496*tzz544; 		/* rule 234 */
	 tzz508 = tx113*tzz502; 		/* rule 235 */
	 tzz533 = -2.*tzz508; 		/* rule 236 */
	 tx119 = tx5*tzz533; 		/* rule 237 */
	 tzz477 = tx46*tzz475; 		/* rule 238 */
	 tzz520 = tx5*tzz477; 		/* rule 239 */
	 tx120 = 2.*tx4*tzz520; 		/* rule 240 */
	 tx121 = tzz494*tzz544; 		/* rule 241 */
	 tzz521 = tzz482*tzz507; 		/* rule 242 */
	 tx122 = tx4*tzz521; 		/* rule 243 */
	 tx123 = tx117 + tx119 + tx122; 		/* rule 244 */
	 tx124 = power2(tx51); 		/* rule 245 */
	 tzz547 = tx116 + tx79; 		/* rule 246 */
	 tzz548 = tx114 + tzz547; 		/* rule 247 */
	 tx125 = tx115 + tx118 + tx120 + tx121 + tx82 + tzz548; 		/* rule 248 */
	 tx126 = tx125*tzz471; 		/* rule 249 */
	 tzz505 = tx51*tzz470; 		/* rule 250 */
	 tx127 = tx123*tzz505; 		/* rule 251 */
	 tx128 = tx124*tzz473; 		/* rule 252 */
	 dhx2x2 = tx126 + tx127 + tx128; 		/* rule 253 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 tx129 = tx13*tzz509; 		/* rule 255 */
	 tx130 = tx55*tzz490; 		/* rule 256 */
	 tzz543 = 2.*tx55; 		/* rule 257 */
	 tx131 = tzz486*tzz543; 		/* rule 258 */
	 tx132 = tx7*tzz533; 		/* rule 259 */
	 tx133 = tzz485*tzz543; 		/* rule 260 */
	 tzz517 = tx7*tzz477; 		/* rule 261 */
	 tx134 = 2.*tx6*tzz517; 		/* rule 262 */
	 tx135 = tx6*tzz521; 		/* rule 263 */
	 tx136 = tx130 + tx132 + tx135; 		/* rule 264 */
	 tx137 = power2(tx59); 		/* rule 265 */
	 tx138 = tx129 + tx131 + tx133 + tx134 + tx94 + tzz548; 		/* rule 266 */
	 tx139 = tx138*tzz471; 		/* rule 267 */
	 tzz514 = tx59*tzz470; 		/* rule 268 */
	 tx140 = tx136*tzz514; 		/* rule 269 */
	 tx141 = tx137*tzz473; 		/* rule 270 */
	 dhy2y2 = tx139 + tx140 + tx141; 		/* rule 271 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 tx142 = tx15*tzz509; 		/* rule 273 */
	 tx143 = tx63*tzz490; 		/* rule 274 */
	 tzz542 = 2.*tx63; 		/* rule 275 */
	 tx144 = tzz484*tzz542; 		/* rule 276 */
	 tx145 = tx8*tzz521; 		/* rule 277 */
	 tx146 = tx9*tzz533; 		/* rule 278 */
	 tx147 = tzz483*tzz542; 		/* rule 279 */
	 tzz516 = tx9*tzz477; 		/* rule 280 */
	 tx148 = 2.*tx8*tzz516; 		/* rule 281 */
	 tx149 = tx143 + tx145 + tx146; 		/* rule 282 */
	 tx150 = power2(tx67); 		/* rule 283 */
	 tx151 = tx102 + tx142 + tx144 + tx147 + tx148 + tzz548; 		/* rule 284 */
	 tx152 = tx151*tzz471; 		/* rule 285 */
	 tzz518 = tx67*tzz470; 		/* rule 286 */
	 tx153 = tx149*tzz518; 		/* rule 287 */
	 tx154 = tx150*tzz473; 		/* rule 288 */
	 dhz2z2 = tx152 + tx153 + tx154; 		/* rule 289 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	 tx155 = tx4*tzz490; 		/* rule 291 */
	 tzz531 = 2.*tzz508; 		/* rule 292 */
	 tx156 = tx5*tzz531; 		/* rule 293 */
	 tx157 = -2.*tx4*tzz494; 		/* rule 294 */
	 tx158 = tx155 + tx156; 		/* rule 295 */
	 tx159 = tx115 + tx116 + tx157; 		/* rule 296 */
	 tx160 = power2(tx70); 		/* rule 297 */
	 tx161 = tx159*tzz471; 		/* rule 298 */
	 tzz523 = tx70*tzz470; 		/* rule 299 */
	 tx162 = tx158*tzz523; 		/* rule 300 */
	 tx163 = tx160*tzz473; 		/* rule 301 */
	 dhx3x3 = tx161 + tx162 + tx163; 		/* rule 302 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 0, dhx3x3);
	 tx164 = tx6*tzz490; 		/* rule 304 */
	 tx165 = tx7*tzz531; 		/* rule 305 */
	 tx166 = -2.*tx6*tzz485; 		/* rule 306 */
	 tx167 = tx164 + tx165; 		/* rule 307 */
	 tx168 = tx116 + tx129 + tx166; 		/* rule 308 */
	 tx169 = power2(tx73); 		/* rule 309 */
	 tx170 = tx168*tzz471; 		/* rule 310 */
	 tzz538 = tx73*tzz470; 		/* rule 311 */
	 tx171 = tx167*tzz538; 		/* rule 312 */
	 tx172 = tx169*tzz473; 		/* rule 313 */
	 dhy3y3 = tx170 + tx171 + tx172; 		/* rule 314 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 1, dhy3y3);
	 tx173 = tx8*tzz490; 		/* rule 316 */
	 tx174 = tx9*tzz531; 		/* rule 317 */
	 tx175 = -2.*tx8*tzz483; 		/* rule 318 */
	 tx176 = tx173 + tx174; 		/* rule 319 */
	 tx177 = tx116 + tx142 + tx175; 		/* rule 320 */
	 tx178 = power2(tx76); 		/* rule 321 */
	 tx179 = tx177*tzz471; 		/* rule 322 */
	 tx180 = tx176*tx76*tzz470; 		/* rule 323 */
	 tx181 = tx178*tzz473; 		/* rule 324 */
	 dhz3z3 = tx179 + tx180 + tx181; 		/* rule 325 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I3, 2, dhz3z3);
	#ifdef ANGLE_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tzz504 = tx6*tzz476; 		/* rule 329 */
	 tx182 = tx5*tzz504; 		/* rule 330 */
	 tzz515 = tx4*tzz476; 		/* rule 331 */
	 tx183 = tx7*tzz515; 		/* rule 332 */
	 tzz541 = tx4*tx6; 		/* rule 333 */
	 tx184 = tzz510*tzz541; 		/* rule 334 */
	 tx185 = tx182 + tx183 + tx184; 		/* rule 335 */
	 tx186 = tx185*tzz471; 		/* rule 336 */
	 tx187 = tx36*tx39*tzz473; 		/* rule 337 */
	 tx188 = tx98*tzz487; 		/* rule 338 */
	 ohx1y1 = tx186 + tx187 + tx188; 		/* rule 339 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tzz498 = tx8*tzz476; 		/* rule 341 */
	 tx189 = tx5*tzz498; 		/* rule 342 */
	 tzz537 = tx8*tzz510; 		/* rule 343 */
	 tx190 = tx4*tzz537; 		/* rule 344 */
	 tx191 = tx9*tzz515; 		/* rule 345 */
	 tx192 = tx189 + tx190 + tx191; 		/* rule 346 */
	 tx193 = tx192*tzz471; 		/* rule 347 */
	 tx194 = tx106*tzz487; 		/* rule 348 */
	 tzz540 = tx42*tzz473; 		/* rule 349 */
	 tx195 = tx36*tzz540; 		/* rule 350 */
	 ohx1z1 = tx193 + tx194 + tx195; 		/* rule 351 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx196 = -tzz495; 		/* rule 353 */
	 tx197 = tx23*tzz469; 		/* rule 354 */
	 tx198 = tx11*tzz467; 		/* rule 355 */
	 tx199 = tx47*tzz515; 		/* rule 356 */
	 tx200 = tx5*tzz496; 		/* rule 357 */
	 tzz491 = -tzz477; 		/* rule 358 */
	 tzz526 = tx5*tzz491; 		/* rule 359 */
	 tx201 = tx4*tzz526; 		/* rule 360 */
	 tzz512 = -3*tzz479; 		/* rule 361 */
	 tx202 = tx10*tzz512; 		/* rule 362 */
	 tzz546 = tx196 + tx197; 		/* rule 363 */
	 tx203 = tx198 + tx199 + tx200 + tx201 + tx202 + tzz546; 		/* rule 364 */
	 tx204 = tx203*tzz471; 		/* rule 365 */
	 tx205 = tx123*tzz487; 		/* rule 366 */
	 tzz525 = tx51*tzz473; 		/* rule 367 */
	 tx206 = tx36*tzz525; 		/* rule 368 */
	 ohx1x2 = tx204 + tx205 + tx206; 		/* rule 369 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx207 = tx55*tzz515; 		/* rule 371 */
	 tx208 = tx5*tzz486; 		/* rule 372 */
	 tzz524 = tx7*tzz491; 		/* rule 373 */
	 tx209 = tx4*tzz524; 		/* rule 374 */
	 tx210 = tx5*tzz485; 		/* rule 375 */
	 tx211 = tzz512*tzz541; 		/* rule 376 */
	 tzz557 = tx210 + tx211; 		/* rule 377 */
	 tx212 = tx207 + tx208 + tx209 + tzz557; 		/* rule 378 */
	 tx213 = tx212*tzz471; 		/* rule 379 */
	 tx214 = tx136*tzz487; 		/* rule 380 */
	 tzz519 = tx59*tzz473; 		/* rule 381 */
	 tx215 = tx36*tzz519; 		/* rule 382 */
	 ohx1y2 = tx213 + tx214 + tx215; 		/* rule 383 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx216 = tx63*tzz515; 		/* rule 385 */
	 tx217 = tx5*tzz484; 		/* rule 386 */
	 tzz536 = tx8*tzz512; 		/* rule 387 */
	 tx218 = tx4*tzz536; 		/* rule 388 */
	 tzz522 = tx9*tzz491; 		/* rule 389 */
	 tx219 = tx4*tzz522; 		/* rule 390 */
	 tx220 = tx5*tzz483; 		/* rule 391 */
	 tzz556 = tx218 + tx220; 		/* rule 392 */
	 tx221 = tx216 + tx217 + tx219 + tzz556; 		/* rule 393 */
	 tx222 = tx221*tzz471; 		/* rule 394 */
	 tx223 = tx149*tzz487; 		/* rule 395 */
	 tzz513 = tx67*tzz473; 		/* rule 396 */
	 tx224 = tx36*tzz513; 		/* rule 397 */
	 ohx1z2 = tx222 + tx223 + tx224; 		/* rule 398 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tx225 = tzz495; 		/* rule 400 */
	 tx226 = tx10*tzz476; 		/* rule 401 */
	 tx227 = tx11*tzz474; 		/* rule 402 */
	 tx228 = tx4*tzz520; 		/* rule 403 */
	 tx229 = tx225 + tx226 + tx227 + tx228; 		/* rule 404 */
	 tx230 = tx229*tzz471; 		/* rule 405 */
	 tx231 = tx158*tzz487; 		/* rule 406 */
	 tzz503 = tx70*tzz473; 		/* rule 407 */
	 tx232 = tx36*tzz503; 		/* rule 408 */
	 ohx1x3 = tx230 + tx231 + tx232; 		/* rule 409 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 0, ohx1x3);
	 tx233 = tx4*tzz504; 		/* rule 411 */
	 tx234 = tx4*tzz517; 		/* rule 412 */
	 tx235 = tx5*tzz500; 		/* rule 413 */
	 tzz555 = tx233 + tx235; 		/* rule 414 */
	 tx236 = tx234 + tzz555; 		/* rule 415 */
	 tx237 = tx236*tzz471; 		/* rule 416 */
	 tx238 = tx167*tzz487; 		/* rule 417 */
	 tzz499 = tx73*tzz473; 		/* rule 418 */
	 tx239 = tx36*tzz499; 		/* rule 419 */
	 ohx1y3 = tx237 + tx238 + tx239; 		/* rule 420 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 1, ohx1y3);
	 tx240 = tx4*tzz498; 		/* rule 422 */
	 tx241 = tx4*tzz516; 		/* rule 423 */
	 tx242 = tx5*tzz492; 		/* rule 424 */
	 tzz554 = tx240 + tx242; 		/* rule 425 */
	 tx243 = tx241 + tzz554; 		/* rule 426 */
	 tx244 = tx243*tzz471; 		/* rule 427 */
	 tx245 = tx176*tzz487; 		/* rule 428 */
	 tzz493 = tx76*tzz473; 		/* rule 429 */
	 tx246 = tx36*tzz493; 		/* rule 430 */
	 ohx1z3 = tx244 + tx245 + tx246; 		/* rule 431 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 2, ohx1z3);
	 tx247 = tx7*tzz498; 		/* rule 433 */
	 tx248 = tx6*tzz537; 		/* rule 434 */
	 tx249 = tx9*tzz504; 		/* rule 435 */
	 tx250 = tx247 + tx248 + tx249; 		/* rule 436 */
	 tx251 = tx250*tzz471; 		/* rule 437 */
	 tx252 = tx106*tzz497; 		/* rule 438 */
	 tx253 = tx39*tzz540; 		/* rule 439 */
	 ohy1z1 = tx251 + tx252 + tx253; 		/* rule 440 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 tx254 = tx47*tzz504; 		/* rule 442 */
	 tx255 = tx6*tzz526; 		/* rule 443 */
	 tx256 = tx7*tzz496; 		/* rule 444 */
	 tx257 = tx254 + tx255 + tx256 + tzz557; 		/* rule 445 */
	 tx258 = tx257*tzz471; 		/* rule 446 */
	 tx259 = tx123*tzz497; 		/* rule 447 */
	 tx260 = tx39*tzz525; 		/* rule 448 */
	 ohy1x2 = tx258 + tx259 + tx260; 		/* rule 449 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx261 = tx13*tzz467; 		/* rule 451 */
	 tx262 = tx55*tzz504; 		/* rule 452 */
	 tx263 = tx7*tzz486; 		/* rule 453 */
	 tx264 = tx6*tzz524; 		/* rule 454 */
	 tx265 = tx12*tzz512; 		/* rule 455 */
	 tx266 = tx261 + tx262 + tx263 + tx264 + tx265 + tzz546; 		/* rule 456 */
	 tx267 = tx266*tzz471; 		/* rule 457 */
	 tx268 = tx136*tzz497; 		/* rule 458 */
	 tx269 = tx39*tzz519; 		/* rule 459 */
	 ohy1y2 = tx267 + tx268 + tx269; 		/* rule 460 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx270 = tx63*tzz504; 		/* rule 462 */
	 tx271 = tx7*tzz484; 		/* rule 463 */
	 tx272 = tx6*tzz536; 		/* rule 464 */
	 tx273 = tx6*tzz522; 		/* rule 465 */
	 tx274 = tx7*tzz483; 		/* rule 466 */
	 tzz553 = tx272 + tx274; 		/* rule 467 */
	 tx275 = tx270 + tx271 + tx273 + tzz553; 		/* rule 468 */
	 tx276 = tx275*tzz471; 		/* rule 469 */
	 tx277 = tx149*tzz497; 		/* rule 470 */
	 tx278 = tx39*tzz513; 		/* rule 471 */
	 ohy1z2 = tx276 + tx277 + tx278; 		/* rule 472 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 tx279 = tx6*tzz520; 		/* rule 474 */
	 tx280 = tx279 + tzz555; 		/* rule 475 */
	 tx281 = tx280*tzz471; 		/* rule 476 */
	 tx282 = tx158*tzz497; 		/* rule 477 */
	 tx283 = tx39*tzz503; 		/* rule 478 */
	 ohy1x3 = tx281 + tx282 + tx283; 		/* rule 479 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 0, ohy1x3);
	 tx284 = tx12*tzz476; 		/* rule 481 */
	 tx285 = tx13*tzz474; 		/* rule 482 */
	 tx286 = tx6*tzz517; 		/* rule 483 */
	 tx287 = tx225 + tx284 + tx285 + tx286; 		/* rule 484 */
	 tx288 = tx287*tzz471; 		/* rule 485 */
	 tx289 = tx167*tzz497; 		/* rule 486 */
	 tx290 = tx39*tzz499; 		/* rule 487 */
	 ohy1y3 = tx288 + tx289 + tx290; 		/* rule 488 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 1, ohy1y3);
	 tx291 = tx6*tzz498; 		/* rule 490 */
	 tx292 = tx6*tzz516; 		/* rule 491 */
	 tx293 = tx7*tzz492; 		/* rule 492 */
	 tzz552 = tx291 + tx293; 		/* rule 493 */
	 tx294 = tx292 + tzz552; 		/* rule 494 */
	 tx295 = tx294*tzz471; 		/* rule 495 */
	 tx296 = tx176*tzz497; 		/* rule 496 */
	 tx297 = tx39*tzz493; 		/* rule 497 */
	 ohy1z3 = tx295 + tx296 + tx297; 		/* rule 498 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 2, ohy1z3);
	 tx298 = tx47*tzz498; 		/* rule 500 */
	 tx299 = tx8*tzz526; 		/* rule 501 */
	 tx300 = tx9*tzz496; 		/* rule 502 */
	 tx301 = tx298 + tx299 + tx300 + tzz556; 		/* rule 503 */
	 tx302 = tx301*tzz471; 		/* rule 504 */
	 tx303 = tx123*tzz501; 		/* rule 505 */
	 tx304 = tx42*tzz525; 		/* rule 506 */
	 ohz1x2 = tx302 + tx303 + tx304; 		/* rule 507 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 tx305 = tx55*tzz498; 		/* rule 509 */
	 tx306 = tx8*tzz524; 		/* rule 510 */
	 tx307 = tx9*tzz486; 		/* rule 511 */
	 tx308 = tx305 + tx306 + tx307 + tzz553; 		/* rule 512 */
	 tx309 = tx308*tzz471; 		/* rule 513 */
	 tx310 = tx136*tzz501; 		/* rule 514 */
	 tx311 = tx42*tzz519; 		/* rule 515 */
	 ohz1y2 = tx309 + tx310 + tx311; 		/* rule 516 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx312 = tx15*tzz467; 		/* rule 518 */
	 tx313 = tx14*tzz512; 		/* rule 519 */
	 tx314 = tx63*tzz498; 		/* rule 520 */
	 tx315 = tx9*tzz484; 		/* rule 521 */
	 tx316 = tx8*tzz522; 		/* rule 522 */
	 tx317 = tx312 + tx313 + tx314 + tx315 + tx316 + tzz546; 		/* rule 523 */
	 tx318 = tx317*tzz471; 		/* rule 524 */
	 tx319 = tx149*tzz501; 		/* rule 525 */
	 tx320 = tx42*tzz513; 		/* rule 526 */
	 ohz1z2 = tx318 + tx319 + tx320; 		/* rule 527 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 tx321 = tx8*tzz520; 		/* rule 529 */
	 tx322 = tx321 + tzz554; 		/* rule 530 */
	 tx323 = tx322*tzz471; 		/* rule 531 */
	 tx324 = tx158*tzz501; 		/* rule 532 */
	 tx325 = tx42*tzz503; 		/* rule 533 */
	 ohz1x3 = tx323 + tx324 + tx325; 		/* rule 534 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 0, ohz1x3);
	 tx326 = tx8*tzz517; 		/* rule 536 */
	 tx327 = tx326 + tzz552; 		/* rule 537 */
	 tx328 = tx327*tzz471; 		/* rule 538 */
	 tx329 = tx167*tzz501; 		/* rule 539 */
	 tx330 = tx42*tzz499; 		/* rule 540 */
	 ohz1y3 = tx328 + tx329 + tx330; 		/* rule 541 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 1, ohz1y3);
	 tx331 = tx14*tzz476; 		/* rule 543 */
	 tx332 = tx15*tzz474; 		/* rule 544 */
	 tx333 = tx8*tzz516; 		/* rule 545 */
	 tx334 = tx225 + tx331 + tx332 + tx333; 		/* rule 546 */
	 tx335 = tx334*tzz471; 		/* rule 547 */
	 tx336 = tx176*tzz501; 		/* rule 548 */
	 tx337 = tx42*tzz493; 		/* rule 549 */
	 ohz1z3 = tx335 + tx336 + tx337; 		/* rule 550 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 2, ohz1z3);
	 tx338 = tx55*tzz496; 		/* rule 552 */
	 tx339 = tx55*tzz494; 		/* rule 553 */
	 tx340 = tx47*tzz486; 		/* rule 554 */
	 tx341 = tx47*tzz485; 		/* rule 555 */
	 tzz539 = tx5*tx7; 		/* rule 556 */
	 tx342 = tzz509*tzz539; 		/* rule 557 */
	 tx343 = tx184 + tx234 + tx279 + tx338 + tx339 + tx340 + tx341 + tx342; 		/* rule 558 */
	 tx344 = tx343*tzz471; 		/* rule 559 */
	 tx345 = tx136*tzz505; 		/* rule 560 */
	 tx346 = tx51*tzz519; 		/* rule 561 */
	 ohx2y2 = tx344 + tx345 + tx346; 		/* rule 562 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 tx347 = tx63*tzz496; 		/* rule 564 */
	 tx348 = tx63*tzz494; 		/* rule 565 */
	 tx349 = tx47*tzz484; 		/* rule 566 */
	 tx350 = tx47*tzz483; 		/* rule 567 */
	 tzz535 = tx9*tzz509; 		/* rule 568 */
	 tx351 = tx5*tzz535; 		/* rule 569 */
	 tx352 = tx190 + tx241 + tx321 + tx347 + tx348 + tx349 + tx350 + tx351; 		/* rule 570 */
	 tx353 = tx352*tzz471; 		/* rule 571 */
	 tx354 = tx149*tzz505; 		/* rule 572 */
	 tx355 = tx51*tzz513; 		/* rule 573 */
	 ohx2z2 = tx353 + tx354 + tx355; 		/* rule 574 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 tzz511 = -3*tzz481; 		/* rule 576 */
	 tx356 = tx11*tzz511; 		/* rule 577 */
	 tx357 = tx10*tzz469; 		/* rule 578 */
	 tx358 = tx23*tzz467; 		/* rule 579 */
	 tx359 = tx4*tzz494; 		/* rule 580 */
	 tx360 = tx47*tzz506; 		/* rule 581 */
	 tzz545 = tx196 + tx358; 		/* rule 582 */
	 tx361 = tx201 + tx356 + tx357 + tx359 + tx360 + tzz545; 		/* rule 583 */
	 tx362 = tx361*tzz471; 		/* rule 584 */
	 tx363 = tx158*tzz505; 		/* rule 585 */
	 tx364 = tx51*tzz503; 		/* rule 586 */
	 ohx2x3 = tx362 + tx363 + tx364; 		/* rule 587 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 0, ohx2x3);
	 tx365 = tx4*tzz486; 		/* rule 589 */
	 tx366 = tx6*tzz494; 		/* rule 590 */
	 tx367 = tx47*tzz500; 		/* rule 591 */
	 tx368 = tzz511*tzz539; 		/* rule 592 */
	 tzz551 = tx365 + tx368; 		/* rule 593 */
	 tx369 = tx209 + tx366 + tx367 + tzz551; 		/* rule 594 */
	 tx370 = tx369*tzz471; 		/* rule 595 */
	 tx371 = tx167*tzz505; 		/* rule 596 */
	 tx372 = tx51*tzz499; 		/* rule 597 */
	 ohx2y3 = tx370 + tx371 + tx372; 		/* rule 598 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 1, ohx2y3);
	 tx373 = tx4*tzz484; 		/* rule 600 */
	 tx374 = tx8*tzz494; 		/* rule 601 */
	 tx375 = tx47*tzz492; 		/* rule 602 */
	 tzz534 = tx9*tzz511; 		/* rule 603 */
	 tx376 = tx5*tzz534; 		/* rule 604 */
	 tzz550 = tx373 + tx376; 		/* rule 605 */
	 tx377 = tx219 + tx374 + tx375 + tzz550; 		/* rule 606 */
	 tx378 = tx377*tzz471; 		/* rule 607 */
	 tx379 = tx176*tzz505; 		/* rule 608 */
	 tx380 = tx51*tzz493; 		/* rule 609 */
	 ohx2z3 = tx378 + tx379 + tx380; 		/* rule 610 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 2, ohx2z3);
	 tx381 = tx63*tzz486; 		/* rule 612 */
	 tx382 = tx63*tzz485; 		/* rule 613 */
	 tx383 = tx55*tzz484; 		/* rule 614 */
	 tx384 = tx55*tzz483; 		/* rule 615 */
	 tx385 = tx7*tzz535; 		/* rule 616 */
	 tx386 = tx248 + tx292 + tx326 + tx381 + tx382 + tx383 + tx384 + tx385; 		/* rule 617 */
	 tx387 = tx386*tzz471; 		/* rule 618 */
	 tx388 = tx149*tzz514; 		/* rule 619 */
	 tx389 = tx59*tzz513; 		/* rule 620 */
	 ohy2z2 = tx387 + tx388 + tx389; 		/* rule 621 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	 tx390 = tx55*tzz506; 		/* rule 623 */
	 tx391 = tx4*tzz485; 		/* rule 624 */
	 tx392 = tx255 + tx390 + tx391 + tzz551; 		/* rule 625 */
	 tx393 = tx392*tzz471; 		/* rule 626 */
	 tx394 = tx158*tzz514; 		/* rule 627 */
	 tx395 = tx59*tzz503; 		/* rule 628 */
	 ohy2x3 = tx393 + tx394 + tx395; 		/* rule 629 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 0, ohy2x3);
	 tx396 = tx13*tzz511; 		/* rule 631 */
	 tx397 = tx12*tzz469; 		/* rule 632 */
	 tx398 = tx55*tzz500; 		/* rule 633 */
	 tx399 = tx6*tzz485; 		/* rule 634 */
	 tx400 = tx264 + tx396 + tx397 + tx398 + tx399 + tzz545; 		/* rule 635 */
	 tx401 = tx400*tzz471; 		/* rule 636 */
	 tx402 = tx167*tzz514; 		/* rule 637 */
	 tx403 = tx59*tzz499; 		/* rule 638 */
	 ohy2y3 = tx401 + tx402 + tx403; 		/* rule 639 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 1, ohy2y3);
	 tx404 = tx6*tzz484; 		/* rule 641 */
	 tx405 = tx8*tzz485; 		/* rule 642 */
	 tx406 = tx55*tzz492; 		/* rule 643 */
	 tx407 = tx7*tzz534; 		/* rule 644 */
	 tzz549 = tx404 + tx407; 		/* rule 645 */
	 tx408 = tx273 + tx405 + tx406 + tzz549; 		/* rule 646 */
	 tx409 = tx408*tzz471; 		/* rule 647 */
	 tx410 = tx176*tzz514; 		/* rule 648 */
	 tx411 = tx59*tzz493; 		/* rule 649 */
	 ohy2z3 = tx409 + tx410 + tx411; 		/* rule 650 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 2, ohy2z3);
	 tx412 = tx63*tzz506; 		/* rule 652 */
	 tx413 = tx4*tzz483; 		/* rule 653 */
	 tx414 = tx299 + tx412 + tx413 + tzz550; 		/* rule 654 */
	 tx415 = tx414*tzz471; 		/* rule 655 */
	 tx416 = tx158*tzz518; 		/* rule 656 */
	 tx417 = tx67*tzz503; 		/* rule 657 */
	 ohz2x3 = tx415 + tx416 + tx417; 		/* rule 658 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 0, ohz2x3);
	 tx418 = tx63*tzz500; 		/* rule 660 */
	 tx419 = tx6*tzz483; 		/* rule 661 */
	 tx420 = tx306 + tx418 + tx419 + tzz549; 		/* rule 662 */
	 tx421 = tx420*tzz471; 		/* rule 663 */
	 tx422 = tx167*tzz518; 		/* rule 664 */
	 tx423 = tx67*tzz499; 		/* rule 665 */
	 ohz2y3 = tx421 + tx422 + tx423; 		/* rule 666 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 1, ohz2y3);
	 tx424 = tx15*tzz511; 		/* rule 668 */
	 tx425 = tx14*tzz469; 		/* rule 669 */
	 tx426 = tx63*tzz492; 		/* rule 670 */
	 tx427 = tx8*tzz483; 		/* rule 671 */
	 tx428 = tx316 + tx424 + tx425 + tx426 + tx427 + tzz545; 		/* rule 672 */
	 tx429 = tx428*tzz471; 		/* rule 673 */
	 tx430 = tx176*tzz518; 		/* rule 674 */
	 tx431 = tx67*tzz493; 		/* rule 675 */
	 ohz2z3 = tx429 + tx430 + tx431; 		/* rule 676 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 2, ohz2z3);
	 tx432 = tx6*tzz506; 		/* rule 678 */
	 tx433 = tx4*tzz500; 		/* rule 679 */
	 tx434 = tx342 + tx432 + tx433; 		/* rule 680 */
	 tx435 = tx434*tzz471; 		/* rule 681 */
	 tx436 = tx167*tzz523; 		/* rule 682 */
	 tx437 = tx70*tzz499; 		/* rule 683 */
	 ohx3y3 = tx435 + tx436 + tx437; 		/* rule 684 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 1, ohx3y3);
	 tx438 = tx8*tzz506; 		/* rule 686 */
	 tx439 = tx4*tzz492; 		/* rule 687 */
	 tx440 = tx351 + tx438 + tx439; 		/* rule 688 */
	 tx441 = tx440*tzz471; 		/* rule 689 */
	 tx442 = tx176*tzz523; 		/* rule 690 */
	 tx443 = tx70*tzz493; 		/* rule 691 */
	 ohx3z3 = tx441 + tx442 + tx443; 		/* rule 692 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 2, ohx3z3);
	 tx444 = tx8*tzz500; 		/* rule 694 */
	 tx445 = tx6*tzz492; 		/* rule 695 */
	 tx446 = tx385 + tx444 + tx445; 		/* rule 696 */
	 tx447 = tx446*tzz471; 		/* rule 697 */
	 tx448 = tx176*tzz538; 		/* rule 698 */
	 tx449 = tx73*tzz493; 		/* rule 699 */
	 ohy3z3 = tx447 + tx448 + tx449; 		/* rule 700 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 2, ohy3z3);
	} /*if calcOffDiagonalHessian */ 
	#endif /* ANGLE_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* ANGLE_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* ANGLE_CALC_FORCE ]*/
	ANGLE_DEBUG_INTERACTIONS(I1,I2,I3);
	SKIP_term: (void)0;
