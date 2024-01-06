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
	 tzz498 = tx27*tx29; 		/* rule 59 */
	 tx31 = -(tx28*tzz498); 		/* rule 60 */
	 tx32 = 1 + tx31; 		/* rule 61 */
	 tzz474 = tx23*tx30; 		/* rule 62 */
	 tzz517 = tx22*tzz474; 		/* rule 63 */
	 tzz522 = -tzz517; 		/* rule 64 */
	 tx33 = tx4*tzz522; 		/* rule 65 */
	 tzz518 = tx22*tx453; 		/* rule 66 */
	 tx34 = tx5*tzz518; 		/* rule 67 */
	 tx455 = mysqrt(tx32); 		/* rule 68 */
	 tx35 = reciprocal(tx455); 		/* rule 69 */
	 tx36 = tx33 + tx34; 		/* rule 70 */
	 tzz466 = AngleDeviation*kt; 		/* rule 71 */
	 tzz467 = -2.*tzz466; 		/* rule 72 */
	 tzz484 = tx35*tzz467; 		/* rule 73 */
	 gx1 = tx36*tzz484; 		/* rule 74 */
	 fx1 = -gx1; 		/* rule 75 */
	ANGLE_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx37 = tx6*tzz522; 		/* rule 77 */
	 tx38 = tx7*tzz518; 		/* rule 78 */
	 tx39 = tx37 + tx38; 		/* rule 79 */
	 gy1 = tx39*tzz484; 		/* rule 80 */
	 fy1 = -gy1; 		/* rule 81 */
	ANGLE_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx40 = tx8*tzz522; 		/* rule 83 */
	 tx41 = tx9*tzz518; 		/* rule 84 */
	 tx42 = tx40 + tx41; 		/* rule 85 */
	 gz1 = tx42*tzz484; 		/* rule 86 */
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
	 tzz468 = tx30*tx457; 		/* rule 97 */
	 tzz502 = tx4*tzz468; 		/* rule 98 */
	 tx48 = tx23*tzz502; 		/* rule 99 */
	 tzz491 = tx453*tx457; 		/* rule 100 */
	 tx49 = tx47*tzz491; 		/* rule 101 */
	 tzz473 = tx453*tx46; 		/* rule 102 */
	 tzz496 = tx5*tzz473; 		/* rule 103 */
	 tx50 = tx23*tzz496; 		/* rule 104 */
	 tx51 = tx48 + tx49 + tx50; 		/* rule 105 */
	 gx2 = tx51*tzz484; 		/* rule 106 */
	 fx2 = -gx2; 		/* rule 107 */
	ANGLE_FORCE_ACCUMULATE(I2, 0, fx2 );
	 tx52 = -y1; 		/* rule 109 */
	 tx53 = 2.*y2; 		/* rule 110 */
	 tx54 = -y3; 		/* rule 111 */
	 tx55 = tx52 + tx53 + tx54; 		/* rule 112 */
	 tx56 = tx55*tzz491; 		/* rule 113 */
	 tzz495 = tx6*tzz468; 		/* rule 114 */
	 tx57 = tx23*tzz495; 		/* rule 115 */
	 tzz490 = tx7*tzz473; 		/* rule 116 */
	 tx58 = tx23*tzz490; 		/* rule 117 */
	 tx59 = tx56 + tx57 + tx58; 		/* rule 118 */
	 gy2 = tx59*tzz484; 		/* rule 119 */
	 fy2 = -gy2; 		/* rule 120 */
	ANGLE_FORCE_ACCUMULATE(I2, 1, fy2 );
	 tx60 = -z1; 		/* rule 122 */
	 tx61 = 2.*z2; 		/* rule 123 */
	 tx62 = -z3; 		/* rule 124 */
	 tx63 = tx60 + tx61 + tx62; 		/* rule 125 */
	 tx64 = tx63*tzz491; 		/* rule 126 */
	 tzz488 = tx8*tzz468; 		/* rule 127 */
	 tx65 = tx23*tzz488; 		/* rule 128 */
	 tzz482 = tx9*tzz473; 		/* rule 129 */
	 tx66 = tx23*tzz482; 		/* rule 130 */
	 tx67 = tx64 + tx65 + tx66; 		/* rule 131 */
	 gz2 = tx67*tzz484; 		/* rule 132 */
	 fz2 = -gz2; 		/* rule 133 */
	ANGLE_FORCE_ACCUMULATE(I2, 2, fz2 );
	 tx68 = tx4*tzz491; 		/* rule 135 */
	 tx69 = -tx50; 		/* rule 136 */
	 tx70 = tx68 + tx69; 		/* rule 137 */
	 gx3 = tx70*tzz484; 		/* rule 138 */
	 fx3 = -gx3; 		/* rule 139 */
	ANGLE_FORCE_ACCUMULATE(I3, 0, fx3 );
	 tx71 = -tx58; 		/* rule 141 */
	 tx72 = tx6*tzz491; 		/* rule 142 */
	 tx73 = tx71 + tx72; 		/* rule 143 */
	 gy3 = tx73*tzz484; 		/* rule 144 */
	 fy3 = -gy3; 		/* rule 145 */
	ANGLE_FORCE_ACCUMULATE(I3, 1, fy3 );
	 tx74 = -tx66; 		/* rule 147 */
	 tx75 = tx8*tzz491; 		/* rule 148 */
	 tx76 = tx74 + tx75; 		/* rule 149 */
	 gz3 = tx76*tzz484; 		/* rule 150 */
	 fz3 = -gz3; 		/* rule 151 */
	ANGLE_FORCE_ACCUMULATE(I3, 2, fz3 );
	#ifdef ANGLE_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx459 = power2(tx452); 		/* rule 155 */
	 tx460 = reciprocal(tx459); 		/* rule 156 */
	 tx77 = tx453*tx460; 		/* rule 157 */
	 tx78 = power2(tx27); 		/* rule 158 */
	 tzz525 = -2.*tzz468; 		/* rule 159 */
	 tx79 = tx18*tzz525; 		/* rule 160 */
	 tzz475 = -tzz468; 		/* rule 161 */
	 tx80 = tx23*tzz475; 		/* rule 162 */
	 tzz481 = -2.*tx28; 		/* rule 163 */
	 tzz485 = tx27*tzz481; 		/* rule 164 */
	 tzz486 = tx23*tzz485; 		/* rule 165 */
	 tx81 = tx5*tzz486; 		/* rule 166 */
	 tzz477 = tx457*tx77; 		/* rule 167 */
	 tzz478 = tx23*tzz477; 		/* rule 168 */
	 tzz506 = 3*tzz478; 		/* rule 169 */
	 tx82 = tx10*tzz506; 		/* rule 170 */
	 tzz503 = tx29*tx78; 		/* rule 171 */
	 tzz516 = tx28*tzz503; 		/* rule 172 */
	 tzz519 = 2.*tzz516; 		/* rule 173 */
	 tx83 = tx4*tzz519; 		/* rule 174 */
	 tx461 = tx32; 		/* rule 175 */
	 tx462 = tx35; 		/* rule 176 */
	 tx463 = reciprocal(tx461); 		/* rule 177 */
	 tx84 = tx462*tx463; 		/* rule 178 */
	 tx85 = tx463; 		/* rule 179 */
	 tx86 = power2(tx36); 		/* rule 180 */
	 tx87 = tx79 + tx80 + tx82; 		/* rule 181 */
	 tx88 = tx81 + tx83; 		/* rule 182 */
	 tzz471 = kt*tx85; 		/* rule 183 */
	 tzz472 = 2.*tzz471; 		/* rule 184 */
	 tx89 = tx86*tzz472; 		/* rule 185 */
	 tzz470 = tx462*tzz467; 		/* rule 186 */
	 tx90 = tx87*tzz470; 		/* rule 187 */
	 tzz469 = tx84*tzz466; 		/* rule 188 */
	 tzz483 = tx36*tzz469; 		/* rule 189 */
	 tx91 = tx88*tzz483; 		/* rule 190 */
	 dhx1x1 = tx89 + tx90 + tx91; 		/* rule 191 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx92 = tx19*tzz525; 		/* rule 193 */
	 tx93 = tx7*tzz486; 		/* rule 194 */
	 tx94 = tx12*tzz506; 		/* rule 195 */
	 tx95 = tx6*tzz519; 		/* rule 196 */
	 tx96 = power2(tx39); 		/* rule 197 */
	 tx97 = tx80 + tx92 + tx94; 		/* rule 198 */
	 tx98 = tx93 + tx95; 		/* rule 199 */
	 tx99 = tx96*tzz472; 		/* rule 200 */
	 tx100 = tx97*tzz470; 		/* rule 201 */
	 tzz492 = tx39*tzz469; 		/* rule 202 */
	 tx101 = tx98*tzz492; 		/* rule 203 */
	 dhy1y1 = tx100 + tx101 + tx99; 		/* rule 204 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx102 = tx20*tzz525; 		/* rule 206 */
	 tx103 = tx14*tzz506; 		/* rule 207 */
	 tx104 = tx8*tzz519; 		/* rule 208 */
	 tx105 = tx9*tzz486; 		/* rule 209 */
	 tx106 = tx104 + tx105; 		/* rule 210 */
	 tx107 = power2(tx42); 		/* rule 211 */
	 tx108 = tx102 + tx103 + tx80; 		/* rule 212 */
	 tx109 = tx108*tzz470; 		/* rule 213 */
	 tzz497 = tx42*tzz469; 		/* rule 214 */
	 tx110 = tx106*tzz497; 		/* rule 215 */
	 tx111 = tx107*tzz472; 		/* rule 216 */
	 dhz1z1 = tx109 + tx110 + tx111; 		/* rule 217 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 tx464 = power2(tx456); 		/* rule 219 */
	 tx465 = reciprocal(tx464); 		/* rule 220 */
	 tx112 = tx457*tx465; 		/* rule 221 */
	 tx113 = power2(tx28); 		/* rule 222 */
	 tx114 = 2.*tzz491; 		/* rule 223 */
	 tzz479 = tx23*tx453; 		/* rule 224 */
	 tzz480 = tx112*tzz479; 		/* rule 225 */
	 tzz505 = 3*tzz480; 		/* rule 226 */
	 tx115 = tx11*tzz505; 		/* rule 227 */
	 tx116 = -(tx23*tzz473); 		/* rule 228 */
	 tzz476 = tx46*tzz474; 		/* rule 229 */
	 tzz521 = 2.*tzz476; 		/* rule 230 */
	 tx117 = tx18*tzz521; 		/* rule 231 */
	 tx118 = tx47*tzz486; 		/* rule 232 */
	 tx119 = 2.*tx47*tzz502; 		/* rule 233 */
	 tzz504 = tx113*tzz498; 		/* rule 234 */
	 tzz523 = -2.*tzz504; 		/* rule 235 */
	 tx120 = tx5*tzz523; 		/* rule 236 */
	 tzz539 = tx47*tzz496; 		/* rule 237 */
	 tx121 = 2.*tzz539; 		/* rule 238 */
	 tx122 = tx4*tzz481*tzz503; 		/* rule 239 */
	 tx123 = tx118 + tx120 + tx122; 		/* rule 240 */
	 tx124 = power2(tx51); 		/* rule 241 */
	 tzz545 = tx116 + tx80; 		/* rule 242 */
	 tzz546 = tx114 + tzz545; 		/* rule 243 */
	 tx125 = tx115 + tx117 + tx119 + tx121 + tx82 + tzz546; 		/* rule 244 */
	 tx126 = tx125*tzz470; 		/* rule 245 */
	 tzz501 = tx51*tzz469; 		/* rule 246 */
	 tx127 = tx123*tzz501; 		/* rule 247 */
	 tx128 = tx124*tzz472; 		/* rule 248 */
	 dhx2x2 = tx126 + tx127 + tx128; 		/* rule 249 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 tx129 = tx13*tzz505; 		/* rule 251 */
	 tx130 = tx19*tzz521; 		/* rule 252 */
	 tx131 = tx55*tzz486; 		/* rule 253 */
	 tx132 = 2.*tx55*tzz495; 		/* rule 254 */
	 tx133 = tx7*tzz523; 		/* rule 255 */
	 tzz536 = tx55*tzz490; 		/* rule 256 */
	 tx134 = 2.*tzz536; 		/* rule 257 */
	 tx135 = tx6*tzz481*tzz503; 		/* rule 258 */
	 tx136 = tx131 + tx133 + tx135; 		/* rule 259 */
	 tx137 = power2(tx59); 		/* rule 260 */
	 tx138 = tx129 + tx130 + tx132 + tx134 + tx94 + tzz546; 		/* rule 261 */
	 tx139 = tx138*tzz470; 		/* rule 262 */
	 tzz510 = tx59*tzz469; 		/* rule 263 */
	 tx140 = tx136*tzz510; 		/* rule 264 */
	 tx141 = tx137*tzz472; 		/* rule 265 */
	 dhy2y2 = tx139 + tx140 + tx141; 		/* rule 266 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 tx142 = tx15*tzz505; 		/* rule 268 */
	 tx143 = tx20*tzz521; 		/* rule 269 */
	 tx144 = tx63*tzz486; 		/* rule 270 */
	 tx145 = 2.*tx63*tzz488; 		/* rule 271 */
	 tx146 = tx8*tzz481*tzz503; 		/* rule 272 */
	 tx147 = tx9*tzz523; 		/* rule 273 */
	 tzz535 = tx63*tzz482; 		/* rule 274 */
	 tx148 = 2.*tzz535; 		/* rule 275 */
	 tx149 = tx144 + tx146 + tx147; 		/* rule 276 */
	 tx150 = power2(tx67); 		/* rule 277 */
	 tx151 = tx103 + tx142 + tx143 + tx145 + tx148 + tzz546; 		/* rule 278 */
	 tx152 = tx151*tzz470; 		/* rule 279 */
	 tzz512 = tx67*tzz469; 		/* rule 280 */
	 tx153 = tx149*tzz512; 		/* rule 281 */
	 tx154 = tx150*tzz472; 		/* rule 282 */
	 dhz2z2 = tx152 + tx153 + tx154; 		/* rule 283 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	 tx155 = tx4*tzz486; 		/* rule 285 */
	 tzz524 = -2.*tzz473; 		/* rule 286 */
	 tx156 = tx18*tzz524; 		/* rule 287 */
	 tzz520 = 2.*tzz504; 		/* rule 288 */
	 tx157 = tx5*tzz520; 		/* rule 289 */
	 tx158 = tx115 + tx116 + tx156; 		/* rule 290 */
	 tx159 = tx155 + tx157; 		/* rule 291 */
	 tx160 = power2(tx70); 		/* rule 292 */
	 tx161 = tx158*tzz470; 		/* rule 293 */
	 tzz514 = tx70*tzz469; 		/* rule 294 */
	 tx162 = tx159*tzz514; 		/* rule 295 */
	 tx163 = tx160*tzz472; 		/* rule 296 */
	 dhx3x3 = tx161 + tx162 + tx163; 		/* rule 297 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 0, dhx3x3);
	 tx164 = tx19*tzz524; 		/* rule 299 */
	 tx165 = tx6*tzz486; 		/* rule 300 */
	 tx166 = tx7*tzz520; 		/* rule 301 */
	 tx167 = tx116 + tx129 + tx164; 		/* rule 302 */
	 tx168 = tx165 + tx166; 		/* rule 303 */
	 tx169 = power2(tx73); 		/* rule 304 */
	 tx170 = tx167*tzz470; 		/* rule 305 */
	 tzz534 = tx73*tzz469; 		/* rule 306 */
	 tx171 = tx168*tzz534; 		/* rule 307 */
	 tx172 = tx169*tzz472; 		/* rule 308 */
	 dhy3y3 = tx170 + tx171 + tx172; 		/* rule 309 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 1, dhy3y3);
	 tx173 = tx20*tzz524; 		/* rule 311 */
	 tx174 = tx8*tzz486; 		/* rule 312 */
	 tx175 = tx9*tzz520; 		/* rule 313 */
	 tx176 = tx116 + tx142 + tx173; 		/* rule 314 */
	 tx177 = tx174 + tx175; 		/* rule 315 */
	 tx178 = power2(tx76); 		/* rule 316 */
	 tx179 = tx176*tzz470; 		/* rule 317 */
	 tx180 = tx177*tx76*tzz469; 		/* rule 318 */
	 tx181 = tx178*tzz472; 		/* rule 319 */
	 dhz3z3 = tx179 + tx180 + tx181; 		/* rule 320 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I3, 2, dhz3z3);
	#ifdef ANGLE_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tzz500 = tx6*tzz475; 		/* rule 324 */
	 tx182 = tx5*tzz500; 		/* rule 325 */
	 tzz511 = tx4*tzz475; 		/* rule 326 */
	 tx183 = tx7*tzz511; 		/* rule 327 */
	 tzz542 = tx4*tx6; 		/* rule 328 */
	 tx184 = tzz506*tzz542; 		/* rule 329 */
	 tx185 = tx182 + tx183 + tx184; 		/* rule 330 */
	 tx186 = tx185*tzz470; 		/* rule 331 */
	 tx187 = tx36*tx39*tzz472; 		/* rule 332 */
	 tx188 = tx98*tzz483; 		/* rule 333 */
	 ohx1y1 = tx186 + tx187 + tx188; 		/* rule 334 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tzz493 = tx8*tzz475; 		/* rule 336 */
	 tx189 = tx5*tzz493; 		/* rule 337 */
	 tzz531 = tx8*tzz506; 		/* rule 338 */
	 tx190 = tx4*tzz531; 		/* rule 339 */
	 tx191 = tx9*tzz511; 		/* rule 340 */
	 tx192 = tx189 + tx190 + tx191; 		/* rule 341 */
	 tx193 = tx192*tzz470; 		/* rule 342 */
	 tx194 = tx106*tzz483; 		/* rule 343 */
	 tzz540 = tx42*tzz472; 		/* rule 344 */
	 tx195 = tx36*tzz540; 		/* rule 345 */
	 ohx1z1 = tx193 + tx194 + tx195; 		/* rule 346 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx196 = -tzz491; 		/* rule 348 */
	 tx197 = tx18*tzz468; 		/* rule 349 */
	 tx198 = tx23*tzz468; 		/* rule 350 */
	 tx199 = tx11*tzz473; 		/* rule 351 */
	 tzz487 = -tzz476; 		/* rule 352 */
	 tx200 = tx18*tzz487; 		/* rule 353 */
	 tx201 = tx47*tzz511; 		/* rule 354 */
	 tzz508 = -3*tzz478; 		/* rule 355 */
	 tx202 = tx10*tzz508; 		/* rule 356 */
	 tzz544 = tx196 + tx198; 		/* rule 357 */
	 tx203 = tx197 + tx199 + tx200 + tx201 + tx202 + tzz544; 		/* rule 358 */
	 tx204 = tx203*tzz470; 		/* rule 359 */
	 tx205 = tx123*tzz483; 		/* rule 360 */
	 tzz515 = tx51*tzz472; 		/* rule 361 */
	 tx206 = tx36*tzz515; 		/* rule 362 */
	 ohx1x2 = tx204 + tx205 + tx206; 		/* rule 363 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx207 = tx55*tzz511; 		/* rule 365 */
	 tx208 = tx5*tzz495; 		/* rule 366 */
	 tzz541 = tx4*tx7; 		/* rule 367 */
	 tx209 = tzz487*tzz541; 		/* rule 368 */
	 tx210 = tx5*tzz490; 		/* rule 369 */
	 tx211 = tzz508*tzz542; 		/* rule 370 */
	 tzz555 = tx210 + tx211; 		/* rule 371 */
	 tx212 = tx207 + tx208 + tx209 + tzz555; 		/* rule 372 */
	 tx213 = tx212*tzz470; 		/* rule 373 */
	 tx214 = tx136*tzz483; 		/* rule 374 */
	 tzz513 = tx59*tzz472; 		/* rule 375 */
	 tx215 = tx36*tzz513; 		/* rule 376 */
	 ohx1y2 = tx213 + tx214 + tx215; 		/* rule 377 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx216 = tx63*tzz511; 		/* rule 379 */
	 tx217 = tx5*tzz488; 		/* rule 380 */
	 tzz530 = tx8*tzz508; 		/* rule 381 */
	 tx218 = tx4*tzz530; 		/* rule 382 */
	 tzz528 = tx9*tzz487; 		/* rule 383 */
	 tx219 = tx4*tzz528; 		/* rule 384 */
	 tx220 = tx5*tzz482; 		/* rule 385 */
	 tzz554 = tx218 + tx220; 		/* rule 386 */
	 tx221 = tx216 + tx217 + tx219 + tzz554; 		/* rule 387 */
	 tx222 = tx221*tzz470; 		/* rule 388 */
	 tx223 = tx149*tzz483; 		/* rule 389 */
	 tzz509 = tx67*tzz472; 		/* rule 390 */
	 tx224 = tx36*tzz509; 		/* rule 391 */
	 ohx1z2 = tx222 + tx223 + tx224; 		/* rule 392 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tx225 = -tx199; 		/* rule 394 */
	 tx226 = tzz491; 		/* rule 395 */
	 tx227 = tx10*tzz475; 		/* rule 396 */
	 tx228 = tx18*tzz476; 		/* rule 397 */
	 tx229 = tx225 + tx226 + tx227 + tx228; 		/* rule 398 */
	 tx230 = tx229*tzz470; 		/* rule 399 */
	 tx231 = tx159*tzz483; 		/* rule 400 */
	 tzz499 = tx70*tzz472; 		/* rule 401 */
	 tx232 = tx36*tzz499; 		/* rule 402 */
	 ohx1x3 = tx230 + tx231 + tx232; 		/* rule 403 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 0, ohx1x3);
	 tx233 = -tx210; 		/* rule 405 */
	 tx234 = tx4*tzz500; 		/* rule 406 */
	 tx235 = tzz476*tzz541; 		/* rule 407 */
	 tzz553 = tx233 + tx234; 		/* rule 408 */
	 tx236 = tx235 + tzz553; 		/* rule 409 */
	 tx237 = tx236*tzz470; 		/* rule 410 */
	 tx238 = tx168*tzz483; 		/* rule 411 */
	 tzz494 = tx73*tzz472; 		/* rule 412 */
	 tx239 = tx36*tzz494; 		/* rule 413 */
	 ohx1y3 = tx237 + tx238 + tx239; 		/* rule 414 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 1, ohx1y3);
	 tx240 = -tx220; 		/* rule 416 */
	 tx241 = tx4*tzz493; 		/* rule 417 */
	 tzz529 = tx9*tzz476; 		/* rule 418 */
	 tx242 = tx4*tzz529; 		/* rule 419 */
	 tzz552 = tx240 + tx241; 		/* rule 420 */
	 tx243 = tx242 + tzz552; 		/* rule 421 */
	 tx244 = tx243*tzz470; 		/* rule 422 */
	 tx245 = tx177*tzz483; 		/* rule 423 */
	 tzz489 = tx76*tzz472; 		/* rule 424 */
	 tx246 = tx36*tzz489; 		/* rule 425 */
	 ohx1z3 = tx244 + tx245 + tx246; 		/* rule 426 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 2, ohx1z3);
	 tx247 = tx7*tzz493; 		/* rule 428 */
	 tx248 = tx6*tzz531; 		/* rule 429 */
	 tx249 = tx9*tzz500; 		/* rule 430 */
	 tx250 = tx247 + tx248 + tx249; 		/* rule 431 */
	 tx251 = tx250*tzz470; 		/* rule 432 */
	 tx252 = tx106*tzz492; 		/* rule 433 */
	 tx253 = tx39*tzz540; 		/* rule 434 */
	 ohy1z1 = tx251 + tx252 + tx253; 		/* rule 435 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 tx254 = tx47*tzz500; 		/* rule 437 */
	 tzz538 = tx5*tx6; 		/* rule 438 */
	 tx255 = tzz487*tzz538; 		/* rule 439 */
	 tx256 = tx7*tzz502; 		/* rule 440 */
	 tx257 = tx254 + tx255 + tx256 + tzz555; 		/* rule 441 */
	 tx258 = tx257*tzz470; 		/* rule 442 */
	 tx259 = tx123*tzz492; 		/* rule 443 */
	 tx260 = tx39*tzz515; 		/* rule 444 */
	 ohy1x2 = tx258 + tx259 + tx260; 		/* rule 445 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx261 = tx19*tzz468; 		/* rule 447 */
	 tx262 = tx13*tzz473; 		/* rule 448 */
	 tx263 = tx19*tzz487; 		/* rule 449 */
	 tx264 = tx55*tzz500; 		/* rule 450 */
	 tx265 = tx12*tzz508; 		/* rule 451 */
	 tx266 = tx261 + tx262 + tx263 + tx264 + tx265 + tzz544; 		/* rule 452 */
	 tx267 = tx266*tzz470; 		/* rule 453 */
	 tx268 = tx136*tzz492; 		/* rule 454 */
	 tx269 = tx39*tzz513; 		/* rule 455 */
	 ohy1y2 = tx267 + tx268 + tx269; 		/* rule 456 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx270 = tx63*tzz500; 		/* rule 458 */
	 tx271 = tx7*tzz488; 		/* rule 459 */
	 tx272 = tx6*tzz530; 		/* rule 460 */
	 tx273 = tx6*tzz528; 		/* rule 461 */
	 tx274 = tx7*tzz482; 		/* rule 462 */
	 tzz551 = tx272 + tx274; 		/* rule 463 */
	 tx275 = tx270 + tx271 + tx273 + tzz551; 		/* rule 464 */
	 tx276 = tx275*tzz470; 		/* rule 465 */
	 tx277 = tx149*tzz492; 		/* rule 466 */
	 tx278 = tx39*tzz509; 		/* rule 467 */
	 ohy1z2 = tx276 + tx277 + tx278; 		/* rule 468 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 tx279 = tzz476*tzz538; 		/* rule 470 */
	 tx280 = tx279 + tzz553; 		/* rule 471 */
	 tx281 = tx280*tzz470; 		/* rule 472 */
	 tx282 = tx159*tzz492; 		/* rule 473 */
	 tx283 = tx39*tzz499; 		/* rule 474 */
	 ohy1x3 = tx281 + tx282 + tx283; 		/* rule 475 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 0, ohy1x3);
	 tx284 = -tx262; 		/* rule 477 */
	 tx285 = tx12*tzz475; 		/* rule 478 */
	 tx286 = tx19*tzz476; 		/* rule 479 */
	 tx287 = tx226 + tx284 + tx285 + tx286; 		/* rule 480 */
	 tx288 = tx287*tzz470; 		/* rule 481 */
	 tx289 = tx168*tzz492; 		/* rule 482 */
	 tx290 = tx39*tzz494; 		/* rule 483 */
	 ohy1y3 = tx288 + tx289 + tx290; 		/* rule 484 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 1, ohy1y3);
	 tx291 = -tx274; 		/* rule 486 */
	 tx292 = tx6*tzz493; 		/* rule 487 */
	 tx293 = tx6*tzz529; 		/* rule 488 */
	 tzz550 = tx291 + tx292; 		/* rule 489 */
	 tx294 = tx293 + tzz550; 		/* rule 490 */
	 tx295 = tx294*tzz470; 		/* rule 491 */
	 tx296 = tx177*tzz492; 		/* rule 492 */
	 tx297 = tx39*tzz489; 		/* rule 493 */
	 ohy1z3 = tx295 + tx296 + tx297; 		/* rule 494 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 2, ohy1z3);
	 tx298 = tx47*tzz493; 		/* rule 496 */
	 tzz532 = tx8*tzz487; 		/* rule 497 */
	 tx299 = tx5*tzz532; 		/* rule 498 */
	 tx300 = tx9*tzz502; 		/* rule 499 */
	 tx301 = tx298 + tx299 + tx300 + tzz554; 		/* rule 500 */
	 tx302 = tx301*tzz470; 		/* rule 501 */
	 tx303 = tx123*tzz497; 		/* rule 502 */
	 tx304 = tx42*tzz515; 		/* rule 503 */
	 ohz1x2 = tx302 + tx303 + tx304; 		/* rule 504 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 tx305 = tx55*tzz493; 		/* rule 506 */
	 tx306 = tx7*tzz532; 		/* rule 507 */
	 tx307 = tx9*tzz495; 		/* rule 508 */
	 tx308 = tx305 + tx306 + tx307 + tzz551; 		/* rule 509 */
	 tx309 = tx308*tzz470; 		/* rule 510 */
	 tx310 = tx136*tzz497; 		/* rule 511 */
	 tx311 = tx42*tzz513; 		/* rule 512 */
	 ohz1y2 = tx309 + tx310 + tx311; 		/* rule 513 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx312 = tx20*tzz468; 		/* rule 515 */
	 tx313 = tx15*tzz473; 		/* rule 516 */
	 tx314 = tx20*tzz487; 		/* rule 517 */
	 tx315 = tx14*tzz508; 		/* rule 518 */
	 tx316 = tx63*tzz493; 		/* rule 519 */
	 tx317 = tx312 + tx313 + tx314 + tx315 + tx316 + tzz544; 		/* rule 520 */
	 tx318 = tx317*tzz470; 		/* rule 521 */
	 tx319 = tx149*tzz497; 		/* rule 522 */
	 tx320 = tx42*tzz509; 		/* rule 523 */
	 ohz1z2 = tx318 + tx319 + tx320; 		/* rule 524 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 tzz533 = tx8*tzz476; 		/* rule 526 */
	 tx321 = tx5*tzz533; 		/* rule 527 */
	 tx322 = tx321 + tzz552; 		/* rule 528 */
	 tx323 = tx322*tzz470; 		/* rule 529 */
	 tx324 = tx159*tzz497; 		/* rule 530 */
	 tx325 = tx42*tzz499; 		/* rule 531 */
	 ohz1x3 = tx323 + tx324 + tx325; 		/* rule 532 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 0, ohz1x3);
	 tx326 = tx7*tzz533; 		/* rule 534 */
	 tx327 = tx326 + tzz550; 		/* rule 535 */
	 tx328 = tx327*tzz470; 		/* rule 536 */
	 tx329 = tx168*tzz497; 		/* rule 537 */
	 tx330 = tx42*tzz494; 		/* rule 538 */
	 ohz1y3 = tx328 + tx329 + tx330; 		/* rule 539 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 1, ohz1y3);
	 tx331 = tx14*tzz475; 		/* rule 541 */
	 tx332 = -tx313; 		/* rule 542 */
	 tx333 = tx20*tzz476; 		/* rule 543 */
	 tx334 = tx226 + tx331 + tx332 + tx333; 		/* rule 544 */
	 tx335 = tx334*tzz470; 		/* rule 545 */
	 tx336 = tx177*tzz497; 		/* rule 546 */
	 tx337 = tx42*tzz489; 		/* rule 547 */
	 ohz1z3 = tx335 + tx336 + tx337; 		/* rule 548 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 2, ohz1z3);
	 tx338 = tx55*tzz502; 		/* rule 550 */
	 tx339 = tx55*tzz496; 		/* rule 551 */
	 tx340 = tx47*tzz495; 		/* rule 552 */
	 tx341 = tx47*tzz490; 		/* rule 553 */
	 tzz537 = tx5*tx7; 		/* rule 554 */
	 tx342 = tzz505*tzz537; 		/* rule 555 */
	 tx343 = tx184 + tx235 + tx279 + tx338 + tx339 + tx340 + tx341 + tx342; 		/* rule 556 */
	 tx344 = tx343*tzz470; 		/* rule 557 */
	 tx345 = tx136*tzz501; 		/* rule 558 */
	 tx346 = tx51*tzz513; 		/* rule 559 */
	 ohx2y2 = tx344 + tx345 + tx346; 		/* rule 560 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 tx347 = tx63*tzz502; 		/* rule 562 */
	 tx348 = tx63*tzz496; 		/* rule 563 */
	 tx349 = tx47*tzz488; 		/* rule 564 */
	 tx350 = tx47*tzz482; 		/* rule 565 */
	 tzz527 = tx9*tzz505; 		/* rule 566 */
	 tx351 = tx5*tzz527; 		/* rule 567 */
	 tx352 = tx190 + tx242 + tx321 + tx347 + tx348 + tx349 + tx350 + tx351; 		/* rule 568 */
	 tx353 = tx352*tzz470; 		/* rule 569 */
	 tx354 = tx149*tzz501; 		/* rule 570 */
	 tx355 = tx51*tzz509; 		/* rule 571 */
	 ohx2z2 = tx353 + tx354 + tx355; 		/* rule 572 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 tzz507 = -3*tzz480; 		/* rule 574 */
	 tx356 = tx11*tzz507; 		/* rule 575 */
	 tx357 = tx10*tzz468; 		/* rule 576 */
	 tx358 = tx18*tzz473; 		/* rule 577 */
	 tx359 = tx23*tzz473; 		/* rule 578 */
	 tx360 = -tzz539; 		/* rule 579 */
	 tzz543 = tx196 + tx359; 		/* rule 580 */
	 tx361 = tx200 + tx356 + tx357 + tx358 + tx360 + tzz543; 		/* rule 581 */
	 tx362 = tx361*tzz470; 		/* rule 582 */
	 tx363 = tx159*tzz501; 		/* rule 583 */
	 tx364 = tx51*tzz499; 		/* rule 584 */
	 ohx2x3 = tx362 + tx363 + tx364; 		/* rule 585 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 0, ohx2x3);
	 tx365 = -tx341; 		/* rule 587 */
	 tx366 = tx4*tzz495; 		/* rule 588 */
	 tx367 = tx6*tzz496; 		/* rule 589 */
	 tx368 = tzz507*tzz537; 		/* rule 590 */
	 tzz549 = tx366 + tx368; 		/* rule 591 */
	 tx369 = tx209 + tx365 + tx367 + tzz549; 		/* rule 592 */
	 tx370 = tx369*tzz470; 		/* rule 593 */
	 tx371 = tx168*tzz501; 		/* rule 594 */
	 tx372 = tx51*tzz494; 		/* rule 595 */
	 ohx2y3 = tx370 + tx371 + tx372; 		/* rule 596 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 1, ohx2y3);
	 tx373 = -tx350; 		/* rule 598 */
	 tx374 = tx4*tzz488; 		/* rule 599 */
	 tx375 = tx8*tzz496; 		/* rule 600 */
	 tzz526 = tx9*tzz507; 		/* rule 601 */
	 tx376 = tx5*tzz526; 		/* rule 602 */
	 tzz548 = tx374 + tx376; 		/* rule 603 */
	 tx377 = tx219 + tx373 + tx375 + tzz548; 		/* rule 604 */
	 tx378 = tx377*tzz470; 		/* rule 605 */
	 tx379 = tx177*tzz501; 		/* rule 606 */
	 tx380 = tx51*tzz489; 		/* rule 607 */
	 ohx2z3 = tx378 + tx379 + tx380; 		/* rule 608 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 2, ohx2z3);
	 tx381 = tx63*tzz495; 		/* rule 610 */
	 tx382 = tx63*tzz490; 		/* rule 611 */
	 tx383 = tx55*tzz488; 		/* rule 612 */
	 tx384 = tx55*tzz482; 		/* rule 613 */
	 tx385 = tx7*tzz527; 		/* rule 614 */
	 tx386 = tx248 + tx293 + tx326 + tx381 + tx382 + tx383 + tx384 + tx385; 		/* rule 615 */
	 tx387 = tx386*tzz470; 		/* rule 616 */
	 tx388 = tx149*tzz510; 		/* rule 617 */
	 tx389 = tx59*tzz509; 		/* rule 618 */
	 ohy2z2 = tx387 + tx388 + tx389; 		/* rule 619 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	 tx390 = -tx339; 		/* rule 621 */
	 tx391 = tx4*tzz490; 		/* rule 622 */
	 tx392 = tx255 + tx390 + tx391 + tzz549; 		/* rule 623 */
	 tx393 = tx392*tzz470; 		/* rule 624 */
	 tx394 = tx159*tzz510; 		/* rule 625 */
	 tx395 = tx59*tzz499; 		/* rule 626 */
	 ohy2x3 = tx393 + tx394 + tx395; 		/* rule 627 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 0, ohy2x3);
	 tx396 = tx13*tzz507; 		/* rule 629 */
	 tx397 = tx12*tzz468; 		/* rule 630 */
	 tx398 = tx19*tzz473; 		/* rule 631 */
	 tx399 = -tzz536; 		/* rule 632 */
	 tx400 = tx263 + tx396 + tx397 + tx398 + tx399 + tzz543; 		/* rule 633 */
	 tx401 = tx400*tzz470; 		/* rule 634 */
	 tx402 = tx168*tzz510; 		/* rule 635 */
	 tx403 = tx59*tzz494; 		/* rule 636 */
	 ohy2y3 = tx401 + tx402 + tx403; 		/* rule 637 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 1, ohy2y3);
	 tx404 = -tx384; 		/* rule 639 */
	 tx405 = tx6*tzz488; 		/* rule 640 */
	 tx406 = tx8*tzz490; 		/* rule 641 */
	 tx407 = tx7*tzz526; 		/* rule 642 */
	 tzz547 = tx405 + tx407; 		/* rule 643 */
	 tx408 = tx273 + tx404 + tx406 + tzz547; 		/* rule 644 */
	 tx409 = tx408*tzz470; 		/* rule 645 */
	 tx410 = tx177*tzz510; 		/* rule 646 */
	 tx411 = tx59*tzz489; 		/* rule 647 */
	 ohy2z3 = tx409 + tx410 + tx411; 		/* rule 648 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 2, ohy2z3);
	 tx412 = -tx348; 		/* rule 650 */
	 tx413 = tx4*tzz482; 		/* rule 651 */
	 tx414 = tx299 + tx412 + tx413 + tzz548; 		/* rule 652 */
	 tx415 = tx414*tzz470; 		/* rule 653 */
	 tx416 = tx159*tzz512; 		/* rule 654 */
	 tx417 = tx67*tzz499; 		/* rule 655 */
	 ohz2x3 = tx415 + tx416 + tx417; 		/* rule 656 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 0, ohz2x3);
	 tx418 = -tx382; 		/* rule 658 */
	 tx419 = tx6*tzz482; 		/* rule 659 */
	 tx420 = tx306 + tx418 + tx419 + tzz547; 		/* rule 660 */
	 tx421 = tx420*tzz470; 		/* rule 661 */
	 tx422 = tx168*tzz512; 		/* rule 662 */
	 tx423 = tx67*tzz494; 		/* rule 663 */
	 ohz2y3 = tx421 + tx422 + tx423; 		/* rule 664 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 1, ohz2y3);
	 tx424 = tx15*tzz507; 		/* rule 666 */
	 tx425 = tx14*tzz468; 		/* rule 667 */
	 tx426 = tx20*tzz473; 		/* rule 668 */
	 tx427 = -tzz535; 		/* rule 669 */
	 tx428 = tx314 + tx424 + tx425 + tx426 + tx427 + tzz543; 		/* rule 670 */
	 tx429 = tx428*tzz470; 		/* rule 671 */
	 tx430 = tx177*tzz512; 		/* rule 672 */
	 tx431 = tx67*tzz489; 		/* rule 673 */
	 ohz2z3 = tx429 + tx430 + tx431; 		/* rule 674 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 2, ohz2z3);
	 tx432 = -tx367; 		/* rule 676 */
	 tx433 = -tx391; 		/* rule 677 */
	 tx434 = tx342 + tx432 + tx433; 		/* rule 678 */
	 tx435 = tx434*tzz470; 		/* rule 679 */
	 tx436 = tx168*tzz514; 		/* rule 680 */
	 tx437 = tx70*tzz494; 		/* rule 681 */
	 ohx3y3 = tx435 + tx436 + tx437; 		/* rule 682 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 1, ohx3y3);
	 tx438 = -tx375; 		/* rule 684 */
	 tx439 = -tx413; 		/* rule 685 */
	 tx440 = tx351 + tx438 + tx439; 		/* rule 686 */
	 tx441 = tx440*tzz470; 		/* rule 687 */
	 tx442 = tx177*tzz514; 		/* rule 688 */
	 tx443 = tx70*tzz489; 		/* rule 689 */
	 ohx3z3 = tx441 + tx442 + tx443; 		/* rule 690 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 2, ohx3z3);
	 tx444 = -tx406; 		/* rule 692 */
	 tx445 = -tx419; 		/* rule 693 */
	 tx446 = tx385 + tx444 + tx445; 		/* rule 694 */
	 tx447 = tx446*tzz470; 		/* rule 695 */
	 tx448 = tx177*tzz534; 		/* rule 696 */
	 tx449 = tx73*tzz489; 		/* rule 697 */
	 ohy3z3 = tx447 + tx448 + tx449; 		/* rule 698 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 2, ohy3z3);
	} /*if calcOffDiagonalHessian */ 
	#endif /* ANGLE_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* ANGLE_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* ANGLE_CALC_FORCE ]*/
	ANGLE_DEBUG_INTERACTIONS(I1,I2,I3);
	SKIP_term: (void)0;
