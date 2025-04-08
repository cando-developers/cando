// C-code
	CHIRAL_RESTRAINT_SET_PARAMETER(K);
	CHIRAL_RESTRAINT_SET_PARAMETER(CO);
	CHIRAL_RESTRAINT_SET_PARAMETER(I1);
	CHIRAL_RESTRAINT_SET_PARAMETER(I2);
	CHIRAL_RESTRAINT_SET_PARAMETER(I3);
	CHIRAL_RESTRAINT_SET_PARAMETER(I4);
	CHIRAL_RESTRAINT_APPLY_ATOM_MASK(I1,I2,I3,I4);
	CHIRAL_RESTRAINT_SET_POSITION(x1,I1,0);
	CHIRAL_RESTRAINT_SET_POSITION(y1,I1,1);
	CHIRAL_RESTRAINT_SET_POSITION(z1,I1,2);
	CHIRAL_RESTRAINT_SET_POSITION(x2,I2,0);
	CHIRAL_RESTRAINT_SET_POSITION(y2,I2,1);
	CHIRAL_RESTRAINT_SET_POSITION(z2,I2,2);
	CHIRAL_RESTRAINT_SET_POSITION(x3,I3,0);
	CHIRAL_RESTRAINT_SET_POSITION(y3,I3,1);
	CHIRAL_RESTRAINT_SET_POSITION(z3,I3,2);
	CHIRAL_RESTRAINT_SET_POSITION(x4,I4,0);
	CHIRAL_RESTRAINT_SET_POSITION(y4,I4,1);
	CHIRAL_RESTRAINT_SET_POSITION(z4,I4,2);
	 tx1 = -x3; 		/* rule 20 */
	 tzz931 = -y3; 		/* rule 21 */
	 tx2 = tzz931; 		/* rule 22 */
	 tzz897 = -z3; 		/* rule 23 */
	 tx3 = tzz897; 		/* rule 24 */
	 tx4 = tx2*x1; 		/* rule 25 */
	 tx5 = tx1 + x1; 		/* rule 26 */
	 tx6 = tx3*x2; 		/* rule 27 */
	 tx7 = tx1 + x2; 		/* rule 28 */
	 tx8 = tx1 + x4; 		/* rule 29 */
	 tx9 = tx3*y1; 		/* rule 30 */
	 tzz933 = -y1; 		/* rule 31 */
	 tx10 = tzz933*x2; 		/* rule 32 */
	 tx11 = x3*y1; 		/* rule 33 */
	 tx12 = tx2 + y1; 		/* rule 34 */
	 tx13 = tx1*y2; 		/* rule 35 */
	 tx14 = x1*y2; 		/* rule 36 */
	 tx15 = tx2 + y2; 		/* rule 37 */
	 tx16 = x2*y3; 		/* rule 38 */
	 tx17 = tx2 + y4; 		/* rule 39 */
	 tx18 = tx1*z1; 		/* rule 40 */
	 tx19 = x2*z1; 		/* rule 41 */
	 tzz899 = -z1; 		/* rule 42 */
	 tx20 = tzz899*y2; 		/* rule 43 */
	 tx21 = y3*z1; 		/* rule 44 */
	 tx22 = tx3 + z1; 		/* rule 45 */
	 tx23 = tx2*z2; 		/* rule 46 */
	 tzz898 = -z2; 		/* rule 47 */
	 tx24 = tzz898*x1; 		/* rule 48 */
	 tx25 = x3*z2; 		/* rule 49 */
	 tx26 = y1*z2; 		/* rule 50 */
	 tx27 = tx3 + z2; 		/* rule 51 */
	 tx28 = x1*z3; 		/* rule 52 */
	 tx29 = y2*z3; 		/* rule 53 */
	 tx30 = tx3 + z4; 		/* rule 54 */
	 tx31 = power2(tx12); 		/* rule 55 */
	 tx32 = power2(tx15); 		/* rule 56 */
	 tx33 = power2(tx17); 		/* rule 57 */
	 tx34 = power2(tx22); 		/* rule 58 */
	 tx35 = power2(tx27); 		/* rule 59 */
	 tx36 = power2(tx30); 		/* rule 60 */
	 tzz942 = tx14 + tx16; 		/* rule 61 */
	 tzz943 = tx11 + tzz942; 		/* rule 62 */
	 tzz944 = tx10 + tzz943; 		/* rule 63 */
	 tx37 = tx13 + tx4 + tzz944; 		/* rule 64 */
	 tx38 = power2(tx5); 		/* rule 65 */
	 tzz937 = tx25 + tx28; 		/* rule 66 */
	 tzz938 = tx24 + tzz937; 		/* rule 67 */
	 tzz941 = tx19 + tzz938; 		/* rule 68 */
	 tx39 = tx18 + tx6 + tzz941; 		/* rule 69 */
	 tx40 = power2(tx7); 		/* rule 70 */
	 tx41 = power2(tx8); 		/* rule 71 */
	 tzz936 = tx26 + tx29; 		/* rule 72 */
	 tzz939 = tx21 + tzz936; 		/* rule 73 */
	 tzz940 = tx20 + tzz939; 		/* rule 74 */
	 tx42 = tx23 + tx9 + tzz940; 		/* rule 75 */
	 tx43 = tx30*tx37; 		/* rule 76 */
	 tx44 = tx31 + tx34 + tx38; 		/* rule 77 */
	 tx45 = tx17*tx39; 		/* rule 78 */
	 tx46 = tx32 + tx35 + tx40; 		/* rule 79 */
	 tx47 = tx33 + tx36 + tx41; 		/* rule 80 */
	 tx48 = tx42*tx8; 		/* rule 81 */
	 tx804 = mysqrt(tx44); 		/* rule 82 */
	 tx49 = reciprocal(tx804); 		/* rule 83 */
	 tx805 = mysqrt(tx46); 		/* rule 84 */
	 tx50 = reciprocal(tx805); 		/* rule 85 */
	 tx806 = mysqrt(tx47); 		/* rule 86 */
	 tx51 = reciprocal(tx806); 		/* rule 87 */
	 tx52 = tx43 + tx45 + tx48; 		/* rule 88 */
	 tzz874 = tx50*tx51; 		/* rule 89 */
	 tx53 = tx49*tx52*tzz874; 		/* rule 90 */
	 tx54 = CO + tx53; 		/* rule 91 */
	 tx807 = power2(tx54); 		/* rule 92 */
	 tx55 = tx54*tx807; 		/* rule 93 */
	 ChiralTest = K*tx55; 		/* rule 94 */
/*If[ChiralTest>0.0,*/
	if (ChiralTest>0.0 ) {
	Energy=ChiralTest;
/*Energy = ChiralTest;*/
	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE(Energy);
	#ifdef CHIRAL_RESTRAINT_CALC_FORCE //[
	if ( calcForce ) {
	 tzz932 = -y2; 		/* rule 102 */
	 tx56 = tzz932*x3; 		/* rule 103 */
	 tx57 = tx4; 		/* rule 104 */
	 tx58 = tzz899*x3; 		/* rule 105 */
	 tx59 = tzz898*y3; 		/* rule 106 */
	 tx60 = tx6; 		/* rule 107 */
	 tx61 = tx9; 		/* rule 108 */
	 tx62 = tx56 + tx57 + tzz944; 		/* rule 109 */
	 tx63 = tx58 + tx60 + tzz941; 		/* rule 110 */
	 tx64 = tx59 + tx61 + tzz940; 		/* rule 111 */
	 tx65 = tzz898; 		/* rule 112 */
	 tx66 = tx30*tx62; 		/* rule 113 */
	 tx67 = tx17*tx63; 		/* rule 114 */
	 tx68 = tx64*tx8; 		/* rule 115 */
	 tx69 = tx65 + z3; 		/* rule 116 */
	 tx70 = tx15*tx30; 		/* rule 117 */
	 tx71 = tx66 + tx67 + tx68; 		/* rule 118 */
	 tx72 = tx17*tx69; 		/* rule 119 */
	 tx808 = tx44; 		/* rule 120 */
	 tx809 = tx49; 		/* rule 121 */
	 tx810 = reciprocal(tx808); 		/* rule 122 */
	 tx73 = tx809*tx810; 		/* rule 123 */
	 tzz886 = tx809*tzz874; 		/* rule 124 */
	 tx74 = tx71*tzz886; 		/* rule 125 */
	 tx75 = tx70 + tx72; 		/* rule 126 */
	 tzz881 = -tx71; 		/* rule 127 */
	 tzz900 = tzz874*tzz881; 		/* rule 128 */
	 tzz904 = tx73*tzz900; 		/* rule 129 */
	 tx76 = tx5*tzz904; 		/* rule 130 */
	 tx77 = CO + tx74; 		/* rule 131 */
	 tx78 = tx75*tzz886; 		/* rule 132 */
	 tx79 = power2(tx77); 		/* rule 133 */
	 tx80 = tx76 + tx78; 		/* rule 134 */
	 tzz826 = K*tx79; 		/* rule 135 */
	 tzz827 = 3*tzz826; 		/* rule 136 */
	 gx1 = tx80*tzz827; 		/* rule 137 */
	 fx1 = -gx1; 		/* rule 138 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx81 = -x2; 		/* rule 140 */
	 tx82 = tx81 + x3; 		/* rule 141 */
	 tx83 = tx27*tx8; 		/* rule 142 */
	 tx84 = tx30*tx82; 		/* rule 143 */
	 tx85 = tx83 + tx84; 		/* rule 144 */
	 tx86 = tx12*tzz904; 		/* rule 145 */
	 tx87 = tx85*tzz886; 		/* rule 146 */
	 tx88 = tx86 + tx87; 		/* rule 147 */
	 gy1 = tx88*tzz827; 		/* rule 148 */
	 fy1 = -gy1; 		/* rule 149 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx89 = tzz932; 		/* rule 151 */
	 tx90 = tx89 + y3; 		/* rule 152 */
	 tx91 = tx17*tx7; 		/* rule 153 */
	 tx92 = tx8*tx90; 		/* rule 154 */
	 tx93 = tx91 + tx92; 		/* rule 155 */
	 tx94 = tx22*tzz904; 		/* rule 156 */
	 tx95 = tx93*tzz886; 		/* rule 157 */
	 tx96 = tx94 + tx95; 		/* rule 158 */
	 gz1 = tx96*tzz827; 		/* rule 159 */
	 fz1 = -gz1; 		/* rule 160 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tx97 = tzz933; 		/* rule 162 */
	 tx98 = tx97 + y3; 		/* rule 163 */
	 tx99 = tx17*tx22; 		/* rule 164 */
	 tx100 = tx30*tx98; 		/* rule 165 */
	 tx811 = tx46; 		/* rule 166 */
	 tx812 = tx50; 		/* rule 167 */
	 tx813 = reciprocal(tx811); 		/* rule 168 */
	 tx101 = tx812*tx813; 		/* rule 169 */
	 tx102 = tx100 + tx99; 		/* rule 170 */
	 tzz823 = tx809*tx812; 		/* rule 171 */
	 tzz907 = tx51*tzz823; 		/* rule 172 */
	 tx103 = tx102*tzz907; 		/* rule 173 */
	 tzz824 = tx101*tx809; 		/* rule 174 */
	 tzz901 = tzz824*tzz881; 		/* rule 175 */
	 tzz906 = tx51*tzz901; 		/* rule 176 */
	 tx104 = tx7*tzz906; 		/* rule 177 */
	 tx105 = tx103 + tx104; 		/* rule 178 */
	 gx2 = tx105*tzz827; 		/* rule 179 */
	 fx2 = -gx2; 		/* rule 180 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I2, 0, fx2 );
	 tx106 = tzz899; 		/* rule 182 */
	 tx107 = tx106 + z3; 		/* rule 183 */
	 tx108 = tx30*tx5; 		/* rule 184 */
	 tx109 = tx107*tx8; 		/* rule 185 */
	 tx110 = tx108 + tx109; 		/* rule 186 */
	 tx111 = tx110*tzz907; 		/* rule 187 */
	 tx112 = tx15*tzz906; 		/* rule 188 */
	 tx113 = tx111 + tx112; 		/* rule 189 */
	 gy2 = tx113*tzz827; 		/* rule 190 */
	 fy2 = -gy2; 		/* rule 191 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I2, 1, fy2 );
	 tx114 = -x1; 		/* rule 193 */
	 tx115 = tx114 + x3; 		/* rule 194 */
	 tx116 = tx115*tx17; 		/* rule 195 */
	 tx117 = tx12*tx8; 		/* rule 196 */
	 tx118 = tx116 + tx117; 		/* rule 197 */
	 tx119 = tx118*tzz907; 		/* rule 198 */
	 tx120 = tx27*tzz906; 		/* rule 199 */
	 tx121 = tx119 + tx120; 		/* rule 200 */
	 gz2 = tx121*tzz827; 		/* rule 201 */
	 fz2 = -gz2; 		/* rule 202 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I2, 2, fz2 );
	 tx122 = tx89 + y1; 		/* rule 204 */
	 tx123 = tx106 + z2; 		/* rule 205 */
	 tx124 = tx123*tx17; 		/* rule 206 */
	 tx125 = tx122*tx30; 		/* rule 207 */
	 tx126 = y2*z1; 		/* rule 208 */
	 tx127 = tx106*y3; 		/* rule 209 */
	 tx128 = tx65*y1; 		/* rule 210 */
	 tx129 = y3*z2; 		/* rule 211 */
	 tx130 = y1*z3; 		/* rule 212 */
	 tx131 = tx3*y2; 		/* rule 213 */
	 tx132 = tx124 + tx125 + tx126 + tx127 + tx128 + tx129 + tx130 + tx131; 		/* rule 214 */
	 tx814 = tx47; 		/* rule 215 */
	 tx815 = tx51; 		/* rule 216 */
	 tx816 = reciprocal(tx814); 		/* rule 217 */
	 tx133 = tx815*tx816; 		/* rule 218 */
	 tzz833 = tx815*tzz823; 		/* rule 219 */
	 tx134 = tx132*tzz833; 		/* rule 220 */
	 tzz831 = tx815*tzz824; 		/* rule 221 */
	 tzz845 = tx7*tzz831; 		/* rule 222 */
	 tx135 = tx71*tzz845; 		/* rule 223 */
	 tzz825 = tx73*tx815; 		/* rule 224 */
	 tzz832 = tx812*tzz825; 		/* rule 225 */
	 tzz846 = tx5*tzz832; 		/* rule 226 */
	 tx136 = tx71*tzz846; 		/* rule 227 */
	 tzz828 = tx133*tzz823; 		/* rule 228 */
	 tzz844 = tx8*tzz828; 		/* rule 229 */
	 tx137 = tx71*tzz844; 		/* rule 230 */
	 tx138 = tx134 + tx135 + tx136 + tx137; 		/* rule 231 */
	 gx3 = tx138*tzz827; 		/* rule 232 */
	 fx3 = -gx3; 		/* rule 233 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I3, 0, fx3 );
	 tx139 = tx114 + x2; 		/* rule 235 */
	 tx140 = tx65 + z1; 		/* rule 236 */
	 tx141 = tx139*tx30; 		/* rule 237 */
	 tx142 = tx140*tx8; 		/* rule 238 */
	 tx143 = tx106*x2; 		/* rule 239 */
	 tx144 = x3*z1; 		/* rule 240 */
	 tx145 = x1*z2; 		/* rule 241 */
	 tx146 = tx65*x3; 		/* rule 242 */
	 tx147 = tx3*x1; 		/* rule 243 */
	 tx148 = x2*z3; 		/* rule 244 */
	 tx149 = tx141 + tx142 + tx143 + tx144 + tx145 + tx146 + tx147 + tx148; 		/* rule 245 */
	 tx150 = tx149*tzz833; 		/* rule 246 */
	 tzz850 = tx17*tzz828; 		/* rule 247 */
	 tx151 = tx71*tzz850; 		/* rule 248 */
	 tzz851 = tx15*tzz831; 		/* rule 249 */
	 tx152 = tx71*tzz851; 		/* rule 250 */
	 tzz852 = tx12*tzz832; 		/* rule 251 */
	 tx153 = tx71*tzz852; 		/* rule 252 */
	 tx154 = tx150 + tx151 + tx152 + tx153; 		/* rule 253 */
	 gy3 = tx154*tzz827; 		/* rule 254 */
	 fy3 = -gy3; 		/* rule 255 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I3, 1, fy3 );
	 tx155 = tx81 + x1; 		/* rule 257 */
	 tx156 = tx97 + y2; 		/* rule 258 */
	 tx157 = tx155*tx17; 		/* rule 259 */
	 tx158 = tx156*tx8; 		/* rule 260 */
	 tx159 = x2*y1; 		/* rule 261 */
	 tx160 = tx97*x3; 		/* rule 262 */
	 tx161 = tx89*x1; 		/* rule 263 */
	 tx162 = x3*y2; 		/* rule 264 */
	 tx163 = x1*y3; 		/* rule 265 */
	 tx164 = tx2*x2; 		/* rule 266 */
	 tx165 = tx157 + tx158 + tx159 + tx160 + tx161 + tx162 + tx163 + tx164; 		/* rule 267 */
	 tx166 = tx165*tzz833; 		/* rule 268 */
	 tzz847 = tx30*tzz828; 		/* rule 269 */
	 tx167 = tx71*tzz847; 		/* rule 270 */
	 tzz848 = tx27*tzz831; 		/* rule 271 */
	 tx168 = tx71*tzz848; 		/* rule 272 */
	 tzz849 = tx22*tzz832; 		/* rule 273 */
	 tx169 = tx71*tzz849; 		/* rule 274 */
	 tx170 = tx166 + tx167 + tx168 + tx169; 		/* rule 275 */
	 gz3 = tx170*tzz827; 		/* rule 276 */
	 fz3 = -gz3; 		/* rule 277 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I3, 2, fz3 );
	 tx171 = tx64*tzz833; 		/* rule 279 */
	 tzz835 = -tzz828; 		/* rule 280 */
	 tzz855 = tx8*tzz835; 		/* rule 281 */
	 tx172 = tx71*tzz855; 		/* rule 282 */
	 tx173 = tx171 + tx172; 		/* rule 283 */
	 gx4 = tx173*tzz827; 		/* rule 284 */
	 fx4 = -gx4; 		/* rule 285 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I4, 0, fx4 );
	 tx174 = tx63*tzz833; 		/* rule 287 */
	 tzz861 = tx17*tzz835; 		/* rule 288 */
	 tx175 = tx71*tzz861; 		/* rule 289 */
	 tx176 = tx174 + tx175; 		/* rule 290 */
	 gy4 = tx176*tzz827; 		/* rule 291 */
	 fy4 = -gy4; 		/* rule 292 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I4, 1, fy4 );
	 tx177 = tx62*tzz833; 		/* rule 294 */
	 tzz857 = tx30*tzz835; 		/* rule 295 */
	 tx178 = tx71*tzz857; 		/* rule 296 */
	 tx179 = tx177 + tx178; 		/* rule 297 */
	 gz4 = tx179*tzz827; 		/* rule 298 */
	 fz4 = -gz4; 		/* rule 299 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I4, 2, fz4 );
	#ifdef CHIRAL_RESTRAINT_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx817 = power2(tx808); 		/* rule 303 */
	 tx818 = reciprocal(tx817); 		/* rule 304 */
	 tx180 = tx809*tx818; 		/* rule 305 */
	 tzz838 = tx71*tx815; 		/* rule 306 */
	 tzz853 = tx812*tzz838; 		/* rule 307 */
	 tzz860 = tx180*tzz853; 		/* rule 308 */
	 tzz878 = 3*tzz860; 		/* rule 309 */
	 tx181 = tx38*tzz878; 		/* rule 310 */
	 tzz836 = -tzz832; 		/* rule 311 */
	 tx182 = tx71*tzz836; 		/* rule 312 */
	 tx183 = -2.*tx75*tzz846; 		/* rule 313 */
	 tx184 = tx181 + tx182 + tx183; 		/* rule 314 */
	 tx185 = power2(tx80); 		/* rule 315 */
	 tzz829 = K*tx77; 		/* rule 316 */
	 tzz830 = 6.*tzz829; 		/* rule 317 */
	 tx186 = tx185*tzz830; 		/* rule 318 */
	 tx187 = tx184*tzz827; 		/* rule 319 */
	 dhx1x1 = tx186 + tx187; 		/* rule 320 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx188 = tx31*tzz878; 		/* rule 322 */
	 tx189 = -2.*tx85*tzz852; 		/* rule 323 */
	 tx190 = tx182 + tx188 + tx189; 		/* rule 324 */
	 tx191 = power2(tx88); 		/* rule 325 */
	 tx192 = tx191*tzz830; 		/* rule 326 */
	 tx193 = tx190*tzz827; 		/* rule 327 */
	 dhy1y1 = tx192 + tx193; 		/* rule 328 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx194 = tx34*tzz878; 		/* rule 330 */
	 tx195 = -2.*tx93*tzz849; 		/* rule 331 */
	 tx196 = tx182 + tx194 + tx195; 		/* rule 332 */
	 tx197 = power2(tx96); 		/* rule 333 */
	 tx198 = tx197*tzz830; 		/* rule 334 */
	 tx199 = tx196*tzz827; 		/* rule 335 */
	 dhz1z1 = tx198 + tx199; 		/* rule 336 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 tx819 = power2(tx811); 		/* rule 338 */
	 tx820 = reciprocal(tx819); 		/* rule 339 */
	 tx200 = tx812*tx820; 		/* rule 340 */
	 tx201 = -2.*tx102*tzz845; 		/* rule 341 */
	 tzz837 = -tzz831; 		/* rule 342 */
	 tx202 = tx71*tzz837; 		/* rule 343 */
	 tzz854 = tx809*tzz838; 		/* rule 344 */
	 tzz859 = tx200*tzz854; 		/* rule 345 */
	 tzz879 = 3*tzz859; 		/* rule 346 */
	 tx203 = tx40*tzz879; 		/* rule 347 */
	 tx204 = power2(tx105); 		/* rule 348 */
	 tx205 = tx201 + tx202 + tx203; 		/* rule 349 */
	 tx206 = tx204*tzz830; 		/* rule 350 */
	 tx207 = tx205*tzz827; 		/* rule 351 */
	 dhx2x2 = tx206 + tx207; 		/* rule 352 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 tx208 = -2.*tx110*tzz851; 		/* rule 354 */
	 tx209 = tx32*tzz879; 		/* rule 355 */
	 tx210 = power2(tx113); 		/* rule 356 */
	 tx211 = tx202 + tx208 + tx209; 		/* rule 357 */
	 tx212 = tx210*tzz830; 		/* rule 358 */
	 tx213 = tx211*tzz827; 		/* rule 359 */
	 dhy2y2 = tx212 + tx213; 		/* rule 360 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 tx214 = -2.*tx118*tzz848; 		/* rule 362 */
	 tx215 = tx35*tzz879; 		/* rule 363 */
	 tx216 = power2(tx121); 		/* rule 364 */
	 tx217 = tx202 + tx214 + tx215; 		/* rule 365 */
	 tx218 = tx216*tzz830; 		/* rule 366 */
	 tx219 = tx217*tzz827; 		/* rule 367 */
	 dhz2z2 = tx218 + tx219; 		/* rule 368 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	 tx821 = power2(tx814); 		/* rule 370 */
	 tx822 = reciprocal(tx821); 		/* rule 371 */
	 tx220 = tx815*tx822; 		/* rule 372 */
	 tzz917 = 2.*tx132; 		/* rule 373 */
	 tx221 = tzz845*tzz917; 		/* rule 374 */
	 tx222 = tx71*tzz835; 		/* rule 375 */
	 tzz856 = tx71*tzz823; 		/* rule 376 */
	 tzz858 = tx220*tzz856; 		/* rule 377 */
	 tzz880 = 3*tzz858; 		/* rule 378 */
	 tx223 = tx41*tzz880; 		/* rule 379 */
	 tx224 = tzz846*tzz917; 		/* rule 380 */
	 tzz842 = tx71*tzz825; 		/* rule 381 */
	 tzz843 = tx101*tzz842; 		/* rule 382 */
	 tzz889 = tx7*tzz843; 		/* rule 383 */
	 tzz928 = 2.*tx5; 		/* rule 384 */
	 tx225 = tzz889*tzz928; 		/* rule 385 */
	 tx226 = tzz844*tzz917; 		/* rule 386 */
	 tzz834 = tx133*tx71; 		/* rule 387 */
	 tzz839 = tzz824*tzz834; 		/* rule 388 */
	 tzz888 = tx8*tzz839; 		/* rule 389 */
	 tx227 = 2.*tx7*tzz888; 		/* rule 390 */
	 tzz840 = tx812*tzz834; 		/* rule 391 */
	 tzz841 = tx73*tzz840; 		/* rule 392 */
	 tzz887 = tx8*tzz841; 		/* rule 393 */
	 tx228 = tzz887*tzz928; 		/* rule 394 */
	 tx229 = power2(tx138); 		/* rule 395 */
	 tzz934 = tx202 + tx222; 		/* rule 396 */
	 tzz935 = tx182 + tzz934; 		/* rule 397 */
	 tx230 = tx181 + tx203 + tx221 + tx223 + tx224 + tx225 + tx226 + tx227 + tx228 + tzz935; 		/* rule 398 */
	 tx231 = tx229*tzz830; 		/* rule 399 */
	 tx232 = tx230*tzz827; 		/* rule 400 */
	 dhx3x3 = tx231 + tx232; 		/* rule 401 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 0, dhx3x3);
	 tzz916 = 2.*tx149; 		/* rule 403 */
	 tx233 = tzz850*tzz916; 		/* rule 404 */
	 tx234 = tzz851*tzz916; 		/* rule 405 */
	 tzz894 = tx17*tzz839; 		/* rule 406 */
	 tx235 = 2.*tx15*tzz894; 		/* rule 407 */
	 tx236 = tx33*tzz880; 		/* rule 408 */
	 tx237 = tzz852*tzz916; 		/* rule 409 */
	 tzz893 = tx17*tzz841; 		/* rule 410 */
	 tzz930 = 2.*tx12; 		/* rule 411 */
	 tx238 = tzz893*tzz930; 		/* rule 412 */
	 tzz896 = tx15*tzz843; 		/* rule 413 */
	 tx239 = tzz896*tzz930; 		/* rule 414 */
	 tx240 = power2(tx154); 		/* rule 415 */
	 tx241 = tx188 + tx209 + tx233 + tx234 + tx235 + tx236 + tx237 + tx238 + tx239 + tzz935; 		/* rule 416 */
	 tx242 = tx240*tzz830; 		/* rule 417 */
	 tx243 = tx241*tzz827; 		/* rule 418 */
	 dhy3y3 = tx242 + tx243; 		/* rule 419 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 1, dhy3y3);
	 tzz915 = 2.*tx165; 		/* rule 421 */
	 tx244 = tzz847*tzz915; 		/* rule 422 */
	 tx245 = tzz848*tzz915; 		/* rule 423 */
	 tzz891 = tx30*tzz839; 		/* rule 424 */
	 tx246 = 2.*tx27*tzz891; 		/* rule 425 */
	 tx247 = tx36*tzz880; 		/* rule 426 */
	 tx248 = tzz849*tzz915; 		/* rule 427 */
	 tzz890 = tx30*tzz841; 		/* rule 428 */
	 tzz929 = 2.*tx22; 		/* rule 429 */
	 tx249 = tzz890*tzz929; 		/* rule 430 */
	 tzz892 = tx27*tzz843; 		/* rule 431 */
	 tx250 = tzz892*tzz929; 		/* rule 432 */
	 tx251 = power2(tx170); 		/* rule 433 */
	 tx252 = tx194 + tx215 + tx244 + tx245 + tx246 + tx247 + tx248 + tx249 + tx250 + tzz935; 		/* rule 434 */
	 tx253 = tx251*tzz830; 		/* rule 435 */
	 tx254 = tx252*tzz827; 		/* rule 436 */
	 dhz3z3 = tx253 + tx254; 		/* rule 437 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I3, 2, dhz3z3);
	 tx255 = -2.*tx64*tzz844; 		/* rule 439 */
	 tx256 = power2(tx173); 		/* rule 440 */
	 tx257 = tx222 + tx223 + tx255; 		/* rule 441 */
	 tx258 = tx256*tzz830; 		/* rule 442 */
	 tx259 = tx257*tzz827; 		/* rule 443 */
	 dhx4x4 = tx258 + tx259; 		/* rule 444 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 0, dhx4x4);
	 tx260 = -2.*tx63*tzz850; 		/* rule 446 */
	 tx261 = power2(tx176); 		/* rule 447 */
	 tx262 = tx222 + tx236 + tx260; 		/* rule 448 */
	 tx263 = tx261*tzz830; 		/* rule 449 */
	 tx264 = tx262*tzz827; 		/* rule 450 */
	 dhy4y4 = tx263 + tx264; 		/* rule 451 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I4, 1, I4, 1, dhy4y4);
	 tx265 = -2.*tx62*tzz847; 		/* rule 453 */
	 tx266 = power2(tx179); 		/* rule 454 */
	 tx267 = tx222 + tx247 + tx265; 		/* rule 455 */
	 tx268 = tx266*tzz830; 		/* rule 456 */
	 tx269 = tx267*tzz827; 		/* rule 457 */
	 dhz4z4 = tx268 + tx269; 		/* rule 458 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I4, 2, I4, 2, dhz4z4);
	#ifdef CHIRAL_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tzz923 = tx5*tzz878; 		/* rule 462 */
	 tx270 = tx12*tzz923; 		/* rule 463 */
	 tzz868 = tx12*tzz836; 		/* rule 464 */
	 tx271 = tx75*tzz868; 		/* rule 465 */
	 tzz864 = tx5*tzz836; 		/* rule 466 */
	 tx272 = tx85*tzz864; 		/* rule 467 */
	 tx273 = tx270 + tx271 + tx272; 		/* rule 468 */
	 tx274 = tx273*tzz827; 		/* rule 469 */
	 tzz869 = tx88*tzz830; 		/* rule 470 */
	 tx275 = tx80*tzz869; 		/* rule 471 */
	 ohx1y1 = tx274 + tx275; 		/* rule 472 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tx276 = tx22*tzz923; 		/* rule 474 */
	 tzz866 = tx22*tzz836; 		/* rule 475 */
	 tx277 = tx75*tzz866; 		/* rule 476 */
	 tx278 = tx93*tzz864; 		/* rule 477 */
	 tx279 = tx276 + tx277 + tx278; 		/* rule 478 */
	 tx280 = tx279*tzz827; 		/* rule 479 */
	 tzz862 = tx96*tzz830; 		/* rule 480 */
	 tx281 = tx80*tzz862; 		/* rule 481 */
	 ohx1z1 = tx280 + tx281; 		/* rule 482 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx282 = tx102*tzz864; 		/* rule 484 */
	 tx283 = tx5*tzz889; 		/* rule 485 */
	 tzz863 = tx7*tzz837; 		/* rule 486 */
	 tx284 = tx75*tzz863; 		/* rule 487 */
	 tx285 = tx282 + tx283 + tx284; 		/* rule 488 */
	 tx286 = tx285*tzz827; 		/* rule 489 */
	 tzz870 = tx80*tzz830; 		/* rule 490 */
	 tx287 = tx105*tzz870; 		/* rule 491 */
	 ohx1x2 = tx286 + tx287; 		/* rule 492 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx288 = tx30*tzz833; 		/* rule 494 */
	 tx289 = tx110*tzz864; 		/* rule 495 */
	 tx290 = tx5*tzz896; 		/* rule 496 */
	 tzz867 = tx15*tzz837; 		/* rule 497 */
	 tx291 = tx75*tzz867; 		/* rule 498 */
	 tx292 = tx288 + tx289 + tx290 + tx291; 		/* rule 499 */
	 tx293 = tx292*tzz827; 		/* rule 500 */
	 tx294 = tx113*tzz870; 		/* rule 501 */
	 ohx1y2 = tx293 + tx294; 		/* rule 502 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx295 = -y4; 		/* rule 504 */
	 tx296 = tx295 + y3; 		/* rule 505 */
	 tx297 = tx296*tzz833; 		/* rule 506 */
	 tx298 = tx118*tzz864; 		/* rule 507 */
	 tx299 = tx5*tzz892; 		/* rule 508 */
	 tzz865 = tx27*tzz837; 		/* rule 509 */
	 tx300 = tx75*tzz865; 		/* rule 510 */
	 tx301 = tx297 + tx298 + tx299 + tx300; 		/* rule 511 */
	 tx302 = tx301*tzz827; 		/* rule 512 */
	 tx303 = tx121*tzz870; 		/* rule 513 */
	 ohx1z2 = tx302 + tx303; 		/* rule 514 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tzz882 = -3*tzz860; 		/* rule 516 */
	 tx304 = tx38*tzz882; 		/* rule 517 */
	 tx305 = tx132*tzz864; 		/* rule 518 */
	 tx306 = tx71*tzz832; 		/* rule 519 */
	 tzz871 = -tzz843; 		/* rule 520 */
	 tzz905 = tx7*tzz871; 		/* rule 521 */
	 tx307 = tx5*tzz905; 		/* rule 522 */
	 tx308 = tx75*tzz845; 		/* rule 523 */
	 tx309 = tx75*tzz846; 		/* rule 524 */
	 tzz872 = -tzz841; 		/* rule 525 */
	 tzz903 = tx8*tzz872; 		/* rule 526 */
	 tx310 = tx5*tzz903; 		/* rule 527 */
	 tx311 = tx75*tzz844; 		/* rule 528 */
	 tx312 = tx304 + tx305 + tx306 + tx307 + tx308 + tx309 + tx310 + tx311; 		/* rule 529 */
	 tx313 = tx312*tzz827; 		/* rule 530 */
	 tx314 = tx138*tzz870; 		/* rule 531 */
	 ohx1x3 = tx313 + tx314; 		/* rule 532 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 0, ohx1x3);
	 tx315 = -z4; 		/* rule 534 */
	 tx316 = tx315 + z2; 		/* rule 535 */
	 tx317 = tx316*tzz833; 		/* rule 536 */
	 tzz922 = tx5*tzz882; 		/* rule 537 */
	 tx318 = tx12*tzz922; 		/* rule 538 */
	 tx319 = tx149*tzz864; 		/* rule 539 */
	 tzz912 = tx17*tzz872; 		/* rule 540 */
	 tx320 = tx5*tzz912; 		/* rule 541 */
	 tzz913 = tx15*tzz871; 		/* rule 542 */
	 tx321 = tx5*tzz913; 		/* rule 543 */
	 tx322 = tx75*tzz850; 		/* rule 544 */
	 tx323 = tx75*tzz851; 		/* rule 545 */
	 tx324 = tx75*tzz852; 		/* rule 546 */
	 tx325 = tx317 + tx318 + tx319 + tx320 + tx321 + tx322 + tx323 + tx324; 		/* rule 547 */
	 tx326 = tx325*tzz827; 		/* rule 548 */
	 tx327 = tx154*tzz870; 		/* rule 549 */
	 ohx1y3 = tx326 + tx327; 		/* rule 550 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 1, ohx1y3);
	 tx328 = tx89 + y4; 		/* rule 552 */
	 tx329 = tx328*tzz833; 		/* rule 553 */
	 tx330 = tx22*tzz922; 		/* rule 554 */
	 tx331 = tx165*tzz864; 		/* rule 555 */
	 tzz909 = tx30*tzz872; 		/* rule 556 */
	 tx332 = tx5*tzz909; 		/* rule 557 */
	 tzz910 = tx27*tzz871; 		/* rule 558 */
	 tx333 = tx5*tzz910; 		/* rule 559 */
	 tx334 = tx75*tzz847; 		/* rule 560 */
	 tx335 = tx75*tzz848; 		/* rule 561 */
	 tx336 = tx75*tzz849; 		/* rule 562 */
	 tx337 = tx329 + tx330 + tx331 + tx332 + tx333 + tx334 + tx335 + tx336; 		/* rule 563 */
	 tx338 = tx337*tzz827; 		/* rule 564 */
	 tx339 = tx170*tzz870; 		/* rule 565 */
	 ohx1z3 = tx338 + tx339; 		/* rule 566 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 2, ohx1z3);
	 tx340 = tx64*tzz864; 		/* rule 568 */
	 tx341 = tx5*tzz887; 		/* rule 569 */
	 tx342 = tx75*tzz855; 		/* rule 570 */
	 tx343 = tx340 + tx341 + tx342; 		/* rule 571 */
	 tx344 = tx343*tzz827; 		/* rule 572 */
	 tx345 = tx173*tzz870; 		/* rule 573 */
	 ohx1x4 = tx344 + tx345; 		/* rule 574 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 0, ohx1x4);
	 tx346 = tx69*tzz833; 		/* rule 576 */
	 tx347 = tx63*tzz864; 		/* rule 577 */
	 tx348 = tx5*tzz893; 		/* rule 578 */
	 tx349 = tx75*tzz861; 		/* rule 579 */
	 tx350 = tx346 + tx347 + tx348 + tx349; 		/* rule 580 */
	 tx351 = tx350*tzz827; 		/* rule 581 */
	 tx352 = tx176*tzz870; 		/* rule 582 */
	 ohx1y4 = tx351 + tx352; 		/* rule 583 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 1, ohx1y4);
	 tx353 = tx15*tzz833; 		/* rule 585 */
	 tx354 = tx62*tzz864; 		/* rule 586 */
	 tx355 = tx5*tzz890; 		/* rule 587 */
	 tx356 = tx75*tzz857; 		/* rule 588 */
	 tx357 = tx353 + tx354 + tx355 + tx356; 		/* rule 589 */
	 tx358 = tx357*tzz827; 		/* rule 590 */
	 tx359 = tx179*tzz870; 		/* rule 591 */
	 ohx1z4 = tx358 + tx359; 		/* rule 592 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 2, ohx1z4);
	 tzz927 = tx12*tx22; 		/* rule 594 */
	 tx360 = tzz878*tzz927; 		/* rule 595 */
	 tx361 = tx85*tzz866; 		/* rule 596 */
	 tx362 = tx93*tzz868; 		/* rule 597 */
	 tx363 = tx360 + tx361 + tx362; 		/* rule 598 */
	 tx364 = tx363*tzz827; 		/* rule 599 */
	 tx365 = tx88*tzz862; 		/* rule 600 */
	 ohy1z1 = tx364 + tx365; 		/* rule 601 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 tx366 = tx315 + z3; 		/* rule 603 */
	 tx367 = tx366*tzz833; 		/* rule 604 */
	 tx368 = tx102*tzz868; 		/* rule 605 */
	 tx369 = tx12*tzz889; 		/* rule 606 */
	 tx370 = tx85*tzz863; 		/* rule 607 */
	 tx371 = tx367 + tx368 + tx369 + tx370; 		/* rule 608 */
	 tx372 = tx371*tzz827; 		/* rule 609 */
	 tx373 = tx105*tzz869; 		/* rule 610 */
	 ohy1x2 = tx372 + tx373; 		/* rule 611 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx374 = tx110*tzz868; 		/* rule 613 */
	 tx375 = tx12*tzz896; 		/* rule 614 */
	 tx376 = tx85*tzz867; 		/* rule 615 */
	 tx377 = tx374 + tx375 + tx376; 		/* rule 616 */
	 tx378 = tx377*tzz827; 		/* rule 617 */
	 tx379 = tx113*tzz869; 		/* rule 618 */
	 ohy1y2 = tx378 + tx379; 		/* rule 619 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx380 = tx118*tzz868; 		/* rule 621 */
	 tx381 = tx12*tzz892; 		/* rule 622 */
	 tx382 = tx8*tzz833; 		/* rule 623 */
	 tx383 = tx85*tzz865; 		/* rule 624 */
	 tx384 = tx380 + tx381 + tx382 + tx383; 		/* rule 625 */
	 tx385 = tx384*tzz827; 		/* rule 626 */
	 tx386 = tx121*tzz869; 		/* rule 627 */
	 ohy1z2 = tx385 + tx386; 		/* rule 628 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 tx387 = tx65 + z4; 		/* rule 630 */
	 tx388 = tx387*tzz833; 		/* rule 631 */
	 tx389 = tx132*tzz868; 		/* rule 632 */
	 tx390 = tx12*tzz905; 		/* rule 633 */
	 tx391 = tx12*tzz903; 		/* rule 634 */
	 tx392 = tx85*tzz845; 		/* rule 635 */
	 tx393 = tx85*tzz846; 		/* rule 636 */
	 tx394 = tx85*tzz844; 		/* rule 637 */
	 tx395 = tx318 + tx388 + tx389 + tx390 + tx391 + tx392 + tx393 + tx394; 		/* rule 638 */
	 tx396 = tx395*tzz827; 		/* rule 639 */
	 tx397 = tx138*tzz869; 		/* rule 640 */
	 ohy1x3 = tx396 + tx397; 		/* rule 641 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 0, ohy1x3);
	 tx398 = tx31*tzz882; 		/* rule 643 */
	 tx399 = tx149*tzz868; 		/* rule 644 */
	 tx400 = tx12*tzz912; 		/* rule 645 */
	 tx401 = tx12*tzz913; 		/* rule 646 */
	 tx402 = tx85*tzz850; 		/* rule 647 */
	 tx403 = tx85*tzz851; 		/* rule 648 */
	 tx404 = tx85*tzz852; 		/* rule 649 */
	 tx405 = tx306 + tx398 + tx399 + tx400 + tx401 + tx402 + tx403 + tx404; 		/* rule 650 */
	 tx406 = tx405*tzz827; 		/* rule 651 */
	 tx407 = tx154*tzz869; 		/* rule 652 */
	 ohy1y3 = tx406 + tx407; 		/* rule 653 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 1, ohy1y3);
	 tx408 = -x4; 		/* rule 655 */
	 tx409 = tx408 + x2; 		/* rule 656 */
	 tx410 = tx409*tzz833; 		/* rule 657 */
	 tx411 = tzz882*tzz927; 		/* rule 658 */
	 tx412 = tx165*tzz868; 		/* rule 659 */
	 tx413 = tx12*tzz909; 		/* rule 660 */
	 tx414 = tx12*tzz910; 		/* rule 661 */
	 tx415 = tx85*tzz847; 		/* rule 662 */
	 tx416 = tx85*tzz848; 		/* rule 663 */
	 tx417 = tx85*tzz849; 		/* rule 664 */
	 tx418 = tx410 + tx411 + tx412 + tx413 + tx414 + tx415 + tx416 + tx417; 		/* rule 665 */
	 tx419 = tx418*tzz827; 		/* rule 666 */
	 tx420 = tx170*tzz869; 		/* rule 667 */
	 ohy1z3 = tx419 + tx420; 		/* rule 668 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 2, ohy1z3);
	 tx421 = tx27*tzz833; 		/* rule 670 */
	 tx422 = tx64*tzz868; 		/* rule 671 */
	 tx423 = tx12*tzz887; 		/* rule 672 */
	 tx424 = tx85*tzz855; 		/* rule 673 */
	 tx425 = tx421 + tx422 + tx423 + tx424; 		/* rule 674 */
	 tx426 = tx425*tzz827; 		/* rule 675 */
	 tx427 = tx173*tzz869; 		/* rule 676 */
	 ohy1x4 = tx426 + tx427; 		/* rule 677 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 0, ohy1x4);
	 tx428 = tx63*tzz868; 		/* rule 679 */
	 tx429 = tx12*tzz893; 		/* rule 680 */
	 tx430 = tx85*tzz861; 		/* rule 681 */
	 tx431 = tx428 + tx429 + tx430; 		/* rule 682 */
	 tx432 = tx431*tzz827; 		/* rule 683 */
	 tx433 = tx176*tzz869; 		/* rule 684 */
	 ohy1y4 = tx432 + tx433; 		/* rule 685 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 1, ohy1y4);
	 tx434 = tx62*tzz868; 		/* rule 687 */
	 tx435 = tx12*tzz890; 		/* rule 688 */
	 tx436 = tx82*tzz833; 		/* rule 689 */
	 tx437 = tx85*tzz857; 		/* rule 690 */
	 tx438 = tx434 + tx435 + tx436 + tx437; 		/* rule 691 */
	 tx439 = tx438*tzz827; 		/* rule 692 */
	 tx440 = tx179*tzz869; 		/* rule 693 */
	 ohy1z4 = tx439 + tx440; 		/* rule 694 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 2, ohy1z4);
	 tx441 = tx17*tzz833; 		/* rule 696 */
	 tx442 = tx102*tzz866; 		/* rule 697 */
	 tx443 = tx22*tzz889; 		/* rule 698 */
	 tx444 = tx93*tzz863; 		/* rule 699 */
	 tx445 = tx441 + tx442 + tx443 + tx444; 		/* rule 700 */
	 tx446 = tx445*tzz827; 		/* rule 701 */
	 tx447 = tx105*tzz862; 		/* rule 702 */
	 ohz1x2 = tx446 + tx447; 		/* rule 703 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 tx448 = tx408 + x3; 		/* rule 705 */
	 tx449 = tx448*tzz833; 		/* rule 706 */
	 tx450 = tx110*tzz866; 		/* rule 707 */
	 tx451 = tx22*tzz896; 		/* rule 708 */
	 tx452 = tx93*tzz867; 		/* rule 709 */
	 tx453 = tx449 + tx450 + tx451 + tx452; 		/* rule 710 */
	 tx454 = tx453*tzz827; 		/* rule 711 */
	 tx455 = tx113*tzz862; 		/* rule 712 */
	 ohz1y2 = tx454 + tx455; 		/* rule 713 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx456 = tx118*tzz866; 		/* rule 715 */
	 tx457 = tx22*tzz892; 		/* rule 716 */
	 tx458 = tx93*tzz865; 		/* rule 717 */
	 tx459 = tx456 + tx457 + tx458; 		/* rule 718 */
	 tx460 = tx459*tzz827; 		/* rule 719 */
	 tx461 = tx121*tzz862; 		/* rule 720 */
	 ohz1z2 = tx460 + tx461; 		/* rule 721 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 tx462 = tx295 + y2; 		/* rule 723 */
	 tx463 = tx462*tzz833; 		/* rule 724 */
	 tx464 = tx132*tzz866; 		/* rule 725 */
	 tx465 = tx22*tzz905; 		/* rule 726 */
	 tx466 = tx22*tzz903; 		/* rule 727 */
	 tx467 = tx93*tzz845; 		/* rule 728 */
	 tx468 = tx93*tzz846; 		/* rule 729 */
	 tx469 = tx93*tzz844; 		/* rule 730 */
	 tx470 = tx330 + tx463 + tx464 + tx465 + tx466 + tx467 + tx468 + tx469; 		/* rule 731 */
	 tx471 = tx470*tzz827; 		/* rule 732 */
	 tx472 = tx138*tzz862; 		/* rule 733 */
	 ohz1x3 = tx471 + tx472; 		/* rule 734 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 0, ohz1x3);
	 tx473 = tx81 + x4; 		/* rule 736 */
	 tx474 = tx473*tzz833; 		/* rule 737 */
	 tx475 = tx149*tzz866; 		/* rule 738 */
	 tx476 = tx22*tzz912; 		/* rule 739 */
	 tx477 = tx22*tzz913; 		/* rule 740 */
	 tx478 = tx93*tzz850; 		/* rule 741 */
	 tx479 = tx93*tzz851; 		/* rule 742 */
	 tx480 = tx93*tzz852; 		/* rule 743 */
	 tx481 = tx411 + tx474 + tx475 + tx476 + tx477 + tx478 + tx479 + tx480; 		/* rule 744 */
	 tx482 = tx481*tzz827; 		/* rule 745 */
	 tx483 = tx154*tzz862; 		/* rule 746 */
	 ohz1y3 = tx482 + tx483; 		/* rule 747 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 1, ohz1y3);
	 tx484 = tx34*tzz882; 		/* rule 749 */
	 tx485 = tx165*tzz866; 		/* rule 750 */
	 tx486 = tx22*tzz909; 		/* rule 751 */
	 tx487 = tx22*tzz910; 		/* rule 752 */
	 tx488 = tx93*tzz847; 		/* rule 753 */
	 tx489 = tx93*tzz848; 		/* rule 754 */
	 tx490 = tx93*tzz849; 		/* rule 755 */
	 tx491 = tx306 + tx484 + tx485 + tx486 + tx487 + tx488 + tx489 + tx490; 		/* rule 756 */
	 tx492 = tx491*tzz827; 		/* rule 757 */
	 tx493 = tx170*tzz862; 		/* rule 758 */
	 ohz1z3 = tx492 + tx493; 		/* rule 759 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 2, ohz1z3);
	 tx494 = tx64*tzz866; 		/* rule 761 */
	 tx495 = tx22*tzz887; 		/* rule 762 */
	 tx496 = tx90*tzz833; 		/* rule 763 */
	 tx497 = tx93*tzz855; 		/* rule 764 */
	 tx498 = tx494 + tx495 + tx496 + tx497; 		/* rule 765 */
	 tx499 = tx498*tzz827; 		/* rule 766 */
	 tx500 = tx173*tzz862; 		/* rule 767 */
	 ohz1x4 = tx499 + tx500; 		/* rule 768 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 0, ohz1x4);
	 tx501 = tx7*tzz833; 		/* rule 770 */
	 tx502 = tx63*tzz866; 		/* rule 771 */
	 tx503 = tx22*tzz893; 		/* rule 772 */
	 tx504 = tx93*tzz861; 		/* rule 773 */
	 tx505 = tx501 + tx502 + tx503 + tx504; 		/* rule 774 */
	 tx506 = tx505*tzz827; 		/* rule 775 */
	 tx507 = tx176*tzz862; 		/* rule 776 */
	 ohz1y4 = tx506 + tx507; 		/* rule 777 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 1, ohz1y4);
	 tx508 = tx62*tzz866; 		/* rule 779 */
	 tx509 = tx22*tzz890; 		/* rule 780 */
	 tx510 = tx93*tzz857; 		/* rule 781 */
	 tx511 = tx508 + tx509 + tx510; 		/* rule 782 */
	 tx512 = tx511*tzz827; 		/* rule 783 */
	 tx513 = tx179*tzz862; 		/* rule 784 */
	 ohz1z4 = tx512 + tx513; 		/* rule 785 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 2, ohz1z4);
	 tx514 = tx102*tzz867; 		/* rule 787 */
	 tx515 = tx110*tzz863; 		/* rule 788 */
	 tzz921 = tx7*tzz879; 		/* rule 789 */
	 tx516 = tx15*tzz921; 		/* rule 790 */
	 tx517 = tx514 + tx515 + tx516; 		/* rule 791 */
	 tx518 = tx105*tx113*tzz830; 		/* rule 792 */
	 tx519 = tx517*tzz827; 		/* rule 793 */
	 ohx2y2 = tx518 + tx519; 		/* rule 794 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 tx520 = tx102*tzz865; 		/* rule 796 */
	 tx521 = tx118*tzz863; 		/* rule 797 */
	 tx522 = tx27*tzz921; 		/* rule 798 */
	 tx523 = tx520 + tx521 + tx522; 		/* rule 799 */
	 tzz926 = tx121*tzz830; 		/* rule 800 */
	 tx524 = tx105*tzz926; 		/* rule 801 */
	 tx525 = tx523*tzz827; 		/* rule 802 */
	 ohx2z2 = tx524 + tx525; 		/* rule 803 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 tx526 = tx102*tzz845; 		/* rule 805 */
	 tx527 = tx132*tzz863; 		/* rule 806 */
	 tx528 = tx71*tzz831; 		/* rule 807 */
	 tzz883 = -3*tzz859; 		/* rule 808 */
	 tx529 = tx40*tzz883; 		/* rule 809 */
	 tx530 = tx102*tzz846; 		/* rule 810 */
	 tx531 = tx102*tzz844; 		/* rule 811 */
	 tzz873 = -tzz839; 		/* rule 812 */
	 tzz902 = tx8*tzz873; 		/* rule 813 */
	 tx532 = tx7*tzz902; 		/* rule 814 */
	 tx533 = tx307 + tx526 + tx527 + tx528 + tx529 + tx530 + tx531 + tx532; 		/* rule 815 */
	 tzz914 = tx138*tzz830; 		/* rule 816 */
	 tx534 = tx105*tzz914; 		/* rule 817 */
	 tx535 = tx533*tzz827; 		/* rule 818 */
	 ohx2x3 = tx534 + tx535; 		/* rule 819 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 0, ohx2x3);
	 tx536 = tx106 + z4; 		/* rule 821 */
	 tx537 = tx102*tzz850; 		/* rule 822 */
	 tx538 = tx102*tzz851; 		/* rule 823 */
	 tx539 = tx536*tzz833; 		/* rule 824 */
	 tx540 = tx149*tzz863; 		/* rule 825 */
	 tzz911 = tx17*tzz873; 		/* rule 826 */
	 tx541 = tx7*tzz911; 		/* rule 827 */
	 tzz920 = tx7*tzz883; 		/* rule 828 */
	 tx542 = tx15*tzz920; 		/* rule 829 */
	 tx543 = tx102*tzz852; 		/* rule 830 */
	 tx544 = tx390 + tx537 + tx538 + tx539 + tx540 + tx541 + tx542 + tx543; 		/* rule 831 */
	 tzz895 = tx154*tzz830; 		/* rule 832 */
	 tx545 = tx105*tzz895; 		/* rule 833 */
	 tx546 = tx544*tzz827; 		/* rule 834 */
	 ohx2y3 = tx545 + tx546; 		/* rule 835 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 1, ohx2y3);
	 tx547 = tx295 + y1; 		/* rule 837 */
	 tx548 = tx102*tzz847; 		/* rule 838 */
	 tx549 = tx102*tzz848; 		/* rule 839 */
	 tx550 = tx547*tzz833; 		/* rule 840 */
	 tx551 = tx165*tzz863; 		/* rule 841 */
	 tzz908 = tx30*tzz873; 		/* rule 842 */
	 tx552 = tx7*tzz908; 		/* rule 843 */
	 tx553 = tx27*tzz920; 		/* rule 844 */
	 tx554 = tx102*tzz849; 		/* rule 845 */
	 tx555 = tx465 + tx548 + tx549 + tx550 + tx551 + tx552 + tx553 + tx554; 		/* rule 846 */
	 tzz885 = tx170*tzz830; 		/* rule 847 */
	 tx556 = tx105*tzz885; 		/* rule 848 */
	 tx557 = tx555*tzz827; 		/* rule 849 */
	 ohx2z3 = tx556 + tx557; 		/* rule 850 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 2, ohx2z3);
	 tx558 = tx64*tzz863; 		/* rule 852 */
	 tx559 = tx102*tzz855; 		/* rule 853 */
	 tx560 = tx7*tzz888; 		/* rule 854 */
	 tx561 = tx558 + tx559 + tx560; 		/* rule 855 */
	 tzz877 = tx173*tzz830; 		/* rule 856 */
	 tx562 = tx105*tzz877; 		/* rule 857 */
	 tx563 = tx561*tzz827; 		/* rule 858 */
	 ohx2x4 = tx562 + tx563; 		/* rule 859 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 0, ohx2x4);
	 tx564 = tx102*tzz861; 		/* rule 861 */
	 tx565 = tx22*tzz833; 		/* rule 862 */
	 tx566 = tx63*tzz863; 		/* rule 863 */
	 tx567 = tx7*tzz894; 		/* rule 864 */
	 tx568 = tx564 + tx565 + tx566 + tx567; 		/* rule 865 */
	 tzz876 = tx176*tzz830; 		/* rule 866 */
	 tx569 = tx105*tzz876; 		/* rule 867 */
	 tx570 = tx568*tzz827; 		/* rule 868 */
	 ohx2y4 = tx569 + tx570; 		/* rule 869 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 1, ohx2y4);
	 tx571 = tx102*tzz857; 		/* rule 871 */
	 tx572 = tx62*tzz863; 		/* rule 872 */
	 tx573 = tx7*tzz891; 		/* rule 873 */
	 tx574 = tx98*tzz833; 		/* rule 874 */
	 tx575 = tx571 + tx572 + tx573 + tx574; 		/* rule 875 */
	 tzz875 = tx179*tzz830; 		/* rule 876 */
	 tx576 = tx105*tzz875; 		/* rule 877 */
	 tx577 = tx575*tzz827; 		/* rule 878 */
	 ohx2z4 = tx576 + tx577; 		/* rule 879 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 2, ohx2z4);
	 tx578 = tx118*tzz867; 		/* rule 881 */
	 tx579 = tx110*tzz865; 		/* rule 882 */
	 tzz925 = tx15*tx27; 		/* rule 883 */
	 tx580 = tzz879*tzz925; 		/* rule 884 */
	 tx581 = tx578 + tx579 + tx580; 		/* rule 885 */
	 tx582 = tx113*tzz926; 		/* rule 886 */
	 tx583 = tx581*tzz827; 		/* rule 887 */
	 ohy2z2 = tx582 + tx583; 		/* rule 888 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	 tx584 = tx315 + z1; 		/* rule 890 */
	 tx585 = tx132*tzz867; 		/* rule 891 */
	 tx586 = tx584*tzz833; 		/* rule 892 */
	 tx587 = tx110*tzz845; 		/* rule 893 */
	 tx588 = tx110*tzz846; 		/* rule 894 */
	 tx589 = tx110*tzz844; 		/* rule 895 */
	 tx590 = tx15*tzz902; 		/* rule 896 */
	 tx591 = tx321 + tx542 + tx585 + tx586 + tx587 + tx588 + tx589 + tx590; 		/* rule 897 */
	 tx592 = tx113*tzz914; 		/* rule 898 */
	 tx593 = tx591*tzz827; 		/* rule 899 */
	 ohy2x3 = tx592 + tx593; 		/* rule 900 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 0, ohy2x3);
	 tx594 = tx110*tzz850; 		/* rule 902 */
	 tx595 = tx110*tzz851; 		/* rule 903 */
	 tx596 = tx149*tzz867; 		/* rule 904 */
	 tx597 = tx15*tzz911; 		/* rule 905 */
	 tx598 = tx32*tzz883; 		/* rule 906 */
	 tx599 = tx110*tzz852; 		/* rule 907 */
	 tx600 = tx401 + tx528 + tx594 + tx595 + tx596 + tx597 + tx598 + tx599; 		/* rule 908 */
	 tx601 = tx113*tzz895; 		/* rule 909 */
	 tx602 = tx600*tzz827; 		/* rule 910 */
	 ohy2y3 = tx601 + tx602; 		/* rule 911 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 1, ohy2y3);
	 tx603 = tx114 + x4; 		/* rule 913 */
	 tx604 = tx110*tzz847; 		/* rule 914 */
	 tx605 = tx165*tzz867; 		/* rule 915 */
	 tx606 = tx110*tzz848; 		/* rule 916 */
	 tx607 = tx603*tzz833; 		/* rule 917 */
	 tx608 = tx15*tzz908; 		/* rule 918 */
	 tx609 = tzz883*tzz925; 		/* rule 919 */
	 tx610 = tx110*tzz849; 		/* rule 920 */
	 tx611 = tx477 + tx604 + tx605 + tx606 + tx607 + tx608 + tx609 + tx610; 		/* rule 921 */
	 tx612 = tx113*tzz885; 		/* rule 922 */
	 tx613 = tx611*tzz827; 		/* rule 923 */
	 ohy2z3 = tx612 + tx613; 		/* rule 924 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 2, ohy2z3);
	 tx614 = tx107*tzz833; 		/* rule 926 */
	 tx615 = tx64*tzz867; 		/* rule 927 */
	 tx616 = tx110*tzz855; 		/* rule 928 */
	 tx617 = tx15*tzz888; 		/* rule 929 */
	 tx618 = tx614 + tx615 + tx616 + tx617; 		/* rule 930 */
	 tx619 = tx113*tzz877; 		/* rule 931 */
	 tx620 = tx618*tzz827; 		/* rule 932 */
	 ohy2x4 = tx619 + tx620; 		/* rule 933 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 0, ohy2x4);
	 tx621 = tx110*tzz861; 		/* rule 935 */
	 tx622 = tx63*tzz867; 		/* rule 936 */
	 tx623 = tx15*tzz894; 		/* rule 937 */
	 tx624 = tx621 + tx622 + tx623; 		/* rule 938 */
	 tx625 = tx113*tzz876; 		/* rule 939 */
	 tx626 = tx624*tzz827; 		/* rule 940 */
	 ohy2y4 = tx625 + tx626; 		/* rule 941 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 1, ohy2y4);
	 tx627 = tx110*tzz857; 		/* rule 943 */
	 tx628 = tx5*tzz833; 		/* rule 944 */
	 tx629 = tx62*tzz867; 		/* rule 945 */
	 tx630 = tx15*tzz891; 		/* rule 946 */
	 tx631 = tx627 + tx628 + tx629 + tx630; 		/* rule 947 */
	 tx632 = tx113*tzz875; 		/* rule 948 */
	 tx633 = tx631*tzz827; 		/* rule 949 */
	 ohy2z4 = tx632 + tx633; 		/* rule 950 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 2, ohy2z4);
	 tx634 = tx97 + y4; 		/* rule 952 */
	 tx635 = tx132*tzz865; 		/* rule 953 */
	 tx636 = tx634*tzz833; 		/* rule 954 */
	 tx637 = tx118*tzz845; 		/* rule 955 */
	 tx638 = tx118*tzz846; 		/* rule 956 */
	 tx639 = tx118*tzz844; 		/* rule 957 */
	 tx640 = tx27*tzz902; 		/* rule 958 */
	 tx641 = tx333 + tx553 + tx635 + tx636 + tx637 + tx638 + tx639 + tx640; 		/* rule 959 */
	 tx642 = tx121*tzz914; 		/* rule 960 */
	 tx643 = tx641*tzz827; 		/* rule 961 */
	 ohz2x3 = tx642 + tx643; 		/* rule 962 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 0, ohz2x3);
	 tx644 = tx408 + x1; 		/* rule 964 */
	 tx645 = tx118*tzz850; 		/* rule 965 */
	 tx646 = tx118*tzz851; 		/* rule 966 */
	 tx647 = tx149*tzz865; 		/* rule 967 */
	 tx648 = tx644*tzz833; 		/* rule 968 */
	 tx649 = tx27*tzz911; 		/* rule 969 */
	 tx650 = tx118*tzz852; 		/* rule 970 */
	 tx651 = tx414 + tx609 + tx645 + tx646 + tx647 + tx648 + tx649 + tx650; 		/* rule 971 */
	 tx652 = tx121*tzz895; 		/* rule 972 */
	 tx653 = tx651*tzz827; 		/* rule 973 */
	 ohz2y3 = tx652 + tx653; 		/* rule 974 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 1, ohz2y3);
	 tx654 = tx118*tzz847; 		/* rule 976 */
	 tx655 = tx118*tzz848; 		/* rule 977 */
	 tx656 = tx165*tzz865; 		/* rule 978 */
	 tx657 = tx27*tzz908; 		/* rule 979 */
	 tx658 = tx35*tzz883; 		/* rule 980 */
	 tx659 = tx118*tzz849; 		/* rule 981 */
	 tx660 = tx487 + tx528 + tx654 + tx655 + tx656 + tx657 + tx658 + tx659; 		/* rule 982 */
	 tx661 = tx121*tzz885; 		/* rule 983 */
	 tx662 = tx660*tzz827; 		/* rule 984 */
	 ohz2z3 = tx661 + tx662; 		/* rule 985 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 2, ohz2z3);
	 tx663 = tx12*tzz833; 		/* rule 987 */
	 tx664 = tx64*tzz865; 		/* rule 988 */
	 tx665 = tx118*tzz855; 		/* rule 989 */
	 tx666 = tx27*tzz888; 		/* rule 990 */
	 tx667 = tx663 + tx664 + tx665 + tx666; 		/* rule 991 */
	 tx668 = tx121*tzz877; 		/* rule 992 */
	 tx669 = tx667*tzz827; 		/* rule 993 */
	 ohz2x4 = tx668 + tx669; 		/* rule 994 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 0, ohz2x4);
	 tx670 = tx118*tzz861; 		/* rule 996 */
	 tx671 = tx115*tzz833; 		/* rule 997 */
	 tx672 = tx63*tzz865; 		/* rule 998 */
	 tx673 = tx27*tzz894; 		/* rule 999 */
	 tx674 = tx670 + tx671 + tx672 + tx673; 		/* rule 1000 */
	 tx675 = tx121*tzz876; 		/* rule 1001 */
	 tx676 = tx674*tzz827; 		/* rule 1002 */
	 ohz2y4 = tx675 + tx676; 		/* rule 1003 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 1, ohz2y4);
	 tx677 = tx118*tzz857; 		/* rule 1005 */
	 tx678 = tx62*tzz865; 		/* rule 1006 */
	 tx679 = tx27*tzz891; 		/* rule 1007 */
	 tx680 = tx677 + tx678 + tx679; 		/* rule 1008 */
	 tx681 = tx121*tzz875; 		/* rule 1009 */
	 tx682 = tx680*tzz827; 		/* rule 1010 */
	 ohz2z4 = tx681 + tx682; 		/* rule 1011 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 2, ohz2z4);
	 tx683 = tx132*tzz850; 		/* rule 1013 */
	 tx684 = tx132*tzz851; 		/* rule 1014 */
	 tx685 = tx149*tzz845; 		/* rule 1015 */
	 tx686 = tx132*tzz852; 		/* rule 1016 */
	 tx687 = tx149*tzz846; 		/* rule 1017 */
	 tx688 = tx149*tzz844; 		/* rule 1018 */
	 tzz919 = tx8*tzz880; 		/* rule 1019 */
	 tx689 = tx17*tzz919; 		/* rule 1020 */
	 tx690 = tx270 + tx290 + tx348 + tx369 + tx423 + tx516 + tx567 + tx617 + tx683 + tx684 + tx685 + tx686 + tx687 + tx688 + tx689; 		/* rule 1021 */
	 tx691 = tx138*tzz895; 		/* rule 1022 */
	 tx692 = tx690*tzz827; 		/* rule 1023 */
	 ohx3y3 = tx691 + tx692; 		/* rule 1024 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 1, ohx3y3);
	 tx693 = tx132*tzz847; 		/* rule 1026 */
	 tx694 = tx132*tzz848; 		/* rule 1027 */
	 tx695 = tx165*tzz845; 		/* rule 1028 */
	 tx696 = tx132*tzz849; 		/* rule 1029 */
	 tx697 = tx165*tzz846; 		/* rule 1030 */
	 tx698 = tx165*tzz844; 		/* rule 1031 */
	 tx699 = tx30*tzz919; 		/* rule 1032 */
	 tx700 = tx276 + tx299 + tx355 + tx443 + tx495 + tx522 + tx573 + tx666 + tx693 + tx694 + tx695 + tx696 + tx697 + tx698 + tx699; 		/* rule 1033 */
	 tx701 = tx138*tzz885; 		/* rule 1034 */
	 tx702 = tx700*tzz827; 		/* rule 1035 */
	 ohx3z3 = tx701 + tx702; 		/* rule 1036 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 2, ohx3z3);
	 tx703 = tx64*tzz845; 		/* rule 1038 */
	 tx704 = tx71*tzz828; 		/* rule 1039 */
	 tzz884 = -3*tzz858; 		/* rule 1040 */
	 tx705 = tx41*tzz884; 		/* rule 1041 */
	 tx706 = tx64*tzz846; 		/* rule 1042 */
	 tx707 = tx132*tzz855; 		/* rule 1043 */
	 tx708 = tx64*tzz844; 		/* rule 1044 */
	 tx709 = tx310 + tx532 + tx703 + tx704 + tx705 + tx706 + tx707 + tx708; 		/* rule 1045 */
	 tx710 = tx138*tzz877; 		/* rule 1046 */
	 tx711 = tx709*tzz827; 		/* rule 1047 */
	 ohx3x4 = tx710 + tx711; 		/* rule 1048 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 0, ohx3x4);
	 tx712 = tx132*tzz861; 		/* rule 1050 */
	 tx713 = tx123*tzz833; 		/* rule 1051 */
	 tx714 = tx63*tzz845; 		/* rule 1052 */
	 tx715 = tx63*tzz846; 		/* rule 1053 */
	 tx716 = tx63*tzz844; 		/* rule 1054 */
	 tzz918 = tx8*tzz884; 		/* rule 1055 */
	 tx717 = tx17*tzz918; 		/* rule 1056 */
	 tx718 = tx320 + tx541 + tx712 + tx713 + tx714 + tx715 + tx716 + tx717; 		/* rule 1057 */
	 tx719 = tx138*tzz876; 		/* rule 1058 */
	 tx720 = tx718*tzz827; 		/* rule 1059 */
	 ohx3y4 = tx719 + tx720; 		/* rule 1060 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 1, ohx3y4);
	 tx721 = tx132*tzz857; 		/* rule 1062 */
	 tx722 = tx122*tzz833; 		/* rule 1063 */
	 tx723 = tx62*tzz845; 		/* rule 1064 */
	 tx724 = tx62*tzz846; 		/* rule 1065 */
	 tx725 = tx62*tzz844; 		/* rule 1066 */
	 tx726 = tx30*tzz918; 		/* rule 1067 */
	 tx727 = tx332 + tx552 + tx721 + tx722 + tx723 + tx724 + tx725 + tx726; 		/* rule 1068 */
	 tx728 = tx138*tzz875; 		/* rule 1069 */
	 tx729 = tx727*tzz827; 		/* rule 1070 */
	 ohx3z4 = tx728 + tx729; 		/* rule 1071 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 2, ohx3z4);
	 tx730 = tx165*tzz850; 		/* rule 1073 */
	 tx731 = tx149*tzz847; 		/* rule 1074 */
	 tx732 = tx165*tzz851; 		/* rule 1075 */
	 tx733 = tx149*tzz848; 		/* rule 1076 */
	 tzz924 = tx17*tx30; 		/* rule 1077 */
	 tx734 = tzz880*tzz924; 		/* rule 1078 */
	 tx735 = tx165*tzz852; 		/* rule 1079 */
	 tx736 = tx149*tzz849; 		/* rule 1080 */
	 tx737 = tx360 + tx381 + tx435 + tx451 + tx503 + tx580 + tx630 + tx673 + tx730 + tx731 + tx732 + tx733 + tx734 + tx735 + tx736; 		/* rule 1081 */
	 tx738 = tx154*tzz885; 		/* rule 1082 */
	 tx739 = tx737*tzz827; 		/* rule 1083 */
	 ohy3z3 = tx738 + tx739; 		/* rule 1084 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 2, ohy3z3);
	 tx740 = tx140*tzz833; 		/* rule 1086 */
	 tx741 = tx64*tzz850; 		/* rule 1087 */
	 tx742 = tx64*tzz851; 		/* rule 1088 */
	 tx743 = tx64*tzz852; 		/* rule 1089 */
	 tx744 = tx149*tzz855; 		/* rule 1090 */
	 tx745 = tx391 + tx590 + tx717 + tx740 + tx741 + tx742 + tx743 + tx744; 		/* rule 1091 */
	 tx746 = tx154*tzz877; 		/* rule 1092 */
	 tx747 = tx745*tzz827; 		/* rule 1093 */
	 ohy3x4 = tx746 + tx747; 		/* rule 1094 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 0, ohy3x4);
	 tx748 = tx149*tzz861; 		/* rule 1096 */
	 tx749 = tx63*tzz850; 		/* rule 1097 */
	 tx750 = tx63*tzz851; 		/* rule 1098 */
	 tx751 = tx33*tzz884; 		/* rule 1099 */
	 tx752 = tx63*tzz852; 		/* rule 1100 */
	 tx753 = tx400 + tx597 + tx704 + tx748 + tx749 + tx750 + tx751 + tx752; 		/* rule 1101 */
	 tx754 = tx154*tzz876; 		/* rule 1102 */
	 tx755 = tx753*tzz827; 		/* rule 1103 */
	 ohy3y4 = tx754 + tx755; 		/* rule 1104 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 1, ohy3y4);
	 tx756 = tx149*tzz857; 		/* rule 1106 */
	 tx757 = tx139*tzz833; 		/* rule 1107 */
	 tx758 = tx62*tzz850; 		/* rule 1108 */
	 tx759 = tx62*tzz851; 		/* rule 1109 */
	 tx760 = tzz884*tzz924; 		/* rule 1110 */
	 tx761 = tx62*tzz852; 		/* rule 1111 */
	 tx762 = tx413 + tx608 + tx756 + tx757 + tx758 + tx759 + tx760 + tx761; 		/* rule 1112 */
	 tx763 = tx154*tzz875; 		/* rule 1113 */
	 tx764 = tx762*tzz827; 		/* rule 1114 */
	 ohy3z4 = tx763 + tx764; 		/* rule 1115 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 2, ohy3z4);
	 tx765 = tx156*tzz833; 		/* rule 1117 */
	 tx766 = tx64*tzz847; 		/* rule 1118 */
	 tx767 = tx64*tzz848; 		/* rule 1119 */
	 tx768 = tx64*tzz849; 		/* rule 1120 */
	 tx769 = tx165*tzz855; 		/* rule 1121 */
	 tx770 = tx466 + tx640 + tx726 + tx765 + tx766 + tx767 + tx768 + tx769; 		/* rule 1122 */
	 tx771 = tx170*tzz877; 		/* rule 1123 */
	 tx772 = tx770*tzz827; 		/* rule 1124 */
	 ohz3x4 = tx771 + tx772; 		/* rule 1125 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 0, ohz3x4);
	 tx773 = tx165*tzz861; 		/* rule 1127 */
	 tx774 = tx155*tzz833; 		/* rule 1128 */
	 tx775 = tx63*tzz847; 		/* rule 1129 */
	 tx776 = tx63*tzz848; 		/* rule 1130 */
	 tx777 = tx63*tzz849; 		/* rule 1131 */
	 tx778 = tx476 + tx649 + tx760 + tx773 + tx774 + tx775 + tx776 + tx777; 		/* rule 1132 */
	 tx779 = tx170*tzz876; 		/* rule 1133 */
	 tx780 = tx778*tzz827; 		/* rule 1134 */
	 ohz3y4 = tx779 + tx780; 		/* rule 1135 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 1, ohz3y4);
	 tx781 = tx165*tzz857; 		/* rule 1137 */
	 tx782 = tx62*tzz847; 		/* rule 1138 */
	 tx783 = tx62*tzz848; 		/* rule 1139 */
	 tx784 = tx36*tzz884; 		/* rule 1140 */
	 tx785 = tx62*tzz849; 		/* rule 1141 */
	 tx786 = tx486 + tx657 + tx704 + tx781 + tx782 + tx783 + tx784 + tx785; 		/* rule 1142 */
	 tx787 = tx170*tzz875; 		/* rule 1143 */
	 tx788 = tx786*tzz827; 		/* rule 1144 */
	 ohz3z4 = tx787 + tx788; 		/* rule 1145 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 2, ohz3z4);
	 tx789 = tx64*tzz861; 		/* rule 1147 */
	 tx790 = tx63*tzz855; 		/* rule 1148 */
	 tx791 = tx689 + tx789 + tx790; 		/* rule 1149 */
	 tx792 = tx173*tzz876; 		/* rule 1150 */
	 tx793 = tx791*tzz827; 		/* rule 1151 */
	 ohx4y4 = tx792 + tx793; 		/* rule 1152 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 1, ohx4y4);
	 tx794 = tx64*tzz857; 		/* rule 1154 */
	 tx795 = tx62*tzz855; 		/* rule 1155 */
	 tx796 = tx699 + tx794 + tx795; 		/* rule 1156 */
	 tx797 = tx173*tzz875; 		/* rule 1157 */
	 tx798 = tx796*tzz827; 		/* rule 1158 */
	 ohx4z4 = tx797 + tx798; 		/* rule 1159 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 2, ohx4z4);
	 tx799 = tx62*tzz861; 		/* rule 1161 */
	 tx800 = tx63*tzz857; 		/* rule 1162 */
	 tx801 = tx734 + tx799 + tx800; 		/* rule 1163 */
	 tx802 = tx176*tzz875; 		/* rule 1164 */
	 tx803 = tx801*tzz827; 		/* rule 1165 */
	 ohy4z4 = tx802 + tx803; 		/* rule 1166 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 1, I4, 2, ohy4z4);
	} /*if calcOffDiagonalHessian */ 
	#endif /* CHIRAL_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* CHIRAL_RESTRAINT_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* CHIRAL_RESTRAINT_CALC_FORCE ]*/
	CHIRAL_RESTRAINT_DEBUG_INTERACTIONS(I1,I2,I3,I4);
	SKIP_term: (void)0;
	} /*ChiralTest>0.0*/
