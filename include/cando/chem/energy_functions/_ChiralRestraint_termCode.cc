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
	 tx2 = -y3; 		/* rule 21 */
	 tx3 = -z3; 		/* rule 22 */
	 tx4 = tx2*x1; 		/* rule 23 */
	 tx5 = tx1 + x1; 		/* rule 24 */
	 tx6 = tx3*x2; 		/* rule 25 */
	 tx7 = tx1 + x2; 		/* rule 26 */
	 tx8 = tx1 + x4; 		/* rule 27 */
	 tx9 = tx3*y1; 		/* rule 28 */
	 tx10 = -(x2*y1); 		/* rule 29 */
	 tx11 = x3*y1; 		/* rule 30 */
	 tx12 = tx2 + y1; 		/* rule 31 */
	 tx13 = tx1*y2; 		/* rule 32 */
	 tx14 = x1*y2; 		/* rule 33 */
	 tx15 = tx2 + y2; 		/* rule 34 */
	 tx16 = x2*y3; 		/* rule 35 */
	 tx17 = tx2 + y4; 		/* rule 36 */
	 tx18 = tx1*z1; 		/* rule 37 */
	 tx19 = x2*z1; 		/* rule 38 */
	 tx20 = -(y2*z1); 		/* rule 39 */
	 tx21 = y3*z1; 		/* rule 40 */
	 tx22 = tx3 + z1; 		/* rule 41 */
	 tx23 = tx2*z2; 		/* rule 42 */
	 tx24 = -(x1*z2); 		/* rule 43 */
	 tx25 = x3*z2; 		/* rule 44 */
	 tx26 = y1*z2; 		/* rule 45 */
	 tx27 = tx3 + z2; 		/* rule 46 */
	 tx28 = x1*z3; 		/* rule 47 */
	 tx29 = y2*z3; 		/* rule 48 */
	 tx30 = tx3 + z4; 		/* rule 49 */
	 tx31 = power2(tx12); 		/* rule 50 */
	 tx32 = power2(tx15); 		/* rule 51 */
	 tx33 = power2(tx17); 		/* rule 52 */
	 tx34 = power2(tx22); 		/* rule 53 */
	 tx35 = power2(tx27); 		/* rule 54 */
	 tx36 = power2(tx30); 		/* rule 55 */
	 tx37 = tx10 + tx11 + tx13 + tx14 + tx16 + tx4; 		/* rule 56 */
	 tx38 = power2(tx5); 		/* rule 57 */
	 tx39 = tx18 + tx19 + tx24 + tx25 + tx28 + tx6; 		/* rule 58 */
	 tx40 = power2(tx7); 		/* rule 59 */
	 tx41 = power2(tx8); 		/* rule 60 */
	 tx42 = tx20 + tx21 + tx23 + tx26 + tx29 + tx9; 		/* rule 61 */
	 tx43 = tx30*tx37; 		/* rule 62 */
	 tx44 = tx31 + tx34 + tx38; 		/* rule 63 */
	 tx45 = tx17*tx39; 		/* rule 64 */
	 tx46 = tx32 + tx35 + tx40; 		/* rule 65 */
	 tx47 = tx33 + tx36 + tx41; 		/* rule 66 */
	 tx48 = tx42*tx8; 		/* rule 67 */
	 tx789 = mysqrt(tx44); 		/* rule 68 */
	 tx49 = reciprocal(tx789); 		/* rule 69 */
	 tx790 = mysqrt(tx46); 		/* rule 70 */
	 tx50 = reciprocal(tx790); 		/* rule 71 */
	 tx791 = mysqrt(tx47); 		/* rule 72 */
	 tx51 = reciprocal(tx791); 		/* rule 73 */
	 tx52 = tx43 + tx45 + tx48; 		/* rule 74 */
	 tx53 = tx49*tx50*tx51*tx52; 		/* rule 75 */
	 tx54 = CO + tx53; 		/* rule 76 */
	 tx792 = power2(tx54); 		/* rule 77 */
	 tx55 = tx54*tx792; 		/* rule 78 */
	 ChiralTest = K*tx55; 		/* rule 79 */
/*If[ChiralTest>0.0,*/
	if (ChiralTest>0.0 ) {
	Energy=ChiralTest;
/*Energy = ChiralTest;*/
	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE(Energy);
	#ifdef CHIRAL_RESTRAINT_CALC_FORCE //[
	if ( calcForce ) {
	 tx56 = -z2; 		/* rule 87 */
	 tx57 = tx56 + z3; 		/* rule 88 */
	 tx58 = tx15*tx30; 		/* rule 89 */
	 tx59 = tx17*tx57; 		/* rule 90 */
	 tx793 = tx44; 		/* rule 91 */
	 tx794 = tx49; 		/* rule 92 */
	 tx795 = reciprocal(tx793); 		/* rule 93 */
	 tx60 = tx794*tx795; 		/* rule 94 */
	 tx61 = tx58 + tx59; 		/* rule 95 */
	 tx62 = -(tx5*tx50*tx51*tx52*tx60); 		/* rule 96 */
	 tx63 = tx50*tx51*tx61*tx794; 		/* rule 97 */
	 tx64 = tx792; 		/* rule 98 */
	 tx65 = tx62 + tx63; 		/* rule 99 */
	 gx1 = 3*K*tx64*tx65; 		/* rule 100 */
	 fx1 = -gx1; 		/* rule 101 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx66 = -x2; 		/* rule 103 */
	 tx67 = tx66 + x3; 		/* rule 104 */
	 tx68 = tx30*tx67; 		/* rule 105 */
	 tx69 = tx27*tx8; 		/* rule 106 */
	 tx70 = tx68 + tx69; 		/* rule 107 */
	 tx71 = -(tx12*tx50*tx51*tx52*tx60); 		/* rule 108 */
	 tx72 = tx50*tx51*tx70*tx794; 		/* rule 109 */
	 tx73 = tx71 + tx72; 		/* rule 110 */
	 gy1 = 3*K*tx64*tx73; 		/* rule 111 */
	 fy1 = -gy1; 		/* rule 112 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx74 = -y2; 		/* rule 114 */
	 tx75 = tx74 + y3; 		/* rule 115 */
	 tx76 = tx17*tx7; 		/* rule 116 */
	 tx77 = tx75*tx8; 		/* rule 117 */
	 tx78 = tx76 + tx77; 		/* rule 118 */
	 tx79 = -(tx22*tx50*tx51*tx52*tx60); 		/* rule 119 */
	 tx80 = tx50*tx51*tx78*tx794; 		/* rule 120 */
	 tx81 = tx79 + tx80; 		/* rule 121 */
	 gz1 = 3*K*tx64*tx81; 		/* rule 122 */
	 fz1 = -gz1; 		/* rule 123 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tx82 = -y1; 		/* rule 125 */
	 tx83 = tx82 + y3; 		/* rule 126 */
	 tx84 = tx17*tx22; 		/* rule 127 */
	 tx85 = tx30*tx83; 		/* rule 128 */
	 tx796 = tx46; 		/* rule 129 */
	 tx797 = tx50; 		/* rule 130 */
	 tx798 = reciprocal(tx796); 		/* rule 131 */
	 tx86 = tx797*tx798; 		/* rule 132 */
	 tx87 = tx84 + tx85; 		/* rule 133 */
	 tx88 = -(tx51*tx52*tx7*tx794*tx86); 		/* rule 134 */
	 tx89 = tx51*tx794*tx797*tx87; 		/* rule 135 */
	 tx90 = tx88 + tx89; 		/* rule 136 */
	 gx2 = 3*K*tx64*tx90; 		/* rule 137 */
	 fx2 = -gx2; 		/* rule 138 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I2, 0, fx2 );
	 tx91 = -z1; 		/* rule 140 */
	 tx92 = tx91 + z3; 		/* rule 141 */
	 tx93 = tx30*tx5; 		/* rule 142 */
	 tx94 = tx8*tx92; 		/* rule 143 */
	 tx95 = tx93 + tx94; 		/* rule 144 */
	 tx96 = -(tx15*tx51*tx52*tx794*tx86); 		/* rule 145 */
	 tx97 = tx51*tx794*tx797*tx95; 		/* rule 146 */
	 tx98 = tx96 + tx97; 		/* rule 147 */
	 gy2 = 3*K*tx64*tx98; 		/* rule 148 */
	 fy2 = -gy2; 		/* rule 149 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I2, 1, fy2 );
	 tx99 = -x1; 		/* rule 151 */
	 tx100 = tx99 + x3; 		/* rule 152 */
	 tx101 = tx100*tx17; 		/* rule 153 */
	 tx102 = tx12*tx8; 		/* rule 154 */
	 tx103 = tx101 + tx102; 		/* rule 155 */
	 tx104 = tx103*tx51*tx794*tx797; 		/* rule 156 */
	 tx105 = -(tx27*tx51*tx52*tx794*tx86); 		/* rule 157 */
	 tx106 = tx104 + tx105; 		/* rule 158 */
	 gz2 = 3*K*tx106*tx64; 		/* rule 159 */
	 fz2 = -gz2; 		/* rule 160 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I2, 2, fz2 );
	 tx107 = tx74 + y1; 		/* rule 162 */
	 tx108 = tx91 + z2; 		/* rule 163 */
	 tx109 = tx108*tx17; 		/* rule 164 */
	 tx110 = -tx26; 		/* rule 165 */
	 tx111 = tx107*tx30; 		/* rule 166 */
	 tx112 = tx3*y2; 		/* rule 167 */
	 tx113 = tx2*z1; 		/* rule 168 */
	 tx114 = y2*z1; 		/* rule 169 */
	 tx115 = y3*z2; 		/* rule 170 */
	 tx116 = y1*z3; 		/* rule 171 */
	 tx117 = tx109 + tx110 + tx111 + tx112 + tx113 + tx114 + tx115 + tx116; 		/* rule 172 */
	 tx799 = tx47; 		/* rule 173 */
	 tx800 = tx51; 		/* rule 174 */
	 tx801 = reciprocal(tx799); 		/* rule 175 */
	 tx118 = tx800*tx801; 		/* rule 176 */
	 tx119 = tx117*tx794*tx797*tx800; 		/* rule 177 */
	 tx120 = tx5*tx52*tx60*tx797*tx800; 		/* rule 178 */
	 tx121 = tx118*tx52*tx794*tx797*tx8; 		/* rule 179 */
	 tx122 = tx52*tx7*tx794*tx800*tx86; 		/* rule 180 */
	 tx123 = tx119 + tx120 + tx121 + tx122; 		/* rule 181 */
	 gx3 = 3*K*tx123*tx64; 		/* rule 182 */
	 fx3 = -gx3; 		/* rule 183 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I3, 0, fx3 );
	 tx124 = tx99 + x2; 		/* rule 185 */
	 tx125 = tx56 + z1; 		/* rule 186 */
	 tx126 = -tx19; 		/* rule 187 */
	 tx127 = tx124*tx30; 		/* rule 188 */
	 tx128 = tx125*tx8; 		/* rule 189 */
	 tx129 = tx3*x1; 		/* rule 190 */
	 tx130 = x3*z1; 		/* rule 191 */
	 tx131 = tx1*z2; 		/* rule 192 */
	 tx132 = x1*z2; 		/* rule 193 */
	 tx133 = x2*z3; 		/* rule 194 */
	 tx134 = tx126 + tx127 + tx128 + tx129 + tx130 + tx131 + tx132 + tx133; 		/* rule 195 */
	 tx135 = tx134*tx794*tx797*tx800; 		/* rule 196 */
	 tx136 = tx118*tx17*tx52*tx794*tx797; 		/* rule 197 */
	 tx137 = tx12*tx52*tx60*tx797*tx800; 		/* rule 198 */
	 tx138 = tx15*tx52*tx794*tx800*tx86; 		/* rule 199 */
	 tx139 = tx135 + tx136 + tx137 + tx138; 		/* rule 200 */
	 gy3 = 3*K*tx139*tx64; 		/* rule 201 */
	 fy3 = -gy3; 		/* rule 202 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I3, 1, fy3 );
	 tx140 = tx66 + x1; 		/* rule 204 */
	 tx141 = tx82 + y2; 		/* rule 205 */
	 tx142 = -tx14; 		/* rule 206 */
	 tx143 = tx140*tx17; 		/* rule 207 */
	 tx144 = tx141*tx8; 		/* rule 208 */
	 tx145 = tx2*x2; 		/* rule 209 */
	 tx146 = tx1*y1; 		/* rule 210 */
	 tx147 = x2*y1; 		/* rule 211 */
	 tx148 = x3*y2; 		/* rule 212 */
	 tx149 = x1*y3; 		/* rule 213 */
	 tx150 = tx142 + tx143 + tx144 + tx145 + tx146 + tx147 + tx148 + tx149; 		/* rule 214 */
	 tx151 = tx150*tx794*tx797*tx800; 		/* rule 215 */
	 tx152 = tx118*tx30*tx52*tx794*tx797; 		/* rule 216 */
	 tx153 = tx22*tx52*tx60*tx797*tx800; 		/* rule 217 */
	 tx154 = tx27*tx52*tx794*tx800*tx86; 		/* rule 218 */
	 tx155 = tx151 + tx152 + tx153 + tx154; 		/* rule 219 */
	 gz3 = 3*K*tx155*tx64; 		/* rule 220 */
	 fz3 = -gz3; 		/* rule 221 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I3, 2, fz3 );
	 tx156 = -tx121; 		/* rule 223 */
	 tx157 = tx42*tx794*tx797*tx800; 		/* rule 224 */
	 tx158 = tx156 + tx157; 		/* rule 225 */
	 gx4 = 3*K*tx158*tx64; 		/* rule 226 */
	 fx4 = -gx4; 		/* rule 227 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I4, 0, fx4 );
	 tx159 = -tx136; 		/* rule 229 */
	 tx160 = tx39*tx794*tx797*tx800; 		/* rule 230 */
	 tx161 = tx159 + tx160; 		/* rule 231 */
	 gy4 = 3*K*tx161*tx64; 		/* rule 232 */
	 fy4 = -gy4; 		/* rule 233 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I4, 1, fy4 );
	 tx162 = -tx152; 		/* rule 235 */
	 tx163 = tx37*tx794*tx797*tx800; 		/* rule 236 */
	 tx164 = tx162 + tx163; 		/* rule 237 */
	 gz4 = 3*K*tx164*tx64; 		/* rule 238 */
	 fz4 = -gz4; 		/* rule 239 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I4, 2, fz4 );
	#ifdef CHIRAL_RESTRAINT_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx802 = power2(tx793); 		/* rule 243 */
	 tx803 = reciprocal(tx802); 		/* rule 244 */
	 tx165 = tx794*tx803; 		/* rule 245 */
	 tx166 = 3*tx165*tx38*tx52*tx797*tx800; 		/* rule 246 */
	 tx167 = -(tx52*tx60*tx797*tx800); 		/* rule 247 */
	 tx168 = -2.*tx5*tx60*tx61*tx797*tx800; 		/* rule 248 */
	 tx169 = tx166 + tx167 + tx168; 		/* rule 249 */
	 tx170 = power2(tx65); 		/* rule 250 */
	 tx171 = 6.*K*tx170*tx54; 		/* rule 251 */
	 tx172 = 3*K*tx169*tx64; 		/* rule 252 */
	 dhx1x1 = tx171 + tx172; 		/* rule 253 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx173 = 3*tx165*tx31*tx52*tx797*tx800; 		/* rule 255 */
	 tx174 = -2.*tx12*tx60*tx70*tx797*tx800; 		/* rule 256 */
	 tx175 = tx167 + tx173 + tx174; 		/* rule 257 */
	 tx176 = power2(tx73); 		/* rule 258 */
	 tx177 = 6.*K*tx176*tx54; 		/* rule 259 */
	 tx178 = 3*K*tx175*tx64; 		/* rule 260 */
	 dhy1y1 = tx177 + tx178; 		/* rule 261 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx179 = 3*tx165*tx34*tx52*tx797*tx800; 		/* rule 263 */
	 tx180 = -2.*tx22*tx60*tx78*tx797*tx800; 		/* rule 264 */
	 tx181 = tx167 + tx179 + tx180; 		/* rule 265 */
	 tx182 = power2(tx81); 		/* rule 266 */
	 tx183 = 6.*K*tx182*tx54; 		/* rule 267 */
	 tx184 = 3*K*tx181*tx64; 		/* rule 268 */
	 dhz1z1 = tx183 + tx184; 		/* rule 269 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 tx804 = power2(tx796); 		/* rule 271 */
	 tx805 = reciprocal(tx804); 		/* rule 272 */
	 tx185 = tx797*tx805; 		/* rule 273 */
	 tx186 = 3*tx185*tx40*tx52*tx794*tx800; 		/* rule 274 */
	 tx187 = -(tx52*tx794*tx800*tx86); 		/* rule 275 */
	 tx188 = -2.*tx7*tx794*tx800*tx86*tx87; 		/* rule 276 */
	 tx189 = tx186 + tx187 + tx188; 		/* rule 277 */
	 tx190 = power2(tx90); 		/* rule 278 */
	 tx191 = 6.*K*tx190*tx54; 		/* rule 279 */
	 tx192 = 3*K*tx189*tx64; 		/* rule 280 */
	 dhx2x2 = tx191 + tx192; 		/* rule 281 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 tx193 = 3*tx185*tx32*tx52*tx794*tx800; 		/* rule 283 */
	 tx194 = -2.*tx15*tx794*tx800*tx86*tx95; 		/* rule 284 */
	 tx195 = tx187 + tx193 + tx194; 		/* rule 285 */
	 tx196 = power2(tx98); 		/* rule 286 */
	 tx197 = 6.*K*tx196*tx54; 		/* rule 287 */
	 tx198 = 3*K*tx195*tx64; 		/* rule 288 */
	 dhy2y2 = tx197 + tx198; 		/* rule 289 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 tx199 = 3*tx185*tx35*tx52*tx794*tx800; 		/* rule 291 */
	 tx200 = -2.*tx103*tx27*tx794*tx800*tx86; 		/* rule 292 */
	 tx201 = power2(tx106); 		/* rule 293 */
	 tx202 = tx187 + tx199 + tx200; 		/* rule 294 */
	 tx203 = 6.*K*tx201*tx54; 		/* rule 295 */
	 tx204 = 3*K*tx202*tx64; 		/* rule 296 */
	 dhz2z2 = tx203 + tx204; 		/* rule 297 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	 tx806 = power2(tx799); 		/* rule 299 */
	 tx807 = reciprocal(tx806); 		/* rule 300 */
	 tx205 = tx800*tx807; 		/* rule 301 */
	 tx206 = -(tx118*tx52*tx794*tx797); 		/* rule 302 */
	 tx207 = 3*tx205*tx41*tx52*tx794*tx797; 		/* rule 303 */
	 tx208 = 2.*tx117*tx5*tx60*tx797*tx800; 		/* rule 304 */
	 tx209 = 2.*tx117*tx118*tx794*tx797*tx8; 		/* rule 305 */
	 tx210 = 2.*tx118*tx5*tx52*tx60*tx797*tx8; 		/* rule 306 */
	 tx211 = 2.*tx117*tx7*tx794*tx800*tx86; 		/* rule 307 */
	 tx212 = 2.*tx5*tx52*tx60*tx7*tx800*tx86; 		/* rule 308 */
	 tx213 = 2.*tx118*tx52*tx7*tx794*tx8*tx86; 		/* rule 309 */
	 tx214 = power2(tx123); 		/* rule 310 */
	 tx215 = tx166 + tx167 + tx186 + tx187 + tx206 + tx207 + tx208 + tx209 + tx210 + tx211 + tx212 + tx213; 		/* rule 311 */
	 tx216 = 6.*K*tx214*tx54; 		/* rule 312 */
	 tx217 = 3*K*tx215*tx64; 		/* rule 313 */
	 dhx3x3 = tx216 + tx217; 		/* rule 314 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 0, dhx3x3);
	 tx218 = 2.*tx118*tx134*tx17*tx794*tx797; 		/* rule 316 */
	 tx219 = 3*tx205*tx33*tx52*tx794*tx797; 		/* rule 317 */
	 tx220 = 2.*tx12*tx134*tx60*tx797*tx800; 		/* rule 318 */
	 tx221 = 2.*tx118*tx12*tx17*tx52*tx60*tx797; 		/* rule 319 */
	 tx222 = 2.*tx134*tx15*tx794*tx800*tx86; 		/* rule 320 */
	 tx223 = 2.*tx118*tx15*tx17*tx52*tx794*tx86; 		/* rule 321 */
	 tx224 = 2.*tx12*tx15*tx52*tx60*tx800*tx86; 		/* rule 322 */
	 tx225 = power2(tx139); 		/* rule 323 */
	 tx226 = tx167 + tx173 + tx187 + tx193 + tx206 + tx218 + tx219 + tx220 + tx221 + tx222 + tx223 + tx224; 		/* rule 324 */
	 tx227 = 6.*K*tx225*tx54; 		/* rule 325 */
	 tx228 = 3*K*tx226*tx64; 		/* rule 326 */
	 dhy3y3 = tx227 + tx228; 		/* rule 327 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 1, dhy3y3);
	 tx229 = 2.*tx118*tx150*tx30*tx794*tx797; 		/* rule 329 */
	 tx230 = 3*tx205*tx36*tx52*tx794*tx797; 		/* rule 330 */
	 tx231 = 2.*tx150*tx22*tx60*tx797*tx800; 		/* rule 331 */
	 tx232 = 2.*tx118*tx22*tx30*tx52*tx60*tx797; 		/* rule 332 */
	 tx233 = 2.*tx150*tx27*tx794*tx800*tx86; 		/* rule 333 */
	 tx234 = 2.*tx118*tx27*tx30*tx52*tx794*tx86; 		/* rule 334 */
	 tx235 = 2.*tx22*tx27*tx52*tx60*tx800*tx86; 		/* rule 335 */
	 tx236 = power2(tx155); 		/* rule 336 */
	 tx237 = tx167 + tx179 + tx187 + tx199 + tx206 + tx229 + tx230 + tx231 + tx232 + tx233 + tx234 + tx235; 		/* rule 337 */
	 tx238 = 6.*K*tx236*tx54; 		/* rule 338 */
	 tx239 = 3*K*tx237*tx64; 		/* rule 339 */
	 dhz3z3 = tx238 + tx239; 		/* rule 340 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I3, 2, dhz3z3);
	 tx240 = -2.*tx118*tx48*tx794*tx797; 		/* rule 342 */
	 tx241 = power2(tx158); 		/* rule 343 */
	 tx242 = tx206 + tx207 + tx240; 		/* rule 344 */
	 tx243 = 6.*K*tx241*tx54; 		/* rule 345 */
	 tx244 = 3*K*tx242*tx64; 		/* rule 346 */
	 dhx4x4 = tx243 + tx244; 		/* rule 347 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 0, dhx4x4);
	 tx245 = -2.*tx118*tx45*tx794*tx797; 		/* rule 349 */
	 tx246 = power2(tx161); 		/* rule 350 */
	 tx247 = tx206 + tx219 + tx245; 		/* rule 351 */
	 tx248 = 6.*K*tx246*tx54; 		/* rule 352 */
	 tx249 = 3*K*tx247*tx64; 		/* rule 353 */
	 dhy4y4 = tx248 + tx249; 		/* rule 354 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I4, 1, I4, 1, dhy4y4);
	 tx250 = -2.*tx118*tx43*tx794*tx797; 		/* rule 356 */
	 tx251 = power2(tx164); 		/* rule 357 */
	 tx252 = tx206 + tx230 + tx250; 		/* rule 358 */
	 tx253 = 6.*K*tx251*tx54; 		/* rule 359 */
	 tx254 = 3*K*tx252*tx64; 		/* rule 360 */
	 dhz4z4 = tx253 + tx254; 		/* rule 361 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I4, 2, I4, 2, dhz4z4);
	#ifdef CHIRAL_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tx255 = 3*tx12*tx165*tx5*tx52*tx797*tx800; 		/* rule 365 */
	 tx256 = -(tx12*tx60*tx61*tx797*tx800); 		/* rule 366 */
	 tx257 = -(tx5*tx60*tx70*tx797*tx800); 		/* rule 367 */
	 tx258 = tx255 + tx256 + tx257; 		/* rule 368 */
	 tx259 = 3*K*tx258*tx64; 		/* rule 369 */
	 tx260 = 6.*K*tx54*tx65*tx73; 		/* rule 370 */
	 ohx1y1 = tx259 + tx260; 		/* rule 371 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tx261 = 3*tx165*tx22*tx5*tx52*tx797*tx800; 		/* rule 373 */
	 tx262 = -(tx22*tx60*tx61*tx797*tx800); 		/* rule 374 */
	 tx263 = -(tx5*tx60*tx78*tx797*tx800); 		/* rule 375 */
	 tx264 = tx261 + tx262 + tx263; 		/* rule 376 */
	 tx265 = 3*K*tx264*tx64; 		/* rule 377 */
	 tx266 = 6.*K*tx54*tx65*tx81; 		/* rule 378 */
	 ohx1z1 = tx265 + tx266; 		/* rule 379 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx267 = tx5*tx52*tx60*tx7*tx800*tx86; 		/* rule 381 */
	 tx268 = -(tx61*tx7*tx794*tx800*tx86); 		/* rule 382 */
	 tx269 = -(tx5*tx60*tx797*tx800*tx87); 		/* rule 383 */
	 tx270 = tx267 + tx268 + tx269; 		/* rule 384 */
	 tx271 = 3*K*tx270*tx64; 		/* rule 385 */
	 tx272 = 6.*K*tx54*tx65*tx90; 		/* rule 386 */
	 ohx1x2 = tx271 + tx272; 		/* rule 387 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx273 = tx30*tx794*tx797*tx800; 		/* rule 389 */
	 tx274 = tx15*tx5*tx52*tx60*tx800*tx86; 		/* rule 390 */
	 tx275 = -(tx15*tx61*tx794*tx800*tx86); 		/* rule 391 */
	 tx276 = -(tx5*tx60*tx797*tx800*tx95); 		/* rule 392 */
	 tx277 = tx273 + tx274 + tx275 + tx276; 		/* rule 393 */
	 tx278 = 3*K*tx277*tx64; 		/* rule 394 */
	 tx279 = 6.*K*tx54*tx65*tx98; 		/* rule 395 */
	 ohx1y2 = tx278 + tx279; 		/* rule 396 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx280 = -y4; 		/* rule 398 */
	 tx281 = tx280 + y3; 		/* rule 399 */
	 tx282 = tx281*tx794*tx797*tx800; 		/* rule 400 */
	 tx283 = -(tx103*tx5*tx60*tx797*tx800); 		/* rule 401 */
	 tx284 = tx27*tx5*tx52*tx60*tx800*tx86; 		/* rule 402 */
	 tx285 = -(tx27*tx61*tx794*tx800*tx86); 		/* rule 403 */
	 tx286 = tx282 + tx283 + tx284 + tx285; 		/* rule 404 */
	 tx287 = 3*K*tx286*tx64; 		/* rule 405 */
	 tx288 = 6.*K*tx106*tx54*tx65; 		/* rule 406 */
	 ohx1z2 = tx287 + tx288; 		/* rule 407 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tx289 = -tx267; 		/* rule 409 */
	 tx290 = -3*tx165*tx38*tx52*tx797*tx800; 		/* rule 410 */
	 tx291 = -(tx117*tx5*tx60*tx797*tx800); 		/* rule 411 */
	 tx292 = tx52*tx60*tx797*tx800; 		/* rule 412 */
	 tx293 = tx5*tx60*tx61*tx797*tx800; 		/* rule 413 */
	 tx294 = -(tx118*tx5*tx52*tx60*tx797*tx8); 		/* rule 414 */
	 tx295 = tx118*tx61*tx794*tx797*tx8; 		/* rule 415 */
	 tx296 = tx61*tx7*tx794*tx800*tx86; 		/* rule 416 */
	 tx297 = tx289 + tx290 + tx291 + tx292 + tx293 + tx294 + tx295 + tx296; 		/* rule 417 */
	 tx298 = 3*K*tx297*tx64; 		/* rule 418 */
	 tx299 = 6.*K*tx123*tx54*tx65; 		/* rule 419 */
	 ohx1x3 = tx298 + tx299; 		/* rule 420 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 0, ohx1x3);
	 tx300 = -z4; 		/* rule 422 */
	 tx301 = tx300 + z2; 		/* rule 423 */
	 tx302 = -tx274; 		/* rule 424 */
	 tx303 = tx301*tx794*tx797*tx800; 		/* rule 425 */
	 tx304 = -3*tx12*tx165*tx5*tx52*tx797*tx800; 		/* rule 426 */
	 tx305 = -(tx134*tx5*tx60*tx797*tx800); 		/* rule 427 */
	 tx306 = -(tx118*tx17*tx5*tx52*tx60*tx797); 		/* rule 428 */
	 tx307 = tx118*tx17*tx61*tx794*tx797; 		/* rule 429 */
	 tx308 = tx12*tx60*tx61*tx797*tx800; 		/* rule 430 */
	 tx309 = tx15*tx61*tx794*tx800*tx86; 		/* rule 431 */
	 tx310 = tx302 + tx303 + tx304 + tx305 + tx306 + tx307 + tx308 + tx309; 		/* rule 432 */
	 tx311 = 3*K*tx310*tx64; 		/* rule 433 */
	 tx312 = 6.*K*tx139*tx54*tx65; 		/* rule 434 */
	 ohx1y3 = tx311 + tx312; 		/* rule 435 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 1, ohx1y3);
	 tx313 = tx74 + y4; 		/* rule 437 */
	 tx314 = -tx284; 		/* rule 438 */
	 tx315 = tx313*tx794*tx797*tx800; 		/* rule 439 */
	 tx316 = -3*tx165*tx22*tx5*tx52*tx797*tx800; 		/* rule 440 */
	 tx317 = -(tx150*tx5*tx60*tx797*tx800); 		/* rule 441 */
	 tx318 = tx118*tx30*tx61*tx794*tx797; 		/* rule 442 */
	 tx319 = tx22*tx60*tx61*tx797*tx800; 		/* rule 443 */
	 tx320 = tx27*tx61*tx794*tx800*tx86; 		/* rule 444 */
	 tx321 = -(tx118*tx52*tx60*tx797*tx93); 		/* rule 445 */
	 tx322 = tx314 + tx315 + tx316 + tx317 + tx318 + tx319 + tx320 + tx321; 		/* rule 446 */
	 tx323 = 3*K*tx322*tx64; 		/* rule 447 */
	 tx324 = 6.*K*tx155*tx54*tx65; 		/* rule 448 */
	 ohx1z3 = tx323 + tx324; 		/* rule 449 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 2, ohx1z3);
	 tx325 = -tx295; 		/* rule 451 */
	 tx326 = -(tx42*tx5*tx60*tx797*tx800); 		/* rule 452 */
	 tx327 = tx118*tx5*tx52*tx60*tx797*tx8; 		/* rule 453 */
	 tx328 = tx325 + tx326 + tx327; 		/* rule 454 */
	 tx329 = 3*K*tx328*tx64; 		/* rule 455 */
	 tx330 = 6.*K*tx158*tx54*tx65; 		/* rule 456 */
	 ohx1x4 = tx329 + tx330; 		/* rule 457 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 0, ohx1x4);
	 tx331 = -tx307; 		/* rule 459 */
	 tx332 = tx57*tx794*tx797*tx800; 		/* rule 460 */
	 tx333 = -(tx39*tx5*tx60*tx797*tx800); 		/* rule 461 */
	 tx334 = tx118*tx17*tx5*tx52*tx60*tx797; 		/* rule 462 */
	 tx335 = tx331 + tx332 + tx333 + tx334; 		/* rule 463 */
	 tx336 = 3*K*tx335*tx64; 		/* rule 464 */
	 tx337 = 6.*K*tx161*tx54*tx65; 		/* rule 465 */
	 ohx1y4 = tx336 + tx337; 		/* rule 466 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 1, ohx1y4);
	 tx338 = -tx318; 		/* rule 468 */
	 tx339 = tx15*tx794*tx797*tx800; 		/* rule 469 */
	 tx340 = -(tx37*tx5*tx60*tx797*tx800); 		/* rule 470 */
	 tx341 = tx118*tx52*tx60*tx797*tx93; 		/* rule 471 */
	 tx342 = tx338 + tx339 + tx340 + tx341; 		/* rule 472 */
	 tx343 = 3*K*tx342*tx64; 		/* rule 473 */
	 tx344 = 6.*K*tx164*tx54*tx65; 		/* rule 474 */
	 ohx1z4 = tx343 + tx344; 		/* rule 475 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 2, ohx1z4);
	 tx345 = 3*tx12*tx165*tx22*tx52*tx797*tx800; 		/* rule 477 */
	 tx346 = -(tx22*tx60*tx70*tx797*tx800); 		/* rule 478 */
	 tx347 = -(tx12*tx60*tx78*tx797*tx800); 		/* rule 479 */
	 tx348 = tx345 + tx346 + tx347; 		/* rule 480 */
	 tx349 = 3*K*tx348*tx64; 		/* rule 481 */
	 tx350 = 6.*K*tx54*tx73*tx81; 		/* rule 482 */
	 ohy1z1 = tx349 + tx350; 		/* rule 483 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 tx351 = tx300 + z3; 		/* rule 485 */
	 tx352 = tx351*tx794*tx797*tx800; 		/* rule 486 */
	 tx353 = tx12*tx52*tx60*tx7*tx800*tx86; 		/* rule 487 */
	 tx354 = -(tx7*tx70*tx794*tx800*tx86); 		/* rule 488 */
	 tx355 = -(tx12*tx60*tx797*tx800*tx87); 		/* rule 489 */
	 tx356 = tx352 + tx353 + tx354 + tx355; 		/* rule 490 */
	 tx357 = 3*K*tx356*tx64; 		/* rule 491 */
	 tx358 = 6.*K*tx54*tx73*tx90; 		/* rule 492 */
	 ohy1x2 = tx357 + tx358; 		/* rule 493 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx359 = tx12*tx15*tx52*tx60*tx800*tx86; 		/* rule 495 */
	 tx360 = -(tx15*tx70*tx794*tx800*tx86); 		/* rule 496 */
	 tx361 = -(tx12*tx60*tx797*tx800*tx95); 		/* rule 497 */
	 tx362 = tx359 + tx360 + tx361; 		/* rule 498 */
	 tx363 = 3*K*tx362*tx64; 		/* rule 499 */
	 tx364 = 6.*K*tx54*tx73*tx98; 		/* rule 500 */
	 ohy1y2 = tx363 + tx364; 		/* rule 501 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx365 = -(tx103*tx12*tx60*tx797*tx800); 		/* rule 503 */
	 tx366 = tx794*tx797*tx8*tx800; 		/* rule 504 */
	 tx367 = tx12*tx27*tx52*tx60*tx800*tx86; 		/* rule 505 */
	 tx368 = -(tx27*tx70*tx794*tx800*tx86); 		/* rule 506 */
	 tx369 = tx365 + tx366 + tx367 + tx368; 		/* rule 507 */
	 tx370 = 3*K*tx369*tx64; 		/* rule 508 */
	 tx371 = 6.*K*tx106*tx54*tx73; 		/* rule 509 */
	 ohy1z2 = tx370 + tx371; 		/* rule 510 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 tx372 = tx56 + z4; 		/* rule 512 */
	 tx373 = -tx353; 		/* rule 513 */
	 tx374 = tx372*tx794*tx797*tx800; 		/* rule 514 */
	 tx375 = -(tx117*tx12*tx60*tx797*tx800); 		/* rule 515 */
	 tx376 = -(tx102*tx118*tx52*tx60*tx797); 		/* rule 516 */
	 tx377 = tx5*tx60*tx70*tx797*tx800; 		/* rule 517 */
	 tx378 = tx118*tx70*tx794*tx797*tx8; 		/* rule 518 */
	 tx379 = tx7*tx70*tx794*tx800*tx86; 		/* rule 519 */
	 tx380 = tx304 + tx373 + tx374 + tx375 + tx376 + tx377 + tx378 + tx379; 		/* rule 520 */
	 tx381 = 3*K*tx380*tx64; 		/* rule 521 */
	 tx382 = 6.*K*tx123*tx54*tx73; 		/* rule 522 */
	 ohy1x3 = tx381 + tx382; 		/* rule 523 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 0, ohy1x3);
	 tx383 = -tx359; 		/* rule 525 */
	 tx384 = -3*tx165*tx31*tx52*tx797*tx800; 		/* rule 526 */
	 tx385 = -(tx12*tx134*tx60*tx797*tx800); 		/* rule 527 */
	 tx386 = -(tx118*tx12*tx17*tx52*tx60*tx797); 		/* rule 528 */
	 tx387 = tx118*tx17*tx70*tx794*tx797; 		/* rule 529 */
	 tx388 = tx12*tx60*tx70*tx797*tx800; 		/* rule 530 */
	 tx389 = tx15*tx70*tx794*tx800*tx86; 		/* rule 531 */
	 tx390 = tx292 + tx383 + tx384 + tx385 + tx386 + tx387 + tx388 + tx389; 		/* rule 532 */
	 tx391 = 3*K*tx390*tx64; 		/* rule 533 */
	 tx392 = 6.*K*tx139*tx54*tx73; 		/* rule 534 */
	 ohy1y3 = tx391 + tx392; 		/* rule 535 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 1, ohy1y3);
	 tx393 = -x4; 		/* rule 537 */
	 tx394 = tx393 + x2; 		/* rule 538 */
	 tx395 = -tx367; 		/* rule 539 */
	 tx396 = tx394*tx794*tx797*tx800; 		/* rule 540 */
	 tx397 = -3*tx12*tx165*tx22*tx52*tx797*tx800; 		/* rule 541 */
	 tx398 = -(tx12*tx150*tx60*tx797*tx800); 		/* rule 542 */
	 tx399 = -(tx118*tx12*tx30*tx52*tx60*tx797); 		/* rule 543 */
	 tx400 = tx118*tx30*tx70*tx794*tx797; 		/* rule 544 */
	 tx401 = tx22*tx60*tx70*tx797*tx800; 		/* rule 545 */
	 tx402 = tx27*tx70*tx794*tx800*tx86; 		/* rule 546 */
	 tx403 = tx395 + tx396 + tx397 + tx398 + tx399 + tx400 + tx401 + tx402; 		/* rule 547 */
	 tx404 = 3*K*tx403*tx64; 		/* rule 548 */
	 tx405 = 6.*K*tx155*tx54*tx73; 		/* rule 549 */
	 ohy1z3 = tx404 + tx405; 		/* rule 550 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 2, ohy1z3);
	 tx406 = -tx378; 		/* rule 552 */
	 tx407 = tx27*tx794*tx797*tx800; 		/* rule 553 */
	 tx408 = -(tx12*tx42*tx60*tx797*tx800); 		/* rule 554 */
	 tx409 = tx102*tx118*tx52*tx60*tx797; 		/* rule 555 */
	 tx410 = tx406 + tx407 + tx408 + tx409; 		/* rule 556 */
	 tx411 = 3*K*tx410*tx64; 		/* rule 557 */
	 tx412 = 6.*K*tx158*tx54*tx73; 		/* rule 558 */
	 ohy1x4 = tx411 + tx412; 		/* rule 559 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 0, ohy1x4);
	 tx413 = -tx387; 		/* rule 561 */
	 tx414 = -(tx12*tx39*tx60*tx797*tx800); 		/* rule 562 */
	 tx415 = tx118*tx12*tx17*tx52*tx60*tx797; 		/* rule 563 */
	 tx416 = tx413 + tx414 + tx415; 		/* rule 564 */
	 tx417 = 3*K*tx416*tx64; 		/* rule 565 */
	 tx418 = 6.*K*tx161*tx54*tx73; 		/* rule 566 */
	 ohy1y4 = tx417 + tx418; 		/* rule 567 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 1, ohy1y4);
	 tx419 = -tx400; 		/* rule 569 */
	 tx420 = -(tx12*tx37*tx60*tx797*tx800); 		/* rule 570 */
	 tx421 = tx118*tx12*tx30*tx52*tx60*tx797; 		/* rule 571 */
	 tx422 = tx67*tx794*tx797*tx800; 		/* rule 572 */
	 tx423 = tx419 + tx420 + tx421 + tx422; 		/* rule 573 */
	 tx424 = 3*K*tx423*tx64; 		/* rule 574 */
	 tx425 = 6.*K*tx164*tx54*tx73; 		/* rule 575 */
	 ohy1z4 = tx424 + tx425; 		/* rule 576 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 2, ohy1z4);
	 tx426 = tx17*tx794*tx797*tx800; 		/* rule 578 */
	 tx427 = tx22*tx52*tx60*tx7*tx800*tx86; 		/* rule 579 */
	 tx428 = -(tx7*tx78*tx794*tx800*tx86); 		/* rule 580 */
	 tx429 = -(tx22*tx60*tx797*tx800*tx87); 		/* rule 581 */
	 tx430 = tx426 + tx427 + tx428 + tx429; 		/* rule 582 */
	 tx431 = 3*K*tx430*tx64; 		/* rule 583 */
	 tx432 = 6.*K*tx54*tx81*tx90; 		/* rule 584 */
	 ohz1x2 = tx431 + tx432; 		/* rule 585 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 tx433 = tx393 + x3; 		/* rule 587 */
	 tx434 = tx433*tx794*tx797*tx800; 		/* rule 588 */
	 tx435 = tx15*tx22*tx52*tx60*tx800*tx86; 		/* rule 589 */
	 tx436 = -(tx15*tx78*tx794*tx800*tx86); 		/* rule 590 */
	 tx437 = -(tx22*tx60*tx797*tx800*tx95); 		/* rule 591 */
	 tx438 = tx434 + tx435 + tx436 + tx437; 		/* rule 592 */
	 tx439 = 3*K*tx438*tx64; 		/* rule 593 */
	 tx440 = 6.*K*tx54*tx81*tx98; 		/* rule 594 */
	 ohz1y2 = tx439 + tx440; 		/* rule 595 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx441 = -(tx103*tx22*tx60*tx797*tx800); 		/* rule 597 */
	 tx442 = tx22*tx27*tx52*tx60*tx800*tx86; 		/* rule 598 */
	 tx443 = -(tx27*tx78*tx794*tx800*tx86); 		/* rule 599 */
	 tx444 = tx441 + tx442 + tx443; 		/* rule 600 */
	 tx445 = 3*K*tx444*tx64; 		/* rule 601 */
	 tx446 = 6.*K*tx106*tx54*tx81; 		/* rule 602 */
	 ohz1z2 = tx445 + tx446; 		/* rule 603 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 tx447 = tx280 + y2; 		/* rule 605 */
	 tx448 = -tx427; 		/* rule 606 */
	 tx449 = tx447*tx794*tx797*tx800; 		/* rule 607 */
	 tx450 = -(tx117*tx22*tx60*tx797*tx800); 		/* rule 608 */
	 tx451 = tx5*tx60*tx78*tx797*tx800; 		/* rule 609 */
	 tx452 = -(tx118*tx22*tx52*tx60*tx797*tx8); 		/* rule 610 */
	 tx453 = tx118*tx78*tx794*tx797*tx8; 		/* rule 611 */
	 tx454 = tx7*tx78*tx794*tx800*tx86; 		/* rule 612 */
	 tx455 = tx316 + tx448 + tx449 + tx450 + tx451 + tx452 + tx453 + tx454; 		/* rule 613 */
	 tx456 = 3*K*tx455*tx64; 		/* rule 614 */
	 tx457 = 6.*K*tx123*tx54*tx81; 		/* rule 615 */
	 ohz1x3 = tx456 + tx457; 		/* rule 616 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 0, ohz1x3);
	 tx458 = tx66 + x4; 		/* rule 618 */
	 tx459 = -tx435; 		/* rule 619 */
	 tx460 = tx458*tx794*tx797*tx800; 		/* rule 620 */
	 tx461 = -(tx134*tx22*tx60*tx797*tx800); 		/* rule 621 */
	 tx462 = tx118*tx17*tx78*tx794*tx797; 		/* rule 622 */
	 tx463 = tx12*tx60*tx78*tx797*tx800; 		/* rule 623 */
	 tx464 = -(tx118*tx52*tx60*tx797*tx84); 		/* rule 624 */
	 tx465 = tx15*tx78*tx794*tx800*tx86; 		/* rule 625 */
	 tx466 = tx397 + tx459 + tx460 + tx461 + tx462 + tx463 + tx464 + tx465; 		/* rule 626 */
	 tx467 = 3*K*tx466*tx64; 		/* rule 627 */
	 tx468 = 6.*K*tx139*tx54*tx81; 		/* rule 628 */
	 ohz1y3 = tx467 + tx468; 		/* rule 629 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 1, ohz1y3);
	 tx469 = -tx442; 		/* rule 631 */
	 tx470 = -3*tx165*tx34*tx52*tx797*tx800; 		/* rule 632 */
	 tx471 = -(tx150*tx22*tx60*tx797*tx800); 		/* rule 633 */
	 tx472 = -(tx118*tx22*tx30*tx52*tx60*tx797); 		/* rule 634 */
	 tx473 = tx118*tx30*tx78*tx794*tx797; 		/* rule 635 */
	 tx474 = tx22*tx60*tx78*tx797*tx800; 		/* rule 636 */
	 tx475 = tx27*tx78*tx794*tx800*tx86; 		/* rule 637 */
	 tx476 = tx292 + tx469 + tx470 + tx471 + tx472 + tx473 + tx474 + tx475; 		/* rule 638 */
	 tx477 = 3*K*tx476*tx64; 		/* rule 639 */
	 tx478 = 6.*K*tx155*tx54*tx81; 		/* rule 640 */
	 ohz1z3 = tx477 + tx478; 		/* rule 641 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 2, ohz1z3);
	 tx479 = -tx453; 		/* rule 643 */
	 tx480 = -(tx22*tx42*tx60*tx797*tx800); 		/* rule 644 */
	 tx481 = tx75*tx794*tx797*tx800; 		/* rule 645 */
	 tx482 = tx118*tx22*tx52*tx60*tx797*tx8; 		/* rule 646 */
	 tx483 = tx479 + tx480 + tx481 + tx482; 		/* rule 647 */
	 tx484 = 3*K*tx483*tx64; 		/* rule 648 */
	 tx485 = 6.*K*tx158*tx54*tx81; 		/* rule 649 */
	 ohz1x4 = tx484 + tx485; 		/* rule 650 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 0, ohz1x4);
	 tx486 = -tx462; 		/* rule 652 */
	 tx487 = -(tx22*tx39*tx60*tx797*tx800); 		/* rule 653 */
	 tx488 = tx7*tx794*tx797*tx800; 		/* rule 654 */
	 tx489 = tx118*tx52*tx60*tx797*tx84; 		/* rule 655 */
	 tx490 = tx486 + tx487 + tx488 + tx489; 		/* rule 656 */
	 tx491 = 3*K*tx490*tx64; 		/* rule 657 */
	 tx492 = 6.*K*tx161*tx54*tx81; 		/* rule 658 */
	 ohz1y4 = tx491 + tx492; 		/* rule 659 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 1, ohz1y4);
	 tx493 = -tx473; 		/* rule 661 */
	 tx494 = -(tx22*tx37*tx60*tx797*tx800); 		/* rule 662 */
	 tx495 = tx118*tx22*tx30*tx52*tx60*tx797; 		/* rule 663 */
	 tx496 = tx493 + tx494 + tx495; 		/* rule 664 */
	 tx497 = 3*K*tx496*tx64; 		/* rule 665 */
	 tx498 = 6.*K*tx164*tx54*tx81; 		/* rule 666 */
	 ohz1z4 = tx497 + tx498; 		/* rule 667 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 2, ohz1z4);
	 tx499 = 3*tx15*tx185*tx52*tx7*tx794*tx800; 		/* rule 669 */
	 tx500 = -(tx15*tx794*tx800*tx86*tx87); 		/* rule 670 */
	 tx501 = -(tx7*tx794*tx800*tx86*tx95); 		/* rule 671 */
	 tx502 = tx499 + tx500 + tx501; 		/* rule 672 */
	 tx503 = 3*K*tx502*tx64; 		/* rule 673 */
	 tx504 = 6.*K*tx54*tx90*tx98; 		/* rule 674 */
	 ohx2y2 = tx503 + tx504; 		/* rule 675 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 tx505 = 3*tx185*tx27*tx52*tx7*tx794*tx800; 		/* rule 677 */
	 tx506 = -(tx103*tx7*tx794*tx800*tx86); 		/* rule 678 */
	 tx507 = -(tx27*tx794*tx800*tx86*tx87); 		/* rule 679 */
	 tx508 = tx505 + tx506 + tx507; 		/* rule 680 */
	 tx509 = 3*K*tx508*tx64; 		/* rule 681 */
	 tx510 = 6.*K*tx106*tx54*tx90; 		/* rule 682 */
	 ohx2z2 = tx509 + tx510; 		/* rule 683 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 tx511 = -3*tx185*tx40*tx52*tx794*tx800; 		/* rule 685 */
	 tx512 = tx52*tx794*tx800*tx86; 		/* rule 686 */
	 tx513 = -(tx117*tx7*tx794*tx800*tx86); 		/* rule 687 */
	 tx514 = -(tx118*tx52*tx7*tx794*tx8*tx86); 		/* rule 688 */
	 tx515 = tx5*tx60*tx797*tx800*tx87; 		/* rule 689 */
	 tx516 = tx118*tx794*tx797*tx8*tx87; 		/* rule 690 */
	 tx517 = tx7*tx794*tx800*tx86*tx87; 		/* rule 691 */
	 tx518 = tx289 + tx511 + tx512 + tx513 + tx514 + tx515 + tx516 + tx517; 		/* rule 692 */
	 tx519 = 3*K*tx518*tx64; 		/* rule 693 */
	 tx520 = 6.*K*tx123*tx54*tx90; 		/* rule 694 */
	 ohx2x3 = tx519 + tx520; 		/* rule 695 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 0, ohx2x3);
	 tx521 = tx91 + z4; 		/* rule 697 */
	 tx522 = tx521*tx794*tx797*tx800; 		/* rule 698 */
	 tx523 = -3*tx15*tx185*tx52*tx7*tx794*tx800; 		/* rule 699 */
	 tx524 = -(tx134*tx7*tx794*tx800*tx86); 		/* rule 700 */
	 tx525 = -(tx118*tx52*tx76*tx794*tx86); 		/* rule 701 */
	 tx526 = tx118*tx17*tx794*tx797*tx87; 		/* rule 702 */
	 tx527 = tx12*tx60*tx797*tx800*tx87; 		/* rule 703 */
	 tx528 = tx15*tx794*tx800*tx86*tx87; 		/* rule 704 */
	 tx529 = tx373 + tx522 + tx523 + tx524 + tx525 + tx526 + tx527 + tx528; 		/* rule 705 */
	 tx530 = 3*K*tx529*tx64; 		/* rule 706 */
	 tx531 = 6.*K*tx139*tx54*tx90; 		/* rule 707 */
	 ohx2y3 = tx530 + tx531; 		/* rule 708 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 1, ohx2y3);
	 tx532 = tx280 + y1; 		/* rule 710 */
	 tx533 = tx532*tx794*tx797*tx800; 		/* rule 711 */
	 tx534 = -3*tx185*tx27*tx52*tx7*tx794*tx800; 		/* rule 712 */
	 tx535 = -(tx150*tx7*tx794*tx800*tx86); 		/* rule 713 */
	 tx536 = -(tx118*tx30*tx52*tx7*tx794*tx86); 		/* rule 714 */
	 tx537 = tx118*tx30*tx794*tx797*tx87; 		/* rule 715 */
	 tx538 = tx22*tx60*tx797*tx800*tx87; 		/* rule 716 */
	 tx539 = tx27*tx794*tx800*tx86*tx87; 		/* rule 717 */
	 tx540 = tx448 + tx533 + tx534 + tx535 + tx536 + tx537 + tx538 + tx539; 		/* rule 718 */
	 tx541 = 3*K*tx540*tx64; 		/* rule 719 */
	 tx542 = 6.*K*tx155*tx54*tx90; 		/* rule 720 */
	 ohx2z3 = tx541 + tx542; 		/* rule 721 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 2, ohx2z3);
	 tx543 = -tx516; 		/* rule 723 */
	 tx544 = -(tx42*tx7*tx794*tx800*tx86); 		/* rule 724 */
	 tx545 = tx118*tx52*tx7*tx794*tx8*tx86; 		/* rule 725 */
	 tx546 = tx543 + tx544 + tx545; 		/* rule 726 */
	 tx547 = 3*K*tx546*tx64; 		/* rule 727 */
	 tx548 = 6.*K*tx158*tx54*tx90; 		/* rule 728 */
	 ohx2x4 = tx547 + tx548; 		/* rule 729 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 0, ohx2x4);
	 tx549 = tx22*tx794*tx797*tx800; 		/* rule 731 */
	 tx550 = -tx526; 		/* rule 732 */
	 tx551 = -(tx39*tx7*tx794*tx800*tx86); 		/* rule 733 */
	 tx552 = tx118*tx52*tx76*tx794*tx86; 		/* rule 734 */
	 tx553 = tx549 + tx550 + tx551 + tx552; 		/* rule 735 */
	 tx554 = 3*K*tx553*tx64; 		/* rule 736 */
	 tx555 = 6.*K*tx161*tx54*tx90; 		/* rule 737 */
	 ohx2y4 = tx554 + tx555; 		/* rule 738 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 1, ohx2y4);
	 tx556 = -tx537; 		/* rule 740 */
	 tx557 = tx794*tx797*tx800*tx83; 		/* rule 741 */
	 tx558 = -(tx37*tx7*tx794*tx800*tx86); 		/* rule 742 */
	 tx559 = tx118*tx30*tx52*tx7*tx794*tx86; 		/* rule 743 */
	 tx560 = tx556 + tx557 + tx558 + tx559; 		/* rule 744 */
	 tx561 = 3*K*tx560*tx64; 		/* rule 745 */
	 tx562 = 6.*K*tx164*tx54*tx90; 		/* rule 746 */
	 ohx2z4 = tx561 + tx562; 		/* rule 747 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 2, ohx2z4);
	 tx563 = 3*tx15*tx185*tx27*tx52*tx794*tx800; 		/* rule 749 */
	 tx564 = -(tx103*tx15*tx794*tx800*tx86); 		/* rule 750 */
	 tx565 = -(tx27*tx794*tx800*tx86*tx95); 		/* rule 751 */
	 tx566 = tx563 + tx564 + tx565; 		/* rule 752 */
	 tx567 = 3*K*tx566*tx64; 		/* rule 753 */
	 tx568 = 6.*K*tx106*tx54*tx98; 		/* rule 754 */
	 ohy2z2 = tx567 + tx568; 		/* rule 755 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	 tx569 = tx300 + z1; 		/* rule 757 */
	 tx570 = tx569*tx794*tx797*tx800; 		/* rule 758 */
	 tx571 = -(tx117*tx15*tx794*tx800*tx86); 		/* rule 759 */
	 tx572 = -(tx118*tx15*tx52*tx794*tx8*tx86); 		/* rule 760 */
	 tx573 = tx5*tx60*tx797*tx800*tx95; 		/* rule 761 */
	 tx574 = tx118*tx794*tx797*tx8*tx95; 		/* rule 762 */
	 tx575 = tx7*tx794*tx800*tx86*tx95; 		/* rule 763 */
	 tx576 = tx302 + tx523 + tx570 + tx571 + tx572 + tx573 + tx574 + tx575; 		/* rule 764 */
	 tx577 = 3*K*tx576*tx64; 		/* rule 765 */
	 tx578 = 6.*K*tx123*tx54*tx98; 		/* rule 766 */
	 ohy2x3 = tx577 + tx578; 		/* rule 767 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 0, ohy2x3);
	 tx579 = -3*tx185*tx32*tx52*tx794*tx800; 		/* rule 769 */
	 tx580 = -(tx134*tx15*tx794*tx800*tx86); 		/* rule 770 */
	 tx581 = -(tx118*tx15*tx17*tx52*tx794*tx86); 		/* rule 771 */
	 tx582 = tx118*tx17*tx794*tx797*tx95; 		/* rule 772 */
	 tx583 = tx12*tx60*tx797*tx800*tx95; 		/* rule 773 */
	 tx584 = tx15*tx794*tx800*tx86*tx95; 		/* rule 774 */
	 tx585 = tx383 + tx512 + tx579 + tx580 + tx581 + tx582 + tx583 + tx584; 		/* rule 775 */
	 tx586 = 3*K*tx585*tx64; 		/* rule 776 */
	 tx587 = 6.*K*tx139*tx54*tx98; 		/* rule 777 */
	 ohy2y3 = tx586 + tx587; 		/* rule 778 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 1, ohy2y3);
	 tx588 = tx99 + x4; 		/* rule 780 */
	 tx589 = -3*tx15*tx185*tx27*tx52*tx794*tx800; 		/* rule 781 */
	 tx590 = tx588*tx794*tx797*tx800; 		/* rule 782 */
	 tx591 = -(tx15*tx150*tx794*tx800*tx86); 		/* rule 783 */
	 tx592 = -(tx118*tx52*tx58*tx794*tx86); 		/* rule 784 */
	 tx593 = tx118*tx30*tx794*tx797*tx95; 		/* rule 785 */
	 tx594 = tx22*tx60*tx797*tx800*tx95; 		/* rule 786 */
	 tx595 = tx27*tx794*tx800*tx86*tx95; 		/* rule 787 */
	 tx596 = tx459 + tx589 + tx590 + tx591 + tx592 + tx593 + tx594 + tx595; 		/* rule 788 */
	 tx597 = 3*K*tx596*tx64; 		/* rule 789 */
	 tx598 = 6.*K*tx155*tx54*tx98; 		/* rule 790 */
	 ohy2z3 = tx597 + tx598; 		/* rule 791 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 2, ohy2z3);
	 tx599 = -tx574; 		/* rule 793 */
	 tx600 = -(tx15*tx42*tx794*tx800*tx86); 		/* rule 794 */
	 tx601 = tx118*tx15*tx52*tx794*tx8*tx86; 		/* rule 795 */
	 tx602 = tx794*tx797*tx800*tx92; 		/* rule 796 */
	 tx603 = tx599 + tx600 + tx601 + tx602; 		/* rule 797 */
	 tx604 = 3*K*tx603*tx64; 		/* rule 798 */
	 tx605 = 6.*K*tx158*tx54*tx98; 		/* rule 799 */
	 ohy2x4 = tx604 + tx605; 		/* rule 800 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 0, ohy2x4);
	 tx606 = -tx582; 		/* rule 802 */
	 tx607 = -(tx15*tx39*tx794*tx800*tx86); 		/* rule 803 */
	 tx608 = tx118*tx15*tx17*tx52*tx794*tx86; 		/* rule 804 */
	 tx609 = tx606 + tx607 + tx608; 		/* rule 805 */
	 tx610 = 3*K*tx609*tx64; 		/* rule 806 */
	 tx611 = 6.*K*tx161*tx54*tx98; 		/* rule 807 */
	 ohy2y4 = tx610 + tx611; 		/* rule 808 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 1, ohy2y4);
	 tx612 = tx5*tx794*tx797*tx800; 		/* rule 810 */
	 tx613 = -tx593; 		/* rule 811 */
	 tx614 = -(tx15*tx37*tx794*tx800*tx86); 		/* rule 812 */
	 tx615 = tx118*tx52*tx58*tx794*tx86; 		/* rule 813 */
	 tx616 = tx612 + tx613 + tx614 + tx615; 		/* rule 814 */
	 tx617 = 3*K*tx616*tx64; 		/* rule 815 */
	 tx618 = 6.*K*tx164*tx54*tx98; 		/* rule 816 */
	 ohy2z4 = tx617 + tx618; 		/* rule 817 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 2, ohy2z4);
	 tx619 = tx82 + y4; 		/* rule 819 */
	 tx620 = tx103*tx5*tx60*tx797*tx800; 		/* rule 820 */
	 tx621 = tx619*tx794*tx797*tx800; 		/* rule 821 */
	 tx622 = tx103*tx118*tx794*tx797*tx8; 		/* rule 822 */
	 tx623 = -(tx117*tx27*tx794*tx800*tx86); 		/* rule 823 */
	 tx624 = -(tx118*tx52*tx69*tx794*tx86); 		/* rule 824 */
	 tx625 = tx103*tx7*tx794*tx800*tx86; 		/* rule 825 */
	 tx626 = tx314 + tx534 + tx620 + tx621 + tx622 + tx623 + tx624 + tx625; 		/* rule 826 */
	 tx627 = 6.*K*tx106*tx123*tx54; 		/* rule 827 */
	 tx628 = 3*K*tx626*tx64; 		/* rule 828 */
	 ohz2x3 = tx627 + tx628; 		/* rule 829 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 0, ohz2x3);
	 tx629 = tx393 + x1; 		/* rule 831 */
	 tx630 = tx103*tx118*tx17*tx794*tx797; 		/* rule 832 */
	 tx631 = tx103*tx12*tx60*tx797*tx800; 		/* rule 833 */
	 tx632 = tx629*tx794*tx797*tx800; 		/* rule 834 */
	 tx633 = tx103*tx15*tx794*tx800*tx86; 		/* rule 835 */
	 tx634 = -(tx134*tx27*tx794*tx800*tx86); 		/* rule 836 */
	 tx635 = -(tx118*tx17*tx27*tx52*tx794*tx86); 		/* rule 837 */
	 tx636 = tx395 + tx589 + tx630 + tx631 + tx632 + tx633 + tx634 + tx635; 		/* rule 838 */
	 tx637 = 6.*K*tx106*tx139*tx54; 		/* rule 839 */
	 tx638 = 3*K*tx636*tx64; 		/* rule 840 */
	 ohz2y3 = tx637 + tx638; 		/* rule 841 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 1, ohz2y3);
	 tx639 = tx103*tx118*tx30*tx794*tx797; 		/* rule 843 */
	 tx640 = -3*tx185*tx35*tx52*tx794*tx800; 		/* rule 844 */
	 tx641 = tx103*tx22*tx60*tx797*tx800; 		/* rule 845 */
	 tx642 = tx103*tx27*tx794*tx800*tx86; 		/* rule 846 */
	 tx643 = -(tx150*tx27*tx794*tx800*tx86); 		/* rule 847 */
	 tx644 = -(tx118*tx27*tx30*tx52*tx794*tx86); 		/* rule 848 */
	 tx645 = tx469 + tx512 + tx639 + tx640 + tx641 + tx642 + tx643 + tx644; 		/* rule 849 */
	 tx646 = 6.*K*tx106*tx155*tx54; 		/* rule 850 */
	 tx647 = 3*K*tx64*tx645; 		/* rule 851 */
	 ohz2z3 = tx646 + tx647; 		/* rule 852 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 2, ohz2z3);
	 tx648 = tx12*tx794*tx797*tx800; 		/* rule 854 */
	 tx649 = -tx622; 		/* rule 855 */
	 tx650 = -(tx27*tx42*tx794*tx800*tx86); 		/* rule 856 */
	 tx651 = tx118*tx52*tx69*tx794*tx86; 		/* rule 857 */
	 tx652 = tx648 + tx649 + tx650 + tx651; 		/* rule 858 */
	 tx653 = 6.*K*tx106*tx158*tx54; 		/* rule 859 */
	 tx654 = 3*K*tx64*tx652; 		/* rule 860 */
	 ohz2x4 = tx653 + tx654; 		/* rule 861 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 0, ohz2x4);
	 tx655 = tx100*tx794*tx797*tx800; 		/* rule 863 */
	 tx656 = -tx630; 		/* rule 864 */
	 tx657 = -(tx27*tx39*tx794*tx800*tx86); 		/* rule 865 */
	 tx658 = tx118*tx17*tx27*tx52*tx794*tx86; 		/* rule 866 */
	 tx659 = tx655 + tx656 + tx657 + tx658; 		/* rule 867 */
	 tx660 = 6.*K*tx106*tx161*tx54; 		/* rule 868 */
	 tx661 = 3*K*tx64*tx659; 		/* rule 869 */
	 ohz2y4 = tx660 + tx661; 		/* rule 870 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 1, ohz2y4);
	 tx662 = -tx639; 		/* rule 872 */
	 tx663 = -(tx27*tx37*tx794*tx800*tx86); 		/* rule 873 */
	 tx664 = tx118*tx27*tx30*tx52*tx794*tx86; 		/* rule 874 */
	 tx665 = tx662 + tx663 + tx664; 		/* rule 875 */
	 tx666 = 6.*K*tx106*tx164*tx54; 		/* rule 876 */
	 tx667 = 3*K*tx64*tx665; 		/* rule 877 */
	 ohz2z4 = tx666 + tx667; 		/* rule 878 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 2, ohz2z4);
	 tx668 = tx117*tx118*tx17*tx794*tx797; 		/* rule 880 */
	 tx669 = tx117*tx12*tx60*tx797*tx800; 		/* rule 881 */
	 tx670 = tx134*tx5*tx60*tx797*tx800; 		/* rule 882 */
	 tx671 = tx118*tx134*tx794*tx797*tx8; 		/* rule 883 */
	 tx672 = 3*tx17*tx205*tx52*tx794*tx797*tx8; 		/* rule 884 */
	 tx673 = tx117*tx15*tx794*tx800*tx86; 		/* rule 885 */
	 tx674 = tx134*tx7*tx794*tx800*tx86; 		/* rule 886 */
	 tx675 = tx255 + tx274 + tx334 + tx353 + tx409 + tx499 + tx552 + tx601 + tx668 + tx669 + tx670 + tx671 + tx672 + tx673 + tx674; 		/* rule 887 */
	 tx676 = 6.*K*tx123*tx139*tx54; 		/* rule 888 */
	 tx677 = 3*K*tx64*tx675; 		/* rule 889 */
	 ohx3y3 = tx676 + tx677; 		/* rule 890 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 1, ohx3y3);
	 tx678 = tx117*tx118*tx30*tx794*tx797; 		/* rule 892 */
	 tx679 = tx117*tx22*tx60*tx797*tx800; 		/* rule 893 */
	 tx680 = tx150*tx5*tx60*tx797*tx800; 		/* rule 894 */
	 tx681 = tx118*tx150*tx794*tx797*tx8; 		/* rule 895 */
	 tx682 = 3*tx205*tx30*tx52*tx794*tx797*tx8; 		/* rule 896 */
	 tx683 = tx117*tx27*tx794*tx800*tx86; 		/* rule 897 */
	 tx684 = tx150*tx7*tx794*tx800*tx86; 		/* rule 898 */
	 tx685 = tx261 + tx284 + tx341 + tx427 + tx482 + tx505 + tx559 + tx651 + tx678 + tx679 + tx680 + tx681 + tx682 + tx683 + tx684; 		/* rule 899 */
	 tx686 = 6.*K*tx123*tx155*tx54; 		/* rule 900 */
	 tx687 = 3*K*tx64*tx685; 		/* rule 901 */
	 ohx3z3 = tx686 + tx687; 		/* rule 902 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 2, ohx3z3);
	 tx688 = tx118*tx48*tx794*tx797; 		/* rule 904 */
	 tx689 = tx118*tx52*tx794*tx797; 		/* rule 905 */
	 tx690 = -3*tx205*tx41*tx52*tx794*tx797; 		/* rule 906 */
	 tx691 = tx42*tx5*tx60*tx797*tx800; 		/* rule 907 */
	 tx692 = -(tx117*tx118*tx794*tx797*tx8); 		/* rule 908 */
	 tx693 = tx42*tx7*tx794*tx800*tx86; 		/* rule 909 */
	 tx694 = tx294 + tx514 + tx688 + tx689 + tx690 + tx691 + tx692 + tx693; 		/* rule 910 */
	 tx695 = 6.*K*tx123*tx158*tx54; 		/* rule 911 */
	 tx696 = 3*K*tx64*tx694; 		/* rule 912 */
	 ohx3x4 = tx695 + tx696; 		/* rule 913 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 0, ohx3x4);
	 tx697 = tx108*tx794*tx797*tx800; 		/* rule 915 */
	 tx698 = tx39*tx5*tx60*tx797*tx800; 		/* rule 916 */
	 tx699 = -tx668; 		/* rule 917 */
	 tx700 = tx118*tx39*tx794*tx797*tx8; 		/* rule 918 */
	 tx701 = -3*tx17*tx205*tx52*tx794*tx797*tx8; 		/* rule 919 */
	 tx702 = tx39*tx7*tx794*tx800*tx86; 		/* rule 920 */
	 tx703 = tx306 + tx525 + tx697 + tx698 + tx699 + tx700 + tx701 + tx702; 		/* rule 921 */
	 tx704 = 6.*K*tx123*tx161*tx54; 		/* rule 922 */
	 tx705 = 3*K*tx64*tx703; 		/* rule 923 */
	 ohx3y4 = tx704 + tx705; 		/* rule 924 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 1, ohx3y4);
	 tx706 = tx107*tx794*tx797*tx800; 		/* rule 926 */
	 tx707 = tx37*tx5*tx60*tx797*tx800; 		/* rule 927 */
	 tx708 = -tx678; 		/* rule 928 */
	 tx709 = tx118*tx37*tx794*tx797*tx8; 		/* rule 929 */
	 tx710 = -3*tx205*tx30*tx52*tx794*tx797*tx8; 		/* rule 930 */
	 tx711 = tx37*tx7*tx794*tx800*tx86; 		/* rule 931 */
	 tx712 = tx321 + tx536 + tx706 + tx707 + tx708 + tx709 + tx710 + tx711; 		/* rule 932 */
	 tx713 = 6.*K*tx123*tx164*tx54; 		/* rule 933 */
	 tx714 = 3*K*tx64*tx712; 		/* rule 934 */
	 ohx3z4 = tx713 + tx714; 		/* rule 935 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 2, ohx3z4);
	 tx715 = tx118*tx150*tx17*tx794*tx797; 		/* rule 937 */
	 tx716 = tx118*tx134*tx30*tx794*tx797; 		/* rule 938 */
	 tx717 = 3*tx17*tx205*tx30*tx52*tx794*tx797; 		/* rule 939 */
	 tx718 = tx12*tx150*tx60*tx797*tx800; 		/* rule 940 */
	 tx719 = tx134*tx22*tx60*tx797*tx800; 		/* rule 941 */
	 tx720 = tx15*tx150*tx794*tx800*tx86; 		/* rule 942 */
	 tx721 = tx134*tx27*tx794*tx800*tx86; 		/* rule 943 */
	 tx722 = tx345 + tx367 + tx421 + tx435 + tx489 + tx563 + tx615 + tx658 + tx715 + tx716 + tx717 + tx718 + tx719 + tx720 + tx721; 		/* rule 944 */
	 tx723 = 6.*K*tx139*tx155*tx54; 		/* rule 945 */
	 tx724 = 3*K*tx64*tx722; 		/* rule 946 */
	 ohy3z3 = tx723 + tx724; 		/* rule 947 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 2, ohy3z3);
	 tx725 = tx118*tx17*tx42*tx794*tx797; 		/* rule 949 */
	 tx726 = tx125*tx794*tx797*tx800; 		/* rule 950 */
	 tx727 = tx12*tx42*tx60*tx797*tx800; 		/* rule 951 */
	 tx728 = -tx671; 		/* rule 952 */
	 tx729 = tx15*tx42*tx794*tx800*tx86; 		/* rule 953 */
	 tx730 = tx376 + tx572 + tx701 + tx725 + tx726 + tx727 + tx728 + tx729; 		/* rule 954 */
	 tx731 = 6.*K*tx139*tx158*tx54; 		/* rule 955 */
	 tx732 = 3*K*tx64*tx730; 		/* rule 956 */
	 ohy3x4 = tx731 + tx732; 		/* rule 957 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 0, ohy3x4);
	 tx733 = -(tx118*tx134*tx17*tx794*tx797); 		/* rule 959 */
	 tx734 = tx118*tx45*tx794*tx797; 		/* rule 960 */
	 tx735 = -3*tx205*tx33*tx52*tx794*tx797; 		/* rule 961 */
	 tx736 = tx12*tx39*tx60*tx797*tx800; 		/* rule 962 */
	 tx737 = tx15*tx39*tx794*tx800*tx86; 		/* rule 963 */
	 tx738 = tx386 + tx581 + tx689 + tx733 + tx734 + tx735 + tx736 + tx737; 		/* rule 964 */
	 tx739 = 6.*K*tx139*tx161*tx54; 		/* rule 965 */
	 tx740 = 3*K*tx64*tx738; 		/* rule 966 */
	 ohy3y4 = tx739 + tx740; 		/* rule 967 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 1, ohy3y4);
	 tx741 = tx118*tx17*tx37*tx794*tx797; 		/* rule 969 */
	 tx742 = tx124*tx794*tx797*tx800; 		/* rule 970 */
	 tx743 = -3*tx17*tx205*tx30*tx52*tx794*tx797; 		/* rule 971 */
	 tx744 = tx12*tx37*tx60*tx797*tx800; 		/* rule 972 */
	 tx745 = -tx716; 		/* rule 973 */
	 tx746 = tx15*tx37*tx794*tx800*tx86; 		/* rule 974 */
	 tx747 = tx399 + tx592 + tx741 + tx742 + tx743 + tx744 + tx745 + tx746; 		/* rule 975 */
	 tx748 = 6.*K*tx139*tx164*tx54; 		/* rule 976 */
	 tx749 = 3*K*tx64*tx747; 		/* rule 977 */
	 ohy3z4 = tx748 + tx749; 		/* rule 978 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 2, ohy3z4);
	 tx750 = tx118*tx30*tx42*tx794*tx797; 		/* rule 980 */
	 tx751 = tx141*tx794*tx797*tx800; 		/* rule 981 */
	 tx752 = tx22*tx42*tx60*tx797*tx800; 		/* rule 982 */
	 tx753 = -tx681; 		/* rule 983 */
	 tx754 = tx27*tx42*tx794*tx800*tx86; 		/* rule 984 */
	 tx755 = tx452 + tx624 + tx710 + tx750 + tx751 + tx752 + tx753 + tx754; 		/* rule 985 */
	 tx756 = 6.*K*tx155*tx158*tx54; 		/* rule 986 */
	 tx757 = 3*K*tx64*tx755; 		/* rule 987 */
	 ohz3x4 = tx756 + tx757; 		/* rule 988 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 0, ohz3x4);
	 tx758 = tx118*tx30*tx39*tx794*tx797; 		/* rule 990 */
	 tx759 = tx140*tx794*tx797*tx800; 		/* rule 991 */
	 tx760 = tx22*tx39*tx60*tx797*tx800; 		/* rule 992 */
	 tx761 = -tx715; 		/* rule 993 */
	 tx762 = tx27*tx39*tx794*tx800*tx86; 		/* rule 994 */
	 tx763 = tx464 + tx635 + tx743 + tx758 + tx759 + tx760 + tx761 + tx762; 		/* rule 995 */
	 tx764 = 6.*K*tx155*tx161*tx54; 		/* rule 996 */
	 tx765 = 3*K*tx64*tx763; 		/* rule 997 */
	 ohz3y4 = tx764 + tx765; 		/* rule 998 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 1, ohz3y4);
	 tx766 = -(tx118*tx150*tx30*tx794*tx797); 		/* rule 1000 */
	 tx767 = tx118*tx43*tx794*tx797; 		/* rule 1001 */
	 tx768 = -3*tx205*tx36*tx52*tx794*tx797; 		/* rule 1002 */
	 tx769 = tx22*tx37*tx60*tx797*tx800; 		/* rule 1003 */
	 tx770 = tx27*tx37*tx794*tx800*tx86; 		/* rule 1004 */
	 tx771 = tx472 + tx644 + tx689 + tx766 + tx767 + tx768 + tx769 + tx770; 		/* rule 1005 */
	 tx772 = 6.*K*tx155*tx164*tx54; 		/* rule 1006 */
	 tx773 = 3*K*tx64*tx771; 		/* rule 1007 */
	 ohz3z4 = tx772 + tx773; 		/* rule 1008 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 2, ohz3z4);
	 tx774 = -tx700; 		/* rule 1010 */
	 tx775 = -tx725; 		/* rule 1011 */
	 tx776 = tx672 + tx774 + tx775; 		/* rule 1012 */
	 tx777 = 6.*K*tx158*tx161*tx54; 		/* rule 1013 */
	 tx778 = 3*K*tx64*tx776; 		/* rule 1014 */
	 ohx4y4 = tx777 + tx778; 		/* rule 1015 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 1, ohx4y4);
	 tx779 = -tx709; 		/* rule 1017 */
	 tx780 = -tx750; 		/* rule 1018 */
	 tx781 = tx682 + tx779 + tx780; 		/* rule 1019 */
	 tx782 = 6.*K*tx158*tx164*tx54; 		/* rule 1020 */
	 tx783 = 3*K*tx64*tx781; 		/* rule 1021 */
	 ohx4z4 = tx782 + tx783; 		/* rule 1022 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 2, ohx4z4);
	 tx784 = -tx741; 		/* rule 1024 */
	 tx785 = -tx758; 		/* rule 1025 */
	 tx786 = tx717 + tx784 + tx785; 		/* rule 1026 */
	 tx787 = 6.*K*tx161*tx164*tx54; 		/* rule 1027 */
	 tx788 = 3*K*tx64*tx786; 		/* rule 1028 */
	 ohy4z4 = tx787 + tx788; 		/* rule 1029 */
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
