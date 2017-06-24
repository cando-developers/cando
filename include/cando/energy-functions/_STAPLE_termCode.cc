// C-code
	STAPLE_SET_PARAMETER(ks,ks);
	STAPLE_SET_PARAMETER(r0,r0);
	STAPLE_SET_PARAMETER(I1,rigidBodyA);
	STAPLE_SET_PARAMETER(I2,rigidBodyB);
	STAPLE_SET_POSITION(ak,I1,0);
	STAPLE_SET_POSITION(bk,I1,1);
	STAPLE_SET_POSITION(ck,I1,2);
	STAPLE_SET_POSITION(dk,I1,3);
	STAPLE_SET_POSITION(xk,I1,4);
	STAPLE_SET_POSITION(yk,I1,5);
	STAPLE_SET_POSITION(zk,I1,6);
	STAPLE_SET_POSITION(al,I2,0);
	STAPLE_SET_POSITION(bl,I2,1);
	STAPLE_SET_POSITION(cl,I2,2);
	STAPLE_SET_POSITION(dl,I2,3);
	STAPLE_SET_POSITION(xl,I2,4);
	STAPLE_SET_POSITION(yl,I2,5);
	STAPLE_SET_POSITION(zl,I2,6);
	STAPLE_SET_POINT(xh1,pointA,getX());
	STAPLE_SET_POINT(yh1,pointA,getY());
	STAPLE_SET_POINT(zh1,pointA,getZ());
	STAPLE_SET_POINT(xh2,pointB,getX());
	STAPLE_SET_POINT(yh2,pointB,getY());
	STAPLE_SET_POINT(zh2,pointB,getZ());
	 tx170 = power2(bk); 		/* rule 25 */
	 tx171 = power2(bl); 		/* rule 26 */
	 tx172 = power2(ck); 		/* rule 27 */
	 tx173 = power2(cl); 		/* rule 28 */
	 tx174 = power2(dk); 		/* rule 29 */
	 tx175 = power2(dl); 		/* rule 30 */
	 tx176 = power2(ak); 		/* rule 31 */
	 tx177 = power2(al); 		/* rule 32 */
	 tzz346 = -2.*ak; 		/* rule 33 */
	 tx178 = bk*tzz346; 		/* rule 34 */
	 tzz345 = -2.*al; 		/* rule 35 */
	 tx179 = bl*tzz345; 		/* rule 36 */
	 tx180 = ck*tzz346; 		/* rule 37 */
	 tzz348 = 2.*ck; 		/* rule 38 */
	 tx181 = ak*tzz348; 		/* rule 39 */
	 tx182 = bk*tzz348; 		/* rule 40 */
	 tx183 = cl*tzz345; 		/* rule 41 */
	 tzz347 = 2.*cl; 		/* rule 42 */
	 tx184 = al*tzz347; 		/* rule 43 */
	 tx185 = bl*tzz347; 		/* rule 44 */
	 tx186 = dk*tzz346; 		/* rule 45 */
	 tzz344 = 2.*dk; 		/* rule 46 */
	 tx187 = ak*tzz344; 		/* rule 47 */
	 tx188 = bk*tzz344; 		/* rule 48 */
	 tx189 = ck*tzz344; 		/* rule 49 */
	 tx190 = dl*tzz345; 		/* rule 50 */
	 tzz349 = 2.*al; 		/* rule 51 */
	 tx191 = dl*tzz349; 		/* rule 52 */
	 tzz350 = 2.*bl; 		/* rule 53 */
	 tx192 = dl*tzz350; 		/* rule 54 */
	 tx193 = dl*tzz347; 		/* rule 55 */
	 tx194 = -tx170; 		/* rule 56 */
	 tx195 = -tx171; 		/* rule 57 */
	 tx196 = -tx172; 		/* rule 58 */
	 tx197 = -tx173; 		/* rule 59 */
	 tx198 = -tx174; 		/* rule 60 */
	 tx199 = -tx175; 		/* rule 61 */
	 tx200 = tx170 + tx172 + tx174 + tx176; 		/* rule 62 */
	 tx201 = tx171 + tx173 + tx175 + tx177; 		/* rule 63 */
	 tx202 = tx182 + tx186; 		/* rule 64 */
	 tx203 = tx182 + tx187; 		/* rule 65 */
	 tx204 = tx180 + tx188; 		/* rule 66 */
	 tx205 = tx181 + tx188; 		/* rule 67 */
	 tx206 = tx178 + tx189; 		/* rule 68 */
	 tx207 = tx185 + tx190; 		/* rule 69 */
	 tx208 = tx185 + tx191; 		/* rule 70 */
	 tx209 = tx183 + tx192; 		/* rule 71 */
	 tx210 = tx184 + tx192; 		/* rule 72 */
	 tx211 = tx179 + tx193; 		/* rule 73 */
	 tzz365 = tx176 + tx194; 		/* rule 74 */
	 tx212 = tx174 + tx196 + tzz365; 		/* rule 75 */
	 tzz364 = tx177 + tx195; 		/* rule 76 */
	 tx213 = tx175 + tx197 + tzz364; 		/* rule 77 */
	 tx214 = tx172 + tx198 + tzz365; 		/* rule 78 */
	 tx215 = tx173 + tx199 + tzz364; 		/* rule 79 */
	 tx216 = tx200*xh1; 		/* rule 80 */
	 tx217 = tx203*xh1; 		/* rule 81 */
	 tx218 = tx204*xh1; 		/* rule 82 */
	 tzz355 = -xh2; 		/* rule 83 */
	 tx219 = tx201*tzz355; 		/* rule 84 */
	 tx220 = tx208*tzz355; 		/* rule 85 */
	 tx221 = tx209*tzz355; 		/* rule 86 */
	 tx222 = -xl; 		/* rule 87 */
	 tx223 = tx202*yh1; 		/* rule 88 */
	 tx224 = tx206*yh1; 		/* rule 89 */
	 tx225 = tx214*yh1; 		/* rule 90 */
	 tzz354 = -yh2; 		/* rule 91 */
	 tx226 = tx207*tzz354; 		/* rule 92 */
	 tx227 = tx211*tzz354; 		/* rule 93 */
	 tx228 = tx215*tzz354; 		/* rule 94 */
	 tx229 = -yl; 		/* rule 95 */
	 tx230 = tx205*zh1; 		/* rule 96 */
	 tx231 = tx206*zh1; 		/* rule 97 */
	 tx232 = tx212*zh1; 		/* rule 98 */
	 tzz353 = -zh2; 		/* rule 99 */
	 tx233 = tx210*tzz353; 		/* rule 100 */
	 tx234 = tx211*tzz353; 		/* rule 101 */
	 tx235 = tx213*tzz353; 		/* rule 102 */
	 tx236 = -zl; 		/* rule 103 */
	 tx237 = tx216 + tx219 + tx222 + tx223 + tx226 + tx230 + tx233 + xk; 		/* rule 104 */
	 tx238 = tx217 + tx220 + tx225 + tx228 + tx229 + tx231 + tx234 + yk; 		/* rule 105 */
	 tx239 = tx218 + tx221 + tx224 + tx227 + tx232 + tx235 + tx236 + zk; 		/* rule 106 */
	 tx240 = power2(tx237); 		/* rule 107 */
	 tx241 = power2(tx238); 		/* rule 108 */
	 tx242 = power2(tx239); 		/* rule 109 */
	 tx243 = tx240 + tx241 + tx242; 		/* rule 110 */
	 tx244 = -r0; 		/* rule 111 */
	 tx245 = mysqrt(tx243); 		/* rule 112 */
	 tx246 = tx244 + tx245; 		/* rule 113 */
	 tx247 = power2(tx246); 		/* rule 114 */
	 Energy = ks*tx247; 		/* rule 115 */
	STAPLE_ENERGY_ACCUMULATE(Energy);
	 tzz352 = 2.*ak; 		/* rule 117 */
	 tx248 = tzz352*xh1; 		/* rule 118 */
	 tx249 = -2.*ck*xh1; 		/* rule 119 */
	 tx250 = tzz344*xh1; 		/* rule 120 */
	 tx251 = tzz352*yh1; 		/* rule 121 */
	 tzz361 = -2.*yh1; 		/* rule 122 */
	 tx252 = bk*tzz361; 		/* rule 123 */
	 tx253 = dk*tzz361; 		/* rule 124 */
	 tx254 = tzz352*zh1; 		/* rule 125 */
	 tzz360 = -2.*zh1; 		/* rule 126 */
	 tx255 = bk*tzz360; 		/* rule 127 */
	 tx256 = tzz348*zh1; 		/* rule 128 */
	 tx257 = tx249 + tx252 + tx254; 		/* rule 129 */
	 tx258 = tx250 + tx251 + tx255; 		/* rule 130 */
	 tx259 = tx248 + tx253 + tx256; 		/* rule 131 */
	 tzz341 = 2.*tx239; 		/* rule 132 */
	 tx260 = tx257*tzz341; 		/* rule 133 */
	 tzz342 = 2.*tx238; 		/* rule 134 */
	 tx261 = tx258*tzz342; 		/* rule 135 */
	 tzz343 = 2.*tx237; 		/* rule 136 */
	 tx262 = tx259*tzz343; 		/* rule 137 */
	 tx263 = reciprocal(tx245); 		/* rule 138 */
	 tx264 = tx260 + tx261 + tx262; 		/* rule 139 */
	 tzz339 = tx246*tx263; 		/* rule 140 */
	 tzz340 = ks*tzz339; 		/* rule 141 */
	 gak = tx264*tzz340; 		/* rule 142 */
	 fak = -gak; 		/* rule 143 */
	STAPLE_FORCE_ACCUMULATE(I1, 0, fak );
	 tzz351 = 2.*bk; 		/* rule 145 */
	 tx265 = tzz351*xh1; 		/* rule 146 */
	 tx266 = tzz348*xh1; 		/* rule 147 */
	 tx267 = tzz346*yh1; 		/* rule 148 */
	 tx268 = tzz348*yh1; 		/* rule 149 */
	 tx269 = tzz346*zh1; 		/* rule 150 */
	 tx270 = tzz344*zh1; 		/* rule 151 */
	 tzz363 = tx250 + tx267; 		/* rule 152 */
	 tx271 = tx255 + tzz363; 		/* rule 153 */
	 tx272 = tx252 + tx266 + tx269; 		/* rule 154 */
	 tx273 = tx265 + tx268 + tx270; 		/* rule 155 */
	 tx274 = tx271*tzz341; 		/* rule 156 */
	 tx275 = tx272*tzz342; 		/* rule 157 */
	 tx276 = tx273*tzz343; 		/* rule 158 */
	 tx277 = tx274 + tx275 + tx276; 		/* rule 159 */
	 gbk = tx277*tzz340; 		/* rule 160 */
	 fbk = -gbk; 		/* rule 161 */
	STAPLE_FORCE_ACCUMULATE(I1, 1, fbk );
	 tx278 = tzz346*xh1; 		/* rule 163 */
	 tx279 = tzz351*yh1; 		/* rule 164 */
	 tx280 = tzz344*yh1; 		/* rule 165 */
	 tx281 = ck*tzz360; 		/* rule 166 */
	 tx282 = tx254 + tx266 + tx279; 		/* rule 167 */
	 tx283 = tx278 + tx280 + tx281; 		/* rule 168 */
	 tx284 = tx273*tzz342; 		/* rule 169 */
	 tx285 = tx282*tzz343; 		/* rule 170 */
	 tx286 = tx283*tzz341; 		/* rule 171 */
	 tx287 = tx284 + tx285 + tx286; 		/* rule 172 */
	 gck = tx287*tzz340; 		/* rule 173 */
	 fck = -gck; 		/* rule 174 */
	STAPLE_FORCE_ACCUMULATE(I1, 2, fck );
	 tx288 = tzz351*zh1; 		/* rule 176 */
	 tx289 = tx288 + tzz363; 		/* rule 177 */
	 tx290 = tx259*tzz342; 		/* rule 178 */
	 tx291 = tx273*tzz341; 		/* rule 179 */
	 tx292 = tx289*tzz343; 		/* rule 180 */
	 tx293 = tx290 + tx291 + tx292; 		/* rule 181 */
	 gdk = tx293*tzz340; 		/* rule 182 */
	 fdk = -gdk; 		/* rule 183 */
	STAPLE_FORCE_ACCUMULATE(I1, 3, fdk );
	 gxk = tzz340*tzz343; 		/* rule 185 */
	 fxk = -gxk; 		/* rule 186 */
	STAPLE_FORCE_ACCUMULATE(I1, 4, fxk );
	 gyk = tzz340*tzz342; 		/* rule 188 */
	 fyk = -gyk; 		/* rule 189 */
	STAPLE_FORCE_ACCUMULATE(I1, 5, fyk );
	 gzk = tzz340*tzz341; 		/* rule 191 */
	 fzk = -gzk; 		/* rule 192 */
	STAPLE_FORCE_ACCUMULATE(I1, 6, fzk );
	 tx294 = tzz345*xh2; 		/* rule 194 */
	 tx295 = tzz347*xh2; 		/* rule 195 */
	 tzz358 = -2.*xh2; 		/* rule 196 */
	 tx296 = dl*tzz358; 		/* rule 197 */
	 tx297 = tzz345*yh2; 		/* rule 198 */
	 tx298 = tzz350*yh2; 		/* rule 199 */
	 tx299 = 2.*dl*yh2; 		/* rule 200 */
	 tx300 = tzz345*zh2; 		/* rule 201 */
	 tx301 = tzz350*zh2; 		/* rule 202 */
	 tzz356 = -2.*zh2; 		/* rule 203 */
	 tx302 = cl*tzz356; 		/* rule 204 */
	 tx303 = tx295 + tx298 + tx300; 		/* rule 205 */
	 tx304 = tx296 + tx297 + tx301; 		/* rule 206 */
	 tx305 = tx294 + tx299 + tx302; 		/* rule 207 */
	 tx306 = tx303*tzz341; 		/* rule 208 */
	 tx307 = tx304*tzz342; 		/* rule 209 */
	 tx308 = tx305*tzz343; 		/* rule 210 */
	 tx309 = tx306 + tx307 + tx308; 		/* rule 211 */
	 gal = tx309*tzz340; 		/* rule 212 */
	 fal = -gal; 		/* rule 213 */
	STAPLE_FORCE_ACCUMULATE(I2, 0, fal );
	 tx310 = bl*tzz358; 		/* rule 215 */
	 tx311 = cl*tzz358; 		/* rule 216 */
	 tx312 = tzz349*yh2; 		/* rule 217 */
	 tzz357 = -2.*yh2; 		/* rule 218 */
	 tx313 = cl*tzz357; 		/* rule 219 */
	 tx314 = tzz349*zh2; 		/* rule 220 */
	 tx315 = dl*tzz356; 		/* rule 221 */
	 tzz362 = tx296 + tx312; 		/* rule 222 */
	 tx316 = tx301 + tzz362; 		/* rule 223 */
	 tx317 = tx298 + tx311 + tx314; 		/* rule 224 */
	 tx318 = tx310 + tx313 + tx315; 		/* rule 225 */
	 tx319 = tx316*tzz341; 		/* rule 226 */
	 tx320 = tx317*tzz342; 		/* rule 227 */
	 tx321 = tx318*tzz343; 		/* rule 228 */
	 tx322 = tx319 + tx320 + tx321; 		/* rule 229 */
	 gbl = tx322*tzz340; 		/* rule 230 */
	 fbl = -gbl; 		/* rule 231 */
	STAPLE_FORCE_ACCUMULATE(I2, 1, fbl );
	 tx323 = tzz349*xh2; 		/* rule 233 */
	 tx324 = bl*tzz357; 		/* rule 234 */
	 tx325 = dl*tzz357; 		/* rule 235 */
	 tx326 = tzz347*zh2; 		/* rule 236 */
	 tx327 = tx300 + tx311 + tx324; 		/* rule 237 */
	 tx328 = tx323 + tx325 + tx326; 		/* rule 238 */
	 tx329 = tx318*tzz342; 		/* rule 239 */
	 tx330 = tx327*tzz343; 		/* rule 240 */
	 tx331 = tx328*tzz341; 		/* rule 241 */
	 tx332 = tx329 + tx330 + tx331; 		/* rule 242 */
	 gcl = tx332*tzz340; 		/* rule 243 */
	 fcl = -gcl; 		/* rule 244 */
	STAPLE_FORCE_ACCUMULATE(I2, 2, fcl );
	 tx333 = bl*tzz356; 		/* rule 246 */
	 tx334 = tx333 + tzz362; 		/* rule 247 */
	 tx335 = tx305*tzz342; 		/* rule 248 */
	 tx336 = tx318*tzz341; 		/* rule 249 */
	 tx337 = tx334*tzz343; 		/* rule 250 */
	 tx338 = tx335 + tx336 + tx337; 		/* rule 251 */
	 gdl = tx338*tzz340; 		/* rule 252 */
	 fdl = -gdl; 		/* rule 253 */
	STAPLE_FORCE_ACCUMULATE(I2, 3, fdl );
	 tzz359 = -2.*tzz340; 		/* rule 255 */
	 gxl = tx237*tzz359; 		/* rule 256 */
	 fxl = -gxl; 		/* rule 257 */
	STAPLE_FORCE_ACCUMULATE(I2, 4, fxl );
	 gyl = tx238*tzz359; 		/* rule 259 */
	 fyl = -gyl; 		/* rule 260 */
	STAPLE_FORCE_ACCUMULATE(I2, 5, fyl );
	 gzl = tx239*tzz359; 		/* rule 262 */
	 fzl = -gzl; 		/* rule 263 */
	STAPLE_FORCE_ACCUMULATE(I2, 6, fzl );
