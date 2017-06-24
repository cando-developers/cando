// C-code
	NONBONDRB_SET_PARAMETER(I1);
	NONBONDRB_SET_PARAMETER(I2);
	NONBONDRB_SET_PARAMETER(dQ1Q2);
	NONBONDRB_SET_PARAMETER(dA);
	NONBONDRB_SET_PARAMETER(dC);
	NONBONDRB_SET_POSITION(ak,I1,0);
	NONBONDRB_SET_POSITION(bk,I1,1);
	NONBONDRB_SET_POSITION(ck,I1,2);
	NONBONDRB_SET_POSITION(dk,I1,3);
	NONBONDRB_SET_POSITION(xk,I1,4);
	NONBONDRB_SET_POSITION(yk,I1,5);
	NONBONDRB_SET_POSITION(zk,I1,6);
	NONBONDRB_SET_POSITION(al,I2,0);
	NONBONDRB_SET_POSITION(bl,I2,1);
	NONBONDRB_SET_POSITION(cl,I2,2);
	NONBONDRB_SET_POSITION(dl,I2,3);
	NONBONDRB_SET_POSITION(xl,I2,4);
	NONBONDRB_SET_POSITION(yl,I2,5);
	NONBONDRB_SET_POSITION(zl,I2,6);
	NONBONDRB_SET_POINT(xh1,ea1,getX());
	NONBONDRB_SET_POINT(yh1,ea1,getY());
	NONBONDRB_SET_POINT(zh1,ea1,getZ());
	NONBONDRB_SET_POINT(xh2,ea2,getX());
	NONBONDRB_SET_POINT(yh2,ea2,getY());
	NONBONDRB_SET_POINT(zh2,ea2,getZ());
	 tx216 = power2(bk); 		/* rule 26 */
	 tx217 = power2(bl); 		/* rule 27 */
	 tx218 = power2(ck); 		/* rule 28 */
	 tx219 = power2(cl); 		/* rule 29 */
	 tx220 = power2(dk); 		/* rule 30 */
	 tx221 = power2(dl); 		/* rule 31 */
	 tx222 = power2(ak); 		/* rule 32 */
	 tx223 = power2(al); 		/* rule 33 */
	 tzz452 = -2.*ak; 		/* rule 34 */
	 tx224 = bk*tzz452; 		/* rule 35 */
	 tzz451 = -2.*al; 		/* rule 36 */
	 tx225 = bl*tzz451; 		/* rule 37 */
	 tx226 = ck*tzz452; 		/* rule 38 */
	 tzz454 = 2.*ck; 		/* rule 39 */
	 tx227 = ak*tzz454; 		/* rule 40 */
	 tx228 = bk*tzz454; 		/* rule 41 */
	 tx229 = cl*tzz451; 		/* rule 42 */
	 tzz453 = 2.*cl; 		/* rule 43 */
	 tx230 = al*tzz453; 		/* rule 44 */
	 tx231 = bl*tzz453; 		/* rule 45 */
	 tx232 = dk*tzz452; 		/* rule 46 */
	 tzz450 = 2.*dk; 		/* rule 47 */
	 tx233 = ak*tzz450; 		/* rule 48 */
	 tx234 = bk*tzz450; 		/* rule 49 */
	 tx235 = ck*tzz450; 		/* rule 50 */
	 tx236 = dl*tzz451; 		/* rule 51 */
	 tzz455 = 2.*al; 		/* rule 52 */
	 tx237 = dl*tzz455; 		/* rule 53 */
	 tzz458 = 2.*bl; 		/* rule 54 */
	 tx238 = dl*tzz458; 		/* rule 55 */
	 tx239 = dl*tzz453; 		/* rule 56 */
	 tx240 = -tx216; 		/* rule 57 */
	 tx241 = -tx217; 		/* rule 58 */
	 tx242 = -tx218; 		/* rule 59 */
	 tx243 = -tx219; 		/* rule 60 */
	 tx244 = -tx220; 		/* rule 61 */
	 tx245 = -tx221; 		/* rule 62 */
	 tx246 = tx216 + tx218 + tx220 + tx222; 		/* rule 63 */
	 tx247 = tx217 + tx219 + tx221 + tx223; 		/* rule 64 */
	 tx248 = tx228 + tx232; 		/* rule 65 */
	 tx249 = tx228 + tx233; 		/* rule 66 */
	 tx250 = tx226 + tx234; 		/* rule 67 */
	 tx251 = tx227 + tx234; 		/* rule 68 */
	 tx252 = tx224 + tx235; 		/* rule 69 */
	 tx253 = tx231 + tx236; 		/* rule 70 */
	 tx254 = tx231 + tx237; 		/* rule 71 */
	 tx255 = tx229 + tx238; 		/* rule 72 */
	 tx256 = tx230 + tx238; 		/* rule 73 */
	 tx257 = tx225 + tx239; 		/* rule 74 */
	 tzz476 = tx222 + tx240; 		/* rule 75 */
	 tx258 = tx220 + tx242 + tzz476; 		/* rule 76 */
	 tzz475 = tx223 + tx241; 		/* rule 77 */
	 tx259 = tx221 + tx243 + tzz475; 		/* rule 78 */
	 tx260 = tx218 + tx244 + tzz476; 		/* rule 79 */
	 tx261 = tx219 + tx245 + tzz475; 		/* rule 80 */
	 tx262 = tx246*xh1; 		/* rule 81 */
	 tx263 = tx249*xh1; 		/* rule 82 */
	 tx264 = tx250*xh1; 		/* rule 83 */
	 tzz463 = -xh2; 		/* rule 84 */
	 tx265 = tx247*tzz463; 		/* rule 85 */
	 tx266 = tx254*tzz463; 		/* rule 86 */
	 tx267 = tx255*tzz463; 		/* rule 87 */
	 tx268 = -xl; 		/* rule 88 */
	 tx269 = tx248*yh1; 		/* rule 89 */
	 tx270 = tx252*yh1; 		/* rule 90 */
	 tx271 = tx260*yh1; 		/* rule 91 */
	 tzz462 = -yh2; 		/* rule 92 */
	 tx272 = tx253*tzz462; 		/* rule 93 */
	 tx273 = tx257*tzz462; 		/* rule 94 */
	 tx274 = tx261*tzz462; 		/* rule 95 */
	 tx275 = -yl; 		/* rule 96 */
	 tx276 = tx251*zh1; 		/* rule 97 */
	 tx277 = tx252*zh1; 		/* rule 98 */
	 tx278 = tx258*zh1; 		/* rule 99 */
	 tzz461 = -zh2; 		/* rule 100 */
	 tx279 = tx256*tzz461; 		/* rule 101 */
	 tx280 = tx257*tzz461; 		/* rule 102 */
	 tx281 = tx259*tzz461; 		/* rule 103 */
	 tx282 = -zl; 		/* rule 104 */
	 tx283 = tx262 + tx265 + tx268 + tx269 + tx272 + tx276 + tx279 + xk; 		/* rule 105 */
	 tx284 = tx263 + tx266 + tx271 + tx274 + tx275 + tx277 + tx280 + yk; 		/* rule 106 */
	 tx285 = tx264 + tx267 + tx270 + tx273 + tx278 + tx281 + tx282 + zk; 		/* rule 107 */
	 tx286 = power2(tx283); 		/* rule 108 */
	 tx287 = power2(tx284); 		/* rule 109 */
	 tx288 = power2(tx285); 		/* rule 110 */
	 tx289 = tx286 + tx287 + tx288; 		/* rule 111 */
	 tx472 = reciprocal(tx289); 		/* rule 112 */
	 tx439 = power2(tx472); 		/* rule 113 */
	 tx440 = tx472; 		/* rule 114 */
	 tx431 = tx439*tx440; 		/* rule 115 */
	 tx290 = power2(tx431); 		/* rule 116 */
	 tx432 = power2(tx440); 		/* rule 117 */
	 tx433 = tx440; 		/* rule 118 */
	 tx291 = tx432*tx433; 		/* rule 119 */
	 tx434 = mysqrt(tx289); 		/* rule 120 */
	 tx292 = reciprocal(tx434); 		/* rule 121 */
	 tx293 = dA*tx290; 		/* rule 122 */
	 tx294 = -(dC*tx291); 		/* rule 123 */
	 tx295 = dQ1Q2*tx292; 		/* rule 124 */
	 Energy = tx293 + tx294 + tx295; 		/* rule 125 */
	NONBONDRB_ENERGY_ACCUMULATE(Energy);
	 tzz460 = 2.*ak; 		/* rule 127 */
	 tx296 = tzz460*xh1; 		/* rule 128 */
	 tx297 = -2.*ck*xh1; 		/* rule 129 */
	 tx298 = tzz450*xh1; 		/* rule 130 */
	 tx299 = tzz460*yh1; 		/* rule 131 */
	 tzz471 = -2.*yh1; 		/* rule 132 */
	 tx300 = bk*tzz471; 		/* rule 133 */
	 tx301 = dk*tzz471; 		/* rule 134 */
	 tx302 = tzz460*zh1; 		/* rule 135 */
	 tzz470 = -2.*zh1; 		/* rule 136 */
	 tx303 = bk*tzz470; 		/* rule 137 */
	 tx304 = tzz454*zh1; 		/* rule 138 */
	 tx305 = tx297 + tx300 + tx302; 		/* rule 139 */
	 tx306 = tx298 + tx299 + tx303; 		/* rule 140 */
	 tx307 = tx296 + tx301 + tx304; 		/* rule 141 */
	 tzz445 = 2.*tx285; 		/* rule 142 */
	 tx308 = tx305*tzz445; 		/* rule 143 */
	 tzz446 = 2.*tx284; 		/* rule 144 */
	 tx309 = tx306*tzz446; 		/* rule 145 */
	 tzz447 = 2.*tx283; 		/* rule 146 */
	 tx310 = tx307*tzz447; 		/* rule 147 */
	 tx435 = power2(tx439); 		/* rule 148 */
	 tx311 = tx431*tx435; 		/* rule 149 */
	 tx312 = power2(tx432); 		/* rule 150 */
	 tx436 = tx289; 		/* rule 151 */
	 tx437 = tx292; 		/* rule 152 */
	 tx438 = reciprocal(tx436); 		/* rule 153 */
	 tx313 = tx437*tx438; 		/* rule 154 */
	 tx314 = tx308 + tx309 + tx310; 		/* rule 155 */
	 tzz443 = dA*tx311; 		/* rule 156 */
	 tzz449 = -6.*tzz443; 		/* rule 157 */
	 tx315 = tx314*tzz449; 		/* rule 158 */
	 tzz442 = dC*tx312; 		/* rule 159 */
	 tzz444 = 3*tzz442; 		/* rule 160 */
	 tx316 = tx314*tzz444; 		/* rule 161 */
	 tzz441 = dQ1Q2*tx313; 		/* rule 162 */
	 tzz448 = -0.5*tzz441; 		/* rule 163 */
	 tx317 = tx314*tzz448; 		/* rule 164 */
	 gak = tx315 + tx316 + tx317; 		/* rule 165 */
	 fak = -gak; 		/* rule 166 */
	NONBONDRB_FORCE_ACCUMULATE(I1, 0, fak );
	 tzz459 = 2.*bk; 		/* rule 168 */
	 tx318 = tzz459*xh1; 		/* rule 169 */
	 tx319 = tzz454*xh1; 		/* rule 170 */
	 tx320 = tzz452*yh1; 		/* rule 171 */
	 tx321 = tzz454*yh1; 		/* rule 172 */
	 tx322 = tzz452*zh1; 		/* rule 173 */
	 tx323 = tzz450*zh1; 		/* rule 174 */
	 tzz474 = tx298 + tx320; 		/* rule 175 */
	 tx324 = tx303 + tzz474; 		/* rule 176 */
	 tx325 = tx300 + tx319 + tx322; 		/* rule 177 */
	 tx326 = tx318 + tx321 + tx323; 		/* rule 178 */
	 tx327 = tx324*tzz445; 		/* rule 179 */
	 tx328 = tx325*tzz446; 		/* rule 180 */
	 tx329 = tx326*tzz447; 		/* rule 181 */
	 tx330 = tx327 + tx328 + tx329; 		/* rule 182 */
	 tx331 = tx330*tzz449; 		/* rule 183 */
	 tx332 = tx330*tzz444; 		/* rule 184 */
	 tx333 = tx330*tzz448; 		/* rule 185 */
	 gbk = tx331 + tx332 + tx333; 		/* rule 186 */
	 fbk = -gbk; 		/* rule 187 */
	NONBONDRB_FORCE_ACCUMULATE(I1, 1, fbk );
	 tx334 = tzz452*xh1; 		/* rule 189 */
	 tx335 = tzz459*yh1; 		/* rule 190 */
	 tx336 = tzz450*yh1; 		/* rule 191 */
	 tx337 = ck*tzz470; 		/* rule 192 */
	 tx338 = tx302 + tx319 + tx335; 		/* rule 193 */
	 tx339 = tx334 + tx336 + tx337; 		/* rule 194 */
	 tx340 = tx326*tzz446; 		/* rule 195 */
	 tx341 = tx338*tzz447; 		/* rule 196 */
	 tx342 = tx339*tzz445; 		/* rule 197 */
	 tx343 = tx340 + tx341 + tx342; 		/* rule 198 */
	 tx344 = tx343*tzz449; 		/* rule 199 */
	 tx345 = tx343*tzz444; 		/* rule 200 */
	 tx346 = tx343*tzz448; 		/* rule 201 */
	 gck = tx344 + tx345 + tx346; 		/* rule 202 */
	 fck = -gck; 		/* rule 203 */
	NONBONDRB_FORCE_ACCUMULATE(I1, 2, fck );
	 tx347 = tzz459*zh1; 		/* rule 205 */
	 tx348 = tx347 + tzz474; 		/* rule 206 */
	 tx349 = tx307*tzz446; 		/* rule 207 */
	 tx350 = tx326*tzz445; 		/* rule 208 */
	 tx351 = tx348*tzz447; 		/* rule 209 */
	 tx352 = tx349 + tx350 + tx351; 		/* rule 210 */
	 tx353 = tx352*tzz449; 		/* rule 211 */
	 tx354 = tx352*tzz444; 		/* rule 212 */
	 tx355 = tx352*tzz448; 		/* rule 213 */
	 gdk = tx353 + tx354 + tx355; 		/* rule 214 */
	 fdk = -gdk; 		/* rule 215 */
	NONBONDRB_FORCE_ACCUMULATE(I1, 3, fdk );
	 tzz469 = -12.*tzz443; 		/* rule 217 */
	 tx356 = tx283*tzz469; 		/* rule 218 */
	 tzz457 = 6.*tzz442; 		/* rule 219 */
	 tx357 = tx283*tzz457; 		/* rule 220 */
	 tzz464 = -tzz441; 		/* rule 221 */
	 tx358 = tx283*tzz464; 		/* rule 222 */
	 gxk = tx356 + tx357 + tx358; 		/* rule 223 */
	 fxk = -gxk; 		/* rule 224 */
	NONBONDRB_FORCE_ACCUMULATE(I1, 4, fxk );
	 tx359 = tx284*tzz469; 		/* rule 226 */
	 tx360 = tx284*tzz457; 		/* rule 227 */
	 tx361 = tx284*tzz464; 		/* rule 228 */
	 gyk = tx359 + tx360 + tx361; 		/* rule 229 */
	 fyk = -gyk; 		/* rule 230 */
	NONBONDRB_FORCE_ACCUMULATE(I1, 5, fyk );
	 tx362 = tx285*tzz469; 		/* rule 232 */
	 tx363 = tx285*tzz457; 		/* rule 233 */
	 tx364 = tx285*tzz464; 		/* rule 234 */
	 gzk = tx362 + tx363 + tx364; 		/* rule 235 */
	 fzk = -gzk; 		/* rule 236 */
	NONBONDRB_FORCE_ACCUMULATE(I1, 6, fzk );
	 tx365 = tzz451*xh2; 		/* rule 238 */
	 tx366 = tzz453*xh2; 		/* rule 239 */
	 tzz467 = -2.*xh2; 		/* rule 240 */
	 tx367 = dl*tzz467; 		/* rule 241 */
	 tx368 = tzz451*yh2; 		/* rule 242 */
	 tx369 = tzz458*yh2; 		/* rule 243 */
	 tx370 = 2.*dl*yh2; 		/* rule 244 */
	 tx371 = tzz451*zh2; 		/* rule 245 */
	 tx372 = tzz458*zh2; 		/* rule 246 */
	 tzz465 = -2.*zh2; 		/* rule 247 */
	 tx373 = cl*tzz465; 		/* rule 248 */
	 tx374 = tx366 + tx369 + tx371; 		/* rule 249 */
	 tx375 = tx367 + tx368 + tx372; 		/* rule 250 */
	 tx376 = tx365 + tx370 + tx373; 		/* rule 251 */
	 tx377 = tx374*tzz445; 		/* rule 252 */
	 tx378 = tx375*tzz446; 		/* rule 253 */
	 tx379 = tx376*tzz447; 		/* rule 254 */
	 tx380 = tx377 + tx378 + tx379; 		/* rule 255 */
	 tx381 = tx380*tzz449; 		/* rule 256 */
	 tx382 = tx380*tzz444; 		/* rule 257 */
	 tx383 = tx380*tzz448; 		/* rule 258 */
	 gal = tx381 + tx382 + tx383; 		/* rule 259 */
	 fal = -gal; 		/* rule 260 */
	NONBONDRB_FORCE_ACCUMULATE(I2, 0, fal );
	 tx384 = bl*tzz467; 		/* rule 262 */
	 tx385 = cl*tzz467; 		/* rule 263 */
	 tx386 = tzz455*yh2; 		/* rule 264 */
	 tzz466 = -2.*yh2; 		/* rule 265 */
	 tx387 = cl*tzz466; 		/* rule 266 */
	 tx388 = tzz455*zh2; 		/* rule 267 */
	 tx389 = dl*tzz465; 		/* rule 268 */
	 tzz473 = tx367 + tx386; 		/* rule 269 */
	 tx390 = tx372 + tzz473; 		/* rule 270 */
	 tx391 = tx369 + tx385 + tx388; 		/* rule 271 */
	 tx392 = tx384 + tx387 + tx389; 		/* rule 272 */
	 tx393 = tx390*tzz445; 		/* rule 273 */
	 tx394 = tx391*tzz446; 		/* rule 274 */
	 tx395 = tx392*tzz447; 		/* rule 275 */
	 tx396 = tx393 + tx394 + tx395; 		/* rule 276 */
	 tx397 = tx396*tzz449; 		/* rule 277 */
	 tx398 = tx396*tzz444; 		/* rule 278 */
	 tx399 = tx396*tzz448; 		/* rule 279 */
	 gbl = tx397 + tx398 + tx399; 		/* rule 280 */
	 fbl = -gbl; 		/* rule 281 */
	NONBONDRB_FORCE_ACCUMULATE(I2, 1, fbl );
	 tx400 = tzz455*xh2; 		/* rule 283 */
	 tx401 = bl*tzz466; 		/* rule 284 */
	 tx402 = dl*tzz466; 		/* rule 285 */
	 tx403 = tzz453*zh2; 		/* rule 286 */
	 tx404 = tx371 + tx385 + tx401; 		/* rule 287 */
	 tx405 = tx400 + tx402 + tx403; 		/* rule 288 */
	 tx406 = tx392*tzz446; 		/* rule 289 */
	 tx407 = tx404*tzz447; 		/* rule 290 */
	 tx408 = tx405*tzz445; 		/* rule 291 */
	 tx409 = tx406 + tx407 + tx408; 		/* rule 292 */
	 tx410 = tx409*tzz449; 		/* rule 293 */
	 tx411 = tx409*tzz444; 		/* rule 294 */
	 tx412 = tx409*tzz448; 		/* rule 295 */
	 gcl = tx410 + tx411 + tx412; 		/* rule 296 */
	 fcl = -gcl; 		/* rule 297 */
	NONBONDRB_FORCE_ACCUMULATE(I2, 2, fcl );
	 tx413 = bl*tzz465; 		/* rule 299 */
	 tx414 = tx413 + tzz473; 		/* rule 300 */
	 tx415 = tx376*tzz446; 		/* rule 301 */
	 tx416 = tx392*tzz445; 		/* rule 302 */
	 tx417 = tx414*tzz447; 		/* rule 303 */
	 tx418 = tx415 + tx416 + tx417; 		/* rule 304 */
	 tx419 = tx418*tzz449; 		/* rule 305 */
	 tx420 = tx418*tzz444; 		/* rule 306 */
	 tx421 = tx418*tzz448; 		/* rule 307 */
	 gdl = tx419 + tx420 + tx421; 		/* rule 308 */
	 fdl = -gdl; 		/* rule 309 */
	NONBONDRB_FORCE_ACCUMULATE(I2, 3, fdl );
	 tzz456 = 12.*tzz443; 		/* rule 311 */
	 tx422 = tx283*tzz456; 		/* rule 312 */
	 tzz468 = -6.*tzz442; 		/* rule 313 */
	 tx423 = tx283*tzz468; 		/* rule 314 */
	 tx424 = tx283*tzz441; 		/* rule 315 */
	 gxl = tx422 + tx423 + tx424; 		/* rule 316 */
	 fxl = -gxl; 		/* rule 317 */
	NONBONDRB_FORCE_ACCUMULATE(I2, 4, fxl );
	 tx425 = tx284*tzz456; 		/* rule 319 */
	 tx426 = tx284*tzz468; 		/* rule 320 */
	 tx427 = tx284*tzz441; 		/* rule 321 */
	 gyl = tx425 + tx426 + tx427; 		/* rule 322 */
	 fyl = -gyl; 		/* rule 323 */
	NONBONDRB_FORCE_ACCUMULATE(I2, 5, fyl );
	 tx428 = tx285*tzz456; 		/* rule 325 */
	 tx429 = tx285*tzz468; 		/* rule 326 */
	 tx430 = tx285*tzz441; 		/* rule 327 */
	 gzl = tx428 + tx429 + tx430; 		/* rule 328 */
	 fzl = -gzl; 		/* rule 329 */
	NONBONDRB_FORCE_ACCUMULATE(I2, 6, fzl );
