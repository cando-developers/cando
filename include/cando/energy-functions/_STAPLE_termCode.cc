// C-code
	STAPLE_SET_PARAMETER(ks,ks);
	STAPLE_SET_PARAMETER(r0,r0);
	STAPLE_SET_PARAMETER(I1,rigidBodyK);
	STAPLE_SET_PARAMETER(I2,rigidBodyL);
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
	STAPLE_SET_POINT(pxk,pointK,getX());
	STAPLE_SET_POINT(pyk,pointK,getY());
	STAPLE_SET_POINT(pzk,pointK,getZ());
	STAPLE_SET_POINT(pxl,pointL,getX());
	STAPLE_SET_POINT(pyl,pointL,getY());
	STAPLE_SET_POINT(pzl,pointL,getZ());
	 tx278 = power2(ak); 		/* rule 25 */
	 tx279 = power2(al); 		/* rule 26 */
	 tx280 = power2(bk); 		/* rule 27 */
	 tx281 = power2(bl); 		/* rule 28 */
	 tx282 = power2(ck); 		/* rule 29 */
	 tx283 = power2(cl); 		/* rule 30 */
	 tx284 = power2(dk); 		/* rule 31 */
	 tx285 = power2(dl); 		/* rule 32 */
	 tzz565 = -pxk; 		/* rule 33 */
	 tzz583 = ak*ck; 		/* rule 34 */
	 tx286 = tzz565*tzz583; 		/* rule 35 */
	 tx287 = bk*ck*pxk; 		/* rule 36 */
	 tzz579 = dk*pxk; 		/* rule 37 */
	 tx288 = ak*tzz579; 		/* rule 38 */
	 tx289 = bk*tzz579; 		/* rule 39 */
	 tzz575 = -pxl; 		/* rule 40 */
	 tx290 = al*cl*tzz575; 		/* rule 41 */
	 tx291 = bl*cl*pxl; 		/* rule 42 */
	 tzz577 = dl*pxl; 		/* rule 43 */
	 tx292 = al*tzz577; 		/* rule 44 */
	 tx293 = bl*tzz577; 		/* rule 45 */
	 tzz584 = ak*bk; 		/* rule 46 */
	 tx294 = pyk*tzz584; 		/* rule 47 */
	 tzz582 = ck*pyk; 		/* rule 48 */
	 tx295 = bk*tzz582; 		/* rule 49 */
	 tzz572 = -pyk; 		/* rule 50 */
	 tx296 = ak*dk*tzz572; 		/* rule 51 */
	 tx297 = dk*tzz582; 		/* rule 52 */
	 tx298 = al*bl*pyl; 		/* rule 53 */
	 tzz581 = cl*pyl; 		/* rule 54 */
	 tx299 = bl*tzz581; 		/* rule 55 */
	 tzz574 = -pyl; 		/* rule 56 */
	 tx300 = al*dl*tzz574; 		/* rule 57 */
	 tx301 = dl*tzz581; 		/* rule 58 */
	 tzz564 = -pzk; 		/* rule 59 */
	 tx302 = tzz564*tzz584; 		/* rule 60 */
	 tx303 = pzk*tzz583; 		/* rule 61 */
	 tzz578 = dk*pzk; 		/* rule 62 */
	 tx304 = bk*tzz578; 		/* rule 63 */
	 tx305 = ck*tzz578; 		/* rule 64 */
	 tzz573 = bl*pzl; 		/* rule 65 */
	 tx306 = -(al*tzz573); 		/* rule 66 */
	 tzz580 = cl*pzl; 		/* rule 67 */
	 tx307 = al*tzz580; 		/* rule 68 */
	 tx308 = dl*tzz573; 		/* rule 69 */
	 tx309 = dl*tzz580; 		/* rule 70 */
	 tx310 = tx280*tzz572; 		/* rule 71 */
	 tx311 = tx280*tzz564; 		/* rule 72 */
	 tx312 = tx281*tzz574; 		/* rule 73 */
	 tzz585 = -pzl; 		/* rule 74 */
	 tx313 = tx281*tzz585; 		/* rule 75 */
	 tx314 = tx282*tzz565; 		/* rule 76 */
	 tx315 = tx282*tzz564; 		/* rule 77 */
	 tx316 = tx283*tzz575; 		/* rule 78 */
	 tx317 = tx283*tzz585; 		/* rule 79 */
	 tx318 = tx284*tzz565; 		/* rule 80 */
	 tx319 = tx284*tzz572; 		/* rule 81 */
	 tx320 = tx278 + tx280 + tx282 + tx284; 		/* rule 82 */
	 tx321 = tx285*tzz575; 		/* rule 83 */
	 tx322 = tx285*tzz574; 		/* rule 84 */
	 tx323 = tx279 + tx281 + tx283 + tx285; 		/* rule 85 */
	 tx324 = tx286 + tx289 + tx294 + tx297 + tx311 + tx315; 		/* rule 86 */
	 tzz599 = tx298 + tx301; 		/* rule 87 */
	 tzz600 = tx293 + tzz599; 		/* rule 88 */
	 tzz603 = tx290 + tzz600; 		/* rule 89 */
	 tx325 = tx313 + tx317 + tzz603; 		/* rule 90 */
	 tx326 = tx295 + tx296 + tx303 + tx304 + tx314 + tx318; 		/* rule 91 */
	 tx327 = tx287 + tx288 + tx302 + tx305 + tx310 + tx319; 		/* rule 92 */
	 tx328 = reciprocal(tx320); 		/* rule 93 */
	 tzz595 = tx307 + tx308; 		/* rule 94 */
	 tzz597 = tx300 + tzz595; 		/* rule 95 */
	 tzz598 = tx299 + tzz597; 		/* rule 96 */
	 tx329 = tx316 + tx321 + tzz598; 		/* rule 97 */
	 tzz596 = tx306 + tx309; 		/* rule 98 */
	 tzz601 = tx292 + tzz596; 		/* rule 99 */
	 tzz602 = tx291 + tzz601; 		/* rule 100 */
	 tx330 = tx312 + tx322 + tzz602; 		/* rule 101 */
	 tx331 = reciprocal(tx323); 		/* rule 102 */
	 tx332 = tzz575; 		/* rule 103 */
	 tx333 = tzz574; 		/* rule 104 */
	 tx334 = tzz585; 		/* rule 105 */
	 tzz556 = 2.*tx328; 		/* rule 106 */
	 tx335 = tx324*tzz556; 		/* rule 107 */
	 tx336 = tx326*tzz556; 		/* rule 108 */
	 tx337 = tx327*tzz556; 		/* rule 109 */
	 tzz555 = -2.*tx331; 		/* rule 110 */
	 tx338 = tx325*tzz555; 		/* rule 111 */
	 tx339 = tx329*tzz555; 		/* rule 112 */
	 tx340 = tx330*tzz555; 		/* rule 113 */
	 tx341 = -xl; 		/* rule 114 */
	 tx342 = -yl; 		/* rule 115 */
	 tx343 = -zl; 		/* rule 116 */
	 tzz588 = tx341 + xk; 		/* rule 117 */
	 tzz590 = tx336 + tzz588; 		/* rule 118 */
	 tzz594 = tx332 + tzz590; 		/* rule 119 */
	 tzz606 = pxk + tzz594; 		/* rule 120 */
	 tx344 = tx339 + tzz606; 		/* rule 121 */
	 tzz587 = tx342 + yk; 		/* rule 122 */
	 tzz589 = tx337 + tzz587; 		/* rule 123 */
	 tzz593 = tx333 + tzz589; 		/* rule 124 */
	 tzz605 = pyk + tzz593; 		/* rule 125 */
	 tx345 = tx340 + tzz605; 		/* rule 126 */
	 tzz586 = tx343 + zk; 		/* rule 127 */
	 tzz591 = tx335 + tzz586; 		/* rule 128 */
	 tzz592 = tx334 + tzz591; 		/* rule 129 */
	 tzz604 = pzk + tzz592; 		/* rule 130 */
	 tx346 = tx338 + tzz604; 		/* rule 131 */
	 tx347 = power2(tx344); 		/* rule 132 */
	 tx348 = power2(tx345); 		/* rule 133 */
	 tx349 = power2(tx346); 		/* rule 134 */
	 tx350 = tx347 + tx348 + tx349; 		/* rule 135 */
	 tx351 = -r0; 		/* rule 136 */
	 tx352 = mysqrt(tx350); 		/* rule 137 */
	 tx353 = tx351 + tx352; 		/* rule 138 */
	 tx354 = power2(tx353); 		/* rule 139 */
	 Energy = ks*tx354; 		/* rule 140 */
	STAPLE_ENERGY_ACCUMULATE(Energy);
	 tx355 = tx283*tx332; 		/* rule 142 */
	 tx356 = tx285*tx332; 		/* rule 143 */
	 tx357 = tx281*tx333; 		/* rule 144 */
	 tx358 = tx285*tx333; 		/* rule 145 */
	 tx359 = tx281*tx334; 		/* rule 146 */
	 tx360 = tx283*tx334; 		/* rule 147 */
	 tx361 = tx355 + tx356 + tzz598; 		/* rule 148 */
	 tx362 = tx357 + tx358 + tzz602; 		/* rule 149 */
	 tx363 = tx359 + tx360 + tzz603; 		/* rule 150 */
	 tx364 = ck*tzz565; 		/* rule 151 */
	 tx365 = tzz579; 		/* rule 152 */
	 tx366 = bk*pyk; 		/* rule 153 */
	 tx367 = dk*tzz572; 		/* rule 154 */
	 tx368 = bk*tzz564; 		/* rule 155 */
	 tx369 = ck*pzk; 		/* rule 156 */
	 tx370 = tx361*tzz555; 		/* rule 157 */
	 tx371 = tx362*tzz555; 		/* rule 158 */
	 tx372 = tx363*tzz555; 		/* rule 159 */
	 tx373 = power2(tx328); 		/* rule 160 */
	 tx374 = tx364 + tx366; 		/* rule 161 */
	 tx375 = tx365 + tx368; 		/* rule 162 */
	 tx376 = tx367 + tx369; 		/* rule 163 */
	 tx377 = tx370 + tzz606; 		/* rule 164 */
	 tx378 = tx371 + tzz605; 		/* rule 165 */
	 tx379 = tx372 + tzz604; 		/* rule 166 */
	 tzz560 = -4.*tx373; 		/* rule 167 */
	 tzz571 = tx324*tzz560; 		/* rule 168 */
	 tx380 = ak*tzz571; 		/* rule 169 */
	 tzz570 = tx326*tzz560; 		/* rule 170 */
	 tx381 = ak*tzz570; 		/* rule 171 */
	 tzz569 = tx327*tzz560; 		/* rule 172 */
	 tx382 = ak*tzz569; 		/* rule 173 */
	 tx383 = tx374*tzz556; 		/* rule 174 */
	 tx384 = tx375*tzz556; 		/* rule 175 */
	 tx385 = tx376*tzz556; 		/* rule 176 */
	 tx386 = power2(tx377); 		/* rule 177 */
	 tx387 = power2(tx378); 		/* rule 178 */
	 tx388 = power2(tx379); 		/* rule 179 */
	 tx389 = tx380 + tx383; 		/* rule 180 */
	 tx390 = tx382 + tx384; 		/* rule 181 */
	 tx391 = tx381 + tx385; 		/* rule 182 */
	 tx392 = tx386 + tx387 + tx388; 		/* rule 183 */
	 tzz561 = 2.*tx379; 		/* rule 184 */
	 tx393 = tx389*tzz561; 		/* rule 185 */
	 tzz562 = 2.*tx378; 		/* rule 186 */
	 tx394 = tx390*tzz562; 		/* rule 187 */
	 tzz563 = 2.*tx377; 		/* rule 188 */
	 tx395 = tx391*tzz563; 		/* rule 189 */
	 tx396 = mysqrt(tx392); 		/* rule 190 */
	 tx397 = reciprocal(tx396); 		/* rule 191 */
	 tx398 = tx393 + tx394 + tx395; 		/* rule 192 */
	 tx399 = tx351 + tx396; 		/* rule 193 */
	 tzz557 = tx397*tx399; 		/* rule 194 */
	 tzz558 = ks*tzz557; 		/* rule 195 */
	 gak = tx398*tzz558; 		/* rule 196 */
	 fak = -gak; 		/* rule 197 */
	STAPLE_FORCE_ACCUMULATE(I1, 0, fak );
	 tx400 = ck*pxk; 		/* rule 199 */
	 tx401 = ak*pyk; 		/* rule 200 */
	 tx402 = tzz582; 		/* rule 201 */
	 tx403 = ak*tzz564; 		/* rule 202 */
	 tx404 = -2.*bk*pzk; 		/* rule 203 */
	 tx405 = tzz578; 		/* rule 204 */
	 tx406 = -2.*tx366; 		/* rule 205 */
	 tx407 = tx365 + tx401 + tx404; 		/* rule 206 */
	 tx408 = tx402 + tx405; 		/* rule 207 */
	 tx409 = tx400 + tx403 + tx406; 		/* rule 208 */
	 tx410 = bk*tzz571; 		/* rule 209 */
	 tx411 = bk*tzz570; 		/* rule 210 */
	 tx412 = bk*tzz569; 		/* rule 211 */
	 tx413 = tx407*tzz556; 		/* rule 212 */
	 tx414 = tx408*tzz556; 		/* rule 213 */
	 tx415 = tx409*tzz556; 		/* rule 214 */
	 tx416 = tx410 + tx413; 		/* rule 215 */
	 tx417 = tx411 + tx414; 		/* rule 216 */
	 tx418 = tx412 + tx415; 		/* rule 217 */
	 tx419 = tx416*tzz561; 		/* rule 218 */
	 tx420 = tx417*tzz563; 		/* rule 219 */
	 tx421 = tx418*tzz562; 		/* rule 220 */
	 tx422 = tx419 + tx420 + tx421; 		/* rule 221 */
	 gbk = tx422*tzz558; 		/* rule 222 */
	 fbk = -gbk; 		/* rule 223 */
	STAPLE_FORCE_ACCUMULATE(I1, 1, fbk );
	 tx423 = ak*tzz565; 		/* rule 225 */
	 tx424 = bk*pxk; 		/* rule 226 */
	 tx425 = dk*pyk; 		/* rule 227 */
	 tx426 = ak*pzk; 		/* rule 228 */
	 tx427 = -2.*tx369; 		/* rule 229 */
	 tx428 = -2.*tx400; 		/* rule 230 */
	 tx429 = tx405 + tx424; 		/* rule 231 */
	 tx430 = tx423 + tx425 + tx427; 		/* rule 232 */
	 tx431 = tx366 + tx426 + tx428; 		/* rule 233 */
	 tx432 = ck*tzz571; 		/* rule 234 */
	 tx433 = ck*tzz570; 		/* rule 235 */
	 tx434 = ck*tzz569; 		/* rule 236 */
	 tx435 = tx429*tzz556; 		/* rule 237 */
	 tx436 = tx430*tzz556; 		/* rule 238 */
	 tx437 = tx431*tzz556; 		/* rule 239 */
	 tx438 = tx434 + tx435; 		/* rule 240 */
	 tx439 = tx432 + tx436; 		/* rule 241 */
	 tx440 = tx433 + tx437; 		/* rule 242 */
	 tx441 = tx438*tzz562; 		/* rule 243 */
	 tx442 = tx439*tzz561; 		/* rule 244 */
	 tx443 = tx440*tzz563; 		/* rule 245 */
	 tx444 = tx441 + tx442 + tx443; 		/* rule 246 */
	 gck = tx444*tzz558; 		/* rule 247 */
	 fck = -gck; 		/* rule 248 */
	STAPLE_FORCE_ACCUMULATE(I1, 2, fck );
	 tx445 = ak*pxk; 		/* rule 250 */
	 tx446 = bk*pzk; 		/* rule 251 */
	 tx447 = -2.*tx365; 		/* rule 252 */
	 tx448 = -tx401; 		/* rule 253 */
	 tx449 = -2.*tx425; 		/* rule 254 */
	 tx450 = tx402 + tx424; 		/* rule 255 */
	 tx451 = tx446 + tx447 + tx448; 		/* rule 256 */
	 tx452 = tx369 + tx445 + tx449; 		/* rule 257 */
	 tx453 = dk*tzz571; 		/* rule 258 */
	 tx454 = dk*tzz570; 		/* rule 259 */
	 tx455 = dk*tzz569; 		/* rule 260 */
	 tx456 = tx450*tzz556; 		/* rule 261 */
	 tx457 = tx451*tzz556; 		/* rule 262 */
	 tx458 = tx452*tzz556; 		/* rule 263 */
	 tx459 = tx453 + tx456; 		/* rule 264 */
	 tx460 = tx454 + tx457; 		/* rule 265 */
	 tx461 = tx455 + tx458; 		/* rule 266 */
	 tx462 = tx459*tzz561; 		/* rule 267 */
	 tx463 = tx460*tzz563; 		/* rule 268 */
	 tx464 = tx461*tzz562; 		/* rule 269 */
	 tx465 = tx462 + tx463 + tx464; 		/* rule 270 */
	 gdk = tx465*tzz558; 		/* rule 271 */
	 fdk = -gdk; 		/* rule 272 */
	STAPLE_FORCE_ACCUMULATE(I1, 3, fdk );
	 gxk = tzz558*tzz563; 		/* rule 274 */
	 fxk = -gxk; 		/* rule 275 */
	STAPLE_FORCE_ACCUMULATE(I1, 4, fxk );
	 gyk = tzz558*tzz562; 		/* rule 277 */
	 fyk = -gyk; 		/* rule 278 */
	STAPLE_FORCE_ACCUMULATE(I1, 5, fyk );
	 gzk = tzz558*tzz561; 		/* rule 280 */
	 fzk = -gzk; 		/* rule 281 */
	STAPLE_FORCE_ACCUMULATE(I1, 6, fzk );
	 tx466 = tzz577; 		/* rule 283 */
	 tx467 = bl*pyl; 		/* rule 284 */
	 tx468 = tzz580; 		/* rule 285 */
	 tx469 = cl*tx332; 		/* rule 286 */
	 tx470 = dl*tx333; 		/* rule 287 */
	 tx471 = bl*tx334; 		/* rule 288 */
	 tx472 = power2(tx331); 		/* rule 289 */
	 tx473 = tx467 + tx469; 		/* rule 290 */
	 tx474 = tx468 + tx470; 		/* rule 291 */
	 tx475 = tx466 + tx471; 		/* rule 292 */
	 tzz559 = 4.*tx472; 		/* rule 293 */
	 tzz568 = tx361*tzz559; 		/* rule 294 */
	 tx476 = al*tzz568; 		/* rule 295 */
	 tzz567 = tx362*tzz559; 		/* rule 296 */
	 tx477 = al*tzz567; 		/* rule 297 */
	 tzz566 = tx363*tzz559; 		/* rule 298 */
	 tx478 = al*tzz566; 		/* rule 299 */
	 tx479 = tx473*tzz555; 		/* rule 300 */
	 tx480 = tx474*tzz555; 		/* rule 301 */
	 tx481 = tx475*tzz555; 		/* rule 302 */
	 tx482 = tx478 + tx479; 		/* rule 303 */
	 tx483 = tx476 + tx480; 		/* rule 304 */
	 tx484 = tx477 + tx481; 		/* rule 305 */
	 tx485 = tx482*tzz561; 		/* rule 306 */
	 tx486 = tx483*tzz563; 		/* rule 307 */
	 tx487 = tx484*tzz562; 		/* rule 308 */
	 tx488 = tx485 + tx486 + tx487; 		/* rule 309 */
	 gal = tx488*tzz558; 		/* rule 310 */
	 fal = -gal; 		/* rule 311 */
	STAPLE_FORCE_ACCUMULATE(I2, 0, fal );
	 tx489 = cl*pxl; 		/* rule 313 */
	 tx490 = al*pyl; 		/* rule 314 */
	 tx491 = tzz581; 		/* rule 315 */
	 tx492 = -2.*tzz573; 		/* rule 316 */
	 tx493 = dl*pzl; 		/* rule 317 */
	 tx494 = al*tx334; 		/* rule 318 */
	 tx495 = -2.*tx467; 		/* rule 319 */
	 tx496 = tx466 + tx490 + tx492; 		/* rule 320 */
	 tx497 = tx491 + tx493; 		/* rule 321 */
	 tx498 = tx489 + tx494 + tx495; 		/* rule 322 */
	 tx499 = bl*tzz568; 		/* rule 323 */
	 tx500 = bl*tzz567; 		/* rule 324 */
	 tx501 = bl*tzz566; 		/* rule 325 */
	 tx502 = tx496*tzz555; 		/* rule 326 */
	 tx503 = tx497*tzz555; 		/* rule 327 */
	 tx504 = tx498*tzz555; 		/* rule 328 */
	 tx505 = tx501 + tx502; 		/* rule 329 */
	 tx506 = tx499 + tx503; 		/* rule 330 */
	 tx507 = tx500 + tx504; 		/* rule 331 */
	 tx508 = tx505*tzz561; 		/* rule 332 */
	 tx509 = tx506*tzz563; 		/* rule 333 */
	 tx510 = tx507*tzz562; 		/* rule 334 */
	 tx511 = tx508 + tx509 + tx510; 		/* rule 335 */
	 gbl = tx511*tzz558; 		/* rule 336 */
	 fbl = -gbl; 		/* rule 337 */
	STAPLE_FORCE_ACCUMULATE(I2, 1, fbl );
	 tx512 = bl*pxl; 		/* rule 339 */
	 tx513 = dl*pyl; 		/* rule 340 */
	 tx514 = al*pzl; 		/* rule 341 */
	 tx515 = al*tx332; 		/* rule 342 */
	 tx516 = -2.*tx468; 		/* rule 343 */
	 tx517 = -2.*tx489; 		/* rule 344 */
	 tx518 = tx493 + tx512; 		/* rule 345 */
	 tx519 = tx513 + tx515 + tx516; 		/* rule 346 */
	 tx520 = tx467 + tx514 + tx517; 		/* rule 347 */
	 tx521 = cl*tzz568; 		/* rule 348 */
	 tx522 = cl*tzz567; 		/* rule 349 */
	 tx523 = cl*tzz566; 		/* rule 350 */
	 tx524 = tx518*tzz555; 		/* rule 351 */
	 tx525 = tx519*tzz555; 		/* rule 352 */
	 tx526 = tx520*tzz555; 		/* rule 353 */
	 tx527 = tx522 + tx524; 		/* rule 354 */
	 tx528 = tx523 + tx525; 		/* rule 355 */
	 tx529 = tx521 + tx526; 		/* rule 356 */
	 tx530 = tx527*tzz562; 		/* rule 357 */
	 tx531 = tx528*tzz561; 		/* rule 358 */
	 tx532 = tx529*tzz563; 		/* rule 359 */
	 tx533 = tx530 + tx531 + tx532; 		/* rule 360 */
	 gcl = tx533*tzz558; 		/* rule 361 */
	 fcl = -gcl; 		/* rule 362 */
	STAPLE_FORCE_ACCUMULATE(I2, 2, fcl );
	 tx534 = al*pxl; 		/* rule 364 */
	 tx535 = tzz573; 		/* rule 365 */
	 tx536 = al*tx333; 		/* rule 366 */
	 tx537 = -2.*tx466; 		/* rule 367 */
	 tx538 = -2.*tx513; 		/* rule 368 */
	 tx539 = tx491 + tx512; 		/* rule 369 */
	 tx540 = tx535 + tx536 + tx537; 		/* rule 370 */
	 tx541 = tx468 + tx534 + tx538; 		/* rule 371 */
	 tx542 = dl*tzz568; 		/* rule 372 */
	 tx543 = dl*tzz567; 		/* rule 373 */
	 tx544 = dl*tzz566; 		/* rule 374 */
	 tx545 = tx539*tzz555; 		/* rule 375 */
	 tx546 = tx540*tzz555; 		/* rule 376 */
	 tx547 = tx541*tzz555; 		/* rule 377 */
	 tx548 = tx544 + tx545; 		/* rule 378 */
	 tx549 = tx542 + tx546; 		/* rule 379 */
	 tx550 = tx543 + tx547; 		/* rule 380 */
	 tx551 = tx548*tzz561; 		/* rule 381 */
	 tx552 = tx549*tzz563; 		/* rule 382 */
	 tx553 = tx550*tzz562; 		/* rule 383 */
	 tx554 = tx551 + tx552 + tx553; 		/* rule 384 */
	 gdl = tx554*tzz558; 		/* rule 385 */
	 fdl = -gdl; 		/* rule 386 */
	STAPLE_FORCE_ACCUMULATE(I2, 3, fdl );
	 tzz576 = -2.*tzz558; 		/* rule 388 */
	 gxl = tx377*tzz576; 		/* rule 389 */
	 fxl = -gxl; 		/* rule 390 */
	STAPLE_FORCE_ACCUMULATE(I2, 4, fxl );
	 gyl = tx378*tzz576; 		/* rule 392 */
	 fyl = -gyl; 		/* rule 393 */
	STAPLE_FORCE_ACCUMULATE(I2, 5, fyl );
	 gzl = tx379*tzz576; 		/* rule 395 */
	 fzl = -gzl; 		/* rule 396 */
	STAPLE_FORCE_ACCUMULATE(I2, 6, fzl );
