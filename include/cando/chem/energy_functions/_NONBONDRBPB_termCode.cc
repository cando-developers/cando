// C-code
	NONBONDRB_SET_PARAMETER(I1);
	NONBONDRB_SET_PARAMETER(I2);
	NONBONDRB_SET_PARAMETER(dQ1Q2);
	NONBONDRB_SET_PARAMETER(dA);
	NONBONDRB_SET_PARAMETER(dC);
	NONBONDRB_SET_POSITION(am,I1,0);
	NONBONDRB_SET_POSITION(bm,I1,1);
	NONBONDRB_SET_POSITION(cm,I1,2);
	NONBONDRB_SET_POSITION(dm,I1,3);
	NONBONDRB_SET_POSITION(xm,I1,4);
	NONBONDRB_SET_POSITION(ym,I1,5);
	NONBONDRB_SET_POSITION(zm,I1,6);
	NONBONDRB_SET_POSITION(an,I2,0);
	NONBONDRB_SET_POSITION(bn,I2,1);
	NONBONDRB_SET_POSITION(cn,I2,2);
	NONBONDRB_SET_POSITION(dn,I2,3);
	NONBONDRB_SET_POSITION(xn,I2,4);
	NONBONDRB_SET_POSITION(yn,I2,5);
	NONBONDRB_SET_POSITION(zn,I2,6);
	NONBONDRB_SET_POINT(pxm,ea1,getX());
	NONBONDRB_SET_POINT(pym,ea1,getY());
	NONBONDRB_SET_POINT(pzm,ea1,getZ());
	NONBONDRB_SET_POINT(pxn,ea2,getX());
	NONBONDRB_SET_POINT(pyn,ea2,getY());
	NONBONDRB_SET_POINT(pzn,ea2,getZ());
	 tx318 = power2(am); 		/* rule 26 */
	 tx319 = power2(an); 		/* rule 27 */
	 tx320 = power2(bm); 		/* rule 28 */
	 tx321 = power2(bn); 		/* rule 29 */
	 tx322 = power2(cm); 		/* rule 30 */
	 tx323 = power2(cn); 		/* rule 31 */
	 tx324 = power2(dm); 		/* rule 32 */
	 tx325 = power2(dn); 		/* rule 33 */
	 tzz679 = cm*pym; 		/* rule 34 */
	 tx326 = bm*tzz679; 		/* rule 35 */
	 tzz666 = -pym; 		/* rule 36 */
	 tx327 = am*dm*tzz666; 		/* rule 37 */
	 tzz678 = cn*pyn; 		/* rule 38 */
	 tx328 = bn*tzz678; 		/* rule 39 */
	 tzz671 = -pyn; 		/* rule 40 */
	 tx329 = an*dn*tzz671; 		/* rule 41 */
	 tzz680 = am*cm; 		/* rule 42 */
	 tx330 = pzm*tzz680; 		/* rule 43 */
	 tzz675 = dm*pzm; 		/* rule 44 */
	 tx331 = bm*tzz675; 		/* rule 45 */
	 tzz677 = cn*pzn; 		/* rule 46 */
	 tx332 = an*tzz677; 		/* rule 47 */
	 tzz667 = bn*pzn; 		/* rule 48 */
	 tx333 = dn*tzz667; 		/* rule 49 */
	 tzz659 = -pxm; 		/* rule 50 */
	 tx334 = tx322*tzz659; 		/* rule 51 */
	 tzz682 = -tx323; 		/* rule 52 */
	 tx335 = pxn*tzz682; 		/* rule 53 */
	 tx336 = tx324*tzz659; 		/* rule 54 */
	 tx337 = tx318 + tx320 + tx322 + tx324; 		/* rule 55 */
	 tx338 = -(pxn*tx325); 		/* rule 56 */
	 tx339 = tx319 + tx321 + tx323 + tx325; 		/* rule 57 */
	 tx340 = tx326 + tx327 + tx330 + tx331 + tx334 + tx336; 		/* rule 58 */
	 tx341 = reciprocal(tx337); 		/* rule 59 */
	 tzz699 = tx332 + tx333; 		/* rule 60 */
	 tzz700 = tx329 + tzz699; 		/* rule 61 */
	 tzz701 = tx328 + tzz700; 		/* rule 62 */
	 tx342 = tx335 + tx338 + tzz701; 		/* rule 63 */
	 tx343 = reciprocal(tx339); 		/* rule 64 */
	 tx344 = -pxn; 		/* rule 65 */
	 tzz646 = 2.*tx341; 		/* rule 66 */
	 tx345 = tx340*tzz646; 		/* rule 67 */
	 tzz645 = -2.*tx343; 		/* rule 68 */
	 tx346 = tx342*tzz645; 		/* rule 69 */
	 tx347 = -xn; 		/* rule 70 */
	 tzz696 = tx347 + xm; 		/* rule 71 */
	 tzz697 = tx345 + tzz696; 		/* rule 72 */
	 tzz698 = tx344 + tzz697; 		/* rule 73 */
	 tzz704 = pxm + tzz698; 		/* rule 74 */
	 DeltaX = tx346 + tzz704; 		/* rule 75 */
	 tx348 = bm*cm*pxm; 		/* rule 76 */
	 tzz676 = dm*pxm; 		/* rule 77 */
	 tx349 = am*tzz676; 		/* rule 78 */
	 tx350 = bn*cn*pxn; 		/* rule 79 */
	 tzz674 = dn*pxn; 		/* rule 80 */
	 tx351 = an*tzz674; 		/* rule 81 */
	 tzz658 = -pzm; 		/* rule 82 */
	 tzz681 = am*bm; 		/* rule 83 */
	 tx352 = tzz658*tzz681; 		/* rule 84 */
	 tx353 = cm*tzz675; 		/* rule 85 */
	 tx354 = -(an*tzz667); 		/* rule 86 */
	 tx355 = dn*tzz677; 		/* rule 87 */
	 tx356 = tx320*tzz666; 		/* rule 88 */
	 tx357 = tx321*tzz671; 		/* rule 89 */
	 tx358 = tx324*tzz666; 		/* rule 90 */
	 tx359 = tx325*tzz671; 		/* rule 91 */
	 tx360 = tx348 + tx349 + tx352 + tx353 + tx356 + tx358; 		/* rule 92 */
	 tzz693 = tx354 + tx355; 		/* rule 93 */
	 tzz694 = tx351 + tzz693; 		/* rule 94 */
	 tzz695 = tx350 + tzz694; 		/* rule 95 */
	 tx361 = tx357 + tx359 + tzz695; 		/* rule 96 */
	 tx362 = tzz671; 		/* rule 97 */
	 tx363 = tx360*tzz646; 		/* rule 98 */
	 tx364 = tx361*tzz645; 		/* rule 99 */
	 tx365 = -yn; 		/* rule 100 */
	 tzz690 = tx365 + ym; 		/* rule 101 */
	 tzz691 = tx363 + tzz690; 		/* rule 102 */
	 tzz692 = tx362 + tzz691; 		/* rule 103 */
	 tzz703 = pym + tzz692; 		/* rule 104 */
	 DeltaY = tx364 + tzz703; 		/* rule 105 */
	 tx366 = tzz659*tzz680; 		/* rule 106 */
	 tx367 = bm*tzz676; 		/* rule 107 */
	 tx368 = bn*tzz674; 		/* rule 108 */
	 tx369 = pym*tzz681; 		/* rule 109 */
	 tx370 = dm*tzz679; 		/* rule 110 */
	 tx371 = an*bn*pyn; 		/* rule 111 */
	 tx372 = dn*tzz678; 		/* rule 112 */
	 tx373 = tx320*tzz658; 		/* rule 113 */
	 tx374 = -(pzn*tx321); 		/* rule 114 */
	 tx375 = tx322*tzz658; 		/* rule 115 */
	 tx376 = pzn*tzz682; 		/* rule 116 */
	 tx377 = an*cn*tx344; 		/* rule 117 */
	 tx378 = tx366 + tx367 + tx369 + tx370 + tx373 + tx375; 		/* rule 118 */
	 tzz687 = tx372 + tx377; 		/* rule 119 */
	 tzz688 = tx371 + tzz687; 		/* rule 120 */
	 tzz689 = tx368 + tzz688; 		/* rule 121 */
	 tx379 = tx374 + tx376 + tzz689; 		/* rule 122 */
	 tx380 = -pzn; 		/* rule 123 */
	 tx381 = tx378*tzz646; 		/* rule 124 */
	 tx382 = tx379*tzz645; 		/* rule 125 */
	 tx383 = -zn; 		/* rule 126 */
	 tzz684 = tx383 + zm; 		/* rule 127 */
	 tzz685 = tx381 + tzz684; 		/* rule 128 */
	 tzz686 = tx380 + tzz685; 		/* rule 129 */
	 tzz702 = pzm + tzz686; 		/* rule 130 */
	 DeltaZ = tx382 + tzz702; 		/* rule 131 */
	 tx384 = tx323*tx344; 		/* rule 132 */
	 tx385 = tx325*tx344; 		/* rule 133 */
	 tx386 = tx321*tx362; 		/* rule 134 */
	 tx387 = tx325*tx362; 		/* rule 135 */
	 tx388 = tx321*tx380; 		/* rule 136 */
	 tx389 = tx323*tx380; 		/* rule 137 */
	 tx390 = tx384 + tx385 + tzz701; 		/* rule 138 */
	 tx391 = tx386 + tx387 + tzz695; 		/* rule 139 */
	 tx392 = tx388 + tx389 + tzz689; 		/* rule 140 */
	 tx393 = tx390*tzz645; 		/* rule 141 */
	 tx394 = tx391*tzz645; 		/* rule 142 */
	 tx395 = tx392*tzz645; 		/* rule 143 */
	 tx396 = tx393 + tzz704; 		/* rule 144 */
	 tx397 = tx394 + tzz703; 		/* rule 145 */
	 tx398 = tx395 + tzz702; 		/* rule 146 */
	 tx399 = PBX(tx396); 		/* rule 147 */
	 tx400 = PBY(tx397); 		/* rule 148 */
	 tx401 = PBZ(tx398); 		/* rule 149 */
	 tx402 = power2(tx399); 		/* rule 150 */
	 tx403 = power2(tx400); 		/* rule 151 */
	 tx404 = power2(tx401); 		/* rule 152 */
	 tx405 = tx402 + tx403 + tx404; 		/* rule 153 */
	 tx683 = reciprocal(tx405); 		/* rule 154 */
	 tx643 = power2(tx683); 		/* rule 155 */
	 tx644 = tx683; 		/* rule 156 */
	 tx635 = tx643*tx644; 		/* rule 157 */
	 tx406 = power2(tx635); 		/* rule 158 */
	 tx636 = power2(tx644); 		/* rule 159 */
	 tx637 = tx644; 		/* rule 160 */
	 tx407 = tx636*tx637; 		/* rule 161 */
	 tx638 = mysqrt(tx405); 		/* rule 162 */
	 tx408 = reciprocal(tx638); 		/* rule 163 */
	 tx409 = dA*tx406; 		/* rule 164 */
	 tx410 = -(dC*tx407); 		/* rule 165 */
	 tx411 = dQ1Q2*tx408; 		/* rule 166 */
	 Energy = tx409 + tx410 + tx411; 		/* rule 167 */
	NONBONDRB_ENERGY_ACCUMULATE(Energy);
	#ifdef NONBONDRB_CALC_FORCE //[
	if ( calcForce ) {
	BAIL_OUT_IF_CUTOFF(DeltaX,DeltaY,DeltaZ)
	 tx412 = cm*tzz659; 		/* rule 172 */
	 tx413 = tzz676; 		/* rule 173 */
	 tx414 = bm*pym; 		/* rule 174 */
	 tx415 = dm*tzz666; 		/* rule 175 */
	 tx416 = bm*tzz658; 		/* rule 176 */
	 tx417 = cm*pzm; 		/* rule 177 */
	 tx418 = power2(tx341); 		/* rule 178 */
	 tx419 = tx412 + tx414; 		/* rule 179 */
	 tx420 = tx413 + tx416; 		/* rule 180 */
	 tx421 = tx415 + tx417; 		/* rule 181 */
	 tzz651 = -4.*tx418; 		/* rule 182 */
	 tzz665 = tx340*tzz651; 		/* rule 183 */
	 tx422 = am*tzz665; 		/* rule 184 */
	 tzz664 = tx360*tzz651; 		/* rule 185 */
	 tx423 = am*tzz664; 		/* rule 186 */
	 tzz663 = tx378*tzz651; 		/* rule 187 */
	 tx424 = am*tzz663; 		/* rule 188 */
	 tx425 = tx419*tzz646; 		/* rule 189 */
	 tx426 = tx420*tzz646; 		/* rule 190 */
	 tx427 = tx421*tzz646; 		/* rule 191 */
	 tx428 = tx424 + tx425; 		/* rule 192 */
	 tx429 = tx423 + tx426; 		/* rule 193 */
	 tx430 = tx422 + tx427; 		/* rule 194 */
	 tzz653 = 2.*tx401; 		/* rule 195 */
	 tx431 = tx428*tzz653; 		/* rule 196 */
	 tzz654 = 2.*tx400; 		/* rule 197 */
	 tx432 = tx429*tzz654; 		/* rule 198 */
	 tzz655 = 2.*tx399; 		/* rule 199 */
	 tx433 = tx430*tzz655; 		/* rule 200 */
	 tx639 = power2(tx643); 		/* rule 201 */
	 tx434 = tx635*tx639; 		/* rule 202 */
	 tx435 = power2(tx636); 		/* rule 203 */
	 tx640 = tx405; 		/* rule 204 */
	 tx641 = tx408; 		/* rule 205 */
	 tx642 = reciprocal(tx640); 		/* rule 206 */
	 tx436 = tx641*tx642; 		/* rule 207 */
	 tx437 = tx431 + tx432 + tx433; 		/* rule 208 */
	 tzz649 = dA*tx434; 		/* rule 209 */
	 tzz657 = -6.*tzz649; 		/* rule 210 */
	 tx438 = tx437*tzz657; 		/* rule 211 */
	 tzz648 = dC*tx435; 		/* rule 212 */
	 tzz652 = 3*tzz648; 		/* rule 213 */
	 tx439 = tx437*tzz652; 		/* rule 214 */
	 tzz647 = dQ1Q2*tx436; 		/* rule 215 */
	 tzz656 = -0.5*tzz647; 		/* rule 216 */
	 tx440 = tx437*tzz656; 		/* rule 217 */
	 gam = tx438 + tx439 + tx440; 		/* rule 218 */
	 fam = -gam; 		/* rule 219 */
	NONBONDRB_FORCE_ACCUMULATE(I1, 0, fam );
	 tx441 = cm*pxm; 		/* rule 221 */
	 tx442 = am*pym; 		/* rule 222 */
	 tx443 = tzz679; 		/* rule 223 */
	 tx444 = am*tzz658; 		/* rule 224 */
	 tx445 = -2.*bm*pzm; 		/* rule 225 */
	 tx446 = tzz675; 		/* rule 226 */
	 tx447 = -2.*tx414; 		/* rule 227 */
	 tx448 = tx413 + tx442 + tx445; 		/* rule 228 */
	 tx449 = tx443 + tx446; 		/* rule 229 */
	 tx450 = tx441 + tx444 + tx447; 		/* rule 230 */
	 tx451 = bm*tzz665; 		/* rule 231 */
	 tx452 = bm*tzz664; 		/* rule 232 */
	 tx453 = bm*tzz663; 		/* rule 233 */
	 tx454 = tx448*tzz646; 		/* rule 234 */
	 tx455 = tx449*tzz646; 		/* rule 235 */
	 tx456 = tx450*tzz646; 		/* rule 236 */
	 tx457 = tx453 + tx454; 		/* rule 237 */
	 tx458 = tx451 + tx455; 		/* rule 238 */
	 tx459 = tx452 + tx456; 		/* rule 239 */
	 tx460 = tx457*tzz653; 		/* rule 240 */
	 tx461 = tx458*tzz655; 		/* rule 241 */
	 tx462 = tx459*tzz654; 		/* rule 242 */
	 tx463 = tx460 + tx461 + tx462; 		/* rule 243 */
	 tx464 = tx463*tzz657; 		/* rule 244 */
	 tx465 = tx463*tzz652; 		/* rule 245 */
	 tx466 = tx463*tzz656; 		/* rule 246 */
	 gbm = tx464 + tx465 + tx466; 		/* rule 247 */
	 fbm = -gbm; 		/* rule 248 */
	NONBONDRB_FORCE_ACCUMULATE(I1, 1, fbm );
	 tx467 = am*tzz659; 		/* rule 250 */
	 tx468 = bm*pxm; 		/* rule 251 */
	 tx469 = dm*pym; 		/* rule 252 */
	 tx470 = am*pzm; 		/* rule 253 */
	 tx471 = -2.*tx417; 		/* rule 254 */
	 tx472 = -2.*tx441; 		/* rule 255 */
	 tx473 = tx446 + tx468; 		/* rule 256 */
	 tx474 = tx467 + tx469 + tx471; 		/* rule 257 */
	 tx475 = tx414 + tx470 + tx472; 		/* rule 258 */
	 tx476 = cm*tzz665; 		/* rule 259 */
	 tx477 = cm*tzz664; 		/* rule 260 */
	 tx478 = cm*tzz663; 		/* rule 261 */
	 tx479 = tx473*tzz646; 		/* rule 262 */
	 tx480 = tx474*tzz646; 		/* rule 263 */
	 tx481 = tx475*tzz646; 		/* rule 264 */
	 tx482 = tx477 + tx479; 		/* rule 265 */
	 tx483 = tx478 + tx480; 		/* rule 266 */
	 tx484 = tx476 + tx481; 		/* rule 267 */
	 tx485 = tx482*tzz654; 		/* rule 268 */
	 tx486 = tx483*tzz653; 		/* rule 269 */
	 tx487 = tx484*tzz655; 		/* rule 270 */
	 tx488 = tx485 + tx486 + tx487; 		/* rule 271 */
	 tx489 = tx488*tzz657; 		/* rule 272 */
	 tx490 = tx488*tzz652; 		/* rule 273 */
	 tx491 = tx488*tzz656; 		/* rule 274 */
	 gcm = tx489 + tx490 + tx491; 		/* rule 275 */
	 fcm = -gcm; 		/* rule 276 */
	NONBONDRB_FORCE_ACCUMULATE(I1, 2, fcm );
	 tx492 = am*pxm; 		/* rule 278 */
	 tx493 = bm*pzm; 		/* rule 279 */
	 tx494 = -2.*tx413; 		/* rule 280 */
	 tx495 = -tx442; 		/* rule 281 */
	 tx496 = -2.*tx469; 		/* rule 282 */
	 tx497 = tx443 + tx468; 		/* rule 283 */
	 tx498 = tx493 + tx494 + tx495; 		/* rule 284 */
	 tx499 = tx417 + tx492 + tx496; 		/* rule 285 */
	 tx500 = dm*tzz665; 		/* rule 286 */
	 tx501 = dm*tzz664; 		/* rule 287 */
	 tx502 = dm*tzz663; 		/* rule 288 */
	 tx503 = tx497*tzz646; 		/* rule 289 */
	 tx504 = tx498*tzz646; 		/* rule 290 */
	 tx505 = tx499*tzz646; 		/* rule 291 */
	 tx506 = tx502 + tx503; 		/* rule 292 */
	 tx507 = tx500 + tx504; 		/* rule 293 */
	 tx508 = tx501 + tx505; 		/* rule 294 */
	 tx509 = tx506*tzz653; 		/* rule 295 */
	 tx510 = tx507*tzz655; 		/* rule 296 */
	 tx511 = tx508*tzz654; 		/* rule 297 */
	 tx512 = tx509 + tx510 + tx511; 		/* rule 298 */
	 tx513 = tx512*tzz657; 		/* rule 299 */
	 tx514 = tx512*tzz652; 		/* rule 300 */
	 tx515 = tx512*tzz656; 		/* rule 301 */
	 gdm = tx513 + tx514 + tx515; 		/* rule 302 */
	 fdm = -gdm; 		/* rule 303 */
	NONBONDRB_FORCE_ACCUMULATE(I1, 3, fdm );
	 tzz673 = -12.*tzz649; 		/* rule 305 */
	 tx516 = tx399*tzz673; 		/* rule 306 */
	 tzz669 = 6.*tzz648; 		/* rule 307 */
	 tx517 = tx399*tzz669; 		/* rule 308 */
	 tzz670 = -tzz647; 		/* rule 309 */
	 tx518 = tx399*tzz670; 		/* rule 310 */
	 gxm = tx516 + tx517 + tx518; 		/* rule 311 */
	 fxm = -gxm; 		/* rule 312 */
	NONBONDRB_FORCE_ACCUMULATE(I1, 4, fxm );
	 tx519 = tx400*tzz673; 		/* rule 314 */
	 tx520 = tx400*tzz669; 		/* rule 315 */
	 tx521 = tx400*tzz670; 		/* rule 316 */
	 gym = tx519 + tx520 + tx521; 		/* rule 317 */
	 fym = -gym; 		/* rule 318 */
	NONBONDRB_FORCE_ACCUMULATE(I1, 5, fym );
	 tx522 = tx401*tzz673; 		/* rule 320 */
	 tx523 = tx401*tzz669; 		/* rule 321 */
	 tx524 = tx401*tzz670; 		/* rule 322 */
	 gzm = tx522 + tx523 + tx524; 		/* rule 323 */
	 fzm = -gzm; 		/* rule 324 */
	NONBONDRB_FORCE_ACCUMULATE(I1, 6, fzm );
	 tx525 = tzz674; 		/* rule 326 */
	 tx526 = bn*pyn; 		/* rule 327 */
	 tx527 = tzz677; 		/* rule 328 */
	 tx528 = cn*tx344; 		/* rule 329 */
	 tx529 = dn*tx362; 		/* rule 330 */
	 tx530 = bn*tx380; 		/* rule 331 */
	 tx531 = power2(tx343); 		/* rule 332 */
	 tx532 = tx526 + tx528; 		/* rule 333 */
	 tx533 = tx527 + tx529; 		/* rule 334 */
	 tx534 = tx525 + tx530; 		/* rule 335 */
	 tzz650 = 4.*tx531; 		/* rule 336 */
	 tzz662 = tx390*tzz650; 		/* rule 337 */
	 tx535 = an*tzz662; 		/* rule 338 */
	 tzz661 = tx391*tzz650; 		/* rule 339 */
	 tx536 = an*tzz661; 		/* rule 340 */
	 tzz660 = tx392*tzz650; 		/* rule 341 */
	 tx537 = an*tzz660; 		/* rule 342 */
	 tx538 = tx532*tzz645; 		/* rule 343 */
	 tx539 = tx533*tzz645; 		/* rule 344 */
	 tx540 = tx534*tzz645; 		/* rule 345 */
	 tx541 = tx537 + tx538; 		/* rule 346 */
	 tx542 = tx535 + tx539; 		/* rule 347 */
	 tx543 = tx536 + tx540; 		/* rule 348 */
	 tx544 = tx541*tzz653; 		/* rule 349 */
	 tx545 = tx542*tzz655; 		/* rule 350 */
	 tx546 = tx543*tzz654; 		/* rule 351 */
	 tx547 = tx544 + tx545 + tx546; 		/* rule 352 */
	 tx548 = tx547*tzz657; 		/* rule 353 */
	 tx549 = tx547*tzz652; 		/* rule 354 */
	 tx550 = tx547*tzz656; 		/* rule 355 */
	 gan = tx548 + tx549 + tx550; 		/* rule 356 */
	 fan = -gan; 		/* rule 357 */
	NONBONDRB_FORCE_ACCUMULATE(I2, 0, fan );
	 tx551 = cn*pxn; 		/* rule 359 */
	 tx552 = an*pyn; 		/* rule 360 */
	 tx553 = tzz678; 		/* rule 361 */
	 tx554 = -2.*tzz667; 		/* rule 362 */
	 tx555 = dn*pzn; 		/* rule 363 */
	 tx556 = an*tx380; 		/* rule 364 */
	 tx557 = -2.*tx526; 		/* rule 365 */
	 tx558 = tx525 + tx552 + tx554; 		/* rule 366 */
	 tx559 = tx553 + tx555; 		/* rule 367 */
	 tx560 = tx551 + tx556 + tx557; 		/* rule 368 */
	 tx561 = bn*tzz662; 		/* rule 369 */
	 tx562 = bn*tzz661; 		/* rule 370 */
	 tx563 = bn*tzz660; 		/* rule 371 */
	 tx564 = tx558*tzz645; 		/* rule 372 */
	 tx565 = tx559*tzz645; 		/* rule 373 */
	 tx566 = tx560*tzz645; 		/* rule 374 */
	 tx567 = tx563 + tx564; 		/* rule 375 */
	 tx568 = tx561 + tx565; 		/* rule 376 */
	 tx569 = tx562 + tx566; 		/* rule 377 */
	 tx570 = tx567*tzz653; 		/* rule 378 */
	 tx571 = tx568*tzz655; 		/* rule 379 */
	 tx572 = tx569*tzz654; 		/* rule 380 */
	 tx573 = tx570 + tx571 + tx572; 		/* rule 381 */
	 tx574 = tx573*tzz657; 		/* rule 382 */
	 tx575 = tx573*tzz652; 		/* rule 383 */
	 tx576 = tx573*tzz656; 		/* rule 384 */
	 gbn = tx574 + tx575 + tx576; 		/* rule 385 */
	 fbn = -gbn; 		/* rule 386 */
	NONBONDRB_FORCE_ACCUMULATE(I2, 1, fbn );
	 tx577 = bn*pxn; 		/* rule 388 */
	 tx578 = dn*pyn; 		/* rule 389 */
	 tx579 = an*pzn; 		/* rule 390 */
	 tx580 = an*tx344; 		/* rule 391 */
	 tx581 = -2.*tx527; 		/* rule 392 */
	 tx582 = -2.*tx551; 		/* rule 393 */
	 tx583 = tx555 + tx577; 		/* rule 394 */
	 tx584 = tx578 + tx580 + tx581; 		/* rule 395 */
	 tx585 = tx526 + tx579 + tx582; 		/* rule 396 */
	 tx586 = cn*tzz662; 		/* rule 397 */
	 tx587 = cn*tzz661; 		/* rule 398 */
	 tx588 = cn*tzz660; 		/* rule 399 */
	 tx589 = tx583*tzz645; 		/* rule 400 */
	 tx590 = tx584*tzz645; 		/* rule 401 */
	 tx591 = tx585*tzz645; 		/* rule 402 */
	 tx592 = tx587 + tx589; 		/* rule 403 */
	 tx593 = tx588 + tx590; 		/* rule 404 */
	 tx594 = tx586 + tx591; 		/* rule 405 */
	 tx595 = tx592*tzz654; 		/* rule 406 */
	 tx596 = tx593*tzz653; 		/* rule 407 */
	 tx597 = tx594*tzz655; 		/* rule 408 */
	 tx598 = tx595 + tx596 + tx597; 		/* rule 409 */
	 tx599 = tx598*tzz657; 		/* rule 410 */
	 tx600 = tx598*tzz652; 		/* rule 411 */
	 tx601 = tx598*tzz656; 		/* rule 412 */
	 gcn = tx599 + tx600 + tx601; 		/* rule 413 */
	 fcn = -gcn; 		/* rule 414 */
	NONBONDRB_FORCE_ACCUMULATE(I2, 2, fcn );
	 tx602 = an*pxn; 		/* rule 416 */
	 tx603 = tzz667; 		/* rule 417 */
	 tx604 = an*tx362; 		/* rule 418 */
	 tx605 = -2.*tx525; 		/* rule 419 */
	 tx606 = -2.*tx578; 		/* rule 420 */
	 tx607 = tx553 + tx577; 		/* rule 421 */
	 tx608 = tx603 + tx604 + tx605; 		/* rule 422 */
	 tx609 = tx527 + tx602 + tx606; 		/* rule 423 */
	 tx610 = dn*tzz662; 		/* rule 424 */
	 tx611 = dn*tzz661; 		/* rule 425 */
	 tx612 = dn*tzz660; 		/* rule 426 */
	 tx613 = tx607*tzz645; 		/* rule 427 */
	 tx614 = tx608*tzz645; 		/* rule 428 */
	 tx615 = tx609*tzz645; 		/* rule 429 */
	 tx616 = tx612 + tx613; 		/* rule 430 */
	 tx617 = tx610 + tx614; 		/* rule 431 */
	 tx618 = tx611 + tx615; 		/* rule 432 */
	 tx619 = tx616*tzz653; 		/* rule 433 */
	 tx620 = tx617*tzz655; 		/* rule 434 */
	 tx621 = tx618*tzz654; 		/* rule 435 */
	 tx622 = tx619 + tx620 + tx621; 		/* rule 436 */
	 tx623 = tx622*tzz657; 		/* rule 437 */
	 tx624 = tx622*tzz652; 		/* rule 438 */
	 tx625 = tx622*tzz656; 		/* rule 439 */
	 gdn = tx623 + tx624 + tx625; 		/* rule 440 */
	 fdn = -gdn; 		/* rule 441 */
	NONBONDRB_FORCE_ACCUMULATE(I2, 3, fdn );
	 tzz668 = 12.*tzz649; 		/* rule 443 */
	 tx626 = tx399*tzz668; 		/* rule 444 */
	 tzz672 = -6.*tzz648; 		/* rule 445 */
	 tx627 = tx399*tzz672; 		/* rule 446 */
	 tx628 = tx399*tzz647; 		/* rule 447 */
	 gxn = tx626 + tx627 + tx628; 		/* rule 448 */
	 fxn = -gxn; 		/* rule 449 */
	NONBONDRB_FORCE_ACCUMULATE(I2, 4, fxn );
	 tx629 = tx400*tzz668; 		/* rule 451 */
	 tx630 = tx400*tzz672; 		/* rule 452 */
	 tx631 = tx400*tzz647; 		/* rule 453 */
	 gyn = tx629 + tx630 + tx631; 		/* rule 454 */
	 fyn = -gyn; 		/* rule 455 */
	NONBONDRB_FORCE_ACCUMULATE(I2, 5, fyn );
	 tx632 = tx401*tzz668; 		/* rule 457 */
	 tx633 = tx401*tzz672; 		/* rule 458 */
	 tx634 = tx401*tzz647; 		/* rule 459 */
	 gzn = tx632 + tx633 + tx634; 		/* rule 460 */
	 fzn = -gzn; 		/* rule 461 */
	NONBONDRB_FORCE_ACCUMULATE(I2, 6, fzn );
	} /*calcForce */
	#endif /* NONBONDRB_CALC_FORCE ]*/
