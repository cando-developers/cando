// C-code
	ANGLE_SET_PARAMETER(kt);
	ANGLE_SET_PARAMETER(t0);
	ANGLE_SET_PARAMETER(I1);
	ANGLE_SET_PARAMETER(I2);
	ANGLE_SET_PARAMETER(I3);
	ANGLE_SET_POSITION(x1,I1,0);
	ANGLE_SET_POSITION(y1,I1,1);
	ANGLE_SET_POSITION(z1,I1,2);
	ANGLE_SET_POSITION(x2,I2,0);
	ANGLE_SET_POSITION(y2,I2,1);
	ANGLE_SET_POSITION(z2,I2,2);
	ANGLE_SET_POSITION(x3,I3,0);
	ANGLE_SET_POSITION(y3,I3,1);
	ANGLE_SET_POSITION(z3,I3,2);
	 tx450 = -x2; 		/* rule 15 */
	 tx451 = -y2; 		/* rule 16 */
	 tx452 = -z2; 		/* rule 17 */
	 tx453 = tx450 + x1; 		/* rule 18 */
	 tx454 = tx450 + x3; 		/* rule 19 */
	 tx455 = tx451 + y1; 		/* rule 20 */
	 tx456 = tx451 + y3; 		/* rule 21 */
	 tx457 = tx452 + z1; 		/* rule 22 */
	 tx458 = tx452 + z3; 		/* rule 23 */
	 tx459 = power2(tx453); 		/* rule 24 */
	 tx460 = power2(tx454); 		/* rule 25 */
	 tx461 = power2(tx455); 		/* rule 26 */
	 tx462 = power2(tx456); 		/* rule 27 */
	 tx463 = power2(tx457); 		/* rule 28 */
	 tx464 = power2(tx458); 		/* rule 29 */
	 tx465 = tx453*tx454; 		/* rule 30 */
	 tx466 = tx455*tx456; 		/* rule 31 */
	 tx467 = tx457*tx458; 		/* rule 32 */
	 tx468 = tx459 + tx461 + tx463; 		/* rule 33 */
	 tx469 = tx460 + tx462 + tx464; 		/* rule 34 */
	 tx470 = tx465 + tx466 + tx467; 		/* rule 35 */
	 tx899 = mysqrt(tx468); 		/* rule 36 */
	 tx471 = reciprocal(tx899); 		/* rule 37 */
	 tx900 = mysqrt(tx469); 		/* rule 38 */
	 tx472 = reciprocal(tx900); 		/* rule 39 */
	 DotNormAbNormCb = tx470*tx471*tx472; 		/* rule 40 */
	 tx473 = -t0; 		/* rule 41 */
	 tx474 = ArcCos(DotNormAbNormCb); 		/* rule 42 */
	 AngleDeviation = tx473 + tx474; 		/* rule 43 */
	bool IllegalAngle=false;
	if(fabs(DotNormAbNormCb)>(1.0-VERYSMALL)) IllegalAngle=true;
	 tx475 = power2(AngleDeviation); 		/* rule 46 */
	 Energy = kt*tx475; 		/* rule 47 */
	ANGLE_ENERGY_ACCUMULATE(Energy);
	#ifdef ANGLE_CALC_FORCE //[
	if ( calcForce ) {
	 tx476 = reciprocal(tx468); 		/* rule 51 */
	 tx477 = reciprocal(tx469); 		/* rule 52 */
	 tx478 = power2(tx470); 		/* rule 53 */
	 tx901 = tx468; 		/* rule 54 */
	 tx902 = tx471; 		/* rule 55 */
	 tx903 = reciprocal(tx901); 		/* rule 56 */
	 tx479 = tx902*tx903; 		/* rule 57 */
	 tzz943 = tx476*tx478; 		/* rule 58 */
	 tx480 = -(tx477*tzz943); 		/* rule 59 */
	 tzz967 = tx472*tx902; 		/* rule 60 */
	 tx481 = tx454*tzz967; 		/* rule 61 */
	 tzz922 = tx470*tx479; 		/* rule 62 */
	 tzz966 = tx472*tzz922; 		/* rule 63 */
	 tzz971 = -tzz966; 		/* rule 64 */
	 tx482 = tx453*tzz971; 		/* rule 65 */
	 tx483 = 1 + tx480; 		/* rule 66 */
	 tx484 = tx481 + tx482; 		/* rule 67 */
	 tx904 = mysqrt(tx483); 		/* rule 68 */
	 tx485 = reciprocal(tx904); 		/* rule 69 */
	 tzz915 = AngleDeviation*kt; 		/* rule 70 */
	 tzz916 = -2.*tzz915; 		/* rule 71 */
	 tzz931 = tx485*tzz916; 		/* rule 72 */
	 gx1 = tx484*tzz931; 		/* rule 73 */
	 fx1 = -gx1; 		/* rule 74 */
	ANGLE_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx486 = tx456*tzz967; 		/* rule 76 */
	 tx487 = tx455*tzz971; 		/* rule 77 */
	 tx488 = tx486 + tx487; 		/* rule 78 */
	 gy1 = tx488*tzz931; 		/* rule 79 */
	 fy1 = -gy1; 		/* rule 80 */
	ANGLE_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx489 = tx458*tzz967; 		/* rule 82 */
	 tx490 = tx457*tzz971; 		/* rule 83 */
	 tx491 = tx489 + tx490; 		/* rule 84 */
	 gz1 = tx491*tzz931; 		/* rule 85 */
	 fz1 = -gz1; 		/* rule 86 */
	ANGLE_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tx492 = -x1; 		/* rule 88 */
	 tx493 = 2.*x2; 		/* rule 89 */
	 tx494 = -x3; 		/* rule 90 */
	 tx905 = tx469; 		/* rule 91 */
	 tx906 = tx472; 		/* rule 92 */
	 tx907 = reciprocal(tx905); 		/* rule 93 */
	 tx495 = tx906*tx907; 		/* rule 94 */
	 tx496 = tx492 + tx493 + tx494; 		/* rule 95 */
	 tzz917 = tx479*tx906; 		/* rule 96 */
	 tzz953 = tx453*tzz917; 		/* rule 97 */
	 tx497 = tx470*tzz953; 		/* rule 98 */
	 tzz921 = tx495*tx902; 		/* rule 99 */
	 tzz946 = tx454*tzz921; 		/* rule 100 */
	 tx498 = tx470*tzz946; 		/* rule 101 */
	 tzz936 = tx902*tx906; 		/* rule 102 */
	 tx499 = tx496*tzz936; 		/* rule 103 */
	 tx500 = tx497 + tx498 + tx499; 		/* rule 104 */
	 gx2 = tx500*tzz931; 		/* rule 105 */
	 fx2 = -gx2; 		/* rule 106 */
	ANGLE_FORCE_ACCUMULATE(I2, 0, fx2 );
	 tx501 = -y1; 		/* rule 108 */
	 tx502 = 2.*y2; 		/* rule 109 */
	 tx503 = -y3; 		/* rule 110 */
	 tx504 = tx501 + tx502 + tx503; 		/* rule 111 */
	 tzz945 = tx455*tzz917; 		/* rule 112 */
	 tx505 = tx470*tzz945; 		/* rule 113 */
	 tzz940 = tx456*tzz921; 		/* rule 114 */
	 tx506 = tx470*tzz940; 		/* rule 115 */
	 tx507 = tx504*tzz936; 		/* rule 116 */
	 tx508 = tx505 + tx506 + tx507; 		/* rule 117 */
	 gy2 = tx508*tzz931; 		/* rule 118 */
	 fy2 = -gy2; 		/* rule 119 */
	ANGLE_FORCE_ACCUMULATE(I2, 1, fy2 );
	 tx509 = -z1; 		/* rule 121 */
	 tx510 = 2.*z2; 		/* rule 122 */
	 tx511 = -z3; 		/* rule 123 */
	 tx512 = tx509 + tx510 + tx511; 		/* rule 124 */
	 tzz939 = tx457*tzz917; 		/* rule 125 */
	 tx513 = tx470*tzz939; 		/* rule 126 */
	 tzz934 = tx458*tzz921; 		/* rule 127 */
	 tx514 = tx470*tzz934; 		/* rule 128 */
	 tx515 = tx512*tzz936; 		/* rule 129 */
	 tx516 = tx513 + tx514 + tx515; 		/* rule 130 */
	 gz2 = tx516*tzz931; 		/* rule 131 */
	 fz2 = -gz2; 		/* rule 132 */
	ANGLE_FORCE_ACCUMULATE(I2, 2, fz2 );
	 tx517 = tx453*tzz936; 		/* rule 134 */
	 tx518 = -tx498; 		/* rule 135 */
	 tx519 = tx517 + tx518; 		/* rule 136 */
	 gx3 = tx519*tzz931; 		/* rule 137 */
	 fx3 = -gx3; 		/* rule 138 */
	ANGLE_FORCE_ACCUMULATE(I3, 0, fx3 );
	 tx520 = tx455*tzz936; 		/* rule 140 */
	 tx521 = -tx506; 		/* rule 141 */
	 tx522 = tx520 + tx521; 		/* rule 142 */
	 gy3 = tx522*tzz931; 		/* rule 143 */
	 fy3 = -gy3; 		/* rule 144 */
	ANGLE_FORCE_ACCUMULATE(I3, 1, fy3 );
	 tx523 = tx457*tzz936; 		/* rule 146 */
	 tx524 = -tx514; 		/* rule 147 */
	 tx525 = tx523 + tx524; 		/* rule 148 */
	 gz3 = tx525*tzz931; 		/* rule 149 */
	 fz3 = -gz3; 		/* rule 150 */
	ANGLE_FORCE_ACCUMULATE(I3, 2, fz3 );
	#ifdef ANGLE_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx908 = power2(tx901); 		/* rule 154 */
	 tx909 = reciprocal(tx908); 		/* rule 155 */
	 tx526 = tx902*tx909; 		/* rule 156 */
	 tx527 = power2(tx476); 		/* rule 157 */
	 tzz929 = -2.*tx477; 		/* rule 158 */
	 tzz932 = tx476*tzz929; 		/* rule 159 */
	 tzz933 = tx470*tzz932; 		/* rule 160 */
	 tx528 = tx454*tzz933; 		/* rule 161 */
	 tzz974 = -2.*tzz917; 		/* rule 162 */
	 tx529 = tx465*tzz974; 		/* rule 163 */
	 tzz923 = -tzz917; 		/* rule 164 */
	 tx530 = tx470*tzz923; 		/* rule 165 */
	 tzz926 = tx526*tx906; 		/* rule 166 */
	 tzz927 = tx470*tzz926; 		/* rule 167 */
	 tzz955 = 3*tzz927; 		/* rule 168 */
	 tx531 = tx459*tzz955; 		/* rule 169 */
	 tzz950 = tx478*tx527; 		/* rule 170 */
	 tzz965 = tx477*tzz950; 		/* rule 171 */
	 tzz968 = 2.*tzz965; 		/* rule 172 */
	 tx532 = tx453*tzz968; 		/* rule 173 */
	 tx910 = tx483; 		/* rule 174 */
	 tx911 = tx485; 		/* rule 175 */
	 tx912 = reciprocal(tx910); 		/* rule 176 */
	 tx533 = tx911*tx912; 		/* rule 177 */
	 tx534 = tx912; 		/* rule 178 */
	 tx535 = power2(tx484); 		/* rule 179 */
	 tx536 = tx529 + tx530 + tx531; 		/* rule 180 */
	 tx537 = tx528 + tx532; 		/* rule 181 */
	 tzz919 = kt*tx534; 		/* rule 182 */
	 tzz920 = 2.*tzz919; 		/* rule 183 */
	 tx538 = tx535*tzz920; 		/* rule 184 */
	 tzz951 = tx911*tzz916; 		/* rule 185 */
	 tx539 = tx536*tzz951; 		/* rule 186 */
	 tzz918 = tx533*tzz915; 		/* rule 187 */
	 tx540 = tx484*tx537*tzz918; 		/* rule 188 */
	 dhx1x1 = tx538 + tx539 + tx540; 		/* rule 189 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx541 = tx456*tzz933; 		/* rule 191 */
	 tx542 = tx466*tzz974; 		/* rule 192 */
	 tx543 = tx461*tzz955; 		/* rule 193 */
	 tx544 = tx455*tzz968; 		/* rule 194 */
	 tx545 = power2(tx488); 		/* rule 195 */
	 tx546 = tx530 + tx542 + tx543; 		/* rule 196 */
	 tx547 = tx541 + tx544; 		/* rule 197 */
	 tx548 = tx545*tzz920; 		/* rule 198 */
	 tx549 = tx546*tzz951; 		/* rule 199 */
	 tzz975 = tx547*tzz918; 		/* rule 200 */
	 tx550 = tx488*tzz975; 		/* rule 201 */
	 dhy1y1 = tx548 + tx549 + tx550; 		/* rule 202 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx551 = tx458*tzz933; 		/* rule 204 */
	 tx552 = tx467*tzz974; 		/* rule 205 */
	 tx553 = tx463*tzz955; 		/* rule 206 */
	 tx554 = tx457*tzz968; 		/* rule 207 */
	 tx555 = power2(tx491); 		/* rule 208 */
	 tx556 = tx530 + tx552 + tx553; 		/* rule 209 */
	 tx557 = tx551 + tx554; 		/* rule 210 */
	 tx558 = tx555*tzz920; 		/* rule 211 */
	 tx559 = tx556*tzz951; 		/* rule 212 */
	 tzz963 = tx557*tzz918; 		/* rule 213 */
	 tx560 = tx491*tzz963; 		/* rule 214 */
	 dhz1z1 = tx558 + tx559 + tx560; 		/* rule 215 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 tx913 = power2(tx905); 		/* rule 217 */
	 tx914 = reciprocal(tx913); 		/* rule 218 */
	 tx561 = tx906*tx914; 		/* rule 219 */
	 tx562 = power2(tx477); 		/* rule 220 */
	 tx563 = 2.*tzz936; 		/* rule 221 */
	 tx564 = -(tx470*tzz921); 		/* rule 222 */
	 tzz924 = tx495*tzz922; 		/* rule 223 */
	 tzz970 = 2.*tzz924; 		/* rule 224 */
	 tx565 = tx465*tzz970; 		/* rule 225 */
	 tx566 = tx496*tzz933; 		/* rule 226 */
	 tx567 = 2.*tx496*tzz953; 		/* rule 227 */
	 tzz978 = tx496*tzz946; 		/* rule 228 */
	 tx568 = 2.*tzz978; 		/* rule 229 */
	 tx569 = tx453*tzz929*tzz950; 		/* rule 230 */
	 tzz925 = tx561*tx902; 		/* rule 231 */
	 tzz928 = tx470*tzz925; 		/* rule 232 */
	 tzz954 = 3*tzz928; 		/* rule 233 */
	 tx570 = tx460*tzz954; 		/* rule 234 */
	 tzz948 = tx562*tzz943; 		/* rule 235 */
	 tzz972 = -2.*tzz948; 		/* rule 236 */
	 tx571 = tx454*tzz972; 		/* rule 237 */
	 tx572 = power2(tx500); 		/* rule 238 */
	 tzz994 = tx563 + tx564; 		/* rule 239 */
	 tzz995 = tx530 + tzz994; 		/* rule 240 */
	 tx573 = tx531 + tx565 + tx567 + tx568 + tx570 + tzz995; 		/* rule 241 */
	 tx574 = tx566 + tx569 + tx571; 		/* rule 242 */
	 tx575 = tx572*tzz920; 		/* rule 243 */
	 tx576 = tx573*tzz951; 		/* rule 244 */
	 tzz961 = tx574*tzz918; 		/* rule 245 */
	 tx577 = tx500*tzz961; 		/* rule 246 */
	 dhx2x2 = tx575 + tx576 + tx577; 		/* rule 247 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 tx578 = tx466*tzz970; 		/* rule 249 */
	 tx579 = tx504*tzz933; 		/* rule 250 */
	 tx580 = 2.*tx504*tzz945; 		/* rule 251 */
	 tzz977 = tx504*tzz940; 		/* rule 252 */
	 tx581 = 2.*tzz977; 		/* rule 253 */
	 tx582 = tx455*tzz929*tzz950; 		/* rule 254 */
	 tx583 = tx462*tzz954; 		/* rule 255 */
	 tx584 = tx456*tzz972; 		/* rule 256 */
	 tx585 = power2(tx508); 		/* rule 257 */
	 tx586 = tx543 + tx578 + tx580 + tx581 + tx583 + tzz995; 		/* rule 258 */
	 tx587 = tx579 + tx582 + tx584; 		/* rule 259 */
	 tx588 = tx585*tzz920; 		/* rule 260 */
	 tx589 = tx586*tzz951; 		/* rule 261 */
	 tzz958 = tx587*tzz918; 		/* rule 262 */
	 tx590 = tx508*tzz958; 		/* rule 263 */
	 dhy2y2 = tx588 + tx589 + tx590; 		/* rule 264 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 tx591 = tx467*tzz970; 		/* rule 266 */
	 tx592 = tx512*tzz933; 		/* rule 267 */
	 tx593 = 2.*tx512*tzz939; 		/* rule 268 */
	 tzz976 = tx512*tzz934; 		/* rule 269 */
	 tx594 = 2.*tzz976; 		/* rule 270 */
	 tx595 = tx457*tzz929*tzz950; 		/* rule 271 */
	 tx596 = tx464*tzz954; 		/* rule 272 */
	 tx597 = tx458*tzz972; 		/* rule 273 */
	 tx598 = power2(tx516); 		/* rule 274 */
	 tx599 = tx553 + tx591 + tx593 + tx594 + tx596 + tzz995; 		/* rule 275 */
	 tx600 = tx592 + tx595 + tx597; 		/* rule 276 */
	 tx601 = tx598*tzz920; 		/* rule 277 */
	 tx602 = tx599*tzz951; 		/* rule 278 */
	 tzz947 = tx600*tzz918; 		/* rule 279 */
	 tx603 = tx516*tzz947; 		/* rule 280 */
	 dhz2z2 = tx601 + tx602 + tx603; 		/* rule 281 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	 tx604 = tx453*tzz933; 		/* rule 283 */
	 tzz973 = -2.*tzz921; 		/* rule 284 */
	 tx605 = tx465*tzz973; 		/* rule 285 */
	 tzz969 = 2.*tzz948; 		/* rule 286 */
	 tx606 = tx454*tzz969; 		/* rule 287 */
	 tx607 = power2(tx519); 		/* rule 288 */
	 tx608 = tx564 + tx570 + tx605; 		/* rule 289 */
	 tx609 = tx604 + tx606; 		/* rule 290 */
	 tx610 = tx607*tzz920; 		/* rule 291 */
	 tx611 = tx608*tzz951; 		/* rule 292 */
	 tzz941 = tx609*tzz918; 		/* rule 293 */
	 tx612 = tx519*tzz941; 		/* rule 294 */
	 dhx3x3 = tx610 + tx611 + tx612; 		/* rule 295 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 0, dhx3x3);
	 tx613 = tx455*tzz933; 		/* rule 297 */
	 tx614 = tx466*tzz973; 		/* rule 298 */
	 tx615 = tx456*tzz969; 		/* rule 299 */
	 tx616 = power2(tx522); 		/* rule 300 */
	 tx617 = tx564 + tx583 + tx614; 		/* rule 301 */
	 tx618 = tx613 + tx615; 		/* rule 302 */
	 tx619 = tx616*tzz920; 		/* rule 303 */
	 tx620 = tx617*tzz951; 		/* rule 304 */
	 tzz937 = tx618*tzz918; 		/* rule 305 */
	 tx621 = tx522*tzz937; 		/* rule 306 */
	 dhy3y3 = tx619 + tx620 + tx621; 		/* rule 307 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 1, dhy3y3);
	 tx622 = tx457*tzz933; 		/* rule 309 */
	 tx623 = tx467*tzz973; 		/* rule 310 */
	 tx624 = tx458*tzz969; 		/* rule 311 */
	 tx625 = power2(tx525); 		/* rule 312 */
	 tx626 = tx564 + tx596 + tx623; 		/* rule 313 */
	 tx627 = tx622 + tx624; 		/* rule 314 */
	 tx628 = tx625*tzz920; 		/* rule 315 */
	 tx629 = tx626*tzz951; 		/* rule 316 */
	 tzz930 = tx627*tzz918; 		/* rule 317 */
	 tx630 = tx525*tzz930; 		/* rule 318 */
	 dhz3z3 = tx628 + tx629 + tx630; 		/* rule 319 */
	ANGLE_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I3, 2, dhz3z3);
	#ifdef ANGLE_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tzz952 = tx455*tzz923; 		/* rule 323 */
	 tx631 = tx454*tzz952; 		/* rule 324 */
	 tzz960 = tx453*tzz923; 		/* rule 325 */
	 tx632 = tx456*tzz960; 		/* rule 326 */
	 tzz991 = tx453*tx455; 		/* rule 327 */
	 tx633 = tzz955*tzz991; 		/* rule 328 */
	 tx634 = tx631 + tx632 + tx633; 		/* rule 329 */
	 tx635 = tx484*tx488*tzz920; 		/* rule 330 */
	 tx636 = tx484*tzz975; 		/* rule 331 */
	 tx637 = tx634*tzz951; 		/* rule 332 */
	 ohx1y1 = tx635 + tx636 + tx637; 		/* rule 333 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tzz944 = tx457*tzz923; 		/* rule 335 */
	 tx638 = tx454*tzz944; 		/* rule 336 */
	 tx639 = tx458*tzz960; 		/* rule 337 */
	 tzz985 = tx457*tzz955; 		/* rule 338 */
	 tx640 = tx453*tzz985; 		/* rule 339 */
	 tx641 = tx638 + tx639 + tx640; 		/* rule 340 */
	 tzz979 = tx491*tzz920; 		/* rule 341 */
	 tx642 = tx484*tzz979; 		/* rule 342 */
	 tx643 = tx484*tzz963; 		/* rule 343 */
	 tx644 = tx641*tzz951; 		/* rule 344 */
	 ohx1z1 = tx642 + tx643 + tx644; 		/* rule 345 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx645 = -tzz936; 		/* rule 347 */
	 tx646 = tx465*tzz917; 		/* rule 348 */
	 tx647 = tx470*tzz917; 		/* rule 349 */
	 tx648 = tx460*tzz921; 		/* rule 350 */
	 tzz935 = -tzz924; 		/* rule 351 */
	 tx649 = tx465*tzz935; 		/* rule 352 */
	 tx650 = tx496*tzz960; 		/* rule 353 */
	 tzz957 = -3*tzz927; 		/* rule 354 */
	 tx651 = tx459*tzz957; 		/* rule 355 */
	 tzz993 = tx645 + tx647; 		/* rule 356 */
	 tx652 = tx646 + tx648 + tx649 + tx650 + tx651 + tzz993; 		/* rule 357 */
	 tzz964 = tx500*tzz920; 		/* rule 358 */
	 tx653 = tx484*tzz964; 		/* rule 359 */
	 tx654 = tx484*tzz961; 		/* rule 360 */
	 tx655 = tx652*tzz951; 		/* rule 361 */
	 ohx1x2 = tx653 + tx654 + tx655; 		/* rule 362 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx656 = tx454*tzz945; 		/* rule 364 */
	 tx657 = tx454*tzz940; 		/* rule 365 */
	 tzz990 = tx453*tx456; 		/* rule 366 */
	 tx658 = tzz935*tzz990; 		/* rule 367 */
	 tx659 = tx504*tzz960; 		/* rule 368 */
	 tx660 = tzz957*tzz991; 		/* rule 369 */
	 tzz1004 = tx657 + tx660; 		/* rule 370 */
	 tx661 = tx656 + tx658 + tx659 + tzz1004; 		/* rule 371 */
	 tzz962 = tx508*tzz920; 		/* rule 372 */
	 tx662 = tx484*tzz962; 		/* rule 373 */
	 tx663 = tx484*tzz958; 		/* rule 374 */
	 tx664 = tx661*tzz951; 		/* rule 375 */
	 ohx1y2 = tx662 + tx663 + tx664; 		/* rule 376 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx665 = tx454*tzz939; 		/* rule 378 */
	 tx666 = tx454*tzz934; 		/* rule 379 */
	 tzz982 = tx458*tzz935; 		/* rule 380 */
	 tx667 = tx453*tzz982; 		/* rule 381 */
	 tx668 = tx512*tzz960; 		/* rule 382 */
	 tzz984 = tx457*tzz957; 		/* rule 383 */
	 tx669 = tx453*tzz984; 		/* rule 384 */
	 tzz1003 = tx666 + tx669; 		/* rule 385 */
	 tx670 = tx665 + tx667 + tx668 + tzz1003; 		/* rule 386 */
	 tzz959 = tx516*tzz920; 		/* rule 387 */
	 tx671 = tx484*tzz959; 		/* rule 388 */
	 tx672 = tx484*tzz947; 		/* rule 389 */
	 tx673 = tx670*tzz951; 		/* rule 390 */
	 ohx1z2 = tx671 + tx672 + tx673; 		/* rule 391 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tx674 = tzz936; 		/* rule 393 */
	 tx675 = tx459*tzz923; 		/* rule 394 */
	 tx676 = tx465*tzz924; 		/* rule 395 */
	 tx677 = -tx648; 		/* rule 396 */
	 tx678 = tx674 + tx675 + tx676 + tx677; 		/* rule 397 */
	 tzz949 = tx519*tzz920; 		/* rule 398 */
	 tx679 = tx484*tzz949; 		/* rule 399 */
	 tx680 = tx484*tzz941; 		/* rule 400 */
	 tx681 = tx678*tzz951; 		/* rule 401 */
	 ohx1x3 = tx679 + tx680 + tx681; 		/* rule 402 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 0, ohx1x3);
	 tx682 = tx453*tzz952; 		/* rule 404 */
	 tx683 = tzz924*tzz990; 		/* rule 405 */
	 tx684 = -tx657; 		/* rule 406 */
	 tzz1002 = tx682 + tx684; 		/* rule 407 */
	 tx685 = tx683 + tzz1002; 		/* rule 408 */
	 tzz942 = tx522*tzz920; 		/* rule 409 */
	 tx686 = tx484*tzz942; 		/* rule 410 */
	 tx687 = tx484*tzz937; 		/* rule 411 */
	 tx688 = tx685*tzz951; 		/* rule 412 */
	 ohx1y3 = tx686 + tx687 + tx688; 		/* rule 413 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 1, ohx1y3);
	 tx689 = tx453*tzz944; 		/* rule 415 */
	 tzz983 = tx458*tzz924; 		/* rule 416 */
	 tx690 = tx453*tzz983; 		/* rule 417 */
	 tx691 = -tx666; 		/* rule 418 */
	 tzz1001 = tx689 + tx691; 		/* rule 419 */
	 tx692 = tx690 + tzz1001; 		/* rule 420 */
	 tzz938 = tx525*tzz920; 		/* rule 421 */
	 tx693 = tx484*tzz938; 		/* rule 422 */
	 tx694 = tx484*tzz930; 		/* rule 423 */
	 tx695 = tx692*tzz951; 		/* rule 424 */
	 ohx1z3 = tx693 + tx694 + tx695; 		/* rule 425 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 2, ohx1z3);
	 tx696 = tx456*tzz944; 		/* rule 427 */
	 tx697 = tx458*tzz952; 		/* rule 428 */
	 tx698 = tx455*tzz985; 		/* rule 429 */
	 tx699 = tx696 + tx697 + tx698; 		/* rule 430 */
	 tx700 = tx488*tzz979; 		/* rule 431 */
	 tx701 = tx488*tzz963; 		/* rule 432 */
	 tx702 = tx699*tzz951; 		/* rule 433 */
	 ohy1z1 = tx700 + tx701 + tx702; 		/* rule 434 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 tx703 = tx456*tzz953; 		/* rule 436 */
	 tzz989 = tx454*tx455; 		/* rule 437 */
	 tx704 = tzz935*tzz989; 		/* rule 438 */
	 tx705 = tx496*tzz952; 		/* rule 439 */
	 tx706 = tx703 + tx704 + tx705 + tzz1004; 		/* rule 440 */
	 tx707 = tx488*tzz964; 		/* rule 441 */
	 tx708 = tx488*tzz961; 		/* rule 442 */
	 tx709 = tx706*tzz951; 		/* rule 443 */
	 ohy1x2 = tx707 + tx708 + tx709; 		/* rule 444 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx710 = tx466*tzz917; 		/* rule 446 */
	 tx711 = tx462*tzz921; 		/* rule 447 */
	 tx712 = tx466*tzz935; 		/* rule 448 */
	 tx713 = tx504*tzz952; 		/* rule 449 */
	 tx714 = tx461*tzz957; 		/* rule 450 */
	 tx715 = tx710 + tx711 + tx712 + tx713 + tx714 + tzz993; 		/* rule 451 */
	 tx716 = tx488*tzz962; 		/* rule 452 */
	 tx717 = tx488*tzz958; 		/* rule 453 */
	 tx718 = tx715*tzz951; 		/* rule 454 */
	 ohy1y2 = tx716 + tx717 + tx718; 		/* rule 455 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx719 = tx456*tzz939; 		/* rule 457 */
	 tx720 = tx456*tzz934; 		/* rule 458 */
	 tx721 = tx455*tzz982; 		/* rule 459 */
	 tx722 = tx512*tzz952; 		/* rule 460 */
	 tx723 = tx455*tzz984; 		/* rule 461 */
	 tzz1000 = tx720 + tx723; 		/* rule 462 */
	 tx724 = tx719 + tx721 + tx722 + tzz1000; 		/* rule 463 */
	 tx725 = tx488*tzz959; 		/* rule 464 */
	 tx726 = tx488*tzz947; 		/* rule 465 */
	 tx727 = tx724*tzz951; 		/* rule 466 */
	 ohy1z2 = tx725 + tx726 + tx727; 		/* rule 467 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 tx728 = tzz924*tzz989; 		/* rule 469 */
	 tx729 = tx728 + tzz1002; 		/* rule 470 */
	 tx730 = tx488*tzz949; 		/* rule 471 */
	 tx731 = tx488*tzz941; 		/* rule 472 */
	 tx732 = tx729*tzz951; 		/* rule 473 */
	 ohy1x3 = tx730 + tx731 + tx732; 		/* rule 474 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 0, ohy1x3);
	 tx733 = tx461*tzz923; 		/* rule 476 */
	 tx734 = tx466*tzz924; 		/* rule 477 */
	 tx735 = -tx711; 		/* rule 478 */
	 tx736 = tx674 + tx733 + tx734 + tx735; 		/* rule 479 */
	 tx737 = tx488*tzz942; 		/* rule 480 */
	 tx738 = tx488*tzz937; 		/* rule 481 */
	 tx739 = tx736*tzz951; 		/* rule 482 */
	 ohy1y3 = tx737 + tx738 + tx739; 		/* rule 483 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 1, ohy1y3);
	 tx740 = tx455*tzz944; 		/* rule 485 */
	 tx741 = tx455*tzz983; 		/* rule 486 */
	 tx742 = -tx720; 		/* rule 487 */
	 tzz999 = tx740 + tx742; 		/* rule 488 */
	 tx743 = tx741 + tzz999; 		/* rule 489 */
	 tx744 = tx488*tzz938; 		/* rule 490 */
	 tx745 = tx488*tzz930; 		/* rule 491 */
	 tx746 = tx743*tzz951; 		/* rule 492 */
	 ohy1z3 = tx744 + tx745 + tx746; 		/* rule 493 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 2, ohy1z3);
	 tx747 = tx458*tzz953; 		/* rule 495 */
	 tzz986 = tx457*tzz935; 		/* rule 496 */
	 tx748 = tx454*tzz986; 		/* rule 497 */
	 tx749 = tx496*tzz944; 		/* rule 498 */
	 tx750 = tx747 + tx748 + tx749 + tzz1003; 		/* rule 499 */
	 tx751 = tx491*tzz964; 		/* rule 500 */
	 tx752 = tx491*tzz961; 		/* rule 501 */
	 tx753 = tx750*tzz951; 		/* rule 502 */
	 ohz1x2 = tx751 + tx752 + tx753; 		/* rule 503 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 tx754 = tx458*tzz945; 		/* rule 505 */
	 tx755 = tx456*tzz986; 		/* rule 506 */
	 tx756 = tx504*tzz944; 		/* rule 507 */
	 tx757 = tx754 + tx755 + tx756 + tzz1000; 		/* rule 508 */
	 tx758 = tx491*tzz962; 		/* rule 509 */
	 tx759 = tx491*tzz958; 		/* rule 510 */
	 tx760 = tx757*tzz951; 		/* rule 511 */
	 ohz1y2 = tx758 + tx759 + tx760; 		/* rule 512 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx761 = tx467*tzz917; 		/* rule 514 */
	 tx762 = tx464*tzz921; 		/* rule 515 */
	 tx763 = tx467*tzz935; 		/* rule 516 */
	 tx764 = tx512*tzz944; 		/* rule 517 */
	 tx765 = tx463*tzz957; 		/* rule 518 */
	 tx766 = tx761 + tx762 + tx763 + tx764 + tx765 + tzz993; 		/* rule 519 */
	 tx767 = tx491*tzz959; 		/* rule 520 */
	 tx768 = tx491*tzz947; 		/* rule 521 */
	 tx769 = tx766*tzz951; 		/* rule 522 */
	 ohz1z2 = tx767 + tx768 + tx769; 		/* rule 523 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 tzz987 = tx457*tzz924; 		/* rule 525 */
	 tx770 = tx454*tzz987; 		/* rule 526 */
	 tx771 = tx770 + tzz1001; 		/* rule 527 */
	 tx772 = tx491*tzz949; 		/* rule 528 */
	 tx773 = tx491*tzz941; 		/* rule 529 */
	 tx774 = tx771*tzz951; 		/* rule 530 */
	 ohz1x3 = tx772 + tx773 + tx774; 		/* rule 531 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 0, ohz1x3);
	 tx775 = tx456*tzz987; 		/* rule 533 */
	 tx776 = tx775 + tzz999; 		/* rule 534 */
	 tx777 = tx491*tzz942; 		/* rule 535 */
	 tx778 = tx491*tzz937; 		/* rule 536 */
	 tx779 = tx776*tzz951; 		/* rule 537 */
	 ohz1y3 = tx777 + tx778 + tx779; 		/* rule 538 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 1, ohz1y3);
	 tx780 = tx463*tzz923; 		/* rule 540 */
	 tx781 = tx467*tzz924; 		/* rule 541 */
	 tx782 = -tx762; 		/* rule 542 */
	 tx783 = tx674 + tx780 + tx781 + tx782; 		/* rule 543 */
	 tx784 = tx491*tzz938; 		/* rule 544 */
	 tx785 = tx491*tzz930; 		/* rule 545 */
	 tx786 = tx783*tzz951; 		/* rule 546 */
	 ohz1z3 = tx784 + tx785 + tx786; 		/* rule 547 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 2, ohz1z3);
	 tx787 = tx496*tzz945; 		/* rule 549 */
	 tx788 = tx496*tzz940; 		/* rule 550 */
	 tx789 = tx504*tzz953; 		/* rule 551 */
	 tx790 = tx504*tzz946; 		/* rule 552 */
	 tzz988 = tx454*tx456; 		/* rule 553 */
	 tx791 = tzz954*tzz988; 		/* rule 554 */
	 tx792 = tx633 + tx683 + tx728 + tx787 + tx788 + tx789 + tx790 + tx791; 		/* rule 555 */
	 tx793 = tx500*tzz962; 		/* rule 556 */
	 tx794 = tx500*tzz958; 		/* rule 557 */
	 tx795 = tx792*tzz951; 		/* rule 558 */
	 ohx2y2 = tx793 + tx794 + tx795; 		/* rule 559 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 tx796 = tx496*tzz939; 		/* rule 561 */
	 tx797 = tx496*tzz934; 		/* rule 562 */
	 tx798 = tx512*tzz953; 		/* rule 563 */
	 tx799 = tx512*tzz946; 		/* rule 564 */
	 tzz981 = tx458*tzz954; 		/* rule 565 */
	 tx800 = tx454*tzz981; 		/* rule 566 */
	 tx801 = tx640 + tx690 + tx770 + tx796 + tx797 + tx798 + tx799 + tx800; 		/* rule 567 */
	 tx802 = tx500*tzz959; 		/* rule 568 */
	 tx803 = tx500*tzz947; 		/* rule 569 */
	 tx804 = tx801*tzz951; 		/* rule 570 */
	 ohx2z2 = tx802 + tx803 + tx804; 		/* rule 571 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 tx805 = tx459*tzz917; 		/* rule 573 */
	 tx806 = tx465*tzz921; 		/* rule 574 */
	 tx807 = tx470*tzz921; 		/* rule 575 */
	 tx808 = -tzz978; 		/* rule 576 */
	 tzz956 = -3*tzz928; 		/* rule 577 */
	 tx809 = tx460*tzz956; 		/* rule 578 */
	 tzz992 = tx645 + tx807; 		/* rule 579 */
	 tx810 = tx649 + tx805 + tx806 + tx808 + tx809 + tzz992; 		/* rule 580 */
	 tx811 = tx500*tzz949; 		/* rule 581 */
	 tx812 = tx500*tzz941; 		/* rule 582 */
	 tx813 = tx810*tzz951; 		/* rule 583 */
	 ohx2x3 = tx811 + tx812 + tx813; 		/* rule 584 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 0, ohx2x3);
	 tx814 = tx453*tzz945; 		/* rule 586 */
	 tx815 = tx455*tzz946; 		/* rule 587 */
	 tx816 = tzz956*tzz988; 		/* rule 588 */
	 tx817 = -tx788; 		/* rule 589 */
	 tzz998 = tx814 + tx816; 		/* rule 590 */
	 tx818 = tx658 + tx815 + tx817 + tzz998; 		/* rule 591 */
	 tx819 = tx500*tzz942; 		/* rule 592 */
	 tx820 = tx500*tzz937; 		/* rule 593 */
	 tx821 = tx818*tzz951; 		/* rule 594 */
	 ohx2y3 = tx819 + tx820 + tx821; 		/* rule 595 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 1, ohx2y3);
	 tx822 = tx453*tzz939; 		/* rule 597 */
	 tx823 = tx457*tzz946; 		/* rule 598 */
	 tzz980 = tx458*tzz956; 		/* rule 599 */
	 tx824 = tx454*tzz980; 		/* rule 600 */
	 tx825 = -tx797; 		/* rule 601 */
	 tzz997 = tx822 + tx824; 		/* rule 602 */
	 tx826 = tx667 + tx823 + tx825 + tzz997; 		/* rule 603 */
	 tx827 = tx500*tzz938; 		/* rule 604 */
	 tx828 = tx500*tzz930; 		/* rule 605 */
	 tx829 = tx826*tzz951; 		/* rule 606 */
	 ohx2z3 = tx827 + tx828 + tx829; 		/* rule 607 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 2, ohx2z3);
	 tx830 = tx504*tzz939; 		/* rule 609 */
	 tx831 = tx504*tzz934; 		/* rule 610 */
	 tx832 = tx512*tzz945; 		/* rule 611 */
	 tx833 = tx512*tzz940; 		/* rule 612 */
	 tx834 = tx456*tzz981; 		/* rule 613 */
	 tx835 = tx698 + tx741 + tx775 + tx830 + tx831 + tx832 + tx833 + tx834; 		/* rule 614 */
	 tx836 = tx508*tzz959; 		/* rule 615 */
	 tx837 = tx508*tzz947; 		/* rule 616 */
	 tx838 = tx835*tzz951; 		/* rule 617 */
	 ohy2z2 = tx836 + tx837 + tx838; 		/* rule 618 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	 tx839 = tx453*tzz940; 		/* rule 620 */
	 tx840 = -tx790; 		/* rule 621 */
	 tx841 = tx704 + tx839 + tx840 + tzz998; 		/* rule 622 */
	 tx842 = tx508*tzz949; 		/* rule 623 */
	 tx843 = tx508*tzz941; 		/* rule 624 */
	 tx844 = tx841*tzz951; 		/* rule 625 */
	 ohy2x3 = tx842 + tx843 + tx844; 		/* rule 626 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 0, ohy2x3);
	 tx845 = tx461*tzz917; 		/* rule 628 */
	 tx846 = tx466*tzz921; 		/* rule 629 */
	 tx847 = -tzz977; 		/* rule 630 */
	 tx848 = tx462*tzz956; 		/* rule 631 */
	 tx849 = tx712 + tx845 + tx846 + tx847 + tx848 + tzz992; 		/* rule 632 */
	 tx850 = tx508*tzz942; 		/* rule 633 */
	 tx851 = tx508*tzz937; 		/* rule 634 */
	 tx852 = tx849*tzz951; 		/* rule 635 */
	 ohy2y3 = tx850 + tx851 + tx852; 		/* rule 636 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 1, ohy2y3);
	 tx853 = tx455*tzz939; 		/* rule 638 */
	 tx854 = tx457*tzz940; 		/* rule 639 */
	 tx855 = tx456*tzz980; 		/* rule 640 */
	 tx856 = -tx831; 		/* rule 641 */
	 tzz996 = tx853 + tx855; 		/* rule 642 */
	 tx857 = tx721 + tx854 + tx856 + tzz996; 		/* rule 643 */
	 tx858 = tx508*tzz938; 		/* rule 644 */
	 tx859 = tx508*tzz930; 		/* rule 645 */
	 tx860 = tx857*tzz951; 		/* rule 646 */
	 ohy2z3 = tx858 + tx859 + tx860; 		/* rule 647 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 2, ohy2z3);
	 tx861 = tx453*tzz934; 		/* rule 649 */
	 tx862 = -tx799; 		/* rule 650 */
	 tx863 = tx748 + tx861 + tx862 + tzz997; 		/* rule 651 */
	 tx864 = tx516*tzz949; 		/* rule 652 */
	 tx865 = tx516*tzz941; 		/* rule 653 */
	 tx866 = tx863*tzz951; 		/* rule 654 */
	 ohz2x3 = tx864 + tx865 + tx866; 		/* rule 655 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 0, ohz2x3);
	 tx867 = tx455*tzz934; 		/* rule 657 */
	 tx868 = -tx833; 		/* rule 658 */
	 tx869 = tx755 + tx867 + tx868 + tzz996; 		/* rule 659 */
	 tx870 = tx516*tzz942; 		/* rule 660 */
	 tx871 = tx516*tzz937; 		/* rule 661 */
	 tx872 = tx869*tzz951; 		/* rule 662 */
	 ohz2y3 = tx870 + tx871 + tx872; 		/* rule 663 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 1, ohz2y3);
	 tx873 = tx463*tzz917; 		/* rule 665 */
	 tx874 = tx467*tzz921; 		/* rule 666 */
	 tx875 = -tzz976; 		/* rule 667 */
	 tx876 = tx464*tzz956; 		/* rule 668 */
	 tx877 = tx763 + tx873 + tx874 + tx875 + tx876 + tzz992; 		/* rule 669 */
	 tx878 = tx516*tzz938; 		/* rule 670 */
	 tx879 = tx516*tzz930; 		/* rule 671 */
	 tx880 = tx877*tzz951; 		/* rule 672 */
	 ohz2z3 = tx878 + tx879 + tx880; 		/* rule 673 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 2, ohz2z3);
	 tx881 = -tx815; 		/* rule 675 */
	 tx882 = -tx839; 		/* rule 676 */
	 tx883 = tx791 + tx881 + tx882; 		/* rule 677 */
	 tx884 = tx519*tzz942; 		/* rule 678 */
	 tx885 = tx519*tzz937; 		/* rule 679 */
	 tx886 = tx883*tzz951; 		/* rule 680 */
	 ohx3y3 = tx884 + tx885 + tx886; 		/* rule 681 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 1, ohx3y3);
	 tx887 = -tx823; 		/* rule 683 */
	 tx888 = -tx861; 		/* rule 684 */
	 tx889 = tx800 + tx887 + tx888; 		/* rule 685 */
	 tx890 = tx519*tzz938; 		/* rule 686 */
	 tx891 = tx519*tzz930; 		/* rule 687 */
	 tx892 = tx889*tzz951; 		/* rule 688 */
	 ohx3z3 = tx890 + tx891 + tx892; 		/* rule 689 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 2, ohx3z3);
	 tx893 = -tx854; 		/* rule 691 */
	 tx894 = -tx867; 		/* rule 692 */
	 tx895 = tx834 + tx893 + tx894; 		/* rule 693 */
	 tx896 = tx522*tzz938; 		/* rule 694 */
	 tx897 = tx522*tzz930; 		/* rule 695 */
	 tx898 = tx895*tzz951; 		/* rule 696 */
	 ohy3z3 = tx896 + tx897 + tx898; 		/* rule 697 */
	ANGLE_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 2, ohy3z3);
	} /*if calcOffDiagonalHessian */ 
	#endif /* ANGLE_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* ANGLE_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* ANGLE_CALC_FORCE ]*/
