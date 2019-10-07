{Name -> "NONBONDRBPB", AdditionalCDeclares -> "", 
 Input -> {ak, bk, ck, dk, xk, yk, zk, al, bl, cl, dl, xl, yl, zl}, 
 Output -> {Energy, DeltaX, DeltaY, DeltaZ, fam, fbm, fcm, fdm, fxm, fym, 
   fzm, fan, fbn, fcn, fdn, fxn, fyn, fzn, nonbondRBDeviation}, 
 EnergyFunction -> EnergyFunction, DerivativeVariables -> 
  {ak, bk, ck, dk, xk, yk, zk, al, bl, cl, dl, xl, yl, zl}, 
 Rules -> {CCode["NONBONDRB_SET_PARAMETER(I1);"], 
   CCode["NONBONDRB_SET_PARAMETER(I2);"], 
   CCode["NONBONDRB_SET_PARAMETER(dQ1Q2);"], 
   CCode["NONBONDRB_SET_PARAMETER(dA);"], 
   CCode["NONBONDRB_SET_PARAMETER(dC);"], 
   CCode["NONBONDRB_SET_POSITION(am,I1,0);"], 
   CCode["NONBONDRB_SET_POSITION(bm,I1,1);"], 
   CCode["NONBONDRB_SET_POSITION(cm,I1,2);"], 
   CCode["NONBONDRB_SET_POSITION(dm,I1,3);"], 
   CCode["NONBONDRB_SET_POSITION(xm,I1,4);"], 
   CCode["NONBONDRB_SET_POSITION(ym,I1,5);"], 
   CCode["NONBONDRB_SET_POSITION(zm,I1,6);"], 
   CCode["NONBONDRB_SET_POSITION(an,I2,0);"], 
   CCode["NONBONDRB_SET_POSITION(bn,I2,1);"], 
   CCode["NONBONDRB_SET_POSITION(cn,I2,2);"], 
   CCode["NONBONDRB_SET_POSITION(dn,I2,3);"], 
   CCode["NONBONDRB_SET_POSITION(xn,I2,4);"], 
   CCode["NONBONDRB_SET_POSITION(yn,I2,5);"], 
   CCode["NONBONDRB_SET_POSITION(zn,I2,6);"], 
   CCode["NONBONDRB_SET_POINT(pxm,ea1,getX());"], 
   CCode["NONBONDRB_SET_POINT(pym,ea1,getY());"], 
   CCode["NONBONDRB_SET_POINT(pzm,ea1,getZ());"], 
   CCode["NONBONDRB_SET_POINT(pxn,ea2,getX());"], 
   CCode["NONBONDRB_SET_POINT(pyn,ea2,getY());"], 
   CCode["NONBONDRB_SET_POINT(pzn,ea2,getZ());"], power2[am] -> tx318, 
   power2[an] -> tx319, power2[bm] -> tx320, power2[bn] -> tx321, 
   power2[cm] -> tx322, power2[cn] -> tx323, power2[dm] -> tx324, 
   power2[dn] -> tx325, cm*pym -> tzz679, bm*tzz679 -> tx326, -pym -> tzz666, 
   am*dm*tzz666 -> tx327, cn*pyn -> tzz678, bn*tzz678 -> tx328, 
   -pyn -> tzz671, an*dn*tzz671 -> tx329, am*cm -> tzz680, 
   pzm*tzz680 -> tx330, dm*pzm -> tzz675, bm*tzz675 -> tx331, 
   cn*pzn -> tzz677, an*tzz677 -> tx332, bn*pzn -> tzz667, 
   dn*tzz667 -> tx333, -pxm -> tzz659, tx322*tzz659 -> tx334, 
   -tx323 -> tzz682, pxn*tzz682 -> tx335, tx324*tzz659 -> tx336, 
   tx318 + tx320 + tx322 + tx324 -> tx337, -(pxn*tx325) -> tx338, 
   tx319 + tx321 + tx323 + tx325 -> tx339, 
   tx326 + tx327 + tx330 + tx331 + tx334 + tx336 -> tx340, 
   reciprocal[tx337] -> tx341, tx332 + tx333 -> tzz699, 
   tx329 + tzz699 -> tzz700, tx328 + tzz700 -> tzz701, 
   tx335 + tx338 + tzz701 -> tx342, reciprocal[tx339] -> tx343, 
   -pxn -> tx344, 2.*tx341 -> tzz646, tx340*tzz646 -> tx345, 
   -2.*tx343 -> tzz645, tx342*tzz645 -> tx346, -xn -> tx347, 
   tx347 + xm -> tzz696, tx345 + tzz696 -> tzz697, tx344 + tzz697 -> tzz698, 
   pxm + tzz698 -> tzz704, tx346 + tzz704 -> DeltaX, bm*cm*pxm -> tx348, 
   dm*pxm -> tzz676, am*tzz676 -> tx349, bn*cn*pxn -> tx350, 
   dn*pxn -> tzz674, an*tzz674 -> tx351, -pzm -> tzz658, am*bm -> tzz681, 
   tzz658*tzz681 -> tx352, cm*tzz675 -> tx353, -(an*tzz667) -> tx354, 
   dn*tzz677 -> tx355, tx320*tzz666 -> tx356, tx321*tzz671 -> tx357, 
   tx324*tzz666 -> tx358, tx325*tzz671 -> tx359, 
   tx348 + tx349 + tx352 + tx353 + tx356 + tx358 -> tx360, 
   tx354 + tx355 -> tzz693, tx351 + tzz693 -> tzz694, 
   tx350 + tzz694 -> tzz695, tx357 + tx359 + tzz695 -> tx361, 
   tzz671 -> tx362, tx360*tzz646 -> tx363, tx361*tzz645 -> tx364, 
   -yn -> tx365, tx365 + ym -> tzz690, tx363 + tzz690 -> tzz691, 
   tx362 + tzz691 -> tzz692, pym + tzz692 -> tzz703, 
   tx364 + tzz703 -> DeltaY, tzz659*tzz680 -> tx366, bm*tzz676 -> tx367, 
   bn*tzz674 -> tx368, pym*tzz681 -> tx369, dm*tzz679 -> tx370, 
   an*bn*pyn -> tx371, dn*tzz678 -> tx372, tx320*tzz658 -> tx373, 
   -(pzn*tx321) -> tx374, tx322*tzz658 -> tx375, pzn*tzz682 -> tx376, 
   an*cn*tx344 -> tx377, tx366 + tx367 + tx369 + tx370 + tx373 + tx375 -> 
    tx378, tx372 + tx377 -> tzz687, tx371 + tzz687 -> tzz688, 
   tx368 + tzz688 -> tzz689, tx374 + tx376 + tzz689 -> tx379, -pzn -> tx380, 
   tx378*tzz646 -> tx381, tx379*tzz645 -> tx382, -zn -> tx383, 
   tx383 + zm -> tzz684, tx381 + tzz684 -> tzz685, tx380 + tzz685 -> tzz686, 
   pzm + tzz686 -> tzz702, tx382 + tzz702 -> DeltaZ, tx323*tx344 -> tx384, 
   tx325*tx344 -> tx385, tx321*tx362 -> tx386, tx325*tx362 -> tx387, 
   tx321*tx380 -> tx388, tx323*tx380 -> tx389, tx384 + tx385 + tzz701 -> 
    tx390, tx386 + tx387 + tzz695 -> tx391, tx388 + tx389 + tzz689 -> tx392, 
   tx390*tzz645 -> tx393, tx391*tzz645 -> tx394, tx392*tzz645 -> tx395, 
   tx393 + tzz704 -> tx396, tx394 + tzz703 -> tx397, tx395 + tzz702 -> tx398, 
   PBX[tx396] -> tx399, PBY[tx397] -> tx400, PBZ[tx398] -> tx401, 
   power2[tx399] -> tx402, power2[tx400] -> tx403, power2[tx401] -> tx404, 
   tx402 + tx403 + tx404 -> tx405, reciprocal[tx405] -> tx683, 
   power2[tx683] -> tx643, tx683 -> tx644, tx643*tx644 -> tx635, 
   power2[tx635] -> tx406, power2[tx644] -> tx636, tx644 -> tx637, 
   tx636*tx637 -> tx407, mysqrt[tx405] -> tx638, reciprocal[tx638] -> tx408, 
   dA*tx406 -> tx409, -(dC*tx407) -> tx410, dQ1Q2*tx408 -> tx411, 
   tx409 + tx410 + tx411 -> Energy, 
   CCode["NONBONDRB_ENERGY_ACCUMULATE(Energy);"], 
   CCode["#ifdef NONBONDRB_CALC_FORCE //["], CCode["if ( calcForce ) {"], 
   CCode["BAIL_OUT_IF_CUTOFF(DeltaX,DeltaY,DeltaZ)"], cm*tzz659 -> tx412, 
   tzz676 -> tx413, bm*pym -> tx414, dm*tzz666 -> tx415, bm*tzz658 -> tx416, 
   cm*pzm -> tx417, power2[tx341] -> tx418, tx412 + tx414 -> tx419, 
   tx413 + tx416 -> tx420, tx415 + tx417 -> tx421, -4.*tx418 -> tzz651, 
   tx340*tzz651 -> tzz665, am*tzz665 -> tx422, tx360*tzz651 -> tzz664, 
   am*tzz664 -> tx423, tx378*tzz651 -> tzz663, am*tzz663 -> tx424, 
   tx419*tzz646 -> tx425, tx420*tzz646 -> tx426, tx421*tzz646 -> tx427, 
   tx424 + tx425 -> tx428, tx423 + tx426 -> tx429, tx422 + tx427 -> tx430, 
   2.*tx401 -> tzz653, tx428*tzz653 -> tx431, 2.*tx400 -> tzz654, 
   tx429*tzz654 -> tx432, 2.*tx399 -> tzz655, tx430*tzz655 -> tx433, 
   power2[tx643] -> tx639, tx635*tx639 -> tx434, power2[tx636] -> tx435, 
   tx405 -> tx640, tx408 -> tx641, reciprocal[tx640] -> tx642, 
   tx641*tx642 -> tx436, tx431 + tx432 + tx433 -> tx437, dA*tx434 -> tzz649, 
   -6.*tzz649 -> tzz657, tx437*tzz657 -> tx438, dC*tx435 -> tzz648, 
   3*tzz648 -> tzz652, tx437*tzz652 -> tx439, dQ1Q2*tx436 -> tzz647, 
   -0.5*tzz647 -> tzz656, tx437*tzz656 -> tx440, 
   tx438 + tx439 + tx440 -> gam, -gam -> fam, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I1, 0, fam );"], cm*pxm -> tx441, 
   am*pym -> tx442, tzz679 -> tx443, am*tzz658 -> tx444, -2.*bm*pzm -> tx445, 
   tzz675 -> tx446, -2.*tx414 -> tx447, tx413 + tx442 + tx445 -> tx448, 
   tx443 + tx446 -> tx449, tx441 + tx444 + tx447 -> tx450, 
   bm*tzz665 -> tx451, bm*tzz664 -> tx452, bm*tzz663 -> tx453, 
   tx448*tzz646 -> tx454, tx449*tzz646 -> tx455, tx450*tzz646 -> tx456, 
   tx453 + tx454 -> tx457, tx451 + tx455 -> tx458, tx452 + tx456 -> tx459, 
   tx457*tzz653 -> tx460, tx458*tzz655 -> tx461, tx459*tzz654 -> tx462, 
   tx460 + tx461 + tx462 -> tx463, tx463*tzz657 -> tx464, 
   tx463*tzz652 -> tx465, tx463*tzz656 -> tx466, 
   tx464 + tx465 + tx466 -> gbm, -gbm -> fbm, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I1, 1, fbm );"], am*tzz659 -> tx467, 
   bm*pxm -> tx468, dm*pym -> tx469, am*pzm -> tx470, -2.*tx417 -> tx471, 
   -2.*tx441 -> tx472, tx446 + tx468 -> tx473, tx467 + tx469 + tx471 -> 
    tx474, tx414 + tx470 + tx472 -> tx475, cm*tzz665 -> tx476, 
   cm*tzz664 -> tx477, cm*tzz663 -> tx478, tx473*tzz646 -> tx479, 
   tx474*tzz646 -> tx480, tx475*tzz646 -> tx481, tx477 + tx479 -> tx482, 
   tx478 + tx480 -> tx483, tx476 + tx481 -> tx484, tx482*tzz654 -> tx485, 
   tx483*tzz653 -> tx486, tx484*tzz655 -> tx487, 
   tx485 + tx486 + tx487 -> tx488, tx488*tzz657 -> tx489, 
   tx488*tzz652 -> tx490, tx488*tzz656 -> tx491, 
   tx489 + tx490 + tx491 -> gcm, -gcm -> fcm, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I1, 2, fcm );"], am*pxm -> tx492, 
   bm*pzm -> tx493, -2.*tx413 -> tx494, -tx442 -> tx495, -2.*tx469 -> tx496, 
   tx443 + tx468 -> tx497, tx493 + tx494 + tx495 -> tx498, 
   tx417 + tx492 + tx496 -> tx499, dm*tzz665 -> tx500, dm*tzz664 -> tx501, 
   dm*tzz663 -> tx502, tx497*tzz646 -> tx503, tx498*tzz646 -> tx504, 
   tx499*tzz646 -> tx505, tx502 + tx503 -> tx506, tx500 + tx504 -> tx507, 
   tx501 + tx505 -> tx508, tx506*tzz653 -> tx509, tx507*tzz655 -> tx510, 
   tx508*tzz654 -> tx511, tx509 + tx510 + tx511 -> tx512, 
   tx512*tzz657 -> tx513, tx512*tzz652 -> tx514, tx512*tzz656 -> tx515, 
   tx513 + tx514 + tx515 -> gdm, -gdm -> fdm, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I1, 3, fdm );"], -12.*tzz649 -> tzz673, 
   tx399*tzz673 -> tx516, 6.*tzz648 -> tzz669, tx399*tzz669 -> tx517, 
   -tzz647 -> tzz670, tx399*tzz670 -> tx518, tx516 + tx517 + tx518 -> gxm, 
   -gxm -> fxm, CCode["NONBONDRB_FORCE_ACCUMULATE(I1, 4, fxm );"], 
   tx400*tzz673 -> tx519, tx400*tzz669 -> tx520, tx400*tzz670 -> tx521, 
   tx519 + tx520 + tx521 -> gym, -gym -> fym, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I1, 5, fym );"], tx401*tzz673 -> tx522, 
   tx401*tzz669 -> tx523, tx401*tzz670 -> tx524, 
   tx522 + tx523 + tx524 -> gzm, -gzm -> fzm, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I1, 6, fzm );"], tzz674 -> tx525, 
   bn*pyn -> tx526, tzz677 -> tx527, cn*tx344 -> tx528, dn*tx362 -> tx529, 
   bn*tx380 -> tx530, power2[tx343] -> tx531, tx526 + tx528 -> tx532, 
   tx527 + tx529 -> tx533, tx525 + tx530 -> tx534, 4.*tx531 -> tzz650, 
   tx390*tzz650 -> tzz662, an*tzz662 -> tx535, tx391*tzz650 -> tzz661, 
   an*tzz661 -> tx536, tx392*tzz650 -> tzz660, an*tzz660 -> tx537, 
   tx532*tzz645 -> tx538, tx533*tzz645 -> tx539, tx534*tzz645 -> tx540, 
   tx537 + tx538 -> tx541, tx535 + tx539 -> tx542, tx536 + tx540 -> tx543, 
   tx541*tzz653 -> tx544, tx542*tzz655 -> tx545, tx543*tzz654 -> tx546, 
   tx544 + tx545 + tx546 -> tx547, tx547*tzz657 -> tx548, 
   tx547*tzz652 -> tx549, tx547*tzz656 -> tx550, 
   tx548 + tx549 + tx550 -> gan, -gan -> fan, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I2, 0, fan );"], cn*pxn -> tx551, 
   an*pyn -> tx552, tzz678 -> tx553, -2.*tzz667 -> tx554, dn*pzn -> tx555, 
   an*tx380 -> tx556, -2.*tx526 -> tx557, tx525 + tx552 + tx554 -> tx558, 
   tx553 + tx555 -> tx559, tx551 + tx556 + tx557 -> tx560, 
   bn*tzz662 -> tx561, bn*tzz661 -> tx562, bn*tzz660 -> tx563, 
   tx558*tzz645 -> tx564, tx559*tzz645 -> tx565, tx560*tzz645 -> tx566, 
   tx563 + tx564 -> tx567, tx561 + tx565 -> tx568, tx562 + tx566 -> tx569, 
   tx567*tzz653 -> tx570, tx568*tzz655 -> tx571, tx569*tzz654 -> tx572, 
   tx570 + tx571 + tx572 -> tx573, tx573*tzz657 -> tx574, 
   tx573*tzz652 -> tx575, tx573*tzz656 -> tx576, 
   tx574 + tx575 + tx576 -> gbn, -gbn -> fbn, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I2, 1, fbn );"], bn*pxn -> tx577, 
   dn*pyn -> tx578, an*pzn -> tx579, an*tx344 -> tx580, -2.*tx527 -> tx581, 
   -2.*tx551 -> tx582, tx555 + tx577 -> tx583, tx578 + tx580 + tx581 -> 
    tx584, tx526 + tx579 + tx582 -> tx585, cn*tzz662 -> tx586, 
   cn*tzz661 -> tx587, cn*tzz660 -> tx588, tx583*tzz645 -> tx589, 
   tx584*tzz645 -> tx590, tx585*tzz645 -> tx591, tx587 + tx589 -> tx592, 
   tx588 + tx590 -> tx593, tx586 + tx591 -> tx594, tx592*tzz654 -> tx595, 
   tx593*tzz653 -> tx596, tx594*tzz655 -> tx597, 
   tx595 + tx596 + tx597 -> tx598, tx598*tzz657 -> tx599, 
   tx598*tzz652 -> tx600, tx598*tzz656 -> tx601, 
   tx599 + tx600 + tx601 -> gcn, -gcn -> fcn, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I2, 2, fcn );"], an*pxn -> tx602, 
   tzz667 -> tx603, an*tx362 -> tx604, -2.*tx525 -> tx605, 
   -2.*tx578 -> tx606, tx553 + tx577 -> tx607, tx603 + tx604 + tx605 -> 
    tx608, tx527 + tx602 + tx606 -> tx609, dn*tzz662 -> tx610, 
   dn*tzz661 -> tx611, dn*tzz660 -> tx612, tx607*tzz645 -> tx613, 
   tx608*tzz645 -> tx614, tx609*tzz645 -> tx615, tx612 + tx613 -> tx616, 
   tx610 + tx614 -> tx617, tx611 + tx615 -> tx618, tx616*tzz653 -> tx619, 
   tx617*tzz655 -> tx620, tx618*tzz654 -> tx621, 
   tx619 + tx620 + tx621 -> tx622, tx622*tzz657 -> tx623, 
   tx622*tzz652 -> tx624, tx622*tzz656 -> tx625, 
   tx623 + tx624 + tx625 -> gdn, -gdn -> fdn, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I2, 3, fdn );"], 12.*tzz649 -> tzz668, 
   tx399*tzz668 -> tx626, -6.*tzz648 -> tzz672, tx399*tzz672 -> tx627, 
   tx399*tzz647 -> tx628, tx626 + tx627 + tx628 -> gxn, -gxn -> fxn, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I2, 4, fxn );"], tx400*tzz668 -> tx629, 
   tx400*tzz672 -> tx630, tx400*tzz647 -> tx631, 
   tx629 + tx630 + tx631 -> gyn, -gyn -> fyn, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I2, 5, fyn );"], tx401*tzz668 -> tx632, 
   tx401*tzz672 -> tx633, tx401*tzz647 -> tx634, 
   tx632 + tx633 + tx634 -> gzn, -gzn -> fzn, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I2, 6, fzn );"], 
   CCode["} /*calcForce */"], CCode["#endif /* NONBONDRB_CALC_FORCE ]*/"]}, 
 HessianStructure -> HessianStructure, CCode -> "NOT USED", 
 MathCode -> "Block[{DeltaX,DeltaY,DeltaZ,Energy,fam,fan,fbm,fbn,fcm,fcn,fdm,\
fdn,fxm,fxn,fym,fyn,fzm,fzn,gam,gan,gbm,gbn,gcm,gcn,gdm,gdn,gxm,gxn,gym,gyn,g\
zm,gzn,nonbondRBDeviation,tx318,tx319,tx320,tx321,tx322,tx323,tx324,tx325,tx3\
26,tx327,tx328,tx329,tx330,tx331,tx332,tx333,tx334,tx335,tx336,tx337,tx338,tx\
339,tx340,tx341,tx342,tx343,tx344,tx345,tx346,tx347,tx348,tx349,tx350,tx351,t\
x352,tx353,tx354,tx355,tx356,tx357,tx358,tx359,tx360,tx361,tx362,tx363,tx364,\
tx365,tx366,tx367,tx368,tx369,tx370,tx371,tx372,tx373,tx374,tx375,tx376,tx377\
,tx378,tx379,tx380,tx381,tx382,tx383,tx384,tx385,tx386,tx387,tx388,tx389,tx39\
0,tx391,tx392,tx393,tx394,tx395,tx396,tx397,tx398,tx399,tx400,tx401,tx402,tx4\
03,tx404,tx405,tx406,tx407,tx408,tx409,tx410,tx411,tx412,tx413,tx414,tx415,tx\
416,tx417,tx418,tx419,tx420,tx421,tx422,tx423,tx424,tx425,tx426,tx427,tx428,t\
x429,tx430,tx431,tx432,tx433,tx434,tx435,tx436,tx437,tx438,tx439,tx440,tx441,\
tx442,tx443,tx444,tx445,tx446,tx447,tx448,tx449,tx450,tx451,tx452,tx453,tx454\
,tx455,tx456,tx457,tx458,tx459,tx460,tx461,tx462,tx463,tx464,tx465,tx466,tx46\
7,tx468,tx469,tx470,tx471,tx472,tx473,tx474,tx475,tx476,tx477,tx478,tx479,tx4\
80,tx481,tx482,tx483,tx484,tx485,tx486,tx487,tx488,tx489,tx490,tx491,tx492,tx\
493,tx494,tx495,tx496,tx497,tx498,tx499,tx500,tx501,tx502,tx503,tx504,tx505,t\
x506,tx507,tx508,tx509,tx510,tx511,tx512,tx513,tx514,tx515,tx516,tx517,tx518,\
tx519,tx520,tx521,tx522,tx523,tx524,tx525,tx526,tx527,tx528,tx529,tx530,tx531\
,tx532,tx533,tx534,tx535,tx536,tx537,tx538,tx539,tx540,tx541,tx542,tx543,tx54\
4,tx545,tx546,tx547,tx548,tx549,tx550,tx551,tx552,tx553,tx554,tx555,tx556,tx5\
57,tx558,tx559,tx560,tx561,tx562,tx563,tx564,tx565,tx566,tx567,tx568,tx569,tx\
570,tx571,tx572,tx573,tx574,tx575,tx576,tx577,tx578,tx579,tx580,tx581,tx582,t\
x583,tx584,tx585,tx586,tx587,tx588,tx589,tx590,tx591,tx592,tx593,tx594,tx595,\
tx596,tx597,tx598,tx599,tx600,tx601,tx602,tx603,tx604,tx605,tx606,tx607,tx608\
,tx609,tx610,tx611,tx612,tx613,tx614,tx615,tx616,tx617,tx618,tx619,tx620,tx62\
1,tx622,tx623,tx624,tx625,tx626,tx627,tx628,tx629,tx630,tx631,tx632,tx633,tx6\
34,tx635,tx636,tx637,tx638,tx639,tx640,tx641,tx642,tx643,tx644,tx683,tzz645,t\
zz646,tzz647,tzz648,tzz649,tzz650,tzz651,tzz652,tzz653,tzz654,tzz655,tzz656,t\
zz657,tzz658,tzz659,tzz660,tzz661,tzz662,tzz663,tzz664,tzz665,tzz666,tzz667,t\
zz668,tzz669,tzz670,tzz671,tzz672,tzz673,tzz674,tzz675,tzz676,tzz677,tzz678,t\
zz679,tzz680,tzz681,tzz682,tzz684,tzz685,tzz686,tzz687,tzz688,tzz689,tzz690,t\
zz691,tzz692,tzz693,tzz694,tzz695,tzz696,tzz697,tzz698,tzz699,tzz700,tzz701,t\
zz702,tzz703,tzz704,xxxDummy},\n\t(*NONBONDRB_SET_PARAMETER(I1);*)\n\t(*NONBO\
NDRB_SET_PARAMETER(I2);*)\n\t(*NONBONDRB_SET_PARAMETER(dQ1Q2);*)\n\t(*NONBOND\
RB_SET_PARAMETER(dA);*)\n\t(*NONBONDRB_SET_PARAMETER(dC);*)\n\t(*NONBONDRB_SE\
T_POSITION(am,I1,0);*)\n\t(*NONBONDRB_SET_POSITION(bm,I1,1);*)\n\t(*NONBONDRB\
_SET_POSITION(cm,I1,2);*)\n\t(*NONBONDRB_SET_POSITION(dm,I1,3);*)\n\t(*NONBON\
DRB_SET_POSITION(xm,I1,4);*)\n\t(*NONBONDRB_SET_POSITION(ym,I1,5);*)\n\t(*NON\
BONDRB_SET_POSITION(zm,I1,6);*)\n\t(*NONBONDRB_SET_POSITION(an,I2,0);*)\n\t(*\
NONBONDRB_SET_POSITION(bn,I2,1);*)\n\t(*NONBONDRB_SET_POSITION(cn,I2,2);*)\n\
\t(*NONBONDRB_SET_POSITION(dn,I2,3);*)\n\t(*NONBONDRB_SET_POSITION(xn,I2,4);*\
)\n\t(*NONBONDRB_SET_POSITION(yn,I2,5);*)\n\t(*NONBONDRB_SET_POSITION(zn,I2,6\
);*)\n\t(*NONBONDRB_SET_POINT(pxm,ea1,getX());*)\n\t(*NONBONDRB_SET_POINT(pym\
,ea1,getY());*)\n\t(*NONBONDRB_SET_POINT(pzm,ea1,getZ());*)\n\t(*NONBONDRB_SE\
T_POINT(pxn,ea2,getX());*)\n\t(*NONBONDRB_SET_POINT(pyn,ea2,getY());*)\n\t(*N\
ONBONDRB_SET_POINT(pzn,ea2,getZ());*)\n\ttx318= mathPower2[am]; (*Rule \
26*)\n\ttx319= mathPower2[an]; (*Rule 27*)\n\ttx320= mathPower2[bm]; (*Rule \
28*)\n\ttx321= mathPower2[bn]; (*Rule 29*)\n\ttx322= mathPower2[cm]; (*Rule \
30*)\n\ttx323= mathPower2[cn]; (*Rule 31*)\n\ttx324= mathPower2[dm]; (*Rule \
32*)\n\ttx325= mathPower2[dn]; (*Rule 33*)\n\ttzz679=cm pym; (*Rule \
34*)\n\ttx326=bm tzz679; (*Rule 35*)\n\ttzz666=-pym; (*Rule 36*)\n\ttx327=am \
dm tzz666; (*Rule 37*)\n\ttzz678=cn pyn; (*Rule 38*)\n\ttx328=bn tzz678; \
(*Rule 39*)\n\ttzz671=-pyn; (*Rule 40*)\n\ttx329=an dn tzz671; (*Rule \
41*)\n\ttzz680=am cm; (*Rule 42*)\n\ttx330=pzm tzz680; (*Rule \
43*)\n\ttzz675=dm pzm; (*Rule 44*)\n\ttx331=bm tzz675; (*Rule \
45*)\n\ttzz677=cn pzn; (*Rule 46*)\n\ttx332=an tzz677; (*Rule \
47*)\n\ttzz667=bn pzn; (*Rule 48*)\n\ttx333=dn tzz667; (*Rule \
49*)\n\ttzz659=-pxm; (*Rule 50*)\n\ttx334=tx322 tzz659; (*Rule \
51*)\n\ttzz682=-tx323; (*Rule 52*)\n\ttx335=pxn tzz682; (*Rule \
53*)\n\ttx336=tx324 tzz659; (*Rule 54*)\n\ttx337=tx318 + tx320 + tx322 + \
tx324; (*Rule 55*)\n\ttx338=-(pxn tx325); (*Rule 56*)\n\ttx339=tx319 + tx321 \
+ tx323 + tx325; (*Rule 57*)\n\ttx340=tx326 + tx327 + tx330 + tx331 + tx334 + \
tx336; (*Rule 58*)\n\ttx341= mathReciprocal[tx337]; (*Rule \
59*)\n\ttzz699=tx332 + tx333; (*Rule 60*)\n\ttzz700=tx329 + tzz699; (*Rule \
61*)\n\ttzz701=tx328 + tzz700; (*Rule 62*)\n\ttx342=tx335 + tx338 + tzz701; \
(*Rule 63*)\n\ttx343= mathReciprocal[tx339]; (*Rule 64*)\n\ttx344=-pxn; \
(*Rule 65*)\n\ttzz646=2. tx341; (*Rule 66*)\n\ttx345=tx340 tzz646; (*Rule \
67*)\n\ttzz645=-2. tx343; (*Rule 68*)\n\ttx346=tx342 tzz645; (*Rule \
69*)\n\ttx347=-xn; (*Rule 70*)\n\ttzz696=tx347 + xm; (*Rule \
71*)\n\ttzz697=tx345 + tzz696; (*Rule 72*)\n\ttzz698=tx344 + tzz697; (*Rule \
73*)\n\ttzz704=pxm + tzz698; (*Rule 74*)\n\tDeltaX=tx346 + tzz704; (*Rule \
75*)\n\ttx348=bm cm pxm; (*Rule 76*)\n\ttzz676=dm pxm; (*Rule \
77*)\n\ttx349=am tzz676; (*Rule 78*)\n\ttx350=bn cn pxn; (*Rule \
79*)\n\ttzz674=dn pxn; (*Rule 80*)\n\ttx351=an tzz674; (*Rule \
81*)\n\ttzz658=-pzm; (*Rule 82*)\n\ttzz681=am bm; (*Rule 83*)\n\ttx352=tzz658 \
tzz681; (*Rule 84*)\n\ttx353=cm tzz675; (*Rule 85*)\n\ttx354=-(an tzz667); \
(*Rule 86*)\n\ttx355=dn tzz677; (*Rule 87*)\n\ttx356=tx320 tzz666; (*Rule \
88*)\n\ttx357=tx321 tzz671; (*Rule 89*)\n\ttx358=tx324 tzz666; (*Rule \
90*)\n\ttx359=tx325 tzz671; (*Rule 91*)\n\ttx360=tx348 + tx349 + tx352 + \
tx353 + tx356 + tx358; (*Rule 92*)\n\ttzz693=tx354 + tx355; (*Rule \
93*)\n\ttzz694=tx351 + tzz693; (*Rule 94*)\n\ttzz695=tx350 + tzz694; (*Rule \
95*)\n\ttx361=tx357 + tx359 + tzz695; (*Rule 96*)\n\ttx362=tzz671; (*Rule \
97*)\n\ttx363=tx360 tzz646; (*Rule 98*)\n\ttx364=tx361 tzz645; (*Rule \
99*)\n\ttx365=-yn; (*Rule 100*)\n\ttzz690=tx365 + ym; (*Rule \
101*)\n\ttzz691=tx363 + tzz690; (*Rule 102*)\n\ttzz692=tx362 + tzz691; (*Rule \
103*)\n\ttzz703=pym + tzz692; (*Rule 104*)\n\tDeltaY=tx364 + tzz703; (*Rule \
105*)\n\ttx366=tzz659 tzz680; (*Rule 106*)\n\ttx367=bm tzz676; (*Rule \
107*)\n\ttx368=bn tzz674; (*Rule 108*)\n\ttx369=pym tzz681; (*Rule \
109*)\n\ttx370=dm tzz679; (*Rule 110*)\n\ttx371=an bn pyn; (*Rule \
111*)\n\ttx372=dn tzz678; (*Rule 112*)\n\ttx373=tx320 tzz658; (*Rule \
113*)\n\ttx374=-(pzn tx321); (*Rule 114*)\n\ttx375=tx322 tzz658; (*Rule \
115*)\n\ttx376=pzn tzz682; (*Rule 116*)\n\ttx377=an cn tx344; (*Rule \
117*)\n\ttx378=tx366 + tx367 + tx369 + tx370 + tx373 + tx375; (*Rule \
118*)\n\ttzz687=tx372 + tx377; (*Rule 119*)\n\ttzz688=tx371 + tzz687; (*Rule \
120*)\n\ttzz689=tx368 + tzz688; (*Rule 121*)\n\ttx379=tx374 + tx376 + tzz689; \
(*Rule 122*)\n\ttx380=-pzn; (*Rule 123*)\n\ttx381=tx378 tzz646; (*Rule \
124*)\n\ttx382=tx379 tzz645; (*Rule 125*)\n\ttx383=-zn; (*Rule \
126*)\n\ttzz684=tx383 + zm; (*Rule 127*)\n\ttzz685=tx381 + tzz684; (*Rule \
128*)\n\ttzz686=tx380 + tzz685; (*Rule 129*)\n\ttzz702=pzm + tzz686; (*Rule \
130*)\n\tDeltaZ=tx382 + tzz702; (*Rule 131*)\n\ttx384=tx323 tx344; (*Rule \
132*)\n\ttx385=tx325 tx344; (*Rule 133*)\n\ttx386=tx321 tx362; (*Rule \
134*)\n\ttx387=tx325 tx362; (*Rule 135*)\n\ttx388=tx321 tx380; (*Rule \
136*)\n\ttx389=tx323 tx380; (*Rule 137*)\n\ttx390=tx384 + tx385 + tzz701; \
(*Rule 138*)\n\ttx391=tx386 + tx387 + tzz695; (*Rule 139*)\n\ttx392=tx388 + \
tx389 + tzz689; (*Rule 140*)\n\ttx393=tx390 tzz645; (*Rule \
141*)\n\ttx394=tx391 tzz645; (*Rule 142*)\n\ttx395=tx392 tzz645; (*Rule \
143*)\n\ttx396=tx393 + tzz704; (*Rule 144*)\n\ttx397=tx394 + tzz703; (*Rule \
145*)\n\ttx398=tx395 + tzz702; (*Rule 146*)\n\ttx399=PBX[tx396]; (*Rule \
147*)\n\ttx400=PBY[tx397]; (*Rule 148*)\n\ttx401=PBZ[tx398]; (*Rule \
149*)\n\ttx402= mathPower2[tx399]; (*Rule 150*)\n\ttx403= mathPower2[tx400]; \
(*Rule 151*)\n\ttx404= mathPower2[tx401]; (*Rule 152*)\n\ttx405=tx402 + tx403 \
+ tx404; (*Rule 153*)\n\ttx683= mathReciprocal[tx405]; (*Rule 154*)\n\ttx643= \
mathPower2[tx683]; (*Rule 155*)\n\ttx644=tx683; (*Rule 156*)\n\ttx635=tx643 \
tx644; (*Rule 157*)\n\ttx406= mathPower2[tx635]; (*Rule 158*)\n\ttx636= \
mathPower2[tx644]; (*Rule 159*)\n\ttx637=tx644; (*Rule 160*)\n\ttx407=tx636 \
tx637; (*Rule 161*)\n\ttx638= mathSqrt[tx405]; (*Rule 162*)\n\ttx408= \
mathReciprocal[tx638]; (*Rule 163*)\n\ttx409=dA tx406; (*Rule \
164*)\n\ttx410=-(dC tx407); (*Rule 165*)\n\ttx411=dQ1Q2 tx408; (*Rule \
166*)\n\tEnergy=tx409 + tx410 + tx411; (*Rule \
167*)\n\t(*NONBONDRB_ENERGY_ACCUMULATE(Energy);*)\n\t(*#ifdef \
NONBONDRB_CALC_FORCE //[*)\n\t(*if ( calcForce ) \
{*)\n\t(*BAIL_OUT_IF_CUTOFF(DeltaX,DeltaY,DeltaZ)*)\n\ttx412=cm tzz659; \
(*Rule 172*)\n\ttx413=tzz676; (*Rule 173*)\n\ttx414=bm pym; (*Rule \
174*)\n\ttx415=dm tzz666; (*Rule 175*)\n\ttx416=bm tzz658; (*Rule \
176*)\n\ttx417=cm pzm; (*Rule 177*)\n\ttx418= mathPower2[tx341]; (*Rule \
178*)\n\ttx419=tx412 + tx414; (*Rule 179*)\n\ttx420=tx413 + tx416; (*Rule \
180*)\n\ttx421=tx415 + tx417; (*Rule 181*)\n\ttzz651=-4. tx418; (*Rule \
182*)\n\ttzz665=tx340 tzz651; (*Rule 183*)\n\ttx422=am tzz665; (*Rule \
184*)\n\ttzz664=tx360 tzz651; (*Rule 185*)\n\ttx423=am tzz664; (*Rule \
186*)\n\ttzz663=tx378 tzz651; (*Rule 187*)\n\ttx424=am tzz663; (*Rule \
188*)\n\ttx425=tx419 tzz646; (*Rule 189*)\n\ttx426=tx420 tzz646; (*Rule \
190*)\n\ttx427=tx421 tzz646; (*Rule 191*)\n\ttx428=tx424 + tx425; (*Rule \
192*)\n\ttx429=tx423 + tx426; (*Rule 193*)\n\ttx430=tx422 + tx427; (*Rule \
194*)\n\ttzz653=2. tx401; (*Rule 195*)\n\ttx431=tx428 tzz653; (*Rule \
196*)\n\ttzz654=2. tx400; (*Rule 197*)\n\ttx432=tx429 tzz654; (*Rule \
198*)\n\ttzz655=2. tx399; (*Rule 199*)\n\ttx433=tx430 tzz655; (*Rule \
200*)\n\ttx639= mathPower2[tx643]; (*Rule 201*)\n\ttx434=tx635 tx639; (*Rule \
202*)\n\ttx435= mathPower2[tx636]; (*Rule 203*)\n\ttx640=tx405; (*Rule \
204*)\n\ttx641=tx408; (*Rule 205*)\n\ttx642= mathReciprocal[tx640]; (*Rule \
206*)\n\ttx436=tx641 tx642; (*Rule 207*)\n\ttx437=tx431 + tx432 + tx433; \
(*Rule 208*)\n\ttzz649=dA tx434; (*Rule 209*)\n\ttzz657=-6. tzz649; (*Rule \
210*)\n\ttx438=tx437 tzz657; (*Rule 211*)\n\ttzz648=dC tx435; (*Rule \
212*)\n\ttzz652=3 tzz648; (*Rule 213*)\n\ttx439=tx437 tzz652; (*Rule \
214*)\n\ttzz647=dQ1Q2 tx436; (*Rule 215*)\n\ttzz656=-0.5 tzz647; (*Rule \
216*)\n\ttx440=tx437 tzz656; (*Rule 217*)\n\tgam=tx438 + tx439 + tx440; \
(*Rule 218*)\n\tfam=-gam; (*Rule 219*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I1, 0, \
fam );*)\n\ttx441=cm pxm; (*Rule 221*)\n\ttx442=am pym; (*Rule \
222*)\n\ttx443=tzz679; (*Rule 223*)\n\ttx444=am tzz658; (*Rule \
224*)\n\ttx445=-2. bm pzm; (*Rule 225*)\n\ttx446=tzz675; (*Rule \
226*)\n\ttx447=-2. tx414; (*Rule 227*)\n\ttx448=tx413 + tx442 + tx445; (*Rule \
228*)\n\ttx449=tx443 + tx446; (*Rule 229*)\n\ttx450=tx441 + tx444 + tx447; \
(*Rule 230*)\n\ttx451=bm tzz665; (*Rule 231*)\n\ttx452=bm tzz664; (*Rule \
232*)\n\ttx453=bm tzz663; (*Rule 233*)\n\ttx454=tx448 tzz646; (*Rule \
234*)\n\ttx455=tx449 tzz646; (*Rule 235*)\n\ttx456=tx450 tzz646; (*Rule \
236*)\n\ttx457=tx453 + tx454; (*Rule 237*)\n\ttx458=tx451 + tx455; (*Rule \
238*)\n\ttx459=tx452 + tx456; (*Rule 239*)\n\ttx460=tx457 tzz653; (*Rule \
240*)\n\ttx461=tx458 tzz655; (*Rule 241*)\n\ttx462=tx459 tzz654; (*Rule \
242*)\n\ttx463=tx460 + tx461 + tx462; (*Rule 243*)\n\ttx464=tx463 tzz657; \
(*Rule 244*)\n\ttx465=tx463 tzz652; (*Rule 245*)\n\ttx466=tx463 tzz656; \
(*Rule 246*)\n\tgbm=tx464 + tx465 + tx466; (*Rule 247*)\n\tfbm=-gbm; (*Rule \
248*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I1, 1, fbm );*)\n\ttx467=am tzz659; \
(*Rule 250*)\n\ttx468=bm pxm; (*Rule 251*)\n\ttx469=dm pym; (*Rule \
252*)\n\ttx470=am pzm; (*Rule 253*)\n\ttx471=-2. tx417; (*Rule \
254*)\n\ttx472=-2. tx441; (*Rule 255*)\n\ttx473=tx446 + tx468; (*Rule \
256*)\n\ttx474=tx467 + tx469 + tx471; (*Rule 257*)\n\ttx475=tx414 + tx470 + \
tx472; (*Rule 258*)\n\ttx476=cm tzz665; (*Rule 259*)\n\ttx477=cm tzz664; \
(*Rule 260*)\n\ttx478=cm tzz663; (*Rule 261*)\n\ttx479=tx473 tzz646; (*Rule \
262*)\n\ttx480=tx474 tzz646; (*Rule 263*)\n\ttx481=tx475 tzz646; (*Rule \
264*)\n\ttx482=tx477 + tx479; (*Rule 265*)\n\ttx483=tx478 + tx480; (*Rule \
266*)\n\ttx484=tx476 + tx481; (*Rule 267*)\n\ttx485=tx482 tzz654; (*Rule \
268*)\n\ttx486=tx483 tzz653; (*Rule 269*)\n\ttx487=tx484 tzz655; (*Rule \
270*)\n\ttx488=tx485 + tx486 + tx487; (*Rule 271*)\n\ttx489=tx488 tzz657; \
(*Rule 272*)\n\ttx490=tx488 tzz652; (*Rule 273*)\n\ttx491=tx488 tzz656; \
(*Rule 274*)\n\tgcm=tx489 + tx490 + tx491; (*Rule 275*)\n\tfcm=-gcm; (*Rule \
276*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I1, 2, fcm );*)\n\ttx492=am pxm; (*Rule \
278*)\n\ttx493=bm pzm; (*Rule 279*)\n\ttx494=-2. tx413; (*Rule \
280*)\n\ttx495=-tx442; (*Rule 281*)\n\ttx496=-2. tx469; (*Rule \
282*)\n\ttx497=tx443 + tx468; (*Rule 283*)\n\ttx498=tx493 + tx494 + tx495; \
(*Rule 284*)\n\ttx499=tx417 + tx492 + tx496; (*Rule 285*)\n\ttx500=dm tzz665; \
(*Rule 286*)\n\ttx501=dm tzz664; (*Rule 287*)\n\ttx502=dm tzz663; (*Rule \
288*)\n\ttx503=tx497 tzz646; (*Rule 289*)\n\ttx504=tx498 tzz646; (*Rule \
290*)\n\ttx505=tx499 tzz646; (*Rule 291*)\n\ttx506=tx502 + tx503; (*Rule \
292*)\n\ttx507=tx500 + tx504; (*Rule 293*)\n\ttx508=tx501 + tx505; (*Rule \
294*)\n\ttx509=tx506 tzz653; (*Rule 295*)\n\ttx510=tx507 tzz655; (*Rule \
296*)\n\ttx511=tx508 tzz654; (*Rule 297*)\n\ttx512=tx509 + tx510 + tx511; \
(*Rule 298*)\n\ttx513=tx512 tzz657; (*Rule 299*)\n\ttx514=tx512 tzz652; \
(*Rule 300*)\n\ttx515=tx512 tzz656; (*Rule 301*)\n\tgdm=tx513 + tx514 + \
tx515; (*Rule 302*)\n\tfdm=-gdm; (*Rule \
303*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I1, 3, fdm );*)\n\ttzz673=-12. tzz649; \
(*Rule 305*)\n\ttx516=tx399 tzz673; (*Rule 306*)\n\ttzz669=6. tzz648; (*Rule \
307*)\n\ttx517=tx399 tzz669; (*Rule 308*)\n\ttzz670=-tzz647; (*Rule \
309*)\n\ttx518=tx399 tzz670; (*Rule 310*)\n\tgxm=tx516 + tx517 + tx518; \
(*Rule 311*)\n\tfxm=-gxm; (*Rule 312*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I1, 4, \
fxm );*)\n\ttx519=tx400 tzz673; (*Rule 314*)\n\ttx520=tx400 tzz669; (*Rule \
315*)\n\ttx521=tx400 tzz670; (*Rule 316*)\n\tgym=tx519 + tx520 + tx521; \
(*Rule 317*)\n\tfym=-gym; (*Rule 318*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I1, 5, \
fym );*)\n\ttx522=tx401 tzz673; (*Rule 320*)\n\ttx523=tx401 tzz669; (*Rule \
321*)\n\ttx524=tx401 tzz670; (*Rule 322*)\n\tgzm=tx522 + tx523 + tx524; \
(*Rule 323*)\n\tfzm=-gzm; (*Rule 324*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I1, 6, \
fzm );*)\n\ttx525=tzz674; (*Rule 326*)\n\ttx526=bn pyn; (*Rule \
327*)\n\ttx527=tzz677; (*Rule 328*)\n\ttx528=cn tx344; (*Rule \
329*)\n\ttx529=dn tx362; (*Rule 330*)\n\ttx530=bn tx380; (*Rule \
331*)\n\ttx531= mathPower2[tx343]; (*Rule 332*)\n\ttx532=tx526 + tx528; \
(*Rule 333*)\n\ttx533=tx527 + tx529; (*Rule 334*)\n\ttx534=tx525 + tx530; \
(*Rule 335*)\n\ttzz650=4. tx531; (*Rule 336*)\n\ttzz662=tx390 tzz650; (*Rule \
337*)\n\ttx535=an tzz662; (*Rule 338*)\n\ttzz661=tx391 tzz650; (*Rule \
339*)\n\ttx536=an tzz661; (*Rule 340*)\n\ttzz660=tx392 tzz650; (*Rule \
341*)\n\ttx537=an tzz660; (*Rule 342*)\n\ttx538=tx532 tzz645; (*Rule \
343*)\n\ttx539=tx533 tzz645; (*Rule 344*)\n\ttx540=tx534 tzz645; (*Rule \
345*)\n\ttx541=tx537 + tx538; (*Rule 346*)\n\ttx542=tx535 + tx539; (*Rule \
347*)\n\ttx543=tx536 + tx540; (*Rule 348*)\n\ttx544=tx541 tzz653; (*Rule \
349*)\n\ttx545=tx542 tzz655; (*Rule 350*)\n\ttx546=tx543 tzz654; (*Rule \
351*)\n\ttx547=tx544 + tx545 + tx546; (*Rule 352*)\n\ttx548=tx547 tzz657; \
(*Rule 353*)\n\ttx549=tx547 tzz652; (*Rule 354*)\n\ttx550=tx547 tzz656; \
(*Rule 355*)\n\tgan=tx548 + tx549 + tx550; (*Rule 356*)\n\tfan=-gan; (*Rule \
357*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I2, 0, fan );*)\n\ttx551=cn pxn; (*Rule \
359*)\n\ttx552=an pyn; (*Rule 360*)\n\ttx553=tzz678; (*Rule \
361*)\n\ttx554=-2. tzz667; (*Rule 362*)\n\ttx555=dn pzn; (*Rule \
363*)\n\ttx556=an tx380; (*Rule 364*)\n\ttx557=-2. tx526; (*Rule \
365*)\n\ttx558=tx525 + tx552 + tx554; (*Rule 366*)\n\ttx559=tx553 + tx555; \
(*Rule 367*)\n\ttx560=tx551 + tx556 + tx557; (*Rule 368*)\n\ttx561=bn tzz662; \
(*Rule 369*)\n\ttx562=bn tzz661; (*Rule 370*)\n\ttx563=bn tzz660; (*Rule \
371*)\n\ttx564=tx558 tzz645; (*Rule 372*)\n\ttx565=tx559 tzz645; (*Rule \
373*)\n\ttx566=tx560 tzz645; (*Rule 374*)\n\ttx567=tx563 + tx564; (*Rule \
375*)\n\ttx568=tx561 + tx565; (*Rule 376*)\n\ttx569=tx562 + tx566; (*Rule \
377*)\n\ttx570=tx567 tzz653; (*Rule 378*)\n\ttx571=tx568 tzz655; (*Rule \
379*)\n\ttx572=tx569 tzz654; (*Rule 380*)\n\ttx573=tx570 + tx571 + tx572; \
(*Rule 381*)\n\ttx574=tx573 tzz657; (*Rule 382*)\n\ttx575=tx573 tzz652; \
(*Rule 383*)\n\ttx576=tx573 tzz656; (*Rule 384*)\n\tgbn=tx574 + tx575 + \
tx576; (*Rule 385*)\n\tfbn=-gbn; (*Rule \
386*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I2, 1, fbn );*)\n\ttx577=bn pxn; (*Rule \
388*)\n\ttx578=dn pyn; (*Rule 389*)\n\ttx579=an pzn; (*Rule 390*)\n\ttx580=an \
tx344; (*Rule 391*)\n\ttx581=-2. tx527; (*Rule 392*)\n\ttx582=-2. tx551; \
(*Rule 393*)\n\ttx583=tx555 + tx577; (*Rule 394*)\n\ttx584=tx578 + tx580 + \
tx581; (*Rule 395*)\n\ttx585=tx526 + tx579 + tx582; (*Rule 396*)\n\ttx586=cn \
tzz662; (*Rule 397*)\n\ttx587=cn tzz661; (*Rule 398*)\n\ttx588=cn tzz660; \
(*Rule 399*)\n\ttx589=tx583 tzz645; (*Rule 400*)\n\ttx590=tx584 tzz645; \
(*Rule 401*)\n\ttx591=tx585 tzz645; (*Rule 402*)\n\ttx592=tx587 + tx589; \
(*Rule 403*)\n\ttx593=tx588 + tx590; (*Rule 404*)\n\ttx594=tx586 + tx591; \
(*Rule 405*)\n\ttx595=tx592 tzz654; (*Rule 406*)\n\ttx596=tx593 tzz653; \
(*Rule 407*)\n\ttx597=tx594 tzz655; (*Rule 408*)\n\ttx598=tx595 + tx596 + \
tx597; (*Rule 409*)\n\ttx599=tx598 tzz657; (*Rule 410*)\n\ttx600=tx598 \
tzz652; (*Rule 411*)\n\ttx601=tx598 tzz656; (*Rule 412*)\n\tgcn=tx599 + tx600 \
+ tx601; (*Rule 413*)\n\tfcn=-gcn; (*Rule \
414*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I2, 2, fcn );*)\n\ttx602=an pxn; (*Rule \
416*)\n\ttx603=tzz667; (*Rule 417*)\n\ttx604=an tx362; (*Rule \
418*)\n\ttx605=-2. tx525; (*Rule 419*)\n\ttx606=-2. tx578; (*Rule \
420*)\n\ttx607=tx553 + tx577; (*Rule 421*)\n\ttx608=tx603 + tx604 + tx605; \
(*Rule 422*)\n\ttx609=tx527 + tx602 + tx606; (*Rule 423*)\n\ttx610=dn tzz662; \
(*Rule 424*)\n\ttx611=dn tzz661; (*Rule 425*)\n\ttx612=dn tzz660; (*Rule \
426*)\n\ttx613=tx607 tzz645; (*Rule 427*)\n\ttx614=tx608 tzz645; (*Rule \
428*)\n\ttx615=tx609 tzz645; (*Rule 429*)\n\ttx616=tx612 + tx613; (*Rule \
430*)\n\ttx617=tx610 + tx614; (*Rule 431*)\n\ttx618=tx611 + tx615; (*Rule \
432*)\n\ttx619=tx616 tzz653; (*Rule 433*)\n\ttx620=tx617 tzz655; (*Rule \
434*)\n\ttx621=tx618 tzz654; (*Rule 435*)\n\ttx622=tx619 + tx620 + tx621; \
(*Rule 436*)\n\ttx623=tx622 tzz657; (*Rule 437*)\n\ttx624=tx622 tzz652; \
(*Rule 438*)\n\ttx625=tx622 tzz656; (*Rule 439*)\n\tgdn=tx623 + tx624 + \
tx625; (*Rule 440*)\n\tfdn=-gdn; (*Rule \
441*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I2, 3, fdn );*)\n\ttzz668=12. tzz649; \
(*Rule 443*)\n\ttx626=tx399 tzz668; (*Rule 444*)\n\ttzz672=-6. tzz648; (*Rule \
445*)\n\ttx627=tx399 tzz672; (*Rule 446*)\n\ttx628=tx399 tzz647; (*Rule \
447*)\n\tgxn=tx626 + tx627 + tx628; (*Rule 448*)\n\tfxn=-gxn; (*Rule \
449*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I2, 4, fxn );*)\n\ttx629=tx400 tzz668; \
(*Rule 451*)\n\ttx630=tx400 tzz672; (*Rule 452*)\n\ttx631=tx400 tzz647; \
(*Rule 453*)\n\tgyn=tx629 + tx630 + tx631; (*Rule 454*)\n\tfyn=-gyn; (*Rule \
455*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I2, 5, fyn );*)\n\ttx632=tx401 tzz668; \
(*Rule 457*)\n\ttx633=tx401 tzz672; (*Rule 458*)\n\ttx634=tx401 tzz647; \
(*Rule 459*)\n\tgzn=tx632 + tx633 + tx634; (*Rule 460*)\n\tfzn=-gzn; (*Rule \
461*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I2, 6, fzn );*)\n\t(*} /*calcForce \
*/*)\n\t(*#endif /* NONBONDRB_CALC_FORCE ]*/*)\n{Energy, DeltaX, DeltaY, \
DeltaZ, fam, fbm, fcm, fdm, fxm, fym, fzm, fan, fbn, fcn, fdn, fxn, fyn, fzn, \
nonbondRBDeviation}]\n"}
