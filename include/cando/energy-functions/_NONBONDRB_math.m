{Name -> "NONBONDRB", AdditionalCDeclares -> "", 
 Input -> {ak, bk, ck, dk, xk, yk, zk, al, bl, cl, dl, xl, yl, zl}, 
 Output -> {Energy, EnergyVdw, EnergyElectrostatic, fam, fbm, fcm, fdm, fxm, 
   fym, fzm, fan, fbn, fcn, fdn, fxn, fyn, fzn, nonbondRBDeviation}, 
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
   CCode["NONBONDRB_SET_POINT(pzn,ea2,getZ());"], power2[am] -> tx325, 
   power2[an] -> tx326, power2[bm] -> tx327, power2[bn] -> tx328, 
   power2[cm] -> tx329, power2[cn] -> tx330, power2[dm] -> tx331, 
   power2[dn] -> tx332, -pxm -> tzz678, am*cm -> tzz700, 
   tzz678*tzz700 -> tx333, bm*cm*pxm -> tx334, dm*pxm -> tzz696, 
   am*tzz696 -> tx335, bm*tzz696 -> tx336, -pxn -> tzz691, 
   an*cn*tzz691 -> tx337, bn*cn*pxn -> tx338, dn*pxn -> tzz694, 
   an*tzz694 -> tx339, bn*tzz694 -> tx340, am*bm -> tzz701, 
   pym*tzz701 -> tx341, cm*pym -> tzz699, bm*tzz699 -> tx342, -pym -> tzz685, 
   am*dm*tzz685 -> tx343, dm*tzz699 -> tx344, an*bn*pyn -> tx345, 
   cn*pyn -> tzz698, bn*tzz698 -> tx346, -pyn -> tzz690, 
   an*dn*tzz690 -> tx347, dn*tzz698 -> tx348, -pzm -> tzz677, 
   tzz677*tzz701 -> tx349, pzm*tzz700 -> tx350, dm*pzm -> tzz695, 
   bm*tzz695 -> tx351, cm*tzz695 -> tx352, bn*pzn -> tzz686, 
   -(an*tzz686) -> tx353, cn*pzn -> tzz697, an*tzz697 -> tx354, 
   dn*tzz686 -> tx355, dn*tzz697 -> tx356, tx327*tzz685 -> tx357, 
   tx327*tzz677 -> tx358, tx328*tzz690 -> tx359, -pzn -> tzz702, 
   tx328*tzz702 -> tx360, tx329*tzz678 -> tx361, tx329*tzz677 -> tx362, 
   tx330*tzz691 -> tx363, tx330*tzz702 -> tx364, tx331*tzz678 -> tx365, 
   tx331*tzz685 -> tx366, tx325 + tx327 + tx329 + tx331 -> tx367, 
   tx332*tzz691 -> tx368, tx332*tzz690 -> tx369, 
   tx326 + tx328 + tx330 + tx332 -> tx370, 
   tx333 + tx336 + tx341 + tx344 + tx358 + tx362 -> tx371, 
   tx345 + tx348 -> tzz719, tx340 + tzz719 -> tzz720, 
   tx337 + tzz720 -> tzz723, tx360 + tx364 + tzz723 -> tx372, 
   tx342 + tx343 + tx350 + tx351 + tx361 + tx365 -> tx373, 
   tx334 + tx335 + tx349 + tx352 + tx357 + tx366 -> tx374, 
   reciprocal[tx367] -> tx375, tx354 + tx355 -> tzz715, 
   tx347 + tzz715 -> tzz717, tx346 + tzz717 -> tzz718, 
   tx363 + tx368 + tzz718 -> tx376, tx353 + tx356 -> tzz716, 
   tx339 + tzz716 -> tzz721, tx338 + tzz721 -> tzz722, 
   tx359 + tx369 + tzz722 -> tx377, reciprocal[tx370] -> tx378, 
   tzz691 -> tx379, tzz690 -> tx380, tzz702 -> tx381, 2.*tx375 -> tzz665, 
   tx371*tzz665 -> tx382, tx373*tzz665 -> tx383, tx374*tzz665 -> tx384, 
   -2.*tx378 -> tzz664, tx372*tzz664 -> tx385, tx376*tzz664 -> tx386, 
   tx377*tzz664 -> tx387, -xn -> tx388, -yn -> tx389, -zn -> tx390, 
   tx388 + xm -> tzz708, tx383 + tzz708 -> tzz710, tx379 + tzz710 -> tzz714, 
   pxm + tzz714 -> tzz726, tx386 + tzz726 -> tx391, tx389 + ym -> tzz707, 
   tx384 + tzz707 -> tzz709, tx380 + tzz709 -> tzz713, 
   pym + tzz713 -> tzz725, tx387 + tzz725 -> tx392, tx390 + zm -> tzz706, 
   tx382 + tzz706 -> tzz711, tx381 + tzz711 -> tzz712, 
   pzm + tzz712 -> tzz724, tx385 + tzz724 -> tx393, power2[tx391] -> tx394, 
   power2[tx392] -> tx395, power2[tx393] -> tx396, 
   tx394 + tx395 + tx396 -> tx397, reciprocal[tx397] -> tx704, 
   power2[tx704] -> tx660, tx704 -> tx661, tx660*tx661 -> tx649, 
   power2[tx649] -> tx398, power2[tx661] -> tx650, tx661 -> tx651, 
   tx650*tx651 -> tx399, dA*tx398 -> tx400, -dC -> tzz703, 
   tx399*tzz703 -> tx401, tx400 + tx401 -> EnergyVdw, tx330*tx379 -> tx402, 
   tx332*tx379 -> tx403, tx328*tx380 -> tx404, tx332*tx380 -> tx405, 
   tx328*tx381 -> tx406, tx330*tx381 -> tx407, tx402 + tx403 + tzz718 -> 
    tx408, tx404 + tx405 + tzz722 -> tx409, tx406 + tx407 + tzz723 -> tx410, 
   tx408*tzz664 -> tx411, tx409*tzz664 -> tx412, tx410*tzz664 -> tx413, 
   tx411 + tzz726 -> tx414, tx412 + tzz725 -> tx415, tx413 + tzz724 -> tx416, 
   power2[tx414] -> tx417, power2[tx415] -> tx418, power2[tx416] -> tx419, 
   tx417 + tx418 + tx419 -> tx420, mysqrt[tx420] -> tx652, 
   reciprocal[tx652] -> tx421, dQ1Q2*tx421 -> EnergyElectrostatic, 
   reciprocal[tx420] -> tx705, power2[tx705] -> tx662, tx705 -> tx663, 
   tx662*tx663 -> tx653, power2[tx653] -> tx422, power2[tx663] -> tx654, 
   tx663 -> tx655, tx654*tx655 -> tx423, dA*tx422 -> tx424, 
   tx423*tzz703 -> tx425, EnergyElectrostatic + tx424 + tx425 -> Energy, 
   CCode["NONBONDRB_EVDW_ENERGY_ACCUMULATE(EnergyVdw);"], 
   CCode["NONBONDRB_EEEL_ENERGY_ACCUMULATE(EnergyElectrostatic);"], 
   CCode["NONBONDRB_ENERGY_ACCUMULATE(Energy);"], 
   CCode["#ifdef NONBONDRB_CALC_FORCE //["], CCode["if ( calcForce ) {"], 
   cm*tzz678 -> tx426, tzz696 -> tx427, bm*pym -> tx428, dm*tzz685 -> tx429, 
   bm*tzz677 -> tx430, cm*pzm -> tx431, power2[tx375] -> tx432, 
   tx426 + tx428 -> tx433, tx427 + tx430 -> tx434, tx429 + tx431 -> tx435, 
   -4.*tx432 -> tzz670, tx371*tzz670 -> tzz684, am*tzz684 -> tx436, 
   tx373*tzz670 -> tzz683, am*tzz683 -> tx437, tx374*tzz670 -> tzz682, 
   am*tzz682 -> tx438, tx433*tzz665 -> tx439, tx434*tzz665 -> tx440, 
   tx435*tzz665 -> tx441, tx436 + tx439 -> tx442, tx438 + tx440 -> tx443, 
   tx437 + tx441 -> tx444, 2.*tx416 -> tzz672, tx442*tzz672 -> tx445, 
   2.*tx415 -> tzz673, tx443*tzz673 -> tx446, 2.*tx414 -> tzz674, 
   tx444*tzz674 -> tx447, power2[tx662] -> tx656, tx653*tx656 -> tx448, 
   power2[tx654] -> tx449, tx420 -> tx657, tx421 -> tx658, 
   reciprocal[tx657] -> tx659, tx658*tx659 -> tx450, 
   tx445 + tx446 + tx447 -> tx451, dA*tx448 -> tzz668, -6.*tzz668 -> tzz676, 
   tx451*tzz676 -> tx452, dC*tx449 -> tzz667, 3*tzz667 -> tzz671, 
   tx451*tzz671 -> tx453, dQ1Q2*tx450 -> tzz666, -0.5*tzz666 -> tzz675, 
   tx451*tzz675 -> tx454, tx452 + tx453 + tx454 -> gam, -gam -> fam, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I1, 0, fam );"], cm*pxm -> tx455, 
   am*pym -> tx456, tzz699 -> tx457, am*tzz677 -> tx458, -2.*bm*pzm -> tx459, 
   tzz695 -> tx460, -2.*tx428 -> tx461, tx427 + tx456 + tx459 -> tx462, 
   tx457 + tx460 -> tx463, tx455 + tx458 + tx461 -> tx464, 
   bm*tzz684 -> tx465, bm*tzz683 -> tx466, bm*tzz682 -> tx467, 
   tx462*tzz665 -> tx468, tx463*tzz665 -> tx469, tx464*tzz665 -> tx470, 
   tx465 + tx468 -> tx471, tx466 + tx469 -> tx472, tx467 + tx470 -> tx473, 
   tx471*tzz672 -> tx474, tx472*tzz674 -> tx475, tx473*tzz673 -> tx476, 
   tx474 + tx475 + tx476 -> tx477, tx477*tzz676 -> tx478, 
   tx477*tzz671 -> tx479, tx477*tzz675 -> tx480, 
   tx478 + tx479 + tx480 -> gbm, -gbm -> fbm, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I1, 1, fbm );"], am*tzz678 -> tx481, 
   bm*pxm -> tx482, dm*pym -> tx483, am*pzm -> tx484, -2.*tx431 -> tx485, 
   -2.*tx455 -> tx486, tx460 + tx482 -> tx487, tx481 + tx483 + tx485 -> 
    tx488, tx428 + tx484 + tx486 -> tx489, cm*tzz684 -> tx490, 
   cm*tzz683 -> tx491, cm*tzz682 -> tx492, tx487*tzz665 -> tx493, 
   tx488*tzz665 -> tx494, tx489*tzz665 -> tx495, tx492 + tx493 -> tx496, 
   tx490 + tx494 -> tx497, tx491 + tx495 -> tx498, tx496*tzz673 -> tx499, 
   tx497*tzz672 -> tx500, tx498*tzz674 -> tx501, 
   tx499 + tx500 + tx501 -> tx502, tx502*tzz676 -> tx503, 
   tx502*tzz671 -> tx504, tx502*tzz675 -> tx505, 
   tx503 + tx504 + tx505 -> gcm, -gcm -> fcm, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I1, 2, fcm );"], am*pxm -> tx506, 
   bm*pzm -> tx507, -2.*tx427 -> tx508, -tx456 -> tx509, -2.*tx483 -> tx510, 
   tx457 + tx482 -> tx511, tx507 + tx508 + tx509 -> tx512, 
   tx431 + tx506 + tx510 -> tx513, dm*tzz684 -> tx514, dm*tzz683 -> tx515, 
   dm*tzz682 -> tx516, tx511*tzz665 -> tx517, tx512*tzz665 -> tx518, 
   tx513*tzz665 -> tx519, tx514 + tx517 -> tx520, tx515 + tx518 -> tx521, 
   tx516 + tx519 -> tx522, tx520*tzz672 -> tx523, tx521*tzz674 -> tx524, 
   tx522*tzz673 -> tx525, tx523 + tx524 + tx525 -> tx526, 
   tx526*tzz676 -> tx527, tx526*tzz671 -> tx528, tx526*tzz675 -> tx529, 
   tx527 + tx528 + tx529 -> gdm, -gdm -> fdm, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I1, 3, fdm );"], -12.*tzz668 -> tzz693, 
   tx414*tzz693 -> tx530, 6.*tzz667 -> tzz688, tx414*tzz688 -> tx531, 
   -tzz666 -> tzz689, tx414*tzz689 -> tx532, tx530 + tx531 + tx532 -> gxm, 
   -gxm -> fxm, CCode["NONBONDRB_FORCE_ACCUMULATE(I1, 4, fxm );"], 
   tx415*tzz693 -> tx533, tx415*tzz688 -> tx534, tx415*tzz689 -> tx535, 
   tx533 + tx534 + tx535 -> gym, -gym -> fym, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I1, 5, fym );"], tx416*tzz693 -> tx536, 
   tx416*tzz688 -> tx537, tx416*tzz689 -> tx538, 
   tx536 + tx537 + tx538 -> gzm, -gzm -> fzm, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I1, 6, fzm );"], tzz694 -> tx539, 
   bn*pyn -> tx540, tzz697 -> tx541, cn*tx379 -> tx542, dn*tx380 -> tx543, 
   bn*tx381 -> tx544, power2[tx378] -> tx545, tx540 + tx542 -> tx546, 
   tx541 + tx543 -> tx547, tx539 + tx544 -> tx548, 4.*tx545 -> tzz669, 
   tx408*tzz669 -> tzz681, an*tzz681 -> tx549, tx409*tzz669 -> tzz680, 
   an*tzz680 -> tx550, tx410*tzz669 -> tzz679, an*tzz679 -> tx551, 
   tx546*tzz664 -> tx552, tx547*tzz664 -> tx553, tx548*tzz664 -> tx554, 
   tx551 + tx552 -> tx555, tx549 + tx553 -> tx556, tx550 + tx554 -> tx557, 
   tx555*tzz672 -> tx558, tx556*tzz674 -> tx559, tx557*tzz673 -> tx560, 
   tx558 + tx559 + tx560 -> tx561, tx561*tzz676 -> tx562, 
   tx561*tzz671 -> tx563, tx561*tzz675 -> tx564, 
   tx562 + tx563 + tx564 -> gan, -gan -> fan, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I2, 0, fan );"], cn*pxn -> tx565, 
   an*pyn -> tx566, tzz698 -> tx567, -2.*tzz686 -> tx568, dn*pzn -> tx569, 
   an*tx381 -> tx570, -2.*tx540 -> tx571, tx539 + tx566 + tx568 -> tx572, 
   tx567 + tx569 -> tx573, tx565 + tx570 + tx571 -> tx574, 
   bn*tzz681 -> tx575, bn*tzz680 -> tx576, bn*tzz679 -> tx577, 
   tx572*tzz664 -> tx578, tx573*tzz664 -> tx579, tx574*tzz664 -> tx580, 
   tx577 + tx578 -> tx581, tx575 + tx579 -> tx582, tx576 + tx580 -> tx583, 
   tx581*tzz672 -> tx584, tx582*tzz674 -> tx585, tx583*tzz673 -> tx586, 
   tx584 + tx585 + tx586 -> tx587, tx587*tzz676 -> tx588, 
   tx587*tzz671 -> tx589, tx587*tzz675 -> tx590, 
   tx588 + tx589 + tx590 -> gbn, -gbn -> fbn, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I2, 1, fbn );"], bn*pxn -> tx591, 
   dn*pyn -> tx592, an*pzn -> tx593, an*tx379 -> tx594, -2.*tx541 -> tx595, 
   -2.*tx565 -> tx596, tx569 + tx591 -> tx597, tx592 + tx594 + tx595 -> 
    tx598, tx540 + tx593 + tx596 -> tx599, cn*tzz681 -> tx600, 
   cn*tzz680 -> tx601, cn*tzz679 -> tx602, tx597*tzz664 -> tx603, 
   tx598*tzz664 -> tx604, tx599*tzz664 -> tx605, tx601 + tx603 -> tx606, 
   tx602 + tx604 -> tx607, tx600 + tx605 -> tx608, tx606*tzz673 -> tx609, 
   tx607*tzz672 -> tx610, tx608*tzz674 -> tx611, 
   tx609 + tx610 + tx611 -> tx612, tx612*tzz676 -> tx613, 
   tx612*tzz671 -> tx614, tx612*tzz675 -> tx615, 
   tx613 + tx614 + tx615 -> gcn, -gcn -> fcn, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I2, 2, fcn );"], an*pxn -> tx616, 
   tzz686 -> tx617, an*tx380 -> tx618, -2.*tx539 -> tx619, 
   -2.*tx592 -> tx620, tx567 + tx591 -> tx621, tx617 + tx618 + tx619 -> 
    tx622, tx541 + tx616 + tx620 -> tx623, dn*tzz681 -> tx624, 
   dn*tzz680 -> tx625, dn*tzz679 -> tx626, tx621*tzz664 -> tx627, 
   tx622*tzz664 -> tx628, tx623*tzz664 -> tx629, tx626 + tx627 -> tx630, 
   tx624 + tx628 -> tx631, tx625 + tx629 -> tx632, tx630*tzz672 -> tx633, 
   tx631*tzz674 -> tx634, tx632*tzz673 -> tx635, 
   tx633 + tx634 + tx635 -> tx636, tx636*tzz676 -> tx637, 
   tx636*tzz671 -> tx638, tx636*tzz675 -> tx639, 
   tx637 + tx638 + tx639 -> gdn, -gdn -> fdn, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I2, 3, fdn );"], 12.*tzz668 -> tzz687, 
   tx414*tzz687 -> tx640, -6.*tzz667 -> tzz692, tx414*tzz692 -> tx641, 
   tx414*tzz666 -> tx642, tx640 + tx641 + tx642 -> gxn, -gxn -> fxn, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I2, 4, fxn );"], tx415*tzz687 -> tx643, 
   tx415*tzz692 -> tx644, tx415*tzz666 -> tx645, 
   tx643 + tx644 + tx645 -> gyn, -gyn -> fyn, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I2, 5, fyn );"], tx416*tzz687 -> tx646, 
   tx416*tzz692 -> tx647, tx416*tzz666 -> tx648, 
   tx646 + tx647 + tx648 -> gzn, -gzn -> fzn, 
   CCode["NONBONDRB_FORCE_ACCUMULATE(I2, 6, fzn );"], 
   CCode["} /*calcForce */"], CCode["#endif /* NONBONDRB_CALC_FORCE ]*/"]}, 
 HessianStructure -> HessianStructure, CCode -> "NOT USED", 
 MathCode -> "Block[{Energy,EnergyElectrostatic,EnergyVdw,fam,fan,fbm,fbn,fcm\
,fcn,fdm,fdn,fxm,fxn,fym,fyn,fzm,fzn,gam,gan,gbm,gbn,gcm,gcn,gdm,gdn,gxm,gxn,\
gym,gyn,gzm,gzn,nonbondRBDeviation,tx325,tx326,tx327,tx328,tx329,tx330,tx331,\
tx332,tx333,tx334,tx335,tx336,tx337,tx338,tx339,tx340,tx341,tx342,tx343,tx344\
,tx345,tx346,tx347,tx348,tx349,tx350,tx351,tx352,tx353,tx354,tx355,tx356,tx35\
7,tx358,tx359,tx360,tx361,tx362,tx363,tx364,tx365,tx366,tx367,tx368,tx369,tx3\
70,tx371,tx372,tx373,tx374,tx375,tx376,tx377,tx378,tx379,tx380,tx381,tx382,tx\
383,tx384,tx385,tx386,tx387,tx388,tx389,tx390,tx391,tx392,tx393,tx394,tx395,t\
x396,tx397,tx398,tx399,tx400,tx401,tx402,tx403,tx404,tx405,tx406,tx407,tx408,\
tx409,tx410,tx411,tx412,tx413,tx414,tx415,tx416,tx417,tx418,tx419,tx420,tx421\
,tx422,tx423,tx424,tx425,tx426,tx427,tx428,tx429,tx430,tx431,tx432,tx433,tx43\
4,tx435,tx436,tx437,tx438,tx439,tx440,tx441,tx442,tx443,tx444,tx445,tx446,tx4\
47,tx448,tx449,tx450,tx451,tx452,tx453,tx454,tx455,tx456,tx457,tx458,tx459,tx\
460,tx461,tx462,tx463,tx464,tx465,tx466,tx467,tx468,tx469,tx470,tx471,tx472,t\
x473,tx474,tx475,tx476,tx477,tx478,tx479,tx480,tx481,tx482,tx483,tx484,tx485,\
tx486,tx487,tx488,tx489,tx490,tx491,tx492,tx493,tx494,tx495,tx496,tx497,tx498\
,tx499,tx500,tx501,tx502,tx503,tx504,tx505,tx506,tx507,tx508,tx509,tx510,tx51\
1,tx512,tx513,tx514,tx515,tx516,tx517,tx518,tx519,tx520,tx521,tx522,tx523,tx5\
24,tx525,tx526,tx527,tx528,tx529,tx530,tx531,tx532,tx533,tx534,tx535,tx536,tx\
537,tx538,tx539,tx540,tx541,tx542,tx543,tx544,tx545,tx546,tx547,tx548,tx549,t\
x550,tx551,tx552,tx553,tx554,tx555,tx556,tx557,tx558,tx559,tx560,tx561,tx562,\
tx563,tx564,tx565,tx566,tx567,tx568,tx569,tx570,tx571,tx572,tx573,tx574,tx575\
,tx576,tx577,tx578,tx579,tx580,tx581,tx582,tx583,tx584,tx585,tx586,tx587,tx58\
8,tx589,tx590,tx591,tx592,tx593,tx594,tx595,tx596,tx597,tx598,tx599,tx600,tx6\
01,tx602,tx603,tx604,tx605,tx606,tx607,tx608,tx609,tx610,tx611,tx612,tx613,tx\
614,tx615,tx616,tx617,tx618,tx619,tx620,tx621,tx622,tx623,tx624,tx625,tx626,t\
x627,tx628,tx629,tx630,tx631,tx632,tx633,tx634,tx635,tx636,tx637,tx638,tx639,\
tx640,tx641,tx642,tx643,tx644,tx645,tx646,tx647,tx648,tx649,tx650,tx651,tx652\
,tx653,tx654,tx655,tx656,tx657,tx658,tx659,tx660,tx661,tx662,tx663,tx704,tx70\
5,tzz664,tzz665,tzz666,tzz667,tzz668,tzz669,tzz670,tzz671,tzz672,tzz673,tzz67\
4,tzz675,tzz676,tzz677,tzz678,tzz679,tzz680,tzz681,tzz682,tzz683,tzz684,tzz68\
5,tzz686,tzz687,tzz688,tzz689,tzz690,tzz691,tzz692,tzz693,tzz694,tzz695,tzz69\
6,tzz697,tzz698,tzz699,tzz700,tzz701,tzz702,tzz703,tzz706,tzz707,tzz708,tzz70\
9,tzz710,tzz711,tzz712,tzz713,tzz714,tzz715,tzz716,tzz717,tzz718,tzz719,tzz72\
0,tzz721,tzz722,tzz723,tzz724,tzz725,tzz726,xxxDummy},\n\t(*NONBONDRB_SET_PAR\
AMETER(I1);*)\n\t(*NONBONDRB_SET_PARAMETER(I2);*)\n\t(*NONBONDRB_SET_PARAMETE\
R(dQ1Q2);*)\n\t(*NONBONDRB_SET_PARAMETER(dA);*)\n\t(*NONBONDRB_SET_PARAMETER(\
dC);*)\n\t(*NONBONDRB_SET_POSITION(am,I1,0);*)\n\t(*NONBONDRB_SET_POSITION(bm\
,I1,1);*)\n\t(*NONBONDRB_SET_POSITION(cm,I1,2);*)\n\t(*NONBONDRB_SET_POSITION\
(dm,I1,3);*)\n\t(*NONBONDRB_SET_POSITION(xm,I1,4);*)\n\t(*NONBONDRB_SET_POSIT\
ION(ym,I1,5);*)\n\t(*NONBONDRB_SET_POSITION(zm,I1,6);*)\n\t(*NONBONDRB_SET_PO\
SITION(an,I2,0);*)\n\t(*NONBONDRB_SET_POSITION(bn,I2,1);*)\n\t(*NONBONDRB_SET\
_POSITION(cn,I2,2);*)\n\t(*NONBONDRB_SET_POSITION(dn,I2,3);*)\n\t(*NONBONDRB_\
SET_POSITION(xn,I2,4);*)\n\t(*NONBONDRB_SET_POSITION(yn,I2,5);*)\n\t(*NONBOND\
RB_SET_POSITION(zn,I2,6);*)\n\t(*NONBONDRB_SET_POINT(pxm,ea1,getX());*)\n\t(*\
NONBONDRB_SET_POINT(pym,ea1,getY());*)\n\t(*NONBONDRB_SET_POINT(pzm,ea1,getZ(\
));*)\n\t(*NONBONDRB_SET_POINT(pxn,ea2,getX());*)\n\t(*NONBONDRB_SET_POINT(py\
n,ea2,getY());*)\n\t(*NONBONDRB_SET_POINT(pzn,ea2,getZ());*)\n\ttx325= \
mathPower2[am]; (*Rule 26*)\n\ttx326= mathPower2[an]; (*Rule 27*)\n\ttx327= \
mathPower2[bm]; (*Rule 28*)\n\ttx328= mathPower2[bn]; (*Rule 29*)\n\ttx329= \
mathPower2[cm]; (*Rule 30*)\n\ttx330= mathPower2[cn]; (*Rule 31*)\n\ttx331= \
mathPower2[dm]; (*Rule 32*)\n\ttx332= mathPower2[dn]; (*Rule \
33*)\n\ttzz678=-pxm; (*Rule 34*)\n\ttzz700=am cm; (*Rule 35*)\n\ttx333=tzz678 \
tzz700; (*Rule 36*)\n\ttx334=bm cm pxm; (*Rule 37*)\n\ttzz696=dm pxm; (*Rule \
38*)\n\ttx335=am tzz696; (*Rule 39*)\n\ttx336=bm tzz696; (*Rule \
40*)\n\ttzz691=-pxn; (*Rule 41*)\n\ttx337=an cn tzz691; (*Rule \
42*)\n\ttx338=bn cn pxn; (*Rule 43*)\n\ttzz694=dn pxn; (*Rule \
44*)\n\ttx339=an tzz694; (*Rule 45*)\n\ttx340=bn tzz694; (*Rule \
46*)\n\ttzz701=am bm; (*Rule 47*)\n\ttx341=pym tzz701; (*Rule \
48*)\n\ttzz699=cm pym; (*Rule 49*)\n\ttx342=bm tzz699; (*Rule \
50*)\n\ttzz685=-pym; (*Rule 51*)\n\ttx343=am dm tzz685; (*Rule \
52*)\n\ttx344=dm tzz699; (*Rule 53*)\n\ttx345=an bn pyn; (*Rule \
54*)\n\ttzz698=cn pyn; (*Rule 55*)\n\ttx346=bn tzz698; (*Rule \
56*)\n\ttzz690=-pyn; (*Rule 57*)\n\ttx347=an dn tzz690; (*Rule \
58*)\n\ttx348=dn tzz698; (*Rule 59*)\n\ttzz677=-pzm; (*Rule \
60*)\n\ttx349=tzz677 tzz701; (*Rule 61*)\n\ttx350=pzm tzz700; (*Rule \
62*)\n\ttzz695=dm pzm; (*Rule 63*)\n\ttx351=bm tzz695; (*Rule \
64*)\n\ttx352=cm tzz695; (*Rule 65*)\n\ttzz686=bn pzn; (*Rule \
66*)\n\ttx353=-(an tzz686); (*Rule 67*)\n\ttzz697=cn pzn; (*Rule \
68*)\n\ttx354=an tzz697; (*Rule 69*)\n\ttx355=dn tzz686; (*Rule \
70*)\n\ttx356=dn tzz697; (*Rule 71*)\n\ttx357=tx327 tzz685; (*Rule \
72*)\n\ttx358=tx327 tzz677; (*Rule 73*)\n\ttx359=tx328 tzz690; (*Rule \
74*)\n\ttzz702=-pzn; (*Rule 75*)\n\ttx360=tx328 tzz702; (*Rule \
76*)\n\ttx361=tx329 tzz678; (*Rule 77*)\n\ttx362=tx329 tzz677; (*Rule \
78*)\n\ttx363=tx330 tzz691; (*Rule 79*)\n\ttx364=tx330 tzz702; (*Rule \
80*)\n\ttx365=tx331 tzz678; (*Rule 81*)\n\ttx366=tx331 tzz685; (*Rule \
82*)\n\ttx367=tx325 + tx327 + tx329 + tx331; (*Rule 83*)\n\ttx368=tx332 \
tzz691; (*Rule 84*)\n\ttx369=tx332 tzz690; (*Rule 85*)\n\ttx370=tx326 + tx328 \
+ tx330 + tx332; (*Rule 86*)\n\ttx371=tx333 + tx336 + tx341 + tx344 + tx358 + \
tx362; (*Rule 87*)\n\ttzz719=tx345 + tx348; (*Rule 88*)\n\ttzz720=tx340 + \
tzz719; (*Rule 89*)\n\ttzz723=tx337 + tzz720; (*Rule 90*)\n\ttx372=tx360 + \
tx364 + tzz723; (*Rule 91*)\n\ttx373=tx342 + tx343 + tx350 + tx351 + tx361 + \
tx365; (*Rule 92*)\n\ttx374=tx334 + tx335 + tx349 + tx352 + tx357 + tx366; \
(*Rule 93*)\n\ttx375= mathReciprocal[tx367]; (*Rule 94*)\n\ttzz715=tx354 + \
tx355; (*Rule 95*)\n\ttzz717=tx347 + tzz715; (*Rule 96*)\n\ttzz718=tx346 + \
tzz717; (*Rule 97*)\n\ttx376=tx363 + tx368 + tzz718; (*Rule \
98*)\n\ttzz716=tx353 + tx356; (*Rule 99*)\n\ttzz721=tx339 + tzz716; (*Rule \
100*)\n\ttzz722=tx338 + tzz721; (*Rule 101*)\n\ttx377=tx359 + tx369 + tzz722; \
(*Rule 102*)\n\ttx378= mathReciprocal[tx370]; (*Rule 103*)\n\ttx379=tzz691; \
(*Rule 104*)\n\ttx380=tzz690; (*Rule 105*)\n\ttx381=tzz702; (*Rule \
106*)\n\ttzz665=2. tx375; (*Rule 107*)\n\ttx382=tx371 tzz665; (*Rule \
108*)\n\ttx383=tx373 tzz665; (*Rule 109*)\n\ttx384=tx374 tzz665; (*Rule \
110*)\n\ttzz664=-2. tx378; (*Rule 111*)\n\ttx385=tx372 tzz664; (*Rule \
112*)\n\ttx386=tx376 tzz664; (*Rule 113*)\n\ttx387=tx377 tzz664; (*Rule \
114*)\n\ttx388=-xn; (*Rule 115*)\n\ttx389=-yn; (*Rule 116*)\n\ttx390=-zn; \
(*Rule 117*)\n\ttzz708=tx388 + xm; (*Rule 118*)\n\ttzz710=tx383 + tzz708; \
(*Rule 119*)\n\ttzz714=tx379 + tzz710; (*Rule 120*)\n\ttzz726=pxm + tzz714; \
(*Rule 121*)\n\ttx391=tx386 + tzz726; (*Rule 122*)\n\ttzz707=tx389 + ym; \
(*Rule 123*)\n\ttzz709=tx384 + tzz707; (*Rule 124*)\n\ttzz713=tx380 + tzz709; \
(*Rule 125*)\n\ttzz725=pym + tzz713; (*Rule 126*)\n\ttx392=tx387 + tzz725; \
(*Rule 127*)\n\ttzz706=tx390 + zm; (*Rule 128*)\n\ttzz711=tx382 + tzz706; \
(*Rule 129*)\n\ttzz712=tx381 + tzz711; (*Rule 130*)\n\ttzz724=pzm + tzz712; \
(*Rule 131*)\n\ttx393=tx385 + tzz724; (*Rule 132*)\n\ttx394= \
mathPower2[tx391]; (*Rule 133*)\n\ttx395= mathPower2[tx392]; (*Rule \
134*)\n\ttx396= mathPower2[tx393]; (*Rule 135*)\n\ttx397=tx394 + tx395 + \
tx396; (*Rule 136*)\n\ttx704= mathReciprocal[tx397]; (*Rule 137*)\n\ttx660= \
mathPower2[tx704]; (*Rule 138*)\n\ttx661=tx704; (*Rule 139*)\n\ttx649=tx660 \
tx661; (*Rule 140*)\n\ttx398= mathPower2[tx649]; (*Rule 141*)\n\ttx650= \
mathPower2[tx661]; (*Rule 142*)\n\ttx651=tx661; (*Rule 143*)\n\ttx399=tx650 \
tx651; (*Rule 144*)\n\ttx400=dA tx398; (*Rule 145*)\n\ttzz703=-dC; (*Rule \
146*)\n\ttx401=tx399 tzz703; (*Rule 147*)\n\tEnergyVdw=tx400 + tx401; (*Rule \
148*)\n\ttx402=tx330 tx379; (*Rule 149*)\n\ttx403=tx332 tx379; (*Rule \
150*)\n\ttx404=tx328 tx380; (*Rule 151*)\n\ttx405=tx332 tx380; (*Rule \
152*)\n\ttx406=tx328 tx381; (*Rule 153*)\n\ttx407=tx330 tx381; (*Rule \
154*)\n\ttx408=tx402 + tx403 + tzz718; (*Rule 155*)\n\ttx409=tx404 + tx405 + \
tzz722; (*Rule 156*)\n\ttx410=tx406 + tx407 + tzz723; (*Rule \
157*)\n\ttx411=tx408 tzz664; (*Rule 158*)\n\ttx412=tx409 tzz664; (*Rule \
159*)\n\ttx413=tx410 tzz664; (*Rule 160*)\n\ttx414=tx411 + tzz726; (*Rule \
161*)\n\ttx415=tx412 + tzz725; (*Rule 162*)\n\ttx416=tx413 + tzz724; (*Rule \
163*)\n\ttx417= mathPower2[tx414]; (*Rule 164*)\n\ttx418= mathPower2[tx415]; \
(*Rule 165*)\n\ttx419= mathPower2[tx416]; (*Rule 166*)\n\ttx420=tx417 + tx418 \
+ tx419; (*Rule 167*)\n\ttx652= mathSqrt[tx420]; (*Rule 168*)\n\ttx421= \
mathReciprocal[tx652]; (*Rule 169*)\n\tEnergyElectrostatic=dQ1Q2 tx421; \
(*Rule 170*)\n\ttx705= mathReciprocal[tx420]; (*Rule 171*)\n\ttx662= \
mathPower2[tx705]; (*Rule 172*)\n\ttx663=tx705; (*Rule 173*)\n\ttx653=tx662 \
tx663; (*Rule 174*)\n\ttx422= mathPower2[tx653]; (*Rule 175*)\n\ttx654= \
mathPower2[tx663]; (*Rule 176*)\n\ttx655=tx663; (*Rule 177*)\n\ttx423=tx654 \
tx655; (*Rule 178*)\n\ttx424=dA tx422; (*Rule 179*)\n\ttx425=tx423 tzz703; \
(*Rule 180*)\n\tEnergy=EnergyElectrostatic + tx424 + tx425; (*Rule \
181*)\n\t(*NONBONDRB_EVDW_ENERGY_ACCUMULATE(EnergyVdw);*)\n\t(*NONBONDRB_EEEL\
_ENERGY_ACCUMULATE(EnergyElectrostatic);*)\n\t(*NONBONDRB_ENERGY_ACCUMULATE(E\
nergy);*)\n\t(*#ifdef NONBONDRB_CALC_FORCE //[*)\n\t(*if ( calcForce ) \
{*)\n\ttx426=cm tzz678; (*Rule 187*)\n\ttx427=tzz696; (*Rule \
188*)\n\ttx428=bm pym; (*Rule 189*)\n\ttx429=dm tzz685; (*Rule \
190*)\n\ttx430=bm tzz677; (*Rule 191*)\n\ttx431=cm pzm; (*Rule \
192*)\n\ttx432= mathPower2[tx375]; (*Rule 193*)\n\ttx433=tx426 + tx428; \
(*Rule 194*)\n\ttx434=tx427 + tx430; (*Rule 195*)\n\ttx435=tx429 + tx431; \
(*Rule 196*)\n\ttzz670=-4. tx432; (*Rule 197*)\n\ttzz684=tx371 tzz670; (*Rule \
198*)\n\ttx436=am tzz684; (*Rule 199*)\n\ttzz683=tx373 tzz670; (*Rule \
200*)\n\ttx437=am tzz683; (*Rule 201*)\n\ttzz682=tx374 tzz670; (*Rule \
202*)\n\ttx438=am tzz682; (*Rule 203*)\n\ttx439=tx433 tzz665; (*Rule \
204*)\n\ttx440=tx434 tzz665; (*Rule 205*)\n\ttx441=tx435 tzz665; (*Rule \
206*)\n\ttx442=tx436 + tx439; (*Rule 207*)\n\ttx443=tx438 + tx440; (*Rule \
208*)\n\ttx444=tx437 + tx441; (*Rule 209*)\n\ttzz672=2. tx416; (*Rule \
210*)\n\ttx445=tx442 tzz672; (*Rule 211*)\n\ttzz673=2. tx415; (*Rule \
212*)\n\ttx446=tx443 tzz673; (*Rule 213*)\n\ttzz674=2. tx414; (*Rule \
214*)\n\ttx447=tx444 tzz674; (*Rule 215*)\n\ttx656= mathPower2[tx662]; (*Rule \
216*)\n\ttx448=tx653 tx656; (*Rule 217*)\n\ttx449= mathPower2[tx654]; (*Rule \
218*)\n\ttx657=tx420; (*Rule 219*)\n\ttx658=tx421; (*Rule 220*)\n\ttx659= \
mathReciprocal[tx657]; (*Rule 221*)\n\ttx450=tx658 tx659; (*Rule \
222*)\n\ttx451=tx445 + tx446 + tx447; (*Rule 223*)\n\ttzz668=dA tx448; (*Rule \
224*)\n\ttzz676=-6. tzz668; (*Rule 225*)\n\ttx452=tx451 tzz676; (*Rule \
226*)\n\ttzz667=dC tx449; (*Rule 227*)\n\ttzz671=3 tzz667; (*Rule \
228*)\n\ttx453=tx451 tzz671; (*Rule 229*)\n\ttzz666=dQ1Q2 tx450; (*Rule \
230*)\n\ttzz675=-0.5 tzz666; (*Rule 231*)\n\ttx454=tx451 tzz675; (*Rule \
232*)\n\tgam=tx452 + tx453 + tx454; (*Rule 233*)\n\tfam=-gam; (*Rule \
234*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I1, 0, fam );*)\n\ttx455=cm pxm; (*Rule \
236*)\n\ttx456=am pym; (*Rule 237*)\n\ttx457=tzz699; (*Rule 238*)\n\ttx458=am \
tzz677; (*Rule 239*)\n\ttx459=-2. bm pzm; (*Rule 240*)\n\ttx460=tzz695; \
(*Rule 241*)\n\ttx461=-2. tx428; (*Rule 242*)\n\ttx462=tx427 + tx456 + tx459; \
(*Rule 243*)\n\ttx463=tx457 + tx460; (*Rule 244*)\n\ttx464=tx455 + tx458 + \
tx461; (*Rule 245*)\n\ttx465=bm tzz684; (*Rule 246*)\n\ttx466=bm tzz683; \
(*Rule 247*)\n\ttx467=bm tzz682; (*Rule 248*)\n\ttx468=tx462 tzz665; (*Rule \
249*)\n\ttx469=tx463 tzz665; (*Rule 250*)\n\ttx470=tx464 tzz665; (*Rule \
251*)\n\ttx471=tx465 + tx468; (*Rule 252*)\n\ttx472=tx466 + tx469; (*Rule \
253*)\n\ttx473=tx467 + tx470; (*Rule 254*)\n\ttx474=tx471 tzz672; (*Rule \
255*)\n\ttx475=tx472 tzz674; (*Rule 256*)\n\ttx476=tx473 tzz673; (*Rule \
257*)\n\ttx477=tx474 + tx475 + tx476; (*Rule 258*)\n\ttx478=tx477 tzz676; \
(*Rule 259*)\n\ttx479=tx477 tzz671; (*Rule 260*)\n\ttx480=tx477 tzz675; \
(*Rule 261*)\n\tgbm=tx478 + tx479 + tx480; (*Rule 262*)\n\tfbm=-gbm; (*Rule \
263*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I1, 1, fbm );*)\n\ttx481=am tzz678; \
(*Rule 265*)\n\ttx482=bm pxm; (*Rule 266*)\n\ttx483=dm pym; (*Rule \
267*)\n\ttx484=am pzm; (*Rule 268*)\n\ttx485=-2. tx431; (*Rule \
269*)\n\ttx486=-2. tx455; (*Rule 270*)\n\ttx487=tx460 + tx482; (*Rule \
271*)\n\ttx488=tx481 + tx483 + tx485; (*Rule 272*)\n\ttx489=tx428 + tx484 + \
tx486; (*Rule 273*)\n\ttx490=cm tzz684; (*Rule 274*)\n\ttx491=cm tzz683; \
(*Rule 275*)\n\ttx492=cm tzz682; (*Rule 276*)\n\ttx493=tx487 tzz665; (*Rule \
277*)\n\ttx494=tx488 tzz665; (*Rule 278*)\n\ttx495=tx489 tzz665; (*Rule \
279*)\n\ttx496=tx492 + tx493; (*Rule 280*)\n\ttx497=tx490 + tx494; (*Rule \
281*)\n\ttx498=tx491 + tx495; (*Rule 282*)\n\ttx499=tx496 tzz673; (*Rule \
283*)\n\ttx500=tx497 tzz672; (*Rule 284*)\n\ttx501=tx498 tzz674; (*Rule \
285*)\n\ttx502=tx499 + tx500 + tx501; (*Rule 286*)\n\ttx503=tx502 tzz676; \
(*Rule 287*)\n\ttx504=tx502 tzz671; (*Rule 288*)\n\ttx505=tx502 tzz675; \
(*Rule 289*)\n\tgcm=tx503 + tx504 + tx505; (*Rule 290*)\n\tfcm=-gcm; (*Rule \
291*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I1, 2, fcm );*)\n\ttx506=am pxm; (*Rule \
293*)\n\ttx507=bm pzm; (*Rule 294*)\n\ttx508=-2. tx427; (*Rule \
295*)\n\ttx509=-tx456; (*Rule 296*)\n\ttx510=-2. tx483; (*Rule \
297*)\n\ttx511=tx457 + tx482; (*Rule 298*)\n\ttx512=tx507 + tx508 + tx509; \
(*Rule 299*)\n\ttx513=tx431 + tx506 + tx510; (*Rule 300*)\n\ttx514=dm tzz684; \
(*Rule 301*)\n\ttx515=dm tzz683; (*Rule 302*)\n\ttx516=dm tzz682; (*Rule \
303*)\n\ttx517=tx511 tzz665; (*Rule 304*)\n\ttx518=tx512 tzz665; (*Rule \
305*)\n\ttx519=tx513 tzz665; (*Rule 306*)\n\ttx520=tx514 + tx517; (*Rule \
307*)\n\ttx521=tx515 + tx518; (*Rule 308*)\n\ttx522=tx516 + tx519; (*Rule \
309*)\n\ttx523=tx520 tzz672; (*Rule 310*)\n\ttx524=tx521 tzz674; (*Rule \
311*)\n\ttx525=tx522 tzz673; (*Rule 312*)\n\ttx526=tx523 + tx524 + tx525; \
(*Rule 313*)\n\ttx527=tx526 tzz676; (*Rule 314*)\n\ttx528=tx526 tzz671; \
(*Rule 315*)\n\ttx529=tx526 tzz675; (*Rule 316*)\n\tgdm=tx527 + tx528 + \
tx529; (*Rule 317*)\n\tfdm=-gdm; (*Rule \
318*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I1, 3, fdm );*)\n\ttzz693=-12. tzz668; \
(*Rule 320*)\n\ttx530=tx414 tzz693; (*Rule 321*)\n\ttzz688=6. tzz667; (*Rule \
322*)\n\ttx531=tx414 tzz688; (*Rule 323*)\n\ttzz689=-tzz666; (*Rule \
324*)\n\ttx532=tx414 tzz689; (*Rule 325*)\n\tgxm=tx530 + tx531 + tx532; \
(*Rule 326*)\n\tfxm=-gxm; (*Rule 327*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I1, 4, \
fxm );*)\n\ttx533=tx415 tzz693; (*Rule 329*)\n\ttx534=tx415 tzz688; (*Rule \
330*)\n\ttx535=tx415 tzz689; (*Rule 331*)\n\tgym=tx533 + tx534 + tx535; \
(*Rule 332*)\n\tfym=-gym; (*Rule 333*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I1, 5, \
fym );*)\n\ttx536=tx416 tzz693; (*Rule 335*)\n\ttx537=tx416 tzz688; (*Rule \
336*)\n\ttx538=tx416 tzz689; (*Rule 337*)\n\tgzm=tx536 + tx537 + tx538; \
(*Rule 338*)\n\tfzm=-gzm; (*Rule 339*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I1, 6, \
fzm );*)\n\ttx539=tzz694; (*Rule 341*)\n\ttx540=bn pyn; (*Rule \
342*)\n\ttx541=tzz697; (*Rule 343*)\n\ttx542=cn tx379; (*Rule \
344*)\n\ttx543=dn tx380; (*Rule 345*)\n\ttx544=bn tx381; (*Rule \
346*)\n\ttx545= mathPower2[tx378]; (*Rule 347*)\n\ttx546=tx540 + tx542; \
(*Rule 348*)\n\ttx547=tx541 + tx543; (*Rule 349*)\n\ttx548=tx539 + tx544; \
(*Rule 350*)\n\ttzz669=4. tx545; (*Rule 351*)\n\ttzz681=tx408 tzz669; (*Rule \
352*)\n\ttx549=an tzz681; (*Rule 353*)\n\ttzz680=tx409 tzz669; (*Rule \
354*)\n\ttx550=an tzz680; (*Rule 355*)\n\ttzz679=tx410 tzz669; (*Rule \
356*)\n\ttx551=an tzz679; (*Rule 357*)\n\ttx552=tx546 tzz664; (*Rule \
358*)\n\ttx553=tx547 tzz664; (*Rule 359*)\n\ttx554=tx548 tzz664; (*Rule \
360*)\n\ttx555=tx551 + tx552; (*Rule 361*)\n\ttx556=tx549 + tx553; (*Rule \
362*)\n\ttx557=tx550 + tx554; (*Rule 363*)\n\ttx558=tx555 tzz672; (*Rule \
364*)\n\ttx559=tx556 tzz674; (*Rule 365*)\n\ttx560=tx557 tzz673; (*Rule \
366*)\n\ttx561=tx558 + tx559 + tx560; (*Rule 367*)\n\ttx562=tx561 tzz676; \
(*Rule 368*)\n\ttx563=tx561 tzz671; (*Rule 369*)\n\ttx564=tx561 tzz675; \
(*Rule 370*)\n\tgan=tx562 + tx563 + tx564; (*Rule 371*)\n\tfan=-gan; (*Rule \
372*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I2, 0, fan );*)\n\ttx565=cn pxn; (*Rule \
374*)\n\ttx566=an pyn; (*Rule 375*)\n\ttx567=tzz698; (*Rule \
376*)\n\ttx568=-2. tzz686; (*Rule 377*)\n\ttx569=dn pzn; (*Rule \
378*)\n\ttx570=an tx381; (*Rule 379*)\n\ttx571=-2. tx540; (*Rule \
380*)\n\ttx572=tx539 + tx566 + tx568; (*Rule 381*)\n\ttx573=tx567 + tx569; \
(*Rule 382*)\n\ttx574=tx565 + tx570 + tx571; (*Rule 383*)\n\ttx575=bn tzz681; \
(*Rule 384*)\n\ttx576=bn tzz680; (*Rule 385*)\n\ttx577=bn tzz679; (*Rule \
386*)\n\ttx578=tx572 tzz664; (*Rule 387*)\n\ttx579=tx573 tzz664; (*Rule \
388*)\n\ttx580=tx574 tzz664; (*Rule 389*)\n\ttx581=tx577 + tx578; (*Rule \
390*)\n\ttx582=tx575 + tx579; (*Rule 391*)\n\ttx583=tx576 + tx580; (*Rule \
392*)\n\ttx584=tx581 tzz672; (*Rule 393*)\n\ttx585=tx582 tzz674; (*Rule \
394*)\n\ttx586=tx583 tzz673; (*Rule 395*)\n\ttx587=tx584 + tx585 + tx586; \
(*Rule 396*)\n\ttx588=tx587 tzz676; (*Rule 397*)\n\ttx589=tx587 tzz671; \
(*Rule 398*)\n\ttx590=tx587 tzz675; (*Rule 399*)\n\tgbn=tx588 + tx589 + \
tx590; (*Rule 400*)\n\tfbn=-gbn; (*Rule \
401*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I2, 1, fbn );*)\n\ttx591=bn pxn; (*Rule \
403*)\n\ttx592=dn pyn; (*Rule 404*)\n\ttx593=an pzn; (*Rule 405*)\n\ttx594=an \
tx379; (*Rule 406*)\n\ttx595=-2. tx541; (*Rule 407*)\n\ttx596=-2. tx565; \
(*Rule 408*)\n\ttx597=tx569 + tx591; (*Rule 409*)\n\ttx598=tx592 + tx594 + \
tx595; (*Rule 410*)\n\ttx599=tx540 + tx593 + tx596; (*Rule 411*)\n\ttx600=cn \
tzz681; (*Rule 412*)\n\ttx601=cn tzz680; (*Rule 413*)\n\ttx602=cn tzz679; \
(*Rule 414*)\n\ttx603=tx597 tzz664; (*Rule 415*)\n\ttx604=tx598 tzz664; \
(*Rule 416*)\n\ttx605=tx599 tzz664; (*Rule 417*)\n\ttx606=tx601 + tx603; \
(*Rule 418*)\n\ttx607=tx602 + tx604; (*Rule 419*)\n\ttx608=tx600 + tx605; \
(*Rule 420*)\n\ttx609=tx606 tzz673; (*Rule 421*)\n\ttx610=tx607 tzz672; \
(*Rule 422*)\n\ttx611=tx608 tzz674; (*Rule 423*)\n\ttx612=tx609 + tx610 + \
tx611; (*Rule 424*)\n\ttx613=tx612 tzz676; (*Rule 425*)\n\ttx614=tx612 \
tzz671; (*Rule 426*)\n\ttx615=tx612 tzz675; (*Rule 427*)\n\tgcn=tx613 + tx614 \
+ tx615; (*Rule 428*)\n\tfcn=-gcn; (*Rule \
429*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I2, 2, fcn );*)\n\ttx616=an pxn; (*Rule \
431*)\n\ttx617=tzz686; (*Rule 432*)\n\ttx618=an tx380; (*Rule \
433*)\n\ttx619=-2. tx539; (*Rule 434*)\n\ttx620=-2. tx592; (*Rule \
435*)\n\ttx621=tx567 + tx591; (*Rule 436*)\n\ttx622=tx617 + tx618 + tx619; \
(*Rule 437*)\n\ttx623=tx541 + tx616 + tx620; (*Rule 438*)\n\ttx624=dn tzz681; \
(*Rule 439*)\n\ttx625=dn tzz680; (*Rule 440*)\n\ttx626=dn tzz679; (*Rule \
441*)\n\ttx627=tx621 tzz664; (*Rule 442*)\n\ttx628=tx622 tzz664; (*Rule \
443*)\n\ttx629=tx623 tzz664; (*Rule 444*)\n\ttx630=tx626 + tx627; (*Rule \
445*)\n\ttx631=tx624 + tx628; (*Rule 446*)\n\ttx632=tx625 + tx629; (*Rule \
447*)\n\ttx633=tx630 tzz672; (*Rule 448*)\n\ttx634=tx631 tzz674; (*Rule \
449*)\n\ttx635=tx632 tzz673; (*Rule 450*)\n\ttx636=tx633 + tx634 + tx635; \
(*Rule 451*)\n\ttx637=tx636 tzz676; (*Rule 452*)\n\ttx638=tx636 tzz671; \
(*Rule 453*)\n\ttx639=tx636 tzz675; (*Rule 454*)\n\tgdn=tx637 + tx638 + \
tx639; (*Rule 455*)\n\tfdn=-gdn; (*Rule \
456*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I2, 3, fdn );*)\n\ttzz687=12. tzz668; \
(*Rule 458*)\n\ttx640=tx414 tzz687; (*Rule 459*)\n\ttzz692=-6. tzz667; (*Rule \
460*)\n\ttx641=tx414 tzz692; (*Rule 461*)\n\ttx642=tx414 tzz666; (*Rule \
462*)\n\tgxn=tx640 + tx641 + tx642; (*Rule 463*)\n\tfxn=-gxn; (*Rule \
464*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I2, 4, fxn );*)\n\ttx643=tx415 tzz687; \
(*Rule 466*)\n\ttx644=tx415 tzz692; (*Rule 467*)\n\ttx645=tx415 tzz666; \
(*Rule 468*)\n\tgyn=tx643 + tx644 + tx645; (*Rule 469*)\n\tfyn=-gyn; (*Rule \
470*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I2, 5, fyn );*)\n\ttx646=tx416 tzz687; \
(*Rule 472*)\n\ttx647=tx416 tzz692; (*Rule 473*)\n\ttx648=tx416 tzz666; \
(*Rule 474*)\n\tgzn=tx646 + tx647 + tx648; (*Rule 475*)\n\tfzn=-gzn; (*Rule \
476*)\n\t(*NONBONDRB_FORCE_ACCUMULATE(I2, 6, fzn );*)\n\t(*} /*calcForce \
*/*)\n\t(*#endif /* NONBONDRB_CALC_FORCE ]*/*)\n{Energy, EnergyVdw, \
EnergyElectrostatic, fam, fbm, fcm, fdm, fxm, fym, fzm, fan, fbn, fcn, fdn, \
fxn, fyn, fzn, nonbondRBDeviation}]\n"}
