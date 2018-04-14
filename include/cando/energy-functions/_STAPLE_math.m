{Name -> "STAPLE", AdditionalCDeclares -> "", 
 Input -> {ak, bk, ck, dk, xk, yk, zk, al, bl, cl, dl, xl, yl, zl}, 
 Output -> {Energy, fak, fbk, fck, fdk, fxk, fyk, fzk, fal, fbl, fcl, fdl, 
   fxl, fyl, fzl}, EnergyFunction -> EnergyFunction, 
 DerivativeVariables -> {ak, bk, ck, dk, xk, yk, zk, al, bl, cl, dl, xl, yl, 
   zl}, Rules -> {CCode["STAPLE_SET_PARAMETER(ks,ks);"], 
   CCode["STAPLE_SET_PARAMETER(r0,r0);"], 
   CCode["STAPLE_SET_PARAMETER(I1,rigidBodyK);"], 
   CCode["STAPLE_SET_PARAMETER(I2,rigidBodyL);"], 
   CCode["STAPLE_SET_POSITION(ak,I1,0);"], 
   CCode["STAPLE_SET_POSITION(bk,I1,1);"], 
   CCode["STAPLE_SET_POSITION(ck,I1,2);"], 
   CCode["STAPLE_SET_POSITION(dk,I1,3);"], 
   CCode["STAPLE_SET_POSITION(xk,I1,4);"], 
   CCode["STAPLE_SET_POSITION(yk,I1,5);"], 
   CCode["STAPLE_SET_POSITION(zk,I1,6);"], 
   CCode["STAPLE_SET_POSITION(al,I2,0);"], 
   CCode["STAPLE_SET_POSITION(bl,I2,1);"], 
   CCode["STAPLE_SET_POSITION(cl,I2,2);"], 
   CCode["STAPLE_SET_POSITION(dl,I2,3);"], 
   CCode["STAPLE_SET_POSITION(xl,I2,4);"], 
   CCode["STAPLE_SET_POSITION(yl,I2,5);"], 
   CCode["STAPLE_SET_POSITION(zl,I2,6);"], 
   CCode["STAPLE_SET_POINT(pxk,pointK,getX());"], 
   CCode["STAPLE_SET_POINT(pyk,pointK,getY());"], 
   CCode["STAPLE_SET_POINT(pzk,pointK,getZ());"], 
   CCode["STAPLE_SET_POINT(pxl,pointL,getX());"], 
   CCode["STAPLE_SET_POINT(pyl,pointL,getY());"], 
   CCode["STAPLE_SET_POINT(pzl,pointL,getZ());"], power2[ak] -> tx278, 
   power2[al] -> tx279, power2[bk] -> tx280, power2[bl] -> tx281, 
   power2[ck] -> tx282, power2[cl] -> tx283, power2[dk] -> tx284, 
   power2[dl] -> tx285, -pxk -> tzz565, ak*ck -> tzz583, 
   tzz565*tzz583 -> tx286, bk*ck*pxk -> tx287, dk*pxk -> tzz579, 
   ak*tzz579 -> tx288, bk*tzz579 -> tx289, -pxl -> tzz575, 
   al*cl*tzz575 -> tx290, bl*cl*pxl -> tx291, dl*pxl -> tzz577, 
   al*tzz577 -> tx292, bl*tzz577 -> tx293, ak*bk -> tzz584, 
   pyk*tzz584 -> tx294, ck*pyk -> tzz582, bk*tzz582 -> tx295, -pyk -> tzz572, 
   ak*dk*tzz572 -> tx296, dk*tzz582 -> tx297, al*bl*pyl -> tx298, 
   cl*pyl -> tzz581, bl*tzz581 -> tx299, -pyl -> tzz574, 
   al*dl*tzz574 -> tx300, dl*tzz581 -> tx301, -pzk -> tzz564, 
   tzz564*tzz584 -> tx302, pzk*tzz583 -> tx303, dk*pzk -> tzz578, 
   bk*tzz578 -> tx304, ck*tzz578 -> tx305, bl*pzl -> tzz573, 
   -(al*tzz573) -> tx306, cl*pzl -> tzz580, al*tzz580 -> tx307, 
   dl*tzz573 -> tx308, dl*tzz580 -> tx309, tx280*tzz572 -> tx310, 
   tx280*tzz564 -> tx311, tx281*tzz574 -> tx312, -pzl -> tzz585, 
   tx281*tzz585 -> tx313, tx282*tzz565 -> tx314, tx282*tzz564 -> tx315, 
   tx283*tzz575 -> tx316, tx283*tzz585 -> tx317, tx284*tzz565 -> tx318, 
   tx284*tzz572 -> tx319, tx278 + tx280 + tx282 + tx284 -> tx320, 
   tx285*tzz575 -> tx321, tx285*tzz574 -> tx322, 
   tx279 + tx281 + tx283 + tx285 -> tx323, 
   tx286 + tx289 + tx294 + tx297 + tx311 + tx315 -> tx324, 
   tx298 + tx301 -> tzz599, tx293 + tzz599 -> tzz600, 
   tx290 + tzz600 -> tzz603, tx313 + tx317 + tzz603 -> tx325, 
   tx295 + tx296 + tx303 + tx304 + tx314 + tx318 -> tx326, 
   tx287 + tx288 + tx302 + tx305 + tx310 + tx319 -> tx327, 
   reciprocal[tx320] -> tx328, tx307 + tx308 -> tzz595, 
   tx300 + tzz595 -> tzz597, tx299 + tzz597 -> tzz598, 
   tx316 + tx321 + tzz598 -> tx329, tx306 + tx309 -> tzz596, 
   tx292 + tzz596 -> tzz601, tx291 + tzz601 -> tzz602, 
   tx312 + tx322 + tzz602 -> tx330, reciprocal[tx323] -> tx331, 
   tzz575 -> tx332, tzz574 -> tx333, tzz585 -> tx334, 2.*tx328 -> tzz556, 
   tx324*tzz556 -> tx335, tx326*tzz556 -> tx336, tx327*tzz556 -> tx337, 
   -2.*tx331 -> tzz555, tx325*tzz555 -> tx338, tx329*tzz555 -> tx339, 
   tx330*tzz555 -> tx340, -xl -> tx341, -yl -> tx342, -zl -> tx343, 
   tx341 + xk -> tzz588, tx336 + tzz588 -> tzz590, tx332 + tzz590 -> tzz594, 
   pxk + tzz594 -> tzz606, tx339 + tzz606 -> tx344, tx342 + yk -> tzz587, 
   tx337 + tzz587 -> tzz589, tx333 + tzz589 -> tzz593, 
   pyk + tzz593 -> tzz605, tx340 + tzz605 -> tx345, tx343 + zk -> tzz586, 
   tx335 + tzz586 -> tzz591, tx334 + tzz591 -> tzz592, 
   pzk + tzz592 -> tzz604, tx338 + tzz604 -> tx346, power2[tx344] -> tx347, 
   power2[tx345] -> tx348, power2[tx346] -> tx349, 
   tx347 + tx348 + tx349 -> tx350, -r0 -> tx351, mysqrt[tx350] -> tx352, 
   tx351 + tx352 -> tx353, power2[tx353] -> tx354, ks*tx354 -> Energy, 
   CCode["STAPLE_ENERGY_ACCUMULATE(Energy);"], 
   CCode["#ifdef STAPLE_CALC_FORCE //["], CCode["if ( calcForce ) {"], 
   tx283*tx332 -> tx355, tx285*tx332 -> tx356, tx281*tx333 -> tx357, 
   tx285*tx333 -> tx358, tx281*tx334 -> tx359, tx283*tx334 -> tx360, 
   tx355 + tx356 + tzz598 -> tx361, tx357 + tx358 + tzz602 -> tx362, 
   tx359 + tx360 + tzz603 -> tx363, ck*tzz565 -> tx364, tzz579 -> tx365, 
   bk*pyk -> tx366, dk*tzz572 -> tx367, bk*tzz564 -> tx368, ck*pzk -> tx369, 
   tx361*tzz555 -> tx370, tx362*tzz555 -> tx371, tx363*tzz555 -> tx372, 
   power2[tx328] -> tx373, tx364 + tx366 -> tx374, tx365 + tx368 -> tx375, 
   tx367 + tx369 -> tx376, tx370 + tzz606 -> tx377, tx371 + tzz605 -> tx378, 
   tx372 + tzz604 -> tx379, -4.*tx373 -> tzz560, tx324*tzz560 -> tzz571, 
   ak*tzz571 -> tx380, tx326*tzz560 -> tzz570, ak*tzz570 -> tx381, 
   tx327*tzz560 -> tzz569, ak*tzz569 -> tx382, tx374*tzz556 -> tx383, 
   tx375*tzz556 -> tx384, tx376*tzz556 -> tx385, power2[tx377] -> tx386, 
   power2[tx378] -> tx387, power2[tx379] -> tx388, tx380 + tx383 -> tx389, 
   tx382 + tx384 -> tx390, tx381 + tx385 -> tx391, 
   tx386 + tx387 + tx388 -> tx392, 2.*tx379 -> tzz561, tx389*tzz561 -> tx393, 
   2.*tx378 -> tzz562, tx390*tzz562 -> tx394, 2.*tx377 -> tzz563, 
   tx391*tzz563 -> tx395, mysqrt[tx392] -> tx396, reciprocal[tx396] -> tx397, 
   tx393 + tx394 + tx395 -> tx398, tx351 + tx396 -> tx399, 
   tx397*tx399 -> tzz557, ks*tzz557 -> tzz558, tx398*tzz558 -> gak, 
   -gak -> fak, CCode["STAPLE_FORCE_ACCUMULATE(I1, 0, fak );"], 
   ck*pxk -> tx400, ak*pyk -> tx401, tzz582 -> tx402, ak*tzz564 -> tx403, 
   -2.*bk*pzk -> tx404, tzz578 -> tx405, -2.*tx366 -> tx406, 
   tx365 + tx401 + tx404 -> tx407, tx402 + tx405 -> tx408, 
   tx400 + tx403 + tx406 -> tx409, bk*tzz571 -> tx410, bk*tzz570 -> tx411, 
   bk*tzz569 -> tx412, tx407*tzz556 -> tx413, tx408*tzz556 -> tx414, 
   tx409*tzz556 -> tx415, tx410 + tx413 -> tx416, tx411 + tx414 -> tx417, 
   tx412 + tx415 -> tx418, tx416*tzz561 -> tx419, tx417*tzz563 -> tx420, 
   tx418*tzz562 -> tx421, tx419 + tx420 + tx421 -> tx422, 
   tx422*tzz558 -> gbk, -gbk -> fbk, 
   CCode["STAPLE_FORCE_ACCUMULATE(I1, 1, fbk );"], ak*tzz565 -> tx423, 
   bk*pxk -> tx424, dk*pyk -> tx425, ak*pzk -> tx426, -2.*tx369 -> tx427, 
   -2.*tx400 -> tx428, tx405 + tx424 -> tx429, tx423 + tx425 + tx427 -> 
    tx430, tx366 + tx426 + tx428 -> tx431, ck*tzz571 -> tx432, 
   ck*tzz570 -> tx433, ck*tzz569 -> tx434, tx429*tzz556 -> tx435, 
   tx430*tzz556 -> tx436, tx431*tzz556 -> tx437, tx434 + tx435 -> tx438, 
   tx432 + tx436 -> tx439, tx433 + tx437 -> tx440, tx438*tzz562 -> tx441, 
   tx439*tzz561 -> tx442, tx440*tzz563 -> tx443, 
   tx441 + tx442 + tx443 -> tx444, tx444*tzz558 -> gck, -gck -> fck, 
   CCode["STAPLE_FORCE_ACCUMULATE(I1, 2, fck );"], ak*pxk -> tx445, 
   bk*pzk -> tx446, -2.*tx365 -> tx447, -tx401 -> tx448, -2.*tx425 -> tx449, 
   tx402 + tx424 -> tx450, tx446 + tx447 + tx448 -> tx451, 
   tx369 + tx445 + tx449 -> tx452, dk*tzz571 -> tx453, dk*tzz570 -> tx454, 
   dk*tzz569 -> tx455, tx450*tzz556 -> tx456, tx451*tzz556 -> tx457, 
   tx452*tzz556 -> tx458, tx453 + tx456 -> tx459, tx454 + tx457 -> tx460, 
   tx455 + tx458 -> tx461, tx459*tzz561 -> tx462, tx460*tzz563 -> tx463, 
   tx461*tzz562 -> tx464, tx462 + tx463 + tx464 -> tx465, 
   tx465*tzz558 -> gdk, -gdk -> fdk, 
   CCode["STAPLE_FORCE_ACCUMULATE(I1, 3, fdk );"], tzz558*tzz563 -> gxk, 
   -gxk -> fxk, CCode["STAPLE_FORCE_ACCUMULATE(I1, 4, fxk );"], 
   tzz558*tzz562 -> gyk, -gyk -> fyk, 
   CCode["STAPLE_FORCE_ACCUMULATE(I1, 5, fyk );"], tzz558*tzz561 -> gzk, 
   -gzk -> fzk, CCode["STAPLE_FORCE_ACCUMULATE(I1, 6, fzk );"], 
   tzz577 -> tx466, bl*pyl -> tx467, tzz580 -> tx468, cl*tx332 -> tx469, 
   dl*tx333 -> tx470, bl*tx334 -> tx471, power2[tx331] -> tx472, 
   tx467 + tx469 -> tx473, tx468 + tx470 -> tx474, tx466 + tx471 -> tx475, 
   4.*tx472 -> tzz559, tx361*tzz559 -> tzz568, al*tzz568 -> tx476, 
   tx362*tzz559 -> tzz567, al*tzz567 -> tx477, tx363*tzz559 -> tzz566, 
   al*tzz566 -> tx478, tx473*tzz555 -> tx479, tx474*tzz555 -> tx480, 
   tx475*tzz555 -> tx481, tx478 + tx479 -> tx482, tx476 + tx480 -> tx483, 
   tx477 + tx481 -> tx484, tx482*tzz561 -> tx485, tx483*tzz563 -> tx486, 
   tx484*tzz562 -> tx487, tx485 + tx486 + tx487 -> tx488, 
   tx488*tzz558 -> gal, -gal -> fal, 
   CCode["STAPLE_FORCE_ACCUMULATE(I2, 0, fal );"], cl*pxl -> tx489, 
   al*pyl -> tx490, tzz581 -> tx491, -2.*tzz573 -> tx492, dl*pzl -> tx493, 
   al*tx334 -> tx494, -2.*tx467 -> tx495, tx466 + tx490 + tx492 -> tx496, 
   tx491 + tx493 -> tx497, tx489 + tx494 + tx495 -> tx498, 
   bl*tzz568 -> tx499, bl*tzz567 -> tx500, bl*tzz566 -> tx501, 
   tx496*tzz555 -> tx502, tx497*tzz555 -> tx503, tx498*tzz555 -> tx504, 
   tx501 + tx502 -> tx505, tx499 + tx503 -> tx506, tx500 + tx504 -> tx507, 
   tx505*tzz561 -> tx508, tx506*tzz563 -> tx509, tx507*tzz562 -> tx510, 
   tx508 + tx509 + tx510 -> tx511, tx511*tzz558 -> gbl, -gbl -> fbl, 
   CCode["STAPLE_FORCE_ACCUMULATE(I2, 1, fbl );"], bl*pxl -> tx512, 
   dl*pyl -> tx513, al*pzl -> tx514, al*tx332 -> tx515, -2.*tx468 -> tx516, 
   -2.*tx489 -> tx517, tx493 + tx512 -> tx518, tx513 + tx515 + tx516 -> 
    tx519, tx467 + tx514 + tx517 -> tx520, cl*tzz568 -> tx521, 
   cl*tzz567 -> tx522, cl*tzz566 -> tx523, tx518*tzz555 -> tx524, 
   tx519*tzz555 -> tx525, tx520*tzz555 -> tx526, tx522 + tx524 -> tx527, 
   tx523 + tx525 -> tx528, tx521 + tx526 -> tx529, tx527*tzz562 -> tx530, 
   tx528*tzz561 -> tx531, tx529*tzz563 -> tx532, 
   tx530 + tx531 + tx532 -> tx533, tx533*tzz558 -> gcl, -gcl -> fcl, 
   CCode["STAPLE_FORCE_ACCUMULATE(I2, 2, fcl );"], al*pxl -> tx534, 
   tzz573 -> tx535, al*tx333 -> tx536, -2.*tx466 -> tx537, 
   -2.*tx513 -> tx538, tx491 + tx512 -> tx539, tx535 + tx536 + tx537 -> 
    tx540, tx468 + tx534 + tx538 -> tx541, dl*tzz568 -> tx542, 
   dl*tzz567 -> tx543, dl*tzz566 -> tx544, tx539*tzz555 -> tx545, 
   tx540*tzz555 -> tx546, tx541*tzz555 -> tx547, tx544 + tx545 -> tx548, 
   tx542 + tx546 -> tx549, tx543 + tx547 -> tx550, tx548*tzz561 -> tx551, 
   tx549*tzz563 -> tx552, tx550*tzz562 -> tx553, 
   tx551 + tx552 + tx553 -> tx554, tx554*tzz558 -> gdl, -gdl -> fdl, 
   CCode["STAPLE_FORCE_ACCUMULATE(I2, 3, fdl );"], -2.*tzz558 -> tzz576, 
   tx377*tzz576 -> gxl, -gxl -> fxl, 
   CCode["STAPLE_FORCE_ACCUMULATE(I2, 4, fxl );"], tx378*tzz576 -> gyl, 
   -gyl -> fyl, CCode["STAPLE_FORCE_ACCUMULATE(I2, 5, fyl );"], 
   tx379*tzz576 -> gzl, -gzl -> fzl, 
   CCode["STAPLE_FORCE_ACCUMULATE(I2, 6, fzl );"], CCode["} /*calcForce */"], 
   CCode["#endif /* STAPLE_CALC_FORCE ]*/"]}, 
 HessianStructure -> HessianStructure, CCode -> "NOT USED", 
 MathCode -> "Block[{Energy,fak,fal,fbk,fbl,fck,fcl,fdk,fdl,fxk,fxl,fyk,fyl,f\
zk,fzl,gak,gal,gbk,gbl,gck,gcl,gdk,gdl,gxk,gxl,gyk,gyl,gzk,gzl,tx278,tx279,tx\
280,tx281,tx282,tx283,tx284,tx285,tx286,tx287,tx288,tx289,tx290,tx291,tx292,t\
x293,tx294,tx295,tx296,tx297,tx298,tx299,tx300,tx301,tx302,tx303,tx304,tx305,\
tx306,tx307,tx308,tx309,tx310,tx311,tx312,tx313,tx314,tx315,tx316,tx317,tx318\
,tx319,tx320,tx321,tx322,tx323,tx324,tx325,tx326,tx327,tx328,tx329,tx330,tx33\
1,tx332,tx333,tx334,tx335,tx336,tx337,tx338,tx339,tx340,tx341,tx342,tx343,tx3\
44,tx345,tx346,tx347,tx348,tx349,tx350,tx351,tx352,tx353,tx354,tx355,tx356,tx\
357,tx358,tx359,tx360,tx361,tx362,tx363,tx364,tx365,tx366,tx367,tx368,tx369,t\
x370,tx371,tx372,tx373,tx374,tx375,tx376,tx377,tx378,tx379,tx380,tx381,tx382,\
tx383,tx384,tx385,tx386,tx387,tx388,tx389,tx390,tx391,tx392,tx393,tx394,tx395\
,tx396,tx397,tx398,tx399,tx400,tx401,tx402,tx403,tx404,tx405,tx406,tx407,tx40\
8,tx409,tx410,tx411,tx412,tx413,tx414,tx415,tx416,tx417,tx418,tx419,tx420,tx4\
21,tx422,tx423,tx424,tx425,tx426,tx427,tx428,tx429,tx430,tx431,tx432,tx433,tx\
434,tx435,tx436,tx437,tx438,tx439,tx440,tx441,tx442,tx443,tx444,tx445,tx446,t\
x447,tx448,tx449,tx450,tx451,tx452,tx453,tx454,tx455,tx456,tx457,tx458,tx459,\
tx460,tx461,tx462,tx463,tx464,tx465,tx466,tx467,tx468,tx469,tx470,tx471,tx472\
,tx473,tx474,tx475,tx476,tx477,tx478,tx479,tx480,tx481,tx482,tx483,tx484,tx48\
5,tx486,tx487,tx488,tx489,tx490,tx491,tx492,tx493,tx494,tx495,tx496,tx497,tx4\
98,tx499,tx500,tx501,tx502,tx503,tx504,tx505,tx506,tx507,tx508,tx509,tx510,tx\
511,tx512,tx513,tx514,tx515,tx516,tx517,tx518,tx519,tx520,tx521,tx522,tx523,t\
x524,tx525,tx526,tx527,tx528,tx529,tx530,tx531,tx532,tx533,tx534,tx535,tx536,\
tx537,tx538,tx539,tx540,tx541,tx542,tx543,tx544,tx545,tx546,tx547,tx548,tx549\
,tx550,tx551,tx552,tx553,tx554,tzz555,tzz556,tzz557,tzz558,tzz559,tzz560,tzz5\
61,tzz562,tzz563,tzz564,tzz565,tzz566,tzz567,tzz568,tzz569,tzz570,tzz571,tzz5\
72,tzz573,tzz574,tzz575,tzz576,tzz577,tzz578,tzz579,tzz580,tzz581,tzz582,tzz5\
83,tzz584,tzz585,tzz586,tzz587,tzz588,tzz589,tzz590,tzz591,tzz592,tzz593,tzz5\
94,tzz595,tzz596,tzz597,tzz598,tzz599,tzz600,tzz601,tzz602,tzz603,tzz604,tzz6\
05,tzz606,xxxDummy},\n\t(*STAPLE_SET_PARAMETER(ks,ks);*)\n\t(*STAPLE_SET_PARA\
METER(r0,r0);*)\n\t(*STAPLE_SET_PARAMETER(I1,rigidBodyK);*)\n\t(*STAPLE_SET_P\
ARAMETER(I2,rigidBodyL);*)\n\t(*STAPLE_SET_POSITION(ak,I1,0);*)\n\t(*STAPLE_S\
ET_POSITION(bk,I1,1);*)\n\t(*STAPLE_SET_POSITION(ck,I1,2);*)\n\t(*STAPLE_SET_\
POSITION(dk,I1,3);*)\n\t(*STAPLE_SET_POSITION(xk,I1,4);*)\n\t(*STAPLE_SET_POS\
ITION(yk,I1,5);*)\n\t(*STAPLE_SET_POSITION(zk,I1,6);*)\n\t(*STAPLE_SET_POSITI\
ON(al,I2,0);*)\n\t(*STAPLE_SET_POSITION(bl,I2,1);*)\n\t(*STAPLE_SET_POSITION(\
cl,I2,2);*)\n\t(*STAPLE_SET_POSITION(dl,I2,3);*)\n\t(*STAPLE_SET_POSITION(xl,\
I2,4);*)\n\t(*STAPLE_SET_POSITION(yl,I2,5);*)\n\t(*STAPLE_SET_POSITION(zl,I2,\
6);*)\n\t(*STAPLE_SET_POINT(pxk,pointK,getX());*)\n\t(*STAPLE_SET_POINT(pyk,p\
ointK,getY());*)\n\t(*STAPLE_SET_POINT(pzk,pointK,getZ());*)\n\t(*STAPLE_SET_\
POINT(pxl,pointL,getX());*)\n\t(*STAPLE_SET_POINT(pyl,pointL,getY());*)\n\t(*\
STAPLE_SET_POINT(pzl,pointL,getZ());*)\n\ttx278= mathPower2[ak]; (*Rule \
25*)\n\ttx279= mathPower2[al]; (*Rule 26*)\n\ttx280= mathPower2[bk]; (*Rule \
27*)\n\ttx281= mathPower2[bl]; (*Rule 28*)\n\ttx282= mathPower2[ck]; (*Rule \
29*)\n\ttx283= mathPower2[cl]; (*Rule 30*)\n\ttx284= mathPower2[dk]; (*Rule \
31*)\n\ttx285= mathPower2[dl]; (*Rule 32*)\n\ttzz565=-pxk; (*Rule \
33*)\n\ttzz583=ak ck; (*Rule 34*)\n\ttx286=tzz565 tzz583; (*Rule \
35*)\n\ttx287=bk ck pxk; (*Rule 36*)\n\ttzz579=dk pxk; (*Rule \
37*)\n\ttx288=ak tzz579; (*Rule 38*)\n\ttx289=bk tzz579; (*Rule \
39*)\n\ttzz575=-pxl; (*Rule 40*)\n\ttx290=al cl tzz575; (*Rule \
41*)\n\ttx291=bl cl pxl; (*Rule 42*)\n\ttzz577=dl pxl; (*Rule \
43*)\n\ttx292=al tzz577; (*Rule 44*)\n\ttx293=bl tzz577; (*Rule \
45*)\n\ttzz584=ak bk; (*Rule 46*)\n\ttx294=pyk tzz584; (*Rule \
47*)\n\ttzz582=ck pyk; (*Rule 48*)\n\ttx295=bk tzz582; (*Rule \
49*)\n\ttzz572=-pyk; (*Rule 50*)\n\ttx296=ak dk tzz572; (*Rule \
51*)\n\ttx297=dk tzz582; (*Rule 52*)\n\ttx298=al bl pyl; (*Rule \
53*)\n\ttzz581=cl pyl; (*Rule 54*)\n\ttx299=bl tzz581; (*Rule \
55*)\n\ttzz574=-pyl; (*Rule 56*)\n\ttx300=al dl tzz574; (*Rule \
57*)\n\ttx301=dl tzz581; (*Rule 58*)\n\ttzz564=-pzk; (*Rule \
59*)\n\ttx302=tzz564 tzz584; (*Rule 60*)\n\ttx303=pzk tzz583; (*Rule \
61*)\n\ttzz578=dk pzk; (*Rule 62*)\n\ttx304=bk tzz578; (*Rule \
63*)\n\ttx305=ck tzz578; (*Rule 64*)\n\ttzz573=bl pzl; (*Rule \
65*)\n\ttx306=-(al tzz573); (*Rule 66*)\n\ttzz580=cl pzl; (*Rule \
67*)\n\ttx307=al tzz580; (*Rule 68*)\n\ttx308=dl tzz573; (*Rule \
69*)\n\ttx309=dl tzz580; (*Rule 70*)\n\ttx310=tx280 tzz572; (*Rule \
71*)\n\ttx311=tx280 tzz564; (*Rule 72*)\n\ttx312=tx281 tzz574; (*Rule \
73*)\n\ttzz585=-pzl; (*Rule 74*)\n\ttx313=tx281 tzz585; (*Rule \
75*)\n\ttx314=tx282 tzz565; (*Rule 76*)\n\ttx315=tx282 tzz564; (*Rule \
77*)\n\ttx316=tx283 tzz575; (*Rule 78*)\n\ttx317=tx283 tzz585; (*Rule \
79*)\n\ttx318=tx284 tzz565; (*Rule 80*)\n\ttx319=tx284 tzz572; (*Rule \
81*)\n\ttx320=tx278 + tx280 + tx282 + tx284; (*Rule 82*)\n\ttx321=tx285 \
tzz575; (*Rule 83*)\n\ttx322=tx285 tzz574; (*Rule 84*)\n\ttx323=tx279 + tx281 \
+ tx283 + tx285; (*Rule 85*)\n\ttx324=tx286 + tx289 + tx294 + tx297 + tx311 + \
tx315; (*Rule 86*)\n\ttzz599=tx298 + tx301; (*Rule 87*)\n\ttzz600=tx293 + \
tzz599; (*Rule 88*)\n\ttzz603=tx290 + tzz600; (*Rule 89*)\n\ttx325=tx313 + \
tx317 + tzz603; (*Rule 90*)\n\ttx326=tx295 + tx296 + tx303 + tx304 + tx314 + \
tx318; (*Rule 91*)\n\ttx327=tx287 + tx288 + tx302 + tx305 + tx310 + tx319; \
(*Rule 92*)\n\ttx328= mathReciprocal[tx320]; (*Rule 93*)\n\ttzz595=tx307 + \
tx308; (*Rule 94*)\n\ttzz597=tx300 + tzz595; (*Rule 95*)\n\ttzz598=tx299 + \
tzz597; (*Rule 96*)\n\ttx329=tx316 + tx321 + tzz598; (*Rule \
97*)\n\ttzz596=tx306 + tx309; (*Rule 98*)\n\ttzz601=tx292 + tzz596; (*Rule \
99*)\n\ttzz602=tx291 + tzz601; (*Rule 100*)\n\ttx330=tx312 + tx322 + tzz602; \
(*Rule 101*)\n\ttx331= mathReciprocal[tx323]; (*Rule 102*)\n\ttx332=tzz575; \
(*Rule 103*)\n\ttx333=tzz574; (*Rule 104*)\n\ttx334=tzz585; (*Rule \
105*)\n\ttzz556=2. tx328; (*Rule 106*)\n\ttx335=tx324 tzz556; (*Rule \
107*)\n\ttx336=tx326 tzz556; (*Rule 108*)\n\ttx337=tx327 tzz556; (*Rule \
109*)\n\ttzz555=-2. tx331; (*Rule 110*)\n\ttx338=tx325 tzz555; (*Rule \
111*)\n\ttx339=tx329 tzz555; (*Rule 112*)\n\ttx340=tx330 tzz555; (*Rule \
113*)\n\ttx341=-xl; (*Rule 114*)\n\ttx342=-yl; (*Rule 115*)\n\ttx343=-zl; \
(*Rule 116*)\n\ttzz588=tx341 + xk; (*Rule 117*)\n\ttzz590=tx336 + tzz588; \
(*Rule 118*)\n\ttzz594=tx332 + tzz590; (*Rule 119*)\n\ttzz606=pxk + tzz594; \
(*Rule 120*)\n\ttx344=tx339 + tzz606; (*Rule 121*)\n\ttzz587=tx342 + yk; \
(*Rule 122*)\n\ttzz589=tx337 + tzz587; (*Rule 123*)\n\ttzz593=tx333 + tzz589; \
(*Rule 124*)\n\ttzz605=pyk + tzz593; (*Rule 125*)\n\ttx345=tx340 + tzz605; \
(*Rule 126*)\n\ttzz586=tx343 + zk; (*Rule 127*)\n\ttzz591=tx335 + tzz586; \
(*Rule 128*)\n\ttzz592=tx334 + tzz591; (*Rule 129*)\n\ttzz604=pzk + tzz592; \
(*Rule 130*)\n\ttx346=tx338 + tzz604; (*Rule 131*)\n\ttx347= \
mathPower2[tx344]; (*Rule 132*)\n\ttx348= mathPower2[tx345]; (*Rule \
133*)\n\ttx349= mathPower2[tx346]; (*Rule 134*)\n\ttx350=tx347 + tx348 + \
tx349; (*Rule 135*)\n\ttx351=-r0; (*Rule 136*)\n\ttx352= mathSqrt[tx350]; \
(*Rule 137*)\n\ttx353=tx351 + tx352; (*Rule 138*)\n\ttx354= \
mathPower2[tx353]; (*Rule 139*)\n\tEnergy=ks tx354; (*Rule \
140*)\n\t(*STAPLE_ENERGY_ACCUMULATE(Energy);*)\n\t(*#ifdef STAPLE_CALC_FORCE \
//[*)\n\t(*if ( calcForce ) {*)\n\ttx355=tx283 tx332; (*Rule \
144*)\n\ttx356=tx285 tx332; (*Rule 145*)\n\ttx357=tx281 tx333; (*Rule \
146*)\n\ttx358=tx285 tx333; (*Rule 147*)\n\ttx359=tx281 tx334; (*Rule \
148*)\n\ttx360=tx283 tx334; (*Rule 149*)\n\ttx361=tx355 + tx356 + tzz598; \
(*Rule 150*)\n\ttx362=tx357 + tx358 + tzz602; (*Rule 151*)\n\ttx363=tx359 + \
tx360 + tzz603; (*Rule 152*)\n\ttx364=ck tzz565; (*Rule \
153*)\n\ttx365=tzz579; (*Rule 154*)\n\ttx366=bk pyk; (*Rule 155*)\n\ttx367=dk \
tzz572; (*Rule 156*)\n\ttx368=bk tzz564; (*Rule 157*)\n\ttx369=ck pzk; (*Rule \
158*)\n\ttx370=tx361 tzz555; (*Rule 159*)\n\ttx371=tx362 tzz555; (*Rule \
160*)\n\ttx372=tx363 tzz555; (*Rule 161*)\n\ttx373= mathPower2[tx328]; (*Rule \
162*)\n\ttx374=tx364 + tx366; (*Rule 163*)\n\ttx375=tx365 + tx368; (*Rule \
164*)\n\ttx376=tx367 + tx369; (*Rule 165*)\n\ttx377=tx370 + tzz606; (*Rule \
166*)\n\ttx378=tx371 + tzz605; (*Rule 167*)\n\ttx379=tx372 + tzz604; (*Rule \
168*)\n\ttzz560=-4. tx373; (*Rule 169*)\n\ttzz571=tx324 tzz560; (*Rule \
170*)\n\ttx380=ak tzz571; (*Rule 171*)\n\ttzz570=tx326 tzz560; (*Rule \
172*)\n\ttx381=ak tzz570; (*Rule 173*)\n\ttzz569=tx327 tzz560; (*Rule \
174*)\n\ttx382=ak tzz569; (*Rule 175*)\n\ttx383=tx374 tzz556; (*Rule \
176*)\n\ttx384=tx375 tzz556; (*Rule 177*)\n\ttx385=tx376 tzz556; (*Rule \
178*)\n\ttx386= mathPower2[tx377]; (*Rule 179*)\n\ttx387= mathPower2[tx378]; \
(*Rule 180*)\n\ttx388= mathPower2[tx379]; (*Rule 181*)\n\ttx389=tx380 + \
tx383; (*Rule 182*)\n\ttx390=tx382 + tx384; (*Rule 183*)\n\ttx391=tx381 + \
tx385; (*Rule 184*)\n\ttx392=tx386 + tx387 + tx388; (*Rule 185*)\n\ttzz561=2. \
tx379; (*Rule 186*)\n\ttx393=tx389 tzz561; (*Rule 187*)\n\ttzz562=2. tx378; \
(*Rule 188*)\n\ttx394=tx390 tzz562; (*Rule 189*)\n\ttzz563=2. tx377; (*Rule \
190*)\n\ttx395=tx391 tzz563; (*Rule 191*)\n\ttx396= mathSqrt[tx392]; (*Rule \
192*)\n\ttx397= mathReciprocal[tx396]; (*Rule 193*)\n\ttx398=tx393 + tx394 + \
tx395; (*Rule 194*)\n\ttx399=tx351 + tx396; (*Rule 195*)\n\ttzz557=tx397 \
tx399; (*Rule 196*)\n\ttzz558=ks tzz557; (*Rule 197*)\n\tgak=tx398 tzz558; \
(*Rule 198*)\n\tfak=-gak; (*Rule 199*)\n\t(*STAPLE_FORCE_ACCUMULATE(I1, 0, \
fak );*)\n\ttx400=ck pxk; (*Rule 201*)\n\ttx401=ak pyk; (*Rule \
202*)\n\ttx402=tzz582; (*Rule 203*)\n\ttx403=ak tzz564; (*Rule \
204*)\n\ttx404=-2. bk pzk; (*Rule 205*)\n\ttx405=tzz578; (*Rule \
206*)\n\ttx406=-2. tx366; (*Rule 207*)\n\ttx407=tx365 + tx401 + tx404; (*Rule \
208*)\n\ttx408=tx402 + tx405; (*Rule 209*)\n\ttx409=tx400 + tx403 + tx406; \
(*Rule 210*)\n\ttx410=bk tzz571; (*Rule 211*)\n\ttx411=bk tzz570; (*Rule \
212*)\n\ttx412=bk tzz569; (*Rule 213*)\n\ttx413=tx407 tzz556; (*Rule \
214*)\n\ttx414=tx408 tzz556; (*Rule 215*)\n\ttx415=tx409 tzz556; (*Rule \
216*)\n\ttx416=tx410 + tx413; (*Rule 217*)\n\ttx417=tx411 + tx414; (*Rule \
218*)\n\ttx418=tx412 + tx415; (*Rule 219*)\n\ttx419=tx416 tzz561; (*Rule \
220*)\n\ttx420=tx417 tzz563; (*Rule 221*)\n\ttx421=tx418 tzz562; (*Rule \
222*)\n\ttx422=tx419 + tx420 + tx421; (*Rule 223*)\n\tgbk=tx422 tzz558; \
(*Rule 224*)\n\tfbk=-gbk; (*Rule 225*)\n\t(*STAPLE_FORCE_ACCUMULATE(I1, 1, \
fbk );*)\n\ttx423=ak tzz565; (*Rule 227*)\n\ttx424=bk pxk; (*Rule \
228*)\n\ttx425=dk pyk; (*Rule 229*)\n\ttx426=ak pzk; (*Rule \
230*)\n\ttx427=-2. tx369; (*Rule 231*)\n\ttx428=-2. tx400; (*Rule \
232*)\n\ttx429=tx405 + tx424; (*Rule 233*)\n\ttx430=tx423 + tx425 + tx427; \
(*Rule 234*)\n\ttx431=tx366 + tx426 + tx428; (*Rule 235*)\n\ttx432=ck tzz571; \
(*Rule 236*)\n\ttx433=ck tzz570; (*Rule 237*)\n\ttx434=ck tzz569; (*Rule \
238*)\n\ttx435=tx429 tzz556; (*Rule 239*)\n\ttx436=tx430 tzz556; (*Rule \
240*)\n\ttx437=tx431 tzz556; (*Rule 241*)\n\ttx438=tx434 + tx435; (*Rule \
242*)\n\ttx439=tx432 + tx436; (*Rule 243*)\n\ttx440=tx433 + tx437; (*Rule \
244*)\n\ttx441=tx438 tzz562; (*Rule 245*)\n\ttx442=tx439 tzz561; (*Rule \
246*)\n\ttx443=tx440 tzz563; (*Rule 247*)\n\ttx444=tx441 + tx442 + tx443; \
(*Rule 248*)\n\tgck=tx444 tzz558; (*Rule 249*)\n\tfck=-gck; (*Rule \
250*)\n\t(*STAPLE_FORCE_ACCUMULATE(I1, 2, fck );*)\n\ttx445=ak pxk; (*Rule \
252*)\n\ttx446=bk pzk; (*Rule 253*)\n\ttx447=-2. tx365; (*Rule \
254*)\n\ttx448=-tx401; (*Rule 255*)\n\ttx449=-2. tx425; (*Rule \
256*)\n\ttx450=tx402 + tx424; (*Rule 257*)\n\ttx451=tx446 + tx447 + tx448; \
(*Rule 258*)\n\ttx452=tx369 + tx445 + tx449; (*Rule 259*)\n\ttx453=dk tzz571; \
(*Rule 260*)\n\ttx454=dk tzz570; (*Rule 261*)\n\ttx455=dk tzz569; (*Rule \
262*)\n\ttx456=tx450 tzz556; (*Rule 263*)\n\ttx457=tx451 tzz556; (*Rule \
264*)\n\ttx458=tx452 tzz556; (*Rule 265*)\n\ttx459=tx453 + tx456; (*Rule \
266*)\n\ttx460=tx454 + tx457; (*Rule 267*)\n\ttx461=tx455 + tx458; (*Rule \
268*)\n\ttx462=tx459 tzz561; (*Rule 269*)\n\ttx463=tx460 tzz563; (*Rule \
270*)\n\ttx464=tx461 tzz562; (*Rule 271*)\n\ttx465=tx462 + tx463 + tx464; \
(*Rule 272*)\n\tgdk=tx465 tzz558; (*Rule 273*)\n\tfdk=-gdk; (*Rule \
274*)\n\t(*STAPLE_FORCE_ACCUMULATE(I1, 3, fdk );*)\n\tgxk=tzz558 tzz563; \
(*Rule 276*)\n\tfxk=-gxk; (*Rule 277*)\n\t(*STAPLE_FORCE_ACCUMULATE(I1, 4, \
fxk );*)\n\tgyk=tzz558 tzz562; (*Rule 279*)\n\tfyk=-gyk; (*Rule \
280*)\n\t(*STAPLE_FORCE_ACCUMULATE(I1, 5, fyk );*)\n\tgzk=tzz558 tzz561; \
(*Rule 282*)\n\tfzk=-gzk; (*Rule 283*)\n\t(*STAPLE_FORCE_ACCUMULATE(I1, 6, \
fzk );*)\n\ttx466=tzz577; (*Rule 285*)\n\ttx467=bl pyl; (*Rule \
286*)\n\ttx468=tzz580; (*Rule 287*)\n\ttx469=cl tx332; (*Rule \
288*)\n\ttx470=dl tx333; (*Rule 289*)\n\ttx471=bl tx334; (*Rule \
290*)\n\ttx472= mathPower2[tx331]; (*Rule 291*)\n\ttx473=tx467 + tx469; \
(*Rule 292*)\n\ttx474=tx468 + tx470; (*Rule 293*)\n\ttx475=tx466 + tx471; \
(*Rule 294*)\n\ttzz559=4. tx472; (*Rule 295*)\n\ttzz568=tx361 tzz559; (*Rule \
296*)\n\ttx476=al tzz568; (*Rule 297*)\n\ttzz567=tx362 tzz559; (*Rule \
298*)\n\ttx477=al tzz567; (*Rule 299*)\n\ttzz566=tx363 tzz559; (*Rule \
300*)\n\ttx478=al tzz566; (*Rule 301*)\n\ttx479=tx473 tzz555; (*Rule \
302*)\n\ttx480=tx474 tzz555; (*Rule 303*)\n\ttx481=tx475 tzz555; (*Rule \
304*)\n\ttx482=tx478 + tx479; (*Rule 305*)\n\ttx483=tx476 + tx480; (*Rule \
306*)\n\ttx484=tx477 + tx481; (*Rule 307*)\n\ttx485=tx482 tzz561; (*Rule \
308*)\n\ttx486=tx483 tzz563; (*Rule 309*)\n\ttx487=tx484 tzz562; (*Rule \
310*)\n\ttx488=tx485 + tx486 + tx487; (*Rule 311*)\n\tgal=tx488 tzz558; \
(*Rule 312*)\n\tfal=-gal; (*Rule 313*)\n\t(*STAPLE_FORCE_ACCUMULATE(I2, 0, \
fal );*)\n\ttx489=cl pxl; (*Rule 315*)\n\ttx490=al pyl; (*Rule \
316*)\n\ttx491=tzz581; (*Rule 317*)\n\ttx492=-2. tzz573; (*Rule \
318*)\n\ttx493=dl pzl; (*Rule 319*)\n\ttx494=al tx334; (*Rule \
320*)\n\ttx495=-2. tx467; (*Rule 321*)\n\ttx496=tx466 + tx490 + tx492; (*Rule \
322*)\n\ttx497=tx491 + tx493; (*Rule 323*)\n\ttx498=tx489 + tx494 + tx495; \
(*Rule 324*)\n\ttx499=bl tzz568; (*Rule 325*)\n\ttx500=bl tzz567; (*Rule \
326*)\n\ttx501=bl tzz566; (*Rule 327*)\n\ttx502=tx496 tzz555; (*Rule \
328*)\n\ttx503=tx497 tzz555; (*Rule 329*)\n\ttx504=tx498 tzz555; (*Rule \
330*)\n\ttx505=tx501 + tx502; (*Rule 331*)\n\ttx506=tx499 + tx503; (*Rule \
332*)\n\ttx507=tx500 + tx504; (*Rule 333*)\n\ttx508=tx505 tzz561; (*Rule \
334*)\n\ttx509=tx506 tzz563; (*Rule 335*)\n\ttx510=tx507 tzz562; (*Rule \
336*)\n\ttx511=tx508 + tx509 + tx510; (*Rule 337*)\n\tgbl=tx511 tzz558; \
(*Rule 338*)\n\tfbl=-gbl; (*Rule 339*)\n\t(*STAPLE_FORCE_ACCUMULATE(I2, 1, \
fbl );*)\n\ttx512=bl pxl; (*Rule 341*)\n\ttx513=dl pyl; (*Rule \
342*)\n\ttx514=al pzl; (*Rule 343*)\n\ttx515=al tx332; (*Rule \
344*)\n\ttx516=-2. tx468; (*Rule 345*)\n\ttx517=-2. tx489; (*Rule \
346*)\n\ttx518=tx493 + tx512; (*Rule 347*)\n\ttx519=tx513 + tx515 + tx516; \
(*Rule 348*)\n\ttx520=tx467 + tx514 + tx517; (*Rule 349*)\n\ttx521=cl tzz568; \
(*Rule 350*)\n\ttx522=cl tzz567; (*Rule 351*)\n\ttx523=cl tzz566; (*Rule \
352*)\n\ttx524=tx518 tzz555; (*Rule 353*)\n\ttx525=tx519 tzz555; (*Rule \
354*)\n\ttx526=tx520 tzz555; (*Rule 355*)\n\ttx527=tx522 + tx524; (*Rule \
356*)\n\ttx528=tx523 + tx525; (*Rule 357*)\n\ttx529=tx521 + tx526; (*Rule \
358*)\n\ttx530=tx527 tzz562; (*Rule 359*)\n\ttx531=tx528 tzz561; (*Rule \
360*)\n\ttx532=tx529 tzz563; (*Rule 361*)\n\ttx533=tx530 + tx531 + tx532; \
(*Rule 362*)\n\tgcl=tx533 tzz558; (*Rule 363*)\n\tfcl=-gcl; (*Rule \
364*)\n\t(*STAPLE_FORCE_ACCUMULATE(I2, 2, fcl );*)\n\ttx534=al pxl; (*Rule \
366*)\n\ttx535=tzz573; (*Rule 367*)\n\ttx536=al tx333; (*Rule \
368*)\n\ttx537=-2. tx466; (*Rule 369*)\n\ttx538=-2. tx513; (*Rule \
370*)\n\ttx539=tx491 + tx512; (*Rule 371*)\n\ttx540=tx535 + tx536 + tx537; \
(*Rule 372*)\n\ttx541=tx468 + tx534 + tx538; (*Rule 373*)\n\ttx542=dl tzz568; \
(*Rule 374*)\n\ttx543=dl tzz567; (*Rule 375*)\n\ttx544=dl tzz566; (*Rule \
376*)\n\ttx545=tx539 tzz555; (*Rule 377*)\n\ttx546=tx540 tzz555; (*Rule \
378*)\n\ttx547=tx541 tzz555; (*Rule 379*)\n\ttx548=tx544 + tx545; (*Rule \
380*)\n\ttx549=tx542 + tx546; (*Rule 381*)\n\ttx550=tx543 + tx547; (*Rule \
382*)\n\ttx551=tx548 tzz561; (*Rule 383*)\n\ttx552=tx549 tzz563; (*Rule \
384*)\n\ttx553=tx550 tzz562; (*Rule 385*)\n\ttx554=tx551 + tx552 + tx553; \
(*Rule 386*)\n\tgdl=tx554 tzz558; (*Rule 387*)\n\tfdl=-gdl; (*Rule \
388*)\n\t(*STAPLE_FORCE_ACCUMULATE(I2, 3, fdl );*)\n\ttzz576=-2. tzz558; \
(*Rule 390*)\n\tgxl=tx377 tzz576; (*Rule 391*)\n\tfxl=-gxl; (*Rule \
392*)\n\t(*STAPLE_FORCE_ACCUMULATE(I2, 4, fxl );*)\n\tgyl=tx378 tzz576; \
(*Rule 394*)\n\tfyl=-gyl; (*Rule 395*)\n\t(*STAPLE_FORCE_ACCUMULATE(I2, 5, \
fyl );*)\n\tgzl=tx379 tzz576; (*Rule 397*)\n\tfzl=-gzl; (*Rule \
398*)\n\t(*STAPLE_FORCE_ACCUMULATE(I2, 6, fzl );*)\n\t(*} /*calcForce \
*/*)\n\t(*#endif /* STAPLE_CALC_FORCE ]*/*)\n{Energy, fak, fbk, fck, fdk, \
fxk, fyk, fzk, fal, fbl, fcl, fdl, fxl, fyl, fzl}]\n"}
