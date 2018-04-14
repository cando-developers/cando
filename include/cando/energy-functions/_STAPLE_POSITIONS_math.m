{Name -> "STAPLE_POSITIONS", AdditionalCDeclares -> "", 
 Input -> staplePositionInput, Output -> {plabkx, plabky, plabkz, plablx, 
   plably, plablz}, EnergyFunction -> EnergyFunction, 
 DerivativeVariables -> staplePositionInput, 
 Rules -> {CCode["STAPLE_SET_PARAMETER(ks,ks);"], 
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
   CCode["STAPLE_SET_POINT(pzl,pointL,getZ());"], power2[ak] -> tx61, 
   power2[bk] -> tx62, power2[ck] -> tx63, power2[dk] -> tx64, 
   ck*pyk -> tzz134, bk*tzz134 -> tx65, -pyk -> tzz126, ak*dk*tzz126 -> tx66, 
   ak*ck -> tzz137, pzk*tzz137 -> tx67, dk*pzk -> tzz131, bk*tzz131 -> tx68, 
   -pxk -> tzz128, tx63*tzz128 -> tx69, tx64*tzz128 -> tx70, 
   tx61 + tx62 + tx63 + tx64 -> tx71, 
   tx65 + tx66 + tx67 + tx68 + tx69 + tx70 -> tx72, reciprocal[tx71] -> tx73, 
   2.*tx73 -> tzz121, tx72*tzz121 -> tx74, pxk + tx74 + xk -> plabkx, 
   bk*ck*pxk -> tx75, dk*pxk -> tzz132, ak*tzz132 -> tx76, -pzk -> tzz124, 
   ak*bk -> tzz138, tzz124*tzz138 -> tx77, ck*tzz131 -> tx78, 
   tx62*tzz126 -> tx79, tx64*tzz126 -> tx80, 
   tx75 + tx76 + tx77 + tx78 + tx79 + tx80 -> tx81, tx81*tzz121 -> tx82, 
   pyk + tx82 + yk -> plabky, tzz128*tzz137 -> tx83, bk*tzz132 -> tx84, 
   pyk*tzz138 -> tx85, dk*tzz134 -> tx86, tx62*tzz124 -> tx87, 
   tx63*tzz124 -> tx88, tx83 + tx84 + tx85 + tx86 + tx87 + tx88 -> tx89, 
   tx89*tzz121 -> tx90, pzk + tx90 + zk -> plabkz, power2[al] -> tx91, 
   power2[bl] -> tx92, power2[cl] -> tx93, power2[dl] -> tx94, 
   cl*pyl -> tzz133, bl*tzz133 -> tx95, -pyl -> tzz125, al*dl*tzz125 -> tx96, 
   al*cl -> tzz135, pzl*tzz135 -> tx97, dl*pzl -> tzz129, bl*tzz129 -> tx98, 
   -pxl -> tzz127, tx93*tzz127 -> tx99, tx94*tzz127 -> tx100, 
   tx91 + tx92 + tx93 + tx94 -> tx101, reciprocal[tx101] -> tx102, 
   tx100 + tx95 + tx96 + tx97 + tx98 + tx99 -> tx103, 2.*tx102 -> tzz122, 
   tx103*tzz122 -> tx104, pxl + tx104 + xl -> plablx, bl*cl*pxl -> tx105, 
   dl*pxl -> tzz130, al*tzz130 -> tx106, -pzl -> tzz123, al*bl -> tzz136, 
   tzz123*tzz136 -> tx107, cl*tzz129 -> tx108, tx92*tzz125 -> tx109, 
   tx94*tzz125 -> tx110, tx105 + tx106 + tx107 + tx108 + tx109 + tx110 -> 
    tx111, tx111*tzz122 -> tx112, pyl + tx112 + yl -> plably, 
   tzz127*tzz135 -> tx113, bl*tzz130 -> tx114, pyl*tzz136 -> tx115, 
   dl*tzz133 -> tx116, tx92*tzz123 -> tx117, tx93*tzz123 -> tx118, 
   tx113 + tx114 + tx115 + tx116 + tx117 + tx118 -> tx119, 
   tx119*tzz122 -> tx120, pzl + tx120 + zl -> plablz}, 
 HessianStructure -> HessianStructure, CCode -> "NOT USED", 
 MathCode -> "Block[{plabkx,plabky,plabkz,plablx,plably,plablz,tx100,tx101,tx\
102,tx103,tx104,tx105,tx106,tx107,tx108,tx109,tx110,tx111,tx112,tx113,tx114,t\
x115,tx116,tx117,tx118,tx119,tx120,tx61,tx62,tx63,tx64,tx65,tx66,tx67,tx68,tx\
69,tx70,tx71,tx72,tx73,tx74,tx75,tx76,tx77,tx78,tx79,tx80,tx81,tx82,tx83,tx84\
,tx85,tx86,tx87,tx88,tx89,tx90,tx91,tx92,tx93,tx94,tx95,tx96,tx97,tx98,tx99,t\
zz121,tzz122,tzz123,tzz124,tzz125,tzz126,tzz127,tzz128,tzz129,tzz130,tzz131,t\
zz132,tzz133,tzz134,tzz135,tzz136,tzz137,tzz138,xxxDummy},\n\t(*STAPLE_SET_PA\
RAMETER(ks,ks);*)\n\t(*STAPLE_SET_PARAMETER(r0,r0);*)\n\t(*STAPLE_SET_PARAMET\
ER(I1,rigidBodyK);*)\n\t(*STAPLE_SET_PARAMETER(I2,rigidBodyL);*)\n\t(*STAPLE_\
SET_POSITION(ak,I1,0);*)\n\t(*STAPLE_SET_POSITION(bk,I1,1);*)\n\t(*STAPLE_SET\
_POSITION(ck,I1,2);*)\n\t(*STAPLE_SET_POSITION(dk,I1,3);*)\n\t(*STAPLE_SET_PO\
SITION(xk,I1,4);*)\n\t(*STAPLE_SET_POSITION(yk,I1,5);*)\n\t(*STAPLE_SET_POSIT\
ION(zk,I1,6);*)\n\t(*STAPLE_SET_POSITION(al,I2,0);*)\n\t(*STAPLE_SET_POSITION\
(bl,I2,1);*)\n\t(*STAPLE_SET_POSITION(cl,I2,2);*)\n\t(*STAPLE_SET_POSITION(dl\
,I2,3);*)\n\t(*STAPLE_SET_POSITION(xl,I2,4);*)\n\t(*STAPLE_SET_POSITION(yl,I2\
,5);*)\n\t(*STAPLE_SET_POSITION(zl,I2,6);*)\n\t(*STAPLE_SET_POINT(pxk,pointK,\
getX());*)\n\t(*STAPLE_SET_POINT(pyk,pointK,getY());*)\n\t(*STAPLE_SET_POINT(\
pzk,pointK,getZ());*)\n\t(*STAPLE_SET_POINT(pxl,pointL,getX());*)\n\t(*STAPLE\
_SET_POINT(pyl,pointL,getY());*)\n\t(*STAPLE_SET_POINT(pzl,pointL,getZ());*)\
\n\ttx61= mathPower2[ak]; (*Rule 25*)\n\ttx62= mathPower2[bk]; (*Rule \
26*)\n\ttx63= mathPower2[ck]; (*Rule 27*)\n\ttx64= mathPower2[dk]; (*Rule \
28*)\n\ttzz134=ck pyk; (*Rule 29*)\n\ttx65=bk tzz134; (*Rule \
30*)\n\ttzz126=-pyk; (*Rule 31*)\n\ttx66=ak dk tzz126; (*Rule \
32*)\n\ttzz137=ak ck; (*Rule 33*)\n\ttx67=pzk tzz137; (*Rule \
34*)\n\ttzz131=dk pzk; (*Rule 35*)\n\ttx68=bk tzz131; (*Rule \
36*)\n\ttzz128=-pxk; (*Rule 37*)\n\ttx69=tx63 tzz128; (*Rule \
38*)\n\ttx70=tx64 tzz128; (*Rule 39*)\n\ttx71=tx61 + tx62 + tx63 + tx64; \
(*Rule 40*)\n\ttx72=tx65 + tx66 + tx67 + tx68 + tx69 + tx70; (*Rule \
41*)\n\ttx73= mathReciprocal[tx71]; (*Rule 42*)\n\ttzz121=2. tx73; (*Rule \
43*)\n\ttx74=tx72 tzz121; (*Rule 44*)\n\tplabkx=pxk + tx74 + xk; (*Rule \
45*)\n\ttx75=bk ck pxk; (*Rule 46*)\n\ttzz132=dk pxk; (*Rule 47*)\n\ttx76=ak \
tzz132; (*Rule 48*)\n\ttzz124=-pzk; (*Rule 49*)\n\ttzz138=ak bk; (*Rule \
50*)\n\ttx77=tzz124 tzz138; (*Rule 51*)\n\ttx78=ck tzz131; (*Rule \
52*)\n\ttx79=tx62 tzz126; (*Rule 53*)\n\ttx80=tx64 tzz126; (*Rule \
54*)\n\ttx81=tx75 + tx76 + tx77 + tx78 + tx79 + tx80; (*Rule \
55*)\n\ttx82=tx81 tzz121; (*Rule 56*)\n\tplabky=pyk + tx82 + yk; (*Rule \
57*)\n\ttx83=tzz128 tzz137; (*Rule 58*)\n\ttx84=bk tzz132; (*Rule \
59*)\n\ttx85=pyk tzz138; (*Rule 60*)\n\ttx86=dk tzz134; (*Rule \
61*)\n\ttx87=tx62 tzz124; (*Rule 62*)\n\ttx88=tx63 tzz124; (*Rule \
63*)\n\ttx89=tx83 + tx84 + tx85 + tx86 + tx87 + tx88; (*Rule \
64*)\n\ttx90=tx89 tzz121; (*Rule 65*)\n\tplabkz=pzk + tx90 + zk; (*Rule \
66*)\n\ttx91= mathPower2[al]; (*Rule 67*)\n\ttx92= mathPower2[bl]; (*Rule \
68*)\n\ttx93= mathPower2[cl]; (*Rule 69*)\n\ttx94= mathPower2[dl]; (*Rule \
70*)\n\ttzz133=cl pyl; (*Rule 71*)\n\ttx95=bl tzz133; (*Rule \
72*)\n\ttzz125=-pyl; (*Rule 73*)\n\ttx96=al dl tzz125; (*Rule \
74*)\n\ttzz135=al cl; (*Rule 75*)\n\ttx97=pzl tzz135; (*Rule \
76*)\n\ttzz129=dl pzl; (*Rule 77*)\n\ttx98=bl tzz129; (*Rule \
78*)\n\ttzz127=-pxl; (*Rule 79*)\n\ttx99=tx93 tzz127; (*Rule \
80*)\n\ttx100=tx94 tzz127; (*Rule 81*)\n\ttx101=tx91 + tx92 + tx93 + tx94; \
(*Rule 82*)\n\ttx102= mathReciprocal[tx101]; (*Rule 83*)\n\ttx103=tx100 + \
tx95 + tx96 + tx97 + tx98 + tx99; (*Rule 84*)\n\ttzz122=2. tx102; (*Rule \
85*)\n\ttx104=tx103 tzz122; (*Rule 86*)\n\tplablx=pxl + tx104 + xl; (*Rule \
87*)\n\ttx105=bl cl pxl; (*Rule 88*)\n\ttzz130=dl pxl; (*Rule \
89*)\n\ttx106=al tzz130; (*Rule 90*)\n\ttzz123=-pzl; (*Rule 91*)\n\ttzz136=al \
bl; (*Rule 92*)\n\ttx107=tzz123 tzz136; (*Rule 93*)\n\ttx108=cl tzz129; \
(*Rule 94*)\n\ttx109=tx92 tzz125; (*Rule 95*)\n\ttx110=tx94 tzz125; (*Rule \
96*)\n\ttx111=tx105 + tx106 + tx107 + tx108 + tx109 + tx110; (*Rule \
97*)\n\ttx112=tx111 tzz122; (*Rule 98*)\n\tplably=pyl + tx112 + yl; (*Rule \
99*)\n\ttx113=tzz127 tzz135; (*Rule 100*)\n\ttx114=bl tzz130; (*Rule \
101*)\n\ttx115=pyl tzz136; (*Rule 102*)\n\ttx116=dl tzz133; (*Rule \
103*)\n\ttx117=tx92 tzz123; (*Rule 104*)\n\ttx118=tx93 tzz123; (*Rule \
105*)\n\ttx119=tx113 + tx114 + tx115 + tx116 + tx117 + tx118; (*Rule \
106*)\n\ttx120=tx119 tzz122; (*Rule 107*)\n\tplablz=pzl + tx120 + zl; (*Rule \
108*)\n{plabkx, plabky, plabkz, plablx, plably, plablz}]\n"}
