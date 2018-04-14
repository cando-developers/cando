{Name -> "NONBOND_POSITIONS", AdditionalCDeclares -> "", 
 Input -> nonbondPositionInput, Output -> {plabmx, plabmy, plabmz}, 
 EnergyFunction -> EnergyFunction, DerivativeVariables -> 
  nonbondPositionInput, 
 Rules -> {CCode["NONBOND_POSITION_RB_SET_PARAMETER(I1);"], 
   CCode["NONBOND_POSITION_RB_SET_POSITION(am,I1,0);"], 
   CCode["NONBOND_POSITION_RB_SET_POSITION(bm,I1,1);"], 
   CCode["NONBOND_POSITION_RB_SET_POSITION(cm,I1,2);"], 
   CCode["NONBOND_POSITION_RB_SET_POSITION(dm,I1,3);"], 
   CCode["NONBOND_POSITION_RB_SET_POSITION(xm,I1,4);"], 
   CCode["NONBOND_POSITION_RB_SET_POSITION(ym,I1,5);"], 
   CCode["NONBOND_POSITION_RB_SET_POSITION(zm,I1,6);"], 
   CCode["NONBOND_POSITION_RB_SET_POINT(pxm,ea1,getX());"], 
   CCode["NONBOND_POSITION_RB_SET_POINT(pym,ea1,getY());"], 
   CCode["NONBOND_POSITION_RB_SET_POINT(pzm,ea1,getZ());"], 
   power2[am] -> tx31, power2[bm] -> tx32, power2[cm] -> tx33, 
   power2[dm] -> tx34, cm*pym -> tzz67, bm*tzz67 -> tx35, -pym -> tzz63, 
   am*dm*tzz63 -> tx36, am*cm -> tzz68, pzm*tzz68 -> tx37, dm*pzm -> tzz65, 
   bm*tzz65 -> tx38, -pxm -> tzz64, tx33*tzz64 -> tx39, tx34*tzz64 -> tx40, 
   tx31 + tx32 + tx33 + tx34 -> tx41, 
   tx35 + tx36 + tx37 + tx38 + tx39 + tx40 -> tx42, reciprocal[tx41] -> tx43, 
   2.*tx43 -> tzz61, tx42*tzz61 -> tx44, pxm + tx44 + xm -> plabmx, 
   bm*cm*pxm -> tx45, dm*pxm -> tzz66, am*tzz66 -> tx46, -pzm -> tzz62, 
   am*bm -> tzz69, tzz62*tzz69 -> tx47, cm*tzz65 -> tx48, tx32*tzz63 -> tx49, 
   tx34*tzz63 -> tx50, tx45 + tx46 + tx47 + tx48 + tx49 + tx50 -> tx51, 
   tx51*tzz61 -> tx52, pym + tx52 + ym -> plabmy, tzz64*tzz68 -> tx53, 
   bm*tzz66 -> tx54, pym*tzz69 -> tx55, dm*tzz67 -> tx56, tx32*tzz62 -> tx57, 
   tx33*tzz62 -> tx58, tx53 + tx54 + tx55 + tx56 + tx57 + tx58 -> tx59, 
   tx59*tzz61 -> tx60, pzm + tx60 + zm -> plabmz}, 
 HessianStructure -> HessianStructure, CCode -> "NOT USED", 
 MathCode -> "Block[{plabmx,plabmy,plabmz,tx31,tx32,tx33,tx34,tx35,tx36,tx37,\
tx38,tx39,tx40,tx41,tx42,tx43,tx44,tx45,tx46,tx47,tx48,tx49,tx50,tx51,tx52,tx\
53,tx54,tx55,tx56,tx57,tx58,tx59,tx60,tzz61,tzz62,tzz63,tzz64,tzz65,tzz66,tzz\
67,tzz68,tzz69,xxxDummy},\n\t(*NONBOND_POSITION_RB_SET_PARAMETER(I1);*)\n\t(*\
NONBOND_POSITION_RB_SET_POSITION(am,I1,0);*)\n\t(*NONBOND_POSITION_RB_SET_POS\
ITION(bm,I1,1);*)\n\t(*NONBOND_POSITION_RB_SET_POSITION(cm,I1,2);*)\n\t(*NONB\
OND_POSITION_RB_SET_POSITION(dm,I1,3);*)\n\t(*NONBOND_POSITION_RB_SET_POSITIO\
N(xm,I1,4);*)\n\t(*NONBOND_POSITION_RB_SET_POSITION(ym,I1,5);*)\n\t(*NONBOND_\
POSITION_RB_SET_POSITION(zm,I1,6);*)\n\t(*NONBOND_POSITION_RB_SET_POINT(pxm,e\
a1,getX());*)\n\t(*NONBOND_POSITION_RB_SET_POINT(pym,ea1,getY());*)\n\t(*NONB\
OND_POSITION_RB_SET_POINT(pzm,ea1,getZ());*)\n\ttx31= mathPower2[am]; (*Rule \
12*)\n\ttx32= mathPower2[bm]; (*Rule 13*)\n\ttx33= mathPower2[cm]; (*Rule \
14*)\n\ttx34= mathPower2[dm]; (*Rule 15*)\n\ttzz67=cm pym; (*Rule \
16*)\n\ttx35=bm tzz67; (*Rule 17*)\n\ttzz63=-pym; (*Rule 18*)\n\ttx36=am dm \
tzz63; (*Rule 19*)\n\ttzz68=am cm; (*Rule 20*)\n\ttx37=pzm tzz68; (*Rule \
21*)\n\ttzz65=dm pzm; (*Rule 22*)\n\ttx38=bm tzz65; (*Rule \
23*)\n\ttzz64=-pxm; (*Rule 24*)\n\ttx39=tx33 tzz64; (*Rule 25*)\n\ttx40=tx34 \
tzz64; (*Rule 26*)\n\ttx41=tx31 + tx32 + tx33 + tx34; (*Rule \
27*)\n\ttx42=tx35 + tx36 + tx37 + tx38 + tx39 + tx40; (*Rule 28*)\n\ttx43= \
mathReciprocal[tx41]; (*Rule 29*)\n\ttzz61=2. tx43; (*Rule 30*)\n\ttx44=tx42 \
tzz61; (*Rule 31*)\n\tplabmx=pxm + tx44 + xm; (*Rule 32*)\n\ttx45=bm cm pxm; \
(*Rule 33*)\n\ttzz66=dm pxm; (*Rule 34*)\n\ttx46=am tzz66; (*Rule \
35*)\n\ttzz62=-pzm; (*Rule 36*)\n\ttzz69=am bm; (*Rule 37*)\n\ttx47=tzz62 \
tzz69; (*Rule 38*)\n\ttx48=cm tzz65; (*Rule 39*)\n\ttx49=tx32 tzz63; (*Rule \
40*)\n\ttx50=tx34 tzz63; (*Rule 41*)\n\ttx51=tx45 + tx46 + tx47 + tx48 + tx49 \
+ tx50; (*Rule 42*)\n\ttx52=tx51 tzz61; (*Rule 43*)\n\tplabmy=pym + tx52 + \
ym; (*Rule 44*)\n\ttx53=tzz64 tzz68; (*Rule 45*)\n\ttx54=bm tzz66; (*Rule \
46*)\n\ttx55=pym tzz69; (*Rule 47*)\n\ttx56=dm tzz67; (*Rule \
48*)\n\ttx57=tx32 tzz62; (*Rule 49*)\n\ttx58=tx33 tzz62; (*Rule \
50*)\n\ttx59=tx53 + tx54 + tx55 + tx56 + tx57 + tx58; (*Rule \
51*)\n\ttx60=tx59 tzz61; (*Rule 52*)\n\tplabmz=pzm + tx60 + zm; (*Rule \
53*)\n{plabmx, plabmy, plabmz}]\n"}
