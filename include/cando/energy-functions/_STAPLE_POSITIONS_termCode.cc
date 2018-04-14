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
	 tx61 = power2(ak); 		/* rule 25 */
	 tx62 = power2(bk); 		/* rule 26 */
	 tx63 = power2(ck); 		/* rule 27 */
	 tx64 = power2(dk); 		/* rule 28 */
	 tzz134 = ck*pyk; 		/* rule 29 */
	 tx65 = bk*tzz134; 		/* rule 30 */
	 tzz126 = -pyk; 		/* rule 31 */
	 tx66 = ak*dk*tzz126; 		/* rule 32 */
	 tzz137 = ak*ck; 		/* rule 33 */
	 tx67 = pzk*tzz137; 		/* rule 34 */
	 tzz131 = dk*pzk; 		/* rule 35 */
	 tx68 = bk*tzz131; 		/* rule 36 */
	 tzz128 = -pxk; 		/* rule 37 */
	 tx69 = tx63*tzz128; 		/* rule 38 */
	 tx70 = tx64*tzz128; 		/* rule 39 */
	 tx71 = tx61 + tx62 + tx63 + tx64; 		/* rule 40 */
	 tx72 = tx65 + tx66 + tx67 + tx68 + tx69 + tx70; 		/* rule 41 */
	 tx73 = reciprocal(tx71); 		/* rule 42 */
	 tzz121 = 2.*tx73; 		/* rule 43 */
	 tx74 = tx72*tzz121; 		/* rule 44 */
	 plabkx = pxk + tx74 + xk; 		/* rule 45 */
	 tx75 = bk*ck*pxk; 		/* rule 46 */
	 tzz132 = dk*pxk; 		/* rule 47 */
	 tx76 = ak*tzz132; 		/* rule 48 */
	 tzz124 = -pzk; 		/* rule 49 */
	 tzz138 = ak*bk; 		/* rule 50 */
	 tx77 = tzz124*tzz138; 		/* rule 51 */
	 tx78 = ck*tzz131; 		/* rule 52 */
	 tx79 = tx62*tzz126; 		/* rule 53 */
	 tx80 = tx64*tzz126; 		/* rule 54 */
	 tx81 = tx75 + tx76 + tx77 + tx78 + tx79 + tx80; 		/* rule 55 */
	 tx82 = tx81*tzz121; 		/* rule 56 */
	 plabky = pyk + tx82 + yk; 		/* rule 57 */
	 tx83 = tzz128*tzz137; 		/* rule 58 */
	 tx84 = bk*tzz132; 		/* rule 59 */
	 tx85 = pyk*tzz138; 		/* rule 60 */
	 tx86 = dk*tzz134; 		/* rule 61 */
	 tx87 = tx62*tzz124; 		/* rule 62 */
	 tx88 = tx63*tzz124; 		/* rule 63 */
	 tx89 = tx83 + tx84 + tx85 + tx86 + tx87 + tx88; 		/* rule 64 */
	 tx90 = tx89*tzz121; 		/* rule 65 */
	 plabkz = pzk + tx90 + zk; 		/* rule 66 */
	 tx91 = power2(al); 		/* rule 67 */
	 tx92 = power2(bl); 		/* rule 68 */
	 tx93 = power2(cl); 		/* rule 69 */
	 tx94 = power2(dl); 		/* rule 70 */
	 tzz133 = cl*pyl; 		/* rule 71 */
	 tx95 = bl*tzz133; 		/* rule 72 */
	 tzz125 = -pyl; 		/* rule 73 */
	 tx96 = al*dl*tzz125; 		/* rule 74 */
	 tzz135 = al*cl; 		/* rule 75 */
	 tx97 = pzl*tzz135; 		/* rule 76 */
	 tzz129 = dl*pzl; 		/* rule 77 */
	 tx98 = bl*tzz129; 		/* rule 78 */
	 tzz127 = -pxl; 		/* rule 79 */
	 tx99 = tx93*tzz127; 		/* rule 80 */
	 tx100 = tx94*tzz127; 		/* rule 81 */
	 tx101 = tx91 + tx92 + tx93 + tx94; 		/* rule 82 */
	 tx102 = reciprocal(tx101); 		/* rule 83 */
	 tx103 = tx100 + tx95 + tx96 + tx97 + tx98 + tx99; 		/* rule 84 */
	 tzz122 = 2.*tx102; 		/* rule 85 */
	 tx104 = tx103*tzz122; 		/* rule 86 */
	 plablx = pxl + tx104 + xl; 		/* rule 87 */
	 tx105 = bl*cl*pxl; 		/* rule 88 */
	 tzz130 = dl*pxl; 		/* rule 89 */
	 tx106 = al*tzz130; 		/* rule 90 */
	 tzz123 = -pzl; 		/* rule 91 */
	 tzz136 = al*bl; 		/* rule 92 */
	 tx107 = tzz123*tzz136; 		/* rule 93 */
	 tx108 = cl*tzz129; 		/* rule 94 */
	 tx109 = tx92*tzz125; 		/* rule 95 */
	 tx110 = tx94*tzz125; 		/* rule 96 */
	 tx111 = tx105 + tx106 + tx107 + tx108 + tx109 + tx110; 		/* rule 97 */
	 tx112 = tx111*tzz122; 		/* rule 98 */
	 plably = pyl + tx112 + yl; 		/* rule 99 */
	 tx113 = tzz127*tzz135; 		/* rule 100 */
	 tx114 = bl*tzz130; 		/* rule 101 */
	 tx115 = pyl*tzz136; 		/* rule 102 */
	 tx116 = dl*tzz133; 		/* rule 103 */
	 tx117 = tx92*tzz123; 		/* rule 104 */
	 tx118 = tx93*tzz123; 		/* rule 105 */
	 tx119 = tx113 + tx114 + tx115 + tx116 + tx117 + tx118; 		/* rule 106 */
	 tx120 = tx119*tzz122; 		/* rule 107 */
	 plablz = pzl + tx120 + zl; 		/* rule 108 */
