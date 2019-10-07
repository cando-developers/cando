// C-code
	NONBOND_POSITION_RB_SET_PARAMETER(I1);
	NONBOND_POSITION_RB_SET_POSITION(am,I1,0);
	NONBOND_POSITION_RB_SET_POSITION(bm,I1,1);
	NONBOND_POSITION_RB_SET_POSITION(cm,I1,2);
	NONBOND_POSITION_RB_SET_POSITION(dm,I1,3);
	NONBOND_POSITION_RB_SET_POSITION(xm,I1,4);
	NONBOND_POSITION_RB_SET_POSITION(ym,I1,5);
	NONBOND_POSITION_RB_SET_POSITION(zm,I1,6);
	NONBOND_POSITION_RB_SET_POINT(pxm,ea1,getX());
	NONBOND_POSITION_RB_SET_POINT(pym,ea1,getY());
	NONBOND_POSITION_RB_SET_POINT(pzm,ea1,getZ());
	 tx31 = power2(am); 		/* rule 12 */
	 tx32 = power2(bm); 		/* rule 13 */
	 tx33 = power2(cm); 		/* rule 14 */
	 tx34 = power2(dm); 		/* rule 15 */
	 tzz67 = cm*pym; 		/* rule 16 */
	 tx35 = bm*tzz67; 		/* rule 17 */
	 tzz63 = -pym; 		/* rule 18 */
	 tx36 = am*dm*tzz63; 		/* rule 19 */
	 tzz68 = am*cm; 		/* rule 20 */
	 tx37 = pzm*tzz68; 		/* rule 21 */
	 tzz65 = dm*pzm; 		/* rule 22 */
	 tx38 = bm*tzz65; 		/* rule 23 */
	 tzz64 = -pxm; 		/* rule 24 */
	 tx39 = tx33*tzz64; 		/* rule 25 */
	 tx40 = tx34*tzz64; 		/* rule 26 */
	 tx41 = tx31 + tx32 + tx33 + tx34; 		/* rule 27 */
	 tx42 = tx35 + tx36 + tx37 + tx38 + tx39 + tx40; 		/* rule 28 */
	 tx43 = reciprocal(tx41); 		/* rule 29 */
	 tzz61 = 2.*tx43; 		/* rule 30 */
	 tx44 = tx42*tzz61; 		/* rule 31 */
	 plabmx = pxm + tx44 + xm; 		/* rule 32 */
	 tx45 = bm*cm*pxm; 		/* rule 33 */
	 tzz66 = dm*pxm; 		/* rule 34 */
	 tx46 = am*tzz66; 		/* rule 35 */
	 tzz62 = -pzm; 		/* rule 36 */
	 tzz69 = am*bm; 		/* rule 37 */
	 tx47 = tzz62*tzz69; 		/* rule 38 */
	 tx48 = cm*tzz65; 		/* rule 39 */
	 tx49 = tx32*tzz63; 		/* rule 40 */
	 tx50 = tx34*tzz63; 		/* rule 41 */
	 tx51 = tx45 + tx46 + tx47 + tx48 + tx49 + tx50; 		/* rule 42 */
	 tx52 = tx51*tzz61; 		/* rule 43 */
	 plabmy = pym + tx52 + ym; 		/* rule 44 */
	 tx53 = tzz64*tzz68; 		/* rule 45 */
	 tx54 = bm*tzz66; 		/* rule 46 */
	 tx55 = pym*tzz69; 		/* rule 47 */
	 tx56 = dm*tzz67; 		/* rule 48 */
	 tx57 = tx32*tzz62; 		/* rule 49 */
	 tx58 = tx33*tzz62; 		/* rule 50 */
	 tx59 = tx53 + tx54 + tx55 + tx56 + tx57 + tx58; 		/* rule 51 */
	 tx60 = tx59*tzz61; 		/* rule 52 */
	 plabmz = pzm + tx60 + zm; 		/* rule 53 */
