// C-code
	CHIRAL_RESTRAINT_SET_PARAMETER(K);
	CHIRAL_RESTRAINT_SET_PARAMETER(CO);
	CHIRAL_RESTRAINT_SET_PARAMETER(I1);
	CHIRAL_RESTRAINT_SET_PARAMETER(I2);
	CHIRAL_RESTRAINT_SET_PARAMETER(I3);
	CHIRAL_RESTRAINT_SET_PARAMETER(I4);
	CHIRAL_RESTRAINT_SET_POSITION(x1,I1,0);
	CHIRAL_RESTRAINT_SET_POSITION(y1,I1,1);
	CHIRAL_RESTRAINT_SET_POSITION(z1,I1,2);
	CHIRAL_RESTRAINT_SET_POSITION(x2,I2,0);
	CHIRAL_RESTRAINT_SET_POSITION(y2,I2,1);
	CHIRAL_RESTRAINT_SET_POSITION(z2,I2,2);
	CHIRAL_RESTRAINT_SET_POSITION(x3,I3,0);
	CHIRAL_RESTRAINT_SET_POSITION(y3,I3,1);
	CHIRAL_RESTRAINT_SET_POSITION(z3,I3,2);
	CHIRAL_RESTRAINT_SET_POSITION(x4,I4,0);
	CHIRAL_RESTRAINT_SET_POSITION(y4,I4,1);
	CHIRAL_RESTRAINT_SET_POSITION(z4,I4,2);
	 tx789 = -x3; 		/* rule 19 */
	 tx790 = -y3; 		/* rule 20 */
	 tx791 = -z3; 		/* rule 21 */
	 tx792 = tx790*x1; 		/* rule 22 */
	 tx793 = tx789 + x1; 		/* rule 23 */
	 tx794 = tx791*x2; 		/* rule 24 */
	 tx795 = tx789 + x2; 		/* rule 25 */
	 tx796 = tx789 + x4; 		/* rule 26 */
	 tx797 = tx791*y1; 		/* rule 27 */
	 tx798 = -(x2*y1); 		/* rule 28 */
	 tx799 = x3*y1; 		/* rule 29 */
	 tx800 = tx790 + y1; 		/* rule 30 */
	 tx801 = tx789*y2; 		/* rule 31 */
	 tx802 = x1*y2; 		/* rule 32 */
	 tx803 = tx790 + y2; 		/* rule 33 */
	 tx804 = x2*y3; 		/* rule 34 */
	 tx805 = tx790 + y4; 		/* rule 35 */
	 tx806 = tx789*z1; 		/* rule 36 */
	 tx807 = x2*z1; 		/* rule 37 */
	 tx808 = -(y2*z1); 		/* rule 38 */
	 tx809 = y3*z1; 		/* rule 39 */
	 tx810 = tx791 + z1; 		/* rule 40 */
	 tx811 = tx790*z2; 		/* rule 41 */
	 tx812 = -(x1*z2); 		/* rule 42 */
	 tx813 = x3*z2; 		/* rule 43 */
	 tx814 = y1*z2; 		/* rule 44 */
	 tx815 = tx791 + z2; 		/* rule 45 */
	 tx816 = x1*z3; 		/* rule 46 */
	 tx817 = y2*z3; 		/* rule 47 */
	 tx818 = tx791 + z4; 		/* rule 48 */
	 tx819 = power2(tx793); 		/* rule 49 */
	 tx820 = power2(tx795); 		/* rule 50 */
	 tx821 = power2(tx796); 		/* rule 51 */
	 tx822 = power2(tx800); 		/* rule 52 */
	 tx823 = power2(tx803); 		/* rule 53 */
	 tx824 = tx792 + tx798 + tx799 + tx801 + tx802 + tx804; 		/* rule 54 */
	 tx825 = power2(tx805); 		/* rule 55 */
	 tx826 = power2(tx810); 		/* rule 56 */
	 tx827 = power2(tx815); 		/* rule 57 */
	 tx828 = tx794 + tx806 + tx807 + tx812 + tx813 + tx816; 		/* rule 58 */
	 tx829 = tx797 + tx808 + tx809 + tx811 + tx814 + tx817; 		/* rule 59 */
	 tx830 = power2(tx818); 		/* rule 60 */
	 tx831 = tx818*tx824; 		/* rule 61 */
	 tx832 = tx819 + tx822 + tx826; 		/* rule 62 */
	 tx833 = tx820 + tx823 + tx827; 		/* rule 63 */
	 tx834 = tx805*tx828; 		/* rule 64 */
	 tx835 = tx796*tx829; 		/* rule 65 */
	 tx836 = tx821 + tx825 + tx830; 		/* rule 66 */
	 tx1577 = mysqrt(tx832); 		/* rule 67 */
	 tx837 = reciprocal(tx1577); 		/* rule 68 */
	 tx1578 = mysqrt(tx833); 		/* rule 69 */
	 tx838 = reciprocal(tx1578); 		/* rule 70 */
	 tx839 = tx831 + tx834 + tx835; 		/* rule 71 */
	 tx1579 = mysqrt(tx836); 		/* rule 72 */
	 tx840 = reciprocal(tx1579); 		/* rule 73 */
	 tx841 = tx837*tx838*tx839*tx840; 		/* rule 74 */
	 tx842 = CO + tx841; 		/* rule 75 */
	 tx1580 = power2(tx842); 		/* rule 76 */
	 tx843 = tx1580*tx842; 		/* rule 77 */
	 ChiralTest = K*tx843; 		/* rule 78 */
/*If[ChiralTest>0.0,*/
	if (ChiralTest>0.0 ) {
	Energy=ChiralTest;
/*Energy = ChiralTest;*/
	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE(Energy);
	#ifdef CHIRAL_RESTRAINT_CALC_FORCE //[
	if ( calcForce ) {
	 tx844 = -z2; 		/* rule 86 */
	 tx845 = tx844 + z3; 		/* rule 87 */
	 tx846 = tx803*tx818; 		/* rule 88 */
	 tx847 = tx805*tx845; 		/* rule 89 */
	 tx1581 = tx832; 		/* rule 90 */
	 tx1582 = tx837; 		/* rule 91 */
	 tx1583 = reciprocal(tx1581); 		/* rule 92 */
	 tx848 = tx1582*tx1583; 		/* rule 93 */
	 tx849 = tx846 + tx847; 		/* rule 94 */
	 tx850 = -(tx793*tx838*tx839*tx840*tx848); 		/* rule 95 */
	 tx851 = tx1582*tx838*tx840*tx849; 		/* rule 96 */
	 tx852 = tx1580; 		/* rule 97 */
	 tx853 = tx850 + tx851; 		/* rule 98 */
	 gx1 = 3*K*tx852*tx853; 		/* rule 99 */
	 fx1 = -gx1; 		/* rule 100 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx854 = -x2; 		/* rule 102 */
	 tx855 = tx854 + x3; 		/* rule 103 */
	 tx856 = tx796*tx815; 		/* rule 104 */
	 tx857 = tx818*tx855; 		/* rule 105 */
	 tx858 = tx856 + tx857; 		/* rule 106 */
	 tx859 = -(tx800*tx838*tx839*tx840*tx848); 		/* rule 107 */
	 tx860 = tx1582*tx838*tx840*tx858; 		/* rule 108 */
	 tx861 = tx859 + tx860; 		/* rule 109 */
	 gy1 = 3*K*tx852*tx861; 		/* rule 110 */
	 fy1 = -gy1; 		/* rule 111 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx862 = -y2; 		/* rule 113 */
	 tx863 = tx862 + y3; 		/* rule 114 */
	 tx864 = tx795*tx805; 		/* rule 115 */
	 tx865 = tx796*tx863; 		/* rule 116 */
	 tx866 = tx864 + tx865; 		/* rule 117 */
	 tx867 = -(tx810*tx838*tx839*tx840*tx848); 		/* rule 118 */
	 tx868 = tx1582*tx838*tx840*tx866; 		/* rule 119 */
	 tx869 = tx867 + tx868; 		/* rule 120 */
	 gz1 = 3*K*tx852*tx869; 		/* rule 121 */
	 fz1 = -gz1; 		/* rule 122 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tx870 = -y1; 		/* rule 124 */
	 tx871 = tx870 + y3; 		/* rule 125 */
	 tx872 = tx805*tx810; 		/* rule 126 */
	 tx873 = tx818*tx871; 		/* rule 127 */
	 tx1584 = tx833; 		/* rule 128 */
	 tx1585 = tx838; 		/* rule 129 */
	 tx1586 = reciprocal(tx1584); 		/* rule 130 */
	 tx874 = tx1585*tx1586; 		/* rule 131 */
	 tx875 = tx872 + tx873; 		/* rule 132 */
	 tx876 = -(tx1582*tx795*tx839*tx840*tx874); 		/* rule 133 */
	 tx877 = tx1582*tx1585*tx840*tx875; 		/* rule 134 */
	 tx878 = tx876 + tx877; 		/* rule 135 */
	 gx2 = 3*K*tx852*tx878; 		/* rule 136 */
	 fx2 = -gx2; 		/* rule 137 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I2, 0, fx2 );
	 tx879 = -z1; 		/* rule 139 */
	 tx880 = tx879 + z3; 		/* rule 140 */
	 tx881 = tx793*tx818; 		/* rule 141 */
	 tx882 = tx796*tx880; 		/* rule 142 */
	 tx883 = tx881 + tx882; 		/* rule 143 */
	 tx884 = -(tx1582*tx803*tx839*tx840*tx874); 		/* rule 144 */
	 tx885 = tx1582*tx1585*tx840*tx883; 		/* rule 145 */
	 tx886 = tx884 + tx885; 		/* rule 146 */
	 gy2 = 3*K*tx852*tx886; 		/* rule 147 */
	 fy2 = -gy2; 		/* rule 148 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I2, 1, fy2 );
	 tx887 = -x1; 		/* rule 150 */
	 tx888 = tx887 + x3; 		/* rule 151 */
	 tx889 = tx796*tx800; 		/* rule 152 */
	 tx890 = tx805*tx888; 		/* rule 153 */
	 tx891 = tx889 + tx890; 		/* rule 154 */
	 tx892 = -(tx1582*tx815*tx839*tx840*tx874); 		/* rule 155 */
	 tx893 = tx1582*tx1585*tx840*tx891; 		/* rule 156 */
	 tx894 = tx892 + tx893; 		/* rule 157 */
	 gz2 = 3*K*tx852*tx894; 		/* rule 158 */
	 fz2 = -gz2; 		/* rule 159 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I2, 2, fz2 );
	 tx895 = tx862 + y1; 		/* rule 161 */
	 tx896 = tx879 + z2; 		/* rule 162 */
	 tx897 = -tx814; 		/* rule 163 */
	 tx898 = tx818*tx895; 		/* rule 164 */
	 tx899 = tx805*tx896; 		/* rule 165 */
	 tx900 = tx791*y2; 		/* rule 166 */
	 tx901 = tx790*z1; 		/* rule 167 */
	 tx902 = y2*z1; 		/* rule 168 */
	 tx903 = y3*z2; 		/* rule 169 */
	 tx904 = y1*z3; 		/* rule 170 */
	 tx1587 = tx836; 		/* rule 171 */
	 tx1588 = tx840; 		/* rule 172 */
	 tx1589 = reciprocal(tx1587); 		/* rule 173 */
	 tx905 = tx1588*tx1589; 		/* rule 174 */
	 tx906 = tx897 + tx898 + tx899 + tx900 + tx901 + tx902 + tx903 + tx904; 		/* rule 175 */
	 tx907 = tx1585*tx1588*tx793*tx839*tx848; 		/* rule 176 */
	 tx908 = tx1582*tx1588*tx795*tx839*tx874; 		/* rule 177 */
	 tx909 = tx1582*tx1585*tx796*tx839*tx905; 		/* rule 178 */
	 tx910 = tx1582*tx1585*tx1588*tx906; 		/* rule 179 */
	 tx911 = tx907 + tx908 + tx909 + tx910; 		/* rule 180 */
	 gx3 = 3*K*tx852*tx911; 		/* rule 181 */
	 fx3 = -gx3; 		/* rule 182 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I3, 0, fx3 );
	 tx912 = tx887 + x2; 		/* rule 184 */
	 tx913 = tx844 + z1; 		/* rule 185 */
	 tx914 = -tx807; 		/* rule 186 */
	 tx915 = tx818*tx912; 		/* rule 187 */
	 tx916 = tx796*tx913; 		/* rule 188 */
	 tx917 = tx791*x1; 		/* rule 189 */
	 tx918 = x3*z1; 		/* rule 190 */
	 tx919 = tx789*z2; 		/* rule 191 */
	 tx920 = x1*z2; 		/* rule 192 */
	 tx921 = x2*z3; 		/* rule 193 */
	 tx922 = tx914 + tx915 + tx916 + tx917 + tx918 + tx919 + tx920 + tx921; 		/* rule 194 */
	 tx923 = tx1585*tx1588*tx800*tx839*tx848; 		/* rule 195 */
	 tx924 = tx1582*tx1588*tx803*tx839*tx874; 		/* rule 196 */
	 tx925 = tx1582*tx1585*tx805*tx839*tx905; 		/* rule 197 */
	 tx926 = tx1582*tx1585*tx1588*tx922; 		/* rule 198 */
	 tx927 = tx923 + tx924 + tx925 + tx926; 		/* rule 199 */
	 gy3 = 3*K*tx852*tx927; 		/* rule 200 */
	 fy3 = -gy3; 		/* rule 201 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I3, 1, fy3 );
	 tx928 = tx854 + x1; 		/* rule 203 */
	 tx929 = tx870 + y2; 		/* rule 204 */
	 tx930 = -tx802; 		/* rule 205 */
	 tx931 = tx805*tx928; 		/* rule 206 */
	 tx932 = tx796*tx929; 		/* rule 207 */
	 tx933 = tx790*x2; 		/* rule 208 */
	 tx934 = tx789*y1; 		/* rule 209 */
	 tx935 = x2*y1; 		/* rule 210 */
	 tx936 = x3*y2; 		/* rule 211 */
	 tx937 = x1*y3; 		/* rule 212 */
	 tx938 = tx930 + tx931 + tx932 + tx933 + tx934 + tx935 + tx936 + tx937; 		/* rule 213 */
	 tx939 = tx1585*tx1588*tx810*tx839*tx848; 		/* rule 214 */
	 tx940 = tx1582*tx1588*tx815*tx839*tx874; 		/* rule 215 */
	 tx941 = tx1582*tx1585*tx818*tx839*tx905; 		/* rule 216 */
	 tx942 = tx1582*tx1585*tx1588*tx938; 		/* rule 217 */
	 tx943 = tx939 + tx940 + tx941 + tx942; 		/* rule 218 */
	 gz3 = 3*K*tx852*tx943; 		/* rule 219 */
	 fz3 = -gz3; 		/* rule 220 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I3, 2, fz3 );
	 tx944 = tx1582*tx1585*tx1588*tx829; 		/* rule 222 */
	 tx945 = -tx909; 		/* rule 223 */
	 tx946 = tx944 + tx945; 		/* rule 224 */
	 gx4 = 3*K*tx852*tx946; 		/* rule 225 */
	 fx4 = -gx4; 		/* rule 226 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I4, 0, fx4 );
	 tx947 = tx1582*tx1585*tx1588*tx828; 		/* rule 228 */
	 tx948 = -tx925; 		/* rule 229 */
	 tx949 = tx947 + tx948; 		/* rule 230 */
	 gy4 = 3*K*tx852*tx949; 		/* rule 231 */
	 fy4 = -gy4; 		/* rule 232 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I4, 1, fy4 );
	 tx950 = tx1582*tx1585*tx1588*tx824; 		/* rule 234 */
	 tx951 = -tx941; 		/* rule 235 */
	 tx952 = tx950 + tx951; 		/* rule 236 */
	 gz4 = 3*K*tx852*tx952; 		/* rule 237 */
	 fz4 = -gz4; 		/* rule 238 */
	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(I4, 2, fz4 );
	#ifdef CHIRAL_RESTRAINT_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx1590 = power2(tx1581); 		/* rule 242 */
	 tx1591 = reciprocal(tx1590); 		/* rule 243 */
	 tx953 = tx1582*tx1591; 		/* rule 244 */
	 tx954 = -(tx1585*tx1588*tx839*tx848); 		/* rule 245 */
	 tx955 = -2.*tx1585*tx1588*tx793*tx848*tx849; 		/* rule 246 */
	 tx956 = 3*tx1585*tx1588*tx819*tx839*tx953; 		/* rule 247 */
	 tx957 = power2(tx853); 		/* rule 248 */
	 tx958 = tx954 + tx955 + tx956; 		/* rule 249 */
	 tx959 = 6.*K*tx842*tx957; 		/* rule 250 */
	 tx960 = 3*K*tx852*tx958; 		/* rule 251 */
	 dhx1x1 = tx959 + tx960; 		/* rule 252 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx961 = -2.*tx1585*tx1588*tx800*tx848*tx858; 		/* rule 254 */
	 tx962 = 3*tx1585*tx1588*tx822*tx839*tx953; 		/* rule 255 */
	 tx963 = power2(tx861); 		/* rule 256 */
	 tx964 = tx954 + tx961 + tx962; 		/* rule 257 */
	 tx965 = 6.*K*tx842*tx963; 		/* rule 258 */
	 tx966 = 3*K*tx852*tx964; 		/* rule 259 */
	 dhy1y1 = tx965 + tx966; 		/* rule 260 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx967 = -2.*tx1585*tx1588*tx810*tx848*tx866; 		/* rule 262 */
	 tx968 = 3*tx1585*tx1588*tx826*tx839*tx953; 		/* rule 263 */
	 tx969 = power2(tx869); 		/* rule 264 */
	 tx970 = tx954 + tx967 + tx968; 		/* rule 265 */
	 tx971 = 6.*K*tx842*tx969; 		/* rule 266 */
	 tx972 = 3*K*tx852*tx970; 		/* rule 267 */
	 dhz1z1 = tx971 + tx972; 		/* rule 268 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 tx1592 = power2(tx1584); 		/* rule 270 */
	 tx1593 = reciprocal(tx1592); 		/* rule 271 */
	 tx973 = tx1585*tx1593; 		/* rule 272 */
	 tx974 = -(tx1582*tx1588*tx839*tx874); 		/* rule 273 */
	 tx975 = -2.*tx1582*tx1588*tx795*tx874*tx875; 		/* rule 274 */
	 tx976 = 3*tx1582*tx1588*tx820*tx839*tx973; 		/* rule 275 */
	 tx977 = power2(tx878); 		/* rule 276 */
	 tx978 = tx974 + tx975 + tx976; 		/* rule 277 */
	 tx979 = 6.*K*tx842*tx977; 		/* rule 278 */
	 tx980 = 3*K*tx852*tx978; 		/* rule 279 */
	 dhx2x2 = tx979 + tx980; 		/* rule 280 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 tx981 = -2.*tx1582*tx1588*tx803*tx874*tx883; 		/* rule 282 */
	 tx982 = 3*tx1582*tx1588*tx823*tx839*tx973; 		/* rule 283 */
	 tx983 = power2(tx886); 		/* rule 284 */
	 tx984 = tx974 + tx981 + tx982; 		/* rule 285 */
	 tx985 = 6.*K*tx842*tx983; 		/* rule 286 */
	 tx986 = 3*K*tx852*tx984; 		/* rule 287 */
	 dhy2y2 = tx985 + tx986; 		/* rule 288 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 tx987 = -2.*tx1582*tx1588*tx815*tx874*tx891; 		/* rule 290 */
	 tx988 = 3*tx1582*tx1588*tx827*tx839*tx973; 		/* rule 291 */
	 tx989 = power2(tx894); 		/* rule 292 */
	 tx990 = tx974 + tx987 + tx988; 		/* rule 293 */
	 tx991 = 6.*K*tx842*tx989; 		/* rule 294 */
	 tx992 = 3*K*tx852*tx990; 		/* rule 295 */
	 dhz2z2 = tx991 + tx992; 		/* rule 296 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	 tx1594 = power2(tx1587); 		/* rule 298 */
	 tx1595 = reciprocal(tx1594); 		/* rule 299 */
	 tx993 = tx1588*tx1595; 		/* rule 300 */
	 tx994 = 2.*tx1588*tx793*tx795*tx839*tx848*tx874; 		/* rule 301 */
	 tx995 = -(tx1582*tx1585*tx839*tx905); 		/* rule 302 */
	 tx996 = 2.*tx1585*tx793*tx796*tx839*tx848*tx905; 		/* rule 303 */
	 tx997 = 2.*tx1582*tx795*tx796*tx839*tx874*tx905; 		/* rule 304 */
	 tx998 = 2.*tx1585*tx1588*tx793*tx848*tx906; 		/* rule 305 */
	 tx999 = 2.*tx1582*tx1588*tx795*tx874*tx906; 		/* rule 306 */
	 tx1000 = 2.*tx1582*tx1585*tx796*tx905*tx906; 		/* rule 307 */
	 tx1001 = 3*tx1582*tx1585*tx821*tx839*tx993; 		/* rule 308 */
	 tx1002 = power2(tx911); 		/* rule 309 */
	 tx1003 = tx1000 + tx1001 + tx954 + tx956 + tx974 + tx976 + tx994 + tx995 + tx996 + tx997 + tx998 + tx999; 		/* rule 310 */
	 tx1004 = 6.*K*tx1002*tx842; 		/* rule 311 */
	 tx1005 = 3*K*tx1003*tx852; 		/* rule 312 */
	 dhx3x3 = tx1004 + tx1005; 		/* rule 313 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 0, dhx3x3);
	 tx1006 = 2.*tx1588*tx800*tx803*tx839*tx848*tx874; 		/* rule 315 */
	 tx1007 = 2.*tx1585*tx800*tx805*tx839*tx848*tx905; 		/* rule 316 */
	 tx1008 = 2.*tx1582*tx803*tx805*tx839*tx874*tx905; 		/* rule 317 */
	 tx1009 = 2.*tx1585*tx1588*tx800*tx848*tx922; 		/* rule 318 */
	 tx1010 = 2.*tx1582*tx1588*tx803*tx874*tx922; 		/* rule 319 */
	 tx1011 = 2.*tx1582*tx1585*tx805*tx905*tx922; 		/* rule 320 */
	 tx1012 = 3*tx1582*tx1585*tx825*tx839*tx993; 		/* rule 321 */
	 tx1013 = power2(tx927); 		/* rule 322 */
	 tx1014 = tx1006 + tx1007 + tx1008 + tx1009 + tx1010 + tx1011 + tx1012 + tx954 + tx962 + tx974 + tx982 + tx995; 		/* rule 323 */
	 tx1015 = 6.*K*tx1013*tx842; 		/* rule 324 */
	 tx1016 = 3*K*tx1014*tx852; 		/* rule 325 */
	 dhy3y3 = tx1015 + tx1016; 		/* rule 326 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 1, dhy3y3);
	 tx1017 = 2.*tx1588*tx810*tx815*tx839*tx848*tx874; 		/* rule 328 */
	 tx1018 = 2.*tx1585*tx810*tx818*tx839*tx848*tx905; 		/* rule 329 */
	 tx1019 = 2.*tx1582*tx815*tx818*tx839*tx874*tx905; 		/* rule 330 */
	 tx1020 = 2.*tx1585*tx1588*tx810*tx848*tx938; 		/* rule 331 */
	 tx1021 = 2.*tx1582*tx1588*tx815*tx874*tx938; 		/* rule 332 */
	 tx1022 = 2.*tx1582*tx1585*tx818*tx905*tx938; 		/* rule 333 */
	 tx1023 = 3*tx1582*tx1585*tx830*tx839*tx993; 		/* rule 334 */
	 tx1024 = power2(tx943); 		/* rule 335 */
	 tx1025 = tx1017 + tx1018 + tx1019 + tx1020 + tx1021 + tx1022 + tx1023 + tx954 + tx968 + tx974 + tx988 + tx995; 		/* rule 336 */
	 tx1026 = 6.*K*tx1024*tx842; 		/* rule 337 */
	 tx1027 = 3*K*tx1025*tx852; 		/* rule 338 */
	 dhz3z3 = tx1026 + tx1027; 		/* rule 339 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I3, 2, dhz3z3);
	 tx1028 = -2.*tx1582*tx1585*tx835*tx905; 		/* rule 341 */
	 tx1029 = power2(tx946); 		/* rule 342 */
	 tx1030 = tx1001 + tx1028 + tx995; 		/* rule 343 */
	 tx1031 = 6.*K*tx1029*tx842; 		/* rule 344 */
	 tx1032 = 3*K*tx1030*tx852; 		/* rule 345 */
	 dhx4x4 = tx1031 + tx1032; 		/* rule 346 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 0, dhx4x4);
	 tx1033 = -2.*tx1582*tx1585*tx834*tx905; 		/* rule 348 */
	 tx1034 = power2(tx949); 		/* rule 349 */
	 tx1035 = tx1012 + tx1033 + tx995; 		/* rule 350 */
	 tx1036 = 6.*K*tx1034*tx842; 		/* rule 351 */
	 tx1037 = 3*K*tx1035*tx852; 		/* rule 352 */
	 dhy4y4 = tx1036 + tx1037; 		/* rule 353 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I4, 1, I4, 1, dhy4y4);
	 tx1038 = -2.*tx1582*tx1585*tx831*tx905; 		/* rule 355 */
	 tx1039 = power2(tx952); 		/* rule 356 */
	 tx1040 = tx1023 + tx1038 + tx995; 		/* rule 357 */
	 tx1041 = 6.*K*tx1039*tx842; 		/* rule 358 */
	 tx1042 = 3*K*tx1040*tx852; 		/* rule 359 */
	 dhz4z4 = tx1041 + tx1042; 		/* rule 360 */
	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I4, 2, I4, 2, dhz4z4);
	#ifdef CHIRAL_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tx1043 = -(tx1585*tx1588*tx800*tx848*tx849); 		/* rule 364 */
	 tx1044 = -(tx1585*tx1588*tx793*tx848*tx858); 		/* rule 365 */
	 tx1045 = 3*tx1585*tx1588*tx793*tx800*tx839*tx953; 		/* rule 366 */
	 tx1046 = tx1043 + tx1044 + tx1045; 		/* rule 367 */
	 tx1047 = 3*K*tx1046*tx852; 		/* rule 368 */
	 tx1048 = 6.*K*tx842*tx853*tx861; 		/* rule 369 */
	 ohx1y1 = tx1047 + tx1048; 		/* rule 370 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tx1049 = -(tx1585*tx1588*tx810*tx848*tx849); 		/* rule 372 */
	 tx1050 = -(tx1585*tx1588*tx793*tx848*tx866); 		/* rule 373 */
	 tx1051 = 3*tx1585*tx1588*tx793*tx810*tx839*tx953; 		/* rule 374 */
	 tx1052 = tx1049 + tx1050 + tx1051; 		/* rule 375 */
	 tx1053 = 3*K*tx1052*tx852; 		/* rule 376 */
	 tx1054 = 6.*K*tx842*tx853*tx869; 		/* rule 377 */
	 ohx1z1 = tx1053 + tx1054; 		/* rule 378 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx1055 = tx1588*tx793*tx795*tx839*tx848*tx874; 		/* rule 380 */
	 tx1056 = -(tx1582*tx1588*tx795*tx849*tx874); 		/* rule 381 */
	 tx1057 = -(tx1585*tx1588*tx793*tx848*tx875); 		/* rule 382 */
	 tx1058 = tx1055 + tx1056 + tx1057; 		/* rule 383 */
	 tx1059 = 3*K*tx1058*tx852; 		/* rule 384 */
	 tx1060 = 6.*K*tx842*tx853*tx878; 		/* rule 385 */
	 ohx1x2 = tx1059 + tx1060; 		/* rule 386 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx1061 = tx1582*tx1585*tx1588*tx818; 		/* rule 388 */
	 tx1062 = tx1588*tx793*tx803*tx839*tx848*tx874; 		/* rule 389 */
	 tx1063 = -(tx1582*tx1588*tx803*tx849*tx874); 		/* rule 390 */
	 tx1064 = -(tx1585*tx1588*tx793*tx848*tx883); 		/* rule 391 */
	 tx1065 = tx1061 + tx1062 + tx1063 + tx1064; 		/* rule 392 */
	 tx1066 = 3*K*tx1065*tx852; 		/* rule 393 */
	 tx1067 = 6.*K*tx842*tx853*tx886; 		/* rule 394 */
	 ohx1y2 = tx1066 + tx1067; 		/* rule 395 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx1068 = -y4; 		/* rule 397 */
	 tx1069 = tx1068 + y3; 		/* rule 398 */
	 tx1070 = tx1069*tx1582*tx1585*tx1588; 		/* rule 399 */
	 tx1071 = tx1588*tx793*tx815*tx839*tx848*tx874; 		/* rule 400 */
	 tx1072 = -(tx1582*tx1588*tx815*tx849*tx874); 		/* rule 401 */
	 tx1073 = -(tx1585*tx1588*tx793*tx848*tx891); 		/* rule 402 */
	 tx1074 = tx1070 + tx1071 + tx1072 + tx1073; 		/* rule 403 */
	 tx1075 = 3*K*tx1074*tx852; 		/* rule 404 */
	 tx1076 = 6.*K*tx842*tx853*tx894; 		/* rule 405 */
	 ohx1z2 = tx1075 + tx1076; 		/* rule 406 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tx1077 = -tx1055; 		/* rule 408 */
	 tx1078 = tx1585*tx1588*tx839*tx848; 		/* rule 409 */
	 tx1079 = tx1585*tx1588*tx793*tx848*tx849; 		/* rule 410 */
	 tx1080 = tx1582*tx1588*tx795*tx849*tx874; 		/* rule 411 */
	 tx1081 = -(tx1585*tx793*tx796*tx839*tx848*tx905); 		/* rule 412 */
	 tx1082 = tx1582*tx1585*tx796*tx849*tx905; 		/* rule 413 */
	 tx1083 = -(tx1585*tx1588*tx793*tx848*tx906); 		/* rule 414 */
	 tx1084 = -3*tx1585*tx1588*tx819*tx839*tx953; 		/* rule 415 */
	 tx1085 = tx1077 + tx1078 + tx1079 + tx1080 + tx1081 + tx1082 + tx1083 + tx1084; 		/* rule 416 */
	 tx1086 = 3*K*tx1085*tx852; 		/* rule 417 */
	 tx1087 = 6.*K*tx842*tx853*tx911; 		/* rule 418 */
	 ohx1x3 = tx1086 + tx1087; 		/* rule 419 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 0, ohx1x3);
	 tx1088 = -z4; 		/* rule 421 */
	 tx1089 = tx1088 + z2; 		/* rule 422 */
	 tx1090 = -tx1062; 		/* rule 423 */
	 tx1091 = tx1089*tx1582*tx1585*tx1588; 		/* rule 424 */
	 tx1092 = tx1585*tx1588*tx800*tx848*tx849; 		/* rule 425 */
	 tx1093 = tx1582*tx1588*tx803*tx849*tx874; 		/* rule 426 */
	 tx1094 = -(tx1585*tx793*tx805*tx839*tx848*tx905); 		/* rule 427 */
	 tx1095 = tx1582*tx1585*tx805*tx849*tx905; 		/* rule 428 */
	 tx1096 = -(tx1585*tx1588*tx793*tx848*tx922); 		/* rule 429 */
	 tx1097 = -3*tx1585*tx1588*tx793*tx800*tx839*tx953; 		/* rule 430 */
	 tx1098 = tx1090 + tx1091 + tx1092 + tx1093 + tx1094 + tx1095 + tx1096 + tx1097; 		/* rule 431 */
	 tx1099 = 3*K*tx1098*tx852; 		/* rule 432 */
	 tx1100 = 6.*K*tx842*tx853*tx927; 		/* rule 433 */
	 ohx1y3 = tx1099 + tx1100; 		/* rule 434 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 1, ohx1y3);
	 tx1101 = tx862 + y4; 		/* rule 436 */
	 tx1102 = -tx1071; 		/* rule 437 */
	 tx1103 = tx1101*tx1582*tx1585*tx1588; 		/* rule 438 */
	 tx1104 = tx1585*tx1588*tx810*tx848*tx849; 		/* rule 439 */
	 tx1105 = tx1582*tx1588*tx815*tx849*tx874; 		/* rule 440 */
	 tx1106 = tx1582*tx1585*tx818*tx849*tx905; 		/* rule 441 */
	 tx1107 = -(tx1585*tx839*tx848*tx881*tx905); 		/* rule 442 */
	 tx1108 = -(tx1585*tx1588*tx793*tx848*tx938); 		/* rule 443 */
	 tx1109 = -3*tx1585*tx1588*tx793*tx810*tx839*tx953; 		/* rule 444 */
	 tx1110 = tx1102 + tx1103 + tx1104 + tx1105 + tx1106 + tx1107 + tx1108 + tx1109; 		/* rule 445 */
	 tx1111 = 3*K*tx1110*tx852; 		/* rule 446 */
	 tx1112 = 6.*K*tx842*tx853*tx943; 		/* rule 447 */
	 ohx1z3 = tx1111 + tx1112; 		/* rule 448 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 2, ohx1z3);
	 tx1113 = -tx1082; 		/* rule 450 */
	 tx1114 = -(tx1585*tx1588*tx793*tx829*tx848); 		/* rule 451 */
	 tx1115 = tx1585*tx793*tx796*tx839*tx848*tx905; 		/* rule 452 */
	 tx1116 = tx1113 + tx1114 + tx1115; 		/* rule 453 */
	 tx1117 = 3*K*tx1116*tx852; 		/* rule 454 */
	 tx1118 = 6.*K*tx842*tx853*tx946; 		/* rule 455 */
	 ohx1x4 = tx1117 + tx1118; 		/* rule 456 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 0, ohx1x4);
	 tx1119 = -tx1095; 		/* rule 458 */
	 tx1120 = tx1582*tx1585*tx1588*tx845; 		/* rule 459 */
	 tx1121 = -(tx1585*tx1588*tx793*tx828*tx848); 		/* rule 460 */
	 tx1122 = tx1585*tx793*tx805*tx839*tx848*tx905; 		/* rule 461 */
	 tx1123 = tx1119 + tx1120 + tx1121 + tx1122; 		/* rule 462 */
	 tx1124 = 3*K*tx1123*tx852; 		/* rule 463 */
	 tx1125 = 6.*K*tx842*tx853*tx949; 		/* rule 464 */
	 ohx1y4 = tx1124 + tx1125; 		/* rule 465 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 1, ohx1y4);
	 tx1126 = -tx1106; 		/* rule 467 */
	 tx1127 = tx1582*tx1585*tx1588*tx803; 		/* rule 468 */
	 tx1128 = -(tx1585*tx1588*tx793*tx824*tx848); 		/* rule 469 */
	 tx1129 = tx1585*tx839*tx848*tx881*tx905; 		/* rule 470 */
	 tx1130 = tx1126 + tx1127 + tx1128 + tx1129; 		/* rule 471 */
	 tx1131 = 3*K*tx1130*tx852; 		/* rule 472 */
	 tx1132 = 6.*K*tx842*tx853*tx952; 		/* rule 473 */
	 ohx1z4 = tx1131 + tx1132; 		/* rule 474 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 2, ohx1z4);
	 tx1133 = -(tx1585*tx1588*tx810*tx848*tx858); 		/* rule 476 */
	 tx1134 = -(tx1585*tx1588*tx800*tx848*tx866); 		/* rule 477 */
	 tx1135 = 3*tx1585*tx1588*tx800*tx810*tx839*tx953; 		/* rule 478 */
	 tx1136 = tx1133 + tx1134 + tx1135; 		/* rule 479 */
	 tx1137 = 3*K*tx1136*tx852; 		/* rule 480 */
	 tx1138 = 6.*K*tx842*tx861*tx869; 		/* rule 481 */
	 ohy1z1 = tx1137 + tx1138; 		/* rule 482 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 tx1139 = tx1088 + z3; 		/* rule 484 */
	 tx1140 = tx1139*tx1582*tx1585*tx1588; 		/* rule 485 */
	 tx1141 = tx1588*tx795*tx800*tx839*tx848*tx874; 		/* rule 486 */
	 tx1142 = -(tx1582*tx1588*tx795*tx858*tx874); 		/* rule 487 */
	 tx1143 = -(tx1585*tx1588*tx800*tx848*tx875); 		/* rule 488 */
	 tx1144 = tx1140 + tx1141 + tx1142 + tx1143; 		/* rule 489 */
	 tx1145 = 3*K*tx1144*tx852; 		/* rule 490 */
	 tx1146 = 6.*K*tx842*tx861*tx878; 		/* rule 491 */
	 ohy1x2 = tx1145 + tx1146; 		/* rule 492 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx1147 = tx1588*tx800*tx803*tx839*tx848*tx874; 		/* rule 494 */
	 tx1148 = -(tx1582*tx1588*tx803*tx858*tx874); 		/* rule 495 */
	 tx1149 = -(tx1585*tx1588*tx800*tx848*tx883); 		/* rule 496 */
	 tx1150 = tx1147 + tx1148 + tx1149; 		/* rule 497 */
	 tx1151 = 3*K*tx1150*tx852; 		/* rule 498 */
	 tx1152 = 6.*K*tx842*tx861*tx886; 		/* rule 499 */
	 ohy1y2 = tx1151 + tx1152; 		/* rule 500 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx1153 = tx1582*tx1585*tx1588*tx796; 		/* rule 502 */
	 tx1154 = tx1588*tx800*tx815*tx839*tx848*tx874; 		/* rule 503 */
	 tx1155 = -(tx1582*tx1588*tx815*tx858*tx874); 		/* rule 504 */
	 tx1156 = -(tx1585*tx1588*tx800*tx848*tx891); 		/* rule 505 */
	 tx1157 = tx1153 + tx1154 + tx1155 + tx1156; 		/* rule 506 */
	 tx1158 = 3*K*tx1157*tx852; 		/* rule 507 */
	 tx1159 = 6.*K*tx842*tx861*tx894; 		/* rule 508 */
	 ohy1z2 = tx1158 + tx1159; 		/* rule 509 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 tx1160 = tx844 + z4; 		/* rule 511 */
	 tx1161 = -tx1141; 		/* rule 512 */
	 tx1162 = tx1160*tx1582*tx1585*tx1588; 		/* rule 513 */
	 tx1163 = tx1585*tx1588*tx793*tx848*tx858; 		/* rule 514 */
	 tx1164 = tx1582*tx1588*tx795*tx858*tx874; 		/* rule 515 */
	 tx1165 = tx1582*tx1585*tx796*tx858*tx905; 		/* rule 516 */
	 tx1166 = -(tx1585*tx839*tx848*tx889*tx905); 		/* rule 517 */
	 tx1167 = -(tx1585*tx1588*tx800*tx848*tx906); 		/* rule 518 */
	 tx1168 = tx1097 + tx1161 + tx1162 + tx1163 + tx1164 + tx1165 + tx1166 + tx1167; 		/* rule 519 */
	 tx1169 = 3*K*tx1168*tx852; 		/* rule 520 */
	 tx1170 = 6.*K*tx842*tx861*tx911; 		/* rule 521 */
	 ohy1x3 = tx1169 + tx1170; 		/* rule 522 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 0, ohy1x3);
	 tx1171 = -tx1147; 		/* rule 524 */
	 tx1172 = tx1585*tx1588*tx800*tx848*tx858; 		/* rule 525 */
	 tx1173 = tx1582*tx1588*tx803*tx858*tx874; 		/* rule 526 */
	 tx1174 = -(tx1585*tx800*tx805*tx839*tx848*tx905); 		/* rule 527 */
	 tx1175 = tx1582*tx1585*tx805*tx858*tx905; 		/* rule 528 */
	 tx1176 = -(tx1585*tx1588*tx800*tx848*tx922); 		/* rule 529 */
	 tx1177 = -3*tx1585*tx1588*tx822*tx839*tx953; 		/* rule 530 */
	 tx1178 = tx1078 + tx1171 + tx1172 + tx1173 + tx1174 + tx1175 + tx1176 + tx1177; 		/* rule 531 */
	 tx1179 = 3*K*tx1178*tx852; 		/* rule 532 */
	 tx1180 = 6.*K*tx842*tx861*tx927; 		/* rule 533 */
	 ohy1y3 = tx1179 + tx1180; 		/* rule 534 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 1, ohy1y3);
	 tx1181 = -x4; 		/* rule 536 */
	 tx1182 = tx1181 + x2; 		/* rule 537 */
	 tx1183 = -tx1154; 		/* rule 538 */
	 tx1184 = tx1182*tx1582*tx1585*tx1588; 		/* rule 539 */
	 tx1185 = tx1585*tx1588*tx810*tx848*tx858; 		/* rule 540 */
	 tx1186 = tx1582*tx1588*tx815*tx858*tx874; 		/* rule 541 */
	 tx1187 = -(tx1585*tx800*tx818*tx839*tx848*tx905); 		/* rule 542 */
	 tx1188 = tx1582*tx1585*tx818*tx858*tx905; 		/* rule 543 */
	 tx1189 = -(tx1585*tx1588*tx800*tx848*tx938); 		/* rule 544 */
	 tx1190 = -3*tx1585*tx1588*tx800*tx810*tx839*tx953; 		/* rule 545 */
	 tx1191 = tx1183 + tx1184 + tx1185 + tx1186 + tx1187 + tx1188 + tx1189 + tx1190; 		/* rule 546 */
	 tx1192 = 3*K*tx1191*tx852; 		/* rule 547 */
	 tx1193 = 6.*K*tx842*tx861*tx943; 		/* rule 548 */
	 ohy1z3 = tx1192 + tx1193; 		/* rule 549 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 2, ohy1z3);
	 tx1194 = -tx1165; 		/* rule 551 */
	 tx1195 = tx1582*tx1585*tx1588*tx815; 		/* rule 552 */
	 tx1196 = -(tx1585*tx1588*tx800*tx829*tx848); 		/* rule 553 */
	 tx1197 = tx1585*tx839*tx848*tx889*tx905; 		/* rule 554 */
	 tx1198 = tx1194 + tx1195 + tx1196 + tx1197; 		/* rule 555 */
	 tx1199 = 3*K*tx1198*tx852; 		/* rule 556 */
	 tx1200 = 6.*K*tx842*tx861*tx946; 		/* rule 557 */
	 ohy1x4 = tx1199 + tx1200; 		/* rule 558 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 0, ohy1x4);
	 tx1201 = -tx1175; 		/* rule 560 */
	 tx1202 = -(tx1585*tx1588*tx800*tx828*tx848); 		/* rule 561 */
	 tx1203 = tx1585*tx800*tx805*tx839*tx848*tx905; 		/* rule 562 */
	 tx1204 = tx1201 + tx1202 + tx1203; 		/* rule 563 */
	 tx1205 = 3*K*tx1204*tx852; 		/* rule 564 */
	 tx1206 = 6.*K*tx842*tx861*tx949; 		/* rule 565 */
	 ohy1y4 = tx1205 + tx1206; 		/* rule 566 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 1, ohy1y4);
	 tx1207 = -tx1188; 		/* rule 568 */
	 tx1208 = -(tx1585*tx1588*tx800*tx824*tx848); 		/* rule 569 */
	 tx1209 = tx1582*tx1585*tx1588*tx855; 		/* rule 570 */
	 tx1210 = tx1585*tx800*tx818*tx839*tx848*tx905; 		/* rule 571 */
	 tx1211 = tx1207 + tx1208 + tx1209 + tx1210; 		/* rule 572 */
	 tx1212 = 3*K*tx1211*tx852; 		/* rule 573 */
	 tx1213 = 6.*K*tx842*tx861*tx952; 		/* rule 574 */
	 ohy1z4 = tx1212 + tx1213; 		/* rule 575 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 2, ohy1z4);
	 tx1214 = tx1582*tx1585*tx1588*tx805; 		/* rule 577 */
	 tx1215 = tx1588*tx795*tx810*tx839*tx848*tx874; 		/* rule 578 */
	 tx1216 = -(tx1582*tx1588*tx795*tx866*tx874); 		/* rule 579 */
	 tx1217 = -(tx1585*tx1588*tx810*tx848*tx875); 		/* rule 580 */
	 tx1218 = tx1214 + tx1215 + tx1216 + tx1217; 		/* rule 581 */
	 tx1219 = 3*K*tx1218*tx852; 		/* rule 582 */
	 tx1220 = 6.*K*tx842*tx869*tx878; 		/* rule 583 */
	 ohz1x2 = tx1219 + tx1220; 		/* rule 584 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 tx1221 = tx1181 + x3; 		/* rule 586 */
	 tx1222 = tx1221*tx1582*tx1585*tx1588; 		/* rule 587 */
	 tx1223 = tx1588*tx803*tx810*tx839*tx848*tx874; 		/* rule 588 */
	 tx1224 = -(tx1582*tx1588*tx803*tx866*tx874); 		/* rule 589 */
	 tx1225 = -(tx1585*tx1588*tx810*tx848*tx883); 		/* rule 590 */
	 tx1226 = tx1222 + tx1223 + tx1224 + tx1225; 		/* rule 591 */
	 tx1227 = 3*K*tx1226*tx852; 		/* rule 592 */
	 tx1228 = 6.*K*tx842*tx869*tx886; 		/* rule 593 */
	 ohz1y2 = tx1227 + tx1228; 		/* rule 594 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx1229 = tx1588*tx810*tx815*tx839*tx848*tx874; 		/* rule 596 */
	 tx1230 = -(tx1582*tx1588*tx815*tx866*tx874); 		/* rule 597 */
	 tx1231 = -(tx1585*tx1588*tx810*tx848*tx891); 		/* rule 598 */
	 tx1232 = tx1229 + tx1230 + tx1231; 		/* rule 599 */
	 tx1233 = 3*K*tx1232*tx852; 		/* rule 600 */
	 tx1234 = 6.*K*tx842*tx869*tx894; 		/* rule 601 */
	 ohz1z2 = tx1233 + tx1234; 		/* rule 602 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 tx1235 = tx1068 + y2; 		/* rule 604 */
	 tx1236 = -tx1215; 		/* rule 605 */
	 tx1237 = tx1235*tx1582*tx1585*tx1588; 		/* rule 606 */
	 tx1238 = tx1585*tx1588*tx793*tx848*tx866; 		/* rule 607 */
	 tx1239 = tx1582*tx1588*tx795*tx866*tx874; 		/* rule 608 */
	 tx1240 = -(tx1585*tx796*tx810*tx839*tx848*tx905); 		/* rule 609 */
	 tx1241 = tx1582*tx1585*tx796*tx866*tx905; 		/* rule 610 */
	 tx1242 = -(tx1585*tx1588*tx810*tx848*tx906); 		/* rule 611 */
	 tx1243 = tx1109 + tx1236 + tx1237 + tx1238 + tx1239 + tx1240 + tx1241 + tx1242; 		/* rule 612 */
	 tx1244 = 3*K*tx1243*tx852; 		/* rule 613 */
	 tx1245 = 6.*K*tx842*tx869*tx911; 		/* rule 614 */
	 ohz1x3 = tx1244 + tx1245; 		/* rule 615 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 0, ohz1x3);
	 tx1246 = tx854 + x4; 		/* rule 617 */
	 tx1247 = -tx1223; 		/* rule 618 */
	 tx1248 = tx1246*tx1582*tx1585*tx1588; 		/* rule 619 */
	 tx1249 = tx1585*tx1588*tx800*tx848*tx866; 		/* rule 620 */
	 tx1250 = tx1582*tx1588*tx803*tx866*tx874; 		/* rule 621 */
	 tx1251 = tx1582*tx1585*tx805*tx866*tx905; 		/* rule 622 */
	 tx1252 = -(tx1585*tx839*tx848*tx872*tx905); 		/* rule 623 */
	 tx1253 = -(tx1585*tx1588*tx810*tx848*tx922); 		/* rule 624 */
	 tx1254 = tx1190 + tx1247 + tx1248 + tx1249 + tx1250 + tx1251 + tx1252 + tx1253; 		/* rule 625 */
	 tx1255 = 3*K*tx1254*tx852; 		/* rule 626 */
	 tx1256 = 6.*K*tx842*tx869*tx927; 		/* rule 627 */
	 ohz1y3 = tx1255 + tx1256; 		/* rule 628 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 1, ohz1y3);
	 tx1257 = -tx1229; 		/* rule 630 */
	 tx1258 = tx1585*tx1588*tx810*tx848*tx866; 		/* rule 631 */
	 tx1259 = tx1582*tx1588*tx815*tx866*tx874; 		/* rule 632 */
	 tx1260 = -(tx1585*tx810*tx818*tx839*tx848*tx905); 		/* rule 633 */
	 tx1261 = tx1582*tx1585*tx818*tx866*tx905; 		/* rule 634 */
	 tx1262 = -(tx1585*tx1588*tx810*tx848*tx938); 		/* rule 635 */
	 tx1263 = -3*tx1585*tx1588*tx826*tx839*tx953; 		/* rule 636 */
	 tx1264 = tx1078 + tx1257 + tx1258 + tx1259 + tx1260 + tx1261 + tx1262 + tx1263; 		/* rule 637 */
	 tx1265 = 3*K*tx1264*tx852; 		/* rule 638 */
	 tx1266 = 6.*K*tx842*tx869*tx943; 		/* rule 639 */
	 ohz1z3 = tx1265 + tx1266; 		/* rule 640 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 2, ohz1z3);
	 tx1267 = -tx1241; 		/* rule 642 */
	 tx1268 = -(tx1585*tx1588*tx810*tx829*tx848); 		/* rule 643 */
	 tx1269 = tx1582*tx1585*tx1588*tx863; 		/* rule 644 */
	 tx1270 = tx1585*tx796*tx810*tx839*tx848*tx905; 		/* rule 645 */
	 tx1271 = tx1267 + tx1268 + tx1269 + tx1270; 		/* rule 646 */
	 tx1272 = 3*K*tx1271*tx852; 		/* rule 647 */
	 tx1273 = 6.*K*tx842*tx869*tx946; 		/* rule 648 */
	 ohz1x4 = tx1272 + tx1273; 		/* rule 649 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 0, ohz1x4);
	 tx1274 = -tx1251; 		/* rule 651 */
	 tx1275 = tx1582*tx1585*tx1588*tx795; 		/* rule 652 */
	 tx1276 = -(tx1585*tx1588*tx810*tx828*tx848); 		/* rule 653 */
	 tx1277 = tx1585*tx839*tx848*tx872*tx905; 		/* rule 654 */
	 tx1278 = tx1274 + tx1275 + tx1276 + tx1277; 		/* rule 655 */
	 tx1279 = 3*K*tx1278*tx852; 		/* rule 656 */
	 tx1280 = 6.*K*tx842*tx869*tx949; 		/* rule 657 */
	 ohz1y4 = tx1279 + tx1280; 		/* rule 658 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 1, ohz1y4);
	 tx1281 = -tx1261; 		/* rule 660 */
	 tx1282 = -(tx1585*tx1588*tx810*tx824*tx848); 		/* rule 661 */
	 tx1283 = tx1585*tx810*tx818*tx839*tx848*tx905; 		/* rule 662 */
	 tx1284 = tx1281 + tx1282 + tx1283; 		/* rule 663 */
	 tx1285 = 3*K*tx1284*tx852; 		/* rule 664 */
	 tx1286 = 6.*K*tx842*tx869*tx952; 		/* rule 665 */
	 ohz1z4 = tx1285 + tx1286; 		/* rule 666 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 2, ohz1z4);
	 tx1287 = -(tx1582*tx1588*tx803*tx874*tx875); 		/* rule 668 */
	 tx1288 = -(tx1582*tx1588*tx795*tx874*tx883); 		/* rule 669 */
	 tx1289 = 3*tx1582*tx1588*tx795*tx803*tx839*tx973; 		/* rule 670 */
	 tx1290 = tx1287 + tx1288 + tx1289; 		/* rule 671 */
	 tx1291 = 3*K*tx1290*tx852; 		/* rule 672 */
	 tx1292 = 6.*K*tx842*tx878*tx886; 		/* rule 673 */
	 ohx2y2 = tx1291 + tx1292; 		/* rule 674 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 tx1293 = -(tx1582*tx1588*tx815*tx874*tx875); 		/* rule 676 */
	 tx1294 = -(tx1582*tx1588*tx795*tx874*tx891); 		/* rule 677 */
	 tx1295 = 3*tx1582*tx1588*tx795*tx815*tx839*tx973; 		/* rule 678 */
	 tx1296 = tx1293 + tx1294 + tx1295; 		/* rule 679 */
	 tx1297 = 3*K*tx1296*tx852; 		/* rule 680 */
	 tx1298 = 6.*K*tx842*tx878*tx894; 		/* rule 681 */
	 ohx2z2 = tx1297 + tx1298; 		/* rule 682 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 tx1299 = tx1582*tx1588*tx839*tx874; 		/* rule 684 */
	 tx1300 = tx1585*tx1588*tx793*tx848*tx875; 		/* rule 685 */
	 tx1301 = tx1582*tx1588*tx795*tx874*tx875; 		/* rule 686 */
	 tx1302 = -(tx1582*tx795*tx796*tx839*tx874*tx905); 		/* rule 687 */
	 tx1303 = tx1582*tx1585*tx796*tx875*tx905; 		/* rule 688 */
	 tx1304 = -(tx1582*tx1588*tx795*tx874*tx906); 		/* rule 689 */
	 tx1305 = -3*tx1582*tx1588*tx820*tx839*tx973; 		/* rule 690 */
	 tx1306 = tx1077 + tx1299 + tx1300 + tx1301 + tx1302 + tx1303 + tx1304 + tx1305; 		/* rule 691 */
	 tx1307 = 3*K*tx1306*tx852; 		/* rule 692 */
	 tx1308 = 6.*K*tx842*tx878*tx911; 		/* rule 693 */
	 ohx2x3 = tx1307 + tx1308; 		/* rule 694 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 0, ohx2x3);
	 tx1309 = tx879 + z4; 		/* rule 696 */
	 tx1310 = tx1309*tx1582*tx1585*tx1588; 		/* rule 697 */
	 tx1311 = tx1585*tx1588*tx800*tx848*tx875; 		/* rule 698 */
	 tx1312 = tx1582*tx1588*tx803*tx874*tx875; 		/* rule 699 */
	 tx1313 = -(tx1582*tx839*tx864*tx874*tx905); 		/* rule 700 */
	 tx1314 = tx1582*tx1585*tx805*tx875*tx905; 		/* rule 701 */
	 tx1315 = -(tx1582*tx1588*tx795*tx874*tx922); 		/* rule 702 */
	 tx1316 = -3*tx1582*tx1588*tx795*tx803*tx839*tx973; 		/* rule 703 */
	 tx1317 = tx1161 + tx1310 + tx1311 + tx1312 + tx1313 + tx1314 + tx1315 + tx1316; 		/* rule 704 */
	 tx1318 = 3*K*tx1317*tx852; 		/* rule 705 */
	 tx1319 = 6.*K*tx842*tx878*tx927; 		/* rule 706 */
	 ohx2y3 = tx1318 + tx1319; 		/* rule 707 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 1, ohx2y3);
	 tx1320 = tx1068 + y1; 		/* rule 709 */
	 tx1321 = tx1320*tx1582*tx1585*tx1588; 		/* rule 710 */
	 tx1322 = tx1585*tx1588*tx810*tx848*tx875; 		/* rule 711 */
	 tx1323 = tx1582*tx1588*tx815*tx874*tx875; 		/* rule 712 */
	 tx1324 = -(tx1582*tx795*tx818*tx839*tx874*tx905); 		/* rule 713 */
	 tx1325 = tx1582*tx1585*tx818*tx875*tx905; 		/* rule 714 */
	 tx1326 = -(tx1582*tx1588*tx795*tx874*tx938); 		/* rule 715 */
	 tx1327 = -3*tx1582*tx1588*tx795*tx815*tx839*tx973; 		/* rule 716 */
	 tx1328 = tx1236 + tx1321 + tx1322 + tx1323 + tx1324 + tx1325 + tx1326 + tx1327; 		/* rule 717 */
	 tx1329 = 3*K*tx1328*tx852; 		/* rule 718 */
	 tx1330 = 6.*K*tx842*tx878*tx943; 		/* rule 719 */
	 ohx2z3 = tx1329 + tx1330; 		/* rule 720 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 2, ohx2z3);
	 tx1331 = -tx1303; 		/* rule 722 */
	 tx1332 = -(tx1582*tx1588*tx795*tx829*tx874); 		/* rule 723 */
	 tx1333 = tx1582*tx795*tx796*tx839*tx874*tx905; 		/* rule 724 */
	 tx1334 = tx1331 + tx1332 + tx1333; 		/* rule 725 */
	 tx1335 = 3*K*tx1334*tx852; 		/* rule 726 */
	 tx1336 = 6.*K*tx842*tx878*tx946; 		/* rule 727 */
	 ohx2x4 = tx1335 + tx1336; 		/* rule 728 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 0, ohx2x4);
	 tx1337 = -tx1314; 		/* rule 730 */
	 tx1338 = tx1582*tx1585*tx1588*tx810; 		/* rule 731 */
	 tx1339 = -(tx1582*tx1588*tx795*tx828*tx874); 		/* rule 732 */
	 tx1340 = tx1582*tx839*tx864*tx874*tx905; 		/* rule 733 */
	 tx1341 = tx1337 + tx1338 + tx1339 + tx1340; 		/* rule 734 */
	 tx1342 = 3*K*tx1341*tx852; 		/* rule 735 */
	 tx1343 = 6.*K*tx842*tx878*tx949; 		/* rule 736 */
	 ohx2y4 = tx1342 + tx1343; 		/* rule 737 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 1, ohx2y4);
	 tx1344 = -tx1325; 		/* rule 739 */
	 tx1345 = tx1582*tx1585*tx1588*tx871; 		/* rule 740 */
	 tx1346 = -(tx1582*tx1588*tx795*tx824*tx874); 		/* rule 741 */
	 tx1347 = tx1582*tx795*tx818*tx839*tx874*tx905; 		/* rule 742 */
	 tx1348 = tx1344 + tx1345 + tx1346 + tx1347; 		/* rule 743 */
	 tx1349 = 3*K*tx1348*tx852; 		/* rule 744 */
	 tx1350 = 6.*K*tx842*tx878*tx952; 		/* rule 745 */
	 ohx2z4 = tx1349 + tx1350; 		/* rule 746 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 2, ohx2z4);
	 tx1351 = -(tx1582*tx1588*tx815*tx874*tx883); 		/* rule 748 */
	 tx1352 = -(tx1582*tx1588*tx803*tx874*tx891); 		/* rule 749 */
	 tx1353 = 3*tx1582*tx1588*tx803*tx815*tx839*tx973; 		/* rule 750 */
	 tx1354 = tx1351 + tx1352 + tx1353; 		/* rule 751 */
	 tx1355 = 3*K*tx1354*tx852; 		/* rule 752 */
	 tx1356 = 6.*K*tx842*tx886*tx894; 		/* rule 753 */
	 ohy2z2 = tx1355 + tx1356; 		/* rule 754 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	 tx1357 = tx1088 + z1; 		/* rule 756 */
	 tx1358 = tx1357*tx1582*tx1585*tx1588; 		/* rule 757 */
	 tx1359 = tx1585*tx1588*tx793*tx848*tx883; 		/* rule 758 */
	 tx1360 = tx1582*tx1588*tx795*tx874*tx883; 		/* rule 759 */
	 tx1361 = -(tx1582*tx796*tx803*tx839*tx874*tx905); 		/* rule 760 */
	 tx1362 = tx1582*tx1585*tx796*tx883*tx905; 		/* rule 761 */
	 tx1363 = -(tx1582*tx1588*tx803*tx874*tx906); 		/* rule 762 */
	 tx1364 = tx1090 + tx1316 + tx1358 + tx1359 + tx1360 + tx1361 + tx1362 + tx1363; 		/* rule 763 */
	 tx1365 = 3*K*tx1364*tx852; 		/* rule 764 */
	 tx1366 = 6.*K*tx842*tx886*tx911; 		/* rule 765 */
	 ohy2x3 = tx1365 + tx1366; 		/* rule 766 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 0, ohy2x3);
	 tx1367 = tx1585*tx1588*tx800*tx848*tx883; 		/* rule 768 */
	 tx1368 = tx1582*tx1588*tx803*tx874*tx883; 		/* rule 769 */
	 tx1369 = -(tx1582*tx803*tx805*tx839*tx874*tx905); 		/* rule 770 */
	 tx1370 = tx1582*tx1585*tx805*tx883*tx905; 		/* rule 771 */
	 tx1371 = -(tx1582*tx1588*tx803*tx874*tx922); 		/* rule 772 */
	 tx1372 = -3*tx1582*tx1588*tx823*tx839*tx973; 		/* rule 773 */
	 tx1373 = tx1171 + tx1299 + tx1367 + tx1368 + tx1369 + tx1370 + tx1371 + tx1372; 		/* rule 774 */
	 tx1374 = 3*K*tx1373*tx852; 		/* rule 775 */
	 tx1375 = 6.*K*tx842*tx886*tx927; 		/* rule 776 */
	 ohy2y3 = tx1374 + tx1375; 		/* rule 777 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 1, ohy2y3);
	 tx1376 = tx887 + x4; 		/* rule 779 */
	 tx1377 = tx1376*tx1582*tx1585*tx1588; 		/* rule 780 */
	 tx1378 = tx1585*tx1588*tx810*tx848*tx883; 		/* rule 781 */
	 tx1379 = tx1582*tx1588*tx815*tx874*tx883; 		/* rule 782 */
	 tx1380 = -(tx1582*tx839*tx846*tx874*tx905); 		/* rule 783 */
	 tx1381 = tx1582*tx1585*tx818*tx883*tx905; 		/* rule 784 */
	 tx1382 = -(tx1582*tx1588*tx803*tx874*tx938); 		/* rule 785 */
	 tx1383 = -3*tx1582*tx1588*tx803*tx815*tx839*tx973; 		/* rule 786 */
	 tx1384 = tx1247 + tx1377 + tx1378 + tx1379 + tx1380 + tx1381 + tx1382 + tx1383; 		/* rule 787 */
	 tx1385 = 3*K*tx1384*tx852; 		/* rule 788 */
	 tx1386 = 6.*K*tx842*tx886*tx943; 		/* rule 789 */
	 ohy2z3 = tx1385 + tx1386; 		/* rule 790 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 2, ohy2z3);
	 tx1387 = -tx1362; 		/* rule 792 */
	 tx1388 = -(tx1582*tx1588*tx803*tx829*tx874); 		/* rule 793 */
	 tx1389 = tx1582*tx1585*tx1588*tx880; 		/* rule 794 */
	 tx1390 = tx1582*tx796*tx803*tx839*tx874*tx905; 		/* rule 795 */
	 tx1391 = tx1387 + tx1388 + tx1389 + tx1390; 		/* rule 796 */
	 tx1392 = 3*K*tx1391*tx852; 		/* rule 797 */
	 tx1393 = 6.*K*tx842*tx886*tx946; 		/* rule 798 */
	 ohy2x4 = tx1392 + tx1393; 		/* rule 799 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 0, ohy2x4);
	 tx1394 = -tx1370; 		/* rule 801 */
	 tx1395 = -(tx1582*tx1588*tx803*tx828*tx874); 		/* rule 802 */
	 tx1396 = tx1582*tx803*tx805*tx839*tx874*tx905; 		/* rule 803 */
	 tx1397 = tx1394 + tx1395 + tx1396; 		/* rule 804 */
	 tx1398 = 3*K*tx1397*tx852; 		/* rule 805 */
	 tx1399 = 6.*K*tx842*tx886*tx949; 		/* rule 806 */
	 ohy2y4 = tx1398 + tx1399; 		/* rule 807 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 1, ohy2y4);
	 tx1400 = -tx1381; 		/* rule 809 */
	 tx1401 = tx1582*tx1585*tx1588*tx793; 		/* rule 810 */
	 tx1402 = -(tx1582*tx1588*tx803*tx824*tx874); 		/* rule 811 */
	 tx1403 = tx1582*tx839*tx846*tx874*tx905; 		/* rule 812 */
	 tx1404 = tx1400 + tx1401 + tx1402 + tx1403; 		/* rule 813 */
	 tx1405 = 3*K*tx1404*tx852; 		/* rule 814 */
	 tx1406 = 6.*K*tx842*tx886*tx952; 		/* rule 815 */
	 ohy2z4 = tx1405 + tx1406; 		/* rule 816 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 2, ohy2z4);
	 tx1407 = tx870 + y4; 		/* rule 818 */
	 tx1408 = tx1407*tx1582*tx1585*tx1588; 		/* rule 819 */
	 tx1409 = tx1585*tx1588*tx793*tx848*tx891; 		/* rule 820 */
	 tx1410 = tx1582*tx1588*tx795*tx874*tx891; 		/* rule 821 */
	 tx1411 = -(tx1582*tx839*tx856*tx874*tx905); 		/* rule 822 */
	 tx1412 = tx1582*tx1585*tx796*tx891*tx905; 		/* rule 823 */
	 tx1413 = -(tx1582*tx1588*tx815*tx874*tx906); 		/* rule 824 */
	 tx1414 = tx1102 + tx1327 + tx1408 + tx1409 + tx1410 + tx1411 + tx1412 + tx1413; 		/* rule 825 */
	 tx1415 = 3*K*tx1414*tx852; 		/* rule 826 */
	 tx1416 = 6.*K*tx842*tx894*tx911; 		/* rule 827 */
	 ohz2x3 = tx1415 + tx1416; 		/* rule 828 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 0, ohz2x3);
	 tx1417 = tx1181 + x1; 		/* rule 830 */
	 tx1418 = tx1417*tx1582*tx1585*tx1588; 		/* rule 831 */
	 tx1419 = tx1585*tx1588*tx800*tx848*tx891; 		/* rule 832 */
	 tx1420 = tx1582*tx1588*tx803*tx874*tx891; 		/* rule 833 */
	 tx1421 = -(tx1582*tx805*tx815*tx839*tx874*tx905); 		/* rule 834 */
	 tx1422 = tx1582*tx1585*tx805*tx891*tx905; 		/* rule 835 */
	 tx1423 = -(tx1582*tx1588*tx815*tx874*tx922); 		/* rule 836 */
	 tx1424 = tx1183 + tx1383 + tx1418 + tx1419 + tx1420 + tx1421 + tx1422 + tx1423; 		/* rule 837 */
	 tx1425 = 3*K*tx1424*tx852; 		/* rule 838 */
	 tx1426 = 6.*K*tx842*tx894*tx927; 		/* rule 839 */
	 ohz2y3 = tx1425 + tx1426; 		/* rule 840 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 1, ohz2y3);
	 tx1427 = tx1585*tx1588*tx810*tx848*tx891; 		/* rule 842 */
	 tx1428 = tx1582*tx1588*tx815*tx874*tx891; 		/* rule 843 */
	 tx1429 = -(tx1582*tx815*tx818*tx839*tx874*tx905); 		/* rule 844 */
	 tx1430 = tx1582*tx1585*tx818*tx891*tx905; 		/* rule 845 */
	 tx1431 = -(tx1582*tx1588*tx815*tx874*tx938); 		/* rule 846 */
	 tx1432 = -3*tx1582*tx1588*tx827*tx839*tx973; 		/* rule 847 */
	 tx1433 = tx1257 + tx1299 + tx1427 + tx1428 + tx1429 + tx1430 + tx1431 + tx1432; 		/* rule 848 */
	 tx1434 = 3*K*tx1433*tx852; 		/* rule 849 */
	 tx1435 = 6.*K*tx842*tx894*tx943; 		/* rule 850 */
	 ohz2z3 = tx1434 + tx1435; 		/* rule 851 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 2, ohz2z3);
	 tx1436 = -tx1412; 		/* rule 853 */
	 tx1437 = tx1582*tx1585*tx1588*tx800; 		/* rule 854 */
	 tx1438 = -(tx1582*tx1588*tx815*tx829*tx874); 		/* rule 855 */
	 tx1439 = tx1582*tx839*tx856*tx874*tx905; 		/* rule 856 */
	 tx1440 = tx1436 + tx1437 + tx1438 + tx1439; 		/* rule 857 */
	 tx1441 = 3*K*tx1440*tx852; 		/* rule 858 */
	 tx1442 = 6.*K*tx842*tx894*tx946; 		/* rule 859 */
	 ohz2x4 = tx1441 + tx1442; 		/* rule 860 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 0, ohz2x4);
	 tx1443 = -tx1422; 		/* rule 862 */
	 tx1444 = -(tx1582*tx1588*tx815*tx828*tx874); 		/* rule 863 */
	 tx1445 = tx1582*tx1585*tx1588*tx888; 		/* rule 864 */
	 tx1446 = tx1582*tx805*tx815*tx839*tx874*tx905; 		/* rule 865 */
	 tx1447 = tx1443 + tx1444 + tx1445 + tx1446; 		/* rule 866 */
	 tx1448 = 3*K*tx1447*tx852; 		/* rule 867 */
	 tx1449 = 6.*K*tx842*tx894*tx949; 		/* rule 868 */
	 ohz2y4 = tx1448 + tx1449; 		/* rule 869 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 1, ohz2y4);
	 tx1450 = -tx1430; 		/* rule 871 */
	 tx1451 = -(tx1582*tx1588*tx815*tx824*tx874); 		/* rule 872 */
	 tx1452 = tx1582*tx815*tx818*tx839*tx874*tx905; 		/* rule 873 */
	 tx1453 = tx1450 + tx1451 + tx1452; 		/* rule 874 */
	 tx1454 = 3*K*tx1453*tx852; 		/* rule 875 */
	 tx1455 = 6.*K*tx842*tx894*tx952; 		/* rule 876 */
	 ohz2z4 = tx1454 + tx1455; 		/* rule 877 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 2, ohz2z4);
	 tx1456 = tx1585*tx1588*tx800*tx848*tx906; 		/* rule 879 */
	 tx1457 = tx1582*tx1588*tx803*tx874*tx906; 		/* rule 880 */
	 tx1458 = tx1582*tx1585*tx805*tx905*tx906; 		/* rule 881 */
	 tx1459 = tx1585*tx1588*tx793*tx848*tx922; 		/* rule 882 */
	 tx1460 = tx1582*tx1588*tx795*tx874*tx922; 		/* rule 883 */
	 tx1461 = tx1582*tx1585*tx796*tx905*tx922; 		/* rule 884 */
	 tx1462 = 3*tx1582*tx1585*tx796*tx805*tx839*tx993; 		/* rule 885 */
	 tx1463 = tx1045 + tx1062 + tx1122 + tx1141 + tx1197 + tx1289 + tx1340 + tx1390 + tx1456 + tx1457 + tx1458 + tx1459 + tx1460 + tx1461 + tx1462; 		/* rule 886 */
	 tx1464 = 3*K*tx1463*tx852; 		/* rule 887 */
	 tx1465 = 6.*K*tx842*tx911*tx927; 		/* rule 888 */
	 ohx3y3 = tx1464 + tx1465; 		/* rule 889 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 1, ohx3y3);
	 tx1466 = tx1585*tx1588*tx810*tx848*tx906; 		/* rule 891 */
	 tx1467 = tx1582*tx1588*tx815*tx874*tx906; 		/* rule 892 */
	 tx1468 = tx1582*tx1585*tx818*tx905*tx906; 		/* rule 893 */
	 tx1469 = tx1585*tx1588*tx793*tx848*tx938; 		/* rule 894 */
	 tx1470 = tx1582*tx1588*tx795*tx874*tx938; 		/* rule 895 */
	 tx1471 = tx1582*tx1585*tx796*tx905*tx938; 		/* rule 896 */
	 tx1472 = 3*tx1582*tx1585*tx796*tx818*tx839*tx993; 		/* rule 897 */
	 tx1473 = tx1051 + tx1071 + tx1129 + tx1215 + tx1270 + tx1295 + tx1347 + tx1439 + tx1466 + tx1467 + tx1468 + tx1469 + tx1470 + tx1471 + tx1472; 		/* rule 898 */
	 tx1474 = 3*K*tx1473*tx852; 		/* rule 899 */
	 tx1475 = 6.*K*tx842*tx911*tx943; 		/* rule 900 */
	 ohx3z3 = tx1474 + tx1475; 		/* rule 901 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 2, ohx3z3);
	 tx1476 = tx1585*tx1588*tx793*tx829*tx848; 		/* rule 903 */
	 tx1477 = tx1582*tx1588*tx795*tx829*tx874; 		/* rule 904 */
	 tx1478 = tx1582*tx1585*tx835*tx905; 		/* rule 905 */
	 tx1479 = tx1582*tx1585*tx839*tx905; 		/* rule 906 */
	 tx1480 = -(tx1582*tx1585*tx796*tx905*tx906); 		/* rule 907 */
	 tx1481 = -3*tx1582*tx1585*tx821*tx839*tx993; 		/* rule 908 */
	 tx1482 = tx1081 + tx1302 + tx1476 + tx1477 + tx1478 + tx1479 + tx1480 + tx1481; 		/* rule 909 */
	 tx1483 = 3*K*tx1482*tx852; 		/* rule 910 */
	 tx1484 = 6.*K*tx842*tx911*tx946; 		/* rule 911 */
	 ohx3x4 = tx1483 + tx1484; 		/* rule 912 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 0, ohx3x4);
	 tx1485 = -tx1458; 		/* rule 914 */
	 tx1486 = tx1585*tx1588*tx793*tx828*tx848; 		/* rule 915 */
	 tx1487 = tx1582*tx1588*tx795*tx828*tx874; 		/* rule 916 */
	 tx1488 = tx1582*tx1585*tx1588*tx896; 		/* rule 917 */
	 tx1489 = tx1582*tx1585*tx796*tx828*tx905; 		/* rule 918 */
	 tx1490 = -3*tx1582*tx1585*tx796*tx805*tx839*tx993; 		/* rule 919 */
	 tx1491 = tx1094 + tx1313 + tx1485 + tx1486 + tx1487 + tx1488 + tx1489 + tx1490; 		/* rule 920 */
	 tx1492 = 3*K*tx1491*tx852; 		/* rule 921 */
	 tx1493 = 6.*K*tx842*tx911*tx949; 		/* rule 922 */
	 ohx3y4 = tx1492 + tx1493; 		/* rule 923 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 1, ohx3y4);
	 tx1494 = -tx1468; 		/* rule 925 */
	 tx1495 = tx1585*tx1588*tx793*tx824*tx848; 		/* rule 926 */
	 tx1496 = tx1582*tx1588*tx795*tx824*tx874; 		/* rule 927 */
	 tx1497 = tx1582*tx1585*tx1588*tx895; 		/* rule 928 */
	 tx1498 = tx1582*tx1585*tx796*tx824*tx905; 		/* rule 929 */
	 tx1499 = -3*tx1582*tx1585*tx796*tx818*tx839*tx993; 		/* rule 930 */
	 tx1500 = tx1107 + tx1324 + tx1494 + tx1495 + tx1496 + tx1497 + tx1498 + tx1499; 		/* rule 931 */
	 tx1501 = 3*K*tx1500*tx852; 		/* rule 932 */
	 tx1502 = 6.*K*tx842*tx911*tx952; 		/* rule 933 */
	 ohx3z4 = tx1501 + tx1502; 		/* rule 934 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 2, ohx3z4);
	 tx1503 = tx1585*tx1588*tx810*tx848*tx922; 		/* rule 936 */
	 tx1504 = tx1582*tx1588*tx815*tx874*tx922; 		/* rule 937 */
	 tx1505 = tx1582*tx1585*tx818*tx905*tx922; 		/* rule 938 */
	 tx1506 = tx1585*tx1588*tx800*tx848*tx938; 		/* rule 939 */
	 tx1507 = tx1582*tx1588*tx803*tx874*tx938; 		/* rule 940 */
	 tx1508 = tx1582*tx1585*tx805*tx905*tx938; 		/* rule 941 */
	 tx1509 = 3*tx1582*tx1585*tx805*tx818*tx839*tx993; 		/* rule 942 */
	 tx1510 = tx1135 + tx1154 + tx1210 + tx1223 + tx1277 + tx1353 + tx1403 + tx1446 + tx1503 + tx1504 + tx1505 + tx1506 + tx1507 + tx1508 + tx1509; 		/* rule 943 */
	 tx1511 = 3*K*tx1510*tx852; 		/* rule 944 */
	 tx1512 = 6.*K*tx842*tx927*tx943; 		/* rule 945 */
	 ohy3z3 = tx1511 + tx1512; 		/* rule 946 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 2, ohy3z3);
	 tx1513 = -tx1461; 		/* rule 948 */
	 tx1514 = tx1585*tx1588*tx800*tx829*tx848; 		/* rule 949 */
	 tx1515 = tx1582*tx1588*tx803*tx829*tx874; 		/* rule 950 */
	 tx1516 = tx1582*tx1585*tx805*tx829*tx905; 		/* rule 951 */
	 tx1517 = tx1582*tx1585*tx1588*tx913; 		/* rule 952 */
	 tx1518 = tx1166 + tx1361 + tx1490 + tx1513 + tx1514 + tx1515 + tx1516 + tx1517; 		/* rule 953 */
	 tx1519 = 3*K*tx1518*tx852; 		/* rule 954 */
	 tx1520 = 6.*K*tx842*tx927*tx946; 		/* rule 955 */
	 ohy3x4 = tx1519 + tx1520; 		/* rule 956 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 0, ohy3x4);
	 tx1521 = tx1585*tx1588*tx800*tx828*tx848; 		/* rule 958 */
	 tx1522 = tx1582*tx1588*tx803*tx828*tx874; 		/* rule 959 */
	 tx1523 = tx1582*tx1585*tx834*tx905; 		/* rule 960 */
	 tx1524 = -(tx1582*tx1585*tx805*tx905*tx922); 		/* rule 961 */
	 tx1525 = -3*tx1582*tx1585*tx825*tx839*tx993; 		/* rule 962 */
	 tx1526 = tx1174 + tx1369 + tx1479 + tx1521 + tx1522 + tx1523 + tx1524 + tx1525; 		/* rule 963 */
	 tx1527 = 3*K*tx1526*tx852; 		/* rule 964 */
	 tx1528 = 6.*K*tx842*tx927*tx949; 		/* rule 965 */
	 ohy3y4 = tx1527 + tx1528; 		/* rule 966 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 1, ohy3y4);
	 tx1529 = -tx1505; 		/* rule 968 */
	 tx1530 = tx1585*tx1588*tx800*tx824*tx848; 		/* rule 969 */
	 tx1531 = tx1582*tx1588*tx803*tx824*tx874; 		/* rule 970 */
	 tx1532 = tx1582*tx1585*tx805*tx824*tx905; 		/* rule 971 */
	 tx1533 = tx1582*tx1585*tx1588*tx912; 		/* rule 972 */
	 tx1534 = -3*tx1582*tx1585*tx805*tx818*tx839*tx993; 		/* rule 973 */
	 tx1535 = tx1187 + tx1380 + tx1529 + tx1530 + tx1531 + tx1532 + tx1533 + tx1534; 		/* rule 974 */
	 tx1536 = 3*K*tx1535*tx852; 		/* rule 975 */
	 tx1537 = 6.*K*tx842*tx927*tx952; 		/* rule 976 */
	 ohy3z4 = tx1536 + tx1537; 		/* rule 977 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 2, ohy3z4);
	 tx1538 = -tx1471; 		/* rule 979 */
	 tx1539 = tx1585*tx1588*tx810*tx829*tx848; 		/* rule 980 */
	 tx1540 = tx1582*tx1588*tx815*tx829*tx874; 		/* rule 981 */
	 tx1541 = tx1582*tx1585*tx818*tx829*tx905; 		/* rule 982 */
	 tx1542 = tx1582*tx1585*tx1588*tx929; 		/* rule 983 */
	 tx1543 = tx1240 + tx1411 + tx1499 + tx1538 + tx1539 + tx1540 + tx1541 + tx1542; 		/* rule 984 */
	 tx1544 = 3*K*tx1543*tx852; 		/* rule 985 */
	 tx1545 = 6.*K*tx842*tx943*tx946; 		/* rule 986 */
	 ohz3x4 = tx1544 + tx1545; 		/* rule 987 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 0, ohz3x4);
	 tx1546 = -tx1508; 		/* rule 989 */
	 tx1547 = tx1585*tx1588*tx810*tx828*tx848; 		/* rule 990 */
	 tx1548 = tx1582*tx1588*tx815*tx828*tx874; 		/* rule 991 */
	 tx1549 = tx1582*tx1585*tx818*tx828*tx905; 		/* rule 992 */
	 tx1550 = tx1582*tx1585*tx1588*tx928; 		/* rule 993 */
	 tx1551 = tx1252 + tx1421 + tx1534 + tx1546 + tx1547 + tx1548 + tx1549 + tx1550; 		/* rule 994 */
	 tx1552 = 3*K*tx1551*tx852; 		/* rule 995 */
	 tx1553 = 6.*K*tx842*tx943*tx949; 		/* rule 996 */
	 ohz3y4 = tx1552 + tx1553; 		/* rule 997 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 1, ohz3y4);
	 tx1554 = tx1585*tx1588*tx810*tx824*tx848; 		/* rule 999 */
	 tx1555 = tx1582*tx1588*tx815*tx824*tx874; 		/* rule 1000 */
	 tx1556 = tx1582*tx1585*tx831*tx905; 		/* rule 1001 */
	 tx1557 = -(tx1582*tx1585*tx818*tx905*tx938); 		/* rule 1002 */
	 tx1558 = -3*tx1582*tx1585*tx830*tx839*tx993; 		/* rule 1003 */
	 tx1559 = tx1260 + tx1429 + tx1479 + tx1554 + tx1555 + tx1556 + tx1557 + tx1558; 		/* rule 1004 */
	 tx1560 = 3*K*tx1559*tx852; 		/* rule 1005 */
	 tx1561 = 6.*K*tx842*tx943*tx952; 		/* rule 1006 */
	 ohz3z4 = tx1560 + tx1561; 		/* rule 1007 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 2, ohz3z4);
	 tx1562 = -tx1489; 		/* rule 1009 */
	 tx1563 = -tx1516; 		/* rule 1010 */
	 tx1564 = tx1462 + tx1562 + tx1563; 		/* rule 1011 */
	 tx1565 = 3*K*tx1564*tx852; 		/* rule 1012 */
	 tx1566 = 6.*K*tx842*tx946*tx949; 		/* rule 1013 */
	 ohx4y4 = tx1565 + tx1566; 		/* rule 1014 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 1, ohx4y4);
	 tx1567 = -tx1498; 		/* rule 1016 */
	 tx1568 = -tx1541; 		/* rule 1017 */
	 tx1569 = tx1472 + tx1567 + tx1568; 		/* rule 1018 */
	 tx1570 = 3*K*tx1569*tx852; 		/* rule 1019 */
	 tx1571 = 6.*K*tx842*tx946*tx952; 		/* rule 1020 */
	 ohx4z4 = tx1570 + tx1571; 		/* rule 1021 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 2, ohx4z4);
	 tx1572 = -tx1532; 		/* rule 1023 */
	 tx1573 = -tx1549; 		/* rule 1024 */
	 tx1574 = tx1509 + tx1572 + tx1573; 		/* rule 1025 */
	 tx1575 = 3*K*tx1574*tx852; 		/* rule 1026 */
	 tx1576 = 6.*K*tx842*tx949*tx952; 		/* rule 1027 */
	 ohy4z4 = tx1575 + tx1576; 		/* rule 1028 */
	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 1, I4, 2, ohy4z4);
	} /*if calcOffDiagonalHessian */ 
	#endif /* CHIRAL_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* CHIRAL_RESTRAINT_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* CHIRAL_RESTRAINT_CALC_FORCE ]*/
	} /*ChiralTest>0.0*/
