// C-code
	DIHEDRAL_SET_PARAMETER(sinPhase);
	DIHEDRAL_SET_PARAMETER(cosPhase);
	DIHEDRAL_SET_PARAMETER(V);
	DIHEDRAL_SET_PARAMETER(DN);
	DIHEDRAL_SET_PARAMETER(IN);
	DIHEDRAL_SET_PARAMETER(I1);
	DIHEDRAL_SET_PARAMETER(I2);
	DIHEDRAL_SET_PARAMETER(I3);
	DIHEDRAL_SET_PARAMETER(I4);
	DIHEDRAL_SET_POSITION(x1,I1,0);
	DIHEDRAL_SET_POSITION(y1,I1,1);
	DIHEDRAL_SET_POSITION(z1,I1,2);
	DIHEDRAL_SET_POSITION(x2,I2,0);
	DIHEDRAL_SET_POSITION(y2,I2,1);
	DIHEDRAL_SET_POSITION(z2,I2,2);
	DIHEDRAL_SET_POSITION(x3,I3,0);
	DIHEDRAL_SET_POSITION(y3,I3,1);
	DIHEDRAL_SET_POSITION(z3,I3,2);
	DIHEDRAL_SET_POSITION(x4,I4,0);
	DIHEDRAL_SET_POSITION(y4,I4,1);
	DIHEDRAL_SET_POSITION(z4,I4,2);
	 tx907 = -(x2*y1); 		/* rule 22 */
	 tx908 = x3*y1; 		/* rule 23 */
	 tx909 = x1*y2; 		/* rule 24 */
	 tx910 = -(x3*y2); 		/* rule 25 */
	 tx911 = -(x1*y3); 		/* rule 26 */
	 tx912 = x2*y3; 		/* rule 27 */
	 tx913 = x2*z1; 		/* rule 28 */
	 tx914 = -(x3*z1); 		/* rule 29 */
	 tx915 = -(y2*z1); 		/* rule 30 */
	 tx916 = y3*z1; 		/* rule 31 */
	 tx917 = -(x1*z2); 		/* rule 32 */
	 tx918 = x3*z2; 		/* rule 33 */
	 tx919 = y1*z2; 		/* rule 34 */
	 tx920 = -(y3*z2); 		/* rule 35 */
	 tx921 = x1*z3; 		/* rule 36 */
	 tx922 = -(x2*z3); 		/* rule 37 */
	 tx923 = -(y1*z3); 		/* rule 38 */
	 tx924 = y2*z3; 		/* rule 39 */
	 tx925 = tx907 + tx908 + tx909 + tx910 + tx911 + tx912; 		/* rule 40 */
	 tx926 = tx913 + tx914 + tx917 + tx918 + tx921 + tx922; 		/* rule 41 */
	 tx927 = tx915 + tx916 + tx919 + tx920 + tx923 + tx924; 		/* rule 42 */
	 tx928 = power2(tx925); 		/* rule 43 */
	 tx929 = power2(tx926); 		/* rule 44 */
	 tx930 = power2(tx927); 		/* rule 45 */
	 tx931 = tx928 + tx929 + tx930; 		/* rule 46 */
	 LenA = mysqrt(tx931); 		/* rule 47 */
	 tx932 = x4*y2; 		/* rule 48 */
	 tx933 = -(x4*y3); 		/* rule 49 */
	 tx934 = -(x2*y4); 		/* rule 50 */
	 tx935 = x3*y4; 		/* rule 51 */
	 tx936 = -(x4*z2); 		/* rule 52 */
	 tx937 = y4*z2; 		/* rule 53 */
	 tx938 = x4*z3; 		/* rule 54 */
	 tx939 = -(y4*z3); 		/* rule 55 */
	 tx940 = x2*z4; 		/* rule 56 */
	 tx941 = -(x3*z4); 		/* rule 57 */
	 tx942 = -(y2*z4); 		/* rule 58 */
	 tx943 = y3*z4; 		/* rule 59 */
	 tx944 = tx910 + tx912 + tx932 + tx933 + tx934 + tx935; 		/* rule 60 */
	 tx945 = tx918 + tx922 + tx936 + tx938 + tx940 + tx941; 		/* rule 61 */
	 tx946 = tx920 + tx924 + tx937 + tx939 + tx942 + tx943; 		/* rule 62 */
	 tx947 = power2(tx944); 		/* rule 63 */
	 tx948 = power2(tx945); 		/* rule 64 */
	 tx949 = power2(tx946); 		/* rule 65 */
	 tx950 = tx947 + tx948 + tx949; 		/* rule 66 */
	 LenB = mysqrt(tx950); 		/* rule 67 */
	 ReciprocalLenA = reciprocal(LenA); 		/* rule 68 */
	 ReciprocalLenB = reciprocal(LenB); 		/* rule 69 */
	if (fabs(LenA)<TENM3) ReciprocalLenA = 0.0;
	if (fabs(LenB)<TENM3) ReciprocalLenB = 0.0;
	 RecLenARecLenB = ReciprocalLenA*ReciprocalLenB; 		/* rule 72 */
	EraseLinearDihedral = 1.0;
	if (RecLenARecLenB==0.0) EraseLinearDihedral = 0.0;
	 tx951 = tx925*tx944; 		/* rule 75 */
	 tx952 = tx926*tx945; 		/* rule 76 */
	 tx953 = tx927*tx946; 		/* rule 77 */
	 tx954 = tx951 + tx952 + tx953; 		/* rule 78 */
	 CosPhi = RecLenARecLenB*tx954; 		/* rule 79 */
	 tx955 = -x3; 		/* rule 80 */
	 tx956 = -y3; 		/* rule 81 */
	 tx957 = -z3; 		/* rule 82 */
	 tx958 = tx955 + x2; 		/* rule 83 */
	 tx959 = tx956 + y2; 		/* rule 84 */
	 tx960 = tx957 + z2; 		/* rule 85 */
	 tx961 = power2(tx958); 		/* rule 86 */
	 tx962 = power2(tx959); 		/* rule 87 */
	 tx963 = power2(tx960); 		/* rule 88 */
	 tx964 = tx955 + x4; 		/* rule 89 */
	 tx965 = tx956 + y4; 		/* rule 90 */
	 tx966 = tx957 + z4; 		/* rule 91 */
	 tx967 = tx961 + tx962 + tx963; 		/* rule 92 */
	 tx968 = tx927*tx964; 		/* rule 93 */
	 tx969 = tx926*tx965; 		/* rule 94 */
	 tx970 = tx925*tx966; 		/* rule 95 */
	 tx971 = mysqrt(tx967); 		/* rule 96 */
	 tx972 = tx968 + tx969 + tx970; 		/* rule 97 */
	 SinPhi = RecLenARecLenB*tx971*tx972; 		/* rule 98 */
	CosPhi=MAX(-1.0,MIN(1.0,CosPhi));
/*CosNPhi = mathCosNPhi[IN,SinPhi,CosPhi];*/
/*SinNPhi = mathSinNPhi[IN,SinPhi,CosPhi];*/
	sinNPhiCosNPhi(IN, &SinNPhi, &CosNPhi, SinPhi, CosPhi);
	 tx973 = CosNPhi*cosPhase; 		/* rule 103 */
	 tx974 = SinNPhi*sinPhase; 		/* rule 104 */
	 DihedralDeviation = 1. + tx973 + tx974; 		/* rule 105 */
	 Energy = DihedralDeviation*EraseLinearDihedral*V; 		/* rule 106 */
	DIHEDRAL_ENERGY_ACCUMULATE(Energy);
	#ifdef DIHEDRAL_CALC_FORCE //[
	if (calcForce ) {
	 tx975 = cosPhase*SinNPhi; 		/* rule 110 */
	 tx976 = -(CosNPhi*sinPhase); 		/* rule 111 */
	 tx977 = tx975 + tx976; 		/* rule 112 */
	 DeDPhi = -(DN*EraseLinearDihedral*tx977*V); 		/* rule 113 */
	 tx978 = reciprocal(tx931); 		/* rule 114 */
	 gx1 = -(DeDPhi*tx927*tx971*tx978); 		/* rule 115 */
	 fx1 = -gx1; 		/* rule 116 */
	DIHEDRAL_FORCE_ACCUMULATE(I1, 0, fx1 );
	 gy1 = -(DeDPhi*tx926*tx971*tx978); 		/* rule 118 */
	 fy1 = -gy1; 		/* rule 119 */
	DIHEDRAL_FORCE_ACCUMULATE(I1, 1, fy1 );
	 gz1 = -(DeDPhi*tx925*tx971*tx978); 		/* rule 121 */
	 fz1 = -gz1; 		/* rule 122 */
	DIHEDRAL_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tx979 = -x2; 		/* rule 124 */
	 tx980 = -y2; 		/* rule 125 */
	 tx981 = -z2; 		/* rule 126 */
	 tx982 = tx979 + x1; 		/* rule 127 */
	 tx983 = tx980 + y1; 		/* rule 128 */
	 tx984 = tx981 + z1; 		/* rule 129 */
	 tx985 = tx958*tx964; 		/* rule 130 */
	 tx986 = tx959*tx965; 		/* rule 131 */
	 tx987 = tx960*tx966; 		/* rule 132 */
	 tx988 = tx958*tx982; 		/* rule 133 */
	 tx989 = tx959*tx983; 		/* rule 134 */
	 tx990 = tx960*tx984; 		/* rule 135 */
	 tx991 = reciprocal(tx950); 		/* rule 136 */
	 tx992 = reciprocal(tx971); 		/* rule 137 */
	 tx993 = tx985 + tx986 + tx987; 		/* rule 138 */
	 tx994 = tx988 + tx989 + tx990; 		/* rule 139 */
	 tx995 = tx927*tx971*tx978; 		/* rule 140 */
	 tx996 = -(tx946*tx991*tx992*tx993); 		/* rule 141 */
	 tx997 = tx927*tx978*tx992*tx994; 		/* rule 142 */
	 tx998 = tx995 + tx996 + tx997; 		/* rule 143 */
	 gx2 = DeDPhi*tx998; 		/* rule 144 */
	 fx2 = -gx2; 		/* rule 145 */
	DIHEDRAL_FORCE_ACCUMULATE(I2, 0, fx2 );
	 tx999 = tx926*tx971*tx978; 		/* rule 147 */
	 tx1000 = -(tx945*tx991*tx992*tx993); 		/* rule 148 */
	 tx1001 = tx926*tx978*tx992*tx994; 		/* rule 149 */
	 tx1002 = tx1000 + tx1001 + tx999; 		/* rule 150 */
	 gy2 = DeDPhi*tx1002; 		/* rule 151 */
	 fy2 = -gy2; 		/* rule 152 */
	DIHEDRAL_FORCE_ACCUMULATE(I2, 1, fy2 );
	 tx1003 = tx925*tx971*tx978; 		/* rule 154 */
	 tx1004 = -(tx944*tx991*tx992*tx993); 		/* rule 155 */
	 tx1005 = tx925*tx978*tx992*tx994; 		/* rule 156 */
	 tx1006 = tx1003 + tx1004 + tx1005; 		/* rule 157 */
	 gz2 = DeDPhi*tx1006; 		/* rule 158 */
	 fz2 = -gz2; 		/* rule 159 */
	DIHEDRAL_FORCE_ACCUMULATE(I2, 2, fz2 );
	 tx1007 = -(tx946*tx971*tx991); 		/* rule 161 */
	 tx1008 = tx946*tx991*tx992*tx993; 		/* rule 162 */
	 tx1009 = -tx997; 		/* rule 163 */
	 tx1010 = tx1007 + tx1008 + tx1009; 		/* rule 164 */
	 gx3 = DeDPhi*tx1010; 		/* rule 165 */
	 fx3 = -gx3; 		/* rule 166 */
	DIHEDRAL_FORCE_ACCUMULATE(I3, 0, fx3 );
	 tx1011 = -tx1001; 		/* rule 168 */
	 tx1012 = -(tx945*tx971*tx991); 		/* rule 169 */
	 tx1013 = tx945*tx991*tx992*tx993; 		/* rule 170 */
	 tx1014 = tx1011 + tx1012 + tx1013; 		/* rule 171 */
	 gy3 = DeDPhi*tx1014; 		/* rule 172 */
	 fy3 = -gy3; 		/* rule 173 */
	DIHEDRAL_FORCE_ACCUMULATE(I3, 1, fy3 );
	 tx1015 = -tx1005; 		/* rule 175 */
	 tx1016 = -(tx944*tx971*tx991); 		/* rule 176 */
	 tx1017 = tx944*tx991*tx992*tx993; 		/* rule 177 */
	 tx1018 = tx1015 + tx1016 + tx1017; 		/* rule 178 */
	 gz3 = DeDPhi*tx1018; 		/* rule 179 */
	 fz3 = -gz3; 		/* rule 180 */
	DIHEDRAL_FORCE_ACCUMULATE(I3, 2, fz3 );
	 gx4 = DeDPhi*tx946*tx971*tx991; 		/* rule 182 */
	 fx4 = -gx4; 		/* rule 183 */
	DIHEDRAL_FORCE_ACCUMULATE(I4, 0, fx4 );
	 gy4 = DeDPhi*tx945*tx971*tx991; 		/* rule 185 */
	 fy4 = -gy4; 		/* rule 186 */
	DIHEDRAL_FORCE_ACCUMULATE(I4, 1, fy4 );
	 gz4 = DeDPhi*tx944*tx971*tx991; 		/* rule 188 */
	 fz4 = -gz4; 		/* rule 189 */
	DIHEDRAL_FORCE_ACCUMULATE(I4, 2, fz4 );
	#ifdef DIHEDRAL_CALC_DIAGONAL_HESSIAN //[
	if (calcDiagonalHessian) {
	 tx1019 = power2(y2); 		/* rule 193 */
	 tx1020 = power2(y3); 		/* rule 194 */
	 tx1021 = power2(z2); 		/* rule 195 */
	 tx1022 = power2(z3); 		/* rule 196 */
	 tx1023 = tx1019*tx955; 		/* rule 197 */
	 tx1024 = tx1021*tx955; 		/* rule 198 */
	 tx1025 = tx908*tx956; 		/* rule 199 */
	 tx1026 = tx1020*tx979; 		/* rule 200 */
	 tx1027 = tx1022*tx979; 		/* rule 201 */
	 tx1028 = tx913*tx981; 		/* rule 202 */
	 tx1029 = tx1019*x1; 		/* rule 203 */
	 tx1030 = tx1020*x1; 		/* rule 204 */
	 tx1031 = tx1021*x1; 		/* rule 205 */
	 tx1032 = tx1022*x1; 		/* rule 206 */
	 tx1033 = tx912*y1; 		/* rule 207 */
	 tx1034 = tx907*y2; 		/* rule 208 */
	 tx1035 = tx908*y2; 		/* rule 209 */
	 tx1036 = tx912*y2; 		/* rule 210 */
	 tx1037 = -2.*tx909*y3; 		/* rule 211 */
	 tx1038 = x3*y2*y3; 		/* rule 212 */
	 tx1039 = tx918*z1; 		/* rule 213 */
	 tx1040 = -2.*tx921*z2; 		/* rule 214 */
	 tx1041 = tx913*z3; 		/* rule 215 */
	 tx1042 = tx914*z3; 		/* rule 216 */
	 tx1043 = tx918*z3; 		/* rule 217 */
	 tx1044 = x2*z2*z3; 		/* rule 218 */
	 tx1045 = power2(DN); 		/* rule 219 */
	 tx1046 = tx1023 + tx1024 + tx1025 + tx1026 + tx1027 + tx1028 + tx1029 + tx1030 + tx1031 + tx1032 + tx1033 + tx1034 + tx1035 + tx1036 + tx1037 + tx1038 + tx1039 + tx1040 + tx1041 + tx1042 + tx1043 + tx1044; 		/* rule 220 */
	 tx1047 = power2(tx978); 		/* rule 221 */
	 tx1048 = tx973 + tx974; 		/* rule 222 */
	 tx1049 = 2.*DeDPhi*tx1046*tx1047*tx927*tx971; 		/* rule 223 */
	 tx1050 = -(tx1045*tx1047*tx1048*tx930*tx967*V); 		/* rule 224 */
	 dhx1x1 = tx1049 + tx1050; 		/* rule 225 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx1051 = power2(x2); 		/* rule 227 */
	 tx1052 = power2(x3); 		/* rule 228 */
	 tx1053 = tx1021*tx956; 		/* rule 229 */
	 tx1054 = tx1051*tx956; 		/* rule 230 */
	 tx1055 = tx916*tx957; 		/* rule 231 */
	 tx1056 = tx909*tx979; 		/* rule 232 */
	 tx1057 = tx1022*tx980; 		/* rule 233 */
	 tx1058 = tx1052*tx980; 		/* rule 234 */
	 tx1059 = tx912*x1; 		/* rule 235 */
	 tx1060 = -2.*tx908*x2; 		/* rule 236 */
	 tx1061 = tx909*x3; 		/* rule 237 */
	 tx1062 = tx911*x3; 		/* rule 238 */
	 tx1063 = tx912*x3; 		/* rule 239 */
	 tx1064 = tx1021*y1; 		/* rule 240 */
	 tx1065 = tx1022*y1; 		/* rule 241 */
	 tx1066 = tx1051*y1; 		/* rule 242 */
	 tx1067 = tx1052*y1; 		/* rule 243 */
	 tx1068 = x2*x3*y2; 		/* rule 244 */
	 tx1069 = tx924*z1; 		/* rule 245 */
	 tx1070 = tx915*z2; 		/* rule 246 */
	 tx1071 = tx916*z2; 		/* rule 247 */
	 tx1072 = tx924*z2; 		/* rule 248 */
	 tx1073 = -2.*tx919*z3; 		/* rule 249 */
	 tx1074 = y3*z2*z3; 		/* rule 250 */
	 tx1075 = tx1053 + tx1054 + tx1055 + tx1056 + tx1057 + tx1058 + tx1059 + tx1060 + tx1061 + tx1062 + tx1063 + tx1064 + tx1065 + tx1066 + tx1067 + tx1068 + tx1069 + tx1070 + tx1071 + tx1072 + tx1073 + tx1074; 		/* rule 251 */
	 tx1076 = 2.*DeDPhi*tx1047*tx1075*tx926*tx971; 		/* rule 252 */
	 tx1077 = -(tx1045*tx1047*tx1048*tx929*tx967*V); 		/* rule 253 */
	 dhy1y1 = tx1076 + tx1077; 		/* rule 254 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx1078 = tx921*tx955; 		/* rule 256 */
	 tx1079 = tx1019*tx957; 		/* rule 257 */
	 tx1080 = tx1051*tx957; 		/* rule 258 */
	 tx1081 = tx919*tx980; 		/* rule 259 */
	 tx1082 = tx1020*tx981; 		/* rule 260 */
	 tx1083 = tx1052*tx981; 		/* rule 261 */
	 tx1084 = tx918*x1; 		/* rule 262 */
	 tx1085 = tx917*x2; 		/* rule 263 */
	 tx1086 = tx918*x2; 		/* rule 264 */
	 tx1087 = tx921*x2; 		/* rule 265 */
	 tx1088 = -2.*tx913*x3; 		/* rule 266 */
	 tx1089 = tx924*y1; 		/* rule 267 */
	 tx1090 = -2.*tx916*y2; 		/* rule 268 */
	 tx1091 = tx919*y3; 		/* rule 269 */
	 tx1092 = tx923*y3; 		/* rule 270 */
	 tx1093 = tx924*y3; 		/* rule 271 */
	 tx1094 = tx1019*z1; 		/* rule 272 */
	 tx1095 = tx1020*z1; 		/* rule 273 */
	 tx1096 = tx1051*z1; 		/* rule 274 */
	 tx1097 = tx1052*z1; 		/* rule 275 */
	 tx1098 = y2*y3*z2; 		/* rule 276 */
	 tx1099 = x2*x3*z3; 		/* rule 277 */
	 tx1100 = tx1078 + tx1079 + tx1080 + tx1081 + tx1082 + tx1083 + tx1084 + tx1085 + tx1086 + tx1087 + tx1088 + tx1089 + tx1090 + tx1091 + tx1092 + tx1093 + tx1094 + tx1095 + tx1096 + tx1097 + tx1098 + tx1099; 		/* rule 278 */
	 tx1101 = 2.*DeDPhi*tx1047*tx1100*tx925*tx971; 		/* rule 279 */
	 tx1102 = -(tx1045*tx1047*tx1048*tx928*tx967*V); 		/* rule 280 */
	 dhz1z1 = tx1101 + tx1102; 		/* rule 281 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 tx1103 = power2(y1); 		/* rule 283 */
	 tx1104 = power2(z1); 		/* rule 284 */
	 tx1105 = -tx1029; 		/* rule 285 */
	 tx1106 = -tx1030; 		/* rule 286 */
	 tx1107 = -tx1031; 		/* rule 287 */
	 tx1108 = -tx1032; 		/* rule 288 */
	 tx1109 = 2.*tx1035; 		/* rule 289 */
	 tx1110 = 2.*tx1039; 		/* rule 290 */
	 tx1111 = tx1103*tx955; 		/* rule 291 */
	 tx1112 = tx1104*tx955; 		/* rule 292 */
	 tx1113 = tx909*tx956; 		/* rule 293 */
	 tx1114 = tx913*tx957; 		/* rule 294 */
	 tx1115 = tx918*tx957; 		/* rule 295 */
	 tx1116 = tx908*tx980; 		/* rule 296 */
	 tx1117 = tx912*tx980; 		/* rule 297 */
	 tx1118 = tx1020*x2; 		/* rule 298 */
	 tx1119 = tx1022*x2; 		/* rule 299 */
	 tx1120 = tx1103*x2; 		/* rule 300 */
	 tx1121 = tx1104*x2; 		/* rule 301 */
	 tx1122 = tx1019*x3; 		/* rule 302 */
	 tx1123 = tx1021*x3; 		/* rule 303 */
	 tx1124 = -(tx909*y1); 		/* rule 304 */
	 tx1125 = x2*y1*y2; 		/* rule 305 */
	 tx1126 = tx907*y3; 		/* rule 306 */
	 tx1127 = tx908*y3; 		/* rule 307 */
	 tx1128 = 2.*tx909*y3; 		/* rule 308 */
	 tx1129 = tx910*y3; 		/* rule 309 */
	 tx1130 = x1*y1*y3; 		/* rule 310 */
	 tx1131 = tx917*z1; 		/* rule 311 */
	 tx1132 = tx921*z1; 		/* rule 312 */
	 tx1133 = tx913*z2; 		/* rule 313 */
	 tx1134 = tx914*z2; 		/* rule 314 */
	 tx1135 = 2.*tx921*z2; 		/* rule 315 */
	 tx1136 = tx922*z2; 		/* rule 316 */
	 tx1137 = tx917*z3; 		/* rule 317 */
	 tx1138 = x3*z1*z3; 		/* rule 318 */
	 tx1139 = tx1023 + tx1024 + tx1028 + tx1029 + tx1031 + tx1034 + tx1036 + tx1044 + tx1109 + tx1110 + tx1111 + tx1112 + tx1113 + tx1114 + tx1120 + tx1121 + tx1124 + tx1126 + tx1130 + tx1131 + tx1132 + tx1137; 		/* rule 319 */
	 tx1140 = tx1105 + tx1106 + tx1107 + tx1108 + tx1114 + tx1115 + tx1116 + tx1117 + tx1118 + tx1119 + tx1122 + tx1123 + tx1125 + tx1126 + tx1127 + tx1128 + tx1129 + tx1133 + tx1134 + tx1135 + tx1136 + tx1138; 		/* rule 320 */
	 tx1141 = power2(y4); 		/* rule 321 */
	 tx1142 = power2(z4); 		/* rule 322 */
	 tx1143 = -2.*tx1035; 		/* rule 323 */
	 tx1144 = -2.*tx1039; 		/* rule 324 */
	 tx1145 = -tx1132; 		/* rule 325 */
	 tx1146 = tx932*tx956; 		/* rule 326 */
	 tx1147 = tx935*tx956; 		/* rule 327 */
	 tx1148 = tx1139*tx927*tx967; 		/* rule 328 */
	 tx1149 = tx1103*tx979; 		/* rule 329 */
	 tx1150 = tx1104*tx979; 		/* rule 330 */
	 tx1151 = tx1141*tx979; 		/* rule 331 */
	 tx1152 = tx1142*tx979; 		/* rule 332 */
	 tx1153 = tx940*tx981; 		/* rule 333 */
	 tx1154 = tx1140*tx927*tx994; 		/* rule 334 */
	 tx1155 = tx1103*x3; 		/* rule 335 */
	 tx1156 = tx1104*x3; 		/* rule 336 */
	 tx1157 = tx1141*x3; 		/* rule 337 */
	 tx1158 = tx1142*x3; 		/* rule 338 */
	 tx1159 = tx1019*x4; 		/* rule 339 */
	 tx1160 = tx1020*x4; 		/* rule 340 */
	 tx1161 = tx1021*x4; 		/* rule 341 */
	 tx1162 = tx1022*x4; 		/* rule 342 */
	 tx1163 = tx909*y1; 		/* rule 343 */
	 tx1164 = tx911*y1; 		/* rule 344 */
	 tx1165 = tx934*y2; 		/* rule 345 */
	 tx1166 = tx935*y2; 		/* rule 346 */
	 tx1167 = tx909*y3; 		/* rule 347 */
	 tx1168 = -2.*tx932*y3; 		/* rule 348 */
	 tx1169 = tx910*y4; 		/* rule 349 */
	 tx1170 = tx912*y4; 		/* rule 350 */
	 tx1171 = 2.*tx1170; 		/* rule 351 */
	 tx1172 = tx932*y4; 		/* rule 352 */
	 tx1173 = tx933*y4; 		/* rule 353 */
	 tx1174 = tx921*z2; 		/* rule 354 */
	 tx1175 = -2.*tx938*z2; 		/* rule 355 */
	 tx1176 = x1*z1*z2; 		/* rule 356 */
	 tx1177 = tx936*z3; 		/* rule 357 */
	 tx1178 = tx940*z3; 		/* rule 358 */
	 tx1179 = 2.*tx1178; 		/* rule 359 */
	 tx1180 = tx941*z3; 		/* rule 360 */
	 tx1181 = -(tx918*z4); 		/* rule 361 */
	 tx1182 = tx918*z4; 		/* rule 362 */
	 tx1183 = -(tx938*z4); 		/* rule 363 */
	 tx1184 = x4*z2*z4; 		/* rule 364 */
	 tx1185 = tx1171; 		/* rule 365 */
	 tx1186 = tx1179; 		/* rule 366 */
	 tx1187 = tx1148 + tx1154; 		/* rule 367 */
	 tx1188 = tx1033 + tx1041 + tx1105 + tx1107 + tx1117 + tx1122 + tx1123 + tx1125 + tx1133 + tx1136 + tx1143 + tx1144 + tx1145 + tx1149 + tx1150 + tx1155 + tx1156 + tx1163 + tx1164 + tx1167 + tx1174 + tx1176; 		/* rule 368 */
	 tx1189 = tx1023 + tx1024 + tx1026 + tx1027 + tx1036 + tx1038 + tx1043 + tx1044 + tx1147 + tx1153 + tx1159 + tx1160 + tx1161 + tx1162 + tx1165 + tx1166 + tx1168 + tx1170 + tx1175 + tx1178 + tx1180 + tx1182; 		/* rule 369 */
	 tx1190 = tx1026 + tx1027 + tx1038 + tx1043 + tx1146 + tx1147 + tx1151 + tx1152 + tx1157 + tx1158 + tx1160 + tx1162 + tx1169 + tx1172 + tx1173 + tx1177 + tx1180 + tx1181 + tx1183 + tx1184 + tx1185 + tx1186; 		/* rule 370 */
	 tx1191 = power2(tx991); 		/* rule 371 */
	 tx1813 = tx967; 		/* rule 372 */
	 tx1814 = reciprocal(tx1813); 		/* rule 373 */
	 tx1815 = tx992; 		/* rule 374 */
	 tx1192 = tx1814*tx1815; 		/* rule 375 */
	 tx1193 = 2.*tx1046*tx1047*tx927*tx971; 		/* rule 376 */
	 tx1194 = tx1046*tx1192*tx927*tx978; 		/* rule 377 */
	 tx1195 = -(tx1189*tx1192*tx946*tx991); 		/* rule 378 */
	 tx1196 = -2.*tx1047*tx1187*tx1815; 		/* rule 379 */
	 tx1197 = -2.*tx1190*tx1191*tx1815*tx946*tx993; 		/* rule 380 */
	 tx1198 = 2.*tx1047*tx1188*tx1815*tx927*tx994; 		/* rule 381 */
	 tx1199 = tx1193 + tx1194 + tx1195 + tx1196 + tx1197 + tx1198; 		/* rule 382 */
	 tx1200 = power2(tx998); 		/* rule 383 */
	 tx1201 = DeDPhi*tx1199; 		/* rule 384 */
	 tx1202 = -(tx1045*tx1048*tx1200*V); 		/* rule 385 */
	 dhx2x2 = tx1201 + tx1202; 		/* rule 386 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 tx1203 = power2(x1); 		/* rule 388 */
	 tx1204 = 2.*tx1059; 		/* rule 389 */
	 tx1205 = -tx1064; 		/* rule 390 */
	 tx1206 = -tx1065; 		/* rule 391 */
	 tx1207 = -tx1066; 		/* rule 392 */
	 tx1208 = -tx1067; 		/* rule 393 */
	 tx1209 = 2.*tx1071; 		/* rule 394 */
	 tx1210 = tx909*tx955; 		/* rule 395 */
	 tx1211 = tx912*tx955; 		/* rule 396 */
	 tx1212 = tx1104*tx956; 		/* rule 397 */
	 tx1213 = tx1203*tx956; 		/* rule 398 */
	 tx1214 = tx919*tx957; 		/* rule 399 */
	 tx1215 = tx916*tx981; 		/* rule 400 */
	 tx1216 = tx924*tx981; 		/* rule 401 */
	 tx1217 = tx907*x1; 		/* rule 402 */
	 tx1218 = tx908*x1; 		/* rule 403 */
	 tx1219 = 2.*tx908*x2; 		/* rule 404 */
	 tx1220 = tx909*x2; 		/* rule 405 */
	 tx1221 = tx910*x2; 		/* rule 406 */
	 tx1222 = tx911*x2; 		/* rule 407 */
	 tx1223 = tx907*x3; 		/* rule 408 */
	 tx1224 = tx1022*y2; 		/* rule 409 */
	 tx1225 = tx1052*y2; 		/* rule 410 */
	 tx1226 = tx1104*y2; 		/* rule 411 */
	 tx1227 = tx1203*y2; 		/* rule 412 */
	 tx1228 = tx1021*y3; 		/* rule 413 */
	 tx1229 = tx1051*y3; 		/* rule 414 */
	 tx1230 = x1*x3*y3; 		/* rule 415 */
	 tx1231 = -(tx919*z1); 		/* rule 416 */
	 tx1232 = y2*z1*z2; 		/* rule 417 */
	 tx1233 = tx915*z3; 		/* rule 418 */
	 tx1234 = tx916*z3; 		/* rule 419 */
	 tx1235 = 2.*tx919*z3; 		/* rule 420 */
	 tx1236 = tx920*z3; 		/* rule 421 */
	 tx1237 = y1*z1*z3; 		/* rule 422 */
	 tx1238 = tx1205 + tx1206 + tx1207 + tx1208 + tx1210 + tx1211 + tx1215 + tx1216 + tx1219 + tx1220 + tx1221 + tx1222 + tx1224 + tx1225 + tx1228 + tx1229 + tx1230 + tx1232 + tx1233 + tx1234 + tx1235 + tx1236; 		/* rule 423 */
	 tx1239 = tx1053 + tx1054 + tx1056 + tx1064 + tx1066 + tx1068 + tx1070 + tx1072 + tx1204 + tx1209 + tx1210 + tx1212 + tx1213 + tx1214 + tx1217 + tx1218 + tx1223 + tx1226 + tx1227 + tx1231 + tx1233 + tx1237; 		/* rule 424 */
	 tx1240 = power2(x4); 		/* rule 425 */
	 tx1241 = -2.*tx1059; 		/* rule 426 */
	 tx1242 = -2.*tx1071; 		/* rule 427 */
	 tx1243 = -tx1218; 		/* rule 428 */
	 tx1244 = tx937*tx957; 		/* rule 429 */
	 tx1245 = tx943*tx957; 		/* rule 430 */
	 tx1246 = tx1239*tx1813*tx926; 		/* rule 431 */
	 tx1247 = tx932*tx979; 		/* rule 432 */
	 tx1248 = tx1104*tx980; 		/* rule 433 */
	 tx1249 = tx1142*tx980; 		/* rule 434 */
	 tx1250 = tx1203*tx980; 		/* rule 435 */
	 tx1251 = tx1240*tx980; 		/* rule 436 */
	 tx1252 = tx1238*tx926*tx994; 		/* rule 437 */
	 tx1253 = tx908*x2; 		/* rule 438 */
	 tx1254 = -2.*tx935*x2; 		/* rule 439 */
	 tx1255 = tx932*x3; 		/* rule 440 */
	 tx1256 = 2.*tx1255; 		/* rule 441 */
	 tx1257 = tx933*x3; 		/* rule 442 */
	 tx1258 = tx934*x3; 		/* rule 443 */
	 tx1259 = -(tx912*x4); 		/* rule 444 */
	 tx1260 = tx912*x4; 		/* rule 445 */
	 tx1261 = -(tx935*x4); 		/* rule 446 */
	 tx1262 = x1*x2*y1; 		/* rule 447 */
	 tx1263 = tx1104*y3; 		/* rule 448 */
	 tx1264 = tx1142*y3; 		/* rule 449 */
	 tx1265 = tx1203*y3; 		/* rule 450 */
	 tx1266 = tx1240*y3; 		/* rule 451 */
	 tx1267 = tx1021*y4; 		/* rule 452 */
	 tx1268 = tx1022*y4; 		/* rule 453 */
	 tx1269 = tx1051*y4; 		/* rule 454 */
	 tx1270 = tx1052*y4; 		/* rule 455 */
	 tx1271 = x2*x4*y4; 		/* rule 456 */
	 tx1272 = tx919*z1; 		/* rule 457 */
	 tx1273 = tx923*z1; 		/* rule 458 */
	 tx1274 = tx942*z2; 		/* rule 459 */
	 tx1275 = tx943*z2; 		/* rule 460 */
	 tx1276 = tx919*z3; 		/* rule 461 */
	 tx1277 = -2.*tx937*z3; 		/* rule 462 */
	 tx1278 = tx920*z4; 		/* rule 463 */
	 tx1279 = tx924*z4; 		/* rule 464 */
	 tx1280 = 2.*tx1279; 		/* rule 465 */
	 tx1281 = tx937*z4; 		/* rule 466 */
	 tx1282 = tx939*z4; 		/* rule 467 */
	 tx1283 = tx1256; 		/* rule 468 */
	 tx1284 = tx1280; 		/* rule 469 */
	 tx1285 = tx1246 + tx1252; 		/* rule 470 */
	 tx1286 = tx1061 + tx1069 + tx1205 + tx1207 + tx1216 + tx1220 + tx1221 + tx1228 + tx1229 + tx1232 + tx1241 + tx1242 + tx1243 + tx1248 + tx1250 + tx1253 + tx1262 + tx1263 + tx1265 + tx1272 + tx1273 + tx1276; 		/* rule 471 */
	 tx1287 = tx1053 + tx1054 + tx1057 + tx1058 + tx1063 + tx1068 + tx1072 + tx1074 + tx1245 + tx1247 + tx1254 + tx1255 + tx1257 + tx1260 + tx1267 + tx1268 + tx1269 + tx1270 + tx1274 + tx1275 + tx1277 + tx1279; 		/* rule 472 */
	 tx1288 = tx1057 + tx1058 + tx1063 + tx1074 + tx1244 + tx1245 + tx1249 + tx1251 + tx1257 + tx1258 + tx1259 + tx1261 + tx1264 + tx1266 + tx1268 + tx1270 + tx1271 + tx1278 + tx1281 + tx1282 + tx1283 + tx1284; 		/* rule 473 */
	 tx1289 = 2.*tx1047*tx1075*tx926*tx971; 		/* rule 474 */
	 tx1290 = tx1075*tx1192*tx926*tx978; 		/* rule 475 */
	 tx1291 = -(tx1192*tx1287*tx945*tx991); 		/* rule 476 */
	 tx1292 = -2.*tx1047*tx1285*tx1815; 		/* rule 477 */
	 tx1293 = -2.*tx1191*tx1288*tx1815*tx945*tx993; 		/* rule 478 */
	 tx1294 = 2.*tx1047*tx1286*tx1815*tx926*tx994; 		/* rule 479 */
	 tx1295 = power2(tx1002); 		/* rule 480 */
	 tx1296 = tx1289 + tx1290 + tx1291 + tx1292 + tx1293 + tx1294; 		/* rule 481 */
	 tx1297 = DeDPhi*tx1296; 		/* rule 482 */
	 tx1298 = -(tx1045*tx1048*tx1295*V); 		/* rule 483 */
	 dhy2y2 = tx1297 + tx1298; 		/* rule 484 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 tx1299 = 2.*tx1087; 		/* rule 486 */
	 tx1300 = 2.*tx1089; 		/* rule 487 */
	 tx1301 = -tx1094; 		/* rule 488 */
	 tx1302 = -tx1095; 		/* rule 489 */
	 tx1303 = -tx1096; 		/* rule 490 */
	 tx1304 = -tx1097; 		/* rule 491 */
	 tx1305 = tx913*tx955; 		/* rule 492 */
	 tx1306 = tx919*tx956; 		/* rule 493 */
	 tx1307 = tx924*tx956; 		/* rule 494 */
	 tx1308 = tx1103*tx957; 		/* rule 495 */
	 tx1309 = tx1203*tx957; 		/* rule 496 */
	 tx1310 = tx918*tx979; 		/* rule 497 */
	 tx1311 = tx921*tx979; 		/* rule 498 */
	 tx1312 = -(tx913*x1); 		/* rule 499 */
	 tx1313 = 2.*tx913*x3; 		/* rule 500 */
	 tx1314 = tx917*x3; 		/* rule 501 */
	 tx1315 = tx921*x3; 		/* rule 502 */
	 tx1316 = tx922*x3; 		/* rule 503 */
	 tx1317 = tx915*y1; 		/* rule 504 */
	 tx1318 = tx916*y1; 		/* rule 505 */
	 tx1319 = 2.*tx916*y2; 		/* rule 506 */
	 tx1320 = tx919*y2; 		/* rule 507 */
	 tx1321 = tx920*y2; 		/* rule 508 */
	 tx1322 = tx923*y2; 		/* rule 509 */
	 tx1323 = tx915*y3; 		/* rule 510 */
	 tx1324 = x1*x3*z1; 		/* rule 511 */
	 tx1325 = tx1020*z2; 		/* rule 512 */
	 tx1326 = tx1052*z2; 		/* rule 513 */
	 tx1327 = tx1103*z2; 		/* rule 514 */
	 tx1328 = tx1203*z2; 		/* rule 515 */
	 tx1329 = x1*x2*z2; 		/* rule 516 */
	 tx1330 = tx1019*z3; 		/* rule 517 */
	 tx1331 = tx1051*z3; 		/* rule 518 */
	 tx1332 = y1*y3*z3; 		/* rule 519 */
	 tx1333 = tx1079 + tx1080 + tx1081 + tx1085 + tx1086 + tx1094 + tx1096 + tx1098 + tx1299 + tx1300 + tx1305 + tx1306 + tx1308 + tx1309 + tx1312 + tx1314 + tx1317 + tx1318 + tx1323 + tx1324 + tx1327 + tx1328; 		/* rule 520 */
	 tx1334 = tx1301 + tx1302 + tx1303 + tx1304 + tx1306 + tx1307 + tx1310 + tx1311 + tx1313 + tx1314 + tx1315 + tx1316 + tx1319 + tx1320 + tx1321 + tx1322 + tx1325 + tx1326 + tx1329 + tx1330 + tx1331 + tx1332; 		/* rule 521 */
	 tx1335 = -2.*tx1087; 		/* rule 522 */
	 tx1336 = -2.*tx1089; 		/* rule 523 */
	 tx1337 = -tx1318; 		/* rule 524 */
	 tx1338 = tx938*tx955; 		/* rule 525 */
	 tx1339 = tx940*tx955; 		/* rule 526 */
	 tx1340 = tx1333*tx1813*tx925; 		/* rule 527 */
	 tx1341 = tx937*tx980; 		/* rule 528 */
	 tx1342 = tx1103*tx981; 		/* rule 529 */
	 tx1343 = tx1141*tx981; 		/* rule 530 */
	 tx1344 = tx1203*tx981; 		/* rule 531 */
	 tx1345 = tx1240*tx981; 		/* rule 532 */
	 tx1346 = tx1334*tx925*tx994; 		/* rule 533 */
	 tx1347 = tx913*x1; 		/* rule 534 */
	 tx1348 = tx914*x1; 		/* rule 535 */
	 tx1349 = tx936*x2; 		/* rule 536 */
	 tx1350 = tx938*x2; 		/* rule 537 */
	 tx1351 = tx913*x3; 		/* rule 538 */
	 tx1352 = -2.*tx940*x3; 		/* rule 539 */
	 tx1353 = tx918*x4; 		/* rule 540 */
	 tx1354 = 2.*tx1353; 		/* rule 541 */
	 tx1355 = tx922*x4; 		/* rule 542 */
	 tx1356 = tx940*x4; 		/* rule 543 */
	 tx1357 = tx941*x4; 		/* rule 544 */
	 tx1358 = tx916*y2; 		/* rule 545 */
	 tx1359 = -2.*tx943*y2; 		/* rule 546 */
	 tx1360 = tx937*y3; 		/* rule 547 */
	 tx1361 = 2.*tx1360; 		/* rule 548 */
	 tx1362 = tx939*y3; 		/* rule 549 */
	 tx1363 = tx942*y3; 		/* rule 550 */
	 tx1364 = -(tx924*y4); 		/* rule 551 */
	 tx1365 = tx924*y4; 		/* rule 552 */
	 tx1366 = -(tx943*y4); 		/* rule 553 */
	 tx1367 = y1*y2*z1; 		/* rule 554 */
	 tx1368 = tx1103*z3; 		/* rule 555 */
	 tx1369 = tx1141*z3; 		/* rule 556 */
	 tx1370 = tx1203*z3; 		/* rule 557 */
	 tx1371 = tx1240*z3; 		/* rule 558 */
	 tx1372 = tx1019*z4; 		/* rule 559 */
	 tx1373 = tx1020*z4; 		/* rule 560 */
	 tx1374 = tx1051*z4; 		/* rule 561 */
	 tx1375 = tx1052*z4; 		/* rule 562 */
	 tx1376 = y2*y4*z4; 		/* rule 563 */
	 tx1377 = tx1354; 		/* rule 564 */
	 tx1378 = tx1361; 		/* rule 565 */
	 tx1379 = tx1340 + tx1346; 		/* rule 566 */
	 tx1380 = tx1084 + tx1091 + tx1301 + tx1303 + tx1310 + tx1320 + tx1321 + tx1329 + tx1330 + tx1331 + tx1335 + tx1336 + tx1337 + tx1342 + tx1344 + tx1347 + tx1348 + tx1351 + tx1358 + tx1367 + tx1368 + tx1370; 		/* rule 567 */
	 tx1381 = tx1079 + tx1080 + tx1082 + tx1083 + tx1086 + tx1093 + tx1098 + tx1099 + tx1338 + tx1341 + tx1349 + tx1350 + tx1352 + tx1353 + tx1359 + tx1360 + tx1362 + tx1365 + tx1372 + tx1373 + tx1374 + tx1375; 		/* rule 568 */
	 tx1382 = tx1082 + tx1083 + tx1093 + tx1099 + tx1338 + tx1339 + tx1343 + tx1345 + tx1355 + tx1356 + tx1357 + tx1362 + tx1363 + tx1364 + tx1366 + tx1369 + tx1371 + tx1373 + tx1375 + tx1376 + tx1377 + tx1378; 		/* rule 569 */
	 tx1383 = 2.*tx1047*tx1100*tx925*tx971; 		/* rule 570 */
	 tx1384 = tx1100*tx1192*tx925*tx978; 		/* rule 571 */
	 tx1385 = -(tx1192*tx1381*tx944*tx991); 		/* rule 572 */
	 tx1386 = -2.*tx1047*tx1379*tx1815; 		/* rule 573 */
	 tx1387 = -2.*tx1191*tx1382*tx1815*tx944*tx993; 		/* rule 574 */
	 tx1388 = 2.*tx1047*tx1380*tx1815*tx925*tx994; 		/* rule 575 */
	 tx1389 = power2(tx1006); 		/* rule 576 */
	 tx1390 = tx1383 + tx1384 + tx1385 + tx1386 + tx1387 + tx1388; 		/* rule 577 */
	 tx1391 = DeDPhi*tx1390; 		/* rule 578 */
	 tx1392 = -(tx1045*tx1048*tx1389*V); 		/* rule 579 */
	 dhz2z2 = tx1391 + tx1392; 		/* rule 580 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	 tx1393 = -tx1159; 		/* rule 582 */
	 tx1394 = -tx1160; 		/* rule 583 */
	 tx1395 = -tx1161; 		/* rule 584 */
	 tx1396 = -tx1162; 		/* rule 585 */
	 tx1397 = -2.*tx1170; 		/* rule 586 */
	 tx1398 = -tx1170; 		/* rule 587 */
	 tx1399 = -tx1172; 		/* rule 588 */
	 tx1400 = -2.*tx1178; 		/* rule 589 */
	 tx1401 = tx1141*tx955; 		/* rule 590 */
	 tx1402 = tx1142*tx955; 		/* rule 591 */
	 tx1403 = tx1141*x2; 		/* rule 592 */
	 tx1404 = tx1142*x2; 		/* rule 593 */
	 tx1405 = tx932*y3; 		/* rule 594 */
	 tx1406 = 2.*tx1405; 		/* rule 595 */
	 tx1407 = tx935*y3; 		/* rule 596 */
	 tx1408 = x2*y2*y4; 		/* rule 597 */
	 tx1409 = x4*y3*y4; 		/* rule 598 */
	 tx1410 = tx938*z2; 		/* rule 599 */
	 tx1411 = 2.*tx1410; 		/* rule 600 */
	 tx1412 = tx940*z2; 		/* rule 601 */
	 tx1413 = tx922*z4; 		/* rule 602 */
	 tx1414 = tx936*z4; 		/* rule 603 */
	 tx1415 = tx938*z4; 		/* rule 604 */
	 tx1416 = x3*z3*z4; 		/* rule 605 */
	 tx1417 = tx1406; 		/* rule 606 */
	 tx1418 = tx1411; 		/* rule 607 */
	 tx1419 = tx1115 + tx1118 + tx1119 + tx1129 + tx1166 + tx1182 + tx1394 + tx1396 + tx1397 + tx1399 + tx1400 + tx1401 + tx1402 + tx1403 + tx1404 + tx1405 + tx1407 + tx1409 + tx1410 + tx1414 + tx1415 + tx1416; 		/* rule 608 */
	 tx1420 = tx1115 + tx1117 + tx1118 + tx1119 + tx1122 + tx1123 + tx1129 + tx1136 + tx1169 + tx1181 + tx1393 + tx1394 + tx1395 + tx1396 + tx1398 + tx1407 + tx1408 + tx1412 + tx1413 + tx1416 + tx1417 + tx1418; 		/* rule 609 */
	 tx1421 = tx1419*tx1813*tx946; 		/* rule 610 */
	 tx1422 = tx1420*tx946*tx993; 		/* rule 611 */
	 tx1423 = tx1421 + tx1422; 		/* rule 612 */
	 tx1424 = -2.*tx1189*tx1191*tx946*tx971; 		/* rule 613 */
	 tx1425 = -2.*tx1191*tx1423*tx1815; 		/* rule 614 */
	 tx1426 = power2(tx1010); 		/* rule 615 */
	 tx1427 = tx1194 + tx1195 + tx1197 + tx1198 + tx1424 + tx1425; 		/* rule 616 */
	 tx1428 = DeDPhi*tx1427; 		/* rule 617 */
	 tx1429 = -(tx1045*tx1048*tx1426*V); 		/* rule 618 */
	 dhx3x3 = tx1428 + tx1429; 		/* rule 619 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 0, dhx3x3);
	 tx1430 = -2.*tx1255; 		/* rule 621 */
	 tx1431 = -tx1267; 		/* rule 622 */
	 tx1432 = -tx1268; 		/* rule 623 */
	 tx1433 = -tx1269; 		/* rule 624 */
	 tx1434 = -tx1270; 		/* rule 625 */
	 tx1435 = -2.*tx1279; 		/* rule 626 */
	 tx1436 = -tx1279; 		/* rule 627 */
	 tx1437 = -tx1281; 		/* rule 628 */
	 tx1438 = tx1142*tx956; 		/* rule 629 */
	 tx1439 = tx1240*tx956; 		/* rule 630 */
	 tx1440 = tx932*x2; 		/* rule 631 */
	 tx1441 = tx935*x2; 		/* rule 632 */
	 tx1442 = 2.*tx1441; 		/* rule 633 */
	 tx1443 = tx910*x4; 		/* rule 634 */
	 tx1444 = tx934*x4; 		/* rule 635 */
	 tx1445 = tx935*x4; 		/* rule 636 */
	 tx1446 = tx1142*y2; 		/* rule 637 */
	 tx1447 = tx1240*y2; 		/* rule 638 */
	 tx1448 = x3*x4*y3; 		/* rule 639 */
	 tx1449 = tx937*z3; 		/* rule 640 */
	 tx1450 = 2.*tx1449; 		/* rule 641 */
	 tx1451 = tx943*z3; 		/* rule 642 */
	 tx1452 = y2*z2*z4; 		/* rule 643 */
	 tx1453 = y4*z3*z4; 		/* rule 644 */
	 tx1454 = tx1442; 		/* rule 645 */
	 tx1455 = tx1450; 		/* rule 646 */
	 tx1456 = tx1211 + tx1224 + tx1225 + tx1236 + tx1260 + tx1275 + tx1430 + tx1432 + tx1434 + tx1435 + tx1437 + tx1438 + tx1439 + tx1441 + tx1444 + tx1445 + tx1446 + tx1447 + tx1448 + tx1449 + tx1451 + tx1453; 		/* rule 647 */
	 tx1457 = tx1211 + tx1216 + tx1221 + tx1224 + tx1225 + tx1228 + tx1229 + tx1236 + tx1259 + tx1278 + tx1431 + tx1432 + tx1433 + tx1434 + tx1436 + tx1440 + tx1443 + tx1448 + tx1451 + tx1452 + tx1454 + tx1455; 		/* rule 648 */
	 tx1458 = tx1456*tx1813*tx945; 		/* rule 649 */
	 tx1459 = tx1457*tx945*tx993; 		/* rule 650 */
	 tx1460 = tx1458 + tx1459; 		/* rule 651 */
	 tx1461 = -2.*tx1191*tx1287*tx945*tx971; 		/* rule 652 */
	 tx1462 = -2.*tx1191*tx1460*tx1815; 		/* rule 653 */
	 tx1463 = power2(tx1014); 		/* rule 654 */
	 tx1464 = tx1290 + tx1291 + tx1293 + tx1294 + tx1461 + tx1462; 		/* rule 655 */
	 tx1465 = DeDPhi*tx1464; 		/* rule 656 */
	 tx1466 = -(tx1045*tx1048*tx1463*V); 		/* rule 657 */
	 dhy3y3 = tx1465 + tx1466; 		/* rule 658 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 1, dhy3y3);
	 tx1467 = -2.*tx1353; 		/* rule 660 */
	 tx1468 = -tx1353; 		/* rule 661 */
	 tx1469 = -tx1356; 		/* rule 662 */
	 tx1470 = -2.*tx1360; 		/* rule 663 */
	 tx1471 = -tx1372; 		/* rule 664 */
	 tx1472 = -tx1373; 		/* rule 665 */
	 tx1473 = -tx1374; 		/* rule 666 */
	 tx1474 = -tx1375; 		/* rule 667 */
	 tx1475 = tx1141*tx957; 		/* rule 668 */
	 tx1476 = tx1240*tx957; 		/* rule 669 */
	 tx1477 = tx938*x3; 		/* rule 670 */
	 tx1478 = tx940*x3; 		/* rule 671 */
	 tx1479 = 2.*tx1478; 		/* rule 672 */
	 tx1480 = tx937*y2; 		/* rule 673 */
	 tx1481 = tx943*y2; 		/* rule 674 */
	 tx1482 = 2.*tx1481; 		/* rule 675 */
	 tx1483 = tx920*y4; 		/* rule 676 */
	 tx1484 = tx942*y4; 		/* rule 677 */
	 tx1485 = tx943*y4; 		/* rule 678 */
	 tx1486 = tx1141*z2; 		/* rule 679 */
	 tx1487 = tx1240*z2; 		/* rule 680 */
	 tx1488 = x2*x4*z2; 		/* rule 681 */
	 tx1489 = y3*y4*z3; 		/* rule 682 */
	 tx1490 = x3*x4*z4; 		/* rule 683 */
	 tx1491 = tx1479; 		/* rule 684 */
	 tx1492 = tx1482; 		/* rule 685 */
	 tx1493 = tx1307 + tx1316 + tx1325 + tx1326 + tx1350 + tx1365 + tx1467 + tx1469 + tx1470 + tx1472 + tx1474 + tx1475 + tx1476 + tx1477 + tx1478 + tx1481 + tx1484 + tx1485 + tx1486 + tx1487 + tx1489 + tx1490; 		/* rule 686 */
	 tx1494 = tx1307 + tx1310 + tx1316 + tx1321 + tx1325 + tx1326 + tx1330 + tx1331 + tx1355 + tx1364 + tx1468 + tx1471 + tx1472 + tx1473 + tx1474 + tx1477 + tx1480 + tx1483 + tx1488 + tx1489 + tx1491 + tx1492; 		/* rule 687 */
	 tx1495 = tx1493*tx1813*tx944; 		/* rule 688 */
	 tx1496 = tx1494*tx944*tx993; 		/* rule 689 */
	 tx1497 = tx1495 + tx1496; 		/* rule 690 */
	 tx1498 = -2.*tx1191*tx1381*tx944*tx971; 		/* rule 691 */
	 tx1499 = -2.*tx1191*tx1497*tx1815; 		/* rule 692 */
	 tx1500 = power2(tx1018); 		/* rule 693 */
	 tx1501 = tx1384 + tx1385 + tx1387 + tx1388 + tx1498 + tx1499; 		/* rule 694 */
	 tx1502 = DeDPhi*tx1501; 		/* rule 695 */
	 tx1503 = -(tx1045*tx1048*tx1500*V); 		/* rule 696 */
	 dhz3z3 = tx1502 + tx1503; 		/* rule 697 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I3, 2, dhz3z3);
	 tx1504 = DeDPhi*tx1424; 		/* rule 699 */
	 tx1505 = -(tx1045*tx1048*tx1191*tx1813*tx949*V); 		/* rule 700 */
	 dhx4x4 = tx1504 + tx1505; 		/* rule 701 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 0, dhx4x4);
	 tx1506 = DeDPhi*tx1461; 		/* rule 703 */
	 tx1507 = -(tx1045*tx1048*tx1191*tx1813*tx948*V); 		/* rule 704 */
	 dhy4y4 = tx1506 + tx1507; 		/* rule 705 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I4, 1, I4, 1, dhy4y4);
	 tx1508 = DeDPhi*tx1498; 		/* rule 707 */
	 tx1509 = -(tx1045*tx1048*tx1191*tx1813*tx947*V); 		/* rule 708 */
	 dhz4z4 = tx1508 + tx1509; 		/* rule 709 */
	DIHEDRAL_DIAGONAL_HESSIAN_ACCUMULATE(I4, 2, I4, 2, dhz4z4);
	#ifdef DIHEDRAL_CALC_OFF_DIAGONAL_HESSIAN //[
	if (calcOffDiagonalHessian) { 
	 tx1510 = tx1046*tx926; 		/* rule 713 */
	 tx1511 = tx1075*tx927; 		/* rule 714 */
	 tx1512 = tx1510 + tx1511; 		/* rule 715 */
	 tx1513 = DeDPhi*tx1047*tx1512*tx971; 		/* rule 716 */
	 tx1514 = -(tx1045*tx1047*tx1048*tx1813*tx926*tx927*V); 		/* rule 717 */
	 ohx1y1 = tx1513 + tx1514; 		/* rule 718 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tx1515 = tx1046*tx925; 		/* rule 720 */
	 tx1516 = tx1100*tx927; 		/* rule 721 */
	 tx1517 = tx1515 + tx1516; 		/* rule 722 */
	 tx1518 = DeDPhi*tx1047*tx1517*tx971; 		/* rule 723 */
	 tx1519 = -(tx1045*tx1047*tx1048*tx1813*tx925*tx927*V); 		/* rule 724 */
	 ohx1z1 = tx1518 + tx1519; 		/* rule 725 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx1520 = -2.*tx1046*tx1047*tx927*tx971; 		/* rule 727 */
	 tx1521 = tx1047*tx1187*tx1815; 		/* rule 728 */
	 tx1522 = tx1520 + tx1521; 		/* rule 729 */
	 tx1523 = DeDPhi*tx1522; 		/* rule 730 */
	 tx1524 = tx1045*tx1048*tx995*tx998*V; 		/* rule 731 */
	 ohx1x2 = tx1523 + tx1524; 		/* rule 732 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx1525 = tx1139*tx1813*tx926; 		/* rule 734 */
	 tx1526 = tx1238*tx927*tx994; 		/* rule 735 */
	 tx1527 = tx1525 + tx1526; 		/* rule 736 */
	 tx1528 = -(tx1047*tx1512*tx971); 		/* rule 737 */
	 tx1529 = tx1047*tx1527*tx1815; 		/* rule 738 */
	 tx1530 = tx1528 + tx1529; 		/* rule 739 */
	 tx1531 = DeDPhi*tx1530; 		/* rule 740 */
	 tx1532 = tx1002*tx1045*tx1048*tx995*V; 		/* rule 741 */
	 ohx1y2 = tx1531 + tx1532; 		/* rule 742 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx1533 = tx1139*tx1813*tx925; 		/* rule 744 */
	 tx1534 = tx1334*tx927*tx994; 		/* rule 745 */
	 tx1535 = tx1533 + tx1534; 		/* rule 746 */
	 tx1536 = -(tx1047*tx1517*tx971); 		/* rule 747 */
	 tx1537 = tx1047*tx1535*tx1815; 		/* rule 748 */
	 tx1538 = tx1536 + tx1537; 		/* rule 749 */
	 tx1539 = DeDPhi*tx1538; 		/* rule 750 */
	 tx1540 = tx1006*tx1045*tx1048*tx995*V; 		/* rule 751 */
	 ohx1z2 = tx1539 + tx1540; 		/* rule 752 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tx1541 = -(DeDPhi*tx1521); 		/* rule 754 */
	 tx1542 = tx1010*tx1045*tx1048*tx995*V; 		/* rule 755 */
	 ohx1x3 = tx1541 + tx1542; 		/* rule 756 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 0, ohx1x3);
	 tx1543 = -(DeDPhi*tx1529); 		/* rule 758 */
	 tx1544 = tx1014*tx1045*tx1048*tx995*V; 		/* rule 759 */
	 ohx1y3 = tx1543 + tx1544; 		/* rule 760 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 1, ohx1y3);
	 tx1545 = -(DeDPhi*tx1537); 		/* rule 762 */
	 tx1546 = tx1018*tx1045*tx1048*tx995*V; 		/* rule 763 */
	 ohx1z3 = tx1545 + tx1546; 		/* rule 764 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 2, ohx1z3);
	 ohx1x4 = tx1045*tx1048*tx1813*tx953*tx978*tx991*V; 		/* rule 766 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 0, ohx1x4);
	 ohx1y4 = tx1045*tx1048*tx1813*tx927*tx945*tx978*tx991*V; 		/* rule 768 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 1, ohx1y4);
	 ohx1z4 = tx1045*tx1048*tx1813*tx927*tx944*tx978*tx991*V; 		/* rule 770 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 2, ohx1z4);
	 tx1547 = tx1075*tx925; 		/* rule 772 */
	 tx1548 = tx1100*tx926; 		/* rule 773 */
	 tx1549 = tx1547 + tx1548; 		/* rule 774 */
	 tx1550 = DeDPhi*tx1047*tx1549*tx971; 		/* rule 775 */
	 tx1551 = -(tx1045*tx1047*tx1048*tx1813*tx925*tx926*V); 		/* rule 776 */
	 ohy1z1 = tx1550 + tx1551; 		/* rule 777 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 tx1552 = tx1239*tx1813*tx927; 		/* rule 779 */
	 tx1553 = tx1140*tx926*tx994; 		/* rule 780 */
	 tx1554 = tx1552 + tx1553; 		/* rule 781 */
	 tx1555 = tx1047*tx1554*tx1815; 		/* rule 782 */
	 tx1556 = tx1528 + tx1555; 		/* rule 783 */
	 tx1557 = DeDPhi*tx1556; 		/* rule 784 */
	 tx1558 = tx1045*tx1048*tx998*tx999*V; 		/* rule 785 */
	 ohy1x2 = tx1557 + tx1558; 		/* rule 786 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx1559 = -2.*tx1047*tx1075*tx926*tx971; 		/* rule 788 */
	 tx1560 = tx1047*tx1285*tx1815; 		/* rule 789 */
	 tx1561 = tx1559 + tx1560; 		/* rule 790 */
	 tx1562 = DeDPhi*tx1561; 		/* rule 791 */
	 tx1563 = tx1002*tx1045*tx1048*tx999*V; 		/* rule 792 */
	 ohy1y2 = tx1562 + tx1563; 		/* rule 793 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx1564 = tx1239*tx1813*tx925; 		/* rule 795 */
	 tx1565 = tx1334*tx926*tx994; 		/* rule 796 */
	 tx1566 = tx1564 + tx1565; 		/* rule 797 */
	 tx1567 = -(tx1047*tx1549*tx971); 		/* rule 798 */
	 tx1568 = tx1047*tx1566*tx1815; 		/* rule 799 */
	 tx1569 = tx1567 + tx1568; 		/* rule 800 */
	 tx1570 = DeDPhi*tx1569; 		/* rule 801 */
	 tx1571 = tx1006*tx1045*tx1048*tx999*V; 		/* rule 802 */
	 ohy1z2 = tx1570 + tx1571; 		/* rule 803 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 tx1572 = -(DeDPhi*tx1555); 		/* rule 805 */
	 tx1573 = tx1010*tx1045*tx1048*tx999*V; 		/* rule 806 */
	 ohy1x3 = tx1572 + tx1573; 		/* rule 807 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 0, ohy1x3);
	 tx1574 = -(DeDPhi*tx1560); 		/* rule 809 */
	 tx1575 = tx1014*tx1045*tx1048*tx999*V; 		/* rule 810 */
	 ohy1y3 = tx1574 + tx1575; 		/* rule 811 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 1, ohy1y3);
	 tx1576 = -(DeDPhi*tx1568); 		/* rule 813 */
	 tx1577 = tx1018*tx1045*tx1048*tx999*V; 		/* rule 814 */
	 ohy1z3 = tx1576 + tx1577; 		/* rule 815 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 2, ohy1z3);
	 ohy1x4 = tx1045*tx1048*tx1813*tx926*tx946*tx978*tx991*V; 		/* rule 817 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 0, ohy1x4);
	 ohy1y4 = tx1045*tx1048*tx1813*tx952*tx978*tx991*V; 		/* rule 819 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 1, ohy1y4);
	 ohy1z4 = tx1045*tx1048*tx1813*tx926*tx944*tx978*tx991*V; 		/* rule 821 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 2, ohy1z4);
	 tx1578 = tx1333*tx1813*tx927; 		/* rule 823 */
	 tx1579 = tx1140*tx925*tx994; 		/* rule 824 */
	 tx1580 = tx1578 + tx1579; 		/* rule 825 */
	 tx1581 = tx1047*tx1580*tx1815; 		/* rule 826 */
	 tx1582 = tx1536 + tx1581; 		/* rule 827 */
	 tx1583 = DeDPhi*tx1582; 		/* rule 828 */
	 tx1584 = tx1003*tx1045*tx1048*tx998*V; 		/* rule 829 */
	 ohz1x2 = tx1583 + tx1584; 		/* rule 830 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 tx1585 = tx1333*tx1813*tx926; 		/* rule 832 */
	 tx1586 = tx1238*tx925*tx994; 		/* rule 833 */
	 tx1587 = tx1585 + tx1586; 		/* rule 834 */
	 tx1588 = tx1047*tx1587*tx1815; 		/* rule 835 */
	 tx1589 = tx1567 + tx1588; 		/* rule 836 */
	 tx1590 = DeDPhi*tx1589; 		/* rule 837 */
	 tx1591 = tx1002*tx1003*tx1045*tx1048*V; 		/* rule 838 */
	 ohz1y2 = tx1590 + tx1591; 		/* rule 839 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx1592 = -2.*tx1047*tx1100*tx925*tx971; 		/* rule 841 */
	 tx1593 = tx1047*tx1379*tx1815; 		/* rule 842 */
	 tx1594 = tx1592 + tx1593; 		/* rule 843 */
	 tx1595 = DeDPhi*tx1594; 		/* rule 844 */
	 tx1596 = tx1003*tx1006*tx1045*tx1048*V; 		/* rule 845 */
	 ohz1z2 = tx1595 + tx1596; 		/* rule 846 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 tx1597 = -(DeDPhi*tx1581); 		/* rule 848 */
	 tx1598 = tx1003*tx1010*tx1045*tx1048*V; 		/* rule 849 */
	 ohz1x3 = tx1597 + tx1598; 		/* rule 850 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 0, ohz1x3);
	 tx1599 = -(DeDPhi*tx1588); 		/* rule 852 */
	 tx1600 = tx1003*tx1014*tx1045*tx1048*V; 		/* rule 853 */
	 ohz1y3 = tx1599 + tx1600; 		/* rule 854 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 1, ohz1y3);
	 tx1601 = -(DeDPhi*tx1593); 		/* rule 856 */
	 tx1602 = tx1003*tx1018*tx1045*tx1048*V; 		/* rule 857 */
	 ohz1z3 = tx1601 + tx1602; 		/* rule 858 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 2, ohz1z3);
	 ohz1x4 = tx1045*tx1048*tx1813*tx925*tx946*tx978*tx991*V; 		/* rule 860 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 0, ohz1x4);
	 ohz1y4 = tx1045*tx1048*tx1813*tx925*tx945*tx978*tx991*V; 		/* rule 862 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 1, ohz1y4);
	 ohz1z4 = tx1045*tx1048*tx1813*tx951*tx978*tx991*V; 		/* rule 864 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 2, ohz1z4);
	 tx1603 = tx1188*tx926; 		/* rule 866 */
	 tx1604 = tx1286*tx927; 		/* rule 867 */
	 tx1605 = tx1189*tx945; 		/* rule 868 */
	 tx1606 = tx1190*tx945; 		/* rule 869 */
	 tx1607 = tx1287*tx946; 		/* rule 870 */
	 tx1608 = tx1288*tx946; 		/* rule 871 */
	 tx1609 = tx1603 + tx1604; 		/* rule 872 */
	 tx1610 = tx1605 + tx1607; 		/* rule 873 */
	 tx1611 = tx1606 + tx1608; 		/* rule 874 */
	 tx1612 = -tx1529; 		/* rule 875 */
	 tx1613 = -tx1555; 		/* rule 876 */
	 tx1614 = tx1047*tx1512*tx971; 		/* rule 877 */
	 tx1615 = 0.5*tx1192*tx1512*tx978; 		/* rule 878 */
	 tx1616 = -0.5*tx1192*tx1610*tx991; 		/* rule 879 */
	 tx1617 = -(tx1191*tx1611*tx1815*tx993); 		/* rule 880 */
	 tx1618 = tx1047*tx1609*tx1815*tx994; 		/* rule 881 */
	 tx1619 = tx1612 + tx1613 + tx1614 + tx1615 + tx1616 + tx1617 + tx1618; 		/* rule 882 */
	 tx1620 = DeDPhi*tx1619; 		/* rule 883 */
	 tx1621 = -(tx1002*tx1045*tx1048*tx998*V); 		/* rule 884 */
	 ohx2y2 = tx1620 + tx1621; 		/* rule 885 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 tx1622 = tx1188*tx925; 		/* rule 887 */
	 tx1623 = tx1380*tx927; 		/* rule 888 */
	 tx1624 = tx1189*tx944; 		/* rule 889 */
	 tx1625 = tx1190*tx944; 		/* rule 890 */
	 tx1626 = tx1381*tx946; 		/* rule 891 */
	 tx1627 = tx1382*tx946; 		/* rule 892 */
	 tx1628 = tx1622 + tx1623; 		/* rule 893 */
	 tx1629 = tx1624 + tx1626; 		/* rule 894 */
	 tx1630 = tx1625 + tx1627; 		/* rule 895 */
	 tx1631 = -tx1537; 		/* rule 896 */
	 tx1632 = -tx1581; 		/* rule 897 */
	 tx1633 = tx1047*tx1517*tx971; 		/* rule 898 */
	 tx1634 = 0.5*tx1192*tx1517*tx978; 		/* rule 899 */
	 tx1635 = -0.5*tx1192*tx1629*tx991; 		/* rule 900 */
	 tx1636 = -(tx1191*tx1630*tx1815*tx993); 		/* rule 901 */
	 tx1637 = tx1047*tx1628*tx1815*tx994; 		/* rule 902 */
	 tx1638 = tx1631 + tx1632 + tx1633 + tx1634 + tx1635 + tx1636 + tx1637; 		/* rule 903 */
	 tx1639 = DeDPhi*tx1638; 		/* rule 904 */
	 tx1640 = -(tx1006*tx1045*tx1048*tx998*V); 		/* rule 905 */
	 ohx2z2 = tx1639 + tx1640; 		/* rule 906 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 tx1641 = -tx1194; 		/* rule 908 */
	 tx1642 = tx1189*tx1192*tx946*tx991; 		/* rule 909 */
	 tx1643 = tx1191*tx1423*tx1815; 		/* rule 910 */
	 tx1644 = 2.*tx1190*tx1191*tx1815*tx946*tx993; 		/* rule 911 */
	 tx1645 = -2.*tx1047*tx1188*tx1815*tx927*tx994; 		/* rule 912 */
	 tx1646 = tx1521 + tx1641 + tx1642 + tx1643 + tx1644 + tx1645; 		/* rule 913 */
	 tx1647 = DeDPhi*tx1646; 		/* rule 914 */
	 tx1648 = -(tx1010*tx1045*tx1048*tx998*V); 		/* rule 915 */
	 ohx2x3 = tx1647 + tx1648; 		/* rule 916 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 0, ohx2x3);
	 tx1649 = tx1456*tx1813*tx946; 		/* rule 918 */
	 tx1650 = tx1420*tx945*tx993; 		/* rule 919 */
	 tx1651 = tx1649 + tx1650; 		/* rule 920 */
	 tx1652 = -tx1618; 		/* rule 921 */
	 tx1653 = -0.5*tx1192*tx1512*tx978; 		/* rule 922 */
	 tx1654 = 0.5*tx1192*tx1610*tx991; 		/* rule 923 */
	 tx1655 = tx1191*tx1651*tx1815; 		/* rule 924 */
	 tx1656 = tx1191*tx1611*tx1815*tx993; 		/* rule 925 */
	 tx1657 = tx1529 + tx1652 + tx1653 + tx1654 + tx1655 + tx1656; 		/* rule 926 */
	 tx1658 = DeDPhi*tx1657; 		/* rule 927 */
	 tx1659 = -(tx1014*tx1045*tx1048*tx998*V); 		/* rule 928 */
	 ohx2y3 = tx1658 + tx1659; 		/* rule 929 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 1, ohx2y3);
	 tx1660 = tx1493*tx1813*tx946; 		/* rule 931 */
	 tx1661 = tx1420*tx944*tx993; 		/* rule 932 */
	 tx1662 = tx1660 + tx1661; 		/* rule 933 */
	 tx1663 = -tx1637; 		/* rule 934 */
	 tx1664 = -0.5*tx1192*tx1517*tx978; 		/* rule 935 */
	 tx1665 = 0.5*tx1192*tx1629*tx991; 		/* rule 936 */
	 tx1666 = tx1191*tx1662*tx1815; 		/* rule 937 */
	 tx1667 = tx1191*tx1630*tx1815*tx993; 		/* rule 938 */
	 tx1668 = tx1537 + tx1663 + tx1664 + tx1665 + tx1666 + tx1667; 		/* rule 939 */
	 tx1669 = DeDPhi*tx1668; 		/* rule 940 */
	 tx1670 = -(tx1018*tx1045*tx1048*tx998*V); 		/* rule 941 */
	 ohx2z3 = tx1669 + tx1670; 		/* rule 942 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 2, ohx2z3);
	 tx1671 = -(DeDPhi*tx1643); 		/* rule 944 */
	 tx1672 = tx1007*tx1045*tx1048*tx998*V; 		/* rule 945 */
	 ohx2x4 = tx1671 + tx1672; 		/* rule 946 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 0, ohx2x4);
	 tx1673 = -(DeDPhi*tx1655); 		/* rule 948 */
	 tx1674 = tx1012*tx1045*tx1048*tx998*V; 		/* rule 949 */
	 ohx2y4 = tx1673 + tx1674; 		/* rule 950 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 1, ohx2y4);
	 tx1675 = -(DeDPhi*tx1666); 		/* rule 952 */
	 tx1676 = tx1016*tx1045*tx1048*tx998*V; 		/* rule 953 */
	 ohx2z4 = tx1675 + tx1676; 		/* rule 954 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 2, ohx2z4);
	 tx1677 = tx1286*tx925; 		/* rule 956 */
	 tx1678 = tx1380*tx926; 		/* rule 957 */
	 tx1679 = tx1287*tx944; 		/* rule 958 */
	 tx1680 = tx1288*tx944; 		/* rule 959 */
	 tx1681 = tx1381*tx945; 		/* rule 960 */
	 tx1682 = tx1382*tx945; 		/* rule 961 */
	 tx1683 = tx1677 + tx1678; 		/* rule 962 */
	 tx1684 = tx1679 + tx1681; 		/* rule 963 */
	 tx1685 = tx1680 + tx1682; 		/* rule 964 */
	 tx1686 = -tx1568; 		/* rule 965 */
	 tx1687 = -tx1588; 		/* rule 966 */
	 tx1688 = tx1047*tx1549*tx971; 		/* rule 967 */
	 tx1689 = 0.5*tx1192*tx1549*tx978; 		/* rule 968 */
	 tx1690 = -0.5*tx1192*tx1684*tx991; 		/* rule 969 */
	 tx1691 = -(tx1191*tx1685*tx1815*tx993); 		/* rule 970 */
	 tx1692 = tx1047*tx1683*tx1815*tx994; 		/* rule 971 */
	 tx1693 = tx1686 + tx1687 + tx1688 + tx1689 + tx1690 + tx1691 + tx1692; 		/* rule 972 */
	 tx1694 = DeDPhi*tx1693; 		/* rule 973 */
	 tx1695 = -(tx1002*tx1006*tx1045*tx1048*V); 		/* rule 974 */
	 ohy2z2 = tx1694 + tx1695; 		/* rule 975 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	 tx1696 = tx1419*tx1813*tx945; 		/* rule 977 */
	 tx1697 = tx1457*tx946*tx993; 		/* rule 978 */
	 tx1698 = tx1696 + tx1697; 		/* rule 979 */
	 tx1699 = tx1191*tx1698*tx1815; 		/* rule 980 */
	 tx1700 = tx1555 + tx1652 + tx1653 + tx1654 + tx1656 + tx1699; 		/* rule 981 */
	 tx1701 = DeDPhi*tx1700; 		/* rule 982 */
	 tx1702 = -(tx1002*tx1010*tx1045*tx1048*V); 		/* rule 983 */
	 ohy2x3 = tx1701 + tx1702; 		/* rule 984 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 0, ohy2x3);
	 tx1703 = -tx1290; 		/* rule 986 */
	 tx1704 = tx1192*tx1287*tx945*tx991; 		/* rule 987 */
	 tx1705 = tx1191*tx1460*tx1815; 		/* rule 988 */
	 tx1706 = 2.*tx1191*tx1288*tx1815*tx945*tx993; 		/* rule 989 */
	 tx1707 = -2.*tx1047*tx1286*tx1815*tx926*tx994; 		/* rule 990 */
	 tx1708 = tx1560 + tx1703 + tx1704 + tx1705 + tx1706 + tx1707; 		/* rule 991 */
	 tx1709 = DeDPhi*tx1708; 		/* rule 992 */
	 tx1710 = -(tx1002*tx1014*tx1045*tx1048*V); 		/* rule 993 */
	 ohy2y3 = tx1709 + tx1710; 		/* rule 994 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 1, ohy2y3);
	 tx1711 = tx1493*tx1813*tx945; 		/* rule 996 */
	 tx1712 = tx1457*tx944*tx993; 		/* rule 997 */
	 tx1713 = tx1711 + tx1712; 		/* rule 998 */
	 tx1714 = -tx1692; 		/* rule 999 */
	 tx1715 = -0.5*tx1192*tx1549*tx978; 		/* rule 1000 */
	 tx1716 = 0.5*tx1192*tx1684*tx991; 		/* rule 1001 */
	 tx1717 = tx1191*tx1713*tx1815; 		/* rule 1002 */
	 tx1718 = tx1191*tx1685*tx1815*tx993; 		/* rule 1003 */
	 tx1719 = tx1568 + tx1714 + tx1715 + tx1716 + tx1717 + tx1718; 		/* rule 1004 */
	 tx1720 = DeDPhi*tx1719; 		/* rule 1005 */
	 tx1721 = -(tx1002*tx1018*tx1045*tx1048*V); 		/* rule 1006 */
	 ohy2z3 = tx1720 + tx1721; 		/* rule 1007 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 2, ohy2z3);
	 tx1722 = -(DeDPhi*tx1699); 		/* rule 1009 */
	 tx1723 = tx1002*tx1007*tx1045*tx1048*V; 		/* rule 1010 */
	 ohy2x4 = tx1722 + tx1723; 		/* rule 1011 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 0, ohy2x4);
	 tx1724 = -(DeDPhi*tx1705); 		/* rule 1013 */
	 tx1725 = tx1002*tx1012*tx1045*tx1048*V; 		/* rule 1014 */
	 ohy2y4 = tx1724 + tx1725; 		/* rule 1015 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 1, ohy2y4);
	 tx1726 = -(DeDPhi*tx1717); 		/* rule 1017 */
	 tx1727 = tx1002*tx1016*tx1045*tx1048*V; 		/* rule 1018 */
	 ohy2z4 = tx1726 + tx1727; 		/* rule 1019 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 2, ohy2z4);
	 tx1728 = tx1419*tx1813*tx944; 		/* rule 1021 */
	 tx1729 = tx1494*tx946*tx993; 		/* rule 1022 */
	 tx1730 = tx1728 + tx1729; 		/* rule 1023 */
	 tx1731 = tx1191*tx1730*tx1815; 		/* rule 1024 */
	 tx1732 = tx1581 + tx1663 + tx1664 + tx1665 + tx1667 + tx1731; 		/* rule 1025 */
	 tx1733 = DeDPhi*tx1732; 		/* rule 1026 */
	 tx1734 = -(tx1006*tx1010*tx1045*tx1048*V); 		/* rule 1027 */
	 ohz2x3 = tx1733 + tx1734; 		/* rule 1028 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 0, ohz2x3);
	 tx1735 = tx1456*tx1813*tx944; 		/* rule 1030 */
	 tx1736 = tx1494*tx945*tx993; 		/* rule 1031 */
	 tx1737 = tx1735 + tx1736; 		/* rule 1032 */
	 tx1738 = tx1191*tx1737*tx1815; 		/* rule 1033 */
	 tx1739 = tx1588 + tx1714 + tx1715 + tx1716 + tx1718 + tx1738; 		/* rule 1034 */
	 tx1740 = DeDPhi*tx1739; 		/* rule 1035 */
	 tx1741 = -(tx1006*tx1014*tx1045*tx1048*V); 		/* rule 1036 */
	 ohz2y3 = tx1740 + tx1741; 		/* rule 1037 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 1, ohz2y3);
	 tx1742 = -tx1384; 		/* rule 1039 */
	 tx1743 = tx1192*tx1381*tx944*tx991; 		/* rule 1040 */
	 tx1744 = tx1191*tx1497*tx1815; 		/* rule 1041 */
	 tx1745 = 2.*tx1191*tx1382*tx1815*tx944*tx993; 		/* rule 1042 */
	 tx1746 = -2.*tx1047*tx1380*tx1815*tx925*tx994; 		/* rule 1043 */
	 tx1747 = tx1593 + tx1742 + tx1743 + tx1744 + tx1745 + tx1746; 		/* rule 1044 */
	 tx1748 = DeDPhi*tx1747; 		/* rule 1045 */
	 tx1749 = -(tx1006*tx1018*tx1045*tx1048*V); 		/* rule 1046 */
	 ohz2z3 = tx1748 + tx1749; 		/* rule 1047 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 2, ohz2z3);
	 tx1750 = -(DeDPhi*tx1731); 		/* rule 1049 */
	 tx1751 = tx1006*tx1007*tx1045*tx1048*V; 		/* rule 1050 */
	 ohz2x4 = tx1750 + tx1751; 		/* rule 1051 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 0, ohz2x4);
	 tx1752 = -(DeDPhi*tx1738); 		/* rule 1053 */
	 tx1753 = tx1006*tx1012*tx1045*tx1048*V; 		/* rule 1054 */
	 ohz2y4 = tx1752 + tx1753; 		/* rule 1055 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 1, ohz2y4);
	 tx1754 = -(DeDPhi*tx1744); 		/* rule 1057 */
	 tx1755 = tx1006*tx1016*tx1045*tx1048*V; 		/* rule 1058 */
	 ohz2z4 = tx1754 + tx1755; 		/* rule 1059 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 2, ohz2z4);
	 tx1756 = -tx1655; 		/* rule 1061 */
	 tx1757 = -tx1699; 		/* rule 1062 */
	 tx1758 = -(tx1191*tx1610*tx971); 		/* rule 1063 */
	 tx1759 = tx1615 + tx1616 + tx1617 + tx1618 + tx1756 + tx1757 + tx1758; 		/* rule 1064 */
	 tx1760 = DeDPhi*tx1759; 		/* rule 1065 */
	 tx1761 = -(tx1010*tx1014*tx1045*tx1048*V); 		/* rule 1066 */
	 ohx3y3 = tx1760 + tx1761; 		/* rule 1067 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 1, ohx3y3);
	 tx1762 = -tx1666; 		/* rule 1069 */
	 tx1763 = -tx1731; 		/* rule 1070 */
	 tx1764 = -(tx1191*tx1629*tx971); 		/* rule 1071 */
	 tx1765 = tx1634 + tx1635 + tx1636 + tx1637 + tx1762 + tx1763 + tx1764; 		/* rule 1072 */
	 tx1766 = DeDPhi*tx1765; 		/* rule 1073 */
	 tx1767 = -(tx1010*tx1018*tx1045*tx1048*V); 		/* rule 1074 */
	 ohx3z3 = tx1766 + tx1767; 		/* rule 1075 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 2, ohx3z3);
	 tx1768 = 2.*tx1189*tx1191*tx946*tx971; 		/* rule 1077 */
	 tx1769 = tx1643 + tx1768; 		/* rule 1078 */
	 tx1770 = DeDPhi*tx1769; 		/* rule 1079 */
	 tx1771 = tx1007*tx1010*tx1045*tx1048*V; 		/* rule 1080 */
	 ohx3x4 = tx1770 + tx1771; 		/* rule 1081 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 0, ohx3x4);
	 tx1772 = tx1191*tx1610*tx971; 		/* rule 1083 */
	 tx1773 = tx1655 + tx1772; 		/* rule 1084 */
	 tx1774 = DeDPhi*tx1773; 		/* rule 1085 */
	 tx1775 = tx1010*tx1012*tx1045*tx1048*V; 		/* rule 1086 */
	 ohx3y4 = tx1774 + tx1775; 		/* rule 1087 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 1, ohx3y4);
	 tx1776 = tx1191*tx1629*tx971; 		/* rule 1089 */
	 tx1777 = tx1666 + tx1776; 		/* rule 1090 */
	 tx1778 = DeDPhi*tx1777; 		/* rule 1091 */
	 tx1779 = tx1010*tx1016*tx1045*tx1048*V; 		/* rule 1092 */
	 ohx3z4 = tx1778 + tx1779; 		/* rule 1093 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 2, ohx3z4);
	 tx1780 = -tx1717; 		/* rule 1095 */
	 tx1781 = -tx1738; 		/* rule 1096 */
	 tx1782 = -(tx1191*tx1684*tx971); 		/* rule 1097 */
	 tx1783 = tx1689 + tx1690 + tx1691 + tx1692 + tx1780 + tx1781 + tx1782; 		/* rule 1098 */
	 tx1784 = DeDPhi*tx1783; 		/* rule 1099 */
	 tx1785 = -(tx1014*tx1018*tx1045*tx1048*V); 		/* rule 1100 */
	 ohy3z3 = tx1784 + tx1785; 		/* rule 1101 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 2, ohy3z3);
	 tx1786 = tx1699 + tx1772; 		/* rule 1103 */
	 tx1787 = DeDPhi*tx1786; 		/* rule 1104 */
	 tx1788 = tx1007*tx1014*tx1045*tx1048*V; 		/* rule 1105 */
	 ohy3x4 = tx1787 + tx1788; 		/* rule 1106 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 0, ohy3x4);
	 tx1789 = 2.*tx1191*tx1287*tx945*tx971; 		/* rule 1108 */
	 tx1790 = tx1705 + tx1789; 		/* rule 1109 */
	 tx1791 = DeDPhi*tx1790; 		/* rule 1110 */
	 tx1792 = tx1012*tx1014*tx1045*tx1048*V; 		/* rule 1111 */
	 ohy3y4 = tx1791 + tx1792; 		/* rule 1112 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 1, ohy3y4);
	 tx1793 = tx1191*tx1684*tx971; 		/* rule 1114 */
	 tx1794 = tx1717 + tx1793; 		/* rule 1115 */
	 tx1795 = DeDPhi*tx1794; 		/* rule 1116 */
	 tx1796 = tx1014*tx1016*tx1045*tx1048*V; 		/* rule 1117 */
	 ohy3z4 = tx1795 + tx1796; 		/* rule 1118 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 2, ohy3z4);
	 tx1797 = tx1731 + tx1776; 		/* rule 1120 */
	 tx1798 = DeDPhi*tx1797; 		/* rule 1121 */
	 tx1799 = tx1007*tx1018*tx1045*tx1048*V; 		/* rule 1122 */
	 ohz3x4 = tx1798 + tx1799; 		/* rule 1123 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 0, ohz3x4);
	 tx1800 = tx1738 + tx1793; 		/* rule 1125 */
	 tx1801 = DeDPhi*tx1800; 		/* rule 1126 */
	 tx1802 = tx1012*tx1018*tx1045*tx1048*V; 		/* rule 1127 */
	 ohz3y4 = tx1801 + tx1802; 		/* rule 1128 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 1, ohz3y4);
	 tx1803 = 2.*tx1191*tx1381*tx944*tx971; 		/* rule 1130 */
	 tx1804 = tx1744 + tx1803; 		/* rule 1131 */
	 tx1805 = DeDPhi*tx1804; 		/* rule 1132 */
	 tx1806 = tx1016*tx1018*tx1045*tx1048*V; 		/* rule 1133 */
	 ohz3z4 = tx1805 + tx1806; 		/* rule 1134 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 2, ohz3z4);
	 tx1807 = DeDPhi*tx1758; 		/* rule 1136 */
	 tx1808 = -(tx1045*tx1048*tx1191*tx1813*tx945*tx946*V); 		/* rule 1137 */
	 ohx4y4 = tx1807 + tx1808; 		/* rule 1138 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 1, ohx4y4);
	 tx1809 = DeDPhi*tx1764; 		/* rule 1140 */
	 tx1810 = -(tx1045*tx1048*tx1191*tx1813*tx944*tx946*V); 		/* rule 1141 */
	 ohx4z4 = tx1809 + tx1810; 		/* rule 1142 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 2, ohx4z4);
	 tx1811 = DeDPhi*tx1782; 		/* rule 1144 */
	 tx1812 = -(tx1045*tx1048*tx1191*tx1813*tx944*tx945*V); 		/* rule 1145 */
	 ohy4z4 = tx1811 + tx1812; 		/* rule 1146 */
	DIHEDRAL_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 1, I4, 2, ohy4z4);
	} /*calcOffDiagonalHessian*/
	#endif // DIHEDRAL_CALC_OFF_DIAGONAL_HESSIAN ]
	} /*calcDiagonalHessian*/
	#endif // DIHEDRAL_CALC_DIAGONAL_HESSIAN ]
	} /*calcForce*/
	#endif // DIHEDRAL_CALC_FORCE ]
