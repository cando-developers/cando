// C-code
	IMPROPER_RESTRAINT_SET_PARAMETER(K);
	IMPROPER_RESTRAINT_SET_PARAMETER(U);
	IMPROPER_RESTRAINT_SET_PARAMETER(L);
	IMPROPER_RESTRAINT_SET_PARAMETER(I1);
	IMPROPER_RESTRAINT_SET_PARAMETER(I2);
	IMPROPER_RESTRAINT_SET_PARAMETER(I3);
	IMPROPER_RESTRAINT_SET_PARAMETER(I4);
	IMPROPER_RESTRAINT_SET_POSITION(x1,I1,0);
	IMPROPER_RESTRAINT_SET_POSITION(y1,I1,1);
	IMPROPER_RESTRAINT_SET_POSITION(z1,I1,2);
	IMPROPER_RESTRAINT_SET_POSITION(x2,I2,0);
	IMPROPER_RESTRAINT_SET_POSITION(y2,I2,1);
	IMPROPER_RESTRAINT_SET_POSITION(z2,I2,2);
	IMPROPER_RESTRAINT_SET_POSITION(x3,I3,0);
	IMPROPER_RESTRAINT_SET_POSITION(y3,I3,1);
	IMPROPER_RESTRAINT_SET_POSITION(z3,I3,2);
	IMPROPER_RESTRAINT_SET_POSITION(x4,I4,0);
	IMPROPER_RESTRAINT_SET_POSITION(y4,I4,1);
	IMPROPER_RESTRAINT_SET_POSITION(z4,I4,2);
/*Energy= 0.0;
fx1= 0.0;
fy1= 0.0;
fz1= 0.0;
fx2= 0.0;
fy2= 0.0;
fz2= 0.0;
fx3= 0.0;
fy3= 0.0;
fz3= 0.0;
fx4= 0.0;
fy4= 0.0;
fz4= 0.0;
dhx1x1= 0.0;
dhy1y1= 0.0;
dhz1z1= 0.0;
dhx2x2= 0.0;
dhy2y2= 0.0;
dhz2z2= 0.0;
dhx3x3= 0.0;
dhy3y3= 0.0;
dhz3z3= 0.0;
dhx4x4= 0.0;
dhy4y4= 0.0;
dhz4z4= 0.0;
ohx1y1= 0.0;
ohx1z1= 0.0;
ohx1x2= 0.0;
ohx1y2= 0.0;
ohx1z2= 0.0;
ohx1x3= 0.0;
ohx1y3= 0.0;
ohx1z3= 0.0;
ohx1x4= 0.0;
ohx1y4= 0.0;
ohx1z4= 0.0;
ohy1z1= 0.0;
ohy1x2= 0.0;
ohy1y2= 0.0;
ohy1z2= 0.0;
ohy1x3= 0.0;
ohy1y3= 0.0;
ohy1z3= 0.0;
ohy1x4= 0.0;
ohy1y4= 0.0;
ohy1z4= 0.0;
ohz1x2= 0.0;
ohz1y2= 0.0;
ohz1z2= 0.0;
ohz1x3= 0.0;
ohz1y3= 0.0;
ohz1z3= 0.0;
ohz1x4= 0.0;
ohz1y4= 0.0;
ohz1z4= 0.0;
ohx2y2= 0.0;
ohx2z2= 0.0;
ohx2x3= 0.0;
ohx2y3= 0.0;
ohx2z3= 0.0;
ohx2x4= 0.0;
ohx2y4= 0.0;
ohx2z4= 0.0;
ohy2z2= 0.0;
ohy2x3= 0.0;
ohy2y3= 0.0;
ohy2z3= 0.0;
ohy2x4= 0.0;
ohy2y4= 0.0;
ohy2z4= 0.0;
ohz2x3= 0.0;
ohz2y3= 0.0;
ohz2z3= 0.0;
ohz2x4= 0.0;
ohz2y4= 0.0;
ohz2z4= 0.0;
ohx3y3= 0.0;
ohx3z3= 0.0;
ohx3x4= 0.0;
ohx3y4= 0.0;
ohx3z4= 0.0;
ohy3z3= 0.0;
ohy3x4= 0.0;
ohy3y4= 0.0;
ohy3z4= 0.0;
ohz3x4= 0.0;
ohz3y4= 0.0;
ohz3z4= 0.0;
ohx4y4= 0.0;
ohx4z4= 0.0;
ohy4z4= 0.0;
*/
	 tx916 = -(x2*y1); 		/* rule 21 */
	 tx917 = x3*y1; 		/* rule 22 */
	 tx918 = x1*y2; 		/* rule 23 */
	 tx919 = -(x3*y2); 		/* rule 24 */
	 tx920 = -(x1*y3); 		/* rule 25 */
	 tx921 = x2*y3; 		/* rule 26 */
	 tx922 = x2*z1; 		/* rule 27 */
	 tx923 = -(x3*z1); 		/* rule 28 */
	 tx924 = -(y2*z1); 		/* rule 29 */
	 tx925 = y3*z1; 		/* rule 30 */
	 tx926 = -(x1*z2); 		/* rule 31 */
	 tx927 = x3*z2; 		/* rule 32 */
	 tx928 = y1*z2; 		/* rule 33 */
	 tx929 = -(y3*z2); 		/* rule 34 */
	 tx930 = x1*z3; 		/* rule 35 */
	 tx931 = -(x2*z3); 		/* rule 36 */
	 tx932 = -(y1*z3); 		/* rule 37 */
	 tx933 = y2*z3; 		/* rule 38 */
	 tx934 = tx916 + tx917 + tx918 + tx919 + tx920 + tx921; 		/* rule 39 */
	 tx935 = tx922 + tx923 + tx926 + tx927 + tx930 + tx931; 		/* rule 40 */
	 tx936 = tx924 + tx925 + tx928 + tx929 + tx932 + tx933; 		/* rule 41 */
	 tx937 = power2(tx934); 		/* rule 42 */
	 tx938 = power2(tx935); 		/* rule 43 */
	 tx939 = power2(tx936); 		/* rule 44 */
	 tx940 = tx937 + tx938 + tx939; 		/* rule 45 */
	 LenA = mysqrt(tx940); 		/* rule 46 */
	 tx941 = x4*y2; 		/* rule 47 */
	 tx942 = -(x4*y3); 		/* rule 48 */
	 tx943 = -(x2*y4); 		/* rule 49 */
	 tx944 = x3*y4; 		/* rule 50 */
	 tx945 = -(x4*z2); 		/* rule 51 */
	 tx946 = y4*z2; 		/* rule 52 */
	 tx947 = x4*z3; 		/* rule 53 */
	 tx948 = -(y4*z3); 		/* rule 54 */
	 tx949 = x2*z4; 		/* rule 55 */
	 tx950 = -(x3*z4); 		/* rule 56 */
	 tx951 = -(y2*z4); 		/* rule 57 */
	 tx952 = y3*z4; 		/* rule 58 */
	 tx953 = tx919 + tx921 + tx941 + tx942 + tx943 + tx944; 		/* rule 59 */
	 tx954 = tx927 + tx931 + tx945 + tx947 + tx949 + tx950; 		/* rule 60 */
	 tx955 = tx929 + tx933 + tx946 + tx948 + tx951 + tx952; 		/* rule 61 */
	 tx956 = power2(tx953); 		/* rule 62 */
	 tx957 = power2(tx954); 		/* rule 63 */
	 tx958 = power2(tx955); 		/* rule 64 */
	 tx959 = tx956 + tx957 + tx958; 		/* rule 65 */
	 LenB = mysqrt(tx959); 		/* rule 66 */
	 ReciprocalLenA = reciprocal(LenA); 		/* rule 67 */
	 ReciprocalLenB = reciprocal(LenB); 		/* rule 68 */
	if (fabs(LenA)<TENM3) ReciprocalLenA = 0.0;
	if (fabs(LenB)<TENM3) ReciprocalLenB = 0.0;
	 RecLenARecLenB = ReciprocalLenA*ReciprocalLenB; 		/* rule 71 */
	EraseLinearDihedral = 1.0;
	if (RecLenARecLenB==0.0) EraseLinearDihedral = 0.0;
	 tx960 = tx934*tx953; 		/* rule 74 */
	 tx961 = tx935*tx954; 		/* rule 75 */
	 tx962 = tx936*tx955; 		/* rule 76 */
	 tx963 = tx960 + tx961 + tx962; 		/* rule 77 */
	 CosPhi = RecLenARecLenB*tx963; 		/* rule 78 */
	 tx964 = -x3; 		/* rule 79 */
	 tx965 = -y3; 		/* rule 80 */
	 tx966 = -z3; 		/* rule 81 */
	 tx967 = tx964 + x2; 		/* rule 82 */
	 tx968 = tx965 + y2; 		/* rule 83 */
	 tx969 = tx966 + z2; 		/* rule 84 */
	 tx970 = power2(tx967); 		/* rule 85 */
	 tx971 = power2(tx968); 		/* rule 86 */
	 tx972 = power2(tx969); 		/* rule 87 */
	 tx973 = tx964 + x4; 		/* rule 88 */
	 tx974 = tx965 + y4; 		/* rule 89 */
	 tx975 = tx966 + z4; 		/* rule 90 */
	 tx976 = tx970 + tx971 + tx972; 		/* rule 91 */
	 tx977 = tx936*tx973; 		/* rule 92 */
	 tx978 = tx935*tx974; 		/* rule 93 */
	 tx979 = tx934*tx975; 		/* rule 94 */
	 tx980 = mysqrt(tx976); 		/* rule 95 */
	 tx981 = tx977 + tx978 + tx979; 		/* rule 96 */
	 SinPhi = RecLenARecLenB*tx980*tx981; 		/* rule 97 */
	CosPhi=MAX(-1.0,MIN(1.0,CosPhi));
/*If[CosPhi>0.1,Phi=ArcSin[SinPhi],
					Phi=ArcCos[CosPhi]*Sign[SinPhi]];*/
	if(CosPhi>0.1){Phi=asin(SinPhi);}else
					{Phi=acos(CosPhi)*SIGN(SinPhi);}
/*RestraintActive = False;*/
	RestraintActive=false;
/*UShift = 0.0; PhiShift=0.0;*/
	UShift = 0.0; PhiShift = 0.0;
/*
If[ L<U && ((L<=Phi)&&(Phi<=U)), 
	RestraintActive = True;
	UShift = 0.0;
	PhiShift = 0.0;
];
*/
	
if( L<U && (L<Phi)&&(Phi<U) ) {
	RestraintActive = true;
	UShift = 0.0;
	PhiShift = 0.0;
}
/*
If[ U<=L && ((Phi<=U)||(Phi>=L)),
	RestraintActive = True;
	UShift = TWOPI;
	If[Phi>=L,
		PhiShift = 0.0;
	,
		PhiShift = TWOPI;
	];
];*/
	
if ( U<=L && ((Phi<=U)||(Phi>=L))) {
	RestraintActive = true;
	UShift = TWOPI;
	if ( Phi>=L ) {
		PhiShift = 0.0;
	} else {
		PhiShift = TWOPI;
	}
}
/*If[RestraintActive,*/
	if (RestraintActive ) {
	 tx982 = -L; 		/* rule 111 */
	 tx983 = tx982 + U + UShift; 		/* rule 112 */
	 tx984 = Phi + PhiShift + tx982; 		/* rule 113 */
	 tx985 = reciprocal(tx983); 		/* rule 114 */
	 Psi = TWOPI*tx984*tx985; 		/* rule 115 */
	 tx986 = Cos(Psi); 		/* rule 116 */
	 tx987 = -tx986; 		/* rule 117 */
	 tx988 = 1 + tx987; 		/* rule 118 */
	 tx1831 = power2(tx988); 		/* rule 119 */
	 tx989 = tx1831*tx988; 		/* rule 120 */
	 Energy = EraseLinearDihedral*K*tx989; 		/* rule 121 */
	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE(Energy);
	#ifdef IMPROPER_RESTRAINT_CALC_FORCE //[
	if (calcForce ) {
	 tx990 = tx1831; 		/* rule 125 */
	 tx991 = Sin(Psi); 		/* rule 126 */
	 DEImproperRestraintDPhi = 3*EraseLinearDihedral*K*TWOPI*tx985*tx990*tx991; 		/* rule 127 */
	 tx992 = reciprocal(tx940); 		/* rule 128 */
	 gx1 = -(DEImproperRestraintDPhi*tx936*tx980*tx992); 		/* rule 129 */
	 fx1 = -gx1; 		/* rule 130 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I1, 0, fx1 );
	 gy1 = -(DEImproperRestraintDPhi*tx935*tx980*tx992); 		/* rule 132 */
	 fy1 = -gy1; 		/* rule 133 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I1, 1, fy1 );
	 gz1 = -(DEImproperRestraintDPhi*tx934*tx980*tx992); 		/* rule 135 */
	 fz1 = -gz1; 		/* rule 136 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tx993 = -x2; 		/* rule 138 */
	 tx994 = -y2; 		/* rule 139 */
	 tx995 = -z2; 		/* rule 140 */
	 tx996 = tx993 + x1; 		/* rule 141 */
	 tx997 = tx994 + y1; 		/* rule 142 */
	 tx998 = tx995 + z1; 		/* rule 143 */
	 tx999 = tx967*tx973; 		/* rule 144 */
	 tx1000 = tx968*tx974; 		/* rule 145 */
	 tx1001 = tx969*tx975; 		/* rule 146 */
	 tx1002 = tx967*tx996; 		/* rule 147 */
	 tx1003 = tx968*tx997; 		/* rule 148 */
	 tx1004 = tx969*tx998; 		/* rule 149 */
	 tx1005 = tx1002 + tx1003 + tx1004; 		/* rule 150 */
	 tx1006 = reciprocal(tx959); 		/* rule 151 */
	 tx1007 = reciprocal(tx980); 		/* rule 152 */
	 tx1008 = tx1000 + tx1001 + tx999; 		/* rule 153 */
	 tx1009 = -(tx1006*tx1007*tx1008*tx955); 		/* rule 154 */
	 tx1010 = tx1005*tx1007*tx936*tx992; 		/* rule 155 */
	 tx1011 = tx936*tx980*tx992; 		/* rule 156 */
	 tx1012 = tx1009 + tx1010 + tx1011; 		/* rule 157 */
	 gx2 = DEImproperRestraintDPhi*tx1012; 		/* rule 158 */
	 fx2 = -gx2; 		/* rule 159 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I2, 0, fx2 );
	 tx1013 = -(tx1006*tx1007*tx1008*tx954); 		/* rule 161 */
	 tx1014 = tx1005*tx1007*tx935*tx992; 		/* rule 162 */
	 tx1015 = tx935*tx980*tx992; 		/* rule 163 */
	 tx1016 = tx1013 + tx1014 + tx1015; 		/* rule 164 */
	 gy2 = DEImproperRestraintDPhi*tx1016; 		/* rule 165 */
	 fy2 = -gy2; 		/* rule 166 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I2, 1, fy2 );
	 tx1017 = -(tx1006*tx1007*tx1008*tx953); 		/* rule 168 */
	 tx1018 = tx1005*tx1007*tx934*tx992; 		/* rule 169 */
	 tx1019 = tx934*tx980*tx992; 		/* rule 170 */
	 tx1020 = tx1017 + tx1018 + tx1019; 		/* rule 171 */
	 gz2 = DEImproperRestraintDPhi*tx1020; 		/* rule 172 */
	 fz2 = -gz2; 		/* rule 173 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I2, 2, fz2 );
	 tx1021 = -tx1010; 		/* rule 175 */
	 tx1022 = tx1006*tx1007*tx1008*tx955; 		/* rule 176 */
	 tx1023 = -(tx1006*tx955*tx980); 		/* rule 177 */
	 tx1024 = tx1021 + tx1022 + tx1023; 		/* rule 178 */
	 gx3 = DEImproperRestraintDPhi*tx1024; 		/* rule 179 */
	 fx3 = -gx3; 		/* rule 180 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I3, 0, fx3 );
	 tx1025 = -tx1014; 		/* rule 182 */
	 tx1026 = tx1006*tx1007*tx1008*tx954; 		/* rule 183 */
	 tx1027 = -(tx1006*tx954*tx980); 		/* rule 184 */
	 tx1028 = tx1025 + tx1026 + tx1027; 		/* rule 185 */
	 gy3 = DEImproperRestraintDPhi*tx1028; 		/* rule 186 */
	 fy3 = -gy3; 		/* rule 187 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I3, 1, fy3 );
	 tx1029 = -tx1018; 		/* rule 189 */
	 tx1030 = tx1006*tx1007*tx1008*tx953; 		/* rule 190 */
	 tx1031 = -(tx1006*tx953*tx980); 		/* rule 191 */
	 tx1032 = tx1029 + tx1030 + tx1031; 		/* rule 192 */
	 gz3 = DEImproperRestraintDPhi*tx1032; 		/* rule 193 */
	 fz3 = -gz3; 		/* rule 194 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I3, 2, fz3 );
	 gx4 = DEImproperRestraintDPhi*tx1006*tx955*tx980; 		/* rule 196 */
	 fx4 = -gx4; 		/* rule 197 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I4, 0, fx4 );
	 gy4 = DEImproperRestraintDPhi*tx1006*tx954*tx980; 		/* rule 199 */
	 fy4 = -gy4; 		/* rule 200 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I4, 1, fy4 );
	 gz4 = DEImproperRestraintDPhi*tx1006*tx953*tx980; 		/* rule 202 */
	 fz4 = -gz4; 		/* rule 203 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I4, 2, fz4 );
	#ifdef IMPROPER_RESTRAINT_CALC_DIAGONAL_HESSIAN //[
	if (calcDiagonalHessian) {
	 tx1033 = power2(TWOPI); 		/* rule 207 */
	 tx1034 = power2(tx985); 		/* rule 208 */
	 tx1035 = power2(tx991); 		/* rule 209 */
	 tx1036 = power2(y2); 		/* rule 210 */
	 tx1037 = power2(y3); 		/* rule 211 */
	 tx1038 = power2(z2); 		/* rule 212 */
	 tx1039 = power2(z3); 		/* rule 213 */
	 tx1040 = tx1036*tx964; 		/* rule 214 */
	 tx1041 = tx1038*tx964; 		/* rule 215 */
	 tx1042 = tx917*tx965; 		/* rule 216 */
	 tx1043 = 6.*EraseLinearDihedral*K*tx1033*tx1034*tx1035*tx988; 		/* rule 217 */
	 tx1044 = 3*EraseLinearDihedral*K*tx1033*tx1034*tx986*tx990; 		/* rule 218 */
	 tx1045 = tx1037*tx993; 		/* rule 219 */
	 tx1046 = tx1039*tx993; 		/* rule 220 */
	 tx1047 = tx922*tx995; 		/* rule 221 */
	 tx1048 = tx1036*x1; 		/* rule 222 */
	 tx1049 = tx1037*x1; 		/* rule 223 */
	 tx1050 = tx1038*x1; 		/* rule 224 */
	 tx1051 = tx1039*x1; 		/* rule 225 */
	 tx1052 = tx921*y1; 		/* rule 226 */
	 tx1053 = tx916*y2; 		/* rule 227 */
	 tx1054 = tx917*y2; 		/* rule 228 */
	 tx1055 = tx921*y2; 		/* rule 229 */
	 tx1056 = -2.*tx918*y3; 		/* rule 230 */
	 tx1057 = x3*y2*y3; 		/* rule 231 */
	 tx1058 = tx927*z1; 		/* rule 232 */
	 tx1059 = -2.*tx930*z2; 		/* rule 233 */
	 tx1060 = tx922*z3; 		/* rule 234 */
	 tx1061 = tx923*z3; 		/* rule 235 */
	 tx1062 = tx927*z3; 		/* rule 236 */
	 tx1063 = x2*z2*z3; 		/* rule 237 */
	 tx1064 = tx1043 + tx1044; 		/* rule 238 */
	 tx1065 = tx1040 + tx1041 + tx1042 + tx1045 + tx1046 + tx1047 + tx1048 + tx1049 + tx1050 + tx1051 + tx1052 + tx1053 + tx1054 + tx1055 + tx1056 + tx1057 + tx1058 + tx1059 + tx1060 + tx1061 + tx1062 + tx1063; 		/* rule 239 */
	 tx1066 = power2(tx992); 		/* rule 240 */
	 tx1067 = tx1064*tx1066*tx939*tx976; 		/* rule 241 */
	 tx1068 = 6.*EraseLinearDihedral*K*TWOPI*tx1065*tx1066*tx936*tx980*tx985*tx990*tx991; 		/* rule 242 */
	 dhx1x1 = tx1067 + tx1068; 		/* rule 243 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx1069 = power2(x2); 		/* rule 245 */
	 tx1070 = power2(x3); 		/* rule 246 */
	 tx1071 = tx1038*tx965; 		/* rule 247 */
	 tx1072 = tx1069*tx965; 		/* rule 248 */
	 tx1073 = tx925*tx966; 		/* rule 249 */
	 tx1074 = tx918*tx993; 		/* rule 250 */
	 tx1075 = tx1039*tx994; 		/* rule 251 */
	 tx1076 = tx1070*tx994; 		/* rule 252 */
	 tx1077 = tx921*x1; 		/* rule 253 */
	 tx1078 = -2.*tx917*x2; 		/* rule 254 */
	 tx1079 = tx918*x3; 		/* rule 255 */
	 tx1080 = tx920*x3; 		/* rule 256 */
	 tx1081 = tx921*x3; 		/* rule 257 */
	 tx1082 = tx1038*y1; 		/* rule 258 */
	 tx1083 = tx1039*y1; 		/* rule 259 */
	 tx1084 = tx1069*y1; 		/* rule 260 */
	 tx1085 = tx1070*y1; 		/* rule 261 */
	 tx1086 = x2*x3*y2; 		/* rule 262 */
	 tx1087 = tx933*z1; 		/* rule 263 */
	 tx1088 = tx924*z2; 		/* rule 264 */
	 tx1089 = tx925*z2; 		/* rule 265 */
	 tx1090 = tx933*z2; 		/* rule 266 */
	 tx1091 = -2.*tx928*z3; 		/* rule 267 */
	 tx1092 = y3*z2*z3; 		/* rule 268 */
	 tx1093 = tx1071 + tx1072 + tx1073 + tx1074 + tx1075 + tx1076 + tx1077 + tx1078 + tx1079 + tx1080 + tx1081 + tx1082 + tx1083 + tx1084 + tx1085 + tx1086 + tx1087 + tx1088 + tx1089 + tx1090 + tx1091 + tx1092; 		/* rule 269 */
	 tx1094 = tx1064*tx1066*tx938*tx976; 		/* rule 270 */
	 tx1095 = 6.*EraseLinearDihedral*K*TWOPI*tx1066*tx1093*tx935*tx980*tx985*tx990*tx991; 		/* rule 271 */
	 dhy1y1 = tx1094 + tx1095; 		/* rule 272 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx1096 = tx930*tx964; 		/* rule 274 */
	 tx1097 = tx1036*tx966; 		/* rule 275 */
	 tx1098 = tx1069*tx966; 		/* rule 276 */
	 tx1099 = tx928*tx994; 		/* rule 277 */
	 tx1100 = tx1037*tx995; 		/* rule 278 */
	 tx1101 = tx1070*tx995; 		/* rule 279 */
	 tx1102 = tx927*x1; 		/* rule 280 */
	 tx1103 = tx926*x2; 		/* rule 281 */
	 tx1104 = tx927*x2; 		/* rule 282 */
	 tx1105 = tx930*x2; 		/* rule 283 */
	 tx1106 = -2.*tx922*x3; 		/* rule 284 */
	 tx1107 = tx933*y1; 		/* rule 285 */
	 tx1108 = -2.*tx925*y2; 		/* rule 286 */
	 tx1109 = tx928*y3; 		/* rule 287 */
	 tx1110 = tx932*y3; 		/* rule 288 */
	 tx1111 = tx933*y3; 		/* rule 289 */
	 tx1112 = tx1036*z1; 		/* rule 290 */
	 tx1113 = tx1037*z1; 		/* rule 291 */
	 tx1114 = tx1069*z1; 		/* rule 292 */
	 tx1115 = tx1070*z1; 		/* rule 293 */
	 tx1116 = y2*y3*z2; 		/* rule 294 */
	 tx1117 = x2*x3*z3; 		/* rule 295 */
	 tx1118 = tx1096 + tx1097 + tx1098 + tx1099 + tx1100 + tx1101 + tx1102 + tx1103 + tx1104 + tx1105 + tx1106 + tx1107 + tx1108 + tx1109 + tx1110 + tx1111 + tx1112 + tx1113 + tx1114 + tx1115 + tx1116 + tx1117; 		/* rule 296 */
	 tx1119 = tx1064*tx1066*tx937*tx976; 		/* rule 297 */
	 tx1120 = 6.*EraseLinearDihedral*K*TWOPI*tx1066*tx1118*tx934*tx980*tx985*tx990*tx991; 		/* rule 298 */
	 dhz1z1 = tx1119 + tx1120; 		/* rule 299 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 tx1121 = power2(y1); 		/* rule 301 */
	 tx1122 = power2(z1); 		/* rule 302 */
	 tx1123 = -tx1048; 		/* rule 303 */
	 tx1124 = -tx1049; 		/* rule 304 */
	 tx1125 = -tx1050; 		/* rule 305 */
	 tx1126 = -tx1051; 		/* rule 306 */
	 tx1127 = 2.*tx1054; 		/* rule 307 */
	 tx1128 = 2.*tx1058; 		/* rule 308 */
	 tx1129 = tx1121*tx964; 		/* rule 309 */
	 tx1130 = tx1122*tx964; 		/* rule 310 */
	 tx1131 = tx918*tx965; 		/* rule 311 */
	 tx1132 = tx922*tx966; 		/* rule 312 */
	 tx1133 = tx927*tx966; 		/* rule 313 */
	 tx1134 = tx917*tx994; 		/* rule 314 */
	 tx1135 = tx921*tx994; 		/* rule 315 */
	 tx1136 = tx1037*x2; 		/* rule 316 */
	 tx1137 = tx1039*x2; 		/* rule 317 */
	 tx1138 = tx1121*x2; 		/* rule 318 */
	 tx1139 = tx1122*x2; 		/* rule 319 */
	 tx1140 = tx1036*x3; 		/* rule 320 */
	 tx1141 = tx1038*x3; 		/* rule 321 */
	 tx1142 = -(tx918*y1); 		/* rule 322 */
	 tx1143 = x2*y1*y2; 		/* rule 323 */
	 tx1144 = tx916*y3; 		/* rule 324 */
	 tx1145 = tx917*y3; 		/* rule 325 */
	 tx1146 = 2.*tx918*y3; 		/* rule 326 */
	 tx1147 = tx919*y3; 		/* rule 327 */
	 tx1148 = x1*y1*y3; 		/* rule 328 */
	 tx1149 = tx926*z1; 		/* rule 329 */
	 tx1150 = tx930*z1; 		/* rule 330 */
	 tx1151 = tx922*z2; 		/* rule 331 */
	 tx1152 = tx923*z2; 		/* rule 332 */
	 tx1153 = 2.*tx930*z2; 		/* rule 333 */
	 tx1154 = tx931*z2; 		/* rule 334 */
	 tx1155 = tx926*z3; 		/* rule 335 */
	 tx1156 = x3*z1*z3; 		/* rule 336 */
	 tx1157 = tx1040 + tx1041 + tx1047 + tx1048 + tx1050 + tx1053 + tx1055 + tx1063 + tx1127 + tx1128 + tx1129 + tx1130 + tx1131 + tx1132 + tx1138 + tx1139 + tx1142 + tx1144 + tx1148 + tx1149 + tx1150 + tx1155; 		/* rule 337 */
	 tx1158 = tx1123 + tx1124 + tx1125 + tx1126 + tx1132 + tx1133 + tx1134 + tx1135 + tx1136 + tx1137 + tx1140 + tx1141 + tx1143 + tx1144 + tx1145 + tx1146 + tx1147 + tx1151 + tx1152 + tx1153 + tx1154 + tx1156; 		/* rule 338 */
	 tx1159 = power2(y4); 		/* rule 339 */
	 tx1160 = power2(z4); 		/* rule 340 */
	 tx1161 = -2.*tx1054; 		/* rule 341 */
	 tx1162 = -2.*tx1058; 		/* rule 342 */
	 tx1163 = -tx1150; 		/* rule 343 */
	 tx1164 = tx1005*tx1158*tx936; 		/* rule 344 */
	 tx1165 = tx941*tx965; 		/* rule 345 */
	 tx1166 = tx944*tx965; 		/* rule 346 */
	 tx1167 = tx1157*tx936*tx976; 		/* rule 347 */
	 tx1168 = tx1121*tx993; 		/* rule 348 */
	 tx1169 = tx1122*tx993; 		/* rule 349 */
	 tx1170 = tx1159*tx993; 		/* rule 350 */
	 tx1171 = tx1160*tx993; 		/* rule 351 */
	 tx1172 = tx949*tx995; 		/* rule 352 */
	 tx1173 = tx1121*x3; 		/* rule 353 */
	 tx1174 = tx1122*x3; 		/* rule 354 */
	 tx1175 = tx1159*x3; 		/* rule 355 */
	 tx1176 = tx1160*x3; 		/* rule 356 */
	 tx1177 = tx1036*x4; 		/* rule 357 */
	 tx1178 = tx1037*x4; 		/* rule 358 */
	 tx1179 = tx1038*x4; 		/* rule 359 */
	 tx1180 = tx1039*x4; 		/* rule 360 */
	 tx1181 = tx918*y1; 		/* rule 361 */
	 tx1182 = tx920*y1; 		/* rule 362 */
	 tx1183 = tx943*y2; 		/* rule 363 */
	 tx1184 = tx944*y2; 		/* rule 364 */
	 tx1185 = tx918*y3; 		/* rule 365 */
	 tx1186 = -2.*tx941*y3; 		/* rule 366 */
	 tx1187 = tx919*y4; 		/* rule 367 */
	 tx1188 = tx921*y4; 		/* rule 368 */
	 tx1189 = 2.*tx1188; 		/* rule 369 */
	 tx1190 = tx941*y4; 		/* rule 370 */
	 tx1191 = tx942*y4; 		/* rule 371 */
	 tx1192 = tx930*z2; 		/* rule 372 */
	 tx1193 = -2.*tx947*z2; 		/* rule 373 */
	 tx1194 = x1*z1*z2; 		/* rule 374 */
	 tx1195 = tx945*z3; 		/* rule 375 */
	 tx1196 = tx949*z3; 		/* rule 376 */
	 tx1197 = 2.*tx1196; 		/* rule 377 */
	 tx1198 = tx950*z3; 		/* rule 378 */
	 tx1199 = -(tx927*z4); 		/* rule 379 */
	 tx1200 = tx927*z4; 		/* rule 380 */
	 tx1201 = -(tx947*z4); 		/* rule 381 */
	 tx1202 = x4*z2*z4; 		/* rule 382 */
	 tx1203 = tx1189; 		/* rule 383 */
	 tx1204 = tx1197; 		/* rule 384 */
	 tx1205 = tx1164 + tx1167; 		/* rule 385 */
	 tx1206 = tx1052 + tx1060 + tx1123 + tx1125 + tx1135 + tx1140 + tx1141 + tx1143 + tx1151 + tx1154 + tx1161 + tx1162 + tx1163 + tx1168 + tx1169 + tx1173 + tx1174 + tx1181 + tx1182 + tx1185 + tx1192 + tx1194; 		/* rule 386 */
	 tx1207 = tx1040 + tx1041 + tx1045 + tx1046 + tx1055 + tx1057 + tx1062 + tx1063 + tx1166 + tx1172 + tx1177 + tx1178 + tx1179 + tx1180 + tx1183 + tx1184 + tx1186 + tx1188 + tx1193 + tx1196 + tx1198 + tx1200; 		/* rule 387 */
	 tx1208 = tx1045 + tx1046 + tx1057 + tx1062 + tx1165 + tx1166 + tx1170 + tx1171 + tx1175 + tx1176 + tx1178 + tx1180 + tx1187 + tx1190 + tx1191 + tx1195 + tx1198 + tx1199 + tx1201 + tx1202 + tx1203 + tx1204; 		/* rule 388 */
	 tx1209 = power2(tx1006); 		/* rule 389 */
	 tx1832 = tx976; 		/* rule 390 */
	 tx1833 = reciprocal(tx1832); 		/* rule 391 */
	 tx1834 = tx1007; 		/* rule 392 */
	 tx1210 = tx1833*tx1834; 		/* rule 393 */
	 tx1211 = -2.*tx1066*tx1205*tx1834; 		/* rule 394 */
	 tx1212 = 2.*tx1005*tx1066*tx1206*tx1834*tx936; 		/* rule 395 */
	 tx1213 = -2.*tx1008*tx1208*tx1209*tx1834*tx955; 		/* rule 396 */
	 tx1214 = -(tx1006*tx1207*tx1210*tx955); 		/* rule 397 */
	 tx1215 = 2.*tx1065*tx1066*tx936*tx980; 		/* rule 398 */
	 tx1216 = tx1065*tx1210*tx936*tx992; 		/* rule 399 */
	 tx1217 = power2(tx1012); 		/* rule 400 */
	 tx1218 = tx1211 + tx1212 + tx1213 + tx1214 + tx1215 + tx1216; 		/* rule 401 */
	 tx1219 = tx1064*tx1217; 		/* rule 402 */
	 tx1220 = DEImproperRestraintDPhi*tx1218; 		/* rule 403 */
	 dhx2x2 = tx1219 + tx1220; 		/* rule 404 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 tx1221 = power2(x1); 		/* rule 406 */
	 tx1222 = 2.*tx1077; 		/* rule 407 */
	 tx1223 = -tx1082; 		/* rule 408 */
	 tx1224 = -tx1083; 		/* rule 409 */
	 tx1225 = -tx1084; 		/* rule 410 */
	 tx1226 = -tx1085; 		/* rule 411 */
	 tx1227 = 2.*tx1089; 		/* rule 412 */
	 tx1228 = tx918*tx964; 		/* rule 413 */
	 tx1229 = tx921*tx964; 		/* rule 414 */
	 tx1230 = tx1122*tx965; 		/* rule 415 */
	 tx1231 = tx1221*tx965; 		/* rule 416 */
	 tx1232 = tx928*tx966; 		/* rule 417 */
	 tx1233 = tx925*tx995; 		/* rule 418 */
	 tx1234 = tx933*tx995; 		/* rule 419 */
	 tx1235 = tx916*x1; 		/* rule 420 */
	 tx1236 = tx917*x1; 		/* rule 421 */
	 tx1237 = 2.*tx917*x2; 		/* rule 422 */
	 tx1238 = tx918*x2; 		/* rule 423 */
	 tx1239 = tx919*x2; 		/* rule 424 */
	 tx1240 = tx920*x2; 		/* rule 425 */
	 tx1241 = tx916*x3; 		/* rule 426 */
	 tx1242 = tx1039*y2; 		/* rule 427 */
	 tx1243 = tx1070*y2; 		/* rule 428 */
	 tx1244 = tx1122*y2; 		/* rule 429 */
	 tx1245 = tx1221*y2; 		/* rule 430 */
	 tx1246 = tx1038*y3; 		/* rule 431 */
	 tx1247 = tx1069*y3; 		/* rule 432 */
	 tx1248 = x1*x3*y3; 		/* rule 433 */
	 tx1249 = -(tx928*z1); 		/* rule 434 */
	 tx1250 = y2*z1*z2; 		/* rule 435 */
	 tx1251 = tx924*z3; 		/* rule 436 */
	 tx1252 = tx925*z3; 		/* rule 437 */
	 tx1253 = 2.*tx928*z3; 		/* rule 438 */
	 tx1254 = tx929*z3; 		/* rule 439 */
	 tx1255 = y1*z1*z3; 		/* rule 440 */
	 tx1256 = tx1223 + tx1224 + tx1225 + tx1226 + tx1228 + tx1229 + tx1233 + tx1234 + tx1237 + tx1238 + tx1239 + tx1240 + tx1242 + tx1243 + tx1246 + tx1247 + tx1248 + tx1250 + tx1251 + tx1252 + tx1253 + tx1254; 		/* rule 441 */
	 tx1257 = tx1071 + tx1072 + tx1074 + tx1082 + tx1084 + tx1086 + tx1088 + tx1090 + tx1222 + tx1227 + tx1228 + tx1230 + tx1231 + tx1232 + tx1235 + tx1236 + tx1241 + tx1244 + tx1245 + tx1249 + tx1251 + tx1255; 		/* rule 442 */
	 tx1258 = power2(x4); 		/* rule 443 */
	 tx1259 = -2.*tx1077; 		/* rule 444 */
	 tx1260 = -2.*tx1089; 		/* rule 445 */
	 tx1261 = -tx1236; 		/* rule 446 */
	 tx1262 = tx1005*tx1256*tx935; 		/* rule 447 */
	 tx1263 = tx946*tx966; 		/* rule 448 */
	 tx1264 = tx952*tx966; 		/* rule 449 */
	 tx1265 = tx1257*tx1832*tx935; 		/* rule 450 */
	 tx1266 = tx941*tx993; 		/* rule 451 */
	 tx1267 = tx1122*tx994; 		/* rule 452 */
	 tx1268 = tx1160*tx994; 		/* rule 453 */
	 tx1269 = tx1221*tx994; 		/* rule 454 */
	 tx1270 = tx1258*tx994; 		/* rule 455 */
	 tx1271 = tx917*x2; 		/* rule 456 */
	 tx1272 = -2.*tx944*x2; 		/* rule 457 */
	 tx1273 = tx941*x3; 		/* rule 458 */
	 tx1274 = 2.*tx1273; 		/* rule 459 */
	 tx1275 = tx942*x3; 		/* rule 460 */
	 tx1276 = tx943*x3; 		/* rule 461 */
	 tx1277 = -(tx921*x4); 		/* rule 462 */
	 tx1278 = tx921*x4; 		/* rule 463 */
	 tx1279 = -(tx944*x4); 		/* rule 464 */
	 tx1280 = x1*x2*y1; 		/* rule 465 */
	 tx1281 = tx1122*y3; 		/* rule 466 */
	 tx1282 = tx1160*y3; 		/* rule 467 */
	 tx1283 = tx1221*y3; 		/* rule 468 */
	 tx1284 = tx1258*y3; 		/* rule 469 */
	 tx1285 = tx1038*y4; 		/* rule 470 */
	 tx1286 = tx1039*y4; 		/* rule 471 */
	 tx1287 = tx1069*y4; 		/* rule 472 */
	 tx1288 = tx1070*y4; 		/* rule 473 */
	 tx1289 = x2*x4*y4; 		/* rule 474 */
	 tx1290 = tx928*z1; 		/* rule 475 */
	 tx1291 = tx932*z1; 		/* rule 476 */
	 tx1292 = tx951*z2; 		/* rule 477 */
	 tx1293 = tx952*z2; 		/* rule 478 */
	 tx1294 = tx928*z3; 		/* rule 479 */
	 tx1295 = -2.*tx946*z3; 		/* rule 480 */
	 tx1296 = tx929*z4; 		/* rule 481 */
	 tx1297 = tx933*z4; 		/* rule 482 */
	 tx1298 = 2.*tx1297; 		/* rule 483 */
	 tx1299 = tx946*z4; 		/* rule 484 */
	 tx1300 = tx948*z4; 		/* rule 485 */
	 tx1301 = tx1274; 		/* rule 486 */
	 tx1302 = tx1298; 		/* rule 487 */
	 tx1303 = tx1262 + tx1265; 		/* rule 488 */
	 tx1304 = tx1079 + tx1087 + tx1223 + tx1225 + tx1234 + tx1238 + tx1239 + tx1246 + tx1247 + tx1250 + tx1259 + tx1260 + tx1261 + tx1267 + tx1269 + tx1271 + tx1280 + tx1281 + tx1283 + tx1290 + tx1291 + tx1294; 		/* rule 489 */
	 tx1305 = tx1071 + tx1072 + tx1075 + tx1076 + tx1081 + tx1086 + tx1090 + tx1092 + tx1264 + tx1266 + tx1272 + tx1273 + tx1275 + tx1278 + tx1285 + tx1286 + tx1287 + tx1288 + tx1292 + tx1293 + tx1295 + tx1297; 		/* rule 490 */
	 tx1306 = tx1075 + tx1076 + tx1081 + tx1092 + tx1263 + tx1264 + tx1268 + tx1270 + tx1275 + tx1276 + tx1277 + tx1279 + tx1282 + tx1284 + tx1286 + tx1288 + tx1289 + tx1296 + tx1299 + tx1300 + tx1301 + tx1302; 		/* rule 491 */
	 tx1307 = -2.*tx1066*tx1303*tx1834; 		/* rule 492 */
	 tx1308 = 2.*tx1005*tx1066*tx1304*tx1834*tx935; 		/* rule 493 */
	 tx1309 = -(tx1006*tx1210*tx1305*tx954); 		/* rule 494 */
	 tx1310 = -2.*tx1008*tx1209*tx1306*tx1834*tx954; 		/* rule 495 */
	 tx1311 = 2.*tx1066*tx1093*tx935*tx980; 		/* rule 496 */
	 tx1312 = tx1093*tx1210*tx935*tx992; 		/* rule 497 */
	 tx1313 = power2(tx1016); 		/* rule 498 */
	 tx1314 = tx1307 + tx1308 + tx1309 + tx1310 + tx1311 + tx1312; 		/* rule 499 */
	 tx1315 = tx1064*tx1313; 		/* rule 500 */
	 tx1316 = DEImproperRestraintDPhi*tx1314; 		/* rule 501 */
	 dhy2y2 = tx1315 + tx1316; 		/* rule 502 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 tx1317 = 2.*tx1105; 		/* rule 504 */
	 tx1318 = 2.*tx1107; 		/* rule 505 */
	 tx1319 = -tx1112; 		/* rule 506 */
	 tx1320 = -tx1113; 		/* rule 507 */
	 tx1321 = -tx1114; 		/* rule 508 */
	 tx1322 = -tx1115; 		/* rule 509 */
	 tx1323 = tx922*tx964; 		/* rule 510 */
	 tx1324 = tx928*tx965; 		/* rule 511 */
	 tx1325 = tx933*tx965; 		/* rule 512 */
	 tx1326 = tx1121*tx966; 		/* rule 513 */
	 tx1327 = tx1221*tx966; 		/* rule 514 */
	 tx1328 = tx927*tx993; 		/* rule 515 */
	 tx1329 = tx930*tx993; 		/* rule 516 */
	 tx1330 = -(tx922*x1); 		/* rule 517 */
	 tx1331 = 2.*tx922*x3; 		/* rule 518 */
	 tx1332 = tx926*x3; 		/* rule 519 */
	 tx1333 = tx930*x3; 		/* rule 520 */
	 tx1334 = tx931*x3; 		/* rule 521 */
	 tx1335 = tx924*y1; 		/* rule 522 */
	 tx1336 = tx925*y1; 		/* rule 523 */
	 tx1337 = 2.*tx925*y2; 		/* rule 524 */
	 tx1338 = tx928*y2; 		/* rule 525 */
	 tx1339 = tx929*y2; 		/* rule 526 */
	 tx1340 = tx932*y2; 		/* rule 527 */
	 tx1341 = tx924*y3; 		/* rule 528 */
	 tx1342 = x1*x3*z1; 		/* rule 529 */
	 tx1343 = tx1037*z2; 		/* rule 530 */
	 tx1344 = tx1070*z2; 		/* rule 531 */
	 tx1345 = tx1121*z2; 		/* rule 532 */
	 tx1346 = tx1221*z2; 		/* rule 533 */
	 tx1347 = x1*x2*z2; 		/* rule 534 */
	 tx1348 = tx1036*z3; 		/* rule 535 */
	 tx1349 = tx1069*z3; 		/* rule 536 */
	 tx1350 = y1*y3*z3; 		/* rule 537 */
	 tx1351 = tx1097 + tx1098 + tx1099 + tx1103 + tx1104 + tx1112 + tx1114 + tx1116 + tx1317 + tx1318 + tx1323 + tx1324 + tx1326 + tx1327 + tx1330 + tx1332 + tx1335 + tx1336 + tx1341 + tx1342 + tx1345 + tx1346; 		/* rule 538 */
	 tx1352 = tx1319 + tx1320 + tx1321 + tx1322 + tx1324 + tx1325 + tx1328 + tx1329 + tx1331 + tx1332 + tx1333 + tx1334 + tx1337 + tx1338 + tx1339 + tx1340 + tx1343 + tx1344 + tx1347 + tx1348 + tx1349 + tx1350; 		/* rule 539 */
	 tx1353 = -2.*tx1105; 		/* rule 540 */
	 tx1354 = -2.*tx1107; 		/* rule 541 */
	 tx1355 = -tx1336; 		/* rule 542 */
	 tx1356 = tx1005*tx1352*tx934; 		/* rule 543 */
	 tx1357 = tx947*tx964; 		/* rule 544 */
	 tx1358 = tx949*tx964; 		/* rule 545 */
	 tx1359 = tx1351*tx1832*tx934; 		/* rule 546 */
	 tx1360 = tx946*tx994; 		/* rule 547 */
	 tx1361 = tx1121*tx995; 		/* rule 548 */
	 tx1362 = tx1159*tx995; 		/* rule 549 */
	 tx1363 = tx1221*tx995; 		/* rule 550 */
	 tx1364 = tx1258*tx995; 		/* rule 551 */
	 tx1365 = tx922*x1; 		/* rule 552 */
	 tx1366 = tx923*x1; 		/* rule 553 */
	 tx1367 = tx945*x2; 		/* rule 554 */
	 tx1368 = tx947*x2; 		/* rule 555 */
	 tx1369 = tx922*x3; 		/* rule 556 */
	 tx1370 = -2.*tx949*x3; 		/* rule 557 */
	 tx1371 = tx927*x4; 		/* rule 558 */
	 tx1372 = 2.*tx1371; 		/* rule 559 */
	 tx1373 = tx931*x4; 		/* rule 560 */
	 tx1374 = tx949*x4; 		/* rule 561 */
	 tx1375 = tx950*x4; 		/* rule 562 */
	 tx1376 = tx925*y2; 		/* rule 563 */
	 tx1377 = -2.*tx952*y2; 		/* rule 564 */
	 tx1378 = tx946*y3; 		/* rule 565 */
	 tx1379 = 2.*tx1378; 		/* rule 566 */
	 tx1380 = tx948*y3; 		/* rule 567 */
	 tx1381 = tx951*y3; 		/* rule 568 */
	 tx1382 = -(tx933*y4); 		/* rule 569 */
	 tx1383 = tx933*y4; 		/* rule 570 */
	 tx1384 = -(tx952*y4); 		/* rule 571 */
	 tx1385 = y1*y2*z1; 		/* rule 572 */
	 tx1386 = tx1121*z3; 		/* rule 573 */
	 tx1387 = tx1159*z3; 		/* rule 574 */
	 tx1388 = tx1221*z3; 		/* rule 575 */
	 tx1389 = tx1258*z3; 		/* rule 576 */
	 tx1390 = tx1036*z4; 		/* rule 577 */
	 tx1391 = tx1037*z4; 		/* rule 578 */
	 tx1392 = tx1069*z4; 		/* rule 579 */
	 tx1393 = tx1070*z4; 		/* rule 580 */
	 tx1394 = y2*y4*z4; 		/* rule 581 */
	 tx1395 = tx1372; 		/* rule 582 */
	 tx1396 = tx1379; 		/* rule 583 */
	 tx1397 = tx1356 + tx1359; 		/* rule 584 */
	 tx1398 = tx1102 + tx1109 + tx1319 + tx1321 + tx1328 + tx1338 + tx1339 + tx1347 + tx1348 + tx1349 + tx1353 + tx1354 + tx1355 + tx1361 + tx1363 + tx1365 + tx1366 + tx1369 + tx1376 + tx1385 + tx1386 + tx1388; 		/* rule 585 */
	 tx1399 = tx1097 + tx1098 + tx1100 + tx1101 + tx1104 + tx1111 + tx1116 + tx1117 + tx1357 + tx1360 + tx1367 + tx1368 + tx1370 + tx1371 + tx1377 + tx1378 + tx1380 + tx1383 + tx1390 + tx1391 + tx1392 + tx1393; 		/* rule 586 */
	 tx1400 = tx1100 + tx1101 + tx1111 + tx1117 + tx1357 + tx1358 + tx1362 + tx1364 + tx1373 + tx1374 + tx1375 + tx1380 + tx1381 + tx1382 + tx1384 + tx1387 + tx1389 + tx1391 + tx1393 + tx1394 + tx1395 + tx1396; 		/* rule 587 */
	 tx1401 = -2.*tx1066*tx1397*tx1834; 		/* rule 588 */
	 tx1402 = 2.*tx1005*tx1066*tx1398*tx1834*tx934; 		/* rule 589 */
	 tx1403 = -(tx1006*tx1210*tx1399*tx953); 		/* rule 590 */
	 tx1404 = -2.*tx1008*tx1209*tx1400*tx1834*tx953; 		/* rule 591 */
	 tx1405 = 2.*tx1066*tx1118*tx934*tx980; 		/* rule 592 */
	 tx1406 = tx1118*tx1210*tx934*tx992; 		/* rule 593 */
	 tx1407 = power2(tx1020); 		/* rule 594 */
	 tx1408 = tx1401 + tx1402 + tx1403 + tx1404 + tx1405 + tx1406; 		/* rule 595 */
	 tx1409 = tx1064*tx1407; 		/* rule 596 */
	 tx1410 = DEImproperRestraintDPhi*tx1408; 		/* rule 597 */
	 dhz2z2 = tx1409 + tx1410; 		/* rule 598 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	 tx1411 = -tx1177; 		/* rule 600 */
	 tx1412 = -tx1178; 		/* rule 601 */
	 tx1413 = -tx1179; 		/* rule 602 */
	 tx1414 = -tx1180; 		/* rule 603 */
	 tx1415 = -2.*tx1188; 		/* rule 604 */
	 tx1416 = -tx1188; 		/* rule 605 */
	 tx1417 = -tx1190; 		/* rule 606 */
	 tx1418 = -2.*tx1196; 		/* rule 607 */
	 tx1419 = tx1159*tx964; 		/* rule 608 */
	 tx1420 = tx1160*tx964; 		/* rule 609 */
	 tx1421 = tx1159*x2; 		/* rule 610 */
	 tx1422 = tx1160*x2; 		/* rule 611 */
	 tx1423 = tx941*y3; 		/* rule 612 */
	 tx1424 = 2.*tx1423; 		/* rule 613 */
	 tx1425 = tx944*y3; 		/* rule 614 */
	 tx1426 = x2*y2*y4; 		/* rule 615 */
	 tx1427 = x4*y3*y4; 		/* rule 616 */
	 tx1428 = tx947*z2; 		/* rule 617 */
	 tx1429 = 2.*tx1428; 		/* rule 618 */
	 tx1430 = tx949*z2; 		/* rule 619 */
	 tx1431 = tx931*z4; 		/* rule 620 */
	 tx1432 = tx945*z4; 		/* rule 621 */
	 tx1433 = tx947*z4; 		/* rule 622 */
	 tx1434 = x3*z3*z4; 		/* rule 623 */
	 tx1435 = tx1424; 		/* rule 624 */
	 tx1436 = tx1429; 		/* rule 625 */
	 tx1437 = tx1133 + tx1136 + tx1137 + tx1147 + tx1184 + tx1200 + tx1412 + tx1414 + tx1415 + tx1417 + tx1418 + tx1419 + tx1420 + tx1421 + tx1422 + tx1423 + tx1425 + tx1427 + tx1428 + tx1432 + tx1433 + tx1434; 		/* rule 626 */
	 tx1438 = tx1133 + tx1135 + tx1136 + tx1137 + tx1140 + tx1141 + tx1147 + tx1154 + tx1187 + tx1199 + tx1411 + tx1412 + tx1413 + tx1414 + tx1416 + tx1425 + tx1426 + tx1430 + tx1431 + tx1434 + tx1435 + tx1436; 		/* rule 627 */
	 tx1439 = tx1008*tx1438*tx955; 		/* rule 628 */
	 tx1440 = tx1437*tx1832*tx955; 		/* rule 629 */
	 tx1441 = tx1439 + tx1440; 		/* rule 630 */
	 tx1442 = -2.*tx1209*tx1441*tx1834; 		/* rule 631 */
	 tx1443 = -2.*tx1207*tx1209*tx955*tx980; 		/* rule 632 */
	 tx1444 = power2(tx1024); 		/* rule 633 */
	 tx1445 = tx1212 + tx1213 + tx1214 + tx1216 + tx1442 + tx1443; 		/* rule 634 */
	 tx1446 = tx1064*tx1444; 		/* rule 635 */
	 tx1447 = DEImproperRestraintDPhi*tx1445; 		/* rule 636 */
	 dhx3x3 = tx1446 + tx1447; 		/* rule 637 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 0, dhx3x3);
	 tx1448 = -2.*tx1273; 		/* rule 639 */
	 tx1449 = -tx1285; 		/* rule 640 */
	 tx1450 = -tx1286; 		/* rule 641 */
	 tx1451 = -tx1287; 		/* rule 642 */
	 tx1452 = -tx1288; 		/* rule 643 */
	 tx1453 = -2.*tx1297; 		/* rule 644 */
	 tx1454 = -tx1297; 		/* rule 645 */
	 tx1455 = -tx1299; 		/* rule 646 */
	 tx1456 = tx1160*tx965; 		/* rule 647 */
	 tx1457 = tx1258*tx965; 		/* rule 648 */
	 tx1458 = tx941*x2; 		/* rule 649 */
	 tx1459 = tx944*x2; 		/* rule 650 */
	 tx1460 = 2.*tx1459; 		/* rule 651 */
	 tx1461 = tx919*x4; 		/* rule 652 */
	 tx1462 = tx943*x4; 		/* rule 653 */
	 tx1463 = tx944*x4; 		/* rule 654 */
	 tx1464 = tx1160*y2; 		/* rule 655 */
	 tx1465 = tx1258*y2; 		/* rule 656 */
	 tx1466 = x3*x4*y3; 		/* rule 657 */
	 tx1467 = tx946*z3; 		/* rule 658 */
	 tx1468 = 2.*tx1467; 		/* rule 659 */
	 tx1469 = tx952*z3; 		/* rule 660 */
	 tx1470 = y2*z2*z4; 		/* rule 661 */
	 tx1471 = y4*z3*z4; 		/* rule 662 */
	 tx1472 = tx1460; 		/* rule 663 */
	 tx1473 = tx1468; 		/* rule 664 */
	 tx1474 = tx1229 + tx1242 + tx1243 + tx1254 + tx1278 + tx1293 + tx1448 + tx1450 + tx1452 + tx1453 + tx1455 + tx1456 + tx1457 + tx1459 + tx1462 + tx1463 + tx1464 + tx1465 + tx1466 + tx1467 + tx1469 + tx1471; 		/* rule 665 */
	 tx1475 = tx1229 + tx1234 + tx1239 + tx1242 + tx1243 + tx1246 + tx1247 + tx1254 + tx1277 + tx1296 + tx1449 + tx1450 + tx1451 + tx1452 + tx1454 + tx1458 + tx1461 + tx1466 + tx1469 + tx1470 + tx1472 + tx1473; 		/* rule 666 */
	 tx1476 = tx1008*tx1475*tx954; 		/* rule 667 */
	 tx1477 = tx1474*tx1832*tx954; 		/* rule 668 */
	 tx1478 = tx1476 + tx1477; 		/* rule 669 */
	 tx1479 = -2.*tx1209*tx1478*tx1834; 		/* rule 670 */
	 tx1480 = -2.*tx1209*tx1305*tx954*tx980; 		/* rule 671 */
	 tx1481 = power2(tx1028); 		/* rule 672 */
	 tx1482 = tx1308 + tx1309 + tx1310 + tx1312 + tx1479 + tx1480; 		/* rule 673 */
	 tx1483 = tx1064*tx1481; 		/* rule 674 */
	 tx1484 = DEImproperRestraintDPhi*tx1482; 		/* rule 675 */
	 dhy3y3 = tx1483 + tx1484; 		/* rule 676 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 1, dhy3y3);
	 tx1485 = -2.*tx1371; 		/* rule 678 */
	 tx1486 = -tx1371; 		/* rule 679 */
	 tx1487 = -tx1374; 		/* rule 680 */
	 tx1488 = -2.*tx1378; 		/* rule 681 */
	 tx1489 = -tx1390; 		/* rule 682 */
	 tx1490 = -tx1391; 		/* rule 683 */
	 tx1491 = -tx1392; 		/* rule 684 */
	 tx1492 = -tx1393; 		/* rule 685 */
	 tx1493 = tx1159*tx966; 		/* rule 686 */
	 tx1494 = tx1258*tx966; 		/* rule 687 */
	 tx1495 = tx947*x3; 		/* rule 688 */
	 tx1496 = tx949*x3; 		/* rule 689 */
	 tx1497 = 2.*tx1496; 		/* rule 690 */
	 tx1498 = tx946*y2; 		/* rule 691 */
	 tx1499 = tx952*y2; 		/* rule 692 */
	 tx1500 = 2.*tx1499; 		/* rule 693 */
	 tx1501 = tx929*y4; 		/* rule 694 */
	 tx1502 = tx951*y4; 		/* rule 695 */
	 tx1503 = tx952*y4; 		/* rule 696 */
	 tx1504 = tx1159*z2; 		/* rule 697 */
	 tx1505 = tx1258*z2; 		/* rule 698 */
	 tx1506 = x2*x4*z2; 		/* rule 699 */
	 tx1507 = y3*y4*z3; 		/* rule 700 */
	 tx1508 = x3*x4*z4; 		/* rule 701 */
	 tx1509 = tx1497; 		/* rule 702 */
	 tx1510 = tx1500; 		/* rule 703 */
	 tx1511 = tx1325 + tx1334 + tx1343 + tx1344 + tx1368 + tx1383 + tx1485 + tx1487 + tx1488 + tx1490 + tx1492 + tx1493 + tx1494 + tx1495 + tx1496 + tx1499 + tx1502 + tx1503 + tx1504 + tx1505 + tx1507 + tx1508; 		/* rule 704 */
	 tx1512 = tx1325 + tx1328 + tx1334 + tx1339 + tx1343 + tx1344 + tx1348 + tx1349 + tx1373 + tx1382 + tx1486 + tx1489 + tx1490 + tx1491 + tx1492 + tx1495 + tx1498 + tx1501 + tx1506 + tx1507 + tx1509 + tx1510; 		/* rule 705 */
	 tx1513 = tx1008*tx1512*tx953; 		/* rule 706 */
	 tx1514 = tx1511*tx1832*tx953; 		/* rule 707 */
	 tx1515 = tx1513 + tx1514; 		/* rule 708 */
	 tx1516 = -2.*tx1209*tx1515*tx1834; 		/* rule 709 */
	 tx1517 = -2.*tx1209*tx1399*tx953*tx980; 		/* rule 710 */
	 tx1518 = power2(tx1032); 		/* rule 711 */
	 tx1519 = tx1402 + tx1403 + tx1404 + tx1406 + tx1516 + tx1517; 		/* rule 712 */
	 tx1520 = tx1064*tx1518; 		/* rule 713 */
	 tx1521 = DEImproperRestraintDPhi*tx1519; 		/* rule 714 */
	 dhz3z3 = tx1520 + tx1521; 		/* rule 715 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I3, 2, dhz3z3);
	 tx1522 = tx1064*tx1209*tx1832*tx958; 		/* rule 717 */
	 tx1523 = -6.*EraseLinearDihedral*K*TWOPI*tx1207*tx1209*tx955*tx980*tx985*tx990*tx991; 		/* rule 718 */
	 dhx4x4 = tx1522 + tx1523; 		/* rule 719 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 0, dhx4x4);
	 tx1524 = tx1064*tx1209*tx1832*tx957; 		/* rule 721 */
	 tx1525 = -6.*EraseLinearDihedral*K*TWOPI*tx1209*tx1305*tx954*tx980*tx985*tx990*tx991; 		/* rule 722 */
	 dhy4y4 = tx1524 + tx1525; 		/* rule 723 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I4, 1, I4, 1, dhy4y4);
	 tx1526 = tx1064*tx1209*tx1832*tx956; 		/* rule 725 */
	 tx1527 = -6.*EraseLinearDihedral*K*TWOPI*tx1209*tx1399*tx953*tx980*tx985*tx990*tx991; 		/* rule 726 */
	 dhz4z4 = tx1526 + tx1527; 		/* rule 727 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I4, 2, I4, 2, dhz4z4);
	#ifdef IMPROPER_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN //[
	if (calcOffDiagonalHessian) { 
	 tx1528 = tx1065*tx935; 		/* rule 731 */
	 tx1529 = tx1093*tx936; 		/* rule 732 */
	 tx1530 = tx1528 + tx1529; 		/* rule 733 */
	 tx1531 = tx1064*tx1066*tx1832*tx935*tx936; 		/* rule 734 */
	 tx1532 = DEImproperRestraintDPhi*tx1066*tx1530*tx980; 		/* rule 735 */
	 ohx1y1 = tx1531 + tx1532; 		/* rule 736 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tx1533 = tx1065*tx934; 		/* rule 738 */
	 tx1534 = tx1118*tx936; 		/* rule 739 */
	 tx1535 = tx1533 + tx1534; 		/* rule 740 */
	 tx1536 = tx1064*tx1066*tx1832*tx934*tx936; 		/* rule 741 */
	 tx1537 = DEImproperRestraintDPhi*tx1066*tx1535*tx980; 		/* rule 742 */
	 ohx1z1 = tx1536 + tx1537; 		/* rule 743 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx1538 = tx1066*tx1205*tx1834; 		/* rule 745 */
	 tx1539 = -2.*tx1065*tx1066*tx936*tx980; 		/* rule 746 */
	 tx1540 = tx1538 + tx1539; 		/* rule 747 */
	 tx1541 = -(tx1011*tx1012*tx1064); 		/* rule 748 */
	 tx1542 = DEImproperRestraintDPhi*tx1540; 		/* rule 749 */
	 ohx1x2 = tx1541 + tx1542; 		/* rule 750 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx1543 = tx1005*tx1256*tx936; 		/* rule 752 */
	 tx1544 = tx1157*tx1832*tx935; 		/* rule 753 */
	 tx1545 = tx1543 + tx1544; 		/* rule 754 */
	 tx1546 = tx1066*tx1545*tx1834; 		/* rule 755 */
	 tx1547 = -(tx1066*tx1530*tx980); 		/* rule 756 */
	 tx1548 = tx1546 + tx1547; 		/* rule 757 */
	 tx1549 = -(tx1011*tx1016*tx1064); 		/* rule 758 */
	 tx1550 = DEImproperRestraintDPhi*tx1548; 		/* rule 759 */
	 ohx1y2 = tx1549 + tx1550; 		/* rule 760 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx1551 = tx1005*tx1352*tx936; 		/* rule 762 */
	 tx1552 = tx1157*tx1832*tx934; 		/* rule 763 */
	 tx1553 = tx1551 + tx1552; 		/* rule 764 */
	 tx1554 = tx1066*tx1553*tx1834; 		/* rule 765 */
	 tx1555 = -(tx1066*tx1535*tx980); 		/* rule 766 */
	 tx1556 = tx1554 + tx1555; 		/* rule 767 */
	 tx1557 = -(tx1011*tx1020*tx1064); 		/* rule 768 */
	 tx1558 = DEImproperRestraintDPhi*tx1556; 		/* rule 769 */
	 ohx1z2 = tx1557 + tx1558; 		/* rule 770 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tx1559 = -(tx1011*tx1024*tx1064); 		/* rule 772 */
	 tx1560 = -3*EraseLinearDihedral*K*TWOPI*tx1538*tx985*tx990*tx991; 		/* rule 773 */
	 ohx1x3 = tx1559 + tx1560; 		/* rule 774 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 0, ohx1x3);
	 tx1561 = -(tx1011*tx1028*tx1064); 		/* rule 776 */
	 tx1562 = -3*EraseLinearDihedral*K*TWOPI*tx1546*tx985*tx990*tx991; 		/* rule 777 */
	 ohx1y3 = tx1561 + tx1562; 		/* rule 778 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 1, ohx1y3);
	 tx1563 = -(tx1011*tx1032*tx1064); 		/* rule 780 */
	 tx1564 = -3*EraseLinearDihedral*K*TWOPI*tx1554*tx985*tx990*tx991; 		/* rule 781 */
	 ohx1z3 = tx1563 + tx1564; 		/* rule 782 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 2, ohx1z3);
	 ohx1x4 = -(tx1006*tx1064*tx1832*tx962*tx992); 		/* rule 784 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 0, ohx1x4);
	 ohx1y4 = -(tx1006*tx1064*tx1832*tx936*tx954*tx992); 		/* rule 786 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 1, ohx1y4);
	 ohx1z4 = -(tx1006*tx1064*tx1832*tx936*tx953*tx992); 		/* rule 788 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 2, ohx1z4);
	 tx1565 = tx1093*tx934; 		/* rule 790 */
	 tx1566 = tx1118*tx935; 		/* rule 791 */
	 tx1567 = tx1565 + tx1566; 		/* rule 792 */
	 tx1568 = tx1064*tx1066*tx1832*tx934*tx935; 		/* rule 793 */
	 tx1569 = DEImproperRestraintDPhi*tx1066*tx1567*tx980; 		/* rule 794 */
	 ohy1z1 = tx1568 + tx1569; 		/* rule 795 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 tx1570 = tx1005*tx1158*tx935; 		/* rule 797 */
	 tx1571 = tx1257*tx1832*tx936; 		/* rule 798 */
	 tx1572 = tx1570 + tx1571; 		/* rule 799 */
	 tx1573 = tx1066*tx1572*tx1834; 		/* rule 800 */
	 tx1574 = tx1547 + tx1573; 		/* rule 801 */
	 tx1575 = -(tx1012*tx1015*tx1064); 		/* rule 802 */
	 tx1576 = DEImproperRestraintDPhi*tx1574; 		/* rule 803 */
	 ohy1x2 = tx1575 + tx1576; 		/* rule 804 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx1577 = tx1066*tx1303*tx1834; 		/* rule 806 */
	 tx1578 = -2.*tx1066*tx1093*tx935*tx980; 		/* rule 807 */
	 tx1579 = tx1577 + tx1578; 		/* rule 808 */
	 tx1580 = -(tx1015*tx1016*tx1064); 		/* rule 809 */
	 tx1581 = DEImproperRestraintDPhi*tx1579; 		/* rule 810 */
	 ohy1y2 = tx1580 + tx1581; 		/* rule 811 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx1582 = tx1005*tx1352*tx935; 		/* rule 813 */
	 tx1583 = tx1257*tx1832*tx934; 		/* rule 814 */
	 tx1584 = tx1582 + tx1583; 		/* rule 815 */
	 tx1585 = tx1066*tx1584*tx1834; 		/* rule 816 */
	 tx1586 = -(tx1066*tx1567*tx980); 		/* rule 817 */
	 tx1587 = tx1585 + tx1586; 		/* rule 818 */
	 tx1588 = -(tx1015*tx1020*tx1064); 		/* rule 819 */
	 tx1589 = DEImproperRestraintDPhi*tx1587; 		/* rule 820 */
	 ohy1z2 = tx1588 + tx1589; 		/* rule 821 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 tx1590 = -(tx1015*tx1024*tx1064); 		/* rule 823 */
	 tx1591 = -3*EraseLinearDihedral*K*TWOPI*tx1573*tx985*tx990*tx991; 		/* rule 824 */
	 ohy1x3 = tx1590 + tx1591; 		/* rule 825 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 0, ohy1x3);
	 tx1592 = -(tx1015*tx1028*tx1064); 		/* rule 827 */
	 tx1593 = -3*EraseLinearDihedral*K*TWOPI*tx1577*tx985*tx990*tx991; 		/* rule 828 */
	 ohy1y3 = tx1592 + tx1593; 		/* rule 829 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 1, ohy1y3);
	 tx1594 = -(tx1015*tx1032*tx1064); 		/* rule 831 */
	 tx1595 = -3*EraseLinearDihedral*K*TWOPI*tx1585*tx985*tx990*tx991; 		/* rule 832 */
	 ohy1z3 = tx1594 + tx1595; 		/* rule 833 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 2, ohy1z3);
	 ohy1x4 = -(tx1006*tx1064*tx1832*tx935*tx955*tx992); 		/* rule 835 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 0, ohy1x4);
	 ohy1y4 = -(tx1006*tx1064*tx1832*tx961*tx992); 		/* rule 837 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 1, ohy1y4);
	 ohy1z4 = -(tx1006*tx1064*tx1832*tx935*tx953*tx992); 		/* rule 839 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 2, ohy1z4);
	 tx1596 = tx1005*tx1158*tx934; 		/* rule 841 */
	 tx1597 = tx1351*tx1832*tx936; 		/* rule 842 */
	 tx1598 = tx1596 + tx1597; 		/* rule 843 */
	 tx1599 = tx1066*tx1598*tx1834; 		/* rule 844 */
	 tx1600 = tx1555 + tx1599; 		/* rule 845 */
	 tx1601 = -(tx1012*tx1019*tx1064); 		/* rule 846 */
	 tx1602 = DEImproperRestraintDPhi*tx1600; 		/* rule 847 */
	 ohz1x2 = tx1601 + tx1602; 		/* rule 848 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 tx1603 = tx1005*tx1256*tx934; 		/* rule 850 */
	 tx1604 = tx1351*tx1832*tx935; 		/* rule 851 */
	 tx1605 = tx1603 + tx1604; 		/* rule 852 */
	 tx1606 = tx1066*tx1605*tx1834; 		/* rule 853 */
	 tx1607 = tx1586 + tx1606; 		/* rule 854 */
	 tx1608 = -(tx1016*tx1019*tx1064); 		/* rule 855 */
	 tx1609 = DEImproperRestraintDPhi*tx1607; 		/* rule 856 */
	 ohz1y2 = tx1608 + tx1609; 		/* rule 857 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx1610 = tx1066*tx1397*tx1834; 		/* rule 859 */
	 tx1611 = -2.*tx1066*tx1118*tx934*tx980; 		/* rule 860 */
	 tx1612 = tx1610 + tx1611; 		/* rule 861 */
	 tx1613 = -(tx1019*tx1020*tx1064); 		/* rule 862 */
	 tx1614 = DEImproperRestraintDPhi*tx1612; 		/* rule 863 */
	 ohz1z2 = tx1613 + tx1614; 		/* rule 864 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 tx1615 = -(tx1019*tx1024*tx1064); 		/* rule 866 */
	 tx1616 = -3*EraseLinearDihedral*K*TWOPI*tx1599*tx985*tx990*tx991; 		/* rule 867 */
	 ohz1x3 = tx1615 + tx1616; 		/* rule 868 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 0, ohz1x3);
	 tx1617 = -(tx1019*tx1028*tx1064); 		/* rule 870 */
	 tx1618 = -3*EraseLinearDihedral*K*TWOPI*tx1606*tx985*tx990*tx991; 		/* rule 871 */
	 ohz1y3 = tx1617 + tx1618; 		/* rule 872 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 1, ohz1y3);
	 tx1619 = -(tx1019*tx1032*tx1064); 		/* rule 874 */
	 tx1620 = -3*EraseLinearDihedral*K*TWOPI*tx1610*tx985*tx990*tx991; 		/* rule 875 */
	 ohz1z3 = tx1619 + tx1620; 		/* rule 876 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 2, ohz1z3);
	 ohz1x4 = -(tx1006*tx1064*tx1832*tx934*tx955*tx992); 		/* rule 878 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 0, ohz1x4);
	 ohz1y4 = -(tx1006*tx1064*tx1832*tx934*tx954*tx992); 		/* rule 880 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 1, ohz1y4);
	 ohz1z4 = -(tx1006*tx1064*tx1832*tx960*tx992); 		/* rule 882 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 2, ohz1z4);
	 tx1621 = tx1206*tx935; 		/* rule 884 */
	 tx1622 = tx1304*tx936; 		/* rule 885 */
	 tx1623 = tx1207*tx954; 		/* rule 886 */
	 tx1624 = tx1208*tx954; 		/* rule 887 */
	 tx1625 = tx1305*tx955; 		/* rule 888 */
	 tx1626 = tx1306*tx955; 		/* rule 889 */
	 tx1627 = tx1621 + tx1622; 		/* rule 890 */
	 tx1628 = tx1623 + tx1625; 		/* rule 891 */
	 tx1629 = tx1624 + tx1626; 		/* rule 892 */
	 tx1630 = -tx1546; 		/* rule 893 */
	 tx1631 = -tx1573; 		/* rule 894 */
	 tx1632 = tx1005*tx1066*tx1627*tx1834; 		/* rule 895 */
	 tx1633 = -0.5*tx1006*tx1210*tx1628; 		/* rule 896 */
	 tx1634 = -(tx1008*tx1209*tx1629*tx1834); 		/* rule 897 */
	 tx1635 = tx1066*tx1530*tx980; 		/* rule 898 */
	 tx1636 = 0.5*tx1210*tx1530*tx992; 		/* rule 899 */
	 tx1637 = tx1630 + tx1631 + tx1632 + tx1633 + tx1634 + tx1635 + tx1636; 		/* rule 900 */
	 tx1638 = tx1012*tx1016*tx1064; 		/* rule 901 */
	 tx1639 = DEImproperRestraintDPhi*tx1637; 		/* rule 902 */
	 ohx2y2 = tx1638 + tx1639; 		/* rule 903 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 tx1640 = tx1206*tx934; 		/* rule 905 */
	 tx1641 = tx1398*tx936; 		/* rule 906 */
	 tx1642 = tx1207*tx953; 		/* rule 907 */
	 tx1643 = tx1208*tx953; 		/* rule 908 */
	 tx1644 = tx1399*tx955; 		/* rule 909 */
	 tx1645 = tx1400*tx955; 		/* rule 910 */
	 tx1646 = tx1640 + tx1641; 		/* rule 911 */
	 tx1647 = tx1642 + tx1644; 		/* rule 912 */
	 tx1648 = tx1643 + tx1645; 		/* rule 913 */
	 tx1649 = -tx1554; 		/* rule 914 */
	 tx1650 = -tx1599; 		/* rule 915 */
	 tx1651 = tx1005*tx1066*tx1646*tx1834; 		/* rule 916 */
	 tx1652 = -0.5*tx1006*tx1210*tx1647; 		/* rule 917 */
	 tx1653 = -(tx1008*tx1209*tx1648*tx1834); 		/* rule 918 */
	 tx1654 = tx1066*tx1535*tx980; 		/* rule 919 */
	 tx1655 = 0.5*tx1210*tx1535*tx992; 		/* rule 920 */
	 tx1656 = tx1649 + tx1650 + tx1651 + tx1652 + tx1653 + tx1654 + tx1655; 		/* rule 921 */
	 tx1657 = tx1012*tx1020*tx1064; 		/* rule 922 */
	 tx1658 = DEImproperRestraintDPhi*tx1656; 		/* rule 923 */
	 ohx2z2 = tx1657 + tx1658; 		/* rule 924 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 tx1659 = -tx1216; 		/* rule 926 */
	 tx1660 = tx1209*tx1441*tx1834; 		/* rule 927 */
	 tx1661 = -2.*tx1005*tx1066*tx1206*tx1834*tx936; 		/* rule 928 */
	 tx1662 = 2.*tx1008*tx1208*tx1209*tx1834*tx955; 		/* rule 929 */
	 tx1663 = tx1006*tx1207*tx1210*tx955; 		/* rule 930 */
	 tx1664 = tx1538 + tx1659 + tx1660 + tx1661 + tx1662 + tx1663; 		/* rule 931 */
	 tx1665 = tx1012*tx1024*tx1064; 		/* rule 932 */
	 tx1666 = DEImproperRestraintDPhi*tx1664; 		/* rule 933 */
	 ohx2x3 = tx1665 + tx1666; 		/* rule 934 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 0, ohx2x3);
	 tx1667 = tx1008*tx1438*tx954; 		/* rule 936 */
	 tx1668 = tx1474*tx1832*tx955; 		/* rule 937 */
	 tx1669 = tx1667 + tx1668; 		/* rule 938 */
	 tx1670 = 0.5*tx1006*tx1210*tx1628; 		/* rule 939 */
	 tx1671 = tx1008*tx1209*tx1629*tx1834; 		/* rule 940 */
	 tx1672 = -tx1632; 		/* rule 941 */
	 tx1673 = tx1209*tx1669*tx1834; 		/* rule 942 */
	 tx1674 = -0.5*tx1210*tx1530*tx992; 		/* rule 943 */
	 tx1675 = tx1546 + tx1670 + tx1671 + tx1672 + tx1673 + tx1674; 		/* rule 944 */
	 tx1676 = tx1012*tx1028*tx1064; 		/* rule 945 */
	 tx1677 = DEImproperRestraintDPhi*tx1675; 		/* rule 946 */
	 ohx2y3 = tx1676 + tx1677; 		/* rule 947 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 1, ohx2y3);
	 tx1678 = tx1008*tx1438*tx953; 		/* rule 949 */
	 tx1679 = tx1511*tx1832*tx955; 		/* rule 950 */
	 tx1680 = tx1678 + tx1679; 		/* rule 951 */
	 tx1681 = 0.5*tx1006*tx1210*tx1647; 		/* rule 952 */
	 tx1682 = tx1008*tx1209*tx1648*tx1834; 		/* rule 953 */
	 tx1683 = -tx1651; 		/* rule 954 */
	 tx1684 = tx1209*tx1680*tx1834; 		/* rule 955 */
	 tx1685 = -0.5*tx1210*tx1535*tx992; 		/* rule 956 */
	 tx1686 = tx1554 + tx1681 + tx1682 + tx1683 + tx1684 + tx1685; 		/* rule 957 */
	 tx1687 = tx1012*tx1032*tx1064; 		/* rule 958 */
	 tx1688 = DEImproperRestraintDPhi*tx1686; 		/* rule 959 */
	 ohx2z3 = tx1687 + tx1688; 		/* rule 960 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 2, ohx2z3);
	 tx1689 = tx1006*tx1012*tx1064*tx955*tx980; 		/* rule 962 */
	 tx1690 = -3*EraseLinearDihedral*K*TWOPI*tx1660*tx985*tx990*tx991; 		/* rule 963 */
	 ohx2x4 = tx1689 + tx1690; 		/* rule 964 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 0, ohx2x4);
	 tx1691 = tx1006*tx1012*tx1064*tx954*tx980; 		/* rule 966 */
	 tx1692 = -3*EraseLinearDihedral*K*TWOPI*tx1673*tx985*tx990*tx991; 		/* rule 967 */
	 ohx2y4 = tx1691 + tx1692; 		/* rule 968 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 1, ohx2y4);
	 tx1693 = tx1006*tx1012*tx1064*tx953*tx980; 		/* rule 970 */
	 tx1694 = -3*EraseLinearDihedral*K*TWOPI*tx1684*tx985*tx990*tx991; 		/* rule 971 */
	 ohx2z4 = tx1693 + tx1694; 		/* rule 972 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 2, ohx2z4);
	 tx1695 = tx1304*tx934; 		/* rule 974 */
	 tx1696 = tx1398*tx935; 		/* rule 975 */
	 tx1697 = tx1305*tx953; 		/* rule 976 */
	 tx1698 = tx1306*tx953; 		/* rule 977 */
	 tx1699 = tx1399*tx954; 		/* rule 978 */
	 tx1700 = tx1400*tx954; 		/* rule 979 */
	 tx1701 = tx1695 + tx1696; 		/* rule 980 */
	 tx1702 = tx1697 + tx1699; 		/* rule 981 */
	 tx1703 = tx1698 + tx1700; 		/* rule 982 */
	 tx1704 = -tx1585; 		/* rule 983 */
	 tx1705 = -tx1606; 		/* rule 984 */
	 tx1706 = tx1005*tx1066*tx1701*tx1834; 		/* rule 985 */
	 tx1707 = -0.5*tx1006*tx1210*tx1702; 		/* rule 986 */
	 tx1708 = -(tx1008*tx1209*tx1703*tx1834); 		/* rule 987 */
	 tx1709 = tx1066*tx1567*tx980; 		/* rule 988 */
	 tx1710 = 0.5*tx1210*tx1567*tx992; 		/* rule 989 */
	 tx1711 = tx1704 + tx1705 + tx1706 + tx1707 + tx1708 + tx1709 + tx1710; 		/* rule 990 */
	 tx1712 = tx1016*tx1020*tx1064; 		/* rule 991 */
	 tx1713 = DEImproperRestraintDPhi*tx1711; 		/* rule 992 */
	 ohy2z2 = tx1712 + tx1713; 		/* rule 993 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	 tx1714 = tx1008*tx1475*tx955; 		/* rule 995 */
	 tx1715 = tx1437*tx1832*tx954; 		/* rule 996 */
	 tx1716 = tx1714 + tx1715; 		/* rule 997 */
	 tx1717 = tx1209*tx1716*tx1834; 		/* rule 998 */
	 tx1718 = tx1573 + tx1670 + tx1671 + tx1672 + tx1674 + tx1717; 		/* rule 999 */
	 tx1719 = tx1016*tx1024*tx1064; 		/* rule 1000 */
	 tx1720 = DEImproperRestraintDPhi*tx1718; 		/* rule 1001 */
	 ohy2x3 = tx1719 + tx1720; 		/* rule 1002 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 0, ohy2x3);
	 tx1721 = -tx1312; 		/* rule 1004 */
	 tx1722 = tx1209*tx1478*tx1834; 		/* rule 1005 */
	 tx1723 = -2.*tx1005*tx1066*tx1304*tx1834*tx935; 		/* rule 1006 */
	 tx1724 = tx1006*tx1210*tx1305*tx954; 		/* rule 1007 */
	 tx1725 = 2.*tx1008*tx1209*tx1306*tx1834*tx954; 		/* rule 1008 */
	 tx1726 = tx1577 + tx1721 + tx1722 + tx1723 + tx1724 + tx1725; 		/* rule 1009 */
	 tx1727 = tx1016*tx1028*tx1064; 		/* rule 1010 */
	 tx1728 = DEImproperRestraintDPhi*tx1726; 		/* rule 1011 */
	 ohy2y3 = tx1727 + tx1728; 		/* rule 1012 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 1, ohy2y3);
	 tx1729 = tx1008*tx1475*tx953; 		/* rule 1014 */
	 tx1730 = tx1511*tx1832*tx954; 		/* rule 1015 */
	 tx1731 = tx1729 + tx1730; 		/* rule 1016 */
	 tx1732 = 0.5*tx1006*tx1210*tx1702; 		/* rule 1017 */
	 tx1733 = tx1008*tx1209*tx1703*tx1834; 		/* rule 1018 */
	 tx1734 = -tx1706; 		/* rule 1019 */
	 tx1735 = tx1209*tx1731*tx1834; 		/* rule 1020 */
	 tx1736 = -0.5*tx1210*tx1567*tx992; 		/* rule 1021 */
	 tx1737 = tx1585 + tx1732 + tx1733 + tx1734 + tx1735 + tx1736; 		/* rule 1022 */
	 tx1738 = tx1016*tx1032*tx1064; 		/* rule 1023 */
	 tx1739 = DEImproperRestraintDPhi*tx1737; 		/* rule 1024 */
	 ohy2z3 = tx1738 + tx1739; 		/* rule 1025 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 2, ohy2z3);
	 tx1740 = tx1006*tx1016*tx1064*tx955*tx980; 		/* rule 1027 */
	 tx1741 = -3*EraseLinearDihedral*K*TWOPI*tx1717*tx985*tx990*tx991; 		/* rule 1028 */
	 ohy2x4 = tx1740 + tx1741; 		/* rule 1029 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 0, ohy2x4);
	 tx1742 = tx1006*tx1016*tx1064*tx954*tx980; 		/* rule 1031 */
	 tx1743 = -3*EraseLinearDihedral*K*TWOPI*tx1722*tx985*tx990*tx991; 		/* rule 1032 */
	 ohy2y4 = tx1742 + tx1743; 		/* rule 1033 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 1, ohy2y4);
	 tx1744 = tx1006*tx1016*tx1064*tx953*tx980; 		/* rule 1035 */
	 tx1745 = -3*EraseLinearDihedral*K*TWOPI*tx1735*tx985*tx990*tx991; 		/* rule 1036 */
	 ohy2z4 = tx1744 + tx1745; 		/* rule 1037 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 2, ohy2z4);
	 tx1746 = tx1008*tx1512*tx955; 		/* rule 1039 */
	 tx1747 = tx1437*tx1832*tx953; 		/* rule 1040 */
	 tx1748 = tx1746 + tx1747; 		/* rule 1041 */
	 tx1749 = tx1209*tx1748*tx1834; 		/* rule 1042 */
	 tx1750 = tx1599 + tx1681 + tx1682 + tx1683 + tx1685 + tx1749; 		/* rule 1043 */
	 tx1751 = tx1020*tx1024*tx1064; 		/* rule 1044 */
	 tx1752 = DEImproperRestraintDPhi*tx1750; 		/* rule 1045 */
	 ohz2x3 = tx1751 + tx1752; 		/* rule 1046 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 0, ohz2x3);
	 tx1753 = tx1008*tx1512*tx954; 		/* rule 1048 */
	 tx1754 = tx1474*tx1832*tx953; 		/* rule 1049 */
	 tx1755 = tx1753 + tx1754; 		/* rule 1050 */
	 tx1756 = tx1209*tx1755*tx1834; 		/* rule 1051 */
	 tx1757 = tx1606 + tx1732 + tx1733 + tx1734 + tx1736 + tx1756; 		/* rule 1052 */
	 tx1758 = tx1020*tx1028*tx1064; 		/* rule 1053 */
	 tx1759 = DEImproperRestraintDPhi*tx1757; 		/* rule 1054 */
	 ohz2y3 = tx1758 + tx1759; 		/* rule 1055 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 1, ohz2y3);
	 tx1760 = -tx1406; 		/* rule 1057 */
	 tx1761 = tx1209*tx1515*tx1834; 		/* rule 1058 */
	 tx1762 = -2.*tx1005*tx1066*tx1398*tx1834*tx934; 		/* rule 1059 */
	 tx1763 = tx1006*tx1210*tx1399*tx953; 		/* rule 1060 */
	 tx1764 = 2.*tx1008*tx1209*tx1400*tx1834*tx953; 		/* rule 1061 */
	 tx1765 = tx1610 + tx1760 + tx1761 + tx1762 + tx1763 + tx1764; 		/* rule 1062 */
	 tx1766 = tx1020*tx1032*tx1064; 		/* rule 1063 */
	 tx1767 = DEImproperRestraintDPhi*tx1765; 		/* rule 1064 */
	 ohz2z3 = tx1766 + tx1767; 		/* rule 1065 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 2, ohz2z3);
	 tx1768 = tx1006*tx1020*tx1064*tx955*tx980; 		/* rule 1067 */
	 tx1769 = -3*EraseLinearDihedral*K*TWOPI*tx1749*tx985*tx990*tx991; 		/* rule 1068 */
	 ohz2x4 = tx1768 + tx1769; 		/* rule 1069 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 0, ohz2x4);
	 tx1770 = tx1006*tx1020*tx1064*tx954*tx980; 		/* rule 1071 */
	 tx1771 = -3*EraseLinearDihedral*K*TWOPI*tx1756*tx985*tx990*tx991; 		/* rule 1072 */
	 ohz2y4 = tx1770 + tx1771; 		/* rule 1073 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 1, ohz2y4);
	 tx1772 = tx1006*tx1020*tx1064*tx953*tx980; 		/* rule 1075 */
	 tx1773 = -3*EraseLinearDihedral*K*TWOPI*tx1761*tx985*tx990*tx991; 		/* rule 1076 */
	 ohz2z4 = tx1772 + tx1773; 		/* rule 1077 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 2, ohz2z4);
	 tx1774 = -tx1673; 		/* rule 1079 */
	 tx1775 = -tx1717; 		/* rule 1080 */
	 tx1776 = -(tx1209*tx1628*tx980); 		/* rule 1081 */
	 tx1777 = tx1632 + tx1633 + tx1634 + tx1636 + tx1774 + tx1775 + tx1776; 		/* rule 1082 */
	 tx1778 = tx1024*tx1028*tx1064; 		/* rule 1083 */
	 tx1779 = DEImproperRestraintDPhi*tx1777; 		/* rule 1084 */
	 ohx3y3 = tx1778 + tx1779; 		/* rule 1085 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 1, ohx3y3);
	 tx1780 = -tx1684; 		/* rule 1087 */
	 tx1781 = -tx1749; 		/* rule 1088 */
	 tx1782 = -(tx1209*tx1647*tx980); 		/* rule 1089 */
	 tx1783 = tx1651 + tx1652 + tx1653 + tx1655 + tx1780 + tx1781 + tx1782; 		/* rule 1090 */
	 tx1784 = tx1024*tx1032*tx1064; 		/* rule 1091 */
	 tx1785 = DEImproperRestraintDPhi*tx1783; 		/* rule 1092 */
	 ohx3z3 = tx1784 + tx1785; 		/* rule 1093 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 2, ohx3z3);
	 tx1786 = 2.*tx1207*tx1209*tx955*tx980; 		/* rule 1095 */
	 tx1787 = tx1660 + tx1786; 		/* rule 1096 */
	 tx1788 = DEImproperRestraintDPhi*tx1787; 		/* rule 1097 */
	 tx1789 = tx1006*tx1024*tx1064*tx955*tx980; 		/* rule 1098 */
	 ohx3x4 = tx1788 + tx1789; 		/* rule 1099 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 0, ohx3x4);
	 tx1790 = tx1209*tx1628*tx980; 		/* rule 1101 */
	 tx1791 = tx1673 + tx1790; 		/* rule 1102 */
	 tx1792 = DEImproperRestraintDPhi*tx1791; 		/* rule 1103 */
	 tx1793 = tx1006*tx1024*tx1064*tx954*tx980; 		/* rule 1104 */
	 ohx3y4 = tx1792 + tx1793; 		/* rule 1105 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 1, ohx3y4);
	 tx1794 = tx1209*tx1647*tx980; 		/* rule 1107 */
	 tx1795 = tx1684 + tx1794; 		/* rule 1108 */
	 tx1796 = DEImproperRestraintDPhi*tx1795; 		/* rule 1109 */
	 tx1797 = tx1006*tx1024*tx1064*tx953*tx980; 		/* rule 1110 */
	 ohx3z4 = tx1796 + tx1797; 		/* rule 1111 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 2, ohx3z4);
	 tx1798 = -tx1735; 		/* rule 1113 */
	 tx1799 = -tx1756; 		/* rule 1114 */
	 tx1800 = -(tx1209*tx1702*tx980); 		/* rule 1115 */
	 tx1801 = tx1706 + tx1707 + tx1708 + tx1710 + tx1798 + tx1799 + tx1800; 		/* rule 1116 */
	 tx1802 = tx1028*tx1032*tx1064; 		/* rule 1117 */
	 tx1803 = DEImproperRestraintDPhi*tx1801; 		/* rule 1118 */
	 ohy3z3 = tx1802 + tx1803; 		/* rule 1119 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 2, ohy3z3);
	 tx1804 = tx1717 + tx1790; 		/* rule 1121 */
	 tx1805 = DEImproperRestraintDPhi*tx1804; 		/* rule 1122 */
	 tx1806 = tx1006*tx1028*tx1064*tx955*tx980; 		/* rule 1123 */
	 ohy3x4 = tx1805 + tx1806; 		/* rule 1124 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 0, ohy3x4);
	 tx1807 = 2.*tx1209*tx1305*tx954*tx980; 		/* rule 1126 */
	 tx1808 = tx1722 + tx1807; 		/* rule 1127 */
	 tx1809 = DEImproperRestraintDPhi*tx1808; 		/* rule 1128 */
	 tx1810 = tx1006*tx1028*tx1064*tx954*tx980; 		/* rule 1129 */
	 ohy3y4 = tx1809 + tx1810; 		/* rule 1130 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 1, ohy3y4);
	 tx1811 = tx1209*tx1702*tx980; 		/* rule 1132 */
	 tx1812 = tx1735 + tx1811; 		/* rule 1133 */
	 tx1813 = DEImproperRestraintDPhi*tx1812; 		/* rule 1134 */
	 tx1814 = tx1006*tx1028*tx1064*tx953*tx980; 		/* rule 1135 */
	 ohy3z4 = tx1813 + tx1814; 		/* rule 1136 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 2, ohy3z4);
	 tx1815 = tx1749 + tx1794; 		/* rule 1138 */
	 tx1816 = DEImproperRestraintDPhi*tx1815; 		/* rule 1139 */
	 tx1817 = tx1006*tx1032*tx1064*tx955*tx980; 		/* rule 1140 */
	 ohz3x4 = tx1816 + tx1817; 		/* rule 1141 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 0, ohz3x4);
	 tx1818 = tx1756 + tx1811; 		/* rule 1143 */
	 tx1819 = DEImproperRestraintDPhi*tx1818; 		/* rule 1144 */
	 tx1820 = tx1006*tx1032*tx1064*tx954*tx980; 		/* rule 1145 */
	 ohz3y4 = tx1819 + tx1820; 		/* rule 1146 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 1, ohz3y4);
	 tx1821 = 2.*tx1209*tx1399*tx953*tx980; 		/* rule 1148 */
	 tx1822 = tx1761 + tx1821; 		/* rule 1149 */
	 tx1823 = DEImproperRestraintDPhi*tx1822; 		/* rule 1150 */
	 tx1824 = tx1006*tx1032*tx1064*tx953*tx980; 		/* rule 1151 */
	 ohz3z4 = tx1823 + tx1824; 		/* rule 1152 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 2, ohz3z4);
	 tx1825 = tx1064*tx1209*tx1832*tx954*tx955; 		/* rule 1154 */
	 tx1826 = -3*EraseLinearDihedral*K*TWOPI*tx1790*tx985*tx990*tx991; 		/* rule 1155 */
	 ohx4y4 = tx1825 + tx1826; 		/* rule 1156 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 1, ohx4y4);
	 tx1827 = tx1064*tx1209*tx1832*tx953*tx955; 		/* rule 1158 */
	 tx1828 = -3*EraseLinearDihedral*K*TWOPI*tx1794*tx985*tx990*tx991; 		/* rule 1159 */
	 ohx4z4 = tx1827 + tx1828; 		/* rule 1160 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 2, ohx4z4);
	 tx1829 = tx1064*tx1209*tx1832*tx953*tx954; 		/* rule 1162 */
	 tx1830 = -3*EraseLinearDihedral*K*TWOPI*tx1811*tx985*tx990*tx991; 		/* rule 1163 */
	 ohy4z4 = tx1829 + tx1830; 		/* rule 1164 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 1, I4, 2, ohy4z4);
	} /*calcOffDiagonalHessian*/
	#endif //IMPROPER_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN ]
	} /*calcDiagonalHessian*/
	#endif //IMPROPER_RESTRAINT_CALC_DIAGONAL_HESSIAN ]
	} /*calcForce*/
	#endif //IMPROPER_RESTRAINT_CALC_FORCE ]
	} /*RestraintActive*/
	IMPROPER_RESTRAINT_PHI_SET(Phi);
/*]; (*RestraintActive*)*/
