// C-code
	IMPROPER_RESTRAINT_SET_PARAMETER(K);
	IMPROPER_RESTRAINT_SET_PARAMETER(U);
	IMPROPER_RESTRAINT_SET_PARAMETER(L);
	IMPROPER_RESTRAINT_SET_PARAMETER(I1);
	IMPROPER_RESTRAINT_SET_PARAMETER(I2);
	IMPROPER_RESTRAINT_SET_PARAMETER(I3);
	IMPROPER_RESTRAINT_SET_PARAMETER(I4);
	IMPROPER_RESTRAINT_APPLY_ATOM_MASK(I1,I2,I3,I4);
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
	 tx1 = -(x2*y1); 		/* rule 22 */
	 tx2 = x3*y1; 		/* rule 23 */
	 tx3 = x1*y2; 		/* rule 24 */
	 tx4 = -(x3*y2); 		/* rule 25 */
	 tx5 = -(x1*y3); 		/* rule 26 */
	 tx6 = x2*y3; 		/* rule 27 */
	 tx7 = x2*z1; 		/* rule 28 */
	 tx8 = -(x3*z1); 		/* rule 29 */
	 tx9 = -(y2*z1); 		/* rule 30 */
	 tx10 = y3*z1; 		/* rule 31 */
	 tx11 = -(x1*z2); 		/* rule 32 */
	 tx12 = x3*z2; 		/* rule 33 */
	 tx13 = y1*z2; 		/* rule 34 */
	 tx14 = -(y3*z2); 		/* rule 35 */
	 tx15 = x1*z3; 		/* rule 36 */
	 tx16 = -(x2*z3); 		/* rule 37 */
	 tx17 = -(y1*z3); 		/* rule 38 */
	 tx18 = y2*z3; 		/* rule 39 */
	 tx19 = tx1 + tx2 + tx3 + tx4 + tx5 + tx6; 		/* rule 40 */
	 tx20 = tx11 + tx12 + tx15 + tx16 + tx7 + tx8; 		/* rule 41 */
	 tx21 = tx10 + tx13 + tx14 + tx17 + tx18 + tx9; 		/* rule 42 */
	 tx22 = power2(tx19); 		/* rule 43 */
	 tx23 = power2(tx20); 		/* rule 44 */
	 tx24 = power2(tx21); 		/* rule 45 */
	 tx25 = tx22 + tx23 + tx24; 		/* rule 46 */
	 LenA = mysqrt(tx25); 		/* rule 47 */
	 tx26 = x4*y2; 		/* rule 48 */
	 tx27 = -(x4*y3); 		/* rule 49 */
	 tx28 = -(x2*y4); 		/* rule 50 */
	 tx29 = x3*y4; 		/* rule 51 */
	 tx30 = -(x4*z2); 		/* rule 52 */
	 tx31 = y4*z2; 		/* rule 53 */
	 tx32 = x4*z3; 		/* rule 54 */
	 tx33 = -(y4*z3); 		/* rule 55 */
	 tx34 = x2*z4; 		/* rule 56 */
	 tx35 = -(x3*z4); 		/* rule 57 */
	 tx36 = -(y2*z4); 		/* rule 58 */
	 tx37 = y3*z4; 		/* rule 59 */
	 tx38 = tx12 + tx16 + tx30 + tx32 + tx34 + tx35; 		/* rule 60 */
	 tx39 = tx14 + tx18 + tx31 + tx33 + tx36 + tx37; 		/* rule 61 */
	 tx40 = tx26 + tx27 + tx28 + tx29 + tx4 + tx6; 		/* rule 62 */
	 tx41 = power2(tx38); 		/* rule 63 */
	 tx42 = power2(tx39); 		/* rule 64 */
	 tx43 = power2(tx40); 		/* rule 65 */
	 tx44 = tx41 + tx42 + tx43; 		/* rule 66 */
	 LenB = mysqrt(tx44); 		/* rule 67 */
	 ReciprocalLenA = reciprocal(LenA); 		/* rule 68 */
	 ReciprocalLenB = reciprocal(LenB); 		/* rule 69 */
	ZERO_SMALL_LEN(LenA,ReciprocalLenA);
	ZERO_SMALL_LEN(LenB,ReciprocalLenB);
	 RecLenARecLenB = ReciprocalLenA*ReciprocalLenB; 		/* rule 72 */
	EraseLinearDihedral = VEC_CONST(1.0);
	EraseLinearDihedral = (RecLenARecLenB==0.0) ? 0.0 : EraseLinearDihedral;
	 tx45 = tx20*tx38; 		/* rule 75 */
	 tx46 = tx21*tx39; 		/* rule 76 */
	 tx47 = tx19*tx40; 		/* rule 77 */
	 tx48 = tx45 + tx46 + tx47; 		/* rule 78 */
	 CosPhi = RecLenARecLenB*tx48; 		/* rule 79 */
	 tx49 = -x3; 		/* rule 80 */
	 tx50 = -y3; 		/* rule 81 */
	 tx51 = -z3; 		/* rule 82 */
	 tx52 = tx49 + x2; 		/* rule 83 */
	 tx53 = tx50 + y2; 		/* rule 84 */
	 tx54 = tx51 + z2; 		/* rule 85 */
	 tx55 = power2(tx52); 		/* rule 86 */
	 tx56 = power2(tx53); 		/* rule 87 */
	 tx57 = power2(tx54); 		/* rule 88 */
	 tx58 = tx49 + x4; 		/* rule 89 */
	 tx59 = tx50 + y4; 		/* rule 90 */
	 tx60 = tx51 + z4; 		/* rule 91 */
	 tx61 = tx55 + tx56 + tx57; 		/* rule 92 */
	 tx62 = tx21*tx58; 		/* rule 93 */
	 tx63 = tx20*tx59; 		/* rule 94 */
	 tx64 = tx19*tx60; 		/* rule 95 */
	 tx65 = mysqrt(tx61); 		/* rule 96 */
	 tx66 = tx62 + tx63 + tx64; 		/* rule 97 */
	 SinPhi = RecLenARecLenB*tx65*tx66; 		/* rule 98 */
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
	 tx67 = -L; 		/* rule 112 */
	 tx68 = tx67 + U + UShift; 		/* rule 113 */
	 tx69 = Phi + PhiShift + tx67; 		/* rule 114 */
	 tx70 = reciprocal(tx68); 		/* rule 115 */
	 Psi = TWOPI*tx69*tx70; 		/* rule 116 */
	 tx71 = Cos(Psi); 		/* rule 117 */
	 tx72 = -tx71; 		/* rule 118 */
	 tx73 = 1 + tx72; 		/* rule 119 */
	 tx916 = power2(tx73); 		/* rule 120 */
	 tx74 = tx73*tx916; 		/* rule 121 */
	 Energy = EraseLinearDihedral*K*tx74; 		/* rule 122 */
	IMPROPER_RESTRAINT_ENERGY_ACCUMULATE(Energy);
	#ifdef IMPROPER_RESTRAINT_CALC_FORCE //[
	if (calcForce ) {
	 tx75 = tx916; 		/* rule 126 */
	 tx76 = Sin(Psi); 		/* rule 127 */
	 DEImproperRestraintDPhi = 3*EraseLinearDihedral*K*TWOPI*tx70*tx75*tx76; 		/* rule 128 */
	 tx77 = reciprocal(tx25); 		/* rule 129 */
	 gx1 = -(DEImproperRestraintDPhi*tx21*tx65*tx77); 		/* rule 130 */
	 fx1 = -gx1; 		/* rule 131 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I1, 0, fx1 );
	 gy1 = -(DEImproperRestraintDPhi*tx20*tx65*tx77); 		/* rule 133 */
	 fy1 = -gy1; 		/* rule 134 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I1, 1, fy1 );
	 gz1 = -(DEImproperRestraintDPhi*tx19*tx65*tx77); 		/* rule 136 */
	 fz1 = -gz1; 		/* rule 137 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tx78 = -x2; 		/* rule 139 */
	 tx79 = -y2; 		/* rule 140 */
	 tx80 = -z2; 		/* rule 141 */
	 tx81 = tx78 + x1; 		/* rule 142 */
	 tx82 = tx79 + y1; 		/* rule 143 */
	 tx83 = tx80 + z1; 		/* rule 144 */
	 tx84 = tx52*tx58; 		/* rule 145 */
	 tx85 = tx53*tx59; 		/* rule 146 */
	 tx86 = tx54*tx60; 		/* rule 147 */
	 tx87 = tx52*tx81; 		/* rule 148 */
	 tx88 = tx53*tx82; 		/* rule 149 */
	 tx89 = tx54*tx83; 		/* rule 150 */
	 tx90 = reciprocal(tx44); 		/* rule 151 */
	 tx91 = reciprocal(tx65); 		/* rule 152 */
	 tx92 = tx84 + tx85 + tx86; 		/* rule 153 */
	 tx93 = tx87 + tx88 + tx89; 		/* rule 154 */
	 tx94 = tx21*tx65*tx77; 		/* rule 155 */
	 tx95 = -(tx39*tx90*tx91*tx92); 		/* rule 156 */
	 tx96 = tx21*tx77*tx91*tx93; 		/* rule 157 */
	 tx97 = tx94 + tx95 + tx96; 		/* rule 158 */
	 gx2 = DEImproperRestraintDPhi*tx97; 		/* rule 159 */
	 fx2 = -gx2; 		/* rule 160 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I2, 0, fx2 );
	 tx98 = tx20*tx65*tx77; 		/* rule 162 */
	 tx99 = -(tx38*tx90*tx91*tx92); 		/* rule 163 */
	 tx100 = tx20*tx77*tx91*tx93; 		/* rule 164 */
	 tx101 = tx100 + tx98 + tx99; 		/* rule 165 */
	 gy2 = DEImproperRestraintDPhi*tx101; 		/* rule 166 */
	 fy2 = -gy2; 		/* rule 167 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I2, 1, fy2 );
	 tx102 = tx19*tx65*tx77; 		/* rule 169 */
	 tx103 = -(tx40*tx90*tx91*tx92); 		/* rule 170 */
	 tx104 = tx19*tx77*tx91*tx93; 		/* rule 171 */
	 tx105 = tx102 + tx103 + tx104; 		/* rule 172 */
	 gz2 = DEImproperRestraintDPhi*tx105; 		/* rule 173 */
	 fz2 = -gz2; 		/* rule 174 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I2, 2, fz2 );
	 tx106 = -(tx39*tx65*tx90); 		/* rule 176 */
	 tx107 = tx39*tx90*tx91*tx92; 		/* rule 177 */
	 tx108 = -tx96; 		/* rule 178 */
	 tx109 = tx106 + tx107 + tx108; 		/* rule 179 */
	 gx3 = DEImproperRestraintDPhi*tx109; 		/* rule 180 */
	 fx3 = -gx3; 		/* rule 181 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I3, 0, fx3 );
	 tx110 = -tx100; 		/* rule 183 */
	 tx111 = -(tx38*tx65*tx90); 		/* rule 184 */
	 tx112 = tx38*tx90*tx91*tx92; 		/* rule 185 */
	 tx113 = tx110 + tx111 + tx112; 		/* rule 186 */
	 gy3 = DEImproperRestraintDPhi*tx113; 		/* rule 187 */
	 fy3 = -gy3; 		/* rule 188 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I3, 1, fy3 );
	 tx114 = -tx104; 		/* rule 190 */
	 tx115 = -(tx40*tx65*tx90); 		/* rule 191 */
	 tx116 = tx40*tx90*tx91*tx92; 		/* rule 192 */
	 tx117 = tx114 + tx115 + tx116; 		/* rule 193 */
	 gz3 = DEImproperRestraintDPhi*tx117; 		/* rule 194 */
	 fz3 = -gz3; 		/* rule 195 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I3, 2, fz3 );
	 gx4 = DEImproperRestraintDPhi*tx39*tx65*tx90; 		/* rule 197 */
	 fx4 = -gx4; 		/* rule 198 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I4, 0, fx4 );
	 gy4 = DEImproperRestraintDPhi*tx38*tx65*tx90; 		/* rule 200 */
	 fy4 = -gy4; 		/* rule 201 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I4, 1, fy4 );
	 gz4 = DEImproperRestraintDPhi*tx40*tx65*tx90; 		/* rule 203 */
	 fz4 = -gz4; 		/* rule 204 */
	IMPROPER_RESTRAINT_FORCE_ACCUMULATE(I4, 2, fz4 );
	#ifdef IMPROPER_RESTRAINT_CALC_DIAGONAL_HESSIAN //[
	if (calcDiagonalHessian) {
	 tx118 = power2(TWOPI); 		/* rule 208 */
	 tx119 = power2(tx70); 		/* rule 209 */
	 tx120 = power2(tx76); 		/* rule 210 */
	 tx121 = power2(y2); 		/* rule 211 */
	 tx122 = power2(y3); 		/* rule 212 */
	 tx123 = power2(z2); 		/* rule 213 */
	 tx124 = power2(z3); 		/* rule 214 */
	 tx125 = tx121*tx49; 		/* rule 215 */
	 tx126 = tx123*tx49; 		/* rule 216 */
	 tx127 = tx2*tx50; 		/* rule 217 */
	 tx128 = 6.*EraseLinearDihedral*K*tx118*tx119*tx120*tx73; 		/* rule 218 */
	 tx129 = 3*EraseLinearDihedral*K*tx118*tx119*tx71*tx75; 		/* rule 219 */
	 tx130 = tx122*tx78; 		/* rule 220 */
	 tx131 = tx124*tx78; 		/* rule 221 */
	 tx132 = tx7*tx80; 		/* rule 222 */
	 tx133 = tx121*x1; 		/* rule 223 */
	 tx134 = tx122*x1; 		/* rule 224 */
	 tx135 = tx123*x1; 		/* rule 225 */
	 tx136 = tx124*x1; 		/* rule 226 */
	 tx137 = tx6*y1; 		/* rule 227 */
	 tx138 = tx1*y2; 		/* rule 228 */
	 tx139 = tx2*y2; 		/* rule 229 */
	 tx140 = tx6*y2; 		/* rule 230 */
	 tx141 = -2.*tx3*y3; 		/* rule 231 */
	 tx142 = x3*y2*y3; 		/* rule 232 */
	 tx143 = tx12*z1; 		/* rule 233 */
	 tx144 = -2.*tx15*z2; 		/* rule 234 */
	 tx145 = tx12*z3; 		/* rule 235 */
	 tx146 = tx7*z3; 		/* rule 236 */
	 tx147 = tx8*z3; 		/* rule 237 */
	 tx148 = x2*z2*z3; 		/* rule 238 */
	 tx149 = tx128 + tx129; 		/* rule 239 */
	 tx150 = tx125 + tx126 + tx127 + tx130 + tx131 + tx132 + tx133 + tx134 + tx135 + tx136 + tx137 + tx138 + tx139 + tx140 + tx141 + tx142 + tx143 + tx144 + tx145 + tx146 + tx147 + tx148; 		/* rule 240 */
	 tx151 = power2(tx77); 		/* rule 241 */
	 tx152 = tx149*tx151*tx24*tx61; 		/* rule 242 */
	 tx153 = 6.*EraseLinearDihedral*K*TWOPI*tx150*tx151*tx21*tx65*tx70*tx75*tx76; 		/* rule 243 */
	 dhx1x1 = tx152 + tx153; 		/* rule 244 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx154 = power2(x2); 		/* rule 246 */
	 tx155 = power2(x3); 		/* rule 247 */
	 tx156 = tx123*tx50; 		/* rule 248 */
	 tx157 = tx154*tx50; 		/* rule 249 */
	 tx158 = tx10*tx51; 		/* rule 250 */
	 tx159 = tx3*tx78; 		/* rule 251 */
	 tx160 = tx124*tx79; 		/* rule 252 */
	 tx161 = tx155*tx79; 		/* rule 253 */
	 tx162 = tx6*x1; 		/* rule 254 */
	 tx163 = -2.*tx2*x2; 		/* rule 255 */
	 tx164 = tx3*x3; 		/* rule 256 */
	 tx165 = tx5*x3; 		/* rule 257 */
	 tx166 = tx6*x3; 		/* rule 258 */
	 tx167 = tx123*y1; 		/* rule 259 */
	 tx168 = tx124*y1; 		/* rule 260 */
	 tx169 = tx154*y1; 		/* rule 261 */
	 tx170 = tx155*y1; 		/* rule 262 */
	 tx171 = x2*x3*y2; 		/* rule 263 */
	 tx172 = tx18*z1; 		/* rule 264 */
	 tx173 = tx10*z2; 		/* rule 265 */
	 tx174 = tx18*z2; 		/* rule 266 */
	 tx175 = tx9*z2; 		/* rule 267 */
	 tx176 = -2.*tx13*z3; 		/* rule 268 */
	 tx177 = y3*z2*z3; 		/* rule 269 */
	 tx178 = tx156 + tx157 + tx158 + tx159 + tx160 + tx161 + tx162 + tx163 + tx164 + tx165 + tx166 + tx167 + tx168 + tx169 + tx170 + tx171 + tx172 + tx173 + tx174 + tx175 + tx176 + tx177; 		/* rule 270 */
	 tx179 = tx149*tx151*tx23*tx61; 		/* rule 271 */
	 tx180 = 6.*EraseLinearDihedral*K*TWOPI*tx151*tx178*tx20*tx65*tx70*tx75*tx76; 		/* rule 272 */
	 dhy1y1 = tx179 + tx180; 		/* rule 273 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx181 = tx15*tx49; 		/* rule 275 */
	 tx182 = tx121*tx51; 		/* rule 276 */
	 tx183 = tx154*tx51; 		/* rule 277 */
	 tx184 = tx13*tx79; 		/* rule 278 */
	 tx185 = tx122*tx80; 		/* rule 279 */
	 tx186 = tx155*tx80; 		/* rule 280 */
	 tx187 = tx12*x1; 		/* rule 281 */
	 tx188 = tx11*x2; 		/* rule 282 */
	 tx189 = tx12*x2; 		/* rule 283 */
	 tx190 = tx15*x2; 		/* rule 284 */
	 tx191 = -2.*tx7*x3; 		/* rule 285 */
	 tx192 = tx18*y1; 		/* rule 286 */
	 tx193 = -2.*tx10*y2; 		/* rule 287 */
	 tx194 = tx13*y3; 		/* rule 288 */
	 tx195 = tx17*y3; 		/* rule 289 */
	 tx196 = tx18*y3; 		/* rule 290 */
	 tx197 = tx121*z1; 		/* rule 291 */
	 tx198 = tx122*z1; 		/* rule 292 */
	 tx199 = tx154*z1; 		/* rule 293 */
	 tx200 = tx155*z1; 		/* rule 294 */
	 tx201 = y2*y3*z2; 		/* rule 295 */
	 tx202 = x2*x3*z3; 		/* rule 296 */
	 tx203 = tx181 + tx182 + tx183 + tx184 + tx185 + tx186 + tx187 + tx188 + tx189 + tx190 + tx191 + tx192 + tx193 + tx194 + tx195 + tx196 + tx197 + tx198 + tx199 + tx200 + tx201 + tx202; 		/* rule 297 */
	 tx204 = tx149*tx151*tx22*tx61; 		/* rule 298 */
	 tx205 = 6.*EraseLinearDihedral*K*TWOPI*tx151*tx19*tx203*tx65*tx70*tx75*tx76; 		/* rule 299 */
	 dhz1z1 = tx204 + tx205; 		/* rule 300 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 tx206 = power2(y1); 		/* rule 302 */
	 tx207 = power2(z1); 		/* rule 303 */
	 tx208 = -tx133; 		/* rule 304 */
	 tx209 = -tx134; 		/* rule 305 */
	 tx210 = -tx135; 		/* rule 306 */
	 tx211 = -tx136; 		/* rule 307 */
	 tx212 = 2.*tx139; 		/* rule 308 */
	 tx213 = 2.*tx143; 		/* rule 309 */
	 tx214 = tx206*tx49; 		/* rule 310 */
	 tx215 = tx207*tx49; 		/* rule 311 */
	 tx216 = tx3*tx50; 		/* rule 312 */
	 tx217 = tx12*tx51; 		/* rule 313 */
	 tx218 = tx51*tx7; 		/* rule 314 */
	 tx219 = tx2*tx79; 		/* rule 315 */
	 tx220 = tx6*tx79; 		/* rule 316 */
	 tx221 = tx122*x2; 		/* rule 317 */
	 tx222 = tx124*x2; 		/* rule 318 */
	 tx223 = tx206*x2; 		/* rule 319 */
	 tx224 = tx207*x2; 		/* rule 320 */
	 tx225 = tx121*x3; 		/* rule 321 */
	 tx226 = tx123*x3; 		/* rule 322 */
	 tx227 = -(tx3*y1); 		/* rule 323 */
	 tx228 = x2*y1*y2; 		/* rule 324 */
	 tx229 = tx1*y3; 		/* rule 325 */
	 tx230 = tx2*y3; 		/* rule 326 */
	 tx231 = 2.*tx3*y3; 		/* rule 327 */
	 tx232 = tx4*y3; 		/* rule 328 */
	 tx233 = x1*y1*y3; 		/* rule 329 */
	 tx234 = tx11*z1; 		/* rule 330 */
	 tx235 = tx15*z1; 		/* rule 331 */
	 tx236 = 2.*tx15*z2; 		/* rule 332 */
	 tx237 = tx16*z2; 		/* rule 333 */
	 tx238 = tx7*z2; 		/* rule 334 */
	 tx239 = tx8*z2; 		/* rule 335 */
	 tx240 = tx11*z3; 		/* rule 336 */
	 tx241 = x3*z1*z3; 		/* rule 337 */
	 tx242 = tx125 + tx126 + tx132 + tx133 + tx135 + tx138 + tx140 + tx148 + tx212 + tx213 + tx214 + tx215 + tx216 + tx218 + tx223 + tx224 + tx227 + tx229 + tx233 + tx234 + tx235 + tx240; 		/* rule 338 */
	 tx243 = tx208 + tx209 + tx210 + tx211 + tx217 + tx218 + tx219 + tx220 + tx221 + tx222 + tx225 + tx226 + tx228 + tx229 + tx230 + tx231 + tx232 + tx236 + tx237 + tx238 + tx239 + tx241; 		/* rule 339 */
	 tx244 = power2(y4); 		/* rule 340 */
	 tx245 = power2(z4); 		/* rule 341 */
	 tx246 = -2.*tx139; 		/* rule 342 */
	 tx247 = -2.*tx143; 		/* rule 343 */
	 tx248 = -tx235; 		/* rule 344 */
	 tx249 = tx26*tx50; 		/* rule 345 */
	 tx250 = tx29*tx50; 		/* rule 346 */
	 tx251 = tx21*tx242*tx61; 		/* rule 347 */
	 tx252 = tx206*tx78; 		/* rule 348 */
	 tx253 = tx207*tx78; 		/* rule 349 */
	 tx254 = tx244*tx78; 		/* rule 350 */
	 tx255 = tx245*tx78; 		/* rule 351 */
	 tx256 = tx34*tx80; 		/* rule 352 */
	 tx257 = tx21*tx243*tx93; 		/* rule 353 */
	 tx258 = tx206*x3; 		/* rule 354 */
	 tx259 = tx207*x3; 		/* rule 355 */
	 tx260 = tx244*x3; 		/* rule 356 */
	 tx261 = tx245*x3; 		/* rule 357 */
	 tx262 = tx121*x4; 		/* rule 358 */
	 tx263 = tx122*x4; 		/* rule 359 */
	 tx264 = tx123*x4; 		/* rule 360 */
	 tx265 = tx124*x4; 		/* rule 361 */
	 tx266 = tx3*y1; 		/* rule 362 */
	 tx267 = tx5*y1; 		/* rule 363 */
	 tx268 = tx28*y2; 		/* rule 364 */
	 tx269 = tx29*y2; 		/* rule 365 */
	 tx270 = -2.*tx26*y3; 		/* rule 366 */
	 tx271 = tx3*y3; 		/* rule 367 */
	 tx272 = tx26*y4; 		/* rule 368 */
	 tx273 = tx27*y4; 		/* rule 369 */
	 tx274 = tx4*y4; 		/* rule 370 */
	 tx275 = tx6*y4; 		/* rule 371 */
	 tx276 = 2.*tx275; 		/* rule 372 */
	 tx277 = tx15*z2; 		/* rule 373 */
	 tx278 = -2.*tx32*z2; 		/* rule 374 */
	 tx279 = x1*z1*z2; 		/* rule 375 */
	 tx280 = tx30*z3; 		/* rule 376 */
	 tx281 = tx34*z3; 		/* rule 377 */
	 tx282 = 2.*tx281; 		/* rule 378 */
	 tx283 = tx35*z3; 		/* rule 379 */
	 tx284 = -(tx12*z4); 		/* rule 380 */
	 tx285 = tx12*z4; 		/* rule 381 */
	 tx286 = -(tx32*z4); 		/* rule 382 */
	 tx287 = x4*z2*z4; 		/* rule 383 */
	 tx288 = tx276; 		/* rule 384 */
	 tx289 = tx282; 		/* rule 385 */
	 tx290 = tx251 + tx257; 		/* rule 386 */
	 tx291 = tx137 + tx146 + tx208 + tx210 + tx220 + tx225 + tx226 + tx228 + tx237 + tx238 + tx246 + tx247 + tx248 + tx252 + tx253 + tx258 + tx259 + tx266 + tx267 + tx271 + tx277 + tx279; 		/* rule 387 */
	 tx292 = tx125 + tx126 + tx130 + tx131 + tx140 + tx142 + tx145 + tx148 + tx250 + tx256 + tx262 + tx263 + tx264 + tx265 + tx268 + tx269 + tx270 + tx275 + tx278 + tx281 + tx283 + tx285; 		/* rule 388 */
	 tx293 = tx130 + tx131 + tx142 + tx145 + tx249 + tx250 + tx254 + tx255 + tx260 + tx261 + tx263 + tx265 + tx272 + tx273 + tx274 + tx280 + tx283 + tx284 + tx286 + tx287 + tx288 + tx289; 		/* rule 389 */
	 tx294 = power2(tx90); 		/* rule 390 */
	 tx917 = tx61; 		/* rule 391 */
	 tx918 = tx91; 		/* rule 392 */
	 tx919 = reciprocal(tx917); 		/* rule 393 */
	 tx295 = tx918*tx919; 		/* rule 394 */
	 tx296 = 2.*tx150*tx151*tx21*tx65; 		/* rule 395 */
	 tx297 = tx150*tx21*tx295*tx77; 		/* rule 396 */
	 tx298 = -(tx292*tx295*tx39*tx90); 		/* rule 397 */
	 tx299 = -2.*tx151*tx290*tx918; 		/* rule 398 */
	 tx300 = -2.*tx293*tx294*tx39*tx918*tx92; 		/* rule 399 */
	 tx301 = 2.*tx151*tx21*tx291*tx918*tx93; 		/* rule 400 */
	 tx302 = tx296 + tx297 + tx298 + tx299 + tx300 + tx301; 		/* rule 401 */
	 tx303 = power2(tx97); 		/* rule 402 */
	 tx304 = DEImproperRestraintDPhi*tx302; 		/* rule 403 */
	 tx305 = tx149*tx303; 		/* rule 404 */
	 dhx2x2 = tx304 + tx305; 		/* rule 405 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 tx306 = power2(x1); 		/* rule 407 */
	 tx307 = 2.*tx162; 		/* rule 408 */
	 tx308 = -tx167; 		/* rule 409 */
	 tx309 = -tx168; 		/* rule 410 */
	 tx310 = -tx169; 		/* rule 411 */
	 tx311 = -tx170; 		/* rule 412 */
	 tx312 = 2.*tx173; 		/* rule 413 */
	 tx313 = tx3*tx49; 		/* rule 414 */
	 tx314 = tx207*tx50; 		/* rule 415 */
	 tx315 = tx306*tx50; 		/* rule 416 */
	 tx316 = tx13*tx51; 		/* rule 417 */
	 tx317 = tx49*tx6; 		/* rule 418 */
	 tx318 = tx10*tx80; 		/* rule 419 */
	 tx319 = tx18*tx80; 		/* rule 420 */
	 tx320 = tx1*x1; 		/* rule 421 */
	 tx321 = tx2*x1; 		/* rule 422 */
	 tx322 = 2.*tx2*x2; 		/* rule 423 */
	 tx323 = tx3*x2; 		/* rule 424 */
	 tx324 = tx4*x2; 		/* rule 425 */
	 tx325 = tx5*x2; 		/* rule 426 */
	 tx326 = tx1*x3; 		/* rule 427 */
	 tx327 = tx124*y2; 		/* rule 428 */
	 tx328 = tx155*y2; 		/* rule 429 */
	 tx329 = tx207*y2; 		/* rule 430 */
	 tx330 = tx306*y2; 		/* rule 431 */
	 tx331 = tx123*y3; 		/* rule 432 */
	 tx332 = tx154*y3; 		/* rule 433 */
	 tx333 = x1*x3*y3; 		/* rule 434 */
	 tx334 = -(tx13*z1); 		/* rule 435 */
	 tx335 = y2*z1*z2; 		/* rule 436 */
	 tx336 = tx10*z3; 		/* rule 437 */
	 tx337 = 2.*tx13*z3; 		/* rule 438 */
	 tx338 = tx14*z3; 		/* rule 439 */
	 tx339 = tx9*z3; 		/* rule 440 */
	 tx340 = y1*z1*z3; 		/* rule 441 */
	 tx341 = tx308 + tx309 + tx310 + tx311 + tx313 + tx317 + tx318 + tx319 + tx322 + tx323 + tx324 + tx325 + tx327 + tx328 + tx331 + tx332 + tx333 + tx335 + tx336 + tx337 + tx338 + tx339; 		/* rule 442 */
	 tx342 = tx156 + tx157 + tx159 + tx167 + tx169 + tx171 + tx174 + tx175 + tx307 + tx312 + tx313 + tx314 + tx315 + tx316 + tx320 + tx321 + tx326 + tx329 + tx330 + tx334 + tx339 + tx340; 		/* rule 443 */
	 tx343 = power2(x4); 		/* rule 444 */
	 tx344 = -2.*tx162; 		/* rule 445 */
	 tx345 = -2.*tx173; 		/* rule 446 */
	 tx346 = -tx321; 		/* rule 447 */
	 tx347 = tx31*tx51; 		/* rule 448 */
	 tx348 = tx37*tx51; 		/* rule 449 */
	 tx349 = tx20*tx342*tx917; 		/* rule 450 */
	 tx350 = tx26*tx78; 		/* rule 451 */
	 tx351 = tx207*tx79; 		/* rule 452 */
	 tx352 = tx245*tx79; 		/* rule 453 */
	 tx353 = tx306*tx79; 		/* rule 454 */
	 tx354 = tx343*tx79; 		/* rule 455 */
	 tx355 = tx20*tx341*tx93; 		/* rule 456 */
	 tx356 = tx2*x2; 		/* rule 457 */
	 tx357 = -2.*tx29*x2; 		/* rule 458 */
	 tx358 = tx26*x3; 		/* rule 459 */
	 tx359 = 2.*tx358; 		/* rule 460 */
	 tx360 = tx27*x3; 		/* rule 461 */
	 tx361 = tx28*x3; 		/* rule 462 */
	 tx362 = -(tx29*x4); 		/* rule 463 */
	 tx363 = -(tx6*x4); 		/* rule 464 */
	 tx364 = tx6*x4; 		/* rule 465 */
	 tx365 = x1*x2*y1; 		/* rule 466 */
	 tx366 = tx207*y3; 		/* rule 467 */
	 tx367 = tx245*y3; 		/* rule 468 */
	 tx368 = tx306*y3; 		/* rule 469 */
	 tx369 = tx343*y3; 		/* rule 470 */
	 tx370 = tx123*y4; 		/* rule 471 */
	 tx371 = tx124*y4; 		/* rule 472 */
	 tx372 = tx154*y4; 		/* rule 473 */
	 tx373 = tx155*y4; 		/* rule 474 */
	 tx374 = x2*x4*y4; 		/* rule 475 */
	 tx375 = tx13*z1; 		/* rule 476 */
	 tx376 = tx17*z1; 		/* rule 477 */
	 tx377 = tx36*z2; 		/* rule 478 */
	 tx378 = tx37*z2; 		/* rule 479 */
	 tx379 = tx13*z3; 		/* rule 480 */
	 tx380 = -2.*tx31*z3; 		/* rule 481 */
	 tx381 = tx14*z4; 		/* rule 482 */
	 tx382 = tx18*z4; 		/* rule 483 */
	 tx383 = 2.*tx382; 		/* rule 484 */
	 tx384 = tx31*z4; 		/* rule 485 */
	 tx385 = tx33*z4; 		/* rule 486 */
	 tx386 = tx359; 		/* rule 487 */
	 tx387 = tx383; 		/* rule 488 */
	 tx388 = tx349 + tx355; 		/* rule 489 */
	 tx389 = tx164 + tx172 + tx308 + tx310 + tx319 + tx323 + tx324 + tx331 + tx332 + tx335 + tx344 + tx345 + tx346 + tx351 + tx353 + tx356 + tx365 + tx366 + tx368 + tx375 + tx376 + tx379; 		/* rule 490 */
	 tx390 = tx156 + tx157 + tx160 + tx161 + tx166 + tx171 + tx174 + tx177 + tx348 + tx350 + tx357 + tx358 + tx360 + tx364 + tx370 + tx371 + tx372 + tx373 + tx377 + tx378 + tx380 + tx382; 		/* rule 491 */
	 tx391 = tx160 + tx161 + tx166 + tx177 + tx347 + tx348 + tx352 + tx354 + tx360 + tx361 + tx362 + tx363 + tx367 + tx369 + tx371 + tx373 + tx374 + tx381 + tx384 + tx385 + tx386 + tx387; 		/* rule 492 */
	 tx392 = 2.*tx151*tx178*tx20*tx65; 		/* rule 493 */
	 tx393 = tx178*tx20*tx295*tx77; 		/* rule 494 */
	 tx394 = -(tx295*tx38*tx390*tx90); 		/* rule 495 */
	 tx395 = -2.*tx151*tx388*tx918; 		/* rule 496 */
	 tx396 = -2.*tx294*tx38*tx391*tx918*tx92; 		/* rule 497 */
	 tx397 = 2.*tx151*tx20*tx389*tx918*tx93; 		/* rule 498 */
	 tx398 = power2(tx101); 		/* rule 499 */
	 tx399 = tx392 + tx393 + tx394 + tx395 + tx396 + tx397; 		/* rule 500 */
	 tx400 = tx149*tx398; 		/* rule 501 */
	 tx401 = DEImproperRestraintDPhi*tx399; 		/* rule 502 */
	 dhy2y2 = tx400 + tx401; 		/* rule 503 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 tx402 = 2.*tx190; 		/* rule 505 */
	 tx403 = 2.*tx192; 		/* rule 506 */
	 tx404 = -tx197; 		/* rule 507 */
	 tx405 = -tx198; 		/* rule 508 */
	 tx406 = -tx199; 		/* rule 509 */
	 tx407 = -tx200; 		/* rule 510 */
	 tx408 = tx13*tx50; 		/* rule 511 */
	 tx409 = tx18*tx50; 		/* rule 512 */
	 tx410 = tx206*tx51; 		/* rule 513 */
	 tx411 = tx306*tx51; 		/* rule 514 */
	 tx412 = tx49*tx7; 		/* rule 515 */
	 tx413 = tx12*tx78; 		/* rule 516 */
	 tx414 = tx15*tx78; 		/* rule 517 */
	 tx415 = -(tx7*x1); 		/* rule 518 */
	 tx416 = tx11*x3; 		/* rule 519 */
	 tx417 = tx15*x3; 		/* rule 520 */
	 tx418 = tx16*x3; 		/* rule 521 */
	 tx419 = 2.*tx7*x3; 		/* rule 522 */
	 tx420 = tx10*y1; 		/* rule 523 */
	 tx421 = tx9*y1; 		/* rule 524 */
	 tx422 = 2.*tx10*y2; 		/* rule 525 */
	 tx423 = tx13*y2; 		/* rule 526 */
	 tx424 = tx14*y2; 		/* rule 527 */
	 tx425 = tx17*y2; 		/* rule 528 */
	 tx426 = tx9*y3; 		/* rule 529 */
	 tx427 = x1*x3*z1; 		/* rule 530 */
	 tx428 = tx122*z2; 		/* rule 531 */
	 tx429 = tx155*z2; 		/* rule 532 */
	 tx430 = tx206*z2; 		/* rule 533 */
	 tx431 = tx306*z2; 		/* rule 534 */
	 tx432 = x1*x2*z2; 		/* rule 535 */
	 tx433 = tx121*z3; 		/* rule 536 */
	 tx434 = tx154*z3; 		/* rule 537 */
	 tx435 = y1*y3*z3; 		/* rule 538 */
	 tx436 = tx182 + tx183 + tx184 + tx188 + tx189 + tx197 + tx199 + tx201 + tx402 + tx403 + tx408 + tx410 + tx411 + tx412 + tx415 + tx416 + tx420 + tx421 + tx426 + tx427 + tx430 + tx431; 		/* rule 539 */
	 tx437 = tx404 + tx405 + tx406 + tx407 + tx408 + tx409 + tx413 + tx414 + tx416 + tx417 + tx418 + tx419 + tx422 + tx423 + tx424 + tx425 + tx428 + tx429 + tx432 + tx433 + tx434 + tx435; 		/* rule 540 */
	 tx438 = -2.*tx190; 		/* rule 541 */
	 tx439 = -2.*tx192; 		/* rule 542 */
	 tx440 = -tx420; 		/* rule 543 */
	 tx441 = tx32*tx49; 		/* rule 544 */
	 tx442 = tx34*tx49; 		/* rule 545 */
	 tx443 = tx19*tx436*tx917; 		/* rule 546 */
	 tx444 = tx31*tx79; 		/* rule 547 */
	 tx445 = tx206*tx80; 		/* rule 548 */
	 tx446 = tx244*tx80; 		/* rule 549 */
	 tx447 = tx306*tx80; 		/* rule 550 */
	 tx448 = tx343*tx80; 		/* rule 551 */
	 tx449 = tx19*tx437*tx93; 		/* rule 552 */
	 tx450 = tx7*x1; 		/* rule 553 */
	 tx451 = tx8*x1; 		/* rule 554 */
	 tx452 = tx30*x2; 		/* rule 555 */
	 tx453 = tx32*x2; 		/* rule 556 */
	 tx454 = -2.*tx34*x3; 		/* rule 557 */
	 tx455 = tx7*x3; 		/* rule 558 */
	 tx456 = tx12*x4; 		/* rule 559 */
	 tx457 = 2.*tx456; 		/* rule 560 */
	 tx458 = tx16*x4; 		/* rule 561 */
	 tx459 = tx34*x4; 		/* rule 562 */
	 tx460 = tx35*x4; 		/* rule 563 */
	 tx461 = tx10*y2; 		/* rule 564 */
	 tx462 = -2.*tx37*y2; 		/* rule 565 */
	 tx463 = tx31*y3; 		/* rule 566 */
	 tx464 = 2.*tx463; 		/* rule 567 */
	 tx465 = tx33*y3; 		/* rule 568 */
	 tx466 = tx36*y3; 		/* rule 569 */
	 tx467 = -(tx18*y4); 		/* rule 570 */
	 tx468 = tx18*y4; 		/* rule 571 */
	 tx469 = -(tx37*y4); 		/* rule 572 */
	 tx470 = y1*y2*z1; 		/* rule 573 */
	 tx471 = tx206*z3; 		/* rule 574 */
	 tx472 = tx244*z3; 		/* rule 575 */
	 tx473 = tx306*z3; 		/* rule 576 */
	 tx474 = tx343*z3; 		/* rule 577 */
	 tx475 = tx121*z4; 		/* rule 578 */
	 tx476 = tx122*z4; 		/* rule 579 */
	 tx477 = tx154*z4; 		/* rule 580 */
	 tx478 = tx155*z4; 		/* rule 581 */
	 tx479 = y2*y4*z4; 		/* rule 582 */
	 tx480 = tx457; 		/* rule 583 */
	 tx481 = tx464; 		/* rule 584 */
	 tx482 = tx443 + tx449; 		/* rule 585 */
	 tx483 = tx187 + tx194 + tx404 + tx406 + tx413 + tx423 + tx424 + tx432 + tx433 + tx434 + tx438 + tx439 + tx440 + tx445 + tx447 + tx450 + tx451 + tx455 + tx461 + tx470 + tx471 + tx473; 		/* rule 586 */
	 tx484 = tx182 + tx183 + tx185 + tx186 + tx189 + tx196 + tx201 + tx202 + tx441 + tx444 + tx452 + tx453 + tx454 + tx456 + tx462 + tx463 + tx465 + tx468 + tx475 + tx476 + tx477 + tx478; 		/* rule 587 */
	 tx485 = tx185 + tx186 + tx196 + tx202 + tx441 + tx442 + tx446 + tx448 + tx458 + tx459 + tx460 + tx465 + tx466 + tx467 + tx469 + tx472 + tx474 + tx476 + tx478 + tx479 + tx480 + tx481; 		/* rule 588 */
	 tx486 = 2.*tx151*tx19*tx203*tx65; 		/* rule 589 */
	 tx487 = tx19*tx203*tx295*tx77; 		/* rule 590 */
	 tx488 = -(tx295*tx40*tx484*tx90); 		/* rule 591 */
	 tx489 = -2.*tx151*tx482*tx918; 		/* rule 592 */
	 tx490 = -2.*tx294*tx40*tx485*tx918*tx92; 		/* rule 593 */
	 tx491 = 2.*tx151*tx19*tx483*tx918*tx93; 		/* rule 594 */
	 tx492 = power2(tx105); 		/* rule 595 */
	 tx493 = tx486 + tx487 + tx488 + tx489 + tx490 + tx491; 		/* rule 596 */
	 tx494 = tx149*tx492; 		/* rule 597 */
	 tx495 = DEImproperRestraintDPhi*tx493; 		/* rule 598 */
	 dhz2z2 = tx494 + tx495; 		/* rule 599 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	 tx496 = -tx262; 		/* rule 601 */
	 tx497 = -tx263; 		/* rule 602 */
	 tx498 = -tx264; 		/* rule 603 */
	 tx499 = -tx265; 		/* rule 604 */
	 tx500 = -tx272; 		/* rule 605 */
	 tx501 = -2.*tx275; 		/* rule 606 */
	 tx502 = -tx275; 		/* rule 607 */
	 tx503 = -2.*tx281; 		/* rule 608 */
	 tx504 = tx244*tx49; 		/* rule 609 */
	 tx505 = tx245*tx49; 		/* rule 610 */
	 tx506 = tx244*x2; 		/* rule 611 */
	 tx507 = tx245*x2; 		/* rule 612 */
	 tx508 = tx26*y3; 		/* rule 613 */
	 tx509 = 2.*tx508; 		/* rule 614 */
	 tx510 = tx29*y3; 		/* rule 615 */
	 tx511 = x2*y2*y4; 		/* rule 616 */
	 tx512 = x4*y3*y4; 		/* rule 617 */
	 tx513 = tx32*z2; 		/* rule 618 */
	 tx514 = 2.*tx513; 		/* rule 619 */
	 tx515 = tx34*z2; 		/* rule 620 */
	 tx516 = tx16*z4; 		/* rule 621 */
	 tx517 = tx30*z4; 		/* rule 622 */
	 tx518 = tx32*z4; 		/* rule 623 */
	 tx519 = x3*z3*z4; 		/* rule 624 */
	 tx520 = tx509; 		/* rule 625 */
	 tx521 = tx514; 		/* rule 626 */
	 tx522 = tx217 + tx221 + tx222 + tx232 + tx269 + tx285 + tx497 + tx499 + tx500 + tx501 + tx503 + tx504 + tx505 + tx506 + tx507 + tx508 + tx510 + tx512 + tx513 + tx517 + tx518 + tx519; 		/* rule 627 */
	 tx523 = tx217 + tx220 + tx221 + tx222 + tx225 + tx226 + tx232 + tx237 + tx274 + tx284 + tx496 + tx497 + tx498 + tx499 + tx502 + tx510 + tx511 + tx515 + tx516 + tx519 + tx520 + tx521; 		/* rule 628 */
	 tx524 = tx39*tx522*tx917; 		/* rule 629 */
	 tx525 = tx39*tx523*tx92; 		/* rule 630 */
	 tx526 = tx524 + tx525; 		/* rule 631 */
	 tx527 = -2.*tx292*tx294*tx39*tx65; 		/* rule 632 */
	 tx528 = -2.*tx294*tx526*tx918; 		/* rule 633 */
	 tx529 = power2(tx109); 		/* rule 634 */
	 tx530 = tx297 + tx298 + tx300 + tx301 + tx527 + tx528; 		/* rule 635 */
	 tx531 = tx149*tx529; 		/* rule 636 */
	 tx532 = DEImproperRestraintDPhi*tx530; 		/* rule 637 */
	 dhx3x3 = tx531 + tx532; 		/* rule 638 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 0, dhx3x3);
	 tx533 = -2.*tx358; 		/* rule 640 */
	 tx534 = -tx370; 		/* rule 641 */
	 tx535 = -tx371; 		/* rule 642 */
	 tx536 = -tx372; 		/* rule 643 */
	 tx537 = -tx373; 		/* rule 644 */
	 tx538 = -2.*tx382; 		/* rule 645 */
	 tx539 = -tx382; 		/* rule 646 */
	 tx540 = -tx384; 		/* rule 647 */
	 tx541 = tx245*tx50; 		/* rule 648 */
	 tx542 = tx343*tx50; 		/* rule 649 */
	 tx543 = tx26*x2; 		/* rule 650 */
	 tx544 = tx29*x2; 		/* rule 651 */
	 tx545 = 2.*tx544; 		/* rule 652 */
	 tx546 = tx28*x4; 		/* rule 653 */
	 tx547 = tx29*x4; 		/* rule 654 */
	 tx548 = tx4*x4; 		/* rule 655 */
	 tx549 = tx245*y2; 		/* rule 656 */
	 tx550 = tx343*y2; 		/* rule 657 */
	 tx551 = x3*x4*y3; 		/* rule 658 */
	 tx552 = tx31*z3; 		/* rule 659 */
	 tx553 = 2.*tx552; 		/* rule 660 */
	 tx554 = tx37*z3; 		/* rule 661 */
	 tx555 = y2*z2*z4; 		/* rule 662 */
	 tx556 = y4*z3*z4; 		/* rule 663 */
	 tx557 = tx545; 		/* rule 664 */
	 tx558 = tx553; 		/* rule 665 */
	 tx559 = tx317 + tx327 + tx328 + tx338 + tx364 + tx378 + tx533 + tx535 + tx537 + tx538 + tx540 + tx541 + tx542 + tx544 + tx546 + tx547 + tx549 + tx550 + tx551 + tx552 + tx554 + tx556; 		/* rule 666 */
	 tx560 = tx317 + tx319 + tx324 + tx327 + tx328 + tx331 + tx332 + tx338 + tx363 + tx381 + tx534 + tx535 + tx536 + tx537 + tx539 + tx543 + tx548 + tx551 + tx554 + tx555 + tx557 + tx558; 		/* rule 667 */
	 tx561 = tx38*tx559*tx917; 		/* rule 668 */
	 tx562 = tx38*tx560*tx92; 		/* rule 669 */
	 tx563 = tx561 + tx562; 		/* rule 670 */
	 tx564 = -2.*tx294*tx38*tx390*tx65; 		/* rule 671 */
	 tx565 = -2.*tx294*tx563*tx918; 		/* rule 672 */
	 tx566 = power2(tx113); 		/* rule 673 */
	 tx567 = tx393 + tx394 + tx396 + tx397 + tx564 + tx565; 		/* rule 674 */
	 tx568 = tx149*tx566; 		/* rule 675 */
	 tx569 = DEImproperRestraintDPhi*tx567; 		/* rule 676 */
	 dhy3y3 = tx568 + tx569; 		/* rule 677 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 1, dhy3y3);
	 tx570 = -2.*tx456; 		/* rule 679 */
	 tx571 = -tx456; 		/* rule 680 */
	 tx572 = -tx459; 		/* rule 681 */
	 tx573 = -2.*tx463; 		/* rule 682 */
	 tx574 = -tx475; 		/* rule 683 */
	 tx575 = -tx476; 		/* rule 684 */
	 tx576 = -tx477; 		/* rule 685 */
	 tx577 = -tx478; 		/* rule 686 */
	 tx578 = tx244*tx51; 		/* rule 687 */
	 tx579 = tx343*tx51; 		/* rule 688 */
	 tx580 = tx32*x3; 		/* rule 689 */
	 tx581 = tx34*x3; 		/* rule 690 */
	 tx582 = 2.*tx581; 		/* rule 691 */
	 tx583 = tx31*y2; 		/* rule 692 */
	 tx584 = tx37*y2; 		/* rule 693 */
	 tx585 = 2.*tx584; 		/* rule 694 */
	 tx586 = tx14*y4; 		/* rule 695 */
	 tx587 = tx36*y4; 		/* rule 696 */
	 tx588 = tx37*y4; 		/* rule 697 */
	 tx589 = tx244*z2; 		/* rule 698 */
	 tx590 = tx343*z2; 		/* rule 699 */
	 tx591 = x2*x4*z2; 		/* rule 700 */
	 tx592 = y3*y4*z3; 		/* rule 701 */
	 tx593 = x3*x4*z4; 		/* rule 702 */
	 tx594 = tx582; 		/* rule 703 */
	 tx595 = tx585; 		/* rule 704 */
	 tx596 = tx409 + tx418 + tx428 + tx429 + tx453 + tx468 + tx570 + tx572 + tx573 + tx575 + tx577 + tx578 + tx579 + tx580 + tx581 + tx584 + tx587 + tx588 + tx589 + tx590 + tx592 + tx593; 		/* rule 705 */
	 tx597 = tx409 + tx413 + tx418 + tx424 + tx428 + tx429 + tx433 + tx434 + tx458 + tx467 + tx571 + tx574 + tx575 + tx576 + tx577 + tx580 + tx583 + tx586 + tx591 + tx592 + tx594 + tx595; 		/* rule 706 */
	 tx598 = tx40*tx596*tx917; 		/* rule 707 */
	 tx599 = tx40*tx597*tx92; 		/* rule 708 */
	 tx600 = tx598 + tx599; 		/* rule 709 */
	 tx601 = -2.*tx294*tx40*tx484*tx65; 		/* rule 710 */
	 tx602 = -2.*tx294*tx600*tx918; 		/* rule 711 */
	 tx603 = power2(tx117); 		/* rule 712 */
	 tx604 = tx487 + tx488 + tx490 + tx491 + tx601 + tx602; 		/* rule 713 */
	 tx605 = tx149*tx603; 		/* rule 714 */
	 tx606 = DEImproperRestraintDPhi*tx604; 		/* rule 715 */
	 dhz3z3 = tx605 + tx606; 		/* rule 716 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I3, 2, dhz3z3);
	 tx607 = tx149*tx294*tx42*tx917; 		/* rule 718 */
	 tx608 = -6.*EraseLinearDihedral*K*TWOPI*tx292*tx294*tx39*tx65*tx70*tx75*tx76; 		/* rule 719 */
	 dhx4x4 = tx607 + tx608; 		/* rule 720 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 0, dhx4x4);
	 tx609 = tx149*tx294*tx41*tx917; 		/* rule 722 */
	 tx610 = -6.*EraseLinearDihedral*K*TWOPI*tx294*tx38*tx390*tx65*tx70*tx75*tx76; 		/* rule 723 */
	 dhy4y4 = tx609 + tx610; 		/* rule 724 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I4, 1, I4, 1, dhy4y4);
	 tx611 = tx149*tx294*tx43*tx917; 		/* rule 726 */
	 tx612 = -6.*EraseLinearDihedral*K*TWOPI*tx294*tx40*tx484*tx65*tx70*tx75*tx76; 		/* rule 727 */
	 dhz4z4 = tx611 + tx612; 		/* rule 728 */
	IMPROPER_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(I4, 2, I4, 2, dhz4z4);
	#ifdef IMPROPER_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN //[
	if (calcOffDiagonalHessian) { 
	 tx613 = tx150*tx20; 		/* rule 732 */
	 tx614 = tx178*tx21; 		/* rule 733 */
	 tx615 = tx613 + tx614; 		/* rule 734 */
	 tx616 = tx149*tx151*tx20*tx21*tx917; 		/* rule 735 */
	 tx617 = DEImproperRestraintDPhi*tx151*tx615*tx65; 		/* rule 736 */
	 ohx1y1 = tx616 + tx617; 		/* rule 737 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tx618 = tx150*tx19; 		/* rule 739 */
	 tx619 = tx203*tx21; 		/* rule 740 */
	 tx620 = tx618 + tx619; 		/* rule 741 */
	 tx621 = tx149*tx151*tx19*tx21*tx917; 		/* rule 742 */
	 tx622 = DEImproperRestraintDPhi*tx151*tx620*tx65; 		/* rule 743 */
	 ohx1z1 = tx621 + tx622; 		/* rule 744 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx623 = -2.*tx150*tx151*tx21*tx65; 		/* rule 746 */
	 tx624 = tx151*tx290*tx918; 		/* rule 747 */
	 tx625 = tx623 + tx624; 		/* rule 748 */
	 tx626 = DEImproperRestraintDPhi*tx625; 		/* rule 749 */
	 tx627 = -(tx149*tx94*tx97); 		/* rule 750 */
	 ohx1x2 = tx626 + tx627; 		/* rule 751 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx628 = tx20*tx242*tx917; 		/* rule 753 */
	 tx629 = tx21*tx341*tx93; 		/* rule 754 */
	 tx630 = tx628 + tx629; 		/* rule 755 */
	 tx631 = -(tx151*tx615*tx65); 		/* rule 756 */
	 tx632 = tx151*tx630*tx918; 		/* rule 757 */
	 tx633 = tx631 + tx632; 		/* rule 758 */
	 tx634 = DEImproperRestraintDPhi*tx633; 		/* rule 759 */
	 tx635 = -(tx101*tx149*tx94); 		/* rule 760 */
	 ohx1y2 = tx634 + tx635; 		/* rule 761 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx636 = tx19*tx242*tx917; 		/* rule 763 */
	 tx637 = tx21*tx437*tx93; 		/* rule 764 */
	 tx638 = tx636 + tx637; 		/* rule 765 */
	 tx639 = -(tx151*tx620*tx65); 		/* rule 766 */
	 tx640 = tx151*tx638*tx918; 		/* rule 767 */
	 tx641 = tx639 + tx640; 		/* rule 768 */
	 tx642 = DEImproperRestraintDPhi*tx641; 		/* rule 769 */
	 tx643 = -(tx105*tx149*tx94); 		/* rule 770 */
	 ohx1z2 = tx642 + tx643; 		/* rule 771 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tx644 = -3*EraseLinearDihedral*K*TWOPI*tx624*tx70*tx75*tx76; 		/* rule 773 */
	 tx645 = -(tx109*tx149*tx94); 		/* rule 774 */
	 ohx1x3 = tx644 + tx645; 		/* rule 775 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 0, ohx1x3);
	 tx646 = -3*EraseLinearDihedral*K*TWOPI*tx632*tx70*tx75*tx76; 		/* rule 777 */
	 tx647 = -(tx113*tx149*tx94); 		/* rule 778 */
	 ohx1y3 = tx646 + tx647; 		/* rule 779 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 1, ohx1y3);
	 tx648 = -3*EraseLinearDihedral*K*TWOPI*tx640*tx70*tx75*tx76; 		/* rule 781 */
	 tx649 = -(tx117*tx149*tx94); 		/* rule 782 */
	 ohx1z3 = tx648 + tx649; 		/* rule 783 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I3, 2, ohx1z3);
	 ohx1x4 = -(tx149*tx46*tx77*tx90*tx917); 		/* rule 785 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 0, ohx1x4);
	 ohx1y4 = -(tx149*tx21*tx38*tx77*tx90*tx917); 		/* rule 787 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 1, ohx1y4);
	 ohx1z4 = -(tx149*tx21*tx40*tx77*tx90*tx917); 		/* rule 789 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I4, 2, ohx1z4);
	 tx650 = tx178*tx19; 		/* rule 791 */
	 tx651 = tx20*tx203; 		/* rule 792 */
	 tx652 = tx650 + tx651; 		/* rule 793 */
	 tx653 = tx149*tx151*tx19*tx20*tx917; 		/* rule 794 */
	 tx654 = DEImproperRestraintDPhi*tx151*tx65*tx652; 		/* rule 795 */
	 ohy1z1 = tx653 + tx654; 		/* rule 796 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 tx655 = tx21*tx342*tx917; 		/* rule 798 */
	 tx656 = tx20*tx243*tx93; 		/* rule 799 */
	 tx657 = tx655 + tx656; 		/* rule 800 */
	 tx658 = tx151*tx657*tx918; 		/* rule 801 */
	 tx659 = tx631 + tx658; 		/* rule 802 */
	 tx660 = DEImproperRestraintDPhi*tx659; 		/* rule 803 */
	 tx661 = -(tx149*tx97*tx98); 		/* rule 804 */
	 ohy1x2 = tx660 + tx661; 		/* rule 805 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx662 = -2.*tx151*tx178*tx20*tx65; 		/* rule 807 */
	 tx663 = tx151*tx388*tx918; 		/* rule 808 */
	 tx664 = tx662 + tx663; 		/* rule 809 */
	 tx665 = DEImproperRestraintDPhi*tx664; 		/* rule 810 */
	 tx666 = -(tx101*tx149*tx98); 		/* rule 811 */
	 ohy1y2 = tx665 + tx666; 		/* rule 812 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx667 = tx19*tx342*tx917; 		/* rule 814 */
	 tx668 = tx20*tx437*tx93; 		/* rule 815 */
	 tx669 = tx667 + tx668; 		/* rule 816 */
	 tx670 = -(tx151*tx65*tx652); 		/* rule 817 */
	 tx671 = tx151*tx669*tx918; 		/* rule 818 */
	 tx672 = tx670 + tx671; 		/* rule 819 */
	 tx673 = DEImproperRestraintDPhi*tx672; 		/* rule 820 */
	 tx674 = -(tx105*tx149*tx98); 		/* rule 821 */
	 ohy1z2 = tx673 + tx674; 		/* rule 822 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 tx675 = -3*EraseLinearDihedral*K*TWOPI*tx658*tx70*tx75*tx76; 		/* rule 824 */
	 tx676 = -(tx109*tx149*tx98); 		/* rule 825 */
	 ohy1x3 = tx675 + tx676; 		/* rule 826 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 0, ohy1x3);
	 tx677 = -3*EraseLinearDihedral*K*TWOPI*tx663*tx70*tx75*tx76; 		/* rule 828 */
	 tx678 = -(tx113*tx149*tx98); 		/* rule 829 */
	 ohy1y3 = tx677 + tx678; 		/* rule 830 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 1, ohy1y3);
	 tx679 = -3*EraseLinearDihedral*K*TWOPI*tx671*tx70*tx75*tx76; 		/* rule 832 */
	 tx680 = -(tx117*tx149*tx98); 		/* rule 833 */
	 ohy1z3 = tx679 + tx680; 		/* rule 834 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I3, 2, ohy1z3);
	 ohy1x4 = -(tx149*tx20*tx39*tx77*tx90*tx917); 		/* rule 836 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 0, ohy1x4);
	 ohy1y4 = -(tx149*tx45*tx77*tx90*tx917); 		/* rule 838 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 1, ohy1y4);
	 ohy1z4 = -(tx149*tx20*tx40*tx77*tx90*tx917); 		/* rule 840 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I4, 2, ohy1z4);
	 tx681 = tx21*tx436*tx917; 		/* rule 842 */
	 tx682 = tx19*tx243*tx93; 		/* rule 843 */
	 tx683 = tx681 + tx682; 		/* rule 844 */
	 tx684 = tx151*tx683*tx918; 		/* rule 845 */
	 tx685 = tx639 + tx684; 		/* rule 846 */
	 tx686 = DEImproperRestraintDPhi*tx685; 		/* rule 847 */
	 tx687 = -(tx102*tx149*tx97); 		/* rule 848 */
	 ohz1x2 = tx686 + tx687; 		/* rule 849 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 tx688 = tx20*tx436*tx917; 		/* rule 851 */
	 tx689 = tx19*tx341*tx93; 		/* rule 852 */
	 tx690 = tx688 + tx689; 		/* rule 853 */
	 tx691 = tx151*tx690*tx918; 		/* rule 854 */
	 tx692 = tx670 + tx691; 		/* rule 855 */
	 tx693 = -(tx101*tx102*tx149); 		/* rule 856 */
	 tx694 = DEImproperRestraintDPhi*tx692; 		/* rule 857 */
	 ohz1y2 = tx693 + tx694; 		/* rule 858 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx695 = -2.*tx151*tx19*tx203*tx65; 		/* rule 860 */
	 tx696 = tx151*tx482*tx918; 		/* rule 861 */
	 tx697 = tx695 + tx696; 		/* rule 862 */
	 tx698 = -(tx102*tx105*tx149); 		/* rule 863 */
	 tx699 = DEImproperRestraintDPhi*tx697; 		/* rule 864 */
	 ohz1z2 = tx698 + tx699; 		/* rule 865 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 tx700 = -(tx102*tx109*tx149); 		/* rule 867 */
	 tx701 = -3*EraseLinearDihedral*K*TWOPI*tx684*tx70*tx75*tx76; 		/* rule 868 */
	 ohz1x3 = tx700 + tx701; 		/* rule 869 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 0, ohz1x3);
	 tx702 = -(tx102*tx113*tx149); 		/* rule 871 */
	 tx703 = -3*EraseLinearDihedral*K*TWOPI*tx691*tx70*tx75*tx76; 		/* rule 872 */
	 ohz1y3 = tx702 + tx703; 		/* rule 873 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 1, ohz1y3);
	 tx704 = -(tx102*tx117*tx149); 		/* rule 875 */
	 tx705 = -3*EraseLinearDihedral*K*TWOPI*tx696*tx70*tx75*tx76; 		/* rule 876 */
	 ohz1z3 = tx704 + tx705; 		/* rule 877 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I3, 2, ohz1z3);
	 ohz1x4 = -(tx149*tx19*tx39*tx77*tx90*tx917); 		/* rule 879 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 0, ohz1x4);
	 ohz1y4 = -(tx149*tx19*tx38*tx77*tx90*tx917); 		/* rule 881 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 1, ohz1y4);
	 ohz1z4 = -(tx149*tx47*tx77*tx90*tx917); 		/* rule 883 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I4, 2, ohz1z4);
	 tx706 = tx20*tx291; 		/* rule 885 */
	 tx707 = tx292*tx38; 		/* rule 886 */
	 tx708 = tx293*tx38; 		/* rule 887 */
	 tx709 = tx21*tx389; 		/* rule 888 */
	 tx710 = tx39*tx390; 		/* rule 889 */
	 tx711 = tx39*tx391; 		/* rule 890 */
	 tx712 = tx706 + tx709; 		/* rule 891 */
	 tx713 = tx707 + tx710; 		/* rule 892 */
	 tx714 = tx708 + tx711; 		/* rule 893 */
	 tx715 = -tx632; 		/* rule 894 */
	 tx716 = tx151*tx615*tx65; 		/* rule 895 */
	 tx717 = -tx658; 		/* rule 896 */
	 tx718 = 0.5*tx295*tx615*tx77; 		/* rule 897 */
	 tx719 = -0.5*tx295*tx713*tx90; 		/* rule 898 */
	 tx720 = -(tx294*tx714*tx918*tx92); 		/* rule 899 */
	 tx721 = tx151*tx712*tx918*tx93; 		/* rule 900 */
	 tx722 = tx715 + tx716 + tx717 + tx718 + tx719 + tx720 + tx721; 		/* rule 901 */
	 tx723 = DEImproperRestraintDPhi*tx722; 		/* rule 902 */
	 tx724 = tx101*tx149*tx97; 		/* rule 903 */
	 ohx2y2 = tx723 + tx724; 		/* rule 904 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 tx725 = tx19*tx291; 		/* rule 906 */
	 tx726 = tx292*tx40; 		/* rule 907 */
	 tx727 = tx293*tx40; 		/* rule 908 */
	 tx728 = tx21*tx483; 		/* rule 909 */
	 tx729 = tx39*tx484; 		/* rule 910 */
	 tx730 = tx39*tx485; 		/* rule 911 */
	 tx731 = tx725 + tx728; 		/* rule 912 */
	 tx732 = tx726 + tx729; 		/* rule 913 */
	 tx733 = tx727 + tx730; 		/* rule 914 */
	 tx734 = -tx640; 		/* rule 915 */
	 tx735 = tx151*tx620*tx65; 		/* rule 916 */
	 tx736 = -tx684; 		/* rule 917 */
	 tx737 = 0.5*tx295*tx620*tx77; 		/* rule 918 */
	 tx738 = -0.5*tx295*tx732*tx90; 		/* rule 919 */
	 tx739 = -(tx294*tx733*tx918*tx92); 		/* rule 920 */
	 tx740 = tx151*tx731*tx918*tx93; 		/* rule 921 */
	 tx741 = tx734 + tx735 + tx736 + tx737 + tx738 + tx739 + tx740; 		/* rule 922 */
	 tx742 = DEImproperRestraintDPhi*tx741; 		/* rule 923 */
	 tx743 = tx105*tx149*tx97; 		/* rule 924 */
	 ohx2z2 = tx742 + tx743; 		/* rule 925 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 tx744 = -tx297; 		/* rule 927 */
	 tx745 = tx292*tx295*tx39*tx90; 		/* rule 928 */
	 tx746 = tx294*tx526*tx918; 		/* rule 929 */
	 tx747 = 2.*tx293*tx294*tx39*tx918*tx92; 		/* rule 930 */
	 tx748 = -2.*tx151*tx21*tx291*tx918*tx93; 		/* rule 931 */
	 tx749 = tx624 + tx744 + tx745 + tx746 + tx747 + tx748; 		/* rule 932 */
	 tx750 = DEImproperRestraintDPhi*tx749; 		/* rule 933 */
	 tx751 = tx109*tx149*tx97; 		/* rule 934 */
	 ohx2x3 = tx750 + tx751; 		/* rule 935 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 0, ohx2x3);
	 tx752 = tx39*tx559*tx917; 		/* rule 937 */
	 tx753 = tx38*tx523*tx92; 		/* rule 938 */
	 tx754 = tx752 + tx753; 		/* rule 939 */
	 tx755 = -tx721; 		/* rule 940 */
	 tx756 = -0.5*tx295*tx615*tx77; 		/* rule 941 */
	 tx757 = 0.5*tx295*tx713*tx90; 		/* rule 942 */
	 tx758 = tx294*tx754*tx918; 		/* rule 943 */
	 tx759 = tx294*tx714*tx918*tx92; 		/* rule 944 */
	 tx760 = tx632 + tx755 + tx756 + tx757 + tx758 + tx759; 		/* rule 945 */
	 tx761 = DEImproperRestraintDPhi*tx760; 		/* rule 946 */
	 tx762 = tx113*tx149*tx97; 		/* rule 947 */
	 ohx2y3 = tx761 + tx762; 		/* rule 948 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 1, ohx2y3);
	 tx763 = tx39*tx596*tx917; 		/* rule 950 */
	 tx764 = tx40*tx523*tx92; 		/* rule 951 */
	 tx765 = tx763 + tx764; 		/* rule 952 */
	 tx766 = -tx740; 		/* rule 953 */
	 tx767 = -0.5*tx295*tx620*tx77; 		/* rule 954 */
	 tx768 = 0.5*tx295*tx732*tx90; 		/* rule 955 */
	 tx769 = tx294*tx765*tx918; 		/* rule 956 */
	 tx770 = tx294*tx733*tx918*tx92; 		/* rule 957 */
	 tx771 = tx640 + tx766 + tx767 + tx768 + tx769 + tx770; 		/* rule 958 */
	 tx772 = DEImproperRestraintDPhi*tx771; 		/* rule 959 */
	 tx773 = tx117*tx149*tx97; 		/* rule 960 */
	 ohx2z3 = tx772 + tx773; 		/* rule 961 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I3, 2, ohx2z3);
	 tx774 = -3*EraseLinearDihedral*K*TWOPI*tx70*tx746*tx75*tx76; 		/* rule 963 */
	 tx775 = tx149*tx39*tx65*tx90*tx97; 		/* rule 964 */
	 ohx2x4 = tx774 + tx775; 		/* rule 965 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 0, ohx2x4);
	 tx776 = -3*EraseLinearDihedral*K*TWOPI*tx70*tx75*tx758*tx76; 		/* rule 967 */
	 tx777 = tx149*tx38*tx65*tx90*tx97; 		/* rule 968 */
	 ohx2y4 = tx776 + tx777; 		/* rule 969 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 1, ohx2y4);
	 tx778 = -3*EraseLinearDihedral*K*TWOPI*tx70*tx75*tx76*tx769; 		/* rule 971 */
	 tx779 = tx149*tx40*tx65*tx90*tx97; 		/* rule 972 */
	 ohx2z4 = tx778 + tx779; 		/* rule 973 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I4, 2, ohx2z4);
	 tx780 = tx19*tx389; 		/* rule 975 */
	 tx781 = tx390*tx40; 		/* rule 976 */
	 tx782 = tx391*tx40; 		/* rule 977 */
	 tx783 = tx20*tx483; 		/* rule 978 */
	 tx784 = tx38*tx484; 		/* rule 979 */
	 tx785 = tx38*tx485; 		/* rule 980 */
	 tx786 = tx780 + tx783; 		/* rule 981 */
	 tx787 = tx781 + tx784; 		/* rule 982 */
	 tx788 = tx782 + tx785; 		/* rule 983 */
	 tx789 = tx151*tx65*tx652; 		/* rule 984 */
	 tx790 = -tx671; 		/* rule 985 */
	 tx791 = -tx691; 		/* rule 986 */
	 tx792 = 0.5*tx295*tx652*tx77; 		/* rule 987 */
	 tx793 = -0.5*tx295*tx787*tx90; 		/* rule 988 */
	 tx794 = -(tx294*tx788*tx918*tx92); 		/* rule 989 */
	 tx795 = tx151*tx786*tx918*tx93; 		/* rule 990 */
	 tx796 = tx789 + tx790 + tx791 + tx792 + tx793 + tx794 + tx795; 		/* rule 991 */
	 tx797 = tx101*tx105*tx149; 		/* rule 992 */
	 tx798 = DEImproperRestraintDPhi*tx796; 		/* rule 993 */
	 ohy2z2 = tx797 + tx798; 		/* rule 994 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	 tx799 = tx38*tx522*tx917; 		/* rule 996 */
	 tx800 = tx39*tx560*tx92; 		/* rule 997 */
	 tx801 = tx799 + tx800; 		/* rule 998 */
	 tx802 = tx294*tx801*tx918; 		/* rule 999 */
	 tx803 = tx658 + tx755 + tx756 + tx757 + tx759 + tx802; 		/* rule 1000 */
	 tx804 = tx101*tx109*tx149; 		/* rule 1001 */
	 tx805 = DEImproperRestraintDPhi*tx803; 		/* rule 1002 */
	 ohy2x3 = tx804 + tx805; 		/* rule 1003 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 0, ohy2x3);
	 tx806 = -tx393; 		/* rule 1005 */
	 tx807 = tx295*tx38*tx390*tx90; 		/* rule 1006 */
	 tx808 = tx294*tx563*tx918; 		/* rule 1007 */
	 tx809 = 2.*tx294*tx38*tx391*tx918*tx92; 		/* rule 1008 */
	 tx810 = -2.*tx151*tx20*tx389*tx918*tx93; 		/* rule 1009 */
	 tx811 = tx663 + tx806 + tx807 + tx808 + tx809 + tx810; 		/* rule 1010 */
	 tx812 = tx101*tx113*tx149; 		/* rule 1011 */
	 tx813 = DEImproperRestraintDPhi*tx811; 		/* rule 1012 */
	 ohy2y3 = tx812 + tx813; 		/* rule 1013 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 1, ohy2y3);
	 tx814 = tx38*tx596*tx917; 		/* rule 1015 */
	 tx815 = tx40*tx560*tx92; 		/* rule 1016 */
	 tx816 = tx814 + tx815; 		/* rule 1017 */
	 tx817 = -0.5*tx295*tx652*tx77; 		/* rule 1018 */
	 tx818 = -tx795; 		/* rule 1019 */
	 tx819 = 0.5*tx295*tx787*tx90; 		/* rule 1020 */
	 tx820 = tx294*tx816*tx918; 		/* rule 1021 */
	 tx821 = tx294*tx788*tx918*tx92; 		/* rule 1022 */
	 tx822 = tx671 + tx817 + tx818 + tx819 + tx820 + tx821; 		/* rule 1023 */
	 tx823 = tx101*tx117*tx149; 		/* rule 1024 */
	 tx824 = DEImproperRestraintDPhi*tx822; 		/* rule 1025 */
	 ohy2z3 = tx823 + tx824; 		/* rule 1026 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I3, 2, ohy2z3);
	 tx825 = -3*EraseLinearDihedral*K*TWOPI*tx70*tx75*tx76*tx802; 		/* rule 1028 */
	 tx826 = tx101*tx149*tx39*tx65*tx90; 		/* rule 1029 */
	 ohy2x4 = tx825 + tx826; 		/* rule 1030 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 0, ohy2x4);
	 tx827 = -3*EraseLinearDihedral*K*TWOPI*tx70*tx75*tx76*tx808; 		/* rule 1032 */
	 tx828 = tx101*tx149*tx38*tx65*tx90; 		/* rule 1033 */
	 ohy2y4 = tx827 + tx828; 		/* rule 1034 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 1, ohy2y4);
	 tx829 = -3*EraseLinearDihedral*K*TWOPI*tx70*tx75*tx76*tx820; 		/* rule 1036 */
	 tx830 = tx101*tx149*tx40*tx65*tx90; 		/* rule 1037 */
	 ohy2z4 = tx829 + tx830; 		/* rule 1038 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I4, 2, ohy2z4);
	 tx831 = tx40*tx522*tx917; 		/* rule 1040 */
	 tx832 = tx39*tx597*tx92; 		/* rule 1041 */
	 tx833 = tx831 + tx832; 		/* rule 1042 */
	 tx834 = tx294*tx833*tx918; 		/* rule 1043 */
	 tx835 = tx684 + tx766 + tx767 + tx768 + tx770 + tx834; 		/* rule 1044 */
	 tx836 = tx105*tx109*tx149; 		/* rule 1045 */
	 tx837 = DEImproperRestraintDPhi*tx835; 		/* rule 1046 */
	 ohz2x3 = tx836 + tx837; 		/* rule 1047 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 0, ohz2x3);
	 tx838 = tx40*tx559*tx917; 		/* rule 1049 */
	 tx839 = tx38*tx597*tx92; 		/* rule 1050 */
	 tx840 = tx838 + tx839; 		/* rule 1051 */
	 tx841 = tx294*tx840*tx918; 		/* rule 1052 */
	 tx842 = tx691 + tx817 + tx818 + tx819 + tx821 + tx841; 		/* rule 1053 */
	 tx843 = tx105*tx113*tx149; 		/* rule 1054 */
	 tx844 = DEImproperRestraintDPhi*tx842; 		/* rule 1055 */
	 ohz2y3 = tx843 + tx844; 		/* rule 1056 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 1, ohz2y3);
	 tx845 = -tx487; 		/* rule 1058 */
	 tx846 = tx295*tx40*tx484*tx90; 		/* rule 1059 */
	 tx847 = tx294*tx600*tx918; 		/* rule 1060 */
	 tx848 = 2.*tx294*tx40*tx485*tx918*tx92; 		/* rule 1061 */
	 tx849 = -2.*tx151*tx19*tx483*tx918*tx93; 		/* rule 1062 */
	 tx850 = tx696 + tx845 + tx846 + tx847 + tx848 + tx849; 		/* rule 1063 */
	 tx851 = tx105*tx117*tx149; 		/* rule 1064 */
	 tx852 = DEImproperRestraintDPhi*tx850; 		/* rule 1065 */
	 ohz2z3 = tx851 + tx852; 		/* rule 1066 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I3, 2, ohz2z3);
	 tx853 = -3*EraseLinearDihedral*K*TWOPI*tx70*tx75*tx76*tx834; 		/* rule 1068 */
	 tx854 = tx105*tx149*tx39*tx65*tx90; 		/* rule 1069 */
	 ohz2x4 = tx853 + tx854; 		/* rule 1070 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 0, ohz2x4);
	 tx855 = -3*EraseLinearDihedral*K*TWOPI*tx70*tx75*tx76*tx841; 		/* rule 1072 */
	 tx856 = tx105*tx149*tx38*tx65*tx90; 		/* rule 1073 */
	 ohz2y4 = tx855 + tx856; 		/* rule 1074 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 1, ohz2y4);
	 tx857 = -3*EraseLinearDihedral*K*TWOPI*tx70*tx75*tx76*tx847; 		/* rule 1076 */
	 tx858 = tx105*tx149*tx40*tx65*tx90; 		/* rule 1077 */
	 ohz2z4 = tx857 + tx858; 		/* rule 1078 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I4, 2, ohz2z4);
	 tx859 = -(tx294*tx65*tx713); 		/* rule 1080 */
	 tx860 = -tx758; 		/* rule 1081 */
	 tx861 = -tx802; 		/* rule 1082 */
	 tx862 = tx718 + tx719 + tx720 + tx721 + tx859 + tx860 + tx861; 		/* rule 1083 */
	 tx863 = tx109*tx113*tx149; 		/* rule 1084 */
	 tx864 = DEImproperRestraintDPhi*tx862; 		/* rule 1085 */
	 ohx3y3 = tx863 + tx864; 		/* rule 1086 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 1, ohx3y3);
	 tx865 = -(tx294*tx65*tx732); 		/* rule 1088 */
	 tx866 = -tx769; 		/* rule 1089 */
	 tx867 = -tx834; 		/* rule 1090 */
	 tx868 = tx737 + tx738 + tx739 + tx740 + tx865 + tx866 + tx867; 		/* rule 1091 */
	 tx869 = tx109*tx117*tx149; 		/* rule 1092 */
	 tx870 = DEImproperRestraintDPhi*tx868; 		/* rule 1093 */
	 ohx3z3 = tx869 + tx870; 		/* rule 1094 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I3, 2, ohx3z3);
	 tx871 = 2.*tx292*tx294*tx39*tx65; 		/* rule 1096 */
	 tx872 = tx746 + tx871; 		/* rule 1097 */
	 tx873 = DEImproperRestraintDPhi*tx872; 		/* rule 1098 */
	 tx874 = tx109*tx149*tx39*tx65*tx90; 		/* rule 1099 */
	 ohx3x4 = tx873 + tx874; 		/* rule 1100 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 0, ohx3x4);
	 tx875 = tx294*tx65*tx713; 		/* rule 1102 */
	 tx876 = tx758 + tx875; 		/* rule 1103 */
	 tx877 = DEImproperRestraintDPhi*tx876; 		/* rule 1104 */
	 tx878 = tx109*tx149*tx38*tx65*tx90; 		/* rule 1105 */
	 ohx3y4 = tx877 + tx878; 		/* rule 1106 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 1, ohx3y4);
	 tx879 = tx294*tx65*tx732; 		/* rule 1108 */
	 tx880 = tx769 + tx879; 		/* rule 1109 */
	 tx881 = DEImproperRestraintDPhi*tx880; 		/* rule 1110 */
	 tx882 = tx109*tx149*tx40*tx65*tx90; 		/* rule 1111 */
	 ohx3z4 = tx881 + tx882; 		/* rule 1112 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 0, I4, 2, ohx3z4);
	 tx883 = -(tx294*tx65*tx787); 		/* rule 1114 */
	 tx884 = -tx820; 		/* rule 1115 */
	 tx885 = -tx841; 		/* rule 1116 */
	 tx886 = tx792 + tx793 + tx794 + tx795 + tx883 + tx884 + tx885; 		/* rule 1117 */
	 tx887 = tx113*tx117*tx149; 		/* rule 1118 */
	 tx888 = DEImproperRestraintDPhi*tx886; 		/* rule 1119 */
	 ohy3z3 = tx887 + tx888; 		/* rule 1120 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I3, 2, ohy3z3);
	 tx889 = tx802 + tx875; 		/* rule 1122 */
	 tx890 = DEImproperRestraintDPhi*tx889; 		/* rule 1123 */
	 tx891 = tx113*tx149*tx39*tx65*tx90; 		/* rule 1124 */
	 ohy3x4 = tx890 + tx891; 		/* rule 1125 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 0, ohy3x4);
	 tx892 = 2.*tx294*tx38*tx390*tx65; 		/* rule 1127 */
	 tx893 = tx808 + tx892; 		/* rule 1128 */
	 tx894 = DEImproperRestraintDPhi*tx893; 		/* rule 1129 */
	 tx895 = tx113*tx149*tx38*tx65*tx90; 		/* rule 1130 */
	 ohy3y4 = tx894 + tx895; 		/* rule 1131 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 1, ohy3y4);
	 tx896 = tx294*tx65*tx787; 		/* rule 1133 */
	 tx897 = tx820 + tx896; 		/* rule 1134 */
	 tx898 = DEImproperRestraintDPhi*tx897; 		/* rule 1135 */
	 tx899 = tx113*tx149*tx40*tx65*tx90; 		/* rule 1136 */
	 ohy3z4 = tx898 + tx899; 		/* rule 1137 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 1, I4, 2, ohy3z4);
	 tx900 = tx834 + tx879; 		/* rule 1139 */
	 tx901 = tx117*tx149*tx39*tx65*tx90; 		/* rule 1140 */
	 tx902 = DEImproperRestraintDPhi*tx900; 		/* rule 1141 */
	 ohz3x4 = tx901 + tx902; 		/* rule 1142 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 0, ohz3x4);
	 tx903 = tx841 + tx896; 		/* rule 1144 */
	 tx904 = tx117*tx149*tx38*tx65*tx90; 		/* rule 1145 */
	 tx905 = DEImproperRestraintDPhi*tx903; 		/* rule 1146 */
	 ohz3y4 = tx904 + tx905; 		/* rule 1147 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 1, ohz3y4);
	 tx906 = 2.*tx294*tx40*tx484*tx65; 		/* rule 1149 */
	 tx907 = tx847 + tx906; 		/* rule 1150 */
	 tx908 = tx117*tx149*tx40*tx65*tx90; 		/* rule 1151 */
	 tx909 = DEImproperRestraintDPhi*tx907; 		/* rule 1152 */
	 ohz3z4 = tx908 + tx909; 		/* rule 1153 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I3, 2, I4, 2, ohz3z4);
	 tx910 = tx149*tx294*tx38*tx39*tx917; 		/* rule 1155 */
	 tx911 = -3*EraseLinearDihedral*K*TWOPI*tx70*tx75*tx76*tx875; 		/* rule 1156 */
	 ohx4y4 = tx910 + tx911; 		/* rule 1157 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 1, ohx4y4);
	 tx912 = tx149*tx294*tx39*tx40*tx917; 		/* rule 1159 */
	 tx913 = -3*EraseLinearDihedral*K*TWOPI*tx70*tx75*tx76*tx879; 		/* rule 1160 */
	 ohx4z4 = tx912 + tx913; 		/* rule 1161 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 0, I4, 2, ohx4z4);
	 tx914 = tx149*tx294*tx38*tx40*tx917; 		/* rule 1163 */
	 tx915 = -3*EraseLinearDihedral*K*TWOPI*tx70*tx75*tx76*tx896; 		/* rule 1164 */
	 ohy4z4 = tx914 + tx915; 		/* rule 1165 */
	IMPROPER_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I4, 1, I4, 2, ohy4z4);
	} /*calcOffDiagonalHessian*/
	#endif //IMPROPER_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN ]
	} /*calcDiagonalHessian*/
	#endif //IMPROPER_RESTRAINT_CALC_DIAGONAL_HESSIAN ]
	} /*calcForce*/
	#endif //IMPROPER_RESTRAINT_CALC_FORCE ]
	} /*RestraintActive*/
	IMPROPER_RESTRAINT_PHI_SET(Phi);
	IMPROPER_RESTRAINT_DEBUG_INTERACTIONS(I1,I2,I3,I4);
	SKIP_term: (void)0;
/*]; (*RestraintActive*)*/
