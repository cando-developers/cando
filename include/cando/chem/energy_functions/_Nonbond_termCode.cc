/*
    File: _Nonbond_termCode.cc
*/
/*
Open Source License
Copyright (c) 2016, Christian E. Schafmeister
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 
This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
// C-code
	NONBOND_SET_PARAMETER(I1);
	NONBOND_SET_PARAMETER(I2);
	NONBOND_SET_PARAMETER(dQ1Q2);
	NONBOND_SET_PARAMETER(dA);
	NONBOND_SET_PARAMETER(dC);
	NONBOND_SET_POSITION(x1,I1,0);
	NONBOND_SET_POSITION(y1,I1,1);
	NONBOND_SET_POSITION(z1,I1,2);
	NONBOND_SET_POSITION(x2,I2,0);
	NONBOND_SET_POSITION(y2,I2,1);
	NONBOND_SET_POSITION(z2,I2,2);
	 tx82 = -x2; 		/* rule 12 */
	 tx83 = -y2; 		/* rule 13 */
	 tx84 = -z2; 		/* rule 14 */
	 tx85 = tx82 + x1; 		/* rule 15 */
	 tx86 = tx83 + y1; 		/* rule 16 */
	 tx87 = tx84 + z1; 		/* rule 17 */
	 tx88 = power2(tx85); 		/* rule 18 */
	 tx89 = power2(tx86); 		/* rule 19 */
	 tx90 = power2(tx87); 		/* rule 20 */
	 tx91 = tx88 + tx89 + tx90; 		/* rule 21 */
	 NonbondDistance = mysqrt(tx91); 		/* rule 22 */
	 tx194 = reciprocal(tx91); 		/* rule 23 */
	 tx172 = power2(tx194); 		/* rule 24 */
	 tx173 = tx194; 		/* rule 25 */
	 tx163 = tx172*tx173; 		/* rule 26 */
	 tx92 = power2(tx163); 		/* rule 27 */
	 tx164 = power2(tx173); 		/* rule 28 */
	 tx165 = tx173; 		/* rule 29 */
	 tx93 = tx164*tx165; 		/* rule 30 */
	 tx94 = dA*tx92; 		/* rule 31 */
	 tx95 = -(dC*tx93); 		/* rule 32 */
	 Evdw = tx94 + tx95; 		/* rule 33 */
	NONBOND_EVDW_ENERGY_ACCUMULATE(Evdw);
	 tx96 = reciprocal(NonbondDistance); 		/* rule 35 */
	 Eeel = dQ1Q2*tx96; 		/* rule 36 */
	NONBOND_EEEL_ENERGY_ACCUMULATE(Eeel);
	 Energy = Eeel + Evdw; 		/* rule 38 */
	NONBOND_ENERGY_ACCUMULATE(Energy);
	#ifdef NONBOND_CALC_FORCE //[
	if ( calcForce ) {
	 tx166 = power2(tx172); 		/* rule 42 */
	 tx97 = tx163*tx166; 		/* rule 43 */
	 tx98 = power2(tx164); 		/* rule 44 */
	 tx167 = tx91; 		/* rule 45 */
	 tx168 = tx96; 		/* rule 46 */
	 tx169 = reciprocal(tx167); 		/* rule 47 */
	 tx99 = tx168*tx169; 		/* rule 48 */
	 tzz178 = dA*tx97; 		/* rule 49 */
	 tzz193 = -12.*tzz178; 		/* rule 50 */
	 tx100 = tx85*tzz193; 		/* rule 51 */
	 tzz177 = dC*tx98; 		/* rule 52 */
	 tzz190 = 6.*tzz177; 		/* rule 53 */
	 tx101 = tx85*tzz190; 		/* rule 54 */
	 tzz179 = dQ1Q2*tx99; 		/* rule 55 */
	 tzz191 = -tzz179; 		/* rule 56 */
	 tx102 = tx85*tzz191; 		/* rule 57 */
	 gx1 = tx100 + tx101 + tx102; 		/* rule 58 */
	 fx1 = -gx1; 		/* rule 59 */
	NONBOND_FORCE_ACCUMULATE(I1, 0, fx1 );
	 tx103 = tx86*tzz193; 		/* rule 61 */
	 tx104 = tx86*tzz190; 		/* rule 62 */
	 tx105 = tx86*tzz191; 		/* rule 63 */
	 gy1 = tx103 + tx104 + tx105; 		/* rule 64 */
	 fy1 = -gy1; 		/* rule 65 */
	NONBOND_FORCE_ACCUMULATE(I1, 1, fy1 );
	 tx106 = tx87*tzz193; 		/* rule 67 */
	 tx107 = tx87*tzz190; 		/* rule 68 */
	 tx108 = tx87*tzz191; 		/* rule 69 */
	 gz1 = tx106 + tx107 + tx108; 		/* rule 70 */
	 fz1 = -gz1; 		/* rule 71 */
	NONBOND_FORCE_ACCUMULATE(I1, 2, fz1 );
	 tzz189 = 12.*tzz178; 		/* rule 73 */
	 tx109 = tx85*tzz189; 		/* rule 74 */
	 tzz192 = -6.*tzz177; 		/* rule 75 */
	 tx110 = tx85*tzz192; 		/* rule 76 */
	 tx111 = tx85*tzz179; 		/* rule 77 */
	 gx2 = tx109 + tx110 + tx111; 		/* rule 78 */
	 fx2 = -gx2; 		/* rule 79 */
	NONBOND_FORCE_ACCUMULATE(I2, 0, fx2 );
	 tx112 = tx86*tzz189; 		/* rule 81 */
	 tx113 = tx86*tzz192; 		/* rule 82 */
	 tx114 = tx86*tzz179; 		/* rule 83 */
	 gy2 = tx112 + tx113 + tx114; 		/* rule 84 */
	 fy2 = -gy2; 		/* rule 85 */
	NONBOND_FORCE_ACCUMULATE(I2, 1, fy2 );
	 tx115 = tx87*tzz189; 		/* rule 87 */
	 tx116 = tx87*tzz192; 		/* rule 88 */
	 tx117 = tx87*tzz179; 		/* rule 89 */
	 gz2 = tx115 + tx116 + tx117; 		/* rule 90 */
	 fz2 = -gz2; 		/* rule 91 */
	NONBOND_FORCE_ACCUMULATE(I2, 2, fz2 );
	#ifdef NONBOND_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx118 = power2(tx166); 		/* rule 95 */
	 tx119 = tx163*tx164; 		/* rule 96 */
	 tx170 = power2(tx167); 		/* rule 97 */
	 tx171 = reciprocal(tx170); 		/* rule 98 */
	 tx120 = tx168*tx171; 		/* rule 99 */
	 tzz176 = dA*tx118; 		/* rule 100 */
	 tzz183 = 168*tzz176; 		/* rule 101 */
	 tx121 = tx88*tzz183; 		/* rule 102 */
	 tzz175 = dC*tx119; 		/* rule 103 */
	 tzz187 = -48*tzz175; 		/* rule 104 */
	 tx122 = tx88*tzz187; 		/* rule 105 */
	 tzz174 = dQ1Q2*tx120; 		/* rule 106 */
	 tzz185 = 3*tzz174; 		/* rule 107 */
	 tx123 = tx88*tzz185; 		/* rule 108 */
	 tx124 = tzz193; 		/* rule 109 */
	 tx125 = tzz190; 		/* rule 110 */
	 tx126 = tzz191; 		/* rule 111 */
	 tzz197 = tx125 + tx126; 		/* rule 112 */
	 tzz198 = tx124 + tzz197; 		/* rule 113 */
	 dhx1x1 = tx121 + tx122 + tx123 + tzz198; 		/* rule 114 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx127 = tx89*tzz183; 		/* rule 116 */
	 tx128 = tx89*tzz187; 		/* rule 117 */
	 tx129 = tx89*tzz185; 		/* rule 118 */
	 dhy1y1 = tx127 + tx128 + tx129 + tzz198; 		/* rule 119 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx130 = tx90*tzz183; 		/* rule 121 */
	 tx131 = tx90*tzz187; 		/* rule 122 */
	 tx132 = tx90*tzz185; 		/* rule 123 */
	 dhz1z1 = tx130 + tx131 + tx132 + tzz198; 		/* rule 124 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 dhx2x2 = dhx1x1; 		/* rule 126 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 dhy2y2 = dhy1y1; 		/* rule 128 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 dhz2z2 = dhz1z1; 		/* rule 130 */
	NONBOND_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	#ifdef NONBOND_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tzz182 = tx85*tx86; 		/* rule 134 */
	 tx133 = tzz182*tzz183; 		/* rule 135 */
	 tx134 = tzz182*tzz187; 		/* rule 136 */
	 tx135 = tzz182*tzz185; 		/* rule 137 */
	 ohx1y1 = tx133 + tx134 + tx135; 		/* rule 138 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tzz181 = tx85*tx87; 		/* rule 140 */
	 tx136 = tzz181*tzz183; 		/* rule 141 */
	 tx137 = tzz181*tzz187; 		/* rule 142 */
	 tx138 = tzz181*tzz185; 		/* rule 143 */
	 ohx1z1 = tx136 + tx137 + tx138; 		/* rule 144 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tzz188 = -168*tzz176; 		/* rule 146 */
	 tx139 = tx88*tzz188; 		/* rule 147 */
	 tzz184 = 48*tzz175; 		/* rule 148 */
	 tx140 = tx88*tzz184; 		/* rule 149 */
	 tzz186 = -3*tzz174; 		/* rule 150 */
	 tx141 = tx88*tzz186; 		/* rule 151 */
	 tx142 = tzz189; 		/* rule 152 */
	 tx143 = tzz192; 		/* rule 153 */
	 tx144 = tzz179; 		/* rule 154 */
	 tzz195 = tx143 + tx144; 		/* rule 155 */
	 tzz196 = tx142 + tzz195; 		/* rule 156 */
	 ohx1x2 = tx139 + tx140 + tx141 + tzz196; 		/* rule 157 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx145 = tzz182*tzz188; 		/* rule 159 */
	 tx146 = tzz182*tzz184; 		/* rule 160 */
	 tx147 = tzz182*tzz186; 		/* rule 161 */
	 ohx1y2 = tx145 + tx146 + tx147; 		/* rule 162 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx148 = tzz181*tzz188; 		/* rule 164 */
	 tx149 = tzz181*tzz184; 		/* rule 165 */
	 tx150 = tzz181*tzz186; 		/* rule 166 */
	 ohx1z2 = tx148 + tx149 + tx150; 		/* rule 167 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tzz180 = tx86*tx87; 		/* rule 169 */
	 tx151 = tzz180*tzz183; 		/* rule 170 */
	 tx152 = tzz180*tzz187; 		/* rule 171 */
	 tx153 = tzz180*tzz185; 		/* rule 172 */
	 ohy1z1 = tx151 + tx152 + tx153; 		/* rule 173 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 ohy1x2 = ohx1y2; 		/* rule 175 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx154 = tx89*tzz188; 		/* rule 177 */
	 tx155 = tx89*tzz184; 		/* rule 178 */
	 tx156 = tx89*tzz186; 		/* rule 179 */
	 ohy1y2 = tx154 + tx155 + tx156 + tzz196; 		/* rule 180 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx157 = tzz180*tzz188; 		/* rule 182 */
	 tx158 = tzz180*tzz184; 		/* rule 183 */
	 tx159 = tzz180*tzz186; 		/* rule 184 */
	 ohy1z2 = tx157 + tx158 + tx159; 		/* rule 185 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 ohz1x2 = ohx1z2; 		/* rule 187 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 ohz1y2 = ohy1z2; 		/* rule 189 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx160 = tx90*tzz188; 		/* rule 191 */
	 tx161 = tx90*tzz184; 		/* rule 192 */
	 tx162 = tx90*tzz186; 		/* rule 193 */
	 ohz1z2 = tx160 + tx161 + tx162 + tzz196; 		/* rule 194 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 ohx2y2 = ohx1y1; 		/* rule 196 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 ohx2z2 = ohx1z1; 		/* rule 198 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 ohy2z2 = ohy1z1; 		/* rule 200 */
	NONBOND_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	} /*if calcOffDiagonalHessian */ 
	#endif /* NONBOND_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* NONBOND_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* NONBOND_CALC_FORCE ]*/
