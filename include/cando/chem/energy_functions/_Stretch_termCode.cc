/*
    File: _Stretch_termCode.cc
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
	STRETCH_SET_PARAMETER(kb);
	STRETCH_SET_PARAMETER(r0);
	STRETCH_SET_PARAMETER(I1);
	STRETCH_SET_PARAMETER(I2);
	STRETCH_SET_POSITION(x1,I1,0);
	STRETCH_SET_POSITION(y1,I1,1);
	STRETCH_SET_POSITION(z1,I1,2);
	STRETCH_SET_POSITION(x2,I2,0);
	STRETCH_SET_POSITION(y2,I2,1);
	STRETCH_SET_POSITION(z2,I2,2);
	 tx44 = -x1; 		/* rule 11 */
	 tx45 = -y1; 		/* rule 12 */
	 tx46 = -z1; 		/* rule 13 */
	 tx47 = tx44 + x2; 		/* rule 14 */
	 tx48 = tx45 + y2; 		/* rule 15 */
	 tx49 = tx46 + z2; 		/* rule 16 */
	 tx50 = power2(tx47); 		/* rule 17 */
	 tx51 = power2(tx48); 		/* rule 18 */
	 tx52 = power2(tx49); 		/* rule 19 */
	 tx53 = tx50 + tx51 + tx52; 		/* rule 20 */
	 tx54 = -r0; 		/* rule 21 */
	 tx55 = mysqrt(tx53); 		/* rule 22 */
	 StretchDeviation = tx54 + tx55; 		/* rule 23 */
	 tx56 = power2(StretchDeviation); 		/* rule 24 */
	 Energy = kb*tx56; 		/* rule 25 */
	STRETCH_ENERGY_ACCUMULATE(Energy);
	#ifdef STRETCH_CALC_FORCE //[
	if ( calcForce ) {
	 tx57 = StretchDeviation; 		/* rule 29 */
	 tx58 = reciprocal(tx55); 		/* rule 30 */
	 gx1 = -2.*kb*tx47*tx57*tx58; 		/* rule 31 */
	 fx1 = -gx1; 		/* rule 32 */
	STRETCH_FORCE_ACCUMULATE(I1, 0, fx1 );
	 gy1 = -2.*kb*tx48*tx57*tx58; 		/* rule 34 */
	 fy1 = -gy1; 		/* rule 35 */
	STRETCH_FORCE_ACCUMULATE(I1, 1, fy1 );
	 gz1 = -2.*kb*tx49*tx57*tx58; 		/* rule 37 */
	 fz1 = -gz1; 		/* rule 38 */
	STRETCH_FORCE_ACCUMULATE(I1, 2, fz1 );
	 gx2 = 2.*kb*tx47*tx57*tx58; 		/* rule 40 */
	 fx2 = -gx2; 		/* rule 41 */
	STRETCH_FORCE_ACCUMULATE(I2, 0, fx2 );
	 gy2 = 2.*kb*tx48*tx57*tx58; 		/* rule 43 */
	 fy2 = -gy2; 		/* rule 44 */
	STRETCH_FORCE_ACCUMULATE(I2, 1, fy2 );
	 gz2 = 2.*kb*tx49*tx57*tx58; 		/* rule 46 */
	 fz2 = -gz2; 		/* rule 47 */
	STRETCH_FORCE_ACCUMULATE(I2, 2, fz2 );
	#ifdef STRETCH_CALC_DIAGONAL_HESSIAN //[
	if ( calcDiagonalHessian ) {
	 tx87 = tx53; 		/* rule 51 */
	 tx88 = tx58; 		/* rule 52 */
	 tx89 = reciprocal(tx87); 		/* rule 53 */
	 tx59 = tx88*tx89; 		/* rule 54 */
	 tx60 = tx89; 		/* rule 55 */
	 tx61 = 2.*kb*tx57*tx88; 		/* rule 56 */
	 tx62 = -2.*kb*tx50*tx57*tx59; 		/* rule 57 */
	 tx63 = 2.*kb*tx50*tx60; 		/* rule 58 */
	 dhx1x1 = tx61 + tx62 + tx63; 		/* rule 59 */
	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 0, dhx1x1);
	 tx64 = -2.*kb*tx51*tx57*tx59; 		/* rule 61 */
	 tx65 = 2.*kb*tx51*tx60; 		/* rule 62 */
	 dhy1y1 = tx61 + tx64 + tx65; 		/* rule 63 */
	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 1, dhy1y1);
	 tx66 = -2.*kb*tx52*tx57*tx59; 		/* rule 65 */
	 tx67 = 2.*kb*tx52*tx60; 		/* rule 66 */
	 dhz1z1 = tx61 + tx66 + tx67; 		/* rule 67 */
	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I1, 2, dhz1z1);
	 dhx2x2 = dhx1x1; 		/* rule 69 */
	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 0, dhx2x2);
	 dhy2y2 = dhy1y1; 		/* rule 71 */
	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 1, dhy2y2);
	 dhz2z2 = dhz1z1; 		/* rule 73 */
	STRETCH_DIAGONAL_HESSIAN_ACCUMULATE(I2, 2, I2, 2, dhz2z2);
	#ifdef STRETCH_CALC_OFF_DIAGONAL_HESSIAN //[
	if ( calcOffDiagonalHessian ) {
	 tx68 = -2.*kb*tx47*tx48*tx57*tx59; 		/* rule 77 */
	 tx69 = 2.*kb*tx47*tx48*tx60; 		/* rule 78 */
	 ohx1y1 = tx68 + tx69; 		/* rule 79 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 1, ohx1y1);
	 tx70 = -2.*kb*tx47*tx49*tx57*tx59; 		/* rule 81 */
	 tx71 = 2.*kb*tx47*tx49*tx60; 		/* rule 82 */
	 ohx1z1 = tx70 + tx71; 		/* rule 83 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I1, 2, ohx1z1);
	 tx72 = -2.*kb*tx57*tx88; 		/* rule 85 */
	 tx73 = 2.*kb*tx50*tx57*tx59; 		/* rule 86 */
	 tx74 = -2.*kb*tx50*tx60; 		/* rule 87 */
	 ohx1x2 = tx72 + tx73 + tx74; 		/* rule 88 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 0, ohx1x2);
	 tx75 = 2.*kb*tx47*tx48*tx57*tx59; 		/* rule 90 */
	 tx76 = -2.*kb*tx47*tx48*tx60; 		/* rule 91 */
	 ohx1y2 = tx75 + tx76; 		/* rule 92 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 1, ohx1y2);
	 tx77 = 2.*kb*tx47*tx49*tx57*tx59; 		/* rule 94 */
	 tx78 = -2.*kb*tx47*tx49*tx60; 		/* rule 95 */
	 ohx1z2 = tx77 + tx78; 		/* rule 96 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 0, I2, 2, ohx1z2);
	 tx79 = -2.*kb*tx48*tx49*tx57*tx59; 		/* rule 98 */
	 tx80 = 2.*kb*tx48*tx49*tx60; 		/* rule 99 */
	 ohy1z1 = tx79 + tx80; 		/* rule 100 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I1, 2, ohy1z1);
	 ohy1x2 = ohx1y2; 		/* rule 102 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 0, ohy1x2);
	 tx81 = 2.*kb*tx51*tx57*tx59; 		/* rule 104 */
	 tx82 = -2.*kb*tx51*tx60; 		/* rule 105 */
	 ohy1y2 = tx72 + tx81 + tx82; 		/* rule 106 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 1, ohy1y2);
	 tx83 = 2.*kb*tx48*tx49*tx57*tx59; 		/* rule 108 */
	 tx84 = -2.*kb*tx48*tx49*tx60; 		/* rule 109 */
	 ohy1z2 = tx83 + tx84; 		/* rule 110 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 1, I2, 2, ohy1z2);
	 ohz1x2 = ohx1z2; 		/* rule 112 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 0, ohz1x2);
	 ohz1y2 = ohy1z2; 		/* rule 114 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 1, ohz1y2);
	 tx85 = 2.*kb*tx52*tx57*tx59; 		/* rule 116 */
	 tx86 = -2.*kb*tx52*tx60; 		/* rule 117 */
	 ohz1z2 = tx72 + tx85 + tx86; 		/* rule 118 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I1, 2, I2, 2, ohz1z2);
	 ohx2y2 = ohx1y1; 		/* rule 120 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 1, ohx2y2);
	 ohx2z2 = ohx1z1; 		/* rule 122 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 0, I2, 2, ohx2z2);
	 ohy2z2 = ohy1z1; 		/* rule 124 */
	STRETCH_OFF_DIAGONAL_HESSIAN_ACCUMULATE(I2, 1, I2, 2, ohy2z2);
	} /*if calcOffDiagonalHessian */ 
	#endif /* STRETCH_CALC_OFF_DIAGONAL_HESSIAN ]*/
	} /*calcDiagonalHessian */
	#endif /* STRETCH_CALC_DIAGONAL_HESSIAN ]*/
	} /*calcForce */
	#endif /* STRETCH_CALC_FORCE ]*/
