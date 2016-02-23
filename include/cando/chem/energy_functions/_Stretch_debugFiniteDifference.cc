/*
    File: _Stretch_debugFiniteDifference.cc
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
// FINITE_DIFFERENCE_TEST 
TEST_FORCE( _evaluateEnergyOnly_Stretch,  delta2*2.0,
(-delta2 + x1, y1, z1, x2, y2, z2, r0, kb),
(delta2 + x1, y1, z1, x2, y2, z2, r0, kb),
fx1,
index);
TEST_FORCE( _evaluateEnergyOnly_Stretch,  delta2*2.0,
(x1, -delta2 + y1, z1, x2, y2, z2, r0, kb),
(x1, delta2 + y1, z1, x2, y2, z2, r0, kb),
fy1,
index);
TEST_FORCE( _evaluateEnergyOnly_Stretch,  delta2*2.0,
(x1, y1, -delta2 + z1, x2, y2, z2, r0, kb),
(x1, y1, delta2 + z1, x2, y2, z2, r0, kb),
fz1,
index);
TEST_FORCE( _evaluateEnergyOnly_Stretch,  delta2*2.0,
(x1, y1, z1, -delta2 + x2, y2, z2, r0, kb),
(x1, y1, z1, delta2 + x2, y2, z2, r0, kb),
fx2,
index);
TEST_FORCE( _evaluateEnergyOnly_Stretch,  delta2*2.0,
(x1, y1, z1, x2, -delta2 + y2, z2, r0, kb),
(x1, y1, z1, x2, delta2 + y2, z2, r0, kb),
fy2,
index);
TEST_FORCE( _evaluateEnergyOnly_Stretch,  delta2*2.0,
(x1, y1, z1, x2, y2, -delta2 + z2, r0, kb),
(x1, y1, z1, x2, y2, delta2 + z2, r0, kb),
fz2,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch, delta2*2.0,
(-delta2 + x1, y1, z1, x2, y2, z2, r0, kb),
(x1, y1, z1, x2, y2, z2, r0, kb),
(delta2 + x1, y1, z1, x2, y2, z2, r0, kb),
dhx1x1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch,
delta2*2.0,
(-delta2 + x1, -delta2 + y1, z1, x2, y2, z2, r0, kb),
(delta2 + x1, -delta2 + y1, z1, x2, y2, z2, r0, kb),
(-delta2 + x1, delta2 + y1, z1, x2, y2, z2, r0, kb),
(delta2 + x1, delta2 + y1, z1, x2, y2, z2, r0, kb),
ohx1y1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch,
delta2*2.0,
(-delta2 + x1, y1, -delta2 + z1, x2, y2, z2, r0, kb),
(delta2 + x1, y1, -delta2 + z1, x2, y2, z2, r0, kb),
(-delta2 + x1, y1, delta2 + z1, x2, y2, z2, r0, kb),
(delta2 + x1, y1, delta2 + z1, x2, y2, z2, r0, kb),
ohx1z1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch,
delta2*2.0,
(-delta2 + x1, y1, z1, -delta2 + x2, y2, z2, r0, kb),
(delta2 + x1, y1, z1, -delta2 + x2, y2, z2, r0, kb),
(-delta2 + x1, y1, z1, delta2 + x2, y2, z2, r0, kb),
(delta2 + x1, y1, z1, delta2 + x2, y2, z2, r0, kb),
ohx1x2,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch,
delta2*2.0,
(-delta2 + x1, y1, z1, x2, -delta2 + y2, z2, r0, kb),
(delta2 + x1, y1, z1, x2, -delta2 + y2, z2, r0, kb),
(-delta2 + x1, y1, z1, x2, delta2 + y2, z2, r0, kb),
(delta2 + x1, y1, z1, x2, delta2 + y2, z2, r0, kb),
ohx1y2,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch,
delta2*2.0,
(-delta2 + x1, y1, z1, x2, y2, -delta2 + z2, r0, kb),
(delta2 + x1, y1, z1, x2, y2, -delta2 + z2, r0, kb),
(-delta2 + x1, y1, z1, x2, y2, delta2 + z2, r0, kb),
(delta2 + x1, y1, z1, x2, y2, delta2 + z2, r0, kb),
ohx1z2,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch, delta2*2.0,
(x1, -delta2 + y1, z1, x2, y2, z2, r0, kb),
(x1, y1, z1, x2, y2, z2, r0, kb),
(x1, delta2 + y1, z1, x2, y2, z2, r0, kb),
dhy1y1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch,
delta2*2.0,
(x1, -delta2 + y1, -delta2 + z1, x2, y2, z2, r0, kb),
(x1, delta2 + y1, -delta2 + z1, x2, y2, z2, r0, kb),
(x1, -delta2 + y1, delta2 + z1, x2, y2, z2, r0, kb),
(x1, delta2 + y1, delta2 + z1, x2, y2, z2, r0, kb),
ohy1z1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch,
delta2*2.0,
(x1, -delta2 + y1, z1, -delta2 + x2, y2, z2, r0, kb),
(x1, delta2 + y1, z1, -delta2 + x2, y2, z2, r0, kb),
(x1, -delta2 + y1, z1, delta2 + x2, y2, z2, r0, kb),
(x1, delta2 + y1, z1, delta2 + x2, y2, z2, r0, kb),
ohy1x2,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch,
delta2*2.0,
(x1, -delta2 + y1, z1, x2, -delta2 + y2, z2, r0, kb),
(x1, delta2 + y1, z1, x2, -delta2 + y2, z2, r0, kb),
(x1, -delta2 + y1, z1, x2, delta2 + y2, z2, r0, kb),
(x1, delta2 + y1, z1, x2, delta2 + y2, z2, r0, kb),
ohy1y2,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch,
delta2*2.0,
(x1, -delta2 + y1, z1, x2, y2, -delta2 + z2, r0, kb),
(x1, delta2 + y1, z1, x2, y2, -delta2 + z2, r0, kb),
(x1, -delta2 + y1, z1, x2, y2, delta2 + z2, r0, kb),
(x1, delta2 + y1, z1, x2, y2, delta2 + z2, r0, kb),
ohy1z2,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch, delta2*2.0,
(x1, y1, -delta2 + z1, x2, y2, z2, r0, kb),
(x1, y1, z1, x2, y2, z2, r0, kb),
(x1, y1, delta2 + z1, x2, y2, z2, r0, kb),
dhz1z1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch,
delta2*2.0,
(x1, y1, -delta2 + z1, -delta2 + x2, y2, z2, r0, kb),
(x1, y1, delta2 + z1, -delta2 + x2, y2, z2, r0, kb),
(x1, y1, -delta2 + z1, delta2 + x2, y2, z2, r0, kb),
(x1, y1, delta2 + z1, delta2 + x2, y2, z2, r0, kb),
ohz1x2,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch,
delta2*2.0,
(x1, y1, -delta2 + z1, x2, -delta2 + y2, z2, r0, kb),
(x1, y1, delta2 + z1, x2, -delta2 + y2, z2, r0, kb),
(x1, y1, -delta2 + z1, x2, delta2 + y2, z2, r0, kb),
(x1, y1, delta2 + z1, x2, delta2 + y2, z2, r0, kb),
ohz1y2,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch,
delta2*2.0,
(x1, y1, -delta2 + z1, x2, y2, -delta2 + z2, r0, kb),
(x1, y1, delta2 + z1, x2, y2, -delta2 + z2, r0, kb),
(x1, y1, -delta2 + z1, x2, y2, delta2 + z2, r0, kb),
(x1, y1, delta2 + z1, x2, y2, delta2 + z2, r0, kb),
ohz1z2,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch, delta2*2.0,
(x1, y1, z1, -delta2 + x2, y2, z2, r0, kb),
(x1, y1, z1, x2, y2, z2, r0, kb),
(x1, y1, z1, delta2 + x2, y2, z2, r0, kb),
dhx2x2,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch,
delta2*2.0,
(x1, y1, z1, -delta2 + x2, -delta2 + y2, z2, r0, kb),
(x1, y1, z1, delta2 + x2, -delta2 + y2, z2, r0, kb),
(x1, y1, z1, -delta2 + x2, delta2 + y2, z2, r0, kb),
(x1, y1, z1, delta2 + x2, delta2 + y2, z2, r0, kb),
ohx2y2,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch,
delta2*2.0,
(x1, y1, z1, -delta2 + x2, y2, -delta2 + z2, r0, kb),
(x1, y1, z1, delta2 + x2, y2, -delta2 + z2, r0, kb),
(x1, y1, z1, -delta2 + x2, y2, delta2 + z2, r0, kb),
(x1, y1, z1, delta2 + x2, y2, delta2 + z2, r0, kb),
ohx2z2,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch, delta2*2.0,
(x1, y1, z1, x2, -delta2 + y2, z2, r0, kb),
(x1, y1, z1, x2, y2, z2, r0, kb),
(x1, y1, z1, x2, delta2 + y2, z2, r0, kb),
dhy2y2,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch,
delta2*2.0,
(x1, y1, z1, x2, -delta2 + y2, -delta2 + z2, r0, kb),
(x1, y1, z1, x2, delta2 + y2, -delta2 + z2, r0, kb),
(x1, y1, z1, x2, -delta2 + y2, delta2 + z2, r0, kb),
(x1, y1, z1, x2, delta2 + y2, delta2 + z2, r0, kb),
ohy2z2,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Stretch, delta2*2.0,
(x1, y1, z1, x2, y2, -delta2 + z2, r0, kb),
(x1, y1, z1, x2, y2, z2, r0, kb),
(x1, y1, z1, x2, y2, delta2 + z2, r0, kb),
dhz2z2,
index);
