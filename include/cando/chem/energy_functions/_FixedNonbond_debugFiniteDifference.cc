/*
    File: _FixedNonbond_debugFiniteDifference.cc
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
TEST_FORCE( _evaluateEnergyOnly_FixedNonbond,  delta2*2.0,
(-delta2 + x1, y1, z1, xf, yf, zf, dA, dC, dQ1Q2),
(delta2 + x1, y1, z1, xf, yf, zf, dA, dC, dQ1Q2),
fx1,
index);
TEST_FORCE( _evaluateEnergyOnly_FixedNonbond,  delta2*2.0,
(x1, -delta2 + y1, z1, xf, yf, zf, dA, dC, dQ1Q2),
(x1, delta2 + y1, z1, xf, yf, zf, dA, dC, dQ1Q2),
fy1,
index);
TEST_FORCE( _evaluateEnergyOnly_FixedNonbond,  delta2*2.0,
(x1, y1, -delta2 + z1, xf, yf, zf, dA, dC, dQ1Q2),
(x1, y1, delta2 + z1, xf, yf, zf, dA, dC, dQ1Q2),
fz1,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_FixedNonbond, delta2*2.0,
(-delta2 + x1, y1, z1, xf, yf, zf, dA, dC, dQ1Q2),
(x1, y1, z1, xf, yf, zf, dA, dC, dQ1Q2),
(delta2 + x1, y1, z1, xf, yf, zf, dA, dC, dQ1Q2),
dhx1x1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_FixedNonbond,
delta2*2.0,
(-delta2 + x1, -delta2 + y1, z1, xf, yf, zf, dA, dC, dQ1Q2),
(delta2 + x1, -delta2 + y1, z1, xf, yf, zf, dA, dC, dQ1Q2),
(-delta2 + x1, delta2 + y1, z1, xf, yf, zf, dA, dC, dQ1Q2),
(delta2 + x1, delta2 + y1, z1, xf, yf, zf, dA, dC, dQ1Q2),
ohx1y1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_FixedNonbond,
delta2*2.0,
(-delta2 + x1, y1, -delta2 + z1, xf, yf, zf, dA, dC, dQ1Q2),
(delta2 + x1, y1, -delta2 + z1, xf, yf, zf, dA, dC, dQ1Q2),
(-delta2 + x1, y1, delta2 + z1, xf, yf, zf, dA, dC, dQ1Q2),
(delta2 + x1, y1, delta2 + z1, xf, yf, zf, dA, dC, dQ1Q2),
ohx1z1,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_FixedNonbond, delta2*2.0,
(x1, -delta2 + y1, z1, xf, yf, zf, dA, dC, dQ1Q2),
(x1, y1, z1, xf, yf, zf, dA, dC, dQ1Q2),
(x1, delta2 + y1, z1, xf, yf, zf, dA, dC, dQ1Q2),
dhy1y1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_FixedNonbond,
delta2*2.0,
(x1, -delta2 + y1, -delta2 + z1, xf, yf, zf, dA, dC, dQ1Q2),
(x1, delta2 + y1, -delta2 + z1, xf, yf, zf, dA, dC, dQ1Q2),
(x1, -delta2 + y1, delta2 + z1, xf, yf, zf, dA, dC, dQ1Q2),
(x1, delta2 + y1, delta2 + z1, xf, yf, zf, dA, dC, dQ1Q2),
ohy1z1,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_FixedNonbond, delta2*2.0,
(x1, y1, -delta2 + z1, xf, yf, zf, dA, dC, dQ1Q2),
(x1, y1, z1, xf, yf, zf, dA, dC, dQ1Q2),
(x1, y1, delta2 + z1, xf, yf, zf, dA, dC, dQ1Q2),
dhz1z1,
index);
