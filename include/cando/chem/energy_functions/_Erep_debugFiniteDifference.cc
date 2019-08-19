// FINITE_DIFFERENCE_TEST 
TEST_FORCE( _evaluateEnergyOnly_Erep,  delta2*2.0,
(-delta2 + x1, y1, z1, x2, y2, z2, crep),
(delta2 + x1, y1, z1, x2, y2, z2, crep),
fx1,
index);
TEST_FORCE( _evaluateEnergyOnly_Erep,  delta2*2.0,
(x1, -delta2 + y1, z1, x2, y2, z2, crep),
(x1, delta2 + y1, z1, x2, y2, z2, crep),
fy1,
index);
TEST_FORCE( _evaluateEnergyOnly_Erep,  delta2*2.0,
(x1, y1, -delta2 + z1, x2, y2, z2, crep),
(x1, y1, delta2 + z1, x2, y2, z2, crep),
fz1,
index);
TEST_FORCE( _evaluateEnergyOnly_Erep,  delta2*2.0,
(x1, y1, z1, -delta2 + x2, y2, z2, crep),
(x1, y1, z1, delta2 + x2, y2, z2, crep),
fx2,
index);
TEST_FORCE( _evaluateEnergyOnly_Erep,  delta2*2.0,
(x1, y1, z1, x2, -delta2 + y2, z2, crep),
(x1, y1, z1, x2, delta2 + y2, z2, crep),
fy2,
index);
TEST_FORCE( _evaluateEnergyOnly_Erep,  delta2*2.0,
(x1, y1, z1, x2, y2, -delta2 + z2, crep),
(x1, y1, z1, x2, y2, delta2 + z2, crep),
fz2,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep, delta2*2.0,
(-delta2 + x1, y1, z1, x2, y2, z2, crep),
(x1, y1, z1, x2, y2, z2, crep),
(delta2 + x1, y1, z1, x2, y2, z2, crep),
dhx1x1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep,
delta2*2.0,
(-delta2 + x1, -delta2 + y1, z1, x2, y2, z2, crep),
(delta2 + x1, -delta2 + y1, z1, x2, y2, z2, crep),
(-delta2 + x1, delta2 + y1, z1, x2, y2, z2, crep),
(delta2 + x1, delta2 + y1, z1, x2, y2, z2, crep),
ohx1y1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep,
delta2*2.0,
(-delta2 + x1, y1, -delta2 + z1, x2, y2, z2, crep),
(delta2 + x1, y1, -delta2 + z1, x2, y2, z2, crep),
(-delta2 + x1, y1, delta2 + z1, x2, y2, z2, crep),
(delta2 + x1, y1, delta2 + z1, x2, y2, z2, crep),
ohx1z1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep,
delta2*2.0,
(-delta2 + x1, y1, z1, -delta2 + x2, y2, z2, crep),
(delta2 + x1, y1, z1, -delta2 + x2, y2, z2, crep),
(-delta2 + x1, y1, z1, delta2 + x2, y2, z2, crep),
(delta2 + x1, y1, z1, delta2 + x2, y2, z2, crep),
ohx1x2,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep,
delta2*2.0,
(-delta2 + x1, y1, z1, x2, -delta2 + y2, z2, crep),
(delta2 + x1, y1, z1, x2, -delta2 + y2, z2, crep),
(-delta2 + x1, y1, z1, x2, delta2 + y2, z2, crep),
(delta2 + x1, y1, z1, x2, delta2 + y2, z2, crep),
ohx1y2,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep,
delta2*2.0,
(-delta2 + x1, y1, z1, x2, y2, -delta2 + z2, crep),
(delta2 + x1, y1, z1, x2, y2, -delta2 + z2, crep),
(-delta2 + x1, y1, z1, x2, y2, delta2 + z2, crep),
(delta2 + x1, y1, z1, x2, y2, delta2 + z2, crep),
ohx1z2,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep, delta2*2.0,
(x1, -delta2 + y1, z1, x2, y2, z2, crep),
(x1, y1, z1, x2, y2, z2, crep),
(x1, delta2 + y1, z1, x2, y2, z2, crep),
dhy1y1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep,
delta2*2.0,
(x1, -delta2 + y1, -delta2 + z1, x2, y2, z2, crep),
(x1, delta2 + y1, -delta2 + z1, x2, y2, z2, crep),
(x1, -delta2 + y1, delta2 + z1, x2, y2, z2, crep),
(x1, delta2 + y1, delta2 + z1, x2, y2, z2, crep),
ohy1z1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep,
delta2*2.0,
(x1, -delta2 + y1, z1, -delta2 + x2, y2, z2, crep),
(x1, delta2 + y1, z1, -delta2 + x2, y2, z2, crep),
(x1, -delta2 + y1, z1, delta2 + x2, y2, z2, crep),
(x1, delta2 + y1, z1, delta2 + x2, y2, z2, crep),
ohy1x2,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep,
delta2*2.0,
(x1, -delta2 + y1, z1, x2, -delta2 + y2, z2, crep),
(x1, delta2 + y1, z1, x2, -delta2 + y2, z2, crep),
(x1, -delta2 + y1, z1, x2, delta2 + y2, z2, crep),
(x1, delta2 + y1, z1, x2, delta2 + y2, z2, crep),
ohy1y2,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep,
delta2*2.0,
(x1, -delta2 + y1, z1, x2, y2, -delta2 + z2, crep),
(x1, delta2 + y1, z1, x2, y2, -delta2 + z2, crep),
(x1, -delta2 + y1, z1, x2, y2, delta2 + z2, crep),
(x1, delta2 + y1, z1, x2, y2, delta2 + z2, crep),
ohy1z2,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep, delta2*2.0,
(x1, y1, -delta2 + z1, x2, y2, z2, crep),
(x1, y1, z1, x2, y2, z2, crep),
(x1, y1, delta2 + z1, x2, y2, z2, crep),
dhz1z1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep,
delta2*2.0,
(x1, y1, -delta2 + z1, -delta2 + x2, y2, z2, crep),
(x1, y1, delta2 + z1, -delta2 + x2, y2, z2, crep),
(x1, y1, -delta2 + z1, delta2 + x2, y2, z2, crep),
(x1, y1, delta2 + z1, delta2 + x2, y2, z2, crep),
ohz1x2,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep,
delta2*2.0,
(x1, y1, -delta2 + z1, x2, -delta2 + y2, z2, crep),
(x1, y1, delta2 + z1, x2, -delta2 + y2, z2, crep),
(x1, y1, -delta2 + z1, x2, delta2 + y2, z2, crep),
(x1, y1, delta2 + z1, x2, delta2 + y2, z2, crep),
ohz1y2,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep,
delta2*2.0,
(x1, y1, -delta2 + z1, x2, y2, -delta2 + z2, crep),
(x1, y1, delta2 + z1, x2, y2, -delta2 + z2, crep),
(x1, y1, -delta2 + z1, x2, y2, delta2 + z2, crep),
(x1, y1, delta2 + z1, x2, y2, delta2 + z2, crep),
ohz1z2,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep, delta2*2.0,
(x1, y1, z1, -delta2 + x2, y2, z2, crep),
(x1, y1, z1, x2, y2, z2, crep),
(x1, y1, z1, delta2 + x2, y2, z2, crep),
dhx2x2,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep,
delta2*2.0,
(x1, y1, z1, -delta2 + x2, -delta2 + y2, z2, crep),
(x1, y1, z1, delta2 + x2, -delta2 + y2, z2, crep),
(x1, y1, z1, -delta2 + x2, delta2 + y2, z2, crep),
(x1, y1, z1, delta2 + x2, delta2 + y2, z2, crep),
ohx2y2,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep,
delta2*2.0,
(x1, y1, z1, -delta2 + x2, y2, -delta2 + z2, crep),
(x1, y1, z1, delta2 + x2, y2, -delta2 + z2, crep),
(x1, y1, z1, -delta2 + x2, y2, delta2 + z2, crep),
(x1, y1, z1, delta2 + x2, y2, delta2 + z2, crep),
ohx2z2,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep, delta2*2.0,
(x1, y1, z1, x2, -delta2 + y2, z2, crep),
(x1, y1, z1, x2, y2, z2, crep),
(x1, y1, z1, x2, delta2 + y2, z2, crep),
dhy2y2,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep,
delta2*2.0,
(x1, y1, z1, x2, -delta2 + y2, -delta2 + z2, crep),
(x1, y1, z1, x2, delta2 + y2, -delta2 + z2, crep),
(x1, y1, z1, x2, -delta2 + y2, delta2 + z2, crep),
(x1, y1, z1, x2, delta2 + y2, delta2 + z2, crep),
ohy2z2,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Erep, delta2*2.0,
(x1, y1, z1, x2, y2, -delta2 + z2, crep),
(x1, y1, z1, x2, y2, z2, crep),
(x1, y1, z1, x2, y2, delta2 + z2, crep),
dhz2z2,
index);
