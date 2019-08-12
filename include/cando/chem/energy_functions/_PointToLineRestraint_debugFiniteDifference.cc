// FINITE_DIFFERENCE_TEST 
TEST_FORCE( _evaluateEnergyOnly_PointToLineRestraint,  delta2*2.0,
(-delta2 + x1, y1, z1, xa, ya, za, xb, yb, zb, ka, ra),
(delta2 + x1, y1, z1, xa, ya, za, xb, yb, zb, ka, ra),
fx1,
index);
TEST_FORCE( _evaluateEnergyOnly_PointToLineRestraint,  delta2*2.0,
(x1, -delta2 + y1, z1, xa, ya, za, xb, yb, zb, ka, ra),
(x1, delta2 + y1, z1, xa, ya, za, xb, yb, zb, ka, ra),
fy1,
index);
TEST_FORCE( _evaluateEnergyOnly_PointToLineRestraint,  delta2*2.0,
(x1, y1, -delta2 + z1, xa, ya, za, xb, yb, zb, ka, ra),
(x1, y1, delta2 + z1, xa, ya, za, xb, yb, zb, ka, ra),
fz1,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_PointToLineRestraint, delta2*2.0,
(-delta2 + x1, y1, z1, xa, ya, za, xb, yb, zb, ka, ra),
(x1, y1, z1, xa, ya, za, xb, yb, zb, ka, ra),
(delta2 + x1, y1, z1, xa, ya, za, xb, yb, zb, ka, ra),
dhx1x1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_PointToLineRestraint,
delta2*2.0,
(-delta2 + x1, -delta2 + y1, z1, xa, ya, za, xb, yb, zb, ka, ra),
(delta2 + x1, -delta2 + y1, z1, xa, ya, za, xb, yb, zb, ka, ra),
(-delta2 + x1, delta2 + y1, z1, xa, ya, za, xb, yb, zb, ka, ra),
(delta2 + x1, delta2 + y1, z1, xa, ya, za, xb, yb, zb, ka, ra),
ohx1y1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_PointToLineRestraint,
delta2*2.0,
(-delta2 + x1, y1, -delta2 + z1, xa, ya, za, xb, yb, zb, ka, ra),
(delta2 + x1, y1, -delta2 + z1, xa, ya, za, xb, yb, zb, ka, ra),
(-delta2 + x1, y1, delta2 + z1, xa, ya, za, xb, yb, zb, ka, ra),
(delta2 + x1, y1, delta2 + z1, xa, ya, za, xb, yb, zb, ka, ra),
ohx1z1,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_PointToLineRestraint, delta2*2.0,
(x1, -delta2 + y1, z1, xa, ya, za, xb, yb, zb, ka, ra),
(x1, y1, z1, xa, ya, za, xb, yb, zb, ka, ra),
(x1, delta2 + y1, z1, xa, ya, za, xb, yb, zb, ka, ra),
dhy1y1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_PointToLineRestraint,
delta2*2.0,
(x1, -delta2 + y1, -delta2 + z1, xa, ya, za, xb, yb, zb, ka, ra),
(x1, delta2 + y1, -delta2 + z1, xa, ya, za, xb, yb, zb, ka, ra),
(x1, -delta2 + y1, delta2 + z1, xa, ya, za, xb, yb, zb, ka, ra),
(x1, delta2 + y1, delta2 + z1, xa, ya, za, xb, yb, zb, ka, ra),
ohy1z1,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_PointToLineRestraint, delta2*2.0,
(x1, y1, -delta2 + z1, xa, ya, za, xb, yb, zb, ka, ra),
(x1, y1, z1, xa, ya, za, xb, yb, zb, ka, ra),
(x1, y1, delta2 + z1, xa, ya, za, xb, yb, zb, ka, ra),
dhz1z1,
index);
