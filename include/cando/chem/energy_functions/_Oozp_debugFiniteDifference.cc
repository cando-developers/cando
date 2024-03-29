// FINITE_DIFFERENCE_TEST 
TEST_FORCE( _evaluateEnergyOnly_Oozp,  delta2*2.0,
(ENERGY_FUNCTION, -delta2 + x1, y1, z1, kb, za),
(ENERGY_FUNCTION, delta2 + x1, y1, z1, kb, za),
fx1,
index);
TEST_FORCE( _evaluateEnergyOnly_Oozp,  delta2*2.0,
(ENERGY_FUNCTION, x1, -delta2 + y1, z1, kb, za),
(ENERGY_FUNCTION, x1, delta2 + y1, z1, kb, za),
fy1,
index);
TEST_FORCE( _evaluateEnergyOnly_Oozp,  delta2*2.0,
(ENERGY_FUNCTION, x1, y1, -delta2 + z1, kb, za),
(ENERGY_FUNCTION, x1, y1, delta2 + z1, kb, za),
fz1,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Oozp, delta2*2.0,
(ENERGY_FUNCTION, -delta2 + x1, y1, z1, kb, za),
(ENERGY_FUNCTION, x1, y1, z1, kb, za),
(ENERGY_FUNCTION, delta2 + x1, y1, z1, kb, za),
dhx1x1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Oozp,
delta2*2.0,
(ENERGY_FUNCTION, -delta2 + x1, -delta2 + y1, z1, kb, za),
(ENERGY_FUNCTION, delta2 + x1, -delta2 + y1, z1, kb, za),
(ENERGY_FUNCTION, -delta2 + x1, delta2 + y1, z1, kb, za),
(ENERGY_FUNCTION, delta2 + x1, delta2 + y1, z1, kb, za),
ohx1y1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Oozp,
delta2*2.0,
(ENERGY_FUNCTION, -delta2 + x1, y1, -delta2 + z1, kb, za),
(ENERGY_FUNCTION, delta2 + x1, y1, -delta2 + z1, kb, za),
(ENERGY_FUNCTION, -delta2 + x1, y1, delta2 + z1, kb, za),
(ENERGY_FUNCTION, delta2 + x1, y1, delta2 + z1, kb, za),
ohx1z1,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Oozp, delta2*2.0,
(ENERGY_FUNCTION, x1, -delta2 + y1, z1, kb, za),
(ENERGY_FUNCTION, x1, y1, z1, kb, za),
(ENERGY_FUNCTION, x1, delta2 + y1, z1, kb, za),
dhy1y1,
index);
TEST_OFF_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Oozp,
delta2*2.0,
(ENERGY_FUNCTION, x1, -delta2 + y1, -delta2 + z1, kb, za),
(ENERGY_FUNCTION, x1, delta2 + y1, -delta2 + z1, kb, za),
(ENERGY_FUNCTION, x1, -delta2 + y1, delta2 + z1, kb, za),
(ENERGY_FUNCTION, x1, delta2 + y1, delta2 + z1, kb, za),
ohy1z1,
index);
TEST_DIAGONAL_HESSIAN( _evaluateEnergyOnly_Oozp, delta2*2.0,
(ENERGY_FUNCTION, x1, y1, -delta2 + z1, kb, za),
(ENERGY_FUNCTION, x1, y1, z1, kb, za),
(ENERGY_FUNCTION, x1, y1, delta2 + z1, kb, za),
dhz1z1,
index);
