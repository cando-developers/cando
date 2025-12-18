void angle_gradient(double kt, double t0, size_t i3x1, size_t i3x2, size_t i3x3, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec) {
  constexpr size_t PositionSize = 9;
  /* !BASE */
  double cos_theta;
  double cse_p101_t1_g353;
  double cse_p1_t10_g342;
  double cse_p1_t11_g343;
  double cse_p1_t12_g344;
  double cse_p1_t1_g333;
  double cse_p1_t1_g345;
  double cse_p1_t1_g346;
  double cse_p1_t1_g347;
  double cse_p1_t1_g349;
  double cse_p1_t1_g350;
  double cse_p1_t2_g334;
  double cse_p1_t3_g335;
  double cse_p1_t4_g336;
  double cse_p1_t5_g337;
  double cse_p1_t6_g338;
  double cse_p1_t7_g339;
  double cse_p1_t8_g340;
  double cse_p1_t9_g341;
  double cse_p2_t1_g351;
  double cse_p51_t1_g352;
  double dot;
  double dtheta;
  double energy;
  double g_x1;
  double g_x2;
  double g_x3;
  double g_y1;
  double g_y2;
  double g_y3;
  double g_z1;
  double g_z2;
  double g_z3;
  double n1;
  double n1_sq;
  double n2;
  double n2_sq;
  double sin_theta;
  double theta;
  double vx1;
  double vx2;
  double vy1;
  double vy2;
  double vz1;
  double vz2;
  DOUBLE x1 = position[i3x1 + 0];
  DOUBLE y1 = position[i3x1 + 1];
  DOUBLE z1 = position[i3x1 + 2];
  DOUBLE x2 = position[i3x2 + 0];
  DOUBLE y2 = position[i3x2 + 1];
  DOUBLE z2 = position[i3x2 + 2];
  DOUBLE x3 = position[i3x3 + 0];
  DOUBLE y3 = position[i3x3 + 1];
  DOUBLE z3 = position[i3x3 + 2];
  {
    /* !BASE */
    cse_p1_t4_g336 = (-(x2));
    cse_p1_t5_g337 = (-(y2));
    cse_p1_t6_g338 = (-(z2));
    vx1 = (cse_p1_t4_g336 + x1);
    vy1 = (cse_p1_t5_g337 + y1);
    vz1 = (cse_p1_t6_g338 + z1);
    vx2 = (cse_p1_t4_g336 + x3);
    vy2 = (cse_p1_t5_g337 + y3);
    vz2 = (cse_p1_t6_g338 + z3);
    dot = ((vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2));
    n1_sq = ((vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1));
    cse_p1_t11_g343 = pow(n1_sq, -0.50000000000000000    );
    n2_sq = ((vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2));
    cse_p1_t12_g344 = pow(n2_sq, -0.50000000000000000    );
    n1 = sqrt(n1_sq);
    n2 = sqrt(n2_sq);
    cse_p1_t1_g333 = (n1 * n2);
    cse_p1_t7_g339 = pow(cse_p1_t1_g333, -1);
    cse_p1_t8_g340 = pow(cse_p1_t1_g333, -2);
    cse_p101_t1_g353 = (cse_p1_t8_g340 * dot);
    cos_theta = (cse_p1_t7_g339 * dot);
    cse_p1_t2_g334 = (1.0000000000000000     + (-(((cos_theta) * (cos_theta)))));
    cse_p1_t3_g335 = (-(((cos_theta) * (cos_theta))));
    cse_p1_t9_g341 = pow(cse_p1_t2_g334, -0.50000000000000000    );
    cse_p1_t10_g342 = ((cos_theta) * (cos_theta));
    if (fabs(cos_theta)>(1.0-VERYSMALL)) linear_angle_error();
    sin_theta = sqrt(cse_p1_t2_g334);
    theta = acos(cos_theta);
    dtheta = (theta + (-(t0)));
    cse_p2_t1_g351 = (dtheta * kt);
    energy = (cse_p2_t1_g351 * dtheta);
    *energy_accumulate += energy;
    cse_p1_t1_g345 = (2.0000000000000000     * cse_p2_t1_g351);
    cse_p1_t1_g346 = (cse_p1_t11_g343 * n2);
    cse_p1_t1_g349 = (cse_p101_t1_g353 * cse_p1_t1_g346);
    cse_p51_t1_g352 = (cse_p1_t1_g345 * cse_p1_t9_g341);
    g_x1 = (-((cse_p51_t1_g352 * ((cse_p1_t7_g339 * vx2) + (-((cse_p1_t1_g349 * vx1)))))));
    KernelGradientAcc(i3x1, 0, g_x1);
    g_y1 = (-((cse_p51_t1_g352 * ((cse_p1_t7_g339 * vy2) + (-((cse_p1_t1_g349 * vy1)))))));
    KernelGradientAcc(i3x1, 1, g_y1);
    g_z1 = (-((cse_p51_t1_g352 * ((cse_p1_t7_g339 * vz2) + (-((cse_p1_t1_g349 * vz1)))))));
    KernelGradientAcc(i3x1, 2, g_z1);
    cse_p1_t1_g347 = (cse_p1_t12_g344 * n1);
    g_x2 = (-((cse_p51_t1_g352 * ((cse_p1_t7_g339 * ((-(vx1)) + (-(vx2)))) + (-((cse_p101_t1_g353 * ((-((cse_p1_t1_g346 * vx1))) + (-((cse_p1_t1_g347 * vx2)))))))))));
    KernelGradientAcc(i3x2, 0, g_x2);
    g_y2 = (-((cse_p51_t1_g352 * ((cse_p1_t7_g339 * ((-(vy1)) + (-(vy2)))) + (-((cse_p101_t1_g353 * ((-((cse_p1_t1_g346 * vy1))) + (-((cse_p1_t1_g347 * vy2)))))))))));
    KernelGradientAcc(i3x2, 1, g_y2);
    g_z2 = (-((cse_p51_t1_g352 * ((cse_p1_t7_g339 * ((-(vz1)) + (-(vz2)))) + (-((cse_p101_t1_g353 * ((-((cse_p1_t1_g346 * vz1))) + (-((cse_p1_t1_g347 * vz2)))))))))));
    KernelGradientAcc(i3x2, 2, g_z2);
    cse_p1_t1_g350 = (cse_p101_t1_g353 * cse_p1_t1_g347);
    g_x3 = (-((cse_p51_t1_g352 * ((cse_p1_t7_g339 * vx1) + (-((cse_p1_t1_g350 * vx2)))))));
    KernelGradientAcc(i3x3, 0, g_x3);
    g_y3 = (-((cse_p51_t1_g352 * ((cse_p1_t7_g339 * vy1) + (-((cse_p1_t1_g350 * vy2)))))));
    KernelGradientAcc(i3x3, 1, g_y3);
    g_z3 = (-((cse_p51_t1_g352 * ((cse_p1_t7_g339 * vz1) + (-((cse_p1_t1_g350 * vz2)))))));
    KernelGradientAcc(i3x3, 2, g_z3);
  }
}
void angle_gradient_fd(double kt, double t0, size_t i3x1, size_t i3x2, size_t i3x3, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 9;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  double e0 = 0.0;
  angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
  {
    double saved = position[i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 0] = saved + h;
    angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved - h;
    angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 0, d);
  }
  {
    double saved = position[i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 1] = saved + h;
    angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved - h;
    angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 1, d);
  }
  {
    double saved = position[i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 2] = saved + h;
    angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved - h;
    angle_energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 2, d);
  }
}

