template <typename HESSIAN>
void nonbond_gradient(double a, double b, double qq, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  constexpr size_t PositionSize = 6;
  /* !BASE */
  double cse_p1_t1_g1067;
  double cse_p1_t1_g1069;
  double cse_p1_t1_g1070;
  double cse_p1_t1_g1071;
  double cse_p1_t1_g1072;
  double cse_p1_t1_g1073;
  double cse_p1_t1_g1074;
  double cse_p1_t1_g1075;
  double cse_p1_t1_g1077;
  double cse_p1_t1_g1078;
  double cse_p1_t2_g1068;
  double cse_p2_t1_g1079;
  double cse_p2_t1_g1080;
  double cse_p2_t1_g1081;
  double dx;
  double dy;
  double dz;
  double energy;
  double e_coul;
  double e_lj;
  double g_x1;
  double g_x2;
  double g_y1;
  double g_y2;
  double g_z1;
  double g_z2;
  double invr;
  double invr2;
  double invr6;
  double r;
  double r2;
  DOUBLE x1 = position[i3x1 + 0];
  DOUBLE y1 = position[i3x1 + 1];
  DOUBLE z1 = position[i3x1 + 2];
  DOUBLE x2 = position[i3x2 + 0];
  DOUBLE y2 = position[i3x2 + 1];
  DOUBLE z2 = position[i3x2 + 2];
  {
    /* !BASE */
    dx = (x1 + (-(x2)));
    dy = (y1 + (-(y2)));
    dz = (z1 + (-(z2)));
    r2 = ((dx * dx) + (dy * dy) + (dz * dz));
    cse_p1_t2_g1068 = pow(r2, -0.50000000000000000    );
    r = sqrt(r2);
    cse_p1_t1_g1067 = pow(r, -2);
    invr = pow(r, -1);
    invr2 = (invr * invr);
    cse_p2_t1_g1081 = (invr2 * invr2);
    invr6 = (cse_p2_t1_g1081 * invr2);
    cse_p2_t1_g1079 = (a * invr6);
    e_lj = ((cse_p2_t1_g1079 * invr6) + (-((b * invr6))));
    e_coul = (invr * qq);
    energy = (e_coul + e_lj);
    *energy_accumulate += energy;
    cse_p2_t1_g1080 = (cse_p1_t1_g1067 * cse_p1_t2_g1068);
    cse_p1_t1_g1069 = (cse_p2_t1_g1080 * cse_p2_t1_g1081 * invr);
    cse_p1_t1_g1070 = (cse_p1_t1_g1069 * cse_p2_t1_g1079);
    cse_p1_t1_g1071 = (cse_p2_t1_g1080 * qq);
    cse_p1_t1_g1072 = (b * cse_p1_t1_g1069);
    cse_p1_t1_g1073 = (-(cse_p1_t1_g1071));
    cse_p1_t1_g1074 = (-12.000000000000000     * cse_p1_t1_g1070);
    cse_p1_t1_g1078 = (6.0000000000000000     * cse_p1_t1_g1072);
    g_x1 = ((cse_p1_t1_g1073 * dx) + (cse_p1_t1_g1074 * dx) + (cse_p1_t1_g1078 * dx));
    KernelGradientAcc(i3x1, 0, g_x1);
    g_y1 = ((cse_p1_t1_g1073 * dy) + (cse_p1_t1_g1074 * dy) + (cse_p1_t1_g1078 * dy));
    KernelGradientAcc(i3x1, 1, g_y1);
    g_z1 = ((cse_p1_t1_g1073 * dz) + (cse_p1_t1_g1074 * dz) + (cse_p1_t1_g1078 * dz));
    KernelGradientAcc(i3x1, 2, g_z1);
    cse_p1_t1_g1075 = (-6.0000000000000000     * cse_p1_t1_g1072);
    cse_p1_t1_g1077 = (12.000000000000000     * cse_p1_t1_g1070);
    g_x2 = ((cse_p1_t1_g1071 * dx) + (cse_p1_t1_g1075 * dx) + (cse_p1_t1_g1077 * dx));
    KernelGradientAcc(i3x2, 0, g_x2);
    g_y2 = ((cse_p1_t1_g1071 * dy) + (cse_p1_t1_g1075 * dy) + (cse_p1_t1_g1077 * dy));
    KernelGradientAcc(i3x2, 1, g_y2);
    g_z2 = ((cse_p1_t1_g1071 * dz) + (cse_p1_t1_g1075 * dz) + (cse_p1_t1_g1077 * dz));
    KernelGradientAcc(i3x2, 2, g_z2);
  }
}
template <typename HESSIAN>
void nonbond_gradient_fd(double a, double b, double qq, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 6;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  double e0 = 0.0;
  nonbond_energy(a, b, qq, i3x1, i3x2, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    nonbond_energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    nonbond_energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    nonbond_energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    nonbond_energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    nonbond_energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    nonbond_energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    nonbond_energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    nonbond_energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    nonbond_energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    nonbond_energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    nonbond_energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    nonbond_energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
}

