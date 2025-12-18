void stretch_gradient(double kb, double r0, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec) {
  constexpr size_t PositionSize = 6;
  /* !BASE */
  double cse_p1_t1_g325;
  double cse_p1_t1_g326;
  double cse_p1_t1_g327;
  double cse_p1_t1_g328;
  double cse_p51_t1_g329;
  double dr;
  double dx;
  double dy;
  double dz;
  double energy;
  double g_x1;
  double g_x2;
  double g_y1;
  double g_y2;
  double g_z1;
  double g_z2;
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
    dx = (x2 + (-(x1)));
    dy = (y2 + (-(y1)));
    dz = (z2 + (-(z1)));
    r2 = ((dx * dx) + (dy * dy) + (dz * dz));
    cse_p1_t1_g325 = pow(r2, -0.50000000000000000    );
    r = sqrt(r2);
    dr = (r + (-(r0)));
    cse_p1_t1_g328 = (dr * kb);
    energy = (cse_p1_t1_g328 * dr);
    *energy_accumulate += energy;
    cse_p1_t1_g326 = (cse_p1_t1_g325 * cse_p1_t1_g328);
    cse_p51_t1_g329 = (-2.0000000000000000     * cse_p1_t1_g326);
    g_x1 = (cse_p51_t1_g329 * dx);
    KernelGradientAcc(i3x1, 0, g_x1);
    g_y1 = (cse_p51_t1_g329 * dy);
    KernelGradientAcc(i3x1, 1, g_y1);
    g_z1 = (cse_p51_t1_g329 * dz);
    KernelGradientAcc(i3x1, 2, g_z1);
    cse_p1_t1_g327 = (2.0000000000000000     * cse_p1_t1_g326);
    g_x2 = (cse_p1_t1_g327 * dx);
    KernelGradientAcc(i3x2, 0, g_x2);
    g_y2 = (cse_p1_t1_g327 * dy);
    KernelGradientAcc(i3x2, 1, g_y2);
    g_z2 = (cse_p1_t1_g327 * dz);
    KernelGradientAcc(i3x2, 2, g_z2);
  }
}
void stretch_gradient_fd(double kb, double r0, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 6;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  double e0 = 0.0;
  stretch_energy(kb, r0, i3x1, i3x2, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    stretch_energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    stretch_energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    stretch_energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    stretch_energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    stretch_energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    stretch_energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    stretch_energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    stretch_energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    stretch_energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    stretch_energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    stretch_energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    stretch_energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
}

