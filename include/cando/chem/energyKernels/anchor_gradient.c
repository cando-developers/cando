void anchor_gradient(double ka, double xa, double ya, double za, size_t i3x1, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec) {
  constexpr size_t PositionSize = 3;
  /* !BASE */
  double cse_p1_t1_g53017;
  double dx;
  double dy;
  double dz;
  double energy;
  double g_x1;
  double g_y1;
  double g_z1;
  double r2;
  DOUBLE x1 = position[i3x1 + 0];
  DOUBLE y1 = position[i3x1 + 1];
  DOUBLE z1 = position[i3x1 + 2];
  {
    /* !BASE */
    dx = (x1 + (-(xa)));
    dy = (y1 + (-(ya)));
    dz = (z1 + (-(za)));
    r2 = ((dx * dx) + (dy * dy) + (dz * dz));
    energy = (ka * r2);
    *energy_accumulate += energy;
    cse_p1_t1_g53017 = (2.0000000000000000     * ka);
    g_x1 = (cse_p1_t1_g53017 * dx);
    KernelGradientAcc(i3x1, 0, g_x1);
    g_y1 = (cse_p1_t1_g53017 * dy);
    KernelGradientAcc(i3x1, 1, g_y1);
    g_z1 = (cse_p1_t1_g53017 * dz);
    KernelGradientAcc(i3x1, 2, g_z1);
  }
}
void anchor_gradient_fd(double ka, double xa, double ya, double za, size_t i3x1, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 3;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  double e0 = 0.0;
  anchor_energy(ka, xa, ya, za, i3x1, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
}

