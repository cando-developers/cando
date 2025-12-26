void anchor_hessian(double ka, double xa, double ya, double za, size_t i3x1, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec) {
  constexpr size_t PositionSize = 3;
  /* !BASE */
  double cse_p1_t1_g53018;
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
    cse_p1_t1_g53018 = (2.0000000000000000     * ka);
    g_x1 = (cse_p1_t1_g53018 * dx);
    KernelGradientAcc(i3x1, 0, g_x1);
    g_y1 = (cse_p1_t1_g53018 * dy);
    KernelGradientAcc(i3x1, 1, g_y1);
    g_z1 = (cse_p1_t1_g53018 * dz);
    KernelGradientAcc(i3x1, 2, g_z1);
  }
}
void anchor_hessian_fd(double ka, double xa, double ya, double za, size_t i3x1, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 3;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  const double invh2 = 1.0/((h*h));
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
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 0, i3x1, 0, hval);
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
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 1, i3x1, 1, hval);
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
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 2, i3x1, 2, hval);
  }
  {
    double saved_i = position[i3x1 + 1];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x1 + 1] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_mm, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 1, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x1 + 2];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x1 + 2] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_mm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 2, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x1 + 2];
    double saved_j = position[i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x1 + 2] = saved_i + h; position[i3x1 + 1] = saved_j + h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    anchor_energy(ka, xa, ya, za, i3x1, position, &e_mm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i; position[i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 2, i3x1, 1, hval);
  }
}

