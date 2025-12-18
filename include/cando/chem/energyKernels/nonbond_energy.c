template <typename HESSIAN>
void nonbond_energy(double a, double b, double qq, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  constexpr size_t PositionSize = 6;
  /* !BASE */
  double dx;
  double dy;
  double dz;
  double energy;
  double e_coul;
  double e_lj;
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
    r = sqrt(r2);
    invr = pow(r, -1);
    invr2 = (invr * invr);
    invr6 = (invr2 * invr2 * invr2);
    e_lj = ((a * invr6 * invr6) + (-((b * invr6))));
    e_coul = (invr * qq);
    energy = (e_coul + e_lj);
    *energy_accumulate += energy;
  }
}
template <typename HESSIAN>
void nonbond_energy_fd(double a, double b, double qq, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  nonbond_energy(a, b, qq, i3x1, i3x2, position, energy_accumulate, force, hessian, dvec, hdvec);
}

