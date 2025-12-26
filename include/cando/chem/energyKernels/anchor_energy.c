void anchor_energy(double ka, double xa, double ya, double za, size_t i3x1, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec) {
  constexpr size_t PositionSize = 3;
  /* !BASE */
  double dx;
  double dy;
  double dz;
  double energy;
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
  }
}
void anchor_energy_fd(double ka, double xa, double ya, double za, size_t i3x1, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec)
{
  anchor_energy(ka, xa, ya, za, i3x1, position, energy_accumulate, force, hessian, dvec, hdvec);
}

