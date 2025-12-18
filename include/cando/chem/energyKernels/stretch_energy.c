void stretch_energy(double kb, double r0, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec) {
  constexpr size_t PositionSize = 6;
  /* !BASE */
  double dr;
  double dx;
  double dy;
  double dz;
  double energy;
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
    r = sqrt(r2);
    dr = (r + (-(r0)));
    energy = (dr * dr * kb);
    *energy_accumulate += energy;
  }
}
void stretch_energy_fd(double kb, double r0, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec)
{
  stretch_energy(kb, r0, i3x1, i3x2, position, energy_accumulate, force, hessian, dvec, hdvec);
}

