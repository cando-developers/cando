void chiral_restraint_energy(double k, double co, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec) {
  constexpr size_t PositionSize = 12;
  /* !BASE */
  double cse_p1_t1_g52956;
  double cse_p1_t2_g52957;
  double cse_p1_t3_g52958;
  double cx;
  double cy;
  double cz;
  double denom;
  double dx13;
  double dx23;
  double dx43;
  double dy13;
  double dy23;
  double dy43;
  double dz13;
  double dz23;
  double dz43;
  double energy;
  double q;
  double r13;
  double r13_2;
  double r23;
  double r23_2;
  double r43;
  double r43_2;
  double v;
  DOUBLE x1 = position[i3x1 + 0];
  DOUBLE y1 = position[i3x1 + 1];
  DOUBLE z1 = position[i3x1 + 2];
  DOUBLE x2 = position[i3x2 + 0];
  DOUBLE y2 = position[i3x2 + 1];
  DOUBLE z2 = position[i3x2 + 2];
  DOUBLE x3 = position[i3x3 + 0];
  DOUBLE y3 = position[i3x3 + 1];
  DOUBLE z3 = position[i3x3 + 2];
  DOUBLE x4 = position[i3x4 + 0];
  DOUBLE y4 = position[i3x4 + 1];
  DOUBLE z4 = position[i3x4 + 2];
  {
    /* !BASE */
    cse_p1_t1_g52956 = (-(x3));
    cse_p1_t2_g52957 = (-(y3));
    cse_p1_t3_g52958 = (-(z3));
    dx13 = (cse_p1_t1_g52956 + x1);
    dy13 = (cse_p1_t2_g52957 + y1);
    dz13 = (cse_p1_t3_g52958 + z1);
    dx23 = (cse_p1_t1_g52956 + x2);
    dy23 = (cse_p1_t2_g52957 + y2);
    dz23 = (cse_p1_t3_g52958 + z2);
    dx43 = (cse_p1_t1_g52956 + x4);
    dy43 = (cse_p1_t2_g52957 + y4);
    dz43 = (cse_p1_t3_g52958 + z4);
    r13_2 = ((dx13 * dx13) + (dy13 * dy13) + (dz13 * dz13));
    r23_2 = ((dx23 * dx23) + (dy23 * dy23) + (dz23 * dz23));
    r43_2 = ((dx43 * dx43) + (dy43 * dy43) + (dz43 * dz43));
    r13 = sqrt(r13_2);
    r23 = sqrt(r23_2);
    r43 = sqrt(r43_2);
    cx = ((dy13 * dz23) + (-((dy23 * dz13))));
    cy = ((dx23 * dz13) + (-((dx13 * dz23))));
    cz = ((dx13 * dy23) + (-((dx23 * dy13))));
    v = ((cx * dx43) + (cy * dy43) + (cz * dz43));
    denom = (r13 * r23 * r43);
    q = (v * pow(denom, -1));
    energy = (k * ((co + q) * (co + q) * (co + q)));
    *energy_accumulate += energy;
  }
}
void chiral_restraint_energy_fd(double k, double co, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec)
{
  chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, energy_accumulate, force, hessian, dvec, hdvec);
}

