void dihedral_restraint_energy(double kdh, double phi0, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec) {
  constexpr size_t PositionSize = 12;
  /* !BASE */
  double cosphi;
  double deltaphi;
  double dot_n1n2;
  double dx12;
  double dx23;
  double dx34;
  double dy12;
  double dy23;
  double dy34;
  double dz12;
  double dz23;
  double dz34;
  double energy;
  double inv_n1n2;
  double n1;
  double n1x;
  double n1y;
  double n1z;
  double n1_2;
  double n2;
  double n2x;
  double n2y;
  double n2z;
  double n2_2;
  double phi;
  double r12;
  double r12_2;
  double r23;
  double r23_2;
  double r34;
  double r34_2;
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
    dx12 = (x2 + (-(x1)));
    dy12 = (y2 + (-(y1)));
    dz12 = (z2 + (-(z1)));
    dx23 = (x3 + (-(x2)));
    dy23 = (y3 + (-(y2)));
    dz23 = (z3 + (-(z2)));
    dx34 = (x4 + (-(x3)));
    dy34 = (y4 + (-(y3)));
    dz34 = (z4 + (-(z3)));
    r12_2 = ((dx12 * dx12) + (dy12 * dy12) + (dz12 * dz12));
    r23_2 = ((dx23 * dx23) + (dy23 * dy23) + (dz23 * dz23));
    r34_2 = ((dx34 * dx34) + (dy34 * dy34) + (dz34 * dz34));
    r12 = sqrt(r12_2);
    r23 = sqrt(r23_2);
    r34 = sqrt(r34_2);
    n1x = ((dy12 * dz23) + (-((dy23 * dz12))));
    n1y = ((dx23 * dz12) + (-((dx12 * dz23))));
    n1z = ((dx12 * dy23) + (-((dx23 * dy12))));
    n2x = ((dy23 * dz34) + (-((dy34 * dz23))));
    n2y = ((dx34 * dz23) + (-((dx23 * dz34))));
    n2z = ((dx23 * dy34) + (-((dx34 * dy23))));
    n1_2 = ((n1x * n1x) + (n1y * n1y) + (n1z * n1z));
    n2_2 = ((n2x * n2x) + (n2y * n2y) + (n2z * n2z));
    n1 = sqrt(n1_2);
    n2 = sqrt(n2_2);
    dot_n1n2 = ((n1x * n2x) + (n1y * n2y) + (n1z * n2z));
    inv_n1n2 = pow((n1 * n2), -1);
    cosphi = (dot_n1n2 * inv_n1n2);
    phi = acos(cosphi);
    deltaphi = (phi + (-(phi0)));
    energy = (deltaphi * deltaphi * kdh);
    *energy_accumulate += energy;
  }
}
void dihedral_restraint_energy_fd(double kdh, double phi0, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec)
{
  dihedral_restraint_energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, energy_accumulate, force, hessian, dvec, hdvec);
}

