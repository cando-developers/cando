void dihedral_energy(double v, double n, double sinphase, double cosphase, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec) {
  constexpr size_t PositionSize = 12;
  /* !BASE */
  double c1x;
  double c1y;
  double c1z;
  double c1_sq;
  double c2x;
  double c2y;
  double c2z;
  double c2_sq;
  double cos_angle;
  double cos_nphi;
  double dot12;
  double dot23;
  double energy;
  double nphi;
  double phi;
  double sin_angle;
  double sin_nphi;
  double t1;
  double t2;
  double v1x;
  double v1y;
  double v1z;
  double v2x;
  double v2y;
  double v2z;
  double v2_len;
  double v2_sq;
  double v3x;
  double v3y;
  double v3z;
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
    v1x = (x2 + (-(x1)));
    v1y = (y2 + (-(y1)));
    v1z = (z2 + (-(z1)));
    v2x = (x3 + (-(x2)));
    v2y = (y3 + (-(y2)));
    v2z = (z3 + (-(z2)));
    v3x = (x4 + (-(x3)));
    v3y = (y4 + (-(y3)));
    v3z = (z4 + (-(z3)));
    c1x = ((v1y * v2z) + (-((v1z * v2y))));
    c1y = ((v1z * v2x) + (-((v1x * v2z))));
    c1z = ((v1x * v2y) + (-((v1y * v2x))));
    c2x = ((v2y * v3z) + (-((v2z * v3y))));
    c2y = ((v2z * v3x) + (-((v2x * v3z))));
    c2z = ((v2x * v3y) + (-((v2y * v3x))));
    c1_sq = ((c1x * c1x) + (c1y * c1y) + (c1z * c1z));
    c2_sq = ((c2x * c2x) + (c2y * c2y) + (c2z * c2z));
    v2_sq = ((v2x * v2x) + (v2y * v2y) + (v2z * v2z));
    v2_len = sqrt(v2_sq);
    dot12 = ((v1x * v2x) + (v1y * v2y) + (v1z * v2z));
    dot23 = ((v2x * v3x) + (v2y * v3y) + (v2z * v3z));
    t1 = (v2_len * ((c2x * v1x) + (c2y * v1y) + (c2z * v1z)));
    t2 = ((c1x * c2x) + (c1y * c2y) + (c1z * c2z));
    phi = atan2(t1, t2);
    nphi = (n * phi);
    sin_nphi = sin(nphi);
    cos_nphi = cos(nphi);
    cos_angle = ((cosphase * cos_nphi) + (sinphase * sin_nphi));
    sin_angle = ((cosphase * sin_nphi) + (-((cos_nphi * sinphase))));
    energy = (v * (1.0000000000000000     + cos_angle));
    *energy_accumulate += energy;
  }
}
void dihedral_energy_fd(double v, double n, double sinphase, double cosphase, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec)
{
  dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, energy_accumulate, force, hessian, dvec, hdvec);
}

