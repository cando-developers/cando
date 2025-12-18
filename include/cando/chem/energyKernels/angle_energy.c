void angle_energy(double kt, double t0, size_t i3x1, size_t i3x2, size_t i3x3, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec) {
  constexpr size_t PositionSize = 9;
  /* !BASE */
  double cos_theta;
  double cse_p1_t1_g330;
  double cse_p1_t2_g331;
  double cse_p1_t3_g332;
  double dot;
  double dtheta;
  double energy;
  double n1;
  double n1_sq;
  double n2;
  double n2_sq;
  double sin_theta;
  double theta;
  double vx1;
  double vx2;
  double vy1;
  double vy2;
  double vz1;
  double vz2;
  DOUBLE x1 = position[i3x1 + 0];
  DOUBLE y1 = position[i3x1 + 1];
  DOUBLE z1 = position[i3x1 + 2];
  DOUBLE x2 = position[i3x2 + 0];
  DOUBLE y2 = position[i3x2 + 1];
  DOUBLE z2 = position[i3x2 + 2];
  DOUBLE x3 = position[i3x3 + 0];
  DOUBLE y3 = position[i3x3 + 1];
  DOUBLE z3 = position[i3x3 + 2];
  {
    /* !BASE */
    cse_p1_t1_g330 = (-(x2));
    cse_p1_t2_g331 = (-(y2));
    cse_p1_t3_g332 = (-(z2));
    vx1 = (cse_p1_t1_g330 + x1);
    vy1 = (cse_p1_t2_g331 + y1);
    vz1 = (cse_p1_t3_g332 + z1);
    vx2 = (cse_p1_t1_g330 + x3);
    vy2 = (cse_p1_t2_g331 + y3);
    vz2 = (cse_p1_t3_g332 + z3);
    dot = ((vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2));
    n1_sq = ((vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1));
    n2_sq = ((vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2));
    n1 = sqrt(n1_sq);
    n2 = sqrt(n2_sq);
    cos_theta = (dot * pow((n1 * n2), -1));
    if (fabs(cos_theta)>(1.0-VERYSMALL)) linear_angle_error();
    sin_theta = sqrt((1.0000000000000000     + (-(((cos_theta) * (cos_theta))))));
    theta = acos(cos_theta);
    dtheta = (theta + (-(t0)));
    energy = (dtheta * dtheta * kt);
    *energy_accumulate += energy;
  }
}
void angle_energy_fd(double kt, double t0, size_t i3x1, size_t i3x2, size_t i3x3, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec)
{
  angle_energy(kt, t0, i3x1, i3x2, i3x3, position, energy_accumulate, force, hessian, dvec, hdvec);
}

