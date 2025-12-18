template <typename HESSIAN>
void nonbond_dd_cutoff_energy(double a, double b, double qq, double invdd, double r_switch, double r_switch2, double r_cut, double r_cut2, double inv_range, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  constexpr size_t PositionSize = 6;
  /* !BASE */
  double cse_p1_t1_g28332;
  double drs;
  double dx;
  double dy;
  double dz;
  double energy;
  double e_base;
  double e_coul;
  double e_lj;
  double invr;
  double invr2;
  double invr3;
  double invr6;
  double r;
  double r2;
  double s;
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
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
    if ((r2 < r_cut2)) {
      /* !CUT */
      invr = pow(r2, -0.50000000000000000    );
      r = (invr * r2);
      invr2 = (invr * invr);
      invr3 = (invr2 * invr);
      cse_p1_t1_g28332 = (invr2 * invr2);
      invr6 = (cse_p1_t1_g28332 * invr2);
      e_lj = ((a * invr6 * invr6) + (-((b * invr6))));
      e_coul = (invdd * invr2 * qq);
      e_base = (e_coul + e_lj);
      if ((r2 < r_switch2)) {
        /* !SWITCH */
        energy = e_base;
        *energy_accumulate += energy;
      } else {
        /* !SKIN */
        drs = (r + (-(r_switch)));
        t1 = (drs * inv_range);
        t2 = (t1 * t1);
        t3 = (t1 * t2);
        t4 = (t2 * t2);
        t5 = (t2 * t3);
        s = (1.0000000000000000     + (-10.000000000000000     * t3) + (-6.0000000000000000     * t5) + (15.000000000000000     * t4));
        energy = (e_base * s);
        *energy_accumulate += energy;
      }
    }
  }
}
template <typename HESSIAN>
void nonbond_dd_cutoff_energy_fd(double a, double b, double qq, double invdd, double r_switch, double r_switch2, double r_cut, double r_cut2, double inv_range, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  nonbond_dd_cutoff_energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, energy_accumulate, force, hessian, dvec, hdvec);
}

