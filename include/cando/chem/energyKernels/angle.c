template <typename HESSIAN>
struct Angle {
  static constexpr size_t PositionSize = 9;
  static std::string description() { return "mathkernel-angle"; };
double energy(double kt, double t0, size_t i3x1, size_t i3x2, size_t i3x3, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double energy_added = 0.0000000000000000    ;
  {
    /* !BASE */
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
      double cse_p11_t1_g27572 = (-(x2));
      double cse_p11_t2_g27573 = (-(y2));
      double cse_p11_t3_g27574 = (-(z2));
      double vx1 = (cse_p11_t1_g27572 + x1);
      double vy1 = (cse_p11_t2_g27573 + y1);
      double vz1 = (cse_p11_t3_g27574 + z1);
      double vx2 = (cse_p11_t1_g27572 + x3);
      double vy2 = (cse_p11_t2_g27573 + y3);
      double vz2 = (cse_p11_t3_g27574 + z3);
      double dot = ((vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2));
      double n1_sq = ((vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1));
      double n2_sq = ((vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2));
      double n1 = sqrt(n1_sq);
      double n2 = sqrt(n2_sq);
      double cos_theta = (dot * (1.0 / ((n1 * n2))));
      if (fabs(cos_theta)>(1.0-VERYSMALL)) linear_angle_error();
      double sin_theta = sqrt((1.0000000000000000     + (-(((cos_theta) * (cos_theta))))));
      double theta = acos(cos_theta);
      double dtheta = (theta + (-(t0)));
      double energy = (dtheta * dtheta * kt);
      energy_added += energy;
      *energy_accumulate += energy;
    }
  }
  return energy_added;
}
void energy_fd(double kt, double t0, size_t i3x1, size_t i3x2, size_t i3x3, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  energy(kt, t0, i3x1, i3x2, i3x3, position, energy_accumulate, force, hessian, dvec, hdvec);
}

double gradient(double kt, double t0, size_t i3x1, size_t i3x2, size_t i3x3, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double energy_added = 0.0000000000000000    ;
  {
    /* !BASE */
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
      double cse_p11_t1_g27583 = (-(x2));
      double cse_p11_t2_g27584 = (-(y2));
      double cse_p11_t3_g27585 = (-(z2));
      double vx1 = (cse_p11_t1_g27583 + x1);
      double vy1 = (cse_p11_t2_g27584 + y1);
      double vz1 = (cse_p11_t3_g27585 + z1);
      double vx2 = (cse_p11_t1_g27583 + x3);
      double vy2 = (cse_p11_t2_g27584 + y3);
      double vz2 = (cse_p11_t3_g27585 + z3);
      double dot = ((vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2));
      double n1_sq = ((vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1));
      double n2_sq = ((vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2));
      double cse_p1_t1_sqrt27575 = sqrt(n1_sq);
      double cse_p1_t2_invsqrt27576 = (1.0 / (cse_p1_t1_sqrt27575));
      double cse_p1_t3_sqrt27577 = sqrt(n2_sq);
      double cse_p1_t4_invsqrt27578 = (1.0 / (cse_p1_t3_sqrt27577));
      double cse_p1_t5_invr27579 = (1.0 / ((cse_p1_t1_sqrt27575 * cse_p1_t3_sqrt27577)));
      double cse_p1_t6_invr227580 = (cse_p1_t5_invr27579 * cse_p1_t5_invr27579);
      double cse_p113_t1_g27594 = (cse_p1_t6_invr227580 * dot);
      double cos_theta = (cse_p1_t5_invr27579 * dot);
      if (fabs(cos_theta)>(1.0-VERYSMALL)) linear_angle_error();
      double cse_p1_t7_sqrt27581 = sqrt((1.0000000000000000     + (-(((cos_theta) * (cos_theta))))));
      double cse_p1_t8_invsqrt27582 = (1.0 / (cse_p1_t7_sqrt27581));
      double theta = acos(cos_theta);
      double dtheta = (theta + (-(t0)));
      double cse_p11_t1_g27592 = (dtheta * kt);
      double energy = (cse_p11_t1_g27592 * dtheta);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t1_g27586 = (2.0000000000000000     * cse_p11_t1_g27592);
      double cse_p11_t1_g27587 = (cse_p1_t2_invsqrt27576 * cse_p1_t3_sqrt27577);
      double cse_p11_t1_g27590 = (cse_p113_t1_g27594 * cse_p11_t1_g27587);
      double cse_p62_t1_g27593 = (cse_p11_t1_g27586 * cse_p1_t8_invsqrt27582);
      double g_x1 = (-((cse_p62_t1_g27593 * ((cse_p1_t5_invr27579 * vx2) + (-((cse_p11_t1_g27590 * vx1)))))));
      KernelGradientAcc(i3x1, 0, g_x1);
      double g_y1 = (-((cse_p62_t1_g27593 * ((cse_p1_t5_invr27579 * vy2) + (-((cse_p11_t1_g27590 * vy1)))))));
      KernelGradientAcc(i3x1, 1, g_y1);
      double g_z1 = (-((cse_p62_t1_g27593 * ((cse_p1_t5_invr27579 * vz2) + (-((cse_p11_t1_g27590 * vz1)))))));
      KernelGradientAcc(i3x1, 2, g_z1);
      double cse_p11_t1_g27588 = (cse_p1_t1_sqrt27575 * cse_p1_t4_invsqrt27578);
      double g_x2 = (-((cse_p62_t1_g27593 * ((cse_p1_t5_invr27579 * ((-(vx1)) + (-(vx2)))) + (-((cse_p113_t1_g27594 * ((-((cse_p11_t1_g27587 * vx1))) + (-((cse_p11_t1_g27588 * vx2)))))))))));
      KernelGradientAcc(i3x2, 0, g_x2);
      double g_y2 = (-((cse_p62_t1_g27593 * ((cse_p1_t5_invr27579 * ((-(vy1)) + (-(vy2)))) + (-((cse_p113_t1_g27594 * ((-((cse_p11_t1_g27587 * vy1))) + (-((cse_p11_t1_g27588 * vy2)))))))))));
      KernelGradientAcc(i3x2, 1, g_y2);
      double g_z2 = (-((cse_p62_t1_g27593 * ((cse_p1_t5_invr27579 * ((-(vz1)) + (-(vz2)))) + (-((cse_p113_t1_g27594 * ((-((cse_p11_t1_g27587 * vz1))) + (-((cse_p11_t1_g27588 * vz2)))))))))));
      KernelGradientAcc(i3x2, 2, g_z2);
      double cse_p11_t1_g27591 = (cse_p113_t1_g27594 * cse_p11_t1_g27588);
      double g_x3 = (-((cse_p62_t1_g27593 * ((cse_p1_t5_invr27579 * vx1) + (-((cse_p11_t1_g27591 * vx2)))))));
      KernelGradientAcc(i3x3, 0, g_x3);
      double g_y3 = (-((cse_p62_t1_g27593 * ((cse_p1_t5_invr27579 * vy1) + (-((cse_p11_t1_g27591 * vy2)))))));
      KernelGradientAcc(i3x3, 1, g_y3);
      double g_z3 = (-((cse_p62_t1_g27593 * ((cse_p1_t5_invr27579 * vz1) + (-((cse_p11_t1_g27591 * vz2)))))));
      KernelGradientAcc(i3x3, 2, g_z3);
    }
  }
  return energy_added;
}
void gradient_fd(double kt, double t0, size_t i3x1, size_t i3x2, size_t i3x3, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 9;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  double e0 = 0.0;
  energy(kt, t0, i3x1, i3x2, i3x3, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
  {
    double saved = position[i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 0] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 0, d);
  }
  {
    double saved = position[i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 1] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 1, d);
  }
  {
    double saved = position[i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 2] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 2, d);
  }
}

double hessian(double kt, double t0, size_t i3x1, size_t i3x2, size_t i3x3, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double energy_added = 0.0000000000000000    ;
  {
    /* !BASE */
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
      double cse_p11_t186_g27788 = (-(x2));
      double cse_p11_t187_g27789 = (-(y2));
      double cse_p11_t188_g27790 = (-(z2));
      double vx1 = (cse_p11_t186_g27788 + x1);
      double cse_p11_t180_g27782 = (-(vx1));
      double vy1 = (cse_p11_t187_g27789 + y1);
      double cse_p11_t182_g27784 = (-(vy1));
      double vz1 = (cse_p11_t188_g27790 + z1);
      double cse_p11_t184_g27786 = (-(vz1));
      double vx2 = (cse_p11_t186_g27788 + x3);
      double cse_p11_t138_g27740 = (cse_p11_t180_g27782 + (-(vx2)));
      double cse_p11_t181_g27783 = (-(vx2));
      double vy2 = (cse_p11_t187_g27789 + y3);
      double cse_p11_t139_g27741 = (cse_p11_t182_g27784 + (-(vy2)));
      double cse_p11_t183_g27785 = (-(vy2));
      double vz2 = (cse_p11_t188_g27790 + z3);
      double cse_p11_t140_g27742 = (cse_p11_t184_g27786 + (-(vz2)));
      double cse_p11_t185_g27787 = (-(vz2));
      double cse_p271_t1_g27908 = (vx1 * vx2);
      double cse_p271_t1_g27912 = (vy1 * vy2);
      double cse_p271_t1_g27915 = (vz1 * vz2);
      double dot = (cse_p271_t1_g27908 + cse_p271_t1_g27912 + cse_p271_t1_g27915);
      double cse_p427_t1_g27941 = (vx1 * vx1);
      double cse_p427_t1_g27946 = (vy1 * vy1);
      double cse_p427_t1_g27949 = (vz1 * vz1);
      double n1_sq = (cse_p427_t1_g27941 + cse_p427_t1_g27946 + cse_p427_t1_g27949);
      double cse_p529_t1_sqrt27961 = sqrt(n1_sq);
      double cse_p529_t2_invr27962 = (1.0 / (n1_sq));
      double cse_p61_t3_invr227808 = (cse_p529_t2_invr27962 * cse_p529_t2_invr27962);
      double cse_p166_t1_g27866 = (-2.0000000000000000     * cse_p61_t3_invr227808);
      double cse_p115_t2_g27857 = (cse_p166_t1_g27866 + cse_p61_t3_invr227808);
      double cse_p271_t1_g27909 = (vx2 * vx2);
      double cse_p271_t1_g27913 = (vy2 * vy2);
      double cse_p271_t1_g27916 = (vz2 * vz2);
      double n2_sq = (cse_p271_t1_g27909 + cse_p271_t1_g27913 + cse_p271_t1_g27916);
      double cse_p529_t3_sqrt27963 = sqrt(n2_sq);
      double cse_p529_t4_invr27964 = (1.0 / (n2_sq));
      double cse_p114_t1_g27855 = (cse_p529_t1_sqrt27961 * cse_p529_t3_sqrt27963);
      double cse_p61_t6_invr227811 = (cse_p529_t4_invr27964 * cse_p529_t4_invr27964);
      double cse_p11_t115_g27717 = (cse_p271_t1_g27908 * cse_p529_t3_sqrt27963 * cse_p529_t4_invr27964);
      double cse_p219_t1_g27890 = (cse_p529_t2_invr27962 * cse_p529_t3_sqrt27963);
      double cse_p64_t1_g27848 = (cse_p219_t1_g27890 * cse_p529_t4_invr27964);
      double cse_p11_t116_g27718 = (cse_p64_t1_g27848 * vx2);
      double cse_p11_t117_g27719 = (cse_p271_t1_g27912 * cse_p529_t3_sqrt27963 * cse_p529_t4_invr27964);
      double cse_p11_t118_g27720 = (cse_p64_t1_g27848 * vy2);
      double cse_p11_t119_g27721 = (cse_p271_t1_g27915 * cse_p529_t3_sqrt27963 * cse_p529_t4_invr27964);
      double cse_p11_t120_g27722 = (cse_p64_t1_g27848 * vz2);
      double cse_p1_t2_invsqrt27596 = (1.0 / (cse_p529_t1_sqrt27961));
      double cse_p167_t1_g27875 = (2.0000000000000000     * cse_p114_t1_g27855 * cse_p61_t6_invr227811);
      double cse_p271_t1_g27910 = (vx2 * vy2);
      double cse_p11_t20_g27622 = (cse_p167_t1_g27875 * cse_p271_t1_g27910);
      double cse_p271_t1_g27911 = (vx2 * vz2);
      double cse_p11_t21_g27623 = (cse_p167_t1_g27875 * cse_p271_t1_g27911);
      double cse_p271_t1_g27914 = (vy2 * vz2);
      double cse_p11_t22_g27624 = (cse_p167_t1_g27875 * cse_p271_t1_g27914);
      double cse_p11_t1_g27800 = (cse_p114_t1_g27855 * cse_p529_t4_invr27964);
      double cse_p11_t91_g27693 = (cse_p11_t1_g27800 * vx2);
      double cse_p168_t1_g27878 = (cse_p11_t1_g27800 * cse_p529_t4_invr27964);
      double cse_p11_t92_g27694 = (cse_p168_t1_g27878 * cse_p271_t1_g27909);
      double cse_p11_t94_g27696 = (cse_p168_t1_g27878 * cse_p271_t1_g27910);
      double cse_p11_t96_g27698 = (cse_p168_t1_g27878 * cse_p271_t1_g27911);
      double cse_p11_t97_g27699 = (cse_p11_t1_g27800 * vy2);
      double cse_p11_t98_g27700 = (cse_p168_t1_g27878 * cse_p271_t1_g27913);
      double cse_p11_t100_g27702 = (cse_p168_t1_g27878 * cse_p271_t1_g27914);
      double cse_p11_t101_g27703 = (cse_p11_t1_g27800 * vz2);
      double cse_p11_t102_g27704 = (cse_p168_t1_g27878 * cse_p271_t1_g27916);
      double cse_p11_t162_g27764 = (-(cse_p11_t91_g27693));
      double cse_p11_t164_g27766 = (-(cse_p11_t94_g27696));
      double cse_p11_t166_g27768 = (-(cse_p11_t96_g27698));
      double cse_p11_t167_g27769 = (-(cse_p11_t97_g27699));
      double cse_p11_t169_g27771 = (-(cse_p11_t100_g27702));
      double cse_p11_t170_g27772 = (-(cse_p11_t101_g27703));
      double cse_p11_t177_g27779 = (-(cse_p11_t1_g27800));
      double cse_p1_t4_invsqrt27598 = (1.0 / (cse_p529_t3_sqrt27963));
      double cse_p11_t1_g27801 = (cse_p114_t1_g27855 * cse_p529_t2_invr27962);
      double cse_p530_t1_g27967 = (cse_p11_t1_g27801 * vz1);
      double cse_p426_t1_g27940 = (cse_p11_t1_g27801 * vy1);
      double cse_p270_t1_g27906 = (cse_p11_t1_g27801 * vx1);
      double cse_p61_t7_invr27812 = (1.0 / (cse_p114_t1_g27855));
      double cse_p61_t8_invr227813 = (cse_p61_t7_invr27812 * cse_p61_t7_invr27812);
      double cse_p374_t1_g27933 = (cse_p11_t1_g27800 * cse_p61_t8_invr227813);
      double cse_p218_t1_g27885 = (cse_p61_t8_invr227813 * dot);
      double cse_p322_t1_g27926 = (cse_p11_t1_g27800 * cse_p218_t1_g27885);
      double cse_p61_t9_invr327814 = (cse_p61_t7_invr27812 * cse_p61_t8_invr227813);
      double cse_p11_t1_g27603 = (cse_p61_t8_invr227813 * (cse_p11_t162_g27764 + (-(cse_p270_t1_g27906))));
      double cse_p11_t2_g27604 = (cse_p61_t8_invr227813 * (cse_p11_t167_g27769 + (-(cse_p426_t1_g27940))));
      double cse_p11_t3_g27605 = (cse_p61_t8_invr227813 * (cse_p11_t170_g27772 + (-(cse_p530_t1_g27967))));
      double cse_p11_t13_g27615 = (cse_p166_t1_g27866 * cse_p529_t3_sqrt27963 * vx1);
      double cse_p11_t14_g27616 = (cse_p166_t1_g27866 * cse_p529_t3_sqrt27963 * vy1);
      double cse_p11_t15_g27617 = (cse_p166_t1_g27866 * cse_p529_t3_sqrt27963 * vz1);
      double cse_p219_t1_g27895 = (2.0000000000000000     * cse_p61_t8_invr227813);
      double cse_p11_t16_g27618 = (cse_p219_t1_g27895 * cse_p61_t3_invr227808);
      double cse_p11_t17_g27619 = (cse_p219_t1_g27895 * cse_p61_t6_invr227811);
      double cse_p63_t1_g27834 = (2.0000000000000000     * cse_p61_t9_invr327814);
      double cse_p11_t19_g27621 = (cse_p168_t1_g27878 * cse_p63_t1_g27834);
      double cse_p11_t23_g27625 = (cse_p11_t1_g27801 * cse_p529_t2_invr27962 * cse_p63_t1_g27834);
      double cse_p11_t33_g27635 = (dot * (-(cse_p11_t1_g27603)));
      double cse_p11_t34_g27636 = (dot * (-(cse_p11_t2_g27604)));
      double cse_p11_t35_g27637 = (dot * (-(cse_p11_t3_g27605)));
      double cse_p11_t37_g27639 = (cse_p322_t1_g27926 * vx2);
      double cse_p11_t38_g27640 = (cse_p322_t1_g27926 * vy2);
      double cse_p11_t39_g27641 = (cse_p322_t1_g27926 * vz2);
      double cse_p11_t40_g27642 = (cse_p11_t1_g27801 * cse_p218_t1_g27885);
      double cse_p11_t41_g27643 = (cse_p218_t1_g27885 * cse_p270_t1_g27906);
      double cse_p11_t42_g27644 = (cse_p218_t1_g27885 * cse_p426_t1_g27940);
      double cse_p11_t43_g27645 = (cse_p218_t1_g27885 * cse_p530_t1_g27967);
      double cse_p11_t93_g27695 = (cse_p271_t1_g27910 * cse_p374_t1_g27933);
      double cse_p11_t95_g27697 = (cse_p271_t1_g27911 * cse_p374_t1_g27933);
      double cse_p11_t99_g27701 = (cse_p271_t1_g27914 * cse_p374_t1_g27933);
      double cse_p63_t1_g27847 = (cse_p219_t1_g27890 * cse_p529_t2_invr27962);
      double cse_p11_t103_g27705 = (cse_p63_t1_g27847 * vx1);
      double cse_p11_t105_g27707 = (cse_p219_t1_g27890 * cse_p427_t1_g27941);
      double cse_p11_t106_g27708 = (cse_p270_t1_g27906 * cse_p61_t8_invr227813 * vy1);
      double cse_p11_t107_g27709 = (cse_p270_t1_g27906 * cse_p61_t8_invr227813 * vz1);
      double cse_p11_t108_g27710 = (cse_p63_t1_g27847 * vy1);
      double cse_p11_t110_g27712 = (cse_p219_t1_g27890 * cse_p427_t1_g27946);
      double cse_p11_t111_g27713 = (cse_p426_t1_g27940 * cse_p61_t8_invr227813 * vz1);
      double cse_p11_t112_g27714 = (cse_p63_t1_g27847 * vz1);
      double cse_p11_t114_g27716 = (cse_p219_t1_g27890 * cse_p427_t1_g27949);
      double cse_p11_t121_g27723 = (cse_p11_t162_g27764 + (-(cse_p270_t1_g27906)));
      double cse_p63_t1_g27841 = (cse_p219_t1_g27890 * cse_p63_t1_g27834);
      double cse_p63_t1_g27819 = (cse_p11_t121_g27723 * cse_p63_t1_g27841);
      double cse_p11_t122_g27724 = (cse_p11_t167_g27769 + (-(cse_p426_t1_g27940)));
      double cse_p63_t2_g27820 = (cse_p11_t122_g27724 * cse_p63_t1_g27841);
      double cse_p11_t123_g27725 = (cse_p11_t170_g27772 + (-(cse_p530_t1_g27967)));
      double cse_p63_t3_g27821 = (cse_p11_t123_g27725 * cse_p63_t1_g27841);
      double cse_p11_t124_g27726 = (cse_p11_t103_g27705 + cse_p11_t116_g27718 + cse_p11_t13_g27615);
      double cse_p63_t5_g27823 = (cse_p11_t124_g27726 * cse_p61_t8_invr227813);
      double cse_p63_t8_g27826 = (cse_p63_t1_g27819 + cse_p63_t5_g27823);
      double cse_p11_t125_g27727 = (cse_p11_t108_g27710 + cse_p11_t118_g27720 + cse_p11_t14_g27616);
      double cse_p63_t6_g27824 = (cse_p11_t125_g27727 * cse_p61_t8_invr227813);
      double cse_p63_t9_g27827 = (cse_p63_t2_g27820 + cse_p63_t6_g27824);
      double cse_p11_t126_g27728 = (cse_p11_t112_g27714 + cse_p11_t120_g27722 + cse_p11_t15_g27617);
      double cse_p63_t7_g27825 = (cse_p11_t126_g27728 * cse_p61_t8_invr227813);
      double cse_p63_t10_g27828 = (cse_p63_t3_g27821 + cse_p63_t7_g27825);
      double cse_p11_t127_g27729 = (cse_p11_t16_g27618 + cse_p11_t23_g27625);
      double cse_p11_t128_g27730 = (cse_p11_t17_g27619 + cse_p11_t19_g27621);
      double cse_p11_t142_g27744 = (-((cse_p11_t121_g27723 * cse_p61_t8_invr227813)));
      double cse_p11_t143_g27745 = (-((cse_p11_t122_g27724 * cse_p61_t8_invr227813)));
      double cse_p11_t144_g27746 = (-((cse_p11_t123_g27725 * cse_p61_t8_invr227813)));
      double cse_p11_t154_g27756 = (-(cse_p322_t1_g27926));
      double cse_p11_t155_g27757 = (-(cse_p11_t37_g27639));
      double cse_p11_t156_g27758 = (-(cse_p11_t38_g27640));
      double cse_p11_t157_g27759 = (-(cse_p11_t39_g27641));
      double cse_p11_t158_g27760 = (-(cse_p11_t40_g27642));
      double cse_p11_t159_g27761 = (-(cse_p11_t41_g27643));
      double cse_p11_t160_g27762 = (-(cse_p11_t42_g27644));
      double cse_p11_t161_g27763 = (-(cse_p11_t43_g27645));
      double cse_p11_t163_g27765 = (-(cse_p11_t93_g27695));
      double cse_p11_t165_g27767 = (-(cse_p11_t95_g27697));
      double cse_p11_t168_g27770 = (-(cse_p11_t99_g27701));
      double cse_p11_t171_g27773 = (-(cse_p270_t1_g27906));
      double cse_p11_t172_g27774 = (-(cse_p11_t106_g27708));
      double cse_p11_t173_g27775 = (-(cse_p11_t107_g27709));
      double cse_p11_t174_g27776 = (-(cse_p426_t1_g27940));
      double cse_p11_t175_g27777 = (-(cse_p11_t111_g27713));
      double cse_p11_t176_g27778 = (-(cse_p530_t1_g27967));
      double cse_p11_t18_g27620 = (2.0000000000000000     * cse_p61_t7_invr27812);
      double cse_p11_t24_g27626 = (cse_p61_t7_invr27812 * (cse_p11_t180_g27782 + cse_p11_t181_g27783));
      double cse_p11_t25_g27627 = (cse_p61_t7_invr27812 * (cse_p11_t182_g27784 + cse_p11_t183_g27785));
      double cse_p11_t26_g27628 = (cse_p61_t7_invr27812 * (cse_p11_t184_g27786 + cse_p11_t185_g27787));
      double cse_p11_t27_g27629 = (cse_p61_t7_invr27812 * vx1);
      double cse_p11_t28_g27630 = (cse_p61_t7_invr27812 * vx2);
      double cse_p11_t29_g27631 = (cse_p61_t7_invr27812 * vy1);
      double cse_p11_t30_g27632 = (cse_p61_t7_invr27812 * vy2);
      double cse_p11_t31_g27633 = (cse_p61_t7_invr27812 * vz1);
      double cse_p11_t32_g27634 = (cse_p61_t7_invr27812 * vz2);
      double cse_p11_t129_g27731 = (cse_p11_t24_g27626 + (-((cse_p218_t1_g27885 * (cse_p11_t162_g27764 + cse_p11_t171_g27773)))));
      double cse_p11_t130_g27732 = (cse_p11_t25_g27627 + (-((cse_p218_t1_g27885 * (cse_p11_t167_g27769 + cse_p11_t174_g27776)))));
      double cse_p11_t131_g27733 = (cse_p11_t26_g27628 + (-((cse_p218_t1_g27885 * (cse_p11_t170_g27772 + cse_p11_t176_g27778)))));
      double cse_p11_t132_g27734 = (cse_p11_t155_g27757 + cse_p11_t27_g27629);
      double cse_p11_t133_g27735 = (cse_p11_t159_g27761 + cse_p11_t28_g27630);
      double cse_p11_t134_g27736 = (cse_p11_t156_g27758 + cse_p11_t29_g27631);
      double cse_p11_t135_g27737 = (cse_p11_t160_g27762 + cse_p11_t30_g27632);
      double cse_p11_t136_g27738 = (cse_p11_t157_g27759 + cse_p11_t31_g27633);
      double cse_p11_t137_g27739 = (cse_p11_t161_g27763 + cse_p11_t32_g27634);
      double cse_p11_t179_g27781 = (-(cse_p61_t7_invr27812));
      double cos_theta = (cse_p61_t7_invr27812 * dot);
      double cse_p61_t10_invsqrt27815 = (1.0 / sqrt((1.0000000000000000     + (-(((cos_theta) * (cos_theta)))))));
      double cse_p11_t4_g27606 = (cse_p11_t129_g27731 * cse_p61_t10_invsqrt27815);
      double cse_p11_t5_g27607 = (cse_p11_t130_g27732 * cse_p61_t10_invsqrt27815);
      double cse_p11_t6_g27608 = (cse_p11_t131_g27733 * cse_p61_t10_invsqrt27815);
      double cse_p11_t7_g27609 = (cse_p11_t132_g27734 * cse_p61_t10_invsqrt27815);
      double cse_p11_t8_g27610 = (cse_p11_t133_g27735 * cse_p61_t10_invsqrt27815);
      double cse_p11_t9_g27611 = (cse_p11_t134_g27736 * cse_p61_t10_invsqrt27815);
      double cse_p11_t10_g27612 = (cse_p11_t135_g27737 * cse_p61_t10_invsqrt27815);
      double cse_p11_t11_g27613 = (cse_p11_t136_g27738 * cse_p61_t10_invsqrt27815);
      double cse_p11_t12_g27614 = (cse_p11_t137_g27739 * cse_p61_t10_invsqrt27815);
      double cse_p63_t1_g27830 = (cse_p61_t10_invsqrt27815 * cse_p61_t10_invsqrt27815 * kt);
      double cse_p167_t1_g27868 = (cse_p11_t129_g27731 * cse_p63_t1_g27830);
      double cse_p11_t44_g27646 = (cse_p11_t129_g27731 * cse_p167_t1_g27868);
      double cse_p11_t45_g27647 = (cse_p11_t130_g27732 * cse_p167_t1_g27868);
      double cse_p11_t46_g27648 = (cse_p11_t131_g27733 * cse_p167_t1_g27868);
      double cse_p11_t47_g27649 = (cse_p11_t132_g27734 * cse_p167_t1_g27868);
      double cse_p11_t48_g27650 = (cse_p11_t133_g27735 * cse_p167_t1_g27868);
      double cse_p11_t49_g27651 = (cse_p11_t134_g27736 * cse_p167_t1_g27868);
      double cse_p11_t50_g27652 = (cse_p11_t135_g27737 * cse_p167_t1_g27868);
      double cse_p11_t51_g27653 = (cse_p11_t136_g27738 * cse_p167_t1_g27868);
      double cse_p11_t52_g27654 = (cse_p11_t137_g27739 * cse_p167_t1_g27868);
      double cse_p167_t1_g27872 = (cse_p11_t130_g27732 * cse_p63_t1_g27830);
      double cse_p11_t53_g27655 = (cse_p11_t130_g27732 * cse_p167_t1_g27872);
      double cse_p11_t54_g27656 = (cse_p11_t131_g27733 * cse_p167_t1_g27872);
      double cse_p11_t55_g27657 = (cse_p11_t132_g27734 * cse_p167_t1_g27872);
      double cse_p11_t56_g27658 = (cse_p11_t133_g27735 * cse_p167_t1_g27872);
      double cse_p11_t57_g27659 = (cse_p11_t134_g27736 * cse_p167_t1_g27872);
      double cse_p11_t58_g27660 = (cse_p11_t135_g27737 * cse_p167_t1_g27872);
      double cse_p11_t59_g27661 = (cse_p11_t136_g27738 * cse_p167_t1_g27872);
      double cse_p11_t60_g27662 = (cse_p11_t137_g27739 * cse_p167_t1_g27872);
      double cse_p167_t1_g27876 = (cse_p11_t131_g27733 * cse_p63_t1_g27830);
      double cse_p11_t61_g27663 = (cse_p11_t131_g27733 * cse_p167_t1_g27876);
      double cse_p11_t62_g27664 = (cse_p11_t132_g27734 * cse_p167_t1_g27876);
      double cse_p11_t63_g27665 = (cse_p11_t133_g27735 * cse_p167_t1_g27876);
      double cse_p11_t64_g27666 = (cse_p11_t134_g27736 * cse_p167_t1_g27876);
      double cse_p11_t65_g27667 = (cse_p11_t135_g27737 * cse_p167_t1_g27876);
      double cse_p11_t66_g27668 = (cse_p11_t136_g27738 * cse_p167_t1_g27876);
      double cse_p11_t67_g27669 = (cse_p11_t137_g27739 * cse_p167_t1_g27876);
      double cse_p219_t1_g27888 = (cse_p11_t132_g27734 * cse_p63_t1_g27830);
      double cse_p11_t68_g27670 = (cse_p11_t132_g27734 * cse_p219_t1_g27888);
      double cse_p11_t69_g27671 = (cse_p11_t133_g27735 * cse_p219_t1_g27888);
      double cse_p11_t70_g27672 = (cse_p11_t134_g27736 * cse_p219_t1_g27888);
      double cse_p11_t71_g27673 = (cse_p11_t135_g27737 * cse_p219_t1_g27888);
      double cse_p11_t72_g27674 = (cse_p11_t136_g27738 * cse_p219_t1_g27888);
      double cse_p11_t73_g27675 = (cse_p11_t137_g27739 * cse_p219_t1_g27888);
      double cse_p220_t1_g27898 = (cse_p11_t133_g27735 * cse_p63_t1_g27830);
      double cse_p11_t74_g27676 = (cse_p11_t133_g27735 * cse_p220_t1_g27898);
      double cse_p11_t75_g27677 = (cse_p11_t134_g27736 * cse_p220_t1_g27898);
      double cse_p11_t76_g27678 = (cse_p11_t135_g27737 * cse_p220_t1_g27898);
      double cse_p11_t77_g27679 = (cse_p11_t136_g27738 * cse_p220_t1_g27898);
      double cse_p11_t78_g27680 = (cse_p11_t137_g27739 * cse_p220_t1_g27898);
      double cse_p272_t1_g27918 = (cse_p11_t134_g27736 * cse_p63_t1_g27830);
      double cse_p11_t79_g27681 = (cse_p11_t134_g27736 * cse_p272_t1_g27918);
      double cse_p11_t80_g27682 = (cse_p11_t135_g27737 * cse_p272_t1_g27918);
      double cse_p11_t81_g27683 = (cse_p11_t136_g27738 * cse_p272_t1_g27918);
      double cse_p11_t82_g27684 = (cse_p11_t137_g27739 * cse_p272_t1_g27918);
      double cse_p428_t1_g27952 = (cse_p11_t135_g27737 * cse_p63_t1_g27830);
      double cse_p11_t83_g27685 = (cse_p11_t135_g27737 * cse_p428_t1_g27952);
      double cse_p11_t84_g27686 = (cse_p11_t136_g27738 * cse_p428_t1_g27952);
      double cse_p11_t85_g27687 = (cse_p11_t137_g27739 * cse_p428_t1_g27952);
      double cse_p11_t86_g27688 = (cse_p11_t136_g27738 * cse_p11_t136_g27738 * cse_p63_t1_g27830);
      double cse_p11_t87_g27689 = (cse_p11_t136_g27738 * cse_p11_t137_g27739 * cse_p63_t1_g27830);
      double cse_p11_t88_g27690 = (cse_p11_t137_g27739 * cse_p11_t137_g27739 * cse_p63_t1_g27830);
      double cse_p11_t141_g27743 = (1.0000000000000000     + (-(((cos_theta) * (cos_theta)))));
      double cse_p11_t145_g27747 = (-(cse_p11_t4_g27606));
      double cse_p11_t146_g27748 = (-(cse_p11_t5_g27607));
      double cse_p11_t147_g27749 = (-(cse_p11_t6_g27608));
      double cse_p11_t148_g27750 = (-(cse_p11_t7_g27609));
      double cse_p11_t149_g27751 = (-(cse_p11_t8_g27610));
      double cse_p11_t150_g27752 = (-(cse_p11_t9_g27611));
      double cse_p11_t151_g27753 = (-(cse_p11_t10_g27612));
      double cse_p11_t152_g27754 = (-(cse_p11_t11_g27613));
      double cse_p11_t153_g27755 = (-(cse_p11_t12_g27614));
      double cse_p11_t178_g27780 = (-(((cos_theta) * (cos_theta))));
      double cse_p529_t5_sqrt27965 = sqrt(cse_p11_t141_g27743);
      double cse_p529_t6_invr27966 = (1.0 / (cse_p11_t141_g27743));
      double cse_p11_t193_g27795 = ((cos_theta) * (cos_theta));
      if (fabs(cos_theta)>(1.0-VERYSMALL)) linear_angle_error();
      double cse_p1_t8_invsqrt27602 = (1.0 / (cse_p529_t5_sqrt27965));
      double theta = acos(cos_theta);
      double dtheta = (theta + (-(t0)));
      double cse_p62_t1_g27818 = (2.0000000000000000     * dtheta * kt);
      double energy = (dtheta * dtheta * kt);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p63_t1_g27833 = (cse_p529_t3_sqrt27963 * dot);
      double cse_p63_t1_g27846 = (cse_p1_t2_invsqrt27596 * cse_p61_t8_invr227813 * cse_p63_t1_g27833);
      double cse_p167_t1_g27870 = (cse_p1_t8_invsqrt27602 * cse_p62_t1_g27818);
      double g_x1 = (-((cse_p167_t1_g27870 * (cse_p11_t28_g27630 + (-((cse_p63_t1_g27846 * vx1)))))));
      KernelGradientAcc(i3x1, 0, g_x1);
      double g_y1 = (-((cse_p167_t1_g27870 * (cse_p11_t30_g27632 + (-((cse_p63_t1_g27846 * vy1)))))));
      KernelGradientAcc(i3x1, 1, g_y1);
      double g_z1 = (-((cse_p167_t1_g27870 * (cse_p11_t32_g27634 + (-((cse_p63_t1_g27846 * vz1)))))));
      KernelGradientAcc(i3x1, 2, g_z1);
      double g_x2 = (-((cse_p167_t1_g27870 * (cse_p11_t24_g27626 + (-((cse_p218_t1_g27885 * ((-((cse_p1_t2_invsqrt27596 * cse_p529_t3_sqrt27963 * vx1))) + (-((cse_p1_t4_invsqrt27598 * cse_p529_t1_sqrt27961 * vx2)))))))))));
      KernelGradientAcc(i3x2, 0, g_x2);
      double g_y2 = (-((cse_p167_t1_g27870 * (cse_p11_t25_g27627 + (-((cse_p218_t1_g27885 * ((-((cse_p1_t2_invsqrt27596 * cse_p529_t3_sqrt27963 * vy1))) + (-((cse_p1_t4_invsqrt27598 * cse_p529_t1_sqrt27961 * vy2)))))))))));
      KernelGradientAcc(i3x2, 1, g_y2);
      double g_z2 = (-((cse_p167_t1_g27870 * (cse_p11_t26_g27628 + (-((cse_p218_t1_g27885 * ((-((cse_p1_t2_invsqrt27596 * cse_p529_t3_sqrt27963 * vz1))) + (-((cse_p1_t4_invsqrt27598 * cse_p529_t1_sqrt27961 * vz2)))))))))));
      KernelGradientAcc(i3x2, 2, g_z2);
      double cse_p115_t1_g27858 = (cse_p529_t1_sqrt27961 * dot);
      double cse_p115_t1_g27859 = (cse_p115_t1_g27858 * cse_p1_t4_invsqrt27598 * cse_p61_t8_invr227813);
      double g_x3 = (-((cse_p167_t1_g27870 * (cse_p11_t27_g27629 + (-((cse_p115_t1_g27859 * vx2)))))));
      KernelGradientAcc(i3x3, 0, g_x3);
      double g_y3 = (-((cse_p167_t1_g27870 * (cse_p11_t29_g27631 + (-((cse_p115_t1_g27859 * vy2)))))));
      KernelGradientAcc(i3x3, 1, g_y3);
      double g_z3 = (-((cse_p167_t1_g27870 * (cse_p11_t31_g27633 + (-((cse_p115_t1_g27859 * vz2)))))));
      KernelGradientAcc(i3x3, 2, g_z3);
      double cse_p271_t1_g27907 = (cse_p529_t1_sqrt27961 * cse_p529_t2_invr27962);
      double cse_p11_t1_g27803 = (cse_p271_t1_g27907 * cse_p61_t8_invr227813);
      double cse_p63_t1_g27831 = (cos_theta * cse_p529_t6_invr27966);
      double cse_p63_t1_g27832 = (cse_p61_t10_invsqrt27815 * cse_p62_t1_g27818);
      double cse_p63_t1_g27836 = (cse_p529_t1_sqrt27961 * cse_p63_t1_g27833);
      double cse_p478_t1_g27960 = (cse_p11_t127_g27729 * cse_p63_t1_g27836);
      double cse_p63_t1_g27840 = (cse_p11_t1_g27803 * cse_p529_t2_invr27962 * cse_p63_t1_g27833);
      double cse_p63_t1_g27842 = (cse_p271_t1_g27907 * cse_p529_t2_invr27962 * cse_p63_t1_g27834 * cse_p63_t1_g27836);
      double cse_p219_t1_g27896 = (cse_p115_t1_g27858 * cse_p219_t1_g27895 * cse_p61_t3_invr227808);
      double cse_p220_t1_g27899 = (cse_p11_t133_g27735 * cse_p63_t1_g27831);
      double cse_p427_t1_g27950 = (-2.0000000000000000     * cse_p11_t1_g27803);
      double h_x1_x1 = (cse_p11_t74_g27676 + cse_p11_t74_g27676 + (-((cse_p63_t1_g27832 * (cse_p11_t158_g27760 + (cse_p11_t133_g27735 * cse_p220_t1_g27899) + (cse_p529_t3_sqrt27963 * vx1 * ((cse_p219_t1_g27896 * vx1) + (cse_p427_t1_g27950 * vx2) + (cse_p63_t1_g27842 * vx1))) + (-((cse_p427_t1_g27941 * cse_p63_t1_g27840))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 0, h_x1_x1);
      double cse_p63_t1_g27837 = (cse_p11_t1_g27803 * cse_p529_t3_sqrt27963);
      double cse_p219_t1_g27892 = (cse_p63_t1_g27837 * vx1);
      double cse_p219_t1_g27893 = (cse_p63_t1_g27837 * vy1);
      double h_x1_y1 = (cse_p11_t76_g27678 + cse_p11_t76_g27678 + (-((cse_p63_t1_g27832 * ((cse_p11_t135_g27737 * cse_p220_t1_g27899) + (cse_p478_t1_g27960 * vx1 * vy1) + (-((cse_p219_t1_g27892 * vy2))) + (-((cse_p219_t1_g27893 * vx2))) + (-((cse_p63_t1_g27840 * vx1 * vy1))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 1, h_x1_y1);
      double cse_p219_t1_g27894 = (cse_p63_t1_g27837 * vz1);
      double h_x1_z1 = (cse_p11_t78_g27680 + cse_p11_t78_g27680 + (-((cse_p63_t1_g27832 * ((cse_p11_t137_g27739 * cse_p220_t1_g27899) + (cse_p478_t1_g27960 * vx1 * vz1) + (-((cse_p219_t1_g27892 * vz2))) + (-((cse_p219_t1_g27894 * vx2))) + (-((cse_p63_t1_g27840 * vx1 * vz1))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 2, h_x1_z1);
      double cse_p167_t1_g27869 = (cse_p11_t129_g27731 * cse_p63_t1_g27831);
      double cse_p220_t1_g27897 = (cse_p11_t121_g27723 * cse_p63_t1_g27834);
      double h_x1_x2 = (cse_p11_t48_g27650 + cse_p11_t48_g27650 + (-((cse_p63_t1_g27832 * (cse_p11_t179_g27781 + (cse_p11_t133_g27735 * cse_p167_t1_g27869) + (cse_p11_t142_g27744 * vx2) + (cse_p11_t1_g27803 * dot * (cse_p11_t105_g27707 + cse_p11_t115_g27717 + cse_p529_t3_sqrt27963)) + (cse_p63_t1_g27836 * vx1 * ((cse_p166_t1_g27866 * cse_p61_t8_invr227813 * vx1) + (cse_p220_t1_g27897 * cse_p529_t2_invr27962))) + (-((cse_p11_t138_g27740 * cse_p219_t1_g27892))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 0, h_x1_x2);
      double cse_p167_t1_g27873 = (cse_p11_t130_g27732 * cse_p63_t1_g27831);
      double h_x1_y2 = (cse_p11_t56_g27658 + cse_p11_t56_g27658 + (-((cse_p63_t1_g27832 * ((cse_p115_t1_g27858 * cse_p63_t9_g27827 * vx1) + (cse_p11_t133_g27735 * cse_p167_t1_g27873) + (cse_p11_t143_g27745 * vx2) + (-((cse_p11_t139_g27741 * cse_p219_t1_g27892))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 1, h_x1_y2);
      double cse_p168_t1_g27877 = (cse_p11_t131_g27733 * cse_p63_t1_g27831);
      double h_x1_z2 = (cse_p11_t63_g27665 + cse_p11_t63_g27665 + (-((cse_p63_t1_g27832 * ((cse_p115_t1_g27858 * cse_p63_t10_g27828 * vx1) + (cse_p11_t133_g27735 * cse_p168_t1_g27877) + (cse_p11_t144_g27746 * vx2) + (-((cse_p11_t140_g27742 * cse_p219_t1_g27892))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 2, h_x1_z2);
      double cse_p11_t1_g27804 = (cse_p11_t1_g27800 * cse_p529_t2_invr27962);
      double cse_p63_t1_g27838 = (cse_p11_t1_g27803 * cse_p529_t4_invr27964 * cse_p63_t1_g27833);
      double cse_p63_t1_g27839 = (cse_p11_t1_g27804 * cse_p63_t1_g27834 * cse_p63_t1_g27836);
      double cse_p167_t1_g27867 = (cse_p374_t1_g27933 * vx2);
      double cse_p219_t1_g27889 = (cse_p11_t132_g27734 * cse_p63_t1_g27831);
      double h_x1_x3 = (cse_p11_t69_g27671 + cse_p11_t69_g27671 + (-((cse_p63_t1_g27832 * (cse_p61_t7_invr27812 + (cse_p11_t133_g27735 * cse_p219_t1_g27889) + (cse_p271_t1_g27908 * cse_p63_t1_g27839) + (-((cse_p167_t1_g27867 * vx2))) + (-((cse_p219_t1_g27892 * vx1))) + (-((cse_p271_t1_g27908 * cse_p63_t1_g27838))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x3, 0, h_x1_x3);
      double cse_p427_t1_g27942 = (vx1 * vy2);
      double h_x1_y3 = (cse_p11_t75_g27677 + cse_p11_t75_g27677 + (-((cse_p63_t1_g27832 * (cse_p11_t163_g27765 + cse_p11_t172_g27774 + (cse_p11_t134_g27736 * cse_p220_t1_g27899) + (cse_p427_t1_g27942 * cse_p63_t1_g27839) + (-((cse_p427_t1_g27942 * cse_p63_t1_g27838))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x3, 1, h_x1_y3);
      double cse_p427_t1_g27943 = (vx1 * vz2);
      double h_x1_z3 = (cse_p11_t77_g27679 + cse_p11_t77_g27679 + (-((cse_p63_t1_g27832 * (cse_p11_t165_g27767 + cse_p11_t173_g27775 + (cse_p11_t136_g27738 * cse_p220_t1_g27899) + (cse_p427_t1_g27943 * cse_p63_t1_g27839) + (-((cse_p427_t1_g27943 * cse_p63_t1_g27838))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x3, 2, h_x1_z3);
      double cse_p428_t1_g27953 = (cse_p11_t135_g27737 * cse_p63_t1_g27831);
      double h_y1_y1 = (cse_p11_t83_g27685 + cse_p11_t83_g27685 + (-((cse_p63_t1_g27832 * (cse_p11_t158_g27760 + (cse_p11_t135_g27737 * cse_p428_t1_g27953) + (cse_p529_t3_sqrt27963 * vy1 * ((cse_p219_t1_g27896 * vy1) + (cse_p427_t1_g27950 * vy2) + (cse_p63_t1_g27842 * vy1))) + (-((cse_p427_t1_g27946 * cse_p63_t1_g27840))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 1, h_y1_y1);
      double h_y1_z1 = (cse_p11_t85_g27687 + cse_p11_t85_g27687 + (-((cse_p63_t1_g27832 * ((cse_p11_t137_g27739 * cse_p428_t1_g27953) + (cse_p478_t1_g27960 * vy1 * vz1) + (-((cse_p219_t1_g27893 * vz2))) + (-((cse_p219_t1_g27894 * vy2))) + (-((cse_p63_t1_g27840 * vy1 * vz1))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 2, h_y1_z1);
      double h_y1_x2 = (cse_p11_t50_g27652 + cse_p11_t50_g27652 + (-((cse_p63_t1_g27832 * ((cse_p115_t1_g27858 * cse_p63_t8_g27826 * vy1) + (cse_p11_t135_g27737 * cse_p167_t1_g27869) + (cse_p11_t142_g27744 * vy2) + (-((cse_p11_t138_g27740 * cse_p219_t1_g27893))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 0, h_y1_x2);
      double cse_p272_t1_g27917 = (cse_p11_t122_g27724 * cse_p63_t1_g27834);
      double h_y1_y2 = (cse_p11_t58_g27660 + cse_p11_t58_g27660 + (-((cse_p63_t1_g27832 * (cse_p11_t179_g27781 + (cse_p11_t135_g27737 * cse_p167_t1_g27873) + (cse_p11_t143_g27745 * vy2) + (cse_p11_t1_g27803 * dot * (cse_p11_t110_g27712 + cse_p11_t117_g27719 + cse_p529_t3_sqrt27963)) + (cse_p63_t1_g27836 * vy1 * ((cse_p166_t1_g27866 * cse_p61_t8_invr227813 * vy1) + (cse_p272_t1_g27917 * cse_p529_t2_invr27962))) + (-((cse_p11_t139_g27741 * cse_p219_t1_g27893))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 1, h_y1_y2);
      double h_y1_z2 = (cse_p11_t65_g27667 + cse_p11_t65_g27667 + (-((cse_p63_t1_g27832 * ((cse_p115_t1_g27858 * cse_p63_t10_g27828 * vy1) + (cse_p11_t135_g27737 * cse_p168_t1_g27877) + (cse_p11_t144_g27746 * vy2) + (-((cse_p11_t140_g27742 * cse_p219_t1_g27893))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 2, h_y1_z2);
      double cse_p427_t1_g27944 = (vx2 * vy1);
      double h_y1_x3 = (cse_p11_t71_g27673 + cse_p11_t71_g27673 + (-((cse_p63_t1_g27832 * (cse_p11_t163_g27765 + cse_p11_t172_g27774 + (cse_p11_t135_g27737 * cse_p219_t1_g27889) + (cse_p427_t1_g27944 * cse_p63_t1_g27839) + (-((cse_p427_t1_g27944 * cse_p63_t1_g27838))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x3, 0, h_y1_x3);
      double cse_p167_t1_g27871 = (cse_p374_t1_g27933 * vy2);
      double cse_p272_t1_g27919 = (cse_p11_t134_g27736 * cse_p63_t1_g27831);
      double h_y1_y3 = (cse_p11_t80_g27682 + cse_p11_t80_g27682 + (-((cse_p63_t1_g27832 * (cse_p61_t7_invr27812 + (cse_p11_t135_g27737 * cse_p272_t1_g27919) + (cse_p271_t1_g27912 * cse_p63_t1_g27839) + (-((cse_p167_t1_g27871 * vy2))) + (-((cse_p219_t1_g27893 * vy1))) + (-((cse_p271_t1_g27912 * cse_p63_t1_g27838))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x3, 1, h_y1_y3);
      double cse_p427_t1_g27947 = (vy1 * vz2);
      double h_y1_z3 = (cse_p11_t84_g27686 + cse_p11_t84_g27686 + (-((cse_p63_t1_g27832 * (cse_p11_t168_g27770 + cse_p11_t175_g27777 + (cse_p11_t136_g27738 * cse_p428_t1_g27953) + (cse_p427_t1_g27947 * cse_p63_t1_g27839) + (-((cse_p427_t1_g27947 * cse_p63_t1_g27838))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x3, 2, h_y1_z3);
      double h_z1_z1 = (cse_p11_t88_g27690 + cse_p11_t88_g27690 + (-((cse_p63_t1_g27832 * (cse_p11_t158_g27760 + (cse_p11_t137_g27739 * cse_p11_t137_g27739 * cse_p63_t1_g27831) + (cse_p529_t3_sqrt27963 * vz1 * ((cse_p219_t1_g27896 * vz1) + (cse_p427_t1_g27950 * vz2) + (cse_p63_t1_g27842 * vz1))) + (-((cse_p427_t1_g27949 * cse_p63_t1_g27840))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x1, 2, h_z1_z1);
      double h_z1_x2 = (cse_p11_t52_g27654 + cse_p11_t52_g27654 + (-((cse_p63_t1_g27832 * ((cse_p115_t1_g27858 * cse_p63_t8_g27826 * vz1) + (cse_p11_t137_g27739 * cse_p167_t1_g27869) + (cse_p11_t142_g27744 * vz2) + (-((cse_p11_t138_g27740 * cse_p219_t1_g27894))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 0, h_z1_x2);
      double h_z1_y2 = (cse_p11_t60_g27662 + cse_p11_t60_g27662 + (-((cse_p63_t1_g27832 * ((cse_p115_t1_g27858 * cse_p63_t9_g27827 * vz1) + (cse_p11_t137_g27739 * cse_p167_t1_g27873) + (cse_p11_t143_g27745 * vz2) + (-((cse_p11_t139_g27741 * cse_p219_t1_g27894))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 1, h_z1_y2);
      double cse_p428_t1_g27951 = (cse_p11_t123_g27725 * cse_p63_t1_g27834);
      double h_z1_z2 = (cse_p11_t67_g27669 + cse_p11_t67_g27669 + (-((cse_p63_t1_g27832 * (cse_p11_t179_g27781 + (cse_p11_t137_g27739 * cse_p168_t1_g27877) + (cse_p11_t144_g27746 * vz2) + (cse_p11_t1_g27803 * dot * (cse_p11_t114_g27716 + cse_p11_t119_g27721 + cse_p529_t3_sqrt27963)) + (cse_p63_t1_g27836 * vz1 * ((cse_p166_t1_g27866 * cse_p61_t8_invr227813 * vz1) + (cse_p428_t1_g27951 * cse_p529_t2_invr27962))) + (-((cse_p11_t140_g27742 * cse_p219_t1_g27894))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 2, h_z1_z2);
      double cse_p427_t1_g27945 = (vx2 * vz1);
      double h_z1_x3 = (cse_p11_t73_g27675 + cse_p11_t73_g27675 + (-((cse_p63_t1_g27832 * (cse_p11_t165_g27767 + cse_p11_t173_g27775 + (cse_p11_t137_g27739 * cse_p219_t1_g27889) + (cse_p427_t1_g27945 * cse_p63_t1_g27839) + (-((cse_p427_t1_g27945 * cse_p63_t1_g27838))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x3, 0, h_z1_x3);
      double cse_p427_t1_g27948 = (vy2 * vz1);
      double h_z1_y3 = (cse_p11_t82_g27684 + cse_p11_t82_g27684 + (-((cse_p63_t1_g27832 * (cse_p11_t168_g27770 + cse_p11_t175_g27777 + (cse_p11_t137_g27739 * cse_p272_t1_g27919) + (cse_p427_t1_g27948 * cse_p63_t1_g27839) + (-((cse_p427_t1_g27948 * cse_p63_t1_g27838))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x3, 1, h_z1_y3);
      double cse_p167_t1_g27874 = (cse_p374_t1_g27933 * vz2);
      double h_z1_z3 = (cse_p11_t87_g27689 + cse_p11_t87_g27689 + (-((cse_p63_t1_g27832 * (cse_p61_t7_invr27812 + (cse_p11_t136_g27738 * cse_p11_t137_g27739 * cse_p63_t1_g27831) + (cse_p271_t1_g27915 * cse_p63_t1_g27839) + (-((cse_p167_t1_g27874 * vz2))) + (-((cse_p219_t1_g27894 * vz1))) + (-((cse_p271_t1_g27915 * cse_p63_t1_g27838))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x3, 2, h_z1_z3);
      double cse_p219_t1_g27886 = (-2.0000000000000000     * cse_p61_t6_invr227811);
      double h_x2_x2 = (cse_p11_t44_g27646 + cse_p11_t44_g27646 + (-((cse_p63_t1_g27832 * (cse_p11_t18_g27620 + (2.0000000000000000     * cse_p11_t138_g27740 * cse_p11_t142_g27744) + (cse_p11_t129_g27731 * cse_p167_t1_g27869) + (dot * ((cse_p11_t121_g27723 * cse_p220_t1_g27897) + (-((cse_p61_t8_invr227813 * (cse_p11_t1_g27800 + cse_p11_t92_g27694 + (cse_p114_t1_g27855 * ((cse_p166_t1_g27866 * cse_p427_t1_g27941) + (cse_p219_t1_g27886 * cse_p271_t1_g27909))) + (cse_p271_t1_g27907 * (cse_p11_t105_g27707 + cse_p11_t115_g27717 + cse_p11_t115_g27717 + cse_p529_t3_sqrt27963)))))))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 0, h_x2_x2);
      double h_x2_y2 = (cse_p11_t45_g27647 + cse_p11_t45_g27647 + (-((cse_p63_t1_g27832 * ((cse_p11_t130_g27732 * cse_p167_t1_g27869) + (cse_p11_t138_g27740 * cse_p11_t143_g27745) + (cse_p11_t139_g27741 * cse_p11_t142_g27744) + (dot * ((cse_p11_t122_g27724 * cse_p220_t1_g27897) + (-((cse_p61_t8_invr227813 * (cse_p11_t94_g27696 + (cse_p114_t1_g27855 * ((cse_p115_t2_g27857 * vx1 * vy1) + (cse_p219_t1_g27886 * cse_p271_t1_g27910))) + (cse_p11_t1_g27804 * (cse_p427_t1_g27942 + cse_p427_t1_g27944)))))))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 1, h_x2_y2);
      double h_x2_z2 = (cse_p11_t46_g27648 + cse_p11_t46_g27648 + (-((cse_p63_t1_g27832 * ((cse_p11_t131_g27733 * cse_p167_t1_g27869) + (cse_p11_t138_g27740 * cse_p11_t144_g27746) + (cse_p11_t140_g27742 * cse_p11_t142_g27744) + (dot * ((cse_p11_t123_g27725 * cse_p220_t1_g27897) + (-((cse_p61_t8_invr227813 * (cse_p11_t96_g27698 + (cse_p114_t1_g27855 * ((cse_p115_t2_g27857 * vx1 * vz1) + (cse_p219_t1_g27886 * cse_p271_t1_g27911))) + (cse_p11_t1_g27804 * (cse_p427_t1_g27943 + cse_p427_t1_g27945)))))))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 2, h_x2_z2);
      double cse_p63_t1_g27843 = (cse_p11_t1_g27800 * cse_p220_t1_g27897);
      double h_x2_x3 = (cse_p11_t47_g27649 + cse_p11_t47_g27649 + (-((cse_p63_t1_g27832 * (cse_p11_t179_g27781 + (cse_p11_t132_g27734 * cse_p167_t1_g27869) + (cse_p11_t142_g27744 * vx1) + (dot * ((cse_p63_t1_g27843 * vx2) + (-((cse_p61_t8_invr227813 * (cse_p11_t177_g27779 + (cse_p167_t1_g27875 * cse_p271_t1_g27909) + (-((cse_p11_t1_g27804 * cse_p271_t1_g27908))) + (-(cse_p11_t92_g27694)))))))) + (-((cse_p11_t138_g27740 * cse_p167_t1_g27867))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x3, 0, h_x2_x3);
      double h_x2_y3 = (cse_p11_t49_g27651 + cse_p11_t49_g27651 + (-((cse_p63_t1_g27832 * ((cse_p11_t134_g27736 * cse_p167_t1_g27869) + (cse_p11_t142_g27744 * vy1) + (dot * ((cse_p63_t1_g27843 * vy2) + (-((cse_p61_t8_invr227813 * (cse_p11_t164_g27766 + cse_p11_t20_g27622 + (-((cse_p11_t1_g27804 * cse_p427_t1_g27942))))))))) + (-((cse_p11_t138_g27740 * cse_p167_t1_g27871))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x3, 1, h_x2_y3);
      double h_x2_z3 = (cse_p11_t51_g27653 + cse_p11_t51_g27653 + (-((cse_p63_t1_g27832 * ((cse_p11_t136_g27738 * cse_p167_t1_g27869) + (cse_p11_t142_g27744 * vz1) + (dot * ((cse_p63_t1_g27843 * vz2) + (-((cse_p61_t8_invr227813 * (cse_p11_t166_g27768 + cse_p11_t21_g27623 + (-((cse_p11_t1_g27804 * cse_p427_t1_g27943))))))))) + (-((cse_p11_t138_g27740 * cse_p167_t1_g27874))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x3, 2, h_x2_z3);
      double h_y2_y2 = (cse_p11_t53_g27655 + cse_p11_t53_g27655 + (-((cse_p63_t1_g27832 * (cse_p11_t18_g27620 + (2.0000000000000000     * cse_p11_t139_g27741 * cse_p11_t143_g27745) + (cse_p11_t130_g27732 * cse_p167_t1_g27873) + (dot * ((cse_p11_t122_g27724 * cse_p272_t1_g27917) + (-((cse_p61_t8_invr227813 * (cse_p11_t1_g27800 + cse_p11_t98_g27700 + (cse_p114_t1_g27855 * ((cse_p166_t1_g27866 * cse_p427_t1_g27946) + (cse_p219_t1_g27886 * cse_p271_t1_g27913))) + (cse_p271_t1_g27907 * (cse_p11_t110_g27712 + cse_p11_t117_g27719 + cse_p11_t117_g27719 + cse_p529_t3_sqrt27963)))))))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 1, h_y2_y2);
      double h_y2_z2 = (cse_p11_t54_g27656 + cse_p11_t54_g27656 + (-((cse_p63_t1_g27832 * ((cse_p11_t131_g27733 * cse_p167_t1_g27873) + (cse_p11_t139_g27741 * cse_p11_t144_g27746) + (cse_p11_t140_g27742 * cse_p11_t143_g27745) + (dot * ((cse_p11_t123_g27725 * cse_p272_t1_g27917) + (-((cse_p61_t8_invr227813 * (cse_p11_t100_g27702 + (cse_p114_t1_g27855 * ((cse_p115_t2_g27857 * vy1 * vz1) + (cse_p219_t1_g27886 * cse_p271_t1_g27914))) + (cse_p11_t1_g27804 * (cse_p427_t1_g27947 + cse_p427_t1_g27948)))))))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 2, h_y2_z2);
      double cse_p63_t1_g27844 = (cse_p11_t1_g27800 * cse_p272_t1_g27917);
      double h_y2_x3 = (cse_p11_t55_g27657 + cse_p11_t55_g27657 + (-((cse_p63_t1_g27832 * ((cse_p11_t132_g27734 * cse_p167_t1_g27873) + (cse_p11_t143_g27745 * vx1) + (dot * ((cse_p63_t1_g27844 * vx2) + (-((cse_p61_t8_invr227813 * (cse_p11_t164_g27766 + cse_p11_t20_g27622 + (-((cse_p11_t1_g27804 * cse_p427_t1_g27944))))))))) + (-((cse_p11_t139_g27741 * cse_p167_t1_g27867))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x3, 0, h_y2_x3);
      double h_y2_y3 = (cse_p11_t57_g27659 + cse_p11_t57_g27659 + (-((cse_p63_t1_g27832 * (cse_p11_t179_g27781 + (cse_p11_t134_g27736 * cse_p167_t1_g27873) + (cse_p11_t143_g27745 * vy1) + (dot * ((cse_p63_t1_g27844 * vy2) + (-((cse_p61_t8_invr227813 * (cse_p11_t177_g27779 + (cse_p167_t1_g27875 * cse_p271_t1_g27913) + (-((cse_p11_t1_g27804 * cse_p271_t1_g27912))) + (-(cse_p11_t98_g27700)))))))) + (-((cse_p11_t139_g27741 * cse_p167_t1_g27871))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x3, 1, h_y2_y3);
      double h_y2_z3 = (cse_p11_t59_g27661 + cse_p11_t59_g27661 + (-((cse_p63_t1_g27832 * ((cse_p11_t136_g27738 * cse_p167_t1_g27873) + (cse_p11_t143_g27745 * vz1) + (dot * ((cse_p63_t1_g27844 * vz2) + (-((cse_p61_t8_invr227813 * (cse_p11_t169_g27771 + cse_p11_t22_g27624 + (-((cse_p11_t1_g27804 * cse_p427_t1_g27947))))))))) + (-((cse_p11_t139_g27741 * cse_p167_t1_g27874))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x3, 2, h_y2_z3);
      double h_z2_z2 = (cse_p11_t61_g27663 + cse_p11_t61_g27663 + (-((cse_p63_t1_g27832 * (cse_p11_t18_g27620 + (2.0000000000000000     * cse_p11_t140_g27742 * cse_p11_t144_g27746) + (cse_p11_t131_g27733 * cse_p168_t1_g27877) + (dot * ((cse_p11_t123_g27725 * cse_p428_t1_g27951) + (-((cse_p61_t8_invr227813 * (cse_p11_t102_g27704 + cse_p11_t1_g27800 + (cse_p114_t1_g27855 * ((cse_p166_t1_g27866 * cse_p427_t1_g27949) + (cse_p219_t1_g27886 * cse_p271_t1_g27916))) + (cse_p271_t1_g27907 * (cse_p11_t114_g27716 + cse_p11_t119_g27721 + cse_p11_t119_g27721 + cse_p529_t3_sqrt27963)))))))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x2, 2, h_z2_z2);
      double cse_p63_t1_g27845 = (cse_p11_t1_g27800 * cse_p428_t1_g27951);
      double h_z2_x3 = (cse_p11_t62_g27664 + cse_p11_t62_g27664 + (-((cse_p63_t1_g27832 * ((cse_p11_t132_g27734 * cse_p168_t1_g27877) + (cse_p11_t144_g27746 * vx1) + (dot * ((cse_p63_t1_g27845 * vx2) + (-((cse_p61_t8_invr227813 * (cse_p11_t166_g27768 + cse_p11_t21_g27623 + (-((cse_p11_t1_g27804 * cse_p427_t1_g27945))))))))) + (-((cse_p11_t140_g27742 * cse_p167_t1_g27867))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x3, 0, h_z2_x3);
      double h_z2_y3 = (cse_p11_t64_g27666 + cse_p11_t64_g27666 + (-((cse_p63_t1_g27832 * ((cse_p11_t134_g27736 * cse_p168_t1_g27877) + (cse_p11_t144_g27746 * vy1) + (dot * ((cse_p63_t1_g27845 * vy2) + (-((cse_p61_t8_invr227813 * (cse_p11_t169_g27771 + cse_p11_t22_g27624 + (-((cse_p11_t1_g27804 * cse_p427_t1_g27948))))))))) + (-((cse_p11_t140_g27742 * cse_p167_t1_g27871))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x3, 1, h_z2_y3);
      double h_z2_z3 = (cse_p11_t66_g27668 + cse_p11_t66_g27668 + (-((cse_p63_t1_g27832 * (cse_p11_t179_g27781 + (cse_p11_t136_g27738 * cse_p168_t1_g27877) + (cse_p11_t144_g27746 * vz1) + (dot * ((cse_p63_t1_g27845 * vz2) + (-((cse_p61_t8_invr227813 * (cse_p11_t177_g27779 + (cse_p167_t1_g27875 * cse_p271_t1_g27916) + (-((cse_p11_t1_g27804 * cse_p271_t1_g27915))) + (-(cse_p11_t102_g27704)))))))) + (-((cse_p11_t140_g27742 * cse_p167_t1_g27874))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x3, 2, h_z2_z3);
      double cse_p219_t1_g27887 = (cse_p11_t128_g27730 * cse_p63_t1_g27836);
      double cse_p219_t1_g27891 = (cse_p529_t4_invr27964 * dot);
      double h_x3_x3 = (cse_p11_t68_g27670 + cse_p11_t68_g27670 + (-((cse_p63_t1_g27832 * (cse_p11_t154_g27756 + (-2.0000000000000000     * cse_p167_t1_g27867 * vx1) + (cse_p11_t132_g27734 * cse_p219_t1_g27889) + (cse_p219_t1_g27887 * cse_p271_t1_g27909) + (-((cse_p167_t1_g27867 * cse_p219_t1_g27891 * vx2))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 0, i3x3, 0, h_x3_x3);
      double h_x3_y3 = (cse_p11_t70_g27672 + cse_p11_t70_g27672 + (-((cse_p63_t1_g27832 * ((cse_p11_t134_g27736 * cse_p219_t1_g27889) + (cse_p219_t1_g27887 * cse_p271_t1_g27910) + (-((cse_p167_t1_g27867 * cse_p219_t1_g27891 * vy2))) + (-((cse_p167_t1_g27867 * vy1))) + (-((cse_p167_t1_g27871 * vx1))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 0, i3x3, 1, h_x3_y3);
      double h_x3_z3 = (cse_p11_t72_g27674 + cse_p11_t72_g27674 + (-((cse_p63_t1_g27832 * ((cse_p11_t136_g27738 * cse_p219_t1_g27889) + (cse_p219_t1_g27887 * cse_p271_t1_g27911) + (-((cse_p167_t1_g27867 * cse_p219_t1_g27891 * vz2))) + (-((cse_p167_t1_g27867 * vz1))) + (-((cse_p167_t1_g27874 * vx1))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 0, i3x3, 2, h_x3_z3);
      double h_y3_y3 = (cse_p11_t79_g27681 + cse_p11_t79_g27681 + (-((cse_p63_t1_g27832 * (cse_p11_t154_g27756 + (-2.0000000000000000     * cse_p167_t1_g27871 * vy1) + (cse_p11_t134_g27736 * cse_p272_t1_g27919) + (cse_p219_t1_g27887 * cse_p271_t1_g27913) + (-((cse_p167_t1_g27871 * cse_p219_t1_g27891 * vy2))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 1, i3x3, 1, h_y3_y3);
      double h_y3_z3 = (cse_p11_t81_g27683 + cse_p11_t81_g27683 + (-((cse_p63_t1_g27832 * ((cse_p11_t136_g27738 * cse_p272_t1_g27919) + (cse_p219_t1_g27887 * cse_p271_t1_g27914) + (-((cse_p167_t1_g27871 * cse_p219_t1_g27891 * vz2))) + (-((cse_p167_t1_g27871 * vz1))) + (-((cse_p167_t1_g27874 * vy1))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 1, i3x3, 2, h_y3_z3);
      double h_z3_z3 = (cse_p11_t86_g27688 + cse_p11_t86_g27688 + (-((cse_p63_t1_g27832 * (cse_p11_t154_g27756 + (-2.0000000000000000     * cse_p167_t1_g27874 * vz1) + (cse_p11_t136_g27738 * cse_p11_t136_g27738 * cse_p63_t1_g27831) + (cse_p219_t1_g27887 * cse_p271_t1_g27916) + (-((cse_p167_t1_g27874 * cse_p219_t1_g27891 * vz2))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 2, i3x3, 2, h_z3_z3);
    }
  }
  return energy_added;
}
void hessian_fd(double kt, double t0, size_t i3x1, size_t i3x2, size_t i3x3, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 9;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  const double invh2 = 1.0/((h*h));
  double e0 = 0.0;
  energy(kt, t0, i3x1, i3x2, i3x3, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
  {
    double saved = position[i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 0] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 0, d);
  }
  {
    double saved = position[i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 1] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 1, d);
  }
  {
    double saved = position[i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 2] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 2, d);
  }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 0, i3x1, 0, hval);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 1, i3x1, 1, hval);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 2, i3x1, 2, hval);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 0, i3x2, 0, hval);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 1, i3x2, 1, hval);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x2, 2, hval);
  }
  {
    double saved = position[i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 0] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 0, i3x3, 0, hval);
  }
  {
    double saved = position[i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 1] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x3, 1, hval);
  }
  {
    double saved = position[i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 2] = saved + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x3, 2, hval);
  }
  {
    double saved_i = position[i3x1 + 1];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x1 + 1] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
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
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
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
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i; position[i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 2, i3x1, 1, hval);
  }
  {
    double saved_i = position[i3x2 + 0];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 0] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 0, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x2 + 0];
    double saved_j = position[i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 0] = saved_i + h; position[i3x1 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i; position[i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 0, i3x1, 1, hval);
  }
  {
    double saved_i = position[i3x2 + 0];
    double saved_j = position[i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 0] = saved_i + h; position[i3x1 + 2] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i; position[i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 0, i3x1, 2, hval);
  }
  {
    double saved_i = position[i3x2 + 1];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 1] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 1, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x2 + 1];
    double saved_j = position[i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 1] = saved_i + h; position[i3x1 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i; position[i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 1, i3x1, 1, hval);
  }
  {
    double saved_i = position[i3x2 + 1];
    double saved_j = position[i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 1] = saved_i + h; position[i3x1 + 2] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i; position[i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 1, i3x1, 2, hval);
  }
  {
    double saved_i = position[i3x2 + 1];
    double saved_j = position[i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 1] = saved_i + h; position[i3x2 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i; position[i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 1, i3x2, 0, hval);
  }
  {
    double saved_i = position[i3x2 + 2];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 2] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x2 + 2];
    double saved_j = position[i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 2] = saved_i + h; position[i3x1 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i; position[i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x1, 1, hval);
  }
  {
    double saved_i = position[i3x2 + 2];
    double saved_j = position[i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 2] = saved_i + h; position[i3x1 + 2] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i; position[i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x1, 2, hval);
  }
  {
    double saved_i = position[i3x2 + 2];
    double saved_j = position[i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 2] = saved_i + h; position[i3x2 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i; position[i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x2, 0, hval);
  }
  {
    double saved_i = position[i3x2 + 2];
    double saved_j = position[i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 2] = saved_i + h; position[i3x2 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i; position[i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x2, 1, hval);
  }
  {
    double saved_i = position[i3x3 + 0];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 0] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 0, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x3 + 0];
    double saved_j = position[i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 0] = saved_i + h; position[i3x1 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i; position[i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 0, i3x1, 1, hval);
  }
  {
    double saved_i = position[i3x3 + 0];
    double saved_j = position[i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 0] = saved_i + h; position[i3x1 + 2] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i; position[i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 0, i3x1, 2, hval);
  }
  {
    double saved_i = position[i3x3 + 0];
    double saved_j = position[i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 0] = saved_i + h; position[i3x2 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i; position[i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 0, i3x2, 0, hval);
  }
  {
    double saved_i = position[i3x3 + 0];
    double saved_j = position[i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 0] = saved_i + h; position[i3x2 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i; position[i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 0, i3x2, 1, hval);
  }
  {
    double saved_i = position[i3x3 + 0];
    double saved_j = position[i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 0] = saved_i + h; position[i3x2 + 2] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i; position[i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 0, i3x2, 2, hval);
  }
  {
    double saved_i = position[i3x3 + 1];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 1] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x3 + 1];
    double saved_j = position[i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 1] = saved_i + h; position[i3x1 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i; position[i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x1, 1, hval);
  }
  {
    double saved_i = position[i3x3 + 1];
    double saved_j = position[i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 1] = saved_i + h; position[i3x1 + 2] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i; position[i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x1, 2, hval);
  }
  {
    double saved_i = position[i3x3 + 1];
    double saved_j = position[i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 1] = saved_i + h; position[i3x2 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i; position[i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x2, 0, hval);
  }
  {
    double saved_i = position[i3x3 + 1];
    double saved_j = position[i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 1] = saved_i + h; position[i3x2 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i; position[i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x2, 1, hval);
  }
  {
    double saved_i = position[i3x3 + 1];
    double saved_j = position[i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 1] = saved_i + h; position[i3x2 + 2] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i; position[i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x2, 2, hval);
  }
  {
    double saved_i = position[i3x3 + 1];
    double saved_j = position[i3x3 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 1] = saved_i + h; position[i3x3 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i; position[i3x3 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x3, 0, hval);
  }
  {
    double saved_i = position[i3x3 + 2];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 2] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x3 + 2];
    double saved_j = position[i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 2] = saved_i + h; position[i3x1 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i; position[i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x1, 1, hval);
  }
  {
    double saved_i = position[i3x3 + 2];
    double saved_j = position[i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 2] = saved_i + h; position[i3x1 + 2] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i; position[i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x1, 2, hval);
  }
  {
    double saved_i = position[i3x3 + 2];
    double saved_j = position[i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 2] = saved_i + h; position[i3x2 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i; position[i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x2, 0, hval);
  }
  {
    double saved_i = position[i3x3 + 2];
    double saved_j = position[i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 2] = saved_i + h; position[i3x2 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i; position[i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x2, 1, hval);
  }
  {
    double saved_i = position[i3x3 + 2];
    double saved_j = position[i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 2] = saved_i + h; position[i3x2 + 2] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i; position[i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x2, 2, hval);
  }
  {
    double saved_i = position[i3x3 + 2];
    double saved_j = position[i3x3 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 2] = saved_i + h; position[i3x3 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i; position[i3x3 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x3, 0, hval);
  }
  {
    double saved_i = position[i3x3 + 2];
    double saved_j = position[i3x3 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 2] = saved_i + h; position[i3x3 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x3 + 1] = saved_j + h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(kt, t0, i3x1, i3x2, i3x3, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i; position[i3x3 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x3, 1, hval);
  }
}

};
