#include "angle.h"

template <typename HESSIAN>
struct Angle {
  static constexpr size_t PositionSize = 9;
  static std::string description() { return "mathkernel-angle"; };
double energy(const angle_term& term, double* position, double* energy_accumulate) {
  double energy_added = 0.00000000000000000e+0;
  {
    /* !BASE */
    DOUBLE x1 = position[term.i3x1 + 0];
    DOUBLE y1 = position[term.i3x1 + 1];
    DOUBLE z1 = position[term.i3x1 + 2];
    DOUBLE x2 = position[term.i3x2 + 0];
    DOUBLE y2 = position[term.i3x2 + 1];
    DOUBLE z2 = position[term.i3x2 + 2];
    DOUBLE x3 = position[term.i3x3 + 0];
    DOUBLE y3 = position[term.i3x3 + 1];
    DOUBLE z3 = position[term.i3x3 + 2];
    {
      /* !BASE */
      double cse_p11_t1_g371772 = (-(x2));
      double cse_p11_t2_g371773 = (-(y2));
      double cse_p11_t3_g371774 = (-(z2));
      double vx1 = (cse_p11_t1_g371772 + x1);
      double vy1 = (cse_p11_t2_g371773 + y1);
      double vz1 = (cse_p11_t3_g371774 + z1);
      double vx2 = (cse_p11_t1_g371772 + x3);
      double vy2 = (cse_p11_t2_g371773 + y3);
      double vz2 = (cse_p11_t3_g371774 + z3);
      double dot = ((vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2));
      double n1_sq = ((vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1));
      double n2_sq = ((vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2));
      double n1 = sqrt(n1_sq);
      double n2 = sqrt(n2_sq);
      double cos_theta = (dot * (1.0 / ((n1 * n2))));
      if (fabs(cos_theta)>(1.0-VERYSMALL)) cos_theta = linear_angle_clamp(cos_theta);
      double sin_theta = sqrt((1.00000000000000000e+0 + (-(((cos_theta) * (cos_theta))))));
      double theta = acos(cos_theta);
      double dtheta = (theta + (-(term.t0)));
      double energy = (dtheta * dtheta * term.kt);
      energy_added += energy;
      *energy_accumulate += energy;
    }
  }
  return energy_added;
}
void energy_fd(const angle_term& term, double* position, double* energy_accumulate)
{
  energy(term, position, energy_accumulate);
}

double gradient(const angle_term& term, double* position, double* energy_accumulate, double* force) {
  double energy_added = 0.00000000000000000e+0;
  {
    /* !BASE */
    DOUBLE x1 = position[term.i3x1 + 0];
    DOUBLE y1 = position[term.i3x1 + 1];
    DOUBLE z1 = position[term.i3x1 + 2];
    DOUBLE x2 = position[term.i3x2 + 0];
    DOUBLE y2 = position[term.i3x2 + 1];
    DOUBLE z2 = position[term.i3x2 + 2];
    DOUBLE x3 = position[term.i3x3 + 0];
    DOUBLE y3 = position[term.i3x3 + 1];
    DOUBLE z3 = position[term.i3x3 + 2];
    {
      /* !BASE */
      double cse_p11_t1_g371783 = (-(x2));
      double cse_p11_t2_g371784 = (-(y2));
      double cse_p11_t3_g371785 = (-(z2));
      double vx1 = (cse_p11_t1_g371783 + x1);
      double vy1 = (cse_p11_t2_g371784 + y1);
      double vz1 = (cse_p11_t3_g371785 + z1);
      double vx2 = (cse_p11_t1_g371783 + x3);
      double vy2 = (cse_p11_t2_g371784 + y3);
      double vz2 = (cse_p11_t3_g371785 + z3);
      double dot = ((vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2));
      double n1_sq = ((vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1));
      double n2_sq = ((vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2));
      double cse_p1_t1_sqrt371775 = sqrt(n1_sq);
      double cse_p1_t2_invsqrt371776 = (1.0 / (cse_p1_t1_sqrt371775));
      double n1 = cse_p1_t1_sqrt371775;
      double cse_p1_t3_sqrt371777 = sqrt(n2_sq);
      double cse_p1_t4_invsqrt371778 = (1.0 / (cse_p1_t3_sqrt371777));
      double n2 = cse_p1_t3_sqrt371777;
      double cse_p1_t5_invr371779 = (1.0 / ((n1 * n2)));
      double cse_p1_t6_invr2371780 = (cse_p1_t5_invr371779 * cse_p1_t5_invr371779);
      double cse_p113_t1_g371794 = (cse_p1_t6_invr2371780 * dot);
      double cos_theta = (cse_p1_t5_invr371779 * dot);
      if (fabs(cos_theta)>(1.0-VERYSMALL)) cos_theta = linear_angle_clamp(cos_theta);
      double cse_p1_t7_sqrt371781 = sqrt((1.00000000000000000e+0 + (-(((cos_theta) * (cos_theta))))));
      double cse_p1_t8_invsqrt371782 = (1.0 / (cse_p1_t7_sqrt371781));
      double sin_theta = cse_p1_t7_sqrt371781;
      double theta = acos(cos_theta);
      double dtheta = (theta + (-(term.t0)));
      double cse_p11_t1_g371792 = (dtheta * term.kt);
      double energy = (cse_p11_t1_g371792 * dtheta);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t1_g371786 = (2.00000000000000000e+0 * cse_p11_t1_g371792);
      double cse_p11_t1_g371787 = (cse_p1_t2_invsqrt371776 * n2);
      double cse_p11_t1_g371790 = (cse_p113_t1_g371794 * cse_p11_t1_g371787);
      double cse_p62_t1_g371793 = (cse_p11_t1_g371786 * cse_p1_t8_invsqrt371782);
      double g_x1 = (-((cse_p62_t1_g371793 * ((cse_p1_t5_invr371779 * vx2) + (-((cse_p11_t1_g371790 * vx1)))))));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = (-((cse_p62_t1_g371793 * ((cse_p1_t5_invr371779 * vy2) + (-((cse_p11_t1_g371790 * vy1)))))));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = (-((cse_p62_t1_g371793 * ((cse_p1_t5_invr371779 * vz2) + (-((cse_p11_t1_g371790 * vz1)))))));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p11_t1_g371788 = (cse_p1_t4_invsqrt371778 * n1);
      double g_x2 = (-((cse_p62_t1_g371793 * ((cse_p1_t5_invr371779 * ((-(vx1)) + (-(vx2)))) + (-((cse_p113_t1_g371794 * ((-((cse_p11_t1_g371787 * vx1))) + (-((cse_p11_t1_g371788 * vx2)))))))))));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (-((cse_p62_t1_g371793 * ((cse_p1_t5_invr371779 * ((-(vy1)) + (-(vy2)))) + (-((cse_p113_t1_g371794 * ((-((cse_p11_t1_g371787 * vy1))) + (-((cse_p11_t1_g371788 * vy2)))))))))));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (-((cse_p62_t1_g371793 * ((cse_p1_t5_invr371779 * ((-(vz1)) + (-(vz2)))) + (-((cse_p113_t1_g371794 * ((-((cse_p11_t1_g371787 * vz1))) + (-((cse_p11_t1_g371788 * vz2)))))))))));
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double cse_p11_t1_g371791 = (cse_p113_t1_g371794 * cse_p11_t1_g371788);
      double g_x3 = (-((cse_p62_t1_g371793 * ((cse_p1_t5_invr371779 * vx1) + (-((cse_p11_t1_g371791 * vx2)))))));
      KernelGradientAcc(term.i3x3, 0, g_x3);
      double g_y3 = (-((cse_p62_t1_g371793 * ((cse_p1_t5_invr371779 * vy1) + (-((cse_p11_t1_g371791 * vy2)))))));
      KernelGradientAcc(term.i3x3, 1, g_y3);
      double g_z3 = (-((cse_p62_t1_g371793 * ((cse_p1_t5_invr371779 * vz1) + (-((cse_p11_t1_g371791 * vz2)))))));
      KernelGradientAcc(term.i3x3, 2, g_z3);
    }
  }
  return energy_added;
}
void gradient_fd(const angle_term& term, double* position, double* energy_accumulate, double* force)
{
  constexpr size_t PositionSize = 9;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  double e0 = 0.0;
  energy(term, position, &e0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[term.i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x1 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x1, 0, d);
  }
  {
    double saved = position[term.i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x1 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x1, 1, d);
  }
  {
    double saved = position[term.i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x1 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x1, 2, d);
  }
  {
    double saved = position[term.i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x2 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x2, 0, d);
  }
  {
    double saved = position[term.i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x2 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x2, 1, d);
  }
  {
    double saved = position[term.i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x2 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x2, 2, d);
  }
  {
    double saved = position[term.i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x3 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x3 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x3 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x3, 0, d);
  }
  {
    double saved = position[term.i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x3 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x3 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x3 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x3, 1, d);
  }
  {
    double saved = position[term.i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x3 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x3 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x3 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x3, 2, d);
  }
}

double hessian(const angle_term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double energy_added = 0.00000000000000000e+0;
  {
    /* !BASE */
    DOUBLE x1 = position[term.i3x1 + 0];
    DOUBLE y1 = position[term.i3x1 + 1];
    DOUBLE z1 = position[term.i3x1 + 2];
    DOUBLE x2 = position[term.i3x2 + 0];
    DOUBLE y2 = position[term.i3x2 + 1];
    DOUBLE z2 = position[term.i3x2 + 2];
    DOUBLE x3 = position[term.i3x3 + 0];
    DOUBLE y3 = position[term.i3x3 + 1];
    DOUBLE z3 = position[term.i3x3 + 2];
    {
      /* !BASE */
      double cse_p11_t186_g371988 = (-(x2));
      double cse_p11_t187_g371989 = (-(y2));
      double cse_p11_t188_g371990 = (-(z2));
      double vx1 = (cse_p11_t186_g371988 + x1);
      double cse_p11_t180_g371982 = (-(vx1));
      double vy1 = (cse_p11_t187_g371989 + y1);
      double cse_p11_t182_g371984 = (-(vy1));
      double vz1 = (cse_p11_t188_g371990 + z1);
      double cse_p11_t184_g371986 = (-(vz1));
      double vx2 = (cse_p11_t186_g371988 + x3);
      double cse_p11_t138_g371940 = (cse_p11_t180_g371982 + (-(vx2)));
      double cse_p11_t181_g371983 = (-(vx2));
      double vy2 = (cse_p11_t187_g371989 + y3);
      double cse_p11_t139_g371941 = (cse_p11_t182_g371984 + (-(vy2)));
      double cse_p11_t183_g371985 = (-(vy2));
      double vz2 = (cse_p11_t188_g371990 + z3);
      double cse_p11_t140_g371942 = (cse_p11_t184_g371986 + (-(vz2)));
      double cse_p11_t185_g371987 = (-(vz2));
      double cse_p167_t1_g372098 = (vx1 * vx2);
      double cse_p167_t1_g372102 = (vy1 * vy2);
      double cse_p167_t1_g372105 = (vz1 * vz2);
      double dot = (cse_p167_t1_g372098 + cse_p167_t1_g372102 + cse_p167_t1_g372105);
      double n1_sq = ((vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1));
      double cse_p529_t1_sqrt372183 = sqrt(n1_sq);
      double cse_p529_t2_invr372184 = (1.0 / (n1_sq));
      double cse_p322_t1_g372129 = (cse_p529_t2_invr372184 * cse_p529_t2_invr372184);
      double cse_p114_t1_g372057 = (-2.00000000000000000e+0 * cse_p322_t1_g372129);
      double cse_p63_t13_g372030 = (cse_p114_t1_g372057 + cse_p322_t1_g372129);
      double cse_p167_t1_g372099 = (vx2 * vx2);
      double cse_p167_t1_g372103 = (vy2 * vy2);
      double cse_p167_t1_g372106 = (vz2 * vz2);
      double n2_sq = (cse_p167_t1_g372099 + cse_p167_t1_g372103 + cse_p167_t1_g372106);
      double cse_p529_t3_sqrt372185 = sqrt(n2_sq);
      double cse_p529_t4_invr372186 = (1.0 / (n2_sq));
      double cse_p61_t6_invr2372010 = (cse_p529_t4_invr372186 * cse_p529_t4_invr372186);
      double cse_p323_t1_g372130 = (cse_p529_t3_sqrt372185 * cse_p529_t4_invr372186);
      double cse_p11_t115_g371917 = (cse_p167_t1_g372098 * cse_p323_t1_g372130);
      double cse_p63_t1_g372038 = (cse_p323_t1_g372130 * cse_p529_t2_invr372184);
      double cse_p11_t116_g371918 = (cse_p63_t1_g372038 * vx2);
      double cse_p11_t117_g371919 = (cse_p167_t1_g372102 * cse_p323_t1_g372130);
      double cse_p11_t118_g371920 = (cse_p63_t1_g372038 * vy2);
      double cse_p11_t119_g371921 = (cse_p167_t1_g372105 * cse_p323_t1_g372130);
      double cse_p11_t120_g371922 = (cse_p63_t1_g372038 * vz2);
      double cse_p1_t2_invsqrt371796 = (1.0 / (cse_p529_t1_sqrt372183));
      double n1 = cse_p529_t1_sqrt372183;
      double cse_p63_t1_g372031 = (cse_p529_t3_sqrt372185 * n1);
      double cse_p478_t1_g372173 = (cse_p529_t4_invr372186 * cse_p63_t1_g372031);
      double cse_p64_t1_g372046 = (cse_p61_t6_invr2372010 * cse_p63_t1_g372031);
      double cse_p116_t1_g372076 = (2.00000000000000000e+0 * cse_p64_t1_g372046);
      double cse_p167_t1_g372100 = (vx2 * vy2);
      double cse_p11_t20_g371822 = (cse_p116_t1_g372076 * cse_p167_t1_g372100);
      double cse_p167_t1_g372101 = (vx2 * vz2);
      double cse_p11_t21_g371823 = (cse_p116_t1_g372076 * cse_p167_t1_g372101);
      double cse_p167_t1_g372104 = (vy2 * vz2);
      double cse_p11_t22_g371824 = (cse_p116_t1_g372076 * cse_p167_t1_g372104);
      double cse_p11_t91_g371893 = (cse_p478_t1_g372173 * vx2);
      double cse_p116_t1_g372070 = (cse_p478_t1_g372173 * cse_p529_t4_invr372186);
      double cse_p11_t92_g371894 = (cse_p116_t1_g372070 * cse_p167_t1_g372099);
      double cse_p11_t94_g371896 = (cse_p116_t1_g372070 * cse_p167_t1_g372100);
      double cse_p11_t96_g371898 = (cse_p116_t1_g372070 * cse_p167_t1_g372101);
      double cse_p11_t97_g371899 = (cse_p478_t1_g372173 * vy2);
      double cse_p11_t98_g371900 = (cse_p116_t1_g372070 * cse_p167_t1_g372103);
      double cse_p11_t100_g371902 = (cse_p116_t1_g372070 * cse_p167_t1_g372104);
      double cse_p11_t101_g371903 = (cse_p478_t1_g372173 * vz2);
      double cse_p11_t102_g371904 = (cse_p116_t1_g372070 * cse_p167_t1_g372106);
      double cse_p11_t162_g371964 = (-(cse_p11_t91_g371893));
      double cse_p11_t164_g371966 = (-(cse_p11_t94_g371896));
      double cse_p11_t166_g371968 = (-(cse_p11_t96_g371898));
      double cse_p11_t167_g371969 = (-(cse_p11_t97_g371899));
      double cse_p11_t169_g371971 = (-(cse_p11_t100_g371902));
      double cse_p11_t170_g371972 = (-(cse_p11_t101_g371903));
      double cse_p11_t177_g371979 = (-(cse_p478_t1_g372173));
      double cse_p1_t4_invsqrt371798 = (1.0 / (cse_p529_t3_sqrt372185));
      double n2 = cse_p529_t3_sqrt372185;
      double cse_p116_t1_g372077 = (cse_p529_t1_sqrt372183 * cse_p529_t2_invr372184);
      double cse_p11_t1_g372001 = (cse_p116_t1_g372077 * n2);
      double cse_p426_t1_g372165 = (cse_p11_t1_g372001 * vz1);
      double cse_p270_t1_g372122 = (cse_p11_t1_g372001 * vy1);
      double cse_p166_t1_g372097 = (cse_p11_t1_g372001 * vx1);
      double cse_p61_t7_invr372011 = (1.0 / ((n1 * n2)));
      double cse_p61_t8_invr2372012 = (cse_p61_t7_invr372011 * cse_p61_t7_invr372011);
      double cse_p374_t1_g372137 = (cse_p61_t8_invr2372012 * dot);
      double cse_p61_t9_invr3372013 = (cse_p61_t7_invr372011 * cse_p61_t8_invr2372012);
      double cse_p63_t1_g372035 = (2.00000000000000000e+0 * cse_p61_t9_invr3372013);
      double cse_p63_t1_g372037 = (cse_p478_t1_g372173 * cse_p63_t1_g372035);
      double cse_p63_t5_g372022 = (cse_p529_t4_invr372186 * cse_p63_t1_g372037);
      double cse_p116_t1_g372075 = (2.00000000000000000e+0 * cse_p61_t8_invr2372012);
      double cse_p63_t12_g372029 = (cse_p63_t5_g372022 + (cse_p116_t1_g372075 * cse_p61_t6_invr2372010));
      double cse_p11_t1_g371803 = (cse_p61_t8_invr2372012 * (cse_p11_t162_g371964 + (-(cse_p166_t1_g372097))));
      double cse_p11_t2_g371804 = (cse_p61_t8_invr2372012 * (cse_p11_t167_g371969 + (-(cse_p270_t1_g372122))));
      double cse_p11_t3_g371805 = (cse_p61_t8_invr2372012 * (cse_p11_t170_g371972 + (-(cse_p426_t1_g372165))));
      double cse_p115_t1_g372067 = (n2 * vx1);
      double cse_p11_t13_g371815 = (cse_p114_t1_g372057 * cse_p115_t1_g372067);
      double cse_p116_t1_g372073 = (n2 * vy1);
      double cse_p11_t14_g371816 = (cse_p114_t1_g372057 * cse_p116_t1_g372073);
      double cse_p117_t1_g372087 = (n2 * vz1);
      double cse_p11_t15_g371817 = (cse_p114_t1_g372057 * cse_p117_t1_g372087);
      double cse_p11_t16_g371818 = (cse_p116_t1_g372075 * cse_p322_t1_g372129);
      double cse_p11_t17_g371819 = (cse_p116_t1_g372075 * cse_p61_t6_invr2372010);
      double cse_p11_t19_g371821 = (cse_p116_t1_g372070 * cse_p63_t1_g372035);
      double cse_p117_t1_g372088 = (cse_p529_t2_invr372184 * cse_p63_t1_g372035);
      double cse_p11_t23_g371825 = (cse_p117_t1_g372088 * cse_p11_t1_g372001);
      double cse_p11_t33_g371835 = (dot * (-(cse_p11_t1_g371803)));
      double cse_p11_t34_g371836 = (dot * (-(cse_p11_t2_g371804)));
      double cse_p11_t35_g371837 = (dot * (-(cse_p11_t3_g371805)));
      double cse_p218_t1_g372115 = (cse_p374_t1_g372137 * cse_p478_t1_g372173);
      double cse_p11_t37_g371839 = (cse_p218_t1_g372115 * vx2);
      double cse_p11_t38_g371840 = (cse_p218_t1_g372115 * vy2);
      double cse_p11_t39_g371841 = (cse_p218_t1_g372115 * vz2);
      double cse_p11_t40_g371842 = (cse_p11_t1_g372001 * cse_p374_t1_g372137);
      double cse_p11_t41_g371843 = (cse_p166_t1_g372097 * cse_p374_t1_g372137);
      double cse_p11_t42_g371844 = (cse_p270_t1_g372122 * cse_p374_t1_g372137);
      double cse_p11_t43_g371845 = (cse_p374_t1_g372137 * cse_p426_t1_g372165);
      double cse_p11_t90_g371892 = (n1 * n2);
      double cse_p377_t1_g372158 = (cse_p478_t1_g372173 * cse_p61_t8_invr2372012);
      double cse_p11_t93_g371895 = (cse_p167_t1_g372100 * cse_p377_t1_g372158);
      double cse_p11_t95_g371897 = (cse_p167_t1_g372101 * cse_p377_t1_g372158);
      double cse_p11_t99_g371901 = (cse_p167_t1_g372104 * cse_p377_t1_g372158);
      double cse_p11_t103_g371905 = (cse_p115_t1_g372067 * cse_p322_t1_g372129);
      double cse_p11_t105_g371907 = (cse_p115_t1_g372067 * cse_p529_t2_invr372184 * vx1);
      double cse_p11_t106_g371908 = (cse_p166_t1_g372097 * cse_p61_t8_invr2372012 * vy1);
      double cse_p11_t107_g371909 = (cse_p166_t1_g372097 * cse_p61_t8_invr2372012 * vz1);
      double cse_p11_t108_g371910 = (cse_p116_t1_g372073 * cse_p322_t1_g372129);
      double cse_p11_t110_g371912 = (cse_p116_t1_g372073 * cse_p529_t2_invr372184 * vy1);
      double cse_p11_t111_g371913 = (cse_p270_t1_g372122 * cse_p61_t8_invr2372012 * vz1);
      double cse_p11_t112_g371914 = (cse_p117_t1_g372087 * cse_p322_t1_g372129);
      double cse_p11_t114_g371916 = (cse_p117_t1_g372087 * cse_p529_t2_invr372184 * vz1);
      double cse_p11_t121_g371923 = (cse_p11_t162_g371964 + (-(cse_p166_t1_g372097)));
      double cse_p63_t1_g372039 = (cse_p117_t1_g372088 * cse_p11_t121_g371923);
      double cse_p63_t2_g372019 = (cse_p63_t1_g372039 * n2);
      double cse_p11_t122_g371924 = (cse_p11_t167_g371969 + (-(cse_p270_t1_g372122)));
      double cse_p63_t1_g372040 = (cse_p117_t1_g372088 * cse_p11_t122_g371924);
      double cse_p63_t3_g372020 = (cse_p63_t1_g372040 * n2);
      double cse_p11_t123_g371925 = (cse_p11_t170_g371972 + (-(cse_p426_t1_g372165)));
      double cse_p63_t1_g372041 = (cse_p117_t1_g372088 * cse_p11_t123_g371925);
      double cse_p63_t4_g372021 = (cse_p63_t1_g372041 * n2);
      double cse_p11_t124_g371926 = (cse_p11_t103_g371905 + cse_p11_t116_g371918 + cse_p11_t13_g371815);
      double cse_p63_t6_g372023 = (cse_p11_t124_g371926 * cse_p61_t8_invr2372012);
      double cse_p63_t9_g372026 = (cse_p63_t2_g372019 + cse_p63_t6_g372023);
      double cse_p11_t125_g371927 = (cse_p11_t108_g371910 + cse_p11_t118_g371920 + cse_p11_t14_g371816);
      double cse_p63_t7_g372024 = (cse_p11_t125_g371927 * cse_p61_t8_invr2372012);
      double cse_p63_t10_g372027 = (cse_p63_t3_g372020 + cse_p63_t7_g372024);
      double cse_p11_t126_g371928 = (cse_p11_t112_g371914 + cse_p11_t120_g371922 + cse_p11_t15_g371817);
      double cse_p63_t8_g372025 = (cse_p11_t126_g371928 * cse_p61_t8_invr2372012);
      double cse_p63_t11_g372028 = (cse_p63_t4_g372021 + cse_p63_t8_g372025);
      double cse_p11_t127_g371929 = (cse_p11_t16_g371818 + cse_p11_t23_g371825);
      double cse_p11_t128_g371930 = (cse_p11_t17_g371819 + cse_p11_t19_g371821);
      double cse_p11_t142_g371944 = (-((cse_p11_t121_g371923 * cse_p61_t8_invr2372012)));
      double cse_p11_t143_g371945 = (-((cse_p11_t122_g371924 * cse_p61_t8_invr2372012)));
      double cse_p11_t144_g371946 = (-((cse_p11_t123_g371925 * cse_p61_t8_invr2372012)));
      double cse_p11_t154_g371956 = (-(cse_p218_t1_g372115));
      double cse_p11_t155_g371957 = (-(cse_p11_t37_g371839));
      double cse_p11_t156_g371958 = (-(cse_p11_t38_g371840));
      double cse_p11_t157_g371959 = (-(cse_p11_t39_g371841));
      double cse_p11_t158_g371960 = (-(cse_p11_t40_g371842));
      double cse_p11_t159_g371961 = (-(cse_p11_t41_g371843));
      double cse_p11_t160_g371962 = (-(cse_p11_t42_g371844));
      double cse_p11_t161_g371963 = (-(cse_p11_t43_g371845));
      double cse_p11_t163_g371965 = (-(cse_p11_t93_g371895));
      double cse_p11_t165_g371967 = (-(cse_p11_t95_g371897));
      double cse_p11_t168_g371970 = (-(cse_p11_t99_g371901));
      double cse_p11_t171_g371973 = (-(cse_p166_t1_g372097));
      double cse_p11_t172_g371974 = (-(cse_p11_t106_g371908));
      double cse_p11_t173_g371975 = (-(cse_p11_t107_g371909));
      double cse_p11_t174_g371976 = (-(cse_p270_t1_g372122));
      double cse_p11_t175_g371977 = (-(cse_p11_t111_g371913));
      double cse_p11_t176_g371978 = (-(cse_p426_t1_g372165));
      double cse_p1_t5_invr371799 = (1.0 / (cse_p11_t90_g371892));
      double cse_p11_t18_g371820 = (2.00000000000000000e+0 * cse_p1_t5_invr371799);
      double cse_p11_t24_g371826 = (cse_p1_t5_invr371799 * (cse_p11_t180_g371982 + cse_p11_t181_g371983));
      double cse_p11_t25_g371827 = (cse_p1_t5_invr371799 * (cse_p11_t182_g371984 + cse_p11_t183_g371985));
      double cse_p11_t26_g371828 = (cse_p1_t5_invr371799 * (cse_p11_t184_g371986 + cse_p11_t185_g371987));
      double cse_p11_t27_g371829 = (cse_p1_t5_invr371799 * vx1);
      double cse_p11_t28_g371830 = (cse_p1_t5_invr371799 * vx2);
      double cse_p11_t29_g371831 = (cse_p1_t5_invr371799 * vy1);
      double cse_p11_t30_g371832 = (cse_p1_t5_invr371799 * vy2);
      double cse_p11_t31_g371833 = (cse_p1_t5_invr371799 * vz1);
      double cse_p11_t32_g371834 = (cse_p1_t5_invr371799 * vz2);
      double cse_p11_t129_g371931 = (cse_p11_t24_g371826 + (-((cse_p374_t1_g372137 * (cse_p11_t162_g371964 + cse_p11_t171_g371973)))));
      double cse_p11_t130_g371932 = (cse_p11_t25_g371827 + (-((cse_p374_t1_g372137 * (cse_p11_t167_g371969 + cse_p11_t174_g371976)))));
      double cse_p11_t131_g371933 = (cse_p11_t26_g371828 + (-((cse_p374_t1_g372137 * (cse_p11_t170_g371972 + cse_p11_t176_g371978)))));
      double cse_p11_t132_g371934 = (cse_p11_t155_g371957 + cse_p11_t27_g371829);
      double cse_p11_t133_g371935 = (cse_p11_t159_g371961 + cse_p11_t28_g371830);
      double cse_p11_t134_g371936 = (cse_p11_t156_g371958 + cse_p11_t29_g371831);
      double cse_p11_t135_g371937 = (cse_p11_t160_g371962 + cse_p11_t30_g371832);
      double cse_p11_t136_g371938 = (cse_p11_t157_g371959 + cse_p11_t31_g371833);
      double cse_p11_t137_g371939 = (cse_p11_t161_g371963 + cse_p11_t32_g371834);
      double cse_p11_t179_g371981 = (-(cse_p1_t5_invr371799));
      double cse_p1_t6_invr2371800 = (cse_p1_t5_invr371799 * cse_p1_t5_invr371799);
      double cos_theta = (cse_p1_t5_invr371799 * dot);
      double cse_p61_t10_invsqrt372014 = (1.0 / sqrt((1.00000000000000000e+0 + (-(((cos_theta) * (cos_theta)))))));
      double cse_p11_t4_g371806 = (cse_p11_t129_g371931 * cse_p61_t10_invsqrt372014);
      double cse_p11_t5_g371807 = (cse_p11_t130_g371932 * cse_p61_t10_invsqrt372014);
      double cse_p11_t6_g371808 = (cse_p11_t131_g371933 * cse_p61_t10_invsqrt372014);
      double cse_p11_t7_g371809 = (cse_p11_t132_g371934 * cse_p61_t10_invsqrt372014);
      double cse_p11_t8_g371810 = (cse_p11_t133_g371935 * cse_p61_t10_invsqrt372014);
      double cse_p11_t9_g371811 = (cse_p11_t134_g371936 * cse_p61_t10_invsqrt372014);
      double cse_p11_t10_g371812 = (cse_p11_t135_g371937 * cse_p61_t10_invsqrt372014);
      double cse_p11_t11_g371813 = (cse_p11_t136_g371938 * cse_p61_t10_invsqrt372014);
      double cse_p11_t12_g371814 = (cse_p11_t137_g371939 * cse_p61_t10_invsqrt372014);
      double cse_p63_t1_g372032 = (cse_p61_t10_invsqrt372014 * cse_p61_t10_invsqrt372014 * term.kt);
      double cse_p115_t1_g372060 = (cse_p11_t129_g371931 * cse_p63_t1_g372032);
      double cse_p11_t44_g371846 = (cse_p115_t1_g372060 * cse_p11_t129_g371931);
      double cse_p11_t45_g371847 = (cse_p115_t1_g372060 * cse_p11_t130_g371932);
      double cse_p11_t46_g371848 = (cse_p115_t1_g372060 * cse_p11_t131_g371933);
      double cse_p11_t47_g371849 = (cse_p115_t1_g372060 * cse_p11_t132_g371934);
      double cse_p11_t48_g371850 = (cse_p115_t1_g372060 * cse_p11_t133_g371935);
      double cse_p11_t49_g371851 = (cse_p115_t1_g372060 * cse_p11_t134_g371936);
      double cse_p11_t50_g371852 = (cse_p115_t1_g372060 * cse_p11_t135_g371937);
      double cse_p11_t51_g371853 = (cse_p115_t1_g372060 * cse_p11_t136_g371938);
      double cse_p11_t52_g371854 = (cse_p115_t1_g372060 * cse_p11_t137_g371939);
      double cse_p115_t1_g372064 = (cse_p11_t130_g371932 * cse_p63_t1_g372032);
      double cse_p11_t53_g371855 = (cse_p115_t1_g372064 * cse_p11_t130_g371932);
      double cse_p11_t54_g371856 = (cse_p115_t1_g372064 * cse_p11_t131_g371933);
      double cse_p11_t55_g371857 = (cse_p115_t1_g372064 * cse_p11_t132_g371934);
      double cse_p11_t56_g371858 = (cse_p115_t1_g372064 * cse_p11_t133_g371935);
      double cse_p11_t57_g371859 = (cse_p115_t1_g372064 * cse_p11_t134_g371936);
      double cse_p11_t58_g371860 = (cse_p115_t1_g372064 * cse_p11_t135_g371937);
      double cse_p11_t59_g371861 = (cse_p115_t1_g372064 * cse_p11_t136_g371938);
      double cse_p11_t60_g371862 = (cse_p115_t1_g372064 * cse_p11_t137_g371939);
      double cse_p116_t1_g372071 = (cse_p11_t131_g371933 * cse_p63_t1_g372032);
      double cse_p11_t61_g371863 = (cse_p116_t1_g372071 * cse_p11_t131_g371933);
      double cse_p11_t62_g371864 = (cse_p116_t1_g372071 * cse_p11_t132_g371934);
      double cse_p11_t63_g371865 = (cse_p116_t1_g372071 * cse_p11_t133_g371935);
      double cse_p11_t64_g371866 = (cse_p116_t1_g372071 * cse_p11_t134_g371936);
      double cse_p11_t65_g371867 = (cse_p116_t1_g372071 * cse_p11_t135_g371937);
      double cse_p11_t66_g371868 = (cse_p116_t1_g372071 * cse_p11_t136_g371938);
      double cse_p11_t67_g371869 = (cse_p116_t1_g372071 * cse_p11_t137_g371939);
      double cse_p116_t1_g372078 = (cse_p11_t132_g371934 * cse_p63_t1_g372032);
      double cse_p11_t68_g371870 = (cse_p116_t1_g372078 * cse_p11_t132_g371934);
      double cse_p11_t69_g371871 = (cse_p116_t1_g372078 * cse_p11_t133_g371935);
      double cse_p11_t70_g371872 = (cse_p116_t1_g372078 * cse_p11_t134_g371936);
      double cse_p11_t71_g371873 = (cse_p116_t1_g372078 * cse_p11_t135_g371937);
      double cse_p11_t72_g371874 = (cse_p116_t1_g372078 * cse_p11_t136_g371938);
      double cse_p11_t73_g371875 = (cse_p116_t1_g372078 * cse_p11_t137_g371939);
      double cse_p117_t1_g372089 = (cse_p11_t133_g371935 * cse_p63_t1_g372032);
      double cse_p11_t74_g371876 = (cse_p117_t1_g372089 * cse_p11_t133_g371935);
      double cse_p11_t75_g371877 = (cse_p117_t1_g372089 * cse_p11_t134_g371936);
      double cse_p11_t76_g371878 = (cse_p117_t1_g372089 * cse_p11_t135_g371937);
      double cse_p11_t77_g371879 = (cse_p117_t1_g372089 * cse_p11_t136_g371938);
      double cse_p11_t78_g371880 = (cse_p117_t1_g372089 * cse_p11_t137_g371939);
      double cse_p167_t1_g372107 = (cse_p11_t134_g371936 * cse_p63_t1_g372032);
      double cse_p11_t79_g371881 = (cse_p11_t134_g371936 * cse_p167_t1_g372107);
      double cse_p11_t80_g371882 = (cse_p11_t135_g371937 * cse_p167_t1_g372107);
      double cse_p11_t81_g371883 = (cse_p11_t136_g371938 * cse_p167_t1_g372107);
      double cse_p11_t82_g371884 = (cse_p11_t137_g371939 * cse_p167_t1_g372107);
      double cse_p376_t1_g372156 = (cse_p11_t135_g371937 * cse_p63_t1_g372032);
      double cse_p11_t83_g371885 = (cse_p11_t135_g371937 * cse_p376_t1_g372156);
      double cse_p11_t84_g371886 = (cse_p11_t136_g371938 * cse_p376_t1_g372156);
      double cse_p11_t85_g371887 = (cse_p11_t137_g371939 * cse_p376_t1_g372156);
      double cse_p11_t86_g371888 = (cse_p11_t136_g371938 * cse_p11_t136_g371938 * cse_p63_t1_g372032);
      double cse_p11_t87_g371889 = (cse_p11_t136_g371938 * cse_p11_t137_g371939 * cse_p63_t1_g372032);
      double cse_p11_t88_g371890 = (cse_p11_t137_g371939 * cse_p11_t137_g371939 * cse_p63_t1_g372032);
      double cse_p11_t141_g371943 = (1.00000000000000000e+0 + (-(((cos_theta) * (cos_theta)))));
      double cse_p11_t145_g371947 = (-(cse_p11_t4_g371806));
      double cse_p11_t146_g371948 = (-(cse_p11_t5_g371807));
      double cse_p11_t147_g371949 = (-(cse_p11_t6_g371808));
      double cse_p11_t148_g371950 = (-(cse_p11_t7_g371809));
      double cse_p11_t149_g371951 = (-(cse_p11_t8_g371810));
      double cse_p11_t150_g371952 = (-(cse_p11_t9_g371811));
      double cse_p11_t151_g371953 = (-(cse_p11_t10_g371812));
      double cse_p11_t152_g371954 = (-(cse_p11_t11_g371813));
      double cse_p11_t153_g371955 = (-(cse_p11_t12_g371814));
      double cse_p11_t178_g371980 = (-(((cos_theta) * (cos_theta))));
      double cse_p529_t5_sqrt372187 = sqrt(cse_p11_t141_g371943);
      double cse_p529_t6_invr372188 = (1.0 / (cse_p11_t141_g371943));
      double cse_p11_t193_g371995 = ((cos_theta) * (cos_theta));
      if (fabs(cos_theta)>(1.0-VERYSMALL)) cos_theta = linear_angle_clamp(cos_theta);
      double cse_p1_t8_invsqrt371802 = (1.0 / (cse_p529_t5_sqrt372187));
      double sin_theta = cse_p529_t5_sqrt372187;
      double theta = acos(cos_theta);
      double dtheta = (theta + (-(term.t0)));
      double cse_p62_t1_g372017 = (2.00000000000000000e+0 * dtheta * term.kt);
      double energy = (dtheta * dtheta * term.kt);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p64_t1_g372042 = (dot * n2);
      double cse_p64_t1_g372050 = (cse_p1_t2_invsqrt371796 * cse_p1_t6_invr2371800 * cse_p64_t1_g372042);
      double cse_p115_t1_g372062 = (cse_p1_t8_invsqrt371802 * cse_p62_t1_g372017);
      double g_x1 = (-((cse_p115_t1_g372062 * (cse_p11_t28_g371830 + (-((cse_p64_t1_g372050 * vx1)))))));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = (-((cse_p115_t1_g372062 * (cse_p11_t30_g371832 + (-((cse_p64_t1_g372050 * vy1)))))));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = (-((cse_p115_t1_g372062 * (cse_p11_t32_g371834 + (-((cse_p64_t1_g372050 * vz1)))))));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p117_t1_g372081 = (cse_p1_t4_invsqrt371798 * n1);
      double cse_p117_t1_g372082 = (cse_p1_t6_invr2371800 * dot);
      double g_x2 = (-((cse_p115_t1_g372062 * (cse_p11_t24_g371826 + (-((cse_p117_t1_g372082 * ((-((cse_p115_t1_g372067 * cse_p1_t2_invsqrt371796))) + (-((cse_p117_t1_g372081 * vx2)))))))))));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (-((cse_p115_t1_g372062 * (cse_p11_t25_g371827 + (-((cse_p117_t1_g372082 * ((-((cse_p116_t1_g372073 * cse_p1_t2_invsqrt371796))) + (-((cse_p117_t1_g372081 * vy2)))))))))));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (-((cse_p115_t1_g372062 * (cse_p11_t26_g371828 + (-((cse_p117_t1_g372082 * ((-((cse_p117_t1_g372081 * vz2))) + (-((cse_p117_t1_g372087 * cse_p1_t2_invsqrt371796)))))))))));
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double cse_p376_t1_g372149 = (cse_p117_t1_g372081 * cse_p117_t1_g372082);
      double g_x3 = (-((cse_p115_t1_g372062 * (cse_p11_t27_g371829 + (-((cse_p376_t1_g372149 * vx2)))))));
      KernelGradientAcc(term.i3x3, 0, g_x3);
      double g_y3 = (-((cse_p115_t1_g372062 * (cse_p11_t29_g371831 + (-((cse_p376_t1_g372149 * vy2)))))));
      KernelGradientAcc(term.i3x3, 1, g_y3);
      double g_z3 = (-((cse_p115_t1_g372062 * (cse_p11_t31_g371833 + (-((cse_p376_t1_g372149 * vz2)))))));
      KernelGradientAcc(term.i3x3, 2, g_z3);
      double cse_p11_t1_g372003 = (cse_p116_t1_g372077 * cse_p61_t8_invr2372012);
      double cse_p63_t1_g372033 = (cos_theta * cse_p529_t6_invr372188);
      double cse_p63_t1_g372034 = (cse_p61_t10_invsqrt372014 * cse_p62_t1_g372017);
      double cse_p64_t1_g372048 = (cse_p11_t1_g372003 * cse_p529_t2_invr372184 * cse_p64_t1_g372042);
      double cse_p530_t1_g372189 = (cse_p64_t1_g372048 * vx1);
      double cse_p116_t1_g372069 = (cse_p529_t1_sqrt372183 * cse_p64_t1_g372042);
      double cse_p64_t1_g372049 = (cse_p116_t1_g372069 * cse_p116_t1_g372077 * cse_p117_t1_g372088);
      double cse_p115_t1_g372059 = (cse_p529_t1_sqrt372183 * dot);
      double cse_p116_t1_g372079 = (cse_p115_t1_g372059 * cse_p116_t1_g372075 * cse_p322_t1_g372129);
      double cse_p118_t1_g372090 = (cse_p11_t133_g371935 * cse_p63_t1_g372033);
      double cse_p375_t1_g372144 = (-2.00000000000000000e+0 * cse_p11_t1_g372003);
      double h_x1_x1 = (cse_p11_t74_g371876 + cse_p11_t74_g371876 + (-((cse_p63_t1_g372034 * (cse_p11_t158_g371960 + (cse_p115_t1_g372067 * ((cse_p116_t1_g372079 * vx1) + (cse_p375_t1_g372144 * vx2) + (cse_p64_t1_g372049 * vx1))) + (cse_p118_t1_g372090 * cse_p11_t133_g371935) + (-((cse_p530_t1_g372189 * vx1))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
      double cse_p64_t1_g372043 = (cse_p11_t1_g372003 * n2);
      double cse_p117_t1_g372084 = (cse_p64_t1_g372043 * vx1);
      double cse_p117_t1_g372085 = (cse_p64_t1_g372043 * vy1);
      double cse_p376_t1_g372148 = (cse_p116_t1_g372069 * vx1);
      double h_x1_y1 = (cse_p11_t76_g371878 + cse_p11_t76_g371878 + (-((cse_p63_t1_g372034 * ((cse_p118_t1_g372090 * cse_p11_t135_g371937) + (cse_p11_t127_g371929 * cse_p376_t1_g372148 * vy1) + (-((cse_p117_t1_g372084 * vy2))) + (-((cse_p117_t1_g372085 * vx2))) + (-((cse_p530_t1_g372189 * vy1))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
      double cse_p117_t1_g372086 = (cse_p64_t1_g372043 * vz1);
      double h_x1_z1 = (cse_p11_t78_g371880 + cse_p11_t78_g371880 + (-((cse_p63_t1_g372034 * ((cse_p118_t1_g372090 * cse_p11_t137_g371939) + (cse_p11_t127_g371929 * cse_p376_t1_g372148 * vz1) + (-((cse_p117_t1_g372084 * vz2))) + (-((cse_p117_t1_g372086 * vx2))) + (-((cse_p530_t1_g372189 * vz1))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
      double cse_p115_t1_g372061 = (cse_p11_t129_g371931 * cse_p63_t1_g372033);
      double cse_p375_t1_g372146 = (cse_p114_t1_g372057 * cse_p61_t8_invr2372012);
      double cse_p427_t1_g372166 = (cse_p11_t1_g372003 * dot);
      double h_x1_x2 = (cse_p11_t48_g371850 + cse_p11_t48_g371850 + (-((cse_p63_t1_g372034 * (cse_p11_t179_g371981 + (cse_p115_t1_g372061 * cse_p11_t133_g371935) + (cse_p11_t142_g371944 * vx2) + (cse_p376_t1_g372148 * (cse_p63_t1_g372039 + (cse_p375_t1_g372146 * vx1))) + (cse_p427_t1_g372166 * (cse_p11_t105_g371907 + cse_p11_t115_g371917 + n2)) + (-((cse_p117_t1_g372084 * cse_p11_t138_g371940))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
      double cse_p115_t1_g372065 = (cse_p11_t130_g371932 * cse_p63_t1_g372033);
      double h_x1_y2 = (cse_p11_t56_g371858 + cse_p11_t56_g371858 + (-((cse_p63_t1_g372034 * ((cse_p115_t1_g372059 * cse_p63_t10_g372027 * vx1) + (cse_p115_t1_g372065 * cse_p11_t133_g371935) + (cse_p11_t143_g371945 * vx2) + (-((cse_p117_t1_g372084 * cse_p11_t139_g371941))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
      double cse_p116_t1_g372072 = (cse_p11_t131_g371933 * cse_p63_t1_g372033);
      double h_x1_z2 = (cse_p11_t63_g371865 + cse_p11_t63_g371865 + (-((cse_p63_t1_g372034 * ((cse_p115_t1_g372059 * cse_p63_t11_g372028 * vx1) + (cse_p116_t1_g372072 * cse_p11_t133_g371935) + (cse_p11_t144_g371946 * vx2) + (-((cse_p117_t1_g372084 * cse_p11_t140_g371942))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
      double cse_p11_t1_g372004 = (cse_p529_t1_sqrt372183 * cse_p63_t1_g372038);
      double cse_p64_t1_g372044 = (cse_p529_t4_invr372186 * dot);
      double cse_p64_t1_g372045 = (cse_p11_t1_g372004 * cse_p63_t1_g372035 * cse_p64_t1_g372042 * n1);
      double cse_p64_t1_g372047 = (cse_p11_t1_g372003 * cse_p529_t3_sqrt372185 * cse_p64_t1_g372044);
      double cse_p115_t1_g372058 = (cse_p377_t1_g372158 * vx2);
      double cse_p117_t1_g372080 = (cse_p11_t132_g371934 * cse_p63_t1_g372033);
      double h_x1_x3 = (cse_p11_t69_g371871 + cse_p11_t69_g371871 + (-((cse_p63_t1_g372034 * (cse_p1_t5_invr371799 + (cse_p117_t1_g372080 * cse_p11_t133_g371935) + (cse_p167_t1_g372098 * cse_p64_t1_g372045) + (-((cse_p115_t1_g372058 * vx2))) + (-((cse_p117_t1_g372084 * vx1))) + (-((cse_p167_t1_g372098 * cse_p64_t1_g372047))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 0, h_x1_x3);
      double cse_p375_t1_g372138 = (vx1 * vy2);
      double h_x1_y3 = (cse_p11_t75_g371877 + cse_p11_t75_g371877 + (-((cse_p63_t1_g372034 * (cse_p11_t163_g371965 + cse_p11_t172_g371974 + (cse_p118_t1_g372090 * cse_p11_t134_g371936) + (cse_p375_t1_g372138 * cse_p64_t1_g372045) + (-((cse_p375_t1_g372138 * cse_p64_t1_g372047))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 1, h_x1_y3);
      double cse_p375_t1_g372139 = (vx1 * vz2);
      double h_x1_z3 = (cse_p11_t77_g371879 + cse_p11_t77_g371879 + (-((cse_p63_t1_g372034 * (cse_p11_t165_g371967 + cse_p11_t173_g371975 + (cse_p118_t1_g372090 * cse_p11_t136_g371938) + (cse_p375_t1_g372139 * cse_p64_t1_g372045) + (-((cse_p375_t1_g372139 * cse_p64_t1_g372047))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 2, h_x1_z3);
      double cse_p376_t1_g372157 = (cse_p11_t135_g371937 * cse_p63_t1_g372033);
      double h_y1_y1 = (cse_p11_t83_g371885 + cse_p11_t83_g371885 + (-((cse_p63_t1_g372034 * (cse_p11_t158_g371960 + (cse_p116_t1_g372073 * ((cse_p116_t1_g372079 * vy1) + (cse_p375_t1_g372144 * vy2) + (cse_p64_t1_g372049 * vy1))) + (cse_p11_t135_g371937 * cse_p376_t1_g372157) + (-((cse_p64_t1_g372048 * vy1 * vy1))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
      double h_y1_z1 = (cse_p11_t85_g371887 + cse_p11_t85_g371887 + (-((cse_p63_t1_g372034 * ((cse_p116_t1_g372069 * cse_p11_t127_g371929 * vy1 * vz1) + (cse_p11_t137_g371939 * cse_p376_t1_g372157) + (-((cse_p117_t1_g372085 * vz2))) + (-((cse_p117_t1_g372086 * vy2))) + (-((cse_p64_t1_g372048 * vy1 * vz1))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
      double h_y1_x2 = (cse_p11_t50_g371852 + cse_p11_t50_g371852 + (-((cse_p63_t1_g372034 * ((cse_p115_t1_g372059 * cse_p63_t9_g372026 * vy1) + (cse_p115_t1_g372061 * cse_p11_t135_g371937) + (cse_p11_t142_g371944 * vy2) + (-((cse_p117_t1_g372085 * cse_p11_t138_g371940))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
      double h_y1_y2 = (cse_p11_t58_g371860 + cse_p11_t58_g371860 + (-((cse_p63_t1_g372034 * (cse_p11_t179_g371981 + (cse_p115_t1_g372065 * cse_p11_t135_g371937) + (cse_p116_t1_g372069 * vy1 * (cse_p63_t1_g372040 + (cse_p375_t1_g372146 * vy1))) + (cse_p11_t143_g371945 * vy2) + (cse_p427_t1_g372166 * (cse_p11_t110_g371912 + cse_p11_t117_g371919 + n2)) + (-((cse_p117_t1_g372085 * cse_p11_t139_g371941))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
      double h_y1_z2 = (cse_p11_t65_g371867 + cse_p11_t65_g371867 + (-((cse_p63_t1_g372034 * ((cse_p115_t1_g372059 * cse_p63_t11_g372028 * vy1) + (cse_p116_t1_g372072 * cse_p11_t135_g371937) + (cse_p11_t144_g371946 * vy2) + (-((cse_p117_t1_g372085 * cse_p11_t140_g371942))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
      double cse_p375_t1_g372140 = (vx2 * vy1);
      double h_y1_x3 = (cse_p11_t71_g371873 + cse_p11_t71_g371873 + (-((cse_p63_t1_g372034 * (cse_p11_t163_g371965 + cse_p11_t172_g371974 + (cse_p117_t1_g372080 * cse_p11_t135_g371937) + (cse_p375_t1_g372140 * cse_p64_t1_g372045) + (-((cse_p375_t1_g372140 * cse_p64_t1_g372047))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 0, h_y1_x3);
      double cse_p115_t1_g372063 = (cse_p377_t1_g372158 * vy2);
      double cse_p168_t1_g372108 = (cse_p11_t134_g371936 * cse_p63_t1_g372033);
      double h_y1_y3 = (cse_p11_t80_g371882 + cse_p11_t80_g371882 + (-((cse_p63_t1_g372034 * (cse_p1_t5_invr371799 + (cse_p11_t135_g371937 * cse_p168_t1_g372108) + (cse_p167_t1_g372102 * cse_p64_t1_g372045) + (-((cse_p115_t1_g372063 * vy2))) + (-((cse_p117_t1_g372085 * vy1))) + (-((cse_p167_t1_g372102 * cse_p64_t1_g372047))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 1, h_y1_y3);
      double cse_p375_t1_g372142 = (vy1 * vz2);
      double h_y1_z3 = (cse_p11_t84_g371886 + cse_p11_t84_g371886 + (-((cse_p63_t1_g372034 * (cse_p11_t168_g371970 + cse_p11_t175_g371977 + (cse_p11_t136_g371938 * cse_p376_t1_g372157) + (cse_p375_t1_g372142 * cse_p64_t1_g372045) + (-((cse_p375_t1_g372142 * cse_p64_t1_g372047))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 2, h_y1_z3);
      double h_z1_z1 = (cse_p11_t88_g371890 + cse_p11_t88_g371890 + (-((cse_p63_t1_g372034 * (cse_p11_t158_g371960 + (cse_p117_t1_g372087 * ((cse_p116_t1_g372079 * vz1) + (cse_p375_t1_g372144 * vz2) + (cse_p64_t1_g372049 * vz1))) + (cse_p11_t137_g371939 * cse_p11_t137_g371939 * cse_p63_t1_g372033) + (-((cse_p64_t1_g372048 * vz1 * vz1))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
      double h_z1_x2 = (cse_p11_t52_g371854 + cse_p11_t52_g371854 + (-((cse_p63_t1_g372034 * ((cse_p115_t1_g372059 * cse_p63_t9_g372026 * vz1) + (cse_p115_t1_g372061 * cse_p11_t137_g371939) + (cse_p11_t142_g371944 * vz2) + (-((cse_p117_t1_g372086 * cse_p11_t138_g371940))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
      double h_z1_y2 = (cse_p11_t60_g371862 + cse_p11_t60_g371862 + (-((cse_p63_t1_g372034 * ((cse_p115_t1_g372059 * cse_p63_t10_g372027 * vz1) + (cse_p115_t1_g372065 * cse_p11_t137_g371939) + (cse_p11_t143_g371945 * vz2) + (-((cse_p117_t1_g372086 * cse_p11_t139_g371941))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
      double h_z1_z2 = (cse_p11_t67_g371869 + cse_p11_t67_g371869 + (-((cse_p63_t1_g372034 * (cse_p11_t179_g371981 + (cse_p116_t1_g372069 * vz1 * (cse_p63_t1_g372041 + (cse_p375_t1_g372146 * vz1))) + (cse_p116_t1_g372072 * cse_p11_t137_g371939) + (cse_p11_t144_g371946 * vz2) + (cse_p427_t1_g372166 * (cse_p11_t114_g371916 + cse_p11_t119_g371921 + n2)) + (-((cse_p117_t1_g372086 * cse_p11_t140_g371942))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
      double cse_p375_t1_g372141 = (vx2 * vz1);
      double h_z1_x3 = (cse_p11_t73_g371875 + cse_p11_t73_g371875 + (-((cse_p63_t1_g372034 * (cse_p11_t165_g371967 + cse_p11_t173_g371975 + (cse_p117_t1_g372080 * cse_p11_t137_g371939) + (cse_p375_t1_g372141 * cse_p64_t1_g372045) + (-((cse_p375_t1_g372141 * cse_p64_t1_g372047))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 0, h_z1_x3);
      double cse_p375_t1_g372143 = (vy2 * vz1);
      double h_z1_y3 = (cse_p11_t82_g371884 + cse_p11_t82_g371884 + (-((cse_p63_t1_g372034 * (cse_p11_t168_g371970 + cse_p11_t175_g371977 + (cse_p11_t137_g371939 * cse_p168_t1_g372108) + (cse_p375_t1_g372143 * cse_p64_t1_g372045) + (-((cse_p375_t1_g372143 * cse_p64_t1_g372047))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 1, h_z1_y3);
      double cse_p115_t1_g372066 = (cse_p377_t1_g372158 * vz2);
      double h_z1_z3 = (cse_p11_t87_g371889 + cse_p11_t87_g371889 + (-((cse_p63_t1_g372034 * (cse_p1_t5_invr371799 + (cse_p11_t136_g371938 * cse_p11_t137_g371939 * cse_p63_t1_g372033) + (cse_p167_t1_g372105 * cse_p64_t1_g372045) + (-((cse_p115_t1_g372066 * vz2))) + (-((cse_p117_t1_g372086 * vz1))) + (-((cse_p167_t1_g372105 * cse_p64_t1_g372047))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 2, h_z1_z3);
      double cse_p116_t1_g372074 = (-2.00000000000000000e+0 * cse_p64_t1_g372046);
      double cse_p375_t1_g372145 = (cse_p114_t1_g372057 * cse_p529_t1_sqrt372183);
      double cse_p376_t1_g372152 = (cse_p11_t121_g371923 * cse_p63_t1_g372035);
      double h_x2_x2 = (cse_p11_t44_g371846 + cse_p11_t44_g371846 + (-((cse_p63_t1_g372034 * (cse_p11_t18_g371820 + (2.00000000000000000e+0 * cse_p11_t138_g371940 * cse_p11_t142_g371944) + (cse_p115_t1_g372061 * cse_p11_t129_g371931) + (dot * ((cse_p11_t121_g371923 * cse_p376_t1_g372152) + (-((cse_p61_t8_invr2372012 * (cse_p11_t92_g371894 + cse_p478_t1_g372173 + (cse_p115_t1_g372067 * cse_p375_t1_g372145 * vx1) + (cse_p116_t1_g372074 * cse_p167_t1_g372099) + (cse_p116_t1_g372077 * (cse_p11_t105_g371907 + cse_p11_t115_g371917 + cse_p11_t115_g371917 + n2)))))))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
      double cse_p115_t1_g372068 = (cse_p115_t1_g372067 * cse_p529_t1_sqrt372183 * cse_p63_t13_g372030);
      double h_x2_y2 = (cse_p11_t45_g371847 + cse_p11_t45_g371847 + (-((cse_p63_t1_g372034 * ((cse_p115_t1_g372061 * cse_p11_t130_g371932) + (cse_p11_t138_g371940 * cse_p11_t143_g371945) + (cse_p11_t139_g371941 * cse_p11_t142_g371944) + (dot * ((cse_p11_t122_g371924 * cse_p376_t1_g372152) + (-((cse_p61_t8_invr2372012 * (cse_p11_t94_g371896 + (cse_p115_t1_g372068 * vy1) + (cse_p116_t1_g372074 * cse_p167_t1_g372100) + (cse_p11_t1_g372004 * (cse_p375_t1_g372138 + cse_p375_t1_g372140)))))))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
      double h_x2_z2 = (cse_p11_t46_g371848 + cse_p11_t46_g371848 + (-((cse_p63_t1_g372034 * ((cse_p115_t1_g372061 * cse_p11_t131_g371933) + (cse_p11_t138_g371940 * cse_p11_t144_g371946) + (cse_p11_t140_g371942 * cse_p11_t142_g371944) + (dot * ((cse_p11_t123_g371925 * cse_p376_t1_g372152) + (-((cse_p61_t8_invr2372012 * (cse_p11_t96_g371898 + (cse_p115_t1_g372068 * vz1) + (cse_p116_t1_g372074 * cse_p167_t1_g372101) + (cse_p11_t1_g372004 * (cse_p375_t1_g372139 + cse_p375_t1_g372141)))))))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
      double cse_p376_t1_g372153 = (cse_p11_t121_g371923 * cse_p63_t1_g372037);
      double h_x2_x3 = (cse_p11_t47_g371849 + cse_p11_t47_g371849 + (-((cse_p63_t1_g372034 * (cse_p11_t179_g371981 + (cse_p115_t1_g372061 * cse_p11_t132_g371934) + (cse_p11_t142_g371944 * vx1) + (dot * ((cse_p376_t1_g372153 * vx2) + (-((cse_p61_t8_invr2372012 * (cse_p11_t177_g371979 + (cse_p116_t1_g372076 * cse_p167_t1_g372099) + (-((cse_p11_t1_g372004 * cse_p167_t1_g372098))) + (-(cse_p11_t92_g371894)))))))) + (-((cse_p115_t1_g372058 * cse_p11_t138_g371940))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 0, h_x2_x3);
      double h_x2_y3 = (cse_p11_t49_g371851 + cse_p11_t49_g371851 + (-((cse_p63_t1_g372034 * ((cse_p115_t1_g372061 * cse_p11_t134_g371936) + (cse_p11_t142_g371944 * vy1) + (dot * ((cse_p376_t1_g372153 * vy2) + (-((cse_p61_t8_invr2372012 * (cse_p11_t164_g371966 + cse_p11_t20_g371822 + (-((cse_p11_t1_g372004 * cse_p375_t1_g372138))))))))) + (-((cse_p115_t1_g372063 * cse_p11_t138_g371940))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 1, h_x2_y3);
      double h_x2_z3 = (cse_p11_t51_g371853 + cse_p11_t51_g371853 + (-((cse_p63_t1_g372034 * ((cse_p115_t1_g372061 * cse_p11_t136_g371938) + (cse_p11_t142_g371944 * vz1) + (dot * ((cse_p376_t1_g372153 * vz2) + (-((cse_p61_t8_invr2372012 * (cse_p11_t166_g371968 + cse_p11_t21_g371823 + (-((cse_p11_t1_g372004 * cse_p375_t1_g372139))))))))) + (-((cse_p115_t1_g372066 * cse_p11_t138_g371940))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 2, h_x2_z3);
      double h_y2_y2 = (cse_p11_t53_g371855 + cse_p11_t53_g371855 + (-((cse_p63_t1_g372034 * (cse_p11_t18_g371820 + (2.00000000000000000e+0 * cse_p11_t139_g371941 * cse_p11_t143_g371945) + (cse_p115_t1_g372065 * cse_p11_t130_g371932) + (dot * ((cse_p11_t122_g371924 * cse_p11_t122_g371924 * cse_p63_t1_g372035) + (-((cse_p61_t8_invr2372012 * (cse_p11_t98_g371900 + cse_p478_t1_g372173 + (cse_p116_t1_g372073 * cse_p375_t1_g372145 * vy1) + (cse_p116_t1_g372074 * cse_p167_t1_g372103) + (cse_p116_t1_g372077 * (cse_p11_t110_g371912 + cse_p11_t117_g371919 + cse_p11_t117_g371919 + n2)))))))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
      double h_y2_z2 = (cse_p11_t54_g371856 + cse_p11_t54_g371856 + (-((cse_p63_t1_g372034 * ((cse_p115_t1_g372065 * cse_p11_t131_g371933) + (cse_p11_t139_g371941 * cse_p11_t144_g371946) + (cse_p11_t140_g371942 * cse_p11_t143_g371945) + (dot * ((cse_p11_t122_g371924 * cse_p11_t123_g371925 * cse_p63_t1_g372035) + (-((cse_p61_t8_invr2372012 * (cse_p11_t100_g371902 + (cse_p116_t1_g372073 * cse_p529_t1_sqrt372183 * cse_p63_t13_g372030 * vz1) + (cse_p116_t1_g372074 * cse_p167_t1_g372104) + (cse_p11_t1_g372004 * (cse_p375_t1_g372142 + cse_p375_t1_g372143)))))))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
      double cse_p376_t1_g372154 = (cse_p11_t122_g371924 * cse_p63_t1_g372037);
      double h_y2_x3 = (cse_p11_t55_g371857 + cse_p11_t55_g371857 + (-((cse_p63_t1_g372034 * ((cse_p115_t1_g372065 * cse_p11_t132_g371934) + (cse_p11_t143_g371945 * vx1) + (dot * ((cse_p376_t1_g372154 * vx2) + (-((cse_p61_t8_invr2372012 * (cse_p11_t164_g371966 + cse_p11_t20_g371822 + (-((cse_p11_t1_g372004 * cse_p375_t1_g372140))))))))) + (-((cse_p115_t1_g372058 * cse_p11_t139_g371941))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 0, h_y2_x3);
      double h_y2_y3 = (cse_p11_t57_g371859 + cse_p11_t57_g371859 + (-((cse_p63_t1_g372034 * (cse_p11_t179_g371981 + (cse_p115_t1_g372065 * cse_p11_t134_g371936) + (cse_p11_t143_g371945 * vy1) + (dot * ((cse_p376_t1_g372154 * vy2) + (-((cse_p61_t8_invr2372012 * (cse_p11_t177_g371979 + (cse_p116_t1_g372076 * cse_p167_t1_g372103) + (-((cse_p11_t1_g372004 * cse_p167_t1_g372102))) + (-(cse_p11_t98_g371900)))))))) + (-((cse_p115_t1_g372063 * cse_p11_t139_g371941))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 1, h_y2_y3);
      double h_y2_z3 = (cse_p11_t59_g371861 + cse_p11_t59_g371861 + (-((cse_p63_t1_g372034 * ((cse_p115_t1_g372065 * cse_p11_t136_g371938) + (cse_p11_t143_g371945 * vz1) + (dot * ((cse_p376_t1_g372154 * vz2) + (-((cse_p61_t8_invr2372012 * (cse_p11_t169_g371971 + cse_p11_t22_g371824 + (-((cse_p11_t1_g372004 * cse_p375_t1_g372142))))))))) + (-((cse_p115_t1_g372066 * cse_p11_t139_g371941))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 2, h_y2_z3);
      double h_z2_z2 = (cse_p11_t61_g371863 + cse_p11_t61_g371863 + (-((cse_p63_t1_g372034 * (cse_p11_t18_g371820 + (2.00000000000000000e+0 * cse_p11_t140_g371942 * cse_p11_t144_g371946) + (cse_p116_t1_g372072 * cse_p11_t131_g371933) + (dot * ((cse_p11_t123_g371925 * cse_p11_t123_g371925 * cse_p63_t1_g372035) + (-((cse_p61_t8_invr2372012 * (cse_p11_t102_g371904 + cse_p478_t1_g372173 + (cse_p116_t1_g372074 * cse_p167_t1_g372106) + (cse_p116_t1_g372077 * (cse_p11_t114_g371916 + cse_p11_t119_g371921 + cse_p11_t119_g371921 + n2)) + (cse_p117_t1_g372087 * cse_p375_t1_g372145 * vz1))))))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      double cse_p376_t1_g372155 = (cse_p11_t123_g371925 * cse_p63_t1_g372037);
      double h_z2_x3 = (cse_p11_t62_g371864 + cse_p11_t62_g371864 + (-((cse_p63_t1_g372034 * ((cse_p116_t1_g372072 * cse_p11_t132_g371934) + (cse_p11_t144_g371946 * vx1) + (dot * ((cse_p376_t1_g372155 * vx2) + (-((cse_p61_t8_invr2372012 * (cse_p11_t166_g371968 + cse_p11_t21_g371823 + (-((cse_p11_t1_g372004 * cse_p375_t1_g372141))))))))) + (-((cse_p115_t1_g372058 * cse_p11_t140_g371942))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 0, h_z2_x3);
      double h_z2_y3 = (cse_p11_t64_g371866 + cse_p11_t64_g371866 + (-((cse_p63_t1_g372034 * ((cse_p116_t1_g372072 * cse_p11_t134_g371936) + (cse_p11_t144_g371946 * vy1) + (dot * ((cse_p376_t1_g372155 * vy2) + (-((cse_p61_t8_invr2372012 * (cse_p11_t169_g371971 + cse_p11_t22_g371824 + (-((cse_p11_t1_g372004 * cse_p375_t1_g372143))))))))) + (-((cse_p115_t1_g372063 * cse_p11_t140_g371942))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 1, h_z2_y3);
      double h_z2_z3 = (cse_p11_t66_g371868 + cse_p11_t66_g371868 + (-((cse_p63_t1_g372034 * (cse_p11_t179_g371981 + (cse_p116_t1_g372072 * cse_p11_t136_g371938) + (cse_p11_t144_g371946 * vz1) + (dot * ((cse_p376_t1_g372155 * vz2) + (-((cse_p61_t8_invr2372012 * (cse_p11_t177_g371979 + (cse_p116_t1_g372076 * cse_p167_t1_g372106) + (-((cse_p11_t1_g372004 * cse_p167_t1_g372105))) + (-(cse_p11_t102_g371904)))))))) + (-((cse_p115_t1_g372066 * cse_p11_t140_g371942))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 2, h_z2_z3);
      double cse_p117_t1_g372083 = (cse_p63_t1_g372031 * dot);
      double cse_p375_t1_g372147 = (cse_p115_t1_g372058 * cse_p64_t1_g372044);
      double cse_p376_t1_g372151 = (cse_p117_t1_g372083 * cse_p63_t12_g372029);
      double h_x3_x3 = (cse_p11_t68_g371870 + cse_p11_t68_g371870 + (-((cse_p63_t1_g372034 * (cse_p11_t154_g371956 + (-2.00000000000000000e+0 * cse_p115_t1_g372058 * vx1) + (cse_p117_t1_g372080 * cse_p11_t132_g371934) + (cse_p167_t1_g372099 * cse_p376_t1_g372151) + (-((cse_p375_t1_g372147 * vx2))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 0, h_x3_x3);
      double cse_p376_t1_g372150 = (cse_p117_t1_g372083 * cse_p11_t128_g371930);
      double h_x3_y3 = (cse_p11_t70_g371872 + cse_p11_t70_g371872 + (-((cse_p63_t1_g372034 * ((cse_p117_t1_g372080 * cse_p11_t134_g371936) + (cse_p167_t1_g372100 * cse_p376_t1_g372150) + (-((cse_p115_t1_g372058 * vy1))) + (-((cse_p115_t1_g372063 * vx1))) + (-((cse_p375_t1_g372147 * vy2))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 1, h_x3_y3);
      double h_x3_z3 = (cse_p11_t72_g371874 + cse_p11_t72_g371874 + (-((cse_p63_t1_g372034 * ((cse_p117_t1_g372080 * cse_p11_t136_g371938) + (cse_p167_t1_g372101 * cse_p376_t1_g372150) + (-((cse_p115_t1_g372058 * vz1))) + (-((cse_p115_t1_g372066 * vx1))) + (-((cse_p375_t1_g372147 * vz2))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 2, h_x3_z3);
      double h_y3_y3 = (cse_p11_t79_g371881 + cse_p11_t79_g371881 + (-((cse_p63_t1_g372034 * (cse_p11_t154_g371956 + (-2.00000000000000000e+0 * cse_p115_t1_g372063 * vy1) + (cse_p11_t134_g371936 * cse_p168_t1_g372108) + (cse_p167_t1_g372103 * cse_p376_t1_g372151) + (-((cse_p115_t1_g372063 * cse_p64_t1_g372044 * vy2))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 1, h_y3_y3);
      double h_y3_z3 = (cse_p11_t81_g371883 + cse_p11_t81_g371883 + (-((cse_p63_t1_g372034 * ((cse_p11_t136_g371938 * cse_p168_t1_g372108) + (cse_p167_t1_g372104 * cse_p376_t1_g372150) + (-((cse_p115_t1_g372063 * cse_p64_t1_g372044 * vz2))) + (-((cse_p115_t1_g372063 * vz1))) + (-((cse_p115_t1_g372066 * vy1))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 2, h_y3_z3);
      double h_z3_z3 = (cse_p11_t86_g371888 + cse_p11_t86_g371888 + (-((cse_p63_t1_g372034 * (cse_p11_t154_g371956 + (-2.00000000000000000e+0 * cse_p115_t1_g372066 * vz1) + (cse_p11_t136_g371938 * cse_p11_t136_g371938 * cse_p63_t1_g372033) + (cse_p167_t1_g372106 * cse_p376_t1_g372151) + (-((cse_p115_t1_g372066 * cse_p64_t1_g372044 * vz2))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x3, 2, h_z3_z3);
    }
  }
  return energy_added;
}
void hessian_fd(const angle_term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 9;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  const double invh2 = 1.0/((h*h));
  double e0 = 0.0;
  energy(term, position, &e0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[term.i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x1 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x1, 0, d);
  }
  {
    double saved = position[term.i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x1 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x1, 1, d);
  }
  {
    double saved = position[term.i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x1 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x1, 2, d);
  }
  {
    double saved = position[term.i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x2 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x2, 0, d);
  }
  {
    double saved = position[term.i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x2 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x2, 1, d);
  }
  {
    double saved = position[term.i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x2 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x2, 2, d);
  }
  {
    double saved = position[term.i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x3 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x3 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x3 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x3, 0, d);
  }
  {
    double saved = position[term.i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x3 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x3 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x3 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x3, 1, d);
  }
  {
    double saved = position[term.i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x3 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x3 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x3 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x3, 2, d);
  }
  {
    double saved = position[term.i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x1 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x1 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 0, term.i3x1, 0, hval);
  }
  {
    double saved = position[term.i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x1 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x1 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 1, term.i3x1, 1, hval);
  }
  {
    double saved = position[term.i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x1 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x1 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 2, term.i3x1, 2, hval);
  }
  {
    double saved = position[term.i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x2 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x2 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 0, term.i3x2, 0, hval);
  }
  {
    double saved = position[term.i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x2 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x2 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 1, term.i3x2, 1, hval);
  }
  {
    double saved = position[term.i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x2 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x2 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 2, term.i3x2, 2, hval);
  }
  {
    double saved = position[term.i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x3 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x3 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x3 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x3, 0, hval);
  }
  {
    double saved = position[term.i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x3 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x3 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x3 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x3, 1, hval);
  }
  {
    double saved = position[term.i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x3 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x3 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x3 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x3, 2, hval);
  }
  {
    double saved_i = position[term.i3x1 + 1];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x1 + 1] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x1 + 1] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x1 + 1] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 1, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x1 + 2];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x1 + 2] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x1 + 2] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x1 + 2] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 2, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x1 + 2];
    double saved_j = position[term.i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x1 + 2] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x1 + 2] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x1 + 2] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 2, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x2 + 0];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x2 + 0] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 0] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x2 + 0] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 0, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x2 + 0];
    double saved_j = position[term.i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x2 + 0] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 0] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x2 + 0] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 0, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x2 + 0];
    double saved_j = position[term.i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x2 + 0] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 0] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x2 + 0] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 0, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x2 + 1];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x2 + 1] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 1] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x2 + 1] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 1, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x2 + 1];
    double saved_j = position[term.i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x2 + 1] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 1] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x2 + 1] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 1, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x2 + 1];
    double saved_j = position[term.i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x2 + 1] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 1] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x2 + 1] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 1, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x2 + 1];
    double saved_j = position[term.i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x2 + 1] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 1] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x2 + 1] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 1, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x2 + 2];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x2 + 2] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x2 + 2] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 2, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x2 + 2];
    double saved_j = position[term.i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x2 + 2] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x2 + 2] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 2, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x2 + 2];
    double saved_j = position[term.i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x2 + 2] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x2 + 2] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 2, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x2 + 2];
    double saved_j = position[term.i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x2 + 2] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x2 + 2] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 2, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x2 + 2];
    double saved_j = position[term.i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x2 + 2] = saved_i + h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x2 + 2] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 2, term.i3x2, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 0];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 0] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 0] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 0] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 0];
    double saved_j = position[term.i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 0] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 0] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 0] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 0];
    double saved_j = position[term.i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 0] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 0] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 0] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x3 + 0];
    double saved_j = position[term.i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 0] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 0] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 0] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 0];
    double saved_j = position[term.i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 0] = saved_i + h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 0] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 0] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x2, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 0];
    double saved_j = position[term.i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 0] = saved_i + h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 0] = saved_i - h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 0] = saved_i; position[term.i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x2, 2, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x2, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x2, 2, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x3 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x3 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x3, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x2, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x2, 2, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x3 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x3 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x3, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x3 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x3 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x3 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x3 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x3, 1, hval);
  }
}

};
