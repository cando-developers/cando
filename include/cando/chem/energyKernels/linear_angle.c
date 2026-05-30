#include "linear_angle.h"

template <typename HESSIAN>
struct Linear_Angle {
  static constexpr size_t PositionSize = 9;
  static std::string description() { return "mathkernel-linear_angle"; };
double energy(const linear_angle_term& term, double* position, double* energy_accumulate) {
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
      double cse_p11_t1_g558934 = (-(x2));
      double cse_p11_t2_g558935 = (-(y2));
      double cse_p11_t3_g558936 = (-(z2));
      double vx1 = (cse_p11_t1_g558934 + x1);
      double vy1 = (cse_p11_t2_g558935 + y1);
      double vz1 = (cse_p11_t3_g558936 + z1);
      double vx2 = (cse_p11_t1_g558934 + x3);
      double vy2 = (cse_p11_t2_g558935 + y3);
      double vz2 = (cse_p11_t3_g558936 + z3);
      double dot = ((vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2));
      double n1_sq = ((vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1));
      double n2_sq = ((vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2));
      double n1 = sqrt(n1_sq);
      double n2 = sqrt(n2_sq);
      double cos_theta = (dot * (1.0 / ((n1 * n2))));
      double energy = (term.kt * (1.00000000000000000e+0 + cos_theta));
      energy_added += energy;
      *energy_accumulate += energy;
    }
  }
  return energy_added;
}
void energy_fd(const linear_angle_term& term, double* position, double* energy_accumulate)
{
  energy(term, position, energy_accumulate);
}

double gradient(const linear_angle_term& term, double* position, double* energy_accumulate, double* force) {
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
      double cse_p11_t1_g558943 = (-(x2));
      double cse_p11_t2_g558944 = (-(y2));
      double cse_p11_t3_g558945 = (-(z2));
      double vx1 = (cse_p11_t1_g558943 + x1);
      double vy1 = (cse_p11_t2_g558944 + y1);
      double vz1 = (cse_p11_t3_g558945 + z1);
      double vx2 = (cse_p11_t1_g558943 + x3);
      double vy2 = (cse_p11_t2_g558944 + y3);
      double vz2 = (cse_p11_t3_g558945 + z3);
      double dot = ((vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2));
      double n1_sq = ((vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1));
      double n2_sq = ((vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2));
      double cse_p1_t1_sqrt558937 = sqrt(n1_sq);
      double cse_p1_t2_invsqrt558938 = (1.0 / (cse_p1_t1_sqrt558937));
      double n1 = cse_p1_t1_sqrt558937;
      double cse_p1_t3_sqrt558939 = sqrt(n2_sq);
      double cse_p1_t4_invsqrt558940 = (1.0 / (cse_p1_t3_sqrt558939));
      double n2 = cse_p1_t3_sqrt558939;
      double cse_p1_t5_invr558941 = (1.0 / ((n1 * n2)));
      double cse_p1_t6_invr2558942 = (cse_p1_t5_invr558941 * cse_p1_t5_invr558941);
      double cse_p113_t1_g558951 = (cse_p1_t6_invr2558942 * dot);
      double cos_theta = (cse_p1_t5_invr558941 * dot);
      double energy = (term.kt * (1.00000000000000000e+0 + cos_theta));
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t1_g558946 = (cse_p1_t2_invsqrt558938 * n2);
      double cse_p11_t1_g558949 = (cse_p113_t1_g558951 * cse_p11_t1_g558946);
      double g_x1 = (term.kt * ((cse_p1_t5_invr558941 * vx2) + (-((cse_p11_t1_g558949 * vx1)))));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = (term.kt * ((cse_p1_t5_invr558941 * vy2) + (-((cse_p11_t1_g558949 * vy1)))));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = (term.kt * ((cse_p1_t5_invr558941 * vz2) + (-((cse_p11_t1_g558949 * vz1)))));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p11_t1_g558947 = (cse_p1_t4_invsqrt558940 * n1);
      double g_x2 = (term.kt * ((cse_p1_t5_invr558941 * ((-(vx1)) + (-(vx2)))) + (-((cse_p113_t1_g558951 * ((-((cse_p11_t1_g558946 * vx1))) + (-((cse_p11_t1_g558947 * vx2)))))))));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (term.kt * ((cse_p1_t5_invr558941 * ((-(vy1)) + (-(vy2)))) + (-((cse_p113_t1_g558951 * ((-((cse_p11_t1_g558946 * vy1))) + (-((cse_p11_t1_g558947 * vy2)))))))));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (term.kt * ((cse_p1_t5_invr558941 * ((-(vz1)) + (-(vz2)))) + (-((cse_p113_t1_g558951 * ((-((cse_p11_t1_g558946 * vz1))) + (-((cse_p11_t1_g558947 * vz2)))))))));
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double cse_p11_t1_g558950 = (cse_p113_t1_g558951 * cse_p11_t1_g558947);
      double g_x3 = (term.kt * ((cse_p1_t5_invr558941 * vx1) + (-((cse_p11_t1_g558950 * vx2)))));
      KernelGradientAcc(term.i3x3, 0, g_x3);
      double g_y3 = (term.kt * ((cse_p1_t5_invr558941 * vy1) + (-((cse_p11_t1_g558950 * vy2)))));
      KernelGradientAcc(term.i3x3, 1, g_y3);
      double g_z3 = (term.kt * ((cse_p1_t5_invr558941 * vz1) + (-((cse_p11_t1_g558950 * vz2)))));
      KernelGradientAcc(term.i3x3, 2, g_z3);
    }
  }
  return energy_added;
}
void gradient_fd(const linear_angle_term& term, double* position, double* energy_accumulate, double* force)
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

double hessian(const linear_angle_term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
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
      double cse_p11_t88_g559045 = (-(x2));
      double cse_p11_t89_g559046 = (-(y2));
      double cse_p11_t90_g559047 = (-(z2));
      double vx1 = (cse_p11_t88_g559045 + x1);
      double cse_p11_t82_g559039 = (-(vx1));
      double vy1 = (cse_p11_t89_g559046 + y1);
      double cse_p11_t84_g559041 = (-(vy1));
      double vz1 = (cse_p11_t90_g559047 + z1);
      double cse_p11_t86_g559043 = (-(vz1));
      double vx2 = (cse_p11_t88_g559045 + x3);
      double cse_p11_t57_g559014 = (cse_p11_t82_g559039 + (-(vx2)));
      double cse_p11_t83_g559040 = (-(vx2));
      double vy2 = (cse_p11_t89_g559046 + y3);
      double cse_p11_t58_g559015 = (cse_p11_t84_g559041 + (-(vy2)));
      double cse_p11_t85_g559042 = (-(vy2));
      double vz2 = (cse_p11_t90_g559047 + z3);
      double cse_p11_t59_g559016 = (cse_p11_t86_g559043 + (-(vz2)));
      double cse_p11_t87_g559044 = (-(vz2));
      double cse_p64_t1_g559107 = (vx1 * vx2);
      double cse_p64_t1_g559110 = (vy1 * vy2);
      double cse_p64_t1_g559113 = (vz1 * vz2);
      double dot = (cse_p64_t1_g559107 + cse_p64_t1_g559110 + cse_p64_t1_g559113);
      double cse_p64_t1_g559106 = (vx1 * vx1);
      double cse_p64_t1_g559109 = (vy1 * vy1);
      double cse_p64_t1_g559112 = (vz1 * vz1);
      double n1_sq = (cse_p64_t1_g559106 + cse_p64_t1_g559109 + cse_p64_t1_g559112);
      double cse_p529_t1_sqrt559200 = sqrt(n1_sq);
      double cse_p529_t2_invr559201 = (1.0 / (n1_sq));
      double cse_p61_t3_invr2559076 = (cse_p529_t2_invr559201 * cse_p529_t2_invr559201);
      double cse_p64_t1_g559108 = (vx2 * vx2);
      double cse_p64_t1_g559111 = (vy2 * vy2);
      double cse_p65_t1_g559114 = (vz2 * vz2);
      double n2_sq = (cse_p64_t1_g559108 + cse_p64_t1_g559111 + cse_p65_t1_g559114);
      double cse_p529_t3_sqrt559202 = sqrt(n2_sq);
      double cse_p529_t4_invr559203 = (1.0 / (n2_sq));
      double cse_p61_t6_invr2559079 = (cse_p529_t4_invr559203 * cse_p529_t4_invr559203);
      double cse_p11_t1_g559054 = (cse_p529_t3_sqrt559202 * cse_p529_t4_invr559203);
      double cse_p11_t43_g559000 = (cse_p11_t1_g559054 * cse_p64_t1_g559107);
      double cse_p116_t1_g559144 = (cse_p11_t1_g559054 * cse_p529_t2_invr559201);
      double cse_p11_t44_g559001 = (cse_p116_t1_g559144 * vx2);
      double cse_p11_t45_g559002 = (cse_p11_t1_g559054 * cse_p64_t1_g559110);
      double cse_p11_t46_g559003 = (cse_p116_t1_g559144 * vy2);
      double cse_p11_t47_g559004 = (cse_p11_t1_g559054 * cse_p64_t1_g559113);
      double cse_p11_t48_g559005 = (cse_p116_t1_g559144 * vz2);
      double cse_p1_t2_invsqrt558953 = (1.0 / (cse_p529_t1_sqrt559200));
      double n1 = cse_p529_t1_sqrt559200;
      double cse_p62_t1_g559083 = (cse_p11_t1_g559054 * n1);
      double cse_p63_t1_g559094 = (cse_p529_t3_sqrt559202 * n1);
      double cse_p63_t1_g559099 = (cse_p61_t6_invr2559079 * cse_p63_t1_g559094);
      double cse_p65_t1_g559116 = (2.00000000000000000e+0 * cse_p63_t1_g559099);
      double cse_p65_t1_g559123 = (vx2 * vy2);
      double cse_p11_t11_g558968 = (cse_p65_t1_g559116 * cse_p65_t1_g559123);
      double cse_p66_t1_g559125 = (vx2 * vz2);
      double cse_p11_t12_g558969 = (cse_p65_t1_g559116 * cse_p66_t1_g559125);
      double cse_p66_t1_g559128 = (vy2 * vz2);
      double cse_p11_t13_g558970 = (cse_p65_t1_g559116 * cse_p66_t1_g559128);
      double cse_p11_t19_g558976 = (cse_p62_t1_g559083 * vx2);
      double cse_p64_t1_g559105 = (cse_p529_t4_invr559203 * cse_p62_t1_g559083);
      double cse_p11_t20_g558977 = (cse_p64_t1_g559105 * cse_p64_t1_g559108);
      double cse_p11_t22_g558979 = (cse_p64_t1_g559105 * cse_p65_t1_g559123);
      double cse_p11_t24_g558981 = (cse_p64_t1_g559105 * cse_p66_t1_g559125);
      double cse_p11_t25_g558982 = (cse_p62_t1_g559083 * vy2);
      double cse_p11_t26_g558983 = (cse_p64_t1_g559105 * cse_p64_t1_g559111);
      double cse_p11_t28_g558985 = (cse_p64_t1_g559105 * cse_p66_t1_g559128);
      double cse_p11_t29_g558986 = (cse_p62_t1_g559083 * vz2);
      double cse_p11_t30_g558987 = (cse_p64_t1_g559105 * cse_p65_t1_g559114);
      double cse_p11_t65_g559022 = (-(cse_p11_t19_g558976));
      double cse_p11_t67_g559024 = (-(cse_p11_t22_g558979));
      double cse_p11_t69_g559026 = (-(cse_p11_t24_g558981));
      double cse_p11_t70_g559027 = (-(cse_p11_t25_g558982));
      double cse_p11_t72_g559029 = (-(cse_p11_t28_g558985));
      double cse_p11_t73_g559030 = (-(cse_p11_t29_g558986));
      double cse_p11_t80_g559037 = (-(cse_p62_t1_g559083));
      double cse_p1_t4_invsqrt558955 = (1.0 / (cse_p529_t3_sqrt559202));
      double n2 = cse_p529_t3_sqrt559202;
      double cse_p63_t1_g559097 = (cse_p529_t1_sqrt559200 * n2);
      double cse_p63_t1_g559100 = (-2.00000000000000000e+0 * cse_p61_t3_invr2559076);
      double cse_p114_t1_g559133 = (cse_p63_t1_g559097 * cse_p63_t1_g559100);
      double cse_p11_t1_g559057 = (cse_p529_t2_invr559201 * cse_p63_t1_g559097);
      double cse_p322_t1_g559180 = (cse_p11_t1_g559057 * vz1);
      double cse_p270_t1_g559175 = (cse_p11_t1_g559057 * vy1);
      double cse_p218_t1_g559170 = (cse_p11_t1_g559057 * vx1);
      double cse_p61_t7_invr559080 = (1.0 / ((n1 * n2)));
      double cse_p61_t8_invr2559081 = (cse_p61_t7_invr559080 * cse_p61_t7_invr559080);
      double cse_p61_t9_invr3559082 = (cse_p61_t7_invr559080 * cse_p61_t8_invr2559081);
      double cse_p11_t1_g558958 = (cse_p61_t8_invr2559081 * (cse_p11_t65_g559022 + (-(cse_p218_t1_g559170))));
      double cse_p11_t2_g558959 = (cse_p61_t8_invr2559081 * (cse_p11_t70_g559027 + (-(cse_p270_t1_g559175))));
      double cse_p11_t3_g558960 = (cse_p61_t8_invr2559081 * (cse_p11_t73_g559030 + (-(cse_p322_t1_g559180))));
      double cse_p117_t1_g559156 = (cse_p63_t1_g559100 * n2);
      double cse_p11_t4_g558961 = (cse_p117_t1_g559156 * vx1);
      double cse_p11_t5_g558962 = (cse_p117_t1_g559156 * vy1);
      double cse_p11_t6_g558963 = (cse_p117_t1_g559156 * vz1);
      double cse_p117_t1_g559160 = (2.00000000000000000e+0 * cse_p61_t8_invr2559081);
      double cse_p11_t7_g558964 = (cse_p117_t1_g559160 * cse_p61_t3_invr2559076);
      double cse_p11_t8_g558965 = (cse_p117_t1_g559160 * cse_p61_t6_invr2559079);
      double cse_p63_t1_g559093 = (2.00000000000000000e+0 * cse_p61_t9_invr3559082);
      double cse_p11_t10_g558967 = (cse_p63_t1_g559093 * cse_p64_t1_g559105);
      double cse_p11_t14_g558971 = (cse_p11_t1_g559057 * cse_p529_t2_invr559201 * cse_p63_t1_g559093);
      double cse_p115_t1_g559135 = (cse_p61_t8_invr2559081 * cse_p62_t1_g559083);
      double cse_p11_t15_g558972 = (cse_p115_t1_g559135 * dot);
      double cse_p115_t1_g559134 = (cse_p11_t1_g559057 * cse_p61_t8_invr2559081);
      double cse_p11_t16_g558973 = (cse_p115_t1_g559134 * dot);
      double cse_p11_t18_g558975 = (n1 * n2);
      double cse_p11_t21_g558978 = (cse_p115_t1_g559135 * cse_p65_t1_g559123);
      double cse_p11_t23_g558980 = (cse_p115_t1_g559135 * cse_p66_t1_g559125);
      double cse_p11_t27_g558984 = (cse_p115_t1_g559135 * cse_p66_t1_g559128);
      double cse_p65_t1_g559119 = (cse_p529_t2_invr559201 * n2);
      double cse_p116_t1_g559145 = (cse_p529_t2_invr559201 * cse_p65_t1_g559119);
      double cse_p11_t31_g558988 = (cse_p116_t1_g559145 * vx1);
      double cse_p11_t33_g558990 = (cse_p64_t1_g559106 * cse_p65_t1_g559119);
      double cse_p115_t1_g559136 = (vx1 * vy1);
      double cse_p11_t34_g558991 = (cse_p115_t1_g559134 * cse_p115_t1_g559136);
      double cse_p115_t1_g559137 = (vx1 * vz1);
      double cse_p11_t35_g558992 = (cse_p115_t1_g559134 * cse_p115_t1_g559137);
      double cse_p11_t36_g558993 = (cse_p116_t1_g559145 * vy1);
      double cse_p11_t38_g558995 = (cse_p64_t1_g559109 * cse_p65_t1_g559119);
      double cse_p115_t1_g559138 = (vy1 * vz1);
      double cse_p11_t39_g558996 = (cse_p115_t1_g559134 * cse_p115_t1_g559138);
      double cse_p11_t40_g558997 = (cse_p116_t1_g559145 * vz1);
      double cse_p11_t42_g558999 = (cse_p64_t1_g559112 * cse_p65_t1_g559119);
      double cse_p11_t49_g559006 = (cse_p11_t65_g559022 + (-(cse_p218_t1_g559170)));
      double cse_p11_t50_g559007 = (cse_p11_t70_g559027 + (-(cse_p270_t1_g559175)));
      double cse_p11_t51_g559008 = (cse_p11_t73_g559030 + (-(cse_p322_t1_g559180)));
      double cse_p11_t52_g559009 = (cse_p11_t31_g558988 + cse_p11_t44_g559001 + cse_p11_t4_g558961);
      double cse_p11_t53_g559010 = (cse_p11_t36_g558993 + cse_p11_t46_g559003 + cse_p11_t5_g558962);
      double cse_p11_t54_g559011 = (cse_p11_t40_g558997 + cse_p11_t48_g559005 + cse_p11_t6_g558963);
      double cse_p11_t55_g559012 = (cse_p11_t14_g558971 + cse_p11_t7_g558964);
      double cse_p11_t56_g559013 = (cse_p11_t10_g558967 + cse_p11_t8_g558965);
      double cse_p11_t60_g559017 = (-(cse_p11_t1_g558958));
      double cse_p11_t61_g559018 = (-(cse_p11_t2_g558959));
      double cse_p11_t62_g559019 = (-(cse_p11_t3_g558960));
      double cse_p11_t63_g559020 = (-(cse_p11_t15_g558972));
      double cse_p11_t64_g559021 = (-(cse_p11_t16_g558973));
      double cse_p11_t66_g559023 = (-(cse_p11_t21_g558978));
      double cse_p11_t68_g559025 = (-(cse_p11_t23_g558980));
      double cse_p11_t71_g559028 = (-(cse_p11_t27_g558984));
      double cse_p11_t74_g559031 = (-(cse_p218_t1_g559170));
      double cse_p11_t75_g559032 = (-(cse_p11_t34_g558991));
      double cse_p11_t76_g559033 = (-(cse_p11_t35_g558992));
      double cse_p11_t77_g559034 = (-(cse_p270_t1_g559175));
      double cse_p11_t78_g559035 = (-(cse_p11_t39_g558996));
      double cse_p11_t79_g559036 = (-(cse_p322_t1_g559180));
      double cse_p1_t5_invr558956 = (1.0 / (cse_p11_t18_g558975));
      double cse_p11_t9_g558966 = (2.00000000000000000e+0 * cse_p1_t5_invr558956);
      double cse_p11_t81_g559038 = (-(cse_p1_t5_invr558956));
      double cse_p1_t6_invr2558957 = (cse_p1_t5_invr558956 * cse_p1_t5_invr558956);
      double cos_theta = (cse_p1_t5_invr558956 * dot);
      double energy = (term.kt * (1.00000000000000000e+0 + cos_theta));
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p64_t1_g559104 = (cse_p1_t6_invr2558957 * dot);
      double cse_p65_t1_g559117 = (cse_p1_t2_invsqrt558953 * n2);
      double cse_p117_t1_g559157 = (cse_p64_t1_g559104 * cse_p65_t1_g559117);
      double g_x1 = (term.kt * ((cse_p1_t5_invr558956 * vx2) + (-((cse_p117_t1_g559157 * vx1)))));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = (term.kt * ((cse_p1_t5_invr558956 * vy2) + (-((cse_p117_t1_g559157 * vy1)))));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = (term.kt * ((cse_p1_t5_invr558956 * vz2) + (-((cse_p117_t1_g559157 * vz1)))));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p65_t1_g559118 = (cse_p1_t4_invsqrt558955 * n1);
      double g_x2 = (term.kt * ((cse_p11_t57_g559014 * cse_p1_t5_invr558956) + (-((cse_p64_t1_g559104 * ((-((cse_p65_t1_g559117 * vx1))) + (-((cse_p65_t1_g559118 * vx2)))))))));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (term.kt * ((cse_p11_t58_g559015 * cse_p1_t5_invr558956) + (-((cse_p64_t1_g559104 * ((-((cse_p65_t1_g559117 * vy1))) + (-((cse_p65_t1_g559118 * vy2)))))))));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (term.kt * ((cse_p11_t59_g559016 * cse_p1_t5_invr558956) + (-((cse_p64_t1_g559104 * ((-((cse_p65_t1_g559117 * vz1))) + (-((cse_p65_t1_g559118 * vz2)))))))));
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double cse_p117_t1_g559158 = (cse_p64_t1_g559104 * cse_p65_t1_g559118);
      double g_x3 = (term.kt * ((cse_p1_t5_invr558956 * vx1) + (-((cse_p117_t1_g559158 * vx2)))));
      KernelGradientAcc(term.i3x3, 0, g_x3);
      double g_y3 = (term.kt * ((cse_p1_t5_invr558956 * vy1) + (-((cse_p117_t1_g559158 * vy2)))));
      KernelGradientAcc(term.i3x3, 1, g_y3);
      double g_z3 = (term.kt * ((cse_p1_t5_invr558956 * vz1) + (-((cse_p117_t1_g559158 * vz2)))));
      KernelGradientAcc(term.i3x3, 2, g_z3);
      double cse_p11_t1_g559055 = (cse_p529_t1_sqrt559200 * cse_p529_t2_invr559201);
      double cse_p11_t1_g559058 = (cse_p11_t1_g559055 * n2);
      double cse_p166_t1_g559165 = (cse_p11_t1_g559058 * cse_p529_t2_invr559201);
      double cse_p63_t7_g559090 = (cse_p114_t1_g559133 + cse_p166_t1_g559165);
      double cse_p11_t1_g559060 = (cse_p61_t8_invr2559081 * dot);
      double cse_p117_t1_g559159 = (2.00000000000000000e+0 * cse_p11_t1_g559060);
      double cse_p63_t3_g559086 = (cse_p117_t1_g559159 * cse_p61_t3_invr2559076);
      double cse_p63_t4_g559087 = (cse_p117_t1_g559159 * cse_p61_t6_invr2559079);
      double cse_p11_t1_g559062 = (cse_p11_t1_g559060 * cse_p166_t1_g559165);
      double cse_p12_t1_g559071 = (cse_p11_t1_g559058 * cse_p63_t1_g559093 * dot);
      double cse_p63_t6_g559089 = (cse_p12_t1_g559071 * cse_p529_t2_invr559201);
      double cse_p63_t9_g559092 = (cse_p63_t3_g559086 + cse_p63_t6_g559089);
      double cse_p63_t1_g559095 = (cse_p11_t1_g559058 * cse_p61_t8_invr2559081);
      double cse_p115_t1_g559141 = (-2.00000000000000000e+0 * cse_p63_t1_g559095);
      double cse_p117_t1_g559155 = (cse_p63_t1_g559097 * cse_p63_t9_g559092);
      double h_x1_x1 = (term.kt * (cse_p11_t64_g559021 + (cse_p115_t1_g559141 * cse_p64_t1_g559107) + (cse_p117_t1_g559155 * cse_p64_t1_g559106) + (-((cse_p11_t1_g559062 * cse_p64_t1_g559106)))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
      double cse_p65_t1_g559120 = (vx1 * vy2);
      double cse_p65_t1_g559122 = (vx2 * vy1);
      double cse_p115_t1_g559139 = (cse_p11_t55_g559012 * cse_p63_t1_g559097 * dot);
      double h_x1_y1 = (term.kt * ((cse_p115_t1_g559136 * cse_p115_t1_g559139) + (-((cse_p115_t1_g559136 * cse_p11_t1_g559062))) + (-((cse_p63_t1_g559095 * cse_p65_t1_g559120))) + (-((cse_p63_t1_g559095 * cse_p65_t1_g559122)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
      double cse_p65_t1_g559121 = (vx1 * vz2);
      double cse_p66_t1_g559124 = (vx2 * vz1);
      double h_x1_z1 = (term.kt * ((cse_p115_t1_g559137 * cse_p115_t1_g559139) + (-((cse_p115_t1_g559137 * cse_p11_t1_g559062))) + (-((cse_p63_t1_g559095 * cse_p65_t1_g559121))) + (-((cse_p63_t1_g559095 * cse_p66_t1_g559124)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
      double cse_p11_t1_g559064 = (cse_p11_t57_g559014 * cse_p63_t1_g559095);
      double cse_p63_t1_g559103 = (cse_p114_t1_g559133 * cse_p11_t1_g559060);
      double cse_p116_t1_g559146 = (cse_p11_t1_g559055 * cse_p11_t1_g559060);
      double cse_p116_t1_g559147 = (cse_p11_t49_g559006 * cse_p12_t1_g559071);
      double h_x1_x2 = (term.kt * (cse_p11_t81_g559038 + (cse_p116_t1_g559146 * (cse_p11_t33_g558990 + cse_p11_t43_g559000 + n2)) + (cse_p116_t1_g559147 * vx1) + (cse_p11_t60_g559017 * vx2) + (cse_p63_t1_g559103 * cse_p64_t1_g559106) + (-((cse_p11_t1_g559064 * vx1)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
      double cse_p11_t1_g559065 = (cse_p11_t58_g559015 * cse_p63_t1_g559095);
      double cse_p115_t1_g559143 = (cse_p11_t1_g559060 * cse_p529_t1_sqrt559200);
      double cse_p12_t1_g559072 = (cse_p115_t1_g559143 * cse_p11_t53_g559010);
      double cse_p116_t1_g559150 = (cse_p11_t50_g559007 * cse_p12_t1_g559071);
      double h_x1_y2 = (term.kt * ((cse_p116_t1_g559150 * vx1) + (cse_p11_t61_g559018 * vx2) + (cse_p12_t1_g559072 * vx1) + (-((cse_p11_t1_g559065 * vx1)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
      double cse_p11_t1_g559066 = (cse_p11_t59_g559016 * cse_p63_t1_g559095);
      double cse_p12_t1_g559073 = (cse_p115_t1_g559143 * cse_p11_t54_g559011);
      double cse_p116_t1_g559152 = (cse_p11_t51_g559008 * cse_p12_t1_g559071);
      double h_x1_z2 = (term.kt * ((cse_p116_t1_g559152 * vx1) + (cse_p11_t62_g559019 * vx2) + (cse_p12_t1_g559073 * vx1) + (-((cse_p11_t1_g559066 * vx1)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
      double cse_p63_t1_g559101 = (cse_p62_t1_g559083 * cse_p63_t1_g559093);
      double cse_p63_t2_g559085 = (cse_p529_t4_invr559203 * cse_p63_t1_g559101 * dot);
      double cse_p63_t8_g559091 = (cse_p63_t2_g559085 + cse_p63_t4_g559087);
      double cse_p63_t1_g559098 = (cse_p11_t1_g559054 * cse_p11_t1_g559055);
      double cse_p11_t1_g559061 = (cse_p11_t1_g559060 * cse_p63_t1_g559098);
      double cse_p63_t1_g559102 = (cse_p12_t1_g559071 * cse_p62_t1_g559083);
      double h_x1_x3 = (term.kt * (cse_p1_t5_invr558956 + (cse_p63_t1_g559102 * cse_p64_t1_g559107) + (-((cse_p115_t1_g559135 * cse_p64_t1_g559108))) + (-((cse_p11_t1_g559061 * cse_p64_t1_g559107))) + (-((cse_p63_t1_g559095 * cse_p64_t1_g559106)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 0, h_x1_x3);
      double h_x1_y3 = (term.kt * (cse_p11_t66_g559023 + cse_p11_t75_g559032 + (cse_p63_t1_g559102 * cse_p65_t1_g559120) + (-((cse_p11_t1_g559061 * cse_p65_t1_g559120)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 1, h_x1_y3);
      double h_x1_z3 = (term.kt * (cse_p11_t68_g559025 + cse_p11_t76_g559033 + (cse_p63_t1_g559102 * cse_p65_t1_g559121) + (-((cse_p11_t1_g559061 * cse_p65_t1_g559121)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 2, h_x1_z3);
      double h_y1_y1 = (term.kt * (cse_p11_t64_g559021 + (cse_p115_t1_g559141 * cse_p64_t1_g559110) + (cse_p117_t1_g559155 * cse_p64_t1_g559109) + (-((cse_p11_t1_g559062 * cse_p64_t1_g559109)))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
      double cse_p66_t1_g559126 = (vy1 * vz2);
      double cse_p66_t1_g559127 = (vy2 * vz1);
      double h_y1_z1 = (term.kt * ((cse_p115_t1_g559138 * cse_p115_t1_g559139) + (-((cse_p115_t1_g559138 * cse_p11_t1_g559062))) + (-((cse_p63_t1_g559095 * cse_p66_t1_g559126))) + (-((cse_p63_t1_g559095 * cse_p66_t1_g559127)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
      double cse_p11_t1_g559070 = (cse_p115_t1_g559143 * cse_p11_t52_g559009);
      double h_y1_x2 = (term.kt * ((cse_p116_t1_g559147 * vy1) + (cse_p11_t1_g559070 * vy1) + (cse_p11_t60_g559017 * vy2) + (-((cse_p11_t1_g559064 * vy1)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
      double h_y1_y2 = (term.kt * (cse_p11_t81_g559038 + (cse_p116_t1_g559146 * (cse_p11_t38_g558995 + cse_p11_t45_g559002 + n2)) + (cse_p116_t1_g559150 * vy1) + (cse_p11_t61_g559018 * vy2) + (cse_p63_t1_g559103 * cse_p64_t1_g559109) + (-((cse_p11_t1_g559065 * vy1)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
      double h_y1_z2 = (term.kt * ((cse_p116_t1_g559152 * vy1) + (cse_p11_t62_g559019 * vy2) + (cse_p12_t1_g559073 * vy1) + (-((cse_p11_t1_g559066 * vy1)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
      double h_y1_x3 = (term.kt * (cse_p11_t66_g559023 + cse_p11_t75_g559032 + (cse_p63_t1_g559102 * cse_p65_t1_g559122) + (-((cse_p11_t1_g559061 * cse_p65_t1_g559122)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 0, h_y1_x3);
      double h_y1_y3 = (term.kt * (cse_p1_t5_invr558956 + (cse_p63_t1_g559102 * cse_p64_t1_g559110) + (-((cse_p115_t1_g559135 * cse_p64_t1_g559111))) + (-((cse_p11_t1_g559061 * cse_p64_t1_g559110))) + (-((cse_p63_t1_g559095 * cse_p64_t1_g559109)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 1, h_y1_y3);
      double h_y1_z3 = (term.kt * (cse_p11_t71_g559028 + cse_p11_t78_g559035 + (cse_p63_t1_g559102 * cse_p66_t1_g559126) + (-((cse_p11_t1_g559061 * cse_p66_t1_g559126)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 2, h_y1_z3);
      double h_z1_z1 = (term.kt * (cse_p11_t64_g559021 + (cse_p115_t1_g559141 * cse_p64_t1_g559113) + (cse_p117_t1_g559155 * cse_p64_t1_g559112) + (-((cse_p11_t1_g559062 * cse_p64_t1_g559112)))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
      double h_z1_x2 = (term.kt * ((cse_p116_t1_g559147 * vz1) + (cse_p11_t1_g559070 * vz1) + (cse_p11_t60_g559017 * vz2) + (-((cse_p11_t1_g559064 * vz1)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
      double h_z1_y2 = (term.kt * ((cse_p116_t1_g559150 * vz1) + (cse_p11_t61_g559018 * vz2) + (cse_p12_t1_g559072 * vz1) + (-((cse_p11_t1_g559065 * vz1)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
      double h_z1_z2 = (term.kt * (cse_p11_t81_g559038 + (cse_p116_t1_g559146 * (cse_p11_t42_g558999 + cse_p11_t47_g559004 + n2)) + (cse_p116_t1_g559152 * vz1) + (cse_p11_t62_g559019 * vz2) + (cse_p63_t1_g559103 * cse_p64_t1_g559112) + (-((cse_p11_t1_g559066 * vz1)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
      double h_z1_x3 = (term.kt * (cse_p11_t68_g559025 + cse_p11_t76_g559033 + (cse_p63_t1_g559102 * cse_p66_t1_g559124) + (-((cse_p11_t1_g559061 * cse_p66_t1_g559124)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 0, h_z1_x3);
      double h_z1_y3 = (term.kt * (cse_p11_t71_g559028 + cse_p11_t78_g559035 + (cse_p63_t1_g559102 * cse_p66_t1_g559127) + (-((cse_p11_t1_g559061 * cse_p66_t1_g559127)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 1, h_z1_y3);
      double h_z1_z3 = (term.kt * (cse_p1_t5_invr558956 + (cse_p63_t1_g559102 * cse_p64_t1_g559113) + (-((cse_p115_t1_g559135 * cse_p65_t1_g559114))) + (-((cse_p11_t1_g559061 * cse_p64_t1_g559113))) + (-((cse_p63_t1_g559095 * cse_p64_t1_g559112)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 2, h_z1_z3);
      double cse_p65_t1_g559115 = (-2.00000000000000000e+0 * cse_p63_t1_g559099);
      double cse_p116_t1_g559148 = (cse_p11_t49_g559006 * cse_p63_t1_g559093);
      double h_x2_x2 = (term.kt * (cse_p11_t9_g558966 + (2.00000000000000000e+0 * cse_p11_t57_g559014 * cse_p11_t60_g559017) + (dot * ((cse_p116_t1_g559148 * cse_p11_t49_g559006) + (-((cse_p61_t8_invr2559081 * (cse_p11_t20_g558977 + cse_p62_t1_g559083 + (cse_p114_t1_g559133 * cse_p64_t1_g559106) + (cse_p11_t1_g559055 * (cse_p11_t33_g558990 + cse_p11_t43_g559000 + cse_p11_t43_g559000 + n2)) + (cse_p64_t1_g559108 * cse_p65_t1_g559115)))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
      double h_x2_y2 = (term.kt * ((cse_p11_t57_g559014 * cse_p11_t61_g559018) + (cse_p11_t58_g559015 * cse_p11_t60_g559017) + (dot * ((cse_p116_t1_g559148 * cse_p11_t50_g559007) + (-((cse_p61_t8_invr2559081 * (cse_p11_t22_g558979 + (cse_p115_t1_g559136 * cse_p63_t7_g559090) + (cse_p63_t1_g559098 * (cse_p65_t1_g559120 + cse_p65_t1_g559122)) + (cse_p65_t1_g559115 * cse_p65_t1_g559123)))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
      double h_x2_z2 = (term.kt * ((cse_p11_t57_g559014 * cse_p11_t62_g559019) + (cse_p11_t59_g559016 * cse_p11_t60_g559017) + (dot * ((cse_p116_t1_g559148 * cse_p11_t51_g559008) + (-((cse_p61_t8_invr2559081 * (cse_p11_t24_g558981 + (cse_p115_t1_g559137 * cse_p63_t7_g559090) + (cse_p63_t1_g559098 * (cse_p65_t1_g559121 + cse_p66_t1_g559124)) + (cse_p65_t1_g559115 * cse_p66_t1_g559125)))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
      double cse_p11_t1_g559067 = (cse_p115_t1_g559135 * cse_p11_t57_g559014);
      double cse_p116_t1_g559149 = (cse_p11_t49_g559006 * cse_p63_t1_g559101);
      double h_x2_x3 = (term.kt * (cse_p11_t81_g559038 + (cse_p11_t60_g559017 * vx1) + (dot * ((cse_p116_t1_g559149 * vx2) + (-((cse_p61_t8_invr2559081 * (cse_p11_t80_g559037 + (cse_p64_t1_g559108 * cse_p65_t1_g559116) + (-((cse_p63_t1_g559098 * cse_p64_t1_g559107))) + (-(cse_p11_t20_g558977)))))))) + (-((cse_p11_t1_g559067 * vx2)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 0, h_x2_x3);
      double h_x2_y3 = (term.kt * ((cse_p11_t60_g559017 * vy1) + (dot * ((cse_p116_t1_g559149 * vy2) + (-((cse_p61_t8_invr2559081 * (cse_p11_t11_g558968 + cse_p11_t67_g559024 + (-((cse_p63_t1_g559098 * cse_p65_t1_g559120))))))))) + (-((cse_p11_t1_g559067 * vy2)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 1, h_x2_y3);
      double h_x2_z3 = (term.kt * ((cse_p11_t60_g559017 * vz1) + (dot * ((cse_p116_t1_g559149 * vz2) + (-((cse_p61_t8_invr2559081 * (cse_p11_t12_g558969 + cse_p11_t69_g559026 + (-((cse_p63_t1_g559098 * cse_p65_t1_g559121))))))))) + (-((cse_p11_t1_g559067 * vz2)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 2, h_x2_z3);
      double h_y2_y2 = (term.kt * (cse_p11_t9_g558966 + (2.00000000000000000e+0 * cse_p11_t58_g559015 * cse_p11_t61_g559018) + (dot * ((cse_p11_t50_g559007 * cse_p11_t50_g559007 * cse_p63_t1_g559093) + (-((cse_p61_t8_invr2559081 * (cse_p11_t26_g558983 + cse_p62_t1_g559083 + (cse_p114_t1_g559133 * cse_p64_t1_g559109) + (cse_p11_t1_g559055 * (cse_p11_t38_g558995 + cse_p11_t45_g559002 + cse_p11_t45_g559002 + n2)) + (cse_p64_t1_g559111 * cse_p65_t1_g559115)))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
      double h_y2_z2 = (term.kt * ((cse_p11_t58_g559015 * cse_p11_t62_g559019) + (cse_p11_t59_g559016 * cse_p11_t61_g559018) + (dot * ((cse_p11_t50_g559007 * cse_p11_t51_g559008 * cse_p63_t1_g559093) + (-((cse_p61_t8_invr2559081 * (cse_p11_t28_g558985 + (cse_p115_t1_g559138 * cse_p63_t7_g559090) + (cse_p63_t1_g559098 * (cse_p66_t1_g559126 + cse_p66_t1_g559127)) + (cse_p65_t1_g559115 * cse_p66_t1_g559128)))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
      double cse_p11_t1_g559068 = (cse_p115_t1_g559135 * cse_p11_t58_g559015);
      double cse_p116_t1_g559151 = (cse_p11_t50_g559007 * cse_p63_t1_g559101);
      double h_y2_x3 = (term.kt * ((cse_p11_t61_g559018 * vx1) + (dot * ((cse_p116_t1_g559151 * vx2) + (-((cse_p61_t8_invr2559081 * (cse_p11_t11_g558968 + cse_p11_t67_g559024 + (-((cse_p63_t1_g559098 * cse_p65_t1_g559122))))))))) + (-((cse_p11_t1_g559068 * vx2)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 0, h_y2_x3);
      double h_y2_y3 = (term.kt * (cse_p11_t81_g559038 + (cse_p11_t61_g559018 * vy1) + (dot * ((cse_p116_t1_g559151 * vy2) + (-((cse_p61_t8_invr2559081 * (cse_p11_t80_g559037 + (cse_p64_t1_g559111 * cse_p65_t1_g559116) + (-((cse_p63_t1_g559098 * cse_p64_t1_g559110))) + (-(cse_p11_t26_g558983)))))))) + (-((cse_p11_t1_g559068 * vy2)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 1, h_y2_y3);
      double h_y2_z3 = (term.kt * ((cse_p11_t61_g559018 * vz1) + (dot * ((cse_p116_t1_g559151 * vz2) + (-((cse_p61_t8_invr2559081 * (cse_p11_t13_g558970 + cse_p11_t72_g559029 + (-((cse_p63_t1_g559098 * cse_p66_t1_g559126))))))))) + (-((cse_p11_t1_g559068 * vz2)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 2, h_y2_z3);
      double h_z2_z2 = (term.kt * (cse_p11_t9_g558966 + (2.00000000000000000e+0 * cse_p11_t59_g559016 * cse_p11_t62_g559019) + (dot * ((cse_p11_t51_g559008 * cse_p11_t51_g559008 * cse_p63_t1_g559093) + (-((cse_p61_t8_invr2559081 * (cse_p11_t30_g558987 + cse_p62_t1_g559083 + (cse_p114_t1_g559133 * cse_p64_t1_g559112) + (cse_p11_t1_g559055 * (cse_p11_t42_g558999 + cse_p11_t47_g559004 + cse_p11_t47_g559004 + n2)) + (cse_p65_t1_g559114 * cse_p65_t1_g559115)))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      double cse_p11_t1_g559069 = (cse_p115_t1_g559135 * cse_p11_t59_g559016);
      double cse_p116_t1_g559153 = (cse_p11_t51_g559008 * cse_p63_t1_g559101);
      double h_z2_x3 = (term.kt * ((cse_p11_t62_g559019 * vx1) + (dot * ((cse_p116_t1_g559153 * vx2) + (-((cse_p61_t8_invr2559081 * (cse_p11_t12_g558969 + cse_p11_t69_g559026 + (-((cse_p63_t1_g559098 * cse_p66_t1_g559124))))))))) + (-((cse_p11_t1_g559069 * vx2)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 0, h_z2_x3);
      double h_z2_y3 = (term.kt * ((cse_p11_t62_g559019 * vy1) + (dot * ((cse_p116_t1_g559153 * vy2) + (-((cse_p61_t8_invr2559081 * (cse_p11_t13_g558970 + cse_p11_t72_g559029 + (-((cse_p63_t1_g559098 * cse_p66_t1_g559127))))))))) + (-((cse_p11_t1_g559069 * vy2)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 1, h_z2_y3);
      double h_z2_z3 = (term.kt * (cse_p11_t81_g559038 + (cse_p11_t62_g559019 * vz1) + (dot * ((cse_p116_t1_g559153 * vz2) + (-((cse_p61_t8_invr2559081 * (cse_p11_t80_g559037 + (cse_p65_t1_g559114 * cse_p65_t1_g559116) + (-((cse_p63_t1_g559098 * cse_p64_t1_g559113))) + (-(cse_p11_t30_g558987)))))))) + (-((cse_p11_t1_g559069 * vz2)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 2, h_z2_z3);
      double cse_p11_t1_g559063 = (cse_p11_t1_g559060 * cse_p529_t4_invr559203 * cse_p62_t1_g559083);
      double cse_p115_t1_g559142 = (-2.00000000000000000e+0 * cse_p115_t1_g559135);
      double cse_p117_t1_g559154 = (cse_p63_t1_g559094 * cse_p63_t8_g559091);
      double h_x3_x3 = (term.kt * (cse_p11_t63_g559020 + (cse_p115_t1_g559142 * cse_p64_t1_g559107) + (cse_p117_t1_g559154 * cse_p64_t1_g559108) + (-((cse_p11_t1_g559063 * cse_p64_t1_g559108)))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 0, h_x3_x3);
      double cse_p115_t1_g559140 = (cse_p11_t56_g559013 * cse_p63_t1_g559094 * dot);
      double h_x3_y3 = (term.kt * ((cse_p115_t1_g559140 * cse_p65_t1_g559123) + (-((cse_p115_t1_g559135 * cse_p65_t1_g559120))) + (-((cse_p115_t1_g559135 * cse_p65_t1_g559122))) + (-((cse_p11_t1_g559063 * cse_p65_t1_g559123)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 1, h_x3_y3);
      double h_x3_z3 = (term.kt * ((cse_p115_t1_g559140 * cse_p66_t1_g559125) + (-((cse_p115_t1_g559135 * cse_p65_t1_g559121))) + (-((cse_p115_t1_g559135 * cse_p66_t1_g559124))) + (-((cse_p11_t1_g559063 * cse_p66_t1_g559125)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 2, h_x3_z3);
      double h_y3_y3 = (term.kt * (cse_p11_t63_g559020 + (cse_p115_t1_g559142 * cse_p64_t1_g559110) + (cse_p117_t1_g559154 * cse_p64_t1_g559111) + (-((cse_p11_t1_g559063 * cse_p64_t1_g559111)))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 1, h_y3_y3);
      double h_y3_z3 = (term.kt * ((cse_p115_t1_g559140 * cse_p66_t1_g559128) + (-((cse_p115_t1_g559135 * cse_p66_t1_g559126))) + (-((cse_p115_t1_g559135 * cse_p66_t1_g559127))) + (-((cse_p11_t1_g559063 * cse_p66_t1_g559128)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 2, h_y3_z3);
      double h_z3_z3 = (term.kt * (cse_p11_t63_g559020 + (cse_p115_t1_g559142 * cse_p64_t1_g559113) + (cse_p117_t1_g559154 * cse_p65_t1_g559114) + (-((cse_p11_t1_g559063 * cse_p65_t1_g559114)))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x3, 2, h_z3_z3);
    }
  }
  return energy_added;
}
void hessian_fd(const linear_angle_term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
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
