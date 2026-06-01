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
      double cse_p11_t1_g1 = (-(x2));
      double cse_p11_t2_g2 = (-(y2));
      double cse_p11_t3_g3 = (-(z2));
      double vx1 = (cse_p11_t1_g1 + x1);
      double vy1 = (cse_p11_t2_g2 + y1);
      double vz1 = (cse_p11_t3_g3 + z1);
      double vx2 = (cse_p11_t1_g1 + x3);
      double vy2 = (cse_p11_t2_g2 + y3);
      double vz2 = (cse_p11_t3_g3 + z3);
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
      double cse_p11_t1_g7 = (-(x2));
      double cse_p11_t2_g8 = (-(y2));
      double cse_p11_t3_g9 = (-(z2));
      double vx1 = (cse_p11_t1_g7 + x1);
      double vy1 = (cse_p11_t2_g8 + y1);
      double vz1 = (cse_p11_t3_g9 + z1);
      double vx2 = (cse_p11_t1_g7 + x3);
      double vy2 = (cse_p11_t2_g8 + y3);
      double vz2 = (cse_p11_t3_g9 + z3);
      double dot = ((vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2));
      double n1_sq = ((vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1));
      double n2_sq = ((vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2));
      double cse_p1_t1_sqrt1 = sqrt(n1_sq);
      double cse_p1_t2_invsqrt2 = (1.0 / (cse_p1_t1_sqrt1));
      double n1 = cse_p1_t1_sqrt1;
      double cse_p1_t3_sqrt3 = sqrt(n2_sq);
      double cse_p1_t4_invsqrt4 = (1.0 / (cse_p1_t3_sqrt3));
      double n2 = cse_p1_t3_sqrt3;
      double cse_p1_t5_invr5 = (1.0 / ((n1 * n2)));
      double cse_p1_t6_invr26 = (cse_p1_t5_invr5 * cse_p1_t5_invr5);
      double cse_p113_t1_g15 = (cse_p1_t6_invr26 * dot);
      double cos_theta = (cse_p1_t5_invr5 * dot);
      double energy = (term.kt * (1.00000000000000000e+0 + cos_theta));
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t1_g10 = (cse_p1_t2_invsqrt2 * n2);
      double cse_p11_t1_g13 = (cse_p113_t1_g15 * cse_p11_t1_g10);
      double g_x1 = (term.kt * ((cse_p1_t5_invr5 * vx2) + (-((cse_p11_t1_g13 * vx1)))));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = (term.kt * ((cse_p1_t5_invr5 * vy2) + (-((cse_p11_t1_g13 * vy1)))));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = (term.kt * ((cse_p1_t5_invr5 * vz2) + (-((cse_p11_t1_g13 * vz1)))));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p11_t1_g11 = (cse_p1_t4_invsqrt4 * n1);
      double g_x2 = (term.kt * ((cse_p1_t5_invr5 * ((-(vx1)) + (-(vx2)))) + (-((cse_p113_t1_g15 * ((-((cse_p11_t1_g10 * vx1))) + (-((cse_p11_t1_g11 * vx2)))))))));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (term.kt * ((cse_p1_t5_invr5 * ((-(vy1)) + (-(vy2)))) + (-((cse_p113_t1_g15 * ((-((cse_p11_t1_g10 * vy1))) + (-((cse_p11_t1_g11 * vy2)))))))));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (term.kt * ((cse_p1_t5_invr5 * ((-(vz1)) + (-(vz2)))) + (-((cse_p113_t1_g15 * ((-((cse_p11_t1_g10 * vz1))) + (-((cse_p11_t1_g11 * vz2)))))))));
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double cse_p11_t1_g14 = (cse_p113_t1_g15 * cse_p11_t1_g11);
      double g_x3 = (term.kt * ((cse_p1_t5_invr5 * vx1) + (-((cse_p11_t1_g14 * vx2)))));
      KernelGradientAcc(term.i3x3, 0, g_x3);
      double g_y3 = (term.kt * ((cse_p1_t5_invr5 * vy1) + (-((cse_p11_t1_g14 * vy2)))));
      KernelGradientAcc(term.i3x3, 1, g_y3);
      double g_z3 = (term.kt * ((cse_p1_t5_invr5 * vz1) + (-((cse_p11_t1_g14 * vz2)))));
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
      double cse_p11_t88_g94 = (-(x2));
      double cse_p11_t89_g95 = (-(y2));
      double cse_p11_t90_g96 = (-(z2));
      double vx1 = (cse_p11_t88_g94 + x1);
      double cse_p11_t82_g88 = (-(vx1));
      double vy1 = (cse_p11_t89_g95 + y1);
      double cse_p11_t84_g90 = (-(vy1));
      double vz1 = (cse_p11_t90_g96 + z1);
      double cse_p11_t86_g92 = (-(vz1));
      double vx2 = (cse_p11_t88_g94 + x3);
      double cse_p11_t57_g63 = (cse_p11_t82_g88 + (-(vx2)));
      double cse_p11_t83_g89 = (-(vx2));
      double vy2 = (cse_p11_t89_g95 + y3);
      double cse_p11_t58_g64 = (cse_p11_t84_g90 + (-(vy2)));
      double cse_p11_t85_g91 = (-(vy2));
      double vz2 = (cse_p11_t90_g96 + z3);
      double cse_p11_t59_g65 = (cse_p11_t86_g92 + (-(vz2)));
      double cse_p11_t87_g93 = (-(vz2));
      double cse_p64_t1_g156 = (vx1 * vx2);
      double cse_p64_t1_g159 = (vy1 * vy2);
      double cse_p64_t1_g162 = (vz1 * vz2);
      double dot = (cse_p64_t1_g156 + cse_p64_t1_g159 + cse_p64_t1_g162);
      double cse_p64_t1_g155 = (vx1 * vx1);
      double cse_p64_t1_g158 = (vy1 * vy1);
      double cse_p64_t1_g161 = (vz1 * vz1);
      double n1_sq = (cse_p64_t1_g155 + cse_p64_t1_g158 + cse_p64_t1_g161);
      double cse_p529_t1_sqrt249 = sqrt(n1_sq);
      double cse_p529_t2_invr250 = (1.0 / (n1_sq));
      double cse_p61_t3_invr2125 = (cse_p529_t2_invr250 * cse_p529_t2_invr250);
      double cse_p64_t1_g157 = (vx2 * vx2);
      double cse_p64_t1_g160 = (vy2 * vy2);
      double cse_p65_t1_g163 = (vz2 * vz2);
      double n2_sq = (cse_p64_t1_g157 + cse_p64_t1_g160 + cse_p65_t1_g163);
      double cse_p529_t3_sqrt251 = sqrt(n2_sq);
      double cse_p529_t4_invr252 = (1.0 / (n2_sq));
      double cse_p61_t6_invr2128 = (cse_p529_t4_invr252 * cse_p529_t4_invr252);
      double cse_p11_t1_g103 = (cse_p529_t3_sqrt251 * cse_p529_t4_invr252);
      double cse_p11_t43_g49 = (cse_p11_t1_g103 * cse_p64_t1_g156);
      double cse_p116_t1_g193 = (cse_p11_t1_g103 * cse_p529_t2_invr250);
      double cse_p11_t44_g50 = (cse_p116_t1_g193 * vx2);
      double cse_p11_t45_g51 = (cse_p11_t1_g103 * cse_p64_t1_g159);
      double cse_p11_t46_g52 = (cse_p116_t1_g193 * vy2);
      double cse_p11_t47_g53 = (cse_p11_t1_g103 * cse_p64_t1_g162);
      double cse_p11_t48_g54 = (cse_p116_t1_g193 * vz2);
      double cse_p1_t2_invsqrt2 = (1.0 / (cse_p529_t1_sqrt249));
      double n1 = cse_p529_t1_sqrt249;
      double cse_p62_t1_g132 = (cse_p11_t1_g103 * n1);
      double cse_p63_t1_g143 = (cse_p529_t3_sqrt251 * n1);
      double cse_p63_t1_g148 = (cse_p61_t6_invr2128 * cse_p63_t1_g143);
      double cse_p65_t1_g165 = (2.00000000000000000e+0 * cse_p63_t1_g148);
      double cse_p65_t1_g172 = (vx2 * vy2);
      double cse_p11_t11_g17 = (cse_p65_t1_g165 * cse_p65_t1_g172);
      double cse_p66_t1_g174 = (vx2 * vz2);
      double cse_p11_t12_g18 = (cse_p65_t1_g165 * cse_p66_t1_g174);
      double cse_p66_t1_g177 = (vy2 * vz2);
      double cse_p11_t13_g19 = (cse_p65_t1_g165 * cse_p66_t1_g177);
      double cse_p11_t19_g25 = (cse_p62_t1_g132 * vx2);
      double cse_p64_t1_g154 = (cse_p529_t4_invr252 * cse_p62_t1_g132);
      double cse_p11_t20_g26 = (cse_p64_t1_g154 * cse_p64_t1_g157);
      double cse_p11_t22_g28 = (cse_p64_t1_g154 * cse_p65_t1_g172);
      double cse_p11_t24_g30 = (cse_p64_t1_g154 * cse_p66_t1_g174);
      double cse_p11_t25_g31 = (cse_p62_t1_g132 * vy2);
      double cse_p11_t26_g32 = (cse_p64_t1_g154 * cse_p64_t1_g160);
      double cse_p11_t28_g34 = (cse_p64_t1_g154 * cse_p66_t1_g177);
      double cse_p11_t29_g35 = (cse_p62_t1_g132 * vz2);
      double cse_p11_t30_g36 = (cse_p64_t1_g154 * cse_p65_t1_g163);
      double cse_p11_t65_g71 = (-(cse_p11_t19_g25));
      double cse_p11_t67_g73 = (-(cse_p11_t22_g28));
      double cse_p11_t69_g75 = (-(cse_p11_t24_g30));
      double cse_p11_t70_g76 = (-(cse_p11_t25_g31));
      double cse_p11_t72_g78 = (-(cse_p11_t28_g34));
      double cse_p11_t73_g79 = (-(cse_p11_t29_g35));
      double cse_p11_t80_g86 = (-(cse_p62_t1_g132));
      double cse_p1_t4_invsqrt4 = (1.0 / (cse_p529_t3_sqrt251));
      double n2 = cse_p529_t3_sqrt251;
      double cse_p63_t1_g146 = (cse_p529_t1_sqrt249 * n2);
      double cse_p63_t1_g149 = (-2.00000000000000000e+0 * cse_p61_t3_invr2125);
      double cse_p114_t1_g182 = (cse_p63_t1_g146 * cse_p63_t1_g149);
      double cse_p11_t1_g106 = (cse_p529_t2_invr250 * cse_p63_t1_g146);
      double cse_p322_t1_g229 = (cse_p11_t1_g106 * vz1);
      double cse_p270_t1_g224 = (cse_p11_t1_g106 * vy1);
      double cse_p218_t1_g219 = (cse_p11_t1_g106 * vx1);
      double cse_p61_t7_invr129 = (1.0 / ((n1 * n2)));
      double cse_p61_t8_invr2130 = (cse_p61_t7_invr129 * cse_p61_t7_invr129);
      double cse_p61_t9_invr3131 = (cse_p61_t7_invr129 * cse_p61_t8_invr2130);
      double cse_p11_t1_g7 = (cse_p61_t8_invr2130 * (cse_p11_t65_g71 + (-(cse_p218_t1_g219))));
      double cse_p11_t2_g8 = (cse_p61_t8_invr2130 * (cse_p11_t70_g76 + (-(cse_p270_t1_g224))));
      double cse_p11_t3_g9 = (cse_p61_t8_invr2130 * (cse_p11_t73_g79 + (-(cse_p322_t1_g229))));
      double cse_p117_t1_g205 = (cse_p63_t1_g149 * n2);
      double cse_p11_t4_g10 = (cse_p117_t1_g205 * vx1);
      double cse_p11_t5_g11 = (cse_p117_t1_g205 * vy1);
      double cse_p11_t6_g12 = (cse_p117_t1_g205 * vz1);
      double cse_p117_t1_g209 = (2.00000000000000000e+0 * cse_p61_t8_invr2130);
      double cse_p11_t7_g13 = (cse_p117_t1_g209 * cse_p61_t3_invr2125);
      double cse_p11_t8_g14 = (cse_p117_t1_g209 * cse_p61_t6_invr2128);
      double cse_p63_t1_g142 = (2.00000000000000000e+0 * cse_p61_t9_invr3131);
      double cse_p11_t10_g16 = (cse_p63_t1_g142 * cse_p64_t1_g154);
      double cse_p11_t14_g20 = (cse_p11_t1_g106 * cse_p529_t2_invr250 * cse_p63_t1_g142);
      double cse_p115_t1_g184 = (cse_p61_t8_invr2130 * cse_p62_t1_g132);
      double cse_p11_t15_g21 = (cse_p115_t1_g184 * dot);
      double cse_p115_t1_g183 = (cse_p11_t1_g106 * cse_p61_t8_invr2130);
      double cse_p11_t16_g22 = (cse_p115_t1_g183 * dot);
      double cse_p11_t18_g24 = (n1 * n2);
      double cse_p11_t21_g27 = (cse_p115_t1_g184 * cse_p65_t1_g172);
      double cse_p11_t23_g29 = (cse_p115_t1_g184 * cse_p66_t1_g174);
      double cse_p11_t27_g33 = (cse_p115_t1_g184 * cse_p66_t1_g177);
      double cse_p65_t1_g168 = (cse_p529_t2_invr250 * n2);
      double cse_p116_t1_g194 = (cse_p529_t2_invr250 * cse_p65_t1_g168);
      double cse_p11_t31_g37 = (cse_p116_t1_g194 * vx1);
      double cse_p11_t33_g39 = (cse_p64_t1_g155 * cse_p65_t1_g168);
      double cse_p115_t1_g185 = (vx1 * vy1);
      double cse_p11_t34_g40 = (cse_p115_t1_g183 * cse_p115_t1_g185);
      double cse_p115_t1_g186 = (vx1 * vz1);
      double cse_p11_t35_g41 = (cse_p115_t1_g183 * cse_p115_t1_g186);
      double cse_p11_t36_g42 = (cse_p116_t1_g194 * vy1);
      double cse_p11_t38_g44 = (cse_p64_t1_g158 * cse_p65_t1_g168);
      double cse_p115_t1_g187 = (vy1 * vz1);
      double cse_p11_t39_g45 = (cse_p115_t1_g183 * cse_p115_t1_g187);
      double cse_p11_t40_g46 = (cse_p116_t1_g194 * vz1);
      double cse_p11_t42_g48 = (cse_p64_t1_g161 * cse_p65_t1_g168);
      double cse_p11_t49_g55 = (cse_p11_t65_g71 + (-(cse_p218_t1_g219)));
      double cse_p11_t50_g56 = (cse_p11_t70_g76 + (-(cse_p270_t1_g224)));
      double cse_p11_t51_g57 = (cse_p11_t73_g79 + (-(cse_p322_t1_g229)));
      double cse_p11_t52_g58 = (cse_p11_t31_g37 + cse_p11_t44_g50 + cse_p11_t4_g10);
      double cse_p11_t53_g59 = (cse_p11_t36_g42 + cse_p11_t46_g52 + cse_p11_t5_g11);
      double cse_p11_t54_g60 = (cse_p11_t40_g46 + cse_p11_t48_g54 + cse_p11_t6_g12);
      double cse_p11_t55_g61 = (cse_p11_t14_g20 + cse_p11_t7_g13);
      double cse_p11_t56_g62 = (cse_p11_t10_g16 + cse_p11_t8_g14);
      double cse_p11_t60_g66 = (-(cse_p11_t1_g7));
      double cse_p11_t61_g67 = (-(cse_p11_t2_g8));
      double cse_p11_t62_g68 = (-(cse_p11_t3_g9));
      double cse_p11_t63_g69 = (-(cse_p11_t15_g21));
      double cse_p11_t64_g70 = (-(cse_p11_t16_g22));
      double cse_p11_t66_g72 = (-(cse_p11_t21_g27));
      double cse_p11_t68_g74 = (-(cse_p11_t23_g29));
      double cse_p11_t71_g77 = (-(cse_p11_t27_g33));
      double cse_p11_t74_g80 = (-(cse_p218_t1_g219));
      double cse_p11_t75_g81 = (-(cse_p11_t34_g40));
      double cse_p11_t76_g82 = (-(cse_p11_t35_g41));
      double cse_p11_t77_g83 = (-(cse_p270_t1_g224));
      double cse_p11_t78_g84 = (-(cse_p11_t39_g45));
      double cse_p11_t79_g85 = (-(cse_p322_t1_g229));
      double cse_p1_t5_invr5 = (1.0 / (cse_p11_t18_g24));
      double cse_p11_t9_g15 = (2.00000000000000000e+0 * cse_p1_t5_invr5);
      double cse_p11_t81_g87 = (-(cse_p1_t5_invr5));
      double cse_p1_t6_invr26 = (cse_p1_t5_invr5 * cse_p1_t5_invr5);
      double cos_theta = (cse_p1_t5_invr5 * dot);
      double energy = (term.kt * (1.00000000000000000e+0 + cos_theta));
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p64_t1_g153 = (cse_p1_t6_invr26 * dot);
      double cse_p65_t1_g166 = (cse_p1_t2_invsqrt2 * n2);
      double cse_p117_t1_g206 = (cse_p64_t1_g153 * cse_p65_t1_g166);
      double g_x1 = (term.kt * ((cse_p1_t5_invr5 * vx2) + (-((cse_p117_t1_g206 * vx1)))));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = (term.kt * ((cse_p1_t5_invr5 * vy2) + (-((cse_p117_t1_g206 * vy1)))));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = (term.kt * ((cse_p1_t5_invr5 * vz2) + (-((cse_p117_t1_g206 * vz1)))));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p65_t1_g167 = (cse_p1_t4_invsqrt4 * n1);
      double g_x2 = (term.kt * ((cse_p11_t57_g63 * cse_p1_t5_invr5) + (-((cse_p64_t1_g153 * ((-((cse_p65_t1_g166 * vx1))) + (-((cse_p65_t1_g167 * vx2)))))))));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (term.kt * ((cse_p11_t58_g64 * cse_p1_t5_invr5) + (-((cse_p64_t1_g153 * ((-((cse_p65_t1_g166 * vy1))) + (-((cse_p65_t1_g167 * vy2)))))))));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (term.kt * ((cse_p11_t59_g65 * cse_p1_t5_invr5) + (-((cse_p64_t1_g153 * ((-((cse_p65_t1_g166 * vz1))) + (-((cse_p65_t1_g167 * vz2)))))))));
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double cse_p117_t1_g207 = (cse_p64_t1_g153 * cse_p65_t1_g167);
      double g_x3 = (term.kt * ((cse_p1_t5_invr5 * vx1) + (-((cse_p117_t1_g207 * vx2)))));
      KernelGradientAcc(term.i3x3, 0, g_x3);
      double g_y3 = (term.kt * ((cse_p1_t5_invr5 * vy1) + (-((cse_p117_t1_g207 * vy2)))));
      KernelGradientAcc(term.i3x3, 1, g_y3);
      double g_z3 = (term.kt * ((cse_p1_t5_invr5 * vz1) + (-((cse_p117_t1_g207 * vz2)))));
      KernelGradientAcc(term.i3x3, 2, g_z3);
      double cse_p11_t1_g104 = (cse_p529_t1_sqrt249 * cse_p529_t2_invr250);
      double cse_p11_t1_g107 = (cse_p11_t1_g104 * n2);
      double cse_p166_t1_g214 = (cse_p11_t1_g107 * cse_p529_t2_invr250);
      double cse_p63_t7_g139 = (cse_p114_t1_g182 + cse_p166_t1_g214);
      double cse_p11_t1_g109 = (cse_p61_t8_invr2130 * dot);
      double cse_p117_t1_g208 = (2.00000000000000000e+0 * cse_p11_t1_g109);
      double cse_p63_t3_g135 = (cse_p117_t1_g208 * cse_p61_t3_invr2125);
      double cse_p63_t4_g136 = (cse_p117_t1_g208 * cse_p61_t6_invr2128);
      double cse_p11_t1_g111 = (cse_p11_t1_g109 * cse_p166_t1_g214);
      double cse_p12_t1_g120 = (cse_p11_t1_g107 * cse_p63_t1_g142 * dot);
      double cse_p63_t6_g138 = (cse_p12_t1_g120 * cse_p529_t2_invr250);
      double cse_p63_t9_g141 = (cse_p63_t3_g135 + cse_p63_t6_g138);
      double cse_p63_t1_g144 = (cse_p11_t1_g107 * cse_p61_t8_invr2130);
      double cse_p115_t1_g190 = (-2.00000000000000000e+0 * cse_p63_t1_g144);
      double cse_p117_t1_g204 = (cse_p63_t1_g146 * cse_p63_t9_g141);
      double h_x1_x1 = (term.kt * (cse_p11_t64_g70 + (cse_p115_t1_g190 * cse_p64_t1_g156) + (cse_p117_t1_g204 * cse_p64_t1_g155) + (-((cse_p11_t1_g111 * cse_p64_t1_g155)))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
      double cse_p65_t1_g169 = (vx1 * vy2);
      double cse_p65_t1_g171 = (vx2 * vy1);
      double cse_p115_t1_g188 = (cse_p11_t55_g61 * cse_p63_t1_g146 * dot);
      double h_x1_y1 = (term.kt * ((cse_p115_t1_g185 * cse_p115_t1_g188) + (-((cse_p115_t1_g185 * cse_p11_t1_g111))) + (-((cse_p63_t1_g144 * cse_p65_t1_g169))) + (-((cse_p63_t1_g144 * cse_p65_t1_g171)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
      double cse_p65_t1_g170 = (vx1 * vz2);
      double cse_p66_t1_g173 = (vx2 * vz1);
      double h_x1_z1 = (term.kt * ((cse_p115_t1_g186 * cse_p115_t1_g188) + (-((cse_p115_t1_g186 * cse_p11_t1_g111))) + (-((cse_p63_t1_g144 * cse_p65_t1_g170))) + (-((cse_p63_t1_g144 * cse_p66_t1_g173)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
      double cse_p11_t1_g113 = (cse_p11_t57_g63 * cse_p63_t1_g144);
      double cse_p63_t1_g152 = (cse_p114_t1_g182 * cse_p11_t1_g109);
      double cse_p116_t1_g195 = (cse_p11_t1_g104 * cse_p11_t1_g109);
      double cse_p116_t1_g196 = (cse_p11_t49_g55 * cse_p12_t1_g120);
      double h_x1_x2 = (term.kt * (cse_p11_t81_g87 + (cse_p116_t1_g195 * (cse_p11_t33_g39 + cse_p11_t43_g49 + n2)) + (cse_p116_t1_g196 * vx1) + (cse_p11_t60_g66 * vx2) + (cse_p63_t1_g152 * cse_p64_t1_g155) + (-((cse_p11_t1_g113 * vx1)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
      double cse_p11_t1_g114 = (cse_p11_t58_g64 * cse_p63_t1_g144);
      double cse_p115_t1_g192 = (cse_p11_t1_g109 * cse_p529_t1_sqrt249);
      double cse_p12_t1_g121 = (cse_p115_t1_g192 * cse_p11_t53_g59);
      double cse_p116_t1_g199 = (cse_p11_t50_g56 * cse_p12_t1_g120);
      double h_x1_y2 = (term.kt * ((cse_p116_t1_g199 * vx1) + (cse_p11_t61_g67 * vx2) + (cse_p12_t1_g121 * vx1) + (-((cse_p11_t1_g114 * vx1)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
      double cse_p11_t1_g115 = (cse_p11_t59_g65 * cse_p63_t1_g144);
      double cse_p12_t1_g122 = (cse_p115_t1_g192 * cse_p11_t54_g60);
      double cse_p116_t1_g201 = (cse_p11_t51_g57 * cse_p12_t1_g120);
      double h_x1_z2 = (term.kt * ((cse_p116_t1_g201 * vx1) + (cse_p11_t62_g68 * vx2) + (cse_p12_t1_g122 * vx1) + (-((cse_p11_t1_g115 * vx1)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
      double cse_p63_t1_g150 = (cse_p62_t1_g132 * cse_p63_t1_g142);
      double cse_p63_t2_g134 = (cse_p529_t4_invr252 * cse_p63_t1_g150 * dot);
      double cse_p63_t8_g140 = (cse_p63_t2_g134 + cse_p63_t4_g136);
      double cse_p63_t1_g147 = (cse_p11_t1_g103 * cse_p11_t1_g104);
      double cse_p11_t1_g110 = (cse_p11_t1_g109 * cse_p63_t1_g147);
      double cse_p63_t1_g151 = (cse_p12_t1_g120 * cse_p62_t1_g132);
      double h_x1_x3 = (term.kt * (cse_p1_t5_invr5 + (cse_p63_t1_g151 * cse_p64_t1_g156) + (-((cse_p115_t1_g184 * cse_p64_t1_g157))) + (-((cse_p11_t1_g110 * cse_p64_t1_g156))) + (-((cse_p63_t1_g144 * cse_p64_t1_g155)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 0, h_x1_x3);
      double h_x1_y3 = (term.kt * (cse_p11_t66_g72 + cse_p11_t75_g81 + (cse_p63_t1_g151 * cse_p65_t1_g169) + (-((cse_p11_t1_g110 * cse_p65_t1_g169)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 1, h_x1_y3);
      double h_x1_z3 = (term.kt * (cse_p11_t68_g74 + cse_p11_t76_g82 + (cse_p63_t1_g151 * cse_p65_t1_g170) + (-((cse_p11_t1_g110 * cse_p65_t1_g170)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 2, h_x1_z3);
      double h_y1_y1 = (term.kt * (cse_p11_t64_g70 + (cse_p115_t1_g190 * cse_p64_t1_g159) + (cse_p117_t1_g204 * cse_p64_t1_g158) + (-((cse_p11_t1_g111 * cse_p64_t1_g158)))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
      double cse_p66_t1_g175 = (vy1 * vz2);
      double cse_p66_t1_g176 = (vy2 * vz1);
      double h_y1_z1 = (term.kt * ((cse_p115_t1_g187 * cse_p115_t1_g188) + (-((cse_p115_t1_g187 * cse_p11_t1_g111))) + (-((cse_p63_t1_g144 * cse_p66_t1_g175))) + (-((cse_p63_t1_g144 * cse_p66_t1_g176)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
      double cse_p11_t1_g119 = (cse_p115_t1_g192 * cse_p11_t52_g58);
      double h_y1_x2 = (term.kt * ((cse_p116_t1_g196 * vy1) + (cse_p11_t1_g119 * vy1) + (cse_p11_t60_g66 * vy2) + (-((cse_p11_t1_g113 * vy1)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
      double h_y1_y2 = (term.kt * (cse_p11_t81_g87 + (cse_p116_t1_g195 * (cse_p11_t38_g44 + cse_p11_t45_g51 + n2)) + (cse_p116_t1_g199 * vy1) + (cse_p11_t61_g67 * vy2) + (cse_p63_t1_g152 * cse_p64_t1_g158) + (-((cse_p11_t1_g114 * vy1)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
      double h_y1_z2 = (term.kt * ((cse_p116_t1_g201 * vy1) + (cse_p11_t62_g68 * vy2) + (cse_p12_t1_g122 * vy1) + (-((cse_p11_t1_g115 * vy1)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
      double h_y1_x3 = (term.kt * (cse_p11_t66_g72 + cse_p11_t75_g81 + (cse_p63_t1_g151 * cse_p65_t1_g171) + (-((cse_p11_t1_g110 * cse_p65_t1_g171)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 0, h_y1_x3);
      double h_y1_y3 = (term.kt * (cse_p1_t5_invr5 + (cse_p63_t1_g151 * cse_p64_t1_g159) + (-((cse_p115_t1_g184 * cse_p64_t1_g160))) + (-((cse_p11_t1_g110 * cse_p64_t1_g159))) + (-((cse_p63_t1_g144 * cse_p64_t1_g158)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 1, h_y1_y3);
      double h_y1_z3 = (term.kt * (cse_p11_t71_g77 + cse_p11_t78_g84 + (cse_p63_t1_g151 * cse_p66_t1_g175) + (-((cse_p11_t1_g110 * cse_p66_t1_g175)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 2, h_y1_z3);
      double h_z1_z1 = (term.kt * (cse_p11_t64_g70 + (cse_p115_t1_g190 * cse_p64_t1_g162) + (cse_p117_t1_g204 * cse_p64_t1_g161) + (-((cse_p11_t1_g111 * cse_p64_t1_g161)))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
      double h_z1_x2 = (term.kt * ((cse_p116_t1_g196 * vz1) + (cse_p11_t1_g119 * vz1) + (cse_p11_t60_g66 * vz2) + (-((cse_p11_t1_g113 * vz1)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
      double h_z1_y2 = (term.kt * ((cse_p116_t1_g199 * vz1) + (cse_p11_t61_g67 * vz2) + (cse_p12_t1_g121 * vz1) + (-((cse_p11_t1_g114 * vz1)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
      double h_z1_z2 = (term.kt * (cse_p11_t81_g87 + (cse_p116_t1_g195 * (cse_p11_t42_g48 + cse_p11_t47_g53 + n2)) + (cse_p116_t1_g201 * vz1) + (cse_p11_t62_g68 * vz2) + (cse_p63_t1_g152 * cse_p64_t1_g161) + (-((cse_p11_t1_g115 * vz1)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
      double h_z1_x3 = (term.kt * (cse_p11_t68_g74 + cse_p11_t76_g82 + (cse_p63_t1_g151 * cse_p66_t1_g173) + (-((cse_p11_t1_g110 * cse_p66_t1_g173)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 0, h_z1_x3);
      double h_z1_y3 = (term.kt * (cse_p11_t71_g77 + cse_p11_t78_g84 + (cse_p63_t1_g151 * cse_p66_t1_g176) + (-((cse_p11_t1_g110 * cse_p66_t1_g176)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 1, h_z1_y3);
      double h_z1_z3 = (term.kt * (cse_p1_t5_invr5 + (cse_p63_t1_g151 * cse_p64_t1_g162) + (-((cse_p115_t1_g184 * cse_p65_t1_g163))) + (-((cse_p11_t1_g110 * cse_p64_t1_g162))) + (-((cse_p63_t1_g144 * cse_p64_t1_g161)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 2, h_z1_z3);
      double cse_p65_t1_g164 = (-2.00000000000000000e+0 * cse_p63_t1_g148);
      double cse_p116_t1_g197 = (cse_p11_t49_g55 * cse_p63_t1_g142);
      double h_x2_x2 = (term.kt * (cse_p11_t9_g15 + (2.00000000000000000e+0 * cse_p11_t57_g63 * cse_p11_t60_g66) + (dot * ((cse_p116_t1_g197 * cse_p11_t49_g55) + (-((cse_p61_t8_invr2130 * (cse_p11_t20_g26 + cse_p62_t1_g132 + (cse_p114_t1_g182 * cse_p64_t1_g155) + (cse_p11_t1_g104 * (cse_p11_t33_g39 + cse_p11_t43_g49 + cse_p11_t43_g49 + n2)) + (cse_p64_t1_g157 * cse_p65_t1_g164)))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
      double h_x2_y2 = (term.kt * ((cse_p11_t57_g63 * cse_p11_t61_g67) + (cse_p11_t58_g64 * cse_p11_t60_g66) + (dot * ((cse_p116_t1_g197 * cse_p11_t50_g56) + (-((cse_p61_t8_invr2130 * (cse_p11_t22_g28 + (cse_p115_t1_g185 * cse_p63_t7_g139) + (cse_p63_t1_g147 * (cse_p65_t1_g169 + cse_p65_t1_g171)) + (cse_p65_t1_g164 * cse_p65_t1_g172)))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
      double h_x2_z2 = (term.kt * ((cse_p11_t57_g63 * cse_p11_t62_g68) + (cse_p11_t59_g65 * cse_p11_t60_g66) + (dot * ((cse_p116_t1_g197 * cse_p11_t51_g57) + (-((cse_p61_t8_invr2130 * (cse_p11_t24_g30 + (cse_p115_t1_g186 * cse_p63_t7_g139) + (cse_p63_t1_g147 * (cse_p65_t1_g170 + cse_p66_t1_g173)) + (cse_p65_t1_g164 * cse_p66_t1_g174)))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
      double cse_p11_t1_g116 = (cse_p115_t1_g184 * cse_p11_t57_g63);
      double cse_p116_t1_g198 = (cse_p11_t49_g55 * cse_p63_t1_g150);
      double h_x2_x3 = (term.kt * (cse_p11_t81_g87 + (cse_p11_t60_g66 * vx1) + (dot * ((cse_p116_t1_g198 * vx2) + (-((cse_p61_t8_invr2130 * (cse_p11_t80_g86 + (cse_p64_t1_g157 * cse_p65_t1_g165) + (-((cse_p63_t1_g147 * cse_p64_t1_g156))) + (-(cse_p11_t20_g26)))))))) + (-((cse_p11_t1_g116 * vx2)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 0, h_x2_x3);
      double h_x2_y3 = (term.kt * ((cse_p11_t60_g66 * vy1) + (dot * ((cse_p116_t1_g198 * vy2) + (-((cse_p61_t8_invr2130 * (cse_p11_t11_g17 + cse_p11_t67_g73 + (-((cse_p63_t1_g147 * cse_p65_t1_g169))))))))) + (-((cse_p11_t1_g116 * vy2)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 1, h_x2_y3);
      double h_x2_z3 = (term.kt * ((cse_p11_t60_g66 * vz1) + (dot * ((cse_p116_t1_g198 * vz2) + (-((cse_p61_t8_invr2130 * (cse_p11_t12_g18 + cse_p11_t69_g75 + (-((cse_p63_t1_g147 * cse_p65_t1_g170))))))))) + (-((cse_p11_t1_g116 * vz2)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 2, h_x2_z3);
      double h_y2_y2 = (term.kt * (cse_p11_t9_g15 + (2.00000000000000000e+0 * cse_p11_t58_g64 * cse_p11_t61_g67) + (dot * ((cse_p11_t50_g56 * cse_p11_t50_g56 * cse_p63_t1_g142) + (-((cse_p61_t8_invr2130 * (cse_p11_t26_g32 + cse_p62_t1_g132 + (cse_p114_t1_g182 * cse_p64_t1_g158) + (cse_p11_t1_g104 * (cse_p11_t38_g44 + cse_p11_t45_g51 + cse_p11_t45_g51 + n2)) + (cse_p64_t1_g160 * cse_p65_t1_g164)))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
      double h_y2_z2 = (term.kt * ((cse_p11_t58_g64 * cse_p11_t62_g68) + (cse_p11_t59_g65 * cse_p11_t61_g67) + (dot * ((cse_p11_t50_g56 * cse_p11_t51_g57 * cse_p63_t1_g142) + (-((cse_p61_t8_invr2130 * (cse_p11_t28_g34 + (cse_p115_t1_g187 * cse_p63_t7_g139) + (cse_p63_t1_g147 * (cse_p66_t1_g175 + cse_p66_t1_g176)) + (cse_p65_t1_g164 * cse_p66_t1_g177)))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
      double cse_p11_t1_g117 = (cse_p115_t1_g184 * cse_p11_t58_g64);
      double cse_p116_t1_g200 = (cse_p11_t50_g56 * cse_p63_t1_g150);
      double h_y2_x3 = (term.kt * ((cse_p11_t61_g67 * vx1) + (dot * ((cse_p116_t1_g200 * vx2) + (-((cse_p61_t8_invr2130 * (cse_p11_t11_g17 + cse_p11_t67_g73 + (-((cse_p63_t1_g147 * cse_p65_t1_g171))))))))) + (-((cse_p11_t1_g117 * vx2)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 0, h_y2_x3);
      double h_y2_y3 = (term.kt * (cse_p11_t81_g87 + (cse_p11_t61_g67 * vy1) + (dot * ((cse_p116_t1_g200 * vy2) + (-((cse_p61_t8_invr2130 * (cse_p11_t80_g86 + (cse_p64_t1_g160 * cse_p65_t1_g165) + (-((cse_p63_t1_g147 * cse_p64_t1_g159))) + (-(cse_p11_t26_g32)))))))) + (-((cse_p11_t1_g117 * vy2)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 1, h_y2_y3);
      double h_y2_z3 = (term.kt * ((cse_p11_t61_g67 * vz1) + (dot * ((cse_p116_t1_g200 * vz2) + (-((cse_p61_t8_invr2130 * (cse_p11_t13_g19 + cse_p11_t72_g78 + (-((cse_p63_t1_g147 * cse_p66_t1_g175))))))))) + (-((cse_p11_t1_g117 * vz2)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 2, h_y2_z3);
      double h_z2_z2 = (term.kt * (cse_p11_t9_g15 + (2.00000000000000000e+0 * cse_p11_t59_g65 * cse_p11_t62_g68) + (dot * ((cse_p11_t51_g57 * cse_p11_t51_g57 * cse_p63_t1_g142) + (-((cse_p61_t8_invr2130 * (cse_p11_t30_g36 + cse_p62_t1_g132 + (cse_p114_t1_g182 * cse_p64_t1_g161) + (cse_p11_t1_g104 * (cse_p11_t42_g48 + cse_p11_t47_g53 + cse_p11_t47_g53 + n2)) + (cse_p65_t1_g163 * cse_p65_t1_g164)))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      double cse_p11_t1_g118 = (cse_p115_t1_g184 * cse_p11_t59_g65);
      double cse_p116_t1_g202 = (cse_p11_t51_g57 * cse_p63_t1_g150);
      double h_z2_x3 = (term.kt * ((cse_p11_t62_g68 * vx1) + (dot * ((cse_p116_t1_g202 * vx2) + (-((cse_p61_t8_invr2130 * (cse_p11_t12_g18 + cse_p11_t69_g75 + (-((cse_p63_t1_g147 * cse_p66_t1_g173))))))))) + (-((cse_p11_t1_g118 * vx2)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 0, h_z2_x3);
      double h_z2_y3 = (term.kt * ((cse_p11_t62_g68 * vy1) + (dot * ((cse_p116_t1_g202 * vy2) + (-((cse_p61_t8_invr2130 * (cse_p11_t13_g19 + cse_p11_t72_g78 + (-((cse_p63_t1_g147 * cse_p66_t1_g176))))))))) + (-((cse_p11_t1_g118 * vy2)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 1, h_z2_y3);
      double h_z2_z3 = (term.kt * (cse_p11_t81_g87 + (cse_p11_t62_g68 * vz1) + (dot * ((cse_p116_t1_g202 * vz2) + (-((cse_p61_t8_invr2130 * (cse_p11_t80_g86 + (cse_p65_t1_g163 * cse_p65_t1_g165) + (-((cse_p63_t1_g147 * cse_p64_t1_g162))) + (-(cse_p11_t30_g36)))))))) + (-((cse_p11_t1_g118 * vz2)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 2, h_z2_z3);
      double cse_p11_t1_g112 = (cse_p11_t1_g109 * cse_p529_t4_invr252 * cse_p62_t1_g132);
      double cse_p115_t1_g191 = (-2.00000000000000000e+0 * cse_p115_t1_g184);
      double cse_p117_t1_g203 = (cse_p63_t1_g143 * cse_p63_t8_g140);
      double h_x3_x3 = (term.kt * (cse_p11_t63_g69 + (cse_p115_t1_g191 * cse_p64_t1_g156) + (cse_p117_t1_g203 * cse_p64_t1_g157) + (-((cse_p11_t1_g112 * cse_p64_t1_g157)))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 0, h_x3_x3);
      double cse_p115_t1_g189 = (cse_p11_t56_g62 * cse_p63_t1_g143 * dot);
      double h_x3_y3 = (term.kt * ((cse_p115_t1_g189 * cse_p65_t1_g172) + (-((cse_p115_t1_g184 * cse_p65_t1_g169))) + (-((cse_p115_t1_g184 * cse_p65_t1_g171))) + (-((cse_p11_t1_g112 * cse_p65_t1_g172)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 1, h_x3_y3);
      double h_x3_z3 = (term.kt * ((cse_p115_t1_g189 * cse_p66_t1_g174) + (-((cse_p115_t1_g184 * cse_p65_t1_g170))) + (-((cse_p115_t1_g184 * cse_p66_t1_g173))) + (-((cse_p11_t1_g112 * cse_p66_t1_g174)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 2, h_x3_z3);
      double h_y3_y3 = (term.kt * (cse_p11_t63_g69 + (cse_p115_t1_g191 * cse_p64_t1_g159) + (cse_p117_t1_g203 * cse_p64_t1_g160) + (-((cse_p11_t1_g112 * cse_p64_t1_g160)))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 1, h_y3_y3);
      double h_y3_z3 = (term.kt * ((cse_p115_t1_g189 * cse_p66_t1_g177) + (-((cse_p115_t1_g184 * cse_p66_t1_g175))) + (-((cse_p115_t1_g184 * cse_p66_t1_g176))) + (-((cse_p11_t1_g112 * cse_p66_t1_g177)))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 2, h_y3_z3);
      double h_z3_z3 = (term.kt * (cse_p11_t63_g69 + (cse_p115_t1_g191 * cse_p64_t1_g162) + (cse_p117_t1_g203 * cse_p65_t1_g163) + (-((cse_p11_t1_g112 * cse_p65_t1_g163)))));
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
