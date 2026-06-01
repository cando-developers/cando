#include "stretch.h"

template <typename HESSIAN>
struct Stretch {
  static constexpr size_t PositionSize = 6;
  static std::string description() { return "mathkernel-stretch"; };
double hessian(const TermStretch& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double energy_added = 0.00000000000000000e+0;
  {
    /* !BASE */
    DOUBLE x1 = position[term.i3x1 + 0];
    DOUBLE y1 = position[term.i3x1 + 1];
    DOUBLE z1 = position[term.i3x1 + 2];
    DOUBLE x2 = position[term.i3x2 + 0];
    DOUBLE y2 = position[term.i3x2 + 1];
    DOUBLE z2 = position[term.i3x2 + 2];
    {
      /* !BASE */
      double dx = (x2 + (-(x1)));
      double dy = (y2 + (-(y1)));
      double dz = (z2 + (-(z1)));
      double cse_p11_t1_g62 = (dy * dy);
      double cse_p11_t1_g63 = (dz * dz);
      double cse_p11_t1_g64 = (dx * dx);
      double r2 = (cse_p11_t1_g62 + cse_p11_t1_g63 + cse_p11_t1_g64);
      double cse_p529_t1_sqrt130 = sqrt(r2);
      double cse_p529_t2_invr131 = (1.0 / (r2));
      double cse_p61_t3_invr281 = (cse_p529_t2_invr131 * cse_p529_t2_invr131);
      double cse_p11_t1_g56 = (cse_p529_t1_sqrt130 * cse_p529_t2_invr131);
      double cse_p63_t1_g88 = (cse_p11_t1_g56 * term.kb);
      double cse_p11_t1_g58 = (cse_p529_t2_invr131 * cse_p63_t1_g88);
      double cse_p11_t1_g61 = (cse_p11_t1_g58 * dx);
      double cse_p12_t1_g68 = (cse_p11_t1_g61 * cse_p529_t1_sqrt130);
      double cse_p11_t7_g9 = (-2.00000000000000000e+0 * cse_p12_t1_g68 * dy);
      double cse_p11_t8_g10 = (-2.00000000000000000e+0 * cse_p12_t1_g68 * dz);
      double cse_p11_t1_g65 = (dy * dz);
      double cse_p12_t1_g69 = (cse_p11_t1_g58 * cse_p11_t1_g65);
      double cse_p13_t1_g78 = (cse_p12_t1_g69 * cse_p529_t1_sqrt130);
      double cse_p11_t9_g11 = (2.00000000000000000e+0 * (-(cse_p13_t1_g78)));
      double cse_p11_t26_g28 = (cse_p11_t1_g56 * cse_p11_t1_g64);
      double cse_p11_t28_g30 = (cse_p12_t1_g68 * dy);
      double cse_p11_t30_g32 = (cse_p12_t1_g68 * dz);
      double cse_p11_t31_g33 = (cse_p11_t1_g56 * cse_p11_t1_g62);
      double cse_p11_t34_g36 = (cse_p11_t1_g56 * cse_p11_t1_g63);
      double cse_p11_t48_g50 = (-(cse_p11_t28_g30));
      double cse_p114_t1_g92 = (2.00000000000000000e+0 * cse_p11_t48_g50);
      double cse_p11_t49_g51 = (-(cse_p11_t30_g32));
      double cse_p63_t2_g84 = (2.00000000000000000e+0 * cse_p11_t49_g51);
      double cse_p11_t50_g52 = (-(cse_p13_t1_g78));
      double cse_p63_t3_g85 = (2.00000000000000000e+0 * cse_p11_t50_g52);
      double cse_p1_t2_invsqrt2 = (1.0 / (cse_p529_t1_sqrt130));
      double r = cse_p529_t1_sqrt130;
      double dr = (r + (-(term.r0)));
      double cse_p62_t1_g82 = (cse_p12_t1_g69 * dr);
      double cse_p13_t1_g76 = (dr * term.kb);
      double cse_p11_t1_g57 = (cse_p13_t1_g76 * cse_p529_t1_sqrt130 * cse_p61_t3_invr281);
      double cse_p11_t1_g60 = (-2.00000000000000000e+0 * cse_p11_t1_g57);
      double cse_p11_t1_g3 = (cse_p11_t1_g60 * cse_p11_t1_g64);
      double cse_p12_t1_g74 = (cse_p11_t1_g60 * dx);
      double cse_p11_t2_g4 = (cse_p12_t1_g74 * dy);
      double cse_p11_t3_g5 = (cse_p12_t1_g74 * dz);
      double cse_p11_t4_g6 = (cse_p11_t1_g60 * cse_p11_t1_g62);
      double cse_p11_t5_g7 = (cse_p11_t1_g60 * cse_p11_t1_g65);
      double cse_p11_t6_g8 = (cse_p11_t1_g60 * cse_p11_t1_g63);
      double cse_p12_t1_g67 = (2.00000000000000000e+0 * cse_p11_t1_g57);
      double cse_p11_t10_g12 = (cse_p11_t1_g64 * cse_p12_t1_g67);
      double cse_p12_t1_g75 = (cse_p12_t1_g67 * dx);
      double cse_p11_t11_g13 = (cse_p12_t1_g75 * dy);
      double cse_p11_t12_g14 = (cse_p12_t1_g75 * dz);
      double cse_p11_t13_g15 = (cse_p11_t1_g62 * cse_p12_t1_g67);
      double cse_p11_t14_g16 = (cse_p11_t1_g65 * cse_p12_t1_g67);
      double cse_p11_t15_g17 = (cse_p11_t1_g63 * cse_p12_t1_g67);
      double cse_p11_t1_g59 = (cse_p529_t2_invr131 * dr);
      double cse_p11_t16_g18 = (cse_p11_t1_g59 * cse_p11_t1_g64);
      double cse_p12_t1_g73 = (cse_p11_t1_g61 * dr);
      double cse_p11_t17_g19 = (cse_p12_t1_g73 * dx);
      double cse_p11_t18_g20 = (cse_p12_t1_g73 * dy);
      double cse_p11_t19_g21 = (cse_p12_t1_g73 * dz);
      double cse_p11_t20_g22 = (cse_p11_t1_g59 * cse_p11_t1_g62);
      double cse_p11_t21_g23 = (cse_p11_t1_g58 * cse_p11_t1_g62 * dr);
      double cse_p11_t23_g25 = (cse_p11_t1_g59 * cse_p11_t1_g63);
      double cse_p11_t24_g26 = (cse_p11_t1_g58 * cse_p11_t1_g63 * dr);
      double cse_p11_t25_g27 = (cse_p11_t1_g56 * cse_p13_t1_g76);
      double cse_p12_t1_g66 = (cse_p11_t1_g61 * (cse_p529_t1_sqrt130 + cse_p529_t1_sqrt130 + dr + dr));
      double cse_p11_t27_g29 = (cse_p12_t1_g66 * dy);
      double cse_p11_t29_g31 = (cse_p12_t1_g66 * dz);
      double cse_p13_t1_g77 = (cse_p12_t1_g69 * (cse_p529_t1_sqrt130 + cse_p529_t1_sqrt130 + dr + dr));
      double cse_p12_t1_g72 = (cse_p63_t1_g88 * (cse_p11_t16_g18 + cse_p11_t16_g18 + cse_p11_t26_g28 + cse_p11_t26_g28 + dr + dr));
      double cse_p12_t1_g70 = (cse_p63_t1_g88 * (cse_p11_t20_g22 + cse_p11_t20_g22 + cse_p11_t31_g33 + cse_p11_t31_g33 + dr + dr));
      double cse_p12_t1_g71 = (cse_p63_t1_g88 * (cse_p11_t23_g25 + cse_p11_t23_g25 + cse_p11_t34_g36 + cse_p11_t34_g36 + dr + dr));
      double cse_p11_t38_g40 = (cse_p11_t16_g18 + cse_p11_t16_g18 + cse_p11_t26_g28 + cse_p11_t26_g28 + dr + dr);
      double cse_p11_t39_g41 = (cse_p11_t20_g22 + cse_p11_t20_g22 + cse_p11_t31_g33 + cse_p11_t31_g33 + dr + dr);
      double cse_p11_t40_g42 = (cse_p11_t23_g25 + cse_p11_t23_g25 + cse_p11_t34_g36 + cse_p11_t34_g36 + dr + dr);
      double cse_p11_t41_g43 = (cse_p529_t1_sqrt130 + cse_p529_t1_sqrt130 + dr + dr);
      double cse_p11_t42_g44 = (-(cse_p11_t17_g19));
      double cse_p11_t43_g45 = (-(cse_p11_t18_g20));
      double cse_p11_t44_g46 = (-(cse_p11_t19_g21));
      double cse_p11_t45_g47 = (-(cse_p11_t21_g23));
      double cse_p11_t46_g48 = (-(cse_p62_t1_g82));
      double cse_p11_t47_g49 = (-(cse_p11_t24_g26));
      double cse_p11_t51_g53 = (-(cse_p11_t25_g27));
      double energy = (cse_p13_t1_g76 * dr);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p63_t1_g86 = (cse_p13_t1_g76 * cse_p1_t2_invsqrt2);
      double cse_p166_t1_g97 = (2.00000000000000000e+0 * cse_p63_t1_g86);
      double cse_p63_t1_g89 = (-2.00000000000000000e+0 * cse_p63_t1_g86);
      double g_x1 = (cse_p63_t1_g89 * dx);
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = (cse_p63_t1_g89 * dy);
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = (cse_p63_t1_g89 * dz);
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double g_x2 = (cse_p166_t1_g97 * dx);
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (cse_p166_t1_g97 * dy);
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (cse_p166_t1_g97 * dz);
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double h_x1_x1 = (cse_p11_t1_g3 + cse_p11_t1_g3 + cse_p12_t1_g72);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
      double h_x1_y1 = (cse_p11_t27_g29 + cse_p11_t2_g4 + cse_p11_t2_g4);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
      double h_x1_z1 = (cse_p11_t29_g31 + cse_p11_t3_g5 + cse_p11_t3_g5);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
      double cse_p63_t1_g87 = (-2.00000000000000000e+0 * cse_p529_t1_sqrt130 * cse_p529_t1_sqrt130 * cse_p529_t2_invr131 * cse_p529_t2_invr131 * term.kb);
      double h_x1_x2 = (cse_p11_t10_g12 + cse_p11_t10_g12 + cse_p11_t42_g44 + cse_p11_t42_g44 + cse_p11_t51_g53 + cse_p11_t51_g53 + (cse_p11_t1_g64 * cse_p63_t1_g87));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
      double h_x1_y2 = (cse_p114_t1_g92 + cse_p11_t11_g13 + cse_p11_t11_g13 + cse_p11_t43_g45 + cse_p11_t43_g45);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
      double h_x1_z2 = (cse_p11_t12_g14 + cse_p11_t12_g14 + cse_p11_t44_g46 + cse_p11_t44_g46 + cse_p63_t2_g84);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
      double h_y1_y1 = (cse_p11_t4_g6 + cse_p11_t4_g6 + cse_p12_t1_g70);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
      double h_y1_z1 = (cse_p11_t5_g7 + cse_p11_t5_g7 + cse_p13_t1_g77);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
      double h_y1_x2 = h_x1_y2;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
      double h_y1_y2 = (cse_p11_t13_g15 + cse_p11_t13_g15 + cse_p11_t45_g47 + cse_p11_t45_g47 + cse_p11_t51_g53 + cse_p11_t51_g53 + (cse_p11_t1_g62 * cse_p63_t1_g87));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
      double h_y1_z2 = (cse_p11_t14_g16 + cse_p11_t14_g16 + cse_p11_t46_g48 + cse_p11_t46_g48 + cse_p63_t3_g85);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
      double h_z1_z1 = (cse_p11_t6_g8 + cse_p11_t6_g8 + cse_p12_t1_g71);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
      double h_z1_x2 = h_x1_z2;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
      double h_z1_y2 = h_y1_z2;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
      double h_z1_z2 = (cse_p11_t15_g17 + cse_p11_t15_g17 + cse_p11_t47_g49 + cse_p11_t47_g49 + cse_p11_t51_g53 + cse_p11_t51_g53 + (cse_p11_t1_g63 * cse_p63_t1_g87));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
      double h_x2_x2 = h_x1_x1;
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
      double h_x2_y2 = h_x1_y1;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
      double h_x2_z2 = h_x1_z1;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
      double h_y2_y2 = h_y1_y1;
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
      double h_y2_z2 = h_y1_z1;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
      double h_z2_z2 = h_z1_z1;
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
    }
  }
  return energy_added;
}
void hessian_fd(const TermStretch& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 6;
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
}

double gradient(const TermStretch& term, double* position, double* energy_accumulate, double* force) {
  double energy_added = 0.00000000000000000e+0;
  {
    /* !BASE */
    DOUBLE x1 = position[term.i3x1 + 0];
    DOUBLE y1 = position[term.i3x1 + 1];
    DOUBLE z1 = position[term.i3x1 + 2];
    DOUBLE x2 = position[term.i3x2 + 0];
    DOUBLE y2 = position[term.i3x2 + 1];
    DOUBLE z2 = position[term.i3x2 + 2];
    {
      /* !BASE */
      double dx = (x2 + (-(x1)));
      double dy = (y2 + (-(y1)));
      double dz = (z2 + (-(z1)));
      double r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      double cse_p1_t1_sqrt1 = sqrt(r2);
      double cse_p1_t2_invsqrt2 = (1.0 / (cse_p1_t1_sqrt1));
      double r = cse_p1_t1_sqrt1;
      double dr = (r + (-(term.r0)));
      double cse_p11_t1_g5 = (dr * term.kb);
      double energy = (cse_p11_t1_g5 * dr);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t1_g3 = (cse_p11_t1_g5 * cse_p1_t2_invsqrt2);
      double cse_p62_t1_g6 = (-2.00000000000000000e+0 * cse_p11_t1_g3);
      double g_x1 = (cse_p62_t1_g6 * dx);
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = (cse_p62_t1_g6 * dy);
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = (cse_p62_t1_g6 * dz);
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p11_t1_g4 = (2.00000000000000000e+0 * cse_p11_t1_g3);
      double g_x2 = (cse_p11_t1_g4 * dx);
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (cse_p11_t1_g4 * dy);
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (cse_p11_t1_g4 * dz);
      KernelGradientAcc(term.i3x2, 2, g_z2);
    }
  }
  return energy_added;
}
void gradient_fd(const TermStretch& term, double* position, double* energy_accumulate, double* force)
{
  constexpr size_t PositionSize = 6;
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
}

double energy(const TermStretch& term, double* position, double* energy_accumulate) {
  double energy_added = 0.00000000000000000e+0;
  {
    /* !BASE */
    DOUBLE x1 = position[term.i3x1 + 0];
    DOUBLE y1 = position[term.i3x1 + 1];
    DOUBLE z1 = position[term.i3x1 + 2];
    DOUBLE x2 = position[term.i3x2 + 0];
    DOUBLE y2 = position[term.i3x2 + 1];
    DOUBLE z2 = position[term.i3x2 + 2];
    {
      /* !BASE */
      double dx = (x2 + (-(x1)));
      double dy = (y2 + (-(y1)));
      double dz = (z2 + (-(z1)));
      double r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      double r = sqrt(r2);
      double dr = (r + (-(term.r0)));
      double energy = (dr * dr * term.kb);
      energy_added += energy;
      *energy_accumulate += energy;
    }
  }
  return energy_added;
}
void energy_fd(const TermStretch& term, double* position, double* energy_accumulate)
{
  energy(term, position, energy_accumulate);
}

};
