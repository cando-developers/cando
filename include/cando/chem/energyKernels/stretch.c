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
      double cse_p11_t1_g371695 = (dy * dy);
      double cse_p11_t1_g371696 = (dz * dz);
      double cse_p11_t1_g371697 = (dx * dx);
      double r2 = (cse_p11_t1_g371695 + cse_p11_t1_g371696 + cse_p11_t1_g371697);
      double cse_p529_t1_sqrt371763 = sqrt(r2);
      double cse_p529_t2_invr371764 = (1.0 / (r2));
      double cse_p61_t3_invr2371714 = (cse_p529_t2_invr371764 * cse_p529_t2_invr371764);
      double cse_p11_t1_g371689 = (cse_p529_t1_sqrt371763 * cse_p529_t2_invr371764);
      double cse_p63_t1_g371721 = (cse_p11_t1_g371689 * term.kb);
      double cse_p11_t1_g371691 = (cse_p529_t2_invr371764 * cse_p63_t1_g371721);
      double cse_p11_t1_g371694 = (cse_p11_t1_g371691 * dx);
      double cse_p12_t1_g371701 = (cse_p11_t1_g371694 * cse_p529_t1_sqrt371763);
      double cse_p11_t7_g371642 = (-2.00000000000000000e+0 * cse_p12_t1_g371701 * dy);
      double cse_p11_t8_g371643 = (-2.00000000000000000e+0 * cse_p12_t1_g371701 * dz);
      double cse_p11_t1_g371698 = (dy * dz);
      double cse_p12_t1_g371702 = (cse_p11_t1_g371691 * cse_p11_t1_g371698);
      double cse_p13_t1_g371711 = (cse_p12_t1_g371702 * cse_p529_t1_sqrt371763);
      double cse_p11_t9_g371644 = (2.00000000000000000e+0 * (-(cse_p13_t1_g371711)));
      double cse_p11_t26_g371661 = (cse_p11_t1_g371689 * cse_p11_t1_g371697);
      double cse_p11_t28_g371663 = (cse_p12_t1_g371701 * dy);
      double cse_p11_t30_g371665 = (cse_p12_t1_g371701 * dz);
      double cse_p11_t31_g371666 = (cse_p11_t1_g371689 * cse_p11_t1_g371695);
      double cse_p11_t34_g371669 = (cse_p11_t1_g371689 * cse_p11_t1_g371696);
      double cse_p11_t48_g371683 = (-(cse_p11_t28_g371663));
      double cse_p114_t1_g371725 = (2.00000000000000000e+0 * cse_p11_t48_g371683);
      double cse_p11_t49_g371684 = (-(cse_p11_t30_g371665));
      double cse_p63_t2_g371717 = (2.00000000000000000e+0 * cse_p11_t49_g371684);
      double cse_p11_t50_g371685 = (-(cse_p13_t1_g371711));
      double cse_p63_t3_g371718 = (2.00000000000000000e+0 * cse_p11_t50_g371685);
      double cse_p1_t2_invsqrt371635 = (1.0 / (cse_p529_t1_sqrt371763));
      double r = cse_p529_t1_sqrt371763;
      double dr = (r + (-(term.r0)));
      double cse_p62_t1_g371715 = (cse_p12_t1_g371702 * dr);
      double cse_p13_t1_g371709 = (dr * term.kb);
      double cse_p11_t1_g371690 = (cse_p13_t1_g371709 * cse_p529_t1_sqrt371763 * cse_p61_t3_invr2371714);
      double cse_p11_t1_g371693 = (-2.00000000000000000e+0 * cse_p11_t1_g371690);
      double cse_p11_t1_g371636 = (cse_p11_t1_g371693 * cse_p11_t1_g371697);
      double cse_p12_t1_g371707 = (cse_p11_t1_g371693 * dx);
      double cse_p11_t2_g371637 = (cse_p12_t1_g371707 * dy);
      double cse_p11_t3_g371638 = (cse_p12_t1_g371707 * dz);
      double cse_p11_t4_g371639 = (cse_p11_t1_g371693 * cse_p11_t1_g371695);
      double cse_p11_t5_g371640 = (cse_p11_t1_g371693 * cse_p11_t1_g371698);
      double cse_p11_t6_g371641 = (cse_p11_t1_g371693 * cse_p11_t1_g371696);
      double cse_p12_t1_g371700 = (2.00000000000000000e+0 * cse_p11_t1_g371690);
      double cse_p11_t10_g371645 = (cse_p11_t1_g371697 * cse_p12_t1_g371700);
      double cse_p12_t1_g371708 = (cse_p12_t1_g371700 * dx);
      double cse_p11_t11_g371646 = (cse_p12_t1_g371708 * dy);
      double cse_p11_t12_g371647 = (cse_p12_t1_g371708 * dz);
      double cse_p11_t13_g371648 = (cse_p11_t1_g371695 * cse_p12_t1_g371700);
      double cse_p11_t14_g371649 = (cse_p11_t1_g371698 * cse_p12_t1_g371700);
      double cse_p11_t15_g371650 = (cse_p11_t1_g371696 * cse_p12_t1_g371700);
      double cse_p11_t1_g371692 = (cse_p529_t2_invr371764 * dr);
      double cse_p11_t16_g371651 = (cse_p11_t1_g371692 * cse_p11_t1_g371697);
      double cse_p12_t1_g371706 = (cse_p11_t1_g371694 * dr);
      double cse_p11_t17_g371652 = (cse_p12_t1_g371706 * dx);
      double cse_p11_t18_g371653 = (cse_p12_t1_g371706 * dy);
      double cse_p11_t19_g371654 = (cse_p12_t1_g371706 * dz);
      double cse_p11_t20_g371655 = (cse_p11_t1_g371692 * cse_p11_t1_g371695);
      double cse_p11_t21_g371656 = (cse_p11_t1_g371691 * cse_p11_t1_g371695 * dr);
      double cse_p11_t23_g371658 = (cse_p11_t1_g371692 * cse_p11_t1_g371696);
      double cse_p11_t24_g371659 = (cse_p11_t1_g371691 * cse_p11_t1_g371696 * dr);
      double cse_p11_t25_g371660 = (cse_p11_t1_g371689 * cse_p13_t1_g371709);
      double cse_p12_t1_g371699 = (cse_p11_t1_g371694 * (cse_p529_t1_sqrt371763 + cse_p529_t1_sqrt371763 + dr + dr));
      double cse_p11_t27_g371662 = (cse_p12_t1_g371699 * dy);
      double cse_p11_t29_g371664 = (cse_p12_t1_g371699 * dz);
      double cse_p13_t1_g371710 = (cse_p12_t1_g371702 * (cse_p529_t1_sqrt371763 + cse_p529_t1_sqrt371763 + dr + dr));
      double cse_p12_t1_g371705 = (cse_p63_t1_g371721 * (cse_p11_t16_g371651 + cse_p11_t16_g371651 + cse_p11_t26_g371661 + cse_p11_t26_g371661 + dr + dr));
      double cse_p12_t1_g371703 = (cse_p63_t1_g371721 * (cse_p11_t20_g371655 + cse_p11_t20_g371655 + cse_p11_t31_g371666 + cse_p11_t31_g371666 + dr + dr));
      double cse_p12_t1_g371704 = (cse_p63_t1_g371721 * (cse_p11_t23_g371658 + cse_p11_t23_g371658 + cse_p11_t34_g371669 + cse_p11_t34_g371669 + dr + dr));
      double cse_p11_t38_g371673 = (cse_p11_t16_g371651 + cse_p11_t16_g371651 + cse_p11_t26_g371661 + cse_p11_t26_g371661 + dr + dr);
      double cse_p11_t39_g371674 = (cse_p11_t20_g371655 + cse_p11_t20_g371655 + cse_p11_t31_g371666 + cse_p11_t31_g371666 + dr + dr);
      double cse_p11_t40_g371675 = (cse_p11_t23_g371658 + cse_p11_t23_g371658 + cse_p11_t34_g371669 + cse_p11_t34_g371669 + dr + dr);
      double cse_p11_t41_g371676 = (cse_p529_t1_sqrt371763 + cse_p529_t1_sqrt371763 + dr + dr);
      double cse_p11_t42_g371677 = (-(cse_p11_t17_g371652));
      double cse_p11_t43_g371678 = (-(cse_p11_t18_g371653));
      double cse_p11_t44_g371679 = (-(cse_p11_t19_g371654));
      double cse_p11_t45_g371680 = (-(cse_p11_t21_g371656));
      double cse_p11_t46_g371681 = (-(cse_p62_t1_g371715));
      double cse_p11_t47_g371682 = (-(cse_p11_t24_g371659));
      double cse_p11_t51_g371686 = (-(cse_p11_t25_g371660));
      double energy = (cse_p13_t1_g371709 * dr);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p63_t1_g371719 = (cse_p13_t1_g371709 * cse_p1_t2_invsqrt371635);
      double cse_p166_t1_g371730 = (2.00000000000000000e+0 * cse_p63_t1_g371719);
      double cse_p63_t1_g371722 = (-2.00000000000000000e+0 * cse_p63_t1_g371719);
      double g_x1 = (cse_p63_t1_g371722 * dx);
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = (cse_p63_t1_g371722 * dy);
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = (cse_p63_t1_g371722 * dz);
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double g_x2 = (cse_p166_t1_g371730 * dx);
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (cse_p166_t1_g371730 * dy);
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (cse_p166_t1_g371730 * dz);
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double h_x1_x1 = (cse_p11_t1_g371636 + cse_p11_t1_g371636 + cse_p12_t1_g371705);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
      double h_x1_y1 = (cse_p11_t27_g371662 + cse_p11_t2_g371637 + cse_p11_t2_g371637);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
      double h_x1_z1 = (cse_p11_t29_g371664 + cse_p11_t3_g371638 + cse_p11_t3_g371638);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
      double cse_p63_t1_g371720 = (-2.00000000000000000e+0 * cse_p529_t1_sqrt371763 * cse_p529_t1_sqrt371763 * cse_p529_t2_invr371764 * cse_p529_t2_invr371764 * term.kb);
      double h_x1_x2 = (cse_p11_t10_g371645 + cse_p11_t10_g371645 + cse_p11_t42_g371677 + cse_p11_t42_g371677 + cse_p11_t51_g371686 + cse_p11_t51_g371686 + (cse_p11_t1_g371697 * cse_p63_t1_g371720));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
      double h_x1_y2 = (cse_p114_t1_g371725 + cse_p11_t11_g371646 + cse_p11_t11_g371646 + cse_p11_t43_g371678 + cse_p11_t43_g371678);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
      double h_x1_z2 = (cse_p11_t12_g371647 + cse_p11_t12_g371647 + cse_p11_t44_g371679 + cse_p11_t44_g371679 + cse_p63_t2_g371717);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
      double h_y1_y1 = (cse_p11_t4_g371639 + cse_p11_t4_g371639 + cse_p12_t1_g371703);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
      double h_y1_z1 = (cse_p11_t5_g371640 + cse_p11_t5_g371640 + cse_p13_t1_g371710);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
      double h_y1_x2 = h_x1_y2;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
      double h_y1_y2 = (cse_p11_t13_g371648 + cse_p11_t13_g371648 + cse_p11_t45_g371680 + cse_p11_t45_g371680 + cse_p11_t51_g371686 + cse_p11_t51_g371686 + (cse_p11_t1_g371695 * cse_p63_t1_g371720));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
      double h_y1_z2 = (cse_p11_t14_g371649 + cse_p11_t14_g371649 + cse_p11_t46_g371681 + cse_p11_t46_g371681 + cse_p63_t3_g371718);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
      double h_z1_z1 = (cse_p11_t6_g371641 + cse_p11_t6_g371641 + cse_p12_t1_g371704);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
      double h_z1_x2 = h_x1_z2;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
      double h_z1_y2 = h_y1_z2;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
      double h_z1_z2 = (cse_p11_t15_g371650 + cse_p11_t15_g371650 + cse_p11_t47_g371682 + cse_p11_t47_g371682 + cse_p11_t51_g371686 + cse_p11_t51_g371686 + (cse_p11_t1_g371696 * cse_p63_t1_g371720));
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
      double cse_p1_t1_sqrt371766 = sqrt(r2);
      double cse_p1_t2_invsqrt371767 = (1.0 / (cse_p1_t1_sqrt371766));
      double r = cse_p1_t1_sqrt371766;
      double dr = (r + (-(term.r0)));
      double cse_p11_t1_g371770 = (dr * term.kb);
      double energy = (cse_p11_t1_g371770 * dr);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t1_g371768 = (cse_p11_t1_g371770 * cse_p1_t2_invsqrt371767);
      double cse_p62_t1_g371771 = (-2.00000000000000000e+0 * cse_p11_t1_g371768);
      double g_x1 = (cse_p62_t1_g371771 * dx);
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = (cse_p62_t1_g371771 * dy);
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = (cse_p62_t1_g371771 * dz);
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p11_t1_g371769 = (2.00000000000000000e+0 * cse_p11_t1_g371768);
      double g_x2 = (cse_p11_t1_g371769 * dx);
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (cse_p11_t1_g371769 * dy);
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (cse_p11_t1_g371769 * dz);
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
