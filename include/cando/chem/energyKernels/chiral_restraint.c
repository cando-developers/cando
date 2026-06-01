#include "chiral_restraint.h"

template <typename HESSIAN>
struct Chiral_Restraint {
  static constexpr size_t PositionSize = 12;
  static std::string description() { return "mathkernel-chiral_restraint"; };
double energy(const chiral_term& term, double* position, double* energy_accumulate) {
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
    DOUBLE x4 = position[term.i3x4 + 0];
    DOUBLE y4 = position[term.i3x4 + 1];
    DOUBLE z4 = position[term.i3x4 + 2];
    {
      /* !BASE */
      double cse_p11_t1_g1 = (-(x3));
      double cse_p11_t2_g2 = (-(y3));
      double cse_p11_t3_g3 = (-(z3));
      double dx13 = (cse_p11_t1_g1 + x1);
      double dy13 = (cse_p11_t2_g2 + y1);
      double dz13 = (cse_p11_t3_g3 + z1);
      double dx23 = (cse_p11_t1_g1 + x2);
      double dy23 = (cse_p11_t2_g2 + y2);
      double dz23 = (cse_p11_t3_g3 + z2);
      double dx43 = (cse_p11_t1_g1 + x4);
      double dy43 = (cse_p11_t2_g2 + y4);
      double dz43 = (cse_p11_t3_g3 + z4);
      double r13_2 = ((dx13 * dx13) + (dy13 * dy13) + (dz13 * dz13));
      double r23_2 = ((dx23 * dx23) + (dy23 * dy23) + (dz23 * dz23));
      double r43_2 = ((dx43 * dx43) + (dy43 * dy43) + (dz43 * dz43));
      double r13 = sqrt(r13_2);
      double r23 = sqrt(r23_2);
      double r43 = sqrt(r43_2);
      double cx = ((dy13 * dz23) + (-((dy23 * dz13))));
      double cy = ((dx23 * dz13) + (-((dx13 * dz23))));
      double cz = ((dx13 * dy23) + (-((dx23 * dy13))));
      double v = ((cx * dx43) + (cy * dy43) + (cz * dz43));
      double denom = (r13 * r23 * r43);
      double q = (v * (1.0 / (denom)));
      double chiral_e = (term.k * ((term.co + q) * (term.co + q) * (term.co + q)));
      if ((chiral_e > 0.00000000000000000e+0)) {
        /* !ACTIVE */
        double energy = chiral_e;
        energy_added += energy;
        *energy_accumulate += energy;
      }
    }
  }
  return energy_added;
}
void energy_fd(const chiral_term& term, double* position, double* energy_accumulate)
{
  energy(term, position, energy_accumulate);
}

double gradient(const chiral_term& term, double* position, double* energy_accumulate, double* force) {
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
    DOUBLE x4 = position[term.i3x4 + 0];
    DOUBLE y4 = position[term.i3x4 + 1];
    DOUBLE z4 = position[term.i3x4 + 2];
    {
      /* !BASE */
      double cse_p11_t8_g16 = (-(x3));
      double cse_p11_t9_g17 = (-(y3));
      double cse_p11_t10_g18 = (-(z3));
      double dx13 = (cse_p11_t8_g16 + x1);
      double dy13 = (cse_p11_t9_g17 + y1);
      double dz13 = (cse_p11_t10_g18 + z1);
      double dx23 = (cse_p11_t8_g16 + x2);
      double dy23 = (cse_p11_t9_g17 + y2);
      double dz23 = (cse_p11_t10_g18 + z2);
      double dx43 = (cse_p11_t8_g16 + x4);
      double dy43 = (cse_p11_t9_g17 + y4);
      double dz43 = (cse_p11_t10_g18 + z4);
      double r13_2 = ((dx13 * dx13) + (dy13 * dy13) + (dz13 * dz13));
      double r23_2 = ((dx23 * dx23) + (dy23 * dy23) + (dz23 * dz23));
      double r43_2 = ((dx43 * dx43) + (dy43 * dy43) + (dz43 * dz43));
      double cse_p1_t1_sqrt1 = sqrt(r13_2);
      double cse_p1_t2_invsqrt2 = (1.0 / (cse_p1_t1_sqrt1));
      double r13 = cse_p1_t1_sqrt1;
      double cse_p1_t3_sqrt3 = sqrt(r23_2);
      double cse_p1_t4_invsqrt4 = (1.0 / (cse_p1_t3_sqrt3));
      double r23 = cse_p1_t3_sqrt3;
      double cse_p1_t5_sqrt5 = sqrt(r43_2);
      double cse_p1_t6_invsqrt6 = (1.0 / (cse_p1_t5_sqrt5));
      double r43 = cse_p1_t5_sqrt5;
      double cx = ((dy13 * dz23) + (-((dy23 * dz13))));
      double cy = ((dx23 * dz13) + (-((dx13 * dz23))));
      double cz = ((dx13 * dy23) + (-((dx23 * dy13))));
      double v = ((cx * dx43) + (cy * dy43) + (cz * dz43));
      double cse_p11_t1_g28 = (r13 * r23);
      double denom = (cse_p11_t1_g28 * r43);
      double cse_p1_t7_invr7 = (1.0 / (denom));
      double cse_p1_t8_invr28 = (cse_p1_t7_invr7 * cse_p1_t7_invr7);
      double q = (cse_p1_t7_invr7 * v);
      double cse_p11_t1_g9 = (term.co + q);
      double chiral_e = (term.k * (cse_p11_t1_g9 * cse_p11_t1_g9 * cse_p11_t1_g9));
      if ((chiral_e > 0.00000000000000000e+0)) {
        /* !ACTIVE */
        double cse_p11_t2_g10 = (-(dx13));
        double cse_p11_t3_g11 = (-(dx23));
        double cse_p11_t4_g12 = (-(dy13));
        double cse_p11_t5_g13 = (-(dy23));
        double cse_p11_t6_g14 = (-(dz13));
        double cse_p11_t7_g15 = (-(dz23));
        double cse_p11_t11_g19 = ((cse_p11_t1_g9) * (cse_p11_t1_g9));
        double energy = chiral_e;
        energy_added += energy;
        *energy_accumulate += energy;
        double cse_p11_t1_g20 = (3.00000000000000000e+0 * cse_p11_t11_g19 * term.k);
        double cse_p11_t1_g21 = (cse_p1_t2_invsqrt2 * r23 * r43);
        double cse_p62_t1_g29 = (cse_p1_t8_invr28 * v);
        double cse_p11_t1_g25 = (cse_p11_t1_g21 * cse_p62_t1_g29);
        double g_x1 = (cse_p11_t1_g20 * ((cse_p1_t7_invr7 * ((cse_p11_t7_g15 * dy43) + (dy23 * dz43))) + (-((cse_p11_t1_g25 * dx13)))));
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = (cse_p11_t1_g20 * ((cse_p1_t7_invr7 * ((cse_p11_t3_g11 * dz43) + (dx43 * dz23))) + (-((cse_p11_t1_g25 * dy13)))));
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = (cse_p11_t1_g20 * ((cse_p1_t7_invr7 * ((cse_p11_t5_g13 * dx43) + (dx23 * dy43))) + (-((cse_p11_t1_g25 * dz13)))));
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double cse_p11_t1_g23 = (cse_p1_t4_invsqrt4 * r13 * r43);
        double cse_p11_t1_g26 = (cse_p11_t1_g23 * cse_p62_t1_g29);
        double g_x2 = (cse_p11_t1_g20 * ((cse_p1_t7_invr7 * ((cse_p11_t4_g12 * dz43) + (dy43 * dz13))) + (-((cse_p11_t1_g26 * dx23)))));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (cse_p11_t1_g20 * ((cse_p1_t7_invr7 * ((cse_p11_t6_g14 * dx43) + (dx13 * dz43))) + (-((cse_p11_t1_g26 * dy23)))));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (cse_p11_t1_g20 * ((cse_p1_t7_invr7 * ((cse_p11_t2_g10 * dy43) + (dx43 * dy13))) + (-((cse_p11_t1_g26 * dz23)))));
        KernelGradientAcc(term.i3x2, 2, g_z2);
        double cse_p11_t1_g24 = (cse_p11_t1_g28 * cse_p1_t6_invsqrt6);
        double g_x3 = (cse_p11_t1_g20 * ((cse_p1_t7_invr7 * ((dy43 * (cse_p11_t6_g14 + dz23)) + (dz43 * (cse_p11_t5_g13 + dy13)) + (-(cx)))) + (-((cse_p62_t1_g29 * ((-((cse_p11_t1_g21 * dx13))) + (-((cse_p11_t1_g23 * dx23))) + (-((cse_p11_t1_g24 * dx43)))))))));
        KernelGradientAcc(term.i3x3, 0, g_x3);
        double g_y3 = (cse_p11_t1_g20 * ((cse_p1_t7_invr7 * ((dx43 * (cse_p11_t7_g15 + dz13)) + (dz43 * (cse_p11_t2_g10 + dx23)) + (-(cy)))) + (-((cse_p62_t1_g29 * ((-((cse_p11_t1_g21 * dy13))) + (-((cse_p11_t1_g23 * dy23))) + (-((cse_p11_t1_g24 * dy43)))))))));
        KernelGradientAcc(term.i3x3, 1, g_y3);
        double g_z3 = (cse_p11_t1_g20 * ((cse_p1_t7_invr7 * ((dx43 * (cse_p11_t4_g12 + dy23)) + (dy43 * (cse_p11_t3_g11 + dx13)) + (-(cz)))) + (-((cse_p62_t1_g29 * ((-((cse_p11_t1_g21 * dz13))) + (-((cse_p11_t1_g23 * dz23))) + (-((cse_p11_t1_g24 * dz43)))))))));
        KernelGradientAcc(term.i3x3, 2, g_z3);
        double cse_p11_t1_g27 = (cse_p11_t1_g24 * cse_p62_t1_g29);
        double g_x4 = (cse_p11_t1_g20 * ((cse_p1_t7_invr7 * cx) + (-((cse_p11_t1_g27 * dx43)))));
        KernelGradientAcc(term.i3x4, 0, g_x4);
        double g_y4 = (cse_p11_t1_g20 * ((cse_p1_t7_invr7 * cy) + (-((cse_p11_t1_g27 * dy43)))));
        KernelGradientAcc(term.i3x4, 1, g_y4);
        double g_z4 = (cse_p11_t1_g20 * ((cse_p1_t7_invr7 * cz) + (-((cse_p11_t1_g27 * dz43)))));
        KernelGradientAcc(term.i3x4, 2, g_z4);
      }
    }
  }
  return energy_added;
}
void gradient_fd(const chiral_term& term, double* position, double* energy_accumulate, double* force)
{
  constexpr size_t PositionSize = 12;
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
  {
    double saved = position[term.i3x4 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x4 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x4 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x4 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x4, 0, d);
  }
  {
    double saved = position[term.i3x4 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x4 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x4 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x4 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x4, 1, d);
  }
  {
    double saved = position[term.i3x4 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x4 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x4 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x4 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x4, 2, d);
  }
}

double hessian(const chiral_term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
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
    DOUBLE x4 = position[term.i3x4 + 0];
    DOUBLE y4 = position[term.i3x4 + 1];
    DOUBLE z4 = position[term.i3x4 + 2];
    {
      /* !BASE */
      double cse_p11_t189_g197 = (-(x3));
      double cse_p11_t190_g198 = (-(y3));
      double cse_p11_t191_g199 = (-(z3));
      double dx13 = (cse_p11_t189_g197 + x1);
      double dy13 = (cse_p11_t190_g198 + y1);
      double dz13 = (cse_p11_t191_g199 + z1);
      double dx23 = (cse_p11_t189_g197 + x2);
      double dy23 = (cse_p11_t190_g198 + y2);
      double dz23 = (cse_p11_t191_g199 + z2);
      double dx43 = (cse_p11_t189_g197 + x4);
      double dy43 = (cse_p11_t190_g198 + y4);
      double dz43 = (cse_p11_t191_g199 + z4);
      double cse_p15_t1_g610 = (dx13 * dx13);
      double cse_p16_t1_g640 = (dy13 * dy13);
      double cse_p16_t1_g642 = (dz13 * dz13);
      double r13_2 = (cse_p15_t1_g610 + cse_p16_t1_g640 + cse_p16_t1_g642);
      double cse_p16_t1_g639 = (dx23 * dx23);
      double cse_p16_t1_g641 = (dy23 * dy23);
      double cse_p16_t1_g643 = (dz23 * dz23);
      double r23_2 = (cse_p16_t1_g639 + cse_p16_t1_g641 + cse_p16_t1_g643);
      double cse_p17_t1_g662 = (dx43 * dx43);
      double cse_p18_t1_g673 = (dy43 * dy43);
      double cse_p18_t1_g674 = (dz43 * dz43);
      double r43_2 = (cse_p17_t1_g662 + cse_p18_t1_g673 + cse_p18_t1_g674);
      double cse_p520_t1_sqrt830 = sqrt(r13_2);
      double cse_p520_t2_invr831 = (1.0 / (r13_2));
      double cse_p61_t3_invr2770 = (cse_p520_t2_invr831 * cse_p520_t2_invr831);
      double cse_p1_t2_invsqrt2 = (1.0 / (cse_p520_t1_sqrt830));
      double r13 = cse_p520_t1_sqrt830;
      double cse_p520_t3_sqrt832 = sqrt(r23_2);
      double cse_p520_t4_invr833 = (1.0 / (r23_2));
      double cse_p61_t6_invr2773 = (cse_p520_t4_invr833 * cse_p520_t4_invr833);
      double cse_p1_t4_invsqrt4 = (1.0 / (cse_p520_t3_sqrt832));
      double r23 = cse_p520_t3_sqrt832;
      double cse_p520_t5_sqrt834 = sqrt(r43_2);
      double cse_p520_t6_invr835 = (1.0 / (r43_2));
      double cse_p61_t9_invr2776 = (cse_p520_t6_invr835 * cse_p520_t6_invr835);
      double cse_p1_t6_invsqrt6 = (1.0 / (cse_p520_t5_sqrt834));
      double r43 = cse_p520_t5_sqrt834;
      double cse_p21_t1_g706 = (dy13 * dz23);
      double cse_p21_t1_g707 = (dy23 * dz13);
      double cx = (cse_p21_t1_g706 + (-(cse_p21_t1_g707)));
      double cse_p20_t1_g700 = (dx13 * dz23);
      double cse_p21_t1_g703 = (dx23 * dz13);
      double cy = (cse_p21_t1_g703 + (-(cse_p20_t1_g700)));
      double cse_p20_t1_g698 = (dx13 * dy23);
      double cse_p20_t1_g701 = (dx23 * dy13);
      double cz = (cse_p20_t1_g698 + (-(cse_p20_t1_g701)));
      double v = ((cx * dx43) + (cy * dy43) + (cz * dz43));
      double cse_p11_t1_g209 = (r13 * r23);
      double denom = (cse_p11_t1_g209 * r43);
      double cse_p61_t10_invr777 = (1.0 / (denom));
      double cse_p61_t11_invr2778 = (cse_p61_t10_invr777 * cse_p61_t10_invr777);
      double cse_p61_t12_invr3779 = (cse_p61_t10_invr777 * cse_p61_t11_invr2778);
      double q = (cse_p61_t10_invr777 * v);
      double cse_p11_t142_g150 = (term.co + q);
      double chiral_e = (term.k * (cse_p11_t142_g150 * cse_p11_t142_g150 * cse_p11_t142_g150));
      if ((chiral_e > 0.00000000000000000e+0)) {
        /* !ACTIVE */
        double cse_p63_t1_g781 = (cse_p61_t11_invr2778 * v);
        double cse_p16_t19_g632 = (-(dx13));
        double cse_p16_t20_g633 = (-(dx23));
        double cse_p16_t21_g634 = (-(dy13));
        double cse_p16_t22_g635 = (-(dy23));
        double cse_p16_t23_g636 = (-(dz13));
        double cse_p16_t24_g637 = (-(dz23));
        double cse_p15_t13_g579 = (cse_p16_t20_g633 + dx13);
        double cse_p15_t14_g580 = (cse_p16_t19_g632 + dx23);
        double cse_p15_t15_g581 = (cse_p16_t22_g635 + dy13);
        double cse_p15_t16_g582 = (cse_p16_t21_g634 + dy23);
        double cse_p15_t17_g583 = (cse_p16_t24_g637 + dz13);
        double cse_p15_t18_g584 = (cse_p16_t23_g636 + dz23);
        double cse_p14_t14_g511 = (dx13 * dz43);
        double cse_p14_t15_g512 = (dx23 * dy43);
        double cse_p14_t16_g513 = (cse_p15_t16_g582 * dx43);
        double cse_p14_t17_g514 = (cse_p15_t17_g583 * dx43);
        double cse_p14_t18_g515 = (cse_p16_t22_g635 * dx43);
        double cse_p14_t19_g516 = (cse_p16_t23_g636 * dx43);
        double cse_p14_t20_g517 = (dx43 * dy13);
        double cse_p14_t21_g518 = (dx43 * dz23);
        double cse_p14_t22_g519 = (dy23 * dz43);
        double cse_p14_t23_g520 = (cse_p15_t13_g579 * dy43);
        double cse_p14_t24_g521 = (cse_p15_t18_g584 * dy43);
        double cse_p14_t25_g522 = (cse_p16_t19_g632 * dy43);
        double cse_p14_t26_g523 = (cse_p16_t24_g637 * dy43);
        double cse_p14_t27_g524 = (dy43 * dz13);
        double cse_p14_t28_g525 = (cse_p15_t14_g580 * dz43);
        double cse_p14_t29_g526 = (cse_p15_t15_g581 * dz43);
        double cse_p14_t30_g527 = (cse_p16_t20_g633 * dz43);
        double cse_p14_t31_g528 = (cse_p16_t21_g634 * dz43);
        double cse_p14_t50_g547 = (-(cx));
        double cse_p14_t51_g548 = (-(cy));
        double cse_p14_t52_g549 = (-(cz));
        double cse_p13_t43_g432 = (cse_p14_t14_g511 + cse_p14_t19_g516);
        double cse_p13_t44_g433 = (cse_p14_t15_g512 + cse_p14_t18_g515);
        double cse_p13_t45_g434 = (cse_p14_t16_g513 + cse_p14_t23_g520 + cse_p14_t52_g549);
        double cse_p13_t46_g435 = (cse_p14_t17_g514 + cse_p14_t28_g525 + cse_p14_t51_g548);
        double cse_p13_t47_g436 = (cse_p14_t20_g517 + cse_p14_t25_g522);
        double cse_p13_t48_g437 = (cse_p14_t21_g518 + cse_p14_t30_g527);
        double cse_p13_t49_g438 = (cse_p14_t22_g519 + cse_p14_t26_g523);
        double cse_p13_t50_g439 = (cse_p14_t24_g521 + cse_p14_t29_g526 + cse_p14_t50_g547);
        double cse_p13_t51_g440 = (cse_p14_t27_g524 + cse_p14_t31_g528);
        double cse_p21_t1_g712 = (2.00000000000000000e+0 * cse_p61_t11_invr2778);
        double cse_p12_t4_g232 = (cse_p21_t1_g712 * cse_p61_t3_invr2770);
        double cse_p12_t5_g233 = (cse_p21_t1_g712 * cse_p61_t6_invr2773);
        double cse_p12_t6_g234 = (cse_p21_t1_g712 * cse_p61_t9_invr2776);
        double cse_p12_t53_g281 = (cse_p13_t43_g432 * cse_p61_t10_invr777);
        double cse_p12_t54_g282 = (cse_p13_t44_g433 * cse_p61_t10_invr777);
        double cse_p12_t55_g283 = (cse_p13_t45_g434 * cse_p61_t10_invr777);
        double cse_p12_t56_g284 = (cse_p13_t46_g435 * cse_p61_t10_invr777);
        double cse_p12_t57_g285 = (cse_p13_t47_g436 * cse_p61_t10_invr777);
        double cse_p12_t58_g286 = (cse_p13_t48_g437 * cse_p61_t10_invr777);
        double cse_p12_t59_g287 = (cse_p13_t49_g438 * cse_p61_t10_invr777);
        double cse_p12_t60_g288 = (cse_p13_t50_g439 * cse_p61_t10_invr777);
        double cse_p12_t61_g289 = (cse_p13_t51_g440 * cse_p61_t10_invr777);
        double cse_p12_t62_g290 = (cse_p61_t10_invr777 * cx);
        double cse_p12_t63_g291 = (cse_p61_t10_invr777 * cy);
        double cse_p12_t64_g292 = (cse_p61_t10_invr777 * cz);
        double cse_p13_t4_g393 = (cse_p520_t1_sqrt830 * cse_p520_t2_invr831);
        double cse_p13_t5_g394 = (cse_p520_t3_sqrt832 * cse_p520_t4_invr833);
        double cse_p13_t6_g395 = (cse_p520_t5_sqrt834 * cse_p520_t6_invr835);
        double cse_p11_t1_g210 = (r13 * r43);
        double cse_p22_t1_g715 = (cse_p11_t1_g210 * cse_p520_t4_invr833);
        double cse_p12_t28_g256 = (cse_p22_t1_g715 * dx23);
        double cse_p12_t29_g257 = (cse_p22_t1_g715 * dy23);
        double cse_p12_t30_g258 = (cse_p22_t1_g715 * dz23);
        double cse_p11_t1_g211 = (r23 * r43);
        double cse_p22_t1_g716 = (cse_p11_t1_g211 * cse_p520_t2_invr831);
        double cse_p12_t38_g266 = (cse_p22_t1_g716 * dx13);
        double cse_p12_t39_g267 = (cse_p22_t1_g716 * dy13);
        double cse_p12_t40_g268 = (cse_p22_t1_g716 * dz13);
        double cse_p18_t1_g676 = (cse_p11_t1_g209 * cse_p13_t6_g395);
        double cse_p25_t1_g748 = (cse_p18_t1_g676 * dx43);
        double cse_p27_t1_g767 = (cse_p18_t1_g676 * dy43);
        double cse_p17_t3_g650 = (cse_p18_t1_g676 * dz43);
        double cse_p16_t10_g623 = (-(cse_p25_t1_g748));
        double cse_p16_t11_g624 = (-(cse_p27_t1_g767));
        double cse_p16_t12_g625 = (-(cse_p17_t3_g650));
        double cse_p21_t1_g711 = (2.00000000000000000e+0 * cse_p61_t12_invr3779);
        double cse_p26_t1_g761 = (cse_p18_t1_g676 * cse_p520_t6_invr835);
        double cse_p12_t7_g235 = (cse_p21_t1_g711 * cse_p26_t1_g761);
        double cse_p17_t1_g659 = (cse_p520_t6_invr835 * dy43);
        double cse_p12_t10_g238 = (cse_p17_t1_g659 * cse_p25_t1_g748);
        double cse_p15_t1_g604 = (cse_p520_t6_invr835 * dz43);
        double cse_p12_t11_g239 = (cse_p15_t1_g604 * cse_p25_t1_g748);
        double cse_p20_t1_g695 = (cse_p18_t1_g676 * cse_p63_t1_g781);
        double cse_p13_t7_g396 = (cse_p20_t1_g695 * dx43);
        double cse_p13_t8_g397 = (cse_p20_t1_g695 * dy43);
        double cse_p13_t9_g398 = (cse_p20_t1_g695 * dz43);
        double cse_p12_t14_g242 = (cse_p15_t1_g604 * cse_p27_t1_g767);
        double cse_p14_t1_g562 = (cse_p13_t6_g395 * dx43);
        double cse_p12_t41_g269 = (cse_p14_t1_g562 * r13);
        double cse_p12_t42_g270 = (cse_p14_t1_g562 * r23);
        double cse_p14_t1_g563 = (cse_p13_t6_g395 * dy43);
        double cse_p12_t43_g271 = (cse_p14_t1_g563 * r13);
        double cse_p12_t44_g272 = (cse_p14_t1_g563 * r23);
        double cse_p14_t1_g564 = (cse_p13_t6_g395 * dz43);
        double cse_p12_t45_g273 = (cse_p14_t1_g564 * r13);
        double cse_p12_t46_g274 = (cse_p14_t1_g564 * r23);
        double cse_p12_t107_g335 = (-(cse_p13_t7_g396));
        double cse_p12_t108_g336 = (-(cse_p13_t8_g397));
        double cse_p12_t109_g337 = (-(cse_p13_t9_g398));
        double cse_p17_t1_g665 = (cse_p11_t1_g211 * cse_p13_t4_g393);
        double cse_p17_t7_g654 = (cse_p17_t1_g665 * dx13);
        double cse_p17_t8_g655 = (cse_p17_t1_g665 * dy13);
        double cse_p17_t9_g656 = (cse_p17_t1_g665 * dz13);
        double cse_p16_t16_g629 = (-(cse_p17_t7_g654));
        double cse_p16_t17_g630 = (-(cse_p17_t8_g655));
        double cse_p16_t18_g631 = (-(cse_p17_t9_g656));
        double cse_p20_t1_g694 = (cse_p17_t1_g665 * cse_p63_t1_g781);
        double cse_p13_t13_g402 = (cse_p20_t1_g694 * dx13);
        double cse_p13_t14_g403 = (cse_p20_t1_g694 * dy13);
        double cse_p13_t15_g404 = (cse_p20_t1_g694 * dz13);
        double cse_p26_t1_g759 = (cse_p17_t1_g665 * cse_p520_t2_invr831);
        double cse_p12_t9_g237 = (cse_p21_t1_g711 * cse_p26_t1_g759);
        double cse_p14_t1_g566 = (cse_p13_t4_g393 * r43);
        double cse_p12_t47_g275 = (cse_p14_t1_g566 * dx13);
        double cse_p12_t48_g276 = (cse_p14_t1_g566 * dy13);
        double cse_p12_t49_g277 = (cse_p14_t1_g566 * dz13);
        double cse_p12_t113_g341 = (-(cse_p13_t13_g402));
        double cse_p12_t114_g342 = (-(cse_p13_t14_g403));
        double cse_p12_t115_g343 = (-(cse_p13_t15_g404));
        double cse_p17_t1_g664 = (cse_p11_t1_g210 * cse_p13_t5_g394);
        double cse_p17_t4_g651 = (cse_p17_t1_g664 * dx23);
        double cse_p17_t5_g652 = (cse_p17_t1_g664 * dy23);
        double cse_p17_t6_g653 = (cse_p17_t1_g664 * dz23);
        double cse_p16_t13_g626 = (-(cse_p17_t4_g651));
        double cse_p16_t14_g627 = (-(cse_p17_t5_g652));
        double cse_p16_t15_g628 = (-(cse_p17_t6_g653));
        double cse_p15_t10_g576 = (cse_p16_t10_g623 + cse_p16_t13_g626 + cse_p16_t16_g629);
        double cse_p15_t11_g577 = (cse_p16_t11_g624 + cse_p16_t14_g627 + cse_p16_t17_g630);
        double cse_p15_t12_g578 = (cse_p16_t12_g625 + cse_p16_t15_g628 + cse_p16_t18_g631);
        double cse_p14_t1_g498 = (cse_p15_t10_g576 * cse_p61_t11_invr2778);
        double cse_p14_t2_g499 = (cse_p15_t11_g577 * cse_p61_t11_invr2778);
        double cse_p14_t3_g500 = (cse_p15_t12_g578 * cse_p61_t11_invr2778);
        double cse_p20_t1_g693 = (cse_p17_t1_g664 * cse_p63_t1_g781);
        double cse_p13_t10_g399 = (cse_p20_t1_g693 * dx23);
        double cse_p13_t11_g400 = (cse_p20_t1_g693 * dy23);
        double cse_p13_t12_g401 = (cse_p20_t1_g693 * dz23);
        double cse_p13_t61_g450 = (-(cse_p14_t1_g498));
        double cse_p13_t62_g451 = (-(cse_p14_t2_g499));
        double cse_p13_t63_g452 = (-(cse_p14_t3_g500));
        double cse_p26_t1_g760 = (cse_p17_t1_g664 * cse_p520_t4_invr833);
        double cse_p12_t8_g236 = (cse_p21_t1_g711 * cse_p26_t1_g760);
        double cse_p24_t1_g741 = (cse_p13_t5_g394 * r43);
        double cse_p12_t50_g278 = (cse_p24_t1_g741 * dx23);
        double cse_p12_t51_g279 = (cse_p24_t1_g741 * dy23);
        double cse_p12_t52_g280 = (cse_p24_t1_g741 * dz23);
        double cse_p12_t83_g311 = (cse_p13_t61_g450 * v);
        double cse_p12_t84_g312 = (cse_p13_t62_g451 * v);
        double cse_p12_t85_g313 = (cse_p13_t63_g452 * v);
        double cse_p12_t110_g338 = (-(cse_p13_t10_g399));
        double cse_p12_t111_g339 = (-(cse_p13_t11_g400));
        double cse_p12_t112_g340 = (-(cse_p13_t12_g401));
        double cse_p19_t1_g679 = (cse_p11_t1_g211 * cse_p520_t1_sqrt830);
        double cse_p13_t1_g492 = (cse_p19_t1_g679 * cse_p61_t3_invr2770);
        double cse_p14_t1_g559 = (-2.00000000000000000e+0 * cse_p13_t1_g492);
        double cse_p17_t1_g661 = (cse_p14_t1_g559 * cse_p63_t1_g781);
        double cse_p20_t1_g697 = (dx13 * dy13);
        double cse_p11_t4_g12 = (cse_p17_t1_g661 * cse_p20_t1_g697);
        double cse_p20_t1_g699 = (dx13 * dz13);
        double cse_p11_t5_g13 = (cse_p17_t1_g661 * cse_p20_t1_g699);
        double cse_p19_t1_g677 = (cse_p11_t1_g210 * cse_p520_t3_sqrt832);
        double cse_p13_t1_g491 = (cse_p19_t1_g677 * cse_p61_t6_invr2773);
        double cse_p14_t1_g558 = (-2.00000000000000000e+0 * cse_p13_t1_g491);
        double cse_p17_t1_g660 = (cse_p14_t1_g558 * cse_p63_t1_g781);
        double cse_p20_t1_g702 = (dx23 * dy23);
        double cse_p11_t6_g14 = (cse_p17_t1_g660 * cse_p20_t1_g702);
        double cse_p21_t1_g704 = (dx23 * dz23);
        double cse_p11_t7_g15 = (cse_p17_t1_g660 * cse_p21_t1_g704);
        double cse_p21_t1_g705 = (dy13 * dz13);
        double cse_p11_t8_g16 = (cse_p17_t1_g661 * cse_p21_t1_g705);
        double cse_p21_t1_g708 = (dy23 * dz23);
        double cse_p11_t9_g17 = (cse_p17_t1_g660 * cse_p21_t1_g708);
        double cse_p13_t1_g490 = (cse_p11_t1_g209 * cse_p520_t5_sqrt834 * cse_p61_t9_invr2776);
        double cse_p16_t1_g646 = (2.00000000000000000e+0 * cse_p13_t1_g490);
        double cse_p27_t1_g766 = (cse_p16_t1_g646 * dx43);
        double cse_p11_t13_g21 = (cse_p27_t1_g766 * dy43);
        double cse_p11_t14_g22 = (cse_p27_t1_g766 * dz43);
        double cse_p25_t1_g751 = (dy43 * dz43);
        double cse_p11_t15_g23 = (cse_p16_t1_g646 * cse_p25_t1_g751);
        double cse_p15_t1_g609 = (cse_p13_t5_g394 * cse_p14_t1_g566);
        double cse_p25_t1_g749 = (dx13 * dx23);
        double cse_p11_t31_g39 = (cse_p15_t1_g609 * cse_p25_t1_g749);
        double cse_p16_t1_g638 = (cse_p15_t1_g609 * cse_p63_t1_g781);
        double cse_p11_t32_g40 = (cse_p16_t1_g638 * cse_p25_t1_g749);
        double cse_p19_t1_g685 = (cse_p13_t4_g393 * r23);
        double cse_p20_t1_g687 = (cse_p19_t1_g685 * dx13);
        double cse_p20_t2_g688 = (cse_p19_t1_g685 * dy13);
        double cse_p20_t3_g689 = (cse_p19_t1_g685 * dz13);
        double cse_p19_t1_g686 = (cse_p13_t5_g394 * r13);
        double cse_p20_t4_g690 = (cse_p19_t1_g686 * dx23);
        double cse_p20_t5_g691 = (cse_p19_t1_g686 * dy23);
        double cse_p20_t6_g692 = (cse_p19_t1_g686 * dz23);
        double cse_p25_t1_g750 = (dy13 * dy23);
        double cse_p11_t57_g65 = (cse_p15_t1_g609 * cse_p25_t1_g750);
        double cse_p11_t58_g66 = (cse_p16_t1_g638 * cse_p25_t1_g750);
        double cse_p26_t1_g753 = (dz13 * dz23);
        double cse_p11_t80_g88 = (cse_p15_t1_g609 * cse_p26_t1_g753);
        double cse_p11_t81_g89 = (cse_p16_t1_g638 * cse_p26_t1_g753);
        double cse_p11_t109_g117 = (cse_p12_t111_g339 + cse_p12_t53_g281);
        double cse_p11_t110_g118 = (cse_p12_t115_g343 + cse_p12_t54_g282);
        double cse_p11_t111_g119 = (cse_p12_t55_g283 + cse_p12_t85_g313);
        double cse_p11_t112_g120 = (cse_p12_t56_g284 + cse_p12_t84_g312);
        double cse_p11_t113_g121 = (cse_p12_t112_g340 + cse_p12_t57_g285);
        double cse_p11_t114_g122 = (cse_p12_t114_g342 + cse_p12_t58_g286);
        double cse_p11_t115_g123 = (cse_p12_t113_g341 + cse_p12_t59_g287);
        double cse_p11_t116_g124 = (cse_p12_t60_g288 + cse_p12_t83_g311);
        double cse_p11_t117_g125 = (cse_p12_t110_g338 + cse_p12_t61_g289);
        double cse_p11_t118_g126 = (cse_p12_t28_g256 + cse_p12_t41_g269 + cse_p12_t47_g275);
        double cse_p11_t119_g127 = (cse_p12_t29_g257 + cse_p12_t43_g271 + cse_p12_t48_g276);
        double cse_p11_t120_g128 = (cse_p12_t30_g258 + cse_p12_t45_g273 + cse_p12_t49_g277);
        double cse_p11_t124_g132 = (cse_p12_t4_g232 + cse_p12_t9_g237);
        double cse_p11_t125_g133 = (cse_p12_t5_g233 + cse_p12_t8_g236);
        double cse_p11_t126_g134 = (cse_p12_t6_g234 + cse_p12_t7_g235);
        double cse_p11_t127_g135 = (cse_p12_t107_g335 + cse_p12_t62_g290);
        double cse_p11_t128_g136 = (cse_p12_t108_g336 + cse_p12_t63_g291);
        double cse_p11_t129_g137 = (cse_p12_t109_g337 + cse_p12_t64_g292);
        double cse_p11_t131_g139 = (cse_p12_t38_g266 + cse_p12_t42_g270 + cse_p12_t50_g278);
        double cse_p11_t137_g145 = (cse_p12_t39_g267 + cse_p12_t44_g272 + cse_p12_t51_g279);
        double cse_p11_t141_g149 = (cse_p12_t40_g268 + cse_p12_t46_g274 + cse_p12_t52_g280);
        double cse_p11_t156_g164 = (-(cse_p12_t10_g238));
        double cse_p11_t157_g165 = (-(cse_p12_t11_g239));
        double cse_p11_t164_g172 = (-(cse_p12_t14_g242));
        double cse_p11_t173_g181 = (-(cse_p18_t1_g676));
        double cse_p11_t174_g182 = (-(cse_p20_t1_g695));
        double cse_p11_t175_g183 = (-(cse_p20_t1_g693));
        double cse_p11_t176_g184 = (-(cse_p20_t1_g694));
        double cse_p11_t182_g190 = (-(dx43));
        double cse_p11_t185_g193 = (-(dy43));
        double cse_p11_t188_g196 = (-(dz43));
        double cse_p11_t192_g200 = ((cse_p11_t142_g150) * (cse_p11_t142_g150));
        double energy = chiral_e;
        energy_added += energy;
        *energy_accumulate += energy;
        double cse_p11_t1_g218 = (3.00000000000000000e+0 * cse_p11_t192_g200 * term.k);
        double cse_p19_t1_g680 = (cse_p11_t1_g211 * cse_p1_t2_invsqrt2);
        double cse_p15_t1_g612 = (cse_p19_t1_g680 * cse_p63_t1_g781);
        double g_x1 = (cse_p11_t1_g218 * (cse_p12_t59_g287 + (-((cse_p15_t1_g612 * dx13)))));
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = (cse_p11_t1_g218 * (cse_p12_t58_g286 + (-((cse_p15_t1_g612 * dy13)))));
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = (cse_p11_t1_g218 * (cse_p12_t54_g282 + (-((cse_p15_t1_g612 * dz13)))));
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double cse_p19_t1_g678 = (cse_p11_t1_g210 * cse_p1_t4_invsqrt4);
        double cse_p15_t1_g611 = (cse_p19_t1_g678 * cse_p63_t1_g781);
        double g_x2 = (cse_p11_t1_g218 * (cse_p12_t61_g289 + (-((cse_p15_t1_g611 * dx23)))));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (cse_p11_t1_g218 * (cse_p12_t53_g281 + (-((cse_p15_t1_g611 * dy23)))));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (cse_p11_t1_g218 * (cse_p12_t57_g285 + (-((cse_p15_t1_g611 * dz23)))));
        KernelGradientAcc(term.i3x2, 2, g_z2);
        double cse_p15_t1_g608 = (cse_p11_t1_g209 * dx43);
        double cse_p27_t1_g763 = (cse_p15_t1_g608 * cse_p1_t6_invsqrt6);
        double g_x3 = (cse_p11_t1_g218 * (cse_p12_t60_g288 + (-((cse_p63_t1_g781 * ((-((cse_p19_t1_g678 * dx23))) + (-((cse_p19_t1_g680 * dx13))) + (-(cse_p27_t1_g763))))))));
        KernelGradientAcc(term.i3x3, 0, g_x3);
        double cse_p22_t1_g714 = (cse_p11_t1_g209 * cse_p1_t6_invsqrt6);
        double g_y3 = (cse_p11_t1_g218 * (cse_p12_t56_g284 + (-((cse_p63_t1_g781 * ((-((cse_p19_t1_g678 * dy23))) + (-((cse_p19_t1_g680 * dy13))) + (-((cse_p22_t1_g714 * dy43)))))))));
        KernelGradientAcc(term.i3x3, 1, g_y3);
        double g_z3 = (cse_p11_t1_g218 * (cse_p12_t55_g283 + (-((cse_p63_t1_g781 * ((-((cse_p19_t1_g678 * dz23))) + (-((cse_p19_t1_g680 * dz13))) + (-((cse_p22_t1_g714 * dz43)))))))));
        KernelGradientAcc(term.i3x3, 2, g_z3);
        double g_x4 = (cse_p11_t1_g218 * (cse_p12_t62_g290 + (-((cse_p27_t1_g763 * cse_p63_t1_g781)))));
        KernelGradientAcc(term.i3x4, 0, g_x4);
        double cse_p15_t1_g613 = (cse_p22_t1_g714 * cse_p63_t1_g781);
        double g_y4 = (cse_p11_t1_g218 * (cse_p12_t63_g291 + (-((cse_p15_t1_g613 * dy43)))));
        KernelGradientAcc(term.i3x4, 1, g_y4);
        double g_z4 = (cse_p11_t1_g218 * (cse_p12_t64_g292 + (-((cse_p15_t1_g613 * dz43)))));
        KernelGradientAcc(term.i3x4, 2, g_z4);
        double cse_p11_t1_g221 = (cse_p17_t1_g665 * cse_p61_t11_invr2778);
        double cse_p12_t1_g370 = (6.00000000000000000e+0 * cse_p11_t142_g150 * term.k);
        double cse_p12_t1_g371 = (2.00000000000000000e+0 * v);
        double cse_p12_t1_g374 = (cse_p12_t1_g371 * cse_p61_t12_invr3779);
        double cse_p22_t1_g717 = (cse_p12_t1_g371 * cse_p61_t11_invr2778);
        double cse_p12_t1_g386 = (cse_p13_t1_g492 * cse_p22_t1_g717);
        double cse_p13_t1_g482 = (cse_p11_t1_g221 * dx13);
        double cse_p13_t1_g483 = (cse_p12_t1_g374 * r43);
        double cse_p16_t1_g647 = (cse_p11_t115_g123 * cse_p12_t1_g370);
        double cse_p17_t1_g657 = (cse_p520_t2_invr831 * v);
        double cse_p19_t1_g681 = (cse_p13_t1_g483 * cse_p26_t1_g759 * cse_p520_t1_sqrt830 * r23);
        double cse_p24_t1_g738 = (cse_p13_t1_g482 * cse_p17_t1_g657);
        double h_x1_x1 = ((cse_p11_t115_g123 * cse_p16_t1_g647) + (cse_p11_t1_g218 * (cse_p11_t176_g184 + (-2.00000000000000000e+0 * cse_p13_t1_g482 * cse_p13_t49_g438) + (cse_p12_t1_g386 * cse_p15_t1_g610) + (cse_p15_t1_g610 * cse_p19_t1_g681) + (-((cse_p24_t1_g738 * dx13))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double cse_p13_t1_g486 = (cse_p11_t1_g221 * dy13);
        double cse_p15_t1_g607 = (cse_p11_t114_g122 * cse_p12_t1_g370);
        double cse_p21_t1_g709 = (cse_p11_t124_g132 * cse_p19_t1_g679 * v);
        double h_x1_y1 = ((cse_p11_t115_g123 * cse_p15_t1_g607) + (cse_p11_t1_g218 * ((cse_p20_t1_g697 * cse_p21_t1_g709) + (-((cse_p13_t1_g482 * cse_p13_t48_g437))) + (-((cse_p13_t1_g486 * cse_p13_t49_g438))) + (-((cse_p24_t1_g738 * dy13))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double cse_p14_t1_g557 = (cse_p11_t1_g221 * dz13);
        double cse_p14_t1_g561 = (cse_p11_t110_g118 * cse_p12_t1_g370);
        double h_x1_z1 = ((cse_p11_t115_g123 * cse_p14_t1_g561) + (cse_p11_t1_g218 * ((cse_p20_t1_g699 * cse_p21_t1_g709) + (-((cse_p13_t1_g482 * cse_p13_t44_g433))) + (-((cse_p13_t49_g438 * cse_p14_t1_g557))) + (-((cse_p24_t1_g738 * dz13))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double cse_p11_t1_g220 = (cse_p17_t1_g664 * cse_p61_t11_invr2778);
        double cse_p13_t1_g481 = (cse_p11_t1_g220 * dx23);
        double cse_p13_t1_g489 = (cse_p11_t1_g209 * cse_p13_t1_g483 * cse_p15_t1_g609);
        double h_x1_x2 = ((cse_p11_t117_g125 * cse_p16_t1_g647) + (cse_p11_t1_g218 * ((cse_p13_t1_g489 * cse_p25_t1_g749) + (-((cse_p13_t1_g481 * cse_p13_t49_g438))) + (-((cse_p13_t1_g482 * cse_p13_t51_g440))) + (-(cse_p11_t32_g40)))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double cse_p13_t1_g485 = (cse_p11_t1_g220 * dy23);
        double cse_p14_t1_g560 = (cse_p11_t109_g117 * cse_p12_t1_g370);
        double h_x1_y2 = ((cse_p11_t115_g123 * cse_p14_t1_g560) + (cse_p11_t1_g218 * ((cse_p13_t1_g489 * cse_p20_t1_g698) + (cse_p61_t10_invr777 * dz43) + (-((cse_p13_t1_g482 * cse_p13_t43_g432))) + (-((cse_p13_t1_g485 * cse_p13_t49_g438))) + (-((cse_p16_t1_g638 * cse_p20_t1_g698))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double cse_p13_t1_g488 = (cse_p11_t1_g220 * dz23);
        double cse_p15_t1_g606 = (cse_p11_t113_g121 * cse_p12_t1_g370);
        double h_x1_z2 = ((cse_p11_t115_g123 * cse_p15_t1_g606) + (cse_p11_t1_g218 * ((cse_p11_t185_g193 * cse_p61_t10_invr777) + (cse_p13_t1_g489 * cse_p20_t1_g700) + (-((cse_p13_t1_g482 * cse_p13_t47_g436))) + (-((cse_p13_t1_g488 * cse_p13_t49_g438))) + (-((cse_p16_t1_g638 * cse_p20_t1_g700))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double cse_p12_t1_g383 = (cse_p12_t1_g374 * cse_p17_t1_g665);
        double cse_p19_t1_g683 = (cse_p19_t1_g685 * cse_p63_t1_g781);
        double cse_p23_t1_g727 = (cse_p12_t1_g383 * cse_p15_t10_g576);
        double cse_p24_t1_g742 = (cse_p520_t2_invr831 * r43);
        double h_x1_x3 = ((cse_p11_t116_g124 * cse_p16_t1_g647) + (cse_p11_t1_g218 * (cse_p11_t32_g40 + (cse_p13_t49_g438 * cse_p13_t61_g450) + (cse_p15_t1_g610 * cse_p17_t1_g661) + (cse_p19_t1_g683 * (r43 + (cse_p14_t1_g562 * dx13) + (cse_p15_t1_g610 * cse_p24_t1_g742))) + (cse_p23_t1_g727 * dx13) + (-((cse_p13_t1_g482 * cse_p13_t50_g439))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 0, h_x1_x3);
        double cse_p24_t1_g739 = (cse_p13_t4_g393 * cse_p63_t1_g781);
        double cse_p13_t1_g496 = (cse_p11_t137_g145 * cse_p24_t1_g739);
        double cse_p15_t1_g603 = (cse_p11_t112_g120 * cse_p12_t1_g370);
        double cse_p23_t1_g728 = (cse_p12_t1_g383 * cse_p15_t11_g577);
        double h_x1_y3 = ((cse_p11_t115_g123 * cse_p15_t1_g603) + (cse_p11_t1_g218 * (cse_p11_t4_g12 + (cse_p13_t1_g496 * dx13) + (cse_p13_t49_g438 * cse_p13_t62_g451) + (cse_p23_t1_g728 * dx13) + (cse_p61_t10_invr777 * (cse_p11_t188_g196 + dz23)) + (-((cse_p13_t1_g482 * cse_p13_t46_g435))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 1, h_x1_y3);
        double cse_p13_t1_g497 = (cse_p11_t141_g149 * cse_p24_t1_g739);
        double cse_p14_t1_g565 = (cse_p11_t111_g119 * cse_p12_t1_g370);
        double cse_p23_t1_g729 = (cse_p12_t1_g383 * cse_p15_t12_g578);
        double h_x1_z3 = ((cse_p11_t115_g123 * cse_p14_t1_g565) + (cse_p11_t1_g218 * (cse_p11_t5_g13 + (cse_p13_t1_g497 * dx13) + (cse_p13_t49_g438 * cse_p13_t63_g452) + (cse_p23_t1_g729 * dx13) + (cse_p61_t10_invr777 * (cse_p16_t22_g635 + dy43)) + (-((cse_p13_t1_g482 * cse_p13_t45_g434))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 2, h_x1_z3);
        double cse_p11_t1_g219 = (cse_p18_t1_g676 * cse_p61_t11_invr2778);
        double cse_p11_t1_g222 = (cse_p18_t1_g676 * cse_p61_t12_invr3779);
        double cse_p17_t1_g666 = (cse_p11_t1_g222 * cse_p12_t1_g371);
        double cse_p12_t1_g381 = (cse_p17_t1_g665 * cse_p17_t1_g666);
        double cse_p12_t1_g384 = (cse_p13_t6_g395 * cse_p19_t1_g683);
        double cse_p13_t1_g480 = (cse_p11_t1_g219 * dx43);
        double cse_p22_t1_g721 = (cse_p12_t1_g381 * dx13);
        double cse_p23_t1_g730 = (cse_p12_t1_g384 * dx13);
        double h_x1_x4 = ((cse_p11_t127_g135 * cse_p16_t1_g647) + (cse_p11_t1_g218 * ((cse_p22_t1_g721 * dx43) + (-((cse_p13_t1_g480 * cse_p13_t49_g438))) + (-((cse_p13_t1_g482 * cx))) + (-((cse_p23_t1_g730 * dx43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x4, 0, h_x1_x4);
        double cse_p13_t1_g484 = (cse_p11_t1_g219 * dy43);
        double h_x1_y4 = ((cse_p11_t128_g136 * cse_p16_t1_g647) + (cse_p11_t1_g218 * ((cse_p16_t24_g637 * cse_p61_t10_invr777) + (cse_p22_t1_g721 * dy43) + (-((cse_p13_t1_g482 * cy))) + (-((cse_p13_t1_g484 * cse_p13_t49_g438))) + (-((cse_p23_t1_g730 * dy43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x4, 1, h_x1_y4);
        double cse_p13_t1_g487 = (cse_p11_t1_g219 * dz43);
        double h_x1_z4 = ((cse_p11_t129_g137 * cse_p16_t1_g647) + (cse_p11_t1_g218 * ((cse_p22_t1_g721 * dz43) + (cse_p61_t10_invr777 * dy23) + (-((cse_p13_t1_g482 * cz))) + (-((cse_p13_t1_g487 * cse_p13_t49_g438))) + (-((cse_p23_t1_g730 * dz43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x4, 2, h_x1_z4);
        double cse_p26_t1_g758 = (cse_p13_t1_g486 * cse_p17_t1_g657);
        double h_y1_y1 = ((cse_p11_t114_g122 * cse_p15_t1_g607) + (cse_p11_t1_g218 * (cse_p11_t176_g184 + (-2.00000000000000000e+0 * cse_p13_t1_g486 * cse_p13_t48_g437) + (cse_p12_t1_g386 * cse_p16_t1_g640) + (cse_p16_t1_g640 * cse_p19_t1_g681) + (-((cse_p26_t1_g758 * dy13))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = ((cse_p11_t114_g122 * cse_p14_t1_g561) + (cse_p11_t1_g218 * ((cse_p21_t1_g705 * cse_p21_t1_g709) + (-((cse_p13_t1_g486 * cse_p13_t44_g433))) + (-((cse_p13_t48_g437 * cse_p14_t1_g557))) + (-((cse_p26_t1_g758 * dz13))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = ((cse_p11_t117_g125 * cse_p15_t1_g607) + (cse_p11_t1_g218 * ((cse_p11_t188_g196 * cse_p61_t10_invr777) + (cse_p13_t1_g489 * cse_p20_t1_g701) + (-((cse_p13_t1_g481 * cse_p13_t48_g437))) + (-((cse_p13_t1_g486 * cse_p13_t51_g440))) + (-((cse_p16_t1_g638 * cse_p20_t1_g701))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = ((cse_p11_t114_g122 * cse_p14_t1_g560) + (cse_p11_t1_g218 * ((cse_p13_t1_g489 * cse_p25_t1_g750) + (-((cse_p13_t1_g485 * cse_p13_t48_g437))) + (-((cse_p13_t1_g486 * cse_p13_t43_g432))) + (-(cse_p11_t58_g66)))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = ((cse_p11_t114_g122 * cse_p15_t1_g606) + (cse_p11_t1_g218 * ((cse_p13_t1_g489 * cse_p21_t1_g706) + (cse_p61_t10_invr777 * dx43) + (-((cse_p13_t1_g486 * cse_p13_t47_g436))) + (-((cse_p13_t1_g488 * cse_p13_t48_g437))) + (-((cse_p16_t1_g638 * cse_p21_t1_g706))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double cse_p13_t1_g495 = (cse_p11_t131_g139 * cse_p24_t1_g739);
        double h_y1_x3 = ((cse_p11_t116_g124 * cse_p15_t1_g607) + (cse_p11_t1_g218 * (cse_p11_t4_g12 + (cse_p13_t1_g495 * dy13) + (cse_p13_t48_g437 * cse_p13_t61_g450) + (cse_p23_t1_g727 * dy13) + (cse_p61_t10_invr777 * (cse_p16_t24_g637 + dz43)) + (-((cse_p13_t1_g486 * cse_p13_t50_g439))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 0, h_y1_x3);
        double h_y1_y3 = ((cse_p11_t114_g122 * cse_p15_t1_g603) + (cse_p11_t1_g218 * (cse_p11_t58_g66 + (cse_p13_t48_g437 * cse_p13_t62_g451) + (cse_p16_t1_g640 * cse_p17_t1_g661) + (cse_p19_t1_g683 * (r43 + (cse_p14_t1_g563 * dy13) + (cse_p16_t1_g640 * cse_p24_t1_g742))) + (cse_p23_t1_g728 * dy13) + (-((cse_p13_t1_g486 * cse_p13_t46_g435))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 1, h_y1_y3);
        double h_y1_z3 = ((cse_p11_t114_g122 * cse_p14_t1_g565) + (cse_p11_t1_g218 * (cse_p11_t8_g16 + (cse_p13_t1_g497 * dy13) + (cse_p13_t48_g437 * cse_p13_t63_g452) + (cse_p23_t1_g729 * dy13) + (cse_p61_t10_invr777 * (cse_p11_t182_g190 + dx23)) + (-((cse_p13_t1_g486 * cse_p13_t45_g434))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 2, h_y1_z3);
        double cse_p22_t1_g722 = (cse_p12_t1_g381 * dy13);
        double cse_p23_t1_g731 = (cse_p12_t1_g384 * dy13);
        double h_y1_x4 = ((cse_p11_t127_g135 * cse_p15_t1_g607) + (cse_p11_t1_g218 * ((cse_p22_t1_g722 * dx43) + (cse_p61_t10_invr777 * dz23) + (-((cse_p13_t1_g480 * cse_p13_t48_g437))) + (-((cse_p13_t1_g486 * cx))) + (-((cse_p23_t1_g731 * dx43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x4, 0, h_y1_x4);
        double h_y1_y4 = ((cse_p11_t128_g136 * cse_p15_t1_g607) + (cse_p11_t1_g218 * ((cse_p22_t1_g722 * dy43) + (-((cse_p13_t1_g484 * cse_p13_t48_g437))) + (-((cse_p13_t1_g486 * cy))) + (-((cse_p23_t1_g731 * dy43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x4, 1, h_y1_y4);
        double h_y1_z4 = ((cse_p11_t129_g137 * cse_p15_t1_g607) + (cse_p11_t1_g218 * ((cse_p16_t20_g633 * cse_p61_t10_invr777) + (cse_p22_t1_g722 * dz43) + (-((cse_p13_t1_g486 * cz))) + (-((cse_p13_t1_g487 * cse_p13_t48_g437))) + (-((cse_p23_t1_g731 * dz43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x4, 2, h_y1_z4);
        double h_z1_z1 = ((cse_p11_t110_g118 * cse_p14_t1_g561) + (cse_p11_t1_g218 * (cse_p11_t176_g184 + (-2.00000000000000000e+0 * cse_p13_t44_g433 * cse_p14_t1_g557) + (cse_p12_t1_g386 * cse_p16_t1_g642) + (cse_p16_t1_g642 * cse_p19_t1_g681) + (-((cse_p14_t1_g557 * cse_p17_t1_g657 * dz13))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = ((cse_p11_t117_g125 * cse_p14_t1_g561) + (cse_p11_t1_g218 * ((cse_p13_t1_g489 * cse_p21_t1_g703) + (cse_p61_t10_invr777 * dy43) + (-((cse_p13_t1_g481 * cse_p13_t44_g433))) + (-((cse_p13_t51_g440 * cse_p14_t1_g557))) + (-((cse_p16_t1_g638 * cse_p21_t1_g703))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = ((cse_p11_t110_g118 * cse_p14_t1_g560) + (cse_p11_t1_g218 * ((cse_p11_t182_g190 * cse_p61_t10_invr777) + (cse_p13_t1_g489 * cse_p21_t1_g707) + (-((cse_p13_t1_g485 * cse_p13_t44_g433))) + (-((cse_p13_t43_g432 * cse_p14_t1_g557))) + (-((cse_p16_t1_g638 * cse_p21_t1_g707))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = ((cse_p11_t113_g121 * cse_p14_t1_g561) + (cse_p11_t1_g218 * ((cse_p13_t1_g489 * cse_p26_t1_g753) + (-((cse_p13_t1_g488 * cse_p13_t44_g433))) + (-((cse_p13_t47_g436 * cse_p14_t1_g557))) + (-(cse_p11_t81_g89)))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_z1_x3 = ((cse_p11_t116_g124 * cse_p14_t1_g561) + (cse_p11_t1_g218 * (cse_p11_t5_g13 + (cse_p13_t1_g495 * dz13) + (cse_p13_t44_g433 * cse_p13_t61_g450) + (cse_p23_t1_g727 * dz13) + (cse_p61_t10_invr777 * (cse_p11_t185_g193 + dy23)) + (-((cse_p13_t50_g439 * cse_p14_t1_g557))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 0, h_z1_x3);
        double h_z1_y3 = ((cse_p11_t112_g120 * cse_p14_t1_g561) + (cse_p11_t1_g218 * (cse_p11_t8_g16 + (cse_p13_t1_g496 * dz13) + (cse_p13_t44_g433 * cse_p13_t62_g451) + (cse_p23_t1_g728 * dz13) + (cse_p61_t10_invr777 * (cse_p16_t20_g633 + dx43)) + (-((cse_p13_t46_g435 * cse_p14_t1_g557))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 1, h_z1_y3);
        double h_z1_z3 = ((cse_p11_t111_g119 * cse_p14_t1_g561) + (cse_p11_t1_g218 * (cse_p11_t81_g89 + (cse_p13_t44_g433 * cse_p13_t63_g452) + (cse_p16_t1_g642 * cse_p17_t1_g661) + (cse_p19_t1_g683 * (r43 + (cse_p14_t1_g564 * dz13) + (cse_p16_t1_g642 * cse_p24_t1_g742))) + (cse_p23_t1_g729 * dz13) + (-((cse_p13_t45_g434 * cse_p14_t1_g557))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 2, h_z1_z3);
        double cse_p23_t1_g723 = (cse_p12_t1_g381 * dz13);
        double cse_p23_t1_g732 = (cse_p12_t1_g384 * dz13);
        double h_z1_x4 = ((cse_p11_t127_g135 * cse_p14_t1_g561) + (cse_p11_t1_g218 * ((cse_p16_t22_g635 * cse_p61_t10_invr777) + (cse_p23_t1_g723 * dx43) + (-((cse_p13_t1_g480 * cse_p13_t44_g433))) + (-((cse_p14_t1_g557 * cx))) + (-((cse_p23_t1_g732 * dx43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x4, 0, h_z1_x4);
        double h_z1_y4 = ((cse_p11_t128_g136 * cse_p14_t1_g561) + (cse_p11_t1_g218 * ((cse_p23_t1_g723 * dy43) + (cse_p61_t10_invr777 * dx23) + (-((cse_p13_t1_g484 * cse_p13_t44_g433))) + (-((cse_p14_t1_g557 * cy))) + (-((cse_p23_t1_g732 * dy43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x4, 1, h_z1_y4);
        double h_z1_z4 = ((cse_p11_t129_g137 * cse_p14_t1_g561) + (cse_p11_t1_g218 * ((cse_p23_t1_g723 * dz43) + (-((cse_p13_t1_g487 * cse_p13_t44_g433))) + (-((cse_p14_t1_g557 * cz))) + (-((cse_p23_t1_g732 * dz43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x4, 2, h_z1_z4);
        double cse_p12_t1_g387 = (cse_p13_t1_g491 * cse_p22_t1_g717);
        double cse_p17_t1_g658 = (cse_p520_t4_invr833 * v);
        double cse_p18_t1_g675 = (cse_p11_t117_g125 * cse_p12_t1_g370);
        double cse_p19_t1_g682 = (cse_p13_t1_g483 * cse_p26_t1_g760 * cse_p520_t3_sqrt832 * r13);
        double cse_p24_t1_g737 = (cse_p13_t1_g481 * cse_p17_t1_g658);
        double h_x2_x2 = ((cse_p11_t117_g125 * cse_p18_t1_g675) + (cse_p11_t1_g218 * (cse_p11_t175_g183 + (-2.00000000000000000e+0 * cse_p13_t1_g481 * cse_p13_t51_g440) + (cse_p12_t1_g387 * cse_p16_t1_g639) + (cse_p16_t1_g639 * cse_p19_t1_g682) + (-((cse_p24_t1_g737 * dx23))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double cse_p21_t1_g710 = (cse_p11_t125_g133 * cse_p19_t1_g677 * v);
        double h_x2_y2 = ((cse_p11_t117_g125 * cse_p14_t1_g560) + (cse_p11_t1_g218 * ((cse_p20_t1_g702 * cse_p21_t1_g710) + (-((cse_p13_t1_g481 * cse_p13_t43_g432))) + (-((cse_p13_t1_g485 * cse_p13_t51_g440))) + (-((cse_p24_t1_g737 * dy23))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = ((cse_p11_t117_g125 * cse_p15_t1_g606) + (cse_p11_t1_g218 * ((cse_p21_t1_g704 * cse_p21_t1_g710) + (-((cse_p13_t1_g481 * cse_p13_t47_g436))) + (-((cse_p13_t1_g488 * cse_p13_t51_g440))) + (-((cse_p24_t1_g737 * dz23))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double cse_p12_t1_g382 = (cse_p12_t1_g374 * cse_p17_t1_g664);
        double cse_p17_t1_g663 = (cse_p11_t116_g124 * cse_p12_t1_g370);
        double cse_p19_t1_g684 = (cse_p19_t1_g686 * cse_p63_t1_g781);
        double cse_p23_t1_g724 = (cse_p12_t1_g382 * cse_p15_t10_g576);
        double cse_p25_t1_g743 = (cse_p520_t4_invr833 * r43);
        double h_x2_x3 = ((cse_p11_t117_g125 * cse_p17_t1_g663) + (cse_p11_t1_g218 * (cse_p11_t32_g40 + (cse_p13_t51_g440 * cse_p13_t61_g450) + (cse_p16_t1_g639 * cse_p17_t1_g660) + (cse_p19_t1_g684 * (r43 + (cse_p14_t1_g562 * dx23) + (cse_p16_t1_g639 * cse_p25_t1_g743))) + (cse_p23_t1_g724 * dx23) + (-((cse_p13_t1_g481 * cse_p13_t50_g439))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 0, h_x2_x3);
        double cse_p24_t1_g740 = (cse_p13_t5_g394 * cse_p63_t1_g781);
        double cse_p13_t1_g493 = (cse_p11_t119_g127 * cse_p24_t1_g740);
        double cse_p23_t1_g725 = (cse_p12_t1_g382 * cse_p15_t11_g577);
        double h_x2_y3 = ((cse_p11_t117_g125 * cse_p15_t1_g603) + (cse_p11_t1_g218 * (cse_p11_t6_g14 + (cse_p13_t1_g493 * dx23) + (cse_p13_t51_g440 * cse_p13_t62_g451) + (cse_p23_t1_g725 * dx23) + (cse_p61_t10_invr777 * (cse_p16_t23_g636 + dz43)) + (-((cse_p13_t1_g481 * cse_p13_t46_g435))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 1, h_x2_y3);
        double cse_p13_t1_g494 = (cse_p11_t120_g128 * cse_p24_t1_g740);
        double cse_p23_t1_g726 = (cse_p12_t1_g382 * cse_p15_t12_g578);
        double h_x2_z3 = ((cse_p11_t117_g125 * cse_p14_t1_g565) + (cse_p11_t1_g218 * (cse_p11_t7_g15 + (cse_p13_t1_g494 * dx23) + (cse_p13_t51_g440 * cse_p13_t63_g452) + (cse_p23_t1_g726 * dx23) + (cse_p61_t10_invr777 * (cse_p11_t185_g193 + dy13)) + (-((cse_p13_t1_g481 * cse_p13_t45_g434))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 2, h_x2_z3);
        double cse_p12_t1_g380 = (cse_p17_t1_g664 * cse_p17_t1_g666);
        double cse_p12_t1_g385 = (cse_p13_t6_g395 * cse_p19_t1_g684);
        double cse_p22_t1_g718 = (cse_p12_t1_g380 * dx23);
        double cse_p24_t1_g733 = (cse_p12_t1_g385 * dx23);
        double h_x2_x4 = ((cse_p11_t127_g135 * cse_p18_t1_g675) + (cse_p11_t1_g218 * ((cse_p22_t1_g718 * dx43) + (-((cse_p13_t1_g480 * cse_p13_t51_g440))) + (-((cse_p13_t1_g481 * cx))) + (-((cse_p24_t1_g733 * dx43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x4, 0, h_x2_x4);
        double h_x2_y4 = ((cse_p11_t128_g136 * cse_p18_t1_g675) + (cse_p11_t1_g218 * ((cse_p22_t1_g718 * dy43) + (cse_p61_t10_invr777 * dz13) + (-((cse_p13_t1_g481 * cy))) + (-((cse_p13_t1_g484 * cse_p13_t51_g440))) + (-((cse_p24_t1_g733 * dy43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x4, 1, h_x2_y4);
        double h_x2_z4 = ((cse_p11_t129_g137 * cse_p18_t1_g675) + (cse_p11_t1_g218 * ((cse_p16_t21_g634 * cse_p61_t10_invr777) + (cse_p22_t1_g718 * dz43) + (-((cse_p13_t1_g481 * cz))) + (-((cse_p13_t1_g487 * cse_p13_t51_g440))) + (-((cse_p24_t1_g733 * dz43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x4, 2, h_x2_z4);
        double cse_p26_t1_g757 = (cse_p13_t1_g485 * cse_p17_t1_g658);
        double h_y2_y2 = ((cse_p11_t109_g117 * cse_p14_t1_g560) + (cse_p11_t1_g218 * (cse_p11_t175_g183 + (-2.00000000000000000e+0 * cse_p13_t1_g485 * cse_p13_t43_g432) + (cse_p12_t1_g387 * cse_p16_t1_g641) + (cse_p16_t1_g641 * cse_p19_t1_g682) + (-((cse_p26_t1_g757 * dy23))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = ((cse_p11_t113_g121 * cse_p14_t1_g560) + (cse_p11_t1_g218 * ((cse_p21_t1_g708 * cse_p21_t1_g710) + (-((cse_p13_t1_g485 * cse_p13_t47_g436))) + (-((cse_p13_t1_g488 * cse_p13_t43_g432))) + (-((cse_p26_t1_g757 * dz23))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double cse_p12_t1_g389 = (cse_p11_t118_g126 * cse_p24_t1_g740);
        double h_y2_x3 = ((cse_p11_t116_g124 * cse_p14_t1_g560) + (cse_p11_t1_g218 * (cse_p11_t6_g14 + (cse_p12_t1_g389 * dy23) + (cse_p13_t43_g432 * cse_p13_t61_g450) + (cse_p23_t1_g724 * dy23) + (cse_p61_t10_invr777 * (cse_p11_t188_g196 + dz13)) + (-((cse_p13_t1_g485 * cse_p13_t50_g439))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 0, h_y2_x3);
        double h_y2_y3 = ((cse_p11_t112_g120 * cse_p14_t1_g560) + (cse_p11_t1_g218 * (cse_p11_t58_g66 + (cse_p13_t43_g432 * cse_p13_t62_g451) + (cse_p16_t1_g641 * cse_p17_t1_g660) + (cse_p19_t1_g684 * (r43 + (cse_p14_t1_g563 * dy23) + (cse_p16_t1_g641 * cse_p25_t1_g743))) + (cse_p23_t1_g725 * dy23) + (-((cse_p13_t1_g485 * cse_p13_t46_g435))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 1, h_y2_y3);
        double h_y2_z3 = ((cse_p11_t111_g119 * cse_p14_t1_g560) + (cse_p11_t1_g218 * (cse_p11_t9_g17 + (cse_p13_t1_g494 * dy23) + (cse_p13_t43_g432 * cse_p13_t63_g452) + (cse_p23_t1_g726 * dy23) + (cse_p61_t10_invr777 * (cse_p16_t19_g632 + dx43)) + (-((cse_p13_t1_g485 * cse_p13_t45_g434))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 2, h_y2_z3);
        double cse_p22_t1_g719 = (cse_p12_t1_g380 * dy23);
        double cse_p24_t1_g734 = (cse_p12_t1_g385 * dy23);
        double h_y2_x4 = ((cse_p11_t127_g135 * cse_p14_t1_g560) + (cse_p11_t1_g218 * ((cse_p16_t23_g636 * cse_p61_t10_invr777) + (cse_p22_t1_g719 * dx43) + (-((cse_p13_t1_g480 * cse_p13_t43_g432))) + (-((cse_p13_t1_g485 * cx))) + (-((cse_p24_t1_g734 * dx43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x4, 0, h_y2_x4);
        double h_y2_y4 = ((cse_p11_t128_g136 * cse_p14_t1_g560) + (cse_p11_t1_g218 * ((cse_p22_t1_g719 * dy43) + (-((cse_p13_t1_g484 * cse_p13_t43_g432))) + (-((cse_p13_t1_g485 * cy))) + (-((cse_p24_t1_g734 * dy43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x4, 1, h_y2_y4);
        double h_y2_z4 = ((cse_p11_t129_g137 * cse_p14_t1_g560) + (cse_p11_t1_g218 * ((cse_p22_t1_g719 * dz43) + (cse_p61_t10_invr777 * dx13) + (-((cse_p13_t1_g485 * cz))) + (-((cse_p13_t1_g487 * cse_p13_t43_g432))) + (-((cse_p24_t1_g734 * dz43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x4, 2, h_y2_z4);
        double h_z2_z2 = ((cse_p11_t113_g121 * cse_p15_t1_g606) + (cse_p11_t1_g218 * (cse_p11_t175_g183 + (-2.00000000000000000e+0 * cse_p13_t1_g488 * cse_p13_t47_g436) + (cse_p12_t1_g387 * cse_p16_t1_g643) + (cse_p16_t1_g643 * cse_p19_t1_g682) + (-((cse_p13_t1_g488 * cse_p17_t1_g658 * dz23))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        double h_z2_x3 = ((cse_p11_t116_g124 * cse_p15_t1_g606) + (cse_p11_t1_g218 * (cse_p11_t7_g15 + (cse_p12_t1_g389 * dz23) + (cse_p13_t47_g436 * cse_p13_t61_g450) + (cse_p23_t1_g724 * dz23) + (cse_p61_t10_invr777 * (cse_p16_t21_g634 + dy43)) + (-((cse_p13_t1_g488 * cse_p13_t50_g439))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 0, h_z2_x3);
        double h_z2_y3 = ((cse_p11_t113_g121 * cse_p15_t1_g603) + (cse_p11_t1_g218 * (cse_p11_t9_g17 + (cse_p13_t1_g493 * dz23) + (cse_p13_t47_g436 * cse_p13_t62_g451) + (cse_p23_t1_g725 * dz23) + (cse_p61_t10_invr777 * (cse_p11_t182_g190 + dx13)) + (-((cse_p13_t1_g488 * cse_p13_t46_g435))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 1, h_z2_y3);
        double h_z2_z3 = ((cse_p11_t113_g121 * cse_p14_t1_g565) + (cse_p11_t1_g218 * (cse_p11_t81_g89 + (cse_p13_t47_g436 * cse_p13_t63_g452) + (cse_p16_t1_g643 * cse_p17_t1_g660) + (cse_p19_t1_g684 * (r43 + (cse_p14_t1_g564 * dz23) + (cse_p16_t1_g643 * cse_p25_t1_g743))) + (cse_p23_t1_g726 * dz23) + (-((cse_p13_t1_g488 * cse_p13_t45_g434))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 2, h_z2_z3);
        double cse_p22_t1_g720 = (cse_p12_t1_g380 * dz23);
        double cse_p24_t1_g735 = (cse_p12_t1_g385 * dz23);
        double h_z2_x4 = ((cse_p11_t127_g135 * cse_p15_t1_g606) + (cse_p11_t1_g218 * ((cse_p22_t1_g720 * dx43) + (cse_p61_t10_invr777 * dy13) + (-((cse_p13_t1_g480 * cse_p13_t47_g436))) + (-((cse_p13_t1_g488 * cx))) + (-((cse_p24_t1_g735 * dx43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x4, 0, h_z2_x4);
        double h_z2_y4 = ((cse_p11_t128_g136 * cse_p15_t1_g606) + (cse_p11_t1_g218 * ((cse_p16_t19_g632 * cse_p61_t10_invr777) + (cse_p22_t1_g720 * dy43) + (-((cse_p13_t1_g484 * cse_p13_t47_g436))) + (-((cse_p13_t1_g488 * cy))) + (-((cse_p24_t1_g735 * dy43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x4, 1, h_z2_y4);
        double h_z2_z4 = ((cse_p11_t129_g137 * cse_p15_t1_g606) + (cse_p11_t1_g218 * ((cse_p22_t1_g720 * dz43) + (-((cse_p13_t1_g487 * cse_p13_t47_g436))) + (-((cse_p13_t1_g488 * cz))) + (-((cse_p24_t1_g735 * dz43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x4, 2, h_z2_z4);
        double cse_p16_t1_g644 = (-2.00000000000000000e+0 * cse_p13_t1_g490);
        double cse_p16_t1_g645 = (-2.00000000000000000e+0 * cse_p61_t12_invr3779);
        double cse_p25_t1_g745 = (cse_p15_t10_g576 * cse_p16_t1_g645);
        double h_x3_x3 = ((cse_p11_t116_g124 * cse_p17_t1_g663) + (cse_p11_t1_g218 * ((2.00000000000000000e+0 * cse_p13_t50_g439 * cse_p13_t61_g450) + (v * ((-((cse_p15_t10_g576 * cse_p25_t1_g745))) + (-((cse_p61_t11_invr2778 * (cse_p11_t31_g39 + cse_p11_t31_g39 + cse_p17_t1_g664 + cse_p17_t1_g665 + cse_p18_t1_g676 + (cse_p14_t1_g558 * cse_p16_t1_g639) + (cse_p14_t1_g559 * cse_p15_t1_g610) + (cse_p14_t1_g562 * (cse_p20_t1_g687 + cse_p20_t1_g687 + cse_p20_t4_g690 + cse_p20_t4_g690 + (cse_p15_t1_g608 * cse_p520_t6_invr835))) + (cse_p15_t1_g610 * cse_p26_t1_g759) + (cse_p16_t1_g639 * cse_p26_t1_g760) + (cse_p16_t1_g644 * cse_p17_t1_g662))))))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 0, h_x3_x3);
        double cse_p27_t1_g765 = (cse_p16_t1_g644 * dx43);
        double h_x3_y3 = ((cse_p11_t116_g124 * cse_p15_t1_g603) + (cse_p11_t1_g218 * ((cse_p13_t46_g435 * cse_p13_t61_g450) + (cse_p13_t50_g439 * cse_p13_t62_g451) + (v * ((-((cse_p15_t11_g577 * cse_p25_t1_g745))) + (-((cse_p61_t11_invr2778 * ((cse_p13_t6_g395 * ((cse_p15_t1_g608 * cse_p17_t1_g659) + (cse_p20_t1_g687 * dy43) + (cse_p20_t2_g688 * dx43) + (cse_p20_t4_g690 * dy43) + (cse_p20_t5_g691 * dx43))) + (cse_p14_t1_g558 * cse_p20_t1_g702) + (cse_p14_t1_g559 * cse_p20_t1_g697) + (cse_p15_t1_g609 * cse_p20_t1_g698) + (cse_p15_t1_g609 * cse_p20_t1_g701) + (cse_p20_t1_g697 * cse_p26_t1_g759) + (cse_p20_t1_g702 * cse_p26_t1_g760) + (cse_p27_t1_g765 * dy43))))))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 1, h_x3_y3);
        double h_x3_z3 = ((cse_p11_t116_g124 * cse_p14_t1_g565) + (cse_p11_t1_g218 * ((cse_p13_t45_g434 * cse_p13_t61_g450) + (cse_p13_t50_g439 * cse_p13_t63_g452) + (v * ((-((cse_p15_t12_g578 * cse_p25_t1_g745))) + (-((cse_p61_t11_invr2778 * ((cse_p13_t6_g395 * ((cse_p15_t1_g604 * cse_p15_t1_g608) + (cse_p20_t1_g687 * dz43) + (cse_p20_t3_g689 * dx43) + (cse_p20_t4_g690 * dz43) + (cse_p20_t6_g692 * dx43))) + (cse_p14_t1_g558 * cse_p21_t1_g704) + (cse_p14_t1_g559 * cse_p20_t1_g699) + (cse_p15_t1_g609 * cse_p20_t1_g700) + (cse_p15_t1_g609 * cse_p21_t1_g703) + (cse_p20_t1_g699 * cse_p26_t1_g759) + (cse_p21_t1_g704 * cse_p26_t1_g760) + (cse_p27_t1_g765 * dz43))))))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 2, h_x3_z3);
        double cse_p15_t1_g602 = (-2.00000000000000000e+0 * cse_p11_t1_g222);
        double cse_p25_t1_g744 = (cse_p15_t10_g576 * cse_p15_t1_g602);
        double h_x3_x4 = ((cse_p11_t127_g135 * cse_p17_t1_g663) + (cse_p11_t1_g218 * ((cse_p13_t61_g450 * cx) + (v * ((-((cse_p25_t1_g744 * dx43))) + (-((cse_p61_t11_invr2778 * (cse_p11_t173_g181 + (cse_p16_t1_g646 * cse_p17_t1_g662) + (-((cse_p14_t1_g562 * cse_p20_t1_g687))) + (-((cse_p14_t1_g562 * cse_p20_t4_g690))) + (-((cse_p17_t1_g662 * cse_p26_t1_g761))))))))) + (-((cse_p13_t1_g480 * cse_p13_t50_g439))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x4, 0, h_x3_x4);
        double h_x3_y4 = ((cse_p11_t128_g136 * cse_p17_t1_g663) + (cse_p11_t1_g218 * ((cse_p13_t61_g450 * cy) + (cse_p15_t18_g584 * cse_p61_t10_invr777) + (v * ((-((cse_p25_t1_g744 * dy43))) + (-((cse_p61_t11_invr2778 * (cse_p11_t13_g21 + cse_p11_t156_g164 + (-((cse_p14_t1_g563 * cse_p20_t1_g687))) + (-((cse_p14_t1_g563 * cse_p20_t4_g690))))))))) + (-((cse_p13_t1_g484 * cse_p13_t50_g439))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x4, 1, h_x3_y4);
        double h_x3_z4 = ((cse_p11_t129_g137 * cse_p17_t1_g663) + (cse_p11_t1_g218 * ((cse_p13_t61_g450 * cz) + (cse_p15_t15_g581 * cse_p61_t10_invr777) + (v * ((-((cse_p25_t1_g744 * dz43))) + (-((cse_p61_t11_invr2778 * (cse_p11_t14_g22 + cse_p11_t157_g165 + (-((cse_p14_t1_g564 * cse_p20_t1_g687))) + (-((cse_p14_t1_g564 * cse_p20_t4_g690))))))))) + (-((cse_p13_t1_g487 * cse_p13_t50_g439))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x4, 2, h_x3_z4);
        double cse_p26_t1_g762 = (cse_p15_t11_g577 * cse_p16_t1_g645);
        double h_y3_y3 = ((cse_p11_t112_g120 * cse_p15_t1_g603) + (cse_p11_t1_g218 * ((2.00000000000000000e+0 * cse_p13_t46_g435 * cse_p13_t62_g451) + (v * ((-((cse_p15_t11_g577 * cse_p26_t1_g762))) + (-((cse_p61_t11_invr2778 * (cse_p11_t57_g65 + cse_p11_t57_g65 + cse_p17_t1_g664 + cse_p17_t1_g665 + cse_p18_t1_g676 + (cse_p14_t1_g558 * cse_p16_t1_g641) + (cse_p14_t1_g559 * cse_p16_t1_g640) + (cse_p14_t1_g563 * (cse_p20_t2_g688 + cse_p20_t2_g688 + cse_p20_t5_g691 + cse_p20_t5_g691 + (cse_p11_t1_g209 * cse_p17_t1_g659))) + (cse_p16_t1_g640 * cse_p26_t1_g759) + (cse_p16_t1_g641 * cse_p26_t1_g760) + (cse_p16_t1_g644 * cse_p18_t1_g673))))))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 1, h_y3_y3);
        double cse_p26_t1_g755 = (cse_p11_t1_g209 * cse_p15_t1_g604);
        double h_y3_z3 = ((cse_p11_t112_g120 * cse_p14_t1_g565) + (cse_p11_t1_g218 * ((cse_p13_t45_g434 * cse_p13_t62_g451) + (cse_p13_t46_g435 * cse_p13_t63_g452) + (v * ((-((cse_p15_t12_g578 * cse_p26_t1_g762))) + (-((cse_p61_t11_invr2778 * ((cse_p13_t6_g395 * ((cse_p20_t2_g688 * dz43) + (cse_p20_t3_g689 * dy43) + (cse_p20_t5_g691 * dz43) + (cse_p20_t6_g692 * dy43) + (cse_p26_t1_g755 * dy43))) + (cse_p14_t1_g558 * cse_p21_t1_g708) + (cse_p14_t1_g559 * cse_p21_t1_g705) + (cse_p15_t1_g609 * cse_p21_t1_g706) + (cse_p15_t1_g609 * cse_p21_t1_g707) + (cse_p16_t1_g644 * cse_p25_t1_g751) + (cse_p21_t1_g705 * cse_p26_t1_g759) + (cse_p21_t1_g708 * cse_p26_t1_g760))))))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 2, h_y3_z3);
        double cse_p25_t1_g746 = (cse_p15_t11_g577 * cse_p15_t1_g602);
        double h_y3_x4 = ((cse_p11_t127_g135 * cse_p15_t1_g603) + (cse_p11_t1_g218 * ((cse_p13_t62_g451 * cx) + (cse_p15_t17_g583 * cse_p61_t10_invr777) + (v * ((-((cse_p25_t1_g746 * dx43))) + (-((cse_p61_t11_invr2778 * (cse_p11_t13_g21 + cse_p11_t156_g164 + (-((cse_p14_t1_g562 * cse_p20_t2_g688))) + (-((cse_p14_t1_g562 * cse_p20_t5_g691))))))))) + (-((cse_p13_t1_g480 * cse_p13_t46_g435))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x4, 0, h_y3_x4);
        double h_y3_y4 = ((cse_p11_t128_g136 * cse_p15_t1_g603) + (cse_p11_t1_g218 * ((cse_p13_t62_g451 * cy) + (v * ((-((cse_p25_t1_g746 * dy43))) + (-((cse_p61_t11_invr2778 * (cse_p11_t173_g181 + (cse_p16_t1_g646 * cse_p18_t1_g673) + (-((cse_p14_t1_g563 * cse_p20_t2_g688))) + (-((cse_p14_t1_g563 * cse_p20_t5_g691))) + (-((cse_p18_t1_g673 * cse_p26_t1_g761))))))))) + (-((cse_p13_t1_g484 * cse_p13_t46_g435))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x4, 1, h_y3_y4);
        double h_y3_z4 = ((cse_p11_t129_g137 * cse_p15_t1_g603) + (cse_p11_t1_g218 * ((cse_p13_t62_g451 * cz) + (cse_p15_t14_g580 * cse_p61_t10_invr777) + (v * ((-((cse_p25_t1_g746 * dz43))) + (-((cse_p61_t11_invr2778 * (cse_p11_t15_g23 + cse_p11_t164_g172 + (-((cse_p14_t1_g564 * cse_p20_t2_g688))) + (-((cse_p14_t1_g564 * cse_p20_t5_g691))))))))) + (-((cse_p13_t1_g487 * cse_p13_t46_g435))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x4, 2, h_y3_z4);
        double h_z3_z3 = ((cse_p11_t111_g119 * cse_p14_t1_g565) + (cse_p11_t1_g218 * ((2.00000000000000000e+0 * cse_p13_t45_g434 * cse_p13_t63_g452) + (v * ((-((cse_p15_t12_g578 * cse_p15_t12_g578 * cse_p16_t1_g645))) + (-((cse_p61_t11_invr2778 * (cse_p11_t80_g88 + cse_p11_t80_g88 + cse_p17_t1_g664 + cse_p17_t1_g665 + cse_p18_t1_g676 + (cse_p14_t1_g558 * cse_p16_t1_g643) + (cse_p14_t1_g559 * cse_p16_t1_g642) + (cse_p14_t1_g564 * (cse_p20_t3_g689 + cse_p20_t3_g689 + cse_p20_t6_g692 + cse_p20_t6_g692 + cse_p26_t1_g755)) + (cse_p16_t1_g642 * cse_p26_t1_g759) + (cse_p16_t1_g643 * cse_p26_t1_g760) + (cse_p16_t1_g644 * cse_p18_t1_g674))))))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x3, 2, h_z3_z3);
        double cse_p25_t1_g747 = (cse_p15_t12_g578 * cse_p15_t1_g602);
        double h_z3_x4 = ((cse_p11_t127_g135 * cse_p14_t1_g565) + (cse_p11_t1_g218 * ((cse_p13_t63_g452 * cx) + (cse_p15_t16_g582 * cse_p61_t10_invr777) + (v * ((-((cse_p25_t1_g747 * dx43))) + (-((cse_p61_t11_invr2778 * (cse_p11_t14_g22 + cse_p11_t157_g165 + (-((cse_p14_t1_g562 * cse_p20_t3_g689))) + (-((cse_p14_t1_g562 * cse_p20_t6_g692))))))))) + (-((cse_p13_t1_g480 * cse_p13_t45_g434))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x4, 0, h_z3_x4);
        double h_z3_y4 = ((cse_p11_t128_g136 * cse_p14_t1_g565) + (cse_p11_t1_g218 * ((cse_p13_t63_g452 * cy) + (cse_p15_t13_g579 * cse_p61_t10_invr777) + (v * ((-((cse_p25_t1_g747 * dy43))) + (-((cse_p61_t11_invr2778 * (cse_p11_t15_g23 + cse_p11_t164_g172 + (-((cse_p14_t1_g563 * cse_p20_t3_g689))) + (-((cse_p14_t1_g563 * cse_p20_t6_g692))))))))) + (-((cse_p13_t1_g484 * cse_p13_t45_g434))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x4, 1, h_z3_y4);
        double h_z3_z4 = ((cse_p11_t129_g137 * cse_p14_t1_g565) + (cse_p11_t1_g218 * ((cse_p13_t63_g452 * cz) + (v * ((-((cse_p25_t1_g747 * dz43))) + (-((cse_p61_t11_invr2778 * (cse_p11_t173_g181 + (cse_p16_t1_g646 * cse_p18_t1_g674) + (-((cse_p14_t1_g564 * cse_p20_t3_g689))) + (-((cse_p14_t1_g564 * cse_p20_t6_g692))) + (-((cse_p18_t1_g674 * cse_p26_t1_g761))))))))) + (-((cse_p13_t1_g487 * cse_p13_t45_g434))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x4, 2, h_z3_z4);
        double cse_p12_t1_g388 = (cse_p13_t1_g490 * cse_p22_t1_g717);
        double cse_p20_t1_g696 = (cse_p17_t1_g666 * cse_p520_t5_sqrt834 * r13 * r23);
        double cse_p22_t1_g713 = (cse_p11_t127_g135 * cse_p12_t1_g370);
        double cse_p24_t1_g736 = (cse_p13_t1_g480 * v);
        double h_x4_x4 = ((cse_p11_t127_g135 * cse_p22_t1_g713) + (cse_p11_t1_g218 * (cse_p11_t174_g182 + (-2.00000000000000000e+0 * cse_p13_t1_g480 * cx) + (cse_p12_t1_g388 * cse_p17_t1_g662) + (cse_p17_t1_g662 * cse_p20_t1_g696 * cse_p520_t6_invr835) + (-((cse_p24_t1_g736 * cse_p520_t6_invr835 * dx43))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 0, term.i3x4, 0, h_x4_x4);
        double cse_p25_t1_g752 = (cse_p11_t126_g134 * cse_p520_t5_sqrt834 * v);
        double cse_p27_t1_g764 = (cse_p15_t1_g608 * cse_p25_t1_g752);
        double h_x4_y4 = ((cse_p11_t128_g136 * cse_p22_t1_g713) + (cse_p11_t1_g218 * ((cse_p27_t1_g764 * dy43) + (-((cse_p13_t1_g480 * cy))) + (-((cse_p13_t1_g484 * cx))) + (-((cse_p17_t1_g659 * cse_p24_t1_g736))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 0, term.i3x4, 1, h_x4_y4);
        double h_x4_z4 = ((cse_p11_t129_g137 * cse_p22_t1_g713) + (cse_p11_t1_g218 * ((cse_p27_t1_g764 * dz43) + (-((cse_p13_t1_g480 * cz))) + (-((cse_p13_t1_g487 * cx))) + (-((cse_p15_t1_g604 * cse_p24_t1_g736))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 0, term.i3x4, 2, h_x4_z4);
        double cse_p26_t1_g754 = (cse_p11_t128_g136 * cse_p12_t1_g370);
        double cse_p26_t1_g756 = (cse_p13_t1_g484 * v);
        double h_y4_y4 = ((cse_p11_t128_g136 * cse_p26_t1_g754) + (cse_p11_t1_g218 * (cse_p11_t174_g182 + (-2.00000000000000000e+0 * cse_p13_t1_g484 * cy) + (cse_p12_t1_g388 * cse_p18_t1_g673) + (cse_p17_t1_g659 * cse_p20_t1_g696 * dy43) + (-((cse_p17_t1_g659 * cse_p26_t1_g756))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 1, term.i3x4, 1, h_y4_y4);
        double h_y4_z4 = ((cse_p11_t129_g137 * cse_p26_t1_g754) + (cse_p11_t1_g218 * ((cse_p11_t1_g209 * cse_p25_t1_g751 * cse_p25_t1_g752) + (-((cse_p13_t1_g484 * cz))) + (-((cse_p13_t1_g487 * cy))) + (-((cse_p15_t1_g604 * cse_p26_t1_g756))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 1, term.i3x4, 2, h_y4_z4);
        double h_z4_z4 = ((cse_p11_t129_g137 * cse_p11_t129_g137 * cse_p12_t1_g370) + (cse_p11_t1_g218 * (cse_p11_t174_g182 + (-2.00000000000000000e+0 * cse_p13_t1_g487 * cz) + (cse_p12_t1_g388 * cse_p18_t1_g674) + (cse_p15_t1_g604 * cse_p20_t1_g696 * dz43) + (-((cse_p13_t1_g487 * cse_p15_t1_g604 * v))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 2, term.i3x4, 2, h_z4_z4);
      }
    }
  }
  return energy_added;
}
void hessian_fd(const chiral_term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 12;
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
    double saved = position[term.i3x4 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x4 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x4 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x4 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x4, 0, d);
  }
  {
    double saved = position[term.i3x4 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x4 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x4 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x4 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x4, 1, d);
  }
  {
    double saved = position[term.i3x4 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x4 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x4 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x4 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x4, 2, d);
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
    double saved = position[term.i3x4 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x4 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x4 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x4 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x4, 0, hval);
  }
  {
    double saved = position[term.i3x4 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x4 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x4 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x4 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x4, 1, hval);
  }
  {
    double saved = position[term.i3x4 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x4 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x4 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x4 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x4, 2, hval);
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
  {
    double saved_i = position[term.i3x4 + 0];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 0] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 0] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 0] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 0];
    double saved_j = position[term.i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 0] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 0] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 0] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 0];
    double saved_j = position[term.i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 0] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 0] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 0] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x4 + 0];
    double saved_j = position[term.i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 0] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 0] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 0] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 0];
    double saved_j = position[term.i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 0] = saved_i + h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 0] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 0] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x2, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 0];
    double saved_j = position[term.i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 0] = saved_i + h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 0] = saved_i - h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 0] = saved_i; position[term.i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x2, 2, hval);
  }
  {
    double saved_i = position[term.i3x4 + 0];
    double saved_j = position[term.i3x3 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 0] = saved_i + h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 0] = saved_i - h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 0] = saved_i; position[term.i3x3 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x3, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 0];
    double saved_j = position[term.i3x3 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 0] = saved_i + h; position[term.i3x3 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 0] = saved_i - h; position[term.i3x3 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 0] = saved_i; position[term.i3x3 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x3, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 0];
    double saved_j = position[term.i3x3 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 0] = saved_i + h; position[term.i3x3 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 0] = saved_i - h; position[term.i3x3 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 0] = saved_i; position[term.i3x3 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x3, 2, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x2, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x2, 2, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x3 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x3 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x3, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x3 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x3 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x3 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x3 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x3, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x3 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x3 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x3 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x3 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x3, 2, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x4 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x4 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x4 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x4 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x4 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x4 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x4, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x2, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x2, 2, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x3 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x3 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x3, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x3 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x3 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x3 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x3 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x3, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x3 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x3 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x3 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x3 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x3, 2, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x4 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x4 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x4 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x4 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x4 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x4 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x4, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x4 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x4 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x4 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x4 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x4 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x4 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x4, 1, hval);
  }
}

};
