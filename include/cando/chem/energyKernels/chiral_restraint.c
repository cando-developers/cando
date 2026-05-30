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
      double cse_p11_t1_g559493 = (-(x3));
      double cse_p11_t2_g559494 = (-(y3));
      double cse_p11_t3_g559495 = (-(z3));
      double dx13 = (cse_p11_t1_g559493 + x1);
      double dy13 = (cse_p11_t2_g559494 + y1);
      double dz13 = (cse_p11_t3_g559495 + z1);
      double dx23 = (cse_p11_t1_g559493 + x2);
      double dy23 = (cse_p11_t2_g559494 + y2);
      double dz23 = (cse_p11_t3_g559495 + z2);
      double dx43 = (cse_p11_t1_g559493 + x4);
      double dy43 = (cse_p11_t2_g559494 + y4);
      double dz43 = (cse_p11_t3_g559495 + z4);
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
      double cse_p11_t8_g559511 = (-(x3));
      double cse_p11_t9_g559512 = (-(y3));
      double cse_p11_t10_g559513 = (-(z3));
      double dx13 = (cse_p11_t8_g559511 + x1);
      double dy13 = (cse_p11_t9_g559512 + y1);
      double dz13 = (cse_p11_t10_g559513 + z1);
      double dx23 = (cse_p11_t8_g559511 + x2);
      double dy23 = (cse_p11_t9_g559512 + y2);
      double dz23 = (cse_p11_t10_g559513 + z2);
      double dx43 = (cse_p11_t8_g559511 + x4);
      double dy43 = (cse_p11_t9_g559512 + y4);
      double dz43 = (cse_p11_t10_g559513 + z4);
      double r13_2 = ((dx13 * dx13) + (dy13 * dy13) + (dz13 * dz13));
      double r23_2 = ((dx23 * dx23) + (dy23 * dy23) + (dz23 * dz23));
      double r43_2 = ((dx43 * dx43) + (dy43 * dy43) + (dz43 * dz43));
      double cse_p1_t1_sqrt559496 = sqrt(r13_2);
      double cse_p1_t2_invsqrt559497 = (1.0 / (cse_p1_t1_sqrt559496));
      double r13 = cse_p1_t1_sqrt559496;
      double cse_p1_t3_sqrt559498 = sqrt(r23_2);
      double cse_p1_t4_invsqrt559499 = (1.0 / (cse_p1_t3_sqrt559498));
      double r23 = cse_p1_t3_sqrt559498;
      double cse_p1_t5_sqrt559500 = sqrt(r43_2);
      double cse_p1_t6_invsqrt559501 = (1.0 / (cse_p1_t5_sqrt559500));
      double r43 = cse_p1_t5_sqrt559500;
      double cx = ((dy13 * dz23) + (-((dy23 * dz13))));
      double cy = ((dx23 * dz13) + (-((dx13 * dz23))));
      double cz = ((dx13 * dy23) + (-((dx23 * dy13))));
      double v = ((cx * dx43) + (cy * dy43) + (cz * dz43));
      double cse_p11_t1_g559523 = (r13 * r23);
      double denom = (cse_p11_t1_g559523 * r43);
      double cse_p1_t7_invr559502 = (1.0 / (denom));
      double cse_p1_t8_invr2559503 = (cse_p1_t7_invr559502 * cse_p1_t7_invr559502);
      double q = (cse_p1_t7_invr559502 * v);
      double cse_p11_t1_g559504 = (term.co + q);
      double chiral_e = (term.k * (cse_p11_t1_g559504 * cse_p11_t1_g559504 * cse_p11_t1_g559504));
      if ((chiral_e > 0.00000000000000000e+0)) {
        /* !ACTIVE */
        double cse_p11_t2_g559505 = (-(dx13));
        double cse_p11_t3_g559506 = (-(dx23));
        double cse_p11_t4_g559507 = (-(dy13));
        double cse_p11_t5_g559508 = (-(dy23));
        double cse_p11_t6_g559509 = (-(dz13));
        double cse_p11_t7_g559510 = (-(dz23));
        double cse_p11_t11_g559514 = ((cse_p11_t1_g559504) * (cse_p11_t1_g559504));
        double energy = chiral_e;
        energy_added += energy;
        *energy_accumulate += energy;
        double cse_p11_t1_g559515 = (3.00000000000000000e+0 * cse_p11_t11_g559514 * term.k);
        double cse_p11_t1_g559516 = (cse_p1_t2_invsqrt559497 * r23 * r43);
        double cse_p62_t1_g559524 = (cse_p1_t8_invr2559503 * v);
        double cse_p11_t1_g559520 = (cse_p11_t1_g559516 * cse_p62_t1_g559524);
        double g_x1 = (cse_p11_t1_g559515 * ((cse_p1_t7_invr559502 * ((cse_p11_t7_g559510 * dy43) + (dy23 * dz43))) + (-((cse_p11_t1_g559520 * dx13)))));
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = (cse_p11_t1_g559515 * ((cse_p1_t7_invr559502 * ((cse_p11_t3_g559506 * dz43) + (dx43 * dz23))) + (-((cse_p11_t1_g559520 * dy13)))));
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = (cse_p11_t1_g559515 * ((cse_p1_t7_invr559502 * ((cse_p11_t5_g559508 * dx43) + (dx23 * dy43))) + (-((cse_p11_t1_g559520 * dz13)))));
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double cse_p11_t1_g559518 = (cse_p1_t4_invsqrt559499 * r13 * r43);
        double cse_p11_t1_g559521 = (cse_p11_t1_g559518 * cse_p62_t1_g559524);
        double g_x2 = (cse_p11_t1_g559515 * ((cse_p1_t7_invr559502 * ((cse_p11_t4_g559507 * dz43) + (dy43 * dz13))) + (-((cse_p11_t1_g559521 * dx23)))));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (cse_p11_t1_g559515 * ((cse_p1_t7_invr559502 * ((cse_p11_t6_g559509 * dx43) + (dx13 * dz43))) + (-((cse_p11_t1_g559521 * dy23)))));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (cse_p11_t1_g559515 * ((cse_p1_t7_invr559502 * ((cse_p11_t2_g559505 * dy43) + (dx43 * dy13))) + (-((cse_p11_t1_g559521 * dz23)))));
        KernelGradientAcc(term.i3x2, 2, g_z2);
        double cse_p11_t1_g559519 = (cse_p11_t1_g559523 * cse_p1_t6_invsqrt559501);
        double g_x3 = (cse_p11_t1_g559515 * ((cse_p1_t7_invr559502 * ((dy43 * (cse_p11_t6_g559509 + dz23)) + (dz43 * (cse_p11_t5_g559508 + dy13)) + (-(cx)))) + (-((cse_p62_t1_g559524 * ((-((cse_p11_t1_g559516 * dx13))) + (-((cse_p11_t1_g559518 * dx23))) + (-((cse_p11_t1_g559519 * dx43)))))))));
        KernelGradientAcc(term.i3x3, 0, g_x3);
        double g_y3 = (cse_p11_t1_g559515 * ((cse_p1_t7_invr559502 * ((dx43 * (cse_p11_t7_g559510 + dz13)) + (dz43 * (cse_p11_t2_g559505 + dx23)) + (-(cy)))) + (-((cse_p62_t1_g559524 * ((-((cse_p11_t1_g559516 * dy13))) + (-((cse_p11_t1_g559518 * dy23))) + (-((cse_p11_t1_g559519 * dy43)))))))));
        KernelGradientAcc(term.i3x3, 1, g_y3);
        double g_z3 = (cse_p11_t1_g559515 * ((cse_p1_t7_invr559502 * ((dx43 * (cse_p11_t4_g559507 + dy23)) + (dy43 * (cse_p11_t3_g559506 + dx13)) + (-(cz)))) + (-((cse_p62_t1_g559524 * ((-((cse_p11_t1_g559516 * dz13))) + (-((cse_p11_t1_g559518 * dz23))) + (-((cse_p11_t1_g559519 * dz43)))))))));
        KernelGradientAcc(term.i3x3, 2, g_z3);
        double cse_p11_t1_g559522 = (cse_p11_t1_g559519 * cse_p62_t1_g559524);
        double g_x4 = (cse_p11_t1_g559515 * ((cse_p1_t7_invr559502 * cx) + (-((cse_p11_t1_g559522 * dx43)))));
        KernelGradientAcc(term.i3x4, 0, g_x4);
        double g_y4 = (cse_p11_t1_g559515 * ((cse_p1_t7_invr559502 * cy) + (-((cse_p11_t1_g559522 * dy43)))));
        KernelGradientAcc(term.i3x4, 1, g_y4);
        double g_z4 = (cse_p11_t1_g559515 * ((cse_p1_t7_invr559502 * cz) + (-((cse_p11_t1_g559522 * dz43)))));
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
      double cse_p11_t189_g559721 = (-(x3));
      double cse_p11_t190_g559722 = (-(y3));
      double cse_p11_t191_g559723 = (-(z3));
      double dx13 = (cse_p11_t189_g559721 + x1);
      double dy13 = (cse_p11_t190_g559722 + y1);
      double dz13 = (cse_p11_t191_g559723 + z1);
      double dx23 = (cse_p11_t189_g559721 + x2);
      double dy23 = (cse_p11_t190_g559722 + y2);
      double dz23 = (cse_p11_t191_g559723 + z2);
      double dx43 = (cse_p11_t189_g559721 + x4);
      double dy43 = (cse_p11_t190_g559722 + y4);
      double dz43 = (cse_p11_t191_g559723 + z4);
      double cse_p15_t1_g560134 = (dx13 * dx13);
      double cse_p16_t1_g560164 = (dy13 * dy13);
      double cse_p16_t1_g560166 = (dz13 * dz13);
      double r13_2 = (cse_p15_t1_g560134 + cse_p16_t1_g560164 + cse_p16_t1_g560166);
      double cse_p16_t1_g560163 = (dx23 * dx23);
      double cse_p16_t1_g560165 = (dy23 * dy23);
      double cse_p16_t1_g560167 = (dz23 * dz23);
      double r23_2 = (cse_p16_t1_g560163 + cse_p16_t1_g560165 + cse_p16_t1_g560167);
      double cse_p17_t1_g560186 = (dx43 * dx43);
      double cse_p18_t1_g560197 = (dy43 * dy43);
      double cse_p18_t1_g560198 = (dz43 * dz43);
      double r43_2 = (cse_p17_t1_g560186 + cse_p18_t1_g560197 + cse_p18_t1_g560198);
      double cse_p520_t1_sqrt560354 = sqrt(r13_2);
      double cse_p520_t2_invr560355 = (1.0 / (r13_2));
      double cse_p61_t3_invr2560294 = (cse_p520_t2_invr560355 * cse_p520_t2_invr560355);
      double cse_p1_t2_invsqrt559526 = (1.0 / (cse_p520_t1_sqrt560354));
      double r13 = cse_p520_t1_sqrt560354;
      double cse_p520_t3_sqrt560356 = sqrt(r23_2);
      double cse_p520_t4_invr560357 = (1.0 / (r23_2));
      double cse_p61_t6_invr2560297 = (cse_p520_t4_invr560357 * cse_p520_t4_invr560357);
      double cse_p1_t4_invsqrt559528 = (1.0 / (cse_p520_t3_sqrt560356));
      double r23 = cse_p520_t3_sqrt560356;
      double cse_p520_t5_sqrt560358 = sqrt(r43_2);
      double cse_p520_t6_invr560359 = (1.0 / (r43_2));
      double cse_p61_t9_invr2560300 = (cse_p520_t6_invr560359 * cse_p520_t6_invr560359);
      double cse_p1_t6_invsqrt559530 = (1.0 / (cse_p520_t5_sqrt560358));
      double r43 = cse_p520_t5_sqrt560358;
      double cse_p21_t1_g560230 = (dy13 * dz23);
      double cse_p21_t1_g560231 = (dy23 * dz13);
      double cx = (cse_p21_t1_g560230 + (-(cse_p21_t1_g560231)));
      double cse_p20_t1_g560224 = (dx13 * dz23);
      double cse_p21_t1_g560227 = (dx23 * dz13);
      double cy = (cse_p21_t1_g560227 + (-(cse_p20_t1_g560224)));
      double cse_p20_t1_g560222 = (dx13 * dy23);
      double cse_p20_t1_g560225 = (dx23 * dy13);
      double cz = (cse_p20_t1_g560222 + (-(cse_p20_t1_g560225)));
      double v = ((cx * dx43) + (cy * dy43) + (cz * dz43));
      double cse_p11_t1_g559733 = (r13 * r23);
      double denom = (cse_p11_t1_g559733 * r43);
      double cse_p61_t10_invr560301 = (1.0 / (denom));
      double cse_p61_t11_invr2560302 = (cse_p61_t10_invr560301 * cse_p61_t10_invr560301);
      double cse_p61_t12_invr3560303 = (cse_p61_t10_invr560301 * cse_p61_t11_invr2560302);
      double q = (cse_p61_t10_invr560301 * v);
      double cse_p11_t142_g559674 = (term.co + q);
      double chiral_e = (term.k * (cse_p11_t142_g559674 * cse_p11_t142_g559674 * cse_p11_t142_g559674));
      if ((chiral_e > 0.00000000000000000e+0)) {
        /* !ACTIVE */
        double cse_p63_t1_g560305 = (cse_p61_t11_invr2560302 * v);
        double cse_p16_t19_g560156 = (-(dx13));
        double cse_p16_t20_g560157 = (-(dx23));
        double cse_p16_t21_g560158 = (-(dy13));
        double cse_p16_t22_g560159 = (-(dy23));
        double cse_p16_t23_g560160 = (-(dz13));
        double cse_p16_t24_g560161 = (-(dz23));
        double cse_p15_t13_g560103 = (cse_p16_t20_g560157 + dx13);
        double cse_p15_t14_g560104 = (cse_p16_t19_g560156 + dx23);
        double cse_p15_t15_g560105 = (cse_p16_t22_g560159 + dy13);
        double cse_p15_t16_g560106 = (cse_p16_t21_g560158 + dy23);
        double cse_p15_t17_g560107 = (cse_p16_t24_g560161 + dz13);
        double cse_p15_t18_g560108 = (cse_p16_t23_g560160 + dz23);
        double cse_p14_t14_g560035 = (dx13 * dz43);
        double cse_p14_t15_g560036 = (dx23 * dy43);
        double cse_p14_t16_g560037 = (cse_p15_t16_g560106 * dx43);
        double cse_p14_t17_g560038 = (cse_p15_t17_g560107 * dx43);
        double cse_p14_t18_g560039 = (cse_p16_t22_g560159 * dx43);
        double cse_p14_t19_g560040 = (cse_p16_t23_g560160 * dx43);
        double cse_p14_t20_g560041 = (dx43 * dy13);
        double cse_p14_t21_g560042 = (dx43 * dz23);
        double cse_p14_t22_g560043 = (dy23 * dz43);
        double cse_p14_t23_g560044 = (cse_p15_t13_g560103 * dy43);
        double cse_p14_t24_g560045 = (cse_p15_t18_g560108 * dy43);
        double cse_p14_t25_g560046 = (cse_p16_t19_g560156 * dy43);
        double cse_p14_t26_g560047 = (cse_p16_t24_g560161 * dy43);
        double cse_p14_t27_g560048 = (dy43 * dz13);
        double cse_p14_t28_g560049 = (cse_p15_t14_g560104 * dz43);
        double cse_p14_t29_g560050 = (cse_p15_t15_g560105 * dz43);
        double cse_p14_t30_g560051 = (cse_p16_t20_g560157 * dz43);
        double cse_p14_t31_g560052 = (cse_p16_t21_g560158 * dz43);
        double cse_p14_t50_g560071 = (-(cx));
        double cse_p14_t51_g560072 = (-(cy));
        double cse_p14_t52_g560073 = (-(cz));
        double cse_p13_t43_g559956 = (cse_p14_t14_g560035 + cse_p14_t19_g560040);
        double cse_p13_t44_g559957 = (cse_p14_t15_g560036 + cse_p14_t18_g560039);
        double cse_p13_t45_g559958 = (cse_p14_t16_g560037 + cse_p14_t23_g560044 + cse_p14_t52_g560073);
        double cse_p13_t46_g559959 = (cse_p14_t17_g560038 + cse_p14_t28_g560049 + cse_p14_t51_g560072);
        double cse_p13_t47_g559960 = (cse_p14_t20_g560041 + cse_p14_t25_g560046);
        double cse_p13_t48_g559961 = (cse_p14_t21_g560042 + cse_p14_t30_g560051);
        double cse_p13_t49_g559962 = (cse_p14_t22_g560043 + cse_p14_t26_g560047);
        double cse_p13_t50_g559963 = (cse_p14_t24_g560045 + cse_p14_t29_g560050 + cse_p14_t50_g560071);
        double cse_p13_t51_g559964 = (cse_p14_t27_g560048 + cse_p14_t31_g560052);
        double cse_p21_t1_g560236 = (2.00000000000000000e+0 * cse_p61_t11_invr2560302);
        double cse_p12_t4_g559756 = (cse_p21_t1_g560236 * cse_p61_t3_invr2560294);
        double cse_p12_t5_g559757 = (cse_p21_t1_g560236 * cse_p61_t6_invr2560297);
        double cse_p12_t6_g559758 = (cse_p21_t1_g560236 * cse_p61_t9_invr2560300);
        double cse_p12_t53_g559805 = (cse_p13_t43_g559956 * cse_p61_t10_invr560301);
        double cse_p12_t54_g559806 = (cse_p13_t44_g559957 * cse_p61_t10_invr560301);
        double cse_p12_t55_g559807 = (cse_p13_t45_g559958 * cse_p61_t10_invr560301);
        double cse_p12_t56_g559808 = (cse_p13_t46_g559959 * cse_p61_t10_invr560301);
        double cse_p12_t57_g559809 = (cse_p13_t47_g559960 * cse_p61_t10_invr560301);
        double cse_p12_t58_g559810 = (cse_p13_t48_g559961 * cse_p61_t10_invr560301);
        double cse_p12_t59_g559811 = (cse_p13_t49_g559962 * cse_p61_t10_invr560301);
        double cse_p12_t60_g559812 = (cse_p13_t50_g559963 * cse_p61_t10_invr560301);
        double cse_p12_t61_g559813 = (cse_p13_t51_g559964 * cse_p61_t10_invr560301);
        double cse_p12_t62_g559814 = (cse_p61_t10_invr560301 * cx);
        double cse_p12_t63_g559815 = (cse_p61_t10_invr560301 * cy);
        double cse_p12_t64_g559816 = (cse_p61_t10_invr560301 * cz);
        double cse_p13_t4_g559917 = (cse_p520_t1_sqrt560354 * cse_p520_t2_invr560355);
        double cse_p13_t5_g559918 = (cse_p520_t3_sqrt560356 * cse_p520_t4_invr560357);
        double cse_p13_t6_g559919 = (cse_p520_t5_sqrt560358 * cse_p520_t6_invr560359);
        double cse_p11_t1_g559734 = (r13 * r43);
        double cse_p22_t1_g560239 = (cse_p11_t1_g559734 * cse_p520_t4_invr560357);
        double cse_p12_t28_g559780 = (cse_p22_t1_g560239 * dx23);
        double cse_p12_t29_g559781 = (cse_p22_t1_g560239 * dy23);
        double cse_p12_t30_g559782 = (cse_p22_t1_g560239 * dz23);
        double cse_p11_t1_g559735 = (r23 * r43);
        double cse_p22_t1_g560240 = (cse_p11_t1_g559735 * cse_p520_t2_invr560355);
        double cse_p12_t38_g559790 = (cse_p22_t1_g560240 * dx13);
        double cse_p12_t39_g559791 = (cse_p22_t1_g560240 * dy13);
        double cse_p12_t40_g559792 = (cse_p22_t1_g560240 * dz13);
        double cse_p18_t1_g560200 = (cse_p11_t1_g559733 * cse_p13_t6_g559919);
        double cse_p25_t1_g560272 = (cse_p18_t1_g560200 * dx43);
        double cse_p27_t1_g560291 = (cse_p18_t1_g560200 * dy43);
        double cse_p17_t3_g560174 = (cse_p18_t1_g560200 * dz43);
        double cse_p16_t10_g560147 = (-(cse_p25_t1_g560272));
        double cse_p16_t11_g560148 = (-(cse_p27_t1_g560291));
        double cse_p16_t12_g560149 = (-(cse_p17_t3_g560174));
        double cse_p21_t1_g560235 = (2.00000000000000000e+0 * cse_p61_t12_invr3560303);
        double cse_p26_t1_g560285 = (cse_p18_t1_g560200 * cse_p520_t6_invr560359);
        double cse_p12_t7_g559759 = (cse_p21_t1_g560235 * cse_p26_t1_g560285);
        double cse_p17_t1_g560183 = (cse_p520_t6_invr560359 * dy43);
        double cse_p12_t10_g559762 = (cse_p17_t1_g560183 * cse_p25_t1_g560272);
        double cse_p15_t1_g560128 = (cse_p520_t6_invr560359 * dz43);
        double cse_p12_t11_g559763 = (cse_p15_t1_g560128 * cse_p25_t1_g560272);
        double cse_p20_t1_g560219 = (cse_p18_t1_g560200 * cse_p63_t1_g560305);
        double cse_p13_t7_g559920 = (cse_p20_t1_g560219 * dx43);
        double cse_p13_t8_g559921 = (cse_p20_t1_g560219 * dy43);
        double cse_p13_t9_g559922 = (cse_p20_t1_g560219 * dz43);
        double cse_p12_t14_g559766 = (cse_p15_t1_g560128 * cse_p27_t1_g560291);
        double cse_p14_t1_g560086 = (cse_p13_t6_g559919 * dx43);
        double cse_p12_t41_g559793 = (cse_p14_t1_g560086 * r13);
        double cse_p12_t42_g559794 = (cse_p14_t1_g560086 * r23);
        double cse_p14_t1_g560087 = (cse_p13_t6_g559919 * dy43);
        double cse_p12_t43_g559795 = (cse_p14_t1_g560087 * r13);
        double cse_p12_t44_g559796 = (cse_p14_t1_g560087 * r23);
        double cse_p14_t1_g560088 = (cse_p13_t6_g559919 * dz43);
        double cse_p12_t45_g559797 = (cse_p14_t1_g560088 * r13);
        double cse_p12_t46_g559798 = (cse_p14_t1_g560088 * r23);
        double cse_p12_t107_g559859 = (-(cse_p13_t7_g559920));
        double cse_p12_t108_g559860 = (-(cse_p13_t8_g559921));
        double cse_p12_t109_g559861 = (-(cse_p13_t9_g559922));
        double cse_p17_t1_g560189 = (cse_p11_t1_g559735 * cse_p13_t4_g559917);
        double cse_p17_t7_g560178 = (cse_p17_t1_g560189 * dx13);
        double cse_p17_t8_g560179 = (cse_p17_t1_g560189 * dy13);
        double cse_p17_t9_g560180 = (cse_p17_t1_g560189 * dz13);
        double cse_p16_t16_g560153 = (-(cse_p17_t7_g560178));
        double cse_p16_t17_g560154 = (-(cse_p17_t8_g560179));
        double cse_p16_t18_g560155 = (-(cse_p17_t9_g560180));
        double cse_p20_t1_g560218 = (cse_p17_t1_g560189 * cse_p63_t1_g560305);
        double cse_p13_t13_g559926 = (cse_p20_t1_g560218 * dx13);
        double cse_p13_t14_g559927 = (cse_p20_t1_g560218 * dy13);
        double cse_p13_t15_g559928 = (cse_p20_t1_g560218 * dz13);
        double cse_p26_t1_g560283 = (cse_p17_t1_g560189 * cse_p520_t2_invr560355);
        double cse_p12_t9_g559761 = (cse_p21_t1_g560235 * cse_p26_t1_g560283);
        double cse_p14_t1_g560090 = (cse_p13_t4_g559917 * r43);
        double cse_p12_t47_g559799 = (cse_p14_t1_g560090 * dx13);
        double cse_p12_t48_g559800 = (cse_p14_t1_g560090 * dy13);
        double cse_p12_t49_g559801 = (cse_p14_t1_g560090 * dz13);
        double cse_p12_t113_g559865 = (-(cse_p13_t13_g559926));
        double cse_p12_t114_g559866 = (-(cse_p13_t14_g559927));
        double cse_p12_t115_g559867 = (-(cse_p13_t15_g559928));
        double cse_p17_t1_g560188 = (cse_p11_t1_g559734 * cse_p13_t5_g559918);
        double cse_p17_t4_g560175 = (cse_p17_t1_g560188 * dx23);
        double cse_p17_t5_g560176 = (cse_p17_t1_g560188 * dy23);
        double cse_p17_t6_g560177 = (cse_p17_t1_g560188 * dz23);
        double cse_p16_t13_g560150 = (-(cse_p17_t4_g560175));
        double cse_p16_t14_g560151 = (-(cse_p17_t5_g560176));
        double cse_p16_t15_g560152 = (-(cse_p17_t6_g560177));
        double cse_p15_t10_g560100 = (cse_p16_t10_g560147 + cse_p16_t13_g560150 + cse_p16_t16_g560153);
        double cse_p15_t11_g560101 = (cse_p16_t11_g560148 + cse_p16_t14_g560151 + cse_p16_t17_g560154);
        double cse_p15_t12_g560102 = (cse_p16_t12_g560149 + cse_p16_t15_g560152 + cse_p16_t18_g560155);
        double cse_p14_t1_g560022 = (cse_p15_t10_g560100 * cse_p61_t11_invr2560302);
        double cse_p14_t2_g560023 = (cse_p15_t11_g560101 * cse_p61_t11_invr2560302);
        double cse_p14_t3_g560024 = (cse_p15_t12_g560102 * cse_p61_t11_invr2560302);
        double cse_p20_t1_g560217 = (cse_p17_t1_g560188 * cse_p63_t1_g560305);
        double cse_p13_t10_g559923 = (cse_p20_t1_g560217 * dx23);
        double cse_p13_t11_g559924 = (cse_p20_t1_g560217 * dy23);
        double cse_p13_t12_g559925 = (cse_p20_t1_g560217 * dz23);
        double cse_p13_t61_g559974 = (-(cse_p14_t1_g560022));
        double cse_p13_t62_g559975 = (-(cse_p14_t2_g560023));
        double cse_p13_t63_g559976 = (-(cse_p14_t3_g560024));
        double cse_p26_t1_g560284 = (cse_p17_t1_g560188 * cse_p520_t4_invr560357);
        double cse_p12_t8_g559760 = (cse_p21_t1_g560235 * cse_p26_t1_g560284);
        double cse_p24_t1_g560265 = (cse_p13_t5_g559918 * r43);
        double cse_p12_t50_g559802 = (cse_p24_t1_g560265 * dx23);
        double cse_p12_t51_g559803 = (cse_p24_t1_g560265 * dy23);
        double cse_p12_t52_g559804 = (cse_p24_t1_g560265 * dz23);
        double cse_p12_t83_g559835 = (cse_p13_t61_g559974 * v);
        double cse_p12_t84_g559836 = (cse_p13_t62_g559975 * v);
        double cse_p12_t85_g559837 = (cse_p13_t63_g559976 * v);
        double cse_p12_t110_g559862 = (-(cse_p13_t10_g559923));
        double cse_p12_t111_g559863 = (-(cse_p13_t11_g559924));
        double cse_p12_t112_g559864 = (-(cse_p13_t12_g559925));
        double cse_p19_t1_g560203 = (cse_p11_t1_g559735 * cse_p520_t1_sqrt560354);
        double cse_p13_t1_g560016 = (cse_p19_t1_g560203 * cse_p61_t3_invr2560294);
        double cse_p14_t1_g560083 = (-2.00000000000000000e+0 * cse_p13_t1_g560016);
        double cse_p17_t1_g560185 = (cse_p14_t1_g560083 * cse_p63_t1_g560305);
        double cse_p20_t1_g560221 = (dx13 * dy13);
        double cse_p11_t4_g559536 = (cse_p17_t1_g560185 * cse_p20_t1_g560221);
        double cse_p20_t1_g560223 = (dx13 * dz13);
        double cse_p11_t5_g559537 = (cse_p17_t1_g560185 * cse_p20_t1_g560223);
        double cse_p19_t1_g560201 = (cse_p11_t1_g559734 * cse_p520_t3_sqrt560356);
        double cse_p13_t1_g560015 = (cse_p19_t1_g560201 * cse_p61_t6_invr2560297);
        double cse_p14_t1_g560082 = (-2.00000000000000000e+0 * cse_p13_t1_g560015);
        double cse_p17_t1_g560184 = (cse_p14_t1_g560082 * cse_p63_t1_g560305);
        double cse_p20_t1_g560226 = (dx23 * dy23);
        double cse_p11_t6_g559538 = (cse_p17_t1_g560184 * cse_p20_t1_g560226);
        double cse_p21_t1_g560228 = (dx23 * dz23);
        double cse_p11_t7_g559539 = (cse_p17_t1_g560184 * cse_p21_t1_g560228);
        double cse_p21_t1_g560229 = (dy13 * dz13);
        double cse_p11_t8_g559540 = (cse_p17_t1_g560185 * cse_p21_t1_g560229);
        double cse_p21_t1_g560232 = (dy23 * dz23);
        double cse_p11_t9_g559541 = (cse_p17_t1_g560184 * cse_p21_t1_g560232);
        double cse_p13_t1_g560014 = (cse_p11_t1_g559733 * cse_p520_t5_sqrt560358 * cse_p61_t9_invr2560300);
        double cse_p16_t1_g560170 = (2.00000000000000000e+0 * cse_p13_t1_g560014);
        double cse_p27_t1_g560290 = (cse_p16_t1_g560170 * dx43);
        double cse_p11_t13_g559545 = (cse_p27_t1_g560290 * dy43);
        double cse_p11_t14_g559546 = (cse_p27_t1_g560290 * dz43);
        double cse_p25_t1_g560275 = (dy43 * dz43);
        double cse_p11_t15_g559547 = (cse_p16_t1_g560170 * cse_p25_t1_g560275);
        double cse_p15_t1_g560133 = (cse_p13_t5_g559918 * cse_p14_t1_g560090);
        double cse_p25_t1_g560273 = (dx13 * dx23);
        double cse_p11_t31_g559563 = (cse_p15_t1_g560133 * cse_p25_t1_g560273);
        double cse_p16_t1_g560162 = (cse_p15_t1_g560133 * cse_p63_t1_g560305);
        double cse_p11_t32_g559564 = (cse_p16_t1_g560162 * cse_p25_t1_g560273);
        double cse_p19_t1_g560209 = (cse_p13_t4_g559917 * r23);
        double cse_p20_t1_g560211 = (cse_p19_t1_g560209 * dx13);
        double cse_p20_t2_g560212 = (cse_p19_t1_g560209 * dy13);
        double cse_p20_t3_g560213 = (cse_p19_t1_g560209 * dz13);
        double cse_p19_t1_g560210 = (cse_p13_t5_g559918 * r13);
        double cse_p20_t4_g560214 = (cse_p19_t1_g560210 * dx23);
        double cse_p20_t5_g560215 = (cse_p19_t1_g560210 * dy23);
        double cse_p20_t6_g560216 = (cse_p19_t1_g560210 * dz23);
        double cse_p25_t1_g560274 = (dy13 * dy23);
        double cse_p11_t57_g559589 = (cse_p15_t1_g560133 * cse_p25_t1_g560274);
        double cse_p11_t58_g559590 = (cse_p16_t1_g560162 * cse_p25_t1_g560274);
        double cse_p26_t1_g560277 = (dz13 * dz23);
        double cse_p11_t80_g559612 = (cse_p15_t1_g560133 * cse_p26_t1_g560277);
        double cse_p11_t81_g559613 = (cse_p16_t1_g560162 * cse_p26_t1_g560277);
        double cse_p11_t109_g559641 = (cse_p12_t111_g559863 + cse_p12_t53_g559805);
        double cse_p11_t110_g559642 = (cse_p12_t115_g559867 + cse_p12_t54_g559806);
        double cse_p11_t111_g559643 = (cse_p12_t55_g559807 + cse_p12_t85_g559837);
        double cse_p11_t112_g559644 = (cse_p12_t56_g559808 + cse_p12_t84_g559836);
        double cse_p11_t113_g559645 = (cse_p12_t112_g559864 + cse_p12_t57_g559809);
        double cse_p11_t114_g559646 = (cse_p12_t114_g559866 + cse_p12_t58_g559810);
        double cse_p11_t115_g559647 = (cse_p12_t113_g559865 + cse_p12_t59_g559811);
        double cse_p11_t116_g559648 = (cse_p12_t60_g559812 + cse_p12_t83_g559835);
        double cse_p11_t117_g559649 = (cse_p12_t110_g559862 + cse_p12_t61_g559813);
        double cse_p11_t118_g559650 = (cse_p12_t28_g559780 + cse_p12_t41_g559793 + cse_p12_t47_g559799);
        double cse_p11_t119_g559651 = (cse_p12_t29_g559781 + cse_p12_t43_g559795 + cse_p12_t48_g559800);
        double cse_p11_t120_g559652 = (cse_p12_t30_g559782 + cse_p12_t45_g559797 + cse_p12_t49_g559801);
        double cse_p11_t124_g559656 = (cse_p12_t4_g559756 + cse_p12_t9_g559761);
        double cse_p11_t125_g559657 = (cse_p12_t5_g559757 + cse_p12_t8_g559760);
        double cse_p11_t126_g559658 = (cse_p12_t6_g559758 + cse_p12_t7_g559759);
        double cse_p11_t127_g559659 = (cse_p12_t107_g559859 + cse_p12_t62_g559814);
        double cse_p11_t128_g559660 = (cse_p12_t108_g559860 + cse_p12_t63_g559815);
        double cse_p11_t129_g559661 = (cse_p12_t109_g559861 + cse_p12_t64_g559816);
        double cse_p11_t131_g559663 = (cse_p12_t38_g559790 + cse_p12_t42_g559794 + cse_p12_t50_g559802);
        double cse_p11_t137_g559669 = (cse_p12_t39_g559791 + cse_p12_t44_g559796 + cse_p12_t51_g559803);
        double cse_p11_t141_g559673 = (cse_p12_t40_g559792 + cse_p12_t46_g559798 + cse_p12_t52_g559804);
        double cse_p11_t156_g559688 = (-(cse_p12_t10_g559762));
        double cse_p11_t157_g559689 = (-(cse_p12_t11_g559763));
        double cse_p11_t164_g559696 = (-(cse_p12_t14_g559766));
        double cse_p11_t173_g559705 = (-(cse_p18_t1_g560200));
        double cse_p11_t174_g559706 = (-(cse_p20_t1_g560219));
        double cse_p11_t175_g559707 = (-(cse_p20_t1_g560217));
        double cse_p11_t176_g559708 = (-(cse_p20_t1_g560218));
        double cse_p11_t182_g559714 = (-(dx43));
        double cse_p11_t185_g559717 = (-(dy43));
        double cse_p11_t188_g559720 = (-(dz43));
        double cse_p11_t192_g559724 = ((cse_p11_t142_g559674) * (cse_p11_t142_g559674));
        double energy = chiral_e;
        energy_added += energy;
        *energy_accumulate += energy;
        double cse_p11_t1_g559742 = (3.00000000000000000e+0 * cse_p11_t192_g559724 * term.k);
        double cse_p19_t1_g560204 = (cse_p11_t1_g559735 * cse_p1_t2_invsqrt559526);
        double cse_p15_t1_g560136 = (cse_p19_t1_g560204 * cse_p63_t1_g560305);
        double g_x1 = (cse_p11_t1_g559742 * (cse_p12_t59_g559811 + (-((cse_p15_t1_g560136 * dx13)))));
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = (cse_p11_t1_g559742 * (cse_p12_t58_g559810 + (-((cse_p15_t1_g560136 * dy13)))));
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = (cse_p11_t1_g559742 * (cse_p12_t54_g559806 + (-((cse_p15_t1_g560136 * dz13)))));
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double cse_p19_t1_g560202 = (cse_p11_t1_g559734 * cse_p1_t4_invsqrt559528);
        double cse_p15_t1_g560135 = (cse_p19_t1_g560202 * cse_p63_t1_g560305);
        double g_x2 = (cse_p11_t1_g559742 * (cse_p12_t61_g559813 + (-((cse_p15_t1_g560135 * dx23)))));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (cse_p11_t1_g559742 * (cse_p12_t53_g559805 + (-((cse_p15_t1_g560135 * dy23)))));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (cse_p11_t1_g559742 * (cse_p12_t57_g559809 + (-((cse_p15_t1_g560135 * dz23)))));
        KernelGradientAcc(term.i3x2, 2, g_z2);
        double cse_p15_t1_g560132 = (cse_p11_t1_g559733 * dx43);
        double cse_p27_t1_g560287 = (cse_p15_t1_g560132 * cse_p1_t6_invsqrt559530);
        double g_x3 = (cse_p11_t1_g559742 * (cse_p12_t60_g559812 + (-((cse_p63_t1_g560305 * ((-((cse_p19_t1_g560202 * dx23))) + (-((cse_p19_t1_g560204 * dx13))) + (-(cse_p27_t1_g560287))))))));
        KernelGradientAcc(term.i3x3, 0, g_x3);
        double cse_p22_t1_g560238 = (cse_p11_t1_g559733 * cse_p1_t6_invsqrt559530);
        double g_y3 = (cse_p11_t1_g559742 * (cse_p12_t56_g559808 + (-((cse_p63_t1_g560305 * ((-((cse_p19_t1_g560202 * dy23))) + (-((cse_p19_t1_g560204 * dy13))) + (-((cse_p22_t1_g560238 * dy43)))))))));
        KernelGradientAcc(term.i3x3, 1, g_y3);
        double g_z3 = (cse_p11_t1_g559742 * (cse_p12_t55_g559807 + (-((cse_p63_t1_g560305 * ((-((cse_p19_t1_g560202 * dz23))) + (-((cse_p19_t1_g560204 * dz13))) + (-((cse_p22_t1_g560238 * dz43)))))))));
        KernelGradientAcc(term.i3x3, 2, g_z3);
        double g_x4 = (cse_p11_t1_g559742 * (cse_p12_t62_g559814 + (-((cse_p27_t1_g560287 * cse_p63_t1_g560305)))));
        KernelGradientAcc(term.i3x4, 0, g_x4);
        double cse_p15_t1_g560137 = (cse_p22_t1_g560238 * cse_p63_t1_g560305);
        double g_y4 = (cse_p11_t1_g559742 * (cse_p12_t63_g559815 + (-((cse_p15_t1_g560137 * dy43)))));
        KernelGradientAcc(term.i3x4, 1, g_y4);
        double g_z4 = (cse_p11_t1_g559742 * (cse_p12_t64_g559816 + (-((cse_p15_t1_g560137 * dz43)))));
        KernelGradientAcc(term.i3x4, 2, g_z4);
        double cse_p11_t1_g559745 = (cse_p17_t1_g560189 * cse_p61_t11_invr2560302);
        double cse_p12_t1_g559894 = (6.00000000000000000e+0 * cse_p11_t142_g559674 * term.k);
        double cse_p12_t1_g559895 = (2.00000000000000000e+0 * v);
        double cse_p12_t1_g559898 = (cse_p12_t1_g559895 * cse_p61_t12_invr3560303);
        double cse_p22_t1_g560241 = (cse_p12_t1_g559895 * cse_p61_t11_invr2560302);
        double cse_p12_t1_g559910 = (cse_p13_t1_g560016 * cse_p22_t1_g560241);
        double cse_p13_t1_g560006 = (cse_p11_t1_g559745 * dx13);
        double cse_p13_t1_g560007 = (cse_p12_t1_g559898 * r43);
        double cse_p16_t1_g560171 = (cse_p11_t115_g559647 * cse_p12_t1_g559894);
        double cse_p17_t1_g560181 = (cse_p520_t2_invr560355 * v);
        double cse_p19_t1_g560205 = (cse_p13_t1_g560007 * cse_p26_t1_g560283 * cse_p520_t1_sqrt560354 * r23);
        double cse_p24_t1_g560262 = (cse_p13_t1_g560006 * cse_p17_t1_g560181);
        double h_x1_x1 = ((cse_p11_t115_g559647 * cse_p16_t1_g560171) + (cse_p11_t1_g559742 * (cse_p11_t176_g559708 + (-2.00000000000000000e+0 * cse_p13_t1_g560006 * cse_p13_t49_g559962) + (cse_p12_t1_g559910 * cse_p15_t1_g560134) + (cse_p15_t1_g560134 * cse_p19_t1_g560205) + (-((cse_p24_t1_g560262 * dx13))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double cse_p13_t1_g560010 = (cse_p11_t1_g559745 * dy13);
        double cse_p15_t1_g560131 = (cse_p11_t114_g559646 * cse_p12_t1_g559894);
        double cse_p21_t1_g560233 = (cse_p11_t124_g559656 * cse_p19_t1_g560203 * v);
        double h_x1_y1 = ((cse_p11_t115_g559647 * cse_p15_t1_g560131) + (cse_p11_t1_g559742 * ((cse_p20_t1_g560221 * cse_p21_t1_g560233) + (-((cse_p13_t1_g560006 * cse_p13_t48_g559961))) + (-((cse_p13_t1_g560010 * cse_p13_t49_g559962))) + (-((cse_p24_t1_g560262 * dy13))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double cse_p14_t1_g560081 = (cse_p11_t1_g559745 * dz13);
        double cse_p14_t1_g560085 = (cse_p11_t110_g559642 * cse_p12_t1_g559894);
        double h_x1_z1 = ((cse_p11_t115_g559647 * cse_p14_t1_g560085) + (cse_p11_t1_g559742 * ((cse_p20_t1_g560223 * cse_p21_t1_g560233) + (-((cse_p13_t1_g560006 * cse_p13_t44_g559957))) + (-((cse_p13_t49_g559962 * cse_p14_t1_g560081))) + (-((cse_p24_t1_g560262 * dz13))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double cse_p11_t1_g559744 = (cse_p17_t1_g560188 * cse_p61_t11_invr2560302);
        double cse_p13_t1_g560005 = (cse_p11_t1_g559744 * dx23);
        double cse_p13_t1_g560013 = (cse_p11_t1_g559733 * cse_p13_t1_g560007 * cse_p15_t1_g560133);
        double h_x1_x2 = ((cse_p11_t117_g559649 * cse_p16_t1_g560171) + (cse_p11_t1_g559742 * ((cse_p13_t1_g560013 * cse_p25_t1_g560273) + (-((cse_p13_t1_g560005 * cse_p13_t49_g559962))) + (-((cse_p13_t1_g560006 * cse_p13_t51_g559964))) + (-(cse_p11_t32_g559564)))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double cse_p13_t1_g560009 = (cse_p11_t1_g559744 * dy23);
        double cse_p14_t1_g560084 = (cse_p11_t109_g559641 * cse_p12_t1_g559894);
        double h_x1_y2 = ((cse_p11_t115_g559647 * cse_p14_t1_g560084) + (cse_p11_t1_g559742 * ((cse_p13_t1_g560013 * cse_p20_t1_g560222) + (cse_p61_t10_invr560301 * dz43) + (-((cse_p13_t1_g560006 * cse_p13_t43_g559956))) + (-((cse_p13_t1_g560009 * cse_p13_t49_g559962))) + (-((cse_p16_t1_g560162 * cse_p20_t1_g560222))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double cse_p13_t1_g560012 = (cse_p11_t1_g559744 * dz23);
        double cse_p15_t1_g560130 = (cse_p11_t113_g559645 * cse_p12_t1_g559894);
        double h_x1_z2 = ((cse_p11_t115_g559647 * cse_p15_t1_g560130) + (cse_p11_t1_g559742 * ((cse_p11_t185_g559717 * cse_p61_t10_invr560301) + (cse_p13_t1_g560013 * cse_p20_t1_g560224) + (-((cse_p13_t1_g560006 * cse_p13_t47_g559960))) + (-((cse_p13_t1_g560012 * cse_p13_t49_g559962))) + (-((cse_p16_t1_g560162 * cse_p20_t1_g560224))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double cse_p12_t1_g559907 = (cse_p12_t1_g559898 * cse_p17_t1_g560189);
        double cse_p19_t1_g560207 = (cse_p19_t1_g560209 * cse_p63_t1_g560305);
        double cse_p23_t1_g560251 = (cse_p12_t1_g559907 * cse_p15_t10_g560100);
        double cse_p24_t1_g560266 = (cse_p520_t2_invr560355 * r43);
        double h_x1_x3 = ((cse_p11_t116_g559648 * cse_p16_t1_g560171) + (cse_p11_t1_g559742 * (cse_p11_t32_g559564 + (cse_p13_t49_g559962 * cse_p13_t61_g559974) + (cse_p15_t1_g560134 * cse_p17_t1_g560185) + (cse_p19_t1_g560207 * (r43 + (cse_p14_t1_g560086 * dx13) + (cse_p15_t1_g560134 * cse_p24_t1_g560266))) + (cse_p23_t1_g560251 * dx13) + (-((cse_p13_t1_g560006 * cse_p13_t50_g559963))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 0, h_x1_x3);
        double cse_p24_t1_g560263 = (cse_p13_t4_g559917 * cse_p63_t1_g560305);
        double cse_p13_t1_g560020 = (cse_p11_t137_g559669 * cse_p24_t1_g560263);
        double cse_p15_t1_g560127 = (cse_p11_t112_g559644 * cse_p12_t1_g559894);
        double cse_p23_t1_g560252 = (cse_p12_t1_g559907 * cse_p15_t11_g560101);
        double h_x1_y3 = ((cse_p11_t115_g559647 * cse_p15_t1_g560127) + (cse_p11_t1_g559742 * (cse_p11_t4_g559536 + (cse_p13_t1_g560020 * dx13) + (cse_p13_t49_g559962 * cse_p13_t62_g559975) + (cse_p23_t1_g560252 * dx13) + (cse_p61_t10_invr560301 * (cse_p11_t188_g559720 + dz23)) + (-((cse_p13_t1_g560006 * cse_p13_t46_g559959))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 1, h_x1_y3);
        double cse_p13_t1_g560021 = (cse_p11_t141_g559673 * cse_p24_t1_g560263);
        double cse_p14_t1_g560089 = (cse_p11_t111_g559643 * cse_p12_t1_g559894);
        double cse_p23_t1_g560253 = (cse_p12_t1_g559907 * cse_p15_t12_g560102);
        double h_x1_z3 = ((cse_p11_t115_g559647 * cse_p14_t1_g560089) + (cse_p11_t1_g559742 * (cse_p11_t5_g559537 + (cse_p13_t1_g560021 * dx13) + (cse_p13_t49_g559962 * cse_p13_t63_g559976) + (cse_p23_t1_g560253 * dx13) + (cse_p61_t10_invr560301 * (cse_p16_t22_g560159 + dy43)) + (-((cse_p13_t1_g560006 * cse_p13_t45_g559958))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 2, h_x1_z3);
        double cse_p11_t1_g559743 = (cse_p18_t1_g560200 * cse_p61_t11_invr2560302);
        double cse_p11_t1_g559746 = (cse_p18_t1_g560200 * cse_p61_t12_invr3560303);
        double cse_p17_t1_g560190 = (cse_p11_t1_g559746 * cse_p12_t1_g559895);
        double cse_p12_t1_g559905 = (cse_p17_t1_g560189 * cse_p17_t1_g560190);
        double cse_p12_t1_g559908 = (cse_p13_t6_g559919 * cse_p19_t1_g560207);
        double cse_p13_t1_g560004 = (cse_p11_t1_g559743 * dx43);
        double cse_p22_t1_g560245 = (cse_p12_t1_g559905 * dx13);
        double cse_p23_t1_g560254 = (cse_p12_t1_g559908 * dx13);
        double h_x1_x4 = ((cse_p11_t127_g559659 * cse_p16_t1_g560171) + (cse_p11_t1_g559742 * ((cse_p22_t1_g560245 * dx43) + (-((cse_p13_t1_g560004 * cse_p13_t49_g559962))) + (-((cse_p13_t1_g560006 * cx))) + (-((cse_p23_t1_g560254 * dx43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x4, 0, h_x1_x4);
        double cse_p13_t1_g560008 = (cse_p11_t1_g559743 * dy43);
        double h_x1_y4 = ((cse_p11_t128_g559660 * cse_p16_t1_g560171) + (cse_p11_t1_g559742 * ((cse_p16_t24_g560161 * cse_p61_t10_invr560301) + (cse_p22_t1_g560245 * dy43) + (-((cse_p13_t1_g560006 * cy))) + (-((cse_p13_t1_g560008 * cse_p13_t49_g559962))) + (-((cse_p23_t1_g560254 * dy43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x4, 1, h_x1_y4);
        double cse_p13_t1_g560011 = (cse_p11_t1_g559743 * dz43);
        double h_x1_z4 = ((cse_p11_t129_g559661 * cse_p16_t1_g560171) + (cse_p11_t1_g559742 * ((cse_p22_t1_g560245 * dz43) + (cse_p61_t10_invr560301 * dy23) + (-((cse_p13_t1_g560006 * cz))) + (-((cse_p13_t1_g560011 * cse_p13_t49_g559962))) + (-((cse_p23_t1_g560254 * dz43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x4, 2, h_x1_z4);
        double cse_p26_t1_g560282 = (cse_p13_t1_g560010 * cse_p17_t1_g560181);
        double h_y1_y1 = ((cse_p11_t114_g559646 * cse_p15_t1_g560131) + (cse_p11_t1_g559742 * (cse_p11_t176_g559708 + (-2.00000000000000000e+0 * cse_p13_t1_g560010 * cse_p13_t48_g559961) + (cse_p12_t1_g559910 * cse_p16_t1_g560164) + (cse_p16_t1_g560164 * cse_p19_t1_g560205) + (-((cse_p26_t1_g560282 * dy13))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = ((cse_p11_t114_g559646 * cse_p14_t1_g560085) + (cse_p11_t1_g559742 * ((cse_p21_t1_g560229 * cse_p21_t1_g560233) + (-((cse_p13_t1_g560010 * cse_p13_t44_g559957))) + (-((cse_p13_t48_g559961 * cse_p14_t1_g560081))) + (-((cse_p26_t1_g560282 * dz13))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = ((cse_p11_t117_g559649 * cse_p15_t1_g560131) + (cse_p11_t1_g559742 * ((cse_p11_t188_g559720 * cse_p61_t10_invr560301) + (cse_p13_t1_g560013 * cse_p20_t1_g560225) + (-((cse_p13_t1_g560005 * cse_p13_t48_g559961))) + (-((cse_p13_t1_g560010 * cse_p13_t51_g559964))) + (-((cse_p16_t1_g560162 * cse_p20_t1_g560225))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = ((cse_p11_t114_g559646 * cse_p14_t1_g560084) + (cse_p11_t1_g559742 * ((cse_p13_t1_g560013 * cse_p25_t1_g560274) + (-((cse_p13_t1_g560009 * cse_p13_t48_g559961))) + (-((cse_p13_t1_g560010 * cse_p13_t43_g559956))) + (-(cse_p11_t58_g559590)))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = ((cse_p11_t114_g559646 * cse_p15_t1_g560130) + (cse_p11_t1_g559742 * ((cse_p13_t1_g560013 * cse_p21_t1_g560230) + (cse_p61_t10_invr560301 * dx43) + (-((cse_p13_t1_g560010 * cse_p13_t47_g559960))) + (-((cse_p13_t1_g560012 * cse_p13_t48_g559961))) + (-((cse_p16_t1_g560162 * cse_p21_t1_g560230))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double cse_p13_t1_g560019 = (cse_p11_t131_g559663 * cse_p24_t1_g560263);
        double h_y1_x3 = ((cse_p11_t116_g559648 * cse_p15_t1_g560131) + (cse_p11_t1_g559742 * (cse_p11_t4_g559536 + (cse_p13_t1_g560019 * dy13) + (cse_p13_t48_g559961 * cse_p13_t61_g559974) + (cse_p23_t1_g560251 * dy13) + (cse_p61_t10_invr560301 * (cse_p16_t24_g560161 + dz43)) + (-((cse_p13_t1_g560010 * cse_p13_t50_g559963))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 0, h_y1_x3);
        double h_y1_y3 = ((cse_p11_t114_g559646 * cse_p15_t1_g560127) + (cse_p11_t1_g559742 * (cse_p11_t58_g559590 + (cse_p13_t48_g559961 * cse_p13_t62_g559975) + (cse_p16_t1_g560164 * cse_p17_t1_g560185) + (cse_p19_t1_g560207 * (r43 + (cse_p14_t1_g560087 * dy13) + (cse_p16_t1_g560164 * cse_p24_t1_g560266))) + (cse_p23_t1_g560252 * dy13) + (-((cse_p13_t1_g560010 * cse_p13_t46_g559959))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 1, h_y1_y3);
        double h_y1_z3 = ((cse_p11_t114_g559646 * cse_p14_t1_g560089) + (cse_p11_t1_g559742 * (cse_p11_t8_g559540 + (cse_p13_t1_g560021 * dy13) + (cse_p13_t48_g559961 * cse_p13_t63_g559976) + (cse_p23_t1_g560253 * dy13) + (cse_p61_t10_invr560301 * (cse_p11_t182_g559714 + dx23)) + (-((cse_p13_t1_g560010 * cse_p13_t45_g559958))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 2, h_y1_z3);
        double cse_p22_t1_g560246 = (cse_p12_t1_g559905 * dy13);
        double cse_p23_t1_g560255 = (cse_p12_t1_g559908 * dy13);
        double h_y1_x4 = ((cse_p11_t127_g559659 * cse_p15_t1_g560131) + (cse_p11_t1_g559742 * ((cse_p22_t1_g560246 * dx43) + (cse_p61_t10_invr560301 * dz23) + (-((cse_p13_t1_g560004 * cse_p13_t48_g559961))) + (-((cse_p13_t1_g560010 * cx))) + (-((cse_p23_t1_g560255 * dx43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x4, 0, h_y1_x4);
        double h_y1_y4 = ((cse_p11_t128_g559660 * cse_p15_t1_g560131) + (cse_p11_t1_g559742 * ((cse_p22_t1_g560246 * dy43) + (-((cse_p13_t1_g560008 * cse_p13_t48_g559961))) + (-((cse_p13_t1_g560010 * cy))) + (-((cse_p23_t1_g560255 * dy43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x4, 1, h_y1_y4);
        double h_y1_z4 = ((cse_p11_t129_g559661 * cse_p15_t1_g560131) + (cse_p11_t1_g559742 * ((cse_p16_t20_g560157 * cse_p61_t10_invr560301) + (cse_p22_t1_g560246 * dz43) + (-((cse_p13_t1_g560010 * cz))) + (-((cse_p13_t1_g560011 * cse_p13_t48_g559961))) + (-((cse_p23_t1_g560255 * dz43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x4, 2, h_y1_z4);
        double h_z1_z1 = ((cse_p11_t110_g559642 * cse_p14_t1_g560085) + (cse_p11_t1_g559742 * (cse_p11_t176_g559708 + (-2.00000000000000000e+0 * cse_p13_t44_g559957 * cse_p14_t1_g560081) + (cse_p12_t1_g559910 * cse_p16_t1_g560166) + (cse_p16_t1_g560166 * cse_p19_t1_g560205) + (-((cse_p14_t1_g560081 * cse_p17_t1_g560181 * dz13))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = ((cse_p11_t117_g559649 * cse_p14_t1_g560085) + (cse_p11_t1_g559742 * ((cse_p13_t1_g560013 * cse_p21_t1_g560227) + (cse_p61_t10_invr560301 * dy43) + (-((cse_p13_t1_g560005 * cse_p13_t44_g559957))) + (-((cse_p13_t51_g559964 * cse_p14_t1_g560081))) + (-((cse_p16_t1_g560162 * cse_p21_t1_g560227))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = ((cse_p11_t110_g559642 * cse_p14_t1_g560084) + (cse_p11_t1_g559742 * ((cse_p11_t182_g559714 * cse_p61_t10_invr560301) + (cse_p13_t1_g560013 * cse_p21_t1_g560231) + (-((cse_p13_t1_g560009 * cse_p13_t44_g559957))) + (-((cse_p13_t43_g559956 * cse_p14_t1_g560081))) + (-((cse_p16_t1_g560162 * cse_p21_t1_g560231))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = ((cse_p11_t113_g559645 * cse_p14_t1_g560085) + (cse_p11_t1_g559742 * ((cse_p13_t1_g560013 * cse_p26_t1_g560277) + (-((cse_p13_t1_g560012 * cse_p13_t44_g559957))) + (-((cse_p13_t47_g559960 * cse_p14_t1_g560081))) + (-(cse_p11_t81_g559613)))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_z1_x3 = ((cse_p11_t116_g559648 * cse_p14_t1_g560085) + (cse_p11_t1_g559742 * (cse_p11_t5_g559537 + (cse_p13_t1_g560019 * dz13) + (cse_p13_t44_g559957 * cse_p13_t61_g559974) + (cse_p23_t1_g560251 * dz13) + (cse_p61_t10_invr560301 * (cse_p11_t185_g559717 + dy23)) + (-((cse_p13_t50_g559963 * cse_p14_t1_g560081))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 0, h_z1_x3);
        double h_z1_y3 = ((cse_p11_t112_g559644 * cse_p14_t1_g560085) + (cse_p11_t1_g559742 * (cse_p11_t8_g559540 + (cse_p13_t1_g560020 * dz13) + (cse_p13_t44_g559957 * cse_p13_t62_g559975) + (cse_p23_t1_g560252 * dz13) + (cse_p61_t10_invr560301 * (cse_p16_t20_g560157 + dx43)) + (-((cse_p13_t46_g559959 * cse_p14_t1_g560081))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 1, h_z1_y3);
        double h_z1_z3 = ((cse_p11_t111_g559643 * cse_p14_t1_g560085) + (cse_p11_t1_g559742 * (cse_p11_t81_g559613 + (cse_p13_t44_g559957 * cse_p13_t63_g559976) + (cse_p16_t1_g560166 * cse_p17_t1_g560185) + (cse_p19_t1_g560207 * (r43 + (cse_p14_t1_g560088 * dz13) + (cse_p16_t1_g560166 * cse_p24_t1_g560266))) + (cse_p23_t1_g560253 * dz13) + (-((cse_p13_t45_g559958 * cse_p14_t1_g560081))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 2, h_z1_z3);
        double cse_p23_t1_g560247 = (cse_p12_t1_g559905 * dz13);
        double cse_p23_t1_g560256 = (cse_p12_t1_g559908 * dz13);
        double h_z1_x4 = ((cse_p11_t127_g559659 * cse_p14_t1_g560085) + (cse_p11_t1_g559742 * ((cse_p16_t22_g560159 * cse_p61_t10_invr560301) + (cse_p23_t1_g560247 * dx43) + (-((cse_p13_t1_g560004 * cse_p13_t44_g559957))) + (-((cse_p14_t1_g560081 * cx))) + (-((cse_p23_t1_g560256 * dx43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x4, 0, h_z1_x4);
        double h_z1_y4 = ((cse_p11_t128_g559660 * cse_p14_t1_g560085) + (cse_p11_t1_g559742 * ((cse_p23_t1_g560247 * dy43) + (cse_p61_t10_invr560301 * dx23) + (-((cse_p13_t1_g560008 * cse_p13_t44_g559957))) + (-((cse_p14_t1_g560081 * cy))) + (-((cse_p23_t1_g560256 * dy43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x4, 1, h_z1_y4);
        double h_z1_z4 = ((cse_p11_t129_g559661 * cse_p14_t1_g560085) + (cse_p11_t1_g559742 * ((cse_p23_t1_g560247 * dz43) + (-((cse_p13_t1_g560011 * cse_p13_t44_g559957))) + (-((cse_p14_t1_g560081 * cz))) + (-((cse_p23_t1_g560256 * dz43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x4, 2, h_z1_z4);
        double cse_p12_t1_g559911 = (cse_p13_t1_g560015 * cse_p22_t1_g560241);
        double cse_p17_t1_g560182 = (cse_p520_t4_invr560357 * v);
        double cse_p18_t1_g560199 = (cse_p11_t117_g559649 * cse_p12_t1_g559894);
        double cse_p19_t1_g560206 = (cse_p13_t1_g560007 * cse_p26_t1_g560284 * cse_p520_t3_sqrt560356 * r13);
        double cse_p24_t1_g560261 = (cse_p13_t1_g560005 * cse_p17_t1_g560182);
        double h_x2_x2 = ((cse_p11_t117_g559649 * cse_p18_t1_g560199) + (cse_p11_t1_g559742 * (cse_p11_t175_g559707 + (-2.00000000000000000e+0 * cse_p13_t1_g560005 * cse_p13_t51_g559964) + (cse_p12_t1_g559911 * cse_p16_t1_g560163) + (cse_p16_t1_g560163 * cse_p19_t1_g560206) + (-((cse_p24_t1_g560261 * dx23))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double cse_p21_t1_g560234 = (cse_p11_t125_g559657 * cse_p19_t1_g560201 * v);
        double h_x2_y2 = ((cse_p11_t117_g559649 * cse_p14_t1_g560084) + (cse_p11_t1_g559742 * ((cse_p20_t1_g560226 * cse_p21_t1_g560234) + (-((cse_p13_t1_g560005 * cse_p13_t43_g559956))) + (-((cse_p13_t1_g560009 * cse_p13_t51_g559964))) + (-((cse_p24_t1_g560261 * dy23))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = ((cse_p11_t117_g559649 * cse_p15_t1_g560130) + (cse_p11_t1_g559742 * ((cse_p21_t1_g560228 * cse_p21_t1_g560234) + (-((cse_p13_t1_g560005 * cse_p13_t47_g559960))) + (-((cse_p13_t1_g560012 * cse_p13_t51_g559964))) + (-((cse_p24_t1_g560261 * dz23))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double cse_p12_t1_g559906 = (cse_p12_t1_g559898 * cse_p17_t1_g560188);
        double cse_p17_t1_g560187 = (cse_p11_t116_g559648 * cse_p12_t1_g559894);
        double cse_p19_t1_g560208 = (cse_p19_t1_g560210 * cse_p63_t1_g560305);
        double cse_p23_t1_g560248 = (cse_p12_t1_g559906 * cse_p15_t10_g560100);
        double cse_p25_t1_g560267 = (cse_p520_t4_invr560357 * r43);
        double h_x2_x3 = ((cse_p11_t117_g559649 * cse_p17_t1_g560187) + (cse_p11_t1_g559742 * (cse_p11_t32_g559564 + (cse_p13_t51_g559964 * cse_p13_t61_g559974) + (cse_p16_t1_g560163 * cse_p17_t1_g560184) + (cse_p19_t1_g560208 * (r43 + (cse_p14_t1_g560086 * dx23) + (cse_p16_t1_g560163 * cse_p25_t1_g560267))) + (cse_p23_t1_g560248 * dx23) + (-((cse_p13_t1_g560005 * cse_p13_t50_g559963))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 0, h_x2_x3);
        double cse_p24_t1_g560264 = (cse_p13_t5_g559918 * cse_p63_t1_g560305);
        double cse_p13_t1_g560017 = (cse_p11_t119_g559651 * cse_p24_t1_g560264);
        double cse_p23_t1_g560249 = (cse_p12_t1_g559906 * cse_p15_t11_g560101);
        double h_x2_y3 = ((cse_p11_t117_g559649 * cse_p15_t1_g560127) + (cse_p11_t1_g559742 * (cse_p11_t6_g559538 + (cse_p13_t1_g560017 * dx23) + (cse_p13_t51_g559964 * cse_p13_t62_g559975) + (cse_p23_t1_g560249 * dx23) + (cse_p61_t10_invr560301 * (cse_p16_t23_g560160 + dz43)) + (-((cse_p13_t1_g560005 * cse_p13_t46_g559959))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 1, h_x2_y3);
        double cse_p13_t1_g560018 = (cse_p11_t120_g559652 * cse_p24_t1_g560264);
        double cse_p23_t1_g560250 = (cse_p12_t1_g559906 * cse_p15_t12_g560102);
        double h_x2_z3 = ((cse_p11_t117_g559649 * cse_p14_t1_g560089) + (cse_p11_t1_g559742 * (cse_p11_t7_g559539 + (cse_p13_t1_g560018 * dx23) + (cse_p13_t51_g559964 * cse_p13_t63_g559976) + (cse_p23_t1_g560250 * dx23) + (cse_p61_t10_invr560301 * (cse_p11_t185_g559717 + dy13)) + (-((cse_p13_t1_g560005 * cse_p13_t45_g559958))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 2, h_x2_z3);
        double cse_p12_t1_g559904 = (cse_p17_t1_g560188 * cse_p17_t1_g560190);
        double cse_p12_t1_g559909 = (cse_p13_t6_g559919 * cse_p19_t1_g560208);
        double cse_p22_t1_g560242 = (cse_p12_t1_g559904 * dx23);
        double cse_p24_t1_g560257 = (cse_p12_t1_g559909 * dx23);
        double h_x2_x4 = ((cse_p11_t127_g559659 * cse_p18_t1_g560199) + (cse_p11_t1_g559742 * ((cse_p22_t1_g560242 * dx43) + (-((cse_p13_t1_g560004 * cse_p13_t51_g559964))) + (-((cse_p13_t1_g560005 * cx))) + (-((cse_p24_t1_g560257 * dx43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x4, 0, h_x2_x4);
        double h_x2_y4 = ((cse_p11_t128_g559660 * cse_p18_t1_g560199) + (cse_p11_t1_g559742 * ((cse_p22_t1_g560242 * dy43) + (cse_p61_t10_invr560301 * dz13) + (-((cse_p13_t1_g560005 * cy))) + (-((cse_p13_t1_g560008 * cse_p13_t51_g559964))) + (-((cse_p24_t1_g560257 * dy43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x4, 1, h_x2_y4);
        double h_x2_z4 = ((cse_p11_t129_g559661 * cse_p18_t1_g560199) + (cse_p11_t1_g559742 * ((cse_p16_t21_g560158 * cse_p61_t10_invr560301) + (cse_p22_t1_g560242 * dz43) + (-((cse_p13_t1_g560005 * cz))) + (-((cse_p13_t1_g560011 * cse_p13_t51_g559964))) + (-((cse_p24_t1_g560257 * dz43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x4, 2, h_x2_z4);
        double cse_p26_t1_g560281 = (cse_p13_t1_g560009 * cse_p17_t1_g560182);
        double h_y2_y2 = ((cse_p11_t109_g559641 * cse_p14_t1_g560084) + (cse_p11_t1_g559742 * (cse_p11_t175_g559707 + (-2.00000000000000000e+0 * cse_p13_t1_g560009 * cse_p13_t43_g559956) + (cse_p12_t1_g559911 * cse_p16_t1_g560165) + (cse_p16_t1_g560165 * cse_p19_t1_g560206) + (-((cse_p26_t1_g560281 * dy23))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = ((cse_p11_t113_g559645 * cse_p14_t1_g560084) + (cse_p11_t1_g559742 * ((cse_p21_t1_g560232 * cse_p21_t1_g560234) + (-((cse_p13_t1_g560009 * cse_p13_t47_g559960))) + (-((cse_p13_t1_g560012 * cse_p13_t43_g559956))) + (-((cse_p26_t1_g560281 * dz23))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double cse_p12_t1_g559913 = (cse_p11_t118_g559650 * cse_p24_t1_g560264);
        double h_y2_x3 = ((cse_p11_t116_g559648 * cse_p14_t1_g560084) + (cse_p11_t1_g559742 * (cse_p11_t6_g559538 + (cse_p12_t1_g559913 * dy23) + (cse_p13_t43_g559956 * cse_p13_t61_g559974) + (cse_p23_t1_g560248 * dy23) + (cse_p61_t10_invr560301 * (cse_p11_t188_g559720 + dz13)) + (-((cse_p13_t1_g560009 * cse_p13_t50_g559963))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 0, h_y2_x3);
        double h_y2_y3 = ((cse_p11_t112_g559644 * cse_p14_t1_g560084) + (cse_p11_t1_g559742 * (cse_p11_t58_g559590 + (cse_p13_t43_g559956 * cse_p13_t62_g559975) + (cse_p16_t1_g560165 * cse_p17_t1_g560184) + (cse_p19_t1_g560208 * (r43 + (cse_p14_t1_g560087 * dy23) + (cse_p16_t1_g560165 * cse_p25_t1_g560267))) + (cse_p23_t1_g560249 * dy23) + (-((cse_p13_t1_g560009 * cse_p13_t46_g559959))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 1, h_y2_y3);
        double h_y2_z3 = ((cse_p11_t111_g559643 * cse_p14_t1_g560084) + (cse_p11_t1_g559742 * (cse_p11_t9_g559541 + (cse_p13_t1_g560018 * dy23) + (cse_p13_t43_g559956 * cse_p13_t63_g559976) + (cse_p23_t1_g560250 * dy23) + (cse_p61_t10_invr560301 * (cse_p16_t19_g560156 + dx43)) + (-((cse_p13_t1_g560009 * cse_p13_t45_g559958))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 2, h_y2_z3);
        double cse_p22_t1_g560243 = (cse_p12_t1_g559904 * dy23);
        double cse_p24_t1_g560258 = (cse_p12_t1_g559909 * dy23);
        double h_y2_x4 = ((cse_p11_t127_g559659 * cse_p14_t1_g560084) + (cse_p11_t1_g559742 * ((cse_p16_t23_g560160 * cse_p61_t10_invr560301) + (cse_p22_t1_g560243 * dx43) + (-((cse_p13_t1_g560004 * cse_p13_t43_g559956))) + (-((cse_p13_t1_g560009 * cx))) + (-((cse_p24_t1_g560258 * dx43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x4, 0, h_y2_x4);
        double h_y2_y4 = ((cse_p11_t128_g559660 * cse_p14_t1_g560084) + (cse_p11_t1_g559742 * ((cse_p22_t1_g560243 * dy43) + (-((cse_p13_t1_g560008 * cse_p13_t43_g559956))) + (-((cse_p13_t1_g560009 * cy))) + (-((cse_p24_t1_g560258 * dy43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x4, 1, h_y2_y4);
        double h_y2_z4 = ((cse_p11_t129_g559661 * cse_p14_t1_g560084) + (cse_p11_t1_g559742 * ((cse_p22_t1_g560243 * dz43) + (cse_p61_t10_invr560301 * dx13) + (-((cse_p13_t1_g560009 * cz))) + (-((cse_p13_t1_g560011 * cse_p13_t43_g559956))) + (-((cse_p24_t1_g560258 * dz43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x4, 2, h_y2_z4);
        double h_z2_z2 = ((cse_p11_t113_g559645 * cse_p15_t1_g560130) + (cse_p11_t1_g559742 * (cse_p11_t175_g559707 + (-2.00000000000000000e+0 * cse_p13_t1_g560012 * cse_p13_t47_g559960) + (cse_p12_t1_g559911 * cse_p16_t1_g560167) + (cse_p16_t1_g560167 * cse_p19_t1_g560206) + (-((cse_p13_t1_g560012 * cse_p17_t1_g560182 * dz23))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        double h_z2_x3 = ((cse_p11_t116_g559648 * cse_p15_t1_g560130) + (cse_p11_t1_g559742 * (cse_p11_t7_g559539 + (cse_p12_t1_g559913 * dz23) + (cse_p13_t47_g559960 * cse_p13_t61_g559974) + (cse_p23_t1_g560248 * dz23) + (cse_p61_t10_invr560301 * (cse_p16_t21_g560158 + dy43)) + (-((cse_p13_t1_g560012 * cse_p13_t50_g559963))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 0, h_z2_x3);
        double h_z2_y3 = ((cse_p11_t113_g559645 * cse_p15_t1_g560127) + (cse_p11_t1_g559742 * (cse_p11_t9_g559541 + (cse_p13_t1_g560017 * dz23) + (cse_p13_t47_g559960 * cse_p13_t62_g559975) + (cse_p23_t1_g560249 * dz23) + (cse_p61_t10_invr560301 * (cse_p11_t182_g559714 + dx13)) + (-((cse_p13_t1_g560012 * cse_p13_t46_g559959))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 1, h_z2_y3);
        double h_z2_z3 = ((cse_p11_t113_g559645 * cse_p14_t1_g560089) + (cse_p11_t1_g559742 * (cse_p11_t81_g559613 + (cse_p13_t47_g559960 * cse_p13_t63_g559976) + (cse_p16_t1_g560167 * cse_p17_t1_g560184) + (cse_p19_t1_g560208 * (r43 + (cse_p14_t1_g560088 * dz23) + (cse_p16_t1_g560167 * cse_p25_t1_g560267))) + (cse_p23_t1_g560250 * dz23) + (-((cse_p13_t1_g560012 * cse_p13_t45_g559958))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 2, h_z2_z3);
        double cse_p22_t1_g560244 = (cse_p12_t1_g559904 * dz23);
        double cse_p24_t1_g560259 = (cse_p12_t1_g559909 * dz23);
        double h_z2_x4 = ((cse_p11_t127_g559659 * cse_p15_t1_g560130) + (cse_p11_t1_g559742 * ((cse_p22_t1_g560244 * dx43) + (cse_p61_t10_invr560301 * dy13) + (-((cse_p13_t1_g560004 * cse_p13_t47_g559960))) + (-((cse_p13_t1_g560012 * cx))) + (-((cse_p24_t1_g560259 * dx43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x4, 0, h_z2_x4);
        double h_z2_y4 = ((cse_p11_t128_g559660 * cse_p15_t1_g560130) + (cse_p11_t1_g559742 * ((cse_p16_t19_g560156 * cse_p61_t10_invr560301) + (cse_p22_t1_g560244 * dy43) + (-((cse_p13_t1_g560008 * cse_p13_t47_g559960))) + (-((cse_p13_t1_g560012 * cy))) + (-((cse_p24_t1_g560259 * dy43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x4, 1, h_z2_y4);
        double h_z2_z4 = ((cse_p11_t129_g559661 * cse_p15_t1_g560130) + (cse_p11_t1_g559742 * ((cse_p22_t1_g560244 * dz43) + (-((cse_p13_t1_g560011 * cse_p13_t47_g559960))) + (-((cse_p13_t1_g560012 * cz))) + (-((cse_p24_t1_g560259 * dz43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x4, 2, h_z2_z4);
        double cse_p16_t1_g560168 = (-2.00000000000000000e+0 * cse_p13_t1_g560014);
        double cse_p16_t1_g560169 = (-2.00000000000000000e+0 * cse_p61_t12_invr3560303);
        double cse_p25_t1_g560269 = (cse_p15_t10_g560100 * cse_p16_t1_g560169);
        double h_x3_x3 = ((cse_p11_t116_g559648 * cse_p17_t1_g560187) + (cse_p11_t1_g559742 * ((2.00000000000000000e+0 * cse_p13_t50_g559963 * cse_p13_t61_g559974) + (v * ((-((cse_p15_t10_g560100 * cse_p25_t1_g560269))) + (-((cse_p61_t11_invr2560302 * (cse_p11_t31_g559563 + cse_p11_t31_g559563 + cse_p17_t1_g560188 + cse_p17_t1_g560189 + cse_p18_t1_g560200 + (cse_p14_t1_g560082 * cse_p16_t1_g560163) + (cse_p14_t1_g560083 * cse_p15_t1_g560134) + (cse_p14_t1_g560086 * (cse_p20_t1_g560211 + cse_p20_t1_g560211 + cse_p20_t4_g560214 + cse_p20_t4_g560214 + (cse_p15_t1_g560132 * cse_p520_t6_invr560359))) + (cse_p15_t1_g560134 * cse_p26_t1_g560283) + (cse_p16_t1_g560163 * cse_p26_t1_g560284) + (cse_p16_t1_g560168 * cse_p17_t1_g560186))))))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 0, h_x3_x3);
        double cse_p27_t1_g560289 = (cse_p16_t1_g560168 * dx43);
        double h_x3_y3 = ((cse_p11_t116_g559648 * cse_p15_t1_g560127) + (cse_p11_t1_g559742 * ((cse_p13_t46_g559959 * cse_p13_t61_g559974) + (cse_p13_t50_g559963 * cse_p13_t62_g559975) + (v * ((-((cse_p15_t11_g560101 * cse_p25_t1_g560269))) + (-((cse_p61_t11_invr2560302 * ((cse_p13_t6_g559919 * ((cse_p15_t1_g560132 * cse_p17_t1_g560183) + (cse_p20_t1_g560211 * dy43) + (cse_p20_t2_g560212 * dx43) + (cse_p20_t4_g560214 * dy43) + (cse_p20_t5_g560215 * dx43))) + (cse_p14_t1_g560082 * cse_p20_t1_g560226) + (cse_p14_t1_g560083 * cse_p20_t1_g560221) + (cse_p15_t1_g560133 * cse_p20_t1_g560222) + (cse_p15_t1_g560133 * cse_p20_t1_g560225) + (cse_p20_t1_g560221 * cse_p26_t1_g560283) + (cse_p20_t1_g560226 * cse_p26_t1_g560284) + (cse_p27_t1_g560289 * dy43))))))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 1, h_x3_y3);
        double h_x3_z3 = ((cse_p11_t116_g559648 * cse_p14_t1_g560089) + (cse_p11_t1_g559742 * ((cse_p13_t45_g559958 * cse_p13_t61_g559974) + (cse_p13_t50_g559963 * cse_p13_t63_g559976) + (v * ((-((cse_p15_t12_g560102 * cse_p25_t1_g560269))) + (-((cse_p61_t11_invr2560302 * ((cse_p13_t6_g559919 * ((cse_p15_t1_g560128 * cse_p15_t1_g560132) + (cse_p20_t1_g560211 * dz43) + (cse_p20_t3_g560213 * dx43) + (cse_p20_t4_g560214 * dz43) + (cse_p20_t6_g560216 * dx43))) + (cse_p14_t1_g560082 * cse_p21_t1_g560228) + (cse_p14_t1_g560083 * cse_p20_t1_g560223) + (cse_p15_t1_g560133 * cse_p20_t1_g560224) + (cse_p15_t1_g560133 * cse_p21_t1_g560227) + (cse_p20_t1_g560223 * cse_p26_t1_g560283) + (cse_p21_t1_g560228 * cse_p26_t1_g560284) + (cse_p27_t1_g560289 * dz43))))))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 2, h_x3_z3);
        double cse_p15_t1_g560126 = (-2.00000000000000000e+0 * cse_p11_t1_g559746);
        double cse_p25_t1_g560268 = (cse_p15_t10_g560100 * cse_p15_t1_g560126);
        double h_x3_x4 = ((cse_p11_t127_g559659 * cse_p17_t1_g560187) + (cse_p11_t1_g559742 * ((cse_p13_t61_g559974 * cx) + (v * ((-((cse_p25_t1_g560268 * dx43))) + (-((cse_p61_t11_invr2560302 * (cse_p11_t173_g559705 + (cse_p16_t1_g560170 * cse_p17_t1_g560186) + (-((cse_p14_t1_g560086 * cse_p20_t1_g560211))) + (-((cse_p14_t1_g560086 * cse_p20_t4_g560214))) + (-((cse_p17_t1_g560186 * cse_p26_t1_g560285))))))))) + (-((cse_p13_t1_g560004 * cse_p13_t50_g559963))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x4, 0, h_x3_x4);
        double h_x3_y4 = ((cse_p11_t128_g559660 * cse_p17_t1_g560187) + (cse_p11_t1_g559742 * ((cse_p13_t61_g559974 * cy) + (cse_p15_t18_g560108 * cse_p61_t10_invr560301) + (v * ((-((cse_p25_t1_g560268 * dy43))) + (-((cse_p61_t11_invr2560302 * (cse_p11_t13_g559545 + cse_p11_t156_g559688 + (-((cse_p14_t1_g560087 * cse_p20_t1_g560211))) + (-((cse_p14_t1_g560087 * cse_p20_t4_g560214))))))))) + (-((cse_p13_t1_g560008 * cse_p13_t50_g559963))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x4, 1, h_x3_y4);
        double h_x3_z4 = ((cse_p11_t129_g559661 * cse_p17_t1_g560187) + (cse_p11_t1_g559742 * ((cse_p13_t61_g559974 * cz) + (cse_p15_t15_g560105 * cse_p61_t10_invr560301) + (v * ((-((cse_p25_t1_g560268 * dz43))) + (-((cse_p61_t11_invr2560302 * (cse_p11_t14_g559546 + cse_p11_t157_g559689 + (-((cse_p14_t1_g560088 * cse_p20_t1_g560211))) + (-((cse_p14_t1_g560088 * cse_p20_t4_g560214))))))))) + (-((cse_p13_t1_g560011 * cse_p13_t50_g559963))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x4, 2, h_x3_z4);
        double cse_p26_t1_g560286 = (cse_p15_t11_g560101 * cse_p16_t1_g560169);
        double h_y3_y3 = ((cse_p11_t112_g559644 * cse_p15_t1_g560127) + (cse_p11_t1_g559742 * ((2.00000000000000000e+0 * cse_p13_t46_g559959 * cse_p13_t62_g559975) + (v * ((-((cse_p15_t11_g560101 * cse_p26_t1_g560286))) + (-((cse_p61_t11_invr2560302 * (cse_p11_t57_g559589 + cse_p11_t57_g559589 + cse_p17_t1_g560188 + cse_p17_t1_g560189 + cse_p18_t1_g560200 + (cse_p14_t1_g560082 * cse_p16_t1_g560165) + (cse_p14_t1_g560083 * cse_p16_t1_g560164) + (cse_p14_t1_g560087 * (cse_p20_t2_g560212 + cse_p20_t2_g560212 + cse_p20_t5_g560215 + cse_p20_t5_g560215 + (cse_p11_t1_g559733 * cse_p17_t1_g560183))) + (cse_p16_t1_g560164 * cse_p26_t1_g560283) + (cse_p16_t1_g560165 * cse_p26_t1_g560284) + (cse_p16_t1_g560168 * cse_p18_t1_g560197))))))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 1, h_y3_y3);
        double cse_p26_t1_g560279 = (cse_p11_t1_g559733 * cse_p15_t1_g560128);
        double h_y3_z3 = ((cse_p11_t112_g559644 * cse_p14_t1_g560089) + (cse_p11_t1_g559742 * ((cse_p13_t45_g559958 * cse_p13_t62_g559975) + (cse_p13_t46_g559959 * cse_p13_t63_g559976) + (v * ((-((cse_p15_t12_g560102 * cse_p26_t1_g560286))) + (-((cse_p61_t11_invr2560302 * ((cse_p13_t6_g559919 * ((cse_p20_t2_g560212 * dz43) + (cse_p20_t3_g560213 * dy43) + (cse_p20_t5_g560215 * dz43) + (cse_p20_t6_g560216 * dy43) + (cse_p26_t1_g560279 * dy43))) + (cse_p14_t1_g560082 * cse_p21_t1_g560232) + (cse_p14_t1_g560083 * cse_p21_t1_g560229) + (cse_p15_t1_g560133 * cse_p21_t1_g560230) + (cse_p15_t1_g560133 * cse_p21_t1_g560231) + (cse_p16_t1_g560168 * cse_p25_t1_g560275) + (cse_p21_t1_g560229 * cse_p26_t1_g560283) + (cse_p21_t1_g560232 * cse_p26_t1_g560284))))))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 2, h_y3_z3);
        double cse_p25_t1_g560270 = (cse_p15_t11_g560101 * cse_p15_t1_g560126);
        double h_y3_x4 = ((cse_p11_t127_g559659 * cse_p15_t1_g560127) + (cse_p11_t1_g559742 * ((cse_p13_t62_g559975 * cx) + (cse_p15_t17_g560107 * cse_p61_t10_invr560301) + (v * ((-((cse_p25_t1_g560270 * dx43))) + (-((cse_p61_t11_invr2560302 * (cse_p11_t13_g559545 + cse_p11_t156_g559688 + (-((cse_p14_t1_g560086 * cse_p20_t2_g560212))) + (-((cse_p14_t1_g560086 * cse_p20_t5_g560215))))))))) + (-((cse_p13_t1_g560004 * cse_p13_t46_g559959))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x4, 0, h_y3_x4);
        double h_y3_y4 = ((cse_p11_t128_g559660 * cse_p15_t1_g560127) + (cse_p11_t1_g559742 * ((cse_p13_t62_g559975 * cy) + (v * ((-((cse_p25_t1_g560270 * dy43))) + (-((cse_p61_t11_invr2560302 * (cse_p11_t173_g559705 + (cse_p16_t1_g560170 * cse_p18_t1_g560197) + (-((cse_p14_t1_g560087 * cse_p20_t2_g560212))) + (-((cse_p14_t1_g560087 * cse_p20_t5_g560215))) + (-((cse_p18_t1_g560197 * cse_p26_t1_g560285))))))))) + (-((cse_p13_t1_g560008 * cse_p13_t46_g559959))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x4, 1, h_y3_y4);
        double h_y3_z4 = ((cse_p11_t129_g559661 * cse_p15_t1_g560127) + (cse_p11_t1_g559742 * ((cse_p13_t62_g559975 * cz) + (cse_p15_t14_g560104 * cse_p61_t10_invr560301) + (v * ((-((cse_p25_t1_g560270 * dz43))) + (-((cse_p61_t11_invr2560302 * (cse_p11_t15_g559547 + cse_p11_t164_g559696 + (-((cse_p14_t1_g560088 * cse_p20_t2_g560212))) + (-((cse_p14_t1_g560088 * cse_p20_t5_g560215))))))))) + (-((cse_p13_t1_g560011 * cse_p13_t46_g559959))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x4, 2, h_y3_z4);
        double h_z3_z3 = ((cse_p11_t111_g559643 * cse_p14_t1_g560089) + (cse_p11_t1_g559742 * ((2.00000000000000000e+0 * cse_p13_t45_g559958 * cse_p13_t63_g559976) + (v * ((-((cse_p15_t12_g560102 * cse_p15_t12_g560102 * cse_p16_t1_g560169))) + (-((cse_p61_t11_invr2560302 * (cse_p11_t80_g559612 + cse_p11_t80_g559612 + cse_p17_t1_g560188 + cse_p17_t1_g560189 + cse_p18_t1_g560200 + (cse_p14_t1_g560082 * cse_p16_t1_g560167) + (cse_p14_t1_g560083 * cse_p16_t1_g560166) + (cse_p14_t1_g560088 * (cse_p20_t3_g560213 + cse_p20_t3_g560213 + cse_p20_t6_g560216 + cse_p20_t6_g560216 + cse_p26_t1_g560279)) + (cse_p16_t1_g560166 * cse_p26_t1_g560283) + (cse_p16_t1_g560167 * cse_p26_t1_g560284) + (cse_p16_t1_g560168 * cse_p18_t1_g560198))))))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x3, 2, h_z3_z3);
        double cse_p25_t1_g560271 = (cse_p15_t12_g560102 * cse_p15_t1_g560126);
        double h_z3_x4 = ((cse_p11_t127_g559659 * cse_p14_t1_g560089) + (cse_p11_t1_g559742 * ((cse_p13_t63_g559976 * cx) + (cse_p15_t16_g560106 * cse_p61_t10_invr560301) + (v * ((-((cse_p25_t1_g560271 * dx43))) + (-((cse_p61_t11_invr2560302 * (cse_p11_t14_g559546 + cse_p11_t157_g559689 + (-((cse_p14_t1_g560086 * cse_p20_t3_g560213))) + (-((cse_p14_t1_g560086 * cse_p20_t6_g560216))))))))) + (-((cse_p13_t1_g560004 * cse_p13_t45_g559958))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x4, 0, h_z3_x4);
        double h_z3_y4 = ((cse_p11_t128_g559660 * cse_p14_t1_g560089) + (cse_p11_t1_g559742 * ((cse_p13_t63_g559976 * cy) + (cse_p15_t13_g560103 * cse_p61_t10_invr560301) + (v * ((-((cse_p25_t1_g560271 * dy43))) + (-((cse_p61_t11_invr2560302 * (cse_p11_t15_g559547 + cse_p11_t164_g559696 + (-((cse_p14_t1_g560087 * cse_p20_t3_g560213))) + (-((cse_p14_t1_g560087 * cse_p20_t6_g560216))))))))) + (-((cse_p13_t1_g560008 * cse_p13_t45_g559958))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x4, 1, h_z3_y4);
        double h_z3_z4 = ((cse_p11_t129_g559661 * cse_p14_t1_g560089) + (cse_p11_t1_g559742 * ((cse_p13_t63_g559976 * cz) + (v * ((-((cse_p25_t1_g560271 * dz43))) + (-((cse_p61_t11_invr2560302 * (cse_p11_t173_g559705 + (cse_p16_t1_g560170 * cse_p18_t1_g560198) + (-((cse_p14_t1_g560088 * cse_p20_t3_g560213))) + (-((cse_p14_t1_g560088 * cse_p20_t6_g560216))) + (-((cse_p18_t1_g560198 * cse_p26_t1_g560285))))))))) + (-((cse_p13_t1_g560011 * cse_p13_t45_g559958))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x4, 2, h_z3_z4);
        double cse_p12_t1_g559912 = (cse_p13_t1_g560014 * cse_p22_t1_g560241);
        double cse_p20_t1_g560220 = (cse_p17_t1_g560190 * cse_p520_t5_sqrt560358 * r13 * r23);
        double cse_p22_t1_g560237 = (cse_p11_t127_g559659 * cse_p12_t1_g559894);
        double cse_p24_t1_g560260 = (cse_p13_t1_g560004 * v);
        double h_x4_x4 = ((cse_p11_t127_g559659 * cse_p22_t1_g560237) + (cse_p11_t1_g559742 * (cse_p11_t174_g559706 + (-2.00000000000000000e+0 * cse_p13_t1_g560004 * cx) + (cse_p12_t1_g559912 * cse_p17_t1_g560186) + (cse_p17_t1_g560186 * cse_p20_t1_g560220 * cse_p520_t6_invr560359) + (-((cse_p24_t1_g560260 * cse_p520_t6_invr560359 * dx43))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 0, term.i3x4, 0, h_x4_x4);
        double cse_p25_t1_g560276 = (cse_p11_t126_g559658 * cse_p520_t5_sqrt560358 * v);
        double cse_p27_t1_g560288 = (cse_p15_t1_g560132 * cse_p25_t1_g560276);
        double h_x4_y4 = ((cse_p11_t128_g559660 * cse_p22_t1_g560237) + (cse_p11_t1_g559742 * ((cse_p27_t1_g560288 * dy43) + (-((cse_p13_t1_g560004 * cy))) + (-((cse_p13_t1_g560008 * cx))) + (-((cse_p17_t1_g560183 * cse_p24_t1_g560260))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 0, term.i3x4, 1, h_x4_y4);
        double h_x4_z4 = ((cse_p11_t129_g559661 * cse_p22_t1_g560237) + (cse_p11_t1_g559742 * ((cse_p27_t1_g560288 * dz43) + (-((cse_p13_t1_g560004 * cz))) + (-((cse_p13_t1_g560011 * cx))) + (-((cse_p15_t1_g560128 * cse_p24_t1_g560260))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 0, term.i3x4, 2, h_x4_z4);
        double cse_p26_t1_g560278 = (cse_p11_t128_g559660 * cse_p12_t1_g559894);
        double cse_p26_t1_g560280 = (cse_p13_t1_g560008 * v);
        double h_y4_y4 = ((cse_p11_t128_g559660 * cse_p26_t1_g560278) + (cse_p11_t1_g559742 * (cse_p11_t174_g559706 + (-2.00000000000000000e+0 * cse_p13_t1_g560008 * cy) + (cse_p12_t1_g559912 * cse_p18_t1_g560197) + (cse_p17_t1_g560183 * cse_p20_t1_g560220 * dy43) + (-((cse_p17_t1_g560183 * cse_p26_t1_g560280))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 1, term.i3x4, 1, h_y4_y4);
        double h_y4_z4 = ((cse_p11_t129_g559661 * cse_p26_t1_g560278) + (cse_p11_t1_g559742 * ((cse_p11_t1_g559733 * cse_p25_t1_g560275 * cse_p25_t1_g560276) + (-((cse_p13_t1_g560008 * cz))) + (-((cse_p13_t1_g560011 * cy))) + (-((cse_p15_t1_g560128 * cse_p26_t1_g560280))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 1, term.i3x4, 2, h_y4_z4);
        double h_z4_z4 = ((cse_p11_t129_g559661 * cse_p11_t129_g559661 * cse_p12_t1_g559894) + (cse_p11_t1_g559742 * (cse_p11_t174_g559706 + (-2.00000000000000000e+0 * cse_p13_t1_g560011 * cz) + (cse_p12_t1_g559912 * cse_p18_t1_g560198) + (cse_p15_t1_g560128 * cse_p20_t1_g560220 * dz43) + (-((cse_p13_t1_g560011 * cse_p15_t1_g560128 * v))))));
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
