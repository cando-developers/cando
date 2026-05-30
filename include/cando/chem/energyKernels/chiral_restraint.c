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
      double cse_p11_t1_g372478 = (-(x3));
      double cse_p11_t2_g372479 = (-(y3));
      double cse_p11_t3_g372480 = (-(z3));
      double dx13 = (cse_p11_t1_g372478 + x1);
      double dy13 = (cse_p11_t2_g372479 + y1);
      double dz13 = (cse_p11_t3_g372480 + z1);
      double dx23 = (cse_p11_t1_g372478 + x2);
      double dy23 = (cse_p11_t2_g372479 + y2);
      double dz23 = (cse_p11_t3_g372480 + z2);
      double dx43 = (cse_p11_t1_g372478 + x4);
      double dy43 = (cse_p11_t2_g372479 + y4);
      double dz43 = (cse_p11_t3_g372480 + z4);
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
      double cse_p11_t8_g372496 = (-(x3));
      double cse_p11_t9_g372497 = (-(y3));
      double cse_p11_t10_g372498 = (-(z3));
      double dx13 = (cse_p11_t8_g372496 + x1);
      double dy13 = (cse_p11_t9_g372497 + y1);
      double dz13 = (cse_p11_t10_g372498 + z1);
      double dx23 = (cse_p11_t8_g372496 + x2);
      double dy23 = (cse_p11_t9_g372497 + y2);
      double dz23 = (cse_p11_t10_g372498 + z2);
      double dx43 = (cse_p11_t8_g372496 + x4);
      double dy43 = (cse_p11_t9_g372497 + y4);
      double dz43 = (cse_p11_t10_g372498 + z4);
      double r13_2 = ((dx13 * dx13) + (dy13 * dy13) + (dz13 * dz13));
      double r23_2 = ((dx23 * dx23) + (dy23 * dy23) + (dz23 * dz23));
      double r43_2 = ((dx43 * dx43) + (dy43 * dy43) + (dz43 * dz43));
      double cse_p1_t1_sqrt372481 = sqrt(r13_2);
      double cse_p1_t2_invsqrt372482 = (1.0 / (cse_p1_t1_sqrt372481));
      double r13 = cse_p1_t1_sqrt372481;
      double cse_p1_t3_sqrt372483 = sqrt(r23_2);
      double cse_p1_t4_invsqrt372484 = (1.0 / (cse_p1_t3_sqrt372483));
      double r23 = cse_p1_t3_sqrt372483;
      double cse_p1_t5_sqrt372485 = sqrt(r43_2);
      double cse_p1_t6_invsqrt372486 = (1.0 / (cse_p1_t5_sqrt372485));
      double r43 = cse_p1_t5_sqrt372485;
      double cx = ((dy13 * dz23) + (-((dy23 * dz13))));
      double cy = ((dx23 * dz13) + (-((dx13 * dz23))));
      double cz = ((dx13 * dy23) + (-((dx23 * dy13))));
      double v = ((cx * dx43) + (cy * dy43) + (cz * dz43));
      double cse_p11_t1_g372508 = (r13 * r23);
      double denom = (cse_p11_t1_g372508 * r43);
      double cse_p1_t7_invr372487 = (1.0 / (denom));
      double cse_p1_t8_invr2372488 = (cse_p1_t7_invr372487 * cse_p1_t7_invr372487);
      double q = (cse_p1_t7_invr372487 * v);
      double cse_p11_t1_g372489 = (term.co + q);
      double chiral_e = (term.k * (cse_p11_t1_g372489 * cse_p11_t1_g372489 * cse_p11_t1_g372489));
      if ((chiral_e > 0.00000000000000000e+0)) {
        /* !ACTIVE */
        double cse_p11_t2_g372490 = (-(dx13));
        double cse_p11_t3_g372491 = (-(dx23));
        double cse_p11_t4_g372492 = (-(dy13));
        double cse_p11_t5_g372493 = (-(dy23));
        double cse_p11_t6_g372494 = (-(dz13));
        double cse_p11_t7_g372495 = (-(dz23));
        double cse_p11_t11_g372499 = ((cse_p11_t1_g372489) * (cse_p11_t1_g372489));
        double energy = chiral_e;
        energy_added += energy;
        *energy_accumulate += energy;
        double cse_p11_t1_g372500 = (3.00000000000000000e+0 * cse_p11_t11_g372499 * term.k);
        double cse_p11_t1_g372501 = (cse_p1_t2_invsqrt372482 * r23 * r43);
        double cse_p62_t1_g372509 = (cse_p1_t8_invr2372488 * v);
        double cse_p11_t1_g372505 = (cse_p11_t1_g372501 * cse_p62_t1_g372509);
        double g_x1 = (cse_p11_t1_g372500 * ((cse_p1_t7_invr372487 * ((cse_p11_t7_g372495 * dy43) + (dy23 * dz43))) + (-((cse_p11_t1_g372505 * dx13)))));
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = (cse_p11_t1_g372500 * ((cse_p1_t7_invr372487 * ((cse_p11_t3_g372491 * dz43) + (dx43 * dz23))) + (-((cse_p11_t1_g372505 * dy13)))));
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = (cse_p11_t1_g372500 * ((cse_p1_t7_invr372487 * ((cse_p11_t5_g372493 * dx43) + (dx23 * dy43))) + (-((cse_p11_t1_g372505 * dz13)))));
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double cse_p11_t1_g372503 = (cse_p1_t4_invsqrt372484 * r13 * r43);
        double cse_p11_t1_g372506 = (cse_p11_t1_g372503 * cse_p62_t1_g372509);
        double g_x2 = (cse_p11_t1_g372500 * ((cse_p1_t7_invr372487 * ((cse_p11_t4_g372492 * dz43) + (dy43 * dz13))) + (-((cse_p11_t1_g372506 * dx23)))));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (cse_p11_t1_g372500 * ((cse_p1_t7_invr372487 * ((cse_p11_t6_g372494 * dx43) + (dx13 * dz43))) + (-((cse_p11_t1_g372506 * dy23)))));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (cse_p11_t1_g372500 * ((cse_p1_t7_invr372487 * ((cse_p11_t2_g372490 * dy43) + (dx43 * dy13))) + (-((cse_p11_t1_g372506 * dz23)))));
        KernelGradientAcc(term.i3x2, 2, g_z2);
        double cse_p11_t1_g372504 = (cse_p11_t1_g372508 * cse_p1_t6_invsqrt372486);
        double g_x3 = (cse_p11_t1_g372500 * ((cse_p1_t7_invr372487 * ((dy43 * (cse_p11_t6_g372494 + dz23)) + (dz43 * (cse_p11_t5_g372493 + dy13)) + (-(cx)))) + (-((cse_p62_t1_g372509 * ((-((cse_p11_t1_g372501 * dx13))) + (-((cse_p11_t1_g372503 * dx23))) + (-((cse_p11_t1_g372504 * dx43)))))))));
        KernelGradientAcc(term.i3x3, 0, g_x3);
        double g_y3 = (cse_p11_t1_g372500 * ((cse_p1_t7_invr372487 * ((dx43 * (cse_p11_t7_g372495 + dz13)) + (dz43 * (cse_p11_t2_g372490 + dx23)) + (-(cy)))) + (-((cse_p62_t1_g372509 * ((-((cse_p11_t1_g372501 * dy13))) + (-((cse_p11_t1_g372503 * dy23))) + (-((cse_p11_t1_g372504 * dy43)))))))));
        KernelGradientAcc(term.i3x3, 1, g_y3);
        double g_z3 = (cse_p11_t1_g372500 * ((cse_p1_t7_invr372487 * ((dx43 * (cse_p11_t4_g372492 + dy23)) + (dy43 * (cse_p11_t3_g372491 + dx13)) + (-(cz)))) + (-((cse_p62_t1_g372509 * ((-((cse_p11_t1_g372501 * dz13))) + (-((cse_p11_t1_g372503 * dz23))) + (-((cse_p11_t1_g372504 * dz43)))))))));
        KernelGradientAcc(term.i3x3, 2, g_z3);
        double cse_p11_t1_g372507 = (cse_p11_t1_g372504 * cse_p62_t1_g372509);
        double g_x4 = (cse_p11_t1_g372500 * ((cse_p1_t7_invr372487 * cx) + (-((cse_p11_t1_g372507 * dx43)))));
        KernelGradientAcc(term.i3x4, 0, g_x4);
        double g_y4 = (cse_p11_t1_g372500 * ((cse_p1_t7_invr372487 * cy) + (-((cse_p11_t1_g372507 * dy43)))));
        KernelGradientAcc(term.i3x4, 1, g_y4);
        double g_z4 = (cse_p11_t1_g372500 * ((cse_p1_t7_invr372487 * cz) + (-((cse_p11_t1_g372507 * dz43)))));
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
      double cse_p11_t189_g372706 = (-(x3));
      double cse_p11_t190_g372707 = (-(y3));
      double cse_p11_t191_g372708 = (-(z3));
      double dx13 = (cse_p11_t189_g372706 + x1);
      double dy13 = (cse_p11_t190_g372707 + y1);
      double dz13 = (cse_p11_t191_g372708 + z1);
      double dx23 = (cse_p11_t189_g372706 + x2);
      double dy23 = (cse_p11_t190_g372707 + y2);
      double dz23 = (cse_p11_t191_g372708 + z2);
      double dx43 = (cse_p11_t189_g372706 + x4);
      double dy43 = (cse_p11_t190_g372707 + y4);
      double dz43 = (cse_p11_t191_g372708 + z4);
      double cse_p15_t1_g373119 = (dx13 * dx13);
      double cse_p16_t1_g373149 = (dy13 * dy13);
      double cse_p16_t1_g373151 = (dz13 * dz13);
      double r13_2 = (cse_p15_t1_g373119 + cse_p16_t1_g373149 + cse_p16_t1_g373151);
      double cse_p16_t1_g373148 = (dx23 * dx23);
      double cse_p16_t1_g373150 = (dy23 * dy23);
      double cse_p16_t1_g373152 = (dz23 * dz23);
      double r23_2 = (cse_p16_t1_g373148 + cse_p16_t1_g373150 + cse_p16_t1_g373152);
      double cse_p17_t1_g373171 = (dx43 * dx43);
      double cse_p18_t1_g373182 = (dy43 * dy43);
      double cse_p18_t1_g373183 = (dz43 * dz43);
      double r43_2 = (cse_p17_t1_g373171 + cse_p18_t1_g373182 + cse_p18_t1_g373183);
      double cse_p520_t1_sqrt373339 = sqrt(r13_2);
      double cse_p520_t2_invr373340 = (1.0 / (r13_2));
      double cse_p61_t3_invr2373279 = (cse_p520_t2_invr373340 * cse_p520_t2_invr373340);
      double cse_p1_t2_invsqrt372511 = (1.0 / (cse_p520_t1_sqrt373339));
      double r13 = cse_p520_t1_sqrt373339;
      double cse_p520_t3_sqrt373341 = sqrt(r23_2);
      double cse_p520_t4_invr373342 = (1.0 / (r23_2));
      double cse_p61_t6_invr2373282 = (cse_p520_t4_invr373342 * cse_p520_t4_invr373342);
      double cse_p1_t4_invsqrt372513 = (1.0 / (cse_p520_t3_sqrt373341));
      double r23 = cse_p520_t3_sqrt373341;
      double cse_p520_t5_sqrt373343 = sqrt(r43_2);
      double cse_p520_t6_invr373344 = (1.0 / (r43_2));
      double cse_p61_t9_invr2373285 = (cse_p520_t6_invr373344 * cse_p520_t6_invr373344);
      double cse_p1_t6_invsqrt372515 = (1.0 / (cse_p520_t5_sqrt373343));
      double r43 = cse_p520_t5_sqrt373343;
      double cse_p21_t1_g373215 = (dy13 * dz23);
      double cse_p21_t1_g373216 = (dy23 * dz13);
      double cx = (cse_p21_t1_g373215 + (-(cse_p21_t1_g373216)));
      double cse_p20_t1_g373209 = (dx13 * dz23);
      double cse_p21_t1_g373212 = (dx23 * dz13);
      double cy = (cse_p21_t1_g373212 + (-(cse_p20_t1_g373209)));
      double cse_p20_t1_g373207 = (dx13 * dy23);
      double cse_p20_t1_g373210 = (dx23 * dy13);
      double cz = (cse_p20_t1_g373207 + (-(cse_p20_t1_g373210)));
      double v = ((cx * dx43) + (cy * dy43) + (cz * dz43));
      double cse_p11_t1_g372718 = (r13 * r23);
      double denom = (cse_p11_t1_g372718 * r43);
      double cse_p61_t10_invr373286 = (1.0 / (denom));
      double cse_p61_t11_invr2373287 = (cse_p61_t10_invr373286 * cse_p61_t10_invr373286);
      double cse_p61_t12_invr3373288 = (cse_p61_t10_invr373286 * cse_p61_t11_invr2373287);
      double q = (cse_p61_t10_invr373286 * v);
      double cse_p11_t142_g372659 = (term.co + q);
      double chiral_e = (term.k * (cse_p11_t142_g372659 * cse_p11_t142_g372659 * cse_p11_t142_g372659));
      if ((chiral_e > 0.00000000000000000e+0)) {
        /* !ACTIVE */
        double cse_p63_t1_g373290 = (cse_p61_t11_invr2373287 * v);
        double cse_p16_t19_g373141 = (-(dx13));
        double cse_p16_t20_g373142 = (-(dx23));
        double cse_p16_t21_g373143 = (-(dy13));
        double cse_p16_t22_g373144 = (-(dy23));
        double cse_p16_t23_g373145 = (-(dz13));
        double cse_p16_t24_g373146 = (-(dz23));
        double cse_p15_t13_g373088 = (cse_p16_t20_g373142 + dx13);
        double cse_p15_t14_g373089 = (cse_p16_t19_g373141 + dx23);
        double cse_p15_t15_g373090 = (cse_p16_t22_g373144 + dy13);
        double cse_p15_t16_g373091 = (cse_p16_t21_g373143 + dy23);
        double cse_p15_t17_g373092 = (cse_p16_t24_g373146 + dz13);
        double cse_p15_t18_g373093 = (cse_p16_t23_g373145 + dz23);
        double cse_p14_t14_g373020 = (dx13 * dz43);
        double cse_p14_t15_g373021 = (dx23 * dy43);
        double cse_p14_t16_g373022 = (cse_p15_t16_g373091 * dx43);
        double cse_p14_t17_g373023 = (cse_p15_t17_g373092 * dx43);
        double cse_p14_t18_g373024 = (cse_p16_t22_g373144 * dx43);
        double cse_p14_t19_g373025 = (cse_p16_t23_g373145 * dx43);
        double cse_p14_t20_g373026 = (dx43 * dy13);
        double cse_p14_t21_g373027 = (dx43 * dz23);
        double cse_p14_t22_g373028 = (dy23 * dz43);
        double cse_p14_t23_g373029 = (cse_p15_t13_g373088 * dy43);
        double cse_p14_t24_g373030 = (cse_p15_t18_g373093 * dy43);
        double cse_p14_t25_g373031 = (cse_p16_t19_g373141 * dy43);
        double cse_p14_t26_g373032 = (cse_p16_t24_g373146 * dy43);
        double cse_p14_t27_g373033 = (dy43 * dz13);
        double cse_p14_t28_g373034 = (cse_p15_t14_g373089 * dz43);
        double cse_p14_t29_g373035 = (cse_p15_t15_g373090 * dz43);
        double cse_p14_t30_g373036 = (cse_p16_t20_g373142 * dz43);
        double cse_p14_t31_g373037 = (cse_p16_t21_g373143 * dz43);
        double cse_p14_t50_g373056 = (-(cx));
        double cse_p14_t51_g373057 = (-(cy));
        double cse_p14_t52_g373058 = (-(cz));
        double cse_p13_t43_g372941 = (cse_p14_t14_g373020 + cse_p14_t19_g373025);
        double cse_p13_t44_g372942 = (cse_p14_t15_g373021 + cse_p14_t18_g373024);
        double cse_p13_t45_g372943 = (cse_p14_t16_g373022 + cse_p14_t23_g373029 + cse_p14_t52_g373058);
        double cse_p13_t46_g372944 = (cse_p14_t17_g373023 + cse_p14_t28_g373034 + cse_p14_t51_g373057);
        double cse_p13_t47_g372945 = (cse_p14_t20_g373026 + cse_p14_t25_g373031);
        double cse_p13_t48_g372946 = (cse_p14_t21_g373027 + cse_p14_t30_g373036);
        double cse_p13_t49_g372947 = (cse_p14_t22_g373028 + cse_p14_t26_g373032);
        double cse_p13_t50_g372948 = (cse_p14_t24_g373030 + cse_p14_t29_g373035 + cse_p14_t50_g373056);
        double cse_p13_t51_g372949 = (cse_p14_t27_g373033 + cse_p14_t31_g373037);
        double cse_p21_t1_g373221 = (2.00000000000000000e+0 * cse_p61_t11_invr2373287);
        double cse_p12_t4_g372741 = (cse_p21_t1_g373221 * cse_p61_t3_invr2373279);
        double cse_p12_t5_g372742 = (cse_p21_t1_g373221 * cse_p61_t6_invr2373282);
        double cse_p12_t6_g372743 = (cse_p21_t1_g373221 * cse_p61_t9_invr2373285);
        double cse_p12_t53_g372790 = (cse_p13_t43_g372941 * cse_p61_t10_invr373286);
        double cse_p12_t54_g372791 = (cse_p13_t44_g372942 * cse_p61_t10_invr373286);
        double cse_p12_t55_g372792 = (cse_p13_t45_g372943 * cse_p61_t10_invr373286);
        double cse_p12_t56_g372793 = (cse_p13_t46_g372944 * cse_p61_t10_invr373286);
        double cse_p12_t57_g372794 = (cse_p13_t47_g372945 * cse_p61_t10_invr373286);
        double cse_p12_t58_g372795 = (cse_p13_t48_g372946 * cse_p61_t10_invr373286);
        double cse_p12_t59_g372796 = (cse_p13_t49_g372947 * cse_p61_t10_invr373286);
        double cse_p12_t60_g372797 = (cse_p13_t50_g372948 * cse_p61_t10_invr373286);
        double cse_p12_t61_g372798 = (cse_p13_t51_g372949 * cse_p61_t10_invr373286);
        double cse_p12_t62_g372799 = (cse_p61_t10_invr373286 * cx);
        double cse_p12_t63_g372800 = (cse_p61_t10_invr373286 * cy);
        double cse_p12_t64_g372801 = (cse_p61_t10_invr373286 * cz);
        double cse_p13_t4_g372902 = (cse_p520_t1_sqrt373339 * cse_p520_t2_invr373340);
        double cse_p13_t5_g372903 = (cse_p520_t3_sqrt373341 * cse_p520_t4_invr373342);
        double cse_p13_t6_g372904 = (cse_p520_t5_sqrt373343 * cse_p520_t6_invr373344);
        double cse_p11_t1_g372719 = (r13 * r43);
        double cse_p22_t1_g373224 = (cse_p11_t1_g372719 * cse_p520_t4_invr373342);
        double cse_p12_t28_g372765 = (cse_p22_t1_g373224 * dx23);
        double cse_p12_t29_g372766 = (cse_p22_t1_g373224 * dy23);
        double cse_p12_t30_g372767 = (cse_p22_t1_g373224 * dz23);
        double cse_p11_t1_g372720 = (r23 * r43);
        double cse_p22_t1_g373225 = (cse_p11_t1_g372720 * cse_p520_t2_invr373340);
        double cse_p12_t38_g372775 = (cse_p22_t1_g373225 * dx13);
        double cse_p12_t39_g372776 = (cse_p22_t1_g373225 * dy13);
        double cse_p12_t40_g372777 = (cse_p22_t1_g373225 * dz13);
        double cse_p18_t1_g373185 = (cse_p11_t1_g372718 * cse_p13_t6_g372904);
        double cse_p25_t1_g373257 = (cse_p18_t1_g373185 * dx43);
        double cse_p27_t1_g373276 = (cse_p18_t1_g373185 * dy43);
        double cse_p17_t3_g373159 = (cse_p18_t1_g373185 * dz43);
        double cse_p16_t10_g373132 = (-(cse_p25_t1_g373257));
        double cse_p16_t11_g373133 = (-(cse_p27_t1_g373276));
        double cse_p16_t12_g373134 = (-(cse_p17_t3_g373159));
        double cse_p21_t1_g373220 = (2.00000000000000000e+0 * cse_p61_t12_invr3373288);
        double cse_p26_t1_g373270 = (cse_p18_t1_g373185 * cse_p520_t6_invr373344);
        double cse_p12_t7_g372744 = (cse_p21_t1_g373220 * cse_p26_t1_g373270);
        double cse_p17_t1_g373168 = (cse_p520_t6_invr373344 * dy43);
        double cse_p12_t10_g372747 = (cse_p17_t1_g373168 * cse_p25_t1_g373257);
        double cse_p15_t1_g373113 = (cse_p520_t6_invr373344 * dz43);
        double cse_p12_t11_g372748 = (cse_p15_t1_g373113 * cse_p25_t1_g373257);
        double cse_p20_t1_g373204 = (cse_p18_t1_g373185 * cse_p63_t1_g373290);
        double cse_p13_t7_g372905 = (cse_p20_t1_g373204 * dx43);
        double cse_p13_t8_g372906 = (cse_p20_t1_g373204 * dy43);
        double cse_p13_t9_g372907 = (cse_p20_t1_g373204 * dz43);
        double cse_p12_t14_g372751 = (cse_p15_t1_g373113 * cse_p27_t1_g373276);
        double cse_p14_t1_g373071 = (cse_p13_t6_g372904 * dx43);
        double cse_p12_t41_g372778 = (cse_p14_t1_g373071 * r13);
        double cse_p12_t42_g372779 = (cse_p14_t1_g373071 * r23);
        double cse_p14_t1_g373072 = (cse_p13_t6_g372904 * dy43);
        double cse_p12_t43_g372780 = (cse_p14_t1_g373072 * r13);
        double cse_p12_t44_g372781 = (cse_p14_t1_g373072 * r23);
        double cse_p14_t1_g373073 = (cse_p13_t6_g372904 * dz43);
        double cse_p12_t45_g372782 = (cse_p14_t1_g373073 * r13);
        double cse_p12_t46_g372783 = (cse_p14_t1_g373073 * r23);
        double cse_p12_t107_g372844 = (-(cse_p13_t7_g372905));
        double cse_p12_t108_g372845 = (-(cse_p13_t8_g372906));
        double cse_p12_t109_g372846 = (-(cse_p13_t9_g372907));
        double cse_p17_t1_g373174 = (cse_p11_t1_g372720 * cse_p13_t4_g372902);
        double cse_p17_t7_g373163 = (cse_p17_t1_g373174 * dx13);
        double cse_p17_t8_g373164 = (cse_p17_t1_g373174 * dy13);
        double cse_p17_t9_g373165 = (cse_p17_t1_g373174 * dz13);
        double cse_p16_t16_g373138 = (-(cse_p17_t7_g373163));
        double cse_p16_t17_g373139 = (-(cse_p17_t8_g373164));
        double cse_p16_t18_g373140 = (-(cse_p17_t9_g373165));
        double cse_p20_t1_g373203 = (cse_p17_t1_g373174 * cse_p63_t1_g373290);
        double cse_p13_t13_g372911 = (cse_p20_t1_g373203 * dx13);
        double cse_p13_t14_g372912 = (cse_p20_t1_g373203 * dy13);
        double cse_p13_t15_g372913 = (cse_p20_t1_g373203 * dz13);
        double cse_p26_t1_g373268 = (cse_p17_t1_g373174 * cse_p520_t2_invr373340);
        double cse_p12_t9_g372746 = (cse_p21_t1_g373220 * cse_p26_t1_g373268);
        double cse_p14_t1_g373075 = (cse_p13_t4_g372902 * r43);
        double cse_p12_t47_g372784 = (cse_p14_t1_g373075 * dx13);
        double cse_p12_t48_g372785 = (cse_p14_t1_g373075 * dy13);
        double cse_p12_t49_g372786 = (cse_p14_t1_g373075 * dz13);
        double cse_p12_t113_g372850 = (-(cse_p13_t13_g372911));
        double cse_p12_t114_g372851 = (-(cse_p13_t14_g372912));
        double cse_p12_t115_g372852 = (-(cse_p13_t15_g372913));
        double cse_p17_t1_g373173 = (cse_p11_t1_g372719 * cse_p13_t5_g372903);
        double cse_p17_t4_g373160 = (cse_p17_t1_g373173 * dx23);
        double cse_p17_t5_g373161 = (cse_p17_t1_g373173 * dy23);
        double cse_p17_t6_g373162 = (cse_p17_t1_g373173 * dz23);
        double cse_p16_t13_g373135 = (-(cse_p17_t4_g373160));
        double cse_p16_t14_g373136 = (-(cse_p17_t5_g373161));
        double cse_p16_t15_g373137 = (-(cse_p17_t6_g373162));
        double cse_p15_t10_g373085 = (cse_p16_t10_g373132 + cse_p16_t13_g373135 + cse_p16_t16_g373138);
        double cse_p15_t11_g373086 = (cse_p16_t11_g373133 + cse_p16_t14_g373136 + cse_p16_t17_g373139);
        double cse_p15_t12_g373087 = (cse_p16_t12_g373134 + cse_p16_t15_g373137 + cse_p16_t18_g373140);
        double cse_p14_t1_g373007 = (cse_p15_t10_g373085 * cse_p61_t11_invr2373287);
        double cse_p14_t2_g373008 = (cse_p15_t11_g373086 * cse_p61_t11_invr2373287);
        double cse_p14_t3_g373009 = (cse_p15_t12_g373087 * cse_p61_t11_invr2373287);
        double cse_p20_t1_g373202 = (cse_p17_t1_g373173 * cse_p63_t1_g373290);
        double cse_p13_t10_g372908 = (cse_p20_t1_g373202 * dx23);
        double cse_p13_t11_g372909 = (cse_p20_t1_g373202 * dy23);
        double cse_p13_t12_g372910 = (cse_p20_t1_g373202 * dz23);
        double cse_p13_t61_g372959 = (-(cse_p14_t1_g373007));
        double cse_p13_t62_g372960 = (-(cse_p14_t2_g373008));
        double cse_p13_t63_g372961 = (-(cse_p14_t3_g373009));
        double cse_p26_t1_g373269 = (cse_p17_t1_g373173 * cse_p520_t4_invr373342);
        double cse_p12_t8_g372745 = (cse_p21_t1_g373220 * cse_p26_t1_g373269);
        double cse_p24_t1_g373250 = (cse_p13_t5_g372903 * r43);
        double cse_p12_t50_g372787 = (cse_p24_t1_g373250 * dx23);
        double cse_p12_t51_g372788 = (cse_p24_t1_g373250 * dy23);
        double cse_p12_t52_g372789 = (cse_p24_t1_g373250 * dz23);
        double cse_p12_t83_g372820 = (cse_p13_t61_g372959 * v);
        double cse_p12_t84_g372821 = (cse_p13_t62_g372960 * v);
        double cse_p12_t85_g372822 = (cse_p13_t63_g372961 * v);
        double cse_p12_t110_g372847 = (-(cse_p13_t10_g372908));
        double cse_p12_t111_g372848 = (-(cse_p13_t11_g372909));
        double cse_p12_t112_g372849 = (-(cse_p13_t12_g372910));
        double cse_p19_t1_g373188 = (cse_p11_t1_g372720 * cse_p520_t1_sqrt373339);
        double cse_p13_t1_g373001 = (cse_p19_t1_g373188 * cse_p61_t3_invr2373279);
        double cse_p14_t1_g373068 = (-2.00000000000000000e+0 * cse_p13_t1_g373001);
        double cse_p17_t1_g373170 = (cse_p14_t1_g373068 * cse_p63_t1_g373290);
        double cse_p20_t1_g373206 = (dx13 * dy13);
        double cse_p11_t4_g372521 = (cse_p17_t1_g373170 * cse_p20_t1_g373206);
        double cse_p20_t1_g373208 = (dx13 * dz13);
        double cse_p11_t5_g372522 = (cse_p17_t1_g373170 * cse_p20_t1_g373208);
        double cse_p19_t1_g373186 = (cse_p11_t1_g372719 * cse_p520_t3_sqrt373341);
        double cse_p13_t1_g373000 = (cse_p19_t1_g373186 * cse_p61_t6_invr2373282);
        double cse_p14_t1_g373067 = (-2.00000000000000000e+0 * cse_p13_t1_g373000);
        double cse_p17_t1_g373169 = (cse_p14_t1_g373067 * cse_p63_t1_g373290);
        double cse_p20_t1_g373211 = (dx23 * dy23);
        double cse_p11_t6_g372523 = (cse_p17_t1_g373169 * cse_p20_t1_g373211);
        double cse_p21_t1_g373213 = (dx23 * dz23);
        double cse_p11_t7_g372524 = (cse_p17_t1_g373169 * cse_p21_t1_g373213);
        double cse_p21_t1_g373214 = (dy13 * dz13);
        double cse_p11_t8_g372525 = (cse_p17_t1_g373170 * cse_p21_t1_g373214);
        double cse_p21_t1_g373217 = (dy23 * dz23);
        double cse_p11_t9_g372526 = (cse_p17_t1_g373169 * cse_p21_t1_g373217);
        double cse_p13_t1_g372999 = (cse_p11_t1_g372718 * cse_p520_t5_sqrt373343 * cse_p61_t9_invr2373285);
        double cse_p16_t1_g373155 = (2.00000000000000000e+0 * cse_p13_t1_g372999);
        double cse_p27_t1_g373275 = (cse_p16_t1_g373155 * dx43);
        double cse_p11_t13_g372530 = (cse_p27_t1_g373275 * dy43);
        double cse_p11_t14_g372531 = (cse_p27_t1_g373275 * dz43);
        double cse_p25_t1_g373260 = (dy43 * dz43);
        double cse_p11_t15_g372532 = (cse_p16_t1_g373155 * cse_p25_t1_g373260);
        double cse_p15_t1_g373118 = (cse_p13_t5_g372903 * cse_p14_t1_g373075);
        double cse_p25_t1_g373258 = (dx13 * dx23);
        double cse_p11_t31_g372548 = (cse_p15_t1_g373118 * cse_p25_t1_g373258);
        double cse_p16_t1_g373147 = (cse_p15_t1_g373118 * cse_p63_t1_g373290);
        double cse_p11_t32_g372549 = (cse_p16_t1_g373147 * cse_p25_t1_g373258);
        double cse_p19_t1_g373194 = (cse_p13_t4_g372902 * r23);
        double cse_p20_t1_g373196 = (cse_p19_t1_g373194 * dx13);
        double cse_p20_t2_g373197 = (cse_p19_t1_g373194 * dy13);
        double cse_p20_t3_g373198 = (cse_p19_t1_g373194 * dz13);
        double cse_p19_t1_g373195 = (cse_p13_t5_g372903 * r13);
        double cse_p20_t4_g373199 = (cse_p19_t1_g373195 * dx23);
        double cse_p20_t5_g373200 = (cse_p19_t1_g373195 * dy23);
        double cse_p20_t6_g373201 = (cse_p19_t1_g373195 * dz23);
        double cse_p25_t1_g373259 = (dy13 * dy23);
        double cse_p11_t57_g372574 = (cse_p15_t1_g373118 * cse_p25_t1_g373259);
        double cse_p11_t58_g372575 = (cse_p16_t1_g373147 * cse_p25_t1_g373259);
        double cse_p26_t1_g373262 = (dz13 * dz23);
        double cse_p11_t80_g372597 = (cse_p15_t1_g373118 * cse_p26_t1_g373262);
        double cse_p11_t81_g372598 = (cse_p16_t1_g373147 * cse_p26_t1_g373262);
        double cse_p11_t109_g372626 = (cse_p12_t111_g372848 + cse_p12_t53_g372790);
        double cse_p11_t110_g372627 = (cse_p12_t115_g372852 + cse_p12_t54_g372791);
        double cse_p11_t111_g372628 = (cse_p12_t55_g372792 + cse_p12_t85_g372822);
        double cse_p11_t112_g372629 = (cse_p12_t56_g372793 + cse_p12_t84_g372821);
        double cse_p11_t113_g372630 = (cse_p12_t112_g372849 + cse_p12_t57_g372794);
        double cse_p11_t114_g372631 = (cse_p12_t114_g372851 + cse_p12_t58_g372795);
        double cse_p11_t115_g372632 = (cse_p12_t113_g372850 + cse_p12_t59_g372796);
        double cse_p11_t116_g372633 = (cse_p12_t60_g372797 + cse_p12_t83_g372820);
        double cse_p11_t117_g372634 = (cse_p12_t110_g372847 + cse_p12_t61_g372798);
        double cse_p11_t118_g372635 = (cse_p12_t28_g372765 + cse_p12_t41_g372778 + cse_p12_t47_g372784);
        double cse_p11_t119_g372636 = (cse_p12_t29_g372766 + cse_p12_t43_g372780 + cse_p12_t48_g372785);
        double cse_p11_t120_g372637 = (cse_p12_t30_g372767 + cse_p12_t45_g372782 + cse_p12_t49_g372786);
        double cse_p11_t124_g372641 = (cse_p12_t4_g372741 + cse_p12_t9_g372746);
        double cse_p11_t125_g372642 = (cse_p12_t5_g372742 + cse_p12_t8_g372745);
        double cse_p11_t126_g372643 = (cse_p12_t6_g372743 + cse_p12_t7_g372744);
        double cse_p11_t127_g372644 = (cse_p12_t107_g372844 + cse_p12_t62_g372799);
        double cse_p11_t128_g372645 = (cse_p12_t108_g372845 + cse_p12_t63_g372800);
        double cse_p11_t129_g372646 = (cse_p12_t109_g372846 + cse_p12_t64_g372801);
        double cse_p11_t131_g372648 = (cse_p12_t38_g372775 + cse_p12_t42_g372779 + cse_p12_t50_g372787);
        double cse_p11_t137_g372654 = (cse_p12_t39_g372776 + cse_p12_t44_g372781 + cse_p12_t51_g372788);
        double cse_p11_t141_g372658 = (cse_p12_t40_g372777 + cse_p12_t46_g372783 + cse_p12_t52_g372789);
        double cse_p11_t156_g372673 = (-(cse_p12_t10_g372747));
        double cse_p11_t157_g372674 = (-(cse_p12_t11_g372748));
        double cse_p11_t164_g372681 = (-(cse_p12_t14_g372751));
        double cse_p11_t173_g372690 = (-(cse_p18_t1_g373185));
        double cse_p11_t174_g372691 = (-(cse_p20_t1_g373204));
        double cse_p11_t175_g372692 = (-(cse_p20_t1_g373202));
        double cse_p11_t176_g372693 = (-(cse_p20_t1_g373203));
        double cse_p11_t182_g372699 = (-(dx43));
        double cse_p11_t185_g372702 = (-(dy43));
        double cse_p11_t188_g372705 = (-(dz43));
        double cse_p11_t192_g372709 = ((cse_p11_t142_g372659) * (cse_p11_t142_g372659));
        double energy = chiral_e;
        energy_added += energy;
        *energy_accumulate += energy;
        double cse_p11_t1_g372727 = (3.00000000000000000e+0 * cse_p11_t192_g372709 * term.k);
        double cse_p19_t1_g373189 = (cse_p11_t1_g372720 * cse_p1_t2_invsqrt372511);
        double cse_p15_t1_g373121 = (cse_p19_t1_g373189 * cse_p63_t1_g373290);
        double g_x1 = (cse_p11_t1_g372727 * (cse_p12_t59_g372796 + (-((cse_p15_t1_g373121 * dx13)))));
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = (cse_p11_t1_g372727 * (cse_p12_t58_g372795 + (-((cse_p15_t1_g373121 * dy13)))));
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = (cse_p11_t1_g372727 * (cse_p12_t54_g372791 + (-((cse_p15_t1_g373121 * dz13)))));
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double cse_p19_t1_g373187 = (cse_p11_t1_g372719 * cse_p1_t4_invsqrt372513);
        double cse_p15_t1_g373120 = (cse_p19_t1_g373187 * cse_p63_t1_g373290);
        double g_x2 = (cse_p11_t1_g372727 * (cse_p12_t61_g372798 + (-((cse_p15_t1_g373120 * dx23)))));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (cse_p11_t1_g372727 * (cse_p12_t53_g372790 + (-((cse_p15_t1_g373120 * dy23)))));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (cse_p11_t1_g372727 * (cse_p12_t57_g372794 + (-((cse_p15_t1_g373120 * dz23)))));
        KernelGradientAcc(term.i3x2, 2, g_z2);
        double cse_p15_t1_g373117 = (cse_p11_t1_g372718 * dx43);
        double cse_p27_t1_g373272 = (cse_p15_t1_g373117 * cse_p1_t6_invsqrt372515);
        double g_x3 = (cse_p11_t1_g372727 * (cse_p12_t60_g372797 + (-((cse_p63_t1_g373290 * ((-((cse_p19_t1_g373187 * dx23))) + (-((cse_p19_t1_g373189 * dx13))) + (-(cse_p27_t1_g373272))))))));
        KernelGradientAcc(term.i3x3, 0, g_x3);
        double cse_p22_t1_g373223 = (cse_p11_t1_g372718 * cse_p1_t6_invsqrt372515);
        double g_y3 = (cse_p11_t1_g372727 * (cse_p12_t56_g372793 + (-((cse_p63_t1_g373290 * ((-((cse_p19_t1_g373187 * dy23))) + (-((cse_p19_t1_g373189 * dy13))) + (-((cse_p22_t1_g373223 * dy43)))))))));
        KernelGradientAcc(term.i3x3, 1, g_y3);
        double g_z3 = (cse_p11_t1_g372727 * (cse_p12_t55_g372792 + (-((cse_p63_t1_g373290 * ((-((cse_p19_t1_g373187 * dz23))) + (-((cse_p19_t1_g373189 * dz13))) + (-((cse_p22_t1_g373223 * dz43)))))))));
        KernelGradientAcc(term.i3x3, 2, g_z3);
        double g_x4 = (cse_p11_t1_g372727 * (cse_p12_t62_g372799 + (-((cse_p27_t1_g373272 * cse_p63_t1_g373290)))));
        KernelGradientAcc(term.i3x4, 0, g_x4);
        double cse_p15_t1_g373122 = (cse_p22_t1_g373223 * cse_p63_t1_g373290);
        double g_y4 = (cse_p11_t1_g372727 * (cse_p12_t63_g372800 + (-((cse_p15_t1_g373122 * dy43)))));
        KernelGradientAcc(term.i3x4, 1, g_y4);
        double g_z4 = (cse_p11_t1_g372727 * (cse_p12_t64_g372801 + (-((cse_p15_t1_g373122 * dz43)))));
        KernelGradientAcc(term.i3x4, 2, g_z4);
        double cse_p11_t1_g372730 = (cse_p17_t1_g373174 * cse_p61_t11_invr2373287);
        double cse_p12_t1_g372879 = (6.00000000000000000e+0 * cse_p11_t142_g372659 * term.k);
        double cse_p12_t1_g372880 = (2.00000000000000000e+0 * v);
        double cse_p12_t1_g372883 = (cse_p12_t1_g372880 * cse_p61_t12_invr3373288);
        double cse_p22_t1_g373226 = (cse_p12_t1_g372880 * cse_p61_t11_invr2373287);
        double cse_p12_t1_g372895 = (cse_p13_t1_g373001 * cse_p22_t1_g373226);
        double cse_p13_t1_g372991 = (cse_p11_t1_g372730 * dx13);
        double cse_p13_t1_g372992 = (cse_p12_t1_g372883 * r43);
        double cse_p16_t1_g373156 = (cse_p11_t115_g372632 * cse_p12_t1_g372879);
        double cse_p17_t1_g373166 = (cse_p520_t2_invr373340 * v);
        double cse_p19_t1_g373190 = (cse_p13_t1_g372992 * cse_p26_t1_g373268 * cse_p520_t1_sqrt373339 * r23);
        double cse_p24_t1_g373247 = (cse_p13_t1_g372991 * cse_p17_t1_g373166);
        double h_x1_x1 = ((cse_p11_t115_g372632 * cse_p16_t1_g373156) + (cse_p11_t1_g372727 * (cse_p11_t176_g372693 + (-2.00000000000000000e+0 * cse_p13_t1_g372991 * cse_p13_t49_g372947) + (cse_p12_t1_g372895 * cse_p15_t1_g373119) + (cse_p15_t1_g373119 * cse_p19_t1_g373190) + (-((cse_p24_t1_g373247 * dx13))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double cse_p13_t1_g372995 = (cse_p11_t1_g372730 * dy13);
        double cse_p15_t1_g373116 = (cse_p11_t114_g372631 * cse_p12_t1_g372879);
        double cse_p21_t1_g373218 = (cse_p11_t124_g372641 * cse_p19_t1_g373188 * v);
        double h_x1_y1 = ((cse_p11_t115_g372632 * cse_p15_t1_g373116) + (cse_p11_t1_g372727 * ((cse_p20_t1_g373206 * cse_p21_t1_g373218) + (-((cse_p13_t1_g372991 * cse_p13_t48_g372946))) + (-((cse_p13_t1_g372995 * cse_p13_t49_g372947))) + (-((cse_p24_t1_g373247 * dy13))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double cse_p14_t1_g373066 = (cse_p11_t1_g372730 * dz13);
        double cse_p14_t1_g373070 = (cse_p11_t110_g372627 * cse_p12_t1_g372879);
        double h_x1_z1 = ((cse_p11_t115_g372632 * cse_p14_t1_g373070) + (cse_p11_t1_g372727 * ((cse_p20_t1_g373208 * cse_p21_t1_g373218) + (-((cse_p13_t1_g372991 * cse_p13_t44_g372942))) + (-((cse_p13_t49_g372947 * cse_p14_t1_g373066))) + (-((cse_p24_t1_g373247 * dz13))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double cse_p11_t1_g372729 = (cse_p17_t1_g373173 * cse_p61_t11_invr2373287);
        double cse_p13_t1_g372990 = (cse_p11_t1_g372729 * dx23);
        double cse_p13_t1_g372998 = (cse_p11_t1_g372718 * cse_p13_t1_g372992 * cse_p15_t1_g373118);
        double h_x1_x2 = ((cse_p11_t117_g372634 * cse_p16_t1_g373156) + (cse_p11_t1_g372727 * ((cse_p13_t1_g372998 * cse_p25_t1_g373258) + (-((cse_p13_t1_g372990 * cse_p13_t49_g372947))) + (-((cse_p13_t1_g372991 * cse_p13_t51_g372949))) + (-(cse_p11_t32_g372549)))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double cse_p13_t1_g372994 = (cse_p11_t1_g372729 * dy23);
        double cse_p14_t1_g373069 = (cse_p11_t109_g372626 * cse_p12_t1_g372879);
        double h_x1_y2 = ((cse_p11_t115_g372632 * cse_p14_t1_g373069) + (cse_p11_t1_g372727 * ((cse_p13_t1_g372998 * cse_p20_t1_g373207) + (cse_p61_t10_invr373286 * dz43) + (-((cse_p13_t1_g372991 * cse_p13_t43_g372941))) + (-((cse_p13_t1_g372994 * cse_p13_t49_g372947))) + (-((cse_p16_t1_g373147 * cse_p20_t1_g373207))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double cse_p13_t1_g372997 = (cse_p11_t1_g372729 * dz23);
        double cse_p15_t1_g373115 = (cse_p11_t113_g372630 * cse_p12_t1_g372879);
        double h_x1_z2 = ((cse_p11_t115_g372632 * cse_p15_t1_g373115) + (cse_p11_t1_g372727 * ((cse_p11_t185_g372702 * cse_p61_t10_invr373286) + (cse_p13_t1_g372998 * cse_p20_t1_g373209) + (-((cse_p13_t1_g372991 * cse_p13_t47_g372945))) + (-((cse_p13_t1_g372997 * cse_p13_t49_g372947))) + (-((cse_p16_t1_g373147 * cse_p20_t1_g373209))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double cse_p12_t1_g372892 = (cse_p12_t1_g372883 * cse_p17_t1_g373174);
        double cse_p19_t1_g373192 = (cse_p19_t1_g373194 * cse_p63_t1_g373290);
        double cse_p23_t1_g373236 = (cse_p12_t1_g372892 * cse_p15_t10_g373085);
        double cse_p24_t1_g373251 = (cse_p520_t2_invr373340 * r43);
        double h_x1_x3 = ((cse_p11_t116_g372633 * cse_p16_t1_g373156) + (cse_p11_t1_g372727 * (cse_p11_t32_g372549 + (cse_p13_t49_g372947 * cse_p13_t61_g372959) + (cse_p15_t1_g373119 * cse_p17_t1_g373170) + (cse_p19_t1_g373192 * (r43 + (cse_p14_t1_g373071 * dx13) + (cse_p15_t1_g373119 * cse_p24_t1_g373251))) + (cse_p23_t1_g373236 * dx13) + (-((cse_p13_t1_g372991 * cse_p13_t50_g372948))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 0, h_x1_x3);
        double cse_p24_t1_g373248 = (cse_p13_t4_g372902 * cse_p63_t1_g373290);
        double cse_p13_t1_g373005 = (cse_p11_t137_g372654 * cse_p24_t1_g373248);
        double cse_p15_t1_g373112 = (cse_p11_t112_g372629 * cse_p12_t1_g372879);
        double cse_p23_t1_g373237 = (cse_p12_t1_g372892 * cse_p15_t11_g373086);
        double h_x1_y3 = ((cse_p11_t115_g372632 * cse_p15_t1_g373112) + (cse_p11_t1_g372727 * (cse_p11_t4_g372521 + (cse_p13_t1_g373005 * dx13) + (cse_p13_t49_g372947 * cse_p13_t62_g372960) + (cse_p23_t1_g373237 * dx13) + (cse_p61_t10_invr373286 * (cse_p11_t188_g372705 + dz23)) + (-((cse_p13_t1_g372991 * cse_p13_t46_g372944))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 1, h_x1_y3);
        double cse_p13_t1_g373006 = (cse_p11_t141_g372658 * cse_p24_t1_g373248);
        double cse_p14_t1_g373074 = (cse_p11_t111_g372628 * cse_p12_t1_g372879);
        double cse_p23_t1_g373238 = (cse_p12_t1_g372892 * cse_p15_t12_g373087);
        double h_x1_z3 = ((cse_p11_t115_g372632 * cse_p14_t1_g373074) + (cse_p11_t1_g372727 * (cse_p11_t5_g372522 + (cse_p13_t1_g373006 * dx13) + (cse_p13_t49_g372947 * cse_p13_t63_g372961) + (cse_p23_t1_g373238 * dx13) + (cse_p61_t10_invr373286 * (cse_p16_t22_g373144 + dy43)) + (-((cse_p13_t1_g372991 * cse_p13_t45_g372943))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 2, h_x1_z3);
        double cse_p11_t1_g372728 = (cse_p18_t1_g373185 * cse_p61_t11_invr2373287);
        double cse_p11_t1_g372731 = (cse_p18_t1_g373185 * cse_p61_t12_invr3373288);
        double cse_p17_t1_g373175 = (cse_p11_t1_g372731 * cse_p12_t1_g372880);
        double cse_p12_t1_g372890 = (cse_p17_t1_g373174 * cse_p17_t1_g373175);
        double cse_p12_t1_g372893 = (cse_p13_t6_g372904 * cse_p19_t1_g373192);
        double cse_p13_t1_g372989 = (cse_p11_t1_g372728 * dx43);
        double cse_p22_t1_g373230 = (cse_p12_t1_g372890 * dx13);
        double cse_p23_t1_g373239 = (cse_p12_t1_g372893 * dx13);
        double h_x1_x4 = ((cse_p11_t127_g372644 * cse_p16_t1_g373156) + (cse_p11_t1_g372727 * ((cse_p22_t1_g373230 * dx43) + (-((cse_p13_t1_g372989 * cse_p13_t49_g372947))) + (-((cse_p13_t1_g372991 * cx))) + (-((cse_p23_t1_g373239 * dx43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x4, 0, h_x1_x4);
        double cse_p13_t1_g372993 = (cse_p11_t1_g372728 * dy43);
        double h_x1_y4 = ((cse_p11_t128_g372645 * cse_p16_t1_g373156) + (cse_p11_t1_g372727 * ((cse_p16_t24_g373146 * cse_p61_t10_invr373286) + (cse_p22_t1_g373230 * dy43) + (-((cse_p13_t1_g372991 * cy))) + (-((cse_p13_t1_g372993 * cse_p13_t49_g372947))) + (-((cse_p23_t1_g373239 * dy43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x4, 1, h_x1_y4);
        double cse_p13_t1_g372996 = (cse_p11_t1_g372728 * dz43);
        double h_x1_z4 = ((cse_p11_t129_g372646 * cse_p16_t1_g373156) + (cse_p11_t1_g372727 * ((cse_p22_t1_g373230 * dz43) + (cse_p61_t10_invr373286 * dy23) + (-((cse_p13_t1_g372991 * cz))) + (-((cse_p13_t1_g372996 * cse_p13_t49_g372947))) + (-((cse_p23_t1_g373239 * dz43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x4, 2, h_x1_z4);
        double cse_p26_t1_g373267 = (cse_p13_t1_g372995 * cse_p17_t1_g373166);
        double h_y1_y1 = ((cse_p11_t114_g372631 * cse_p15_t1_g373116) + (cse_p11_t1_g372727 * (cse_p11_t176_g372693 + (-2.00000000000000000e+0 * cse_p13_t1_g372995 * cse_p13_t48_g372946) + (cse_p12_t1_g372895 * cse_p16_t1_g373149) + (cse_p16_t1_g373149 * cse_p19_t1_g373190) + (-((cse_p26_t1_g373267 * dy13))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = ((cse_p11_t114_g372631 * cse_p14_t1_g373070) + (cse_p11_t1_g372727 * ((cse_p21_t1_g373214 * cse_p21_t1_g373218) + (-((cse_p13_t1_g372995 * cse_p13_t44_g372942))) + (-((cse_p13_t48_g372946 * cse_p14_t1_g373066))) + (-((cse_p26_t1_g373267 * dz13))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = ((cse_p11_t117_g372634 * cse_p15_t1_g373116) + (cse_p11_t1_g372727 * ((cse_p11_t188_g372705 * cse_p61_t10_invr373286) + (cse_p13_t1_g372998 * cse_p20_t1_g373210) + (-((cse_p13_t1_g372990 * cse_p13_t48_g372946))) + (-((cse_p13_t1_g372995 * cse_p13_t51_g372949))) + (-((cse_p16_t1_g373147 * cse_p20_t1_g373210))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = ((cse_p11_t114_g372631 * cse_p14_t1_g373069) + (cse_p11_t1_g372727 * ((cse_p13_t1_g372998 * cse_p25_t1_g373259) + (-((cse_p13_t1_g372994 * cse_p13_t48_g372946))) + (-((cse_p13_t1_g372995 * cse_p13_t43_g372941))) + (-(cse_p11_t58_g372575)))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = ((cse_p11_t114_g372631 * cse_p15_t1_g373115) + (cse_p11_t1_g372727 * ((cse_p13_t1_g372998 * cse_p21_t1_g373215) + (cse_p61_t10_invr373286 * dx43) + (-((cse_p13_t1_g372995 * cse_p13_t47_g372945))) + (-((cse_p13_t1_g372997 * cse_p13_t48_g372946))) + (-((cse_p16_t1_g373147 * cse_p21_t1_g373215))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double cse_p13_t1_g373004 = (cse_p11_t131_g372648 * cse_p24_t1_g373248);
        double h_y1_x3 = ((cse_p11_t116_g372633 * cse_p15_t1_g373116) + (cse_p11_t1_g372727 * (cse_p11_t4_g372521 + (cse_p13_t1_g373004 * dy13) + (cse_p13_t48_g372946 * cse_p13_t61_g372959) + (cse_p23_t1_g373236 * dy13) + (cse_p61_t10_invr373286 * (cse_p16_t24_g373146 + dz43)) + (-((cse_p13_t1_g372995 * cse_p13_t50_g372948))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 0, h_y1_x3);
        double h_y1_y3 = ((cse_p11_t114_g372631 * cse_p15_t1_g373112) + (cse_p11_t1_g372727 * (cse_p11_t58_g372575 + (cse_p13_t48_g372946 * cse_p13_t62_g372960) + (cse_p16_t1_g373149 * cse_p17_t1_g373170) + (cse_p19_t1_g373192 * (r43 + (cse_p14_t1_g373072 * dy13) + (cse_p16_t1_g373149 * cse_p24_t1_g373251))) + (cse_p23_t1_g373237 * dy13) + (-((cse_p13_t1_g372995 * cse_p13_t46_g372944))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 1, h_y1_y3);
        double h_y1_z3 = ((cse_p11_t114_g372631 * cse_p14_t1_g373074) + (cse_p11_t1_g372727 * (cse_p11_t8_g372525 + (cse_p13_t1_g373006 * dy13) + (cse_p13_t48_g372946 * cse_p13_t63_g372961) + (cse_p23_t1_g373238 * dy13) + (cse_p61_t10_invr373286 * (cse_p11_t182_g372699 + dx23)) + (-((cse_p13_t1_g372995 * cse_p13_t45_g372943))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 2, h_y1_z3);
        double cse_p22_t1_g373231 = (cse_p12_t1_g372890 * dy13);
        double cse_p23_t1_g373240 = (cse_p12_t1_g372893 * dy13);
        double h_y1_x4 = ((cse_p11_t127_g372644 * cse_p15_t1_g373116) + (cse_p11_t1_g372727 * ((cse_p22_t1_g373231 * dx43) + (cse_p61_t10_invr373286 * dz23) + (-((cse_p13_t1_g372989 * cse_p13_t48_g372946))) + (-((cse_p13_t1_g372995 * cx))) + (-((cse_p23_t1_g373240 * dx43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x4, 0, h_y1_x4);
        double h_y1_y4 = ((cse_p11_t128_g372645 * cse_p15_t1_g373116) + (cse_p11_t1_g372727 * ((cse_p22_t1_g373231 * dy43) + (-((cse_p13_t1_g372993 * cse_p13_t48_g372946))) + (-((cse_p13_t1_g372995 * cy))) + (-((cse_p23_t1_g373240 * dy43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x4, 1, h_y1_y4);
        double h_y1_z4 = ((cse_p11_t129_g372646 * cse_p15_t1_g373116) + (cse_p11_t1_g372727 * ((cse_p16_t20_g373142 * cse_p61_t10_invr373286) + (cse_p22_t1_g373231 * dz43) + (-((cse_p13_t1_g372995 * cz))) + (-((cse_p13_t1_g372996 * cse_p13_t48_g372946))) + (-((cse_p23_t1_g373240 * dz43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x4, 2, h_y1_z4);
        double h_z1_z1 = ((cse_p11_t110_g372627 * cse_p14_t1_g373070) + (cse_p11_t1_g372727 * (cse_p11_t176_g372693 + (-2.00000000000000000e+0 * cse_p13_t44_g372942 * cse_p14_t1_g373066) + (cse_p12_t1_g372895 * cse_p16_t1_g373151) + (cse_p16_t1_g373151 * cse_p19_t1_g373190) + (-((cse_p14_t1_g373066 * cse_p17_t1_g373166 * dz13))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = ((cse_p11_t117_g372634 * cse_p14_t1_g373070) + (cse_p11_t1_g372727 * ((cse_p13_t1_g372998 * cse_p21_t1_g373212) + (cse_p61_t10_invr373286 * dy43) + (-((cse_p13_t1_g372990 * cse_p13_t44_g372942))) + (-((cse_p13_t51_g372949 * cse_p14_t1_g373066))) + (-((cse_p16_t1_g373147 * cse_p21_t1_g373212))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = ((cse_p11_t110_g372627 * cse_p14_t1_g373069) + (cse_p11_t1_g372727 * ((cse_p11_t182_g372699 * cse_p61_t10_invr373286) + (cse_p13_t1_g372998 * cse_p21_t1_g373216) + (-((cse_p13_t1_g372994 * cse_p13_t44_g372942))) + (-((cse_p13_t43_g372941 * cse_p14_t1_g373066))) + (-((cse_p16_t1_g373147 * cse_p21_t1_g373216))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = ((cse_p11_t113_g372630 * cse_p14_t1_g373070) + (cse_p11_t1_g372727 * ((cse_p13_t1_g372998 * cse_p26_t1_g373262) + (-((cse_p13_t1_g372997 * cse_p13_t44_g372942))) + (-((cse_p13_t47_g372945 * cse_p14_t1_g373066))) + (-(cse_p11_t81_g372598)))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_z1_x3 = ((cse_p11_t116_g372633 * cse_p14_t1_g373070) + (cse_p11_t1_g372727 * (cse_p11_t5_g372522 + (cse_p13_t1_g373004 * dz13) + (cse_p13_t44_g372942 * cse_p13_t61_g372959) + (cse_p23_t1_g373236 * dz13) + (cse_p61_t10_invr373286 * (cse_p11_t185_g372702 + dy23)) + (-((cse_p13_t50_g372948 * cse_p14_t1_g373066))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 0, h_z1_x3);
        double h_z1_y3 = ((cse_p11_t112_g372629 * cse_p14_t1_g373070) + (cse_p11_t1_g372727 * (cse_p11_t8_g372525 + (cse_p13_t1_g373005 * dz13) + (cse_p13_t44_g372942 * cse_p13_t62_g372960) + (cse_p23_t1_g373237 * dz13) + (cse_p61_t10_invr373286 * (cse_p16_t20_g373142 + dx43)) + (-((cse_p13_t46_g372944 * cse_p14_t1_g373066))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 1, h_z1_y3);
        double h_z1_z3 = ((cse_p11_t111_g372628 * cse_p14_t1_g373070) + (cse_p11_t1_g372727 * (cse_p11_t81_g372598 + (cse_p13_t44_g372942 * cse_p13_t63_g372961) + (cse_p16_t1_g373151 * cse_p17_t1_g373170) + (cse_p19_t1_g373192 * (r43 + (cse_p14_t1_g373073 * dz13) + (cse_p16_t1_g373151 * cse_p24_t1_g373251))) + (cse_p23_t1_g373238 * dz13) + (-((cse_p13_t45_g372943 * cse_p14_t1_g373066))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 2, h_z1_z3);
        double cse_p23_t1_g373232 = (cse_p12_t1_g372890 * dz13);
        double cse_p23_t1_g373241 = (cse_p12_t1_g372893 * dz13);
        double h_z1_x4 = ((cse_p11_t127_g372644 * cse_p14_t1_g373070) + (cse_p11_t1_g372727 * ((cse_p16_t22_g373144 * cse_p61_t10_invr373286) + (cse_p23_t1_g373232 * dx43) + (-((cse_p13_t1_g372989 * cse_p13_t44_g372942))) + (-((cse_p14_t1_g373066 * cx))) + (-((cse_p23_t1_g373241 * dx43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x4, 0, h_z1_x4);
        double h_z1_y4 = ((cse_p11_t128_g372645 * cse_p14_t1_g373070) + (cse_p11_t1_g372727 * ((cse_p23_t1_g373232 * dy43) + (cse_p61_t10_invr373286 * dx23) + (-((cse_p13_t1_g372993 * cse_p13_t44_g372942))) + (-((cse_p14_t1_g373066 * cy))) + (-((cse_p23_t1_g373241 * dy43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x4, 1, h_z1_y4);
        double h_z1_z4 = ((cse_p11_t129_g372646 * cse_p14_t1_g373070) + (cse_p11_t1_g372727 * ((cse_p23_t1_g373232 * dz43) + (-((cse_p13_t1_g372996 * cse_p13_t44_g372942))) + (-((cse_p14_t1_g373066 * cz))) + (-((cse_p23_t1_g373241 * dz43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x4, 2, h_z1_z4);
        double cse_p12_t1_g372896 = (cse_p13_t1_g373000 * cse_p22_t1_g373226);
        double cse_p17_t1_g373167 = (cse_p520_t4_invr373342 * v);
        double cse_p18_t1_g373184 = (cse_p11_t117_g372634 * cse_p12_t1_g372879);
        double cse_p19_t1_g373191 = (cse_p13_t1_g372992 * cse_p26_t1_g373269 * cse_p520_t3_sqrt373341 * r13);
        double cse_p24_t1_g373246 = (cse_p13_t1_g372990 * cse_p17_t1_g373167);
        double h_x2_x2 = ((cse_p11_t117_g372634 * cse_p18_t1_g373184) + (cse_p11_t1_g372727 * (cse_p11_t175_g372692 + (-2.00000000000000000e+0 * cse_p13_t1_g372990 * cse_p13_t51_g372949) + (cse_p12_t1_g372896 * cse_p16_t1_g373148) + (cse_p16_t1_g373148 * cse_p19_t1_g373191) + (-((cse_p24_t1_g373246 * dx23))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double cse_p21_t1_g373219 = (cse_p11_t125_g372642 * cse_p19_t1_g373186 * v);
        double h_x2_y2 = ((cse_p11_t117_g372634 * cse_p14_t1_g373069) + (cse_p11_t1_g372727 * ((cse_p20_t1_g373211 * cse_p21_t1_g373219) + (-((cse_p13_t1_g372990 * cse_p13_t43_g372941))) + (-((cse_p13_t1_g372994 * cse_p13_t51_g372949))) + (-((cse_p24_t1_g373246 * dy23))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = ((cse_p11_t117_g372634 * cse_p15_t1_g373115) + (cse_p11_t1_g372727 * ((cse_p21_t1_g373213 * cse_p21_t1_g373219) + (-((cse_p13_t1_g372990 * cse_p13_t47_g372945))) + (-((cse_p13_t1_g372997 * cse_p13_t51_g372949))) + (-((cse_p24_t1_g373246 * dz23))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double cse_p12_t1_g372891 = (cse_p12_t1_g372883 * cse_p17_t1_g373173);
        double cse_p17_t1_g373172 = (cse_p11_t116_g372633 * cse_p12_t1_g372879);
        double cse_p19_t1_g373193 = (cse_p19_t1_g373195 * cse_p63_t1_g373290);
        double cse_p23_t1_g373233 = (cse_p12_t1_g372891 * cse_p15_t10_g373085);
        double cse_p25_t1_g373252 = (cse_p520_t4_invr373342 * r43);
        double h_x2_x3 = ((cse_p11_t117_g372634 * cse_p17_t1_g373172) + (cse_p11_t1_g372727 * (cse_p11_t32_g372549 + (cse_p13_t51_g372949 * cse_p13_t61_g372959) + (cse_p16_t1_g373148 * cse_p17_t1_g373169) + (cse_p19_t1_g373193 * (r43 + (cse_p14_t1_g373071 * dx23) + (cse_p16_t1_g373148 * cse_p25_t1_g373252))) + (cse_p23_t1_g373233 * dx23) + (-((cse_p13_t1_g372990 * cse_p13_t50_g372948))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 0, h_x2_x3);
        double cse_p24_t1_g373249 = (cse_p13_t5_g372903 * cse_p63_t1_g373290);
        double cse_p13_t1_g373002 = (cse_p11_t119_g372636 * cse_p24_t1_g373249);
        double cse_p23_t1_g373234 = (cse_p12_t1_g372891 * cse_p15_t11_g373086);
        double h_x2_y3 = ((cse_p11_t117_g372634 * cse_p15_t1_g373112) + (cse_p11_t1_g372727 * (cse_p11_t6_g372523 + (cse_p13_t1_g373002 * dx23) + (cse_p13_t51_g372949 * cse_p13_t62_g372960) + (cse_p23_t1_g373234 * dx23) + (cse_p61_t10_invr373286 * (cse_p16_t23_g373145 + dz43)) + (-((cse_p13_t1_g372990 * cse_p13_t46_g372944))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 1, h_x2_y3);
        double cse_p13_t1_g373003 = (cse_p11_t120_g372637 * cse_p24_t1_g373249);
        double cse_p23_t1_g373235 = (cse_p12_t1_g372891 * cse_p15_t12_g373087);
        double h_x2_z3 = ((cse_p11_t117_g372634 * cse_p14_t1_g373074) + (cse_p11_t1_g372727 * (cse_p11_t7_g372524 + (cse_p13_t1_g373003 * dx23) + (cse_p13_t51_g372949 * cse_p13_t63_g372961) + (cse_p23_t1_g373235 * dx23) + (cse_p61_t10_invr373286 * (cse_p11_t185_g372702 + dy13)) + (-((cse_p13_t1_g372990 * cse_p13_t45_g372943))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 2, h_x2_z3);
        double cse_p12_t1_g372889 = (cse_p17_t1_g373173 * cse_p17_t1_g373175);
        double cse_p12_t1_g372894 = (cse_p13_t6_g372904 * cse_p19_t1_g373193);
        double cse_p22_t1_g373227 = (cse_p12_t1_g372889 * dx23);
        double cse_p24_t1_g373242 = (cse_p12_t1_g372894 * dx23);
        double h_x2_x4 = ((cse_p11_t127_g372644 * cse_p18_t1_g373184) + (cse_p11_t1_g372727 * ((cse_p22_t1_g373227 * dx43) + (-((cse_p13_t1_g372989 * cse_p13_t51_g372949))) + (-((cse_p13_t1_g372990 * cx))) + (-((cse_p24_t1_g373242 * dx43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x4, 0, h_x2_x4);
        double h_x2_y4 = ((cse_p11_t128_g372645 * cse_p18_t1_g373184) + (cse_p11_t1_g372727 * ((cse_p22_t1_g373227 * dy43) + (cse_p61_t10_invr373286 * dz13) + (-((cse_p13_t1_g372990 * cy))) + (-((cse_p13_t1_g372993 * cse_p13_t51_g372949))) + (-((cse_p24_t1_g373242 * dy43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x4, 1, h_x2_y4);
        double h_x2_z4 = ((cse_p11_t129_g372646 * cse_p18_t1_g373184) + (cse_p11_t1_g372727 * ((cse_p16_t21_g373143 * cse_p61_t10_invr373286) + (cse_p22_t1_g373227 * dz43) + (-((cse_p13_t1_g372990 * cz))) + (-((cse_p13_t1_g372996 * cse_p13_t51_g372949))) + (-((cse_p24_t1_g373242 * dz43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x4, 2, h_x2_z4);
        double cse_p26_t1_g373266 = (cse_p13_t1_g372994 * cse_p17_t1_g373167);
        double h_y2_y2 = ((cse_p11_t109_g372626 * cse_p14_t1_g373069) + (cse_p11_t1_g372727 * (cse_p11_t175_g372692 + (-2.00000000000000000e+0 * cse_p13_t1_g372994 * cse_p13_t43_g372941) + (cse_p12_t1_g372896 * cse_p16_t1_g373150) + (cse_p16_t1_g373150 * cse_p19_t1_g373191) + (-((cse_p26_t1_g373266 * dy23))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = ((cse_p11_t113_g372630 * cse_p14_t1_g373069) + (cse_p11_t1_g372727 * ((cse_p21_t1_g373217 * cse_p21_t1_g373219) + (-((cse_p13_t1_g372994 * cse_p13_t47_g372945))) + (-((cse_p13_t1_g372997 * cse_p13_t43_g372941))) + (-((cse_p26_t1_g373266 * dz23))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double cse_p12_t1_g372898 = (cse_p11_t118_g372635 * cse_p24_t1_g373249);
        double h_y2_x3 = ((cse_p11_t116_g372633 * cse_p14_t1_g373069) + (cse_p11_t1_g372727 * (cse_p11_t6_g372523 + (cse_p12_t1_g372898 * dy23) + (cse_p13_t43_g372941 * cse_p13_t61_g372959) + (cse_p23_t1_g373233 * dy23) + (cse_p61_t10_invr373286 * (cse_p11_t188_g372705 + dz13)) + (-((cse_p13_t1_g372994 * cse_p13_t50_g372948))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 0, h_y2_x3);
        double h_y2_y3 = ((cse_p11_t112_g372629 * cse_p14_t1_g373069) + (cse_p11_t1_g372727 * (cse_p11_t58_g372575 + (cse_p13_t43_g372941 * cse_p13_t62_g372960) + (cse_p16_t1_g373150 * cse_p17_t1_g373169) + (cse_p19_t1_g373193 * (r43 + (cse_p14_t1_g373072 * dy23) + (cse_p16_t1_g373150 * cse_p25_t1_g373252))) + (cse_p23_t1_g373234 * dy23) + (-((cse_p13_t1_g372994 * cse_p13_t46_g372944))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 1, h_y2_y3);
        double h_y2_z3 = ((cse_p11_t111_g372628 * cse_p14_t1_g373069) + (cse_p11_t1_g372727 * (cse_p11_t9_g372526 + (cse_p13_t1_g373003 * dy23) + (cse_p13_t43_g372941 * cse_p13_t63_g372961) + (cse_p23_t1_g373235 * dy23) + (cse_p61_t10_invr373286 * (cse_p16_t19_g373141 + dx43)) + (-((cse_p13_t1_g372994 * cse_p13_t45_g372943))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 2, h_y2_z3);
        double cse_p22_t1_g373228 = (cse_p12_t1_g372889 * dy23);
        double cse_p24_t1_g373243 = (cse_p12_t1_g372894 * dy23);
        double h_y2_x4 = ((cse_p11_t127_g372644 * cse_p14_t1_g373069) + (cse_p11_t1_g372727 * ((cse_p16_t23_g373145 * cse_p61_t10_invr373286) + (cse_p22_t1_g373228 * dx43) + (-((cse_p13_t1_g372989 * cse_p13_t43_g372941))) + (-((cse_p13_t1_g372994 * cx))) + (-((cse_p24_t1_g373243 * dx43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x4, 0, h_y2_x4);
        double h_y2_y4 = ((cse_p11_t128_g372645 * cse_p14_t1_g373069) + (cse_p11_t1_g372727 * ((cse_p22_t1_g373228 * dy43) + (-((cse_p13_t1_g372993 * cse_p13_t43_g372941))) + (-((cse_p13_t1_g372994 * cy))) + (-((cse_p24_t1_g373243 * dy43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x4, 1, h_y2_y4);
        double h_y2_z4 = ((cse_p11_t129_g372646 * cse_p14_t1_g373069) + (cse_p11_t1_g372727 * ((cse_p22_t1_g373228 * dz43) + (cse_p61_t10_invr373286 * dx13) + (-((cse_p13_t1_g372994 * cz))) + (-((cse_p13_t1_g372996 * cse_p13_t43_g372941))) + (-((cse_p24_t1_g373243 * dz43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x4, 2, h_y2_z4);
        double h_z2_z2 = ((cse_p11_t113_g372630 * cse_p15_t1_g373115) + (cse_p11_t1_g372727 * (cse_p11_t175_g372692 + (-2.00000000000000000e+0 * cse_p13_t1_g372997 * cse_p13_t47_g372945) + (cse_p12_t1_g372896 * cse_p16_t1_g373152) + (cse_p16_t1_g373152 * cse_p19_t1_g373191) + (-((cse_p13_t1_g372997 * cse_p17_t1_g373167 * dz23))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        double h_z2_x3 = ((cse_p11_t116_g372633 * cse_p15_t1_g373115) + (cse_p11_t1_g372727 * (cse_p11_t7_g372524 + (cse_p12_t1_g372898 * dz23) + (cse_p13_t47_g372945 * cse_p13_t61_g372959) + (cse_p23_t1_g373233 * dz23) + (cse_p61_t10_invr373286 * (cse_p16_t21_g373143 + dy43)) + (-((cse_p13_t1_g372997 * cse_p13_t50_g372948))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 0, h_z2_x3);
        double h_z2_y3 = ((cse_p11_t113_g372630 * cse_p15_t1_g373112) + (cse_p11_t1_g372727 * (cse_p11_t9_g372526 + (cse_p13_t1_g373002 * dz23) + (cse_p13_t47_g372945 * cse_p13_t62_g372960) + (cse_p23_t1_g373234 * dz23) + (cse_p61_t10_invr373286 * (cse_p11_t182_g372699 + dx13)) + (-((cse_p13_t1_g372997 * cse_p13_t46_g372944))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 1, h_z2_y3);
        double h_z2_z3 = ((cse_p11_t113_g372630 * cse_p14_t1_g373074) + (cse_p11_t1_g372727 * (cse_p11_t81_g372598 + (cse_p13_t47_g372945 * cse_p13_t63_g372961) + (cse_p16_t1_g373152 * cse_p17_t1_g373169) + (cse_p19_t1_g373193 * (r43 + (cse_p14_t1_g373073 * dz23) + (cse_p16_t1_g373152 * cse_p25_t1_g373252))) + (cse_p23_t1_g373235 * dz23) + (-((cse_p13_t1_g372997 * cse_p13_t45_g372943))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 2, h_z2_z3);
        double cse_p22_t1_g373229 = (cse_p12_t1_g372889 * dz23);
        double cse_p24_t1_g373244 = (cse_p12_t1_g372894 * dz23);
        double h_z2_x4 = ((cse_p11_t127_g372644 * cse_p15_t1_g373115) + (cse_p11_t1_g372727 * ((cse_p22_t1_g373229 * dx43) + (cse_p61_t10_invr373286 * dy13) + (-((cse_p13_t1_g372989 * cse_p13_t47_g372945))) + (-((cse_p13_t1_g372997 * cx))) + (-((cse_p24_t1_g373244 * dx43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x4, 0, h_z2_x4);
        double h_z2_y4 = ((cse_p11_t128_g372645 * cse_p15_t1_g373115) + (cse_p11_t1_g372727 * ((cse_p16_t19_g373141 * cse_p61_t10_invr373286) + (cse_p22_t1_g373229 * dy43) + (-((cse_p13_t1_g372993 * cse_p13_t47_g372945))) + (-((cse_p13_t1_g372997 * cy))) + (-((cse_p24_t1_g373244 * dy43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x4, 1, h_z2_y4);
        double h_z2_z4 = ((cse_p11_t129_g372646 * cse_p15_t1_g373115) + (cse_p11_t1_g372727 * ((cse_p22_t1_g373229 * dz43) + (-((cse_p13_t1_g372996 * cse_p13_t47_g372945))) + (-((cse_p13_t1_g372997 * cz))) + (-((cse_p24_t1_g373244 * dz43))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x4, 2, h_z2_z4);
        double cse_p16_t1_g373153 = (-2.00000000000000000e+0 * cse_p13_t1_g372999);
        double cse_p16_t1_g373154 = (-2.00000000000000000e+0 * cse_p61_t12_invr3373288);
        double cse_p25_t1_g373254 = (cse_p15_t10_g373085 * cse_p16_t1_g373154);
        double h_x3_x3 = ((cse_p11_t116_g372633 * cse_p17_t1_g373172) + (cse_p11_t1_g372727 * ((2.00000000000000000e+0 * cse_p13_t50_g372948 * cse_p13_t61_g372959) + (v * ((-((cse_p15_t10_g373085 * cse_p25_t1_g373254))) + (-((cse_p61_t11_invr2373287 * (cse_p11_t31_g372548 + cse_p11_t31_g372548 + cse_p17_t1_g373173 + cse_p17_t1_g373174 + cse_p18_t1_g373185 + (cse_p14_t1_g373067 * cse_p16_t1_g373148) + (cse_p14_t1_g373068 * cse_p15_t1_g373119) + (cse_p14_t1_g373071 * (cse_p20_t1_g373196 + cse_p20_t1_g373196 + cse_p20_t4_g373199 + cse_p20_t4_g373199 + (cse_p15_t1_g373117 * cse_p520_t6_invr373344))) + (cse_p15_t1_g373119 * cse_p26_t1_g373268) + (cse_p16_t1_g373148 * cse_p26_t1_g373269) + (cse_p16_t1_g373153 * cse_p17_t1_g373171))))))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 0, h_x3_x3);
        double cse_p27_t1_g373274 = (cse_p16_t1_g373153 * dx43);
        double h_x3_y3 = ((cse_p11_t116_g372633 * cse_p15_t1_g373112) + (cse_p11_t1_g372727 * ((cse_p13_t46_g372944 * cse_p13_t61_g372959) + (cse_p13_t50_g372948 * cse_p13_t62_g372960) + (v * ((-((cse_p15_t11_g373086 * cse_p25_t1_g373254))) + (-((cse_p61_t11_invr2373287 * ((cse_p13_t6_g372904 * ((cse_p15_t1_g373117 * cse_p17_t1_g373168) + (cse_p20_t1_g373196 * dy43) + (cse_p20_t2_g373197 * dx43) + (cse_p20_t4_g373199 * dy43) + (cse_p20_t5_g373200 * dx43))) + (cse_p14_t1_g373067 * cse_p20_t1_g373211) + (cse_p14_t1_g373068 * cse_p20_t1_g373206) + (cse_p15_t1_g373118 * cse_p20_t1_g373207) + (cse_p15_t1_g373118 * cse_p20_t1_g373210) + (cse_p20_t1_g373206 * cse_p26_t1_g373268) + (cse_p20_t1_g373211 * cse_p26_t1_g373269) + (cse_p27_t1_g373274 * dy43))))))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 1, h_x3_y3);
        double h_x3_z3 = ((cse_p11_t116_g372633 * cse_p14_t1_g373074) + (cse_p11_t1_g372727 * ((cse_p13_t45_g372943 * cse_p13_t61_g372959) + (cse_p13_t50_g372948 * cse_p13_t63_g372961) + (v * ((-((cse_p15_t12_g373087 * cse_p25_t1_g373254))) + (-((cse_p61_t11_invr2373287 * ((cse_p13_t6_g372904 * ((cse_p15_t1_g373113 * cse_p15_t1_g373117) + (cse_p20_t1_g373196 * dz43) + (cse_p20_t3_g373198 * dx43) + (cse_p20_t4_g373199 * dz43) + (cse_p20_t6_g373201 * dx43))) + (cse_p14_t1_g373067 * cse_p21_t1_g373213) + (cse_p14_t1_g373068 * cse_p20_t1_g373208) + (cse_p15_t1_g373118 * cse_p20_t1_g373209) + (cse_p15_t1_g373118 * cse_p21_t1_g373212) + (cse_p20_t1_g373208 * cse_p26_t1_g373268) + (cse_p21_t1_g373213 * cse_p26_t1_g373269) + (cse_p27_t1_g373274 * dz43))))))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 2, h_x3_z3);
        double cse_p15_t1_g373111 = (-2.00000000000000000e+0 * cse_p11_t1_g372731);
        double cse_p25_t1_g373253 = (cse_p15_t10_g373085 * cse_p15_t1_g373111);
        double h_x3_x4 = ((cse_p11_t127_g372644 * cse_p17_t1_g373172) + (cse_p11_t1_g372727 * ((cse_p13_t61_g372959 * cx) + (v * ((-((cse_p25_t1_g373253 * dx43))) + (-((cse_p61_t11_invr2373287 * (cse_p11_t173_g372690 + (cse_p16_t1_g373155 * cse_p17_t1_g373171) + (-((cse_p14_t1_g373071 * cse_p20_t1_g373196))) + (-((cse_p14_t1_g373071 * cse_p20_t4_g373199))) + (-((cse_p17_t1_g373171 * cse_p26_t1_g373270))))))))) + (-((cse_p13_t1_g372989 * cse_p13_t50_g372948))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x4, 0, h_x3_x4);
        double h_x3_y4 = ((cse_p11_t128_g372645 * cse_p17_t1_g373172) + (cse_p11_t1_g372727 * ((cse_p13_t61_g372959 * cy) + (cse_p15_t18_g373093 * cse_p61_t10_invr373286) + (v * ((-((cse_p25_t1_g373253 * dy43))) + (-((cse_p61_t11_invr2373287 * (cse_p11_t13_g372530 + cse_p11_t156_g372673 + (-((cse_p14_t1_g373072 * cse_p20_t1_g373196))) + (-((cse_p14_t1_g373072 * cse_p20_t4_g373199))))))))) + (-((cse_p13_t1_g372993 * cse_p13_t50_g372948))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x4, 1, h_x3_y4);
        double h_x3_z4 = ((cse_p11_t129_g372646 * cse_p17_t1_g373172) + (cse_p11_t1_g372727 * ((cse_p13_t61_g372959 * cz) + (cse_p15_t15_g373090 * cse_p61_t10_invr373286) + (v * ((-((cse_p25_t1_g373253 * dz43))) + (-((cse_p61_t11_invr2373287 * (cse_p11_t14_g372531 + cse_p11_t157_g372674 + (-((cse_p14_t1_g373073 * cse_p20_t1_g373196))) + (-((cse_p14_t1_g373073 * cse_p20_t4_g373199))))))))) + (-((cse_p13_t1_g372996 * cse_p13_t50_g372948))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x4, 2, h_x3_z4);
        double cse_p26_t1_g373271 = (cse_p15_t11_g373086 * cse_p16_t1_g373154);
        double h_y3_y3 = ((cse_p11_t112_g372629 * cse_p15_t1_g373112) + (cse_p11_t1_g372727 * ((2.00000000000000000e+0 * cse_p13_t46_g372944 * cse_p13_t62_g372960) + (v * ((-((cse_p15_t11_g373086 * cse_p26_t1_g373271))) + (-((cse_p61_t11_invr2373287 * (cse_p11_t57_g372574 + cse_p11_t57_g372574 + cse_p17_t1_g373173 + cse_p17_t1_g373174 + cse_p18_t1_g373185 + (cse_p14_t1_g373067 * cse_p16_t1_g373150) + (cse_p14_t1_g373068 * cse_p16_t1_g373149) + (cse_p14_t1_g373072 * (cse_p20_t2_g373197 + cse_p20_t2_g373197 + cse_p20_t5_g373200 + cse_p20_t5_g373200 + (cse_p11_t1_g372718 * cse_p17_t1_g373168))) + (cse_p16_t1_g373149 * cse_p26_t1_g373268) + (cse_p16_t1_g373150 * cse_p26_t1_g373269) + (cse_p16_t1_g373153 * cse_p18_t1_g373182))))))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 1, h_y3_y3);
        double cse_p26_t1_g373264 = (cse_p11_t1_g372718 * cse_p15_t1_g373113);
        double h_y3_z3 = ((cse_p11_t112_g372629 * cse_p14_t1_g373074) + (cse_p11_t1_g372727 * ((cse_p13_t45_g372943 * cse_p13_t62_g372960) + (cse_p13_t46_g372944 * cse_p13_t63_g372961) + (v * ((-((cse_p15_t12_g373087 * cse_p26_t1_g373271))) + (-((cse_p61_t11_invr2373287 * ((cse_p13_t6_g372904 * ((cse_p20_t2_g373197 * dz43) + (cse_p20_t3_g373198 * dy43) + (cse_p20_t5_g373200 * dz43) + (cse_p20_t6_g373201 * dy43) + (cse_p26_t1_g373264 * dy43))) + (cse_p14_t1_g373067 * cse_p21_t1_g373217) + (cse_p14_t1_g373068 * cse_p21_t1_g373214) + (cse_p15_t1_g373118 * cse_p21_t1_g373215) + (cse_p15_t1_g373118 * cse_p21_t1_g373216) + (cse_p16_t1_g373153 * cse_p25_t1_g373260) + (cse_p21_t1_g373214 * cse_p26_t1_g373268) + (cse_p21_t1_g373217 * cse_p26_t1_g373269))))))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 2, h_y3_z3);
        double cse_p25_t1_g373255 = (cse_p15_t11_g373086 * cse_p15_t1_g373111);
        double h_y3_x4 = ((cse_p11_t127_g372644 * cse_p15_t1_g373112) + (cse_p11_t1_g372727 * ((cse_p13_t62_g372960 * cx) + (cse_p15_t17_g373092 * cse_p61_t10_invr373286) + (v * ((-((cse_p25_t1_g373255 * dx43))) + (-((cse_p61_t11_invr2373287 * (cse_p11_t13_g372530 + cse_p11_t156_g372673 + (-((cse_p14_t1_g373071 * cse_p20_t2_g373197))) + (-((cse_p14_t1_g373071 * cse_p20_t5_g373200))))))))) + (-((cse_p13_t1_g372989 * cse_p13_t46_g372944))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x4, 0, h_y3_x4);
        double h_y3_y4 = ((cse_p11_t128_g372645 * cse_p15_t1_g373112) + (cse_p11_t1_g372727 * ((cse_p13_t62_g372960 * cy) + (v * ((-((cse_p25_t1_g373255 * dy43))) + (-((cse_p61_t11_invr2373287 * (cse_p11_t173_g372690 + (cse_p16_t1_g373155 * cse_p18_t1_g373182) + (-((cse_p14_t1_g373072 * cse_p20_t2_g373197))) + (-((cse_p14_t1_g373072 * cse_p20_t5_g373200))) + (-((cse_p18_t1_g373182 * cse_p26_t1_g373270))))))))) + (-((cse_p13_t1_g372993 * cse_p13_t46_g372944))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x4, 1, h_y3_y4);
        double h_y3_z4 = ((cse_p11_t129_g372646 * cse_p15_t1_g373112) + (cse_p11_t1_g372727 * ((cse_p13_t62_g372960 * cz) + (cse_p15_t14_g373089 * cse_p61_t10_invr373286) + (v * ((-((cse_p25_t1_g373255 * dz43))) + (-((cse_p61_t11_invr2373287 * (cse_p11_t15_g372532 + cse_p11_t164_g372681 + (-((cse_p14_t1_g373073 * cse_p20_t2_g373197))) + (-((cse_p14_t1_g373073 * cse_p20_t5_g373200))))))))) + (-((cse_p13_t1_g372996 * cse_p13_t46_g372944))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x4, 2, h_y3_z4);
        double h_z3_z3 = ((cse_p11_t111_g372628 * cse_p14_t1_g373074) + (cse_p11_t1_g372727 * ((2.00000000000000000e+0 * cse_p13_t45_g372943 * cse_p13_t63_g372961) + (v * ((-((cse_p15_t12_g373087 * cse_p15_t12_g373087 * cse_p16_t1_g373154))) + (-((cse_p61_t11_invr2373287 * (cse_p11_t80_g372597 + cse_p11_t80_g372597 + cse_p17_t1_g373173 + cse_p17_t1_g373174 + cse_p18_t1_g373185 + (cse_p14_t1_g373067 * cse_p16_t1_g373152) + (cse_p14_t1_g373068 * cse_p16_t1_g373151) + (cse_p14_t1_g373073 * (cse_p20_t3_g373198 + cse_p20_t3_g373198 + cse_p20_t6_g373201 + cse_p20_t6_g373201 + cse_p26_t1_g373264)) + (cse_p16_t1_g373151 * cse_p26_t1_g373268) + (cse_p16_t1_g373152 * cse_p26_t1_g373269) + (cse_p16_t1_g373153 * cse_p18_t1_g373183))))))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x3, 2, h_z3_z3);
        double cse_p25_t1_g373256 = (cse_p15_t12_g373087 * cse_p15_t1_g373111);
        double h_z3_x4 = ((cse_p11_t127_g372644 * cse_p14_t1_g373074) + (cse_p11_t1_g372727 * ((cse_p13_t63_g372961 * cx) + (cse_p15_t16_g373091 * cse_p61_t10_invr373286) + (v * ((-((cse_p25_t1_g373256 * dx43))) + (-((cse_p61_t11_invr2373287 * (cse_p11_t14_g372531 + cse_p11_t157_g372674 + (-((cse_p14_t1_g373071 * cse_p20_t3_g373198))) + (-((cse_p14_t1_g373071 * cse_p20_t6_g373201))))))))) + (-((cse_p13_t1_g372989 * cse_p13_t45_g372943))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x4, 0, h_z3_x4);
        double h_z3_y4 = ((cse_p11_t128_g372645 * cse_p14_t1_g373074) + (cse_p11_t1_g372727 * ((cse_p13_t63_g372961 * cy) + (cse_p15_t13_g373088 * cse_p61_t10_invr373286) + (v * ((-((cse_p25_t1_g373256 * dy43))) + (-((cse_p61_t11_invr2373287 * (cse_p11_t15_g372532 + cse_p11_t164_g372681 + (-((cse_p14_t1_g373072 * cse_p20_t3_g373198))) + (-((cse_p14_t1_g373072 * cse_p20_t6_g373201))))))))) + (-((cse_p13_t1_g372993 * cse_p13_t45_g372943))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x4, 1, h_z3_y4);
        double h_z3_z4 = ((cse_p11_t129_g372646 * cse_p14_t1_g373074) + (cse_p11_t1_g372727 * ((cse_p13_t63_g372961 * cz) + (v * ((-((cse_p25_t1_g373256 * dz43))) + (-((cse_p61_t11_invr2373287 * (cse_p11_t173_g372690 + (cse_p16_t1_g373155 * cse_p18_t1_g373183) + (-((cse_p14_t1_g373073 * cse_p20_t3_g373198))) + (-((cse_p14_t1_g373073 * cse_p20_t6_g373201))) + (-((cse_p18_t1_g373183 * cse_p26_t1_g373270))))))))) + (-((cse_p13_t1_g372996 * cse_p13_t45_g372943))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x4, 2, h_z3_z4);
        double cse_p12_t1_g372897 = (cse_p13_t1_g372999 * cse_p22_t1_g373226);
        double cse_p20_t1_g373205 = (cse_p17_t1_g373175 * cse_p520_t5_sqrt373343 * r13 * r23);
        double cse_p22_t1_g373222 = (cse_p11_t127_g372644 * cse_p12_t1_g372879);
        double cse_p24_t1_g373245 = (cse_p13_t1_g372989 * v);
        double h_x4_x4 = ((cse_p11_t127_g372644 * cse_p22_t1_g373222) + (cse_p11_t1_g372727 * (cse_p11_t174_g372691 + (-2.00000000000000000e+0 * cse_p13_t1_g372989 * cx) + (cse_p12_t1_g372897 * cse_p17_t1_g373171) + (cse_p17_t1_g373171 * cse_p20_t1_g373205 * cse_p520_t6_invr373344) + (-((cse_p24_t1_g373245 * cse_p520_t6_invr373344 * dx43))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 0, term.i3x4, 0, h_x4_x4);
        double cse_p25_t1_g373261 = (cse_p11_t126_g372643 * cse_p520_t5_sqrt373343 * v);
        double cse_p27_t1_g373273 = (cse_p15_t1_g373117 * cse_p25_t1_g373261);
        double h_x4_y4 = ((cse_p11_t128_g372645 * cse_p22_t1_g373222) + (cse_p11_t1_g372727 * ((cse_p27_t1_g373273 * dy43) + (-((cse_p13_t1_g372989 * cy))) + (-((cse_p13_t1_g372993 * cx))) + (-((cse_p17_t1_g373168 * cse_p24_t1_g373245))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 0, term.i3x4, 1, h_x4_y4);
        double h_x4_z4 = ((cse_p11_t129_g372646 * cse_p22_t1_g373222) + (cse_p11_t1_g372727 * ((cse_p27_t1_g373273 * dz43) + (-((cse_p13_t1_g372989 * cz))) + (-((cse_p13_t1_g372996 * cx))) + (-((cse_p15_t1_g373113 * cse_p24_t1_g373245))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 0, term.i3x4, 2, h_x4_z4);
        double cse_p26_t1_g373263 = (cse_p11_t128_g372645 * cse_p12_t1_g372879);
        double cse_p26_t1_g373265 = (cse_p13_t1_g372993 * v);
        double h_y4_y4 = ((cse_p11_t128_g372645 * cse_p26_t1_g373263) + (cse_p11_t1_g372727 * (cse_p11_t174_g372691 + (-2.00000000000000000e+0 * cse_p13_t1_g372993 * cy) + (cse_p12_t1_g372897 * cse_p18_t1_g373182) + (cse_p17_t1_g373168 * cse_p20_t1_g373205 * dy43) + (-((cse_p17_t1_g373168 * cse_p26_t1_g373265))))));
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 1, term.i3x4, 1, h_y4_y4);
        double h_y4_z4 = ((cse_p11_t129_g372646 * cse_p26_t1_g373263) + (cse_p11_t1_g372727 * ((cse_p11_t1_g372718 * cse_p25_t1_g373260 * cse_p25_t1_g373261) + (-((cse_p13_t1_g372993 * cz))) + (-((cse_p13_t1_g372996 * cy))) + (-((cse_p15_t1_g373113 * cse_p26_t1_g373265))))));
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 1, term.i3x4, 2, h_y4_z4);
        double h_z4_z4 = ((cse_p11_t129_g372646 * cse_p11_t129_g372646 * cse_p12_t1_g372879) + (cse_p11_t1_g372727 * (cse_p11_t174_g372691 + (-2.00000000000000000e+0 * cse_p13_t1_g372996 * cz) + (cse_p12_t1_g372897 * cse_p18_t1_g373183) + (cse_p15_t1_g373113 * cse_p20_t1_g373205 * dz43) + (-((cse_p13_t1_g372996 * cse_p15_t1_g373113 * v))))));
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
