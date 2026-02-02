#include "amber_nonbond14.h"

template <typename HESSIAN>
struct Amber_Nonbond14 {
  static constexpr size_t PositionSize = 6;
  static std::string description() { return "mathkernel-amber_nonbond14"; };
double gradient(const AmberNonbond14Term& term, double* position, double* energy_accumulate, double* force) {
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
      double cse_p11_t1_g115553 = (-3.00000000000000000e+0 * term.c);
      double dx = (x1 + (-(x2)));
      double dy = (y1 + (-(y2)));
      double dz = (z1 + (-(z2)));
      double r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      double cse_p1_t1_sqrt115549 = sqrt(r2);
      double cse_p1_t2_invsqrt115550 = (1.0 / (cse_p1_t1_sqrt115549));
      double r = cse_p1_t1_sqrt115549;
      double cse_p1_t3_invr115551 = (1.0 / (r));
      double cse_p1_t4_invr2115552 = (cse_p1_t3_invr115551 * cse_p1_t3_invr115551);
      double invr = cse_p1_t3_invr115551;
      double invr2 = (invr * invr);
      double cse_p12_t1_g115564 = (invr2 * invr2);
      double invr6 = (cse_p12_t1_g115564 * invr2);
      double cse_p11_t1_g115563 = (term.a * invr6);
      double cse_p11_t2_g115554 = (6.00000000000000000e+0 * cse_p11_t1_g115563);
      double cse_p11_t3_g115555 = (cse_p11_t1_g115553 + cse_p11_t2_g115554);
      double e_lj = ((cse_p11_t1_g115563 * invr6) + (-((term.c * invr6))));
      double e_coul = (invr * term.kqq);
      double energy = (e_coul + e_lj);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t1_g115556 = (cse_p1_t2_invsqrt115550 * cse_p1_t4_invr2115552);
      double cse_p11_t1_g115557 = (cse_p11_t1_g115556 * cse_p11_t3_g115555 * cse_p12_t1_g115564 * invr);
      double cse_p11_t1_g115558 = (cse_p11_t1_g115556 * term.kqq);
      double cse_p11_t1_g115559 = (-(cse_p11_t1_g115558));
      double cse_p11_t1_g115560 = (-2.00000000000000000e+0 * cse_p11_t1_g115557);
      double g_x1 = ((cse_p11_t1_g115559 * dx) + (cse_p11_t1_g115560 * dx));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = ((cse_p11_t1_g115559 * dy) + (cse_p11_t1_g115560 * dy));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = ((cse_p11_t1_g115559 * dz) + (cse_p11_t1_g115560 * dz));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p11_t1_g115562 = (2.00000000000000000e+0 * cse_p11_t1_g115557);
      double g_x2 = ((cse_p11_t1_g115558 * dx) + (cse_p11_t1_g115562 * dx));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = ((cse_p11_t1_g115558 * dy) + (cse_p11_t1_g115562 * dy));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = ((cse_p11_t1_g115558 * dz) + (cse_p11_t1_g115562 * dz));
      KernelGradientAcc(term.i3x2, 2, g_z2);
    }
  }
  return energy_added;
}
void gradient_fd(const AmberNonbond14Term& term, double* position, double* energy_accumulate, double* force)
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

double energy(const AmberNonbond14Term& term, double* position, double* energy_accumulate) {
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
      double dx = (x1 + (-(x2)));
      double dy = (y1 + (-(y2)));
      double dz = (z1 + (-(z2)));
      double r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      double r = sqrt(r2);
      double invr = (1.0 / (r));
      double invr2 = (invr * invr);
      double invr6 = (invr2 * invr2 * invr2);
      double e_lj = ((term.a * invr6 * invr6) + (-((term.c * invr6))));
      double e_coul = (invr * term.kqq);
      double energy = (e_coul + e_lj);
      energy_added += energy;
      *energy_accumulate += energy;
    }
  }
  return energy_added;
}
void energy_fd(const AmberNonbond14Term& term, double* position, double* energy_accumulate)
{
  energy(term, position, energy_accumulate);
}

double hessian(const AmberNonbond14Term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
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
      double cse_p11_t5_g115573 = (-3.00000000000000000e+0 * term.c);
      double dx = (x1 + (-(x2)));
      double cse_p11_t11_g115579 = (dx * dx);
      double cse_p11_t28_g115596 = (-(cse_p11_t11_g115579));
      double cse_p11_t33_g115601 = (-(dx));
      double dy = (y1 + (-(y2)));
      double cse_p11_t15_g115583 = (dy * dy);
      double cse_p11_t29_g115597 = (-(cse_p11_t15_g115583));
      double cse_p11_t34_g115602 = (-(dy));
      double dz = (z1 + (-(z2)));
      double cse_p11_t19_g115587 = (dz * dz);
      double cse_p11_t30_g115598 = (-(cse_p11_t19_g115587));
      double cse_p11_t35_g115603 = (-(dz));
      double r2 = (cse_p11_t11_g115579 + cse_p11_t15_g115583 + cse_p11_t19_g115587);
      double cse_p1_t1_sqrt115565 = sqrt(r2);
      double cse_p1_t2_invsqrt115566 = (1.0 / (cse_p1_t1_sqrt115565));
      double r = cse_p1_t1_sqrt115565;
      double cse_p62_t1_g115627 = (r * r);
      double cse_p61_t1_invr115621 = (1.0 / (((r) * (r))));
      double cse_p11_t1_g115620 = (cse_p61_t1_invr115621 * term.kqq);
      double cse_p61_t4_invr115624 = (1.0 / ((r * r * r)));
      double cse_p61_t3_invr115623 = (1.0 / (pow(r, 7)));
      double cse_p166_t1_g115641 = (6.00000000000000000e+0 * term.c * cse_p61_t3_invr115623);
      double cse_p61_t2_invr115622 = (1.0 / (pow(r, 13)));
      double cse_p114_t1_g115640 = (-1.20000000000000000e+1 * term.a * cse_p61_t2_invr115622);
      double cse_p11_t1_g115616 = (cse_p61_t4_invr115624 * (cse_p114_t1_g115640 + cse_p166_t1_g115641 + (-(cse_p11_t1_g115620))));
      double cse_p11_t1_g115569 = (cse_p11_t1_g115616 * (cse_p11_t28_g115596 + cse_p62_t1_g115627));
      double cse_p11_t2_g115570 = (cse_p11_t1_g115616 * (cse_p11_t29_g115597 + cse_p62_t1_g115627));
      double cse_p11_t3_g115571 = (cse_p11_t1_g115616 * (cse_p11_t30_g115598 + cse_p62_t1_g115627));
      double cse_p61_t5_invr115625 = (1.0 / (pow(r, 8)));
      double cse_p11_t6_g115574 = (-4.20000000000000040e+1 * term.c * cse_p61_t5_invr115625);
      double cse_p61_t6_invr115626 = (1.0 / (pow(r, 14)));
      double cse_p11_t7_g115575 = (1.56000000000000000e+2 * term.a * cse_p61_t6_invr115626);
      double cse_p11_t8_g115576 = (2.00000000000000000e+0 * cse_p61_t4_invr115624 * term.kqq);
      double cse_p11_t12_g115580 = (cse_p11_t1_g115616 * dx * dy);
      double cse_p63_t1_g115633 = (cse_p11_t1_g115616 * dz);
      double cse_p11_t13_g115581 = (cse_p63_t1_g115633 * dx);
      double cse_p11_t14_g115582 = (cse_p11_t1_g115616 * cse_p11_t33_g115601 * dy);
      double cse_p11_t16_g115584 = (cse_p63_t1_g115633 * dy);
      double cse_p11_t17_g115585 = (cse_p11_t33_g115601 * cse_p63_t1_g115633);
      double cse_p11_t18_g115586 = (cse_p11_t34_g115602 * cse_p63_t1_g115633);
      double cse_p11_t22_g115590 = (cse_p114_t1_g115640 + cse_p166_t1_g115641 + (-(cse_p11_t1_g115620)));
      double cse_p11_t24_g115592 = (cse_p11_t6_g115574 + cse_p11_t7_g115575 + cse_p11_t8_g115576);
      double cse_p11_t25_g115593 = (cse_p11_t28_g115596 + cse_p62_t1_g115627);
      double cse_p11_t26_g115594 = (cse_p11_t29_g115597 + cse_p62_t1_g115627);
      double cse_p11_t27_g115595 = (cse_p11_t30_g115598 + cse_p62_t1_g115627);
      double cse_p11_t31_g115599 = (-(cse_p11_t1_g115620));
      double cse_p11_t32_g115600 = (-(cse_p62_t1_g115627));
      double cse_p11_t42_g115610 = pow(r, 13);
      double cse_p11_t43_g115611 = pow(r, 14);
      double cse_p11_t44_g115612 = ((r) * (r));
      double cse_p11_t45_g115613 = (r * r * r);
      double cse_p11_t46_g115614 = pow(r, 7);
      double cse_p11_t47_g115615 = pow(r, 8);
      double cse_p1_t3_invr115567 = (1.0 / (r));
      double cse_p1_t4_invr2115568 = (cse_p1_t3_invr115567 * cse_p1_t3_invr115567);
      double invr = cse_p1_t3_invr115567;
      double invr2 = (invr * invr);
      double invr6 = (invr2 * invr2 * invr2);
      double cse_p11_t9_g115577 = (6.00000000000000000e+0 * term.a * invr6);
      double cse_p11_t23_g115591 = (cse_p11_t5_g115573 + cse_p11_t9_g115577);
      double e_lj = ((term.a * invr6 * invr6) + (-((term.c * invr6))));
      double e_coul = (invr * term.kqq);
      double energy = (e_coul + e_lj);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t1_g115618 = (cse_p1_t2_invsqrt115566 * cse_p1_t4_invr2115568);
      double cse_p63_t1_g115638 = (cse_p11_t1_g115618 * term.kqq);
      double cse_p63_t1_g115628 = (cse_p63_t1_g115638 * dx);
      double cse_p63_t2_g115629 = (cse_p63_t1_g115638 * dy);
      double cse_p63_t3_g115630 = (cse_p63_t1_g115638 * dz);
      double cse_p11_t1_g115619 = (cse_p11_t1_g115618 * cse_p11_t23_g115591 * invr2 * invr2 * invr);
      double cse_p63_t1_g115635 = (-2.00000000000000000e+0 * cse_p11_t1_g115619);
      double g_x1 = ((cse_p63_t1_g115635 * dx) + (-(cse_p63_t1_g115628)));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = ((cse_p63_t1_g115635 * dy) + (-(cse_p63_t2_g115629)));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = ((cse_p63_t1_g115635 * dz) + (-(cse_p63_t3_g115630)));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p63_t1_g115636 = (2.00000000000000000e+0 * cse_p11_t1_g115619);
      double g_x2 = (cse_p63_t1_g115628 + (cse_p63_t1_g115636 * dx));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (cse_p63_t2_g115629 + (cse_p63_t1_g115636 * dy));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (cse_p63_t3_g115630 + (cse_p63_t1_g115636 * dz));
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double cse_p11_t1_g115617 = (cse_p11_t24_g115592 * cse_p1_t3_invr115567 * cse_p1_t3_invr115567);
      double cse_p63_t1_g115637 = (cse_p11_t1_g115617 * dx);
      double h_x1_x1 = (cse_p11_t1_g115569 + (cse_p63_t1_g115637 * dx));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
      double h_x1_y1 = (cse_p11_t14_g115582 + (cse_p63_t1_g115637 * dy));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
      double h_x1_z1 = (cse_p11_t17_g115585 + (cse_p63_t1_g115637 * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
      double cse_p63_t1_g115631 = (cse_p11_t1_g115617 * cse_p11_t33_g115601);
      double cse_p63_t1_g115639 = (cse_p11_t22_g115590 * cse_p61_t4_invr115624);
      double h_x1_x2 = ((cse_p63_t1_g115631 * dx) + (cse_p63_t1_g115639 * (cse_p11_t11_g115579 + cse_p11_t32_g115600)));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
      double cse_p63_t1_g115632 = (cse_p11_t1_g115617 * cse_p11_t34_g115602);
      double h_x1_y2 = (cse_p11_t12_g115580 + (cse_p63_t1_g115632 * dx));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
      double cse_p63_t1_g115634 = (cse_p11_t1_g115617 * cse_p11_t35_g115603);
      double h_x1_z2 = (cse_p11_t13_g115581 + (cse_p63_t1_g115634 * dx));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
      double h_y1_y1 = (cse_p11_t2_g115570 + (cse_p11_t1_g115617 * dy * dy));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
      double h_y1_z1 = (cse_p11_t18_g115586 + (cse_p11_t1_g115617 * dy * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
      double h_y1_x2 = (cse_p11_t12_g115580 + (cse_p63_t1_g115631 * dy));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
      double h_y1_y2 = ((cse_p63_t1_g115632 * dy) + (cse_p63_t1_g115639 * (cse_p11_t15_g115583 + cse_p11_t32_g115600)));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
      double h_y1_z2 = (cse_p11_t16_g115584 + (cse_p63_t1_g115634 * dy));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
      double h_z1_z1 = (cse_p11_t3_g115571 + (cse_p11_t1_g115617 * dz * dz));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
      double h_z1_x2 = (cse_p11_t13_g115581 + (cse_p63_t1_g115631 * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
      double h_z1_y2 = (cse_p11_t16_g115584 + (cse_p63_t1_g115632 * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
      double h_z1_z2 = ((cse_p63_t1_g115634 * dz) + (cse_p63_t1_g115639 * (cse_p11_t19_g115587 + cse_p11_t32_g115600)));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
      double h_x2_x2 = (cse_p11_t1_g115569 + (cse_p11_t33_g115601 * cse_p63_t1_g115631));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
      double h_x2_y2 = (cse_p11_t14_g115582 + (cse_p11_t34_g115602 * cse_p63_t1_g115631));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
      double h_x2_z2 = (cse_p11_t17_g115585 + (cse_p11_t35_g115603 * cse_p63_t1_g115631));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
      double h_y2_y2 = (cse_p11_t2_g115570 + (cse_p11_t34_g115602 * cse_p63_t1_g115632));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
      double h_y2_z2 = (cse_p11_t18_g115586 + (cse_p11_t35_g115603 * cse_p63_t1_g115632));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
      double h_z2_z2 = (cse_p11_t3_g115571 + (cse_p11_t35_g115603 * cse_p63_t1_g115634));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
    }
  }
  return energy_added;
}
void hessian_fd(const AmberNonbond14Term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
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

};
