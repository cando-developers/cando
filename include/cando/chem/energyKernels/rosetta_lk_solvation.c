#include "rosetta_lk_solvation.h"

template <typename HESSIAN>
struct Rosetta_Lk_Solvation {
  static constexpr size_t PositionSize = 6;
  static std::string description() { return "mathkernel-rosetta_lk_solvation"; };
double gradient(const rosetta_lk_solvation_parameters& params, const rosetta_lk_solvation_term& term, double* position, double* energy_accumulate, double* force) {
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
      double cse_p1_t1_sqrt285547 = sqrt(r2);
      double cse_p1_t2_invsqrt285548 = (1.0 / (cse_p1_t1_sqrt285547));
      double rr = cse_p1_t1_sqrt285547;
      if ((rr <= term.r0_low)) {
        /* !COND-1 */
        double energy = (term.f_contact_i + term.f_contact_j);
        energy_added += energy;
        *energy_accumulate += energy;
        double g_x1 = 0.00000000000000000e+0;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = 0.00000000000000000e+0;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = 0.00000000000000000e+0;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = 0.00000000000000000e+0;
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = 0.00000000000000000e+0;
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = 0.00000000000000000e+0;
        KernelGradientAcc(term.i3x2, 2, g_z2);
        double de_dr = 0.00000000000000000e+0;
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.r1_low)) {
          /* !COND-2 */
          double cse_p528_t4_g294965 = (2.00000000000000000e+0 * rr);
          double cse_p528_t1_g294980 = (cse_p1_t2_invsqrt285548 * dx);
          double cse_p528_t6_g294967 = (cse_p1_t2_invsqrt285548 * dy);
          double cse_p528_t7_g294968 = (cse_p1_t2_invsqrt285548 * dz);
          double cse_p528_t8_g294969 = (cse_p1_t2_invsqrt285548 * rr);
          double cse_p528_t9_g294970 = (rr * rr);
          double cse_p14_t2_g285665 = (3.00000000000000000e+0 * term.aa_low_j);
          double cse_p12_t1_g285640 = (3.00000000000000000e+0 * term.aa_low_i);
          double cse_p13_t11_g285651 = (cse_p12_t1_g285640 + cse_p14_t2_g285665);
          double cse_p12_t10_g285609 = (term.bb_low_i * rr);
          double cse_p12_t11_g285610 = (term.bb_low_j * rr);
          double cse_p12_t27_g285626 = (cse_p13_t11_g285651 * cse_p528_t9_g294970);
          double cse_p13_t1_g285655 = (term.bb_low_i * cse_p528_t8_g294969);
          double cse_p12_t7_g285606 = (cse_p13_t1_g285655 * dx);
          double cse_p12_t8_g285607 = (cse_p13_t1_g285655 * dy);
          double cse_p12_t9_g285608 = (cse_p13_t1_g285655 * dz);
          double cse_p12_t12_g285611 = (cse_p528_t8_g294969 * rr);
          double cse_p11_t27_g285579 = (term.cc_low_i + term.cc_low_j + cse_p12_t10_g285609 + cse_p12_t10_g285609 + cse_p12_t11_g285610 + cse_p12_t11_g285610 + cse_p12_t27_g285626);
          double cse_p11_t31_g285583 = (-(cse_p12_t7_g285606));
          double cse_p11_t32_g285584 = (-(cse_p12_t8_g285607));
          double cse_p11_t33_g285585 = (-(cse_p12_t9_g285608));
          double energy = (term.dd_low_i + term.dd_low_j + (term.cc_low_i * rr) + (term.cc_low_j * rr) + (cse_p528_t9_g294970 * (term.bb_low_i + term.bb_low_j + (term.aa_low_i * rr) + (term.aa_low_j * rr))));
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (cse_p11_t27_g285579 * cse_p528_t1_g294980);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t27_g285579 * cse_p528_t6_g294967);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t27_g285579 * cse_p528_t7_g294968);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double cse_p12_t1_g285632 = (cse_p12_t12_g285611 * dx);
          double cse_p528_t1_g294962 = (-3.00000000000000000e+0 * cse_p12_t1_g285632);
          double cse_p13_t1_g285656 = (term.bb_low_j * cse_p528_t8_g294969);
          double cse_p63_t1_g286432 = (-2.00000000000000000e+0 * cse_p13_t1_g285656);
          double g_x2 = (cse_p11_t31_g285583 + cse_p11_t31_g285583 + (term.aa_low_i * cse_p528_t1_g294962) + (term.aa_low_j * cse_p528_t1_g294962) + (cse_p63_t1_g286432 * dx) + (-((term.cc_low_i * cse_p528_t1_g294980))) + (-((term.cc_low_j * cse_p528_t1_g294980))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double cse_p12_t1_g285633 = (cse_p12_t12_g285611 * dy);
          double cse_p528_t2_g294963 = (-3.00000000000000000e+0 * cse_p12_t1_g285633);
          double g_y2 = (cse_p11_t32_g285584 + cse_p11_t32_g285584 + (term.aa_low_i * cse_p528_t2_g294963) + (term.aa_low_j * cse_p528_t2_g294963) + (cse_p63_t1_g286432 * dy) + (-((term.cc_low_i * cse_p528_t6_g294967))) + (-((term.cc_low_j * cse_p528_t6_g294967))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double cse_p12_t1_g285634 = (cse_p12_t12_g285611 * dz);
          double cse_p528_t3_g294964 = (-3.00000000000000000e+0 * cse_p12_t1_g285634);
          double g_z2 = (cse_p11_t33_g285585 + cse_p11_t33_g285585 + (term.aa_low_i * cse_p528_t3_g294964) + (term.aa_low_j * cse_p528_t3_g294964) + (cse_p63_t1_g286432 * dz) + (-((term.cc_low_i * cse_p528_t7_g294968))) + (-((term.cc_low_j * cse_p528_t7_g294968))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (term.cc_low_i + term.cc_low_j + cse_p12_t27_g285626 + (term.bb_low_i * cse_p528_t4_g294965) + (term.bb_low_j * cse_p528_t4_g294965));
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= term.r2_high)) {
            /* !COND-3 */
            double cse_p1_t3_invr285549 = (1.0 / (rr));
            double cse_p1_t4_invr2285550 = (cse_p1_t3_invr285549 * cse_p1_t3_invr285549);
            double inv_rr = cse_p1_t3_invr285549;
            double inv_rr2 = (inv_rr * inv_rr);
            double di = (rr + (-(term.lj_radius_i)));
            double cse_p14_t7_g285670 = (di * di);
            double cse_p13_t12_g285652 = (-(cse_p14_t7_g285670));
            double cse_p1_t5_invr285551 = (1.0 / (term.lambda2_i));
            double cse_p12_t20_g285619 = (cse_p13_t12_g285652 * cse_p1_t5_invr285551);
            double cse_p11_t36_g285588 = exp(cse_p12_t20_g285619);
            double exp_i = cse_p11_t36_g285588;
            double cse_p12_t24_g285623 = (exp_i * term.pref_i_j);
            double fi = (exp_i * inv_rr2 * term.pref_i_j);
            double hi = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t5_invr285551 * di)));
            double dfi = (fi * hi);
            double dj = (rr + (-(term.lj_radius_j)));
            double cse_p14_t8_g285671 = (dj * dj);
            double cse_p13_t13_g285653 = (-(cse_p14_t8_g285671));
            double cse_p1_t6_invr285552 = (1.0 / (term.lambda2_j));
            double cse_p12_t21_g285620 = (cse_p13_t13_g285653 * cse_p1_t6_invr285552);
            double cse_p11_t37_g285589 = exp(cse_p12_t21_g285620);
            double exp_j = cse_p11_t37_g285589;
            double cse_p12_t25_g285624 = (exp_j * term.pref_j_i);
            double cse_p11_t25_g285577 = (cse_p12_t24_g285623 + cse_p12_t25_g285624);
            double fj = (exp_j * inv_rr2 * term.pref_j_i);
            double hj = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t6_invr285552 * dj)));
            double dfj = (fj * hj);
            double energy = (fi + fj);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p11_t1_g285591 = (-2.00000000000000000e+0 * inv_rr2);
            double cse_p528_t1_g294980 = (cse_p1_t2_invsqrt285548 * dx);
            double cse_p11_t1_g285596 = (cse_p11_t1_g285591 * cse_p11_t36_g285588 * cse_p1_t5_invr285551 * term.pref_i_j);
            double cse_p11_t1_g285597 = (cse_p11_t1_g285591 * cse_p11_t37_g285589 * cse_p1_t6_invr285552 * term.pref_j_i);
            double cse_p11_t1_g285598 = (cse_p11_t25_g285577 * cse_p1_t4_invr2285550 * inv_rr);
            double cse_p12_t1_g285636 = (-2.00000000000000000e+0 * cse_p11_t1_g285598);
            double cse_p63_t1_g286429 = (cse_p11_t1_g285596 * di);
            double cse_p64_t13_g286448 = (cse_p528_t1_g294980 * cse_p63_t1_g286429);
            double cse_p63_t1_g286430 = (cse_p11_t1_g285597 * dj);
            double cse_p64_t14_g286449 = (cse_p528_t1_g294980 * cse_p63_t1_g286430);
            double g_x1 = (cse_p64_t13_g286448 + cse_p64_t14_g286449 + (cse_p12_t1_g285636 * cse_p528_t1_g294980));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p11_t1_g285594 = (cse_p1_t2_invsqrt285548 * dy);
            double cse_p64_t5_g286440 = (cse_p11_t1_g285594 * cse_p63_t1_g286429);
            double cse_p64_t6_g286441 = (cse_p11_t1_g285594 * cse_p63_t1_g286430);
            double g_y1 = (cse_p64_t5_g286440 + cse_p64_t6_g286441 + (cse_p11_t1_g285594 * cse_p12_t1_g285636));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p11_t1_g285595 = (cse_p1_t2_invsqrt285548 * dz);
            double cse_p64_t7_g286442 = (cse_p11_t1_g285595 * cse_p63_t1_g286429);
            double cse_p64_t8_g286443 = (cse_p11_t1_g285595 * cse_p63_t1_g286430);
            double g_z1 = (cse_p64_t7_g286442 + cse_p64_t8_g286443 + (cse_p11_t1_g285595 * cse_p12_t1_g285636));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double cse_p12_t1_g285637 = (2.00000000000000000e+0 * cse_p11_t1_g285598);
            double g_x2 = ((cse_p12_t1_g285637 * cse_p528_t1_g294980) + (-(cse_p64_t13_g286448)) + (-(cse_p64_t14_g286449)));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = ((cse_p11_t1_g285594 * cse_p12_t1_g285637) + (-(cse_p64_t5_g286440)) + (-(cse_p64_t6_g286441)));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = ((cse_p11_t1_g285595 * cse_p12_t1_g285637) + (-(cse_p64_t7_g286442)) + (-(cse_p64_t8_g286443)));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (dfi + dfj);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= term.r3_high)) {
              /* !COND-4 */
              double cse_p12_t1_g285638 = (3.00000000000000000e+0 * term.aa_high_i);
              double cse_p12_t1_g285639 = (3.00000000000000000e+0 * term.aa_high_j);
              double cse_p13_t10_g285650 = (cse_p12_t1_g285638 + cse_p12_t1_g285639);
              double cse_p12_t5_g285604 = (term.bb_high_i * rr);
              double cse_p12_t6_g285605 = (term.bb_high_j * rr);
              double cse_p12_t1_g285635 = (rr * rr);
              double cse_p12_t26_g285625 = (cse_p12_t1_g285635 * cse_p13_t10_g285650);
              double cse_p11_t1_g285590 = (cse_p1_t2_invsqrt285548 * rr);
              double cse_p11_t1_g285599 = (term.bb_high_i * cse_p11_t1_g285590);
              double cse_p12_t13_g285612 = (cse_p11_t1_g285599 * dx);
              double cse_p12_t14_g285613 = (cse_p11_t1_g285599 * dy);
              double cse_p12_t15_g285614 = (cse_p11_t1_g285599 * dz);
              double cse_p11_t26_g285578 = (term.cc_high_i + term.cc_high_j + cse_p12_t26_g285625 + cse_p12_t5_g285604 + cse_p12_t5_g285604 + cse_p12_t6_g285605 + cse_p12_t6_g285605);
              double cse_p11_t28_g285580 = (-(cse_p12_t13_g285612));
              double cse_p11_t29_g285581 = (-(cse_p12_t14_g285613));
              double cse_p11_t30_g285582 = (-(cse_p12_t15_g285614));
              double energy = (term.dd_high_i + term.dd_high_j + (term.cc_high_i * rr) + (term.cc_high_j * rr) + (cse_p12_t1_g285635 * (term.bb_high_i + term.bb_high_j + (term.aa_high_i * rr) + (term.aa_high_j * rr))));
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p528_t1_g294980 = (cse_p1_t2_invsqrt285548 * dx);
              double g_x1 = (cse_p11_t26_g285578 * cse_p528_t1_g294980);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double cse_p11_t1_g285594 = (cse_p1_t2_invsqrt285548 * dy);
              double g_y1 = (cse_p11_t1_g285594 * cse_p11_t26_g285578);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double cse_p11_t1_g285595 = (cse_p1_t2_invsqrt285548 * dz);
              double g_z1 = (cse_p11_t1_g285595 * cse_p11_t26_g285578);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p11_t1_g285592 = (cse_p11_t1_g285590 * rr);
              double cse_p12_t1_g285632 = (cse_p11_t1_g285592 * dx);
              double cse_p13_t1_g285654 = (term.bb_high_j * cse_p11_t1_g285590);
              double cse_p63_t1_g286426 = (-3.00000000000000000e+0 * cse_p12_t1_g285632);
              double cse_p63_t1_g286431 = (-2.00000000000000000e+0 * cse_p13_t1_g285654);
              double g_x2 = (cse_p11_t28_g285580 + cse_p11_t28_g285580 + (term.aa_high_i * cse_p63_t1_g286426) + (term.aa_high_j * cse_p63_t1_g286426) + (cse_p63_t1_g286431 * dx) + (-((term.cc_high_i * cse_p528_t1_g294980))) + (-((term.cc_high_j * cse_p528_t1_g294980))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double cse_p12_t1_g285633 = (cse_p11_t1_g285592 * dy);
              double cse_p63_t1_g286427 = (-3.00000000000000000e+0 * cse_p12_t1_g285633);
              double g_y2 = (cse_p11_t29_g285581 + cse_p11_t29_g285581 + (term.aa_high_i * cse_p63_t1_g286427) + (term.aa_high_j * cse_p63_t1_g286427) + (cse_p63_t1_g286431 * dy) + (-((term.cc_high_i * cse_p11_t1_g285594))) + (-((term.cc_high_j * cse_p11_t1_g285594))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double cse_p12_t1_g285634 = (cse_p11_t1_g285592 * dz);
              double cse_p63_t1_g286428 = (-3.00000000000000000e+0 * cse_p12_t1_g285634);
              double g_z2 = (cse_p11_t30_g285582 + cse_p11_t30_g285582 + (term.aa_high_i * cse_p63_t1_g286428) + (term.aa_high_j * cse_p63_t1_g286428) + (cse_p63_t1_g286431 * dz) + (-((term.cc_high_i * cse_p11_t1_g285595))) + (-((term.cc_high_j * cse_p11_t1_g285595))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double cse_p12_t1_g285631 = (2.00000000000000000e+0 * rr);
              double de_dr = (term.cc_high_i + term.cc_high_j + cse_p12_t26_g285625 + (term.bb_high_i * cse_p12_t1_g285631) + (term.bb_high_j * cse_p12_t1_g285631));
            } else {
              /* !COND-ELSE-4 */
              {
                /* !COND-ELSE-5 */
                double energy = 0.00000000000000000e+0;
                energy_added += energy;
                *energy_accumulate += energy;
                double g_x1 = 0.00000000000000000e+0;
                KernelGradientAcc(term.i3x1, 0, g_x1);
                double g_y1 = 0.00000000000000000e+0;
                KernelGradientAcc(term.i3x1, 1, g_y1);
                double g_z1 = 0.00000000000000000e+0;
                KernelGradientAcc(term.i3x1, 2, g_z1);
                double g_x2 = 0.00000000000000000e+0;
                KernelGradientAcc(term.i3x2, 0, g_x2);
                double g_y2 = 0.00000000000000000e+0;
                KernelGradientAcc(term.i3x2, 1, g_y2);
                double g_z2 = 0.00000000000000000e+0;
                KernelGradientAcc(term.i3x2, 2, g_z2);
                double de_dr = 0.00000000000000000e+0;
              }
            }
          }
        }
      }
    }
  }
  return energy_added;
}
void gradient_fd(const rosetta_lk_solvation_parameters& params, const rosetta_lk_solvation_term& term, double* position, double* energy_accumulate, double* force)
{
  constexpr size_t PositionSize = 6;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  double e0 = 0.0;
  energy(params, term, position, &e0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[term.i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 0] = saved + h;
    energy(params, term, position, &e_plus);
    position[term.i3x1 + 0] = saved - h;
    energy(params, term, position, &e_minus);
    position[term.i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x1, 0, d);
  }
  {
    double saved = position[term.i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 1] = saved + h;
    energy(params, term, position, &e_plus);
    position[term.i3x1 + 1] = saved - h;
    energy(params, term, position, &e_minus);
    position[term.i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x1, 1, d);
  }
  {
    double saved = position[term.i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 2] = saved + h;
    energy(params, term, position, &e_plus);
    position[term.i3x1 + 2] = saved - h;
    energy(params, term, position, &e_minus);
    position[term.i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x1, 2, d);
  }
  {
    double saved = position[term.i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 0] = saved + h;
    energy(params, term, position, &e_plus);
    position[term.i3x2 + 0] = saved - h;
    energy(params, term, position, &e_minus);
    position[term.i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x2, 0, d);
  }
  {
    double saved = position[term.i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 1] = saved + h;
    energy(params, term, position, &e_plus);
    position[term.i3x2 + 1] = saved - h;
    energy(params, term, position, &e_minus);
    position[term.i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x2, 1, d);
  }
  {
    double saved = position[term.i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 2] = saved + h;
    energy(params, term, position, &e_plus);
    position[term.i3x2 + 2] = saved - h;
    energy(params, term, position, &e_minus);
    position[term.i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x2, 2, d);
  }
}

double energy(const rosetta_lk_solvation_parameters& params, const rosetta_lk_solvation_term& term, double* position, double* energy_accumulate) {
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
      double rr = sqrt(r2);
      if ((rr <= term.r0_low)) {
        /* !COND-1 */
        double energy = (term.f_contact_i + term.f_contact_j);
        energy_added += energy;
        *energy_accumulate += energy;
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.r1_low)) {
          /* !COND-2 */
          double energy = (term.dd_low_i + term.dd_low_j + (term.cc_low_i * rr) + (term.cc_low_j * rr) + (rr * rr * (term.bb_low_i + term.bb_low_j + (term.aa_low_i * rr) + (term.aa_low_j * rr))));
          energy_added += energy;
          *energy_accumulate += energy;
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= term.r2_high)) {
            /* !COND-3 */
            double inv_rr = (1.0 / (rr));
            double inv_rr2 = (inv_rr * inv_rr);
            double di = (rr + (-(term.lj_radius_i)));
            double cse_p61_t1_invr295533 = (1.0 / (term.lambda2_i));
            double exp_i = exp((-((cse_p61_t1_invr295533 * di * di))));
            double fi = (exp_i * inv_rr2 * term.pref_i_j);
            double hi = (-2.00000000000000000e+0 * (inv_rr + (cse_p61_t1_invr295533 * di)));
            double dj = (rr + (-(term.lj_radius_j)));
            double cse_p61_t2_invr295534 = (1.0 / (term.lambda2_j));
            double exp_j = exp((-((cse_p61_t2_invr295534 * dj * dj))));
            double fj = (exp_j * inv_rr2 * term.pref_j_i);
            double hj = (-2.00000000000000000e+0 * (inv_rr + (cse_p61_t2_invr295534 * dj)));
            double energy = (fi + fj);
            energy_added += energy;
            *energy_accumulate += energy;
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= term.r3_high)) {
              /* !COND-4 */
              double energy = (term.dd_high_i + term.dd_high_j + (term.cc_high_i * rr) + (term.cc_high_j * rr) + (rr * rr * (term.bb_high_i + term.bb_high_j + (term.aa_high_i * rr) + (term.aa_high_j * rr))));
              energy_added += energy;
              *energy_accumulate += energy;
            } else {
              /* !COND-ELSE-4 */
              {
                /* !COND-ELSE-5 */
                double energy = 0.00000000000000000e+0;
                energy_added += energy;
                *energy_accumulate += energy;
              }
            }
          }
        }
      }
    }
  }
  return energy_added;
}
void energy_fd(const rosetta_lk_solvation_parameters& params, const rosetta_lk_solvation_term& term, double* position, double* energy_accumulate)
{
  energy(params, term, position, energy_accumulate);
}

double hessian(const rosetta_lk_solvation_parameters& params, const rosetta_lk_solvation_term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
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
      double cse_p11_t53_g300543 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t54_g300544 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t55_g300545 = (dz * dz);
      double r2 = (cse_p11_t53_g300543 + cse_p11_t54_g300544 + cse_p11_t55_g300545);
      double cse_p1_t1_sqrt300485 = sqrt(r2);
      double cse_p1_t2_invsqrt300486 = (1.0 / (cse_p1_t1_sqrt300485));
      double rr = cse_p1_t1_sqrt300485;
      if ((rr <= term.r0_low)) {
        /* !COND-1 */
        double cse_p60_t43_g304592 = (rr * rr);
        double cse_p60_t71_g304620 = (-(cse_p11_t53_g300543));
        double cse_p60_t72_g304621 = (-(cse_p11_t54_g300544));
        double cse_p60_t73_g304622 = (-(cse_p11_t55_g300545));
        double cse_p60_t75_g304624 = (-(dx));
        double cse_p60_t76_g304625 = (-(dy));
        double cse_p60_t77_g304626 = (-(dz));
        double cse_p565_t1_invr340180 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr340181 = (1.0 / (rr));
        double cse_p60_t80_g304629 = (rr * rr * rr);
        double cse_p15_t75_g301023 = (-(cse_p60_t43_g304592));
        double cse_p14_t82_g300921 = (cse_p11_t53_g300543 + cse_p15_t75_g301023);
        double cse_p14_t84_g300923 = (cse_p11_t54_g300544 + cse_p15_t75_g301023);
        double cse_p14_t86_g300925 = (cse_p11_t55_g300545 + cse_p15_t75_g301023);
        double cse_p13_t109_g300824 = (cse_p60_t43_g304592 + cse_p60_t71_g304620);
        double cse_p13_t110_g300825 = (cse_p60_t43_g304592 + cse_p60_t72_g304621);
        double cse_p13_t111_g300826 = (cse_p60_t43_g304592 + cse_p60_t73_g304622);
        double energy = (term.f_contact_i + term.f_contact_j);
        energy_added += energy;
        *energy_accumulate += energy;
        double g_x1 = 0.00000000000000000e+0;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = 0.00000000000000000e+0;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = 0.00000000000000000e+0;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = 0.00000000000000000e+0;
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = 0.00000000000000000e+0;
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = 0.00000000000000000e+0;
        KernelGradientAcc(term.i3x2, 2, g_z2);
        double de_dr = 0.00000000000000000e+0;
        double cse_p60_t42_g304591 = (cse_p565_t1_invr340180 * de_dr);
        double cse_p13_t67_g300782 = (cse_p14_t82_g300921 * cse_p60_t42_g304591);
        double cse_p13_t68_g300783 = (cse_p14_t84_g300923 * cse_p60_t42_g304591);
        double cse_p13_t69_g300784 = (cse_p14_t86_g300925 * cse_p60_t42_g304591);
        double cse_p13_t70_g300785 = (cse_p60_t42_g304591 * dx * dy);
        double cse_p13_t71_g300786 = (cse_p60_t42_g304591 * cse_p60_t75_g304624 * dy);
        double cse_p12_t27_g300639 = (cse_p13_t109_g300824 * cse_p60_t42_g304591);
        double cse_p12_t28_g300640 = (cse_p13_t110_g300825 * cse_p60_t42_g304591);
        double cse_p12_t29_g300641 = (cse_p13_t111_g300826 * cse_p60_t42_g304591);
        double cse_p12_t32_g300644 = (cse_p60_t42_g304591 * dz);
        double cse_p13_t57_g300772 = (cse_p12_t32_g300644 * cse_p60_t75_g304624);
        double cse_p13_t58_g300773 = (cse_p12_t32_g300644 * cse_p60_t76_g304625);
        double cse_p13_t59_g300774 = (cse_p12_t32_g300644 * dx);
        double cse_p13_t60_g300775 = (cse_p12_t32_g300644 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p60_t41_g304590 = (cse_p565_t2_invr340181 * cse_p565_t2_invr340181 * d2e_dr2);
        double cse_p13_t61_g300776 = (cse_p60_t41_g304590 * dx * dx);
        double cse_p13_t62_g300777 = (cse_p60_t41_g304590 * dx * dy);
        double cse_p13_t63_g300778 = (cse_p60_t41_g304590 * dx * dz);
        double cse_p13_t64_g300779 = (cse_p60_t41_g304590 * dy * dy);
        double cse_p13_t65_g300780 = (cse_p60_t41_g304590 * dy * dz);
        double cse_p13_t66_g300781 = (cse_p60_t41_g304590 * dz * dz);
        double cse_p12_t15_g300627 = (cse_p60_t41_g304590 * cse_p60_t75_g304624);
        double cse_p13_t42_g300757 = (cse_p12_t15_g300627 * cse_p60_t75_g304624);
        double cse_p13_t43_g300758 = (cse_p12_t15_g300627 * cse_p60_t76_g304625);
        double cse_p13_t44_g300759 = (cse_p12_t15_g300627 * cse_p60_t77_g304626);
        double cse_p13_t45_g300760 = (cse_p12_t15_g300627 * dx);
        double cse_p13_t46_g300761 = (cse_p12_t15_g300627 * dy);
        double cse_p13_t47_g300762 = (cse_p12_t15_g300627 * dz);
        double cse_p12_t16_g300628 = (cse_p60_t41_g304590 * cse_p60_t76_g304625);
        double cse_p13_t48_g300763 = (cse_p12_t16_g300628 * cse_p60_t76_g304625);
        double cse_p13_t49_g300764 = (cse_p12_t16_g300628 * cse_p60_t77_g304626);
        double cse_p13_t50_g300765 = (cse_p12_t16_g300628 * dx);
        double cse_p13_t51_g300766 = (cse_p12_t16_g300628 * dy);
        double cse_p13_t52_g300767 = (cse_p12_t16_g300628 * dz);
        double cse_p12_t17_g300629 = (cse_p60_t41_g304590 * cse_p60_t77_g304626);
        double cse_p13_t53_g300768 = (cse_p12_t17_g300629 * cse_p60_t77_g304626);
        double cse_p13_t54_g300769 = (cse_p12_t17_g300629 * dx);
        double cse_p13_t55_g300770 = (cse_p12_t17_g300629 * dy);
        double cse_p13_t56_g300771 = (cse_p12_t17_g300629 * dz);
        double cse_p12_t68_g300680 = (cse_p13_t62_g300777 + cse_p13_t71_g300786);
        double cse_p12_t69_g300681 = (cse_p13_t57_g300772 + cse_p13_t63_g300778);
        double cse_p12_t70_g300682 = (cse_p13_t58_g300773 + cse_p13_t65_g300780);
        double cse_p12_t86_g300698 = (cse_p12_t28_g300640 + cse_p13_t64_g300779);
        double cse_p12_t88_g300700 = (cse_p12_t29_g300641 + cse_p13_t66_g300781);
        double cse_p12_t90_g300702 = (cse_p12_t27_g300639 + cse_p13_t61_g300776);
        double cse_p12_t63_g300675 = (cse_p13_t45_g300760 + cse_p13_t67_g300782);
        double cse_p12_t71_g300683 = (cse_p13_t46_g300761 + cse_p13_t70_g300785);
        double cse_p12_t73_g300685 = (cse_p13_t43_g300758 + cse_p13_t71_g300786);
        double cse_p12_t74_g300686 = (cse_p13_t44_g300759 + cse_p13_t57_g300772);
        double cse_p12_t75_g300687 = (cse_p13_t47_g300762 + cse_p13_t59_g300774);
        double cse_p12_t91_g300703 = (cse_p12_t27_g300639 + cse_p13_t42_g300757);
        double cse_p12_t64_g300676 = (cse_p13_t51_g300766 + cse_p13_t68_g300783);
        double cse_p12_t72_g300684 = (cse_p13_t50_g300765 + cse_p13_t70_g300785);
        double cse_p12_t76_g300688 = (cse_p13_t49_g300764 + cse_p13_t58_g300773);
        double cse_p12_t77_g300689 = (cse_p13_t52_g300767 + cse_p13_t60_g300775);
        double cse_p12_t87_g300699 = (cse_p12_t28_g300640 + cse_p13_t48_g300763);
        double cse_p12_t65_g300677 = (cse_p13_t56_g300771 + cse_p13_t69_g300784);
        double cse_p12_t78_g300690 = (cse_p13_t54_g300769 + cse_p13_t59_g300774);
        double cse_p12_t79_g300691 = (cse_p13_t55_g300770 + cse_p13_t60_g300775);
        double cse_p12_t89_g300701 = (cse_p12_t29_g300641 + cse_p13_t53_g300768);
        double h_x1_x1 = cse_p12_t90_g300702;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t68_g300680;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t69_g300681;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t63_g300675;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t72_g300684;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t78_g300690;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t86_g300698;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t70_g300682;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t71_g300683;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t64_g300676;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t79_g300691;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t88_g300700;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t75_g300687;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t77_g300689;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t65_g300677;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t91_g300703;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t73_g300685;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t74_g300686;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t87_g300699;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t76_g300688;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t89_g300701;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.r1_low)) {
          /* !COND-2 */
          double cse_p60_t38_g304587 = (cse_p1_t2_invsqrt300486 * dx);
          double cse_p60_t39_g304588 = (cse_p1_t2_invsqrt300486 * dy);
          double cse_p60_t40_g304589 = (cse_p1_t2_invsqrt300486 * dz);
          double cse_p14_t3_g300842 = (3.00000000000000000e+0 * term.aa_low_i);
          double cse_p14_t4_g300843 = (3.00000000000000000e+0 * term.aa_low_j);
          double cse_p13_t81_g300796 = (cse_p14_t3_g300842 + cse_p14_t4_g300843);
          double cse_p12_t10_g300622 = (term.bb_low_i * cse_p60_t38_g304587 * rr);
          double cse_p12_t11_g300623 = (term.bb_low_i * cse_p60_t39_g304588 * rr);
          double cse_p12_t12_g300624 = (term.bb_low_i * cse_p60_t40_g304589 * rr);
          double cse_p12_t13_g300625 = (term.bb_low_i * rr);
          double cse_p12_t14_g300626 = (term.bb_low_j * rr);
          double cse_p11_t1_g300605 = (rr * rr);
          double cse_p60_t59_g304608 = (cse_p11_t53_g300543 + (-(cse_p11_t1_g300605)));
          double cse_p60_t60_g304609 = (cse_p11_t54_g300544 + (-(cse_p11_t1_g300605)));
          double cse_p60_t61_g304610 = (cse_p11_t55_g300545 + (-(cse_p11_t1_g300605)));
          double cse_p60_t62_g304611 = (cse_p11_t1_g300605 + (-(cse_p11_t53_g300543)));
          double cse_p60_t63_g304612 = (cse_p11_t1_g300605 + (-(cse_p11_t54_g300544)));
          double cse_p60_t64_g304613 = (cse_p11_t1_g300605 + (-(cse_p11_t55_g300545)));
          double cse_p60_t74_g304623 = (-(cse_p11_t1_g300605));
          double cse_p12_t34_g300646 = (cse_p11_t1_g300605 * cse_p13_t81_g300796);
          double cse_p11_t92_g300582 = (term.cc_low_i + term.cc_low_j + cse_p12_t13_g300625 + cse_p12_t13_g300625 + cse_p12_t14_g300626 + cse_p12_t14_g300626 + cse_p12_t34_g300646);
          double cse_p11_t96_g300586 = (-(cse_p12_t10_g300622));
          double cse_p11_t97_g300587 = (-(cse_p12_t11_g300623));
          double cse_p11_t98_g300588 = (-(cse_p12_t12_g300624));
          double energy = (term.dd_low_i + term.dd_low_j + (term.cc_low_i * rr) + (term.cc_low_j * rr) + (cse_p11_t1_g300605 * (term.bb_low_i + term.bb_low_j + (term.aa_low_i * rr) + (term.aa_low_j * rr))));
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (cse_p11_t92_g300582 * cse_p60_t38_g304587);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t92_g300582 * cse_p60_t39_g304588);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t92_g300582 * cse_p60_t40_g304589);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p11_t96_g300586 + cse_p11_t96_g300586 + (-2.00000000000000000e+0 * term.bb_low_j * cse_p60_t38_g304587 * rr) + (cse_p11_t1_g300605 * cse_p60_t38_g304587 * ((-3.00000000000000000e+0 * term.aa_low_i) + (-3.00000000000000000e+0 * term.aa_low_j))) + (-((term.cc_low_i * cse_p60_t38_g304587))) + (-((term.cc_low_j * cse_p60_t38_g304587))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p11_t97_g300587 + cse_p11_t97_g300587 + (-2.00000000000000000e+0 * term.bb_low_j * cse_p60_t39_g304588 * rr) + (cse_p11_t1_g300605 * cse_p60_t39_g304588 * ((-3.00000000000000000e+0 * term.aa_low_i) + (-3.00000000000000000e+0 * term.aa_low_j))) + (-((term.cc_low_i * cse_p60_t39_g304588))) + (-((term.cc_low_j * cse_p60_t39_g304588))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p11_t98_g300588 + cse_p11_t98_g300588 + (-2.00000000000000000e+0 * term.bb_low_j * cse_p60_t40_g304589 * rr) + (cse_p11_t1_g300605 * cse_p60_t40_g304589 * ((-3.00000000000000000e+0 * term.aa_low_i) + (-3.00000000000000000e+0 * term.aa_low_j))) + (-((term.cc_low_i * cse_p60_t40_g304589))) + (-((term.cc_low_j * cse_p60_t40_g304589))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (term.cc_low_i + term.cc_low_j + cse_p12_t34_g300646 + (2.00000000000000000e+0 * term.bb_low_i * rr) + (2.00000000000000000e+0 * term.bb_low_j * rr));
          double d2e_dr2 = ((2.00000000000000000e+0 * term.bb_low_i) + (2.00000000000000000e+0 * term.bb_low_j) + (6.00000000000000000e+0 * term.aa_low_i * rr) + (6.00000000000000000e+0 * term.aa_low_j * rr));
          double cse_p565_t2_invr340181 = (1.0 / (rr));
          double cse_p11_t1_g300603 = (cse_p565_t2_invr340181 * cse_p565_t2_invr340181 * d2e_dr2);
          double cse_p60_t1_g304550 = (cse_p11_t1_g300603 * (-(dx)));
          double cse_p60_t2_g304551 = (cse_p11_t1_g300603 * (-(dy)));
          double cse_p60_t3_g304552 = (cse_p11_t1_g300603 * (-(dz)));
          double cse_p60_t4_g304553 = (cse_p11_t1_g300603 * dx * dx);
          double cse_p60_t5_g304554 = (cse_p11_t1_g300603 * dx * dy);
          double cse_p60_t6_g304555 = (cse_p11_t1_g300603 * dx * dz);
          double cse_p60_t7_g304556 = (cse_p11_t1_g300603 * dy * dy);
          double cse_p60_t8_g304557 = (cse_p11_t1_g300603 * dy * dz);
          double cse_p60_t9_g304558 = (cse_p11_t1_g300603 * dz * dz);
          double cse_p14_t46_g300885 = (cse_p60_t1_g304550 * (-(dx)));
          double cse_p14_t47_g300886 = (cse_p60_t1_g304550 * (-(dy)));
          double cse_p14_t48_g300887 = (cse_p60_t1_g304550 * (-(dz)));
          double cse_p14_t49_g300888 = (cse_p60_t1_g304550 * dx);
          double cse_p14_t50_g300889 = (cse_p60_t1_g304550 * dy);
          double cse_p14_t51_g300890 = (cse_p60_t1_g304550 * dz);
          double cse_p14_t52_g300891 = (cse_p60_t2_g304551 * (-(dy)));
          double cse_p14_t53_g300892 = (cse_p60_t2_g304551 * (-(dz)));
          double cse_p14_t54_g300893 = (cse_p60_t2_g304551 * dx);
          double cse_p14_t55_g300894 = (cse_p60_t2_g304551 * dy);
          double cse_p14_t56_g300895 = (cse_p60_t2_g304551 * dz);
          double cse_p14_t57_g300896 = (cse_p60_t3_g304552 * (-(dz)));
          double cse_p14_t58_g300897 = (cse_p60_t3_g304552 * dx);
          double cse_p14_t59_g300898 = (cse_p60_t3_g304552 * dy);
          double cse_p14_t60_g300899 = (cse_p60_t3_g304552 * dz);
          double cse_p565_t1_invr340180 = (1.0 / ((rr * rr * rr)));
          double cse_p11_t1_g300604 = (cse_p565_t1_invr340180 * de_dr);
          double cse_p60_t10_g304559 = (cse_p11_t1_g300604 * (cse_p11_t53_g300543 + cse_p60_t74_g304623));
          double cse_p60_t11_g304560 = (cse_p11_t1_g300604 * (cse_p11_t54_g300544 + cse_p60_t74_g304623));
          double cse_p60_t12_g304561 = (cse_p11_t1_g300604 * (cse_p11_t55_g300545 + cse_p60_t74_g304623));
          double cse_p60_t13_g304562 = (cse_p11_t1_g300604 * cse_p60_t62_g304611);
          double cse_p60_t14_g304563 = (cse_p11_t1_g300604 * cse_p60_t63_g304612);
          double cse_p60_t15_g304564 = (cse_p11_t1_g300604 * cse_p60_t64_g304613);
          double cse_p60_t16_g304565 = (cse_p11_t1_g300604 * dx * dy);
          double cse_p60_t17_g304566 = (cse_p11_t1_g300604 * dy * (-(dx)));
          double cse_p60_t18_g304567 = (cse_p11_t1_g300604 * dz);
          double cse_p60_t47_g304596 = (cse_p60_t17_g304566 + cse_p60_t5_g304554);
          double cse_p14_t42_g300881 = (cse_p60_t18_g304567 * (-(dx)));
          double cse_p14_t43_g300882 = (cse_p60_t18_g304567 * (-(dy)));
          double cse_p14_t44_g300883 = (cse_p60_t18_g304567 * dx);
          double cse_p14_t45_g300884 = (cse_p60_t18_g304567 * dy);
          double cse_p13_t107_g300822 = (cse_p60_t13_g304562 + cse_p60_t4_g304553);
          double h_x1_x1 = cse_p13_t107_g300822;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t47_g304596;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t83_g300798 = (cse_p14_t42_g300881 + cse_p60_t6_g304555);
          double cse_p13_t84_g300799 = (cse_p14_t43_g300882 + cse_p60_t8_g304557);
          double h_x1_z1 = cse_p13_t83_g300798;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t77_g300792 = (cse_p14_t49_g300888 + cse_p60_t10_g304559);
          double cse_p13_t85_g300800 = (cse_p14_t50_g300889 + cse_p60_t16_g304565);
          double cse_p13_t87_g300802 = (cse_p14_t47_g300886 + cse_p60_t17_g304566);
          double cse_p13_t88_g300803 = (cse_p14_t42_g300881 + cse_p14_t48_g300887);
          double cse_p13_t89_g300804 = (cse_p14_t44_g300883 + cse_p14_t51_g300890);
          double cse_p13_t108_g300823 = (cse_p14_t46_g300885 + cse_p60_t13_g304562);
          double h_x1_x2 = cse_p13_t77_g300792;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t78_g300793 = (cse_p14_t55_g300894 + cse_p60_t11_g304560);
          double cse_p13_t86_g300801 = (cse_p14_t54_g300893 + cse_p60_t16_g304565);
          double cse_p13_t90_g300805 = (cse_p14_t43_g300882 + cse_p14_t53_g300892);
          double cse_p13_t91_g300806 = (cse_p14_t45_g300884 + cse_p14_t56_g300895);
          double h_x1_y2 = cse_p13_t86_g300801;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t79_g300794 = (cse_p14_t60_g300899 + cse_p60_t12_g304561);
          double cse_p13_t92_g300807 = (cse_p14_t44_g300883 + cse_p14_t58_g300897);
          double cse_p13_t93_g300808 = (cse_p14_t45_g300884 + cse_p14_t59_g300898);
          double h_x1_z2 = cse_p13_t92_g300807;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t103_g300818 = (cse_p60_t14_g304563 + cse_p60_t7_g304556);
          double cse_p13_t104_g300819 = (cse_p14_t52_g300891 + cse_p60_t14_g304563);
          double h_y1_y1 = cse_p13_t103_g300818;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t84_g300799;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t85_g300800;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t78_g300793;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t93_g300808;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t105_g300820 = (cse_p60_t15_g304564 + cse_p60_t9_g304558);
          double cse_p13_t106_g300821 = (cse_p14_t57_g300896 + cse_p60_t15_g304564);
          double h_z1_z1 = cse_p13_t105_g300820;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t89_g300804;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t91_g300806;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t79_g300794;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t108_g300823;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t87_g300802;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t88_g300803;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t104_g300819;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t90_g300805;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t106_g300821;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= term.r2_high)) {
            /* !COND-3 */
            double cse_p565_t2_invr340181 = (1.0 / (rr));
            double cse_p1_t4_invr2300488 = (cse_p565_t2_invr340181 * cse_p565_t2_invr340181);
            double inv_rr = cse_p565_t2_invr340181;
            double inv_rr2 = (inv_rr * inv_rr);
            double cse_p11_t1_g300491 = (2.00000000000000000e+0 * inv_rr2);
            double di = (rr + (-(term.lj_radius_i)));
            double cse_p14_t63_g300902 = (di * di);
            double cse_p13_t112_g300827 = (-(cse_p14_t63_g300902));
            double cse_p1_t5_invr300489 = (1.0 / (term.lambda2_i));
            double cse_p12_t54_g300666 = (cse_p13_t112_g300827 * cse_p1_t5_invr300489);
            double cse_p11_t108_g300598 = exp(cse_p12_t54_g300666);
            double exp_i = cse_p11_t108_g300598;
            double cse_p12_t60_g300672 = (exp_i * term.pref_i_j);
            double fi = (exp_i * inv_rr2 * term.pref_i_j);
            double hi = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t5_invr300489 * di)));
            double dfi = (fi * hi);
            double d2fi = (fi * (cse_p11_t1_g300491 + (-2.00000000000000000e+0 * cse_p1_t5_invr300489) + (hi * hi)));
            double dj = (rr + (-(term.lj_radius_j)));
            double cse_p14_t64_g300903 = (dj * dj);
            double cse_p13_t113_g300828 = (-(cse_p14_t64_g300903));
            double cse_p1_t6_invr300490 = (1.0 / (term.lambda2_j));
            double cse_p12_t55_g300667 = (cse_p13_t113_g300828 * cse_p1_t6_invr300490);
            double cse_p11_t109_g300599 = exp(cse_p12_t55_g300667);
            double exp_j = cse_p11_t109_g300599;
            double cse_p12_t61_g300673 = (exp_j * term.pref_j_i);
            double cse_p11_t87_g300577 = (cse_p12_t60_g300672 + cse_p12_t61_g300673);
            double fj = (exp_j * inv_rr2 * term.pref_j_i);
            double hj = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t6_invr300490 * dj)));
            double dfj = (fj * hj);
            double d2fj = (fj * (cse_p11_t1_g300491 + (-2.00000000000000000e+0 * cse_p1_t6_invr300490) + (hj * hj)));
            double energy = (fi + fj);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p14_t1_g300946 = (cse_p1_t2_invsqrt300486 * dx);
            double g_x1 = ((-2.00000000000000000e+0 * cse_p11_t87_g300577 * cse_p14_t1_g300946 * cse_p1_t4_invr2300488 * inv_rr) + (cse_p14_t1_g300946 * inv_rr2 * ((-2.00000000000000000e+0 * cse_p11_t108_g300598 * cse_p1_t5_invr300489 * di * term.pref_i_j) + (-2.00000000000000000e+0 * cse_p11_t109_g300599 * cse_p1_t6_invr300490 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p14_t1_g300947 = (cse_p1_t2_invsqrt300486 * dy);
            double g_y1 = ((-2.00000000000000000e+0 * cse_p11_t87_g300577 * cse_p14_t1_g300947 * cse_p1_t4_invr2300488 * inv_rr) + (cse_p14_t1_g300947 * inv_rr2 * ((-2.00000000000000000e+0 * cse_p11_t108_g300598 * cse_p1_t5_invr300489 * di * term.pref_i_j) + (-2.00000000000000000e+0 * cse_p11_t109_g300599 * cse_p1_t6_invr300490 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p14_t1_g300948 = (cse_p1_t2_invsqrt300486 * dz);
            double g_z1 = ((-2.00000000000000000e+0 * cse_p11_t87_g300577 * cse_p14_t1_g300948 * cse_p1_t4_invr2300488 * inv_rr) + (cse_p14_t1_g300948 * inv_rr2 * ((-2.00000000000000000e+0 * cse_p11_t108_g300598 * cse_p1_t5_invr300489 * di * term.pref_i_j) + (-2.00000000000000000e+0 * cse_p11_t109_g300599 * cse_p1_t6_invr300490 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = ((2.00000000000000000e+0 * cse_p11_t87_g300577 * cse_p14_t1_g300946 * cse_p1_t4_invr2300488 * inv_rr) + (cse_p14_t1_g300946 * inv_rr2 * ((2.00000000000000000e+0 * cse_p11_t108_g300598 * cse_p1_t5_invr300489 * di * term.pref_i_j) + (2.00000000000000000e+0 * cse_p11_t109_g300599 * cse_p1_t6_invr300490 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = ((2.00000000000000000e+0 * cse_p11_t87_g300577 * cse_p14_t1_g300947 * cse_p1_t4_invr2300488 * inv_rr) + (cse_p14_t1_g300947 * inv_rr2 * ((2.00000000000000000e+0 * cse_p11_t108_g300598 * cse_p1_t5_invr300489 * di * term.pref_i_j) + (2.00000000000000000e+0 * cse_p11_t109_g300599 * cse_p1_t6_invr300490 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = ((2.00000000000000000e+0 * cse_p11_t87_g300577 * cse_p14_t1_g300948 * cse_p1_t4_invr2300488 * inv_rr) + (cse_p14_t1_g300948 * inv_rr2 * ((2.00000000000000000e+0 * cse_p11_t108_g300598 * cse_p1_t5_invr300489 * di * term.pref_i_j) + (2.00000000000000000e+0 * cse_p11_t109_g300599 * cse_p1_t6_invr300490 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (dfi + dfj);
            double d2e_dr2 = (d2fi + d2fj);
            double cse_p11_t1_g300603 = (cse_p565_t2_invr340181 * cse_p565_t2_invr340181 * d2e_dr2);
            double cse_p565_t1_invr340180 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g300604 = (cse_p565_t1_invr340180 * de_dr);
            double cse_p11_t1_g300605 = (rr * rr);
            double cse_p11_t1_g300612 = (cse_p11_t1_g300604 * (cse_p11_t1_g300605 + (-(cse_p11_t53_g300543))));
            double cse_p60_t69_g304618 = (cse_p11_t1_g300612 + (cse_p11_t1_g300603 * dx * dx));
            double h_x1_x1 = cse_p60_t69_g304618;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = ((cse_p11_t1_g300603 * dx * dy) + (cse_p11_t1_g300604 * dy * (-(dx))));
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p11_t1_g300608 = (cse_p11_t1_g300604 * dz);
            double cse_p60_t30_g304579 = (cse_p11_t1_g300608 * (-(dx)));
            double cse_p60_t31_g304580 = (cse_p11_t1_g300608 * (-(dy)));
            double cse_p60_t32_g304581 = (cse_p11_t1_g300608 * dx);
            double cse_p60_t33_g304582 = (cse_p11_t1_g300608 * dy);
            double cse_p60_t48_g304597 = (cse_p60_t30_g304579 + (cse_p11_t1_g300603 * dx * dz));
            double cse_p60_t49_g304598 = (cse_p60_t31_g304580 + (cse_p11_t1_g300603 * dy * dz));
            double h_x1_z1 = cse_p60_t48_g304597;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p11_t1_g300606 = (cse_p11_t1_g300603 * (-(dx)));
            double cse_p60_t19_g304568 = (cse_p11_t1_g300606 * (-(dx)));
            double cse_p60_t20_g304569 = (cse_p11_t1_g300606 * (-(dy)));
            double cse_p60_t21_g304570 = (cse_p11_t1_g300606 * (-(dz)));
            double cse_p60_t22_g304571 = (cse_p11_t1_g300606 * dx);
            double cse_p60_t23_g304572 = (cse_p11_t1_g300606 * dy);
            double cse_p60_t24_g304573 = (cse_p11_t1_g300606 * dz);
            double cse_p60_t44_g304593 = (cse_p60_t22_g304571 + (cse_p11_t1_g300604 * (cse_p11_t53_g300543 + (-(cse_p11_t1_g300605)))));
            double cse_p60_t50_g304599 = (cse_p60_t23_g304572 + (cse_p11_t1_g300604 * dx * dy));
            double cse_p60_t52_g304601 = (cse_p60_t20_g304569 + (cse_p11_t1_g300604 * dy * (-(dx))));
            double cse_p60_t53_g304602 = (cse_p60_t21_g304570 + cse_p60_t30_g304579);
            double cse_p60_t54_g304603 = (cse_p60_t24_g304573 + cse_p60_t32_g304581);
            double cse_p60_t70_g304619 = (cse_p11_t1_g300612 + cse_p60_t19_g304568);
            double h_x1_x2 = cse_p60_t44_g304593;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p11_t1_g300607 = (cse_p11_t1_g300603 * (-(dy)));
            double cse_p60_t25_g304574 = (cse_p11_t1_g300607 * (-(dy)));
            double cse_p60_t26_g304575 = (cse_p11_t1_g300607 * (-(dz)));
            double cse_p60_t27_g304576 = (cse_p11_t1_g300607 * dx);
            double cse_p60_t28_g304577 = (cse_p11_t1_g300607 * dy);
            double cse_p60_t29_g304578 = (cse_p11_t1_g300607 * dz);
            double cse_p60_t45_g304594 = (cse_p60_t28_g304577 + (cse_p11_t1_g300604 * (cse_p11_t54_g300544 + (-(cse_p11_t1_g300605)))));
            double cse_p60_t51_g304600 = (cse_p60_t27_g304576 + (cse_p11_t1_g300604 * dx * dy));
            double cse_p60_t55_g304604 = (cse_p60_t26_g304575 + cse_p60_t31_g304580);
            double cse_p60_t56_g304605 = (cse_p60_t29_g304578 + cse_p60_t33_g304582);
            double h_x1_y2 = cse_p60_t51_g304600;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p11_t1_g300609 = (cse_p11_t1_g300603 * (-(dz)));
            double cse_p60_t34_g304583 = (cse_p11_t1_g300609 * (-(dz)));
            double cse_p60_t35_g304584 = (cse_p11_t1_g300609 * dx);
            double cse_p60_t36_g304585 = (cse_p11_t1_g300609 * dy);
            double cse_p60_t37_g304586 = (cse_p11_t1_g300609 * dz);
            double cse_p60_t46_g304595 = (cse_p60_t37_g304586 + (cse_p11_t1_g300604 * (cse_p11_t55_g300545 + (-(cse_p11_t1_g300605)))));
            double cse_p60_t57_g304606 = (cse_p60_t32_g304581 + cse_p60_t35_g304584);
            double cse_p60_t58_g304607 = (cse_p60_t33_g304582 + cse_p60_t36_g304585);
            double h_x1_z2 = cse_p60_t57_g304606;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p11_t1_g300610 = (cse_p11_t1_g300604 * (cse_p11_t1_g300605 + (-(cse_p11_t54_g300544))));
            double cse_p60_t65_g304614 = (cse_p11_t1_g300610 + (cse_p11_t1_g300603 * dy * dy));
            double cse_p60_t66_g304615 = (cse_p11_t1_g300610 + cse_p60_t25_g304574);
            double h_y1_y1 = cse_p60_t65_g304614;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t49_g304598;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t50_g304599;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t45_g304594;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t58_g304607;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p11_t1_g300611 = (cse_p11_t1_g300604 * (cse_p11_t1_g300605 + (-(cse_p11_t55_g300545))));
            double cse_p60_t67_g304616 = (cse_p11_t1_g300611 + (cse_p11_t1_g300603 * dz * dz));
            double cse_p60_t68_g304617 = (cse_p11_t1_g300611 + cse_p60_t34_g304583);
            double h_z1_z1 = cse_p60_t67_g304616;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t54_g304603;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t56_g304605;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t46_g304595;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t70_g304619;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t52_g304601;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t53_g304602;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t66_g304615;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t55_g304604;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t68_g304617;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= term.r3_high)) {
              /* !COND-4 */
              double cse_p14_t1_g300840 = (3.00000000000000000e+0 * term.aa_high_i);
              double cse_p14_t2_g300841 = (3.00000000000000000e+0 * term.aa_high_j);
              double cse_p13_t80_g300795 = (cse_p14_t1_g300840 + cse_p14_t2_g300841);
              double cse_p14_t1_g300946 = (cse_p1_t2_invsqrt300486 * dx);
              double cse_p12_t5_g300617 = (term.bb_high_i * cse_p14_t1_g300946 * rr);
              double cse_p14_t1_g300947 = (cse_p1_t2_invsqrt300486 * dy);
              double cse_p12_t6_g300618 = (term.bb_high_i * cse_p14_t1_g300947 * rr);
              double cse_p14_t1_g300948 = (cse_p1_t2_invsqrt300486 * dz);
              double cse_p12_t7_g300619 = (term.bb_high_i * cse_p14_t1_g300948 * rr);
              double cse_p12_t8_g300620 = (term.bb_high_i * rr);
              double cse_p12_t9_g300621 = (term.bb_high_j * rr);
              double cse_p11_t1_g300605 = (rr * rr);
              double cse_p12_t33_g300645 = (cse_p11_t1_g300605 * cse_p13_t80_g300795);
              double cse_p11_t91_g300581 = (term.cc_high_i + term.cc_high_j + cse_p12_t33_g300645 + cse_p12_t8_g300620 + cse_p12_t8_g300620 + cse_p12_t9_g300621 + cse_p12_t9_g300621);
              double cse_p11_t93_g300583 = (-(cse_p12_t5_g300617));
              double cse_p11_t94_g300584 = (-(cse_p12_t6_g300618));
              double cse_p11_t95_g300585 = (-(cse_p12_t7_g300619));
              double energy = (term.dd_high_i + term.dd_high_j + (term.cc_high_i * rr) + (term.cc_high_j * rr) + (cse_p11_t1_g300605 * (term.bb_high_i + term.bb_high_j + (term.aa_high_i * rr) + (term.aa_high_j * rr))));
              energy_added += energy;
              *energy_accumulate += energy;
              double g_x1 = (cse_p11_t91_g300581 * cse_p14_t1_g300946);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t91_g300581 * cse_p14_t1_g300947);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t91_g300581 * cse_p14_t1_g300948);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = (cse_p11_t93_g300583 + cse_p11_t93_g300583 + (-2.00000000000000000e+0 * term.bb_high_j * cse_p14_t1_g300946 * rr) + (cse_p11_t1_g300605 * cse_p14_t1_g300946 * ((-3.00000000000000000e+0 * term.aa_high_i) + (-3.00000000000000000e+0 * term.aa_high_j))) + (-((term.cc_high_i * cse_p14_t1_g300946))) + (-((term.cc_high_j * cse_p14_t1_g300946))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (cse_p11_t94_g300584 + cse_p11_t94_g300584 + (-2.00000000000000000e+0 * term.bb_high_j * cse_p14_t1_g300947 * rr) + (cse_p11_t1_g300605 * cse_p14_t1_g300947 * ((-3.00000000000000000e+0 * term.aa_high_i) + (-3.00000000000000000e+0 * term.aa_high_j))) + (-((term.cc_high_i * cse_p14_t1_g300947))) + (-((term.cc_high_j * cse_p14_t1_g300947))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (cse_p11_t95_g300585 + cse_p11_t95_g300585 + (-2.00000000000000000e+0 * term.bb_high_j * cse_p14_t1_g300948 * rr) + (cse_p11_t1_g300605 * cse_p14_t1_g300948 * ((-3.00000000000000000e+0 * term.aa_high_i) + (-3.00000000000000000e+0 * term.aa_high_j))) + (-((term.cc_high_i * cse_p14_t1_g300948))) + (-((term.cc_high_j * cse_p14_t1_g300948))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double de_dr = (term.cc_high_i + term.cc_high_j + cse_p12_t33_g300645 + (2.00000000000000000e+0 * term.bb_high_i * rr) + (2.00000000000000000e+0 * term.bb_high_j * rr));
              double d2e_dr2 = ((2.00000000000000000e+0 * term.bb_high_i) + (2.00000000000000000e+0 * term.bb_high_j) + (6.00000000000000000e+0 * term.aa_high_i * rr) + (6.00000000000000000e+0 * term.aa_high_j * rr));
              double cse_p565_t2_invr340181 = (1.0 / (rr));
              double cse_p11_t1_g300603 = (cse_p565_t2_invr340181 * cse_p565_t2_invr340181 * d2e_dr2);
              double cse_p565_t1_invr340180 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g300604 = (cse_p565_t1_invr340180 * de_dr);
              double cse_p11_t1_g300612 = (cse_p11_t1_g300604 * (cse_p11_t1_g300605 + (-(cse_p11_t53_g300543))));
              double h_x1_x1 = (cse_p11_t1_g300612 + (cse_p11_t1_g300603 * dx * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double h_x1_y1 = ((cse_p11_t1_g300603 * dx * dy) + (cse_p11_t1_g300604 * dy * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g300608 = (cse_p11_t1_g300604 * dz);
              double h_x1_z1 = ((cse_p11_t1_g300603 * dx * dz) + (cse_p11_t1_g300608 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g300606 = (cse_p11_t1_g300603 * (-(dx)));
              double h_x1_x2 = ((cse_p11_t1_g300604 * (cse_p11_t53_g300543 + (-(cse_p11_t1_g300605)))) + (cse_p11_t1_g300606 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g300607 = (cse_p11_t1_g300603 * (-(dy)));
              double h_x1_y2 = ((cse_p11_t1_g300604 * dx * dy) + (cse_p11_t1_g300607 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g300609 = (cse_p11_t1_g300603 * (-(dz)));
              double h_x1_z2 = ((cse_p11_t1_g300608 * dx) + (cse_p11_t1_g300609 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g300610 = (cse_p11_t1_g300604 * (cse_p11_t1_g300605 + (-(cse_p11_t54_g300544))));
              double h_y1_y1 = (cse_p11_t1_g300610 + (cse_p11_t1_g300603 * dy * dy));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = ((cse_p11_t1_g300603 * dy * dz) + (cse_p11_t1_g300608 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = ((cse_p11_t1_g300604 * dx * dy) + (cse_p11_t1_g300606 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = ((cse_p11_t1_g300604 * (cse_p11_t54_g300544 + (-(cse_p11_t1_g300605)))) + (cse_p11_t1_g300607 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = ((cse_p11_t1_g300608 * dy) + (cse_p11_t1_g300609 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g300611 = (cse_p11_t1_g300604 * (cse_p11_t1_g300605 + (-(cse_p11_t55_g300545))));
              double h_z1_z1 = (cse_p11_t1_g300611 + (cse_p11_t1_g300603 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = ((cse_p11_t1_g300606 * dz) + (cse_p11_t1_g300608 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = ((cse_p11_t1_g300607 * dz) + (cse_p11_t1_g300608 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = ((cse_p11_t1_g300604 * (cse_p11_t55_g300545 + (-(cse_p11_t1_g300605)))) + (cse_p11_t1_g300609 * dz));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p11_t1_g300612 + (cse_p11_t1_g300606 * (-(dx))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = ((cse_p11_t1_g300604 * dy * (-(dx))) + (cse_p11_t1_g300606 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = ((cse_p11_t1_g300606 * (-(dz))) + (cse_p11_t1_g300608 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = (cse_p11_t1_g300610 + (cse_p11_t1_g300607 * (-(dy))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = ((cse_p11_t1_g300607 * (-(dz))) + (cse_p11_t1_g300608 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = (cse_p11_t1_g300611 + (cse_p11_t1_g300609 * (-(dz))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
            } else {
              /* !COND-ELSE-4 */
              {
                /* !COND-ELSE-5 */
                double energy = 0.00000000000000000e+0;
                energy_added += energy;
                *energy_accumulate += energy;
                double g_x1 = 0.00000000000000000e+0;
                KernelGradientAcc(term.i3x1, 0, g_x1);
                double g_y1 = 0.00000000000000000e+0;
                KernelGradientAcc(term.i3x1, 1, g_y1);
                double g_z1 = 0.00000000000000000e+0;
                KernelGradientAcc(term.i3x1, 2, g_z1);
                double g_x2 = 0.00000000000000000e+0;
                KernelGradientAcc(term.i3x2, 0, g_x2);
                double g_y2 = 0.00000000000000000e+0;
                KernelGradientAcc(term.i3x2, 1, g_y2);
                double g_z2 = 0.00000000000000000e+0;
                KernelGradientAcc(term.i3x2, 2, g_z2);
                double de_dr = 0.00000000000000000e+0;
                double d2e_dr2 = 0.00000000000000000e+0;
                double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t53_g300543))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t2_invr340181 = (1.0 / (rr));
                double cse_p565_t1_invr340180 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g300605 = (rr * rr);
                double cse_p11_t1_g300604 = (cse_p565_t1_invr340180 * de_dr);
                double cse_p11_t1_g300603 = (cse_p565_t2_invr340181 * cse_p565_t2_invr340181 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g300603 * dx * dy) + (cse_p11_t1_g300604 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double h_x1_z1 = ((cse_p11_t1_g300603 * dx * dz) + (cse_p11_t1_g300604 * dz * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g300608 = (cse_p11_t1_g300604 * dz);
                double h_x1_x2 = ((cse_p11_t1_g300603 * dx * (-(dx))) + (cse_p11_t1_g300604 * (cse_p11_t53_g300543 + (-(cse_p11_t1_g300605)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g300606 = (cse_p11_t1_g300603 * (-(dx)));
                double h_x1_y2 = ((cse_p11_t1_g300603 * dx * (-(dy))) + (cse_p11_t1_g300604 * dx * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g300607 = (cse_p11_t1_g300603 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g300603 * dx * (-(dz))) + (cse_p11_t1_g300608 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g300609 = (cse_p11_t1_g300603 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g300603 * dy * dy) + (cse_p11_t1_g300604 * (cse_p11_t1_g300605 + (-(cse_p11_t54_g300544)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g300610 = (cse_p11_t1_g300604 * (cse_p11_t1_g300605 + (-(cse_p11_t54_g300544))));
                double h_y1_z1 = ((cse_p11_t1_g300603 * dy * dz) + (cse_p11_t1_g300608 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g300604 * dx * dy) + (cse_p11_t1_g300606 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g300604 * (cse_p11_t54_g300544 + (-(cse_p11_t1_g300605)))) + (cse_p11_t1_g300607 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g300608 * dy) + (cse_p11_t1_g300609 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g300603 * dz * dz) + (cse_p11_t1_g300604 * (cse_p11_t1_g300605 + (-(cse_p11_t55_g300545)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g300611 = (cse_p11_t1_g300604 * (cse_p11_t1_g300605 + (-(cse_p11_t55_g300545))));
                double h_z1_x2 = ((cse_p11_t1_g300606 * dz) + (cse_p11_t1_g300608 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g300607 * dz) + (cse_p11_t1_g300608 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g300604 * (cse_p11_t55_g300545 + (-(cse_p11_t1_g300605)))) + (cse_p11_t1_g300609 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g300604 * (cse_p11_t1_g300605 + (-(cse_p11_t53_g300543)))) + (cse_p11_t1_g300606 * (-(dx))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p11_t1_g300612 = (cse_p11_t1_g300604 * (cse_p11_t1_g300605 + (-(cse_p11_t53_g300543))));
                double h_x2_y2 = ((cse_p11_t1_g300604 * dy * (-(dx))) + (cse_p11_t1_g300606 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double h_x2_z2 = ((cse_p11_t1_g300606 * (-(dz))) + (cse_p11_t1_g300608 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g300610 + (cse_p11_t1_g300607 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g300607 * (-(dz))) + (cse_p11_t1_g300608 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g300611 + (cse_p11_t1_g300609 * (-(dz))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
              }
            }
          }
        }
      }
    }
  }
  return energy_added;
}
void hessian_fd(const rosetta_lk_solvation_parameters& params, const rosetta_lk_solvation_term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 6;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  const double invh2 = 1.0/((h*h));
  double e0 = 0.0;
  energy(params, term, position, &e0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[term.i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 0] = saved + h;
    energy(params, term, position, &e_plus);
    position[term.i3x1 + 0] = saved - h;
    energy(params, term, position, &e_minus);
    position[term.i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x1, 0, d);
  }
  {
    double saved = position[term.i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 1] = saved + h;
    energy(params, term, position, &e_plus);
    position[term.i3x1 + 1] = saved - h;
    energy(params, term, position, &e_minus);
    position[term.i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x1, 1, d);
  }
  {
    double saved = position[term.i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 2] = saved + h;
    energy(params, term, position, &e_plus);
    position[term.i3x1 + 2] = saved - h;
    energy(params, term, position, &e_minus);
    position[term.i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x1, 2, d);
  }
  {
    double saved = position[term.i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 0] = saved + h;
    energy(params, term, position, &e_plus);
    position[term.i3x2 + 0] = saved - h;
    energy(params, term, position, &e_minus);
    position[term.i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x2, 0, d);
  }
  {
    double saved = position[term.i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 1] = saved + h;
    energy(params, term, position, &e_plus);
    position[term.i3x2 + 1] = saved - h;
    energy(params, term, position, &e_minus);
    position[term.i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x2, 1, d);
  }
  {
    double saved = position[term.i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 2] = saved + h;
    energy(params, term, position, &e_plus);
    position[term.i3x2 + 2] = saved - h;
    energy(params, term, position, &e_minus);
    position[term.i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x2, 2, d);
  }
  {
    double saved = position[term.i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 0] = saved + h;
    energy(params, term, position, &e_plus);
    position[term.i3x1 + 0] = saved - h;
    energy(params, term, position, &e_minus);
    position[term.i3x1 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 0, term.i3x1, 0, hval);
  }
  {
    double saved = position[term.i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 1] = saved + h;
    energy(params, term, position, &e_plus);
    position[term.i3x1 + 1] = saved - h;
    energy(params, term, position, &e_minus);
    position[term.i3x1 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 1, term.i3x1, 1, hval);
  }
  {
    double saved = position[term.i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 2] = saved + h;
    energy(params, term, position, &e_plus);
    position[term.i3x1 + 2] = saved - h;
    energy(params, term, position, &e_minus);
    position[term.i3x1 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 2, term.i3x1, 2, hval);
  }
  {
    double saved = position[term.i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 0] = saved + h;
    energy(params, term, position, &e_plus);
    position[term.i3x2 + 0] = saved - h;
    energy(params, term, position, &e_minus);
    position[term.i3x2 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 0, term.i3x2, 0, hval);
  }
  {
    double saved = position[term.i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 1] = saved + h;
    energy(params, term, position, &e_plus);
    position[term.i3x2 + 1] = saved - h;
    energy(params, term, position, &e_minus);
    position[term.i3x2 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 1, term.i3x2, 1, hval);
  }
  {
    double saved = position[term.i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 2] = saved + h;
    energy(params, term, position, &e_plus);
    position[term.i3x2 + 2] = saved - h;
    energy(params, term, position, &e_minus);
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
    energy(params, term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(params, term, position, &e_pm);
    position[term.i3x1 + 1] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(params, term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(params, term, position, &e_mm);
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
    energy(params, term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(params, term, position, &e_pm);
    position[term.i3x1 + 2] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(params, term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(params, term, position, &e_mm);
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
    energy(params, term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(params, term, position, &e_pm);
    position[term.i3x1 + 2] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(params, term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(params, term, position, &e_mm);
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
    energy(params, term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(params, term, position, &e_pm);
    position[term.i3x2 + 0] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(params, term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(params, term, position, &e_mm);
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
    energy(params, term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(params, term, position, &e_pm);
    position[term.i3x2 + 0] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(params, term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(params, term, position, &e_mm);
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
    energy(params, term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(params, term, position, &e_pm);
    position[term.i3x2 + 0] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(params, term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(params, term, position, &e_mm);
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
    energy(params, term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(params, term, position, &e_pm);
    position[term.i3x2 + 1] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(params, term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(params, term, position, &e_mm);
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
    energy(params, term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(params, term, position, &e_pm);
    position[term.i3x2 + 1] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(params, term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(params, term, position, &e_mm);
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
    energy(params, term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(params, term, position, &e_pm);
    position[term.i3x2 + 1] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(params, term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(params, term, position, &e_mm);
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
    energy(params, term, position, &e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(params, term, position, &e_pm);
    position[term.i3x2 + 1] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(params, term, position, &e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(params, term, position, &e_mm);
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
    energy(params, term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(params, term, position, &e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(params, term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(params, term, position, &e_mm);
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
    energy(params, term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(params, term, position, &e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(params, term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(params, term, position, &e_mm);
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
    energy(params, term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(params, term, position, &e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(params, term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(params, term, position, &e_mm);
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
    energy(params, term, position, &e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(params, term, position, &e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(params, term, position, &e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(params, term, position, &e_mm);
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
    energy(params, term, position, &e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(params, term, position, &e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(params, term, position, &e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(params, term, position, &e_mm);
    position[term.i3x2 + 2] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 2, term.i3x2, 1, hval);
  }
}

};
