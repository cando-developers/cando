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
      double cse_p1_t1_sqrt130533 = sqrt(r2);
      double cse_p1_t2_invsqrt130534 = (1.0 / (cse_p1_t1_sqrt130533));
      double rr = cse_p1_t1_sqrt130533;
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
          double cse_p528_t4_g139951 = (2.00000000000000000e+0 * rr);
          double cse_p528_t1_g139966 = (cse_p1_t2_invsqrt130534 * dx);
          double cse_p528_t6_g139953 = (cse_p1_t2_invsqrt130534 * dy);
          double cse_p528_t7_g139954 = (cse_p1_t2_invsqrt130534 * dz);
          double cse_p528_t8_g139955 = (cse_p1_t2_invsqrt130534 * rr);
          double cse_p528_t9_g139956 = (rr * rr);
          double cse_p14_t2_g130651 = (3.00000000000000000e+0 * term.aa_low_j);
          double cse_p12_t1_g130626 = (3.00000000000000000e+0 * term.aa_low_i);
          double cse_p13_t11_g130637 = (cse_p12_t1_g130626 + cse_p14_t2_g130651);
          double cse_p12_t10_g130595 = (term.bb_low_i * rr);
          double cse_p12_t11_g130596 = (term.bb_low_j * rr);
          double cse_p12_t27_g130612 = (cse_p13_t11_g130637 * cse_p528_t9_g139956);
          double cse_p13_t1_g130641 = (term.bb_low_i * cse_p528_t8_g139955);
          double cse_p12_t7_g130592 = (cse_p13_t1_g130641 * dx);
          double cse_p12_t8_g130593 = (cse_p13_t1_g130641 * dy);
          double cse_p12_t9_g130594 = (cse_p13_t1_g130641 * dz);
          double cse_p12_t12_g130597 = (cse_p528_t8_g139955 * rr);
          double cse_p11_t27_g130565 = (term.cc_low_i + term.cc_low_j + cse_p12_t10_g130595 + cse_p12_t10_g130595 + cse_p12_t11_g130596 + cse_p12_t11_g130596 + cse_p12_t27_g130612);
          double cse_p11_t31_g130569 = (-(cse_p12_t7_g130592));
          double cse_p11_t32_g130570 = (-(cse_p12_t8_g130593));
          double cse_p11_t33_g130571 = (-(cse_p12_t9_g130594));
          double energy = (term.dd_low_i + term.dd_low_j + (term.cc_low_i * rr) + (term.cc_low_j * rr) + (cse_p528_t9_g139956 * (term.bb_low_i + term.bb_low_j + (term.aa_low_i * rr) + (term.aa_low_j * rr))));
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (cse_p11_t27_g130565 * cse_p528_t1_g139966);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t27_g130565 * cse_p528_t6_g139953);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t27_g130565 * cse_p528_t7_g139954);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double cse_p12_t1_g130618 = (cse_p12_t12_g130597 * dx);
          double cse_p528_t1_g139948 = (-3.00000000000000000e+0 * cse_p12_t1_g130618);
          double cse_p13_t1_g130642 = (term.bb_low_j * cse_p528_t8_g139955);
          double cse_p63_t1_g131418 = (-2.00000000000000000e+0 * cse_p13_t1_g130642);
          double g_x2 = (cse_p11_t31_g130569 + cse_p11_t31_g130569 + (term.aa_low_i * cse_p528_t1_g139948) + (term.aa_low_j * cse_p528_t1_g139948) + (cse_p63_t1_g131418 * dx) + (-((term.cc_low_i * cse_p528_t1_g139966))) + (-((term.cc_low_j * cse_p528_t1_g139966))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double cse_p12_t1_g130619 = (cse_p12_t12_g130597 * dy);
          double cse_p528_t2_g139949 = (-3.00000000000000000e+0 * cse_p12_t1_g130619);
          double g_y2 = (cse_p11_t32_g130570 + cse_p11_t32_g130570 + (term.aa_low_i * cse_p528_t2_g139949) + (term.aa_low_j * cse_p528_t2_g139949) + (cse_p63_t1_g131418 * dy) + (-((term.cc_low_i * cse_p528_t6_g139953))) + (-((term.cc_low_j * cse_p528_t6_g139953))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double cse_p12_t1_g130620 = (cse_p12_t12_g130597 * dz);
          double cse_p528_t3_g139950 = (-3.00000000000000000e+0 * cse_p12_t1_g130620);
          double g_z2 = (cse_p11_t33_g130571 + cse_p11_t33_g130571 + (term.aa_low_i * cse_p528_t3_g139950) + (term.aa_low_j * cse_p528_t3_g139950) + (cse_p63_t1_g131418 * dz) + (-((term.cc_low_i * cse_p528_t7_g139954))) + (-((term.cc_low_j * cse_p528_t7_g139954))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (term.cc_low_i + term.cc_low_j + cse_p12_t27_g130612 + (term.bb_low_i * cse_p528_t4_g139951) + (term.bb_low_j * cse_p528_t4_g139951));
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= term.r2_high)) {
            /* !COND-3 */
            double cse_p1_t3_invr130535 = (1.0 / (rr));
            double cse_p1_t4_invr2130536 = (cse_p1_t3_invr130535 * cse_p1_t3_invr130535);
            double inv_rr = cse_p1_t3_invr130535;
            double inv_rr2 = (inv_rr * inv_rr);
            double di = (rr + (-(term.lj_radius_i)));
            double cse_p14_t7_g130656 = (di * di);
            double cse_p13_t12_g130638 = (-(cse_p14_t7_g130656));
            double cse_p1_t5_invr130537 = (1.0 / (term.lambda2_i));
            double cse_p12_t20_g130605 = (cse_p13_t12_g130638 * cse_p1_t5_invr130537);
            double cse_p11_t36_g130574 = exp(cse_p12_t20_g130605);
            double exp_i = cse_p11_t36_g130574;
            double cse_p12_t24_g130609 = (exp_i * term.pref_i_j);
            double fi = (exp_i * inv_rr2 * term.pref_i_j);
            double hi = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t5_invr130537 * di)));
            double dfi = (fi * hi);
            double dj = (rr + (-(term.lj_radius_j)));
            double cse_p14_t8_g130657 = (dj * dj);
            double cse_p13_t13_g130639 = (-(cse_p14_t8_g130657));
            double cse_p1_t6_invr130538 = (1.0 / (term.lambda2_j));
            double cse_p12_t21_g130606 = (cse_p13_t13_g130639 * cse_p1_t6_invr130538);
            double cse_p11_t37_g130575 = exp(cse_p12_t21_g130606);
            double exp_j = cse_p11_t37_g130575;
            double cse_p12_t25_g130610 = (exp_j * term.pref_j_i);
            double cse_p11_t25_g130563 = (cse_p12_t24_g130609 + cse_p12_t25_g130610);
            double fj = (exp_j * inv_rr2 * term.pref_j_i);
            double hj = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t6_invr130538 * dj)));
            double dfj = (fj * hj);
            double energy = (fi + fj);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p11_t1_g130577 = (-2.00000000000000000e+0 * inv_rr2);
            double cse_p528_t1_g139966 = (cse_p1_t2_invsqrt130534 * dx);
            double cse_p11_t1_g130582 = (cse_p11_t1_g130577 * cse_p11_t36_g130574 * cse_p1_t5_invr130537 * term.pref_i_j);
            double cse_p11_t1_g130583 = (cse_p11_t1_g130577 * cse_p11_t37_g130575 * cse_p1_t6_invr130538 * term.pref_j_i);
            double cse_p11_t1_g130584 = (cse_p11_t25_g130563 * cse_p1_t4_invr2130536 * inv_rr);
            double cse_p12_t1_g130622 = (-2.00000000000000000e+0 * cse_p11_t1_g130584);
            double cse_p63_t1_g131415 = (cse_p11_t1_g130582 * di);
            double cse_p64_t13_g131434 = (cse_p528_t1_g139966 * cse_p63_t1_g131415);
            double cse_p63_t1_g131416 = (cse_p11_t1_g130583 * dj);
            double cse_p64_t14_g131435 = (cse_p528_t1_g139966 * cse_p63_t1_g131416);
            double g_x1 = (cse_p64_t13_g131434 + cse_p64_t14_g131435 + (cse_p12_t1_g130622 * cse_p528_t1_g139966));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p11_t1_g130580 = (cse_p1_t2_invsqrt130534 * dy);
            double cse_p64_t5_g131426 = (cse_p11_t1_g130580 * cse_p63_t1_g131415);
            double cse_p64_t6_g131427 = (cse_p11_t1_g130580 * cse_p63_t1_g131416);
            double g_y1 = (cse_p64_t5_g131426 + cse_p64_t6_g131427 + (cse_p11_t1_g130580 * cse_p12_t1_g130622));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p11_t1_g130581 = (cse_p1_t2_invsqrt130534 * dz);
            double cse_p64_t7_g131428 = (cse_p11_t1_g130581 * cse_p63_t1_g131415);
            double cse_p64_t8_g131429 = (cse_p11_t1_g130581 * cse_p63_t1_g131416);
            double g_z1 = (cse_p64_t7_g131428 + cse_p64_t8_g131429 + (cse_p11_t1_g130581 * cse_p12_t1_g130622));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double cse_p12_t1_g130623 = (2.00000000000000000e+0 * cse_p11_t1_g130584);
            double g_x2 = ((cse_p12_t1_g130623 * cse_p528_t1_g139966) + (-(cse_p64_t13_g131434)) + (-(cse_p64_t14_g131435)));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = ((cse_p11_t1_g130580 * cse_p12_t1_g130623) + (-(cse_p64_t5_g131426)) + (-(cse_p64_t6_g131427)));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = ((cse_p11_t1_g130581 * cse_p12_t1_g130623) + (-(cse_p64_t7_g131428)) + (-(cse_p64_t8_g131429)));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (dfi + dfj);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= term.r3_high)) {
              /* !COND-4 */
              double cse_p12_t1_g130624 = (3.00000000000000000e+0 * term.aa_high_i);
              double cse_p12_t1_g130625 = (3.00000000000000000e+0 * term.aa_high_j);
              double cse_p13_t10_g130636 = (cse_p12_t1_g130624 + cse_p12_t1_g130625);
              double cse_p12_t5_g130590 = (term.bb_high_i * rr);
              double cse_p12_t6_g130591 = (term.bb_high_j * rr);
              double cse_p12_t1_g130621 = (rr * rr);
              double cse_p12_t26_g130611 = (cse_p12_t1_g130621 * cse_p13_t10_g130636);
              double cse_p11_t1_g130576 = (cse_p1_t2_invsqrt130534 * rr);
              double cse_p11_t1_g130585 = (term.bb_high_i * cse_p11_t1_g130576);
              double cse_p12_t13_g130598 = (cse_p11_t1_g130585 * dx);
              double cse_p12_t14_g130599 = (cse_p11_t1_g130585 * dy);
              double cse_p12_t15_g130600 = (cse_p11_t1_g130585 * dz);
              double cse_p11_t26_g130564 = (term.cc_high_i + term.cc_high_j + cse_p12_t26_g130611 + cse_p12_t5_g130590 + cse_p12_t5_g130590 + cse_p12_t6_g130591 + cse_p12_t6_g130591);
              double cse_p11_t28_g130566 = (-(cse_p12_t13_g130598));
              double cse_p11_t29_g130567 = (-(cse_p12_t14_g130599));
              double cse_p11_t30_g130568 = (-(cse_p12_t15_g130600));
              double energy = (term.dd_high_i + term.dd_high_j + (term.cc_high_i * rr) + (term.cc_high_j * rr) + (cse_p12_t1_g130621 * (term.bb_high_i + term.bb_high_j + (term.aa_high_i * rr) + (term.aa_high_j * rr))));
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p528_t1_g139966 = (cse_p1_t2_invsqrt130534 * dx);
              double g_x1 = (cse_p11_t26_g130564 * cse_p528_t1_g139966);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double cse_p11_t1_g130580 = (cse_p1_t2_invsqrt130534 * dy);
              double g_y1 = (cse_p11_t1_g130580 * cse_p11_t26_g130564);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double cse_p11_t1_g130581 = (cse_p1_t2_invsqrt130534 * dz);
              double g_z1 = (cse_p11_t1_g130581 * cse_p11_t26_g130564);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p11_t1_g130578 = (cse_p11_t1_g130576 * rr);
              double cse_p12_t1_g130618 = (cse_p11_t1_g130578 * dx);
              double cse_p13_t1_g130640 = (term.bb_high_j * cse_p11_t1_g130576);
              double cse_p63_t1_g131412 = (-3.00000000000000000e+0 * cse_p12_t1_g130618);
              double cse_p63_t1_g131417 = (-2.00000000000000000e+0 * cse_p13_t1_g130640);
              double g_x2 = (cse_p11_t28_g130566 + cse_p11_t28_g130566 + (term.aa_high_i * cse_p63_t1_g131412) + (term.aa_high_j * cse_p63_t1_g131412) + (cse_p63_t1_g131417 * dx) + (-((term.cc_high_i * cse_p528_t1_g139966))) + (-((term.cc_high_j * cse_p528_t1_g139966))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double cse_p12_t1_g130619 = (cse_p11_t1_g130578 * dy);
              double cse_p63_t1_g131413 = (-3.00000000000000000e+0 * cse_p12_t1_g130619);
              double g_y2 = (cse_p11_t29_g130567 + cse_p11_t29_g130567 + (term.aa_high_i * cse_p63_t1_g131413) + (term.aa_high_j * cse_p63_t1_g131413) + (cse_p63_t1_g131417 * dy) + (-((term.cc_high_i * cse_p11_t1_g130580))) + (-((term.cc_high_j * cse_p11_t1_g130580))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double cse_p12_t1_g130620 = (cse_p11_t1_g130578 * dz);
              double cse_p63_t1_g131414 = (-3.00000000000000000e+0 * cse_p12_t1_g130620);
              double g_z2 = (cse_p11_t30_g130568 + cse_p11_t30_g130568 + (term.aa_high_i * cse_p63_t1_g131414) + (term.aa_high_j * cse_p63_t1_g131414) + (cse_p63_t1_g131417 * dz) + (-((term.cc_high_i * cse_p11_t1_g130581))) + (-((term.cc_high_j * cse_p11_t1_g130581))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double cse_p12_t1_g130617 = (2.00000000000000000e+0 * rr);
              double de_dr = (term.cc_high_i + term.cc_high_j + cse_p12_t26_g130611 + (term.bb_high_i * cse_p12_t1_g130617) + (term.bb_high_j * cse_p12_t1_g130617));
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
            double cse_p61_t1_invr140519 = (1.0 / (term.lambda2_i));
            double exp_i = exp((-((cse_p61_t1_invr140519 * di * di))));
            double fi = (exp_i * inv_rr2 * term.pref_i_j);
            double hi = (-2.00000000000000000e+0 * (inv_rr + (cse_p61_t1_invr140519 * di)));
            double dj = (rr + (-(term.lj_radius_j)));
            double cse_p61_t2_invr140520 = (1.0 / (term.lambda2_j));
            double exp_j = exp((-((cse_p61_t2_invr140520 * dj * dj))));
            double fj = (exp_j * inv_rr2 * term.pref_j_i);
            double hj = (-2.00000000000000000e+0 * (inv_rr + (cse_p61_t2_invr140520 * dj)));
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
      double cse_p11_t53_g145529 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t54_g145530 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t55_g145531 = (dz * dz);
      double r2 = (cse_p11_t53_g145529 + cse_p11_t54_g145530 + cse_p11_t55_g145531);
      double cse_p1_t1_sqrt145471 = sqrt(r2);
      double cse_p1_t2_invsqrt145472 = (1.0 / (cse_p1_t1_sqrt145471));
      double rr = cse_p1_t1_sqrt145471;
      if ((rr <= term.r0_low)) {
        /* !COND-1 */
        double cse_p60_t43_g149578 = (rr * rr);
        double cse_p60_t71_g149606 = (-(cse_p11_t53_g145529));
        double cse_p60_t72_g149607 = (-(cse_p11_t54_g145530));
        double cse_p60_t73_g149608 = (-(cse_p11_t55_g145531));
        double cse_p60_t75_g149610 = (-(dx));
        double cse_p60_t76_g149611 = (-(dy));
        double cse_p60_t77_g149612 = (-(dz));
        double cse_p565_t1_invr185166 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr185167 = (1.0 / (rr));
        double cse_p60_t80_g149615 = (rr * rr * rr);
        double cse_p15_t75_g146009 = (-(cse_p60_t43_g149578));
        double cse_p14_t82_g145907 = (cse_p11_t53_g145529 + cse_p15_t75_g146009);
        double cse_p14_t84_g145909 = (cse_p11_t54_g145530 + cse_p15_t75_g146009);
        double cse_p14_t86_g145911 = (cse_p11_t55_g145531 + cse_p15_t75_g146009);
        double cse_p13_t109_g145810 = (cse_p60_t43_g149578 + cse_p60_t71_g149606);
        double cse_p13_t110_g145811 = (cse_p60_t43_g149578 + cse_p60_t72_g149607);
        double cse_p13_t111_g145812 = (cse_p60_t43_g149578 + cse_p60_t73_g149608);
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
        double cse_p60_t42_g149577 = (cse_p565_t1_invr185166 * de_dr);
        double cse_p13_t67_g145768 = (cse_p14_t82_g145907 * cse_p60_t42_g149577);
        double cse_p13_t68_g145769 = (cse_p14_t84_g145909 * cse_p60_t42_g149577);
        double cse_p13_t69_g145770 = (cse_p14_t86_g145911 * cse_p60_t42_g149577);
        double cse_p13_t70_g145771 = (cse_p60_t42_g149577 * dx * dy);
        double cse_p13_t71_g145772 = (cse_p60_t42_g149577 * cse_p60_t75_g149610 * dy);
        double cse_p12_t27_g145625 = (cse_p13_t109_g145810 * cse_p60_t42_g149577);
        double cse_p12_t28_g145626 = (cse_p13_t110_g145811 * cse_p60_t42_g149577);
        double cse_p12_t29_g145627 = (cse_p13_t111_g145812 * cse_p60_t42_g149577);
        double cse_p12_t32_g145630 = (cse_p60_t42_g149577 * dz);
        double cse_p13_t57_g145758 = (cse_p12_t32_g145630 * cse_p60_t75_g149610);
        double cse_p13_t58_g145759 = (cse_p12_t32_g145630 * cse_p60_t76_g149611);
        double cse_p13_t59_g145760 = (cse_p12_t32_g145630 * dx);
        double cse_p13_t60_g145761 = (cse_p12_t32_g145630 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p60_t41_g149576 = (cse_p565_t2_invr185167 * cse_p565_t2_invr185167 * d2e_dr2);
        double cse_p13_t61_g145762 = (cse_p60_t41_g149576 * dx * dx);
        double cse_p13_t62_g145763 = (cse_p60_t41_g149576 * dx * dy);
        double cse_p13_t63_g145764 = (cse_p60_t41_g149576 * dx * dz);
        double cse_p13_t64_g145765 = (cse_p60_t41_g149576 * dy * dy);
        double cse_p13_t65_g145766 = (cse_p60_t41_g149576 * dy * dz);
        double cse_p13_t66_g145767 = (cse_p60_t41_g149576 * dz * dz);
        double cse_p12_t15_g145613 = (cse_p60_t41_g149576 * cse_p60_t75_g149610);
        double cse_p13_t42_g145743 = (cse_p12_t15_g145613 * cse_p60_t75_g149610);
        double cse_p13_t43_g145744 = (cse_p12_t15_g145613 * cse_p60_t76_g149611);
        double cse_p13_t44_g145745 = (cse_p12_t15_g145613 * cse_p60_t77_g149612);
        double cse_p13_t45_g145746 = (cse_p12_t15_g145613 * dx);
        double cse_p13_t46_g145747 = (cse_p12_t15_g145613 * dy);
        double cse_p13_t47_g145748 = (cse_p12_t15_g145613 * dz);
        double cse_p12_t16_g145614 = (cse_p60_t41_g149576 * cse_p60_t76_g149611);
        double cse_p13_t48_g145749 = (cse_p12_t16_g145614 * cse_p60_t76_g149611);
        double cse_p13_t49_g145750 = (cse_p12_t16_g145614 * cse_p60_t77_g149612);
        double cse_p13_t50_g145751 = (cse_p12_t16_g145614 * dx);
        double cse_p13_t51_g145752 = (cse_p12_t16_g145614 * dy);
        double cse_p13_t52_g145753 = (cse_p12_t16_g145614 * dz);
        double cse_p12_t17_g145615 = (cse_p60_t41_g149576 * cse_p60_t77_g149612);
        double cse_p13_t53_g145754 = (cse_p12_t17_g145615 * cse_p60_t77_g149612);
        double cse_p13_t54_g145755 = (cse_p12_t17_g145615 * dx);
        double cse_p13_t55_g145756 = (cse_p12_t17_g145615 * dy);
        double cse_p13_t56_g145757 = (cse_p12_t17_g145615 * dz);
        double cse_p12_t68_g145666 = (cse_p13_t62_g145763 + cse_p13_t71_g145772);
        double cse_p12_t69_g145667 = (cse_p13_t57_g145758 + cse_p13_t63_g145764);
        double cse_p12_t70_g145668 = (cse_p13_t58_g145759 + cse_p13_t65_g145766);
        double cse_p12_t86_g145684 = (cse_p12_t28_g145626 + cse_p13_t64_g145765);
        double cse_p12_t88_g145686 = (cse_p12_t29_g145627 + cse_p13_t66_g145767);
        double cse_p12_t90_g145688 = (cse_p12_t27_g145625 + cse_p13_t61_g145762);
        double cse_p12_t63_g145661 = (cse_p13_t45_g145746 + cse_p13_t67_g145768);
        double cse_p12_t71_g145669 = (cse_p13_t46_g145747 + cse_p13_t70_g145771);
        double cse_p12_t73_g145671 = (cse_p13_t43_g145744 + cse_p13_t71_g145772);
        double cse_p12_t74_g145672 = (cse_p13_t44_g145745 + cse_p13_t57_g145758);
        double cse_p12_t75_g145673 = (cse_p13_t47_g145748 + cse_p13_t59_g145760);
        double cse_p12_t91_g145689 = (cse_p12_t27_g145625 + cse_p13_t42_g145743);
        double cse_p12_t64_g145662 = (cse_p13_t51_g145752 + cse_p13_t68_g145769);
        double cse_p12_t72_g145670 = (cse_p13_t50_g145751 + cse_p13_t70_g145771);
        double cse_p12_t76_g145674 = (cse_p13_t49_g145750 + cse_p13_t58_g145759);
        double cse_p12_t77_g145675 = (cse_p13_t52_g145753 + cse_p13_t60_g145761);
        double cse_p12_t87_g145685 = (cse_p12_t28_g145626 + cse_p13_t48_g145749);
        double cse_p12_t65_g145663 = (cse_p13_t56_g145757 + cse_p13_t69_g145770);
        double cse_p12_t78_g145676 = (cse_p13_t54_g145755 + cse_p13_t59_g145760);
        double cse_p12_t79_g145677 = (cse_p13_t55_g145756 + cse_p13_t60_g145761);
        double cse_p12_t89_g145687 = (cse_p12_t29_g145627 + cse_p13_t53_g145754);
        double h_x1_x1 = cse_p12_t90_g145688;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t68_g145666;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t69_g145667;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t63_g145661;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t72_g145670;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t78_g145676;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t86_g145684;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t70_g145668;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t71_g145669;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t64_g145662;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t79_g145677;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t88_g145686;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t75_g145673;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t77_g145675;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t65_g145663;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t91_g145689;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t73_g145671;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t74_g145672;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t87_g145685;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t76_g145674;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t89_g145687;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.r1_low)) {
          /* !COND-2 */
          double cse_p60_t38_g149573 = (cse_p1_t2_invsqrt145472 * dx);
          double cse_p60_t39_g149574 = (cse_p1_t2_invsqrt145472 * dy);
          double cse_p60_t40_g149575 = (cse_p1_t2_invsqrt145472 * dz);
          double cse_p14_t3_g145828 = (3.00000000000000000e+0 * term.aa_low_i);
          double cse_p14_t4_g145829 = (3.00000000000000000e+0 * term.aa_low_j);
          double cse_p13_t81_g145782 = (cse_p14_t3_g145828 + cse_p14_t4_g145829);
          double cse_p12_t10_g145608 = (term.bb_low_i * cse_p60_t38_g149573 * rr);
          double cse_p12_t11_g145609 = (term.bb_low_i * cse_p60_t39_g149574 * rr);
          double cse_p12_t12_g145610 = (term.bb_low_i * cse_p60_t40_g149575 * rr);
          double cse_p12_t13_g145611 = (term.bb_low_i * rr);
          double cse_p12_t14_g145612 = (term.bb_low_j * rr);
          double cse_p11_t1_g145591 = (rr * rr);
          double cse_p60_t59_g149594 = (cse_p11_t53_g145529 + (-(cse_p11_t1_g145591)));
          double cse_p60_t60_g149595 = (cse_p11_t54_g145530 + (-(cse_p11_t1_g145591)));
          double cse_p60_t61_g149596 = (cse_p11_t55_g145531 + (-(cse_p11_t1_g145591)));
          double cse_p60_t62_g149597 = (cse_p11_t1_g145591 + (-(cse_p11_t53_g145529)));
          double cse_p60_t63_g149598 = (cse_p11_t1_g145591 + (-(cse_p11_t54_g145530)));
          double cse_p60_t64_g149599 = (cse_p11_t1_g145591 + (-(cse_p11_t55_g145531)));
          double cse_p60_t74_g149609 = (-(cse_p11_t1_g145591));
          double cse_p12_t34_g145632 = (cse_p11_t1_g145591 * cse_p13_t81_g145782);
          double cse_p11_t92_g145568 = (term.cc_low_i + term.cc_low_j + cse_p12_t13_g145611 + cse_p12_t13_g145611 + cse_p12_t14_g145612 + cse_p12_t14_g145612 + cse_p12_t34_g145632);
          double cse_p11_t96_g145572 = (-(cse_p12_t10_g145608));
          double cse_p11_t97_g145573 = (-(cse_p12_t11_g145609));
          double cse_p11_t98_g145574 = (-(cse_p12_t12_g145610));
          double energy = (term.dd_low_i + term.dd_low_j + (term.cc_low_i * rr) + (term.cc_low_j * rr) + (cse_p11_t1_g145591 * (term.bb_low_i + term.bb_low_j + (term.aa_low_i * rr) + (term.aa_low_j * rr))));
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (cse_p11_t92_g145568 * cse_p60_t38_g149573);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t92_g145568 * cse_p60_t39_g149574);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t92_g145568 * cse_p60_t40_g149575);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p11_t96_g145572 + cse_p11_t96_g145572 + (-2.00000000000000000e+0 * term.bb_low_j * cse_p60_t38_g149573 * rr) + (cse_p11_t1_g145591 * cse_p60_t38_g149573 * ((-3.00000000000000000e+0 * term.aa_low_i) + (-3.00000000000000000e+0 * term.aa_low_j))) + (-((term.cc_low_i * cse_p60_t38_g149573))) + (-((term.cc_low_j * cse_p60_t38_g149573))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p11_t97_g145573 + cse_p11_t97_g145573 + (-2.00000000000000000e+0 * term.bb_low_j * cse_p60_t39_g149574 * rr) + (cse_p11_t1_g145591 * cse_p60_t39_g149574 * ((-3.00000000000000000e+0 * term.aa_low_i) + (-3.00000000000000000e+0 * term.aa_low_j))) + (-((term.cc_low_i * cse_p60_t39_g149574))) + (-((term.cc_low_j * cse_p60_t39_g149574))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p11_t98_g145574 + cse_p11_t98_g145574 + (-2.00000000000000000e+0 * term.bb_low_j * cse_p60_t40_g149575 * rr) + (cse_p11_t1_g145591 * cse_p60_t40_g149575 * ((-3.00000000000000000e+0 * term.aa_low_i) + (-3.00000000000000000e+0 * term.aa_low_j))) + (-((term.cc_low_i * cse_p60_t40_g149575))) + (-((term.cc_low_j * cse_p60_t40_g149575))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (term.cc_low_i + term.cc_low_j + cse_p12_t34_g145632 + (2.00000000000000000e+0 * term.bb_low_i * rr) + (2.00000000000000000e+0 * term.bb_low_j * rr));
          double d2e_dr2 = ((2.00000000000000000e+0 * term.bb_low_i) + (2.00000000000000000e+0 * term.bb_low_j) + (6.00000000000000000e+0 * term.aa_low_i * rr) + (6.00000000000000000e+0 * term.aa_low_j * rr));
          double cse_p565_t2_invr185167 = (1.0 / (rr));
          double cse_p11_t1_g145589 = (cse_p565_t2_invr185167 * cse_p565_t2_invr185167 * d2e_dr2);
          double cse_p60_t1_g149536 = (cse_p11_t1_g145589 * (-(dx)));
          double cse_p60_t2_g149537 = (cse_p11_t1_g145589 * (-(dy)));
          double cse_p60_t3_g149538 = (cse_p11_t1_g145589 * (-(dz)));
          double cse_p60_t4_g149539 = (cse_p11_t1_g145589 * dx * dx);
          double cse_p60_t5_g149540 = (cse_p11_t1_g145589 * dx * dy);
          double cse_p60_t6_g149541 = (cse_p11_t1_g145589 * dx * dz);
          double cse_p60_t7_g149542 = (cse_p11_t1_g145589 * dy * dy);
          double cse_p60_t8_g149543 = (cse_p11_t1_g145589 * dy * dz);
          double cse_p60_t9_g149544 = (cse_p11_t1_g145589 * dz * dz);
          double cse_p14_t46_g145871 = (cse_p60_t1_g149536 * (-(dx)));
          double cse_p14_t47_g145872 = (cse_p60_t1_g149536 * (-(dy)));
          double cse_p14_t48_g145873 = (cse_p60_t1_g149536 * (-(dz)));
          double cse_p14_t49_g145874 = (cse_p60_t1_g149536 * dx);
          double cse_p14_t50_g145875 = (cse_p60_t1_g149536 * dy);
          double cse_p14_t51_g145876 = (cse_p60_t1_g149536 * dz);
          double cse_p14_t52_g145877 = (cse_p60_t2_g149537 * (-(dy)));
          double cse_p14_t53_g145878 = (cse_p60_t2_g149537 * (-(dz)));
          double cse_p14_t54_g145879 = (cse_p60_t2_g149537 * dx);
          double cse_p14_t55_g145880 = (cse_p60_t2_g149537 * dy);
          double cse_p14_t56_g145881 = (cse_p60_t2_g149537 * dz);
          double cse_p14_t57_g145882 = (cse_p60_t3_g149538 * (-(dz)));
          double cse_p14_t58_g145883 = (cse_p60_t3_g149538 * dx);
          double cse_p14_t59_g145884 = (cse_p60_t3_g149538 * dy);
          double cse_p14_t60_g145885 = (cse_p60_t3_g149538 * dz);
          double cse_p565_t1_invr185166 = (1.0 / ((rr * rr * rr)));
          double cse_p11_t1_g145590 = (cse_p565_t1_invr185166 * de_dr);
          double cse_p60_t10_g149545 = (cse_p11_t1_g145590 * (cse_p11_t53_g145529 + cse_p60_t74_g149609));
          double cse_p60_t11_g149546 = (cse_p11_t1_g145590 * (cse_p11_t54_g145530 + cse_p60_t74_g149609));
          double cse_p60_t12_g149547 = (cse_p11_t1_g145590 * (cse_p11_t55_g145531 + cse_p60_t74_g149609));
          double cse_p60_t13_g149548 = (cse_p11_t1_g145590 * cse_p60_t62_g149597);
          double cse_p60_t14_g149549 = (cse_p11_t1_g145590 * cse_p60_t63_g149598);
          double cse_p60_t15_g149550 = (cse_p11_t1_g145590 * cse_p60_t64_g149599);
          double cse_p60_t16_g149551 = (cse_p11_t1_g145590 * dx * dy);
          double cse_p60_t17_g149552 = (cse_p11_t1_g145590 * dy * (-(dx)));
          double cse_p60_t18_g149553 = (cse_p11_t1_g145590 * dz);
          double cse_p60_t47_g149582 = (cse_p60_t17_g149552 + cse_p60_t5_g149540);
          double cse_p14_t42_g145867 = (cse_p60_t18_g149553 * (-(dx)));
          double cse_p14_t43_g145868 = (cse_p60_t18_g149553 * (-(dy)));
          double cse_p14_t44_g145869 = (cse_p60_t18_g149553 * dx);
          double cse_p14_t45_g145870 = (cse_p60_t18_g149553 * dy);
          double cse_p13_t107_g145808 = (cse_p60_t13_g149548 + cse_p60_t4_g149539);
          double h_x1_x1 = cse_p13_t107_g145808;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t47_g149582;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t83_g145784 = (cse_p14_t42_g145867 + cse_p60_t6_g149541);
          double cse_p13_t84_g145785 = (cse_p14_t43_g145868 + cse_p60_t8_g149543);
          double h_x1_z1 = cse_p13_t83_g145784;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t77_g145778 = (cse_p14_t49_g145874 + cse_p60_t10_g149545);
          double cse_p13_t85_g145786 = (cse_p14_t50_g145875 + cse_p60_t16_g149551);
          double cse_p13_t87_g145788 = (cse_p14_t47_g145872 + cse_p60_t17_g149552);
          double cse_p13_t88_g145789 = (cse_p14_t42_g145867 + cse_p14_t48_g145873);
          double cse_p13_t89_g145790 = (cse_p14_t44_g145869 + cse_p14_t51_g145876);
          double cse_p13_t108_g145809 = (cse_p14_t46_g145871 + cse_p60_t13_g149548);
          double h_x1_x2 = cse_p13_t77_g145778;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t78_g145779 = (cse_p14_t55_g145880 + cse_p60_t11_g149546);
          double cse_p13_t86_g145787 = (cse_p14_t54_g145879 + cse_p60_t16_g149551);
          double cse_p13_t90_g145791 = (cse_p14_t43_g145868 + cse_p14_t53_g145878);
          double cse_p13_t91_g145792 = (cse_p14_t45_g145870 + cse_p14_t56_g145881);
          double h_x1_y2 = cse_p13_t86_g145787;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t79_g145780 = (cse_p14_t60_g145885 + cse_p60_t12_g149547);
          double cse_p13_t92_g145793 = (cse_p14_t44_g145869 + cse_p14_t58_g145883);
          double cse_p13_t93_g145794 = (cse_p14_t45_g145870 + cse_p14_t59_g145884);
          double h_x1_z2 = cse_p13_t92_g145793;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t103_g145804 = (cse_p60_t14_g149549 + cse_p60_t7_g149542);
          double cse_p13_t104_g145805 = (cse_p14_t52_g145877 + cse_p60_t14_g149549);
          double h_y1_y1 = cse_p13_t103_g145804;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t84_g145785;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t85_g145786;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t78_g145779;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t93_g145794;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t105_g145806 = (cse_p60_t15_g149550 + cse_p60_t9_g149544);
          double cse_p13_t106_g145807 = (cse_p14_t57_g145882 + cse_p60_t15_g149550);
          double h_z1_z1 = cse_p13_t105_g145806;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t89_g145790;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t91_g145792;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t79_g145780;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t108_g145809;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t87_g145788;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t88_g145789;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t104_g145805;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t90_g145791;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t106_g145807;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= term.r2_high)) {
            /* !COND-3 */
            double cse_p565_t2_invr185167 = (1.0 / (rr));
            double cse_p1_t4_invr2145474 = (cse_p565_t2_invr185167 * cse_p565_t2_invr185167);
            double inv_rr = cse_p565_t2_invr185167;
            double inv_rr2 = (inv_rr * inv_rr);
            double cse_p11_t1_g145477 = (2.00000000000000000e+0 * inv_rr2);
            double di = (rr + (-(term.lj_radius_i)));
            double cse_p14_t63_g145888 = (di * di);
            double cse_p13_t112_g145813 = (-(cse_p14_t63_g145888));
            double cse_p1_t5_invr145475 = (1.0 / (term.lambda2_i));
            double cse_p12_t54_g145652 = (cse_p13_t112_g145813 * cse_p1_t5_invr145475);
            double cse_p11_t108_g145584 = exp(cse_p12_t54_g145652);
            double exp_i = cse_p11_t108_g145584;
            double cse_p12_t60_g145658 = (exp_i * term.pref_i_j);
            double fi = (exp_i * inv_rr2 * term.pref_i_j);
            double hi = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t5_invr145475 * di)));
            double dfi = (fi * hi);
            double d2fi = (fi * (cse_p11_t1_g145477 + (-2.00000000000000000e+0 * cse_p1_t5_invr145475) + (hi * hi)));
            double dj = (rr + (-(term.lj_radius_j)));
            double cse_p14_t64_g145889 = (dj * dj);
            double cse_p13_t113_g145814 = (-(cse_p14_t64_g145889));
            double cse_p1_t6_invr145476 = (1.0 / (term.lambda2_j));
            double cse_p12_t55_g145653 = (cse_p13_t113_g145814 * cse_p1_t6_invr145476);
            double cse_p11_t109_g145585 = exp(cse_p12_t55_g145653);
            double exp_j = cse_p11_t109_g145585;
            double cse_p12_t61_g145659 = (exp_j * term.pref_j_i);
            double cse_p11_t87_g145563 = (cse_p12_t60_g145658 + cse_p12_t61_g145659);
            double fj = (exp_j * inv_rr2 * term.pref_j_i);
            double hj = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t6_invr145476 * dj)));
            double dfj = (fj * hj);
            double d2fj = (fj * (cse_p11_t1_g145477 + (-2.00000000000000000e+0 * cse_p1_t6_invr145476) + (hj * hj)));
            double energy = (fi + fj);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p14_t1_g145932 = (cse_p1_t2_invsqrt145472 * dx);
            double g_x1 = ((-2.00000000000000000e+0 * cse_p11_t87_g145563 * cse_p14_t1_g145932 * cse_p1_t4_invr2145474 * inv_rr) + (cse_p14_t1_g145932 * inv_rr2 * ((-2.00000000000000000e+0 * cse_p11_t108_g145584 * cse_p1_t5_invr145475 * di * term.pref_i_j) + (-2.00000000000000000e+0 * cse_p11_t109_g145585 * cse_p1_t6_invr145476 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p14_t1_g145933 = (cse_p1_t2_invsqrt145472 * dy);
            double g_y1 = ((-2.00000000000000000e+0 * cse_p11_t87_g145563 * cse_p14_t1_g145933 * cse_p1_t4_invr2145474 * inv_rr) + (cse_p14_t1_g145933 * inv_rr2 * ((-2.00000000000000000e+0 * cse_p11_t108_g145584 * cse_p1_t5_invr145475 * di * term.pref_i_j) + (-2.00000000000000000e+0 * cse_p11_t109_g145585 * cse_p1_t6_invr145476 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p14_t1_g145934 = (cse_p1_t2_invsqrt145472 * dz);
            double g_z1 = ((-2.00000000000000000e+0 * cse_p11_t87_g145563 * cse_p14_t1_g145934 * cse_p1_t4_invr2145474 * inv_rr) + (cse_p14_t1_g145934 * inv_rr2 * ((-2.00000000000000000e+0 * cse_p11_t108_g145584 * cse_p1_t5_invr145475 * di * term.pref_i_j) + (-2.00000000000000000e+0 * cse_p11_t109_g145585 * cse_p1_t6_invr145476 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = ((2.00000000000000000e+0 * cse_p11_t87_g145563 * cse_p14_t1_g145932 * cse_p1_t4_invr2145474 * inv_rr) + (cse_p14_t1_g145932 * inv_rr2 * ((2.00000000000000000e+0 * cse_p11_t108_g145584 * cse_p1_t5_invr145475 * di * term.pref_i_j) + (2.00000000000000000e+0 * cse_p11_t109_g145585 * cse_p1_t6_invr145476 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = ((2.00000000000000000e+0 * cse_p11_t87_g145563 * cse_p14_t1_g145933 * cse_p1_t4_invr2145474 * inv_rr) + (cse_p14_t1_g145933 * inv_rr2 * ((2.00000000000000000e+0 * cse_p11_t108_g145584 * cse_p1_t5_invr145475 * di * term.pref_i_j) + (2.00000000000000000e+0 * cse_p11_t109_g145585 * cse_p1_t6_invr145476 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = ((2.00000000000000000e+0 * cse_p11_t87_g145563 * cse_p14_t1_g145934 * cse_p1_t4_invr2145474 * inv_rr) + (cse_p14_t1_g145934 * inv_rr2 * ((2.00000000000000000e+0 * cse_p11_t108_g145584 * cse_p1_t5_invr145475 * di * term.pref_i_j) + (2.00000000000000000e+0 * cse_p11_t109_g145585 * cse_p1_t6_invr145476 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (dfi + dfj);
            double d2e_dr2 = (d2fi + d2fj);
            double cse_p11_t1_g145589 = (cse_p565_t2_invr185167 * cse_p565_t2_invr185167 * d2e_dr2);
            double cse_p565_t1_invr185166 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g145590 = (cse_p565_t1_invr185166 * de_dr);
            double cse_p11_t1_g145591 = (rr * rr);
            double cse_p11_t1_g145598 = (cse_p11_t1_g145590 * (cse_p11_t1_g145591 + (-(cse_p11_t53_g145529))));
            double cse_p60_t69_g149604 = (cse_p11_t1_g145598 + (cse_p11_t1_g145589 * dx * dx));
            double h_x1_x1 = cse_p60_t69_g149604;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = ((cse_p11_t1_g145589 * dx * dy) + (cse_p11_t1_g145590 * dy * (-(dx))));
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p11_t1_g145594 = (cse_p11_t1_g145590 * dz);
            double cse_p60_t30_g149565 = (cse_p11_t1_g145594 * (-(dx)));
            double cse_p60_t31_g149566 = (cse_p11_t1_g145594 * (-(dy)));
            double cse_p60_t32_g149567 = (cse_p11_t1_g145594 * dx);
            double cse_p60_t33_g149568 = (cse_p11_t1_g145594 * dy);
            double cse_p60_t48_g149583 = (cse_p60_t30_g149565 + (cse_p11_t1_g145589 * dx * dz));
            double cse_p60_t49_g149584 = (cse_p60_t31_g149566 + (cse_p11_t1_g145589 * dy * dz));
            double h_x1_z1 = cse_p60_t48_g149583;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p11_t1_g145592 = (cse_p11_t1_g145589 * (-(dx)));
            double cse_p60_t19_g149554 = (cse_p11_t1_g145592 * (-(dx)));
            double cse_p60_t20_g149555 = (cse_p11_t1_g145592 * (-(dy)));
            double cse_p60_t21_g149556 = (cse_p11_t1_g145592 * (-(dz)));
            double cse_p60_t22_g149557 = (cse_p11_t1_g145592 * dx);
            double cse_p60_t23_g149558 = (cse_p11_t1_g145592 * dy);
            double cse_p60_t24_g149559 = (cse_p11_t1_g145592 * dz);
            double cse_p60_t44_g149579 = (cse_p60_t22_g149557 + (cse_p11_t1_g145590 * (cse_p11_t53_g145529 + (-(cse_p11_t1_g145591)))));
            double cse_p60_t50_g149585 = (cse_p60_t23_g149558 + (cse_p11_t1_g145590 * dx * dy));
            double cse_p60_t52_g149587 = (cse_p60_t20_g149555 + (cse_p11_t1_g145590 * dy * (-(dx))));
            double cse_p60_t53_g149588 = (cse_p60_t21_g149556 + cse_p60_t30_g149565);
            double cse_p60_t54_g149589 = (cse_p60_t24_g149559 + cse_p60_t32_g149567);
            double cse_p60_t70_g149605 = (cse_p11_t1_g145598 + cse_p60_t19_g149554);
            double h_x1_x2 = cse_p60_t44_g149579;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p11_t1_g145593 = (cse_p11_t1_g145589 * (-(dy)));
            double cse_p60_t25_g149560 = (cse_p11_t1_g145593 * (-(dy)));
            double cse_p60_t26_g149561 = (cse_p11_t1_g145593 * (-(dz)));
            double cse_p60_t27_g149562 = (cse_p11_t1_g145593 * dx);
            double cse_p60_t28_g149563 = (cse_p11_t1_g145593 * dy);
            double cse_p60_t29_g149564 = (cse_p11_t1_g145593 * dz);
            double cse_p60_t45_g149580 = (cse_p60_t28_g149563 + (cse_p11_t1_g145590 * (cse_p11_t54_g145530 + (-(cse_p11_t1_g145591)))));
            double cse_p60_t51_g149586 = (cse_p60_t27_g149562 + (cse_p11_t1_g145590 * dx * dy));
            double cse_p60_t55_g149590 = (cse_p60_t26_g149561 + cse_p60_t31_g149566);
            double cse_p60_t56_g149591 = (cse_p60_t29_g149564 + cse_p60_t33_g149568);
            double h_x1_y2 = cse_p60_t51_g149586;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p11_t1_g145595 = (cse_p11_t1_g145589 * (-(dz)));
            double cse_p60_t34_g149569 = (cse_p11_t1_g145595 * (-(dz)));
            double cse_p60_t35_g149570 = (cse_p11_t1_g145595 * dx);
            double cse_p60_t36_g149571 = (cse_p11_t1_g145595 * dy);
            double cse_p60_t37_g149572 = (cse_p11_t1_g145595 * dz);
            double cse_p60_t46_g149581 = (cse_p60_t37_g149572 + (cse_p11_t1_g145590 * (cse_p11_t55_g145531 + (-(cse_p11_t1_g145591)))));
            double cse_p60_t57_g149592 = (cse_p60_t32_g149567 + cse_p60_t35_g149570);
            double cse_p60_t58_g149593 = (cse_p60_t33_g149568 + cse_p60_t36_g149571);
            double h_x1_z2 = cse_p60_t57_g149592;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p11_t1_g145596 = (cse_p11_t1_g145590 * (cse_p11_t1_g145591 + (-(cse_p11_t54_g145530))));
            double cse_p60_t65_g149600 = (cse_p11_t1_g145596 + (cse_p11_t1_g145589 * dy * dy));
            double cse_p60_t66_g149601 = (cse_p11_t1_g145596 + cse_p60_t25_g149560);
            double h_y1_y1 = cse_p60_t65_g149600;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t49_g149584;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t50_g149585;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t45_g149580;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t58_g149593;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p11_t1_g145597 = (cse_p11_t1_g145590 * (cse_p11_t1_g145591 + (-(cse_p11_t55_g145531))));
            double cse_p60_t67_g149602 = (cse_p11_t1_g145597 + (cse_p11_t1_g145589 * dz * dz));
            double cse_p60_t68_g149603 = (cse_p11_t1_g145597 + cse_p60_t34_g149569);
            double h_z1_z1 = cse_p60_t67_g149602;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t54_g149589;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t56_g149591;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t46_g149581;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t70_g149605;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t52_g149587;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t53_g149588;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t66_g149601;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t55_g149590;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t68_g149603;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= term.r3_high)) {
              /* !COND-4 */
              double cse_p14_t1_g145826 = (3.00000000000000000e+0 * term.aa_high_i);
              double cse_p14_t2_g145827 = (3.00000000000000000e+0 * term.aa_high_j);
              double cse_p13_t80_g145781 = (cse_p14_t1_g145826 + cse_p14_t2_g145827);
              double cse_p14_t1_g145932 = (cse_p1_t2_invsqrt145472 * dx);
              double cse_p12_t5_g145603 = (term.bb_high_i * cse_p14_t1_g145932 * rr);
              double cse_p14_t1_g145933 = (cse_p1_t2_invsqrt145472 * dy);
              double cse_p12_t6_g145604 = (term.bb_high_i * cse_p14_t1_g145933 * rr);
              double cse_p14_t1_g145934 = (cse_p1_t2_invsqrt145472 * dz);
              double cse_p12_t7_g145605 = (term.bb_high_i * cse_p14_t1_g145934 * rr);
              double cse_p12_t8_g145606 = (term.bb_high_i * rr);
              double cse_p12_t9_g145607 = (term.bb_high_j * rr);
              double cse_p11_t1_g145591 = (rr * rr);
              double cse_p12_t33_g145631 = (cse_p11_t1_g145591 * cse_p13_t80_g145781);
              double cse_p11_t91_g145567 = (term.cc_high_i + term.cc_high_j + cse_p12_t33_g145631 + cse_p12_t8_g145606 + cse_p12_t8_g145606 + cse_p12_t9_g145607 + cse_p12_t9_g145607);
              double cse_p11_t93_g145569 = (-(cse_p12_t5_g145603));
              double cse_p11_t94_g145570 = (-(cse_p12_t6_g145604));
              double cse_p11_t95_g145571 = (-(cse_p12_t7_g145605));
              double energy = (term.dd_high_i + term.dd_high_j + (term.cc_high_i * rr) + (term.cc_high_j * rr) + (cse_p11_t1_g145591 * (term.bb_high_i + term.bb_high_j + (term.aa_high_i * rr) + (term.aa_high_j * rr))));
              energy_added += energy;
              *energy_accumulate += energy;
              double g_x1 = (cse_p11_t91_g145567 * cse_p14_t1_g145932);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t91_g145567 * cse_p14_t1_g145933);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t91_g145567 * cse_p14_t1_g145934);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = (cse_p11_t93_g145569 + cse_p11_t93_g145569 + (-2.00000000000000000e+0 * term.bb_high_j * cse_p14_t1_g145932 * rr) + (cse_p11_t1_g145591 * cse_p14_t1_g145932 * ((-3.00000000000000000e+0 * term.aa_high_i) + (-3.00000000000000000e+0 * term.aa_high_j))) + (-((term.cc_high_i * cse_p14_t1_g145932))) + (-((term.cc_high_j * cse_p14_t1_g145932))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (cse_p11_t94_g145570 + cse_p11_t94_g145570 + (-2.00000000000000000e+0 * term.bb_high_j * cse_p14_t1_g145933 * rr) + (cse_p11_t1_g145591 * cse_p14_t1_g145933 * ((-3.00000000000000000e+0 * term.aa_high_i) + (-3.00000000000000000e+0 * term.aa_high_j))) + (-((term.cc_high_i * cse_p14_t1_g145933))) + (-((term.cc_high_j * cse_p14_t1_g145933))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (cse_p11_t95_g145571 + cse_p11_t95_g145571 + (-2.00000000000000000e+0 * term.bb_high_j * cse_p14_t1_g145934 * rr) + (cse_p11_t1_g145591 * cse_p14_t1_g145934 * ((-3.00000000000000000e+0 * term.aa_high_i) + (-3.00000000000000000e+0 * term.aa_high_j))) + (-((term.cc_high_i * cse_p14_t1_g145934))) + (-((term.cc_high_j * cse_p14_t1_g145934))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double de_dr = (term.cc_high_i + term.cc_high_j + cse_p12_t33_g145631 + (2.00000000000000000e+0 * term.bb_high_i * rr) + (2.00000000000000000e+0 * term.bb_high_j * rr));
              double d2e_dr2 = ((2.00000000000000000e+0 * term.bb_high_i) + (2.00000000000000000e+0 * term.bb_high_j) + (6.00000000000000000e+0 * term.aa_high_i * rr) + (6.00000000000000000e+0 * term.aa_high_j * rr));
              double cse_p565_t2_invr185167 = (1.0 / (rr));
              double cse_p11_t1_g145589 = (cse_p565_t2_invr185167 * cse_p565_t2_invr185167 * d2e_dr2);
              double cse_p565_t1_invr185166 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g145590 = (cse_p565_t1_invr185166 * de_dr);
              double cse_p11_t1_g145598 = (cse_p11_t1_g145590 * (cse_p11_t1_g145591 + (-(cse_p11_t53_g145529))));
              double h_x1_x1 = (cse_p11_t1_g145598 + (cse_p11_t1_g145589 * dx * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double h_x1_y1 = ((cse_p11_t1_g145589 * dx * dy) + (cse_p11_t1_g145590 * dy * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g145594 = (cse_p11_t1_g145590 * dz);
              double h_x1_z1 = ((cse_p11_t1_g145589 * dx * dz) + (cse_p11_t1_g145594 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g145592 = (cse_p11_t1_g145589 * (-(dx)));
              double h_x1_x2 = ((cse_p11_t1_g145590 * (cse_p11_t53_g145529 + (-(cse_p11_t1_g145591)))) + (cse_p11_t1_g145592 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g145593 = (cse_p11_t1_g145589 * (-(dy)));
              double h_x1_y2 = ((cse_p11_t1_g145590 * dx * dy) + (cse_p11_t1_g145593 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g145595 = (cse_p11_t1_g145589 * (-(dz)));
              double h_x1_z2 = ((cse_p11_t1_g145594 * dx) + (cse_p11_t1_g145595 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g145596 = (cse_p11_t1_g145590 * (cse_p11_t1_g145591 + (-(cse_p11_t54_g145530))));
              double h_y1_y1 = (cse_p11_t1_g145596 + (cse_p11_t1_g145589 * dy * dy));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = ((cse_p11_t1_g145589 * dy * dz) + (cse_p11_t1_g145594 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = ((cse_p11_t1_g145590 * dx * dy) + (cse_p11_t1_g145592 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = ((cse_p11_t1_g145590 * (cse_p11_t54_g145530 + (-(cse_p11_t1_g145591)))) + (cse_p11_t1_g145593 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = ((cse_p11_t1_g145594 * dy) + (cse_p11_t1_g145595 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g145597 = (cse_p11_t1_g145590 * (cse_p11_t1_g145591 + (-(cse_p11_t55_g145531))));
              double h_z1_z1 = (cse_p11_t1_g145597 + (cse_p11_t1_g145589 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = ((cse_p11_t1_g145592 * dz) + (cse_p11_t1_g145594 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = ((cse_p11_t1_g145593 * dz) + (cse_p11_t1_g145594 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = ((cse_p11_t1_g145590 * (cse_p11_t55_g145531 + (-(cse_p11_t1_g145591)))) + (cse_p11_t1_g145595 * dz));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p11_t1_g145598 + (cse_p11_t1_g145592 * (-(dx))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = ((cse_p11_t1_g145590 * dy * (-(dx))) + (cse_p11_t1_g145592 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = ((cse_p11_t1_g145592 * (-(dz))) + (cse_p11_t1_g145594 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = (cse_p11_t1_g145596 + (cse_p11_t1_g145593 * (-(dy))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = ((cse_p11_t1_g145593 * (-(dz))) + (cse_p11_t1_g145594 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = (cse_p11_t1_g145597 + (cse_p11_t1_g145595 * (-(dz))));
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
                double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t53_g145529))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t2_invr185167 = (1.0 / (rr));
                double cse_p565_t1_invr185166 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g145591 = (rr * rr);
                double cse_p11_t1_g145590 = (cse_p565_t1_invr185166 * de_dr);
                double cse_p11_t1_g145589 = (cse_p565_t2_invr185167 * cse_p565_t2_invr185167 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g145589 * dx * dy) + (cse_p11_t1_g145590 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double h_x1_z1 = ((cse_p11_t1_g145589 * dx * dz) + (cse_p11_t1_g145590 * dz * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g145594 = (cse_p11_t1_g145590 * dz);
                double h_x1_x2 = ((cse_p11_t1_g145589 * dx * (-(dx))) + (cse_p11_t1_g145590 * (cse_p11_t53_g145529 + (-(cse_p11_t1_g145591)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g145592 = (cse_p11_t1_g145589 * (-(dx)));
                double h_x1_y2 = ((cse_p11_t1_g145589 * dx * (-(dy))) + (cse_p11_t1_g145590 * dx * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g145593 = (cse_p11_t1_g145589 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g145589 * dx * (-(dz))) + (cse_p11_t1_g145594 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g145595 = (cse_p11_t1_g145589 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g145589 * dy * dy) + (cse_p11_t1_g145590 * (cse_p11_t1_g145591 + (-(cse_p11_t54_g145530)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g145596 = (cse_p11_t1_g145590 * (cse_p11_t1_g145591 + (-(cse_p11_t54_g145530))));
                double h_y1_z1 = ((cse_p11_t1_g145589 * dy * dz) + (cse_p11_t1_g145594 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g145590 * dx * dy) + (cse_p11_t1_g145592 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g145590 * (cse_p11_t54_g145530 + (-(cse_p11_t1_g145591)))) + (cse_p11_t1_g145593 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g145594 * dy) + (cse_p11_t1_g145595 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g145589 * dz * dz) + (cse_p11_t1_g145590 * (cse_p11_t1_g145591 + (-(cse_p11_t55_g145531)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g145597 = (cse_p11_t1_g145590 * (cse_p11_t1_g145591 + (-(cse_p11_t55_g145531))));
                double h_z1_x2 = ((cse_p11_t1_g145592 * dz) + (cse_p11_t1_g145594 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g145593 * dz) + (cse_p11_t1_g145594 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g145590 * (cse_p11_t55_g145531 + (-(cse_p11_t1_g145591)))) + (cse_p11_t1_g145595 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g145590 * (cse_p11_t1_g145591 + (-(cse_p11_t53_g145529)))) + (cse_p11_t1_g145592 * (-(dx))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p11_t1_g145598 = (cse_p11_t1_g145590 * (cse_p11_t1_g145591 + (-(cse_p11_t53_g145529))));
                double h_x2_y2 = ((cse_p11_t1_g145590 * dy * (-(dx))) + (cse_p11_t1_g145592 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double h_x2_z2 = ((cse_p11_t1_g145592 * (-(dz))) + (cse_p11_t1_g145594 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g145596 + (cse_p11_t1_g145593 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g145593 * (-(dz))) + (cse_p11_t1_g145594 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g145597 + (cse_p11_t1_g145595 * (-(dz))));
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
