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
      double cse_p1_t1_sqrt503617 = sqrt(r2);
      double cse_p1_t2_invsqrt503618 = (1.0 / (cse_p1_t1_sqrt503617));
      double rr = cse_p1_t1_sqrt503617;
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
          double cse_p528_t4_g513035 = (2.00000000000000000e+0 * rr);
          double cse_p528_t1_g513050 = (cse_p1_t2_invsqrt503618 * dx);
          double cse_p528_t6_g513037 = (cse_p1_t2_invsqrt503618 * dy);
          double cse_p528_t7_g513038 = (cse_p1_t2_invsqrt503618 * dz);
          double cse_p528_t8_g513039 = (cse_p1_t2_invsqrt503618 * rr);
          double cse_p528_t9_g513040 = (rr * rr);
          double cse_p14_t2_g503735 = (3.00000000000000000e+0 * term.aa_low_j);
          double cse_p12_t1_g503710 = (3.00000000000000000e+0 * term.aa_low_i);
          double cse_p13_t11_g503721 = (cse_p12_t1_g503710 + cse_p14_t2_g503735);
          double cse_p12_t10_g503679 = (term.bb_low_i * rr);
          double cse_p12_t11_g503680 = (term.bb_low_j * rr);
          double cse_p12_t27_g503696 = (cse_p13_t11_g503721 * cse_p528_t9_g513040);
          double cse_p13_t1_g503725 = (term.bb_low_i * cse_p528_t8_g513039);
          double cse_p12_t7_g503676 = (cse_p13_t1_g503725 * dx);
          double cse_p12_t8_g503677 = (cse_p13_t1_g503725 * dy);
          double cse_p12_t9_g503678 = (cse_p13_t1_g503725 * dz);
          double cse_p12_t12_g503681 = (cse_p528_t8_g513039 * rr);
          double cse_p11_t27_g503649 = (term.cc_low_i + term.cc_low_j + cse_p12_t10_g503679 + cse_p12_t10_g503679 + cse_p12_t11_g503680 + cse_p12_t11_g503680 + cse_p12_t27_g503696);
          double cse_p11_t31_g503653 = (-(cse_p12_t7_g503676));
          double cse_p11_t32_g503654 = (-(cse_p12_t8_g503677));
          double cse_p11_t33_g503655 = (-(cse_p12_t9_g503678));
          double energy = (term.dd_low_i + term.dd_low_j + (term.cc_low_i * rr) + (term.cc_low_j * rr) + (cse_p528_t9_g513040 * (term.bb_low_i + term.bb_low_j + (term.aa_low_i * rr) + (term.aa_low_j * rr))));
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (cse_p11_t27_g503649 * cse_p528_t1_g513050);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t27_g503649 * cse_p528_t6_g513037);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t27_g503649 * cse_p528_t7_g513038);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double cse_p12_t1_g503702 = (cse_p12_t12_g503681 * dx);
          double cse_p528_t1_g513032 = (-3.00000000000000000e+0 * cse_p12_t1_g503702);
          double cse_p13_t1_g503726 = (term.bb_low_j * cse_p528_t8_g513039);
          double cse_p63_t1_g504502 = (-2.00000000000000000e+0 * cse_p13_t1_g503726);
          double g_x2 = (cse_p11_t31_g503653 + cse_p11_t31_g503653 + (term.aa_low_i * cse_p528_t1_g513032) + (term.aa_low_j * cse_p528_t1_g513032) + (cse_p63_t1_g504502 * dx) + (-((term.cc_low_i * cse_p528_t1_g513050))) + (-((term.cc_low_j * cse_p528_t1_g513050))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double cse_p12_t1_g503703 = (cse_p12_t12_g503681 * dy);
          double cse_p528_t2_g513033 = (-3.00000000000000000e+0 * cse_p12_t1_g503703);
          double g_y2 = (cse_p11_t32_g503654 + cse_p11_t32_g503654 + (term.aa_low_i * cse_p528_t2_g513033) + (term.aa_low_j * cse_p528_t2_g513033) + (cse_p63_t1_g504502 * dy) + (-((term.cc_low_i * cse_p528_t6_g513037))) + (-((term.cc_low_j * cse_p528_t6_g513037))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double cse_p12_t1_g503704 = (cse_p12_t12_g503681 * dz);
          double cse_p528_t3_g513034 = (-3.00000000000000000e+0 * cse_p12_t1_g503704);
          double g_z2 = (cse_p11_t33_g503655 + cse_p11_t33_g503655 + (term.aa_low_i * cse_p528_t3_g513034) + (term.aa_low_j * cse_p528_t3_g513034) + (cse_p63_t1_g504502 * dz) + (-((term.cc_low_i * cse_p528_t7_g513038))) + (-((term.cc_low_j * cse_p528_t7_g513038))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (term.cc_low_i + term.cc_low_j + cse_p12_t27_g503696 + (term.bb_low_i * cse_p528_t4_g513035) + (term.bb_low_j * cse_p528_t4_g513035));
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= term.r2_high)) {
            /* !COND-3 */
            double cse_p1_t3_invr503619 = (1.0 / (rr));
            double cse_p1_t4_invr2503620 = (cse_p1_t3_invr503619 * cse_p1_t3_invr503619);
            double inv_rr = cse_p1_t3_invr503619;
            double inv_rr2 = (inv_rr * inv_rr);
            double di = (rr + (-(term.lj_radius_i)));
            double cse_p14_t7_g503740 = (di * di);
            double cse_p13_t12_g503722 = (-(cse_p14_t7_g503740));
            double cse_p1_t5_invr503621 = (1.0 / (term.lambda2_i));
            double cse_p12_t20_g503689 = (cse_p13_t12_g503722 * cse_p1_t5_invr503621);
            double cse_p11_t36_g503658 = exp(cse_p12_t20_g503689);
            double exp_i = cse_p11_t36_g503658;
            double cse_p12_t24_g503693 = (exp_i * term.pref_i_j);
            double fi = (exp_i * inv_rr2 * term.pref_i_j);
            double hi = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t5_invr503621 * di)));
            double dfi = (fi * hi);
            double dj = (rr + (-(term.lj_radius_j)));
            double cse_p14_t8_g503741 = (dj * dj);
            double cse_p13_t13_g503723 = (-(cse_p14_t8_g503741));
            double cse_p1_t6_invr503622 = (1.0 / (term.lambda2_j));
            double cse_p12_t21_g503690 = (cse_p13_t13_g503723 * cse_p1_t6_invr503622);
            double cse_p11_t37_g503659 = exp(cse_p12_t21_g503690);
            double exp_j = cse_p11_t37_g503659;
            double cse_p12_t25_g503694 = (exp_j * term.pref_j_i);
            double cse_p11_t25_g503647 = (cse_p12_t24_g503693 + cse_p12_t25_g503694);
            double fj = (exp_j * inv_rr2 * term.pref_j_i);
            double hj = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t6_invr503622 * dj)));
            double dfj = (fj * hj);
            double energy = (fi + fj);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p11_t1_g503661 = (-2.00000000000000000e+0 * inv_rr2);
            double cse_p528_t1_g513050 = (cse_p1_t2_invsqrt503618 * dx);
            double cse_p11_t1_g503666 = (cse_p11_t1_g503661 * cse_p11_t36_g503658 * cse_p1_t5_invr503621 * term.pref_i_j);
            double cse_p11_t1_g503667 = (cse_p11_t1_g503661 * cse_p11_t37_g503659 * cse_p1_t6_invr503622 * term.pref_j_i);
            double cse_p11_t1_g503668 = (cse_p11_t25_g503647 * cse_p1_t4_invr2503620 * inv_rr);
            double cse_p12_t1_g503706 = (-2.00000000000000000e+0 * cse_p11_t1_g503668);
            double cse_p63_t1_g504499 = (cse_p11_t1_g503666 * di);
            double cse_p64_t13_g504518 = (cse_p528_t1_g513050 * cse_p63_t1_g504499);
            double cse_p63_t1_g504500 = (cse_p11_t1_g503667 * dj);
            double cse_p64_t14_g504519 = (cse_p528_t1_g513050 * cse_p63_t1_g504500);
            double g_x1 = (cse_p64_t13_g504518 + cse_p64_t14_g504519 + (cse_p12_t1_g503706 * cse_p528_t1_g513050));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p11_t1_g503664 = (cse_p1_t2_invsqrt503618 * dy);
            double cse_p64_t5_g504510 = (cse_p11_t1_g503664 * cse_p63_t1_g504499);
            double cse_p64_t6_g504511 = (cse_p11_t1_g503664 * cse_p63_t1_g504500);
            double g_y1 = (cse_p64_t5_g504510 + cse_p64_t6_g504511 + (cse_p11_t1_g503664 * cse_p12_t1_g503706));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p11_t1_g503665 = (cse_p1_t2_invsqrt503618 * dz);
            double cse_p64_t7_g504512 = (cse_p11_t1_g503665 * cse_p63_t1_g504499);
            double cse_p64_t8_g504513 = (cse_p11_t1_g503665 * cse_p63_t1_g504500);
            double g_z1 = (cse_p64_t7_g504512 + cse_p64_t8_g504513 + (cse_p11_t1_g503665 * cse_p12_t1_g503706));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double cse_p12_t1_g503707 = (2.00000000000000000e+0 * cse_p11_t1_g503668);
            double g_x2 = ((cse_p12_t1_g503707 * cse_p528_t1_g513050) + (-(cse_p64_t13_g504518)) + (-(cse_p64_t14_g504519)));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = ((cse_p11_t1_g503664 * cse_p12_t1_g503707) + (-(cse_p64_t5_g504510)) + (-(cse_p64_t6_g504511)));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = ((cse_p11_t1_g503665 * cse_p12_t1_g503707) + (-(cse_p64_t7_g504512)) + (-(cse_p64_t8_g504513)));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (dfi + dfj);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= term.r3_high)) {
              /* !COND-4 */
              double cse_p12_t1_g503708 = (3.00000000000000000e+0 * term.aa_high_i);
              double cse_p12_t1_g503709 = (3.00000000000000000e+0 * term.aa_high_j);
              double cse_p13_t10_g503720 = (cse_p12_t1_g503708 + cse_p12_t1_g503709);
              double cse_p12_t5_g503674 = (term.bb_high_i * rr);
              double cse_p12_t6_g503675 = (term.bb_high_j * rr);
              double cse_p12_t1_g503705 = (rr * rr);
              double cse_p12_t26_g503695 = (cse_p12_t1_g503705 * cse_p13_t10_g503720);
              double cse_p11_t1_g503660 = (cse_p1_t2_invsqrt503618 * rr);
              double cse_p11_t1_g503669 = (term.bb_high_i * cse_p11_t1_g503660);
              double cse_p12_t13_g503682 = (cse_p11_t1_g503669 * dx);
              double cse_p12_t14_g503683 = (cse_p11_t1_g503669 * dy);
              double cse_p12_t15_g503684 = (cse_p11_t1_g503669 * dz);
              double cse_p11_t26_g503648 = (term.cc_high_i + term.cc_high_j + cse_p12_t26_g503695 + cse_p12_t5_g503674 + cse_p12_t5_g503674 + cse_p12_t6_g503675 + cse_p12_t6_g503675);
              double cse_p11_t28_g503650 = (-(cse_p12_t13_g503682));
              double cse_p11_t29_g503651 = (-(cse_p12_t14_g503683));
              double cse_p11_t30_g503652 = (-(cse_p12_t15_g503684));
              double energy = (term.dd_high_i + term.dd_high_j + (term.cc_high_i * rr) + (term.cc_high_j * rr) + (cse_p12_t1_g503705 * (term.bb_high_i + term.bb_high_j + (term.aa_high_i * rr) + (term.aa_high_j * rr))));
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p528_t1_g513050 = (cse_p1_t2_invsqrt503618 * dx);
              double g_x1 = (cse_p11_t26_g503648 * cse_p528_t1_g513050);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double cse_p11_t1_g503664 = (cse_p1_t2_invsqrt503618 * dy);
              double g_y1 = (cse_p11_t1_g503664 * cse_p11_t26_g503648);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double cse_p11_t1_g503665 = (cse_p1_t2_invsqrt503618 * dz);
              double g_z1 = (cse_p11_t1_g503665 * cse_p11_t26_g503648);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p11_t1_g503662 = (cse_p11_t1_g503660 * rr);
              double cse_p12_t1_g503702 = (cse_p11_t1_g503662 * dx);
              double cse_p13_t1_g503724 = (term.bb_high_j * cse_p11_t1_g503660);
              double cse_p63_t1_g504496 = (-3.00000000000000000e+0 * cse_p12_t1_g503702);
              double cse_p63_t1_g504501 = (-2.00000000000000000e+0 * cse_p13_t1_g503724);
              double g_x2 = (cse_p11_t28_g503650 + cse_p11_t28_g503650 + (term.aa_high_i * cse_p63_t1_g504496) + (term.aa_high_j * cse_p63_t1_g504496) + (cse_p63_t1_g504501 * dx) + (-((term.cc_high_i * cse_p528_t1_g513050))) + (-((term.cc_high_j * cse_p528_t1_g513050))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double cse_p12_t1_g503703 = (cse_p11_t1_g503662 * dy);
              double cse_p63_t1_g504497 = (-3.00000000000000000e+0 * cse_p12_t1_g503703);
              double g_y2 = (cse_p11_t29_g503651 + cse_p11_t29_g503651 + (term.aa_high_i * cse_p63_t1_g504497) + (term.aa_high_j * cse_p63_t1_g504497) + (cse_p63_t1_g504501 * dy) + (-((term.cc_high_i * cse_p11_t1_g503664))) + (-((term.cc_high_j * cse_p11_t1_g503664))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double cse_p12_t1_g503704 = (cse_p11_t1_g503662 * dz);
              double cse_p63_t1_g504498 = (-3.00000000000000000e+0 * cse_p12_t1_g503704);
              double g_z2 = (cse_p11_t30_g503652 + cse_p11_t30_g503652 + (term.aa_high_i * cse_p63_t1_g504498) + (term.aa_high_j * cse_p63_t1_g504498) + (cse_p63_t1_g504501 * dz) + (-((term.cc_high_i * cse_p11_t1_g503665))) + (-((term.cc_high_j * cse_p11_t1_g503665))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double cse_p12_t1_g503701 = (2.00000000000000000e+0 * rr);
              double de_dr = (term.cc_high_i + term.cc_high_j + cse_p12_t26_g503695 + (term.bb_high_i * cse_p12_t1_g503701) + (term.bb_high_j * cse_p12_t1_g503701));
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
            double cse_p61_t1_invr513603 = (1.0 / (term.lambda2_i));
            double exp_i = exp((-((cse_p61_t1_invr513603 * di * di))));
            double fi = (exp_i * inv_rr2 * term.pref_i_j);
            double hi = (-2.00000000000000000e+0 * (inv_rr + (cse_p61_t1_invr513603 * di)));
            double dj = (rr + (-(term.lj_radius_j)));
            double cse_p61_t2_invr513604 = (1.0 / (term.lambda2_j));
            double exp_j = exp((-((cse_p61_t2_invr513604 * dj * dj))));
            double fj = (exp_j * inv_rr2 * term.pref_j_i);
            double hj = (-2.00000000000000000e+0 * (inv_rr + (cse_p61_t2_invr513604 * dj)));
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
      double cse_p11_t53_g518613 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t54_g518614 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t55_g518615 = (dz * dz);
      double r2 = (cse_p11_t53_g518613 + cse_p11_t54_g518614 + cse_p11_t55_g518615);
      double cse_p1_t1_sqrt518555 = sqrt(r2);
      double cse_p1_t2_invsqrt518556 = (1.0 / (cse_p1_t1_sqrt518555));
      double rr = cse_p1_t1_sqrt518555;
      if ((rr <= term.r0_low)) {
        /* !COND-1 */
        double cse_p60_t43_g522662 = (rr * rr);
        double cse_p60_t71_g522690 = (-(cse_p11_t53_g518613));
        double cse_p60_t72_g522691 = (-(cse_p11_t54_g518614));
        double cse_p60_t73_g522692 = (-(cse_p11_t55_g518615));
        double cse_p60_t75_g522694 = (-(dx));
        double cse_p60_t76_g522695 = (-(dy));
        double cse_p60_t77_g522696 = (-(dz));
        double cse_p565_t1_invr558250 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr558251 = (1.0 / (rr));
        double cse_p60_t80_g522699 = (rr * rr * rr);
        double cse_p15_t75_g519093 = (-(cse_p60_t43_g522662));
        double cse_p14_t82_g518991 = (cse_p11_t53_g518613 + cse_p15_t75_g519093);
        double cse_p14_t84_g518993 = (cse_p11_t54_g518614 + cse_p15_t75_g519093);
        double cse_p14_t86_g518995 = (cse_p11_t55_g518615 + cse_p15_t75_g519093);
        double cse_p13_t109_g518894 = (cse_p60_t43_g522662 + cse_p60_t71_g522690);
        double cse_p13_t110_g518895 = (cse_p60_t43_g522662 + cse_p60_t72_g522691);
        double cse_p13_t111_g518896 = (cse_p60_t43_g522662 + cse_p60_t73_g522692);
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
        double cse_p60_t42_g522661 = (cse_p565_t1_invr558250 * de_dr);
        double cse_p13_t67_g518852 = (cse_p14_t82_g518991 * cse_p60_t42_g522661);
        double cse_p13_t68_g518853 = (cse_p14_t84_g518993 * cse_p60_t42_g522661);
        double cse_p13_t69_g518854 = (cse_p14_t86_g518995 * cse_p60_t42_g522661);
        double cse_p13_t70_g518855 = (cse_p60_t42_g522661 * dx * dy);
        double cse_p13_t71_g518856 = (cse_p60_t42_g522661 * cse_p60_t75_g522694 * dy);
        double cse_p12_t27_g518709 = (cse_p13_t109_g518894 * cse_p60_t42_g522661);
        double cse_p12_t28_g518710 = (cse_p13_t110_g518895 * cse_p60_t42_g522661);
        double cse_p12_t29_g518711 = (cse_p13_t111_g518896 * cse_p60_t42_g522661);
        double cse_p12_t32_g518714 = (cse_p60_t42_g522661 * dz);
        double cse_p13_t57_g518842 = (cse_p12_t32_g518714 * cse_p60_t75_g522694);
        double cse_p13_t58_g518843 = (cse_p12_t32_g518714 * cse_p60_t76_g522695);
        double cse_p13_t59_g518844 = (cse_p12_t32_g518714 * dx);
        double cse_p13_t60_g518845 = (cse_p12_t32_g518714 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p60_t41_g522660 = (cse_p565_t2_invr558251 * cse_p565_t2_invr558251 * d2e_dr2);
        double cse_p13_t61_g518846 = (cse_p60_t41_g522660 * dx * dx);
        double cse_p13_t62_g518847 = (cse_p60_t41_g522660 * dx * dy);
        double cse_p13_t63_g518848 = (cse_p60_t41_g522660 * dx * dz);
        double cse_p13_t64_g518849 = (cse_p60_t41_g522660 * dy * dy);
        double cse_p13_t65_g518850 = (cse_p60_t41_g522660 * dy * dz);
        double cse_p13_t66_g518851 = (cse_p60_t41_g522660 * dz * dz);
        double cse_p12_t15_g518697 = (cse_p60_t41_g522660 * cse_p60_t75_g522694);
        double cse_p13_t42_g518827 = (cse_p12_t15_g518697 * cse_p60_t75_g522694);
        double cse_p13_t43_g518828 = (cse_p12_t15_g518697 * cse_p60_t76_g522695);
        double cse_p13_t44_g518829 = (cse_p12_t15_g518697 * cse_p60_t77_g522696);
        double cse_p13_t45_g518830 = (cse_p12_t15_g518697 * dx);
        double cse_p13_t46_g518831 = (cse_p12_t15_g518697 * dy);
        double cse_p13_t47_g518832 = (cse_p12_t15_g518697 * dz);
        double cse_p12_t16_g518698 = (cse_p60_t41_g522660 * cse_p60_t76_g522695);
        double cse_p13_t48_g518833 = (cse_p12_t16_g518698 * cse_p60_t76_g522695);
        double cse_p13_t49_g518834 = (cse_p12_t16_g518698 * cse_p60_t77_g522696);
        double cse_p13_t50_g518835 = (cse_p12_t16_g518698 * dx);
        double cse_p13_t51_g518836 = (cse_p12_t16_g518698 * dy);
        double cse_p13_t52_g518837 = (cse_p12_t16_g518698 * dz);
        double cse_p12_t17_g518699 = (cse_p60_t41_g522660 * cse_p60_t77_g522696);
        double cse_p13_t53_g518838 = (cse_p12_t17_g518699 * cse_p60_t77_g522696);
        double cse_p13_t54_g518839 = (cse_p12_t17_g518699 * dx);
        double cse_p13_t55_g518840 = (cse_p12_t17_g518699 * dy);
        double cse_p13_t56_g518841 = (cse_p12_t17_g518699 * dz);
        double cse_p12_t68_g518750 = (cse_p13_t62_g518847 + cse_p13_t71_g518856);
        double cse_p12_t69_g518751 = (cse_p13_t57_g518842 + cse_p13_t63_g518848);
        double cse_p12_t70_g518752 = (cse_p13_t58_g518843 + cse_p13_t65_g518850);
        double cse_p12_t86_g518768 = (cse_p12_t28_g518710 + cse_p13_t64_g518849);
        double cse_p12_t88_g518770 = (cse_p12_t29_g518711 + cse_p13_t66_g518851);
        double cse_p12_t90_g518772 = (cse_p12_t27_g518709 + cse_p13_t61_g518846);
        double cse_p12_t63_g518745 = (cse_p13_t45_g518830 + cse_p13_t67_g518852);
        double cse_p12_t71_g518753 = (cse_p13_t46_g518831 + cse_p13_t70_g518855);
        double cse_p12_t73_g518755 = (cse_p13_t43_g518828 + cse_p13_t71_g518856);
        double cse_p12_t74_g518756 = (cse_p13_t44_g518829 + cse_p13_t57_g518842);
        double cse_p12_t75_g518757 = (cse_p13_t47_g518832 + cse_p13_t59_g518844);
        double cse_p12_t91_g518773 = (cse_p12_t27_g518709 + cse_p13_t42_g518827);
        double cse_p12_t64_g518746 = (cse_p13_t51_g518836 + cse_p13_t68_g518853);
        double cse_p12_t72_g518754 = (cse_p13_t50_g518835 + cse_p13_t70_g518855);
        double cse_p12_t76_g518758 = (cse_p13_t49_g518834 + cse_p13_t58_g518843);
        double cse_p12_t77_g518759 = (cse_p13_t52_g518837 + cse_p13_t60_g518845);
        double cse_p12_t87_g518769 = (cse_p12_t28_g518710 + cse_p13_t48_g518833);
        double cse_p12_t65_g518747 = (cse_p13_t56_g518841 + cse_p13_t69_g518854);
        double cse_p12_t78_g518760 = (cse_p13_t54_g518839 + cse_p13_t59_g518844);
        double cse_p12_t79_g518761 = (cse_p13_t55_g518840 + cse_p13_t60_g518845);
        double cse_p12_t89_g518771 = (cse_p12_t29_g518711 + cse_p13_t53_g518838);
        double h_x1_x1 = cse_p12_t90_g518772;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t68_g518750;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t69_g518751;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t63_g518745;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t72_g518754;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t78_g518760;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t86_g518768;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t70_g518752;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t71_g518753;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t64_g518746;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t79_g518761;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t88_g518770;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t75_g518757;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t77_g518759;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t65_g518747;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t91_g518773;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t73_g518755;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t74_g518756;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t87_g518769;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t76_g518758;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t89_g518771;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.r1_low)) {
          /* !COND-2 */
          double cse_p60_t38_g522657 = (cse_p1_t2_invsqrt518556 * dx);
          double cse_p60_t39_g522658 = (cse_p1_t2_invsqrt518556 * dy);
          double cse_p60_t40_g522659 = (cse_p1_t2_invsqrt518556 * dz);
          double cse_p14_t3_g518912 = (3.00000000000000000e+0 * term.aa_low_i);
          double cse_p14_t4_g518913 = (3.00000000000000000e+0 * term.aa_low_j);
          double cse_p13_t81_g518866 = (cse_p14_t3_g518912 + cse_p14_t4_g518913);
          double cse_p12_t10_g518692 = (term.bb_low_i * cse_p60_t38_g522657 * rr);
          double cse_p12_t11_g518693 = (term.bb_low_i * cse_p60_t39_g522658 * rr);
          double cse_p12_t12_g518694 = (term.bb_low_i * cse_p60_t40_g522659 * rr);
          double cse_p12_t13_g518695 = (term.bb_low_i * rr);
          double cse_p12_t14_g518696 = (term.bb_low_j * rr);
          double cse_p11_t1_g518675 = (rr * rr);
          double cse_p60_t59_g522678 = (cse_p11_t53_g518613 + (-(cse_p11_t1_g518675)));
          double cse_p60_t60_g522679 = (cse_p11_t54_g518614 + (-(cse_p11_t1_g518675)));
          double cse_p60_t61_g522680 = (cse_p11_t55_g518615 + (-(cse_p11_t1_g518675)));
          double cse_p60_t62_g522681 = (cse_p11_t1_g518675 + (-(cse_p11_t53_g518613)));
          double cse_p60_t63_g522682 = (cse_p11_t1_g518675 + (-(cse_p11_t54_g518614)));
          double cse_p60_t64_g522683 = (cse_p11_t1_g518675 + (-(cse_p11_t55_g518615)));
          double cse_p60_t74_g522693 = (-(cse_p11_t1_g518675));
          double cse_p12_t34_g518716 = (cse_p11_t1_g518675 * cse_p13_t81_g518866);
          double cse_p11_t92_g518652 = (term.cc_low_i + term.cc_low_j + cse_p12_t13_g518695 + cse_p12_t13_g518695 + cse_p12_t14_g518696 + cse_p12_t14_g518696 + cse_p12_t34_g518716);
          double cse_p11_t96_g518656 = (-(cse_p12_t10_g518692));
          double cse_p11_t97_g518657 = (-(cse_p12_t11_g518693));
          double cse_p11_t98_g518658 = (-(cse_p12_t12_g518694));
          double energy = (term.dd_low_i + term.dd_low_j + (term.cc_low_i * rr) + (term.cc_low_j * rr) + (cse_p11_t1_g518675 * (term.bb_low_i + term.bb_low_j + (term.aa_low_i * rr) + (term.aa_low_j * rr))));
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (cse_p11_t92_g518652 * cse_p60_t38_g522657);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t92_g518652 * cse_p60_t39_g522658);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t92_g518652 * cse_p60_t40_g522659);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p11_t96_g518656 + cse_p11_t96_g518656 + (-2.00000000000000000e+0 * term.bb_low_j * cse_p60_t38_g522657 * rr) + (cse_p11_t1_g518675 * cse_p60_t38_g522657 * ((-3.00000000000000000e+0 * term.aa_low_i) + (-3.00000000000000000e+0 * term.aa_low_j))) + (-((term.cc_low_i * cse_p60_t38_g522657))) + (-((term.cc_low_j * cse_p60_t38_g522657))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p11_t97_g518657 + cse_p11_t97_g518657 + (-2.00000000000000000e+0 * term.bb_low_j * cse_p60_t39_g522658 * rr) + (cse_p11_t1_g518675 * cse_p60_t39_g522658 * ((-3.00000000000000000e+0 * term.aa_low_i) + (-3.00000000000000000e+0 * term.aa_low_j))) + (-((term.cc_low_i * cse_p60_t39_g522658))) + (-((term.cc_low_j * cse_p60_t39_g522658))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p11_t98_g518658 + cse_p11_t98_g518658 + (-2.00000000000000000e+0 * term.bb_low_j * cse_p60_t40_g522659 * rr) + (cse_p11_t1_g518675 * cse_p60_t40_g522659 * ((-3.00000000000000000e+0 * term.aa_low_i) + (-3.00000000000000000e+0 * term.aa_low_j))) + (-((term.cc_low_i * cse_p60_t40_g522659))) + (-((term.cc_low_j * cse_p60_t40_g522659))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (term.cc_low_i + term.cc_low_j + cse_p12_t34_g518716 + (2.00000000000000000e+0 * term.bb_low_i * rr) + (2.00000000000000000e+0 * term.bb_low_j * rr));
          double d2e_dr2 = ((2.00000000000000000e+0 * term.bb_low_i) + (2.00000000000000000e+0 * term.bb_low_j) + (6.00000000000000000e+0 * term.aa_low_i * rr) + (6.00000000000000000e+0 * term.aa_low_j * rr));
          double cse_p565_t2_invr558251 = (1.0 / (rr));
          double cse_p11_t1_g518673 = (cse_p565_t2_invr558251 * cse_p565_t2_invr558251 * d2e_dr2);
          double cse_p60_t1_g522620 = (cse_p11_t1_g518673 * (-(dx)));
          double cse_p60_t2_g522621 = (cse_p11_t1_g518673 * (-(dy)));
          double cse_p60_t3_g522622 = (cse_p11_t1_g518673 * (-(dz)));
          double cse_p60_t4_g522623 = (cse_p11_t1_g518673 * dx * dx);
          double cse_p60_t5_g522624 = (cse_p11_t1_g518673 * dx * dy);
          double cse_p60_t6_g522625 = (cse_p11_t1_g518673 * dx * dz);
          double cse_p60_t7_g522626 = (cse_p11_t1_g518673 * dy * dy);
          double cse_p60_t8_g522627 = (cse_p11_t1_g518673 * dy * dz);
          double cse_p60_t9_g522628 = (cse_p11_t1_g518673 * dz * dz);
          double cse_p14_t46_g518955 = (cse_p60_t1_g522620 * (-(dx)));
          double cse_p14_t47_g518956 = (cse_p60_t1_g522620 * (-(dy)));
          double cse_p14_t48_g518957 = (cse_p60_t1_g522620 * (-(dz)));
          double cse_p14_t49_g518958 = (cse_p60_t1_g522620 * dx);
          double cse_p14_t50_g518959 = (cse_p60_t1_g522620 * dy);
          double cse_p14_t51_g518960 = (cse_p60_t1_g522620 * dz);
          double cse_p14_t52_g518961 = (cse_p60_t2_g522621 * (-(dy)));
          double cse_p14_t53_g518962 = (cse_p60_t2_g522621 * (-(dz)));
          double cse_p14_t54_g518963 = (cse_p60_t2_g522621 * dx);
          double cse_p14_t55_g518964 = (cse_p60_t2_g522621 * dy);
          double cse_p14_t56_g518965 = (cse_p60_t2_g522621 * dz);
          double cse_p14_t57_g518966 = (cse_p60_t3_g522622 * (-(dz)));
          double cse_p14_t58_g518967 = (cse_p60_t3_g522622 * dx);
          double cse_p14_t59_g518968 = (cse_p60_t3_g522622 * dy);
          double cse_p14_t60_g518969 = (cse_p60_t3_g522622 * dz);
          double cse_p565_t1_invr558250 = (1.0 / ((rr * rr * rr)));
          double cse_p11_t1_g518674 = (cse_p565_t1_invr558250 * de_dr);
          double cse_p60_t10_g522629 = (cse_p11_t1_g518674 * (cse_p11_t53_g518613 + cse_p60_t74_g522693));
          double cse_p60_t11_g522630 = (cse_p11_t1_g518674 * (cse_p11_t54_g518614 + cse_p60_t74_g522693));
          double cse_p60_t12_g522631 = (cse_p11_t1_g518674 * (cse_p11_t55_g518615 + cse_p60_t74_g522693));
          double cse_p60_t13_g522632 = (cse_p11_t1_g518674 * cse_p60_t62_g522681);
          double cse_p60_t14_g522633 = (cse_p11_t1_g518674 * cse_p60_t63_g522682);
          double cse_p60_t15_g522634 = (cse_p11_t1_g518674 * cse_p60_t64_g522683);
          double cse_p60_t16_g522635 = (cse_p11_t1_g518674 * dx * dy);
          double cse_p60_t17_g522636 = (cse_p11_t1_g518674 * dy * (-(dx)));
          double cse_p60_t18_g522637 = (cse_p11_t1_g518674 * dz);
          double cse_p60_t47_g522666 = (cse_p60_t17_g522636 + cse_p60_t5_g522624);
          double cse_p14_t42_g518951 = (cse_p60_t18_g522637 * (-(dx)));
          double cse_p14_t43_g518952 = (cse_p60_t18_g522637 * (-(dy)));
          double cse_p14_t44_g518953 = (cse_p60_t18_g522637 * dx);
          double cse_p14_t45_g518954 = (cse_p60_t18_g522637 * dy);
          double cse_p13_t107_g518892 = (cse_p60_t13_g522632 + cse_p60_t4_g522623);
          double h_x1_x1 = cse_p13_t107_g518892;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t47_g522666;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t83_g518868 = (cse_p14_t42_g518951 + cse_p60_t6_g522625);
          double cse_p13_t84_g518869 = (cse_p14_t43_g518952 + cse_p60_t8_g522627);
          double h_x1_z1 = cse_p13_t83_g518868;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t77_g518862 = (cse_p14_t49_g518958 + cse_p60_t10_g522629);
          double cse_p13_t85_g518870 = (cse_p14_t50_g518959 + cse_p60_t16_g522635);
          double cse_p13_t87_g518872 = (cse_p14_t47_g518956 + cse_p60_t17_g522636);
          double cse_p13_t88_g518873 = (cse_p14_t42_g518951 + cse_p14_t48_g518957);
          double cse_p13_t89_g518874 = (cse_p14_t44_g518953 + cse_p14_t51_g518960);
          double cse_p13_t108_g518893 = (cse_p14_t46_g518955 + cse_p60_t13_g522632);
          double h_x1_x2 = cse_p13_t77_g518862;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t78_g518863 = (cse_p14_t55_g518964 + cse_p60_t11_g522630);
          double cse_p13_t86_g518871 = (cse_p14_t54_g518963 + cse_p60_t16_g522635);
          double cse_p13_t90_g518875 = (cse_p14_t43_g518952 + cse_p14_t53_g518962);
          double cse_p13_t91_g518876 = (cse_p14_t45_g518954 + cse_p14_t56_g518965);
          double h_x1_y2 = cse_p13_t86_g518871;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t79_g518864 = (cse_p14_t60_g518969 + cse_p60_t12_g522631);
          double cse_p13_t92_g518877 = (cse_p14_t44_g518953 + cse_p14_t58_g518967);
          double cse_p13_t93_g518878 = (cse_p14_t45_g518954 + cse_p14_t59_g518968);
          double h_x1_z2 = cse_p13_t92_g518877;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t103_g518888 = (cse_p60_t14_g522633 + cse_p60_t7_g522626);
          double cse_p13_t104_g518889 = (cse_p14_t52_g518961 + cse_p60_t14_g522633);
          double h_y1_y1 = cse_p13_t103_g518888;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t84_g518869;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t85_g518870;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t78_g518863;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t93_g518878;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t105_g518890 = (cse_p60_t15_g522634 + cse_p60_t9_g522628);
          double cse_p13_t106_g518891 = (cse_p14_t57_g518966 + cse_p60_t15_g522634);
          double h_z1_z1 = cse_p13_t105_g518890;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t89_g518874;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t91_g518876;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t79_g518864;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t108_g518893;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t87_g518872;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t88_g518873;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t104_g518889;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t90_g518875;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t106_g518891;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= term.r2_high)) {
            /* !COND-3 */
            double cse_p565_t2_invr558251 = (1.0 / (rr));
            double cse_p1_t4_invr2518558 = (cse_p565_t2_invr558251 * cse_p565_t2_invr558251);
            double inv_rr = cse_p565_t2_invr558251;
            double inv_rr2 = (inv_rr * inv_rr);
            double cse_p11_t1_g518561 = (2.00000000000000000e+0 * inv_rr2);
            double di = (rr + (-(term.lj_radius_i)));
            double cse_p14_t63_g518972 = (di * di);
            double cse_p13_t112_g518897 = (-(cse_p14_t63_g518972));
            double cse_p1_t5_invr518559 = (1.0 / (term.lambda2_i));
            double cse_p12_t54_g518736 = (cse_p13_t112_g518897 * cse_p1_t5_invr518559);
            double cse_p11_t108_g518668 = exp(cse_p12_t54_g518736);
            double exp_i = cse_p11_t108_g518668;
            double cse_p12_t60_g518742 = (exp_i * term.pref_i_j);
            double fi = (exp_i * inv_rr2 * term.pref_i_j);
            double hi = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t5_invr518559 * di)));
            double dfi = (fi * hi);
            double d2fi = (fi * (cse_p11_t1_g518561 + (-2.00000000000000000e+0 * cse_p1_t5_invr518559) + (hi * hi)));
            double dj = (rr + (-(term.lj_radius_j)));
            double cse_p14_t64_g518973 = (dj * dj);
            double cse_p13_t113_g518898 = (-(cse_p14_t64_g518973));
            double cse_p1_t6_invr518560 = (1.0 / (term.lambda2_j));
            double cse_p12_t55_g518737 = (cse_p13_t113_g518898 * cse_p1_t6_invr518560);
            double cse_p11_t109_g518669 = exp(cse_p12_t55_g518737);
            double exp_j = cse_p11_t109_g518669;
            double cse_p12_t61_g518743 = (exp_j * term.pref_j_i);
            double cse_p11_t87_g518647 = (cse_p12_t60_g518742 + cse_p12_t61_g518743);
            double fj = (exp_j * inv_rr2 * term.pref_j_i);
            double hj = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t6_invr518560 * dj)));
            double dfj = (fj * hj);
            double d2fj = (fj * (cse_p11_t1_g518561 + (-2.00000000000000000e+0 * cse_p1_t6_invr518560) + (hj * hj)));
            double energy = (fi + fj);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p14_t1_g519016 = (cse_p1_t2_invsqrt518556 * dx);
            double g_x1 = ((-2.00000000000000000e+0 * cse_p11_t87_g518647 * cse_p14_t1_g519016 * cse_p1_t4_invr2518558 * inv_rr) + (cse_p14_t1_g519016 * inv_rr2 * ((-2.00000000000000000e+0 * cse_p11_t108_g518668 * cse_p1_t5_invr518559 * di * term.pref_i_j) + (-2.00000000000000000e+0 * cse_p11_t109_g518669 * cse_p1_t6_invr518560 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p14_t1_g519017 = (cse_p1_t2_invsqrt518556 * dy);
            double g_y1 = ((-2.00000000000000000e+0 * cse_p11_t87_g518647 * cse_p14_t1_g519017 * cse_p1_t4_invr2518558 * inv_rr) + (cse_p14_t1_g519017 * inv_rr2 * ((-2.00000000000000000e+0 * cse_p11_t108_g518668 * cse_p1_t5_invr518559 * di * term.pref_i_j) + (-2.00000000000000000e+0 * cse_p11_t109_g518669 * cse_p1_t6_invr518560 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p14_t1_g519018 = (cse_p1_t2_invsqrt518556 * dz);
            double g_z1 = ((-2.00000000000000000e+0 * cse_p11_t87_g518647 * cse_p14_t1_g519018 * cse_p1_t4_invr2518558 * inv_rr) + (cse_p14_t1_g519018 * inv_rr2 * ((-2.00000000000000000e+0 * cse_p11_t108_g518668 * cse_p1_t5_invr518559 * di * term.pref_i_j) + (-2.00000000000000000e+0 * cse_p11_t109_g518669 * cse_p1_t6_invr518560 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = ((2.00000000000000000e+0 * cse_p11_t87_g518647 * cse_p14_t1_g519016 * cse_p1_t4_invr2518558 * inv_rr) + (cse_p14_t1_g519016 * inv_rr2 * ((2.00000000000000000e+0 * cse_p11_t108_g518668 * cse_p1_t5_invr518559 * di * term.pref_i_j) + (2.00000000000000000e+0 * cse_p11_t109_g518669 * cse_p1_t6_invr518560 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = ((2.00000000000000000e+0 * cse_p11_t87_g518647 * cse_p14_t1_g519017 * cse_p1_t4_invr2518558 * inv_rr) + (cse_p14_t1_g519017 * inv_rr2 * ((2.00000000000000000e+0 * cse_p11_t108_g518668 * cse_p1_t5_invr518559 * di * term.pref_i_j) + (2.00000000000000000e+0 * cse_p11_t109_g518669 * cse_p1_t6_invr518560 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = ((2.00000000000000000e+0 * cse_p11_t87_g518647 * cse_p14_t1_g519018 * cse_p1_t4_invr2518558 * inv_rr) + (cse_p14_t1_g519018 * inv_rr2 * ((2.00000000000000000e+0 * cse_p11_t108_g518668 * cse_p1_t5_invr518559 * di * term.pref_i_j) + (2.00000000000000000e+0 * cse_p11_t109_g518669 * cse_p1_t6_invr518560 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (dfi + dfj);
            double d2e_dr2 = (d2fi + d2fj);
            double cse_p11_t1_g518673 = (cse_p565_t2_invr558251 * cse_p565_t2_invr558251 * d2e_dr2);
            double cse_p565_t1_invr558250 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g518674 = (cse_p565_t1_invr558250 * de_dr);
            double cse_p11_t1_g518675 = (rr * rr);
            double cse_p11_t1_g518682 = (cse_p11_t1_g518674 * (cse_p11_t1_g518675 + (-(cse_p11_t53_g518613))));
            double cse_p60_t69_g522688 = (cse_p11_t1_g518682 + (cse_p11_t1_g518673 * dx * dx));
            double h_x1_x1 = cse_p60_t69_g522688;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = ((cse_p11_t1_g518673 * dx * dy) + (cse_p11_t1_g518674 * dy * (-(dx))));
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p11_t1_g518678 = (cse_p11_t1_g518674 * dz);
            double cse_p60_t30_g522649 = (cse_p11_t1_g518678 * (-(dx)));
            double cse_p60_t31_g522650 = (cse_p11_t1_g518678 * (-(dy)));
            double cse_p60_t32_g522651 = (cse_p11_t1_g518678 * dx);
            double cse_p60_t33_g522652 = (cse_p11_t1_g518678 * dy);
            double cse_p60_t48_g522667 = (cse_p60_t30_g522649 + (cse_p11_t1_g518673 * dx * dz));
            double cse_p60_t49_g522668 = (cse_p60_t31_g522650 + (cse_p11_t1_g518673 * dy * dz));
            double h_x1_z1 = cse_p60_t48_g522667;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p11_t1_g518676 = (cse_p11_t1_g518673 * (-(dx)));
            double cse_p60_t19_g522638 = (cse_p11_t1_g518676 * (-(dx)));
            double cse_p60_t20_g522639 = (cse_p11_t1_g518676 * (-(dy)));
            double cse_p60_t21_g522640 = (cse_p11_t1_g518676 * (-(dz)));
            double cse_p60_t22_g522641 = (cse_p11_t1_g518676 * dx);
            double cse_p60_t23_g522642 = (cse_p11_t1_g518676 * dy);
            double cse_p60_t24_g522643 = (cse_p11_t1_g518676 * dz);
            double cse_p60_t44_g522663 = (cse_p60_t22_g522641 + (cse_p11_t1_g518674 * (cse_p11_t53_g518613 + (-(cse_p11_t1_g518675)))));
            double cse_p60_t50_g522669 = (cse_p60_t23_g522642 + (cse_p11_t1_g518674 * dx * dy));
            double cse_p60_t52_g522671 = (cse_p60_t20_g522639 + (cse_p11_t1_g518674 * dy * (-(dx))));
            double cse_p60_t53_g522672 = (cse_p60_t21_g522640 + cse_p60_t30_g522649);
            double cse_p60_t54_g522673 = (cse_p60_t24_g522643 + cse_p60_t32_g522651);
            double cse_p60_t70_g522689 = (cse_p11_t1_g518682 + cse_p60_t19_g522638);
            double h_x1_x2 = cse_p60_t44_g522663;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p11_t1_g518677 = (cse_p11_t1_g518673 * (-(dy)));
            double cse_p60_t25_g522644 = (cse_p11_t1_g518677 * (-(dy)));
            double cse_p60_t26_g522645 = (cse_p11_t1_g518677 * (-(dz)));
            double cse_p60_t27_g522646 = (cse_p11_t1_g518677 * dx);
            double cse_p60_t28_g522647 = (cse_p11_t1_g518677 * dy);
            double cse_p60_t29_g522648 = (cse_p11_t1_g518677 * dz);
            double cse_p60_t45_g522664 = (cse_p60_t28_g522647 + (cse_p11_t1_g518674 * (cse_p11_t54_g518614 + (-(cse_p11_t1_g518675)))));
            double cse_p60_t51_g522670 = (cse_p60_t27_g522646 + (cse_p11_t1_g518674 * dx * dy));
            double cse_p60_t55_g522674 = (cse_p60_t26_g522645 + cse_p60_t31_g522650);
            double cse_p60_t56_g522675 = (cse_p60_t29_g522648 + cse_p60_t33_g522652);
            double h_x1_y2 = cse_p60_t51_g522670;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p11_t1_g518679 = (cse_p11_t1_g518673 * (-(dz)));
            double cse_p60_t34_g522653 = (cse_p11_t1_g518679 * (-(dz)));
            double cse_p60_t35_g522654 = (cse_p11_t1_g518679 * dx);
            double cse_p60_t36_g522655 = (cse_p11_t1_g518679 * dy);
            double cse_p60_t37_g522656 = (cse_p11_t1_g518679 * dz);
            double cse_p60_t46_g522665 = (cse_p60_t37_g522656 + (cse_p11_t1_g518674 * (cse_p11_t55_g518615 + (-(cse_p11_t1_g518675)))));
            double cse_p60_t57_g522676 = (cse_p60_t32_g522651 + cse_p60_t35_g522654);
            double cse_p60_t58_g522677 = (cse_p60_t33_g522652 + cse_p60_t36_g522655);
            double h_x1_z2 = cse_p60_t57_g522676;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p11_t1_g518680 = (cse_p11_t1_g518674 * (cse_p11_t1_g518675 + (-(cse_p11_t54_g518614))));
            double cse_p60_t65_g522684 = (cse_p11_t1_g518680 + (cse_p11_t1_g518673 * dy * dy));
            double cse_p60_t66_g522685 = (cse_p11_t1_g518680 + cse_p60_t25_g522644);
            double h_y1_y1 = cse_p60_t65_g522684;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t49_g522668;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t50_g522669;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t45_g522664;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t58_g522677;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p11_t1_g518681 = (cse_p11_t1_g518674 * (cse_p11_t1_g518675 + (-(cse_p11_t55_g518615))));
            double cse_p60_t67_g522686 = (cse_p11_t1_g518681 + (cse_p11_t1_g518673 * dz * dz));
            double cse_p60_t68_g522687 = (cse_p11_t1_g518681 + cse_p60_t34_g522653);
            double h_z1_z1 = cse_p60_t67_g522686;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t54_g522673;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t56_g522675;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t46_g522665;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t70_g522689;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t52_g522671;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t53_g522672;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t66_g522685;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t55_g522674;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t68_g522687;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= term.r3_high)) {
              /* !COND-4 */
              double cse_p14_t1_g518910 = (3.00000000000000000e+0 * term.aa_high_i);
              double cse_p14_t2_g518911 = (3.00000000000000000e+0 * term.aa_high_j);
              double cse_p13_t80_g518865 = (cse_p14_t1_g518910 + cse_p14_t2_g518911);
              double cse_p14_t1_g519016 = (cse_p1_t2_invsqrt518556 * dx);
              double cse_p12_t5_g518687 = (term.bb_high_i * cse_p14_t1_g519016 * rr);
              double cse_p14_t1_g519017 = (cse_p1_t2_invsqrt518556 * dy);
              double cse_p12_t6_g518688 = (term.bb_high_i * cse_p14_t1_g519017 * rr);
              double cse_p14_t1_g519018 = (cse_p1_t2_invsqrt518556 * dz);
              double cse_p12_t7_g518689 = (term.bb_high_i * cse_p14_t1_g519018 * rr);
              double cse_p12_t8_g518690 = (term.bb_high_i * rr);
              double cse_p12_t9_g518691 = (term.bb_high_j * rr);
              double cse_p11_t1_g518675 = (rr * rr);
              double cse_p12_t33_g518715 = (cse_p11_t1_g518675 * cse_p13_t80_g518865);
              double cse_p11_t91_g518651 = (term.cc_high_i + term.cc_high_j + cse_p12_t33_g518715 + cse_p12_t8_g518690 + cse_p12_t8_g518690 + cse_p12_t9_g518691 + cse_p12_t9_g518691);
              double cse_p11_t93_g518653 = (-(cse_p12_t5_g518687));
              double cse_p11_t94_g518654 = (-(cse_p12_t6_g518688));
              double cse_p11_t95_g518655 = (-(cse_p12_t7_g518689));
              double energy = (term.dd_high_i + term.dd_high_j + (term.cc_high_i * rr) + (term.cc_high_j * rr) + (cse_p11_t1_g518675 * (term.bb_high_i + term.bb_high_j + (term.aa_high_i * rr) + (term.aa_high_j * rr))));
              energy_added += energy;
              *energy_accumulate += energy;
              double g_x1 = (cse_p11_t91_g518651 * cse_p14_t1_g519016);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t91_g518651 * cse_p14_t1_g519017);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t91_g518651 * cse_p14_t1_g519018);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = (cse_p11_t93_g518653 + cse_p11_t93_g518653 + (-2.00000000000000000e+0 * term.bb_high_j * cse_p14_t1_g519016 * rr) + (cse_p11_t1_g518675 * cse_p14_t1_g519016 * ((-3.00000000000000000e+0 * term.aa_high_i) + (-3.00000000000000000e+0 * term.aa_high_j))) + (-((term.cc_high_i * cse_p14_t1_g519016))) + (-((term.cc_high_j * cse_p14_t1_g519016))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (cse_p11_t94_g518654 + cse_p11_t94_g518654 + (-2.00000000000000000e+0 * term.bb_high_j * cse_p14_t1_g519017 * rr) + (cse_p11_t1_g518675 * cse_p14_t1_g519017 * ((-3.00000000000000000e+0 * term.aa_high_i) + (-3.00000000000000000e+0 * term.aa_high_j))) + (-((term.cc_high_i * cse_p14_t1_g519017))) + (-((term.cc_high_j * cse_p14_t1_g519017))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (cse_p11_t95_g518655 + cse_p11_t95_g518655 + (-2.00000000000000000e+0 * term.bb_high_j * cse_p14_t1_g519018 * rr) + (cse_p11_t1_g518675 * cse_p14_t1_g519018 * ((-3.00000000000000000e+0 * term.aa_high_i) + (-3.00000000000000000e+0 * term.aa_high_j))) + (-((term.cc_high_i * cse_p14_t1_g519018))) + (-((term.cc_high_j * cse_p14_t1_g519018))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double de_dr = (term.cc_high_i + term.cc_high_j + cse_p12_t33_g518715 + (2.00000000000000000e+0 * term.bb_high_i * rr) + (2.00000000000000000e+0 * term.bb_high_j * rr));
              double d2e_dr2 = ((2.00000000000000000e+0 * term.bb_high_i) + (2.00000000000000000e+0 * term.bb_high_j) + (6.00000000000000000e+0 * term.aa_high_i * rr) + (6.00000000000000000e+0 * term.aa_high_j * rr));
              double cse_p565_t2_invr558251 = (1.0 / (rr));
              double cse_p11_t1_g518673 = (cse_p565_t2_invr558251 * cse_p565_t2_invr558251 * d2e_dr2);
              double cse_p565_t1_invr558250 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g518674 = (cse_p565_t1_invr558250 * de_dr);
              double cse_p11_t1_g518682 = (cse_p11_t1_g518674 * (cse_p11_t1_g518675 + (-(cse_p11_t53_g518613))));
              double h_x1_x1 = (cse_p11_t1_g518682 + (cse_p11_t1_g518673 * dx * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double h_x1_y1 = ((cse_p11_t1_g518673 * dx * dy) + (cse_p11_t1_g518674 * dy * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g518678 = (cse_p11_t1_g518674 * dz);
              double h_x1_z1 = ((cse_p11_t1_g518673 * dx * dz) + (cse_p11_t1_g518678 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g518676 = (cse_p11_t1_g518673 * (-(dx)));
              double h_x1_x2 = ((cse_p11_t1_g518674 * (cse_p11_t53_g518613 + (-(cse_p11_t1_g518675)))) + (cse_p11_t1_g518676 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g518677 = (cse_p11_t1_g518673 * (-(dy)));
              double h_x1_y2 = ((cse_p11_t1_g518674 * dx * dy) + (cse_p11_t1_g518677 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g518679 = (cse_p11_t1_g518673 * (-(dz)));
              double h_x1_z2 = ((cse_p11_t1_g518678 * dx) + (cse_p11_t1_g518679 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g518680 = (cse_p11_t1_g518674 * (cse_p11_t1_g518675 + (-(cse_p11_t54_g518614))));
              double h_y1_y1 = (cse_p11_t1_g518680 + (cse_p11_t1_g518673 * dy * dy));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = ((cse_p11_t1_g518673 * dy * dz) + (cse_p11_t1_g518678 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = ((cse_p11_t1_g518674 * dx * dy) + (cse_p11_t1_g518676 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = ((cse_p11_t1_g518674 * (cse_p11_t54_g518614 + (-(cse_p11_t1_g518675)))) + (cse_p11_t1_g518677 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = ((cse_p11_t1_g518678 * dy) + (cse_p11_t1_g518679 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g518681 = (cse_p11_t1_g518674 * (cse_p11_t1_g518675 + (-(cse_p11_t55_g518615))));
              double h_z1_z1 = (cse_p11_t1_g518681 + (cse_p11_t1_g518673 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = ((cse_p11_t1_g518676 * dz) + (cse_p11_t1_g518678 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = ((cse_p11_t1_g518677 * dz) + (cse_p11_t1_g518678 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = ((cse_p11_t1_g518674 * (cse_p11_t55_g518615 + (-(cse_p11_t1_g518675)))) + (cse_p11_t1_g518679 * dz));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p11_t1_g518682 + (cse_p11_t1_g518676 * (-(dx))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = ((cse_p11_t1_g518674 * dy * (-(dx))) + (cse_p11_t1_g518676 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = ((cse_p11_t1_g518676 * (-(dz))) + (cse_p11_t1_g518678 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = (cse_p11_t1_g518680 + (cse_p11_t1_g518677 * (-(dy))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = ((cse_p11_t1_g518677 * (-(dz))) + (cse_p11_t1_g518678 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = (cse_p11_t1_g518681 + (cse_p11_t1_g518679 * (-(dz))));
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
                double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t53_g518613))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t2_invr558251 = (1.0 / (rr));
                double cse_p565_t1_invr558250 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g518675 = (rr * rr);
                double cse_p11_t1_g518674 = (cse_p565_t1_invr558250 * de_dr);
                double cse_p11_t1_g518673 = (cse_p565_t2_invr558251 * cse_p565_t2_invr558251 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g518673 * dx * dy) + (cse_p11_t1_g518674 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double h_x1_z1 = ((cse_p11_t1_g518673 * dx * dz) + (cse_p11_t1_g518674 * dz * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g518678 = (cse_p11_t1_g518674 * dz);
                double h_x1_x2 = ((cse_p11_t1_g518673 * dx * (-(dx))) + (cse_p11_t1_g518674 * (cse_p11_t53_g518613 + (-(cse_p11_t1_g518675)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g518676 = (cse_p11_t1_g518673 * (-(dx)));
                double h_x1_y2 = ((cse_p11_t1_g518673 * dx * (-(dy))) + (cse_p11_t1_g518674 * dx * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g518677 = (cse_p11_t1_g518673 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g518673 * dx * (-(dz))) + (cse_p11_t1_g518678 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g518679 = (cse_p11_t1_g518673 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g518673 * dy * dy) + (cse_p11_t1_g518674 * (cse_p11_t1_g518675 + (-(cse_p11_t54_g518614)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g518680 = (cse_p11_t1_g518674 * (cse_p11_t1_g518675 + (-(cse_p11_t54_g518614))));
                double h_y1_z1 = ((cse_p11_t1_g518673 * dy * dz) + (cse_p11_t1_g518678 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g518674 * dx * dy) + (cse_p11_t1_g518676 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g518674 * (cse_p11_t54_g518614 + (-(cse_p11_t1_g518675)))) + (cse_p11_t1_g518677 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g518678 * dy) + (cse_p11_t1_g518679 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g518673 * dz * dz) + (cse_p11_t1_g518674 * (cse_p11_t1_g518675 + (-(cse_p11_t55_g518615)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g518681 = (cse_p11_t1_g518674 * (cse_p11_t1_g518675 + (-(cse_p11_t55_g518615))));
                double h_z1_x2 = ((cse_p11_t1_g518676 * dz) + (cse_p11_t1_g518678 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g518677 * dz) + (cse_p11_t1_g518678 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g518674 * (cse_p11_t55_g518615 + (-(cse_p11_t1_g518675)))) + (cse_p11_t1_g518679 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g518674 * (cse_p11_t1_g518675 + (-(cse_p11_t53_g518613)))) + (cse_p11_t1_g518676 * (-(dx))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p11_t1_g518682 = (cse_p11_t1_g518674 * (cse_p11_t1_g518675 + (-(cse_p11_t53_g518613))));
                double h_x2_y2 = ((cse_p11_t1_g518674 * dy * (-(dx))) + (cse_p11_t1_g518676 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double h_x2_z2 = ((cse_p11_t1_g518676 * (-(dz))) + (cse_p11_t1_g518678 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g518680 + (cse_p11_t1_g518677 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g518677 * (-(dz))) + (cse_p11_t1_g518678 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g518681 + (cse_p11_t1_g518679 * (-(dz))));
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
