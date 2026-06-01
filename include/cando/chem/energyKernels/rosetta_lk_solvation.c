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
      double cse_p1_t1_sqrt1 = sqrt(r2);
      double cse_p1_t2_invsqrt2 = (1.0 / (cse_p1_t1_sqrt1));
      double rr = cse_p1_t1_sqrt1;
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
          double cse_p528_t4_g9419 = (2.00000000000000000e+0 * rr);
          double cse_p528_t1_g9434 = (cse_p1_t2_invsqrt2 * dx);
          double cse_p528_t6_g9421 = (cse_p1_t2_invsqrt2 * dy);
          double cse_p528_t7_g9422 = (cse_p1_t2_invsqrt2 * dz);
          double cse_p528_t8_g9423 = (cse_p1_t2_invsqrt2 * rr);
          double cse_p528_t9_g9424 = (rr * rr);
          double cse_p14_t2_g119 = (3.00000000000000000e+0 * term.aa_low_j);
          double cse_p12_t1_g94 = (3.00000000000000000e+0 * term.aa_low_i);
          double cse_p13_t11_g105 = (cse_p12_t1_g94 + cse_p14_t2_g119);
          double cse_p12_t10_g63 = (term.bb_low_i * rr);
          double cse_p12_t11_g64 = (term.bb_low_j * rr);
          double cse_p12_t27_g80 = (cse_p13_t11_g105 * cse_p528_t9_g9424);
          double cse_p13_t1_g109 = (term.bb_low_i * cse_p528_t8_g9423);
          double cse_p12_t7_g60 = (cse_p13_t1_g109 * dx);
          double cse_p12_t8_g61 = (cse_p13_t1_g109 * dy);
          double cse_p12_t9_g62 = (cse_p13_t1_g109 * dz);
          double cse_p12_t12_g65 = (cse_p528_t8_g9423 * rr);
          double cse_p11_t27_g33 = (term.cc_low_i + term.cc_low_j + cse_p12_t10_g63 + cse_p12_t10_g63 + cse_p12_t11_g64 + cse_p12_t11_g64 + cse_p12_t27_g80);
          double cse_p11_t31_g37 = (-(cse_p12_t7_g60));
          double cse_p11_t32_g38 = (-(cse_p12_t8_g61));
          double cse_p11_t33_g39 = (-(cse_p12_t9_g62));
          double energy = (term.dd_low_i + term.dd_low_j + (term.cc_low_i * rr) + (term.cc_low_j * rr) + (cse_p528_t9_g9424 * (term.bb_low_i + term.bb_low_j + (term.aa_low_i * rr) + (term.aa_low_j * rr))));
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (cse_p11_t27_g33 * cse_p528_t1_g9434);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t27_g33 * cse_p528_t6_g9421);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t27_g33 * cse_p528_t7_g9422);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double cse_p12_t1_g86 = (cse_p12_t12_g65 * dx);
          double cse_p528_t1_g9416 = (-3.00000000000000000e+0 * cse_p12_t1_g86);
          double cse_p13_t1_g110 = (term.bb_low_j * cse_p528_t8_g9423);
          double cse_p63_t1_g886 = (-2.00000000000000000e+0 * cse_p13_t1_g110);
          double g_x2 = (cse_p11_t31_g37 + cse_p11_t31_g37 + (term.aa_low_i * cse_p528_t1_g9416) + (term.aa_low_j * cse_p528_t1_g9416) + (cse_p63_t1_g886 * dx) + (-((term.cc_low_i * cse_p528_t1_g9434))) + (-((term.cc_low_j * cse_p528_t1_g9434))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double cse_p12_t1_g87 = (cse_p12_t12_g65 * dy);
          double cse_p528_t2_g9417 = (-3.00000000000000000e+0 * cse_p12_t1_g87);
          double g_y2 = (cse_p11_t32_g38 + cse_p11_t32_g38 + (term.aa_low_i * cse_p528_t2_g9417) + (term.aa_low_j * cse_p528_t2_g9417) + (cse_p63_t1_g886 * dy) + (-((term.cc_low_i * cse_p528_t6_g9421))) + (-((term.cc_low_j * cse_p528_t6_g9421))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double cse_p12_t1_g88 = (cse_p12_t12_g65 * dz);
          double cse_p528_t3_g9418 = (-3.00000000000000000e+0 * cse_p12_t1_g88);
          double g_z2 = (cse_p11_t33_g39 + cse_p11_t33_g39 + (term.aa_low_i * cse_p528_t3_g9418) + (term.aa_low_j * cse_p528_t3_g9418) + (cse_p63_t1_g886 * dz) + (-((term.cc_low_i * cse_p528_t7_g9422))) + (-((term.cc_low_j * cse_p528_t7_g9422))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (term.cc_low_i + term.cc_low_j + cse_p12_t27_g80 + (term.bb_low_i * cse_p528_t4_g9419) + (term.bb_low_j * cse_p528_t4_g9419));
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= term.r2_high)) {
            /* !COND-3 */
            double cse_p1_t3_invr3 = (1.0 / (rr));
            double cse_p1_t4_invr24 = (cse_p1_t3_invr3 * cse_p1_t3_invr3);
            double inv_rr = cse_p1_t3_invr3;
            double inv_rr2 = (inv_rr * inv_rr);
            double di = (rr + (-(term.lj_radius_i)));
            double cse_p14_t7_g124 = (di * di);
            double cse_p13_t12_g106 = (-(cse_p14_t7_g124));
            double cse_p1_t5_invr5 = (1.0 / (term.lambda2_i));
            double cse_p12_t20_g73 = (cse_p13_t12_g106 * cse_p1_t5_invr5);
            double cse_p11_t36_g42 = exp(cse_p12_t20_g73);
            double exp_i = cse_p11_t36_g42;
            double cse_p12_t24_g77 = (exp_i * term.pref_i_j);
            double fi = (exp_i * inv_rr2 * term.pref_i_j);
            double hi = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t5_invr5 * di)));
            double dfi = (fi * hi);
            double dj = (rr + (-(term.lj_radius_j)));
            double cse_p14_t8_g125 = (dj * dj);
            double cse_p13_t13_g107 = (-(cse_p14_t8_g125));
            double cse_p1_t6_invr6 = (1.0 / (term.lambda2_j));
            double cse_p12_t21_g74 = (cse_p13_t13_g107 * cse_p1_t6_invr6);
            double cse_p11_t37_g43 = exp(cse_p12_t21_g74);
            double exp_j = cse_p11_t37_g43;
            double cse_p12_t25_g78 = (exp_j * term.pref_j_i);
            double cse_p11_t25_g31 = (cse_p12_t24_g77 + cse_p12_t25_g78);
            double fj = (exp_j * inv_rr2 * term.pref_j_i);
            double hj = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t6_invr6 * dj)));
            double dfj = (fj * hj);
            double energy = (fi + fj);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p11_t1_g45 = (-2.00000000000000000e+0 * inv_rr2);
            double cse_p528_t1_g9434 = (cse_p1_t2_invsqrt2 * dx);
            double cse_p11_t1_g50 = (cse_p11_t1_g45 * cse_p11_t36_g42 * cse_p1_t5_invr5 * term.pref_i_j);
            double cse_p11_t1_g51 = (cse_p11_t1_g45 * cse_p11_t37_g43 * cse_p1_t6_invr6 * term.pref_j_i);
            double cse_p11_t1_g52 = (cse_p11_t25_g31 * cse_p1_t4_invr24 * inv_rr);
            double cse_p12_t1_g90 = (-2.00000000000000000e+0 * cse_p11_t1_g52);
            double cse_p63_t1_g883 = (cse_p11_t1_g50 * di);
            double cse_p64_t13_g902 = (cse_p528_t1_g9434 * cse_p63_t1_g883);
            double cse_p63_t1_g884 = (cse_p11_t1_g51 * dj);
            double cse_p64_t14_g903 = (cse_p528_t1_g9434 * cse_p63_t1_g884);
            double g_x1 = (cse_p64_t13_g902 + cse_p64_t14_g903 + (cse_p12_t1_g90 * cse_p528_t1_g9434));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p11_t1_g48 = (cse_p1_t2_invsqrt2 * dy);
            double cse_p64_t5_g894 = (cse_p11_t1_g48 * cse_p63_t1_g883);
            double cse_p64_t6_g895 = (cse_p11_t1_g48 * cse_p63_t1_g884);
            double g_y1 = (cse_p64_t5_g894 + cse_p64_t6_g895 + (cse_p11_t1_g48 * cse_p12_t1_g90));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p11_t1_g49 = (cse_p1_t2_invsqrt2 * dz);
            double cse_p64_t7_g896 = (cse_p11_t1_g49 * cse_p63_t1_g883);
            double cse_p64_t8_g897 = (cse_p11_t1_g49 * cse_p63_t1_g884);
            double g_z1 = (cse_p64_t7_g896 + cse_p64_t8_g897 + (cse_p11_t1_g49 * cse_p12_t1_g90));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double cse_p12_t1_g91 = (2.00000000000000000e+0 * cse_p11_t1_g52);
            double g_x2 = ((cse_p12_t1_g91 * cse_p528_t1_g9434) + (-(cse_p64_t13_g902)) + (-(cse_p64_t14_g903)));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = ((cse_p11_t1_g48 * cse_p12_t1_g91) + (-(cse_p64_t5_g894)) + (-(cse_p64_t6_g895)));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = ((cse_p11_t1_g49 * cse_p12_t1_g91) + (-(cse_p64_t7_g896)) + (-(cse_p64_t8_g897)));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (dfi + dfj);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= term.r3_high)) {
              /* !COND-4 */
              double cse_p12_t1_g92 = (3.00000000000000000e+0 * term.aa_high_i);
              double cse_p12_t1_g93 = (3.00000000000000000e+0 * term.aa_high_j);
              double cse_p13_t10_g104 = (cse_p12_t1_g92 + cse_p12_t1_g93);
              double cse_p12_t5_g58 = (term.bb_high_i * rr);
              double cse_p12_t6_g59 = (term.bb_high_j * rr);
              double cse_p12_t1_g89 = (rr * rr);
              double cse_p12_t26_g79 = (cse_p12_t1_g89 * cse_p13_t10_g104);
              double cse_p11_t1_g44 = (cse_p1_t2_invsqrt2 * rr);
              double cse_p11_t1_g53 = (term.bb_high_i * cse_p11_t1_g44);
              double cse_p12_t13_g66 = (cse_p11_t1_g53 * dx);
              double cse_p12_t14_g67 = (cse_p11_t1_g53 * dy);
              double cse_p12_t15_g68 = (cse_p11_t1_g53 * dz);
              double cse_p11_t26_g32 = (term.cc_high_i + term.cc_high_j + cse_p12_t26_g79 + cse_p12_t5_g58 + cse_p12_t5_g58 + cse_p12_t6_g59 + cse_p12_t6_g59);
              double cse_p11_t28_g34 = (-(cse_p12_t13_g66));
              double cse_p11_t29_g35 = (-(cse_p12_t14_g67));
              double cse_p11_t30_g36 = (-(cse_p12_t15_g68));
              double energy = (term.dd_high_i + term.dd_high_j + (term.cc_high_i * rr) + (term.cc_high_j * rr) + (cse_p12_t1_g89 * (term.bb_high_i + term.bb_high_j + (term.aa_high_i * rr) + (term.aa_high_j * rr))));
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p528_t1_g9434 = (cse_p1_t2_invsqrt2 * dx);
              double g_x1 = (cse_p11_t26_g32 * cse_p528_t1_g9434);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double cse_p11_t1_g48 = (cse_p1_t2_invsqrt2 * dy);
              double g_y1 = (cse_p11_t1_g48 * cse_p11_t26_g32);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double cse_p11_t1_g49 = (cse_p1_t2_invsqrt2 * dz);
              double g_z1 = (cse_p11_t1_g49 * cse_p11_t26_g32);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p11_t1_g46 = (cse_p11_t1_g44 * rr);
              double cse_p12_t1_g86 = (cse_p11_t1_g46 * dx);
              double cse_p13_t1_g108 = (term.bb_high_j * cse_p11_t1_g44);
              double cse_p63_t1_g880 = (-3.00000000000000000e+0 * cse_p12_t1_g86);
              double cse_p63_t1_g885 = (-2.00000000000000000e+0 * cse_p13_t1_g108);
              double g_x2 = (cse_p11_t28_g34 + cse_p11_t28_g34 + (term.aa_high_i * cse_p63_t1_g880) + (term.aa_high_j * cse_p63_t1_g880) + (cse_p63_t1_g885 * dx) + (-((term.cc_high_i * cse_p528_t1_g9434))) + (-((term.cc_high_j * cse_p528_t1_g9434))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double cse_p12_t1_g87 = (cse_p11_t1_g46 * dy);
              double cse_p63_t1_g881 = (-3.00000000000000000e+0 * cse_p12_t1_g87);
              double g_y2 = (cse_p11_t29_g35 + cse_p11_t29_g35 + (term.aa_high_i * cse_p63_t1_g881) + (term.aa_high_j * cse_p63_t1_g881) + (cse_p63_t1_g885 * dy) + (-((term.cc_high_i * cse_p11_t1_g48))) + (-((term.cc_high_j * cse_p11_t1_g48))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double cse_p12_t1_g88 = (cse_p11_t1_g46 * dz);
              double cse_p63_t1_g882 = (-3.00000000000000000e+0 * cse_p12_t1_g88);
              double g_z2 = (cse_p11_t30_g36 + cse_p11_t30_g36 + (term.aa_high_i * cse_p63_t1_g882) + (term.aa_high_j * cse_p63_t1_g882) + (cse_p63_t1_g885 * dz) + (-((term.cc_high_i * cse_p11_t1_g49))) + (-((term.cc_high_j * cse_p11_t1_g49))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double cse_p12_t1_g85 = (2.00000000000000000e+0 * rr);
              double de_dr = (term.cc_high_i + term.cc_high_j + cse_p12_t26_g79 + (term.bb_high_i * cse_p12_t1_g85) + (term.bb_high_j * cse_p12_t1_g85));
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
            double cse_p61_t1_invr552 = (1.0 / (term.lambda2_i));
            double exp_i = exp((-((cse_p61_t1_invr552 * di * di))));
            double fi = (exp_i * inv_rr2 * term.pref_i_j);
            double hi = (-2.00000000000000000e+0 * (inv_rr + (cse_p61_t1_invr552 * di)));
            double dj = (rr + (-(term.lj_radius_j)));
            double cse_p61_t2_invr553 = (1.0 / (term.lambda2_j));
            double exp_j = exp((-((cse_p61_t2_invr553 * dj * dj))));
            double fj = (exp_j * inv_rr2 * term.pref_j_i);
            double hj = (-2.00000000000000000e+0 * (inv_rr + (cse_p61_t2_invr553 * dj)));
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
      double cse_p11_t53_g59 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t54_g60 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t55_g61 = (dz * dz);
      double r2 = (cse_p11_t53_g59 + cse_p11_t54_g60 + cse_p11_t55_g61);
      double cse_p1_t1_sqrt1 = sqrt(r2);
      double cse_p1_t2_invsqrt2 = (1.0 / (cse_p1_t1_sqrt1));
      double rr = cse_p1_t1_sqrt1;
      if ((rr <= term.r0_low)) {
        /* !COND-1 */
        double cse_p60_t43_g4108 = (rr * rr);
        double cse_p60_t71_g4136 = (-(cse_p11_t53_g59));
        double cse_p60_t72_g4137 = (-(cse_p11_t54_g60));
        double cse_p60_t73_g4138 = (-(cse_p11_t55_g61));
        double cse_p60_t75_g4140 = (-(dx));
        double cse_p60_t76_g4141 = (-(dy));
        double cse_p60_t77_g4142 = (-(dz));
        double cse_p565_t1_invr39696 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr39697 = (1.0 / (rr));
        double cse_p60_t80_g4145 = (rr * rr * rr);
        double cse_p15_t75_g539 = (-(cse_p60_t43_g4108));
        double cse_p14_t82_g437 = (cse_p11_t53_g59 + cse_p15_t75_g539);
        double cse_p14_t84_g439 = (cse_p11_t54_g60 + cse_p15_t75_g539);
        double cse_p14_t86_g441 = (cse_p11_t55_g61 + cse_p15_t75_g539);
        double cse_p13_t109_g340 = (cse_p60_t43_g4108 + cse_p60_t71_g4136);
        double cse_p13_t110_g341 = (cse_p60_t43_g4108 + cse_p60_t72_g4137);
        double cse_p13_t111_g342 = (cse_p60_t43_g4108 + cse_p60_t73_g4138);
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
        double cse_p60_t42_g4107 = (cse_p565_t1_invr39696 * de_dr);
        double cse_p13_t67_g298 = (cse_p14_t82_g437 * cse_p60_t42_g4107);
        double cse_p13_t68_g299 = (cse_p14_t84_g439 * cse_p60_t42_g4107);
        double cse_p13_t69_g300 = (cse_p14_t86_g441 * cse_p60_t42_g4107);
        double cse_p13_t70_g301 = (cse_p60_t42_g4107 * dx * dy);
        double cse_p13_t71_g302 = (cse_p60_t42_g4107 * cse_p60_t75_g4140 * dy);
        double cse_p12_t27_g155 = (cse_p13_t109_g340 * cse_p60_t42_g4107);
        double cse_p12_t28_g156 = (cse_p13_t110_g341 * cse_p60_t42_g4107);
        double cse_p12_t29_g157 = (cse_p13_t111_g342 * cse_p60_t42_g4107);
        double cse_p12_t32_g160 = (cse_p60_t42_g4107 * dz);
        double cse_p13_t57_g288 = (cse_p12_t32_g160 * cse_p60_t75_g4140);
        double cse_p13_t58_g289 = (cse_p12_t32_g160 * cse_p60_t76_g4141);
        double cse_p13_t59_g290 = (cse_p12_t32_g160 * dx);
        double cse_p13_t60_g291 = (cse_p12_t32_g160 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p60_t41_g4106 = (cse_p565_t2_invr39697 * cse_p565_t2_invr39697 * d2e_dr2);
        double cse_p13_t61_g292 = (cse_p60_t41_g4106 * dx * dx);
        double cse_p13_t62_g293 = (cse_p60_t41_g4106 * dx * dy);
        double cse_p13_t63_g294 = (cse_p60_t41_g4106 * dx * dz);
        double cse_p13_t64_g295 = (cse_p60_t41_g4106 * dy * dy);
        double cse_p13_t65_g296 = (cse_p60_t41_g4106 * dy * dz);
        double cse_p13_t66_g297 = (cse_p60_t41_g4106 * dz * dz);
        double cse_p12_t15_g143 = (cse_p60_t41_g4106 * cse_p60_t75_g4140);
        double cse_p13_t42_g273 = (cse_p12_t15_g143 * cse_p60_t75_g4140);
        double cse_p13_t43_g274 = (cse_p12_t15_g143 * cse_p60_t76_g4141);
        double cse_p13_t44_g275 = (cse_p12_t15_g143 * cse_p60_t77_g4142);
        double cse_p13_t45_g276 = (cse_p12_t15_g143 * dx);
        double cse_p13_t46_g277 = (cse_p12_t15_g143 * dy);
        double cse_p13_t47_g278 = (cse_p12_t15_g143 * dz);
        double cse_p12_t16_g144 = (cse_p60_t41_g4106 * cse_p60_t76_g4141);
        double cse_p13_t48_g279 = (cse_p12_t16_g144 * cse_p60_t76_g4141);
        double cse_p13_t49_g280 = (cse_p12_t16_g144 * cse_p60_t77_g4142);
        double cse_p13_t50_g281 = (cse_p12_t16_g144 * dx);
        double cse_p13_t51_g282 = (cse_p12_t16_g144 * dy);
        double cse_p13_t52_g283 = (cse_p12_t16_g144 * dz);
        double cse_p12_t17_g145 = (cse_p60_t41_g4106 * cse_p60_t77_g4142);
        double cse_p13_t53_g284 = (cse_p12_t17_g145 * cse_p60_t77_g4142);
        double cse_p13_t54_g285 = (cse_p12_t17_g145 * dx);
        double cse_p13_t55_g286 = (cse_p12_t17_g145 * dy);
        double cse_p13_t56_g287 = (cse_p12_t17_g145 * dz);
        double cse_p12_t68_g196 = (cse_p13_t62_g293 + cse_p13_t71_g302);
        double cse_p12_t69_g197 = (cse_p13_t57_g288 + cse_p13_t63_g294);
        double cse_p12_t70_g198 = (cse_p13_t58_g289 + cse_p13_t65_g296);
        double cse_p12_t86_g214 = (cse_p12_t28_g156 + cse_p13_t64_g295);
        double cse_p12_t88_g216 = (cse_p12_t29_g157 + cse_p13_t66_g297);
        double cse_p12_t90_g218 = (cse_p12_t27_g155 + cse_p13_t61_g292);
        double cse_p12_t63_g191 = (cse_p13_t45_g276 + cse_p13_t67_g298);
        double cse_p12_t71_g199 = (cse_p13_t46_g277 + cse_p13_t70_g301);
        double cse_p12_t73_g201 = (cse_p13_t43_g274 + cse_p13_t71_g302);
        double cse_p12_t74_g202 = (cse_p13_t44_g275 + cse_p13_t57_g288);
        double cse_p12_t75_g203 = (cse_p13_t47_g278 + cse_p13_t59_g290);
        double cse_p12_t91_g219 = (cse_p12_t27_g155 + cse_p13_t42_g273);
        double cse_p12_t64_g192 = (cse_p13_t51_g282 + cse_p13_t68_g299);
        double cse_p12_t72_g200 = (cse_p13_t50_g281 + cse_p13_t70_g301);
        double cse_p12_t76_g204 = (cse_p13_t49_g280 + cse_p13_t58_g289);
        double cse_p12_t77_g205 = (cse_p13_t52_g283 + cse_p13_t60_g291);
        double cse_p12_t87_g215 = (cse_p12_t28_g156 + cse_p13_t48_g279);
        double cse_p12_t65_g193 = (cse_p13_t56_g287 + cse_p13_t69_g300);
        double cse_p12_t78_g206 = (cse_p13_t54_g285 + cse_p13_t59_g290);
        double cse_p12_t79_g207 = (cse_p13_t55_g286 + cse_p13_t60_g291);
        double cse_p12_t89_g217 = (cse_p12_t29_g157 + cse_p13_t53_g284);
        double h_x1_x1 = cse_p12_t90_g218;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t68_g196;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t69_g197;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t63_g191;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t72_g200;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t78_g206;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t86_g214;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t70_g198;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t71_g199;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t64_g192;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t79_g207;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t88_g216;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t75_g203;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t77_g205;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t65_g193;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t91_g219;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t73_g201;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t74_g202;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t87_g215;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t76_g204;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t89_g217;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.r1_low)) {
          /* !COND-2 */
          double cse_p60_t38_g4103 = (cse_p1_t2_invsqrt2 * dx);
          double cse_p60_t39_g4104 = (cse_p1_t2_invsqrt2 * dy);
          double cse_p60_t40_g4105 = (cse_p1_t2_invsqrt2 * dz);
          double cse_p14_t3_g358 = (3.00000000000000000e+0 * term.aa_low_i);
          double cse_p14_t4_g359 = (3.00000000000000000e+0 * term.aa_low_j);
          double cse_p13_t81_g312 = (cse_p14_t3_g358 + cse_p14_t4_g359);
          double cse_p12_t10_g138 = (term.bb_low_i * cse_p60_t38_g4103 * rr);
          double cse_p12_t11_g139 = (term.bb_low_i * cse_p60_t39_g4104 * rr);
          double cse_p12_t12_g140 = (term.bb_low_i * cse_p60_t40_g4105 * rr);
          double cse_p12_t13_g141 = (term.bb_low_i * rr);
          double cse_p12_t14_g142 = (term.bb_low_j * rr);
          double cse_p11_t1_g121 = (rr * rr);
          double cse_p60_t59_g4124 = (cse_p11_t53_g59 + (-(cse_p11_t1_g121)));
          double cse_p60_t60_g4125 = (cse_p11_t54_g60 + (-(cse_p11_t1_g121)));
          double cse_p60_t61_g4126 = (cse_p11_t55_g61 + (-(cse_p11_t1_g121)));
          double cse_p60_t62_g4127 = (cse_p11_t1_g121 + (-(cse_p11_t53_g59)));
          double cse_p60_t63_g4128 = (cse_p11_t1_g121 + (-(cse_p11_t54_g60)));
          double cse_p60_t64_g4129 = (cse_p11_t1_g121 + (-(cse_p11_t55_g61)));
          double cse_p60_t74_g4139 = (-(cse_p11_t1_g121));
          double cse_p12_t34_g162 = (cse_p11_t1_g121 * cse_p13_t81_g312);
          double cse_p11_t92_g98 = (term.cc_low_i + term.cc_low_j + cse_p12_t13_g141 + cse_p12_t13_g141 + cse_p12_t14_g142 + cse_p12_t14_g142 + cse_p12_t34_g162);
          double cse_p11_t96_g102 = (-(cse_p12_t10_g138));
          double cse_p11_t97_g103 = (-(cse_p12_t11_g139));
          double cse_p11_t98_g104 = (-(cse_p12_t12_g140));
          double energy = (term.dd_low_i + term.dd_low_j + (term.cc_low_i * rr) + (term.cc_low_j * rr) + (cse_p11_t1_g121 * (term.bb_low_i + term.bb_low_j + (term.aa_low_i * rr) + (term.aa_low_j * rr))));
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (cse_p11_t92_g98 * cse_p60_t38_g4103);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t92_g98 * cse_p60_t39_g4104);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t92_g98 * cse_p60_t40_g4105);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p11_t96_g102 + cse_p11_t96_g102 + (-2.00000000000000000e+0 * term.bb_low_j * cse_p60_t38_g4103 * rr) + (cse_p11_t1_g121 * cse_p60_t38_g4103 * ((-3.00000000000000000e+0 * term.aa_low_i) + (-3.00000000000000000e+0 * term.aa_low_j))) + (-((term.cc_low_i * cse_p60_t38_g4103))) + (-((term.cc_low_j * cse_p60_t38_g4103))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p11_t97_g103 + cse_p11_t97_g103 + (-2.00000000000000000e+0 * term.bb_low_j * cse_p60_t39_g4104 * rr) + (cse_p11_t1_g121 * cse_p60_t39_g4104 * ((-3.00000000000000000e+0 * term.aa_low_i) + (-3.00000000000000000e+0 * term.aa_low_j))) + (-((term.cc_low_i * cse_p60_t39_g4104))) + (-((term.cc_low_j * cse_p60_t39_g4104))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p11_t98_g104 + cse_p11_t98_g104 + (-2.00000000000000000e+0 * term.bb_low_j * cse_p60_t40_g4105 * rr) + (cse_p11_t1_g121 * cse_p60_t40_g4105 * ((-3.00000000000000000e+0 * term.aa_low_i) + (-3.00000000000000000e+0 * term.aa_low_j))) + (-((term.cc_low_i * cse_p60_t40_g4105))) + (-((term.cc_low_j * cse_p60_t40_g4105))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (term.cc_low_i + term.cc_low_j + cse_p12_t34_g162 + (2.00000000000000000e+0 * term.bb_low_i * rr) + (2.00000000000000000e+0 * term.bb_low_j * rr));
          double d2e_dr2 = ((2.00000000000000000e+0 * term.bb_low_i) + (2.00000000000000000e+0 * term.bb_low_j) + (6.00000000000000000e+0 * term.aa_low_i * rr) + (6.00000000000000000e+0 * term.aa_low_j * rr));
          double cse_p565_t2_invr39697 = (1.0 / (rr));
          double cse_p11_t1_g119 = (cse_p565_t2_invr39697 * cse_p565_t2_invr39697 * d2e_dr2);
          double cse_p60_t1_g4066 = (cse_p11_t1_g119 * (-(dx)));
          double cse_p60_t2_g4067 = (cse_p11_t1_g119 * (-(dy)));
          double cse_p60_t3_g4068 = (cse_p11_t1_g119 * (-(dz)));
          double cse_p60_t4_g4069 = (cse_p11_t1_g119 * dx * dx);
          double cse_p60_t5_g4070 = (cse_p11_t1_g119 * dx * dy);
          double cse_p60_t6_g4071 = (cse_p11_t1_g119 * dx * dz);
          double cse_p60_t7_g4072 = (cse_p11_t1_g119 * dy * dy);
          double cse_p60_t8_g4073 = (cse_p11_t1_g119 * dy * dz);
          double cse_p60_t9_g4074 = (cse_p11_t1_g119 * dz * dz);
          double cse_p14_t46_g401 = (cse_p60_t1_g4066 * (-(dx)));
          double cse_p14_t47_g402 = (cse_p60_t1_g4066 * (-(dy)));
          double cse_p14_t48_g403 = (cse_p60_t1_g4066 * (-(dz)));
          double cse_p14_t49_g404 = (cse_p60_t1_g4066 * dx);
          double cse_p14_t50_g405 = (cse_p60_t1_g4066 * dy);
          double cse_p14_t51_g406 = (cse_p60_t1_g4066 * dz);
          double cse_p14_t52_g407 = (cse_p60_t2_g4067 * (-(dy)));
          double cse_p14_t53_g408 = (cse_p60_t2_g4067 * (-(dz)));
          double cse_p14_t54_g409 = (cse_p60_t2_g4067 * dx);
          double cse_p14_t55_g410 = (cse_p60_t2_g4067 * dy);
          double cse_p14_t56_g411 = (cse_p60_t2_g4067 * dz);
          double cse_p14_t57_g412 = (cse_p60_t3_g4068 * (-(dz)));
          double cse_p14_t58_g413 = (cse_p60_t3_g4068 * dx);
          double cse_p14_t59_g414 = (cse_p60_t3_g4068 * dy);
          double cse_p14_t60_g415 = (cse_p60_t3_g4068 * dz);
          double cse_p565_t1_invr39696 = (1.0 / ((rr * rr * rr)));
          double cse_p11_t1_g120 = (cse_p565_t1_invr39696 * de_dr);
          double cse_p60_t10_g4075 = (cse_p11_t1_g120 * (cse_p11_t53_g59 + cse_p60_t74_g4139));
          double cse_p60_t11_g4076 = (cse_p11_t1_g120 * (cse_p11_t54_g60 + cse_p60_t74_g4139));
          double cse_p60_t12_g4077 = (cse_p11_t1_g120 * (cse_p11_t55_g61 + cse_p60_t74_g4139));
          double cse_p60_t13_g4078 = (cse_p11_t1_g120 * cse_p60_t62_g4127);
          double cse_p60_t14_g4079 = (cse_p11_t1_g120 * cse_p60_t63_g4128);
          double cse_p60_t15_g4080 = (cse_p11_t1_g120 * cse_p60_t64_g4129);
          double cse_p60_t16_g4081 = (cse_p11_t1_g120 * dx * dy);
          double cse_p60_t17_g4082 = (cse_p11_t1_g120 * dy * (-(dx)));
          double cse_p60_t18_g4083 = (cse_p11_t1_g120 * dz);
          double cse_p60_t47_g4112 = (cse_p60_t17_g4082 + cse_p60_t5_g4070);
          double cse_p14_t42_g397 = (cse_p60_t18_g4083 * (-(dx)));
          double cse_p14_t43_g398 = (cse_p60_t18_g4083 * (-(dy)));
          double cse_p14_t44_g399 = (cse_p60_t18_g4083 * dx);
          double cse_p14_t45_g400 = (cse_p60_t18_g4083 * dy);
          double cse_p13_t107_g338 = (cse_p60_t13_g4078 + cse_p60_t4_g4069);
          double h_x1_x1 = cse_p13_t107_g338;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t47_g4112;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t83_g314 = (cse_p14_t42_g397 + cse_p60_t6_g4071);
          double cse_p13_t84_g315 = (cse_p14_t43_g398 + cse_p60_t8_g4073);
          double h_x1_z1 = cse_p13_t83_g314;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t77_g308 = (cse_p14_t49_g404 + cse_p60_t10_g4075);
          double cse_p13_t85_g316 = (cse_p14_t50_g405 + cse_p60_t16_g4081);
          double cse_p13_t87_g318 = (cse_p14_t47_g402 + cse_p60_t17_g4082);
          double cse_p13_t88_g319 = (cse_p14_t42_g397 + cse_p14_t48_g403);
          double cse_p13_t89_g320 = (cse_p14_t44_g399 + cse_p14_t51_g406);
          double cse_p13_t108_g339 = (cse_p14_t46_g401 + cse_p60_t13_g4078);
          double h_x1_x2 = cse_p13_t77_g308;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t78_g309 = (cse_p14_t55_g410 + cse_p60_t11_g4076);
          double cse_p13_t86_g317 = (cse_p14_t54_g409 + cse_p60_t16_g4081);
          double cse_p13_t90_g321 = (cse_p14_t43_g398 + cse_p14_t53_g408);
          double cse_p13_t91_g322 = (cse_p14_t45_g400 + cse_p14_t56_g411);
          double h_x1_y2 = cse_p13_t86_g317;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t79_g310 = (cse_p14_t60_g415 + cse_p60_t12_g4077);
          double cse_p13_t92_g323 = (cse_p14_t44_g399 + cse_p14_t58_g413);
          double cse_p13_t93_g324 = (cse_p14_t45_g400 + cse_p14_t59_g414);
          double h_x1_z2 = cse_p13_t92_g323;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t103_g334 = (cse_p60_t14_g4079 + cse_p60_t7_g4072);
          double cse_p13_t104_g335 = (cse_p14_t52_g407 + cse_p60_t14_g4079);
          double h_y1_y1 = cse_p13_t103_g334;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t84_g315;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t85_g316;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t78_g309;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t93_g324;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t105_g336 = (cse_p60_t15_g4080 + cse_p60_t9_g4074);
          double cse_p13_t106_g337 = (cse_p14_t57_g412 + cse_p60_t15_g4080);
          double h_z1_z1 = cse_p13_t105_g336;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t89_g320;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t91_g322;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t79_g310;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t108_g339;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t87_g318;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t88_g319;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t104_g335;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t90_g321;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t106_g337;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= term.r2_high)) {
            /* !COND-3 */
            double cse_p565_t2_invr39697 = (1.0 / (rr));
            double cse_p1_t4_invr24 = (cse_p565_t2_invr39697 * cse_p565_t2_invr39697);
            double inv_rr = cse_p565_t2_invr39697;
            double inv_rr2 = (inv_rr * inv_rr);
            double cse_p11_t1_g7 = (2.00000000000000000e+0 * inv_rr2);
            double di = (rr + (-(term.lj_radius_i)));
            double cse_p14_t63_g418 = (di * di);
            double cse_p13_t112_g343 = (-(cse_p14_t63_g418));
            double cse_p1_t5_invr5 = (1.0 / (term.lambda2_i));
            double cse_p12_t54_g182 = (cse_p13_t112_g343 * cse_p1_t5_invr5);
            double cse_p11_t108_g114 = exp(cse_p12_t54_g182);
            double exp_i = cse_p11_t108_g114;
            double cse_p12_t60_g188 = (exp_i * term.pref_i_j);
            double fi = (exp_i * inv_rr2 * term.pref_i_j);
            double hi = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t5_invr5 * di)));
            double dfi = (fi * hi);
            double d2fi = (fi * (cse_p11_t1_g7 + (-2.00000000000000000e+0 * cse_p1_t5_invr5) + (hi * hi)));
            double dj = (rr + (-(term.lj_radius_j)));
            double cse_p14_t64_g419 = (dj * dj);
            double cse_p13_t113_g344 = (-(cse_p14_t64_g419));
            double cse_p1_t6_invr6 = (1.0 / (term.lambda2_j));
            double cse_p12_t55_g183 = (cse_p13_t113_g344 * cse_p1_t6_invr6);
            double cse_p11_t109_g115 = exp(cse_p12_t55_g183);
            double exp_j = cse_p11_t109_g115;
            double cse_p12_t61_g189 = (exp_j * term.pref_j_i);
            double cse_p11_t87_g93 = (cse_p12_t60_g188 + cse_p12_t61_g189);
            double fj = (exp_j * inv_rr2 * term.pref_j_i);
            double hj = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t6_invr6 * dj)));
            double dfj = (fj * hj);
            double d2fj = (fj * (cse_p11_t1_g7 + (-2.00000000000000000e+0 * cse_p1_t6_invr6) + (hj * hj)));
            double energy = (fi + fj);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p14_t1_g462 = (cse_p1_t2_invsqrt2 * dx);
            double g_x1 = ((-2.00000000000000000e+0 * cse_p11_t87_g93 * cse_p14_t1_g462 * cse_p1_t4_invr24 * inv_rr) + (cse_p14_t1_g462 * inv_rr2 * ((-2.00000000000000000e+0 * cse_p11_t108_g114 * cse_p1_t5_invr5 * di * term.pref_i_j) + (-2.00000000000000000e+0 * cse_p11_t109_g115 * cse_p1_t6_invr6 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p14_t1_g463 = (cse_p1_t2_invsqrt2 * dy);
            double g_y1 = ((-2.00000000000000000e+0 * cse_p11_t87_g93 * cse_p14_t1_g463 * cse_p1_t4_invr24 * inv_rr) + (cse_p14_t1_g463 * inv_rr2 * ((-2.00000000000000000e+0 * cse_p11_t108_g114 * cse_p1_t5_invr5 * di * term.pref_i_j) + (-2.00000000000000000e+0 * cse_p11_t109_g115 * cse_p1_t6_invr6 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p14_t1_g464 = (cse_p1_t2_invsqrt2 * dz);
            double g_z1 = ((-2.00000000000000000e+0 * cse_p11_t87_g93 * cse_p14_t1_g464 * cse_p1_t4_invr24 * inv_rr) + (cse_p14_t1_g464 * inv_rr2 * ((-2.00000000000000000e+0 * cse_p11_t108_g114 * cse_p1_t5_invr5 * di * term.pref_i_j) + (-2.00000000000000000e+0 * cse_p11_t109_g115 * cse_p1_t6_invr6 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = ((2.00000000000000000e+0 * cse_p11_t87_g93 * cse_p14_t1_g462 * cse_p1_t4_invr24 * inv_rr) + (cse_p14_t1_g462 * inv_rr2 * ((2.00000000000000000e+0 * cse_p11_t108_g114 * cse_p1_t5_invr5 * di * term.pref_i_j) + (2.00000000000000000e+0 * cse_p11_t109_g115 * cse_p1_t6_invr6 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = ((2.00000000000000000e+0 * cse_p11_t87_g93 * cse_p14_t1_g463 * cse_p1_t4_invr24 * inv_rr) + (cse_p14_t1_g463 * inv_rr2 * ((2.00000000000000000e+0 * cse_p11_t108_g114 * cse_p1_t5_invr5 * di * term.pref_i_j) + (2.00000000000000000e+0 * cse_p11_t109_g115 * cse_p1_t6_invr6 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = ((2.00000000000000000e+0 * cse_p11_t87_g93 * cse_p14_t1_g464 * cse_p1_t4_invr24 * inv_rr) + (cse_p14_t1_g464 * inv_rr2 * ((2.00000000000000000e+0 * cse_p11_t108_g114 * cse_p1_t5_invr5 * di * term.pref_i_j) + (2.00000000000000000e+0 * cse_p11_t109_g115 * cse_p1_t6_invr6 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (dfi + dfj);
            double d2e_dr2 = (d2fi + d2fj);
            double cse_p11_t1_g119 = (cse_p565_t2_invr39697 * cse_p565_t2_invr39697 * d2e_dr2);
            double cse_p565_t1_invr39696 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g120 = (cse_p565_t1_invr39696 * de_dr);
            double cse_p11_t1_g121 = (rr * rr);
            double cse_p11_t1_g128 = (cse_p11_t1_g120 * (cse_p11_t1_g121 + (-(cse_p11_t53_g59))));
            double cse_p60_t69_g4134 = (cse_p11_t1_g128 + (cse_p11_t1_g119 * dx * dx));
            double h_x1_x1 = cse_p60_t69_g4134;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = ((cse_p11_t1_g119 * dx * dy) + (cse_p11_t1_g120 * dy * (-(dx))));
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p11_t1_g124 = (cse_p11_t1_g120 * dz);
            double cse_p60_t30_g4095 = (cse_p11_t1_g124 * (-(dx)));
            double cse_p60_t31_g4096 = (cse_p11_t1_g124 * (-(dy)));
            double cse_p60_t32_g4097 = (cse_p11_t1_g124 * dx);
            double cse_p60_t33_g4098 = (cse_p11_t1_g124 * dy);
            double cse_p60_t48_g4113 = (cse_p60_t30_g4095 + (cse_p11_t1_g119 * dx * dz));
            double cse_p60_t49_g4114 = (cse_p60_t31_g4096 + (cse_p11_t1_g119 * dy * dz));
            double h_x1_z1 = cse_p60_t48_g4113;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p11_t1_g122 = (cse_p11_t1_g119 * (-(dx)));
            double cse_p60_t19_g4084 = (cse_p11_t1_g122 * (-(dx)));
            double cse_p60_t20_g4085 = (cse_p11_t1_g122 * (-(dy)));
            double cse_p60_t21_g4086 = (cse_p11_t1_g122 * (-(dz)));
            double cse_p60_t22_g4087 = (cse_p11_t1_g122 * dx);
            double cse_p60_t23_g4088 = (cse_p11_t1_g122 * dy);
            double cse_p60_t24_g4089 = (cse_p11_t1_g122 * dz);
            double cse_p60_t44_g4109 = (cse_p60_t22_g4087 + (cse_p11_t1_g120 * (cse_p11_t53_g59 + (-(cse_p11_t1_g121)))));
            double cse_p60_t50_g4115 = (cse_p60_t23_g4088 + (cse_p11_t1_g120 * dx * dy));
            double cse_p60_t52_g4117 = (cse_p60_t20_g4085 + (cse_p11_t1_g120 * dy * (-(dx))));
            double cse_p60_t53_g4118 = (cse_p60_t21_g4086 + cse_p60_t30_g4095);
            double cse_p60_t54_g4119 = (cse_p60_t24_g4089 + cse_p60_t32_g4097);
            double cse_p60_t70_g4135 = (cse_p11_t1_g128 + cse_p60_t19_g4084);
            double h_x1_x2 = cse_p60_t44_g4109;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p11_t1_g123 = (cse_p11_t1_g119 * (-(dy)));
            double cse_p60_t25_g4090 = (cse_p11_t1_g123 * (-(dy)));
            double cse_p60_t26_g4091 = (cse_p11_t1_g123 * (-(dz)));
            double cse_p60_t27_g4092 = (cse_p11_t1_g123 * dx);
            double cse_p60_t28_g4093 = (cse_p11_t1_g123 * dy);
            double cse_p60_t29_g4094 = (cse_p11_t1_g123 * dz);
            double cse_p60_t45_g4110 = (cse_p60_t28_g4093 + (cse_p11_t1_g120 * (cse_p11_t54_g60 + (-(cse_p11_t1_g121)))));
            double cse_p60_t51_g4116 = (cse_p60_t27_g4092 + (cse_p11_t1_g120 * dx * dy));
            double cse_p60_t55_g4120 = (cse_p60_t26_g4091 + cse_p60_t31_g4096);
            double cse_p60_t56_g4121 = (cse_p60_t29_g4094 + cse_p60_t33_g4098);
            double h_x1_y2 = cse_p60_t51_g4116;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p11_t1_g125 = (cse_p11_t1_g119 * (-(dz)));
            double cse_p60_t34_g4099 = (cse_p11_t1_g125 * (-(dz)));
            double cse_p60_t35_g4100 = (cse_p11_t1_g125 * dx);
            double cse_p60_t36_g4101 = (cse_p11_t1_g125 * dy);
            double cse_p60_t37_g4102 = (cse_p11_t1_g125 * dz);
            double cse_p60_t46_g4111 = (cse_p60_t37_g4102 + (cse_p11_t1_g120 * (cse_p11_t55_g61 + (-(cse_p11_t1_g121)))));
            double cse_p60_t57_g4122 = (cse_p60_t32_g4097 + cse_p60_t35_g4100);
            double cse_p60_t58_g4123 = (cse_p60_t33_g4098 + cse_p60_t36_g4101);
            double h_x1_z2 = cse_p60_t57_g4122;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p11_t1_g126 = (cse_p11_t1_g120 * (cse_p11_t1_g121 + (-(cse_p11_t54_g60))));
            double cse_p60_t65_g4130 = (cse_p11_t1_g126 + (cse_p11_t1_g119 * dy * dy));
            double cse_p60_t66_g4131 = (cse_p11_t1_g126 + cse_p60_t25_g4090);
            double h_y1_y1 = cse_p60_t65_g4130;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t49_g4114;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t50_g4115;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t45_g4110;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t58_g4123;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p11_t1_g127 = (cse_p11_t1_g120 * (cse_p11_t1_g121 + (-(cse_p11_t55_g61))));
            double cse_p60_t67_g4132 = (cse_p11_t1_g127 + (cse_p11_t1_g119 * dz * dz));
            double cse_p60_t68_g4133 = (cse_p11_t1_g127 + cse_p60_t34_g4099);
            double h_z1_z1 = cse_p60_t67_g4132;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t54_g4119;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t56_g4121;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t46_g4111;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t70_g4135;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t52_g4117;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t53_g4118;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t66_g4131;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t55_g4120;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t68_g4133;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= term.r3_high)) {
              /* !COND-4 */
              double cse_p14_t1_g356 = (3.00000000000000000e+0 * term.aa_high_i);
              double cse_p14_t2_g357 = (3.00000000000000000e+0 * term.aa_high_j);
              double cse_p13_t80_g311 = (cse_p14_t1_g356 + cse_p14_t2_g357);
              double cse_p14_t1_g462 = (cse_p1_t2_invsqrt2 * dx);
              double cse_p12_t5_g133 = (term.bb_high_i * cse_p14_t1_g462 * rr);
              double cse_p14_t1_g463 = (cse_p1_t2_invsqrt2 * dy);
              double cse_p12_t6_g134 = (term.bb_high_i * cse_p14_t1_g463 * rr);
              double cse_p14_t1_g464 = (cse_p1_t2_invsqrt2 * dz);
              double cse_p12_t7_g135 = (term.bb_high_i * cse_p14_t1_g464 * rr);
              double cse_p12_t8_g136 = (term.bb_high_i * rr);
              double cse_p12_t9_g137 = (term.bb_high_j * rr);
              double cse_p11_t1_g121 = (rr * rr);
              double cse_p12_t33_g161 = (cse_p11_t1_g121 * cse_p13_t80_g311);
              double cse_p11_t91_g97 = (term.cc_high_i + term.cc_high_j + cse_p12_t33_g161 + cse_p12_t8_g136 + cse_p12_t8_g136 + cse_p12_t9_g137 + cse_p12_t9_g137);
              double cse_p11_t93_g99 = (-(cse_p12_t5_g133));
              double cse_p11_t94_g100 = (-(cse_p12_t6_g134));
              double cse_p11_t95_g101 = (-(cse_p12_t7_g135));
              double energy = (term.dd_high_i + term.dd_high_j + (term.cc_high_i * rr) + (term.cc_high_j * rr) + (cse_p11_t1_g121 * (term.bb_high_i + term.bb_high_j + (term.aa_high_i * rr) + (term.aa_high_j * rr))));
              energy_added += energy;
              *energy_accumulate += energy;
              double g_x1 = (cse_p11_t91_g97 * cse_p14_t1_g462);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t91_g97 * cse_p14_t1_g463);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t91_g97 * cse_p14_t1_g464);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = (cse_p11_t93_g99 + cse_p11_t93_g99 + (-2.00000000000000000e+0 * term.bb_high_j * cse_p14_t1_g462 * rr) + (cse_p11_t1_g121 * cse_p14_t1_g462 * ((-3.00000000000000000e+0 * term.aa_high_i) + (-3.00000000000000000e+0 * term.aa_high_j))) + (-((term.cc_high_i * cse_p14_t1_g462))) + (-((term.cc_high_j * cse_p14_t1_g462))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (cse_p11_t94_g100 + cse_p11_t94_g100 + (-2.00000000000000000e+0 * term.bb_high_j * cse_p14_t1_g463 * rr) + (cse_p11_t1_g121 * cse_p14_t1_g463 * ((-3.00000000000000000e+0 * term.aa_high_i) + (-3.00000000000000000e+0 * term.aa_high_j))) + (-((term.cc_high_i * cse_p14_t1_g463))) + (-((term.cc_high_j * cse_p14_t1_g463))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (cse_p11_t95_g101 + cse_p11_t95_g101 + (-2.00000000000000000e+0 * term.bb_high_j * cse_p14_t1_g464 * rr) + (cse_p11_t1_g121 * cse_p14_t1_g464 * ((-3.00000000000000000e+0 * term.aa_high_i) + (-3.00000000000000000e+0 * term.aa_high_j))) + (-((term.cc_high_i * cse_p14_t1_g464))) + (-((term.cc_high_j * cse_p14_t1_g464))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double de_dr = (term.cc_high_i + term.cc_high_j + cse_p12_t33_g161 + (2.00000000000000000e+0 * term.bb_high_i * rr) + (2.00000000000000000e+0 * term.bb_high_j * rr));
              double d2e_dr2 = ((2.00000000000000000e+0 * term.bb_high_i) + (2.00000000000000000e+0 * term.bb_high_j) + (6.00000000000000000e+0 * term.aa_high_i * rr) + (6.00000000000000000e+0 * term.aa_high_j * rr));
              double cse_p565_t2_invr39697 = (1.0 / (rr));
              double cse_p11_t1_g119 = (cse_p565_t2_invr39697 * cse_p565_t2_invr39697 * d2e_dr2);
              double cse_p565_t1_invr39696 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g120 = (cse_p565_t1_invr39696 * de_dr);
              double cse_p11_t1_g128 = (cse_p11_t1_g120 * (cse_p11_t1_g121 + (-(cse_p11_t53_g59))));
              double h_x1_x1 = (cse_p11_t1_g128 + (cse_p11_t1_g119 * dx * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double h_x1_y1 = ((cse_p11_t1_g119 * dx * dy) + (cse_p11_t1_g120 * dy * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g124 = (cse_p11_t1_g120 * dz);
              double h_x1_z1 = ((cse_p11_t1_g119 * dx * dz) + (cse_p11_t1_g124 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g122 = (cse_p11_t1_g119 * (-(dx)));
              double h_x1_x2 = ((cse_p11_t1_g120 * (cse_p11_t53_g59 + (-(cse_p11_t1_g121)))) + (cse_p11_t1_g122 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g123 = (cse_p11_t1_g119 * (-(dy)));
              double h_x1_y2 = ((cse_p11_t1_g120 * dx * dy) + (cse_p11_t1_g123 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g125 = (cse_p11_t1_g119 * (-(dz)));
              double h_x1_z2 = ((cse_p11_t1_g124 * dx) + (cse_p11_t1_g125 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g126 = (cse_p11_t1_g120 * (cse_p11_t1_g121 + (-(cse_p11_t54_g60))));
              double h_y1_y1 = (cse_p11_t1_g126 + (cse_p11_t1_g119 * dy * dy));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = ((cse_p11_t1_g119 * dy * dz) + (cse_p11_t1_g124 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = ((cse_p11_t1_g120 * dx * dy) + (cse_p11_t1_g122 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = ((cse_p11_t1_g120 * (cse_p11_t54_g60 + (-(cse_p11_t1_g121)))) + (cse_p11_t1_g123 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = ((cse_p11_t1_g124 * dy) + (cse_p11_t1_g125 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g127 = (cse_p11_t1_g120 * (cse_p11_t1_g121 + (-(cse_p11_t55_g61))));
              double h_z1_z1 = (cse_p11_t1_g127 + (cse_p11_t1_g119 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = ((cse_p11_t1_g122 * dz) + (cse_p11_t1_g124 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = ((cse_p11_t1_g123 * dz) + (cse_p11_t1_g124 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = ((cse_p11_t1_g120 * (cse_p11_t55_g61 + (-(cse_p11_t1_g121)))) + (cse_p11_t1_g125 * dz));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p11_t1_g128 + (cse_p11_t1_g122 * (-(dx))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = ((cse_p11_t1_g120 * dy * (-(dx))) + (cse_p11_t1_g122 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = ((cse_p11_t1_g122 * (-(dz))) + (cse_p11_t1_g124 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = (cse_p11_t1_g126 + (cse_p11_t1_g123 * (-(dy))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = ((cse_p11_t1_g123 * (-(dz))) + (cse_p11_t1_g124 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = (cse_p11_t1_g127 + (cse_p11_t1_g125 * (-(dz))));
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
                double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t53_g59))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t2_invr39697 = (1.0 / (rr));
                double cse_p565_t1_invr39696 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g121 = (rr * rr);
                double cse_p11_t1_g120 = (cse_p565_t1_invr39696 * de_dr);
                double cse_p11_t1_g119 = (cse_p565_t2_invr39697 * cse_p565_t2_invr39697 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g119 * dx * dy) + (cse_p11_t1_g120 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double h_x1_z1 = ((cse_p11_t1_g119 * dx * dz) + (cse_p11_t1_g120 * dz * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g124 = (cse_p11_t1_g120 * dz);
                double h_x1_x2 = ((cse_p11_t1_g119 * dx * (-(dx))) + (cse_p11_t1_g120 * (cse_p11_t53_g59 + (-(cse_p11_t1_g121)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g122 = (cse_p11_t1_g119 * (-(dx)));
                double h_x1_y2 = ((cse_p11_t1_g119 * dx * (-(dy))) + (cse_p11_t1_g120 * dx * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g123 = (cse_p11_t1_g119 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g119 * dx * (-(dz))) + (cse_p11_t1_g124 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g125 = (cse_p11_t1_g119 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g119 * dy * dy) + (cse_p11_t1_g120 * (cse_p11_t1_g121 + (-(cse_p11_t54_g60)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g126 = (cse_p11_t1_g120 * (cse_p11_t1_g121 + (-(cse_p11_t54_g60))));
                double h_y1_z1 = ((cse_p11_t1_g119 * dy * dz) + (cse_p11_t1_g124 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g120 * dx * dy) + (cse_p11_t1_g122 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g120 * (cse_p11_t54_g60 + (-(cse_p11_t1_g121)))) + (cse_p11_t1_g123 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g124 * dy) + (cse_p11_t1_g125 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g119 * dz * dz) + (cse_p11_t1_g120 * (cse_p11_t1_g121 + (-(cse_p11_t55_g61)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g127 = (cse_p11_t1_g120 * (cse_p11_t1_g121 + (-(cse_p11_t55_g61))));
                double h_z1_x2 = ((cse_p11_t1_g122 * dz) + (cse_p11_t1_g124 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g123 * dz) + (cse_p11_t1_g124 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g120 * (cse_p11_t55_g61 + (-(cse_p11_t1_g121)))) + (cse_p11_t1_g125 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g120 * (cse_p11_t1_g121 + (-(cse_p11_t53_g59)))) + (cse_p11_t1_g122 * (-(dx))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p11_t1_g128 = (cse_p11_t1_g120 * (cse_p11_t1_g121 + (-(cse_p11_t53_g59))));
                double h_x2_y2 = ((cse_p11_t1_g120 * dy * (-(dx))) + (cse_p11_t1_g122 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double h_x2_z2 = ((cse_p11_t1_g122 * (-(dz))) + (cse_p11_t1_g124 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g126 + (cse_p11_t1_g123 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g123 * (-(dz))) + (cse_p11_t1_g124 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g127 + (cse_p11_t1_g125 * (-(dz))));
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
