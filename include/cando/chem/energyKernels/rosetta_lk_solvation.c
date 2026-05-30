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
      double cse_p1_t1_sqrt316873 = sqrt(r2);
      double cse_p1_t2_invsqrt316874 = (1.0 / (cse_p1_t1_sqrt316873));
      double rr = cse_p1_t1_sqrt316873;
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
          double cse_p528_t4_g326291 = (2.00000000000000000e+0 * rr);
          double cse_p528_t1_g326306 = (cse_p1_t2_invsqrt316874 * dx);
          double cse_p528_t6_g326293 = (cse_p1_t2_invsqrt316874 * dy);
          double cse_p528_t7_g326294 = (cse_p1_t2_invsqrt316874 * dz);
          double cse_p528_t8_g326295 = (cse_p1_t2_invsqrt316874 * rr);
          double cse_p528_t9_g326296 = (rr * rr);
          double cse_p14_t2_g316991 = (3.00000000000000000e+0 * term.aa_low_j);
          double cse_p12_t1_g316966 = (3.00000000000000000e+0 * term.aa_low_i);
          double cse_p13_t11_g316977 = (cse_p12_t1_g316966 + cse_p14_t2_g316991);
          double cse_p12_t10_g316935 = (term.bb_low_i * rr);
          double cse_p12_t11_g316936 = (term.bb_low_j * rr);
          double cse_p12_t27_g316952 = (cse_p13_t11_g316977 * cse_p528_t9_g326296);
          double cse_p13_t1_g316981 = (term.bb_low_i * cse_p528_t8_g326295);
          double cse_p12_t7_g316932 = (cse_p13_t1_g316981 * dx);
          double cse_p12_t8_g316933 = (cse_p13_t1_g316981 * dy);
          double cse_p12_t9_g316934 = (cse_p13_t1_g316981 * dz);
          double cse_p12_t12_g316937 = (cse_p528_t8_g326295 * rr);
          double cse_p11_t27_g316905 = (term.cc_low_i + term.cc_low_j + cse_p12_t10_g316935 + cse_p12_t10_g316935 + cse_p12_t11_g316936 + cse_p12_t11_g316936 + cse_p12_t27_g316952);
          double cse_p11_t31_g316909 = (-(cse_p12_t7_g316932));
          double cse_p11_t32_g316910 = (-(cse_p12_t8_g316933));
          double cse_p11_t33_g316911 = (-(cse_p12_t9_g316934));
          double energy = (term.dd_low_i + term.dd_low_j + (term.cc_low_i * rr) + (term.cc_low_j * rr) + (cse_p528_t9_g326296 * (term.bb_low_i + term.bb_low_j + (term.aa_low_i * rr) + (term.aa_low_j * rr))));
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (cse_p11_t27_g316905 * cse_p528_t1_g326306);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t27_g316905 * cse_p528_t6_g326293);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t27_g316905 * cse_p528_t7_g326294);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double cse_p12_t1_g316958 = (cse_p12_t12_g316937 * dx);
          double cse_p528_t1_g326288 = (-3.00000000000000000e+0 * cse_p12_t1_g316958);
          double cse_p13_t1_g316982 = (term.bb_low_j * cse_p528_t8_g326295);
          double cse_p63_t1_g317758 = (-2.00000000000000000e+0 * cse_p13_t1_g316982);
          double g_x2 = (cse_p11_t31_g316909 + cse_p11_t31_g316909 + (term.aa_low_i * cse_p528_t1_g326288) + (term.aa_low_j * cse_p528_t1_g326288) + (cse_p63_t1_g317758 * dx) + (-((term.cc_low_i * cse_p528_t1_g326306))) + (-((term.cc_low_j * cse_p528_t1_g326306))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double cse_p12_t1_g316959 = (cse_p12_t12_g316937 * dy);
          double cse_p528_t2_g326289 = (-3.00000000000000000e+0 * cse_p12_t1_g316959);
          double g_y2 = (cse_p11_t32_g316910 + cse_p11_t32_g316910 + (term.aa_low_i * cse_p528_t2_g326289) + (term.aa_low_j * cse_p528_t2_g326289) + (cse_p63_t1_g317758 * dy) + (-((term.cc_low_i * cse_p528_t6_g326293))) + (-((term.cc_low_j * cse_p528_t6_g326293))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double cse_p12_t1_g316960 = (cse_p12_t12_g316937 * dz);
          double cse_p528_t3_g326290 = (-3.00000000000000000e+0 * cse_p12_t1_g316960);
          double g_z2 = (cse_p11_t33_g316911 + cse_p11_t33_g316911 + (term.aa_low_i * cse_p528_t3_g326290) + (term.aa_low_j * cse_p528_t3_g326290) + (cse_p63_t1_g317758 * dz) + (-((term.cc_low_i * cse_p528_t7_g326294))) + (-((term.cc_low_j * cse_p528_t7_g326294))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (term.cc_low_i + term.cc_low_j + cse_p12_t27_g316952 + (term.bb_low_i * cse_p528_t4_g326291) + (term.bb_low_j * cse_p528_t4_g326291));
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= term.r2_high)) {
            /* !COND-3 */
            double cse_p1_t3_invr316875 = (1.0 / (rr));
            double cse_p1_t4_invr2316876 = (cse_p1_t3_invr316875 * cse_p1_t3_invr316875);
            double inv_rr = cse_p1_t3_invr316875;
            double inv_rr2 = (inv_rr * inv_rr);
            double di = (rr + (-(term.lj_radius_i)));
            double cse_p14_t7_g316996 = (di * di);
            double cse_p13_t12_g316978 = (-(cse_p14_t7_g316996));
            double cse_p1_t5_invr316877 = (1.0 / (term.lambda2_i));
            double cse_p12_t20_g316945 = (cse_p13_t12_g316978 * cse_p1_t5_invr316877);
            double cse_p11_t36_g316914 = exp(cse_p12_t20_g316945);
            double exp_i = cse_p11_t36_g316914;
            double cse_p12_t24_g316949 = (exp_i * term.pref_i_j);
            double fi = (exp_i * inv_rr2 * term.pref_i_j);
            double hi = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t5_invr316877 * di)));
            double dfi = (fi * hi);
            double dj = (rr + (-(term.lj_radius_j)));
            double cse_p14_t8_g316997 = (dj * dj);
            double cse_p13_t13_g316979 = (-(cse_p14_t8_g316997));
            double cse_p1_t6_invr316878 = (1.0 / (term.lambda2_j));
            double cse_p12_t21_g316946 = (cse_p13_t13_g316979 * cse_p1_t6_invr316878);
            double cse_p11_t37_g316915 = exp(cse_p12_t21_g316946);
            double exp_j = cse_p11_t37_g316915;
            double cse_p12_t25_g316950 = (exp_j * term.pref_j_i);
            double cse_p11_t25_g316903 = (cse_p12_t24_g316949 + cse_p12_t25_g316950);
            double fj = (exp_j * inv_rr2 * term.pref_j_i);
            double hj = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t6_invr316878 * dj)));
            double dfj = (fj * hj);
            double energy = (fi + fj);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p11_t1_g316917 = (-2.00000000000000000e+0 * inv_rr2);
            double cse_p528_t1_g326306 = (cse_p1_t2_invsqrt316874 * dx);
            double cse_p11_t1_g316922 = (cse_p11_t1_g316917 * cse_p11_t36_g316914 * cse_p1_t5_invr316877 * term.pref_i_j);
            double cse_p11_t1_g316923 = (cse_p11_t1_g316917 * cse_p11_t37_g316915 * cse_p1_t6_invr316878 * term.pref_j_i);
            double cse_p11_t1_g316924 = (cse_p11_t25_g316903 * cse_p1_t4_invr2316876 * inv_rr);
            double cse_p12_t1_g316962 = (-2.00000000000000000e+0 * cse_p11_t1_g316924);
            double cse_p63_t1_g317755 = (cse_p11_t1_g316922 * di);
            double cse_p64_t13_g317774 = (cse_p528_t1_g326306 * cse_p63_t1_g317755);
            double cse_p63_t1_g317756 = (cse_p11_t1_g316923 * dj);
            double cse_p64_t14_g317775 = (cse_p528_t1_g326306 * cse_p63_t1_g317756);
            double g_x1 = (cse_p64_t13_g317774 + cse_p64_t14_g317775 + (cse_p12_t1_g316962 * cse_p528_t1_g326306));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p11_t1_g316920 = (cse_p1_t2_invsqrt316874 * dy);
            double cse_p64_t5_g317766 = (cse_p11_t1_g316920 * cse_p63_t1_g317755);
            double cse_p64_t6_g317767 = (cse_p11_t1_g316920 * cse_p63_t1_g317756);
            double g_y1 = (cse_p64_t5_g317766 + cse_p64_t6_g317767 + (cse_p11_t1_g316920 * cse_p12_t1_g316962));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p11_t1_g316921 = (cse_p1_t2_invsqrt316874 * dz);
            double cse_p64_t7_g317768 = (cse_p11_t1_g316921 * cse_p63_t1_g317755);
            double cse_p64_t8_g317769 = (cse_p11_t1_g316921 * cse_p63_t1_g317756);
            double g_z1 = (cse_p64_t7_g317768 + cse_p64_t8_g317769 + (cse_p11_t1_g316921 * cse_p12_t1_g316962));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double cse_p12_t1_g316963 = (2.00000000000000000e+0 * cse_p11_t1_g316924);
            double g_x2 = ((cse_p12_t1_g316963 * cse_p528_t1_g326306) + (-(cse_p64_t13_g317774)) + (-(cse_p64_t14_g317775)));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = ((cse_p11_t1_g316920 * cse_p12_t1_g316963) + (-(cse_p64_t5_g317766)) + (-(cse_p64_t6_g317767)));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = ((cse_p11_t1_g316921 * cse_p12_t1_g316963) + (-(cse_p64_t7_g317768)) + (-(cse_p64_t8_g317769)));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (dfi + dfj);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= term.r3_high)) {
              /* !COND-4 */
              double cse_p12_t1_g316964 = (3.00000000000000000e+0 * term.aa_high_i);
              double cse_p12_t1_g316965 = (3.00000000000000000e+0 * term.aa_high_j);
              double cse_p13_t10_g316976 = (cse_p12_t1_g316964 + cse_p12_t1_g316965);
              double cse_p12_t5_g316930 = (term.bb_high_i * rr);
              double cse_p12_t6_g316931 = (term.bb_high_j * rr);
              double cse_p12_t1_g316961 = (rr * rr);
              double cse_p12_t26_g316951 = (cse_p12_t1_g316961 * cse_p13_t10_g316976);
              double cse_p11_t1_g316916 = (cse_p1_t2_invsqrt316874 * rr);
              double cse_p11_t1_g316925 = (term.bb_high_i * cse_p11_t1_g316916);
              double cse_p12_t13_g316938 = (cse_p11_t1_g316925 * dx);
              double cse_p12_t14_g316939 = (cse_p11_t1_g316925 * dy);
              double cse_p12_t15_g316940 = (cse_p11_t1_g316925 * dz);
              double cse_p11_t26_g316904 = (term.cc_high_i + term.cc_high_j + cse_p12_t26_g316951 + cse_p12_t5_g316930 + cse_p12_t5_g316930 + cse_p12_t6_g316931 + cse_p12_t6_g316931);
              double cse_p11_t28_g316906 = (-(cse_p12_t13_g316938));
              double cse_p11_t29_g316907 = (-(cse_p12_t14_g316939));
              double cse_p11_t30_g316908 = (-(cse_p12_t15_g316940));
              double energy = (term.dd_high_i + term.dd_high_j + (term.cc_high_i * rr) + (term.cc_high_j * rr) + (cse_p12_t1_g316961 * (term.bb_high_i + term.bb_high_j + (term.aa_high_i * rr) + (term.aa_high_j * rr))));
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p528_t1_g326306 = (cse_p1_t2_invsqrt316874 * dx);
              double g_x1 = (cse_p11_t26_g316904 * cse_p528_t1_g326306);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double cse_p11_t1_g316920 = (cse_p1_t2_invsqrt316874 * dy);
              double g_y1 = (cse_p11_t1_g316920 * cse_p11_t26_g316904);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double cse_p11_t1_g316921 = (cse_p1_t2_invsqrt316874 * dz);
              double g_z1 = (cse_p11_t1_g316921 * cse_p11_t26_g316904);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p11_t1_g316918 = (cse_p11_t1_g316916 * rr);
              double cse_p12_t1_g316958 = (cse_p11_t1_g316918 * dx);
              double cse_p13_t1_g316980 = (term.bb_high_j * cse_p11_t1_g316916);
              double cse_p63_t1_g317752 = (-3.00000000000000000e+0 * cse_p12_t1_g316958);
              double cse_p63_t1_g317757 = (-2.00000000000000000e+0 * cse_p13_t1_g316980);
              double g_x2 = (cse_p11_t28_g316906 + cse_p11_t28_g316906 + (term.aa_high_i * cse_p63_t1_g317752) + (term.aa_high_j * cse_p63_t1_g317752) + (cse_p63_t1_g317757 * dx) + (-((term.cc_high_i * cse_p528_t1_g326306))) + (-((term.cc_high_j * cse_p528_t1_g326306))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double cse_p12_t1_g316959 = (cse_p11_t1_g316918 * dy);
              double cse_p63_t1_g317753 = (-3.00000000000000000e+0 * cse_p12_t1_g316959);
              double g_y2 = (cse_p11_t29_g316907 + cse_p11_t29_g316907 + (term.aa_high_i * cse_p63_t1_g317753) + (term.aa_high_j * cse_p63_t1_g317753) + (cse_p63_t1_g317757 * dy) + (-((term.cc_high_i * cse_p11_t1_g316920))) + (-((term.cc_high_j * cse_p11_t1_g316920))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double cse_p12_t1_g316960 = (cse_p11_t1_g316918 * dz);
              double cse_p63_t1_g317754 = (-3.00000000000000000e+0 * cse_p12_t1_g316960);
              double g_z2 = (cse_p11_t30_g316908 + cse_p11_t30_g316908 + (term.aa_high_i * cse_p63_t1_g317754) + (term.aa_high_j * cse_p63_t1_g317754) + (cse_p63_t1_g317757 * dz) + (-((term.cc_high_i * cse_p11_t1_g316921))) + (-((term.cc_high_j * cse_p11_t1_g316921))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double cse_p12_t1_g316957 = (2.00000000000000000e+0 * rr);
              double de_dr = (term.cc_high_i + term.cc_high_j + cse_p12_t26_g316951 + (term.bb_high_i * cse_p12_t1_g316957) + (term.bb_high_j * cse_p12_t1_g316957));
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
            double cse_p61_t1_invr326859 = (1.0 / (term.lambda2_i));
            double exp_i = exp((-((cse_p61_t1_invr326859 * di * di))));
            double fi = (exp_i * inv_rr2 * term.pref_i_j);
            double hi = (-2.00000000000000000e+0 * (inv_rr + (cse_p61_t1_invr326859 * di)));
            double dj = (rr + (-(term.lj_radius_j)));
            double cse_p61_t2_invr326860 = (1.0 / (term.lambda2_j));
            double exp_j = exp((-((cse_p61_t2_invr326860 * dj * dj))));
            double fj = (exp_j * inv_rr2 * term.pref_j_i);
            double hj = (-2.00000000000000000e+0 * (inv_rr + (cse_p61_t2_invr326860 * dj)));
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
      double cse_p11_t53_g331869 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t54_g331870 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t55_g331871 = (dz * dz);
      double r2 = (cse_p11_t53_g331869 + cse_p11_t54_g331870 + cse_p11_t55_g331871);
      double cse_p1_t1_sqrt331811 = sqrt(r2);
      double cse_p1_t2_invsqrt331812 = (1.0 / (cse_p1_t1_sqrt331811));
      double rr = cse_p1_t1_sqrt331811;
      if ((rr <= term.r0_low)) {
        /* !COND-1 */
        double cse_p60_t43_g335918 = (rr * rr);
        double cse_p60_t71_g335946 = (-(cse_p11_t53_g331869));
        double cse_p60_t72_g335947 = (-(cse_p11_t54_g331870));
        double cse_p60_t73_g335948 = (-(cse_p11_t55_g331871));
        double cse_p60_t75_g335950 = (-(dx));
        double cse_p60_t76_g335951 = (-(dy));
        double cse_p60_t77_g335952 = (-(dz));
        double cse_p565_t1_invr371506 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr371507 = (1.0 / (rr));
        double cse_p60_t80_g335955 = (rr * rr * rr);
        double cse_p15_t75_g332349 = (-(cse_p60_t43_g335918));
        double cse_p14_t82_g332247 = (cse_p11_t53_g331869 + cse_p15_t75_g332349);
        double cse_p14_t84_g332249 = (cse_p11_t54_g331870 + cse_p15_t75_g332349);
        double cse_p14_t86_g332251 = (cse_p11_t55_g331871 + cse_p15_t75_g332349);
        double cse_p13_t109_g332150 = (cse_p60_t43_g335918 + cse_p60_t71_g335946);
        double cse_p13_t110_g332151 = (cse_p60_t43_g335918 + cse_p60_t72_g335947);
        double cse_p13_t111_g332152 = (cse_p60_t43_g335918 + cse_p60_t73_g335948);
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
        double cse_p60_t42_g335917 = (cse_p565_t1_invr371506 * de_dr);
        double cse_p13_t67_g332108 = (cse_p14_t82_g332247 * cse_p60_t42_g335917);
        double cse_p13_t68_g332109 = (cse_p14_t84_g332249 * cse_p60_t42_g335917);
        double cse_p13_t69_g332110 = (cse_p14_t86_g332251 * cse_p60_t42_g335917);
        double cse_p13_t70_g332111 = (cse_p60_t42_g335917 * dx * dy);
        double cse_p13_t71_g332112 = (cse_p60_t42_g335917 * cse_p60_t75_g335950 * dy);
        double cse_p12_t27_g331965 = (cse_p13_t109_g332150 * cse_p60_t42_g335917);
        double cse_p12_t28_g331966 = (cse_p13_t110_g332151 * cse_p60_t42_g335917);
        double cse_p12_t29_g331967 = (cse_p13_t111_g332152 * cse_p60_t42_g335917);
        double cse_p12_t32_g331970 = (cse_p60_t42_g335917 * dz);
        double cse_p13_t57_g332098 = (cse_p12_t32_g331970 * cse_p60_t75_g335950);
        double cse_p13_t58_g332099 = (cse_p12_t32_g331970 * cse_p60_t76_g335951);
        double cse_p13_t59_g332100 = (cse_p12_t32_g331970 * dx);
        double cse_p13_t60_g332101 = (cse_p12_t32_g331970 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p60_t41_g335916 = (cse_p565_t2_invr371507 * cse_p565_t2_invr371507 * d2e_dr2);
        double cse_p13_t61_g332102 = (cse_p60_t41_g335916 * dx * dx);
        double cse_p13_t62_g332103 = (cse_p60_t41_g335916 * dx * dy);
        double cse_p13_t63_g332104 = (cse_p60_t41_g335916 * dx * dz);
        double cse_p13_t64_g332105 = (cse_p60_t41_g335916 * dy * dy);
        double cse_p13_t65_g332106 = (cse_p60_t41_g335916 * dy * dz);
        double cse_p13_t66_g332107 = (cse_p60_t41_g335916 * dz * dz);
        double cse_p12_t15_g331953 = (cse_p60_t41_g335916 * cse_p60_t75_g335950);
        double cse_p13_t42_g332083 = (cse_p12_t15_g331953 * cse_p60_t75_g335950);
        double cse_p13_t43_g332084 = (cse_p12_t15_g331953 * cse_p60_t76_g335951);
        double cse_p13_t44_g332085 = (cse_p12_t15_g331953 * cse_p60_t77_g335952);
        double cse_p13_t45_g332086 = (cse_p12_t15_g331953 * dx);
        double cse_p13_t46_g332087 = (cse_p12_t15_g331953 * dy);
        double cse_p13_t47_g332088 = (cse_p12_t15_g331953 * dz);
        double cse_p12_t16_g331954 = (cse_p60_t41_g335916 * cse_p60_t76_g335951);
        double cse_p13_t48_g332089 = (cse_p12_t16_g331954 * cse_p60_t76_g335951);
        double cse_p13_t49_g332090 = (cse_p12_t16_g331954 * cse_p60_t77_g335952);
        double cse_p13_t50_g332091 = (cse_p12_t16_g331954 * dx);
        double cse_p13_t51_g332092 = (cse_p12_t16_g331954 * dy);
        double cse_p13_t52_g332093 = (cse_p12_t16_g331954 * dz);
        double cse_p12_t17_g331955 = (cse_p60_t41_g335916 * cse_p60_t77_g335952);
        double cse_p13_t53_g332094 = (cse_p12_t17_g331955 * cse_p60_t77_g335952);
        double cse_p13_t54_g332095 = (cse_p12_t17_g331955 * dx);
        double cse_p13_t55_g332096 = (cse_p12_t17_g331955 * dy);
        double cse_p13_t56_g332097 = (cse_p12_t17_g331955 * dz);
        double cse_p12_t68_g332006 = (cse_p13_t62_g332103 + cse_p13_t71_g332112);
        double cse_p12_t69_g332007 = (cse_p13_t57_g332098 + cse_p13_t63_g332104);
        double cse_p12_t70_g332008 = (cse_p13_t58_g332099 + cse_p13_t65_g332106);
        double cse_p12_t86_g332024 = (cse_p12_t28_g331966 + cse_p13_t64_g332105);
        double cse_p12_t88_g332026 = (cse_p12_t29_g331967 + cse_p13_t66_g332107);
        double cse_p12_t90_g332028 = (cse_p12_t27_g331965 + cse_p13_t61_g332102);
        double cse_p12_t63_g332001 = (cse_p13_t45_g332086 + cse_p13_t67_g332108);
        double cse_p12_t71_g332009 = (cse_p13_t46_g332087 + cse_p13_t70_g332111);
        double cse_p12_t73_g332011 = (cse_p13_t43_g332084 + cse_p13_t71_g332112);
        double cse_p12_t74_g332012 = (cse_p13_t44_g332085 + cse_p13_t57_g332098);
        double cse_p12_t75_g332013 = (cse_p13_t47_g332088 + cse_p13_t59_g332100);
        double cse_p12_t91_g332029 = (cse_p12_t27_g331965 + cse_p13_t42_g332083);
        double cse_p12_t64_g332002 = (cse_p13_t51_g332092 + cse_p13_t68_g332109);
        double cse_p12_t72_g332010 = (cse_p13_t50_g332091 + cse_p13_t70_g332111);
        double cse_p12_t76_g332014 = (cse_p13_t49_g332090 + cse_p13_t58_g332099);
        double cse_p12_t77_g332015 = (cse_p13_t52_g332093 + cse_p13_t60_g332101);
        double cse_p12_t87_g332025 = (cse_p12_t28_g331966 + cse_p13_t48_g332089);
        double cse_p12_t65_g332003 = (cse_p13_t56_g332097 + cse_p13_t69_g332110);
        double cse_p12_t78_g332016 = (cse_p13_t54_g332095 + cse_p13_t59_g332100);
        double cse_p12_t79_g332017 = (cse_p13_t55_g332096 + cse_p13_t60_g332101);
        double cse_p12_t89_g332027 = (cse_p12_t29_g331967 + cse_p13_t53_g332094);
        double h_x1_x1 = cse_p12_t90_g332028;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t68_g332006;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t69_g332007;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t63_g332001;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t72_g332010;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t78_g332016;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t86_g332024;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t70_g332008;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t71_g332009;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t64_g332002;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t79_g332017;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t88_g332026;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t75_g332013;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t77_g332015;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t65_g332003;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t91_g332029;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t73_g332011;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t74_g332012;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t87_g332025;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t76_g332014;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t89_g332027;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.r1_low)) {
          /* !COND-2 */
          double cse_p60_t38_g335913 = (cse_p1_t2_invsqrt331812 * dx);
          double cse_p60_t39_g335914 = (cse_p1_t2_invsqrt331812 * dy);
          double cse_p60_t40_g335915 = (cse_p1_t2_invsqrt331812 * dz);
          double cse_p14_t3_g332168 = (3.00000000000000000e+0 * term.aa_low_i);
          double cse_p14_t4_g332169 = (3.00000000000000000e+0 * term.aa_low_j);
          double cse_p13_t81_g332122 = (cse_p14_t3_g332168 + cse_p14_t4_g332169);
          double cse_p12_t10_g331948 = (term.bb_low_i * cse_p60_t38_g335913 * rr);
          double cse_p12_t11_g331949 = (term.bb_low_i * cse_p60_t39_g335914 * rr);
          double cse_p12_t12_g331950 = (term.bb_low_i * cse_p60_t40_g335915 * rr);
          double cse_p12_t13_g331951 = (term.bb_low_i * rr);
          double cse_p12_t14_g331952 = (term.bb_low_j * rr);
          double cse_p11_t1_g331931 = (rr * rr);
          double cse_p60_t59_g335934 = (cse_p11_t53_g331869 + (-(cse_p11_t1_g331931)));
          double cse_p60_t60_g335935 = (cse_p11_t54_g331870 + (-(cse_p11_t1_g331931)));
          double cse_p60_t61_g335936 = (cse_p11_t55_g331871 + (-(cse_p11_t1_g331931)));
          double cse_p60_t62_g335937 = (cse_p11_t1_g331931 + (-(cse_p11_t53_g331869)));
          double cse_p60_t63_g335938 = (cse_p11_t1_g331931 + (-(cse_p11_t54_g331870)));
          double cse_p60_t64_g335939 = (cse_p11_t1_g331931 + (-(cse_p11_t55_g331871)));
          double cse_p60_t74_g335949 = (-(cse_p11_t1_g331931));
          double cse_p12_t34_g331972 = (cse_p11_t1_g331931 * cse_p13_t81_g332122);
          double cse_p11_t92_g331908 = (term.cc_low_i + term.cc_low_j + cse_p12_t13_g331951 + cse_p12_t13_g331951 + cse_p12_t14_g331952 + cse_p12_t14_g331952 + cse_p12_t34_g331972);
          double cse_p11_t96_g331912 = (-(cse_p12_t10_g331948));
          double cse_p11_t97_g331913 = (-(cse_p12_t11_g331949));
          double cse_p11_t98_g331914 = (-(cse_p12_t12_g331950));
          double energy = (term.dd_low_i + term.dd_low_j + (term.cc_low_i * rr) + (term.cc_low_j * rr) + (cse_p11_t1_g331931 * (term.bb_low_i + term.bb_low_j + (term.aa_low_i * rr) + (term.aa_low_j * rr))));
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (cse_p11_t92_g331908 * cse_p60_t38_g335913);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t92_g331908 * cse_p60_t39_g335914);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t92_g331908 * cse_p60_t40_g335915);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p11_t96_g331912 + cse_p11_t96_g331912 + (-2.00000000000000000e+0 * term.bb_low_j * cse_p60_t38_g335913 * rr) + (cse_p11_t1_g331931 * cse_p60_t38_g335913 * ((-3.00000000000000000e+0 * term.aa_low_i) + (-3.00000000000000000e+0 * term.aa_low_j))) + (-((term.cc_low_i * cse_p60_t38_g335913))) + (-((term.cc_low_j * cse_p60_t38_g335913))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p11_t97_g331913 + cse_p11_t97_g331913 + (-2.00000000000000000e+0 * term.bb_low_j * cse_p60_t39_g335914 * rr) + (cse_p11_t1_g331931 * cse_p60_t39_g335914 * ((-3.00000000000000000e+0 * term.aa_low_i) + (-3.00000000000000000e+0 * term.aa_low_j))) + (-((term.cc_low_i * cse_p60_t39_g335914))) + (-((term.cc_low_j * cse_p60_t39_g335914))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p11_t98_g331914 + cse_p11_t98_g331914 + (-2.00000000000000000e+0 * term.bb_low_j * cse_p60_t40_g335915 * rr) + (cse_p11_t1_g331931 * cse_p60_t40_g335915 * ((-3.00000000000000000e+0 * term.aa_low_i) + (-3.00000000000000000e+0 * term.aa_low_j))) + (-((term.cc_low_i * cse_p60_t40_g335915))) + (-((term.cc_low_j * cse_p60_t40_g335915))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (term.cc_low_i + term.cc_low_j + cse_p12_t34_g331972 + (2.00000000000000000e+0 * term.bb_low_i * rr) + (2.00000000000000000e+0 * term.bb_low_j * rr));
          double d2e_dr2 = ((2.00000000000000000e+0 * term.bb_low_i) + (2.00000000000000000e+0 * term.bb_low_j) + (6.00000000000000000e+0 * term.aa_low_i * rr) + (6.00000000000000000e+0 * term.aa_low_j * rr));
          double cse_p565_t2_invr371507 = (1.0 / (rr));
          double cse_p11_t1_g331929 = (cse_p565_t2_invr371507 * cse_p565_t2_invr371507 * d2e_dr2);
          double cse_p60_t1_g335876 = (cse_p11_t1_g331929 * (-(dx)));
          double cse_p60_t2_g335877 = (cse_p11_t1_g331929 * (-(dy)));
          double cse_p60_t3_g335878 = (cse_p11_t1_g331929 * (-(dz)));
          double cse_p60_t4_g335879 = (cse_p11_t1_g331929 * dx * dx);
          double cse_p60_t5_g335880 = (cse_p11_t1_g331929 * dx * dy);
          double cse_p60_t6_g335881 = (cse_p11_t1_g331929 * dx * dz);
          double cse_p60_t7_g335882 = (cse_p11_t1_g331929 * dy * dy);
          double cse_p60_t8_g335883 = (cse_p11_t1_g331929 * dy * dz);
          double cse_p60_t9_g335884 = (cse_p11_t1_g331929 * dz * dz);
          double cse_p14_t46_g332211 = (cse_p60_t1_g335876 * (-(dx)));
          double cse_p14_t47_g332212 = (cse_p60_t1_g335876 * (-(dy)));
          double cse_p14_t48_g332213 = (cse_p60_t1_g335876 * (-(dz)));
          double cse_p14_t49_g332214 = (cse_p60_t1_g335876 * dx);
          double cse_p14_t50_g332215 = (cse_p60_t1_g335876 * dy);
          double cse_p14_t51_g332216 = (cse_p60_t1_g335876 * dz);
          double cse_p14_t52_g332217 = (cse_p60_t2_g335877 * (-(dy)));
          double cse_p14_t53_g332218 = (cse_p60_t2_g335877 * (-(dz)));
          double cse_p14_t54_g332219 = (cse_p60_t2_g335877 * dx);
          double cse_p14_t55_g332220 = (cse_p60_t2_g335877 * dy);
          double cse_p14_t56_g332221 = (cse_p60_t2_g335877 * dz);
          double cse_p14_t57_g332222 = (cse_p60_t3_g335878 * (-(dz)));
          double cse_p14_t58_g332223 = (cse_p60_t3_g335878 * dx);
          double cse_p14_t59_g332224 = (cse_p60_t3_g335878 * dy);
          double cse_p14_t60_g332225 = (cse_p60_t3_g335878 * dz);
          double cse_p565_t1_invr371506 = (1.0 / ((rr * rr * rr)));
          double cse_p11_t1_g331930 = (cse_p565_t1_invr371506 * de_dr);
          double cse_p60_t10_g335885 = (cse_p11_t1_g331930 * (cse_p11_t53_g331869 + cse_p60_t74_g335949));
          double cse_p60_t11_g335886 = (cse_p11_t1_g331930 * (cse_p11_t54_g331870 + cse_p60_t74_g335949));
          double cse_p60_t12_g335887 = (cse_p11_t1_g331930 * (cse_p11_t55_g331871 + cse_p60_t74_g335949));
          double cse_p60_t13_g335888 = (cse_p11_t1_g331930 * cse_p60_t62_g335937);
          double cse_p60_t14_g335889 = (cse_p11_t1_g331930 * cse_p60_t63_g335938);
          double cse_p60_t15_g335890 = (cse_p11_t1_g331930 * cse_p60_t64_g335939);
          double cse_p60_t16_g335891 = (cse_p11_t1_g331930 * dx * dy);
          double cse_p60_t17_g335892 = (cse_p11_t1_g331930 * dy * (-(dx)));
          double cse_p60_t18_g335893 = (cse_p11_t1_g331930 * dz);
          double cse_p60_t47_g335922 = (cse_p60_t17_g335892 + cse_p60_t5_g335880);
          double cse_p14_t42_g332207 = (cse_p60_t18_g335893 * (-(dx)));
          double cse_p14_t43_g332208 = (cse_p60_t18_g335893 * (-(dy)));
          double cse_p14_t44_g332209 = (cse_p60_t18_g335893 * dx);
          double cse_p14_t45_g332210 = (cse_p60_t18_g335893 * dy);
          double cse_p13_t107_g332148 = (cse_p60_t13_g335888 + cse_p60_t4_g335879);
          double h_x1_x1 = cse_p13_t107_g332148;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t47_g335922;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t83_g332124 = (cse_p14_t42_g332207 + cse_p60_t6_g335881);
          double cse_p13_t84_g332125 = (cse_p14_t43_g332208 + cse_p60_t8_g335883);
          double h_x1_z1 = cse_p13_t83_g332124;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t77_g332118 = (cse_p14_t49_g332214 + cse_p60_t10_g335885);
          double cse_p13_t85_g332126 = (cse_p14_t50_g332215 + cse_p60_t16_g335891);
          double cse_p13_t87_g332128 = (cse_p14_t47_g332212 + cse_p60_t17_g335892);
          double cse_p13_t88_g332129 = (cse_p14_t42_g332207 + cse_p14_t48_g332213);
          double cse_p13_t89_g332130 = (cse_p14_t44_g332209 + cse_p14_t51_g332216);
          double cse_p13_t108_g332149 = (cse_p14_t46_g332211 + cse_p60_t13_g335888);
          double h_x1_x2 = cse_p13_t77_g332118;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t78_g332119 = (cse_p14_t55_g332220 + cse_p60_t11_g335886);
          double cse_p13_t86_g332127 = (cse_p14_t54_g332219 + cse_p60_t16_g335891);
          double cse_p13_t90_g332131 = (cse_p14_t43_g332208 + cse_p14_t53_g332218);
          double cse_p13_t91_g332132 = (cse_p14_t45_g332210 + cse_p14_t56_g332221);
          double h_x1_y2 = cse_p13_t86_g332127;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t79_g332120 = (cse_p14_t60_g332225 + cse_p60_t12_g335887);
          double cse_p13_t92_g332133 = (cse_p14_t44_g332209 + cse_p14_t58_g332223);
          double cse_p13_t93_g332134 = (cse_p14_t45_g332210 + cse_p14_t59_g332224);
          double h_x1_z2 = cse_p13_t92_g332133;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t103_g332144 = (cse_p60_t14_g335889 + cse_p60_t7_g335882);
          double cse_p13_t104_g332145 = (cse_p14_t52_g332217 + cse_p60_t14_g335889);
          double h_y1_y1 = cse_p13_t103_g332144;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t84_g332125;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t85_g332126;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t78_g332119;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t93_g332134;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t105_g332146 = (cse_p60_t15_g335890 + cse_p60_t9_g335884);
          double cse_p13_t106_g332147 = (cse_p14_t57_g332222 + cse_p60_t15_g335890);
          double h_z1_z1 = cse_p13_t105_g332146;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t89_g332130;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t91_g332132;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t79_g332120;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t108_g332149;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t87_g332128;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t88_g332129;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t104_g332145;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t90_g332131;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t106_g332147;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= term.r2_high)) {
            /* !COND-3 */
            double cse_p565_t2_invr371507 = (1.0 / (rr));
            double cse_p1_t4_invr2331814 = (cse_p565_t2_invr371507 * cse_p565_t2_invr371507);
            double inv_rr = cse_p565_t2_invr371507;
            double inv_rr2 = (inv_rr * inv_rr);
            double cse_p11_t1_g331817 = (2.00000000000000000e+0 * inv_rr2);
            double di = (rr + (-(term.lj_radius_i)));
            double cse_p14_t63_g332228 = (di * di);
            double cse_p13_t112_g332153 = (-(cse_p14_t63_g332228));
            double cse_p1_t5_invr331815 = (1.0 / (term.lambda2_i));
            double cse_p12_t54_g331992 = (cse_p13_t112_g332153 * cse_p1_t5_invr331815);
            double cse_p11_t108_g331924 = exp(cse_p12_t54_g331992);
            double exp_i = cse_p11_t108_g331924;
            double cse_p12_t60_g331998 = (exp_i * term.pref_i_j);
            double fi = (exp_i * inv_rr2 * term.pref_i_j);
            double hi = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t5_invr331815 * di)));
            double dfi = (fi * hi);
            double d2fi = (fi * (cse_p11_t1_g331817 + (-2.00000000000000000e+0 * cse_p1_t5_invr331815) + (hi * hi)));
            double dj = (rr + (-(term.lj_radius_j)));
            double cse_p14_t64_g332229 = (dj * dj);
            double cse_p13_t113_g332154 = (-(cse_p14_t64_g332229));
            double cse_p1_t6_invr331816 = (1.0 / (term.lambda2_j));
            double cse_p12_t55_g331993 = (cse_p13_t113_g332154 * cse_p1_t6_invr331816);
            double cse_p11_t109_g331925 = exp(cse_p12_t55_g331993);
            double exp_j = cse_p11_t109_g331925;
            double cse_p12_t61_g331999 = (exp_j * term.pref_j_i);
            double cse_p11_t87_g331903 = (cse_p12_t60_g331998 + cse_p12_t61_g331999);
            double fj = (exp_j * inv_rr2 * term.pref_j_i);
            double hj = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t6_invr331816 * dj)));
            double dfj = (fj * hj);
            double d2fj = (fj * (cse_p11_t1_g331817 + (-2.00000000000000000e+0 * cse_p1_t6_invr331816) + (hj * hj)));
            double energy = (fi + fj);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p14_t1_g332272 = (cse_p1_t2_invsqrt331812 * dx);
            double g_x1 = ((-2.00000000000000000e+0 * cse_p11_t87_g331903 * cse_p14_t1_g332272 * cse_p1_t4_invr2331814 * inv_rr) + (cse_p14_t1_g332272 * inv_rr2 * ((-2.00000000000000000e+0 * cse_p11_t108_g331924 * cse_p1_t5_invr331815 * di * term.pref_i_j) + (-2.00000000000000000e+0 * cse_p11_t109_g331925 * cse_p1_t6_invr331816 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p14_t1_g332273 = (cse_p1_t2_invsqrt331812 * dy);
            double g_y1 = ((-2.00000000000000000e+0 * cse_p11_t87_g331903 * cse_p14_t1_g332273 * cse_p1_t4_invr2331814 * inv_rr) + (cse_p14_t1_g332273 * inv_rr2 * ((-2.00000000000000000e+0 * cse_p11_t108_g331924 * cse_p1_t5_invr331815 * di * term.pref_i_j) + (-2.00000000000000000e+0 * cse_p11_t109_g331925 * cse_p1_t6_invr331816 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p14_t1_g332274 = (cse_p1_t2_invsqrt331812 * dz);
            double g_z1 = ((-2.00000000000000000e+0 * cse_p11_t87_g331903 * cse_p14_t1_g332274 * cse_p1_t4_invr2331814 * inv_rr) + (cse_p14_t1_g332274 * inv_rr2 * ((-2.00000000000000000e+0 * cse_p11_t108_g331924 * cse_p1_t5_invr331815 * di * term.pref_i_j) + (-2.00000000000000000e+0 * cse_p11_t109_g331925 * cse_p1_t6_invr331816 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = ((2.00000000000000000e+0 * cse_p11_t87_g331903 * cse_p14_t1_g332272 * cse_p1_t4_invr2331814 * inv_rr) + (cse_p14_t1_g332272 * inv_rr2 * ((2.00000000000000000e+0 * cse_p11_t108_g331924 * cse_p1_t5_invr331815 * di * term.pref_i_j) + (2.00000000000000000e+0 * cse_p11_t109_g331925 * cse_p1_t6_invr331816 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = ((2.00000000000000000e+0 * cse_p11_t87_g331903 * cse_p14_t1_g332273 * cse_p1_t4_invr2331814 * inv_rr) + (cse_p14_t1_g332273 * inv_rr2 * ((2.00000000000000000e+0 * cse_p11_t108_g331924 * cse_p1_t5_invr331815 * di * term.pref_i_j) + (2.00000000000000000e+0 * cse_p11_t109_g331925 * cse_p1_t6_invr331816 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = ((2.00000000000000000e+0 * cse_p11_t87_g331903 * cse_p14_t1_g332274 * cse_p1_t4_invr2331814 * inv_rr) + (cse_p14_t1_g332274 * inv_rr2 * ((2.00000000000000000e+0 * cse_p11_t108_g331924 * cse_p1_t5_invr331815 * di * term.pref_i_j) + (2.00000000000000000e+0 * cse_p11_t109_g331925 * cse_p1_t6_invr331816 * dj * term.pref_j_i))));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (dfi + dfj);
            double d2e_dr2 = (d2fi + d2fj);
            double cse_p11_t1_g331929 = (cse_p565_t2_invr371507 * cse_p565_t2_invr371507 * d2e_dr2);
            double cse_p565_t1_invr371506 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g331930 = (cse_p565_t1_invr371506 * de_dr);
            double cse_p11_t1_g331931 = (rr * rr);
            double cse_p11_t1_g331938 = (cse_p11_t1_g331930 * (cse_p11_t1_g331931 + (-(cse_p11_t53_g331869))));
            double cse_p60_t69_g335944 = (cse_p11_t1_g331938 + (cse_p11_t1_g331929 * dx * dx));
            double h_x1_x1 = cse_p60_t69_g335944;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = ((cse_p11_t1_g331929 * dx * dy) + (cse_p11_t1_g331930 * dy * (-(dx))));
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p11_t1_g331934 = (cse_p11_t1_g331930 * dz);
            double cse_p60_t30_g335905 = (cse_p11_t1_g331934 * (-(dx)));
            double cse_p60_t31_g335906 = (cse_p11_t1_g331934 * (-(dy)));
            double cse_p60_t32_g335907 = (cse_p11_t1_g331934 * dx);
            double cse_p60_t33_g335908 = (cse_p11_t1_g331934 * dy);
            double cse_p60_t48_g335923 = (cse_p60_t30_g335905 + (cse_p11_t1_g331929 * dx * dz));
            double cse_p60_t49_g335924 = (cse_p60_t31_g335906 + (cse_p11_t1_g331929 * dy * dz));
            double h_x1_z1 = cse_p60_t48_g335923;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p11_t1_g331932 = (cse_p11_t1_g331929 * (-(dx)));
            double cse_p60_t19_g335894 = (cse_p11_t1_g331932 * (-(dx)));
            double cse_p60_t20_g335895 = (cse_p11_t1_g331932 * (-(dy)));
            double cse_p60_t21_g335896 = (cse_p11_t1_g331932 * (-(dz)));
            double cse_p60_t22_g335897 = (cse_p11_t1_g331932 * dx);
            double cse_p60_t23_g335898 = (cse_p11_t1_g331932 * dy);
            double cse_p60_t24_g335899 = (cse_p11_t1_g331932 * dz);
            double cse_p60_t44_g335919 = (cse_p60_t22_g335897 + (cse_p11_t1_g331930 * (cse_p11_t53_g331869 + (-(cse_p11_t1_g331931)))));
            double cse_p60_t50_g335925 = (cse_p60_t23_g335898 + (cse_p11_t1_g331930 * dx * dy));
            double cse_p60_t52_g335927 = (cse_p60_t20_g335895 + (cse_p11_t1_g331930 * dy * (-(dx))));
            double cse_p60_t53_g335928 = (cse_p60_t21_g335896 + cse_p60_t30_g335905);
            double cse_p60_t54_g335929 = (cse_p60_t24_g335899 + cse_p60_t32_g335907);
            double cse_p60_t70_g335945 = (cse_p11_t1_g331938 + cse_p60_t19_g335894);
            double h_x1_x2 = cse_p60_t44_g335919;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p11_t1_g331933 = (cse_p11_t1_g331929 * (-(dy)));
            double cse_p60_t25_g335900 = (cse_p11_t1_g331933 * (-(dy)));
            double cse_p60_t26_g335901 = (cse_p11_t1_g331933 * (-(dz)));
            double cse_p60_t27_g335902 = (cse_p11_t1_g331933 * dx);
            double cse_p60_t28_g335903 = (cse_p11_t1_g331933 * dy);
            double cse_p60_t29_g335904 = (cse_p11_t1_g331933 * dz);
            double cse_p60_t45_g335920 = (cse_p60_t28_g335903 + (cse_p11_t1_g331930 * (cse_p11_t54_g331870 + (-(cse_p11_t1_g331931)))));
            double cse_p60_t51_g335926 = (cse_p60_t27_g335902 + (cse_p11_t1_g331930 * dx * dy));
            double cse_p60_t55_g335930 = (cse_p60_t26_g335901 + cse_p60_t31_g335906);
            double cse_p60_t56_g335931 = (cse_p60_t29_g335904 + cse_p60_t33_g335908);
            double h_x1_y2 = cse_p60_t51_g335926;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p11_t1_g331935 = (cse_p11_t1_g331929 * (-(dz)));
            double cse_p60_t34_g335909 = (cse_p11_t1_g331935 * (-(dz)));
            double cse_p60_t35_g335910 = (cse_p11_t1_g331935 * dx);
            double cse_p60_t36_g335911 = (cse_p11_t1_g331935 * dy);
            double cse_p60_t37_g335912 = (cse_p11_t1_g331935 * dz);
            double cse_p60_t46_g335921 = (cse_p60_t37_g335912 + (cse_p11_t1_g331930 * (cse_p11_t55_g331871 + (-(cse_p11_t1_g331931)))));
            double cse_p60_t57_g335932 = (cse_p60_t32_g335907 + cse_p60_t35_g335910);
            double cse_p60_t58_g335933 = (cse_p60_t33_g335908 + cse_p60_t36_g335911);
            double h_x1_z2 = cse_p60_t57_g335932;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p11_t1_g331936 = (cse_p11_t1_g331930 * (cse_p11_t1_g331931 + (-(cse_p11_t54_g331870))));
            double cse_p60_t65_g335940 = (cse_p11_t1_g331936 + (cse_p11_t1_g331929 * dy * dy));
            double cse_p60_t66_g335941 = (cse_p11_t1_g331936 + cse_p60_t25_g335900);
            double h_y1_y1 = cse_p60_t65_g335940;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t49_g335924;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t50_g335925;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t45_g335920;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t58_g335933;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p11_t1_g331937 = (cse_p11_t1_g331930 * (cse_p11_t1_g331931 + (-(cse_p11_t55_g331871))));
            double cse_p60_t67_g335942 = (cse_p11_t1_g331937 + (cse_p11_t1_g331929 * dz * dz));
            double cse_p60_t68_g335943 = (cse_p11_t1_g331937 + cse_p60_t34_g335909);
            double h_z1_z1 = cse_p60_t67_g335942;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t54_g335929;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t56_g335931;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t46_g335921;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t70_g335945;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t52_g335927;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t53_g335928;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t66_g335941;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t55_g335930;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t68_g335943;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= term.r3_high)) {
              /* !COND-4 */
              double cse_p14_t1_g332166 = (3.00000000000000000e+0 * term.aa_high_i);
              double cse_p14_t2_g332167 = (3.00000000000000000e+0 * term.aa_high_j);
              double cse_p13_t80_g332121 = (cse_p14_t1_g332166 + cse_p14_t2_g332167);
              double cse_p14_t1_g332272 = (cse_p1_t2_invsqrt331812 * dx);
              double cse_p12_t5_g331943 = (term.bb_high_i * cse_p14_t1_g332272 * rr);
              double cse_p14_t1_g332273 = (cse_p1_t2_invsqrt331812 * dy);
              double cse_p12_t6_g331944 = (term.bb_high_i * cse_p14_t1_g332273 * rr);
              double cse_p14_t1_g332274 = (cse_p1_t2_invsqrt331812 * dz);
              double cse_p12_t7_g331945 = (term.bb_high_i * cse_p14_t1_g332274 * rr);
              double cse_p12_t8_g331946 = (term.bb_high_i * rr);
              double cse_p12_t9_g331947 = (term.bb_high_j * rr);
              double cse_p11_t1_g331931 = (rr * rr);
              double cse_p12_t33_g331971 = (cse_p11_t1_g331931 * cse_p13_t80_g332121);
              double cse_p11_t91_g331907 = (term.cc_high_i + term.cc_high_j + cse_p12_t33_g331971 + cse_p12_t8_g331946 + cse_p12_t8_g331946 + cse_p12_t9_g331947 + cse_p12_t9_g331947);
              double cse_p11_t93_g331909 = (-(cse_p12_t5_g331943));
              double cse_p11_t94_g331910 = (-(cse_p12_t6_g331944));
              double cse_p11_t95_g331911 = (-(cse_p12_t7_g331945));
              double energy = (term.dd_high_i + term.dd_high_j + (term.cc_high_i * rr) + (term.cc_high_j * rr) + (cse_p11_t1_g331931 * (term.bb_high_i + term.bb_high_j + (term.aa_high_i * rr) + (term.aa_high_j * rr))));
              energy_added += energy;
              *energy_accumulate += energy;
              double g_x1 = (cse_p11_t91_g331907 * cse_p14_t1_g332272);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t91_g331907 * cse_p14_t1_g332273);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t91_g331907 * cse_p14_t1_g332274);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = (cse_p11_t93_g331909 + cse_p11_t93_g331909 + (-2.00000000000000000e+0 * term.bb_high_j * cse_p14_t1_g332272 * rr) + (cse_p11_t1_g331931 * cse_p14_t1_g332272 * ((-3.00000000000000000e+0 * term.aa_high_i) + (-3.00000000000000000e+0 * term.aa_high_j))) + (-((term.cc_high_i * cse_p14_t1_g332272))) + (-((term.cc_high_j * cse_p14_t1_g332272))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (cse_p11_t94_g331910 + cse_p11_t94_g331910 + (-2.00000000000000000e+0 * term.bb_high_j * cse_p14_t1_g332273 * rr) + (cse_p11_t1_g331931 * cse_p14_t1_g332273 * ((-3.00000000000000000e+0 * term.aa_high_i) + (-3.00000000000000000e+0 * term.aa_high_j))) + (-((term.cc_high_i * cse_p14_t1_g332273))) + (-((term.cc_high_j * cse_p14_t1_g332273))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (cse_p11_t95_g331911 + cse_p11_t95_g331911 + (-2.00000000000000000e+0 * term.bb_high_j * cse_p14_t1_g332274 * rr) + (cse_p11_t1_g331931 * cse_p14_t1_g332274 * ((-3.00000000000000000e+0 * term.aa_high_i) + (-3.00000000000000000e+0 * term.aa_high_j))) + (-((term.cc_high_i * cse_p14_t1_g332274))) + (-((term.cc_high_j * cse_p14_t1_g332274))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double de_dr = (term.cc_high_i + term.cc_high_j + cse_p12_t33_g331971 + (2.00000000000000000e+0 * term.bb_high_i * rr) + (2.00000000000000000e+0 * term.bb_high_j * rr));
              double d2e_dr2 = ((2.00000000000000000e+0 * term.bb_high_i) + (2.00000000000000000e+0 * term.bb_high_j) + (6.00000000000000000e+0 * term.aa_high_i * rr) + (6.00000000000000000e+0 * term.aa_high_j * rr));
              double cse_p565_t2_invr371507 = (1.0 / (rr));
              double cse_p11_t1_g331929 = (cse_p565_t2_invr371507 * cse_p565_t2_invr371507 * d2e_dr2);
              double cse_p565_t1_invr371506 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g331930 = (cse_p565_t1_invr371506 * de_dr);
              double cse_p11_t1_g331938 = (cse_p11_t1_g331930 * (cse_p11_t1_g331931 + (-(cse_p11_t53_g331869))));
              double h_x1_x1 = (cse_p11_t1_g331938 + (cse_p11_t1_g331929 * dx * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double h_x1_y1 = ((cse_p11_t1_g331929 * dx * dy) + (cse_p11_t1_g331930 * dy * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g331934 = (cse_p11_t1_g331930 * dz);
              double h_x1_z1 = ((cse_p11_t1_g331929 * dx * dz) + (cse_p11_t1_g331934 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g331932 = (cse_p11_t1_g331929 * (-(dx)));
              double h_x1_x2 = ((cse_p11_t1_g331930 * (cse_p11_t53_g331869 + (-(cse_p11_t1_g331931)))) + (cse_p11_t1_g331932 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g331933 = (cse_p11_t1_g331929 * (-(dy)));
              double h_x1_y2 = ((cse_p11_t1_g331930 * dx * dy) + (cse_p11_t1_g331933 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g331935 = (cse_p11_t1_g331929 * (-(dz)));
              double h_x1_z2 = ((cse_p11_t1_g331934 * dx) + (cse_p11_t1_g331935 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g331936 = (cse_p11_t1_g331930 * (cse_p11_t1_g331931 + (-(cse_p11_t54_g331870))));
              double h_y1_y1 = (cse_p11_t1_g331936 + (cse_p11_t1_g331929 * dy * dy));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = ((cse_p11_t1_g331929 * dy * dz) + (cse_p11_t1_g331934 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = ((cse_p11_t1_g331930 * dx * dy) + (cse_p11_t1_g331932 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = ((cse_p11_t1_g331930 * (cse_p11_t54_g331870 + (-(cse_p11_t1_g331931)))) + (cse_p11_t1_g331933 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = ((cse_p11_t1_g331934 * dy) + (cse_p11_t1_g331935 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g331937 = (cse_p11_t1_g331930 * (cse_p11_t1_g331931 + (-(cse_p11_t55_g331871))));
              double h_z1_z1 = (cse_p11_t1_g331937 + (cse_p11_t1_g331929 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = ((cse_p11_t1_g331932 * dz) + (cse_p11_t1_g331934 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = ((cse_p11_t1_g331933 * dz) + (cse_p11_t1_g331934 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = ((cse_p11_t1_g331930 * (cse_p11_t55_g331871 + (-(cse_p11_t1_g331931)))) + (cse_p11_t1_g331935 * dz));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p11_t1_g331938 + (cse_p11_t1_g331932 * (-(dx))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = ((cse_p11_t1_g331930 * dy * (-(dx))) + (cse_p11_t1_g331932 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = ((cse_p11_t1_g331932 * (-(dz))) + (cse_p11_t1_g331934 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = (cse_p11_t1_g331936 + (cse_p11_t1_g331933 * (-(dy))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = ((cse_p11_t1_g331933 * (-(dz))) + (cse_p11_t1_g331934 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = (cse_p11_t1_g331937 + (cse_p11_t1_g331935 * (-(dz))));
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
                double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t53_g331869))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t2_invr371507 = (1.0 / (rr));
                double cse_p565_t1_invr371506 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g331931 = (rr * rr);
                double cse_p11_t1_g331930 = (cse_p565_t1_invr371506 * de_dr);
                double cse_p11_t1_g331929 = (cse_p565_t2_invr371507 * cse_p565_t2_invr371507 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g331929 * dx * dy) + (cse_p11_t1_g331930 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double h_x1_z1 = ((cse_p11_t1_g331929 * dx * dz) + (cse_p11_t1_g331930 * dz * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g331934 = (cse_p11_t1_g331930 * dz);
                double h_x1_x2 = ((cse_p11_t1_g331929 * dx * (-(dx))) + (cse_p11_t1_g331930 * (cse_p11_t53_g331869 + (-(cse_p11_t1_g331931)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g331932 = (cse_p11_t1_g331929 * (-(dx)));
                double h_x1_y2 = ((cse_p11_t1_g331929 * dx * (-(dy))) + (cse_p11_t1_g331930 * dx * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g331933 = (cse_p11_t1_g331929 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g331929 * dx * (-(dz))) + (cse_p11_t1_g331934 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g331935 = (cse_p11_t1_g331929 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g331929 * dy * dy) + (cse_p11_t1_g331930 * (cse_p11_t1_g331931 + (-(cse_p11_t54_g331870)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g331936 = (cse_p11_t1_g331930 * (cse_p11_t1_g331931 + (-(cse_p11_t54_g331870))));
                double h_y1_z1 = ((cse_p11_t1_g331929 * dy * dz) + (cse_p11_t1_g331934 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g331930 * dx * dy) + (cse_p11_t1_g331932 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g331930 * (cse_p11_t54_g331870 + (-(cse_p11_t1_g331931)))) + (cse_p11_t1_g331933 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g331934 * dy) + (cse_p11_t1_g331935 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g331929 * dz * dz) + (cse_p11_t1_g331930 * (cse_p11_t1_g331931 + (-(cse_p11_t55_g331871)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g331937 = (cse_p11_t1_g331930 * (cse_p11_t1_g331931 + (-(cse_p11_t55_g331871))));
                double h_z1_x2 = ((cse_p11_t1_g331932 * dz) + (cse_p11_t1_g331934 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g331933 * dz) + (cse_p11_t1_g331934 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g331930 * (cse_p11_t55_g331871 + (-(cse_p11_t1_g331931)))) + (cse_p11_t1_g331935 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g331930 * (cse_p11_t1_g331931 + (-(cse_p11_t53_g331869)))) + (cse_p11_t1_g331932 * (-(dx))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p11_t1_g331938 = (cse_p11_t1_g331930 * (cse_p11_t1_g331931 + (-(cse_p11_t53_g331869))));
                double h_x2_y2 = ((cse_p11_t1_g331930 * dy * (-(dx))) + (cse_p11_t1_g331932 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double h_x2_z2 = ((cse_p11_t1_g331932 * (-(dz))) + (cse_p11_t1_g331934 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g331936 + (cse_p11_t1_g331933 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g331933 * (-(dz))) + (cse_p11_t1_g331934 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g331937 + (cse_p11_t1_g331935 * (-(dz))));
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
