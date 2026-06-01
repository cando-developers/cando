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
      double rr = sqrt(r2);
      if ((rr <= term.r0_low)) {
        /* !COND-1 */
        double cse_p60_t8_g1048 = (-(dx));
        double cse_p60_t9_g1049 = (-(dy));
        double cse_p60_t10_g1050 = (-(dz));
        double cse_p555_t1_invr10556 = (1.0 / (rr));
        double energy = (term.f_contact_i + term.f_contact_j);
        energy_added += energy;
        *energy_accumulate += energy;
        double de_dr = 0.00000000000000000e+0;
        double cse_p11_t1_g23 = (cse_p555_t1_invr10556 * de_dr);
        double cse_p60_t1_g1061 = (cse_p11_t1_g23 * cse_p60_t8_g1048);
        double cse_p60_t3_g1043 = (cse_p11_t1_g23 * cse_p60_t9_g1049);
        double cse_p60_t4_g1044 = (cse_p11_t1_g23 * cse_p60_t10_g1050);
        double cse_p60_t5_g1045 = (cse_p11_t1_g23 * dx);
        double cse_p60_t6_g1046 = (cse_p11_t1_g23 * dy);
        double cse_p60_t7_g1047 = (cse_p11_t1_g23 * dz);
        double g_x1 = cse_p60_t5_g1045;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = cse_p60_t6_g1046;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = cse_p60_t7_g1047;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = cse_p60_t1_g1061;
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = cse_p60_t3_g1043;
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = cse_p60_t4_g1044;
        KernelGradientAcc(term.i3x2, 2, g_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.r1_low)) {
          /* !COND-2 */
          double energy = (term.dd_low_i + term.dd_low_j + (term.cc_low_i * rr) + (term.cc_low_j * rr) + (rr * rr * (term.bb_low_i + term.bb_low_j + (term.aa_low_i * rr) + (term.aa_low_j * rr))));
          energy_added += energy;
          *energy_accumulate += energy;
          double de_dr = (term.cc_low_i + term.cc_low_j + (2.00000000000000000e+0 * term.bb_low_i * rr) + (2.00000000000000000e+0 * term.bb_low_j * rr) + (rr * rr * ((3.00000000000000000e+0 * term.aa_low_i) + (3.00000000000000000e+0 * term.aa_low_j))));
          double cse_p555_t1_invr10556 = (1.0 / (rr));
          double cse_p60_t1_g1041 = (cse_p555_t1_invr10556 * de_dr);
          double g_x1 = (cse_p60_t1_g1041 * dx);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p60_t1_g1041 * dy);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p60_t1_g1041 * dz);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double cse_p12_t1_g53 = (cse_p60_t1_g1041 * (-(dx)));
          double g_x2 = cse_p12_t1_g53;
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p60_t1_g1041 * (-(dy)));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p60_t1_g1041 * (-(dz)));
          KernelGradientAcc(term.i3x2, 2, g_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= term.r2_high)) {
            /* !COND-3 */
            double cse_p11_t10_g20 = (1.0 / (term.lambda2_i));
            double cse_p11_t11_g21 = (1.0 / (term.lambda2_j));
            double cse_p555_t1_invr10556 = (1.0 / (rr));
            double inv_rr = cse_p555_t1_invr10556;
            double inv_rr2 = (inv_rr * inv_rr);
            double di = (rr + (-(term.lj_radius_i)));
            double exp_i = exp((-((cse_p11_t10_g20 * di * di))));
            double fi = (exp_i * inv_rr2 * term.pref_i_j);
            double hi = (-2.00000000000000000e+0 * (inv_rr + (cse_p11_t10_g20 * di)));
            double dfi = (fi * hi);
            double dj = (rr + (-(term.lj_radius_j)));
            double exp_j = exp((-((cse_p11_t11_g21 * dj * dj))));
            double fj = (exp_j * inv_rr2 * term.pref_j_i);
            double hj = (-2.00000000000000000e+0 * (inv_rr + (cse_p11_t11_g21 * dj)));
            double dfj = (fj * hj);
            double energy = (fi + fj);
            energy_added += energy;
            *energy_accumulate += energy;
            double de_dr = (dfi + dfj);
            double cse_p11_t1_g23 = (cse_p555_t1_invr10556 * de_dr);
            double g_x1 = (cse_p11_t1_g23 * dx);
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double g_y1 = (cse_p11_t1_g23 * dy);
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double g_z1 = (cse_p11_t1_g23 * dz);
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double cse_p60_t1_g1061 = (cse_p11_t1_g23 * (-(dx)));
            double g_x2 = cse_p60_t1_g1061;
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (cse_p11_t1_g23 * (-(dy)));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (cse_p11_t1_g23 * (-(dz)));
            KernelGradientAcc(term.i3x2, 2, g_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= term.r3_high)) {
              /* !COND-4 */
              double energy = (term.dd_high_i + term.dd_high_j + (term.cc_high_i * rr) + (term.cc_high_j * rr) + (rr * rr * (term.bb_high_i + term.bb_high_j + (term.aa_high_i * rr) + (term.aa_high_j * rr))));
              energy_added += energy;
              *energy_accumulate += energy;
              double de_dr = (term.cc_high_i + term.cc_high_j + (2.00000000000000000e+0 * term.bb_high_i * rr) + (2.00000000000000000e+0 * term.bb_high_j * rr) + (rr * rr * ((3.00000000000000000e+0 * term.aa_high_i) + (3.00000000000000000e+0 * term.aa_high_j))));
              double cse_p555_t1_invr10556 = (1.0 / (rr));
              double cse_p11_t1_g23 = (cse_p555_t1_invr10556 * de_dr);
              double g_x1 = (cse_p11_t1_g23 * dx);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t1_g23 * dy);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t1_g23 * dz);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p60_t1_g1061 = (cse_p11_t1_g23 * (-(dx)));
              double g_x2 = cse_p60_t1_g1061;
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (cse_p11_t1_g23 * (-(dy)));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (cse_p11_t1_g23 * (-(dz)));
              KernelGradientAcc(term.i3x2, 2, g_z2);
            } else {
              /* !COND-ELSE-4 */
              {
                /* !COND-ELSE-5 */
                double energy = 0.00000000000000000e+0;
                energy_added += energy;
                *energy_accumulate += energy;
                double de_dr = 0.00000000000000000e+0;
                double g_x1 = (de_dr * dx * (1.0 / (rr)));
                KernelGradientAcc(term.i3x1, 0, g_x1);
                double cse_p555_t1_invr10556 = (1.0 / (rr));
                double g_y1 = (cse_p555_t1_invr10556 * de_dr * dy);
                KernelGradientAcc(term.i3x1, 1, g_y1);
                double cse_p11_t1_g23 = (cse_p555_t1_invr10556 * de_dr);
                double g_z1 = (cse_p11_t1_g23 * dz);
                KernelGradientAcc(term.i3x1, 2, g_z1);
                double g_x2 = (cse_p11_t1_g23 * (-(dx)));
                KernelGradientAcc(term.i3x2, 0, g_x2);
                double g_y2 = (cse_p11_t1_g23 * (-(dy)));
                KernelGradientAcc(term.i3x2, 1, g_y2);
                double g_z2 = (cse_p11_t1_g23 * (-(dz)));
                KernelGradientAcc(term.i3x2, 2, g_z2);
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
      double cse_p11_t41_g53 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t42_g54 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t43_g55 = (dz * dz);
      double r2 = (cse_p11_t41_g53 + cse_p11_t42_g54 + cse_p11_t43_g55);
      double rr = sqrt(r2);
      if ((rr <= term.r0_low)) {
        /* !COND-1 */
        double cse_p60_t47_g4244 = (rr * rr);
        double cse_p60_t75_g4272 = (-(cse_p11_t41_g53));
        double cse_p60_t76_g4273 = (-(cse_p11_t42_g54));
        double cse_p60_t77_g4274 = (-(cse_p11_t43_g55));
        double cse_p60_t79_g4276 = (-(dx));
        double cse_p60_t80_g4277 = (-(dy));
        double cse_p60_t81_g4278 = (-(dz));
        double cse_p565_t1_invr41452 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr41453 = (1.0 / (rr));
        double cse_p60_t84_g4281 = (rr * rr * rr);
        double cse_p15_t79_g495 = (-(cse_p60_t47_g4244));
        double cse_p14_t83_g392 = (cse_p11_t41_g53 + cse_p15_t79_g495);
        double cse_p14_t85_g394 = (cse_p11_t42_g54 + cse_p15_t79_g495);
        double cse_p14_t87_g396 = (cse_p11_t43_g55 + cse_p15_t79_g495);
        double cse_p13_t108_g296 = (cse_p60_t47_g4244 + cse_p60_t75_g4272);
        double cse_p13_t109_g297 = (cse_p60_t47_g4244 + cse_p60_t76_g4273);
        double cse_p13_t110_g298 = (cse_p60_t47_g4244 + cse_p60_t77_g4274);
        double energy = (term.f_contact_i + term.f_contact_j);
        energy_added += energy;
        *energy_accumulate += energy;
        double de_dr = 0.00000000000000000e+0;
        double cse_p60_t46_g4243 = (cse_p565_t1_invr41452 * de_dr);
        double cse_p13_t70_g258 = (cse_p14_t83_g392 * cse_p60_t46_g4243);
        double cse_p13_t71_g259 = (cse_p14_t85_g394 * cse_p60_t46_g4243);
        double cse_p13_t72_g260 = (cse_p14_t87_g396 * cse_p60_t46_g4243);
        double cse_p13_t73_g261 = (cse_p60_t46_g4243 * dx * dy);
        double cse_p13_t74_g262 = (cse_p60_t46_g4243 * cse_p60_t79_g4276 * dy);
        double cse_p12_t1_g188 = (cse_p13_t108_g296 * cse_p60_t46_g4243);
        double cse_p12_t25_g128 = (cse_p13_t109_g297 * cse_p60_t46_g4243);
        double cse_p12_t26_g129 = (cse_p13_t110_g298 * cse_p60_t46_g4243);
        double cse_p12_t29_g132 = (cse_p60_t46_g4243 * dz);
        double cse_p13_t60_g248 = (cse_p12_t29_g132 * cse_p60_t79_g4276);
        double cse_p13_t61_g249 = (cse_p12_t29_g132 * cse_p60_t80_g4277);
        double cse_p13_t62_g250 = (cse_p12_t29_g132 * dx);
        double cse_p13_t63_g251 = (cse_p12_t29_g132 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p60_t45_g4242 = (cse_p565_t2_invr41453 * cse_p565_t2_invr41453 * d2e_dr2);
        double cse_p13_t64_g252 = (cse_p60_t45_g4242 * dx * dx);
        double cse_p13_t65_g253 = (cse_p60_t45_g4242 * dx * dy);
        double cse_p13_t66_g254 = (cse_p60_t45_g4242 * dx * dz);
        double cse_p13_t67_g255 = (cse_p60_t45_g4242 * dy * dy);
        double cse_p13_t68_g256 = (cse_p60_t45_g4242 * dy * dz);
        double cse_p13_t69_g257 = (cse_p60_t45_g4242 * dz * dz);
        double cse_p12_t12_g115 = (cse_p60_t45_g4242 * cse_p60_t79_g4276);
        double cse_p13_t45_g233 = (cse_p12_t12_g115 * cse_p60_t79_g4276);
        double cse_p13_t46_g234 = (cse_p12_t12_g115 * cse_p60_t80_g4277);
        double cse_p13_t47_g235 = (cse_p12_t12_g115 * cse_p60_t81_g4278);
        double cse_p13_t48_g236 = (cse_p12_t12_g115 * dx);
        double cse_p13_t49_g237 = (cse_p12_t12_g115 * dy);
        double cse_p13_t50_g238 = (cse_p12_t12_g115 * dz);
        double cse_p12_t13_g116 = (cse_p60_t45_g4242 * cse_p60_t80_g4277);
        double cse_p13_t51_g239 = (cse_p12_t13_g116 * cse_p60_t80_g4277);
        double cse_p13_t52_g240 = (cse_p12_t13_g116 * cse_p60_t81_g4278);
        double cse_p13_t53_g241 = (cse_p12_t13_g116 * dx);
        double cse_p13_t54_g242 = (cse_p12_t13_g116 * dy);
        double cse_p13_t55_g243 = (cse_p12_t13_g116 * dz);
        double cse_p12_t14_g117 = (cse_p60_t45_g4242 * cse_p60_t81_g4278);
        double cse_p13_t56_g244 = (cse_p12_t14_g117 * cse_p60_t81_g4278);
        double cse_p13_t57_g245 = (cse_p12_t14_g117 * dx);
        double cse_p13_t58_g246 = (cse_p12_t14_g117 * dy);
        double cse_p13_t59_g247 = (cse_p12_t14_g117 * dz);
        double cse_p12_t54_g157 = (cse_p12_t1_g188 + cse_p13_t64_g252);
        double cse_p12_t55_g158 = (cse_p13_t65_g253 + cse_p13_t74_g262);
        double cse_p12_t56_g159 = (cse_p13_t60_g248 + cse_p13_t66_g254);
        double cse_p12_t57_g160 = (cse_p13_t61_g249 + cse_p13_t68_g256);
        double cse_p12_t68_g171 = (cse_p12_t25_g128 + cse_p13_t67_g255);
        double cse_p12_t70_g173 = (cse_p12_t26_g129 + cse_p13_t69_g257);
        double cse_p12_t48_g151 = (cse_p13_t48_g236 + cse_p13_t70_g258);
        double cse_p12_t50_g153 = (cse_p12_t1_g188 + cse_p13_t45_g233);
        double cse_p12_t58_g161 = (cse_p13_t49_g237 + cse_p13_t73_g261);
        double cse_p12_t60_g163 = (cse_p13_t46_g234 + cse_p13_t74_g262);
        double cse_p12_t61_g164 = (cse_p13_t47_g235 + cse_p13_t60_g248);
        double cse_p12_t62_g165 = (cse_p13_t50_g238 + cse_p13_t62_g250);
        double cse_p12_t49_g152 = (cse_p13_t54_g242 + cse_p13_t71_g259);
        double cse_p12_t59_g162 = (cse_p13_t53_g241 + cse_p13_t73_g261);
        double cse_p12_t63_g166 = (cse_p13_t52_g240 + cse_p13_t61_g249);
        double cse_p12_t64_g167 = (cse_p13_t55_g243 + cse_p13_t63_g251);
        double cse_p12_t69_g172 = (cse_p12_t25_g128 + cse_p13_t51_g239);
        double cse_p12_t51_g154 = (cse_p13_t57_g245 + cse_p13_t62_g250);
        double cse_p12_t52_g155 = (cse_p13_t58_g246 + cse_p13_t63_g251);
        double cse_p12_t53_g156 = (cse_p13_t59_g247 + cse_p13_t72_g260);
        double cse_p12_t71_g174 = (cse_p12_t26_g129 + cse_p13_t56_g244);
        double cse_p11_t1_g103 = (cse_p565_t2_invr41453 * de_dr);
        double cse_p60_t6_g4203 = (cse_p11_t1_g103 * cse_p60_t79_g4276);
        double cse_p60_t7_g4204 = (cse_p11_t1_g103 * cse_p60_t80_g4277);
        double cse_p60_t8_g4205 = (cse_p11_t1_g103 * cse_p60_t81_g4278);
        double cse_p60_t9_g4206 = (cse_p11_t1_g103 * dx);
        double cse_p60_t10_g4207 = (cse_p11_t1_g103 * dy);
        double cse_p60_t11_g4208 = (cse_p11_t1_g103 * dz);
        double g_x1 = cse_p60_t9_g4206;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = cse_p60_t10_g4207;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = cse_p60_t11_g4208;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = cse_p60_t6_g4203;
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = cse_p60_t7_g4204;
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = cse_p60_t8_g4205;
        KernelGradientAcc(term.i3x2, 2, g_z2);
        double h_x1_x1 = cse_p12_t54_g157;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t55_g158;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t56_g159;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t48_g151;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t59_g162;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t51_g154;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t68_g171;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t57_g160;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t58_g161;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t49_g152;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t52_g155;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t70_g173;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t62_g165;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t64_g167;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t53_g156;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t50_g153;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t60_g163;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t61_g164;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t69_g172;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t63_g166;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t71_g174;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.r1_low)) {
          /* !COND-2 */
          double cse_p11_t1_g96 = (rr * rr);
          double cse_p60_t63_g4260 = (cse_p11_t41_g53 + (-(cse_p11_t1_g96)));
          double cse_p60_t64_g4261 = (cse_p11_t42_g54 + (-(cse_p11_t1_g96)));
          double cse_p60_t65_g4262 = (cse_p11_t43_g55 + (-(cse_p11_t1_g96)));
          double cse_p60_t70_g4267 = (cse_p11_t1_g96 + (-(cse_p11_t41_g53)));
          double cse_p60_t71_g4268 = (cse_p11_t1_g96 + (-(cse_p11_t42_g54)));
          double cse_p60_t72_g4269 = (cse_p11_t1_g96 + (-(cse_p11_t43_g55)));
          double cse_p60_t78_g4275 = (-(cse_p11_t1_g96));
          double energy = (term.dd_low_i + term.dd_low_j + (term.cc_low_i * rr) + (term.cc_low_j * rr) + (cse_p11_t1_g96 * (term.bb_low_i + term.bb_low_j + (term.aa_low_i * rr) + (term.aa_low_j * rr))));
          energy_added += energy;
          *energy_accumulate += energy;
          double de_dr = (term.cc_low_i + term.cc_low_j + (2.00000000000000000e+0 * term.bb_low_i * rr) + (2.00000000000000000e+0 * term.bb_low_j * rr) + (cse_p11_t1_g96 * ((3.00000000000000000e+0 * term.aa_low_i) + (3.00000000000000000e+0 * term.aa_low_j))));
          double d2e_dr2 = ((2.00000000000000000e+0 * term.bb_low_i) + (2.00000000000000000e+0 * term.bb_low_j) + (6.00000000000000000e+0 * term.aa_low_i * rr) + (6.00000000000000000e+0 * term.aa_low_j * rr));
          double cse_p565_t2_invr41453 = (1.0 / (rr));
          double cse_p60_t1_g4198 = (cse_p565_t2_invr41453 * de_dr);
          double g_x1 = (cse_p60_t1_g4198 * dx);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p60_t1_g4198 * dy);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p60_t1_g4198 * dz);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p60_t1_g4198 * (-(dx)));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p60_t1_g4198 * (-(dy)));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p60_t1_g4198 * (-(dz)));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double cse_p11_t1_g94 = (cse_p565_t2_invr41453 * cse_p565_t2_invr41453 * d2e_dr2);
          double cse_p60_t12_g4209 = (cse_p11_t1_g94 * (-(dx)));
          double cse_p60_t13_g4210 = (cse_p11_t1_g94 * (-(dy)));
          double cse_p60_t14_g4211 = (cse_p11_t1_g94 * (-(dz)));
          double cse_p60_t15_g4212 = (cse_p11_t1_g94 * dx * dx);
          double cse_p60_t16_g4213 = (cse_p11_t1_g94 * dx * dy);
          double cse_p60_t17_g4214 = (cse_p11_t1_g94 * dx * dz);
          double cse_p60_t18_g4215 = (cse_p11_t1_g94 * dy * dy);
          double cse_p60_t19_g4216 = (cse_p11_t1_g94 * dy * dz);
          double cse_p60_t20_g4217 = (cse_p11_t1_g94 * dz * dz);
          double cse_p14_t45_g354 = (cse_p60_t12_g4209 * (-(dx)));
          double cse_p14_t46_g355 = (cse_p60_t12_g4209 * (-(dy)));
          double cse_p14_t47_g356 = (cse_p60_t12_g4209 * (-(dz)));
          double cse_p14_t48_g357 = (cse_p60_t12_g4209 * dx);
          double cse_p14_t49_g358 = (cse_p60_t12_g4209 * dy);
          double cse_p14_t50_g359 = (cse_p60_t12_g4209 * dz);
          double cse_p14_t51_g360 = (cse_p60_t13_g4210 * (-(dy)));
          double cse_p14_t52_g361 = (cse_p60_t13_g4210 * (-(dz)));
          double cse_p14_t53_g362 = (cse_p60_t13_g4210 * dx);
          double cse_p14_t54_g363 = (cse_p60_t13_g4210 * dy);
          double cse_p14_t55_g364 = (cse_p60_t13_g4210 * dz);
          double cse_p14_t56_g365 = (cse_p60_t14_g4211 * (-(dz)));
          double cse_p14_t57_g366 = (cse_p60_t14_g4211 * dx);
          double cse_p14_t58_g367 = (cse_p60_t14_g4211 * dy);
          double cse_p14_t59_g368 = (cse_p60_t14_g4211 * dz);
          double cse_p565_t1_invr41452 = (1.0 / ((rr * rr * rr)));
          double cse_p11_t1_g95 = (cse_p565_t1_invr41452 * de_dr);
          double cse_p60_t21_g4218 = (cse_p11_t1_g95 * (cse_p11_t41_g53 + cse_p60_t78_g4275));
          double cse_p60_t22_g4219 = (cse_p11_t1_g95 * (cse_p11_t42_g54 + cse_p60_t78_g4275));
          double cse_p60_t23_g4220 = (cse_p11_t1_g95 * (cse_p11_t43_g55 + cse_p60_t78_g4275));
          double cse_p60_t24_g4221 = (cse_p11_t1_g95 * cse_p60_t70_g4267);
          double cse_p60_t25_g4222 = (cse_p11_t1_g95 * cse_p60_t71_g4268);
          double cse_p60_t26_g4223 = (cse_p11_t1_g95 * cse_p60_t72_g4269);
          double cse_p60_t27_g4224 = (cse_p11_t1_g95 * dx * dy);
          double cse_p60_t28_g4225 = (cse_p11_t1_g95 * dy * (-(dx)));
          double cse_p60_t29_g4226 = (cse_p11_t1_g95 * dz);
          double cse_p60_t53_g4250 = (cse_p60_t16_g4213 + cse_p60_t28_g4225);
          double cse_p14_t60_g369 = (cse_p60_t29_g4226 * (-(dx)));
          double cse_p14_t61_g370 = (cse_p60_t29_g4226 * (-(dy)));
          double cse_p14_t62_g371 = (cse_p60_t29_g4226 * dx);
          double cse_p14_t63_g372 = (cse_p60_t29_g4226 * dy);
          double cse_p13_t106_g294 = (cse_p60_t15_g4212 + cse_p60_t24_g4221);
          double h_x1_x1 = cse_p13_t106_g294;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t53_g4250;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t84_g272 = (cse_p14_t60_g369 + cse_p60_t17_g4214);
          double cse_p13_t85_g273 = (cse_p14_t61_g370 + cse_p60_t19_g4216);
          double h_x1_z1 = cse_p13_t84_g272;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t78_g266 = (cse_p14_t48_g357 + cse_p60_t21_g4218);
          double cse_p13_t86_g274 = (cse_p14_t49_g358 + cse_p60_t27_g4224);
          double cse_p13_t88_g276 = (cse_p14_t46_g355 + cse_p60_t28_g4225);
          double cse_p13_t89_g277 = (cse_p14_t47_g356 + cse_p14_t60_g369);
          double cse_p13_t90_g278 = (cse_p14_t50_g359 + cse_p14_t62_g371);
          double cse_p13_t107_g295 = (cse_p14_t45_g354 + cse_p60_t24_g4221);
          double h_x1_x2 = cse_p13_t78_g266;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t79_g267 = (cse_p14_t54_g363 + cse_p60_t22_g4219);
          double cse_p13_t87_g275 = (cse_p14_t53_g362 + cse_p60_t27_g4224);
          double cse_p13_t91_g279 = (cse_p14_t52_g361 + cse_p14_t61_g370);
          double cse_p13_t92_g280 = (cse_p14_t55_g364 + cse_p14_t63_g372);
          double h_x1_y2 = cse_p13_t87_g275;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t80_g268 = (cse_p14_t57_g366 + cse_p14_t62_g371);
          double cse_p13_t81_g269 = (cse_p14_t58_g367 + cse_p14_t63_g372);
          double cse_p13_t82_g270 = (cse_p14_t59_g368 + cse_p60_t23_g4220);
          double h_x1_z2 = cse_p13_t80_g268;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t99_g287 = (cse_p60_t18_g4215 + cse_p60_t25_g4222);
          double cse_p13_t100_g288 = (cse_p14_t51_g360 + cse_p60_t25_g4222);
          double h_y1_y1 = cse_p13_t99_g287;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t85_g273;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t86_g274;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t79_g267;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t81_g269;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t101_g289 = (cse_p60_t20_g4217 + cse_p60_t26_g4223);
          double cse_p13_t102_g290 = (cse_p14_t56_g365 + cse_p60_t26_g4223);
          double h_z1_z1 = cse_p13_t101_g289;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t90_g278;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t92_g280;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t82_g270;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t107_g295;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t88_g276;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t89_g277;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t100_g288;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t91_g279;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t102_g290;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= term.r2_high)) {
            /* !COND-3 */
            double cse_p565_t2_invr41453 = (1.0 / (rr));
            double inv_rr = cse_p565_t2_invr41453;
            double inv_rr2 = (inv_rr * inv_rr);
            double cse_p11_t1_g13 = (2.00000000000000000e+0 * inv_rr2);
            double di = (rr + (-(term.lj_radius_i)));
            double cse_p1_t2_invr2 = (1.0 / (term.lambda2_i));
            double exp_i = exp((-((cse_p1_t2_invr2 * di * di))));
            double fi = (exp_i * inv_rr2 * term.pref_i_j);
            double hi = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t2_invr2 * di)));
            double dfi = (fi * hi);
            double d2fi = (fi * (cse_p11_t1_g13 + (-2.00000000000000000e+0 * cse_p1_t2_invr2) + (hi * hi)));
            double dj = (rr + (-(term.lj_radius_j)));
            double cse_p1_t3_invr3 = (1.0 / (term.lambda2_j));
            double exp_j = exp((-((cse_p1_t3_invr3 * dj * dj))));
            double fj = (exp_j * inv_rr2 * term.pref_j_i);
            double hj = (-2.00000000000000000e+0 * (inv_rr + (cse_p1_t3_invr3 * dj)));
            double dfj = (fj * hj);
            double d2fj = (fj * (cse_p11_t1_g13 + (-2.00000000000000000e+0 * cse_p1_t3_invr3) + (hj * hj)));
            double energy = (fi + fj);
            energy_added += energy;
            *energy_accumulate += energy;
            double de_dr = (dfi + dfj);
            double d2e_dr2 = (d2fi + d2fj);
            double cse_p11_t1_g103 = (cse_p565_t2_invr41453 * de_dr);
            double g_x1 = (cse_p11_t1_g103 * dx);
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double g_y1 = (cse_p11_t1_g103 * dy);
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double g_z1 = (cse_p11_t1_g103 * dz);
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (cse_p11_t1_g103 * (-(dx)));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (cse_p11_t1_g103 * (-(dy)));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (cse_p11_t1_g103 * (-(dz)));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double cse_p11_t1_g94 = (cse_p565_t2_invr41453 * cse_p565_t2_invr41453 * d2e_dr2);
            double cse_p565_t1_invr41452 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g95 = (cse_p565_t1_invr41452 * de_dr);
            double cse_p11_t1_g96 = (rr * rr);
            double cse_p12_t1_g188 = (cse_p11_t1_g95 * (cse_p11_t1_g96 + (-(cse_p11_t41_g53))));
            double cse_p60_t73_g4270 = (cse_p12_t1_g188 + (cse_p11_t1_g94 * dx * dx));
            double h_x1_x1 = cse_p60_t73_g4270;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = ((cse_p11_t1_g94 * dx * dy) + (cse_p11_t1_g95 * dy * (-(dx))));
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p11_t1_g99 = (cse_p11_t1_g95 * dz);
            double cse_p60_t41_g4238 = (cse_p11_t1_g99 * (-(dx)));
            double cse_p60_t42_g4239 = (cse_p11_t1_g99 * (-(dy)));
            double cse_p60_t43_g4240 = (cse_p11_t1_g99 * dx);
            double cse_p60_t44_g4241 = (cse_p11_t1_g99 * dy);
            double cse_p60_t54_g4251 = (cse_p60_t41_g4238 + (cse_p11_t1_g94 * dx * dz));
            double cse_p60_t55_g4252 = (cse_p60_t42_g4239 + (cse_p11_t1_g94 * dy * dz));
            double h_x1_z1 = cse_p60_t54_g4251;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p11_t1_g97 = (cse_p11_t1_g94 * (-(dx)));
            double cse_p60_t30_g4227 = (cse_p11_t1_g97 * (-(dx)));
            double cse_p60_t31_g4228 = (cse_p11_t1_g97 * (-(dy)));
            double cse_p60_t32_g4229 = (cse_p11_t1_g97 * (-(dz)));
            double cse_p60_t33_g4230 = (cse_p11_t1_g97 * dx);
            double cse_p60_t34_g4231 = (cse_p11_t1_g97 * dy);
            double cse_p60_t35_g4232 = (cse_p11_t1_g97 * dz);
            double cse_p60_t48_g4245 = (cse_p60_t33_g4230 + (cse_p11_t1_g95 * (cse_p11_t41_g53 + (-(cse_p11_t1_g96)))));
            double cse_p60_t56_g4253 = (cse_p60_t34_g4231 + (cse_p11_t1_g95 * dx * dy));
            double cse_p60_t58_g4255 = (cse_p60_t31_g4228 + (cse_p11_t1_g95 * dy * (-(dx))));
            double cse_p60_t59_g4256 = (cse_p60_t32_g4229 + cse_p60_t41_g4238);
            double cse_p60_t60_g4257 = (cse_p60_t35_g4232 + cse_p60_t43_g4240);
            double cse_p60_t74_g4271 = (cse_p12_t1_g188 + cse_p60_t30_g4227);
            double h_x1_x2 = cse_p60_t48_g4245;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p11_t1_g98 = (cse_p11_t1_g94 * (-(dy)));
            double cse_p60_t36_g4233 = (cse_p11_t1_g98 * (-(dy)));
            double cse_p60_t37_g4234 = (cse_p11_t1_g98 * (-(dz)));
            double cse_p60_t38_g4235 = (cse_p11_t1_g98 * dx);
            double cse_p60_t39_g4236 = (cse_p11_t1_g98 * dy);
            double cse_p60_t40_g4237 = (cse_p11_t1_g98 * dz);
            double cse_p60_t49_g4246 = (cse_p60_t39_g4236 + (cse_p11_t1_g95 * (cse_p11_t42_g54 + (-(cse_p11_t1_g96)))));
            double cse_p60_t57_g4254 = (cse_p60_t38_g4235 + (cse_p11_t1_g95 * dx * dy));
            double cse_p60_t61_g4258 = (cse_p60_t37_g4234 + cse_p60_t42_g4239);
            double cse_p60_t62_g4259 = (cse_p60_t40_g4237 + cse_p60_t44_g4241);
            double h_x1_y2 = cse_p60_t57_g4254;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p11_t1_g100 = (cse_p11_t1_g94 * (-(dz)));
            double cse_p60_t2_g4199 = (cse_p11_t1_g100 * (-(dz)));
            double cse_p60_t3_g4200 = (cse_p11_t1_g100 * dx);
            double cse_p60_t4_g4201 = (cse_p11_t1_g100 * dy);
            double cse_p60_t5_g4202 = (cse_p11_t1_g100 * dz);
            double cse_p60_t50_g4247 = (cse_p60_t3_g4200 + cse_p60_t43_g4240);
            double cse_p60_t51_g4248 = (cse_p60_t44_g4241 + cse_p60_t4_g4201);
            double cse_p60_t52_g4249 = (cse_p60_t5_g4202 + (cse_p11_t1_g95 * (cse_p11_t43_g55 + (-(cse_p11_t1_g96)))));
            double h_x1_z2 = cse_p60_t50_g4247;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p11_t1_g101 = (cse_p11_t1_g95 * (cse_p11_t1_g96 + (-(cse_p11_t42_g54))));
            double cse_p60_t66_g4263 = (cse_p11_t1_g101 + (cse_p11_t1_g94 * dy * dy));
            double cse_p60_t67_g4264 = (cse_p11_t1_g101 + cse_p60_t36_g4233);
            double h_y1_y1 = cse_p60_t66_g4263;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t55_g4252;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t56_g4253;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t49_g4246;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t51_g4248;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p11_t1_g102 = (cse_p11_t1_g95 * (cse_p11_t1_g96 + (-(cse_p11_t43_g55))));
            double cse_p60_t68_g4265 = (cse_p11_t1_g102 + (cse_p11_t1_g94 * dz * dz));
            double cse_p60_t69_g4266 = (cse_p11_t1_g102 + cse_p60_t2_g4199);
            double h_z1_z1 = cse_p60_t68_g4265;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t60_g4257;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t62_g4259;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t52_g4249;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t74_g4271;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t58_g4255;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t59_g4256;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t67_g4264;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t61_g4258;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t69_g4266;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= term.r3_high)) {
              /* !COND-4 */
              double cse_p11_t1_g96 = (rr * rr);
              double energy = (term.dd_high_i + term.dd_high_j + (term.cc_high_i * rr) + (term.cc_high_j * rr) + (cse_p11_t1_g96 * (term.bb_high_i + term.bb_high_j + (term.aa_high_i * rr) + (term.aa_high_j * rr))));
              energy_added += energy;
              *energy_accumulate += energy;
              double de_dr = (term.cc_high_i + term.cc_high_j + (2.00000000000000000e+0 * term.bb_high_i * rr) + (2.00000000000000000e+0 * term.bb_high_j * rr) + (cse_p11_t1_g96 * ((3.00000000000000000e+0 * term.aa_high_i) + (3.00000000000000000e+0 * term.aa_high_j))));
              double d2e_dr2 = ((2.00000000000000000e+0 * term.bb_high_i) + (2.00000000000000000e+0 * term.bb_high_j) + (6.00000000000000000e+0 * term.aa_high_i * rr) + (6.00000000000000000e+0 * term.aa_high_j * rr));
              double cse_p565_t2_invr41453 = (1.0 / (rr));
              double cse_p11_t1_g103 = (cse_p565_t2_invr41453 * de_dr);
              double g_x1 = (cse_p11_t1_g103 * dx);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t1_g103 * dy);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t1_g103 * dz);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = (cse_p11_t1_g103 * (-(dx)));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (cse_p11_t1_g103 * (-(dy)));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (cse_p11_t1_g103 * (-(dz)));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double cse_p11_t1_g94 = (cse_p565_t2_invr41453 * cse_p565_t2_invr41453 * d2e_dr2);
              double cse_p565_t1_invr41452 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g95 = (cse_p565_t1_invr41452 * de_dr);
              double cse_p12_t1_g188 = (cse_p11_t1_g95 * (cse_p11_t1_g96 + (-(cse_p11_t41_g53))));
              double h_x1_x1 = (cse_p12_t1_g188 + (cse_p11_t1_g94 * dx * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double h_x1_y1 = ((cse_p11_t1_g94 * dx * dy) + (cse_p11_t1_g95 * dy * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g99 = (cse_p11_t1_g95 * dz);
              double h_x1_z1 = ((cse_p11_t1_g94 * dx * dz) + (cse_p11_t1_g99 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g97 = (cse_p11_t1_g94 * (-(dx)));
              double h_x1_x2 = ((cse_p11_t1_g95 * (cse_p11_t41_g53 + (-(cse_p11_t1_g96)))) + (cse_p11_t1_g97 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g98 = (cse_p11_t1_g94 * (-(dy)));
              double h_x1_y2 = ((cse_p11_t1_g95 * dx * dy) + (cse_p11_t1_g98 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g100 = (cse_p11_t1_g94 * (-(dz)));
              double h_x1_z2 = ((cse_p11_t1_g100 * dx) + (cse_p11_t1_g99 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g101 = (cse_p11_t1_g95 * (cse_p11_t1_g96 + (-(cse_p11_t42_g54))));
              double h_y1_y1 = (cse_p11_t1_g101 + (cse_p11_t1_g94 * dy * dy));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = ((cse_p11_t1_g94 * dy * dz) + (cse_p11_t1_g99 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = ((cse_p11_t1_g95 * dx * dy) + (cse_p11_t1_g97 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = ((cse_p11_t1_g95 * (cse_p11_t42_g54 + (-(cse_p11_t1_g96)))) + (cse_p11_t1_g98 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = ((cse_p11_t1_g100 * dy) + (cse_p11_t1_g99 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g102 = (cse_p11_t1_g95 * (cse_p11_t1_g96 + (-(cse_p11_t43_g55))));
              double h_z1_z1 = (cse_p11_t1_g102 + (cse_p11_t1_g94 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = ((cse_p11_t1_g97 * dz) + (cse_p11_t1_g99 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = ((cse_p11_t1_g98 * dz) + (cse_p11_t1_g99 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = ((cse_p11_t1_g100 * dz) + (cse_p11_t1_g95 * (cse_p11_t43_g55 + (-(cse_p11_t1_g96)))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p12_t1_g188 + (cse_p11_t1_g97 * (-(dx))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = ((cse_p11_t1_g95 * dy * (-(dx))) + (cse_p11_t1_g97 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = ((cse_p11_t1_g97 * (-(dz))) + (cse_p11_t1_g99 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = (cse_p11_t1_g101 + (cse_p11_t1_g98 * (-(dy))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = ((cse_p11_t1_g98 * (-(dz))) + (cse_p11_t1_g99 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = (cse_p11_t1_g102 + (cse_p11_t1_g100 * (-(dz))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
            } else {
              /* !COND-ELSE-4 */
              {
                /* !COND-ELSE-5 */
                double energy = 0.00000000000000000e+0;
                energy_added += energy;
                *energy_accumulate += energy;
                double de_dr = 0.00000000000000000e+0;
                double d2e_dr2 = 0.00000000000000000e+0;
                double g_x1 = (de_dr * dx * (1.0 / (rr)));
                KernelGradientAcc(term.i3x1, 0, g_x1);
                double cse_p565_t2_invr41453 = (1.0 / (rr));
                double g_y1 = (cse_p565_t2_invr41453 * de_dr * dy);
                KernelGradientAcc(term.i3x1, 1, g_y1);
                double cse_p11_t1_g103 = (cse_p565_t2_invr41453 * de_dr);
                double g_z1 = (cse_p11_t1_g103 * dz);
                KernelGradientAcc(term.i3x1, 2, g_z1);
                double g_x2 = (cse_p11_t1_g103 * (-(dx)));
                KernelGradientAcc(term.i3x2, 0, g_x2);
                double g_y2 = (cse_p11_t1_g103 * (-(dy)));
                KernelGradientAcc(term.i3x2, 1, g_y2);
                double g_z2 = (cse_p11_t1_g103 * (-(dz)));
                KernelGradientAcc(term.i3x2, 2, g_z2);
                double h_x1_x1 = ((cse_p565_t2_invr41453 * cse_p565_t2_invr41453 * d2e_dr2 * dx * dx) + (de_dr * ((rr * rr) + (-(cse_p11_t41_g53))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t1_invr41452 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g96 = (rr * rr);
                double cse_p11_t1_g95 = (cse_p565_t1_invr41452 * de_dr);
                double cse_p11_t1_g94 = (cse_p565_t2_invr41453 * cse_p565_t2_invr41453 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g94 * dx * dy) + (cse_p11_t1_g95 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double h_x1_z1 = ((cse_p11_t1_g94 * dx * dz) + (cse_p11_t1_g95 * dz * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g99 = (cse_p11_t1_g95 * dz);
                double h_x1_x2 = ((cse_p11_t1_g94 * dx * (-(dx))) + (cse_p11_t1_g95 * (cse_p11_t41_g53 + (-(cse_p11_t1_g96)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g97 = (cse_p11_t1_g94 * (-(dx)));
                double h_x1_y2 = ((cse_p11_t1_g94 * dx * (-(dy))) + (cse_p11_t1_g95 * dx * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g98 = (cse_p11_t1_g94 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g94 * dx * (-(dz))) + (cse_p11_t1_g99 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g100 = (cse_p11_t1_g94 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g94 * dy * dy) + (cse_p11_t1_g95 * (cse_p11_t1_g96 + (-(cse_p11_t42_g54)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g101 = (cse_p11_t1_g95 * (cse_p11_t1_g96 + (-(cse_p11_t42_g54))));
                double h_y1_z1 = ((cse_p11_t1_g94 * dy * dz) + (cse_p11_t1_g99 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g95 * dx * dy) + (cse_p11_t1_g97 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g95 * (cse_p11_t42_g54 + (-(cse_p11_t1_g96)))) + (cse_p11_t1_g98 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g100 * dy) + (cse_p11_t1_g99 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g94 * dz * dz) + (cse_p11_t1_g95 * (cse_p11_t1_g96 + (-(cse_p11_t43_g55)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g102 = (cse_p11_t1_g95 * (cse_p11_t1_g96 + (-(cse_p11_t43_g55))));
                double h_z1_x2 = ((cse_p11_t1_g97 * dz) + (cse_p11_t1_g99 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g98 * dz) + (cse_p11_t1_g99 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g100 * dz) + (cse_p11_t1_g95 * (cse_p11_t43_g55 + (-(cse_p11_t1_g96)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g95 * (cse_p11_t1_g96 + (-(cse_p11_t41_g53)))) + (cse_p11_t1_g97 * (-(dx))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p12_t1_g188 = (cse_p11_t1_g95 * (cse_p11_t1_g96 + (-(cse_p11_t41_g53))));
                double h_x2_y2 = ((cse_p11_t1_g95 * dy * (-(dx))) + (cse_p11_t1_g97 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double h_x2_z2 = ((cse_p11_t1_g97 * (-(dz))) + (cse_p11_t1_g99 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g101 + (cse_p11_t1_g98 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g98 * (-(dz))) + (cse_p11_t1_g99 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g102 + (cse_p11_t1_g100 * (-(dz))));
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
