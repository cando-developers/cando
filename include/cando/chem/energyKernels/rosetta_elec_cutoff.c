#include "rosetta_elec_cutoff.h"

template <typename HESSIAN>
struct Rosetta_Elec_Cutoff {
  static constexpr size_t PositionSize = 6;
  static std::string description() { return "mathkernel-rosetta_elec_cutoff"; };
double gradient(const rosetta_elec_parameters& params, const rosetta_elec_term& term, double* position, double* energy_accumulate, double* force) {
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
      if ((rr <= params.rmin)) {
        /* !COND-1 */
        double cse_p60_t8_g1051 = (-(dx));
        double cse_p60_t9_g1052 = (-(dy));
        double cse_p60_t10_g1053 = (-(dz));
        double cse_p556_t1_invr10560 = (1.0 / (rr));
        double energy = (params.elec_weight * term.e_rmin);
        energy_added += energy;
        *energy_accumulate += energy;
        double de_dr = 0.00000000000000000e+0;
        double cse_p11_t1_g25 = (cse_p556_t1_invr10560 * de_dr);
        double cse_p60_t1_g1064 = (cse_p11_t1_g25 * cse_p60_t8_g1051);
        double cse_p60_t3_g1046 = (cse_p11_t1_g25 * cse_p60_t9_g1052);
        double cse_p60_t4_g1047 = (cse_p11_t1_g25 * cse_p60_t10_g1053);
        double cse_p60_t5_g1048 = (cse_p11_t1_g25 * dx);
        double cse_p60_t6_g1049 = (cse_p11_t1_g25 * dy);
        double cse_p60_t7_g1050 = (cse_p11_t1_g25 * dz);
        double g_x1 = cse_p60_t5_g1048;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = cse_p60_t6_g1049;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = cse_p60_t7_g1050;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = cse_p60_t1_g1064;
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = cse_p60_t3_g1046;
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = cse_p60_t4_g1047;
        KernelGradientAcc(term.i3x2, 2, g_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= params.rlow)) {
          /* !COND-2 */
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (rr * rr * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
          double de_dr = (params.elec_weight * (term.cc_low + (2.00000000000000000e+0 * term.bb_low * rr) + (3.00000000000000000e+0 * term.aa_low * rr * rr)));
          double cse_p556_t1_invr10560 = (1.0 / (rr));
          double cse_p60_t1_g1044 = (cse_p556_t1_invr10560 * de_dr);
          double g_x1 = (cse_p60_t1_g1044 * dx);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p60_t1_g1044 * dy);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p60_t1_g1044 * dz);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double cse_p12_t1_g56 = (cse_p60_t1_g1044 * (-(dx)));
          double g_x2 = cse_p12_t1_g56;
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p60_t1_g1044 * (-(dy)));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p60_t1_g1044 * (-(dz)));
          KernelGradientAcc(term.i3x2, 2, g_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t11_g45 = (-(params.eps_solvent));
            double cse_p11_t8_g18 = (cse_p12_t11_g45 + params.eps_core);
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p11_t1_g11 = (5.00000000000000000e-1 * rx * rx);
            double cse_p11_t13_g23 = (-(rx));
            double ex = exp(cse_p11_t13_g23);
            double gx = (ex * (1.00000000000000000e+0 + cse_p11_t1_g11 + rx));
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p11_t13_g23 + cse_p11_t1_g11));
            double ep = (params.eps_solvent + (cse_p11_t8_g18 * gx));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t8_g18 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t8_g18 * gx_pp);
            double ee = (term.kqq * (1.0 / (ep)));
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
            double de_dr = (params.elec_weight * ep_p * (1.0 / ((ep * ep))) * (-(term.kqq)));
            double cse_p556_t1_invr10560 = (1.0 / (rr));
            double cse_p11_t1_g25 = (cse_p556_t1_invr10560 * de_dr);
            double g_x1 = (cse_p11_t1_g25 * dx);
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double g_y1 = (cse_p11_t1_g25 * dy);
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double g_z1 = (cse_p11_t1_g25 * dz);
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double cse_p60_t1_g1064 = (cse_p11_t1_g25 * (-(dx)));
            double g_x2 = cse_p60_t1_g1064;
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (cse_p11_t1_g25 * (-(dy)));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (cse_p11_t1_g25 * (-(dz)));
            KernelGradientAcc(term.i3x2, 2, g_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (rr * rr * (term.bb_high + (term.aa_high * rr)))));
              energy_added += energy;
              *energy_accumulate += energy;
              double de_dr = (params.elec_weight * (term.cc_high + (2.00000000000000000e+0 * term.bb_high * rr) + (3.00000000000000000e+0 * term.aa_high * rr * rr)));
              double cse_p556_t1_invr10560 = (1.0 / (rr));
              double cse_p11_t1_g25 = (cse_p556_t1_invr10560 * de_dr);
              double g_x1 = (cse_p11_t1_g25 * dx);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t1_g25 * dy);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t1_g25 * dz);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p60_t1_g1064 = (cse_p11_t1_g25 * (-(dx)));
              double g_x2 = cse_p60_t1_g1064;
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (cse_p11_t1_g25 * (-(dy)));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (cse_p11_t1_g25 * (-(dz)));
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
                double cse_p556_t1_invr10560 = (1.0 / (rr));
                double g_y1 = (cse_p556_t1_invr10560 * de_dr * dy);
                KernelGradientAcc(term.i3x1, 1, g_y1);
                double cse_p11_t1_g25 = (cse_p556_t1_invr10560 * de_dr);
                double g_z1 = (cse_p11_t1_g25 * dz);
                KernelGradientAcc(term.i3x1, 2, g_z1);
                double g_x2 = (cse_p11_t1_g25 * (-(dx)));
                KernelGradientAcc(term.i3x2, 0, g_x2);
                double g_y2 = (cse_p11_t1_g25 * (-(dy)));
                KernelGradientAcc(term.i3x2, 1, g_y2);
                double g_z2 = (cse_p11_t1_g25 * (-(dz)));
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
void gradient_fd(const rosetta_elec_parameters& params, const rosetta_elec_term& term, double* position, double* energy_accumulate, double* force)
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

double energy(const rosetta_elec_parameters& params, const rosetta_elec_term& term, double* position, double* energy_accumulate) {
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
      if ((rr <= params.rmin)) {
        /* !COND-1 */
        double energy = (params.elec_weight * term.e_rmin);
        energy_added += energy;
        *energy_accumulate += energy;
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= params.rlow)) {
          /* !COND-2 */
          double cse_p528_t1_g5513 = (rr * rr);
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p528_t1_g5513 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t2_g16 = (-(params.eps_solvent));
            double cse_p11_t2_g2 = (cse_p12_t2_g16 + params.eps_core);
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p11_t1_g1 = (5.00000000000000000e-1 * rx * rx);
            double cse_p11_t4_g4 = (-(rx));
            double ex = exp(cse_p11_t4_g4);
            double gx = (ex * (1.00000000000000000e+0 + cse_p11_t1_g1 + rx));
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p11_t1_g1 + cse_p11_t4_g4));
            double ep = (params.eps_solvent + (cse_p11_t2_g2 * gx));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t2_g2 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t2_g2 * gx_pp);
            double ee = (term.kqq * (1.0 / (ep)));
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p528_t1_g5513 = (rr * rr);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p528_t1_g5513 * (term.bb_high + (term.aa_high * rr)))));
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
void energy_fd(const rosetta_elec_parameters& params, const rosetta_elec_term& term, double* position, double* energy_accumulate)
{
  energy(params, term, position, energy_accumulate);
}

double hessian(const rosetta_elec_parameters& params, const rosetta_elec_term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
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
      double cse_p11_t41_g51 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t42_g52 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t43_g53 = (dz * dz);
      double r2 = (cse_p11_t41_g51 + cse_p11_t42_g52 + cse_p11_t43_g53);
      double rr = sqrt(r2);
      if ((rr <= params.rmin)) {
        /* !COND-1 */
        double cse_p60_t47_g4250 = (rr * rr);
        double cse_p60_t75_g4278 = (-(cse_p11_t41_g51));
        double cse_p60_t76_g4279 = (-(cse_p11_t42_g52));
        double cse_p60_t77_g4280 = (-(cse_p11_t43_g53));
        double cse_p60_t79_g4282 = (-(dx));
        double cse_p60_t80_g4283 = (-(dy));
        double cse_p60_t81_g4284 = (-(dz));
        double cse_p565_t1_invr41458 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr41459 = (1.0 / (rr));
        double cse_p60_t84_g4287 = (rr * rr * rr);
        double cse_p15_t79_g501 = (-(cse_p60_t47_g4250));
        double cse_p14_t83_g398 = (cse_p11_t41_g51 + cse_p15_t79_g501);
        double cse_p14_t85_g400 = (cse_p11_t42_g52 + cse_p15_t79_g501);
        double cse_p14_t87_g402 = (cse_p11_t43_g53 + cse_p15_t79_g501);
        double cse_p13_t108_g302 = (cse_p60_t47_g4250 + cse_p60_t75_g4278);
        double cse_p13_t109_g303 = (cse_p60_t47_g4250 + cse_p60_t76_g4279);
        double cse_p13_t110_g304 = (cse_p60_t47_g4250 + cse_p60_t77_g4280);
        double energy = (params.elec_weight * term.e_rmin);
        energy_added += energy;
        *energy_accumulate += energy;
        double de_dr = 0.00000000000000000e+0;
        double cse_p60_t46_g4249 = (cse_p565_t1_invr41458 * de_dr);
        double cse_p13_t70_g264 = (cse_p14_t83_g398 * cse_p60_t46_g4249);
        double cse_p13_t71_g265 = (cse_p14_t85_g400 * cse_p60_t46_g4249);
        double cse_p13_t72_g266 = (cse_p14_t87_g402 * cse_p60_t46_g4249);
        double cse_p13_t73_g267 = (cse_p60_t46_g4249 * dx * dy);
        double cse_p13_t74_g268 = (cse_p60_t46_g4249 * cse_p60_t79_g4282 * dy);
        double cse_p12_t1_g194 = (cse_p13_t108_g302 * cse_p60_t46_g4249);
        double cse_p12_t40_g147 = (cse_p13_t109_g303 * cse_p60_t46_g4249);
        double cse_p12_t41_g148 = (cse_p13_t110_g304 * cse_p60_t46_g4249);
        double cse_p12_t44_g151 = (cse_p60_t46_g4249 * dz);
        double cse_p13_t60_g254 = (cse_p12_t44_g151 * cse_p60_t79_g4282);
        double cse_p13_t61_g255 = (cse_p12_t44_g151 * cse_p60_t80_g4283);
        double cse_p13_t62_g256 = (cse_p12_t44_g151 * dx);
        double cse_p13_t63_g257 = (cse_p12_t44_g151 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p60_t45_g4248 = (cse_p565_t2_invr41459 * cse_p565_t2_invr41459 * d2e_dr2);
        double cse_p13_t64_g258 = (cse_p60_t45_g4248 * dx * dx);
        double cse_p13_t65_g259 = (cse_p60_t45_g4248 * dx * dy);
        double cse_p13_t66_g260 = (cse_p60_t45_g4248 * dx * dz);
        double cse_p13_t67_g261 = (cse_p60_t45_g4248 * dy * dy);
        double cse_p13_t68_g262 = (cse_p60_t45_g4248 * dy * dz);
        double cse_p13_t69_g263 = (cse_p60_t45_g4248 * dz * dz);
        double cse_p12_t27_g134 = (cse_p60_t45_g4248 * cse_p60_t79_g4282);
        double cse_p13_t45_g239 = (cse_p12_t27_g134 * cse_p60_t79_g4282);
        double cse_p13_t46_g240 = (cse_p12_t27_g134 * cse_p60_t80_g4283);
        double cse_p13_t47_g241 = (cse_p12_t27_g134 * cse_p60_t81_g4284);
        double cse_p13_t48_g242 = (cse_p12_t27_g134 * dx);
        double cse_p13_t49_g243 = (cse_p12_t27_g134 * dy);
        double cse_p13_t50_g244 = (cse_p12_t27_g134 * dz);
        double cse_p12_t28_g135 = (cse_p60_t45_g4248 * cse_p60_t80_g4283);
        double cse_p13_t51_g245 = (cse_p12_t28_g135 * cse_p60_t80_g4283);
        double cse_p13_t52_g246 = (cse_p12_t28_g135 * cse_p60_t81_g4284);
        double cse_p13_t53_g247 = (cse_p12_t28_g135 * dx);
        double cse_p13_t54_g248 = (cse_p12_t28_g135 * dy);
        double cse_p13_t55_g249 = (cse_p12_t28_g135 * dz);
        double cse_p12_t29_g136 = (cse_p60_t45_g4248 * cse_p60_t81_g4284);
        double cse_p13_t56_g250 = (cse_p12_t29_g136 * cse_p60_t81_g4284);
        double cse_p13_t57_g251 = (cse_p12_t29_g136 * dx);
        double cse_p13_t58_g252 = (cse_p12_t29_g136 * dy);
        double cse_p13_t59_g253 = (cse_p12_t29_g136 * dz);
        double cse_p12_t59_g166 = (cse_p13_t60_g254 + cse_p13_t66_g260);
        double cse_p12_t60_g167 = (cse_p13_t61_g255 + cse_p13_t68_g262);
        double cse_p12_t64_g171 = (cse_p12_t1_g194 + cse_p13_t64_g258);
        double cse_p12_t65_g172 = (cse_p13_t65_g259 + cse_p13_t74_g268);
        double cse_p12_t72_g179 = (cse_p12_t40_g147 + cse_p13_t67_g261);
        double cse_p12_t74_g181 = (cse_p12_t41_g148 + cse_p13_t69_g263);
        double cse_p12_t49_g156 = (cse_p12_t1_g194 + cse_p13_t45_g239);
        double cse_p12_t50_g157 = (cse_p13_t46_g240 + cse_p13_t74_g268);
        double cse_p12_t51_g158 = (cse_p13_t47_g241 + cse_p13_t60_g254);
        double cse_p12_t52_g159 = (cse_p13_t48_g242 + cse_p13_t70_g264);
        double cse_p12_t53_g160 = (cse_p13_t49_g243 + cse_p13_t73_g267);
        double cse_p12_t54_g161 = (cse_p13_t50_g244 + cse_p13_t62_g256);
        double cse_p12_t55_g162 = (cse_p13_t52_g246 + cse_p13_t61_g255);
        double cse_p12_t56_g163 = (cse_p13_t53_g247 + cse_p13_t73_g267);
        double cse_p12_t57_g164 = (cse_p13_t54_g248 + cse_p13_t71_g265);
        double cse_p12_t58_g165 = (cse_p13_t55_g249 + cse_p13_t63_g257);
        double cse_p12_t73_g180 = (cse_p12_t40_g147 + cse_p13_t51_g245);
        double cse_p12_t61_g168 = (cse_p13_t57_g251 + cse_p13_t62_g256);
        double cse_p12_t62_g169 = (cse_p13_t58_g252 + cse_p13_t63_g257);
        double cse_p12_t63_g170 = (cse_p13_t59_g253 + cse_p13_t72_g266);
        double cse_p12_t75_g182 = (cse_p12_t41_g148 + cse_p13_t56_g250);
        double cse_p11_t1_g107 = (cse_p565_t2_invr41459 * de_dr);
        double cse_p60_t21_g4224 = (cse_p11_t1_g107 * cse_p60_t79_g4282);
        double cse_p60_t22_g4225 = (cse_p11_t1_g107 * cse_p60_t80_g4283);
        double cse_p60_t23_g4226 = (cse_p11_t1_g107 * cse_p60_t81_g4284);
        double cse_p60_t24_g4227 = (cse_p11_t1_g107 * dx);
        double cse_p60_t25_g4228 = (cse_p11_t1_g107 * dy);
        double cse_p60_t26_g4229 = (cse_p11_t1_g107 * dz);
        double g_x1 = cse_p60_t24_g4227;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = cse_p60_t25_g4228;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = cse_p60_t26_g4229;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = cse_p60_t21_g4224;
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = cse_p60_t22_g4225;
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = cse_p60_t23_g4226;
        KernelGradientAcc(term.i3x2, 2, g_z2);
        double h_x1_x1 = cse_p12_t64_g171;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t65_g172;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t59_g166;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t52_g159;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t56_g163;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t61_g168;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t72_g179;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t60_g167;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t53_g160;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t57_g164;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t62_g169;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t74_g181;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t54_g161;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t58_g165;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t63_g170;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t49_g156;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t50_g157;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t51_g158;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t73_g180;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t55_g162;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t75_g182;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= params.rlow)) {
          /* !COND-2 */
          double cse_p11_t1_g100 = (rr * rr);
          double cse_p60_t63_g4266 = (cse_p11_t41_g51 + (-(cse_p11_t1_g100)));
          double cse_p60_t64_g4267 = (cse_p11_t42_g52 + (-(cse_p11_t1_g100)));
          double cse_p60_t65_g4268 = (cse_p11_t43_g53 + (-(cse_p11_t1_g100)));
          double cse_p60_t66_g4269 = (cse_p11_t1_g100 + (-(cse_p11_t41_g51)));
          double cse_p60_t67_g4270 = (cse_p11_t1_g100 + (-(cse_p11_t42_g52)));
          double cse_p60_t68_g4271 = (cse_p11_t1_g100 + (-(cse_p11_t43_g53)));
          double cse_p60_t78_g4281 = (-(cse_p11_t1_g100));
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p11_t1_g100 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
          double de_dr = (params.elec_weight * (term.cc_low + (2.00000000000000000e+0 * term.bb_low * rr) + (3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g100)));
          double d2e_dr2 = (params.elec_weight * ((2.00000000000000000e+0 * term.bb_low) + (6.00000000000000000e+0 * term.aa_low * rr)));
          double cse_p565_t2_invr41459 = (1.0 / (rr));
          double cse_p60_t1_g4204 = (cse_p565_t2_invr41459 * de_dr);
          double g_x1 = (cse_p60_t1_g4204 * dx);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p60_t1_g4204 * dy);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p60_t1_g4204 * dz);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p60_t1_g4204 * (-(dx)));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p60_t1_g4204 * (-(dy)));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p60_t1_g4204 * (-(dz)));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double cse_p11_t1_g98 = (cse_p565_t2_invr41459 * cse_p565_t2_invr41459 * d2e_dr2);
          double cse_p60_t27_g4230 = (cse_p11_t1_g98 * (-(dx)));
          double cse_p60_t28_g4231 = (cse_p11_t1_g98 * (-(dy)));
          double cse_p60_t29_g4232 = (cse_p11_t1_g98 * (-(dz)));
          double cse_p60_t30_g4233 = (cse_p11_t1_g98 * dx * dx);
          double cse_p60_t31_g4234 = (cse_p11_t1_g98 * dx * dy);
          double cse_p60_t32_g4235 = (cse_p11_t1_g98 * dx * dz);
          double cse_p60_t33_g4236 = (cse_p11_t1_g98 * dy * dy);
          double cse_p60_t34_g4237 = (cse_p11_t1_g98 * dy * dz);
          double cse_p60_t35_g4238 = (cse_p11_t1_g98 * dz * dz);
          double cse_p14_t45_g360 = (cse_p60_t27_g4230 * (-(dx)));
          double cse_p14_t46_g361 = (cse_p60_t27_g4230 * (-(dy)));
          double cse_p14_t47_g362 = (cse_p60_t27_g4230 * (-(dz)));
          double cse_p14_t48_g363 = (cse_p60_t27_g4230 * dx);
          double cse_p14_t49_g364 = (cse_p60_t27_g4230 * dy);
          double cse_p14_t50_g365 = (cse_p60_t27_g4230 * dz);
          double cse_p14_t51_g366 = (cse_p60_t28_g4231 * (-(dy)));
          double cse_p14_t52_g367 = (cse_p60_t28_g4231 * (-(dz)));
          double cse_p14_t53_g368 = (cse_p60_t28_g4231 * dx);
          double cse_p14_t54_g369 = (cse_p60_t28_g4231 * dy);
          double cse_p14_t55_g370 = (cse_p60_t28_g4231 * dz);
          double cse_p14_t56_g371 = (cse_p60_t29_g4232 * (-(dz)));
          double cse_p14_t57_g372 = (cse_p60_t29_g4232 * dx);
          double cse_p14_t58_g373 = (cse_p60_t29_g4232 * dy);
          double cse_p14_t59_g374 = (cse_p60_t29_g4232 * dz);
          double cse_p565_t1_invr41458 = (1.0 / ((rr * rr * rr)));
          double cse_p11_t1_g99 = (cse_p565_t1_invr41458 * de_dr);
          double cse_p60_t36_g4239 = (cse_p11_t1_g99 * (cse_p11_t41_g51 + cse_p60_t78_g4281));
          double cse_p60_t37_g4240 = (cse_p11_t1_g99 * (cse_p11_t42_g52 + cse_p60_t78_g4281));
          double cse_p60_t38_g4241 = (cse_p11_t1_g99 * (cse_p11_t43_g53 + cse_p60_t78_g4281));
          double cse_p60_t39_g4242 = (cse_p11_t1_g99 * cse_p60_t66_g4269);
          double cse_p60_t40_g4243 = (cse_p11_t1_g99 * cse_p60_t67_g4270);
          double cse_p60_t41_g4244 = (cse_p11_t1_g99 * cse_p60_t68_g4271);
          double cse_p60_t42_g4245 = (cse_p11_t1_g99 * dx * dy);
          double cse_p60_t43_g4246 = (cse_p11_t1_g99 * dy * (-(dx)));
          double cse_p60_t44_g4247 = (cse_p11_t1_g99 * dz);
          double cse_p60_t62_g4265 = (cse_p60_t31_g4234 + cse_p60_t43_g4246);
          double cse_p14_t60_g375 = (cse_p60_t44_g4247 * (-(dx)));
          double cse_p14_t61_g376 = (cse_p60_t44_g4247 * (-(dy)));
          double cse_p14_t62_g377 = (cse_p60_t44_g4247 * dx);
          double cse_p14_t63_g378 = (cse_p60_t44_g4247 * dy);
          double cse_p13_t106_g300 = (cse_p60_t30_g4233 + cse_p60_t39_g4242);
          double h_x1_x1 = cse_p13_t106_g300;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t62_g4265;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t87_g281 = (cse_p14_t60_g375 + cse_p60_t32_g4235);
          double cse_p13_t88_g282 = (cse_p14_t61_g376 + cse_p60_t34_g4237);
          double h_x1_z1 = cse_p13_t87_g281;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t78_g272 = (cse_p14_t46_g361 + cse_p60_t43_g4246);
          double cse_p13_t79_g273 = (cse_p14_t47_g362 + cse_p14_t60_g375);
          double cse_p13_t80_g274 = (cse_p14_t48_g363 + cse_p60_t36_g4239);
          double cse_p13_t81_g275 = (cse_p14_t49_g364 + cse_p60_t42_g4245);
          double cse_p13_t82_g276 = (cse_p14_t50_g365 + cse_p14_t62_g377);
          double cse_p13_t107_g301 = (cse_p14_t45_g360 + cse_p60_t39_g4242);
          double h_x1_x2 = cse_p13_t80_g274;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t83_g277 = (cse_p14_t52_g367 + cse_p14_t61_g376);
          double cse_p13_t84_g278 = (cse_p14_t53_g368 + cse_p60_t42_g4245);
          double cse_p13_t85_g279 = (cse_p14_t54_g369 + cse_p60_t37_g4240);
          double cse_p13_t86_g280 = (cse_p14_t55_g370 + cse_p14_t63_g378);
          double h_x1_y2 = cse_p13_t84_g278;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t89_g283 = (cse_p14_t57_g372 + cse_p14_t62_g377);
          double cse_p13_t90_g284 = (cse_p14_t58_g373 + cse_p14_t63_g378);
          double cse_p13_t91_g285 = (cse_p14_t59_g374 + cse_p60_t38_g4241);
          double h_x1_z2 = cse_p13_t89_g283;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t102_g296 = (cse_p60_t33_g4236 + cse_p60_t40_g4243);
          double cse_p13_t103_g297 = (cse_p14_t51_g366 + cse_p60_t40_g4243);
          double h_y1_y1 = cse_p13_t102_g296;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t88_g282;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t81_g275;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t85_g279;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t90_g284;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t104_g298 = (cse_p60_t35_g4238 + cse_p60_t41_g4244);
          double cse_p13_t105_g299 = (cse_p14_t56_g371 + cse_p60_t41_g4244);
          double h_z1_z1 = cse_p13_t104_g298;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t82_g276;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t86_g280;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t91_g285;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t107_g301;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t78_g272;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t79_g273;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t103_g297;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t83_g277;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t105_g299;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t83_g190 = (-(params.eps_solvent));
            double cse_p11_t73_g83 = (cse_p12_t83_g190 + params.eps_core);
            double cse_p11_t82_g92 = (-(term.kqq));
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p11_t1_g11 = (5.00000000000000000e-1 * rx * rx);
            double cse_p11_t83_g93 = (-(rx));
            double ex = exp(cse_p11_t83_g93);
            double gx = (ex * (1.00000000000000000e+0 + cse_p11_t1_g11 + rx));
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p11_t1_g11 + cse_p11_t83_g93));
            double ep = (params.eps_solvent + (cse_p11_t73_g83 * gx));
            double cse_p12_t47_g154 = (ep * ep);
            double cse_p11_t84_g94 = (1.0 / (cse_p12_t47_g154));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t73_g83 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t73_g83 * gx_pp);
            double ee = (term.kqq * (1.0 / (ep)));
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
            double de_dr = (cse_p11_t82_g92 * cse_p11_t84_g94 * params.elec_weight * ep_p);
            double d2e_dr2 = (cse_p11_t82_g92 * params.elec_weight * ((-2.00000000000000000e+0 * ep_p * ep_p * (1.0 / ((ep * ep * ep)))) + (cse_p11_t84_g94 * ep_pp)));
            double cse_p565_t2_invr41459 = (1.0 / (rr));
            double cse_p11_t1_g107 = (cse_p565_t2_invr41459 * de_dr);
            double g_x1 = (cse_p11_t1_g107 * dx);
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double g_y1 = (cse_p11_t1_g107 * dy);
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double g_z1 = (cse_p11_t1_g107 * dz);
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (cse_p11_t1_g107 * (-(dx)));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (cse_p11_t1_g107 * (-(dy)));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (cse_p11_t1_g107 * (-(dz)));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double cse_p11_t1_g98 = (cse_p565_t2_invr41459 * cse_p565_t2_invr41459 * d2e_dr2);
            double cse_p565_t1_invr41458 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g99 = (cse_p565_t1_invr41458 * de_dr);
            double cse_p11_t1_g100 = (rr * rr);
            double cse_p12_t1_g194 = (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t41_g51))));
            double cse_p60_t73_g4276 = (cse_p12_t1_g194 + (cse_p11_t1_g98 * dx * dx));
            double h_x1_x1 = cse_p60_t73_g4276;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = ((cse_p11_t1_g98 * dx * dy) + (cse_p11_t1_g99 * dy * (-(dx))));
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p11_t1_g103 = (cse_p11_t1_g99 * dz);
            double cse_p60_t13_g4216 = (cse_p11_t1_g103 * (-(dx)));
            double cse_p60_t14_g4217 = (cse_p11_t1_g103 * (-(dy)));
            double cse_p60_t15_g4218 = (cse_p11_t1_g103 * dx);
            double cse_p60_t16_g4219 = (cse_p11_t1_g103 * dy);
            double cse_p60_t57_g4260 = (cse_p60_t13_g4216 + (cse_p11_t1_g98 * dx * dz));
            double cse_p60_t58_g4261 = (cse_p60_t14_g4217 + (cse_p11_t1_g98 * dy * dz));
            double h_x1_z1 = cse_p60_t57_g4260;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p11_t1_g101 = (cse_p11_t1_g98 * (-(dx)));
            double cse_p60_t2_g4205 = (cse_p11_t1_g101 * (-(dx)));
            double cse_p60_t3_g4206 = (cse_p11_t1_g101 * (-(dy)));
            double cse_p60_t4_g4207 = (cse_p11_t1_g101 * (-(dz)));
            double cse_p60_t5_g4208 = (cse_p11_t1_g101 * dx);
            double cse_p60_t6_g4209 = (cse_p11_t1_g101 * dy);
            double cse_p60_t7_g4210 = (cse_p11_t1_g101 * dz);
            double cse_p60_t48_g4251 = (cse_p60_t3_g4206 + (cse_p11_t1_g99 * dy * (-(dx))));
            double cse_p60_t49_g4252 = (cse_p60_t13_g4216 + cse_p60_t4_g4207);
            double cse_p60_t50_g4253 = (cse_p60_t5_g4208 + (cse_p11_t1_g99 * (cse_p11_t41_g51 + (-(cse_p11_t1_g100)))));
            double cse_p60_t51_g4254 = (cse_p60_t6_g4209 + (cse_p11_t1_g99 * dx * dy));
            double cse_p60_t52_g4255 = (cse_p60_t15_g4218 + cse_p60_t7_g4210);
            double cse_p60_t74_g4277 = (cse_p12_t1_g194 + cse_p60_t2_g4205);
            double h_x1_x2 = cse_p60_t50_g4253;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p11_t1_g102 = (cse_p11_t1_g98 * (-(dy)));
            double cse_p60_t8_g4211 = (cse_p11_t1_g102 * (-(dy)));
            double cse_p60_t9_g4212 = (cse_p11_t1_g102 * (-(dz)));
            double cse_p60_t10_g4213 = (cse_p11_t1_g102 * dx);
            double cse_p60_t11_g4214 = (cse_p11_t1_g102 * dy);
            double cse_p60_t12_g4215 = (cse_p11_t1_g102 * dz);
            double cse_p60_t53_g4256 = (cse_p60_t14_g4217 + cse_p60_t9_g4212);
            double cse_p60_t54_g4257 = (cse_p60_t10_g4213 + (cse_p11_t1_g99 * dx * dy));
            double cse_p60_t55_g4258 = (cse_p60_t11_g4214 + (cse_p11_t1_g99 * (cse_p11_t42_g52 + (-(cse_p11_t1_g100)))));
            double cse_p60_t56_g4259 = (cse_p60_t12_g4215 + cse_p60_t16_g4219);
            double h_x1_y2 = cse_p60_t54_g4257;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p11_t1_g104 = (cse_p11_t1_g98 * (-(dz)));
            double cse_p60_t17_g4220 = (cse_p11_t1_g104 * (-(dz)));
            double cse_p60_t18_g4221 = (cse_p11_t1_g104 * dx);
            double cse_p60_t19_g4222 = (cse_p11_t1_g104 * dy);
            double cse_p60_t20_g4223 = (cse_p11_t1_g104 * dz);
            double cse_p60_t59_g4262 = (cse_p60_t15_g4218 + cse_p60_t18_g4221);
            double cse_p60_t60_g4263 = (cse_p60_t16_g4219 + cse_p60_t19_g4222);
            double cse_p60_t61_g4264 = (cse_p60_t20_g4223 + (cse_p11_t1_g99 * (cse_p11_t43_g53 + (-(cse_p11_t1_g100)))));
            double h_x1_z2 = cse_p60_t59_g4262;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p11_t1_g105 = (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t42_g52))));
            double cse_p60_t69_g4272 = (cse_p11_t1_g105 + (cse_p11_t1_g98 * dy * dy));
            double cse_p60_t70_g4273 = (cse_p11_t1_g105 + cse_p60_t8_g4211);
            double h_y1_y1 = cse_p60_t69_g4272;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t58_g4261;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t51_g4254;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t55_g4258;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t60_g4263;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p11_t1_g106 = (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t43_g53))));
            double cse_p60_t71_g4274 = (cse_p11_t1_g106 + (cse_p11_t1_g98 * dz * dz));
            double cse_p60_t72_g4275 = (cse_p11_t1_g106 + cse_p60_t17_g4220);
            double h_z1_z1 = cse_p60_t71_g4274;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t52_g4255;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t56_g4259;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t61_g4264;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t74_g4277;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t48_g4251;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t49_g4252;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t70_g4273;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t53_g4256;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t72_g4275;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p11_t1_g100 = (rr * rr);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p11_t1_g100 * (term.bb_high + (term.aa_high * rr)))));
              energy_added += energy;
              *energy_accumulate += energy;
              double de_dr = (params.elec_weight * (term.cc_high + (2.00000000000000000e+0 * term.bb_high * rr) + (3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g100)));
              double d2e_dr2 = (params.elec_weight * ((2.00000000000000000e+0 * term.bb_high) + (6.00000000000000000e+0 * term.aa_high * rr)));
              double cse_p565_t2_invr41459 = (1.0 / (rr));
              double cse_p11_t1_g107 = (cse_p565_t2_invr41459 * de_dr);
              double g_x1 = (cse_p11_t1_g107 * dx);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t1_g107 * dy);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t1_g107 * dz);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = (cse_p11_t1_g107 * (-(dx)));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (cse_p11_t1_g107 * (-(dy)));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (cse_p11_t1_g107 * (-(dz)));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double cse_p11_t1_g98 = (cse_p565_t2_invr41459 * cse_p565_t2_invr41459 * d2e_dr2);
              double cse_p565_t1_invr41458 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g99 = (cse_p565_t1_invr41458 * de_dr);
              double cse_p12_t1_g194 = (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t41_g51))));
              double h_x1_x1 = (cse_p12_t1_g194 + (cse_p11_t1_g98 * dx * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double h_x1_y1 = ((cse_p11_t1_g98 * dx * dy) + (cse_p11_t1_g99 * dy * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g103 = (cse_p11_t1_g99 * dz);
              double h_x1_z1 = ((cse_p11_t1_g103 * (-(dx))) + (cse_p11_t1_g98 * dx * dz));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g101 = (cse_p11_t1_g98 * (-(dx)));
              double h_x1_x2 = ((cse_p11_t1_g101 * dx) + (cse_p11_t1_g99 * (cse_p11_t41_g51 + (-(cse_p11_t1_g100)))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g102 = (cse_p11_t1_g98 * (-(dy)));
              double h_x1_y2 = ((cse_p11_t1_g102 * dx) + (cse_p11_t1_g99 * dx * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g104 = (cse_p11_t1_g98 * (-(dz)));
              double h_x1_z2 = ((cse_p11_t1_g103 * dx) + (cse_p11_t1_g104 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g105 = (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t42_g52))));
              double h_y1_y1 = (cse_p11_t1_g105 + (cse_p11_t1_g98 * dy * dy));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = ((cse_p11_t1_g103 * (-(dy))) + (cse_p11_t1_g98 * dy * dz));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = ((cse_p11_t1_g101 * dy) + (cse_p11_t1_g99 * dx * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = ((cse_p11_t1_g102 * dy) + (cse_p11_t1_g99 * (cse_p11_t42_g52 + (-(cse_p11_t1_g100)))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = ((cse_p11_t1_g103 * dy) + (cse_p11_t1_g104 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g106 = (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t43_g53))));
              double h_z1_z1 = (cse_p11_t1_g106 + (cse_p11_t1_g98 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = ((cse_p11_t1_g101 * dz) + (cse_p11_t1_g103 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = ((cse_p11_t1_g102 * dz) + (cse_p11_t1_g103 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = ((cse_p11_t1_g104 * dz) + (cse_p11_t1_g99 * (cse_p11_t43_g53 + (-(cse_p11_t1_g100)))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p12_t1_g194 + (cse_p11_t1_g101 * (-(dx))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = ((cse_p11_t1_g101 * (-(dy))) + (cse_p11_t1_g99 * dy * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = ((cse_p11_t1_g101 * (-(dz))) + (cse_p11_t1_g103 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = (cse_p11_t1_g105 + (cse_p11_t1_g102 * (-(dy))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = ((cse_p11_t1_g102 * (-(dz))) + (cse_p11_t1_g103 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = (cse_p11_t1_g106 + (cse_p11_t1_g104 * (-(dz))));
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
                double cse_p565_t2_invr41459 = (1.0 / (rr));
                double g_y1 = (cse_p565_t2_invr41459 * de_dr * dy);
                KernelGradientAcc(term.i3x1, 1, g_y1);
                double cse_p11_t1_g107 = (cse_p565_t2_invr41459 * de_dr);
                double g_z1 = (cse_p11_t1_g107 * dz);
                KernelGradientAcc(term.i3x1, 2, g_z1);
                double g_x2 = (cse_p11_t1_g107 * (-(dx)));
                KernelGradientAcc(term.i3x2, 0, g_x2);
                double g_y2 = (cse_p11_t1_g107 * (-(dy)));
                KernelGradientAcc(term.i3x2, 1, g_y2);
                double g_z2 = (cse_p11_t1_g107 * (-(dz)));
                KernelGradientAcc(term.i3x2, 2, g_z2);
                double h_x1_x1 = ((cse_p565_t2_invr41459 * cse_p565_t2_invr41459 * d2e_dr2 * dx * dx) + (de_dr * ((rr * rr) + (-(cse_p11_t41_g51))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t1_invr41458 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g100 = (rr * rr);
                double cse_p11_t1_g99 = (cse_p565_t1_invr41458 * de_dr);
                double cse_p11_t1_g98 = (cse_p565_t2_invr41459 * cse_p565_t2_invr41459 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g98 * dx * dy) + (cse_p11_t1_g99 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double h_x1_z1 = ((cse_p11_t1_g98 * dx * dz) + (cse_p11_t1_g99 * dz * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g103 = (cse_p11_t1_g99 * dz);
                double h_x1_x2 = ((cse_p11_t1_g98 * dx * (-(dx))) + (cse_p11_t1_g99 * (cse_p11_t41_g51 + (-(cse_p11_t1_g100)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g101 = (cse_p11_t1_g98 * (-(dx)));
                double h_x1_y2 = ((cse_p11_t1_g98 * dx * (-(dy))) + (cse_p11_t1_g99 * dx * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g102 = (cse_p11_t1_g98 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g103 * dx) + (cse_p11_t1_g98 * dx * (-(dz))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g104 = (cse_p11_t1_g98 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g98 * dy * dy) + (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t42_g52)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g105 = (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t42_g52))));
                double h_y1_z1 = ((cse_p11_t1_g103 * (-(dy))) + (cse_p11_t1_g98 * dy * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g101 * dy) + (cse_p11_t1_g99 * dx * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g102 * dy) + (cse_p11_t1_g99 * (cse_p11_t42_g52 + (-(cse_p11_t1_g100)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g103 * dy) + (cse_p11_t1_g104 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g98 * dz * dz) + (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t43_g53)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g106 = (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t43_g53))));
                double h_z1_x2 = ((cse_p11_t1_g101 * dz) + (cse_p11_t1_g103 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g102 * dz) + (cse_p11_t1_g103 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g104 * dz) + (cse_p11_t1_g99 * (cse_p11_t43_g53 + (-(cse_p11_t1_g100)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g101 * (-(dx))) + (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t41_g51)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p12_t1_g194 = (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t41_g51))));
                double h_x2_y2 = ((cse_p11_t1_g101 * (-(dy))) + (cse_p11_t1_g99 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double h_x2_z2 = ((cse_p11_t1_g101 * (-(dz))) + (cse_p11_t1_g103 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g105 + (cse_p11_t1_g102 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g102 * (-(dz))) + (cse_p11_t1_g103 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g106 + (cse_p11_t1_g104 * (-(dz))));
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
void hessian_fd(const rosetta_elec_parameters& params, const rosetta_elec_term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
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
