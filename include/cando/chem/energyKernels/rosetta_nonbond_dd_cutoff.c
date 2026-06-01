#include "rosetta_nonbond_dd_cutoff.h"

template <typename HESSIAN>
struct Rosetta_Nonbond_Dd_Cutoff {
  static constexpr size_t PositionSize = 6;
  static std::string description() { return "mathkernel-rosetta_nonbond_dd_cutoff"; };
double gradient(const rosetta_nonbond_parameters& params, const rosetta_nonbond_term& term, double* position, double* energy_accumulate, double* force) {
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
      if ((rr <= (6.00000023841857900e-1 * term.sigma))) {
        /* !COND-1 */
        double fa_rep = (term.bij + (term.mij * rr));
        double fa_atr = (-(term.epsilon));
        double de_dr = (term.mij * params.rep_weight);
        double energy = (fa_atr + (fa_rep * params.rep_weight));
        energy_added += energy;
        *energy_accumulate += energy;
        double cse_p565_t1_invr17055 = (1.0 / (rr));
        double g_x1 = (cse_p565_t1_invr17055 * de_dr * dx);
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = (cse_p565_t1_invr17055 * de_dr * dy);
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = (cse_p565_t1_invr17055 * de_dr * dz);
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = (cse_p565_t1_invr17055 * de_dr * (-(dx)));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (cse_p565_t1_invr17055 * de_dr * (-(dy)));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (cse_p565_t1_invr17055 * de_dr * (-(dz)));
        KernelGradientAcc(term.i3x2, 2, g_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.sigma)) {
          /* !COND-2 */
          double cse_p565_t1_invr17055 = (1.0 / (rr));
          double sr = (cse_p565_t1_invr17055 * term.sigma);
          double sr2 = (sr * sr);
          double sr6 = (sr2 * sr2 * sr2);
          double sr12 = (sr6 * sr6);
          double fa_rep_epsilon = (term.epsilon * (sr12 + (-2.00000000000000000e+0 * sr6)));
          double fa_rep = (term.epsilon + fa_rep_epsilon);
          double fa_atr = (-(term.epsilon));
          double de_dr = (cse_p565_t1_invr17055 * term.epsilon * params.rep_weight * ((-1.20000000000000000e+1 * sr12) + (1.20000000000000000e+1 * sr6)));
          double energy = (fa_atr + (fa_rep * params.rep_weight));
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (cse_p565_t1_invr17055 * de_dr * dx);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p565_t1_invr17055 * de_dr * dy);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p565_t1_invr17055 * de_dr * dz);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p565_t1_invr17055 * de_dr * (-(dx)));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p565_t1_invr17055 * de_dr * (-(dy)));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p565_t1_invr17055 * de_dr * (-(dz)));
          KernelGradientAcc(term.i3x2, 2, g_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rswitch)) {
            /* !COND-3 */
            double cse_p565_t1_invr17055 = (1.0 / (rr));
            double sr = (cse_p565_t1_invr17055 * term.sigma);
            double sr2 = (sr * sr);
            double sr6 = (sr2 * sr2 * sr2);
            double sr12 = (sr6 * sr6);
            double fa_atr = (term.epsilon * (sr12 + (-2.00000000000000000e+0 * sr6)));
            double de_dr = (cse_p565_t1_invr17055 * term.epsilon * ((-1.20000000000000000e+1 * sr12) + (1.20000000000000000e+1 * sr6)));
            double energy = fa_atr;
            energy_added += energy;
            *energy_accumulate += energy;
            double g_x1 = (cse_p565_t1_invr17055 * de_dr * dx);
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double g_y1 = (cse_p565_t1_invr17055 * de_dr * dy);
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double g_z1 = (cse_p565_t1_invr17055 * de_dr * dz);
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (cse_p565_t1_invr17055 * de_dr * (-(dx)));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (cse_p565_t1_invr17055 * de_dr * (-(dy)));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (cse_p565_t1_invr17055 * de_dr * (-(dz)));
            KernelGradientAcc(term.i3x2, 2, g_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double fa_atr = (term.dd + (term.cc * rr) + (rr * rr * (term.bb + (term.aa * rr))));
              double de_dr = (term.cc + (2.00000000000000000e+0 * term.bb * rr) + (3.00000000000000000e+0 * term.aa * rr * rr));
              double energy = fa_atr;
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p565_t1_invr17055 = (1.0 / (rr));
              double g_x1 = (cse_p565_t1_invr17055 * de_dr * dx);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p565_t1_invr17055 * de_dr * dy);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p565_t1_invr17055 * de_dr * dz);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = (cse_p565_t1_invr17055 * de_dr * (-(dx)));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (cse_p565_t1_invr17055 * de_dr * (-(dy)));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (cse_p565_t1_invr17055 * de_dr * (-(dz)));
              KernelGradientAcc(term.i3x2, 2, g_z2);
            } else {
              /* !COND-ELSE-4 */
              {
                /* !COND-ELSE-5 */
                double de_dr = 0.00000000000000000e+0;
                double energy = 0.00000000000000000e+0;
                energy_added += energy;
                *energy_accumulate += energy;
                double g_x1 = (de_dr * dx * (1.0 / (rr)));
                KernelGradientAcc(term.i3x1, 0, g_x1);
                double cse_p565_t1_invr17055 = (1.0 / (rr));
                double g_y1 = (cse_p565_t1_invr17055 * de_dr * dy);
                KernelGradientAcc(term.i3x1, 1, g_y1);
                double g_z1 = (cse_p565_t1_invr17055 * de_dr * dz);
                KernelGradientAcc(term.i3x1, 2, g_z1);
                double g_x2 = (cse_p565_t1_invr17055 * de_dr * (-(dx)));
                KernelGradientAcc(term.i3x2, 0, g_x2);
                double g_y2 = (cse_p565_t1_invr17055 * de_dr * (-(dy)));
                KernelGradientAcc(term.i3x2, 1, g_y2);
                double g_z2 = (cse_p565_t1_invr17055 * de_dr * (-(dz)));
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
void gradient_fd(const rosetta_nonbond_parameters& params, const rosetta_nonbond_term& term, double* position, double* energy_accumulate, double* force)
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

double energy(const rosetta_nonbond_parameters& params, const rosetta_nonbond_term& term, double* position, double* energy_accumulate) {
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
      if ((rr <= (6.00000023841857900e-1 * term.sigma))) {
        /* !COND-1 */
        double fa_rep = (term.bij + (term.mij * rr));
        double fa_atr = (-(term.epsilon));
        double energy = (fa_atr + (fa_rep * params.rep_weight));
        energy_added += energy;
        *energy_accumulate += energy;
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.sigma)) {
          /* !COND-2 */
          double sr = (term.sigma * (1.0 / (rr)));
          double sr2 = (sr * sr);
          double sr6 = (sr2 * sr2 * sr2);
          double sr12 = (sr6 * sr6);
          double fa_rep_epsilon = (term.epsilon * (sr12 + (-2.00000000000000000e+0 * sr6)));
          double fa_rep = (term.epsilon + fa_rep_epsilon);
          double fa_atr = (-(term.epsilon));
          double energy = (fa_atr + (fa_rep * params.rep_weight));
          energy_added += energy;
          *energy_accumulate += energy;
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rswitch)) {
            /* !COND-3 */
            double sr = (term.sigma * (1.0 / (rr)));
            double sr2 = (sr * sr);
            double sr6 = (sr2 * sr2 * sr2);
            double sr12 = (sr6 * sr6);
            double fa_atr = (term.epsilon * (sr12 + (-2.00000000000000000e+0 * sr6)));
            double energy = fa_atr;
            energy_added += energy;
            *energy_accumulate += energy;
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double fa_atr = (term.dd + (term.cc * rr) + (rr * rr * (term.bb + (term.aa * rr))));
              double energy = fa_atr;
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
void energy_fd(const rosetta_nonbond_parameters& params, const rosetta_nonbond_term& term, double* position, double* energy_accumulate)
{
  energy(params, term, position, energy_accumulate);
}

double hessian(const rosetta_nonbond_parameters& params, const rosetta_nonbond_term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
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
      double cse_p11_t46_g56 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t47_g57 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t48_g58 = (dz * dz);
      double r2 = (cse_p11_t46_g56 + cse_p11_t47_g57 + cse_p11_t48_g58);
      double rr = sqrt(r2);
      if ((rr <= (6.00000023841857900e-1 * term.sigma))) {
        /* !COND-1 */
        double cse_p60_t55_g5082 = (rr * rr);
        double cse_p60_t90_g5117 = (-(cse_p11_t46_g56));
        double cse_p60_t91_g5118 = (-(cse_p11_t47_g57));
        double cse_p60_t92_g5119 = (-(cse_p11_t48_g58));
        double cse_p60_t94_g5121 = (-(dx));
        double cse_p60_t95_g5122 = (-(dy));
        double cse_p60_t96_g5123 = (-(dz));
        double cse_p60_t97_g5124 = (-(term.epsilon));
        double cse_p565_t1_invr47717 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr47718 = (1.0 / (rr));
        double cse_p60_t101_g5128 = (rr * rr * rr);
        double cse_p15_t94_g575 = (-(cse_p60_t55_g5082));
        double cse_p14_t96_g453 = (cse_p11_t46_g56 + cse_p15_t94_g575);
        double cse_p14_t98_g455 = (cse_p11_t47_g57 + cse_p15_t94_g575);
        double cse_p14_t100_g457 = (cse_p11_t48_g58 + cse_p15_t94_g575);
        double cse_p13_t121_g340 = (cse_p60_t55_g5082 + cse_p60_t90_g5117);
        double cse_p13_t122_g341 = (cse_p60_t55_g5082 + cse_p60_t91_g5118);
        double cse_p13_t123_g342 = (cse_p60_t55_g5082 + cse_p60_t92_g5119);
        double fa_rep = (term.bij + (term.mij * rr));
        double cse_p60_t54_g5081 = (fa_rep * params.rep_weight);
        double fa_atr = cse_p60_t97_g5124;
        double cse_p60_t88_g5115 = (cse_p60_t54_g5081 + fa_atr);
        double de_dr = (term.mij * params.rep_weight);
        double cse_p60_t52_g5079 = (cse_p565_t1_invr47717 * de_dr);
        double cse_p13_t65_g284 = (cse_p14_t96_g453 * cse_p60_t52_g5079);
        double cse_p13_t66_g285 = (cse_p14_t98_g455 * cse_p60_t52_g5079);
        double cse_p13_t67_g286 = (cse_p14_t100_g457 * cse_p60_t52_g5079);
        double cse_p13_t68_g287 = (cse_p60_t52_g5079 * dx * dy);
        double cse_p13_t69_g288 = (cse_p60_t52_g5079 * cse_p60_t94_g5121 * dy);
        double cse_p12_t1_g219 = (cse_p13_t121_g340 * cse_p60_t52_g5079);
        double cse_p12_t21_g138 = (cse_p13_t122_g341 * cse_p60_t52_g5079);
        double cse_p12_t22_g139 = (cse_p13_t123_g342 * cse_p60_t52_g5079);
        double cse_p12_t25_g142 = (cse_p60_t52_g5079 * dz);
        double cse_p13_t55_g274 = (cse_p12_t25_g142 * cse_p60_t94_g5121);
        double cse_p13_t56_g275 = (cse_p12_t25_g142 * cse_p60_t95_g5122);
        double cse_p13_t57_g276 = (cse_p12_t25_g142 * dx);
        double cse_p13_t58_g277 = (cse_p12_t25_g142 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p60_t51_g5078 = (cse_p565_t2_invr47718 * cse_p565_t2_invr47718 * d2e_dr2);
        double cse_p13_t59_g278 = (cse_p60_t51_g5078 * dx * dx);
        double cse_p13_t60_g279 = (cse_p60_t51_g5078 * dx * dy);
        double cse_p13_t61_g280 = (cse_p60_t51_g5078 * dx * dz);
        double cse_p13_t62_g281 = (cse_p60_t51_g5078 * dy * dy);
        double cse_p13_t63_g282 = (cse_p60_t51_g5078 * dy * dz);
        double cse_p13_t64_g283 = (cse_p60_t51_g5078 * dz * dz);
        double cse_p12_t8_g125 = (cse_p60_t51_g5078 * cse_p60_t94_g5121);
        double cse_p13_t70_g289 = (cse_p12_t8_g125 * cse_p60_t94_g5121);
        double cse_p13_t71_g290 = (cse_p12_t8_g125 * cse_p60_t95_g5122);
        double cse_p13_t72_g291 = (cse_p12_t8_g125 * cse_p60_t96_g5123);
        double cse_p13_t73_g292 = (cse_p12_t8_g125 * dx);
        double cse_p13_t74_g293 = (cse_p12_t8_g125 * dy);
        double cse_p13_t75_g294 = (cse_p12_t8_g125 * dz);
        double cse_p12_t9_g126 = (cse_p60_t51_g5078 * cse_p60_t95_g5122);
        double cse_p13_t76_g295 = (cse_p12_t9_g126 * cse_p60_t95_g5122);
        double cse_p13_t77_g296 = (cse_p12_t9_g126 * cse_p60_t96_g5123);
        double cse_p13_t78_g297 = (cse_p12_t9_g126 * dx);
        double cse_p13_t79_g298 = (cse_p12_t9_g126 * dy);
        double cse_p13_t80_g299 = (cse_p12_t9_g126 * dz);
        double cse_p12_t10_g127 = (cse_p60_t51_g5078 * cse_p60_t96_g5123);
        double cse_p13_t51_g270 = (cse_p12_t10_g127 * cse_p60_t96_g5123);
        double cse_p13_t52_g271 = (cse_p12_t10_g127 * dx);
        double cse_p13_t53_g272 = (cse_p12_t10_g127 * dy);
        double cse_p13_t54_g273 = (cse_p12_t10_g127 * dz);
        double cse_p12_t65_g182 = (cse_p12_t1_g219 + cse_p13_t59_g278);
        double cse_p12_t66_g183 = (cse_p13_t60_g279 + cse_p13_t69_g288);
        double cse_p12_t67_g184 = (cse_p13_t55_g274 + cse_p13_t61_g280);
        double cse_p12_t68_g185 = (cse_p13_t56_g275 + cse_p13_t63_g282);
        double cse_p12_t84_g201 = (cse_p12_t21_g138 + cse_p13_t62_g281);
        double cse_p12_t86_g203 = (cse_p12_t22_g139 + cse_p13_t64_g283);
        double cse_p12_t59_g176 = (cse_p13_t65_g284 + cse_p13_t73_g292);
        double cse_p12_t62_g179 = (cse_p12_t1_g219 + cse_p13_t70_g289);
        double cse_p12_t69_g186 = (cse_p13_t68_g287 + cse_p13_t74_g293);
        double cse_p12_t71_g188 = (cse_p13_t69_g288 + cse_p13_t71_g290);
        double cse_p12_t72_g189 = (cse_p13_t55_g274 + cse_p13_t72_g291);
        double cse_p12_t73_g190 = (cse_p13_t57_g276 + cse_p13_t75_g294);
        double cse_p12_t60_g177 = (cse_p13_t66_g285 + cse_p13_t79_g298);
        double cse_p12_t70_g187 = (cse_p13_t68_g287 + cse_p13_t78_g297);
        double cse_p12_t74_g191 = (cse_p13_t56_g275 + cse_p13_t77_g296);
        double cse_p12_t75_g192 = (cse_p13_t58_g277 + cse_p13_t80_g299);
        double cse_p12_t85_g202 = (cse_p12_t21_g138 + cse_p13_t76_g295);
        double cse_p12_t61_g178 = (cse_p13_t54_g273 + cse_p13_t67_g286);
        double cse_p12_t76_g193 = (cse_p13_t52_g271 + cse_p13_t57_g276);
        double cse_p12_t77_g194 = (cse_p13_t53_g272 + cse_p13_t58_g277);
        double cse_p12_t87_g204 = (cse_p12_t22_g139 + cse_p13_t51_g270);
        double energy = cse_p60_t88_g5115;
        energy_added += energy;
        *energy_accumulate += energy;
        double cse_p11_t1_g117 = (cse_p565_t2_invr47718 * de_dr);
        double cse_p60_t45_g5072 = (cse_p11_t1_g117 * cse_p60_t94_g5121);
        double cse_p60_t46_g5073 = (cse_p11_t1_g117 * cse_p60_t95_g5122);
        double cse_p60_t47_g5074 = (cse_p11_t1_g117 * cse_p60_t96_g5123);
        double cse_p60_t48_g5075 = (cse_p11_t1_g117 * dx);
        double cse_p60_t49_g5076 = (cse_p11_t1_g117 * dy);
        double cse_p60_t50_g5077 = (cse_p11_t1_g117 * dz);
        double g_x1 = cse_p60_t48_g5075;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = cse_p60_t49_g5076;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = cse_p60_t50_g5077;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = cse_p60_t45_g5072;
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = cse_p60_t46_g5073;
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = cse_p60_t47_g5074;
        KernelGradientAcc(term.i3x2, 2, g_z2);
        double h_x1_x1 = cse_p12_t65_g182;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t66_g183;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t67_g184;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t59_g176;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t70_g187;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t76_g193;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t84_g201;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t68_g185;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t69_g186;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t60_g177;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t77_g194;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t86_g203;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t73_g190;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t75_g192;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t61_g178;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t62_g179;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t71_g188;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t72_g189;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t85_g202;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t74_g191;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t87_g204;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.sigma)) {
          /* !COND-2 */
          double cse_p11_t1_g110 = (rr * rr);
          double cse_p60_t76_g5103 = (cse_p11_t46_g56 + (-(cse_p11_t1_g110)));
          double cse_p60_t77_g5104 = (cse_p11_t47_g57 + (-(cse_p11_t1_g110)));
          double cse_p60_t78_g5105 = (cse_p11_t48_g58 + (-(cse_p11_t1_g110)));
          double cse_p60_t79_g5106 = (cse_p11_t1_g110 + (-(cse_p11_t46_g56)));
          double cse_p60_t80_g5107 = (cse_p11_t1_g110 + (-(cse_p11_t47_g57)));
          double cse_p60_t81_g5108 = (cse_p11_t1_g110 + (-(cse_p11_t48_g58)));
          double cse_p60_t93_g5120 = (-(cse_p11_t1_g110));
          double cse_p565_t3_invr47719 = (1.0 / (cse_p11_t1_g110));
          double cse_p565_t2_invr47718 = (1.0 / (rr));
          double cse_p60_t7_g5034 = (cse_p565_t2_invr47718 * term.sigma);
          double sr = cse_p60_t7_g5034;
          double cse_p60_t56_g5083 = (sr * sr);
          double sr2 = cse_p60_t56_g5083;
          double cse_p60_t57_g5084 = (sr2 * sr2 * sr2);
          double sr6 = cse_p60_t57_g5084;
          double cse_p60_t2_g5029 = (-2.00000000000000000e+0 * sr6);
          double cse_p60_t3_g5030 = (-4.20000000000000040e+1 * sr6);
          double cse_p60_t4_g5031 = (1.20000000000000000e+1 * sr6);
          double cse_p60_t58_g5085 = (sr6 * sr6);
          double sr12 = cse_p60_t58_g5085;
          double cse_p60_t1_g5028 = (-1.20000000000000000e+1 * sr12);
          double cse_p60_t5_g5032 = (1.56000000000000000e+2 * sr12);
          double cse_p60_t53_g5080 = (term.epsilon * (cse_p60_t2_g5029 + sr12));
          double cse_p60_t62_g5089 = (cse_p60_t1_g5028 + cse_p60_t4_g5031);
          double cse_p60_t63_g5090 = (cse_p60_t3_g5030 + cse_p60_t5_g5032);
          double cse_p60_t89_g5116 = (cse_p60_t2_g5029 + sr12);
          double fa_rep_epsilon = cse_p60_t53_g5080;
          double fa_rep = (term.epsilon + fa_rep_epsilon);
          double fa_atr = (-(term.epsilon));
          double de_dr = (cse_p565_t2_invr47718 * cse_p60_t62_g5089 * term.epsilon * params.rep_weight);
          double cse_p60_t6_g5033 = (cse_p565_t2_invr47718 * de_dr);
          double d2e_dr2 = (cse_p565_t3_invr47719 * cse_p60_t63_g5090 * term.epsilon * params.rep_weight);
          double energy = (fa_atr + (fa_rep * params.rep_weight));
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (cse_p60_t6_g5033 * dx);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p60_t6_g5033 * dy);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p60_t6_g5033 * dz);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p60_t6_g5033 * (-(dx)));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p60_t6_g5033 * (-(dy)));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p60_t6_g5033 * (-(dz)));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double cse_p11_t1_g108 = (cse_p565_t2_invr47718 * cse_p565_t2_invr47718 * d2e_dr2);
          double cse_p60_t8_g5035 = (cse_p11_t1_g108 * (-(dx)));
          double cse_p60_t9_g5036 = (cse_p11_t1_g108 * (-(dy)));
          double cse_p60_t10_g5037 = (cse_p11_t1_g108 * (-(dz)));
          double cse_p60_t11_g5038 = (cse_p11_t1_g108 * dx * dx);
          double cse_p60_t12_g5039 = (cse_p11_t1_g108 * dx * dy);
          double cse_p60_t13_g5040 = (cse_p11_t1_g108 * dx * dz);
          double cse_p60_t14_g5041 = (cse_p11_t1_g108 * dy * dy);
          double cse_p60_t15_g5042 = (cse_p11_t1_g108 * dy * dz);
          double cse_p60_t16_g5043 = (cse_p11_t1_g108 * dz * dz);
          double cse_p14_t59_g416 = (cse_p60_t8_g5035 * (-(dx)));
          double cse_p14_t60_g417 = (cse_p60_t8_g5035 * (-(dy)));
          double cse_p14_t61_g418 = (cse_p60_t8_g5035 * (-(dz)));
          double cse_p14_t62_g419 = (cse_p60_t8_g5035 * dx);
          double cse_p14_t63_g420 = (cse_p60_t8_g5035 * dy);
          double cse_p14_t64_g421 = (cse_p60_t8_g5035 * dz);
          double cse_p14_t65_g422 = (cse_p60_t9_g5036 * (-(dy)));
          double cse_p14_t66_g423 = (cse_p60_t9_g5036 * (-(dz)));
          double cse_p14_t67_g424 = (cse_p60_t9_g5036 * dx);
          double cse_p14_t68_g425 = (cse_p60_t9_g5036 * dy);
          double cse_p14_t69_g426 = (cse_p60_t9_g5036 * dz);
          double cse_p14_t51_g408 = (cse_p60_t10_g5037 * (-(dz)));
          double cse_p14_t52_g409 = (cse_p60_t10_g5037 * dx);
          double cse_p14_t53_g410 = (cse_p60_t10_g5037 * dy);
          double cse_p14_t54_g411 = (cse_p60_t10_g5037 * dz);
          double cse_p565_t1_invr47717 = (1.0 / ((rr * rr * rr)));
          double cse_p11_t1_g109 = (cse_p565_t1_invr47717 * de_dr);
          double cse_p60_t17_g5044 = (cse_p11_t1_g109 * (cse_p11_t46_g56 + cse_p60_t93_g5120));
          double cse_p60_t18_g5045 = (cse_p11_t1_g109 * (cse_p11_t47_g57 + cse_p60_t93_g5120));
          double cse_p60_t19_g5046 = (cse_p11_t1_g109 * (cse_p11_t48_g58 + cse_p60_t93_g5120));
          double cse_p60_t20_g5047 = (cse_p11_t1_g109 * cse_p60_t79_g5106);
          double cse_p60_t21_g5048 = (cse_p11_t1_g109 * cse_p60_t80_g5107);
          double cse_p60_t22_g5049 = (cse_p11_t1_g109 * cse_p60_t81_g5108);
          double cse_p60_t23_g5050 = (cse_p11_t1_g109 * dx * dy);
          double cse_p60_t24_g5051 = (cse_p11_t1_g109 * dy * (-(dx)));
          double cse_p60_t25_g5052 = (cse_p11_t1_g109 * dz);
          double cse_p60_t64_g5091 = (cse_p60_t12_g5039 + cse_p60_t24_g5051);
          double cse_p14_t55_g412 = (cse_p60_t25_g5052 * (-(dx)));
          double cse_p14_t56_g413 = (cse_p60_t25_g5052 * (-(dy)));
          double cse_p14_t57_g414 = (cse_p60_t25_g5052 * dx);
          double cse_p14_t58_g415 = (cse_p60_t25_g5052 * dy);
          double cse_p13_t119_g338 = (cse_p60_t11_g5038 + cse_p60_t20_g5047);
          double h_x1_x1 = cse_p13_t119_g338;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t64_g5091;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t95_g314 = (cse_p14_t55_g412 + cse_p60_t13_g5040);
          double cse_p13_t96_g315 = (cse_p14_t56_g413 + cse_p60_t15_g5042);
          double h_x1_z1 = cse_p13_t95_g314;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t89_g308 = (cse_p14_t62_g419 + cse_p60_t17_g5044);
          double cse_p13_t97_g316 = (cse_p14_t63_g420 + cse_p60_t23_g5050);
          double cse_p13_t99_g318 = (cse_p14_t60_g417 + cse_p60_t24_g5051);
          double cse_p13_t100_g319 = (cse_p14_t55_g412 + cse_p14_t61_g418);
          double cse_p13_t101_g320 = (cse_p14_t57_g414 + cse_p14_t64_g421);
          double cse_p13_t120_g339 = (cse_p14_t59_g416 + cse_p60_t20_g5047);
          double h_x1_x2 = cse_p13_t89_g308;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t90_g309 = (cse_p14_t68_g425 + cse_p60_t18_g5045);
          double cse_p13_t98_g317 = (cse_p14_t67_g424 + cse_p60_t23_g5050);
          double cse_p13_t102_g321 = (cse_p14_t56_g413 + cse_p14_t66_g423);
          double cse_p13_t103_g322 = (cse_p14_t58_g415 + cse_p14_t69_g426);
          double h_x1_y2 = cse_p13_t98_g317;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t91_g310 = (cse_p14_t54_g411 + cse_p60_t19_g5046);
          double cse_p13_t104_g323 = (cse_p14_t52_g409 + cse_p14_t57_g414);
          double cse_p13_t105_g324 = (cse_p14_t53_g410 + cse_p14_t58_g415);
          double h_x1_z2 = cse_p13_t104_g323;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t115_g334 = (cse_p60_t14_g5041 + cse_p60_t21_g5048);
          double cse_p13_t116_g335 = (cse_p14_t65_g422 + cse_p60_t21_g5048);
          double h_y1_y1 = cse_p13_t115_g334;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t96_g315;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t97_g316;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t90_g309;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t105_g324;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t117_g336 = (cse_p60_t16_g5043 + cse_p60_t22_g5049);
          double cse_p13_t118_g337 = (cse_p14_t51_g408 + cse_p60_t22_g5049);
          double h_z1_z1 = cse_p13_t117_g336;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t101_g320;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t103_g322;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t91_g310;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t120_g339;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t99_g318;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t100_g319;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t116_g335;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t102_g321;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t118_g337;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rswitch)) {
            /* !COND-3 */
            double cse_p565_t2_invr47718 = (1.0 / (rr));
            double sr = (cse_p565_t2_invr47718 * term.sigma);
            double sr2 = (sr * sr);
            double sr6 = (sr2 * sr2 * sr2);
            double sr12 = (sr6 * sr6);
            double fa_atr = (term.epsilon * (sr12 + (-2.00000000000000000e+0 * sr6)));
            double de_dr = (cse_p565_t2_invr47718 * term.epsilon * ((-1.20000000000000000e+1 * sr12) + (1.20000000000000000e+1 * sr6)));
            double cse_p11_t1_g110 = (rr * rr);
            double cse_p565_t3_invr47719 = (1.0 / (cse_p11_t1_g110));
            double d2e_dr2 = (cse_p565_t3_invr47719 * term.epsilon * ((-4.20000000000000040e+1 * sr6) + (1.56000000000000000e+2 * sr12)));
            double energy = fa_atr;
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p11_t1_g117 = (cse_p565_t2_invr47718 * de_dr);
            double g_x1 = (cse_p11_t1_g117 * dx);
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double g_y1 = (cse_p11_t1_g117 * dy);
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double g_z1 = (cse_p11_t1_g117 * dz);
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (cse_p11_t1_g117 * (-(dx)));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (cse_p11_t1_g117 * (-(dy)));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (cse_p11_t1_g117 * (-(dz)));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double cse_p11_t1_g108 = (cse_p565_t2_invr47718 * cse_p565_t2_invr47718 * d2e_dr2);
            double cse_p565_t1_invr47717 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g109 = (cse_p565_t1_invr47717 * de_dr);
            double cse_p12_t1_g219 = (cse_p11_t1_g109 * (cse_p11_t1_g110 + (-(cse_p11_t46_g56))));
            double cse_p60_t86_g5113 = (cse_p12_t1_g219 + (cse_p11_t1_g108 * dx * dx));
            double h_x1_x1 = cse_p60_t86_g5113;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = ((cse_p11_t1_g108 * dx * dy) + (cse_p11_t1_g109 * dy * (-(dx))));
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p11_t1_g113 = (cse_p11_t1_g109 * dz);
            double cse_p60_t37_g5064 = (cse_p11_t1_g113 * (-(dx)));
            double cse_p60_t38_g5065 = (cse_p11_t1_g113 * (-(dy)));
            double cse_p60_t39_g5066 = (cse_p11_t1_g113 * dx);
            double cse_p60_t40_g5067 = (cse_p11_t1_g113 * dy);
            double cse_p60_t65_g5092 = (cse_p60_t37_g5064 + (cse_p11_t1_g108 * dx * dz));
            double cse_p60_t66_g5093 = (cse_p60_t38_g5065 + (cse_p11_t1_g108 * dy * dz));
            double h_x1_z1 = cse_p60_t65_g5092;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p11_t1_g111 = (cse_p11_t1_g108 * (-(dx)));
            double cse_p60_t26_g5053 = (cse_p11_t1_g111 * (-(dx)));
            double cse_p60_t27_g5054 = (cse_p11_t1_g111 * (-(dy)));
            double cse_p60_t28_g5055 = (cse_p11_t1_g111 * (-(dz)));
            double cse_p60_t29_g5056 = (cse_p11_t1_g111 * dx);
            double cse_p60_t30_g5057 = (cse_p11_t1_g111 * dy);
            double cse_p60_t31_g5058 = (cse_p11_t1_g111 * dz);
            double cse_p60_t59_g5086 = (cse_p60_t29_g5056 + (cse_p11_t1_g109 * (cse_p11_t46_g56 + (-(cse_p11_t1_g110)))));
            double cse_p60_t67_g5094 = (cse_p60_t30_g5057 + (cse_p11_t1_g109 * dx * dy));
            double cse_p60_t69_g5096 = (cse_p60_t27_g5054 + (cse_p11_t1_g109 * dy * (-(dx))));
            double cse_p60_t70_g5097 = (cse_p60_t28_g5055 + cse_p60_t37_g5064);
            double cse_p60_t71_g5098 = (cse_p60_t31_g5058 + cse_p60_t39_g5066);
            double cse_p60_t87_g5114 = (cse_p12_t1_g219 + cse_p60_t26_g5053);
            double h_x1_x2 = cse_p60_t59_g5086;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p11_t1_g112 = (cse_p11_t1_g108 * (-(dy)));
            double cse_p60_t32_g5059 = (cse_p11_t1_g112 * (-(dy)));
            double cse_p60_t33_g5060 = (cse_p11_t1_g112 * (-(dz)));
            double cse_p60_t34_g5061 = (cse_p11_t1_g112 * dx);
            double cse_p60_t35_g5062 = (cse_p11_t1_g112 * dy);
            double cse_p60_t36_g5063 = (cse_p11_t1_g112 * dz);
            double cse_p60_t60_g5087 = (cse_p60_t35_g5062 + (cse_p11_t1_g109 * (cse_p11_t47_g57 + (-(cse_p11_t1_g110)))));
            double cse_p60_t68_g5095 = (cse_p60_t34_g5061 + (cse_p11_t1_g109 * dx * dy));
            double cse_p60_t72_g5099 = (cse_p60_t33_g5060 + cse_p60_t38_g5065);
            double cse_p60_t73_g5100 = (cse_p60_t36_g5063 + cse_p60_t40_g5067);
            double h_x1_y2 = cse_p60_t68_g5095;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p11_t1_g114 = (cse_p11_t1_g108 * (-(dz)));
            double cse_p60_t41_g5068 = (cse_p11_t1_g114 * (-(dz)));
            double cse_p60_t42_g5069 = (cse_p11_t1_g114 * dx);
            double cse_p60_t43_g5070 = (cse_p11_t1_g114 * dy);
            double cse_p60_t44_g5071 = (cse_p11_t1_g114 * dz);
            double cse_p60_t61_g5088 = (cse_p60_t44_g5071 + (cse_p11_t1_g109 * (cse_p11_t48_g58 + (-(cse_p11_t1_g110)))));
            double cse_p60_t74_g5101 = (cse_p60_t39_g5066 + cse_p60_t42_g5069);
            double cse_p60_t75_g5102 = (cse_p60_t40_g5067 + cse_p60_t43_g5070);
            double h_x1_z2 = cse_p60_t74_g5101;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p11_t1_g115 = (cse_p11_t1_g109 * (cse_p11_t1_g110 + (-(cse_p11_t47_g57))));
            double cse_p60_t82_g5109 = (cse_p11_t1_g115 + (cse_p11_t1_g108 * dy * dy));
            double cse_p60_t83_g5110 = (cse_p11_t1_g115 + cse_p60_t32_g5059);
            double h_y1_y1 = cse_p60_t82_g5109;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t66_g5093;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t67_g5094;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t60_g5087;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t75_g5102;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p11_t1_g116 = (cse_p11_t1_g109 * (cse_p11_t1_g110 + (-(cse_p11_t48_g58))));
            double cse_p60_t84_g5111 = (cse_p11_t1_g116 + (cse_p11_t1_g108 * dz * dz));
            double cse_p60_t85_g5112 = (cse_p11_t1_g116 + cse_p60_t41_g5068);
            double h_z1_z1 = cse_p60_t84_g5111;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t71_g5098;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t73_g5100;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t61_g5088;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t87_g5114;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t69_g5096;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t70_g5097;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t83_g5110;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t72_g5099;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t85_g5112;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p11_t1_g110 = (rr * rr);
              double fa_atr = (term.dd + (term.cc * rr) + (cse_p11_t1_g110 * (term.bb + (term.aa * rr))));
              double de_dr = (term.cc + (2.00000000000000000e+0 * term.bb * rr) + (3.00000000000000000e+0 * term.aa * cse_p11_t1_g110));
              double d2e_dr2 = ((2.00000000000000000e+0 * term.bb) + (6.00000000000000000e+0 * term.aa * rr));
              double energy = fa_atr;
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p565_t2_invr47718 = (1.0 / (rr));
              double cse_p11_t1_g117 = (cse_p565_t2_invr47718 * de_dr);
              double g_x1 = (cse_p11_t1_g117 * dx);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t1_g117 * dy);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t1_g117 * dz);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = (cse_p11_t1_g117 * (-(dx)));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (cse_p11_t1_g117 * (-(dy)));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (cse_p11_t1_g117 * (-(dz)));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double cse_p11_t1_g108 = (cse_p565_t2_invr47718 * cse_p565_t2_invr47718 * d2e_dr2);
              double cse_p565_t1_invr47717 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g109 = (cse_p565_t1_invr47717 * de_dr);
              double cse_p12_t1_g219 = (cse_p11_t1_g109 * (cse_p11_t1_g110 + (-(cse_p11_t46_g56))));
              double h_x1_x1 = (cse_p12_t1_g219 + (cse_p11_t1_g108 * dx * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double h_x1_y1 = ((cse_p11_t1_g108 * dx * dy) + (cse_p11_t1_g109 * dy * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g113 = (cse_p11_t1_g109 * dz);
              double h_x1_z1 = ((cse_p11_t1_g108 * dx * dz) + (cse_p11_t1_g113 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g111 = (cse_p11_t1_g108 * (-(dx)));
              double h_x1_x2 = ((cse_p11_t1_g109 * (cse_p11_t46_g56 + (-(cse_p11_t1_g110)))) + (cse_p11_t1_g111 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g112 = (cse_p11_t1_g108 * (-(dy)));
              double h_x1_y2 = ((cse_p11_t1_g109 * dx * dy) + (cse_p11_t1_g112 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g114 = (cse_p11_t1_g108 * (-(dz)));
              double h_x1_z2 = ((cse_p11_t1_g113 * dx) + (cse_p11_t1_g114 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g115 = (cse_p11_t1_g109 * (cse_p11_t1_g110 + (-(cse_p11_t47_g57))));
              double h_y1_y1 = (cse_p11_t1_g115 + (cse_p11_t1_g108 * dy * dy));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = ((cse_p11_t1_g108 * dy * dz) + (cse_p11_t1_g113 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = ((cse_p11_t1_g109 * dx * dy) + (cse_p11_t1_g111 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = ((cse_p11_t1_g109 * (cse_p11_t47_g57 + (-(cse_p11_t1_g110)))) + (cse_p11_t1_g112 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = ((cse_p11_t1_g113 * dy) + (cse_p11_t1_g114 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g116 = (cse_p11_t1_g109 * (cse_p11_t1_g110 + (-(cse_p11_t48_g58))));
              double h_z1_z1 = (cse_p11_t1_g116 + (cse_p11_t1_g108 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = ((cse_p11_t1_g111 * dz) + (cse_p11_t1_g113 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = ((cse_p11_t1_g112 * dz) + (cse_p11_t1_g113 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = ((cse_p11_t1_g109 * (cse_p11_t48_g58 + (-(cse_p11_t1_g110)))) + (cse_p11_t1_g114 * dz));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p12_t1_g219 + (cse_p11_t1_g111 * (-(dx))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = ((cse_p11_t1_g109 * dy * (-(dx))) + (cse_p11_t1_g111 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = ((cse_p11_t1_g111 * (-(dz))) + (cse_p11_t1_g113 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = (cse_p11_t1_g115 + (cse_p11_t1_g112 * (-(dy))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = ((cse_p11_t1_g112 * (-(dz))) + (cse_p11_t1_g113 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = (cse_p11_t1_g116 + (cse_p11_t1_g114 * (-(dz))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
            } else {
              /* !COND-ELSE-4 */
              {
                /* !COND-ELSE-5 */
                double de_dr = 0.00000000000000000e+0;
                double d2e_dr2 = 0.00000000000000000e+0;
                double energy = 0.00000000000000000e+0;
                energy_added += energy;
                *energy_accumulate += energy;
                double g_x1 = (de_dr * dx * (1.0 / (rr)));
                KernelGradientAcc(term.i3x1, 0, g_x1);
                double cse_p565_t2_invr47718 = (1.0 / (rr));
                double g_y1 = (cse_p565_t2_invr47718 * de_dr * dy);
                KernelGradientAcc(term.i3x1, 1, g_y1);
                double cse_p11_t1_g117 = (cse_p565_t2_invr47718 * de_dr);
                double g_z1 = (cse_p11_t1_g117 * dz);
                KernelGradientAcc(term.i3x1, 2, g_z1);
                double g_x2 = (cse_p11_t1_g117 * (-(dx)));
                KernelGradientAcc(term.i3x2, 0, g_x2);
                double g_y2 = (cse_p11_t1_g117 * (-(dy)));
                KernelGradientAcc(term.i3x2, 1, g_y2);
                double g_z2 = (cse_p11_t1_g117 * (-(dz)));
                KernelGradientAcc(term.i3x2, 2, g_z2);
                double h_x1_x1 = ((cse_p565_t2_invr47718 * cse_p565_t2_invr47718 * d2e_dr2 * dx * dx) + (de_dr * ((rr * rr) + (-(cse_p11_t46_g56))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t1_invr47717 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g110 = (rr * rr);
                double cse_p11_t1_g109 = (cse_p565_t1_invr47717 * de_dr);
                double cse_p11_t1_g108 = (cse_p565_t2_invr47718 * cse_p565_t2_invr47718 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g108 * dx * dy) + (cse_p11_t1_g109 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double h_x1_z1 = ((cse_p11_t1_g108 * dx * dz) + (cse_p11_t1_g109 * dz * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g113 = (cse_p11_t1_g109 * dz);
                double h_x1_x2 = ((cse_p11_t1_g108 * dx * (-(dx))) + (cse_p11_t1_g109 * (cse_p11_t46_g56 + (-(cse_p11_t1_g110)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g111 = (cse_p11_t1_g108 * (-(dx)));
                double h_x1_y2 = ((cse_p11_t1_g108 * dx * (-(dy))) + (cse_p11_t1_g109 * dx * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g112 = (cse_p11_t1_g108 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g108 * dx * (-(dz))) + (cse_p11_t1_g113 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g114 = (cse_p11_t1_g108 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g108 * dy * dy) + (cse_p11_t1_g109 * (cse_p11_t1_g110 + (-(cse_p11_t47_g57)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g115 = (cse_p11_t1_g109 * (cse_p11_t1_g110 + (-(cse_p11_t47_g57))));
                double h_y1_z1 = ((cse_p11_t1_g108 * dy * dz) + (cse_p11_t1_g113 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g109 * dx * dy) + (cse_p11_t1_g111 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g109 * (cse_p11_t47_g57 + (-(cse_p11_t1_g110)))) + (cse_p11_t1_g112 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g113 * dy) + (cse_p11_t1_g114 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g108 * dz * dz) + (cse_p11_t1_g109 * (cse_p11_t1_g110 + (-(cse_p11_t48_g58)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g116 = (cse_p11_t1_g109 * (cse_p11_t1_g110 + (-(cse_p11_t48_g58))));
                double h_z1_x2 = ((cse_p11_t1_g111 * dz) + (cse_p11_t1_g113 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g112 * dz) + (cse_p11_t1_g113 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g109 * (cse_p11_t48_g58 + (-(cse_p11_t1_g110)))) + (cse_p11_t1_g114 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g109 * (cse_p11_t1_g110 + (-(cse_p11_t46_g56)))) + (cse_p11_t1_g111 * (-(dx))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p12_t1_g219 = (cse_p11_t1_g109 * (cse_p11_t1_g110 + (-(cse_p11_t46_g56))));
                double h_x2_y2 = ((cse_p11_t1_g109 * dy * (-(dx))) + (cse_p11_t1_g111 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double h_x2_z2 = ((cse_p11_t1_g111 * (-(dz))) + (cse_p11_t1_g113 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g115 + (cse_p11_t1_g112 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g112 * (-(dz))) + (cse_p11_t1_g113 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g116 + (cse_p11_t1_g114 * (-(dz))));
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
void hessian_fd(const rosetta_nonbond_parameters& params, const rosetta_nonbond_term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
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
