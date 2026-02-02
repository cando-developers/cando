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
      double cse_p1_t1_sqrt52257 = sqrt(r2);
      double cse_p1_t2_invsqrt52258 = (1.0 / (cse_p1_t1_sqrt52257));
      double rr = cse_p1_t1_sqrt52257;
      if ((rr <= (6.00000023841857900e-1 * term.sigma))) {
        /* !COND-1 */
        double cse_p528_t7_g65798 = (cse_p1_t2_invsqrt52258 * dx);
        double cse_p528_t8_g65799 = (cse_p1_t2_invsqrt52258 * dy);
        double cse_p528_t9_g65800 = (cse_p1_t2_invsqrt52258 * dz);
        double cse_p528_t16_g65807 = (-(term.epsilon));
        double cse_p12_t1_g52325 = (term.mij * params.rep_weight);
        double cse_p11_t9_g52269 = (cse_p12_t1_g52325 * cse_p528_t7_g65798);
        double cse_p11_t10_g52270 = (cse_p12_t1_g52325 * cse_p528_t8_g65799);
        double cse_p11_t11_g52271 = (cse_p12_t1_g52325 * cse_p528_t9_g65800);
        double fa_rep = (term.bij + (term.mij * rr));
        double cse_p528_t11_g65802 = (fa_rep * params.rep_weight);
        double fa_atr = cse_p528_t16_g65807;
        double cse_p13_t20_g52351 = (cse_p528_t11_g65802 + fa_atr);
        double de_dr = cse_p12_t1_g52325;
        double energy = cse_p13_t20_g52351;
        energy_added += energy;
        *energy_accumulate += energy;
        double g_x1 = cse_p11_t9_g52269;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = cse_p11_t10_g52270;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = cse_p11_t11_g52271;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = (-(cse_p11_t9_g52269));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (-(cse_p11_t10_g52270));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (-(cse_p11_t11_g52271));
        KernelGradientAcc(term.i3x2, 2, g_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.sigma)) {
          /* !COND-2 */
          double cse_p529_t1_invr65819 = (1.0 / (rr));
          double cse_p479_t11_g64477 = (cse_p529_t1_invr65819 * cse_p529_t1_invr65819);
          double cse_p479_t12_g64478 = (cse_p529_t1_invr65819 * term.sigma);
          double sr = cse_p479_t12_g64478;
          double cse_p528_t1_g65818 = (cse_p1_t2_invsqrt52258 * term.sigma * sr);
          double cse_p528_t12_g65803 = (sr * sr);
          double cse_p11_t1_g52289 = (cse_p479_t11_g64477 * cse_p528_t1_g65818);
          double cse_p528_t4_g65795 = (cse_p11_t1_g52289 * dx);
          double cse_p528_t5_g65796 = (cse_p11_t1_g52289 * dy);
          double cse_p528_t6_g65797 = (cse_p11_t1_g52289 * dz);
          double cse_p12_t19_g52317 = (-(cse_p528_t4_g65795));
          double cse_p12_t20_g52318 = (-(cse_p528_t5_g65796));
          double cse_p12_t21_g52319 = (-(cse_p528_t6_g65797));
          double sr2 = cse_p528_t12_g65803;
          double cse_p528_t13_g65804 = (sr2 * sr2);
          double cse_p13_t4_g52335 = (6.00000000000000000e+0 * cse_p528_t13_g65804 * term.epsilon);
          double cse_p13_t8_g52339 = (cse_p528_t13_g65804 * sr2);
          double sr6 = cse_p13_t8_g52339;
          double cse_p528_t2_g65793 = (-2.00000000000000000e+0 * sr6);
          double cse_p528_t3_g65794 = (1.20000000000000000e+1 * sr6);
          double cse_p528_t14_g65805 = (sr6 * sr6);
          double cse_p528_t15_g65806 = (-2.00000000000000000e+0 + sr6 + sr6);
          double sr12 = cse_p528_t14_g65805;
          double cse_p528_t1_g65792 = (-1.20000000000000000e+1 * sr12);
          double cse_p14_t15_g52378 = (cse_p528_t2_g65793 + sr12);
          double cse_p13_t1_g52363 = (cse_p14_t15_g52378 * term.epsilon);
          double cse_p13_t19_g52350 = (cse_p528_t1_g65792 + cse_p528_t3_g65794);
          double fa_rep_epsilon = cse_p13_t1_g52363;
          double fa_rep = (term.epsilon + fa_rep_epsilon);
          double fa_atr = (-(term.epsilon));
          double de_dr = (cse_p13_t19_g52350 * cse_p529_t1_invr65819 * term.epsilon * params.rep_weight);
          double cse_p12_t1_g52329 = (fa_rep * params.rep_weight);
          double energy = (cse_p12_t1_g52329 + fa_atr);
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p11_t1_g52292 = (cse_p13_t4_g52335 * cse_p528_t15_g65806 * params.rep_weight);
          double g_x1 = (cse_p11_t1_g52292 * cse_p12_t19_g52317);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t1_g52292 * cse_p12_t20_g52318);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t1_g52292 * cse_p12_t21_g52319);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p11_t1_g52292 * cse_p528_t4_g65795);
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p11_t1_g52292 * cse_p528_t5_g65796);
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p11_t1_g52292 * cse_p528_t6_g65797);
          KernelGradientAcc(term.i3x2, 2, g_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.r_switch)) {
            /* !COND-3 */
            double cse_p529_t1_invr65819 = (1.0 / (rr));
            double cse_p1_t4_invr252260 = (cse_p529_t1_invr65819 * cse_p529_t1_invr65819);
            double sr = (cse_p529_t1_invr65819 * term.sigma);
            double sr2 = (sr * sr);
            double cse_p12_t1_g52326 = (sr2 * sr2);
            double sr6 = (cse_p12_t1_g52326 * sr2);
            double sr12 = (sr6 * sr6);
            double cse_p12_t1_g52324 = (-2.00000000000000000e+0 * sr6);
            double cse_p13_t1_g52363 = (term.epsilon * (cse_p12_t1_g52324 + sr12));
            double fa_atr = cse_p13_t1_g52363;
            double cse_p12_t1_g52327 = (-1.20000000000000000e+1 * sr12);
            double cse_p12_t1_g52328 = (1.20000000000000000e+1 * sr6);
            double de_dr = (cse_p529_t1_invr65819 * term.epsilon * (cse_p12_t1_g52327 + cse_p12_t1_g52328));
            double energy = fa_atr;
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p528_t1_g65818 = (cse_p1_t2_invsqrt52258 * term.sigma * sr);
            double cse_p11_t1_g52289 = (cse_p1_t4_invr252260 * cse_p528_t1_g65818);
            double cse_p11_t1_g52290 = (6.00000000000000000e+0 * cse_p12_t1_g52326 * term.epsilon);
            double cse_p11_t1_g52291 = (cse_p11_t1_g52290 * (-2.00000000000000000e+0 + sr6 + sr6));
            double cse_p11_t1_g52295 = (cse_p11_t1_g52289 * dx);
            double g_x1 = (cse_p11_t1_g52291 * (-(cse_p11_t1_g52295)));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p11_t1_g52296 = (cse_p11_t1_g52289 * dy);
            double g_y1 = (cse_p11_t1_g52291 * (-(cse_p11_t1_g52296)));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p11_t1_g52297 = (cse_p11_t1_g52289 * dz);
            double g_z1 = (cse_p11_t1_g52291 * (-(cse_p11_t1_g52297)));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (cse_p11_t1_g52291 * cse_p11_t1_g52295);
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (cse_p11_t1_g52291 * cse_p11_t1_g52296);
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (cse_p11_t1_g52291 * cse_p11_t1_g52297);
            KernelGradientAcc(term.i3x2, 2, g_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.r_cut)) {
              /* !COND-4 */
              double cse_p11_t1_g52294 = (term.aa * rr * rr);
              double cse_p12_t4_g52302 = (3.00000000000000000e+0 * cse_p11_t1_g52294);
              double cse_p11_t1_g52293 = (term.bb * rr);
              double cse_p11_t21_g52281 = (term.cc + cse_p11_t1_g52293 + cse_p11_t1_g52293 + cse_p12_t4_g52302);
              double fa_atr = (term.dd + (term.cc * rr) + (rr * rr * (term.bb + (term.aa * rr))));
              double de_dr = (term.cc + cse_p12_t4_g52302 + (2.00000000000000000e+0 * cse_p11_t1_g52293));
              double energy = fa_atr;
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p11_t1_g52298 = (cse_p1_t2_invsqrt52258 * dx);
              double g_x1 = (cse_p11_t1_g52298 * cse_p11_t21_g52281);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double cse_p12_t1_g52322 = (cse_p1_t2_invsqrt52258 * dy);
              double g_y1 = (cse_p11_t21_g52281 * cse_p12_t1_g52322);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double cse_p12_t1_g52323 = (cse_p1_t2_invsqrt52258 * dz);
              double g_z1 = (cse_p11_t21_g52281 * cse_p12_t1_g52323);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p63_t1_g53654 = (-2.00000000000000000e+0 * cse_p11_t1_g52293);
              double cse_p63_t1_g53655 = (-3.00000000000000000e+0 * cse_p11_t1_g52294);
              double g_x2 = ((cse_p11_t1_g52298 * cse_p63_t1_g53654) + (cse_p11_t1_g52298 * cse_p63_t1_g53655) + (-((term.cc * cse_p11_t1_g52298))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = ((cse_p12_t1_g52322 * cse_p63_t1_g53654) + (cse_p12_t1_g52322 * cse_p63_t1_g53655) + (-((term.cc * cse_p12_t1_g52322))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = ((cse_p12_t1_g52323 * cse_p63_t1_g53654) + (cse_p12_t1_g52323 * cse_p63_t1_g53655) + (-((term.cc * cse_p12_t1_g52323))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
            } else {
              /* !COND-ELSE-4 */
              {
                /* !COND-ELSE-5 */
                double de_dr = 0.00000000000000000e+0;
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
          if ((rr <= params.r_switch)) {
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
            if ((rr <= params.r_cut)) {
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
      double cse_p11_t49_g74953 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t50_g74954 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t51_g74955 = (dz * dz);
      double r2 = (cse_p11_t49_g74953 + cse_p11_t50_g74954 + cse_p11_t51_g74955);
      double cse_p1_t1_sqrt74901 = sqrt(r2);
      double cse_p1_t2_invsqrt74902 = (1.0 / (cse_p1_t1_sqrt74901));
      double rr = cse_p1_t1_sqrt74901;
      if ((rr <= (6.00000023841857900e-1 * term.sigma))) {
        /* !COND-1 */
        double cse_p564_t46_g115508 = (rr * rr);
        double cse_p564_t58_g115520 = (-(cse_p11_t49_g74953));
        double cse_p564_t59_g115521 = (-(cse_p11_t50_g74954));
        double cse_p564_t60_g115522 = (-(cse_p11_t51_g74955));
        double cse_p564_t61_g115523 = (-(dx));
        double cse_p564_t62_g115524 = (-(dy));
        double cse_p564_t63_g115525 = (-(dz));
        double cse_p564_t64_g115526 = (-(term.epsilon));
        double cse_p565_t1_invr115541 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr115542 = (1.0 / (rr));
        double cse_p564_t68_g115530 = (rr * rr * rr);
        double cse_p15_t89_g75461 = (-(cse_p564_t46_g115508));
        double cse_p14_t90_g75343 = (cse_p11_t49_g74953 + cse_p15_t89_g75461);
        double cse_p14_t92_g75345 = (cse_p11_t50_g74954 + cse_p15_t89_g75461);
        double cse_p14_t94_g75347 = (cse_p11_t51_g74955 + cse_p15_t89_g75461);
        double cse_p13_t116_g75236 = (cse_p564_t46_g115508 + cse_p564_t58_g115520);
        double cse_p13_t117_g75237 = (cse_p564_t46_g115508 + cse_p564_t59_g115521);
        double cse_p13_t118_g75238 = (cse_p564_t46_g115508 + cse_p564_t60_g115522);
        double cse_p11_t11_g74915 = (cse_p1_t2_invsqrt74902 * dx * term.mij * params.rep_weight);
        double cse_p11_t12_g74916 = (cse_p1_t2_invsqrt74902 * dy * term.mij * params.rep_weight);
        double cse_p11_t13_g74917 = (cse_p1_t2_invsqrt74902 * dz * term.mij * params.rep_weight);
        double fa_rep = (term.bij + (term.mij * rr));
        double cse_p564_t45_g115507 = (fa_rep * params.rep_weight);
        double fa_atr = cse_p564_t64_g115526;
        double cse_p60_t83_g79776 = (cse_p564_t45_g115507 + fa_atr);
        double de_dr = (term.mij * params.rep_weight);
        double cse_p515_t36_g111676 = (cse_p565_t1_invr115541 * de_dr);
        double cse_p13_t68_g75188 = (cse_p14_t90_g75343 * cse_p515_t36_g111676);
        double cse_p13_t69_g75189 = (cse_p14_t92_g75345 * cse_p515_t36_g111676);
        double cse_p13_t70_g75190 = (cse_p14_t94_g75347 * cse_p515_t36_g111676);
        double cse_p13_t71_g75191 = (cse_p515_t36_g111676 * dx * dy);
        double cse_p13_t72_g75192 = (cse_p515_t36_g111676 * cse_p564_t61_g115523 * dy);
        double cse_p12_t23_g75042 = (cse_p13_t116_g75236 * cse_p515_t36_g111676);
        double cse_p12_t24_g75043 = (cse_p13_t117_g75237 * cse_p515_t36_g111676);
        double cse_p12_t25_g75044 = (cse_p13_t118_g75238 * cse_p515_t36_g111676);
        double cse_p12_t28_g75047 = (cse_p515_t36_g111676 * dz);
        double cse_p13_t58_g75178 = (cse_p12_t28_g75047 * cse_p564_t61_g115523);
        double cse_p13_t59_g75179 = (cse_p12_t28_g75047 * cse_p564_t62_g115524);
        double cse_p13_t60_g75180 = (cse_p12_t28_g75047 * dx);
        double cse_p13_t61_g75181 = (cse_p12_t28_g75047 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p515_t37_g111677 = (cse_p565_t2_invr115542 * cse_p565_t2_invr115542 * d2e_dr2);
        double cse_p13_t62_g75182 = (cse_p515_t37_g111677 * dx * dx);
        double cse_p13_t63_g75183 = (cse_p515_t37_g111677 * dx * dy);
        double cse_p13_t64_g75184 = (cse_p515_t37_g111677 * dx * dz);
        double cse_p13_t65_g75185 = (cse_p515_t37_g111677 * dy * dy);
        double cse_p13_t66_g75186 = (cse_p515_t37_g111677 * dy * dz);
        double cse_p13_t67_g75187 = (cse_p515_t37_g111677 * dz * dz);
        double cse_p12_t11_g75030 = (cse_p515_t37_g111677 * cse_p564_t61_g115523);
        double cse_p13_t43_g75163 = (cse_p12_t11_g75030 * cse_p564_t61_g115523);
        double cse_p13_t44_g75164 = (cse_p12_t11_g75030 * cse_p564_t62_g115524);
        double cse_p13_t45_g75165 = (cse_p12_t11_g75030 * cse_p564_t63_g115525);
        double cse_p13_t46_g75166 = (cse_p12_t11_g75030 * dx);
        double cse_p13_t47_g75167 = (cse_p12_t11_g75030 * dy);
        double cse_p13_t48_g75168 = (cse_p12_t11_g75030 * dz);
        double cse_p12_t12_g75031 = (cse_p515_t37_g111677 * cse_p564_t62_g115524);
        double cse_p13_t49_g75169 = (cse_p12_t12_g75031 * cse_p564_t62_g115524);
        double cse_p13_t50_g75170 = (cse_p12_t12_g75031 * cse_p564_t63_g115525);
        double cse_p13_t51_g75171 = (cse_p12_t12_g75031 * dx);
        double cse_p13_t52_g75172 = (cse_p12_t12_g75031 * dy);
        double cse_p13_t53_g75173 = (cse_p12_t12_g75031 * dz);
        double cse_p12_t13_g75032 = (cse_p515_t37_g111677 * cse_p564_t63_g115525);
        double cse_p13_t54_g75174 = (cse_p12_t13_g75032 * cse_p564_t63_g115525);
        double cse_p13_t55_g75175 = (cse_p12_t13_g75032 * dx);
        double cse_p13_t56_g75176 = (cse_p12_t13_g75032 * dy);
        double cse_p13_t57_g75177 = (cse_p12_t13_g75032 * dz);
        double cse_p12_t63_g75082 = (cse_p13_t63_g75183 + cse_p13_t72_g75192);
        double cse_p12_t64_g75083 = (cse_p13_t58_g75178 + cse_p13_t64_g75184);
        double cse_p12_t65_g75084 = (cse_p13_t59_g75179 + cse_p13_t66_g75186);
        double cse_p12_t82_g75101 = (cse_p12_t24_g75043 + cse_p13_t65_g75185);
        double cse_p12_t84_g75103 = (cse_p12_t25_g75044 + cse_p13_t67_g75187);
        double cse_p12_t86_g75105 = (cse_p12_t23_g75042 + cse_p13_t62_g75182);
        double cse_p12_t58_g75077 = (cse_p13_t46_g75166 + cse_p13_t68_g75188);
        double cse_p12_t66_g75085 = (cse_p13_t47_g75167 + cse_p13_t71_g75191);
        double cse_p12_t68_g75087 = (cse_p13_t44_g75164 + cse_p13_t72_g75192);
        double cse_p12_t69_g75088 = (cse_p13_t45_g75165 + cse_p13_t58_g75178);
        double cse_p12_t70_g75089 = (cse_p13_t48_g75168 + cse_p13_t60_g75180);
        double cse_p12_t87_g75106 = (cse_p12_t23_g75042 + cse_p13_t43_g75163);
        double cse_p12_t59_g75078 = (cse_p13_t52_g75172 + cse_p13_t69_g75189);
        double cse_p12_t67_g75086 = (cse_p13_t51_g75171 + cse_p13_t71_g75191);
        double cse_p12_t71_g75090 = (cse_p13_t50_g75170 + cse_p13_t59_g75179);
        double cse_p12_t72_g75091 = (cse_p13_t53_g75173 + cse_p13_t61_g75181);
        double cse_p12_t83_g75102 = (cse_p12_t24_g75043 + cse_p13_t49_g75169);
        double cse_p12_t60_g75079 = (cse_p13_t57_g75177 + cse_p13_t70_g75190);
        double cse_p12_t73_g75092 = (cse_p13_t55_g75175 + cse_p13_t60_g75180);
        double cse_p12_t74_g75093 = (cse_p13_t56_g75176 + cse_p13_t61_g75181);
        double cse_p12_t85_g75104 = (cse_p12_t25_g75044 + cse_p13_t54_g75174);
        double energy = cse_p60_t83_g79776;
        energy_added += energy;
        *energy_accumulate += energy;
        double g_x1 = cse_p11_t11_g74915;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = cse_p11_t12_g74916;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = cse_p11_t13_g74917;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = (-(cse_p11_t11_g74915));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (-(cse_p11_t12_g74916));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (-(cse_p11_t13_g74917));
        KernelGradientAcc(term.i3x2, 2, g_z2);
        double h_x1_x1 = cse_p12_t86_g75105;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t63_g75082;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t64_g75083;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t58_g75077;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t67_g75086;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t73_g75092;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t82_g75101;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t65_g75084;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t66_g75085;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t59_g75078;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t74_g75093;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t84_g75103;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t70_g75089;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t72_g75091;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t60_g75079;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t87_g75106;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t68_g75087;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t69_g75088;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t83_g75102;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t71_g75090;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t85_g75104;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.sigma)) {
          /* !COND-2 */
          double cse_p11_t1_g75012 = (rr * rr);
          double cse_p564_t51_g115513 = (cse_p11_t1_g75012 + (-(cse_p11_t49_g74953)));
          double cse_p564_t52_g115514 = (cse_p11_t1_g75012 + (-(cse_p11_t50_g74954)));
          double cse_p564_t53_g115515 = (cse_p11_t1_g75012 + (-(cse_p11_t51_g74955)));
          double cse_p564_t54_g115516 = (cse_p11_t49_g74953 + (-(cse_p11_t1_g75012)));
          double cse_p564_t55_g115517 = (cse_p11_t50_g74954 + (-(cse_p11_t1_g75012)));
          double cse_p564_t56_g115518 = (cse_p11_t51_g74955 + (-(cse_p11_t1_g75012)));
          double cse_p564_t57_g115519 = (-(cse_p11_t1_g75012));
          double cse_p565_t3_invr115543 = (1.0 / (cse_p11_t1_g75012));
          double cse_p565_t2_invr115542 = (1.0 / (rr));
          double cse_p564_t38_g115500 = (cse_p565_t2_invr115542 * cse_p565_t2_invr115542);
          double cse_p564_t39_g115501 = (cse_p565_t2_invr115542 * term.sigma);
          double sr = cse_p564_t39_g115501;
          double cse_p564_t47_g115509 = (sr * sr);
          double cse_p67_t1_g79894 = (cse_p1_t2_invsqrt74902 * cse_p564_t38_g115500 * term.sigma * sr);
          double cse_p12_t8_g75027 = (cse_p67_t1_g79894 * dx);
          double cse_p12_t9_g75028 = (cse_p67_t1_g79894 * dy);
          double cse_p12_t10_g75029 = (cse_p67_t1_g79894 * dz);
          double cse_p11_t91_g74995 = (-(cse_p12_t8_g75027));
          double cse_p11_t92_g74996 = (-(cse_p12_t9_g75028));
          double cse_p11_t93_g74997 = (-(cse_p12_t10_g75029));
          double sr2 = cse_p564_t47_g115509;
          double cse_p564_t48_g115510 = (sr2 * sr2 * sr2);
          double sr6 = cse_p564_t48_g115510;
          double cse_p564_t2_g115464 = (-2.00000000000000000e+0 * sr6);
          double cse_p564_t3_g115465 = (-4.20000000000000040e+1 * sr6);
          double cse_p564_t4_g115466 = (1.20000000000000000e+1 * sr6);
          double cse_p564_t49_g115511 = (sr6 * sr6);
          double cse_p564_t50_g115512 = (-2.00000000000000000e+0 + sr6 + sr6);
          double sr12 = cse_p564_t49_g115511;
          double cse_p564_t1_g115463 = (-1.20000000000000000e+1 * sr12);
          double cse_p564_t5_g115467 = (1.56000000000000000e+2 * sr12);
          double cse_p67_t78_g79875 = (cse_p564_t2_g115464 + sr12);
          double cse_p60_t47_g79740 = (cse_p67_t78_g79875 * term.epsilon);
          double cse_p60_t56_g79749 = (cse_p564_t1_g115463 + cse_p564_t4_g115466);
          double cse_p60_t57_g79750 = (cse_p564_t3_g115465 + cse_p564_t5_g115467);
          double fa_rep_epsilon = cse_p60_t47_g79740;
          double fa_rep = (term.epsilon + fa_rep_epsilon);
          double fa_atr = (-(term.epsilon));
          double de_dr = (cse_p565_t2_invr115542 * cse_p60_t56_g79749 * term.epsilon * params.rep_weight);
          double d2e_dr2 = (cse_p565_t3_invr115543 * cse_p60_t57_g79750 * term.epsilon * params.rep_weight);
          double cse_p564_t37_g115499 = (cse_p565_t2_invr115542 * cse_p565_t2_invr115542 * d2e_dr2);
          double energy = (fa_atr + (fa_rep * params.rep_weight));
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p67_t1_g79891 = (6.00000000000000000e+0 * cse_p564_t50_g115512 * term.epsilon * params.rep_weight * sr2 * sr2);
          double g_x1 = (cse_p11_t91_g74995 * cse_p67_t1_g79891);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t92_g74996 * cse_p67_t1_g79891);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t93_g74997 * cse_p67_t1_g79891);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p67_t1_g79891 * cse_p67_t1_g79894 * dx);
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p67_t1_g79891 * cse_p67_t1_g79894 * dy);
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p67_t1_g79891 * cse_p67_t1_g79894 * dz);
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double cse_p68_t8_g79905 = (cse_p564_t37_g115499 * (-(dx)));
          double cse_p68_t9_g79906 = (cse_p564_t37_g115499 * (-(dy)));
          double cse_p68_t10_g79907 = (cse_p564_t37_g115499 * (-(dz)));
          double cse_p68_t11_g79908 = (cse_p564_t37_g115499 * dx);
          double cse_p68_t12_g79909 = (cse_p564_t37_g115499 * dy);
          double cse_p68_t13_g79910 = (cse_p564_t37_g115499 * dz * dz);
          double cse_p68_t45_g79942 = (cse_p68_t11_g79908 * dx);
          double cse_p68_t46_g79943 = (cse_p68_t11_g79908 * dy);
          double cse_p68_t47_g79944 = (cse_p68_t11_g79908 * dz);
          double cse_p68_t50_g79947 = (cse_p68_t12_g79909 * dy);
          double cse_p68_t51_g79948 = (cse_p68_t12_g79909 * dz);
          double cse_p14_t47_g75300 = (cse_p68_t8_g79905 * (-(dx)));
          double cse_p14_t48_g75301 = (cse_p68_t8_g79905 * (-(dy)));
          double cse_p14_t49_g75302 = (cse_p68_t8_g79905 * (-(dz)));
          double cse_p14_t50_g75303 = (cse_p68_t8_g79905 * dx);
          double cse_p14_t51_g75304 = (cse_p68_t8_g79905 * dy);
          double cse_p14_t52_g75305 = (cse_p68_t8_g79905 * dz);
          double cse_p14_t53_g75306 = (cse_p68_t9_g79906 * (-(dy)));
          double cse_p14_t54_g75307 = (cse_p68_t9_g79906 * (-(dz)));
          double cse_p14_t55_g75308 = (cse_p68_t9_g79906 * dx);
          double cse_p14_t56_g75309 = (cse_p68_t9_g79906 * dy);
          double cse_p14_t57_g75310 = (cse_p68_t9_g79906 * dz);
          double cse_p14_t58_g75311 = (cse_p68_t10_g79907 * (-(dz)));
          double cse_p14_t59_g75312 = (cse_p68_t10_g79907 * dx);
          double cse_p14_t60_g75313 = (cse_p68_t10_g79907 * dy);
          double cse_p14_t61_g75314 = (cse_p68_t10_g79907 * dz);
          double cse_p565_t1_invr115541 = (1.0 / ((rr * rr * rr)));
          double cse_p564_t36_g115498 = (cse_p565_t1_invr115541 * de_dr);
          double cse_p68_t20_g79917 = (cse_p564_t36_g115498 * dy);
          double cse_p68_t21_g79918 = (cse_p564_t36_g115498 * dz);
          double cse_p68_t48_g79945 = (cse_p68_t20_g79917 * (-(dx)));
          double cse_p68_t49_g79946 = (cse_p68_t20_g79917 * dx);
          double cse_p60_t15_g79708 = (cse_p564_t36_g115498 * (cse_p11_t49_g74953 + cse_p564_t57_g115519));
          double cse_p60_t16_g79709 = (cse_p564_t36_g115498 * (cse_p11_t50_g74954 + cse_p564_t57_g115519));
          double cse_p60_t17_g79710 = (cse_p564_t36_g115498 * (cse_p11_t51_g74955 + cse_p564_t57_g115519));
          double cse_p60_t18_g79711 = (cse_p564_t36_g115498 * cse_p564_t51_g115513);
          double cse_p60_t19_g79712 = (cse_p564_t36_g115498 * cse_p564_t52_g115514);
          double cse_p60_t20_g79713 = (cse_p564_t36_g115498 * cse_p564_t53_g115515);
          double cse_p60_t58_g79751 = (cse_p68_t46_g79943 + cse_p68_t48_g79945);
          double cse_p14_t43_g75296 = (cse_p68_t21_g79918 * (-(dx)));
          double cse_p14_t44_g75297 = (cse_p68_t21_g79918 * (-(dy)));
          double cse_p14_t45_g75298 = (cse_p68_t21_g79918 * dx);
          double cse_p14_t46_g75299 = (cse_p68_t21_g79918 * dy);
          double cse_p13_t114_g75234 = (cse_p60_t18_g79711 + cse_p68_t45_g79942);
          double h_x1_x1 = cse_p13_t114_g75234;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t58_g79751;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t89_g75209 = (cse_p14_t43_g75296 + cse_p68_t47_g79944);
          double cse_p13_t90_g75210 = (cse_p14_t44_g75297 + cse_p68_t51_g79948);
          double h_x1_z1 = cse_p13_t89_g75209;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t83_g75203 = (cse_p14_t50_g75303 + cse_p60_t15_g79708);
          double cse_p13_t91_g75211 = (cse_p14_t51_g75304 + cse_p68_t49_g79946);
          double cse_p13_t93_g75213 = (cse_p14_t48_g75301 + cse_p68_t48_g79945);
          double cse_p13_t94_g75214 = (cse_p14_t43_g75296 + cse_p14_t49_g75302);
          double cse_p13_t95_g75215 = (cse_p14_t45_g75298 + cse_p14_t52_g75305);
          double cse_p13_t115_g75235 = (cse_p14_t47_g75300 + cse_p60_t18_g79711);
          double h_x1_x2 = cse_p13_t83_g75203;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t84_g75204 = (cse_p14_t56_g75309 + cse_p60_t16_g79709);
          double cse_p13_t92_g75212 = (cse_p14_t55_g75308 + cse_p68_t49_g79946);
          double cse_p13_t96_g75216 = (cse_p14_t44_g75297 + cse_p14_t54_g75307);
          double cse_p13_t97_g75217 = (cse_p14_t46_g75299 + cse_p14_t57_g75310);
          double h_x1_y2 = cse_p13_t92_g75212;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t85_g75205 = (cse_p14_t61_g75314 + cse_p60_t17_g79710);
          double cse_p13_t98_g75218 = (cse_p14_t45_g75298 + cse_p14_t59_g75312);
          double cse_p13_t99_g75219 = (cse_p14_t46_g75299 + cse_p14_t60_g75313);
          double h_x1_z2 = cse_p13_t98_g75218;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t110_g75230 = (cse_p60_t19_g79712 + cse_p68_t50_g79947);
          double cse_p13_t111_g75231 = (cse_p14_t53_g75306 + cse_p60_t19_g79712);
          double h_y1_y1 = cse_p13_t110_g75230;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t90_g75210;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t91_g75211;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t84_g75204;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t99_g75219;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t112_g75232 = (cse_p60_t20_g79713 + cse_p68_t13_g79910);
          double cse_p13_t113_g75233 = (cse_p14_t58_g75311 + cse_p60_t20_g79713);
          double h_z1_z1 = cse_p13_t112_g75232;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t95_g75215;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t97_g75217;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t85_g75205;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t115_g75235;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t93_g75213;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t94_g75214;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t111_g75231;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t96_g75216;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t113_g75233;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.r_switch)) {
            /* !COND-3 */
            double cse_p565_t2_invr115542 = (1.0 / (rr));
            double cse_p1_t4_invr274904 = (cse_p565_t2_invr115542 * cse_p565_t2_invr115542);
            double sr = (cse_p565_t2_invr115542 * term.sigma);
            double sr2 = (sr * sr);
            double sr6 = (sr2 * sr2 * sr2);
            double sr12 = (sr6 * sr6);
            double fa_atr = (term.epsilon * (sr12 + (-2.00000000000000000e+0 * sr6)));
            double de_dr = (cse_p565_t2_invr115542 * term.epsilon * ((-1.20000000000000000e+1 * sr12) + (1.20000000000000000e+1 * sr6)));
            double cse_p11_t1_g75012 = (rr * rr);
            double cse_p565_t3_invr115543 = (1.0 / (cse_p11_t1_g75012));
            double d2e_dr2 = (cse_p565_t3_invr115543 * term.epsilon * ((-4.20000000000000040e+1 * sr6) + (1.56000000000000000e+2 * sr12)));
            double energy = fa_atr;
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p67_t1_g79888 = (cse_p1_t2_invsqrt74902 * cse_p1_t4_invr274904 * term.epsilon * term.sigma * sr2 * sr2 * sr * (-2.00000000000000000e+0 + sr6 + sr6));
            double g_x1 = (-6.00000000000000000e+0 * cse_p67_t1_g79888 * dx);
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double g_y1 = (-6.00000000000000000e+0 * cse_p67_t1_g79888 * dy);
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double g_z1 = (-6.00000000000000000e+0 * cse_p67_t1_g79888 * dz);
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (6.00000000000000000e+0 * cse_p67_t1_g79888 * dx);
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (6.00000000000000000e+0 * cse_p67_t1_g79888 * dy);
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (6.00000000000000000e+0 * cse_p67_t1_g79888 * dz);
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double cse_p11_t1_g75010 = (cse_p565_t2_invr115542 * cse_p565_t2_invr115542 * d2e_dr2);
            double cse_p564_t6_g115468 = (cse_p11_t1_g75010 * (-(dx)));
            double cse_p564_t7_g115469 = (cse_p11_t1_g75010 * (-(dy)));
            double cse_p564_t8_g115470 = (cse_p11_t1_g75010 * (-(dz)));
            double cse_p564_t9_g115471 = (cse_p11_t1_g75010 * dx);
            double cse_p564_t10_g115472 = (cse_p11_t1_g75010 * dy);
            double cse_p564_t11_g115473 = (cse_p11_t1_g75010 * dz * dz);
            double cse_p70_t45_g80136 = (cse_p564_t9_g115471 * dx);
            double cse_p70_t46_g80137 = (cse_p564_t9_g115471 * dy);
            double cse_p565_t1_invr115541 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g75011 = (cse_p565_t1_invr115541 * de_dr);
            double cse_p564_t1_g115540 = (cse_p11_t1_g75011 * (cse_p11_t1_g75012 + (-(cse_p11_t49_g74953))));
            double cse_p564_t13_g115475 = (cse_p11_t1_g75011 * (cse_p11_t1_g75012 + (-(cse_p11_t50_g74954))));
            double cse_p564_t14_g115476 = (cse_p11_t1_g75011 * (cse_p11_t1_g75012 + (-(cse_p11_t51_g74955))));
            double cse_p564_t15_g115477 = (cse_p11_t1_g75011 * (cse_p11_t49_g74953 + (-(cse_p11_t1_g75012))));
            double cse_p564_t16_g115478 = (cse_p11_t1_g75011 * (cse_p11_t50_g74954 + (-(cse_p11_t1_g75012))));
            double cse_p564_t17_g115479 = (cse_p11_t1_g75011 * (cse_p11_t51_g74955 + (-(cse_p11_t1_g75012))));
            double cse_p564_t18_g115480 = (cse_p11_t1_g75011 * dy);
            double cse_p564_t19_g115481 = (cse_p11_t1_g75011 * dz);
            double cse_p69_t45_g80040 = (cse_p564_t9_g115471 * dz);
            double cse_p69_t74_g80069 = (cse_p564_t1_g115540 + cse_p70_t45_g80136);
            double cse_p69_t46_g80041 = (cse_p564_t18_g115480 * (-(dx)));
            double cse_p69_t47_g80042 = (cse_p564_t18_g115480 * dx);
            double cse_p69_t73_g80068 = (cse_p69_t46_g80041 + cse_p70_t46_g80137);
            double h_x1_x1 = cse_p69_t74_g80069;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = cse_p69_t73_g80068;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p69_t31_g80026 = (cse_p564_t19_g115481 * (-(dx)));
            double cse_p69_t32_g80027 = (cse_p564_t19_g115481 * (-(dy)));
            double cse_p69_t33_g80028 = (cse_p564_t19_g115481 * dx);
            double cse_p69_t34_g80029 = (cse_p564_t19_g115481 * dy);
            double cse_p60_t59_g79752 = (cse_p69_t31_g80026 + cse_p69_t45_g80040);
            double cse_p69_t48_g80043 = (cse_p564_t10_g115472 * dy);
            double cse_p69_t49_g80044 = (cse_p564_t10_g115472 * dz);
            double cse_p60_t60_g79753 = (cse_p69_t32_g80027 + cse_p69_t49_g80044);
            double h_x1_z1 = cse_p60_t59_g79752;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p69_t20_g80015 = (cse_p564_t6_g115468 * (-(dx)));
            double cse_p69_t21_g80016 = (cse_p564_t6_g115468 * (-(dy)));
            double cse_p69_t22_g80017 = (cse_p564_t6_g115468 * (-(dz)));
            double cse_p69_t23_g80018 = (cse_p564_t6_g115468 * dx);
            double cse_p69_t24_g80019 = (cse_p564_t6_g115468 * dy);
            double cse_p69_t25_g80020 = (cse_p564_t6_g115468 * dz);
            double cse_p60_t53_g79746 = (cse_p564_t15_g115477 + cse_p69_t23_g80018);
            double cse_p60_t61_g79754 = (cse_p69_t24_g80019 + cse_p69_t47_g80042);
            double cse_p60_t63_g79756 = (cse_p69_t21_g80016 + cse_p69_t46_g80041);
            double cse_p60_t64_g79757 = (cse_p69_t22_g80017 + cse_p69_t31_g80026);
            double cse_p60_t65_g79758 = (cse_p69_t25_g80020 + cse_p69_t33_g80028);
            double cse_p60_t82_g79775 = (cse_p564_t1_g115540 + cse_p69_t20_g80015);
            double h_x1_x2 = cse_p60_t53_g79746;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p69_t26_g80021 = (cse_p564_t7_g115469 * (-(dy)));
            double cse_p69_t27_g80022 = (cse_p564_t7_g115469 * (-(dz)));
            double cse_p69_t28_g80023 = (cse_p564_t7_g115469 * dx);
            double cse_p69_t29_g80024 = (cse_p564_t7_g115469 * dy);
            double cse_p69_t30_g80025 = (cse_p564_t7_g115469 * dz);
            double cse_p60_t54_g79747 = (cse_p564_t16_g115478 + cse_p69_t29_g80024);
            double cse_p60_t62_g79755 = (cse_p69_t28_g80023 + cse_p69_t47_g80042);
            double cse_p60_t66_g79759 = (cse_p69_t27_g80022 + cse_p69_t32_g80027);
            double cse_p60_t67_g79760 = (cse_p69_t30_g80025 + cse_p69_t34_g80029);
            double h_x1_y2 = cse_p60_t62_g79755;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p69_t35_g80030 = (cse_p564_t8_g115470 * (-(dz)));
            double cse_p69_t36_g80031 = (cse_p564_t8_g115470 * dx);
            double cse_p69_t37_g80032 = (cse_p564_t8_g115470 * dy);
            double cse_p69_t38_g80033 = (cse_p564_t8_g115470 * dz);
            double cse_p60_t55_g79748 = (cse_p564_t17_g115479 + cse_p69_t38_g80033);
            double cse_p60_t68_g79761 = (cse_p69_t33_g80028 + cse_p69_t36_g80031);
            double cse_p60_t69_g79762 = (cse_p69_t34_g80029 + cse_p69_t37_g80032);
            double h_x1_z2 = cse_p60_t68_g79761;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p67_t70_g79867 = (cse_p564_t13_g115475 + cse_p69_t48_g80043);
            double cse_p60_t78_g79771 = (cse_p564_t13_g115475 + cse_p69_t26_g80021);
            double h_y1_y1 = cse_p67_t70_g79867;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t60_g79753;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t61_g79754;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t54_g79747;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t69_g79762;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p67_t72_g79869 = (cse_p564_t11_g115473 + cse_p564_t14_g115476);
            double cse_p60_t80_g79773 = (cse_p564_t14_g115476 + cse_p69_t35_g80030);
            double h_z1_z1 = cse_p67_t72_g79869;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t65_g79758;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t67_g79760;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t55_g79748;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t82_g79775;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t63_g79756;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t64_g79757;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t78_g79771;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t66_g79759;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t80_g79773;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.r_cut)) {
              /* !COND-4 */
              double cse_p68_t2_g79899 = (-2.00000000000000000e+0 * term.bb * rr);
              double cse_p12_t7_g75026 = (term.bb * rr);
              double cse_p11_t1_g75012 = (rr * rr);
              double cse_p68_t4_g79901 = (-3.00000000000000000e+0 * term.aa * cse_p11_t1_g75012);
              double cse_p67_t56_g79853 = (cse_p68_t2_g79899 + cse_p68_t4_g79901);
              double cse_p12_t6_g75025 = (3.00000000000000000e+0 * term.aa * cse_p11_t1_g75012);
              double cse_p11_t88_g74992 = (term.cc + cse_p12_t6_g75025 + cse_p12_t7_g75026 + cse_p12_t7_g75026);
              double fa_atr = (term.dd + (term.cc * rr) + (cse_p11_t1_g75012 * (term.bb + (term.aa * rr))));
              double de_dr = (term.cc + cse_p12_t6_g75025 + (2.00000000000000000e+0 * term.bb * rr));
              double d2e_dr2 = ((2.00000000000000000e+0 * term.bb) + (6.00000000000000000e+0 * term.aa * rr));
              double energy = fa_atr;
              energy_added += energy;
              *energy_accumulate += energy;
              double g_x1 = (cse_p11_t88_g74992 * cse_p1_t2_invsqrt74902 * dx);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t88_g74992 * cse_p1_t2_invsqrt74902 * dy);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t88_g74992 * cse_p1_t2_invsqrt74902 * dz);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = ((cse_p1_t2_invsqrt74902 * cse_p67_t56_g79853 * dx) + (-((term.cc * cse_p1_t2_invsqrt74902 * dx))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = ((cse_p1_t2_invsqrt74902 * cse_p67_t56_g79853 * dy) + (-((term.cc * cse_p1_t2_invsqrt74902 * dy))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = ((cse_p1_t2_invsqrt74902 * cse_p67_t56_g79853 * dz) + (-((term.cc * cse_p1_t2_invsqrt74902 * dz))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double cse_p565_t2_invr115542 = (1.0 / (rr));
              double cse_p11_t1_g75010 = (cse_p565_t2_invr115542 * cse_p565_t2_invr115542 * d2e_dr2);
              double cse_p565_t1_invr115541 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g75011 = (cse_p565_t1_invr115541 * de_dr);
              double cse_p564_t1_g115540 = (cse_p11_t1_g75011 * (cse_p11_t1_g75012 + (-(cse_p11_t49_g74953))));
              double cse_p67_t1_g79889 = (cse_p11_t1_g75010 * dx);
              double cse_p564_t40_g115502 = (cse_p67_t1_g79889 * dz);
              double h_x1_x1 = (cse_p564_t1_g115540 + (cse_p67_t1_g79889 * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double cse_p67_t1_g79890 = (cse_p11_t1_g75011 * dy);
              double cse_p564_t41_g115503 = (cse_p67_t1_g79890 * (-(dx)));
              double cse_p564_t42_g115504 = (cse_p67_t1_g79890 * dx);
              double h_x1_y1 = (cse_p564_t41_g115503 + (cse_p67_t1_g79889 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g75015 = (cse_p11_t1_g75011 * dz);
              double cse_p564_t29_g115491 = (cse_p11_t1_g75015 * (-(dx)));
              double cse_p564_t30_g115492 = (cse_p11_t1_g75015 * (-(dy)));
              double cse_p564_t31_g115493 = (cse_p11_t1_g75015 * dx);
              double cse_p564_t32_g115494 = (cse_p11_t1_g75015 * dy);
              double cse_p67_t1_g79892 = (cse_p11_t1_g75010 * dy);
              double cse_p564_t43_g115505 = (cse_p67_t1_g79892 * dy);
              double cse_p564_t44_g115506 = (cse_p67_t1_g79892 * dz);
              double cse_p116_t57_g84244 = (cse_p564_t30_g115492 + cse_p564_t44_g115506);
              double h_x1_z1 = (cse_p564_t29_g115491 + cse_p564_t40_g115502);
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g75013 = (cse_p11_t1_g75010 * (-(dx)));
              double cse_p564_t20_g115482 = (cse_p11_t1_g75013 * (-(dx)));
              double cse_p564_t21_g115483 = (cse_p11_t1_g75013 * (-(dy)));
              double cse_p564_t22_g115484 = (cse_p11_t1_g75013 * (-(dz)));
              double cse_p564_t23_g115485 = (cse_p11_t1_g75013 * dy);
              double cse_p564_t24_g115486 = (cse_p11_t1_g75013 * dz);
              double cse_p116_t52_g84239 = (cse_p564_t22_g115484 + cse_p564_t29_g115491);
              double cse_p116_t53_g84240 = (cse_p564_t23_g115485 + cse_p564_t42_g115504);
              double cse_p116_t54_g84241 = (cse_p564_t24_g115486 + cse_p564_t31_g115493);
              double cse_p67_t60_g79857 = (cse_p564_t21_g115483 + cse_p564_t41_g115503);
              double h_x1_x2 = ((cse_p11_t1_g75011 * (cse_p11_t49_g74953 + (-(cse_p11_t1_g75012)))) + (cse_p11_t1_g75013 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g75014 = (cse_p11_t1_g75010 * (-(dy)));
              double cse_p564_t25_g115487 = (cse_p11_t1_g75014 * (-(dy)));
              double cse_p564_t26_g115488 = (cse_p11_t1_g75014 * (-(dz)));
              double cse_p564_t27_g115489 = (cse_p11_t1_g75014 * dy);
              double cse_p564_t28_g115490 = (cse_p11_t1_g75014 * dz);
              double cse_p116_t50_g84237 = (cse_p564_t27_g115489 + (cse_p11_t1_g75011 * (cse_p11_t50_g74954 + (-(cse_p11_t1_g75012)))));
              double cse_p116_t55_g84242 = (cse_p564_t26_g115488 + cse_p564_t30_g115492);
              double cse_p116_t56_g84243 = (cse_p564_t28_g115490 + cse_p564_t32_g115494);
              double h_x1_y2 = (cse_p564_t42_g115504 + (cse_p11_t1_g75014 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g75016 = (cse_p11_t1_g75010 * (-(dz)));
              double cse_p564_t33_g115495 = (cse_p11_t1_g75016 * (-(dz)));
              double cse_p564_t34_g115496 = (cse_p11_t1_g75016 * dy);
              double cse_p564_t35_g115497 = (cse_p11_t1_g75016 * dz);
              double cse_p116_t51_g84238 = (cse_p564_t35_g115497 + (cse_p11_t1_g75011 * (cse_p11_t51_g74955 + (-(cse_p11_t1_g75012)))));
              double cse_p116_t58_g84245 = (cse_p564_t32_g115494 + cse_p564_t34_g115496);
              double h_x1_z2 = (cse_p564_t31_g115493 + (cse_p11_t1_g75016 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g75017 = (cse_p11_t1_g75011 * (cse_p11_t1_g75012 + (-(cse_p11_t50_g74954))));
              double cse_p116_t63_g84250 = (cse_p11_t1_g75017 + cse_p564_t25_g115487);
              double h_y1_y1 = (cse_p11_t1_g75017 + cse_p564_t43_g115505);
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = cse_p116_t57_g84244;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = cse_p116_t53_g84240;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = cse_p116_t50_g84237;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = cse_p116_t58_g84245;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g75018 = (cse_p11_t1_g75011 * (cse_p11_t1_g75012 + (-(cse_p11_t51_g74955))));
              double cse_p116_t64_g84251 = (cse_p11_t1_g75018 + cse_p564_t33_g115495);
              double h_z1_z1 = (cse_p11_t1_g75018 + (cse_p11_t1_g75010 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = cse_p116_t54_g84241;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = cse_p116_t56_g84243;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = cse_p116_t51_g84238;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p564_t1_g115540 + cse_p564_t20_g115482);
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = cse_p67_t60_g79857;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = cse_p116_t52_g84239;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = cse_p116_t63_g84250;
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = cse_p116_t55_g84242;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = cse_p116_t64_g84251;
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
                double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t49_g74953))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t2_invr115542 = (1.0 / (rr));
                double cse_p565_t1_invr115541 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g75012 = (rr * rr);
                double cse_p11_t1_g75011 = (cse_p565_t1_invr115541 * de_dr);
                double cse_p11_t1_g75010 = (cse_p565_t2_invr115542 * cse_p565_t2_invr115542 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g75010 * dx * dy) + (cse_p11_t1_g75011 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double cse_p67_t1_g79892 = (cse_p11_t1_g75010 * dy);
                double cse_p67_t1_g79890 = (cse_p11_t1_g75011 * dy);
                double cse_p67_t1_g79889 = (cse_p11_t1_g75010 * dx);
                double h_x1_z1 = ((cse_p11_t1_g75011 * dz * (-(dx))) + (cse_p67_t1_g79889 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g75015 = (cse_p11_t1_g75011 * dz);
                double h_x1_x2 = ((cse_p11_t1_g75011 * (cse_p11_t49_g74953 + (-(cse_p11_t1_g75012)))) + (cse_p67_t1_g79889 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g75013 = (cse_p11_t1_g75010 * (-(dx)));
                double h_x1_y2 = ((cse_p67_t1_g79889 * (-(dy))) + (cse_p67_t1_g79890 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g75014 = (cse_p11_t1_g75010 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g75015 * dx) + (cse_p67_t1_g79889 * (-(dz))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g75016 = (cse_p11_t1_g75010 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g75011 * (cse_p11_t1_g75012 + (-(cse_p11_t50_g74954)))) + (cse_p67_t1_g79892 * dy));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g75017 = (cse_p11_t1_g75011 * (cse_p11_t1_g75012 + (-(cse_p11_t50_g74954))));
                double h_y1_z1 = ((cse_p11_t1_g75015 * (-(dy))) + (cse_p67_t1_g79892 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g75013 * dy) + (cse_p67_t1_g79890 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g75011 * (cse_p11_t50_g74954 + (-(cse_p11_t1_g75012)))) + (cse_p11_t1_g75014 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g75015 * dy) + (cse_p11_t1_g75016 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g75010 * dz * dz) + (cse_p11_t1_g75011 * (cse_p11_t1_g75012 + (-(cse_p11_t51_g74955)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g75018 = (cse_p11_t1_g75011 * (cse_p11_t1_g75012 + (-(cse_p11_t51_g74955))));
                double h_z1_x2 = ((cse_p11_t1_g75013 * dz) + (cse_p11_t1_g75015 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g75014 * dz) + (cse_p11_t1_g75015 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g75011 * (cse_p11_t51_g74955 + (-(cse_p11_t1_g75012)))) + (cse_p11_t1_g75016 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g75011 * (cse_p11_t1_g75012 + (-(cse_p11_t49_g74953)))) + (cse_p11_t1_g75013 * (-(dx))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p564_t1_g115540 = (cse_p11_t1_g75011 * (cse_p11_t1_g75012 + (-(cse_p11_t49_g74953))));
                double h_x2_y2 = ((cse_p11_t1_g75013 * (-(dy))) + (cse_p67_t1_g79890 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double cse_p67_t1_g79893 = (cse_p67_t1_g79890 * (-(dx)));
                double h_x2_z2 = ((cse_p11_t1_g75013 * (-(dz))) + (cse_p11_t1_g75015 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g75017 + (cse_p11_t1_g75014 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g75014 * (-(dz))) + (cse_p11_t1_g75015 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g75018 + (cse_p11_t1_g75016 * (-(dz))));
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
