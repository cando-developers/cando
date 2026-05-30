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
      double cse_p1_t1_sqrt440325 = sqrt(r2);
      double cse_p1_t2_invsqrt440326 = (1.0 / (cse_p1_t1_sqrt440325));
      double rr = cse_p1_t1_sqrt440325;
      if ((rr <= (6.00000023841857900e-1 * term.sigma))) {
        /* !COND-1 */
        double cse_p528_t7_g453866 = (cse_p1_t2_invsqrt440326 * dx);
        double cse_p528_t8_g453867 = (cse_p1_t2_invsqrt440326 * dy);
        double cse_p528_t9_g453868 = (cse_p1_t2_invsqrt440326 * dz);
        double cse_p528_t16_g453875 = (-(term.epsilon));
        double cse_p12_t1_g440393 = (term.mij * params.rep_weight);
        double cse_p11_t9_g440337 = (cse_p12_t1_g440393 * cse_p528_t7_g453866);
        double cse_p11_t10_g440338 = (cse_p12_t1_g440393 * cse_p528_t8_g453867);
        double cse_p11_t11_g440339 = (cse_p12_t1_g440393 * cse_p528_t9_g453868);
        double fa_rep = (term.bij + (term.mij * rr));
        double cse_p528_t11_g453870 = (fa_rep * params.rep_weight);
        double fa_atr = cse_p528_t16_g453875;
        double cse_p13_t20_g440419 = (cse_p528_t11_g453870 + fa_atr);
        double de_dr = cse_p12_t1_g440393;
        double energy = cse_p13_t20_g440419;
        energy_added += energy;
        *energy_accumulate += energy;
        double g_x1 = cse_p11_t9_g440337;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = cse_p11_t10_g440338;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = cse_p11_t11_g440339;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = (-(cse_p11_t9_g440337));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (-(cse_p11_t10_g440338));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (-(cse_p11_t11_g440339));
        KernelGradientAcc(term.i3x2, 2, g_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.sigma)) {
          /* !COND-2 */
          double cse_p529_t1_invr453887 = (1.0 / (rr));
          double cse_p479_t11_g452545 = (cse_p529_t1_invr453887 * cse_p529_t1_invr453887);
          double cse_p479_t12_g452546 = (cse_p529_t1_invr453887 * term.sigma);
          double sr = cse_p479_t12_g452546;
          double cse_p528_t1_g453886 = (cse_p1_t2_invsqrt440326 * term.sigma * sr);
          double cse_p528_t12_g453871 = (sr * sr);
          double cse_p11_t1_g440357 = (cse_p479_t11_g452545 * cse_p528_t1_g453886);
          double cse_p528_t4_g453863 = (cse_p11_t1_g440357 * dx);
          double cse_p528_t5_g453864 = (cse_p11_t1_g440357 * dy);
          double cse_p528_t6_g453865 = (cse_p11_t1_g440357 * dz);
          double cse_p12_t19_g440385 = (-(cse_p528_t4_g453863));
          double cse_p12_t20_g440386 = (-(cse_p528_t5_g453864));
          double cse_p12_t21_g440387 = (-(cse_p528_t6_g453865));
          double sr2 = cse_p528_t12_g453871;
          double cse_p528_t13_g453872 = (sr2 * sr2);
          double cse_p13_t4_g440403 = (6.00000000000000000e+0 * cse_p528_t13_g453872 * term.epsilon);
          double cse_p13_t8_g440407 = (cse_p528_t13_g453872 * sr2);
          double sr6 = cse_p13_t8_g440407;
          double cse_p528_t2_g453861 = (-2.00000000000000000e+0 * sr6);
          double cse_p528_t3_g453862 = (1.20000000000000000e+1 * sr6);
          double cse_p528_t14_g453873 = (sr6 * sr6);
          double cse_p528_t15_g453874 = (-2.00000000000000000e+0 + sr6 + sr6);
          double sr12 = cse_p528_t14_g453873;
          double cse_p528_t1_g453860 = (-1.20000000000000000e+1 * sr12);
          double cse_p14_t15_g440446 = (cse_p528_t2_g453861 + sr12);
          double cse_p13_t1_g440431 = (cse_p14_t15_g440446 * term.epsilon);
          double cse_p13_t19_g440418 = (cse_p528_t1_g453860 + cse_p528_t3_g453862);
          double fa_rep_epsilon = cse_p13_t1_g440431;
          double fa_rep = (term.epsilon + fa_rep_epsilon);
          double fa_atr = (-(term.epsilon));
          double de_dr = (cse_p13_t19_g440418 * cse_p529_t1_invr453887 * term.epsilon * params.rep_weight);
          double cse_p12_t1_g440397 = (fa_rep * params.rep_weight);
          double energy = (cse_p12_t1_g440397 + fa_atr);
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p11_t1_g440360 = (cse_p13_t4_g440403 * cse_p528_t15_g453874 * params.rep_weight);
          double g_x1 = (cse_p11_t1_g440360 * cse_p12_t19_g440385);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t1_g440360 * cse_p12_t20_g440386);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t1_g440360 * cse_p12_t21_g440387);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p11_t1_g440360 * cse_p528_t4_g453863);
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p11_t1_g440360 * cse_p528_t5_g453864);
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p11_t1_g440360 * cse_p528_t6_g453865);
          KernelGradientAcc(term.i3x2, 2, g_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rswitch)) {
            /* !COND-3 */
            double cse_p529_t1_invr453887 = (1.0 / (rr));
            double cse_p1_t4_invr2440328 = (cse_p529_t1_invr453887 * cse_p529_t1_invr453887);
            double sr = (cse_p529_t1_invr453887 * term.sigma);
            double sr2 = (sr * sr);
            double cse_p12_t1_g440394 = (sr2 * sr2);
            double sr6 = (cse_p12_t1_g440394 * sr2);
            double sr12 = (sr6 * sr6);
            double cse_p12_t1_g440392 = (-2.00000000000000000e+0 * sr6);
            double cse_p13_t1_g440431 = (term.epsilon * (cse_p12_t1_g440392 + sr12));
            double fa_atr = cse_p13_t1_g440431;
            double cse_p12_t1_g440395 = (-1.20000000000000000e+1 * sr12);
            double cse_p12_t1_g440396 = (1.20000000000000000e+1 * sr6);
            double de_dr = (cse_p529_t1_invr453887 * term.epsilon * (cse_p12_t1_g440395 + cse_p12_t1_g440396));
            double energy = fa_atr;
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p528_t1_g453886 = (cse_p1_t2_invsqrt440326 * term.sigma * sr);
            double cse_p11_t1_g440357 = (cse_p1_t4_invr2440328 * cse_p528_t1_g453886);
            double cse_p11_t1_g440358 = (6.00000000000000000e+0 * cse_p12_t1_g440394 * term.epsilon);
            double cse_p11_t1_g440359 = (cse_p11_t1_g440358 * (-2.00000000000000000e+0 + sr6 + sr6));
            double cse_p11_t1_g440363 = (cse_p11_t1_g440357 * dx);
            double g_x1 = (cse_p11_t1_g440359 * (-(cse_p11_t1_g440363)));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p11_t1_g440364 = (cse_p11_t1_g440357 * dy);
            double g_y1 = (cse_p11_t1_g440359 * (-(cse_p11_t1_g440364)));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p11_t1_g440365 = (cse_p11_t1_g440357 * dz);
            double g_z1 = (cse_p11_t1_g440359 * (-(cse_p11_t1_g440365)));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (cse_p11_t1_g440359 * cse_p11_t1_g440363);
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (cse_p11_t1_g440359 * cse_p11_t1_g440364);
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (cse_p11_t1_g440359 * cse_p11_t1_g440365);
            KernelGradientAcc(term.i3x2, 2, g_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p11_t1_g440362 = (term.aa * rr * rr);
              double cse_p12_t4_g440370 = (3.00000000000000000e+0 * cse_p11_t1_g440362);
              double cse_p11_t1_g440361 = (term.bb * rr);
              double cse_p11_t21_g440349 = (term.cc + cse_p11_t1_g440361 + cse_p11_t1_g440361 + cse_p12_t4_g440370);
              double fa_atr = (term.dd + (term.cc * rr) + (rr * rr * (term.bb + (term.aa * rr))));
              double de_dr = (term.cc + cse_p12_t4_g440370 + (2.00000000000000000e+0 * cse_p11_t1_g440361));
              double energy = fa_atr;
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p11_t1_g440366 = (cse_p1_t2_invsqrt440326 * dx);
              double g_x1 = (cse_p11_t1_g440366 * cse_p11_t21_g440349);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double cse_p12_t1_g440390 = (cse_p1_t2_invsqrt440326 * dy);
              double g_y1 = (cse_p11_t21_g440349 * cse_p12_t1_g440390);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double cse_p12_t1_g440391 = (cse_p1_t2_invsqrt440326 * dz);
              double g_z1 = (cse_p11_t21_g440349 * cse_p12_t1_g440391);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p63_t1_g441722 = (-2.00000000000000000e+0 * cse_p11_t1_g440361);
              double cse_p63_t1_g441723 = (-3.00000000000000000e+0 * cse_p11_t1_g440362);
              double g_x2 = ((cse_p11_t1_g440366 * cse_p63_t1_g441722) + (cse_p11_t1_g440366 * cse_p63_t1_g441723) + (-((term.cc * cse_p11_t1_g440366))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = ((cse_p12_t1_g440390 * cse_p63_t1_g441722) + (cse_p12_t1_g440390 * cse_p63_t1_g441723) + (-((term.cc * cse_p12_t1_g440390))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = ((cse_p12_t1_g440391 * cse_p63_t1_g441722) + (cse_p12_t1_g440391 * cse_p63_t1_g441723) + (-((term.cc * cse_p12_t1_g440391))));
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
      double cse_p11_t49_g463021 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t50_g463022 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t51_g463023 = (dz * dz);
      double r2 = (cse_p11_t49_g463021 + cse_p11_t50_g463022 + cse_p11_t51_g463023);
      double cse_p1_t1_sqrt462969 = sqrt(r2);
      double cse_p1_t2_invsqrt462970 = (1.0 / (cse_p1_t1_sqrt462969));
      double rr = cse_p1_t1_sqrt462969;
      if ((rr <= (6.00000023841857900e-1 * term.sigma))) {
        /* !COND-1 */
        double cse_p564_t46_g503576 = (rr * rr);
        double cse_p564_t58_g503588 = (-(cse_p11_t49_g463021));
        double cse_p564_t59_g503589 = (-(cse_p11_t50_g463022));
        double cse_p564_t60_g503590 = (-(cse_p11_t51_g463023));
        double cse_p564_t61_g503591 = (-(dx));
        double cse_p564_t62_g503592 = (-(dy));
        double cse_p564_t63_g503593 = (-(dz));
        double cse_p564_t64_g503594 = (-(term.epsilon));
        double cse_p565_t1_invr503609 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr503610 = (1.0 / (rr));
        double cse_p564_t68_g503598 = (rr * rr * rr);
        double cse_p15_t89_g463529 = (-(cse_p564_t46_g503576));
        double cse_p14_t90_g463411 = (cse_p11_t49_g463021 + cse_p15_t89_g463529);
        double cse_p14_t92_g463413 = (cse_p11_t50_g463022 + cse_p15_t89_g463529);
        double cse_p14_t94_g463415 = (cse_p11_t51_g463023 + cse_p15_t89_g463529);
        double cse_p13_t116_g463304 = (cse_p564_t46_g503576 + cse_p564_t58_g503588);
        double cse_p13_t117_g463305 = (cse_p564_t46_g503576 + cse_p564_t59_g503589);
        double cse_p13_t118_g463306 = (cse_p564_t46_g503576 + cse_p564_t60_g503590);
        double cse_p11_t11_g462983 = (cse_p1_t2_invsqrt462970 * dx * term.mij * params.rep_weight);
        double cse_p11_t12_g462984 = (cse_p1_t2_invsqrt462970 * dy * term.mij * params.rep_weight);
        double cse_p11_t13_g462985 = (cse_p1_t2_invsqrt462970 * dz * term.mij * params.rep_weight);
        double fa_rep = (term.bij + (term.mij * rr));
        double cse_p564_t45_g503575 = (fa_rep * params.rep_weight);
        double fa_atr = cse_p564_t64_g503594;
        double cse_p60_t83_g467844 = (cse_p564_t45_g503575 + fa_atr);
        double de_dr = (term.mij * params.rep_weight);
        double cse_p515_t36_g499744 = (cse_p565_t1_invr503609 * de_dr);
        double cse_p13_t68_g463256 = (cse_p14_t90_g463411 * cse_p515_t36_g499744);
        double cse_p13_t69_g463257 = (cse_p14_t92_g463413 * cse_p515_t36_g499744);
        double cse_p13_t70_g463258 = (cse_p14_t94_g463415 * cse_p515_t36_g499744);
        double cse_p13_t71_g463259 = (cse_p515_t36_g499744 * dx * dy);
        double cse_p13_t72_g463260 = (cse_p515_t36_g499744 * cse_p564_t61_g503591 * dy);
        double cse_p12_t23_g463110 = (cse_p13_t116_g463304 * cse_p515_t36_g499744);
        double cse_p12_t24_g463111 = (cse_p13_t117_g463305 * cse_p515_t36_g499744);
        double cse_p12_t25_g463112 = (cse_p13_t118_g463306 * cse_p515_t36_g499744);
        double cse_p12_t28_g463115 = (cse_p515_t36_g499744 * dz);
        double cse_p13_t58_g463246 = (cse_p12_t28_g463115 * cse_p564_t61_g503591);
        double cse_p13_t59_g463247 = (cse_p12_t28_g463115 * cse_p564_t62_g503592);
        double cse_p13_t60_g463248 = (cse_p12_t28_g463115 * dx);
        double cse_p13_t61_g463249 = (cse_p12_t28_g463115 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p515_t37_g499745 = (cse_p565_t2_invr503610 * cse_p565_t2_invr503610 * d2e_dr2);
        double cse_p13_t62_g463250 = (cse_p515_t37_g499745 * dx * dx);
        double cse_p13_t63_g463251 = (cse_p515_t37_g499745 * dx * dy);
        double cse_p13_t64_g463252 = (cse_p515_t37_g499745 * dx * dz);
        double cse_p13_t65_g463253 = (cse_p515_t37_g499745 * dy * dy);
        double cse_p13_t66_g463254 = (cse_p515_t37_g499745 * dy * dz);
        double cse_p13_t67_g463255 = (cse_p515_t37_g499745 * dz * dz);
        double cse_p12_t11_g463098 = (cse_p515_t37_g499745 * cse_p564_t61_g503591);
        double cse_p13_t43_g463231 = (cse_p12_t11_g463098 * cse_p564_t61_g503591);
        double cse_p13_t44_g463232 = (cse_p12_t11_g463098 * cse_p564_t62_g503592);
        double cse_p13_t45_g463233 = (cse_p12_t11_g463098 * cse_p564_t63_g503593);
        double cse_p13_t46_g463234 = (cse_p12_t11_g463098 * dx);
        double cse_p13_t47_g463235 = (cse_p12_t11_g463098 * dy);
        double cse_p13_t48_g463236 = (cse_p12_t11_g463098 * dz);
        double cse_p12_t12_g463099 = (cse_p515_t37_g499745 * cse_p564_t62_g503592);
        double cse_p13_t49_g463237 = (cse_p12_t12_g463099 * cse_p564_t62_g503592);
        double cse_p13_t50_g463238 = (cse_p12_t12_g463099 * cse_p564_t63_g503593);
        double cse_p13_t51_g463239 = (cse_p12_t12_g463099 * dx);
        double cse_p13_t52_g463240 = (cse_p12_t12_g463099 * dy);
        double cse_p13_t53_g463241 = (cse_p12_t12_g463099 * dz);
        double cse_p12_t13_g463100 = (cse_p515_t37_g499745 * cse_p564_t63_g503593);
        double cse_p13_t54_g463242 = (cse_p12_t13_g463100 * cse_p564_t63_g503593);
        double cse_p13_t55_g463243 = (cse_p12_t13_g463100 * dx);
        double cse_p13_t56_g463244 = (cse_p12_t13_g463100 * dy);
        double cse_p13_t57_g463245 = (cse_p12_t13_g463100 * dz);
        double cse_p12_t63_g463150 = (cse_p13_t63_g463251 + cse_p13_t72_g463260);
        double cse_p12_t64_g463151 = (cse_p13_t58_g463246 + cse_p13_t64_g463252);
        double cse_p12_t65_g463152 = (cse_p13_t59_g463247 + cse_p13_t66_g463254);
        double cse_p12_t82_g463169 = (cse_p12_t24_g463111 + cse_p13_t65_g463253);
        double cse_p12_t84_g463171 = (cse_p12_t25_g463112 + cse_p13_t67_g463255);
        double cse_p12_t86_g463173 = (cse_p12_t23_g463110 + cse_p13_t62_g463250);
        double cse_p12_t58_g463145 = (cse_p13_t46_g463234 + cse_p13_t68_g463256);
        double cse_p12_t66_g463153 = (cse_p13_t47_g463235 + cse_p13_t71_g463259);
        double cse_p12_t68_g463155 = (cse_p13_t44_g463232 + cse_p13_t72_g463260);
        double cse_p12_t69_g463156 = (cse_p13_t45_g463233 + cse_p13_t58_g463246);
        double cse_p12_t70_g463157 = (cse_p13_t48_g463236 + cse_p13_t60_g463248);
        double cse_p12_t87_g463174 = (cse_p12_t23_g463110 + cse_p13_t43_g463231);
        double cse_p12_t59_g463146 = (cse_p13_t52_g463240 + cse_p13_t69_g463257);
        double cse_p12_t67_g463154 = (cse_p13_t51_g463239 + cse_p13_t71_g463259);
        double cse_p12_t71_g463158 = (cse_p13_t50_g463238 + cse_p13_t59_g463247);
        double cse_p12_t72_g463159 = (cse_p13_t53_g463241 + cse_p13_t61_g463249);
        double cse_p12_t83_g463170 = (cse_p12_t24_g463111 + cse_p13_t49_g463237);
        double cse_p12_t60_g463147 = (cse_p13_t57_g463245 + cse_p13_t70_g463258);
        double cse_p12_t73_g463160 = (cse_p13_t55_g463243 + cse_p13_t60_g463248);
        double cse_p12_t74_g463161 = (cse_p13_t56_g463244 + cse_p13_t61_g463249);
        double cse_p12_t85_g463172 = (cse_p12_t25_g463112 + cse_p13_t54_g463242);
        double energy = cse_p60_t83_g467844;
        energy_added += energy;
        *energy_accumulate += energy;
        double g_x1 = cse_p11_t11_g462983;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = cse_p11_t12_g462984;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = cse_p11_t13_g462985;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = (-(cse_p11_t11_g462983));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (-(cse_p11_t12_g462984));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (-(cse_p11_t13_g462985));
        KernelGradientAcc(term.i3x2, 2, g_z2);
        double h_x1_x1 = cse_p12_t86_g463173;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t63_g463150;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t64_g463151;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t58_g463145;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t67_g463154;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t73_g463160;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t82_g463169;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t65_g463152;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t66_g463153;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t59_g463146;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t74_g463161;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t84_g463171;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t70_g463157;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t72_g463159;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t60_g463147;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t87_g463174;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t68_g463155;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t69_g463156;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t83_g463170;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t71_g463158;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t85_g463172;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.sigma)) {
          /* !COND-2 */
          double cse_p11_t1_g463080 = (rr * rr);
          double cse_p564_t51_g503581 = (cse_p11_t1_g463080 + (-(cse_p11_t49_g463021)));
          double cse_p564_t52_g503582 = (cse_p11_t1_g463080 + (-(cse_p11_t50_g463022)));
          double cse_p564_t53_g503583 = (cse_p11_t1_g463080 + (-(cse_p11_t51_g463023)));
          double cse_p564_t54_g503584 = (cse_p11_t49_g463021 + (-(cse_p11_t1_g463080)));
          double cse_p564_t55_g503585 = (cse_p11_t50_g463022 + (-(cse_p11_t1_g463080)));
          double cse_p564_t56_g503586 = (cse_p11_t51_g463023 + (-(cse_p11_t1_g463080)));
          double cse_p564_t57_g503587 = (-(cse_p11_t1_g463080));
          double cse_p565_t3_invr503611 = (1.0 / (cse_p11_t1_g463080));
          double cse_p565_t2_invr503610 = (1.0 / (rr));
          double cse_p564_t38_g503568 = (cse_p565_t2_invr503610 * cse_p565_t2_invr503610);
          double cse_p564_t39_g503569 = (cse_p565_t2_invr503610 * term.sigma);
          double sr = cse_p564_t39_g503569;
          double cse_p564_t47_g503577 = (sr * sr);
          double cse_p67_t1_g467962 = (cse_p1_t2_invsqrt462970 * cse_p564_t38_g503568 * term.sigma * sr);
          double cse_p12_t8_g463095 = (cse_p67_t1_g467962 * dx);
          double cse_p12_t9_g463096 = (cse_p67_t1_g467962 * dy);
          double cse_p12_t10_g463097 = (cse_p67_t1_g467962 * dz);
          double cse_p11_t91_g463063 = (-(cse_p12_t8_g463095));
          double cse_p11_t92_g463064 = (-(cse_p12_t9_g463096));
          double cse_p11_t93_g463065 = (-(cse_p12_t10_g463097));
          double sr2 = cse_p564_t47_g503577;
          double cse_p564_t48_g503578 = (sr2 * sr2 * sr2);
          double sr6 = cse_p564_t48_g503578;
          double cse_p564_t2_g503532 = (-2.00000000000000000e+0 * sr6);
          double cse_p564_t3_g503533 = (-4.20000000000000040e+1 * sr6);
          double cse_p564_t4_g503534 = (1.20000000000000000e+1 * sr6);
          double cse_p564_t49_g503579 = (sr6 * sr6);
          double cse_p564_t50_g503580 = (-2.00000000000000000e+0 + sr6 + sr6);
          double sr12 = cse_p564_t49_g503579;
          double cse_p564_t1_g503531 = (-1.20000000000000000e+1 * sr12);
          double cse_p564_t5_g503535 = (1.56000000000000000e+2 * sr12);
          double cse_p67_t78_g467943 = (cse_p564_t2_g503532 + sr12);
          double cse_p60_t47_g467808 = (cse_p67_t78_g467943 * term.epsilon);
          double cse_p60_t56_g467817 = (cse_p564_t1_g503531 + cse_p564_t4_g503534);
          double cse_p60_t57_g467818 = (cse_p564_t3_g503533 + cse_p564_t5_g503535);
          double fa_rep_epsilon = cse_p60_t47_g467808;
          double fa_rep = (term.epsilon + fa_rep_epsilon);
          double fa_atr = (-(term.epsilon));
          double de_dr = (cse_p565_t2_invr503610 * cse_p60_t56_g467817 * term.epsilon * params.rep_weight);
          double d2e_dr2 = (cse_p565_t3_invr503611 * cse_p60_t57_g467818 * term.epsilon * params.rep_weight);
          double cse_p564_t37_g503567 = (cse_p565_t2_invr503610 * cse_p565_t2_invr503610 * d2e_dr2);
          double energy = (fa_atr + (fa_rep * params.rep_weight));
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p67_t1_g467959 = (6.00000000000000000e+0 * cse_p564_t50_g503580 * term.epsilon * params.rep_weight * sr2 * sr2);
          double g_x1 = (cse_p11_t91_g463063 * cse_p67_t1_g467959);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t92_g463064 * cse_p67_t1_g467959);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t93_g463065 * cse_p67_t1_g467959);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p67_t1_g467959 * cse_p67_t1_g467962 * dx);
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p67_t1_g467959 * cse_p67_t1_g467962 * dy);
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p67_t1_g467959 * cse_p67_t1_g467962 * dz);
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double cse_p68_t8_g467973 = (cse_p564_t37_g503567 * (-(dx)));
          double cse_p68_t9_g467974 = (cse_p564_t37_g503567 * (-(dy)));
          double cse_p68_t10_g467975 = (cse_p564_t37_g503567 * (-(dz)));
          double cse_p68_t11_g467976 = (cse_p564_t37_g503567 * dx);
          double cse_p68_t12_g467977 = (cse_p564_t37_g503567 * dy);
          double cse_p68_t13_g467978 = (cse_p564_t37_g503567 * dz * dz);
          double cse_p68_t45_g468010 = (cse_p68_t11_g467976 * dx);
          double cse_p68_t46_g468011 = (cse_p68_t11_g467976 * dy);
          double cse_p68_t47_g468012 = (cse_p68_t11_g467976 * dz);
          double cse_p68_t50_g468015 = (cse_p68_t12_g467977 * dy);
          double cse_p68_t51_g468016 = (cse_p68_t12_g467977 * dz);
          double cse_p14_t47_g463368 = (cse_p68_t8_g467973 * (-(dx)));
          double cse_p14_t48_g463369 = (cse_p68_t8_g467973 * (-(dy)));
          double cse_p14_t49_g463370 = (cse_p68_t8_g467973 * (-(dz)));
          double cse_p14_t50_g463371 = (cse_p68_t8_g467973 * dx);
          double cse_p14_t51_g463372 = (cse_p68_t8_g467973 * dy);
          double cse_p14_t52_g463373 = (cse_p68_t8_g467973 * dz);
          double cse_p14_t53_g463374 = (cse_p68_t9_g467974 * (-(dy)));
          double cse_p14_t54_g463375 = (cse_p68_t9_g467974 * (-(dz)));
          double cse_p14_t55_g463376 = (cse_p68_t9_g467974 * dx);
          double cse_p14_t56_g463377 = (cse_p68_t9_g467974 * dy);
          double cse_p14_t57_g463378 = (cse_p68_t9_g467974 * dz);
          double cse_p14_t58_g463379 = (cse_p68_t10_g467975 * (-(dz)));
          double cse_p14_t59_g463380 = (cse_p68_t10_g467975 * dx);
          double cse_p14_t60_g463381 = (cse_p68_t10_g467975 * dy);
          double cse_p14_t61_g463382 = (cse_p68_t10_g467975 * dz);
          double cse_p565_t1_invr503609 = (1.0 / ((rr * rr * rr)));
          double cse_p564_t36_g503566 = (cse_p565_t1_invr503609 * de_dr);
          double cse_p68_t20_g467985 = (cse_p564_t36_g503566 * dy);
          double cse_p68_t21_g467986 = (cse_p564_t36_g503566 * dz);
          double cse_p68_t48_g468013 = (cse_p68_t20_g467985 * (-(dx)));
          double cse_p68_t49_g468014 = (cse_p68_t20_g467985 * dx);
          double cse_p60_t15_g467776 = (cse_p564_t36_g503566 * (cse_p11_t49_g463021 + cse_p564_t57_g503587));
          double cse_p60_t16_g467777 = (cse_p564_t36_g503566 * (cse_p11_t50_g463022 + cse_p564_t57_g503587));
          double cse_p60_t17_g467778 = (cse_p564_t36_g503566 * (cse_p11_t51_g463023 + cse_p564_t57_g503587));
          double cse_p60_t18_g467779 = (cse_p564_t36_g503566 * cse_p564_t51_g503581);
          double cse_p60_t19_g467780 = (cse_p564_t36_g503566 * cse_p564_t52_g503582);
          double cse_p60_t20_g467781 = (cse_p564_t36_g503566 * cse_p564_t53_g503583);
          double cse_p60_t58_g467819 = (cse_p68_t46_g468011 + cse_p68_t48_g468013);
          double cse_p14_t43_g463364 = (cse_p68_t21_g467986 * (-(dx)));
          double cse_p14_t44_g463365 = (cse_p68_t21_g467986 * (-(dy)));
          double cse_p14_t45_g463366 = (cse_p68_t21_g467986 * dx);
          double cse_p14_t46_g463367 = (cse_p68_t21_g467986 * dy);
          double cse_p13_t114_g463302 = (cse_p60_t18_g467779 + cse_p68_t45_g468010);
          double h_x1_x1 = cse_p13_t114_g463302;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t58_g467819;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t89_g463277 = (cse_p14_t43_g463364 + cse_p68_t47_g468012);
          double cse_p13_t90_g463278 = (cse_p14_t44_g463365 + cse_p68_t51_g468016);
          double h_x1_z1 = cse_p13_t89_g463277;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t83_g463271 = (cse_p14_t50_g463371 + cse_p60_t15_g467776);
          double cse_p13_t91_g463279 = (cse_p14_t51_g463372 + cse_p68_t49_g468014);
          double cse_p13_t93_g463281 = (cse_p14_t48_g463369 + cse_p68_t48_g468013);
          double cse_p13_t94_g463282 = (cse_p14_t43_g463364 + cse_p14_t49_g463370);
          double cse_p13_t95_g463283 = (cse_p14_t45_g463366 + cse_p14_t52_g463373);
          double cse_p13_t115_g463303 = (cse_p14_t47_g463368 + cse_p60_t18_g467779);
          double h_x1_x2 = cse_p13_t83_g463271;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t84_g463272 = (cse_p14_t56_g463377 + cse_p60_t16_g467777);
          double cse_p13_t92_g463280 = (cse_p14_t55_g463376 + cse_p68_t49_g468014);
          double cse_p13_t96_g463284 = (cse_p14_t44_g463365 + cse_p14_t54_g463375);
          double cse_p13_t97_g463285 = (cse_p14_t46_g463367 + cse_p14_t57_g463378);
          double h_x1_y2 = cse_p13_t92_g463280;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t85_g463273 = (cse_p14_t61_g463382 + cse_p60_t17_g467778);
          double cse_p13_t98_g463286 = (cse_p14_t45_g463366 + cse_p14_t59_g463380);
          double cse_p13_t99_g463287 = (cse_p14_t46_g463367 + cse_p14_t60_g463381);
          double h_x1_z2 = cse_p13_t98_g463286;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t110_g463298 = (cse_p60_t19_g467780 + cse_p68_t50_g468015);
          double cse_p13_t111_g463299 = (cse_p14_t53_g463374 + cse_p60_t19_g467780);
          double h_y1_y1 = cse_p13_t110_g463298;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t90_g463278;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t91_g463279;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t84_g463272;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t99_g463287;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t112_g463300 = (cse_p60_t20_g467781 + cse_p68_t13_g467978);
          double cse_p13_t113_g463301 = (cse_p14_t58_g463379 + cse_p60_t20_g467781);
          double h_z1_z1 = cse_p13_t112_g463300;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t95_g463283;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t97_g463285;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t85_g463273;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t115_g463303;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t93_g463281;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t94_g463282;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t111_g463299;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t96_g463284;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t113_g463301;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rswitch)) {
            /* !COND-3 */
            double cse_p565_t2_invr503610 = (1.0 / (rr));
            double cse_p1_t4_invr2462972 = (cse_p565_t2_invr503610 * cse_p565_t2_invr503610);
            double sr = (cse_p565_t2_invr503610 * term.sigma);
            double sr2 = (sr * sr);
            double sr6 = (sr2 * sr2 * sr2);
            double sr12 = (sr6 * sr6);
            double fa_atr = (term.epsilon * (sr12 + (-2.00000000000000000e+0 * sr6)));
            double de_dr = (cse_p565_t2_invr503610 * term.epsilon * ((-1.20000000000000000e+1 * sr12) + (1.20000000000000000e+1 * sr6)));
            double cse_p11_t1_g463080 = (rr * rr);
            double cse_p565_t3_invr503611 = (1.0 / (cse_p11_t1_g463080));
            double d2e_dr2 = (cse_p565_t3_invr503611 * term.epsilon * ((-4.20000000000000040e+1 * sr6) + (1.56000000000000000e+2 * sr12)));
            double energy = fa_atr;
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p67_t1_g467956 = (cse_p1_t2_invsqrt462970 * cse_p1_t4_invr2462972 * term.epsilon * term.sigma * sr2 * sr2 * sr * (-2.00000000000000000e+0 + sr6 + sr6));
            double g_x1 = (-6.00000000000000000e+0 * cse_p67_t1_g467956 * dx);
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double g_y1 = (-6.00000000000000000e+0 * cse_p67_t1_g467956 * dy);
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double g_z1 = (-6.00000000000000000e+0 * cse_p67_t1_g467956 * dz);
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (6.00000000000000000e+0 * cse_p67_t1_g467956 * dx);
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (6.00000000000000000e+0 * cse_p67_t1_g467956 * dy);
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (6.00000000000000000e+0 * cse_p67_t1_g467956 * dz);
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double cse_p11_t1_g463078 = (cse_p565_t2_invr503610 * cse_p565_t2_invr503610 * d2e_dr2);
            double cse_p564_t6_g503536 = (cse_p11_t1_g463078 * (-(dx)));
            double cse_p564_t7_g503537 = (cse_p11_t1_g463078 * (-(dy)));
            double cse_p564_t8_g503538 = (cse_p11_t1_g463078 * (-(dz)));
            double cse_p564_t9_g503539 = (cse_p11_t1_g463078 * dx);
            double cse_p564_t10_g503540 = (cse_p11_t1_g463078 * dy);
            double cse_p564_t11_g503541 = (cse_p11_t1_g463078 * dz * dz);
            double cse_p70_t45_g468204 = (cse_p564_t9_g503539 * dx);
            double cse_p70_t46_g468205 = (cse_p564_t9_g503539 * dy);
            double cse_p565_t1_invr503609 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g463079 = (cse_p565_t1_invr503609 * de_dr);
            double cse_p564_t1_g503608 = (cse_p11_t1_g463079 * (cse_p11_t1_g463080 + (-(cse_p11_t49_g463021))));
            double cse_p564_t13_g503543 = (cse_p11_t1_g463079 * (cse_p11_t1_g463080 + (-(cse_p11_t50_g463022))));
            double cse_p564_t14_g503544 = (cse_p11_t1_g463079 * (cse_p11_t1_g463080 + (-(cse_p11_t51_g463023))));
            double cse_p564_t15_g503545 = (cse_p11_t1_g463079 * (cse_p11_t49_g463021 + (-(cse_p11_t1_g463080))));
            double cse_p564_t16_g503546 = (cse_p11_t1_g463079 * (cse_p11_t50_g463022 + (-(cse_p11_t1_g463080))));
            double cse_p564_t17_g503547 = (cse_p11_t1_g463079 * (cse_p11_t51_g463023 + (-(cse_p11_t1_g463080))));
            double cse_p564_t18_g503548 = (cse_p11_t1_g463079 * dy);
            double cse_p564_t19_g503549 = (cse_p11_t1_g463079 * dz);
            double cse_p69_t45_g468108 = (cse_p564_t9_g503539 * dz);
            double cse_p69_t74_g468137 = (cse_p564_t1_g503608 + cse_p70_t45_g468204);
            double cse_p69_t46_g468109 = (cse_p564_t18_g503548 * (-(dx)));
            double cse_p69_t47_g468110 = (cse_p564_t18_g503548 * dx);
            double cse_p69_t73_g468136 = (cse_p69_t46_g468109 + cse_p70_t46_g468205);
            double h_x1_x1 = cse_p69_t74_g468137;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = cse_p69_t73_g468136;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p69_t31_g468094 = (cse_p564_t19_g503549 * (-(dx)));
            double cse_p69_t32_g468095 = (cse_p564_t19_g503549 * (-(dy)));
            double cse_p69_t33_g468096 = (cse_p564_t19_g503549 * dx);
            double cse_p69_t34_g468097 = (cse_p564_t19_g503549 * dy);
            double cse_p60_t59_g467820 = (cse_p69_t31_g468094 + cse_p69_t45_g468108);
            double cse_p69_t48_g468111 = (cse_p564_t10_g503540 * dy);
            double cse_p69_t49_g468112 = (cse_p564_t10_g503540 * dz);
            double cse_p60_t60_g467821 = (cse_p69_t32_g468095 + cse_p69_t49_g468112);
            double h_x1_z1 = cse_p60_t59_g467820;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p69_t20_g468083 = (cse_p564_t6_g503536 * (-(dx)));
            double cse_p69_t21_g468084 = (cse_p564_t6_g503536 * (-(dy)));
            double cse_p69_t22_g468085 = (cse_p564_t6_g503536 * (-(dz)));
            double cse_p69_t23_g468086 = (cse_p564_t6_g503536 * dx);
            double cse_p69_t24_g468087 = (cse_p564_t6_g503536 * dy);
            double cse_p69_t25_g468088 = (cse_p564_t6_g503536 * dz);
            double cse_p60_t53_g467814 = (cse_p564_t15_g503545 + cse_p69_t23_g468086);
            double cse_p60_t61_g467822 = (cse_p69_t24_g468087 + cse_p69_t47_g468110);
            double cse_p60_t63_g467824 = (cse_p69_t21_g468084 + cse_p69_t46_g468109);
            double cse_p60_t64_g467825 = (cse_p69_t22_g468085 + cse_p69_t31_g468094);
            double cse_p60_t65_g467826 = (cse_p69_t25_g468088 + cse_p69_t33_g468096);
            double cse_p60_t82_g467843 = (cse_p564_t1_g503608 + cse_p69_t20_g468083);
            double h_x1_x2 = cse_p60_t53_g467814;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p69_t26_g468089 = (cse_p564_t7_g503537 * (-(dy)));
            double cse_p69_t27_g468090 = (cse_p564_t7_g503537 * (-(dz)));
            double cse_p69_t28_g468091 = (cse_p564_t7_g503537 * dx);
            double cse_p69_t29_g468092 = (cse_p564_t7_g503537 * dy);
            double cse_p69_t30_g468093 = (cse_p564_t7_g503537 * dz);
            double cse_p60_t54_g467815 = (cse_p564_t16_g503546 + cse_p69_t29_g468092);
            double cse_p60_t62_g467823 = (cse_p69_t28_g468091 + cse_p69_t47_g468110);
            double cse_p60_t66_g467827 = (cse_p69_t27_g468090 + cse_p69_t32_g468095);
            double cse_p60_t67_g467828 = (cse_p69_t30_g468093 + cse_p69_t34_g468097);
            double h_x1_y2 = cse_p60_t62_g467823;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p69_t35_g468098 = (cse_p564_t8_g503538 * (-(dz)));
            double cse_p69_t36_g468099 = (cse_p564_t8_g503538 * dx);
            double cse_p69_t37_g468100 = (cse_p564_t8_g503538 * dy);
            double cse_p69_t38_g468101 = (cse_p564_t8_g503538 * dz);
            double cse_p60_t55_g467816 = (cse_p564_t17_g503547 + cse_p69_t38_g468101);
            double cse_p60_t68_g467829 = (cse_p69_t33_g468096 + cse_p69_t36_g468099);
            double cse_p60_t69_g467830 = (cse_p69_t34_g468097 + cse_p69_t37_g468100);
            double h_x1_z2 = cse_p60_t68_g467829;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p67_t70_g467935 = (cse_p564_t13_g503543 + cse_p69_t48_g468111);
            double cse_p60_t78_g467839 = (cse_p564_t13_g503543 + cse_p69_t26_g468089);
            double h_y1_y1 = cse_p67_t70_g467935;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t60_g467821;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t61_g467822;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t54_g467815;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t69_g467830;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p67_t72_g467937 = (cse_p564_t11_g503541 + cse_p564_t14_g503544);
            double cse_p60_t80_g467841 = (cse_p564_t14_g503544 + cse_p69_t35_g468098);
            double h_z1_z1 = cse_p67_t72_g467937;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t65_g467826;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t67_g467828;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t55_g467816;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t82_g467843;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t63_g467824;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t64_g467825;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t78_g467839;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t66_g467827;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t80_g467841;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p68_t2_g467967 = (-2.00000000000000000e+0 * term.bb * rr);
              double cse_p12_t7_g463094 = (term.bb * rr);
              double cse_p11_t1_g463080 = (rr * rr);
              double cse_p68_t4_g467969 = (-3.00000000000000000e+0 * term.aa * cse_p11_t1_g463080);
              double cse_p67_t56_g467921 = (cse_p68_t2_g467967 + cse_p68_t4_g467969);
              double cse_p12_t6_g463093 = (3.00000000000000000e+0 * term.aa * cse_p11_t1_g463080);
              double cse_p11_t88_g463060 = (term.cc + cse_p12_t6_g463093 + cse_p12_t7_g463094 + cse_p12_t7_g463094);
              double fa_atr = (term.dd + (term.cc * rr) + (cse_p11_t1_g463080 * (term.bb + (term.aa * rr))));
              double de_dr = (term.cc + cse_p12_t6_g463093 + (2.00000000000000000e+0 * term.bb * rr));
              double d2e_dr2 = ((2.00000000000000000e+0 * term.bb) + (6.00000000000000000e+0 * term.aa * rr));
              double energy = fa_atr;
              energy_added += energy;
              *energy_accumulate += energy;
              double g_x1 = (cse_p11_t88_g463060 * cse_p1_t2_invsqrt462970 * dx);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t88_g463060 * cse_p1_t2_invsqrt462970 * dy);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t88_g463060 * cse_p1_t2_invsqrt462970 * dz);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = ((cse_p1_t2_invsqrt462970 * cse_p67_t56_g467921 * dx) + (-((term.cc * cse_p1_t2_invsqrt462970 * dx))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = ((cse_p1_t2_invsqrt462970 * cse_p67_t56_g467921 * dy) + (-((term.cc * cse_p1_t2_invsqrt462970 * dy))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = ((cse_p1_t2_invsqrt462970 * cse_p67_t56_g467921 * dz) + (-((term.cc * cse_p1_t2_invsqrt462970 * dz))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double cse_p565_t2_invr503610 = (1.0 / (rr));
              double cse_p11_t1_g463078 = (cse_p565_t2_invr503610 * cse_p565_t2_invr503610 * d2e_dr2);
              double cse_p565_t1_invr503609 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g463079 = (cse_p565_t1_invr503609 * de_dr);
              double cse_p564_t1_g503608 = (cse_p11_t1_g463079 * (cse_p11_t1_g463080 + (-(cse_p11_t49_g463021))));
              double cse_p67_t1_g467957 = (cse_p11_t1_g463078 * dx);
              double cse_p564_t40_g503570 = (cse_p67_t1_g467957 * dz);
              double h_x1_x1 = (cse_p564_t1_g503608 + (cse_p67_t1_g467957 * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double cse_p67_t1_g467958 = (cse_p11_t1_g463079 * dy);
              double cse_p564_t41_g503571 = (cse_p67_t1_g467958 * (-(dx)));
              double cse_p564_t42_g503572 = (cse_p67_t1_g467958 * dx);
              double h_x1_y1 = (cse_p564_t41_g503571 + (cse_p67_t1_g467957 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g463083 = (cse_p11_t1_g463079 * dz);
              double cse_p564_t29_g503559 = (cse_p11_t1_g463083 * (-(dx)));
              double cse_p564_t30_g503560 = (cse_p11_t1_g463083 * (-(dy)));
              double cse_p564_t31_g503561 = (cse_p11_t1_g463083 * dx);
              double cse_p564_t32_g503562 = (cse_p11_t1_g463083 * dy);
              double cse_p67_t1_g467960 = (cse_p11_t1_g463078 * dy);
              double cse_p564_t43_g503573 = (cse_p67_t1_g467960 * dy);
              double cse_p564_t44_g503574 = (cse_p67_t1_g467960 * dz);
              double cse_p116_t57_g472312 = (cse_p564_t30_g503560 + cse_p564_t44_g503574);
              double h_x1_z1 = (cse_p564_t29_g503559 + cse_p564_t40_g503570);
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g463081 = (cse_p11_t1_g463078 * (-(dx)));
              double cse_p564_t20_g503550 = (cse_p11_t1_g463081 * (-(dx)));
              double cse_p564_t21_g503551 = (cse_p11_t1_g463081 * (-(dy)));
              double cse_p564_t22_g503552 = (cse_p11_t1_g463081 * (-(dz)));
              double cse_p564_t23_g503553 = (cse_p11_t1_g463081 * dy);
              double cse_p564_t24_g503554 = (cse_p11_t1_g463081 * dz);
              double cse_p116_t52_g472307 = (cse_p564_t22_g503552 + cse_p564_t29_g503559);
              double cse_p116_t53_g472308 = (cse_p564_t23_g503553 + cse_p564_t42_g503572);
              double cse_p116_t54_g472309 = (cse_p564_t24_g503554 + cse_p564_t31_g503561);
              double cse_p67_t60_g467925 = (cse_p564_t21_g503551 + cse_p564_t41_g503571);
              double h_x1_x2 = ((cse_p11_t1_g463079 * (cse_p11_t49_g463021 + (-(cse_p11_t1_g463080)))) + (cse_p11_t1_g463081 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g463082 = (cse_p11_t1_g463078 * (-(dy)));
              double cse_p564_t25_g503555 = (cse_p11_t1_g463082 * (-(dy)));
              double cse_p564_t26_g503556 = (cse_p11_t1_g463082 * (-(dz)));
              double cse_p564_t27_g503557 = (cse_p11_t1_g463082 * dy);
              double cse_p564_t28_g503558 = (cse_p11_t1_g463082 * dz);
              double cse_p116_t50_g472305 = (cse_p564_t27_g503557 + (cse_p11_t1_g463079 * (cse_p11_t50_g463022 + (-(cse_p11_t1_g463080)))));
              double cse_p116_t55_g472310 = (cse_p564_t26_g503556 + cse_p564_t30_g503560);
              double cse_p116_t56_g472311 = (cse_p564_t28_g503558 + cse_p564_t32_g503562);
              double h_x1_y2 = (cse_p564_t42_g503572 + (cse_p11_t1_g463082 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g463084 = (cse_p11_t1_g463078 * (-(dz)));
              double cse_p564_t33_g503563 = (cse_p11_t1_g463084 * (-(dz)));
              double cse_p564_t34_g503564 = (cse_p11_t1_g463084 * dy);
              double cse_p564_t35_g503565 = (cse_p11_t1_g463084 * dz);
              double cse_p116_t51_g472306 = (cse_p564_t35_g503565 + (cse_p11_t1_g463079 * (cse_p11_t51_g463023 + (-(cse_p11_t1_g463080)))));
              double cse_p116_t58_g472313 = (cse_p564_t32_g503562 + cse_p564_t34_g503564);
              double h_x1_z2 = (cse_p564_t31_g503561 + (cse_p11_t1_g463084 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g463085 = (cse_p11_t1_g463079 * (cse_p11_t1_g463080 + (-(cse_p11_t50_g463022))));
              double cse_p116_t63_g472318 = (cse_p11_t1_g463085 + cse_p564_t25_g503555);
              double h_y1_y1 = (cse_p11_t1_g463085 + cse_p564_t43_g503573);
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = cse_p116_t57_g472312;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = cse_p116_t53_g472308;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = cse_p116_t50_g472305;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = cse_p116_t58_g472313;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g463086 = (cse_p11_t1_g463079 * (cse_p11_t1_g463080 + (-(cse_p11_t51_g463023))));
              double cse_p116_t64_g472319 = (cse_p11_t1_g463086 + cse_p564_t33_g503563);
              double h_z1_z1 = (cse_p11_t1_g463086 + (cse_p11_t1_g463078 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = cse_p116_t54_g472309;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = cse_p116_t56_g472311;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = cse_p116_t51_g472306;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p564_t1_g503608 + cse_p564_t20_g503550);
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = cse_p67_t60_g467925;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = cse_p116_t52_g472307;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = cse_p116_t63_g472318;
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = cse_p116_t55_g472310;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = cse_p116_t64_g472319;
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
                double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t49_g463021))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t2_invr503610 = (1.0 / (rr));
                double cse_p565_t1_invr503609 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g463080 = (rr * rr);
                double cse_p11_t1_g463079 = (cse_p565_t1_invr503609 * de_dr);
                double cse_p11_t1_g463078 = (cse_p565_t2_invr503610 * cse_p565_t2_invr503610 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g463078 * dx * dy) + (cse_p11_t1_g463079 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double cse_p67_t1_g467960 = (cse_p11_t1_g463078 * dy);
                double cse_p67_t1_g467958 = (cse_p11_t1_g463079 * dy);
                double cse_p67_t1_g467957 = (cse_p11_t1_g463078 * dx);
                double h_x1_z1 = ((cse_p11_t1_g463079 * dz * (-(dx))) + (cse_p67_t1_g467957 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g463083 = (cse_p11_t1_g463079 * dz);
                double h_x1_x2 = ((cse_p11_t1_g463079 * (cse_p11_t49_g463021 + (-(cse_p11_t1_g463080)))) + (cse_p67_t1_g467957 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g463081 = (cse_p11_t1_g463078 * (-(dx)));
                double h_x1_y2 = ((cse_p67_t1_g467957 * (-(dy))) + (cse_p67_t1_g467958 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g463082 = (cse_p11_t1_g463078 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g463083 * dx) + (cse_p67_t1_g467957 * (-(dz))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g463084 = (cse_p11_t1_g463078 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g463079 * (cse_p11_t1_g463080 + (-(cse_p11_t50_g463022)))) + (cse_p67_t1_g467960 * dy));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g463085 = (cse_p11_t1_g463079 * (cse_p11_t1_g463080 + (-(cse_p11_t50_g463022))));
                double h_y1_z1 = ((cse_p11_t1_g463083 * (-(dy))) + (cse_p67_t1_g467960 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g463081 * dy) + (cse_p67_t1_g467958 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g463079 * (cse_p11_t50_g463022 + (-(cse_p11_t1_g463080)))) + (cse_p11_t1_g463082 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g463083 * dy) + (cse_p11_t1_g463084 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g463078 * dz * dz) + (cse_p11_t1_g463079 * (cse_p11_t1_g463080 + (-(cse_p11_t51_g463023)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g463086 = (cse_p11_t1_g463079 * (cse_p11_t1_g463080 + (-(cse_p11_t51_g463023))));
                double h_z1_x2 = ((cse_p11_t1_g463081 * dz) + (cse_p11_t1_g463083 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g463082 * dz) + (cse_p11_t1_g463083 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g463079 * (cse_p11_t51_g463023 + (-(cse_p11_t1_g463080)))) + (cse_p11_t1_g463084 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g463079 * (cse_p11_t1_g463080 + (-(cse_p11_t49_g463021)))) + (cse_p11_t1_g463081 * (-(dx))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p564_t1_g503608 = (cse_p11_t1_g463079 * (cse_p11_t1_g463080 + (-(cse_p11_t49_g463021))));
                double h_x2_y2 = ((cse_p11_t1_g463081 * (-(dy))) + (cse_p67_t1_g467958 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double cse_p67_t1_g467961 = (cse_p67_t1_g467958 * (-(dx)));
                double h_x2_z2 = ((cse_p11_t1_g463081 * (-(dz))) + (cse_p11_t1_g463083 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g463085 + (cse_p11_t1_g463082 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g463082 * (-(dz))) + (cse_p11_t1_g463083 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g463086 + (cse_p11_t1_g463084 * (-(dz))));
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
