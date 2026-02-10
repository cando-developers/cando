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
      double cse_p1_t1_sqrt222255 = sqrt(r2);
      double cse_p1_t2_invsqrt222256 = (1.0 / (cse_p1_t1_sqrt222255));
      double rr = cse_p1_t1_sqrt222255;
      if ((rr <= (6.00000023841857900e-1 * term.sigma))) {
        /* !COND-1 */
        double cse_p528_t7_g235796 = (cse_p1_t2_invsqrt222256 * dx);
        double cse_p528_t8_g235797 = (cse_p1_t2_invsqrt222256 * dy);
        double cse_p528_t9_g235798 = (cse_p1_t2_invsqrt222256 * dz);
        double cse_p528_t16_g235805 = (-(term.epsilon));
        double cse_p12_t1_g222323 = (term.mij * params.rep_weight);
        double cse_p11_t9_g222267 = (cse_p12_t1_g222323 * cse_p528_t7_g235796);
        double cse_p11_t10_g222268 = (cse_p12_t1_g222323 * cse_p528_t8_g235797);
        double cse_p11_t11_g222269 = (cse_p12_t1_g222323 * cse_p528_t9_g235798);
        double fa_rep = (term.bij + (term.mij * rr));
        double cse_p528_t11_g235800 = (fa_rep * params.rep_weight);
        double fa_atr = cse_p528_t16_g235805;
        double cse_p13_t20_g222349 = (cse_p528_t11_g235800 + fa_atr);
        double de_dr = cse_p12_t1_g222323;
        double energy = cse_p13_t20_g222349;
        energy_added += energy;
        *energy_accumulate += energy;
        double g_x1 = cse_p11_t9_g222267;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = cse_p11_t10_g222268;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = cse_p11_t11_g222269;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = (-(cse_p11_t9_g222267));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (-(cse_p11_t10_g222268));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (-(cse_p11_t11_g222269));
        KernelGradientAcc(term.i3x2, 2, g_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.sigma)) {
          /* !COND-2 */
          double cse_p529_t1_invr235817 = (1.0 / (rr));
          double cse_p479_t11_g234475 = (cse_p529_t1_invr235817 * cse_p529_t1_invr235817);
          double cse_p479_t12_g234476 = (cse_p529_t1_invr235817 * term.sigma);
          double sr = cse_p479_t12_g234476;
          double cse_p528_t1_g235816 = (cse_p1_t2_invsqrt222256 * term.sigma * sr);
          double cse_p528_t12_g235801 = (sr * sr);
          double cse_p11_t1_g222287 = (cse_p479_t11_g234475 * cse_p528_t1_g235816);
          double cse_p528_t4_g235793 = (cse_p11_t1_g222287 * dx);
          double cse_p528_t5_g235794 = (cse_p11_t1_g222287 * dy);
          double cse_p528_t6_g235795 = (cse_p11_t1_g222287 * dz);
          double cse_p12_t19_g222315 = (-(cse_p528_t4_g235793));
          double cse_p12_t20_g222316 = (-(cse_p528_t5_g235794));
          double cse_p12_t21_g222317 = (-(cse_p528_t6_g235795));
          double sr2 = cse_p528_t12_g235801;
          double cse_p528_t13_g235802 = (sr2 * sr2);
          double cse_p13_t4_g222333 = (6.00000000000000000e+0 * cse_p528_t13_g235802 * term.epsilon);
          double cse_p13_t8_g222337 = (cse_p528_t13_g235802 * sr2);
          double sr6 = cse_p13_t8_g222337;
          double cse_p528_t2_g235791 = (-2.00000000000000000e+0 * sr6);
          double cse_p528_t3_g235792 = (1.20000000000000000e+1 * sr6);
          double cse_p528_t14_g235803 = (sr6 * sr6);
          double cse_p528_t15_g235804 = (-2.00000000000000000e+0 + sr6 + sr6);
          double sr12 = cse_p528_t14_g235803;
          double cse_p528_t1_g235790 = (-1.20000000000000000e+1 * sr12);
          double cse_p14_t15_g222376 = (cse_p528_t2_g235791 + sr12);
          double cse_p13_t1_g222361 = (cse_p14_t15_g222376 * term.epsilon);
          double cse_p13_t19_g222348 = (cse_p528_t1_g235790 + cse_p528_t3_g235792);
          double fa_rep_epsilon = cse_p13_t1_g222361;
          double fa_rep = (term.epsilon + fa_rep_epsilon);
          double fa_atr = (-(term.epsilon));
          double de_dr = (cse_p13_t19_g222348 * cse_p529_t1_invr235817 * term.epsilon * params.rep_weight);
          double cse_p12_t1_g222327 = (fa_rep * params.rep_weight);
          double energy = (cse_p12_t1_g222327 + fa_atr);
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p11_t1_g222290 = (cse_p13_t4_g222333 * cse_p528_t15_g235804 * params.rep_weight);
          double g_x1 = (cse_p11_t1_g222290 * cse_p12_t19_g222315);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t1_g222290 * cse_p12_t20_g222316);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t1_g222290 * cse_p12_t21_g222317);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p11_t1_g222290 * cse_p528_t4_g235793);
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p11_t1_g222290 * cse_p528_t5_g235794);
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p11_t1_g222290 * cse_p528_t6_g235795);
          KernelGradientAcc(term.i3x2, 2, g_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rswitch)) {
            /* !COND-3 */
            double cse_p529_t1_invr235817 = (1.0 / (rr));
            double cse_p1_t4_invr2222258 = (cse_p529_t1_invr235817 * cse_p529_t1_invr235817);
            double sr = (cse_p529_t1_invr235817 * term.sigma);
            double sr2 = (sr * sr);
            double cse_p12_t1_g222324 = (sr2 * sr2);
            double sr6 = (cse_p12_t1_g222324 * sr2);
            double sr12 = (sr6 * sr6);
            double cse_p12_t1_g222322 = (-2.00000000000000000e+0 * sr6);
            double cse_p13_t1_g222361 = (term.epsilon * (cse_p12_t1_g222322 + sr12));
            double fa_atr = cse_p13_t1_g222361;
            double cse_p12_t1_g222325 = (-1.20000000000000000e+1 * sr12);
            double cse_p12_t1_g222326 = (1.20000000000000000e+1 * sr6);
            double de_dr = (cse_p529_t1_invr235817 * term.epsilon * (cse_p12_t1_g222325 + cse_p12_t1_g222326));
            double energy = fa_atr;
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p528_t1_g235816 = (cse_p1_t2_invsqrt222256 * term.sigma * sr);
            double cse_p11_t1_g222287 = (cse_p1_t4_invr2222258 * cse_p528_t1_g235816);
            double cse_p11_t1_g222288 = (6.00000000000000000e+0 * cse_p12_t1_g222324 * term.epsilon);
            double cse_p11_t1_g222289 = (cse_p11_t1_g222288 * (-2.00000000000000000e+0 + sr6 + sr6));
            double cse_p11_t1_g222293 = (cse_p11_t1_g222287 * dx);
            double g_x1 = (cse_p11_t1_g222289 * (-(cse_p11_t1_g222293)));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p11_t1_g222294 = (cse_p11_t1_g222287 * dy);
            double g_y1 = (cse_p11_t1_g222289 * (-(cse_p11_t1_g222294)));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p11_t1_g222295 = (cse_p11_t1_g222287 * dz);
            double g_z1 = (cse_p11_t1_g222289 * (-(cse_p11_t1_g222295)));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (cse_p11_t1_g222289 * cse_p11_t1_g222293);
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (cse_p11_t1_g222289 * cse_p11_t1_g222294);
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (cse_p11_t1_g222289 * cse_p11_t1_g222295);
            KernelGradientAcc(term.i3x2, 2, g_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p11_t1_g222292 = (term.aa * rr * rr);
              double cse_p12_t4_g222300 = (3.00000000000000000e+0 * cse_p11_t1_g222292);
              double cse_p11_t1_g222291 = (term.bb * rr);
              double cse_p11_t21_g222279 = (term.cc + cse_p11_t1_g222291 + cse_p11_t1_g222291 + cse_p12_t4_g222300);
              double fa_atr = (term.dd + (term.cc * rr) + (rr * rr * (term.bb + (term.aa * rr))));
              double de_dr = (term.cc + cse_p12_t4_g222300 + (2.00000000000000000e+0 * cse_p11_t1_g222291));
              double energy = fa_atr;
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p11_t1_g222296 = (cse_p1_t2_invsqrt222256 * dx);
              double g_x1 = (cse_p11_t1_g222296 * cse_p11_t21_g222279);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double cse_p12_t1_g222320 = (cse_p1_t2_invsqrt222256 * dy);
              double g_y1 = (cse_p11_t21_g222279 * cse_p12_t1_g222320);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double cse_p12_t1_g222321 = (cse_p1_t2_invsqrt222256 * dz);
              double g_z1 = (cse_p11_t21_g222279 * cse_p12_t1_g222321);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p63_t1_g223652 = (-2.00000000000000000e+0 * cse_p11_t1_g222291);
              double cse_p63_t1_g223653 = (-3.00000000000000000e+0 * cse_p11_t1_g222292);
              double g_x2 = ((cse_p11_t1_g222296 * cse_p63_t1_g223652) + (cse_p11_t1_g222296 * cse_p63_t1_g223653) + (-((term.cc * cse_p11_t1_g222296))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = ((cse_p12_t1_g222320 * cse_p63_t1_g223652) + (cse_p12_t1_g222320 * cse_p63_t1_g223653) + (-((term.cc * cse_p12_t1_g222320))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = ((cse_p12_t1_g222321 * cse_p63_t1_g223652) + (cse_p12_t1_g222321 * cse_p63_t1_g223653) + (-((term.cc * cse_p12_t1_g222321))));
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
      double cse_p11_t49_g244951 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t50_g244952 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t51_g244953 = (dz * dz);
      double r2 = (cse_p11_t49_g244951 + cse_p11_t50_g244952 + cse_p11_t51_g244953);
      double cse_p1_t1_sqrt244899 = sqrt(r2);
      double cse_p1_t2_invsqrt244900 = (1.0 / (cse_p1_t1_sqrt244899));
      double rr = cse_p1_t1_sqrt244899;
      if ((rr <= (6.00000023841857900e-1 * term.sigma))) {
        /* !COND-1 */
        double cse_p564_t46_g285506 = (rr * rr);
        double cse_p564_t58_g285518 = (-(cse_p11_t49_g244951));
        double cse_p564_t59_g285519 = (-(cse_p11_t50_g244952));
        double cse_p564_t60_g285520 = (-(cse_p11_t51_g244953));
        double cse_p564_t61_g285521 = (-(dx));
        double cse_p564_t62_g285522 = (-(dy));
        double cse_p564_t63_g285523 = (-(dz));
        double cse_p564_t64_g285524 = (-(term.epsilon));
        double cse_p565_t1_invr285539 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr285540 = (1.0 / (rr));
        double cse_p564_t68_g285528 = (rr * rr * rr);
        double cse_p15_t89_g245459 = (-(cse_p564_t46_g285506));
        double cse_p14_t90_g245341 = (cse_p11_t49_g244951 + cse_p15_t89_g245459);
        double cse_p14_t92_g245343 = (cse_p11_t50_g244952 + cse_p15_t89_g245459);
        double cse_p14_t94_g245345 = (cse_p11_t51_g244953 + cse_p15_t89_g245459);
        double cse_p13_t116_g245234 = (cse_p564_t46_g285506 + cse_p564_t58_g285518);
        double cse_p13_t117_g245235 = (cse_p564_t46_g285506 + cse_p564_t59_g285519);
        double cse_p13_t118_g245236 = (cse_p564_t46_g285506 + cse_p564_t60_g285520);
        double cse_p11_t11_g244913 = (cse_p1_t2_invsqrt244900 * dx * term.mij * params.rep_weight);
        double cse_p11_t12_g244914 = (cse_p1_t2_invsqrt244900 * dy * term.mij * params.rep_weight);
        double cse_p11_t13_g244915 = (cse_p1_t2_invsqrt244900 * dz * term.mij * params.rep_weight);
        double fa_rep = (term.bij + (term.mij * rr));
        double cse_p564_t45_g285505 = (fa_rep * params.rep_weight);
        double fa_atr = cse_p564_t64_g285524;
        double cse_p60_t83_g249774 = (cse_p564_t45_g285505 + fa_atr);
        double de_dr = (term.mij * params.rep_weight);
        double cse_p515_t36_g281674 = (cse_p565_t1_invr285539 * de_dr);
        double cse_p13_t68_g245186 = (cse_p14_t90_g245341 * cse_p515_t36_g281674);
        double cse_p13_t69_g245187 = (cse_p14_t92_g245343 * cse_p515_t36_g281674);
        double cse_p13_t70_g245188 = (cse_p14_t94_g245345 * cse_p515_t36_g281674);
        double cse_p13_t71_g245189 = (cse_p515_t36_g281674 * dx * dy);
        double cse_p13_t72_g245190 = (cse_p515_t36_g281674 * cse_p564_t61_g285521 * dy);
        double cse_p12_t23_g245040 = (cse_p13_t116_g245234 * cse_p515_t36_g281674);
        double cse_p12_t24_g245041 = (cse_p13_t117_g245235 * cse_p515_t36_g281674);
        double cse_p12_t25_g245042 = (cse_p13_t118_g245236 * cse_p515_t36_g281674);
        double cse_p12_t28_g245045 = (cse_p515_t36_g281674 * dz);
        double cse_p13_t58_g245176 = (cse_p12_t28_g245045 * cse_p564_t61_g285521);
        double cse_p13_t59_g245177 = (cse_p12_t28_g245045 * cse_p564_t62_g285522);
        double cse_p13_t60_g245178 = (cse_p12_t28_g245045 * dx);
        double cse_p13_t61_g245179 = (cse_p12_t28_g245045 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p515_t37_g281675 = (cse_p565_t2_invr285540 * cse_p565_t2_invr285540 * d2e_dr2);
        double cse_p13_t62_g245180 = (cse_p515_t37_g281675 * dx * dx);
        double cse_p13_t63_g245181 = (cse_p515_t37_g281675 * dx * dy);
        double cse_p13_t64_g245182 = (cse_p515_t37_g281675 * dx * dz);
        double cse_p13_t65_g245183 = (cse_p515_t37_g281675 * dy * dy);
        double cse_p13_t66_g245184 = (cse_p515_t37_g281675 * dy * dz);
        double cse_p13_t67_g245185 = (cse_p515_t37_g281675 * dz * dz);
        double cse_p12_t11_g245028 = (cse_p515_t37_g281675 * cse_p564_t61_g285521);
        double cse_p13_t43_g245161 = (cse_p12_t11_g245028 * cse_p564_t61_g285521);
        double cse_p13_t44_g245162 = (cse_p12_t11_g245028 * cse_p564_t62_g285522);
        double cse_p13_t45_g245163 = (cse_p12_t11_g245028 * cse_p564_t63_g285523);
        double cse_p13_t46_g245164 = (cse_p12_t11_g245028 * dx);
        double cse_p13_t47_g245165 = (cse_p12_t11_g245028 * dy);
        double cse_p13_t48_g245166 = (cse_p12_t11_g245028 * dz);
        double cse_p12_t12_g245029 = (cse_p515_t37_g281675 * cse_p564_t62_g285522);
        double cse_p13_t49_g245167 = (cse_p12_t12_g245029 * cse_p564_t62_g285522);
        double cse_p13_t50_g245168 = (cse_p12_t12_g245029 * cse_p564_t63_g285523);
        double cse_p13_t51_g245169 = (cse_p12_t12_g245029 * dx);
        double cse_p13_t52_g245170 = (cse_p12_t12_g245029 * dy);
        double cse_p13_t53_g245171 = (cse_p12_t12_g245029 * dz);
        double cse_p12_t13_g245030 = (cse_p515_t37_g281675 * cse_p564_t63_g285523);
        double cse_p13_t54_g245172 = (cse_p12_t13_g245030 * cse_p564_t63_g285523);
        double cse_p13_t55_g245173 = (cse_p12_t13_g245030 * dx);
        double cse_p13_t56_g245174 = (cse_p12_t13_g245030 * dy);
        double cse_p13_t57_g245175 = (cse_p12_t13_g245030 * dz);
        double cse_p12_t63_g245080 = (cse_p13_t63_g245181 + cse_p13_t72_g245190);
        double cse_p12_t64_g245081 = (cse_p13_t58_g245176 + cse_p13_t64_g245182);
        double cse_p12_t65_g245082 = (cse_p13_t59_g245177 + cse_p13_t66_g245184);
        double cse_p12_t82_g245099 = (cse_p12_t24_g245041 + cse_p13_t65_g245183);
        double cse_p12_t84_g245101 = (cse_p12_t25_g245042 + cse_p13_t67_g245185);
        double cse_p12_t86_g245103 = (cse_p12_t23_g245040 + cse_p13_t62_g245180);
        double cse_p12_t58_g245075 = (cse_p13_t46_g245164 + cse_p13_t68_g245186);
        double cse_p12_t66_g245083 = (cse_p13_t47_g245165 + cse_p13_t71_g245189);
        double cse_p12_t68_g245085 = (cse_p13_t44_g245162 + cse_p13_t72_g245190);
        double cse_p12_t69_g245086 = (cse_p13_t45_g245163 + cse_p13_t58_g245176);
        double cse_p12_t70_g245087 = (cse_p13_t48_g245166 + cse_p13_t60_g245178);
        double cse_p12_t87_g245104 = (cse_p12_t23_g245040 + cse_p13_t43_g245161);
        double cse_p12_t59_g245076 = (cse_p13_t52_g245170 + cse_p13_t69_g245187);
        double cse_p12_t67_g245084 = (cse_p13_t51_g245169 + cse_p13_t71_g245189);
        double cse_p12_t71_g245088 = (cse_p13_t50_g245168 + cse_p13_t59_g245177);
        double cse_p12_t72_g245089 = (cse_p13_t53_g245171 + cse_p13_t61_g245179);
        double cse_p12_t83_g245100 = (cse_p12_t24_g245041 + cse_p13_t49_g245167);
        double cse_p12_t60_g245077 = (cse_p13_t57_g245175 + cse_p13_t70_g245188);
        double cse_p12_t73_g245090 = (cse_p13_t55_g245173 + cse_p13_t60_g245178);
        double cse_p12_t74_g245091 = (cse_p13_t56_g245174 + cse_p13_t61_g245179);
        double cse_p12_t85_g245102 = (cse_p12_t25_g245042 + cse_p13_t54_g245172);
        double energy = cse_p60_t83_g249774;
        energy_added += energy;
        *energy_accumulate += energy;
        double g_x1 = cse_p11_t11_g244913;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = cse_p11_t12_g244914;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = cse_p11_t13_g244915;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = (-(cse_p11_t11_g244913));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (-(cse_p11_t12_g244914));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (-(cse_p11_t13_g244915));
        KernelGradientAcc(term.i3x2, 2, g_z2);
        double h_x1_x1 = cse_p12_t86_g245103;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t63_g245080;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t64_g245081;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t58_g245075;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t67_g245084;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t73_g245090;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t82_g245099;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t65_g245082;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t66_g245083;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t59_g245076;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t74_g245091;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t84_g245101;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t70_g245087;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t72_g245089;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t60_g245077;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t87_g245104;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t68_g245085;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t69_g245086;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t83_g245100;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t71_g245088;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t85_g245102;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.sigma)) {
          /* !COND-2 */
          double cse_p11_t1_g245010 = (rr * rr);
          double cse_p564_t51_g285511 = (cse_p11_t1_g245010 + (-(cse_p11_t49_g244951)));
          double cse_p564_t52_g285512 = (cse_p11_t1_g245010 + (-(cse_p11_t50_g244952)));
          double cse_p564_t53_g285513 = (cse_p11_t1_g245010 + (-(cse_p11_t51_g244953)));
          double cse_p564_t54_g285514 = (cse_p11_t49_g244951 + (-(cse_p11_t1_g245010)));
          double cse_p564_t55_g285515 = (cse_p11_t50_g244952 + (-(cse_p11_t1_g245010)));
          double cse_p564_t56_g285516 = (cse_p11_t51_g244953 + (-(cse_p11_t1_g245010)));
          double cse_p564_t57_g285517 = (-(cse_p11_t1_g245010));
          double cse_p565_t3_invr285541 = (1.0 / (cse_p11_t1_g245010));
          double cse_p565_t2_invr285540 = (1.0 / (rr));
          double cse_p564_t38_g285498 = (cse_p565_t2_invr285540 * cse_p565_t2_invr285540);
          double cse_p564_t39_g285499 = (cse_p565_t2_invr285540 * term.sigma);
          double sr = cse_p564_t39_g285499;
          double cse_p564_t47_g285507 = (sr * sr);
          double cse_p67_t1_g249892 = (cse_p1_t2_invsqrt244900 * cse_p564_t38_g285498 * term.sigma * sr);
          double cse_p12_t8_g245025 = (cse_p67_t1_g249892 * dx);
          double cse_p12_t9_g245026 = (cse_p67_t1_g249892 * dy);
          double cse_p12_t10_g245027 = (cse_p67_t1_g249892 * dz);
          double cse_p11_t91_g244993 = (-(cse_p12_t8_g245025));
          double cse_p11_t92_g244994 = (-(cse_p12_t9_g245026));
          double cse_p11_t93_g244995 = (-(cse_p12_t10_g245027));
          double sr2 = cse_p564_t47_g285507;
          double cse_p564_t48_g285508 = (sr2 * sr2 * sr2);
          double sr6 = cse_p564_t48_g285508;
          double cse_p564_t2_g285462 = (-2.00000000000000000e+0 * sr6);
          double cse_p564_t3_g285463 = (-4.20000000000000040e+1 * sr6);
          double cse_p564_t4_g285464 = (1.20000000000000000e+1 * sr6);
          double cse_p564_t49_g285509 = (sr6 * sr6);
          double cse_p564_t50_g285510 = (-2.00000000000000000e+0 + sr6 + sr6);
          double sr12 = cse_p564_t49_g285509;
          double cse_p564_t1_g285461 = (-1.20000000000000000e+1 * sr12);
          double cse_p564_t5_g285465 = (1.56000000000000000e+2 * sr12);
          double cse_p67_t78_g249873 = (cse_p564_t2_g285462 + sr12);
          double cse_p60_t47_g249738 = (cse_p67_t78_g249873 * term.epsilon);
          double cse_p60_t56_g249747 = (cse_p564_t1_g285461 + cse_p564_t4_g285464);
          double cse_p60_t57_g249748 = (cse_p564_t3_g285463 + cse_p564_t5_g285465);
          double fa_rep_epsilon = cse_p60_t47_g249738;
          double fa_rep = (term.epsilon + fa_rep_epsilon);
          double fa_atr = (-(term.epsilon));
          double de_dr = (cse_p565_t2_invr285540 * cse_p60_t56_g249747 * term.epsilon * params.rep_weight);
          double d2e_dr2 = (cse_p565_t3_invr285541 * cse_p60_t57_g249748 * term.epsilon * params.rep_weight);
          double cse_p564_t37_g285497 = (cse_p565_t2_invr285540 * cse_p565_t2_invr285540 * d2e_dr2);
          double energy = (fa_atr + (fa_rep * params.rep_weight));
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p67_t1_g249889 = (6.00000000000000000e+0 * cse_p564_t50_g285510 * term.epsilon * params.rep_weight * sr2 * sr2);
          double g_x1 = (cse_p11_t91_g244993 * cse_p67_t1_g249889);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t92_g244994 * cse_p67_t1_g249889);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t93_g244995 * cse_p67_t1_g249889);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p67_t1_g249889 * cse_p67_t1_g249892 * dx);
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p67_t1_g249889 * cse_p67_t1_g249892 * dy);
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p67_t1_g249889 * cse_p67_t1_g249892 * dz);
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double cse_p68_t8_g249903 = (cse_p564_t37_g285497 * (-(dx)));
          double cse_p68_t9_g249904 = (cse_p564_t37_g285497 * (-(dy)));
          double cse_p68_t10_g249905 = (cse_p564_t37_g285497 * (-(dz)));
          double cse_p68_t11_g249906 = (cse_p564_t37_g285497 * dx);
          double cse_p68_t12_g249907 = (cse_p564_t37_g285497 * dy);
          double cse_p68_t13_g249908 = (cse_p564_t37_g285497 * dz * dz);
          double cse_p68_t45_g249940 = (cse_p68_t11_g249906 * dx);
          double cse_p68_t46_g249941 = (cse_p68_t11_g249906 * dy);
          double cse_p68_t47_g249942 = (cse_p68_t11_g249906 * dz);
          double cse_p68_t50_g249945 = (cse_p68_t12_g249907 * dy);
          double cse_p68_t51_g249946 = (cse_p68_t12_g249907 * dz);
          double cse_p14_t47_g245298 = (cse_p68_t8_g249903 * (-(dx)));
          double cse_p14_t48_g245299 = (cse_p68_t8_g249903 * (-(dy)));
          double cse_p14_t49_g245300 = (cse_p68_t8_g249903 * (-(dz)));
          double cse_p14_t50_g245301 = (cse_p68_t8_g249903 * dx);
          double cse_p14_t51_g245302 = (cse_p68_t8_g249903 * dy);
          double cse_p14_t52_g245303 = (cse_p68_t8_g249903 * dz);
          double cse_p14_t53_g245304 = (cse_p68_t9_g249904 * (-(dy)));
          double cse_p14_t54_g245305 = (cse_p68_t9_g249904 * (-(dz)));
          double cse_p14_t55_g245306 = (cse_p68_t9_g249904 * dx);
          double cse_p14_t56_g245307 = (cse_p68_t9_g249904 * dy);
          double cse_p14_t57_g245308 = (cse_p68_t9_g249904 * dz);
          double cse_p14_t58_g245309 = (cse_p68_t10_g249905 * (-(dz)));
          double cse_p14_t59_g245310 = (cse_p68_t10_g249905 * dx);
          double cse_p14_t60_g245311 = (cse_p68_t10_g249905 * dy);
          double cse_p14_t61_g245312 = (cse_p68_t10_g249905 * dz);
          double cse_p565_t1_invr285539 = (1.0 / ((rr * rr * rr)));
          double cse_p564_t36_g285496 = (cse_p565_t1_invr285539 * de_dr);
          double cse_p68_t20_g249915 = (cse_p564_t36_g285496 * dy);
          double cse_p68_t21_g249916 = (cse_p564_t36_g285496 * dz);
          double cse_p68_t48_g249943 = (cse_p68_t20_g249915 * (-(dx)));
          double cse_p68_t49_g249944 = (cse_p68_t20_g249915 * dx);
          double cse_p60_t15_g249706 = (cse_p564_t36_g285496 * (cse_p11_t49_g244951 + cse_p564_t57_g285517));
          double cse_p60_t16_g249707 = (cse_p564_t36_g285496 * (cse_p11_t50_g244952 + cse_p564_t57_g285517));
          double cse_p60_t17_g249708 = (cse_p564_t36_g285496 * (cse_p11_t51_g244953 + cse_p564_t57_g285517));
          double cse_p60_t18_g249709 = (cse_p564_t36_g285496 * cse_p564_t51_g285511);
          double cse_p60_t19_g249710 = (cse_p564_t36_g285496 * cse_p564_t52_g285512);
          double cse_p60_t20_g249711 = (cse_p564_t36_g285496 * cse_p564_t53_g285513);
          double cse_p60_t58_g249749 = (cse_p68_t46_g249941 + cse_p68_t48_g249943);
          double cse_p14_t43_g245294 = (cse_p68_t21_g249916 * (-(dx)));
          double cse_p14_t44_g245295 = (cse_p68_t21_g249916 * (-(dy)));
          double cse_p14_t45_g245296 = (cse_p68_t21_g249916 * dx);
          double cse_p14_t46_g245297 = (cse_p68_t21_g249916 * dy);
          double cse_p13_t114_g245232 = (cse_p60_t18_g249709 + cse_p68_t45_g249940);
          double h_x1_x1 = cse_p13_t114_g245232;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t58_g249749;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t89_g245207 = (cse_p14_t43_g245294 + cse_p68_t47_g249942);
          double cse_p13_t90_g245208 = (cse_p14_t44_g245295 + cse_p68_t51_g249946);
          double h_x1_z1 = cse_p13_t89_g245207;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t83_g245201 = (cse_p14_t50_g245301 + cse_p60_t15_g249706);
          double cse_p13_t91_g245209 = (cse_p14_t51_g245302 + cse_p68_t49_g249944);
          double cse_p13_t93_g245211 = (cse_p14_t48_g245299 + cse_p68_t48_g249943);
          double cse_p13_t94_g245212 = (cse_p14_t43_g245294 + cse_p14_t49_g245300);
          double cse_p13_t95_g245213 = (cse_p14_t45_g245296 + cse_p14_t52_g245303);
          double cse_p13_t115_g245233 = (cse_p14_t47_g245298 + cse_p60_t18_g249709);
          double h_x1_x2 = cse_p13_t83_g245201;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t84_g245202 = (cse_p14_t56_g245307 + cse_p60_t16_g249707);
          double cse_p13_t92_g245210 = (cse_p14_t55_g245306 + cse_p68_t49_g249944);
          double cse_p13_t96_g245214 = (cse_p14_t44_g245295 + cse_p14_t54_g245305);
          double cse_p13_t97_g245215 = (cse_p14_t46_g245297 + cse_p14_t57_g245308);
          double h_x1_y2 = cse_p13_t92_g245210;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t85_g245203 = (cse_p14_t61_g245312 + cse_p60_t17_g249708);
          double cse_p13_t98_g245216 = (cse_p14_t45_g245296 + cse_p14_t59_g245310);
          double cse_p13_t99_g245217 = (cse_p14_t46_g245297 + cse_p14_t60_g245311);
          double h_x1_z2 = cse_p13_t98_g245216;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t110_g245228 = (cse_p60_t19_g249710 + cse_p68_t50_g249945);
          double cse_p13_t111_g245229 = (cse_p14_t53_g245304 + cse_p60_t19_g249710);
          double h_y1_y1 = cse_p13_t110_g245228;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t90_g245208;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t91_g245209;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t84_g245202;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t99_g245217;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t112_g245230 = (cse_p60_t20_g249711 + cse_p68_t13_g249908);
          double cse_p13_t113_g245231 = (cse_p14_t58_g245309 + cse_p60_t20_g249711);
          double h_z1_z1 = cse_p13_t112_g245230;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t95_g245213;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t97_g245215;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t85_g245203;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t115_g245233;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t93_g245211;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t94_g245212;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t111_g245229;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t96_g245214;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t113_g245231;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rswitch)) {
            /* !COND-3 */
            double cse_p565_t2_invr285540 = (1.0 / (rr));
            double cse_p1_t4_invr2244902 = (cse_p565_t2_invr285540 * cse_p565_t2_invr285540);
            double sr = (cse_p565_t2_invr285540 * term.sigma);
            double sr2 = (sr * sr);
            double sr6 = (sr2 * sr2 * sr2);
            double sr12 = (sr6 * sr6);
            double fa_atr = (term.epsilon * (sr12 + (-2.00000000000000000e+0 * sr6)));
            double de_dr = (cse_p565_t2_invr285540 * term.epsilon * ((-1.20000000000000000e+1 * sr12) + (1.20000000000000000e+1 * sr6)));
            double cse_p11_t1_g245010 = (rr * rr);
            double cse_p565_t3_invr285541 = (1.0 / (cse_p11_t1_g245010));
            double d2e_dr2 = (cse_p565_t3_invr285541 * term.epsilon * ((-4.20000000000000040e+1 * sr6) + (1.56000000000000000e+2 * sr12)));
            double energy = fa_atr;
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p67_t1_g249886 = (cse_p1_t2_invsqrt244900 * cse_p1_t4_invr2244902 * term.epsilon * term.sigma * sr2 * sr2 * sr * (-2.00000000000000000e+0 + sr6 + sr6));
            double g_x1 = (-6.00000000000000000e+0 * cse_p67_t1_g249886 * dx);
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double g_y1 = (-6.00000000000000000e+0 * cse_p67_t1_g249886 * dy);
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double g_z1 = (-6.00000000000000000e+0 * cse_p67_t1_g249886 * dz);
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (6.00000000000000000e+0 * cse_p67_t1_g249886 * dx);
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (6.00000000000000000e+0 * cse_p67_t1_g249886 * dy);
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (6.00000000000000000e+0 * cse_p67_t1_g249886 * dz);
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double cse_p11_t1_g245008 = (cse_p565_t2_invr285540 * cse_p565_t2_invr285540 * d2e_dr2);
            double cse_p564_t6_g285466 = (cse_p11_t1_g245008 * (-(dx)));
            double cse_p564_t7_g285467 = (cse_p11_t1_g245008 * (-(dy)));
            double cse_p564_t8_g285468 = (cse_p11_t1_g245008 * (-(dz)));
            double cse_p564_t9_g285469 = (cse_p11_t1_g245008 * dx);
            double cse_p564_t10_g285470 = (cse_p11_t1_g245008 * dy);
            double cse_p564_t11_g285471 = (cse_p11_t1_g245008 * dz * dz);
            double cse_p70_t45_g250134 = (cse_p564_t9_g285469 * dx);
            double cse_p70_t46_g250135 = (cse_p564_t9_g285469 * dy);
            double cse_p565_t1_invr285539 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g245009 = (cse_p565_t1_invr285539 * de_dr);
            double cse_p564_t1_g285538 = (cse_p11_t1_g245009 * (cse_p11_t1_g245010 + (-(cse_p11_t49_g244951))));
            double cse_p564_t13_g285473 = (cse_p11_t1_g245009 * (cse_p11_t1_g245010 + (-(cse_p11_t50_g244952))));
            double cse_p564_t14_g285474 = (cse_p11_t1_g245009 * (cse_p11_t1_g245010 + (-(cse_p11_t51_g244953))));
            double cse_p564_t15_g285475 = (cse_p11_t1_g245009 * (cse_p11_t49_g244951 + (-(cse_p11_t1_g245010))));
            double cse_p564_t16_g285476 = (cse_p11_t1_g245009 * (cse_p11_t50_g244952 + (-(cse_p11_t1_g245010))));
            double cse_p564_t17_g285477 = (cse_p11_t1_g245009 * (cse_p11_t51_g244953 + (-(cse_p11_t1_g245010))));
            double cse_p564_t18_g285478 = (cse_p11_t1_g245009 * dy);
            double cse_p564_t19_g285479 = (cse_p11_t1_g245009 * dz);
            double cse_p69_t45_g250038 = (cse_p564_t9_g285469 * dz);
            double cse_p69_t74_g250067 = (cse_p564_t1_g285538 + cse_p70_t45_g250134);
            double cse_p69_t46_g250039 = (cse_p564_t18_g285478 * (-(dx)));
            double cse_p69_t47_g250040 = (cse_p564_t18_g285478 * dx);
            double cse_p69_t73_g250066 = (cse_p69_t46_g250039 + cse_p70_t46_g250135);
            double h_x1_x1 = cse_p69_t74_g250067;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = cse_p69_t73_g250066;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p69_t31_g250024 = (cse_p564_t19_g285479 * (-(dx)));
            double cse_p69_t32_g250025 = (cse_p564_t19_g285479 * (-(dy)));
            double cse_p69_t33_g250026 = (cse_p564_t19_g285479 * dx);
            double cse_p69_t34_g250027 = (cse_p564_t19_g285479 * dy);
            double cse_p60_t59_g249750 = (cse_p69_t31_g250024 + cse_p69_t45_g250038);
            double cse_p69_t48_g250041 = (cse_p564_t10_g285470 * dy);
            double cse_p69_t49_g250042 = (cse_p564_t10_g285470 * dz);
            double cse_p60_t60_g249751 = (cse_p69_t32_g250025 + cse_p69_t49_g250042);
            double h_x1_z1 = cse_p60_t59_g249750;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p69_t20_g250013 = (cse_p564_t6_g285466 * (-(dx)));
            double cse_p69_t21_g250014 = (cse_p564_t6_g285466 * (-(dy)));
            double cse_p69_t22_g250015 = (cse_p564_t6_g285466 * (-(dz)));
            double cse_p69_t23_g250016 = (cse_p564_t6_g285466 * dx);
            double cse_p69_t24_g250017 = (cse_p564_t6_g285466 * dy);
            double cse_p69_t25_g250018 = (cse_p564_t6_g285466 * dz);
            double cse_p60_t53_g249744 = (cse_p564_t15_g285475 + cse_p69_t23_g250016);
            double cse_p60_t61_g249752 = (cse_p69_t24_g250017 + cse_p69_t47_g250040);
            double cse_p60_t63_g249754 = (cse_p69_t21_g250014 + cse_p69_t46_g250039);
            double cse_p60_t64_g249755 = (cse_p69_t22_g250015 + cse_p69_t31_g250024);
            double cse_p60_t65_g249756 = (cse_p69_t25_g250018 + cse_p69_t33_g250026);
            double cse_p60_t82_g249773 = (cse_p564_t1_g285538 + cse_p69_t20_g250013);
            double h_x1_x2 = cse_p60_t53_g249744;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p69_t26_g250019 = (cse_p564_t7_g285467 * (-(dy)));
            double cse_p69_t27_g250020 = (cse_p564_t7_g285467 * (-(dz)));
            double cse_p69_t28_g250021 = (cse_p564_t7_g285467 * dx);
            double cse_p69_t29_g250022 = (cse_p564_t7_g285467 * dy);
            double cse_p69_t30_g250023 = (cse_p564_t7_g285467 * dz);
            double cse_p60_t54_g249745 = (cse_p564_t16_g285476 + cse_p69_t29_g250022);
            double cse_p60_t62_g249753 = (cse_p69_t28_g250021 + cse_p69_t47_g250040);
            double cse_p60_t66_g249757 = (cse_p69_t27_g250020 + cse_p69_t32_g250025);
            double cse_p60_t67_g249758 = (cse_p69_t30_g250023 + cse_p69_t34_g250027);
            double h_x1_y2 = cse_p60_t62_g249753;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p69_t35_g250028 = (cse_p564_t8_g285468 * (-(dz)));
            double cse_p69_t36_g250029 = (cse_p564_t8_g285468 * dx);
            double cse_p69_t37_g250030 = (cse_p564_t8_g285468 * dy);
            double cse_p69_t38_g250031 = (cse_p564_t8_g285468 * dz);
            double cse_p60_t55_g249746 = (cse_p564_t17_g285477 + cse_p69_t38_g250031);
            double cse_p60_t68_g249759 = (cse_p69_t33_g250026 + cse_p69_t36_g250029);
            double cse_p60_t69_g249760 = (cse_p69_t34_g250027 + cse_p69_t37_g250030);
            double h_x1_z2 = cse_p60_t68_g249759;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p67_t70_g249865 = (cse_p564_t13_g285473 + cse_p69_t48_g250041);
            double cse_p60_t78_g249769 = (cse_p564_t13_g285473 + cse_p69_t26_g250019);
            double h_y1_y1 = cse_p67_t70_g249865;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t60_g249751;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t61_g249752;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t54_g249745;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t69_g249760;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p67_t72_g249867 = (cse_p564_t11_g285471 + cse_p564_t14_g285474);
            double cse_p60_t80_g249771 = (cse_p564_t14_g285474 + cse_p69_t35_g250028);
            double h_z1_z1 = cse_p67_t72_g249867;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t65_g249756;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t67_g249758;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t55_g249746;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t82_g249773;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t63_g249754;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t64_g249755;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t78_g249769;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t66_g249757;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t80_g249771;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p68_t2_g249897 = (-2.00000000000000000e+0 * term.bb * rr);
              double cse_p12_t7_g245024 = (term.bb * rr);
              double cse_p11_t1_g245010 = (rr * rr);
              double cse_p68_t4_g249899 = (-3.00000000000000000e+0 * term.aa * cse_p11_t1_g245010);
              double cse_p67_t56_g249851 = (cse_p68_t2_g249897 + cse_p68_t4_g249899);
              double cse_p12_t6_g245023 = (3.00000000000000000e+0 * term.aa * cse_p11_t1_g245010);
              double cse_p11_t88_g244990 = (term.cc + cse_p12_t6_g245023 + cse_p12_t7_g245024 + cse_p12_t7_g245024);
              double fa_atr = (term.dd + (term.cc * rr) + (cse_p11_t1_g245010 * (term.bb + (term.aa * rr))));
              double de_dr = (term.cc + cse_p12_t6_g245023 + (2.00000000000000000e+0 * term.bb * rr));
              double d2e_dr2 = ((2.00000000000000000e+0 * term.bb) + (6.00000000000000000e+0 * term.aa * rr));
              double energy = fa_atr;
              energy_added += energy;
              *energy_accumulate += energy;
              double g_x1 = (cse_p11_t88_g244990 * cse_p1_t2_invsqrt244900 * dx);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t88_g244990 * cse_p1_t2_invsqrt244900 * dy);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t88_g244990 * cse_p1_t2_invsqrt244900 * dz);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = ((cse_p1_t2_invsqrt244900 * cse_p67_t56_g249851 * dx) + (-((term.cc * cse_p1_t2_invsqrt244900 * dx))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = ((cse_p1_t2_invsqrt244900 * cse_p67_t56_g249851 * dy) + (-((term.cc * cse_p1_t2_invsqrt244900 * dy))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = ((cse_p1_t2_invsqrt244900 * cse_p67_t56_g249851 * dz) + (-((term.cc * cse_p1_t2_invsqrt244900 * dz))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double cse_p565_t2_invr285540 = (1.0 / (rr));
              double cse_p11_t1_g245008 = (cse_p565_t2_invr285540 * cse_p565_t2_invr285540 * d2e_dr2);
              double cse_p565_t1_invr285539 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g245009 = (cse_p565_t1_invr285539 * de_dr);
              double cse_p564_t1_g285538 = (cse_p11_t1_g245009 * (cse_p11_t1_g245010 + (-(cse_p11_t49_g244951))));
              double cse_p67_t1_g249887 = (cse_p11_t1_g245008 * dx);
              double cse_p564_t40_g285500 = (cse_p67_t1_g249887 * dz);
              double h_x1_x1 = (cse_p564_t1_g285538 + (cse_p67_t1_g249887 * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double cse_p67_t1_g249888 = (cse_p11_t1_g245009 * dy);
              double cse_p564_t41_g285501 = (cse_p67_t1_g249888 * (-(dx)));
              double cse_p564_t42_g285502 = (cse_p67_t1_g249888 * dx);
              double h_x1_y1 = (cse_p564_t41_g285501 + (cse_p67_t1_g249887 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g245013 = (cse_p11_t1_g245009 * dz);
              double cse_p564_t29_g285489 = (cse_p11_t1_g245013 * (-(dx)));
              double cse_p564_t30_g285490 = (cse_p11_t1_g245013 * (-(dy)));
              double cse_p564_t31_g285491 = (cse_p11_t1_g245013 * dx);
              double cse_p564_t32_g285492 = (cse_p11_t1_g245013 * dy);
              double cse_p67_t1_g249890 = (cse_p11_t1_g245008 * dy);
              double cse_p564_t43_g285503 = (cse_p67_t1_g249890 * dy);
              double cse_p564_t44_g285504 = (cse_p67_t1_g249890 * dz);
              double cse_p116_t57_g254242 = (cse_p564_t30_g285490 + cse_p564_t44_g285504);
              double h_x1_z1 = (cse_p564_t29_g285489 + cse_p564_t40_g285500);
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g245011 = (cse_p11_t1_g245008 * (-(dx)));
              double cse_p564_t20_g285480 = (cse_p11_t1_g245011 * (-(dx)));
              double cse_p564_t21_g285481 = (cse_p11_t1_g245011 * (-(dy)));
              double cse_p564_t22_g285482 = (cse_p11_t1_g245011 * (-(dz)));
              double cse_p564_t23_g285483 = (cse_p11_t1_g245011 * dy);
              double cse_p564_t24_g285484 = (cse_p11_t1_g245011 * dz);
              double cse_p116_t52_g254237 = (cse_p564_t22_g285482 + cse_p564_t29_g285489);
              double cse_p116_t53_g254238 = (cse_p564_t23_g285483 + cse_p564_t42_g285502);
              double cse_p116_t54_g254239 = (cse_p564_t24_g285484 + cse_p564_t31_g285491);
              double cse_p67_t60_g249855 = (cse_p564_t21_g285481 + cse_p564_t41_g285501);
              double h_x1_x2 = ((cse_p11_t1_g245009 * (cse_p11_t49_g244951 + (-(cse_p11_t1_g245010)))) + (cse_p11_t1_g245011 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g245012 = (cse_p11_t1_g245008 * (-(dy)));
              double cse_p564_t25_g285485 = (cse_p11_t1_g245012 * (-(dy)));
              double cse_p564_t26_g285486 = (cse_p11_t1_g245012 * (-(dz)));
              double cse_p564_t27_g285487 = (cse_p11_t1_g245012 * dy);
              double cse_p564_t28_g285488 = (cse_p11_t1_g245012 * dz);
              double cse_p116_t50_g254235 = (cse_p564_t27_g285487 + (cse_p11_t1_g245009 * (cse_p11_t50_g244952 + (-(cse_p11_t1_g245010)))));
              double cse_p116_t55_g254240 = (cse_p564_t26_g285486 + cse_p564_t30_g285490);
              double cse_p116_t56_g254241 = (cse_p564_t28_g285488 + cse_p564_t32_g285492);
              double h_x1_y2 = (cse_p564_t42_g285502 + (cse_p11_t1_g245012 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g245014 = (cse_p11_t1_g245008 * (-(dz)));
              double cse_p564_t33_g285493 = (cse_p11_t1_g245014 * (-(dz)));
              double cse_p564_t34_g285494 = (cse_p11_t1_g245014 * dy);
              double cse_p564_t35_g285495 = (cse_p11_t1_g245014 * dz);
              double cse_p116_t51_g254236 = (cse_p564_t35_g285495 + (cse_p11_t1_g245009 * (cse_p11_t51_g244953 + (-(cse_p11_t1_g245010)))));
              double cse_p116_t58_g254243 = (cse_p564_t32_g285492 + cse_p564_t34_g285494);
              double h_x1_z2 = (cse_p564_t31_g285491 + (cse_p11_t1_g245014 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g245015 = (cse_p11_t1_g245009 * (cse_p11_t1_g245010 + (-(cse_p11_t50_g244952))));
              double cse_p116_t63_g254248 = (cse_p11_t1_g245015 + cse_p564_t25_g285485);
              double h_y1_y1 = (cse_p11_t1_g245015 + cse_p564_t43_g285503);
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = cse_p116_t57_g254242;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = cse_p116_t53_g254238;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = cse_p116_t50_g254235;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = cse_p116_t58_g254243;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g245016 = (cse_p11_t1_g245009 * (cse_p11_t1_g245010 + (-(cse_p11_t51_g244953))));
              double cse_p116_t64_g254249 = (cse_p11_t1_g245016 + cse_p564_t33_g285493);
              double h_z1_z1 = (cse_p11_t1_g245016 + (cse_p11_t1_g245008 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = cse_p116_t54_g254239;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = cse_p116_t56_g254241;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = cse_p116_t51_g254236;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p564_t1_g285538 + cse_p564_t20_g285480);
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = cse_p67_t60_g249855;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = cse_p116_t52_g254237;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = cse_p116_t63_g254248;
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = cse_p116_t55_g254240;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = cse_p116_t64_g254249;
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
                double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t49_g244951))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t2_invr285540 = (1.0 / (rr));
                double cse_p565_t1_invr285539 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g245010 = (rr * rr);
                double cse_p11_t1_g245009 = (cse_p565_t1_invr285539 * de_dr);
                double cse_p11_t1_g245008 = (cse_p565_t2_invr285540 * cse_p565_t2_invr285540 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g245008 * dx * dy) + (cse_p11_t1_g245009 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double cse_p67_t1_g249890 = (cse_p11_t1_g245008 * dy);
                double cse_p67_t1_g249888 = (cse_p11_t1_g245009 * dy);
                double cse_p67_t1_g249887 = (cse_p11_t1_g245008 * dx);
                double h_x1_z1 = ((cse_p11_t1_g245009 * dz * (-(dx))) + (cse_p67_t1_g249887 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g245013 = (cse_p11_t1_g245009 * dz);
                double h_x1_x2 = ((cse_p11_t1_g245009 * (cse_p11_t49_g244951 + (-(cse_p11_t1_g245010)))) + (cse_p67_t1_g249887 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g245011 = (cse_p11_t1_g245008 * (-(dx)));
                double h_x1_y2 = ((cse_p67_t1_g249887 * (-(dy))) + (cse_p67_t1_g249888 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g245012 = (cse_p11_t1_g245008 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g245013 * dx) + (cse_p67_t1_g249887 * (-(dz))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g245014 = (cse_p11_t1_g245008 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g245009 * (cse_p11_t1_g245010 + (-(cse_p11_t50_g244952)))) + (cse_p67_t1_g249890 * dy));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g245015 = (cse_p11_t1_g245009 * (cse_p11_t1_g245010 + (-(cse_p11_t50_g244952))));
                double h_y1_z1 = ((cse_p11_t1_g245013 * (-(dy))) + (cse_p67_t1_g249890 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g245011 * dy) + (cse_p67_t1_g249888 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g245009 * (cse_p11_t50_g244952 + (-(cse_p11_t1_g245010)))) + (cse_p11_t1_g245012 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g245013 * dy) + (cse_p11_t1_g245014 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g245008 * dz * dz) + (cse_p11_t1_g245009 * (cse_p11_t1_g245010 + (-(cse_p11_t51_g244953)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g245016 = (cse_p11_t1_g245009 * (cse_p11_t1_g245010 + (-(cse_p11_t51_g244953))));
                double h_z1_x2 = ((cse_p11_t1_g245011 * dz) + (cse_p11_t1_g245013 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g245012 * dz) + (cse_p11_t1_g245013 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g245009 * (cse_p11_t51_g244953 + (-(cse_p11_t1_g245010)))) + (cse_p11_t1_g245014 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g245009 * (cse_p11_t1_g245010 + (-(cse_p11_t49_g244951)))) + (cse_p11_t1_g245011 * (-(dx))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p564_t1_g285538 = (cse_p11_t1_g245009 * (cse_p11_t1_g245010 + (-(cse_p11_t49_g244951))));
                double h_x2_y2 = ((cse_p11_t1_g245011 * (-(dy))) + (cse_p67_t1_g249888 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double cse_p67_t1_g249891 = (cse_p67_t1_g249888 * (-(dx)));
                double h_x2_z2 = ((cse_p11_t1_g245011 * (-(dz))) + (cse_p11_t1_g245013 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g245015 + (cse_p11_t1_g245012 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g245012 * (-(dz))) + (cse_p11_t1_g245013 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g245016 + (cse_p11_t1_g245014 * (-(dz))));
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
