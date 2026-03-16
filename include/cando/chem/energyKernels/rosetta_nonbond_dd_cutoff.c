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
      double cse_p1_t1_sqrt67241 = sqrt(r2);
      double cse_p1_t2_invsqrt67242 = (1.0 / (cse_p1_t1_sqrt67241));
      double rr = cse_p1_t1_sqrt67241;
      if ((rr <= (6.00000023841857900e-1 * term.sigma))) {
        /* !COND-1 */
        double cse_p528_t7_g80782 = (cse_p1_t2_invsqrt67242 * dx);
        double cse_p528_t8_g80783 = (cse_p1_t2_invsqrt67242 * dy);
        double cse_p528_t9_g80784 = (cse_p1_t2_invsqrt67242 * dz);
        double cse_p528_t16_g80791 = (-(term.epsilon));
        double cse_p12_t1_g67309 = (term.mij * params.rep_weight);
        double cse_p11_t9_g67253 = (cse_p12_t1_g67309 * cse_p528_t7_g80782);
        double cse_p11_t10_g67254 = (cse_p12_t1_g67309 * cse_p528_t8_g80783);
        double cse_p11_t11_g67255 = (cse_p12_t1_g67309 * cse_p528_t9_g80784);
        double fa_rep = (term.bij + (term.mij * rr));
        double cse_p528_t11_g80786 = (fa_rep * params.rep_weight);
        double fa_atr = cse_p528_t16_g80791;
        double cse_p13_t20_g67335 = (cse_p528_t11_g80786 + fa_atr);
        double de_dr = cse_p12_t1_g67309;
        double energy = cse_p13_t20_g67335;
        energy_added += energy;
        *energy_accumulate += energy;
        double g_x1 = cse_p11_t9_g67253;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = cse_p11_t10_g67254;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = cse_p11_t11_g67255;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = (-(cse_p11_t9_g67253));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (-(cse_p11_t10_g67254));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (-(cse_p11_t11_g67255));
        KernelGradientAcc(term.i3x2, 2, g_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.sigma)) {
          /* !COND-2 */
          double cse_p529_t1_invr80803 = (1.0 / (rr));
          double cse_p479_t11_g79461 = (cse_p529_t1_invr80803 * cse_p529_t1_invr80803);
          double cse_p479_t12_g79462 = (cse_p529_t1_invr80803 * term.sigma);
          double sr = cse_p479_t12_g79462;
          double cse_p528_t1_g80802 = (cse_p1_t2_invsqrt67242 * term.sigma * sr);
          double cse_p528_t12_g80787 = (sr * sr);
          double cse_p11_t1_g67273 = (cse_p479_t11_g79461 * cse_p528_t1_g80802);
          double cse_p528_t4_g80779 = (cse_p11_t1_g67273 * dx);
          double cse_p528_t5_g80780 = (cse_p11_t1_g67273 * dy);
          double cse_p528_t6_g80781 = (cse_p11_t1_g67273 * dz);
          double cse_p12_t19_g67301 = (-(cse_p528_t4_g80779));
          double cse_p12_t20_g67302 = (-(cse_p528_t5_g80780));
          double cse_p12_t21_g67303 = (-(cse_p528_t6_g80781));
          double sr2 = cse_p528_t12_g80787;
          double cse_p528_t13_g80788 = (sr2 * sr2);
          double cse_p13_t4_g67319 = (6.00000000000000000e+0 * cse_p528_t13_g80788 * term.epsilon);
          double cse_p13_t8_g67323 = (cse_p528_t13_g80788 * sr2);
          double sr6 = cse_p13_t8_g67323;
          double cse_p528_t2_g80777 = (-2.00000000000000000e+0 * sr6);
          double cse_p528_t3_g80778 = (1.20000000000000000e+1 * sr6);
          double cse_p528_t14_g80789 = (sr6 * sr6);
          double cse_p528_t15_g80790 = (-2.00000000000000000e+0 + sr6 + sr6);
          double sr12 = cse_p528_t14_g80789;
          double cse_p528_t1_g80776 = (-1.20000000000000000e+1 * sr12);
          double cse_p14_t15_g67362 = (cse_p528_t2_g80777 + sr12);
          double cse_p13_t1_g67347 = (cse_p14_t15_g67362 * term.epsilon);
          double cse_p13_t19_g67334 = (cse_p528_t1_g80776 + cse_p528_t3_g80778);
          double fa_rep_epsilon = cse_p13_t1_g67347;
          double fa_rep = (term.epsilon + fa_rep_epsilon);
          double fa_atr = (-(term.epsilon));
          double de_dr = (cse_p13_t19_g67334 * cse_p529_t1_invr80803 * term.epsilon * params.rep_weight);
          double cse_p12_t1_g67313 = (fa_rep * params.rep_weight);
          double energy = (cse_p12_t1_g67313 + fa_atr);
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p11_t1_g67276 = (cse_p13_t4_g67319 * cse_p528_t15_g80790 * params.rep_weight);
          double g_x1 = (cse_p11_t1_g67276 * cse_p12_t19_g67301);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t1_g67276 * cse_p12_t20_g67302);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t1_g67276 * cse_p12_t21_g67303);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p11_t1_g67276 * cse_p528_t4_g80779);
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p11_t1_g67276 * cse_p528_t5_g80780);
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p11_t1_g67276 * cse_p528_t6_g80781);
          KernelGradientAcc(term.i3x2, 2, g_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rswitch)) {
            /* !COND-3 */
            double cse_p529_t1_invr80803 = (1.0 / (rr));
            double cse_p1_t4_invr267244 = (cse_p529_t1_invr80803 * cse_p529_t1_invr80803);
            double sr = (cse_p529_t1_invr80803 * term.sigma);
            double sr2 = (sr * sr);
            double cse_p12_t1_g67310 = (sr2 * sr2);
            double sr6 = (cse_p12_t1_g67310 * sr2);
            double sr12 = (sr6 * sr6);
            double cse_p12_t1_g67308 = (-2.00000000000000000e+0 * sr6);
            double cse_p13_t1_g67347 = (term.epsilon * (cse_p12_t1_g67308 + sr12));
            double fa_atr = cse_p13_t1_g67347;
            double cse_p12_t1_g67311 = (-1.20000000000000000e+1 * sr12);
            double cse_p12_t1_g67312 = (1.20000000000000000e+1 * sr6);
            double de_dr = (cse_p529_t1_invr80803 * term.epsilon * (cse_p12_t1_g67311 + cse_p12_t1_g67312));
            double energy = fa_atr;
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p528_t1_g80802 = (cse_p1_t2_invsqrt67242 * term.sigma * sr);
            double cse_p11_t1_g67273 = (cse_p1_t4_invr267244 * cse_p528_t1_g80802);
            double cse_p11_t1_g67274 = (6.00000000000000000e+0 * cse_p12_t1_g67310 * term.epsilon);
            double cse_p11_t1_g67275 = (cse_p11_t1_g67274 * (-2.00000000000000000e+0 + sr6 + sr6));
            double cse_p11_t1_g67279 = (cse_p11_t1_g67273 * dx);
            double g_x1 = (cse_p11_t1_g67275 * (-(cse_p11_t1_g67279)));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p11_t1_g67280 = (cse_p11_t1_g67273 * dy);
            double g_y1 = (cse_p11_t1_g67275 * (-(cse_p11_t1_g67280)));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p11_t1_g67281 = (cse_p11_t1_g67273 * dz);
            double g_z1 = (cse_p11_t1_g67275 * (-(cse_p11_t1_g67281)));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (cse_p11_t1_g67275 * cse_p11_t1_g67279);
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (cse_p11_t1_g67275 * cse_p11_t1_g67280);
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (cse_p11_t1_g67275 * cse_p11_t1_g67281);
            KernelGradientAcc(term.i3x2, 2, g_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p11_t1_g67278 = (term.aa * rr * rr);
              double cse_p12_t4_g67286 = (3.00000000000000000e+0 * cse_p11_t1_g67278);
              double cse_p11_t1_g67277 = (term.bb * rr);
              double cse_p11_t21_g67265 = (term.cc + cse_p11_t1_g67277 + cse_p11_t1_g67277 + cse_p12_t4_g67286);
              double fa_atr = (term.dd + (term.cc * rr) + (rr * rr * (term.bb + (term.aa * rr))));
              double de_dr = (term.cc + cse_p12_t4_g67286 + (2.00000000000000000e+0 * cse_p11_t1_g67277));
              double energy = fa_atr;
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p11_t1_g67282 = (cse_p1_t2_invsqrt67242 * dx);
              double g_x1 = (cse_p11_t1_g67282 * cse_p11_t21_g67265);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double cse_p12_t1_g67306 = (cse_p1_t2_invsqrt67242 * dy);
              double g_y1 = (cse_p11_t21_g67265 * cse_p12_t1_g67306);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double cse_p12_t1_g67307 = (cse_p1_t2_invsqrt67242 * dz);
              double g_z1 = (cse_p11_t21_g67265 * cse_p12_t1_g67307);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p63_t1_g68638 = (-2.00000000000000000e+0 * cse_p11_t1_g67277);
              double cse_p63_t1_g68639 = (-3.00000000000000000e+0 * cse_p11_t1_g67278);
              double g_x2 = ((cse_p11_t1_g67282 * cse_p63_t1_g68638) + (cse_p11_t1_g67282 * cse_p63_t1_g68639) + (-((term.cc * cse_p11_t1_g67282))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = ((cse_p12_t1_g67306 * cse_p63_t1_g68638) + (cse_p12_t1_g67306 * cse_p63_t1_g68639) + (-((term.cc * cse_p12_t1_g67306))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = ((cse_p12_t1_g67307 * cse_p63_t1_g68638) + (cse_p12_t1_g67307 * cse_p63_t1_g68639) + (-((term.cc * cse_p12_t1_g67307))));
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
      double cse_p11_t49_g89937 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t50_g89938 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t51_g89939 = (dz * dz);
      double r2 = (cse_p11_t49_g89937 + cse_p11_t50_g89938 + cse_p11_t51_g89939);
      double cse_p1_t1_sqrt89885 = sqrt(r2);
      double cse_p1_t2_invsqrt89886 = (1.0 / (cse_p1_t1_sqrt89885));
      double rr = cse_p1_t1_sqrt89885;
      if ((rr <= (6.00000023841857900e-1 * term.sigma))) {
        /* !COND-1 */
        double cse_p564_t46_g130492 = (rr * rr);
        double cse_p564_t58_g130504 = (-(cse_p11_t49_g89937));
        double cse_p564_t59_g130505 = (-(cse_p11_t50_g89938));
        double cse_p564_t60_g130506 = (-(cse_p11_t51_g89939));
        double cse_p564_t61_g130507 = (-(dx));
        double cse_p564_t62_g130508 = (-(dy));
        double cse_p564_t63_g130509 = (-(dz));
        double cse_p564_t64_g130510 = (-(term.epsilon));
        double cse_p565_t1_invr130525 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr130526 = (1.0 / (rr));
        double cse_p564_t68_g130514 = (rr * rr * rr);
        double cse_p15_t89_g90445 = (-(cse_p564_t46_g130492));
        double cse_p14_t90_g90327 = (cse_p11_t49_g89937 + cse_p15_t89_g90445);
        double cse_p14_t92_g90329 = (cse_p11_t50_g89938 + cse_p15_t89_g90445);
        double cse_p14_t94_g90331 = (cse_p11_t51_g89939 + cse_p15_t89_g90445);
        double cse_p13_t116_g90220 = (cse_p564_t46_g130492 + cse_p564_t58_g130504);
        double cse_p13_t117_g90221 = (cse_p564_t46_g130492 + cse_p564_t59_g130505);
        double cse_p13_t118_g90222 = (cse_p564_t46_g130492 + cse_p564_t60_g130506);
        double cse_p11_t11_g89899 = (cse_p1_t2_invsqrt89886 * dx * term.mij * params.rep_weight);
        double cse_p11_t12_g89900 = (cse_p1_t2_invsqrt89886 * dy * term.mij * params.rep_weight);
        double cse_p11_t13_g89901 = (cse_p1_t2_invsqrt89886 * dz * term.mij * params.rep_weight);
        double fa_rep = (term.bij + (term.mij * rr));
        double cse_p564_t45_g130491 = (fa_rep * params.rep_weight);
        double fa_atr = cse_p564_t64_g130510;
        double cse_p60_t83_g94760 = (cse_p564_t45_g130491 + fa_atr);
        double de_dr = (term.mij * params.rep_weight);
        double cse_p515_t36_g126660 = (cse_p565_t1_invr130525 * de_dr);
        double cse_p13_t68_g90172 = (cse_p14_t90_g90327 * cse_p515_t36_g126660);
        double cse_p13_t69_g90173 = (cse_p14_t92_g90329 * cse_p515_t36_g126660);
        double cse_p13_t70_g90174 = (cse_p14_t94_g90331 * cse_p515_t36_g126660);
        double cse_p13_t71_g90175 = (cse_p515_t36_g126660 * dx * dy);
        double cse_p13_t72_g90176 = (cse_p515_t36_g126660 * cse_p564_t61_g130507 * dy);
        double cse_p12_t23_g90026 = (cse_p13_t116_g90220 * cse_p515_t36_g126660);
        double cse_p12_t24_g90027 = (cse_p13_t117_g90221 * cse_p515_t36_g126660);
        double cse_p12_t25_g90028 = (cse_p13_t118_g90222 * cse_p515_t36_g126660);
        double cse_p12_t28_g90031 = (cse_p515_t36_g126660 * dz);
        double cse_p13_t58_g90162 = (cse_p12_t28_g90031 * cse_p564_t61_g130507);
        double cse_p13_t59_g90163 = (cse_p12_t28_g90031 * cse_p564_t62_g130508);
        double cse_p13_t60_g90164 = (cse_p12_t28_g90031 * dx);
        double cse_p13_t61_g90165 = (cse_p12_t28_g90031 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p515_t37_g126661 = (cse_p565_t2_invr130526 * cse_p565_t2_invr130526 * d2e_dr2);
        double cse_p13_t62_g90166 = (cse_p515_t37_g126661 * dx * dx);
        double cse_p13_t63_g90167 = (cse_p515_t37_g126661 * dx * dy);
        double cse_p13_t64_g90168 = (cse_p515_t37_g126661 * dx * dz);
        double cse_p13_t65_g90169 = (cse_p515_t37_g126661 * dy * dy);
        double cse_p13_t66_g90170 = (cse_p515_t37_g126661 * dy * dz);
        double cse_p13_t67_g90171 = (cse_p515_t37_g126661 * dz * dz);
        double cse_p12_t11_g90014 = (cse_p515_t37_g126661 * cse_p564_t61_g130507);
        double cse_p13_t43_g90147 = (cse_p12_t11_g90014 * cse_p564_t61_g130507);
        double cse_p13_t44_g90148 = (cse_p12_t11_g90014 * cse_p564_t62_g130508);
        double cse_p13_t45_g90149 = (cse_p12_t11_g90014 * cse_p564_t63_g130509);
        double cse_p13_t46_g90150 = (cse_p12_t11_g90014 * dx);
        double cse_p13_t47_g90151 = (cse_p12_t11_g90014 * dy);
        double cse_p13_t48_g90152 = (cse_p12_t11_g90014 * dz);
        double cse_p12_t12_g90015 = (cse_p515_t37_g126661 * cse_p564_t62_g130508);
        double cse_p13_t49_g90153 = (cse_p12_t12_g90015 * cse_p564_t62_g130508);
        double cse_p13_t50_g90154 = (cse_p12_t12_g90015 * cse_p564_t63_g130509);
        double cse_p13_t51_g90155 = (cse_p12_t12_g90015 * dx);
        double cse_p13_t52_g90156 = (cse_p12_t12_g90015 * dy);
        double cse_p13_t53_g90157 = (cse_p12_t12_g90015 * dz);
        double cse_p12_t13_g90016 = (cse_p515_t37_g126661 * cse_p564_t63_g130509);
        double cse_p13_t54_g90158 = (cse_p12_t13_g90016 * cse_p564_t63_g130509);
        double cse_p13_t55_g90159 = (cse_p12_t13_g90016 * dx);
        double cse_p13_t56_g90160 = (cse_p12_t13_g90016 * dy);
        double cse_p13_t57_g90161 = (cse_p12_t13_g90016 * dz);
        double cse_p12_t63_g90066 = (cse_p13_t63_g90167 + cse_p13_t72_g90176);
        double cse_p12_t64_g90067 = (cse_p13_t58_g90162 + cse_p13_t64_g90168);
        double cse_p12_t65_g90068 = (cse_p13_t59_g90163 + cse_p13_t66_g90170);
        double cse_p12_t82_g90085 = (cse_p12_t24_g90027 + cse_p13_t65_g90169);
        double cse_p12_t84_g90087 = (cse_p12_t25_g90028 + cse_p13_t67_g90171);
        double cse_p12_t86_g90089 = (cse_p12_t23_g90026 + cse_p13_t62_g90166);
        double cse_p12_t58_g90061 = (cse_p13_t46_g90150 + cse_p13_t68_g90172);
        double cse_p12_t66_g90069 = (cse_p13_t47_g90151 + cse_p13_t71_g90175);
        double cse_p12_t68_g90071 = (cse_p13_t44_g90148 + cse_p13_t72_g90176);
        double cse_p12_t69_g90072 = (cse_p13_t45_g90149 + cse_p13_t58_g90162);
        double cse_p12_t70_g90073 = (cse_p13_t48_g90152 + cse_p13_t60_g90164);
        double cse_p12_t87_g90090 = (cse_p12_t23_g90026 + cse_p13_t43_g90147);
        double cse_p12_t59_g90062 = (cse_p13_t52_g90156 + cse_p13_t69_g90173);
        double cse_p12_t67_g90070 = (cse_p13_t51_g90155 + cse_p13_t71_g90175);
        double cse_p12_t71_g90074 = (cse_p13_t50_g90154 + cse_p13_t59_g90163);
        double cse_p12_t72_g90075 = (cse_p13_t53_g90157 + cse_p13_t61_g90165);
        double cse_p12_t83_g90086 = (cse_p12_t24_g90027 + cse_p13_t49_g90153);
        double cse_p12_t60_g90063 = (cse_p13_t57_g90161 + cse_p13_t70_g90174);
        double cse_p12_t73_g90076 = (cse_p13_t55_g90159 + cse_p13_t60_g90164);
        double cse_p12_t74_g90077 = (cse_p13_t56_g90160 + cse_p13_t61_g90165);
        double cse_p12_t85_g90088 = (cse_p12_t25_g90028 + cse_p13_t54_g90158);
        double energy = cse_p60_t83_g94760;
        energy_added += energy;
        *energy_accumulate += energy;
        double g_x1 = cse_p11_t11_g89899;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = cse_p11_t12_g89900;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = cse_p11_t13_g89901;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = (-(cse_p11_t11_g89899));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (-(cse_p11_t12_g89900));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (-(cse_p11_t13_g89901));
        KernelGradientAcc(term.i3x2, 2, g_z2);
        double h_x1_x1 = cse_p12_t86_g90089;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t63_g90066;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t64_g90067;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t58_g90061;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t67_g90070;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t73_g90076;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t82_g90085;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t65_g90068;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t66_g90069;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t59_g90062;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t74_g90077;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t84_g90087;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t70_g90073;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t72_g90075;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t60_g90063;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t87_g90090;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t68_g90071;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t69_g90072;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t83_g90086;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t71_g90074;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t85_g90088;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.sigma)) {
          /* !COND-2 */
          double cse_p11_t1_g89996 = (rr * rr);
          double cse_p564_t51_g130497 = (cse_p11_t1_g89996 + (-(cse_p11_t49_g89937)));
          double cse_p564_t52_g130498 = (cse_p11_t1_g89996 + (-(cse_p11_t50_g89938)));
          double cse_p564_t53_g130499 = (cse_p11_t1_g89996 + (-(cse_p11_t51_g89939)));
          double cse_p564_t54_g130500 = (cse_p11_t49_g89937 + (-(cse_p11_t1_g89996)));
          double cse_p564_t55_g130501 = (cse_p11_t50_g89938 + (-(cse_p11_t1_g89996)));
          double cse_p564_t56_g130502 = (cse_p11_t51_g89939 + (-(cse_p11_t1_g89996)));
          double cse_p564_t57_g130503 = (-(cse_p11_t1_g89996));
          double cse_p565_t3_invr130527 = (1.0 / (cse_p11_t1_g89996));
          double cse_p565_t2_invr130526 = (1.0 / (rr));
          double cse_p564_t38_g130484 = (cse_p565_t2_invr130526 * cse_p565_t2_invr130526);
          double cse_p564_t39_g130485 = (cse_p565_t2_invr130526 * term.sigma);
          double sr = cse_p564_t39_g130485;
          double cse_p564_t47_g130493 = (sr * sr);
          double cse_p67_t1_g94878 = (cse_p1_t2_invsqrt89886 * cse_p564_t38_g130484 * term.sigma * sr);
          double cse_p12_t8_g90011 = (cse_p67_t1_g94878 * dx);
          double cse_p12_t9_g90012 = (cse_p67_t1_g94878 * dy);
          double cse_p12_t10_g90013 = (cse_p67_t1_g94878 * dz);
          double cse_p11_t91_g89979 = (-(cse_p12_t8_g90011));
          double cse_p11_t92_g89980 = (-(cse_p12_t9_g90012));
          double cse_p11_t93_g89981 = (-(cse_p12_t10_g90013));
          double sr2 = cse_p564_t47_g130493;
          double cse_p564_t48_g130494 = (sr2 * sr2 * sr2);
          double sr6 = cse_p564_t48_g130494;
          double cse_p564_t2_g130448 = (-2.00000000000000000e+0 * sr6);
          double cse_p564_t3_g130449 = (-4.20000000000000040e+1 * sr6);
          double cse_p564_t4_g130450 = (1.20000000000000000e+1 * sr6);
          double cse_p564_t49_g130495 = (sr6 * sr6);
          double cse_p564_t50_g130496 = (-2.00000000000000000e+0 + sr6 + sr6);
          double sr12 = cse_p564_t49_g130495;
          double cse_p564_t1_g130447 = (-1.20000000000000000e+1 * sr12);
          double cse_p564_t5_g130451 = (1.56000000000000000e+2 * sr12);
          double cse_p67_t78_g94859 = (cse_p564_t2_g130448 + sr12);
          double cse_p60_t47_g94724 = (cse_p67_t78_g94859 * term.epsilon);
          double cse_p60_t56_g94733 = (cse_p564_t1_g130447 + cse_p564_t4_g130450);
          double cse_p60_t57_g94734 = (cse_p564_t3_g130449 + cse_p564_t5_g130451);
          double fa_rep_epsilon = cse_p60_t47_g94724;
          double fa_rep = (term.epsilon + fa_rep_epsilon);
          double fa_atr = (-(term.epsilon));
          double de_dr = (cse_p565_t2_invr130526 * cse_p60_t56_g94733 * term.epsilon * params.rep_weight);
          double d2e_dr2 = (cse_p565_t3_invr130527 * cse_p60_t57_g94734 * term.epsilon * params.rep_weight);
          double cse_p564_t37_g130483 = (cse_p565_t2_invr130526 * cse_p565_t2_invr130526 * d2e_dr2);
          double energy = (fa_atr + (fa_rep * params.rep_weight));
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p67_t1_g94875 = (6.00000000000000000e+0 * cse_p564_t50_g130496 * term.epsilon * params.rep_weight * sr2 * sr2);
          double g_x1 = (cse_p11_t91_g89979 * cse_p67_t1_g94875);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t92_g89980 * cse_p67_t1_g94875);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t93_g89981 * cse_p67_t1_g94875);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p67_t1_g94875 * cse_p67_t1_g94878 * dx);
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p67_t1_g94875 * cse_p67_t1_g94878 * dy);
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p67_t1_g94875 * cse_p67_t1_g94878 * dz);
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double cse_p68_t8_g94889 = (cse_p564_t37_g130483 * (-(dx)));
          double cse_p68_t9_g94890 = (cse_p564_t37_g130483 * (-(dy)));
          double cse_p68_t10_g94891 = (cse_p564_t37_g130483 * (-(dz)));
          double cse_p68_t11_g94892 = (cse_p564_t37_g130483 * dx);
          double cse_p68_t12_g94893 = (cse_p564_t37_g130483 * dy);
          double cse_p68_t13_g94894 = (cse_p564_t37_g130483 * dz * dz);
          double cse_p68_t45_g94926 = (cse_p68_t11_g94892 * dx);
          double cse_p68_t46_g94927 = (cse_p68_t11_g94892 * dy);
          double cse_p68_t47_g94928 = (cse_p68_t11_g94892 * dz);
          double cse_p68_t50_g94931 = (cse_p68_t12_g94893 * dy);
          double cse_p68_t51_g94932 = (cse_p68_t12_g94893 * dz);
          double cse_p14_t47_g90284 = (cse_p68_t8_g94889 * (-(dx)));
          double cse_p14_t48_g90285 = (cse_p68_t8_g94889 * (-(dy)));
          double cse_p14_t49_g90286 = (cse_p68_t8_g94889 * (-(dz)));
          double cse_p14_t50_g90287 = (cse_p68_t8_g94889 * dx);
          double cse_p14_t51_g90288 = (cse_p68_t8_g94889 * dy);
          double cse_p14_t52_g90289 = (cse_p68_t8_g94889 * dz);
          double cse_p14_t53_g90290 = (cse_p68_t9_g94890 * (-(dy)));
          double cse_p14_t54_g90291 = (cse_p68_t9_g94890 * (-(dz)));
          double cse_p14_t55_g90292 = (cse_p68_t9_g94890 * dx);
          double cse_p14_t56_g90293 = (cse_p68_t9_g94890 * dy);
          double cse_p14_t57_g90294 = (cse_p68_t9_g94890 * dz);
          double cse_p14_t58_g90295 = (cse_p68_t10_g94891 * (-(dz)));
          double cse_p14_t59_g90296 = (cse_p68_t10_g94891 * dx);
          double cse_p14_t60_g90297 = (cse_p68_t10_g94891 * dy);
          double cse_p14_t61_g90298 = (cse_p68_t10_g94891 * dz);
          double cse_p565_t1_invr130525 = (1.0 / ((rr * rr * rr)));
          double cse_p564_t36_g130482 = (cse_p565_t1_invr130525 * de_dr);
          double cse_p68_t20_g94901 = (cse_p564_t36_g130482 * dy);
          double cse_p68_t21_g94902 = (cse_p564_t36_g130482 * dz);
          double cse_p68_t48_g94929 = (cse_p68_t20_g94901 * (-(dx)));
          double cse_p68_t49_g94930 = (cse_p68_t20_g94901 * dx);
          double cse_p60_t15_g94692 = (cse_p564_t36_g130482 * (cse_p11_t49_g89937 + cse_p564_t57_g130503));
          double cse_p60_t16_g94693 = (cse_p564_t36_g130482 * (cse_p11_t50_g89938 + cse_p564_t57_g130503));
          double cse_p60_t17_g94694 = (cse_p564_t36_g130482 * (cse_p11_t51_g89939 + cse_p564_t57_g130503));
          double cse_p60_t18_g94695 = (cse_p564_t36_g130482 * cse_p564_t51_g130497);
          double cse_p60_t19_g94696 = (cse_p564_t36_g130482 * cse_p564_t52_g130498);
          double cse_p60_t20_g94697 = (cse_p564_t36_g130482 * cse_p564_t53_g130499);
          double cse_p60_t58_g94735 = (cse_p68_t46_g94927 + cse_p68_t48_g94929);
          double cse_p14_t43_g90280 = (cse_p68_t21_g94902 * (-(dx)));
          double cse_p14_t44_g90281 = (cse_p68_t21_g94902 * (-(dy)));
          double cse_p14_t45_g90282 = (cse_p68_t21_g94902 * dx);
          double cse_p14_t46_g90283 = (cse_p68_t21_g94902 * dy);
          double cse_p13_t114_g90218 = (cse_p60_t18_g94695 + cse_p68_t45_g94926);
          double h_x1_x1 = cse_p13_t114_g90218;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t58_g94735;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t89_g90193 = (cse_p14_t43_g90280 + cse_p68_t47_g94928);
          double cse_p13_t90_g90194 = (cse_p14_t44_g90281 + cse_p68_t51_g94932);
          double h_x1_z1 = cse_p13_t89_g90193;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t83_g90187 = (cse_p14_t50_g90287 + cse_p60_t15_g94692);
          double cse_p13_t91_g90195 = (cse_p14_t51_g90288 + cse_p68_t49_g94930);
          double cse_p13_t93_g90197 = (cse_p14_t48_g90285 + cse_p68_t48_g94929);
          double cse_p13_t94_g90198 = (cse_p14_t43_g90280 + cse_p14_t49_g90286);
          double cse_p13_t95_g90199 = (cse_p14_t45_g90282 + cse_p14_t52_g90289);
          double cse_p13_t115_g90219 = (cse_p14_t47_g90284 + cse_p60_t18_g94695);
          double h_x1_x2 = cse_p13_t83_g90187;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t84_g90188 = (cse_p14_t56_g90293 + cse_p60_t16_g94693);
          double cse_p13_t92_g90196 = (cse_p14_t55_g90292 + cse_p68_t49_g94930);
          double cse_p13_t96_g90200 = (cse_p14_t44_g90281 + cse_p14_t54_g90291);
          double cse_p13_t97_g90201 = (cse_p14_t46_g90283 + cse_p14_t57_g90294);
          double h_x1_y2 = cse_p13_t92_g90196;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t85_g90189 = (cse_p14_t61_g90298 + cse_p60_t17_g94694);
          double cse_p13_t98_g90202 = (cse_p14_t45_g90282 + cse_p14_t59_g90296);
          double cse_p13_t99_g90203 = (cse_p14_t46_g90283 + cse_p14_t60_g90297);
          double h_x1_z2 = cse_p13_t98_g90202;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t110_g90214 = (cse_p60_t19_g94696 + cse_p68_t50_g94931);
          double cse_p13_t111_g90215 = (cse_p14_t53_g90290 + cse_p60_t19_g94696);
          double h_y1_y1 = cse_p13_t110_g90214;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t90_g90194;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t91_g90195;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t84_g90188;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t99_g90203;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t112_g90216 = (cse_p60_t20_g94697 + cse_p68_t13_g94894);
          double cse_p13_t113_g90217 = (cse_p14_t58_g90295 + cse_p60_t20_g94697);
          double h_z1_z1 = cse_p13_t112_g90216;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t95_g90199;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t97_g90201;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t85_g90189;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t115_g90219;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t93_g90197;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t94_g90198;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t111_g90215;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t96_g90200;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t113_g90217;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rswitch)) {
            /* !COND-3 */
            double cse_p565_t2_invr130526 = (1.0 / (rr));
            double cse_p1_t4_invr289888 = (cse_p565_t2_invr130526 * cse_p565_t2_invr130526);
            double sr = (cse_p565_t2_invr130526 * term.sigma);
            double sr2 = (sr * sr);
            double sr6 = (sr2 * sr2 * sr2);
            double sr12 = (sr6 * sr6);
            double fa_atr = (term.epsilon * (sr12 + (-2.00000000000000000e+0 * sr6)));
            double de_dr = (cse_p565_t2_invr130526 * term.epsilon * ((-1.20000000000000000e+1 * sr12) + (1.20000000000000000e+1 * sr6)));
            double cse_p11_t1_g89996 = (rr * rr);
            double cse_p565_t3_invr130527 = (1.0 / (cse_p11_t1_g89996));
            double d2e_dr2 = (cse_p565_t3_invr130527 * term.epsilon * ((-4.20000000000000040e+1 * sr6) + (1.56000000000000000e+2 * sr12)));
            double energy = fa_atr;
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p67_t1_g94872 = (cse_p1_t2_invsqrt89886 * cse_p1_t4_invr289888 * term.epsilon * term.sigma * sr2 * sr2 * sr * (-2.00000000000000000e+0 + sr6 + sr6));
            double g_x1 = (-6.00000000000000000e+0 * cse_p67_t1_g94872 * dx);
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double g_y1 = (-6.00000000000000000e+0 * cse_p67_t1_g94872 * dy);
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double g_z1 = (-6.00000000000000000e+0 * cse_p67_t1_g94872 * dz);
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (6.00000000000000000e+0 * cse_p67_t1_g94872 * dx);
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (6.00000000000000000e+0 * cse_p67_t1_g94872 * dy);
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (6.00000000000000000e+0 * cse_p67_t1_g94872 * dz);
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double cse_p11_t1_g89994 = (cse_p565_t2_invr130526 * cse_p565_t2_invr130526 * d2e_dr2);
            double cse_p564_t6_g130452 = (cse_p11_t1_g89994 * (-(dx)));
            double cse_p564_t7_g130453 = (cse_p11_t1_g89994 * (-(dy)));
            double cse_p564_t8_g130454 = (cse_p11_t1_g89994 * (-(dz)));
            double cse_p564_t9_g130455 = (cse_p11_t1_g89994 * dx);
            double cse_p564_t10_g130456 = (cse_p11_t1_g89994 * dy);
            double cse_p564_t11_g130457 = (cse_p11_t1_g89994 * dz * dz);
            double cse_p70_t45_g95120 = (cse_p564_t9_g130455 * dx);
            double cse_p70_t46_g95121 = (cse_p564_t9_g130455 * dy);
            double cse_p565_t1_invr130525 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g89995 = (cse_p565_t1_invr130525 * de_dr);
            double cse_p564_t1_g130524 = (cse_p11_t1_g89995 * (cse_p11_t1_g89996 + (-(cse_p11_t49_g89937))));
            double cse_p564_t13_g130459 = (cse_p11_t1_g89995 * (cse_p11_t1_g89996 + (-(cse_p11_t50_g89938))));
            double cse_p564_t14_g130460 = (cse_p11_t1_g89995 * (cse_p11_t1_g89996 + (-(cse_p11_t51_g89939))));
            double cse_p564_t15_g130461 = (cse_p11_t1_g89995 * (cse_p11_t49_g89937 + (-(cse_p11_t1_g89996))));
            double cse_p564_t16_g130462 = (cse_p11_t1_g89995 * (cse_p11_t50_g89938 + (-(cse_p11_t1_g89996))));
            double cse_p564_t17_g130463 = (cse_p11_t1_g89995 * (cse_p11_t51_g89939 + (-(cse_p11_t1_g89996))));
            double cse_p564_t18_g130464 = (cse_p11_t1_g89995 * dy);
            double cse_p564_t19_g130465 = (cse_p11_t1_g89995 * dz);
            double cse_p69_t45_g95024 = (cse_p564_t9_g130455 * dz);
            double cse_p69_t74_g95053 = (cse_p564_t1_g130524 + cse_p70_t45_g95120);
            double cse_p69_t46_g95025 = (cse_p564_t18_g130464 * (-(dx)));
            double cse_p69_t47_g95026 = (cse_p564_t18_g130464 * dx);
            double cse_p69_t73_g95052 = (cse_p69_t46_g95025 + cse_p70_t46_g95121);
            double h_x1_x1 = cse_p69_t74_g95053;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = cse_p69_t73_g95052;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p69_t31_g95010 = (cse_p564_t19_g130465 * (-(dx)));
            double cse_p69_t32_g95011 = (cse_p564_t19_g130465 * (-(dy)));
            double cse_p69_t33_g95012 = (cse_p564_t19_g130465 * dx);
            double cse_p69_t34_g95013 = (cse_p564_t19_g130465 * dy);
            double cse_p60_t59_g94736 = (cse_p69_t31_g95010 + cse_p69_t45_g95024);
            double cse_p69_t48_g95027 = (cse_p564_t10_g130456 * dy);
            double cse_p69_t49_g95028 = (cse_p564_t10_g130456 * dz);
            double cse_p60_t60_g94737 = (cse_p69_t32_g95011 + cse_p69_t49_g95028);
            double h_x1_z1 = cse_p60_t59_g94736;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p69_t20_g94999 = (cse_p564_t6_g130452 * (-(dx)));
            double cse_p69_t21_g95000 = (cse_p564_t6_g130452 * (-(dy)));
            double cse_p69_t22_g95001 = (cse_p564_t6_g130452 * (-(dz)));
            double cse_p69_t23_g95002 = (cse_p564_t6_g130452 * dx);
            double cse_p69_t24_g95003 = (cse_p564_t6_g130452 * dy);
            double cse_p69_t25_g95004 = (cse_p564_t6_g130452 * dz);
            double cse_p60_t53_g94730 = (cse_p564_t15_g130461 + cse_p69_t23_g95002);
            double cse_p60_t61_g94738 = (cse_p69_t24_g95003 + cse_p69_t47_g95026);
            double cse_p60_t63_g94740 = (cse_p69_t21_g95000 + cse_p69_t46_g95025);
            double cse_p60_t64_g94741 = (cse_p69_t22_g95001 + cse_p69_t31_g95010);
            double cse_p60_t65_g94742 = (cse_p69_t25_g95004 + cse_p69_t33_g95012);
            double cse_p60_t82_g94759 = (cse_p564_t1_g130524 + cse_p69_t20_g94999);
            double h_x1_x2 = cse_p60_t53_g94730;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p69_t26_g95005 = (cse_p564_t7_g130453 * (-(dy)));
            double cse_p69_t27_g95006 = (cse_p564_t7_g130453 * (-(dz)));
            double cse_p69_t28_g95007 = (cse_p564_t7_g130453 * dx);
            double cse_p69_t29_g95008 = (cse_p564_t7_g130453 * dy);
            double cse_p69_t30_g95009 = (cse_p564_t7_g130453 * dz);
            double cse_p60_t54_g94731 = (cse_p564_t16_g130462 + cse_p69_t29_g95008);
            double cse_p60_t62_g94739 = (cse_p69_t28_g95007 + cse_p69_t47_g95026);
            double cse_p60_t66_g94743 = (cse_p69_t27_g95006 + cse_p69_t32_g95011);
            double cse_p60_t67_g94744 = (cse_p69_t30_g95009 + cse_p69_t34_g95013);
            double h_x1_y2 = cse_p60_t62_g94739;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p69_t35_g95014 = (cse_p564_t8_g130454 * (-(dz)));
            double cse_p69_t36_g95015 = (cse_p564_t8_g130454 * dx);
            double cse_p69_t37_g95016 = (cse_p564_t8_g130454 * dy);
            double cse_p69_t38_g95017 = (cse_p564_t8_g130454 * dz);
            double cse_p60_t55_g94732 = (cse_p564_t17_g130463 + cse_p69_t38_g95017);
            double cse_p60_t68_g94745 = (cse_p69_t33_g95012 + cse_p69_t36_g95015);
            double cse_p60_t69_g94746 = (cse_p69_t34_g95013 + cse_p69_t37_g95016);
            double h_x1_z2 = cse_p60_t68_g94745;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p67_t70_g94851 = (cse_p564_t13_g130459 + cse_p69_t48_g95027);
            double cse_p60_t78_g94755 = (cse_p564_t13_g130459 + cse_p69_t26_g95005);
            double h_y1_y1 = cse_p67_t70_g94851;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t60_g94737;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t61_g94738;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t54_g94731;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t69_g94746;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p67_t72_g94853 = (cse_p564_t11_g130457 + cse_p564_t14_g130460);
            double cse_p60_t80_g94757 = (cse_p564_t14_g130460 + cse_p69_t35_g95014);
            double h_z1_z1 = cse_p67_t72_g94853;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t65_g94742;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t67_g94744;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t55_g94732;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t82_g94759;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t63_g94740;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t64_g94741;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t78_g94755;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t66_g94743;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t80_g94757;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p68_t2_g94883 = (-2.00000000000000000e+0 * term.bb * rr);
              double cse_p12_t7_g90010 = (term.bb * rr);
              double cse_p11_t1_g89996 = (rr * rr);
              double cse_p68_t4_g94885 = (-3.00000000000000000e+0 * term.aa * cse_p11_t1_g89996);
              double cse_p67_t56_g94837 = (cse_p68_t2_g94883 + cse_p68_t4_g94885);
              double cse_p12_t6_g90009 = (3.00000000000000000e+0 * term.aa * cse_p11_t1_g89996);
              double cse_p11_t88_g89976 = (term.cc + cse_p12_t6_g90009 + cse_p12_t7_g90010 + cse_p12_t7_g90010);
              double fa_atr = (term.dd + (term.cc * rr) + (cse_p11_t1_g89996 * (term.bb + (term.aa * rr))));
              double de_dr = (term.cc + cse_p12_t6_g90009 + (2.00000000000000000e+0 * term.bb * rr));
              double d2e_dr2 = ((2.00000000000000000e+0 * term.bb) + (6.00000000000000000e+0 * term.aa * rr));
              double energy = fa_atr;
              energy_added += energy;
              *energy_accumulate += energy;
              double g_x1 = (cse_p11_t88_g89976 * cse_p1_t2_invsqrt89886 * dx);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t88_g89976 * cse_p1_t2_invsqrt89886 * dy);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t88_g89976 * cse_p1_t2_invsqrt89886 * dz);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = ((cse_p1_t2_invsqrt89886 * cse_p67_t56_g94837 * dx) + (-((term.cc * cse_p1_t2_invsqrt89886 * dx))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = ((cse_p1_t2_invsqrt89886 * cse_p67_t56_g94837 * dy) + (-((term.cc * cse_p1_t2_invsqrt89886 * dy))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = ((cse_p1_t2_invsqrt89886 * cse_p67_t56_g94837 * dz) + (-((term.cc * cse_p1_t2_invsqrt89886 * dz))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double cse_p565_t2_invr130526 = (1.0 / (rr));
              double cse_p11_t1_g89994 = (cse_p565_t2_invr130526 * cse_p565_t2_invr130526 * d2e_dr2);
              double cse_p565_t1_invr130525 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g89995 = (cse_p565_t1_invr130525 * de_dr);
              double cse_p564_t1_g130524 = (cse_p11_t1_g89995 * (cse_p11_t1_g89996 + (-(cse_p11_t49_g89937))));
              double cse_p67_t1_g94873 = (cse_p11_t1_g89994 * dx);
              double cse_p564_t40_g130486 = (cse_p67_t1_g94873 * dz);
              double h_x1_x1 = (cse_p564_t1_g130524 + (cse_p67_t1_g94873 * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double cse_p67_t1_g94874 = (cse_p11_t1_g89995 * dy);
              double cse_p564_t41_g130487 = (cse_p67_t1_g94874 * (-(dx)));
              double cse_p564_t42_g130488 = (cse_p67_t1_g94874 * dx);
              double h_x1_y1 = (cse_p564_t41_g130487 + (cse_p67_t1_g94873 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g89999 = (cse_p11_t1_g89995 * dz);
              double cse_p564_t29_g130475 = (cse_p11_t1_g89999 * (-(dx)));
              double cse_p564_t30_g130476 = (cse_p11_t1_g89999 * (-(dy)));
              double cse_p564_t31_g130477 = (cse_p11_t1_g89999 * dx);
              double cse_p564_t32_g130478 = (cse_p11_t1_g89999 * dy);
              double cse_p67_t1_g94876 = (cse_p11_t1_g89994 * dy);
              double cse_p564_t43_g130489 = (cse_p67_t1_g94876 * dy);
              double cse_p564_t44_g130490 = (cse_p67_t1_g94876 * dz);
              double cse_p116_t57_g99228 = (cse_p564_t30_g130476 + cse_p564_t44_g130490);
              double h_x1_z1 = (cse_p564_t29_g130475 + cse_p564_t40_g130486);
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g89997 = (cse_p11_t1_g89994 * (-(dx)));
              double cse_p564_t20_g130466 = (cse_p11_t1_g89997 * (-(dx)));
              double cse_p564_t21_g130467 = (cse_p11_t1_g89997 * (-(dy)));
              double cse_p564_t22_g130468 = (cse_p11_t1_g89997 * (-(dz)));
              double cse_p564_t23_g130469 = (cse_p11_t1_g89997 * dy);
              double cse_p564_t24_g130470 = (cse_p11_t1_g89997 * dz);
              double cse_p116_t52_g99223 = (cse_p564_t22_g130468 + cse_p564_t29_g130475);
              double cse_p116_t53_g99224 = (cse_p564_t23_g130469 + cse_p564_t42_g130488);
              double cse_p116_t54_g99225 = (cse_p564_t24_g130470 + cse_p564_t31_g130477);
              double cse_p67_t60_g94841 = (cse_p564_t21_g130467 + cse_p564_t41_g130487);
              double h_x1_x2 = ((cse_p11_t1_g89995 * (cse_p11_t49_g89937 + (-(cse_p11_t1_g89996)))) + (cse_p11_t1_g89997 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g89998 = (cse_p11_t1_g89994 * (-(dy)));
              double cse_p564_t25_g130471 = (cse_p11_t1_g89998 * (-(dy)));
              double cse_p564_t26_g130472 = (cse_p11_t1_g89998 * (-(dz)));
              double cse_p564_t27_g130473 = (cse_p11_t1_g89998 * dy);
              double cse_p564_t28_g130474 = (cse_p11_t1_g89998 * dz);
              double cse_p116_t50_g99221 = (cse_p564_t27_g130473 + (cse_p11_t1_g89995 * (cse_p11_t50_g89938 + (-(cse_p11_t1_g89996)))));
              double cse_p116_t55_g99226 = (cse_p564_t26_g130472 + cse_p564_t30_g130476);
              double cse_p116_t56_g99227 = (cse_p564_t28_g130474 + cse_p564_t32_g130478);
              double h_x1_y2 = (cse_p564_t42_g130488 + (cse_p11_t1_g89998 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g90000 = (cse_p11_t1_g89994 * (-(dz)));
              double cse_p564_t33_g130479 = (cse_p11_t1_g90000 * (-(dz)));
              double cse_p564_t34_g130480 = (cse_p11_t1_g90000 * dy);
              double cse_p564_t35_g130481 = (cse_p11_t1_g90000 * dz);
              double cse_p116_t51_g99222 = (cse_p564_t35_g130481 + (cse_p11_t1_g89995 * (cse_p11_t51_g89939 + (-(cse_p11_t1_g89996)))));
              double cse_p116_t58_g99229 = (cse_p564_t32_g130478 + cse_p564_t34_g130480);
              double h_x1_z2 = (cse_p564_t31_g130477 + (cse_p11_t1_g90000 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g90001 = (cse_p11_t1_g89995 * (cse_p11_t1_g89996 + (-(cse_p11_t50_g89938))));
              double cse_p116_t63_g99234 = (cse_p11_t1_g90001 + cse_p564_t25_g130471);
              double h_y1_y1 = (cse_p11_t1_g90001 + cse_p564_t43_g130489);
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = cse_p116_t57_g99228;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = cse_p116_t53_g99224;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = cse_p116_t50_g99221;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = cse_p116_t58_g99229;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g90002 = (cse_p11_t1_g89995 * (cse_p11_t1_g89996 + (-(cse_p11_t51_g89939))));
              double cse_p116_t64_g99235 = (cse_p11_t1_g90002 + cse_p564_t33_g130479);
              double h_z1_z1 = (cse_p11_t1_g90002 + (cse_p11_t1_g89994 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = cse_p116_t54_g99225;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = cse_p116_t56_g99227;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = cse_p116_t51_g99222;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p564_t1_g130524 + cse_p564_t20_g130466);
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = cse_p67_t60_g94841;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = cse_p116_t52_g99223;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = cse_p116_t63_g99234;
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = cse_p116_t55_g99226;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = cse_p116_t64_g99235;
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
                double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t49_g89937))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t2_invr130526 = (1.0 / (rr));
                double cse_p565_t1_invr130525 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g89996 = (rr * rr);
                double cse_p11_t1_g89995 = (cse_p565_t1_invr130525 * de_dr);
                double cse_p11_t1_g89994 = (cse_p565_t2_invr130526 * cse_p565_t2_invr130526 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g89994 * dx * dy) + (cse_p11_t1_g89995 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double cse_p67_t1_g94876 = (cse_p11_t1_g89994 * dy);
                double cse_p67_t1_g94874 = (cse_p11_t1_g89995 * dy);
                double cse_p67_t1_g94873 = (cse_p11_t1_g89994 * dx);
                double h_x1_z1 = ((cse_p11_t1_g89995 * dz * (-(dx))) + (cse_p67_t1_g94873 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g89999 = (cse_p11_t1_g89995 * dz);
                double h_x1_x2 = ((cse_p11_t1_g89995 * (cse_p11_t49_g89937 + (-(cse_p11_t1_g89996)))) + (cse_p67_t1_g94873 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g89997 = (cse_p11_t1_g89994 * (-(dx)));
                double h_x1_y2 = ((cse_p67_t1_g94873 * (-(dy))) + (cse_p67_t1_g94874 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g89998 = (cse_p11_t1_g89994 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g89999 * dx) + (cse_p67_t1_g94873 * (-(dz))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g90000 = (cse_p11_t1_g89994 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g89995 * (cse_p11_t1_g89996 + (-(cse_p11_t50_g89938)))) + (cse_p67_t1_g94876 * dy));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g90001 = (cse_p11_t1_g89995 * (cse_p11_t1_g89996 + (-(cse_p11_t50_g89938))));
                double h_y1_z1 = ((cse_p11_t1_g89999 * (-(dy))) + (cse_p67_t1_g94876 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g89997 * dy) + (cse_p67_t1_g94874 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g89995 * (cse_p11_t50_g89938 + (-(cse_p11_t1_g89996)))) + (cse_p11_t1_g89998 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g89999 * dy) + (cse_p11_t1_g90000 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g89994 * dz * dz) + (cse_p11_t1_g89995 * (cse_p11_t1_g89996 + (-(cse_p11_t51_g89939)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g90002 = (cse_p11_t1_g89995 * (cse_p11_t1_g89996 + (-(cse_p11_t51_g89939))));
                double h_z1_x2 = ((cse_p11_t1_g89997 * dz) + (cse_p11_t1_g89999 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g89998 * dz) + (cse_p11_t1_g89999 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g89995 * (cse_p11_t51_g89939 + (-(cse_p11_t1_g89996)))) + (cse_p11_t1_g90000 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g89995 * (cse_p11_t1_g89996 + (-(cse_p11_t49_g89937)))) + (cse_p11_t1_g89997 * (-(dx))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p564_t1_g130524 = (cse_p11_t1_g89995 * (cse_p11_t1_g89996 + (-(cse_p11_t49_g89937))));
                double h_x2_y2 = ((cse_p11_t1_g89997 * (-(dy))) + (cse_p67_t1_g94874 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double cse_p67_t1_g94877 = (cse_p67_t1_g94874 * (-(dx)));
                double h_x2_z2 = ((cse_p11_t1_g89997 * (-(dz))) + (cse_p11_t1_g89999 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g90001 + (cse_p11_t1_g89998 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g89998 * (-(dz))) + (cse_p11_t1_g89999 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g90002 + (cse_p11_t1_g90000 * (-(dz))));
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
