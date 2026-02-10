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
      double cse_p1_t1_sqrt170316 = sqrt(r2);
      double cse_p1_t2_invsqrt170317 = (1.0 / (cse_p1_t1_sqrt170316));
      double rr = cse_p1_t1_sqrt170316;
      if ((rr <= params.rmin)) {
        /* !COND-1 */
        double energy = (params.elec_weight * term.e_rmin);
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
        if ((rr <= params.rlow)) {
          /* !COND-2 */
          double cse_p528_t1_g177347 = (cse_p1_t2_invsqrt170317 * dx);
          double cse_p528_t2_g177335 = (cse_p1_t2_invsqrt170317 * dy);
          double cse_p528_t3_g177336 = (cse_p1_t2_invsqrt170317 * dz);
          double cse_p528_t4_g177337 = (rr * rr);
          double cse_p12_t1_g170357 = (term.aa_low * cse_p528_t4_g177337);
          double cse_p12_t4_g170349 = (3.00000000000000000e+0 * cse_p12_t1_g170357);
          double cse_p11_t1_g170341 = (term.bb_low * rr);
          double cse_p11_t12_g170331 = (term.cc_low + cse_p11_t1_g170341 + cse_p11_t1_g170341 + cse_p12_t4_g170349);
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p528_t4_g177337 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p12_t1_g170361 = (cse_p11_t12_g170331 * params.elec_weight);
          double g_x1 = (cse_p12_t1_g170361 * cse_p528_t1_g177347);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p12_t1_g170361 * cse_p528_t2_g177335);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p12_t1_g170361 * cse_p528_t3_g177336);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double cse_p63_t1_g171046 = (-2.00000000000000000e+0 * cse_p11_t1_g170341);
          double cse_p63_t1_g171048 = (-3.00000000000000000e+0 * cse_p12_t1_g170357);
          double g_x2 = (params.elec_weight * ((cse_p528_t1_g177347 * cse_p63_t1_g171046) + (cse_p528_t1_g177347 * cse_p63_t1_g171048) + (-((term.cc_low * cse_p528_t1_g177347)))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (params.elec_weight * ((cse_p528_t2_g177335 * cse_p63_t1_g171046) + (cse_p528_t2_g177335 * cse_p63_t1_g171048) + (-((term.cc_low * cse_p528_t2_g177335)))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (params.elec_weight * ((cse_p528_t3_g177336 * cse_p63_t1_g171046) + (cse_p528_t3_g177336 * cse_p63_t1_g171048) + (-((term.cc_low * cse_p528_t3_g177336)))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (params.elec_weight * (term.cc_low + cse_p12_t4_g170349 + (2.00000000000000000e+0 * cse_p11_t1_g170341)));
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t10_g170355 = (-(params.eps_solvent));
            double cse_p11_t13_g170332 = (cse_p12_t10_g170355 + params.eps_core);
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p12_t1_g170346 = (-1.25000000000000000e-1 * rx);
            double cse_p12_t2_g170347 = (1.25000000000000000e-1 * rx);
            double cse_p12_t3_g170348 = (5.00000000000000000e-1 * rx * rx);
            double cse_p12_t11_g170356 = (-(rx));
            double cse_p11_t8_g170327 = (-2.50000000000000000e-1 + cse_p12_t1_g170346 + cse_p12_t1_g170346);
            double cse_p11_t9_g170328 = (2.50000000000000000e-1 + cse_p12_t2_g170347 + cse_p12_t2_g170347);
            double cse_p11_t10_g170329 = (1.00000000000000000e+0 + cse_p12_t3_g170348 + rx);
            double cse_p11_t16_g170335 = exp(cse_p12_t11_g170356);
            double ex = cse_p11_t16_g170335;
            double gx = (cse_p11_t10_g170329 * ex);
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p12_t11_g170356 + cse_p12_t3_g170348));
            double ep = (params.eps_solvent + (cse_p11_t13_g170332 * gx));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t13_g170332 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t13_g170332 * gx_pp);
            double cse_p1_t3_invr170318 = (1.0 / (ep));
            double cse_p1_t4_invr2170319 = (cse_p1_t3_invr170318 * cse_p1_t3_invr170318);
            double ee = (cse_p1_t3_invr170318 * term.kqq);
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p528_t1_g177347 = (cse_p1_t2_invsqrt170317 * dx);
            double cse_p11_t1_g170342 = (cse_p11_t10_g170329 * cse_p11_t16_g170335);
            double cse_p11_t1_g170343 = (cse_p11_t13_g170332 * cse_p1_t4_invr2170319);
            double cse_p11_t1_g170344 = (params.elec_weight * term.kqq);
            double cse_p12_t1_g170358 = (-2.50000000000000000e-1 * cse_p11_t1_g170342);
            double cse_p12_t1_g170363 = (cse_p11_t9_g170328 * ex);
            double cse_p63_t1_g171044 = (cse_p11_t1_g170343 * cse_p11_t1_g170344);
            double g_x1 = (-((cse_p63_t1_g171044 * ((cse_p12_t1_g170358 * cse_p528_t1_g177347) + (cse_p12_t1_g170363 * cse_p528_t1_g177347)))));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p11_t1_g170337 = (cse_p1_t2_invsqrt170317 * dy);
            double g_y1 = (-((cse_p63_t1_g171044 * ((cse_p11_t1_g170337 * cse_p12_t1_g170358) + (cse_p11_t1_g170337 * cse_p12_t1_g170363)))));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p11_t1_g170338 = (cse_p1_t2_invsqrt170317 * dz);
            double g_z1 = (-((cse_p63_t1_g171044 * ((cse_p11_t1_g170338 * cse_p12_t1_g170358) + (cse_p11_t1_g170338 * cse_p12_t1_g170363)))));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double cse_p12_t1_g170359 = (2.50000000000000000e-1 * cse_p11_t1_g170342);
            double cse_p12_t1_g170362 = (cse_p11_t8_g170327 * ex);
            double g_x2 = (-((cse_p63_t1_g171044 * ((cse_p12_t1_g170359 * cse_p528_t1_g177347) + (cse_p12_t1_g170362 * cse_p528_t1_g177347)))));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (-((cse_p63_t1_g171044 * ((cse_p11_t1_g170337 * cse_p12_t1_g170359) + (cse_p11_t1_g170337 * cse_p12_t1_g170362)))));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (-((cse_p63_t1_g171044 * ((cse_p11_t1_g170338 * cse_p12_t1_g170359) + (cse_p11_t1_g170338 * cse_p12_t1_g170362)))));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (params.elec_weight * ep_p * (1.0 / ((ep * ep))) * (-(term.kqq)));
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p11_t1_g170339 = (rr * rr);
              double cse_p11_t1_g170345 = (term.aa_high * cse_p11_t1_g170339);
              double cse_p12_t5_g170350 = (3.00000000000000000e+0 * cse_p11_t1_g170345);
              double cse_p11_t1_g170340 = (term.bb_high * rr);
              double cse_p11_t11_g170330 = (term.cc_high + cse_p11_t1_g170340 + cse_p11_t1_g170340 + cse_p12_t5_g170350);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p11_t1_g170339 * (term.bb_high + (term.aa_high * rr)))));
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p528_t1_g177347 = (cse_p1_t2_invsqrt170317 * dx);
              double cse_p12_t1_g170360 = (cse_p11_t11_g170330 * params.elec_weight);
              double g_x1 = (cse_p12_t1_g170360 * cse_p528_t1_g177347);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double cse_p11_t1_g170337 = (cse_p1_t2_invsqrt170317 * dy);
              double g_y1 = (cse_p11_t1_g170337 * cse_p12_t1_g170360);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double cse_p11_t1_g170338 = (cse_p1_t2_invsqrt170317 * dz);
              double g_z1 = (cse_p11_t1_g170338 * cse_p12_t1_g170360);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p63_t1_g171045 = (-2.00000000000000000e+0 * cse_p11_t1_g170340);
              double cse_p63_t1_g171047 = (-3.00000000000000000e+0 * cse_p11_t1_g170345);
              double g_x2 = (params.elec_weight * ((cse_p528_t1_g177347 * cse_p63_t1_g171045) + (cse_p528_t1_g177347 * cse_p63_t1_g171047) + (-((term.cc_high * cse_p528_t1_g177347)))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (params.elec_weight * ((cse_p11_t1_g170337 * cse_p63_t1_g171045) + (cse_p11_t1_g170337 * cse_p63_t1_g171047) + (-((term.cc_high * cse_p11_t1_g170337)))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (params.elec_weight * ((cse_p11_t1_g170338 * cse_p63_t1_g171045) + (cse_p11_t1_g170338 * cse_p63_t1_g171047) + (-((term.cc_high * cse_p11_t1_g170338)))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double de_dr = (params.elec_weight * (term.cc_high + cse_p12_t5_g170350 + (2.00000000000000000e+0 * cse_p11_t1_g170340)));
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
          double cse_p528_t1_g182861 = (rr * rr);
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p528_t1_g182861 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t2_g177364 = (-(params.eps_solvent));
            double cse_p11_t2_g177350 = (cse_p12_t2_g177364 + params.eps_core);
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p11_t1_g177349 = (5.00000000000000000e-1 * rx * rx);
            double cse_p11_t4_g177352 = (-(rx));
            double ex = exp(cse_p11_t4_g177352);
            double gx = (ex * (1.00000000000000000e+0 + cse_p11_t1_g177349 + rx));
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p11_t1_g177349 + cse_p11_t4_g177352));
            double ep = (params.eps_solvent + (cse_p11_t2_g177350 * gx));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t2_g177350 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t2_g177350 * gx_pp);
            double ee = (term.kqq * (1.0 / (ep)));
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p528_t1_g182861 = (rr * rr);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p528_t1_g182861 * (term.bb_high + (term.aa_high * rr)))));
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
      double cse_p11_t41_g182907 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t42_g182908 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t43_g182909 = (dz * dz);
      double r2 = (cse_p11_t41_g182907 + cse_p11_t42_g182908 + cse_p11_t43_g182909);
      double cse_p1_t1_sqrt182863 = sqrt(r2);
      double cse_p1_t2_invsqrt182864 = (1.0 / (cse_p1_t1_sqrt182863));
      double rr = cse_p1_t1_sqrt182863;
      if ((rr <= params.rmin)) {
        /* !COND-1 */
        double cse_p60_t40_g186777 = (rr * rr);
        double cse_p60_t68_g186805 = (-(cse_p11_t41_g182907));
        double cse_p60_t69_g186806 = (-(cse_p11_t42_g182908));
        double cse_p60_t70_g186807 = (-(cse_p11_t43_g182909));
        double cse_p60_t72_g186809 = (-(dx));
        double cse_p60_t73_g186810 = (-(dy));
        double cse_p60_t74_g186811 = (-(dz));
        double cse_p565_t1_invr222248 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr222249 = (1.0 / (rr));
        double cse_p60_t77_g186814 = (rr * rr * rr);
        double cse_p15_t72_g183343 = (-(cse_p60_t40_g186777));
        double cse_p14_t76_g183247 = (cse_p11_t41_g182907 + cse_p15_t72_g183343);
        double cse_p14_t78_g183249 = (cse_p11_t42_g182908 + cse_p15_t72_g183343);
        double cse_p14_t80_g183251 = (cse_p11_t43_g182909 + cse_p15_t72_g183343);
        double cse_p13_t101_g183158 = (cse_p60_t40_g186777 + cse_p60_t68_g186805);
        double cse_p13_t102_g183159 = (cse_p60_t40_g186777 + cse_p60_t69_g186806);
        double cse_p13_t103_g183160 = (cse_p60_t40_g186777 + cse_p60_t70_g186807);
        double energy = (params.elec_weight * term.e_rmin);
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
        double cse_p60_t39_g186776 = (cse_p565_t1_invr222248 * de_dr);
        double cse_p13_t52_g183109 = (cse_p14_t76_g183247 * cse_p60_t39_g186776);
        double cse_p13_t53_g183110 = (cse_p14_t78_g183249 * cse_p60_t39_g186776);
        double cse_p13_t54_g183111 = (cse_p14_t80_g183251 * cse_p60_t39_g186776);
        double cse_p13_t55_g183112 = (cse_p60_t39_g186776 * dx * dy);
        double cse_p13_t56_g183113 = (cse_p60_t39_g186776 * cse_p60_t72_g186809 * dy);
        double cse_p12_t20_g182990 = (cse_p13_t101_g183158 * cse_p60_t39_g186776);
        double cse_p12_t21_g182991 = (cse_p13_t102_g183159 * cse_p60_t39_g186776);
        double cse_p12_t22_g182992 = (cse_p13_t103_g183160 * cse_p60_t39_g186776);
        double cse_p12_t25_g182995 = (cse_p60_t39_g186776 * dz);
        double cse_p13_t42_g183099 = (cse_p12_t25_g182995 * cse_p60_t72_g186809);
        double cse_p13_t43_g183100 = (cse_p12_t25_g182995 * cse_p60_t73_g186810);
        double cse_p13_t44_g183101 = (cse_p12_t25_g182995 * dx);
        double cse_p13_t45_g183102 = (cse_p12_t25_g182995 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p60_t38_g186775 = (cse_p565_t2_invr222249 * cse_p565_t2_invr222249 * d2e_dr2);
        double cse_p13_t46_g183103 = (cse_p60_t38_g186775 * dx * dx);
        double cse_p13_t47_g183104 = (cse_p60_t38_g186775 * dx * dy);
        double cse_p13_t48_g183105 = (cse_p60_t38_g186775 * dx * dz);
        double cse_p13_t49_g183106 = (cse_p60_t38_g186775 * dy * dy);
        double cse_p13_t50_g183107 = (cse_p60_t38_g186775 * dy * dz);
        double cse_p13_t51_g183108 = (cse_p60_t38_g186775 * dz * dz);
        double cse_p12_t8_g182978 = (cse_p60_t38_g186775 * cse_p60_t72_g186809);
        double cse_p13_t57_g183114 = (cse_p12_t8_g182978 * cse_p60_t72_g186809);
        double cse_p13_t58_g183115 = (cse_p12_t8_g182978 * cse_p60_t73_g186810);
        double cse_p13_t59_g183116 = (cse_p12_t8_g182978 * cse_p60_t74_g186811);
        double cse_p13_t60_g183117 = (cse_p12_t8_g182978 * dx);
        double cse_p13_t61_g183118 = (cse_p12_t8_g182978 * dy);
        double cse_p13_t62_g183119 = (cse_p12_t8_g182978 * dz);
        double cse_p12_t9_g182979 = (cse_p60_t38_g186775 * cse_p60_t73_g186810);
        double cse_p13_t63_g183120 = (cse_p12_t9_g182979 * cse_p60_t73_g186810);
        double cse_p13_t64_g183121 = (cse_p12_t9_g182979 * cse_p60_t74_g186811);
        double cse_p13_t65_g183122 = (cse_p12_t9_g182979 * dx);
        double cse_p13_t66_g183123 = (cse_p12_t9_g182979 * dy);
        double cse_p13_t67_g183124 = (cse_p12_t9_g182979 * dz);
        double cse_p12_t10_g182980 = (cse_p60_t38_g186775 * cse_p60_t74_g186811);
        double cse_p13_t38_g183095 = (cse_p12_t10_g182980 * cse_p60_t74_g186811);
        double cse_p13_t39_g183096 = (cse_p12_t10_g182980 * dx);
        double cse_p13_t40_g183097 = (cse_p12_t10_g182980 * dy);
        double cse_p13_t41_g183098 = (cse_p12_t10_g182980 * dz);
        double cse_p12_t52_g183022 = (cse_p13_t47_g183104 + cse_p13_t56_g183113);
        double cse_p12_t53_g183023 = (cse_p13_t42_g183099 + cse_p13_t48_g183105);
        double cse_p12_t54_g183024 = (cse_p13_t43_g183100 + cse_p13_t50_g183107);
        double cse_p12_t70_g183040 = (cse_p12_t21_g182991 + cse_p13_t49_g183106);
        double cse_p12_t72_g183042 = (cse_p12_t22_g182992 + cse_p13_t51_g183108);
        double cse_p12_t74_g183044 = (cse_p12_t20_g182990 + cse_p13_t46_g183103);
        double cse_p12_t49_g183019 = (cse_p13_t52_g183109 + cse_p13_t60_g183117);
        double cse_p12_t55_g183025 = (cse_p13_t55_g183112 + cse_p13_t61_g183118);
        double cse_p12_t57_g183027 = (cse_p13_t56_g183113 + cse_p13_t58_g183115);
        double cse_p12_t58_g183028 = (cse_p13_t42_g183099 + cse_p13_t59_g183116);
        double cse_p12_t59_g183029 = (cse_p13_t44_g183101 + cse_p13_t62_g183119);
        double cse_p12_t75_g183045 = (cse_p12_t20_g182990 + cse_p13_t57_g183114);
        double cse_p12_t50_g183020 = (cse_p13_t53_g183110 + cse_p13_t66_g183123);
        double cse_p12_t56_g183026 = (cse_p13_t55_g183112 + cse_p13_t65_g183122);
        double cse_p12_t60_g183030 = (cse_p13_t43_g183100 + cse_p13_t64_g183121);
        double cse_p12_t61_g183031 = (cse_p13_t45_g183102 + cse_p13_t67_g183124);
        double cse_p12_t71_g183041 = (cse_p12_t21_g182991 + cse_p13_t63_g183120);
        double cse_p12_t51_g183021 = (cse_p13_t41_g183098 + cse_p13_t54_g183111);
        double cse_p12_t62_g183032 = (cse_p13_t39_g183096 + cse_p13_t44_g183101);
        double cse_p12_t63_g183033 = (cse_p13_t40_g183097 + cse_p13_t45_g183102);
        double cse_p12_t73_g183043 = (cse_p12_t22_g182992 + cse_p13_t38_g183095);
        double h_x1_x1 = cse_p12_t74_g183044;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t52_g183022;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t53_g183023;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t49_g183019;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t56_g183026;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t62_g183032;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t70_g183040;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t54_g183024;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t55_g183025;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t50_g183020;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t63_g183033;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t72_g183042;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t59_g183029;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t61_g183031;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t51_g183021;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t75_g183045;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t57_g183027;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t58_g183028;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t71_g183041;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t60_g183030;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t73_g183043;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= params.rlow)) {
          /* !COND-2 */
          double cse_p12_t7_g182977 = (term.bb_low * rr);
          double cse_p11_t1_g182962 = (rr * rr);
          double cse_p60_t56_g186793 = (cse_p11_t41_g182907 + (-(cse_p11_t1_g182962)));
          double cse_p60_t57_g186794 = (cse_p11_t42_g182908 + (-(cse_p11_t1_g182962)));
          double cse_p60_t58_g186795 = (cse_p11_t43_g182909 + (-(cse_p11_t1_g182962)));
          double cse_p60_t59_g186796 = (cse_p11_t1_g182962 + (-(cse_p11_t41_g182907)));
          double cse_p60_t60_g186797 = (cse_p11_t1_g182962 + (-(cse_p11_t42_g182908)));
          double cse_p60_t61_g186798 = (cse_p11_t1_g182962 + (-(cse_p11_t43_g182909)));
          double cse_p60_t71_g186808 = (-(cse_p11_t1_g182962));
          double cse_p12_t5_g182975 = (3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g182962);
          double cse_p11_t77_g182943 = (term.cc_low + cse_p12_t5_g182975 + cse_p12_t7_g182977 + cse_p12_t7_g182977);
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p11_t1_g182962 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (cse_p11_t77_g182943 * cse_p1_t2_invsqrt182864 * dx * params.elec_weight);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t77_g182943 * cse_p1_t2_invsqrt182864 * dy * params.elec_weight);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t77_g182943 * cse_p1_t2_invsqrt182864 * dz * params.elec_weight);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (params.elec_weight * ((cse_p1_t2_invsqrt182864 * dx * ((-2.00000000000000000e+0 * term.bb_low * rr) + (-3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g182962))) + (-((term.cc_low * cse_p1_t2_invsqrt182864 * dx)))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (params.elec_weight * ((cse_p1_t2_invsqrt182864 * dy * ((-2.00000000000000000e+0 * term.bb_low * rr) + (-3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g182962))) + (-((term.cc_low * cse_p1_t2_invsqrt182864 * dy)))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (params.elec_weight * ((cse_p1_t2_invsqrt182864 * dz * ((-2.00000000000000000e+0 * term.bb_low * rr) + (-3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g182962))) + (-((term.cc_low * cse_p1_t2_invsqrt182864 * dz)))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (params.elec_weight * (term.cc_low + cse_p12_t5_g182975 + (2.00000000000000000e+0 * term.bb_low * rr)));
          double d2e_dr2 = (params.elec_weight * ((2.00000000000000000e+0 * term.bb_low) + (6.00000000000000000e+0 * term.aa_low * rr)));
          double cse_p565_t2_invr222249 = (1.0 / (rr));
          double cse_p11_t1_g182960 = (cse_p565_t2_invr222249 * cse_p565_t2_invr222249 * d2e_dr2);
          double cse_p60_t1_g186738 = (cse_p11_t1_g182960 * (-(dx)));
          double cse_p60_t2_g186739 = (cse_p11_t1_g182960 * (-(dy)));
          double cse_p60_t3_g186740 = (cse_p11_t1_g182960 * (-(dz)));
          double cse_p60_t4_g186741 = (cse_p11_t1_g182960 * dx * dx);
          double cse_p60_t5_g186742 = (cse_p11_t1_g182960 * dx * dy);
          double cse_p60_t6_g186743 = (cse_p11_t1_g182960 * dx * dz);
          double cse_p60_t7_g186744 = (cse_p11_t1_g182960 * dy * dy);
          double cse_p60_t8_g186745 = (cse_p11_t1_g182960 * dy * dz);
          double cse_p60_t9_g186746 = (cse_p11_t1_g182960 * dz * dz);
          double cse_p14_t42_g183213 = (cse_p60_t1_g186738 * (-(dx)));
          double cse_p14_t43_g183214 = (cse_p60_t1_g186738 * (-(dy)));
          double cse_p14_t44_g183215 = (cse_p60_t1_g186738 * (-(dz)));
          double cse_p14_t45_g183216 = (cse_p60_t1_g186738 * dx);
          double cse_p14_t46_g183217 = (cse_p60_t1_g186738 * dy);
          double cse_p14_t47_g183218 = (cse_p60_t1_g186738 * dz);
          double cse_p14_t48_g183219 = (cse_p60_t2_g186739 * (-(dy)));
          double cse_p14_t49_g183220 = (cse_p60_t2_g186739 * (-(dz)));
          double cse_p14_t50_g183221 = (cse_p60_t2_g186739 * dx);
          double cse_p14_t51_g183222 = (cse_p60_t2_g186739 * dy);
          double cse_p14_t52_g183223 = (cse_p60_t2_g186739 * dz);
          double cse_p14_t53_g183224 = (cse_p60_t3_g186740 * (-(dz)));
          double cse_p14_t54_g183225 = (cse_p60_t3_g186740 * dx);
          double cse_p14_t55_g183226 = (cse_p60_t3_g186740 * dy);
          double cse_p14_t56_g183227 = (cse_p60_t3_g186740 * dz);
          double cse_p565_t1_invr222248 = (1.0 / ((rr * rr * rr)));
          double cse_p11_t1_g182961 = (cse_p565_t1_invr222248 * de_dr);
          double cse_p60_t10_g186747 = (cse_p11_t1_g182961 * (cse_p11_t41_g182907 + cse_p60_t71_g186808));
          double cse_p60_t11_g186748 = (cse_p11_t1_g182961 * (cse_p11_t42_g182908 + cse_p60_t71_g186808));
          double cse_p60_t12_g186749 = (cse_p11_t1_g182961 * (cse_p11_t43_g182909 + cse_p60_t71_g186808));
          double cse_p60_t13_g186750 = (cse_p11_t1_g182961 * cse_p60_t59_g186796);
          double cse_p60_t14_g186751 = (cse_p11_t1_g182961 * cse_p60_t60_g186797);
          double cse_p60_t15_g186752 = (cse_p11_t1_g182961 * cse_p60_t61_g186798);
          double cse_p60_t16_g186753 = (cse_p11_t1_g182961 * dx * dy);
          double cse_p60_t17_g186754 = (cse_p11_t1_g182961 * dy * (-(dx)));
          double cse_p60_t18_g186755 = (cse_p11_t1_g182961 * dz);
          double cse_p60_t44_g186781 = (cse_p60_t17_g186754 + cse_p60_t5_g186742);
          double cse_p14_t38_g183209 = (cse_p60_t18_g186755 * (-(dx)));
          double cse_p14_t39_g183210 = (cse_p60_t18_g186755 * (-(dy)));
          double cse_p14_t40_g183211 = (cse_p60_t18_g186755 * dx);
          double cse_p14_t41_g183212 = (cse_p60_t18_g186755 * dy);
          double cse_p13_t99_g183156 = (cse_p60_t13_g186750 + cse_p60_t4_g186741);
          double h_x1_x1 = cse_p13_t99_g183156;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t44_g186781;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t75_g183132 = (cse_p14_t38_g183209 + cse_p60_t6_g186743);
          double cse_p13_t76_g183133 = (cse_p14_t39_g183210 + cse_p60_t8_g186745);
          double h_x1_z1 = cse_p13_t75_g183132;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t71_g183128 = (cse_p14_t45_g183216 + cse_p60_t10_g186747);
          double cse_p13_t77_g183134 = (cse_p14_t46_g183217 + cse_p60_t16_g186753);
          double cse_p13_t79_g183136 = (cse_p14_t43_g183214 + cse_p60_t17_g186754);
          double cse_p13_t80_g183137 = (cse_p14_t38_g183209 + cse_p14_t44_g183215);
          double cse_p13_t81_g183138 = (cse_p14_t40_g183211 + cse_p14_t47_g183218);
          double cse_p13_t100_g183157 = (cse_p14_t42_g183213 + cse_p60_t13_g186750);
          double h_x1_x2 = cse_p13_t71_g183128;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t72_g183129 = (cse_p14_t51_g183222 + cse_p60_t11_g186748);
          double cse_p13_t78_g183135 = (cse_p14_t50_g183221 + cse_p60_t16_g186753);
          double cse_p13_t82_g183139 = (cse_p14_t39_g183210 + cse_p14_t49_g183220);
          double cse_p13_t83_g183140 = (cse_p14_t41_g183212 + cse_p14_t52_g183223);
          double h_x1_y2 = cse_p13_t78_g183135;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t73_g183130 = (cse_p14_t56_g183227 + cse_p60_t12_g186749);
          double cse_p13_t84_g183141 = (cse_p14_t40_g183211 + cse_p14_t54_g183225);
          double cse_p13_t85_g183142 = (cse_p14_t41_g183212 + cse_p14_t55_g183226);
          double h_x1_z2 = cse_p13_t84_g183141;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t95_g183152 = (cse_p60_t14_g186751 + cse_p60_t7_g186744);
          double cse_p13_t96_g183153 = (cse_p14_t48_g183219 + cse_p60_t14_g186751);
          double h_y1_y1 = cse_p13_t95_g183152;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t76_g183133;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t77_g183134;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t72_g183129;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t85_g183142;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t97_g183154 = (cse_p60_t15_g186752 + cse_p60_t9_g186746);
          double cse_p13_t98_g183155 = (cse_p14_t53_g183224 + cse_p60_t15_g186752);
          double h_z1_z1 = cse_p13_t97_g183154;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t81_g183138;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t83_g183140;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t73_g183130;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t100_g183157;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t79_g183136;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t80_g183137;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t96_g183153;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t82_g183139;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t98_g183155;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t83_g183053 = (-(params.eps_solvent));
            double cse_p11_t78_g182944 = (cse_p12_t83_g183053 + params.eps_core);
            double cse_p11_t87_g182953 = (-(term.kqq));
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p12_t1_g182971 = (-1.25000000000000000e-1 * rx);
            double cse_p12_t2_g182972 = (1.25000000000000000e-1 * rx);
            double cse_p12_t3_g182973 = (5.00000000000000000e-1 * rx * rx);
            double cse_p12_t84_g183054 = (-(rx));
            double cse_p11_t73_g182939 = (-2.50000000000000000e-1 + cse_p12_t1_g182971 + cse_p12_t1_g182971);
            double cse_p11_t74_g182940 = (2.50000000000000000e-1 + cse_p12_t2_g182972 + cse_p12_t2_g182972);
            double cse_p11_t75_g182941 = (1.00000000000000000e+0 + cse_p12_t3_g182973 + rx);
            double cse_p11_t89_g182955 = exp(cse_p12_t84_g183054);
            double ex = cse_p11_t89_g182955;
            double gx = (cse_p11_t75_g182941 * ex);
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p12_t3_g182973 + cse_p12_t84_g183054));
            double ep = (params.eps_solvent + (cse_p11_t78_g182944 * gx));
            double cse_p12_t47_g183017 = (ep * ep);
            double cse_p11_t90_g182956 = (1.0 / (cse_p12_t47_g183017));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t78_g182944 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t78_g182944 * gx_pp);
            double cse_p1_t3_invr182865 = (1.0 / (ep));
            double cse_p1_t4_invr2182866 = (cse_p1_t3_invr182865 * cse_p1_t3_invr182865);
            double ee = (cse_p1_t3_invr182865 * term.kqq);
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p11_t1_g182970 = (cse_p11_t75_g182941 * cse_p11_t89_g182955 * cse_p1_t2_invsqrt182864);
            double g_x1 = (-((cse_p11_t78_g182944 * cse_p1_t4_invr2182866 * params.elec_weight * term.kqq * ((-2.50000000000000000e-1 * cse_p11_t1_g182970 * dx) + (cse_p11_t74_g182940 * cse_p1_t2_invsqrt182864 * dx * ex)))));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double g_y1 = (-((cse_p11_t78_g182944 * cse_p1_t4_invr2182866 * params.elec_weight * term.kqq * ((-2.50000000000000000e-1 * cse_p11_t1_g182970 * dy) + (cse_p11_t74_g182940 * cse_p1_t2_invsqrt182864 * dy * ex)))));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double g_z1 = (-((cse_p11_t78_g182944 * cse_p1_t4_invr2182866 * params.elec_weight * term.kqq * ((-2.50000000000000000e-1 * cse_p11_t1_g182970 * dz) + (cse_p11_t74_g182940 * cse_p1_t2_invsqrt182864 * dz * ex)))));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (-((cse_p11_t78_g182944 * cse_p1_t4_invr2182866 * params.elec_weight * term.kqq * ((2.50000000000000000e-1 * cse_p11_t1_g182970 * dx) + (cse_p11_t73_g182939 * cse_p1_t2_invsqrt182864 * dx * ex)))));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (-((cse_p11_t78_g182944 * cse_p1_t4_invr2182866 * params.elec_weight * term.kqq * ((2.50000000000000000e-1 * cse_p11_t1_g182970 * dy) + (cse_p11_t73_g182939 * cse_p1_t2_invsqrt182864 * dy * ex)))));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (-((cse_p11_t78_g182944 * cse_p1_t4_invr2182866 * params.elec_weight * term.kqq * ((2.50000000000000000e-1 * cse_p11_t1_g182970 * dz) + (cse_p11_t73_g182939 * cse_p1_t2_invsqrt182864 * dz * ex)))));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (cse_p11_t87_g182953 * cse_p11_t90_g182956 * params.elec_weight * ep_p);
            double d2e_dr2 = (cse_p11_t87_g182953 * params.elec_weight * ((-2.00000000000000000e+0 * ep_p * ep_p * (1.0 / ((ep * ep * ep)))) + (cse_p11_t90_g182956 * ep_pp)));
            double cse_p565_t2_invr222249 = (1.0 / (rr));
            double cse_p11_t1_g182960 = (cse_p565_t2_invr222249 * cse_p565_t2_invr222249 * d2e_dr2);
            double cse_p565_t1_invr222248 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g182961 = (cse_p565_t1_invr222248 * de_dr);
            double cse_p11_t1_g182962 = (rr * rr);
            double cse_p11_t1_g182969 = (cse_p11_t1_g182961 * (cse_p11_t1_g182962 + (-(cse_p11_t41_g182907))));
            double cse_p60_t66_g186803 = (cse_p11_t1_g182969 + (cse_p11_t1_g182960 * dx * dx));
            double h_x1_x1 = cse_p60_t66_g186803;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = ((cse_p11_t1_g182960 * dx * dy) + (cse_p11_t1_g182961 * dy * (-(dx))));
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p11_t1_g182965 = (cse_p11_t1_g182961 * dz);
            double cse_p60_t30_g186767 = (cse_p11_t1_g182965 * (-(dx)));
            double cse_p60_t31_g186768 = (cse_p11_t1_g182965 * (-(dy)));
            double cse_p60_t32_g186769 = (cse_p11_t1_g182965 * dx);
            double cse_p60_t33_g186770 = (cse_p11_t1_g182965 * dy);
            double cse_p60_t45_g186782 = (cse_p60_t30_g186767 + (cse_p11_t1_g182960 * dx * dz));
            double cse_p60_t46_g186783 = (cse_p60_t31_g186768 + (cse_p11_t1_g182960 * dy * dz));
            double h_x1_z1 = cse_p60_t45_g186782;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p11_t1_g182963 = (cse_p11_t1_g182960 * (-(dx)));
            double cse_p60_t19_g186756 = (cse_p11_t1_g182963 * (-(dx)));
            double cse_p60_t20_g186757 = (cse_p11_t1_g182963 * (-(dy)));
            double cse_p60_t21_g186758 = (cse_p11_t1_g182963 * (-(dz)));
            double cse_p60_t22_g186759 = (cse_p11_t1_g182963 * dx);
            double cse_p60_t23_g186760 = (cse_p11_t1_g182963 * dy);
            double cse_p60_t24_g186761 = (cse_p11_t1_g182963 * dz);
            double cse_p60_t41_g186778 = (cse_p60_t22_g186759 + (cse_p11_t1_g182961 * (cse_p11_t41_g182907 + (-(cse_p11_t1_g182962)))));
            double cse_p60_t47_g186784 = (cse_p60_t23_g186760 + (cse_p11_t1_g182961 * dx * dy));
            double cse_p60_t49_g186786 = (cse_p60_t20_g186757 + (cse_p11_t1_g182961 * dy * (-(dx))));
            double cse_p60_t50_g186787 = (cse_p60_t21_g186758 + cse_p60_t30_g186767);
            double cse_p60_t51_g186788 = (cse_p60_t24_g186761 + cse_p60_t32_g186769);
            double cse_p60_t67_g186804 = (cse_p11_t1_g182969 + cse_p60_t19_g186756);
            double h_x1_x2 = cse_p60_t41_g186778;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p11_t1_g182964 = (cse_p11_t1_g182960 * (-(dy)));
            double cse_p60_t25_g186762 = (cse_p11_t1_g182964 * (-(dy)));
            double cse_p60_t26_g186763 = (cse_p11_t1_g182964 * (-(dz)));
            double cse_p60_t27_g186764 = (cse_p11_t1_g182964 * dx);
            double cse_p60_t28_g186765 = (cse_p11_t1_g182964 * dy);
            double cse_p60_t29_g186766 = (cse_p11_t1_g182964 * dz);
            double cse_p60_t42_g186779 = (cse_p60_t28_g186765 + (cse_p11_t1_g182961 * (cse_p11_t42_g182908 + (-(cse_p11_t1_g182962)))));
            double cse_p60_t48_g186785 = (cse_p60_t27_g186764 + (cse_p11_t1_g182961 * dx * dy));
            double cse_p60_t52_g186789 = (cse_p60_t26_g186763 + cse_p60_t31_g186768);
            double cse_p60_t53_g186790 = (cse_p60_t29_g186766 + cse_p60_t33_g186770);
            double h_x1_y2 = cse_p60_t48_g186785;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p11_t1_g182966 = (cse_p11_t1_g182960 * (-(dz)));
            double cse_p60_t34_g186771 = (cse_p11_t1_g182966 * (-(dz)));
            double cse_p60_t35_g186772 = (cse_p11_t1_g182966 * dx);
            double cse_p60_t36_g186773 = (cse_p11_t1_g182966 * dy);
            double cse_p60_t37_g186774 = (cse_p11_t1_g182966 * dz);
            double cse_p60_t43_g186780 = (cse_p60_t37_g186774 + (cse_p11_t1_g182961 * (cse_p11_t43_g182909 + (-(cse_p11_t1_g182962)))));
            double cse_p60_t54_g186791 = (cse_p60_t32_g186769 + cse_p60_t35_g186772);
            double cse_p60_t55_g186792 = (cse_p60_t33_g186770 + cse_p60_t36_g186773);
            double h_x1_z2 = cse_p60_t54_g186791;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p11_t1_g182967 = (cse_p11_t1_g182961 * (cse_p11_t1_g182962 + (-(cse_p11_t42_g182908))));
            double cse_p60_t62_g186799 = (cse_p11_t1_g182967 + (cse_p11_t1_g182960 * dy * dy));
            double cse_p60_t63_g186800 = (cse_p11_t1_g182967 + cse_p60_t25_g186762);
            double h_y1_y1 = cse_p60_t62_g186799;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t46_g186783;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t47_g186784;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t42_g186779;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t55_g186792;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p11_t1_g182968 = (cse_p11_t1_g182961 * (cse_p11_t1_g182962 + (-(cse_p11_t43_g182909))));
            double cse_p60_t64_g186801 = (cse_p11_t1_g182968 + (cse_p11_t1_g182960 * dz * dz));
            double cse_p60_t65_g186802 = (cse_p11_t1_g182968 + cse_p60_t34_g186771);
            double h_z1_z1 = cse_p60_t64_g186801;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t51_g186788;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t53_g186790;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t43_g186780;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t67_g186804;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t49_g186786;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t50_g186787;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t63_g186800;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t52_g186789;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t65_g186802;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p12_t6_g182976 = (term.bb_high * rr);
              double cse_p11_t1_g182962 = (rr * rr);
              double cse_p12_t4_g182974 = (3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g182962);
              double cse_p11_t76_g182942 = (term.cc_high + cse_p12_t4_g182974 + cse_p12_t6_g182976 + cse_p12_t6_g182976);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p11_t1_g182962 * (term.bb_high + (term.aa_high * rr)))));
              energy_added += energy;
              *energy_accumulate += energy;
              double g_x1 = (cse_p11_t76_g182942 * cse_p1_t2_invsqrt182864 * dx * params.elec_weight);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t76_g182942 * cse_p1_t2_invsqrt182864 * dy * params.elec_weight);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t76_g182942 * cse_p1_t2_invsqrt182864 * dz * params.elec_weight);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = (params.elec_weight * ((cse_p1_t2_invsqrt182864 * dx * ((-2.00000000000000000e+0 * term.bb_high * rr) + (-3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g182962))) + (-((term.cc_high * cse_p1_t2_invsqrt182864 * dx)))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (params.elec_weight * ((cse_p1_t2_invsqrt182864 * dy * ((-2.00000000000000000e+0 * term.bb_high * rr) + (-3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g182962))) + (-((term.cc_high * cse_p1_t2_invsqrt182864 * dy)))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (params.elec_weight * ((cse_p1_t2_invsqrt182864 * dz * ((-2.00000000000000000e+0 * term.bb_high * rr) + (-3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g182962))) + (-((term.cc_high * cse_p1_t2_invsqrt182864 * dz)))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double de_dr = (params.elec_weight * (term.cc_high + cse_p12_t4_g182974 + (2.00000000000000000e+0 * term.bb_high * rr)));
              double d2e_dr2 = (params.elec_weight * ((2.00000000000000000e+0 * term.bb_high) + (6.00000000000000000e+0 * term.aa_high * rr)));
              double cse_p565_t2_invr222249 = (1.0 / (rr));
              double cse_p11_t1_g182960 = (cse_p565_t2_invr222249 * cse_p565_t2_invr222249 * d2e_dr2);
              double cse_p565_t1_invr222248 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g182961 = (cse_p565_t1_invr222248 * de_dr);
              double cse_p11_t1_g182969 = (cse_p11_t1_g182961 * (cse_p11_t1_g182962 + (-(cse_p11_t41_g182907))));
              double h_x1_x1 = (cse_p11_t1_g182969 + (cse_p11_t1_g182960 * dx * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double h_x1_y1 = ((cse_p11_t1_g182960 * dx * dy) + (cse_p11_t1_g182961 * dy * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g182965 = (cse_p11_t1_g182961 * dz);
              double h_x1_z1 = ((cse_p11_t1_g182960 * dx * dz) + (cse_p11_t1_g182965 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g182963 = (cse_p11_t1_g182960 * (-(dx)));
              double h_x1_x2 = ((cse_p11_t1_g182961 * (cse_p11_t41_g182907 + (-(cse_p11_t1_g182962)))) + (cse_p11_t1_g182963 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g182964 = (cse_p11_t1_g182960 * (-(dy)));
              double h_x1_y2 = ((cse_p11_t1_g182961 * dx * dy) + (cse_p11_t1_g182964 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g182966 = (cse_p11_t1_g182960 * (-(dz)));
              double h_x1_z2 = ((cse_p11_t1_g182965 * dx) + (cse_p11_t1_g182966 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g182967 = (cse_p11_t1_g182961 * (cse_p11_t1_g182962 + (-(cse_p11_t42_g182908))));
              double h_y1_y1 = (cse_p11_t1_g182967 + (cse_p11_t1_g182960 * dy * dy));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = ((cse_p11_t1_g182960 * dy * dz) + (cse_p11_t1_g182965 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = ((cse_p11_t1_g182961 * dx * dy) + (cse_p11_t1_g182963 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = ((cse_p11_t1_g182961 * (cse_p11_t42_g182908 + (-(cse_p11_t1_g182962)))) + (cse_p11_t1_g182964 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = ((cse_p11_t1_g182965 * dy) + (cse_p11_t1_g182966 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g182968 = (cse_p11_t1_g182961 * (cse_p11_t1_g182962 + (-(cse_p11_t43_g182909))));
              double h_z1_z1 = (cse_p11_t1_g182968 + (cse_p11_t1_g182960 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = ((cse_p11_t1_g182963 * dz) + (cse_p11_t1_g182965 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = ((cse_p11_t1_g182964 * dz) + (cse_p11_t1_g182965 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = ((cse_p11_t1_g182961 * (cse_p11_t43_g182909 + (-(cse_p11_t1_g182962)))) + (cse_p11_t1_g182966 * dz));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p11_t1_g182969 + (cse_p11_t1_g182963 * (-(dx))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = ((cse_p11_t1_g182961 * dy * (-(dx))) + (cse_p11_t1_g182963 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = ((cse_p11_t1_g182963 * (-(dz))) + (cse_p11_t1_g182965 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = (cse_p11_t1_g182967 + (cse_p11_t1_g182964 * (-(dy))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = ((cse_p11_t1_g182964 * (-(dz))) + (cse_p11_t1_g182965 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = (cse_p11_t1_g182968 + (cse_p11_t1_g182966 * (-(dz))));
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
                double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t41_g182907))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t2_invr222249 = (1.0 / (rr));
                double cse_p565_t1_invr222248 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g182962 = (rr * rr);
                double cse_p11_t1_g182961 = (cse_p565_t1_invr222248 * de_dr);
                double cse_p11_t1_g182960 = (cse_p565_t2_invr222249 * cse_p565_t2_invr222249 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g182960 * dx * dy) + (cse_p11_t1_g182961 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double h_x1_z1 = ((cse_p11_t1_g182960 * dx * dz) + (cse_p11_t1_g182961 * dz * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g182965 = (cse_p11_t1_g182961 * dz);
                double h_x1_x2 = ((cse_p11_t1_g182960 * dx * (-(dx))) + (cse_p11_t1_g182961 * (cse_p11_t41_g182907 + (-(cse_p11_t1_g182962)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g182963 = (cse_p11_t1_g182960 * (-(dx)));
                double h_x1_y2 = ((cse_p11_t1_g182960 * dx * (-(dy))) + (cse_p11_t1_g182961 * dx * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g182964 = (cse_p11_t1_g182960 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g182960 * dx * (-(dz))) + (cse_p11_t1_g182965 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g182966 = (cse_p11_t1_g182960 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g182960 * dy * dy) + (cse_p11_t1_g182961 * (cse_p11_t1_g182962 + (-(cse_p11_t42_g182908)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g182967 = (cse_p11_t1_g182961 * (cse_p11_t1_g182962 + (-(cse_p11_t42_g182908))));
                double h_y1_z1 = ((cse_p11_t1_g182960 * dy * dz) + (cse_p11_t1_g182965 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g182961 * dx * dy) + (cse_p11_t1_g182963 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g182961 * (cse_p11_t42_g182908 + (-(cse_p11_t1_g182962)))) + (cse_p11_t1_g182964 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g182965 * dy) + (cse_p11_t1_g182966 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g182960 * dz * dz) + (cse_p11_t1_g182961 * (cse_p11_t1_g182962 + (-(cse_p11_t43_g182909)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g182968 = (cse_p11_t1_g182961 * (cse_p11_t1_g182962 + (-(cse_p11_t43_g182909))));
                double h_z1_x2 = ((cse_p11_t1_g182963 * dz) + (cse_p11_t1_g182965 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g182964 * dz) + (cse_p11_t1_g182965 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g182961 * (cse_p11_t43_g182909 + (-(cse_p11_t1_g182962)))) + (cse_p11_t1_g182966 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g182961 * (cse_p11_t1_g182962 + (-(cse_p11_t41_g182907)))) + (cse_p11_t1_g182963 * (-(dx))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p11_t1_g182969 = (cse_p11_t1_g182961 * (cse_p11_t1_g182962 + (-(cse_p11_t41_g182907))));
                double h_x2_y2 = ((cse_p11_t1_g182961 * dy * (-(dx))) + (cse_p11_t1_g182963 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double h_x2_z2 = ((cse_p11_t1_g182963 * (-(dz))) + (cse_p11_t1_g182965 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g182967 + (cse_p11_t1_g182964 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g182964 * (-(dz))) + (cse_p11_t1_g182965 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g182968 + (cse_p11_t1_g182966 * (-(dz))));
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
