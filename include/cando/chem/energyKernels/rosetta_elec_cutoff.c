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
      double cse_p1_t1_sqrt15302 = sqrt(r2);
      double cse_p1_t2_invsqrt15303 = (1.0 / (cse_p1_t1_sqrt15302));
      double rr = cse_p1_t1_sqrt15302;
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
          double cse_p528_t1_g22333 = (cse_p1_t2_invsqrt15303 * dx);
          double cse_p528_t2_g22321 = (cse_p1_t2_invsqrt15303 * dy);
          double cse_p528_t3_g22322 = (cse_p1_t2_invsqrt15303 * dz);
          double cse_p528_t4_g22323 = (rr * rr);
          double cse_p12_t1_g15343 = (term.aa_low * cse_p528_t4_g22323);
          double cse_p12_t4_g15335 = (3.00000000000000000e+0 * cse_p12_t1_g15343);
          double cse_p11_t1_g15327 = (term.bb_low * rr);
          double cse_p11_t12_g15317 = (term.cc_low + cse_p11_t1_g15327 + cse_p11_t1_g15327 + cse_p12_t4_g15335);
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p528_t4_g22323 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p12_t1_g15347 = (cse_p11_t12_g15317 * params.elec_weight);
          double g_x1 = (cse_p12_t1_g15347 * cse_p528_t1_g22333);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p12_t1_g15347 * cse_p528_t2_g22321);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p12_t1_g15347 * cse_p528_t3_g22322);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double cse_p63_t1_g16032 = (-2.00000000000000000e+0 * cse_p11_t1_g15327);
          double cse_p63_t1_g16034 = (-3.00000000000000000e+0 * cse_p12_t1_g15343);
          double g_x2 = (params.elec_weight * ((cse_p528_t1_g22333 * cse_p63_t1_g16032) + (cse_p528_t1_g22333 * cse_p63_t1_g16034) + (-((term.cc_low * cse_p528_t1_g22333)))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (params.elec_weight * ((cse_p528_t2_g22321 * cse_p63_t1_g16032) + (cse_p528_t2_g22321 * cse_p63_t1_g16034) + (-((term.cc_low * cse_p528_t2_g22321)))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (params.elec_weight * ((cse_p528_t3_g22322 * cse_p63_t1_g16032) + (cse_p528_t3_g22322 * cse_p63_t1_g16034) + (-((term.cc_low * cse_p528_t3_g22322)))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (params.elec_weight * (term.cc_low + cse_p12_t4_g15335 + (2.00000000000000000e+0 * cse_p11_t1_g15327)));
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t10_g15341 = (-(params.eps_solvent));
            double cse_p11_t13_g15318 = (cse_p12_t10_g15341 + params.eps_core);
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p12_t1_g15332 = (-1.25000000000000000e-1 * rx);
            double cse_p12_t2_g15333 = (1.25000000000000000e-1 * rx);
            double cse_p12_t3_g15334 = (5.00000000000000000e-1 * rx * rx);
            double cse_p12_t11_g15342 = (-(rx));
            double cse_p11_t8_g15313 = (-2.50000000000000000e-1 + cse_p12_t1_g15332 + cse_p12_t1_g15332);
            double cse_p11_t9_g15314 = (2.50000000000000000e-1 + cse_p12_t2_g15333 + cse_p12_t2_g15333);
            double cse_p11_t10_g15315 = (1.00000000000000000e+0 + cse_p12_t3_g15334 + rx);
            double cse_p11_t16_g15321 = exp(cse_p12_t11_g15342);
            double ex = cse_p11_t16_g15321;
            double gx = (cse_p11_t10_g15315 * ex);
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p12_t11_g15342 + cse_p12_t3_g15334));
            double ep = (params.eps_solvent + (cse_p11_t13_g15318 * gx));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t13_g15318 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t13_g15318 * gx_pp);
            double cse_p1_t3_invr15304 = (1.0 / (ep));
            double cse_p1_t4_invr215305 = (cse_p1_t3_invr15304 * cse_p1_t3_invr15304);
            double ee = (cse_p1_t3_invr15304 * term.kqq);
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p528_t1_g22333 = (cse_p1_t2_invsqrt15303 * dx);
            double cse_p11_t1_g15328 = (cse_p11_t10_g15315 * cse_p11_t16_g15321);
            double cse_p11_t1_g15329 = (cse_p11_t13_g15318 * cse_p1_t4_invr215305);
            double cse_p11_t1_g15330 = (params.elec_weight * term.kqq);
            double cse_p12_t1_g15344 = (-2.50000000000000000e-1 * cse_p11_t1_g15328);
            double cse_p12_t1_g15349 = (cse_p11_t9_g15314 * ex);
            double cse_p63_t1_g16030 = (cse_p11_t1_g15329 * cse_p11_t1_g15330);
            double g_x1 = (-((cse_p63_t1_g16030 * ((cse_p12_t1_g15344 * cse_p528_t1_g22333) + (cse_p12_t1_g15349 * cse_p528_t1_g22333)))));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p11_t1_g15323 = (cse_p1_t2_invsqrt15303 * dy);
            double g_y1 = (-((cse_p63_t1_g16030 * ((cse_p11_t1_g15323 * cse_p12_t1_g15344) + (cse_p11_t1_g15323 * cse_p12_t1_g15349)))));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p11_t1_g15324 = (cse_p1_t2_invsqrt15303 * dz);
            double g_z1 = (-((cse_p63_t1_g16030 * ((cse_p11_t1_g15324 * cse_p12_t1_g15344) + (cse_p11_t1_g15324 * cse_p12_t1_g15349)))));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double cse_p12_t1_g15345 = (2.50000000000000000e-1 * cse_p11_t1_g15328);
            double cse_p12_t1_g15348 = (cse_p11_t8_g15313 * ex);
            double g_x2 = (-((cse_p63_t1_g16030 * ((cse_p12_t1_g15345 * cse_p528_t1_g22333) + (cse_p12_t1_g15348 * cse_p528_t1_g22333)))));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (-((cse_p63_t1_g16030 * ((cse_p11_t1_g15323 * cse_p12_t1_g15345) + (cse_p11_t1_g15323 * cse_p12_t1_g15348)))));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (-((cse_p63_t1_g16030 * ((cse_p11_t1_g15324 * cse_p12_t1_g15345) + (cse_p11_t1_g15324 * cse_p12_t1_g15348)))));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (params.elec_weight * ep_p * (1.0 / ((ep * ep))) * (-(term.kqq)));
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p11_t1_g15325 = (rr * rr);
              double cse_p11_t1_g15331 = (term.aa_high * cse_p11_t1_g15325);
              double cse_p12_t5_g15336 = (3.00000000000000000e+0 * cse_p11_t1_g15331);
              double cse_p11_t1_g15326 = (term.bb_high * rr);
              double cse_p11_t11_g15316 = (term.cc_high + cse_p11_t1_g15326 + cse_p11_t1_g15326 + cse_p12_t5_g15336);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p11_t1_g15325 * (term.bb_high + (term.aa_high * rr)))));
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p528_t1_g22333 = (cse_p1_t2_invsqrt15303 * dx);
              double cse_p12_t1_g15346 = (cse_p11_t11_g15316 * params.elec_weight);
              double g_x1 = (cse_p12_t1_g15346 * cse_p528_t1_g22333);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double cse_p11_t1_g15323 = (cse_p1_t2_invsqrt15303 * dy);
              double g_y1 = (cse_p11_t1_g15323 * cse_p12_t1_g15346);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double cse_p11_t1_g15324 = (cse_p1_t2_invsqrt15303 * dz);
              double g_z1 = (cse_p11_t1_g15324 * cse_p12_t1_g15346);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p63_t1_g16031 = (-2.00000000000000000e+0 * cse_p11_t1_g15326);
              double cse_p63_t1_g16033 = (-3.00000000000000000e+0 * cse_p11_t1_g15331);
              double g_x2 = (params.elec_weight * ((cse_p528_t1_g22333 * cse_p63_t1_g16031) + (cse_p528_t1_g22333 * cse_p63_t1_g16033) + (-((term.cc_high * cse_p528_t1_g22333)))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (params.elec_weight * ((cse_p11_t1_g15323 * cse_p63_t1_g16031) + (cse_p11_t1_g15323 * cse_p63_t1_g16033) + (-((term.cc_high * cse_p11_t1_g15323)))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (params.elec_weight * ((cse_p11_t1_g15324 * cse_p63_t1_g16031) + (cse_p11_t1_g15324 * cse_p63_t1_g16033) + (-((term.cc_high * cse_p11_t1_g15324)))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double de_dr = (params.elec_weight * (term.cc_high + cse_p12_t5_g15336 + (2.00000000000000000e+0 * cse_p11_t1_g15326)));
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
          double cse_p528_t1_g27847 = (rr * rr);
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p528_t1_g27847 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t2_g22350 = (-(params.eps_solvent));
            double cse_p11_t2_g22336 = (cse_p12_t2_g22350 + params.eps_core);
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p11_t1_g22335 = (5.00000000000000000e-1 * rx * rx);
            double cse_p11_t4_g22338 = (-(rx));
            double ex = exp(cse_p11_t4_g22338);
            double gx = (ex * (1.00000000000000000e+0 + cse_p11_t1_g22335 + rx));
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p11_t1_g22335 + cse_p11_t4_g22338));
            double ep = (params.eps_solvent + (cse_p11_t2_g22336 * gx));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t2_g22336 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t2_g22336 * gx_pp);
            double ee = (term.kqq * (1.0 / (ep)));
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p528_t1_g27847 = (rr * rr);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p528_t1_g27847 * (term.bb_high + (term.aa_high * rr)))));
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
      double cse_p11_t41_g27893 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t42_g27894 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t43_g27895 = (dz * dz);
      double r2 = (cse_p11_t41_g27893 + cse_p11_t42_g27894 + cse_p11_t43_g27895);
      double cse_p1_t1_sqrt27849 = sqrt(r2);
      double cse_p1_t2_invsqrt27850 = (1.0 / (cse_p1_t1_sqrt27849));
      double rr = cse_p1_t1_sqrt27849;
      if ((rr <= params.rmin)) {
        /* !COND-1 */
        double cse_p60_t40_g31763 = (rr * rr);
        double cse_p60_t68_g31791 = (-(cse_p11_t41_g27893));
        double cse_p60_t69_g31792 = (-(cse_p11_t42_g27894));
        double cse_p60_t70_g31793 = (-(cse_p11_t43_g27895));
        double cse_p60_t72_g31795 = (-(dx));
        double cse_p60_t73_g31796 = (-(dy));
        double cse_p60_t74_g31797 = (-(dz));
        double cse_p565_t1_invr67234 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr67235 = (1.0 / (rr));
        double cse_p60_t77_g31800 = (rr * rr * rr);
        double cse_p15_t72_g28329 = (-(cse_p60_t40_g31763));
        double cse_p14_t76_g28233 = (cse_p11_t41_g27893 + cse_p15_t72_g28329);
        double cse_p14_t78_g28235 = (cse_p11_t42_g27894 + cse_p15_t72_g28329);
        double cse_p14_t80_g28237 = (cse_p11_t43_g27895 + cse_p15_t72_g28329);
        double cse_p13_t101_g28144 = (cse_p60_t40_g31763 + cse_p60_t68_g31791);
        double cse_p13_t102_g28145 = (cse_p60_t40_g31763 + cse_p60_t69_g31792);
        double cse_p13_t103_g28146 = (cse_p60_t40_g31763 + cse_p60_t70_g31793);
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
        double cse_p60_t39_g31762 = (cse_p565_t1_invr67234 * de_dr);
        double cse_p13_t52_g28095 = (cse_p14_t76_g28233 * cse_p60_t39_g31762);
        double cse_p13_t53_g28096 = (cse_p14_t78_g28235 * cse_p60_t39_g31762);
        double cse_p13_t54_g28097 = (cse_p14_t80_g28237 * cse_p60_t39_g31762);
        double cse_p13_t55_g28098 = (cse_p60_t39_g31762 * dx * dy);
        double cse_p13_t56_g28099 = (cse_p60_t39_g31762 * cse_p60_t72_g31795 * dy);
        double cse_p12_t20_g27976 = (cse_p13_t101_g28144 * cse_p60_t39_g31762);
        double cse_p12_t21_g27977 = (cse_p13_t102_g28145 * cse_p60_t39_g31762);
        double cse_p12_t22_g27978 = (cse_p13_t103_g28146 * cse_p60_t39_g31762);
        double cse_p12_t25_g27981 = (cse_p60_t39_g31762 * dz);
        double cse_p13_t42_g28085 = (cse_p12_t25_g27981 * cse_p60_t72_g31795);
        double cse_p13_t43_g28086 = (cse_p12_t25_g27981 * cse_p60_t73_g31796);
        double cse_p13_t44_g28087 = (cse_p12_t25_g27981 * dx);
        double cse_p13_t45_g28088 = (cse_p12_t25_g27981 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p60_t38_g31761 = (cse_p565_t2_invr67235 * cse_p565_t2_invr67235 * d2e_dr2);
        double cse_p13_t46_g28089 = (cse_p60_t38_g31761 * dx * dx);
        double cse_p13_t47_g28090 = (cse_p60_t38_g31761 * dx * dy);
        double cse_p13_t48_g28091 = (cse_p60_t38_g31761 * dx * dz);
        double cse_p13_t49_g28092 = (cse_p60_t38_g31761 * dy * dy);
        double cse_p13_t50_g28093 = (cse_p60_t38_g31761 * dy * dz);
        double cse_p13_t51_g28094 = (cse_p60_t38_g31761 * dz * dz);
        double cse_p12_t8_g27964 = (cse_p60_t38_g31761 * cse_p60_t72_g31795);
        double cse_p13_t57_g28100 = (cse_p12_t8_g27964 * cse_p60_t72_g31795);
        double cse_p13_t58_g28101 = (cse_p12_t8_g27964 * cse_p60_t73_g31796);
        double cse_p13_t59_g28102 = (cse_p12_t8_g27964 * cse_p60_t74_g31797);
        double cse_p13_t60_g28103 = (cse_p12_t8_g27964 * dx);
        double cse_p13_t61_g28104 = (cse_p12_t8_g27964 * dy);
        double cse_p13_t62_g28105 = (cse_p12_t8_g27964 * dz);
        double cse_p12_t9_g27965 = (cse_p60_t38_g31761 * cse_p60_t73_g31796);
        double cse_p13_t63_g28106 = (cse_p12_t9_g27965 * cse_p60_t73_g31796);
        double cse_p13_t64_g28107 = (cse_p12_t9_g27965 * cse_p60_t74_g31797);
        double cse_p13_t65_g28108 = (cse_p12_t9_g27965 * dx);
        double cse_p13_t66_g28109 = (cse_p12_t9_g27965 * dy);
        double cse_p13_t67_g28110 = (cse_p12_t9_g27965 * dz);
        double cse_p12_t10_g27966 = (cse_p60_t38_g31761 * cse_p60_t74_g31797);
        double cse_p13_t38_g28081 = (cse_p12_t10_g27966 * cse_p60_t74_g31797);
        double cse_p13_t39_g28082 = (cse_p12_t10_g27966 * dx);
        double cse_p13_t40_g28083 = (cse_p12_t10_g27966 * dy);
        double cse_p13_t41_g28084 = (cse_p12_t10_g27966 * dz);
        double cse_p12_t52_g28008 = (cse_p13_t47_g28090 + cse_p13_t56_g28099);
        double cse_p12_t53_g28009 = (cse_p13_t42_g28085 + cse_p13_t48_g28091);
        double cse_p12_t54_g28010 = (cse_p13_t43_g28086 + cse_p13_t50_g28093);
        double cse_p12_t70_g28026 = (cse_p12_t21_g27977 + cse_p13_t49_g28092);
        double cse_p12_t72_g28028 = (cse_p12_t22_g27978 + cse_p13_t51_g28094);
        double cse_p12_t74_g28030 = (cse_p12_t20_g27976 + cse_p13_t46_g28089);
        double cse_p12_t49_g28005 = (cse_p13_t52_g28095 + cse_p13_t60_g28103);
        double cse_p12_t55_g28011 = (cse_p13_t55_g28098 + cse_p13_t61_g28104);
        double cse_p12_t57_g28013 = (cse_p13_t56_g28099 + cse_p13_t58_g28101);
        double cse_p12_t58_g28014 = (cse_p13_t42_g28085 + cse_p13_t59_g28102);
        double cse_p12_t59_g28015 = (cse_p13_t44_g28087 + cse_p13_t62_g28105);
        double cse_p12_t75_g28031 = (cse_p12_t20_g27976 + cse_p13_t57_g28100);
        double cse_p12_t50_g28006 = (cse_p13_t53_g28096 + cse_p13_t66_g28109);
        double cse_p12_t56_g28012 = (cse_p13_t55_g28098 + cse_p13_t65_g28108);
        double cse_p12_t60_g28016 = (cse_p13_t43_g28086 + cse_p13_t64_g28107);
        double cse_p12_t61_g28017 = (cse_p13_t45_g28088 + cse_p13_t67_g28110);
        double cse_p12_t71_g28027 = (cse_p12_t21_g27977 + cse_p13_t63_g28106);
        double cse_p12_t51_g28007 = (cse_p13_t41_g28084 + cse_p13_t54_g28097);
        double cse_p12_t62_g28018 = (cse_p13_t39_g28082 + cse_p13_t44_g28087);
        double cse_p12_t63_g28019 = (cse_p13_t40_g28083 + cse_p13_t45_g28088);
        double cse_p12_t73_g28029 = (cse_p12_t22_g27978 + cse_p13_t38_g28081);
        double h_x1_x1 = cse_p12_t74_g28030;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t52_g28008;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t53_g28009;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t49_g28005;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t56_g28012;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t62_g28018;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t70_g28026;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t54_g28010;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t55_g28011;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t50_g28006;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t63_g28019;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t72_g28028;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t59_g28015;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t61_g28017;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t51_g28007;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t75_g28031;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t57_g28013;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t58_g28014;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t71_g28027;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t60_g28016;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t73_g28029;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= params.rlow)) {
          /* !COND-2 */
          double cse_p12_t7_g27963 = (term.bb_low * rr);
          double cse_p11_t1_g27948 = (rr * rr);
          double cse_p60_t56_g31779 = (cse_p11_t41_g27893 + (-(cse_p11_t1_g27948)));
          double cse_p60_t57_g31780 = (cse_p11_t42_g27894 + (-(cse_p11_t1_g27948)));
          double cse_p60_t58_g31781 = (cse_p11_t43_g27895 + (-(cse_p11_t1_g27948)));
          double cse_p60_t59_g31782 = (cse_p11_t1_g27948 + (-(cse_p11_t41_g27893)));
          double cse_p60_t60_g31783 = (cse_p11_t1_g27948 + (-(cse_p11_t42_g27894)));
          double cse_p60_t61_g31784 = (cse_p11_t1_g27948 + (-(cse_p11_t43_g27895)));
          double cse_p60_t71_g31794 = (-(cse_p11_t1_g27948));
          double cse_p12_t5_g27961 = (3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g27948);
          double cse_p11_t77_g27929 = (term.cc_low + cse_p12_t5_g27961 + cse_p12_t7_g27963 + cse_p12_t7_g27963);
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p11_t1_g27948 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (cse_p11_t77_g27929 * cse_p1_t2_invsqrt27850 * dx * params.elec_weight);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t77_g27929 * cse_p1_t2_invsqrt27850 * dy * params.elec_weight);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t77_g27929 * cse_p1_t2_invsqrt27850 * dz * params.elec_weight);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (params.elec_weight * ((cse_p1_t2_invsqrt27850 * dx * ((-2.00000000000000000e+0 * term.bb_low * rr) + (-3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g27948))) + (-((term.cc_low * cse_p1_t2_invsqrt27850 * dx)))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (params.elec_weight * ((cse_p1_t2_invsqrt27850 * dy * ((-2.00000000000000000e+0 * term.bb_low * rr) + (-3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g27948))) + (-((term.cc_low * cse_p1_t2_invsqrt27850 * dy)))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (params.elec_weight * ((cse_p1_t2_invsqrt27850 * dz * ((-2.00000000000000000e+0 * term.bb_low * rr) + (-3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g27948))) + (-((term.cc_low * cse_p1_t2_invsqrt27850 * dz)))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (params.elec_weight * (term.cc_low + cse_p12_t5_g27961 + (2.00000000000000000e+0 * term.bb_low * rr)));
          double d2e_dr2 = (params.elec_weight * ((2.00000000000000000e+0 * term.bb_low) + (6.00000000000000000e+0 * term.aa_low * rr)));
          double cse_p565_t2_invr67235 = (1.0 / (rr));
          double cse_p11_t1_g27946 = (cse_p565_t2_invr67235 * cse_p565_t2_invr67235 * d2e_dr2);
          double cse_p60_t1_g31724 = (cse_p11_t1_g27946 * (-(dx)));
          double cse_p60_t2_g31725 = (cse_p11_t1_g27946 * (-(dy)));
          double cse_p60_t3_g31726 = (cse_p11_t1_g27946 * (-(dz)));
          double cse_p60_t4_g31727 = (cse_p11_t1_g27946 * dx * dx);
          double cse_p60_t5_g31728 = (cse_p11_t1_g27946 * dx * dy);
          double cse_p60_t6_g31729 = (cse_p11_t1_g27946 * dx * dz);
          double cse_p60_t7_g31730 = (cse_p11_t1_g27946 * dy * dy);
          double cse_p60_t8_g31731 = (cse_p11_t1_g27946 * dy * dz);
          double cse_p60_t9_g31732 = (cse_p11_t1_g27946 * dz * dz);
          double cse_p14_t42_g28199 = (cse_p60_t1_g31724 * (-(dx)));
          double cse_p14_t43_g28200 = (cse_p60_t1_g31724 * (-(dy)));
          double cse_p14_t44_g28201 = (cse_p60_t1_g31724 * (-(dz)));
          double cse_p14_t45_g28202 = (cse_p60_t1_g31724 * dx);
          double cse_p14_t46_g28203 = (cse_p60_t1_g31724 * dy);
          double cse_p14_t47_g28204 = (cse_p60_t1_g31724 * dz);
          double cse_p14_t48_g28205 = (cse_p60_t2_g31725 * (-(dy)));
          double cse_p14_t49_g28206 = (cse_p60_t2_g31725 * (-(dz)));
          double cse_p14_t50_g28207 = (cse_p60_t2_g31725 * dx);
          double cse_p14_t51_g28208 = (cse_p60_t2_g31725 * dy);
          double cse_p14_t52_g28209 = (cse_p60_t2_g31725 * dz);
          double cse_p14_t53_g28210 = (cse_p60_t3_g31726 * (-(dz)));
          double cse_p14_t54_g28211 = (cse_p60_t3_g31726 * dx);
          double cse_p14_t55_g28212 = (cse_p60_t3_g31726 * dy);
          double cse_p14_t56_g28213 = (cse_p60_t3_g31726 * dz);
          double cse_p565_t1_invr67234 = (1.0 / ((rr * rr * rr)));
          double cse_p11_t1_g27947 = (cse_p565_t1_invr67234 * de_dr);
          double cse_p60_t10_g31733 = (cse_p11_t1_g27947 * (cse_p11_t41_g27893 + cse_p60_t71_g31794));
          double cse_p60_t11_g31734 = (cse_p11_t1_g27947 * (cse_p11_t42_g27894 + cse_p60_t71_g31794));
          double cse_p60_t12_g31735 = (cse_p11_t1_g27947 * (cse_p11_t43_g27895 + cse_p60_t71_g31794));
          double cse_p60_t13_g31736 = (cse_p11_t1_g27947 * cse_p60_t59_g31782);
          double cse_p60_t14_g31737 = (cse_p11_t1_g27947 * cse_p60_t60_g31783);
          double cse_p60_t15_g31738 = (cse_p11_t1_g27947 * cse_p60_t61_g31784);
          double cse_p60_t16_g31739 = (cse_p11_t1_g27947 * dx * dy);
          double cse_p60_t17_g31740 = (cse_p11_t1_g27947 * dy * (-(dx)));
          double cse_p60_t18_g31741 = (cse_p11_t1_g27947 * dz);
          double cse_p60_t44_g31767 = (cse_p60_t17_g31740 + cse_p60_t5_g31728);
          double cse_p14_t38_g28195 = (cse_p60_t18_g31741 * (-(dx)));
          double cse_p14_t39_g28196 = (cse_p60_t18_g31741 * (-(dy)));
          double cse_p14_t40_g28197 = (cse_p60_t18_g31741 * dx);
          double cse_p14_t41_g28198 = (cse_p60_t18_g31741 * dy);
          double cse_p13_t99_g28142 = (cse_p60_t13_g31736 + cse_p60_t4_g31727);
          double h_x1_x1 = cse_p13_t99_g28142;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t44_g31767;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t75_g28118 = (cse_p14_t38_g28195 + cse_p60_t6_g31729);
          double cse_p13_t76_g28119 = (cse_p14_t39_g28196 + cse_p60_t8_g31731);
          double h_x1_z1 = cse_p13_t75_g28118;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t71_g28114 = (cse_p14_t45_g28202 + cse_p60_t10_g31733);
          double cse_p13_t77_g28120 = (cse_p14_t46_g28203 + cse_p60_t16_g31739);
          double cse_p13_t79_g28122 = (cse_p14_t43_g28200 + cse_p60_t17_g31740);
          double cse_p13_t80_g28123 = (cse_p14_t38_g28195 + cse_p14_t44_g28201);
          double cse_p13_t81_g28124 = (cse_p14_t40_g28197 + cse_p14_t47_g28204);
          double cse_p13_t100_g28143 = (cse_p14_t42_g28199 + cse_p60_t13_g31736);
          double h_x1_x2 = cse_p13_t71_g28114;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t72_g28115 = (cse_p14_t51_g28208 + cse_p60_t11_g31734);
          double cse_p13_t78_g28121 = (cse_p14_t50_g28207 + cse_p60_t16_g31739);
          double cse_p13_t82_g28125 = (cse_p14_t39_g28196 + cse_p14_t49_g28206);
          double cse_p13_t83_g28126 = (cse_p14_t41_g28198 + cse_p14_t52_g28209);
          double h_x1_y2 = cse_p13_t78_g28121;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t73_g28116 = (cse_p14_t56_g28213 + cse_p60_t12_g31735);
          double cse_p13_t84_g28127 = (cse_p14_t40_g28197 + cse_p14_t54_g28211);
          double cse_p13_t85_g28128 = (cse_p14_t41_g28198 + cse_p14_t55_g28212);
          double h_x1_z2 = cse_p13_t84_g28127;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t95_g28138 = (cse_p60_t14_g31737 + cse_p60_t7_g31730);
          double cse_p13_t96_g28139 = (cse_p14_t48_g28205 + cse_p60_t14_g31737);
          double h_y1_y1 = cse_p13_t95_g28138;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t76_g28119;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t77_g28120;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t72_g28115;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t85_g28128;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t97_g28140 = (cse_p60_t15_g31738 + cse_p60_t9_g31732);
          double cse_p13_t98_g28141 = (cse_p14_t53_g28210 + cse_p60_t15_g31738);
          double h_z1_z1 = cse_p13_t97_g28140;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t81_g28124;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t83_g28126;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t73_g28116;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t100_g28143;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t79_g28122;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t80_g28123;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t96_g28139;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t82_g28125;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t98_g28141;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t83_g28039 = (-(params.eps_solvent));
            double cse_p11_t78_g27930 = (cse_p12_t83_g28039 + params.eps_core);
            double cse_p11_t87_g27939 = (-(term.kqq));
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p12_t1_g27957 = (-1.25000000000000000e-1 * rx);
            double cse_p12_t2_g27958 = (1.25000000000000000e-1 * rx);
            double cse_p12_t3_g27959 = (5.00000000000000000e-1 * rx * rx);
            double cse_p12_t84_g28040 = (-(rx));
            double cse_p11_t73_g27925 = (-2.50000000000000000e-1 + cse_p12_t1_g27957 + cse_p12_t1_g27957);
            double cse_p11_t74_g27926 = (2.50000000000000000e-1 + cse_p12_t2_g27958 + cse_p12_t2_g27958);
            double cse_p11_t75_g27927 = (1.00000000000000000e+0 + cse_p12_t3_g27959 + rx);
            double cse_p11_t89_g27941 = exp(cse_p12_t84_g28040);
            double ex = cse_p11_t89_g27941;
            double gx = (cse_p11_t75_g27927 * ex);
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p12_t3_g27959 + cse_p12_t84_g28040));
            double ep = (params.eps_solvent + (cse_p11_t78_g27930 * gx));
            double cse_p12_t47_g28003 = (ep * ep);
            double cse_p11_t90_g27942 = (1.0 / (cse_p12_t47_g28003));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t78_g27930 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t78_g27930 * gx_pp);
            double cse_p1_t3_invr27851 = (1.0 / (ep));
            double cse_p1_t4_invr227852 = (cse_p1_t3_invr27851 * cse_p1_t3_invr27851);
            double ee = (cse_p1_t3_invr27851 * term.kqq);
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p11_t1_g27956 = (cse_p11_t75_g27927 * cse_p11_t89_g27941 * cse_p1_t2_invsqrt27850);
            double g_x1 = (-((cse_p11_t78_g27930 * cse_p1_t4_invr227852 * params.elec_weight * term.kqq * ((-2.50000000000000000e-1 * cse_p11_t1_g27956 * dx) + (cse_p11_t74_g27926 * cse_p1_t2_invsqrt27850 * dx * ex)))));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double g_y1 = (-((cse_p11_t78_g27930 * cse_p1_t4_invr227852 * params.elec_weight * term.kqq * ((-2.50000000000000000e-1 * cse_p11_t1_g27956 * dy) + (cse_p11_t74_g27926 * cse_p1_t2_invsqrt27850 * dy * ex)))));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double g_z1 = (-((cse_p11_t78_g27930 * cse_p1_t4_invr227852 * params.elec_weight * term.kqq * ((-2.50000000000000000e-1 * cse_p11_t1_g27956 * dz) + (cse_p11_t74_g27926 * cse_p1_t2_invsqrt27850 * dz * ex)))));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (-((cse_p11_t78_g27930 * cse_p1_t4_invr227852 * params.elec_weight * term.kqq * ((2.50000000000000000e-1 * cse_p11_t1_g27956 * dx) + (cse_p11_t73_g27925 * cse_p1_t2_invsqrt27850 * dx * ex)))));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (-((cse_p11_t78_g27930 * cse_p1_t4_invr227852 * params.elec_weight * term.kqq * ((2.50000000000000000e-1 * cse_p11_t1_g27956 * dy) + (cse_p11_t73_g27925 * cse_p1_t2_invsqrt27850 * dy * ex)))));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (-((cse_p11_t78_g27930 * cse_p1_t4_invr227852 * params.elec_weight * term.kqq * ((2.50000000000000000e-1 * cse_p11_t1_g27956 * dz) + (cse_p11_t73_g27925 * cse_p1_t2_invsqrt27850 * dz * ex)))));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (cse_p11_t87_g27939 * cse_p11_t90_g27942 * params.elec_weight * ep_p);
            double d2e_dr2 = (cse_p11_t87_g27939 * params.elec_weight * ((-2.00000000000000000e+0 * ep_p * ep_p * (1.0 / ((ep * ep * ep)))) + (cse_p11_t90_g27942 * ep_pp)));
            double cse_p565_t2_invr67235 = (1.0 / (rr));
            double cse_p11_t1_g27946 = (cse_p565_t2_invr67235 * cse_p565_t2_invr67235 * d2e_dr2);
            double cse_p565_t1_invr67234 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g27947 = (cse_p565_t1_invr67234 * de_dr);
            double cse_p11_t1_g27948 = (rr * rr);
            double cse_p11_t1_g27955 = (cse_p11_t1_g27947 * (cse_p11_t1_g27948 + (-(cse_p11_t41_g27893))));
            double cse_p60_t66_g31789 = (cse_p11_t1_g27955 + (cse_p11_t1_g27946 * dx * dx));
            double h_x1_x1 = cse_p60_t66_g31789;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = ((cse_p11_t1_g27946 * dx * dy) + (cse_p11_t1_g27947 * dy * (-(dx))));
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p11_t1_g27951 = (cse_p11_t1_g27947 * dz);
            double cse_p60_t30_g31753 = (cse_p11_t1_g27951 * (-(dx)));
            double cse_p60_t31_g31754 = (cse_p11_t1_g27951 * (-(dy)));
            double cse_p60_t32_g31755 = (cse_p11_t1_g27951 * dx);
            double cse_p60_t33_g31756 = (cse_p11_t1_g27951 * dy);
            double cse_p60_t45_g31768 = (cse_p60_t30_g31753 + (cse_p11_t1_g27946 * dx * dz));
            double cse_p60_t46_g31769 = (cse_p60_t31_g31754 + (cse_p11_t1_g27946 * dy * dz));
            double h_x1_z1 = cse_p60_t45_g31768;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p11_t1_g27949 = (cse_p11_t1_g27946 * (-(dx)));
            double cse_p60_t19_g31742 = (cse_p11_t1_g27949 * (-(dx)));
            double cse_p60_t20_g31743 = (cse_p11_t1_g27949 * (-(dy)));
            double cse_p60_t21_g31744 = (cse_p11_t1_g27949 * (-(dz)));
            double cse_p60_t22_g31745 = (cse_p11_t1_g27949 * dx);
            double cse_p60_t23_g31746 = (cse_p11_t1_g27949 * dy);
            double cse_p60_t24_g31747 = (cse_p11_t1_g27949 * dz);
            double cse_p60_t41_g31764 = (cse_p60_t22_g31745 + (cse_p11_t1_g27947 * (cse_p11_t41_g27893 + (-(cse_p11_t1_g27948)))));
            double cse_p60_t47_g31770 = (cse_p60_t23_g31746 + (cse_p11_t1_g27947 * dx * dy));
            double cse_p60_t49_g31772 = (cse_p60_t20_g31743 + (cse_p11_t1_g27947 * dy * (-(dx))));
            double cse_p60_t50_g31773 = (cse_p60_t21_g31744 + cse_p60_t30_g31753);
            double cse_p60_t51_g31774 = (cse_p60_t24_g31747 + cse_p60_t32_g31755);
            double cse_p60_t67_g31790 = (cse_p11_t1_g27955 + cse_p60_t19_g31742);
            double h_x1_x2 = cse_p60_t41_g31764;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p11_t1_g27950 = (cse_p11_t1_g27946 * (-(dy)));
            double cse_p60_t25_g31748 = (cse_p11_t1_g27950 * (-(dy)));
            double cse_p60_t26_g31749 = (cse_p11_t1_g27950 * (-(dz)));
            double cse_p60_t27_g31750 = (cse_p11_t1_g27950 * dx);
            double cse_p60_t28_g31751 = (cse_p11_t1_g27950 * dy);
            double cse_p60_t29_g31752 = (cse_p11_t1_g27950 * dz);
            double cse_p60_t42_g31765 = (cse_p60_t28_g31751 + (cse_p11_t1_g27947 * (cse_p11_t42_g27894 + (-(cse_p11_t1_g27948)))));
            double cse_p60_t48_g31771 = (cse_p60_t27_g31750 + (cse_p11_t1_g27947 * dx * dy));
            double cse_p60_t52_g31775 = (cse_p60_t26_g31749 + cse_p60_t31_g31754);
            double cse_p60_t53_g31776 = (cse_p60_t29_g31752 + cse_p60_t33_g31756);
            double h_x1_y2 = cse_p60_t48_g31771;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p11_t1_g27952 = (cse_p11_t1_g27946 * (-(dz)));
            double cse_p60_t34_g31757 = (cse_p11_t1_g27952 * (-(dz)));
            double cse_p60_t35_g31758 = (cse_p11_t1_g27952 * dx);
            double cse_p60_t36_g31759 = (cse_p11_t1_g27952 * dy);
            double cse_p60_t37_g31760 = (cse_p11_t1_g27952 * dz);
            double cse_p60_t43_g31766 = (cse_p60_t37_g31760 + (cse_p11_t1_g27947 * (cse_p11_t43_g27895 + (-(cse_p11_t1_g27948)))));
            double cse_p60_t54_g31777 = (cse_p60_t32_g31755 + cse_p60_t35_g31758);
            double cse_p60_t55_g31778 = (cse_p60_t33_g31756 + cse_p60_t36_g31759);
            double h_x1_z2 = cse_p60_t54_g31777;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p11_t1_g27953 = (cse_p11_t1_g27947 * (cse_p11_t1_g27948 + (-(cse_p11_t42_g27894))));
            double cse_p60_t62_g31785 = (cse_p11_t1_g27953 + (cse_p11_t1_g27946 * dy * dy));
            double cse_p60_t63_g31786 = (cse_p11_t1_g27953 + cse_p60_t25_g31748);
            double h_y1_y1 = cse_p60_t62_g31785;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t46_g31769;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t47_g31770;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t42_g31765;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t55_g31778;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p11_t1_g27954 = (cse_p11_t1_g27947 * (cse_p11_t1_g27948 + (-(cse_p11_t43_g27895))));
            double cse_p60_t64_g31787 = (cse_p11_t1_g27954 + (cse_p11_t1_g27946 * dz * dz));
            double cse_p60_t65_g31788 = (cse_p11_t1_g27954 + cse_p60_t34_g31757);
            double h_z1_z1 = cse_p60_t64_g31787;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t51_g31774;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t53_g31776;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t43_g31766;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t67_g31790;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t49_g31772;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t50_g31773;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t63_g31786;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t52_g31775;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t65_g31788;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p12_t6_g27962 = (term.bb_high * rr);
              double cse_p11_t1_g27948 = (rr * rr);
              double cse_p12_t4_g27960 = (3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g27948);
              double cse_p11_t76_g27928 = (term.cc_high + cse_p12_t4_g27960 + cse_p12_t6_g27962 + cse_p12_t6_g27962);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p11_t1_g27948 * (term.bb_high + (term.aa_high * rr)))));
              energy_added += energy;
              *energy_accumulate += energy;
              double g_x1 = (cse_p11_t76_g27928 * cse_p1_t2_invsqrt27850 * dx * params.elec_weight);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t76_g27928 * cse_p1_t2_invsqrt27850 * dy * params.elec_weight);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t76_g27928 * cse_p1_t2_invsqrt27850 * dz * params.elec_weight);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = (params.elec_weight * ((cse_p1_t2_invsqrt27850 * dx * ((-2.00000000000000000e+0 * term.bb_high * rr) + (-3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g27948))) + (-((term.cc_high * cse_p1_t2_invsqrt27850 * dx)))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (params.elec_weight * ((cse_p1_t2_invsqrt27850 * dy * ((-2.00000000000000000e+0 * term.bb_high * rr) + (-3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g27948))) + (-((term.cc_high * cse_p1_t2_invsqrt27850 * dy)))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (params.elec_weight * ((cse_p1_t2_invsqrt27850 * dz * ((-2.00000000000000000e+0 * term.bb_high * rr) + (-3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g27948))) + (-((term.cc_high * cse_p1_t2_invsqrt27850 * dz)))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double de_dr = (params.elec_weight * (term.cc_high + cse_p12_t4_g27960 + (2.00000000000000000e+0 * term.bb_high * rr)));
              double d2e_dr2 = (params.elec_weight * ((2.00000000000000000e+0 * term.bb_high) + (6.00000000000000000e+0 * term.aa_high * rr)));
              double cse_p565_t2_invr67235 = (1.0 / (rr));
              double cse_p11_t1_g27946 = (cse_p565_t2_invr67235 * cse_p565_t2_invr67235 * d2e_dr2);
              double cse_p565_t1_invr67234 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g27947 = (cse_p565_t1_invr67234 * de_dr);
              double cse_p11_t1_g27955 = (cse_p11_t1_g27947 * (cse_p11_t1_g27948 + (-(cse_p11_t41_g27893))));
              double h_x1_x1 = (cse_p11_t1_g27955 + (cse_p11_t1_g27946 * dx * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double h_x1_y1 = ((cse_p11_t1_g27946 * dx * dy) + (cse_p11_t1_g27947 * dy * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g27951 = (cse_p11_t1_g27947 * dz);
              double h_x1_z1 = ((cse_p11_t1_g27946 * dx * dz) + (cse_p11_t1_g27951 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g27949 = (cse_p11_t1_g27946 * (-(dx)));
              double h_x1_x2 = ((cse_p11_t1_g27947 * (cse_p11_t41_g27893 + (-(cse_p11_t1_g27948)))) + (cse_p11_t1_g27949 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g27950 = (cse_p11_t1_g27946 * (-(dy)));
              double h_x1_y2 = ((cse_p11_t1_g27947 * dx * dy) + (cse_p11_t1_g27950 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g27952 = (cse_p11_t1_g27946 * (-(dz)));
              double h_x1_z2 = ((cse_p11_t1_g27951 * dx) + (cse_p11_t1_g27952 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g27953 = (cse_p11_t1_g27947 * (cse_p11_t1_g27948 + (-(cse_p11_t42_g27894))));
              double h_y1_y1 = (cse_p11_t1_g27953 + (cse_p11_t1_g27946 * dy * dy));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = ((cse_p11_t1_g27946 * dy * dz) + (cse_p11_t1_g27951 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = ((cse_p11_t1_g27947 * dx * dy) + (cse_p11_t1_g27949 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = ((cse_p11_t1_g27947 * (cse_p11_t42_g27894 + (-(cse_p11_t1_g27948)))) + (cse_p11_t1_g27950 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = ((cse_p11_t1_g27951 * dy) + (cse_p11_t1_g27952 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g27954 = (cse_p11_t1_g27947 * (cse_p11_t1_g27948 + (-(cse_p11_t43_g27895))));
              double h_z1_z1 = (cse_p11_t1_g27954 + (cse_p11_t1_g27946 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = ((cse_p11_t1_g27949 * dz) + (cse_p11_t1_g27951 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = ((cse_p11_t1_g27950 * dz) + (cse_p11_t1_g27951 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = ((cse_p11_t1_g27947 * (cse_p11_t43_g27895 + (-(cse_p11_t1_g27948)))) + (cse_p11_t1_g27952 * dz));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p11_t1_g27955 + (cse_p11_t1_g27949 * (-(dx))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = ((cse_p11_t1_g27947 * dy * (-(dx))) + (cse_p11_t1_g27949 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = ((cse_p11_t1_g27949 * (-(dz))) + (cse_p11_t1_g27951 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = (cse_p11_t1_g27953 + (cse_p11_t1_g27950 * (-(dy))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = ((cse_p11_t1_g27950 * (-(dz))) + (cse_p11_t1_g27951 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = (cse_p11_t1_g27954 + (cse_p11_t1_g27952 * (-(dz))));
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
                double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t41_g27893))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t2_invr67235 = (1.0 / (rr));
                double cse_p565_t1_invr67234 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g27948 = (rr * rr);
                double cse_p11_t1_g27947 = (cse_p565_t1_invr67234 * de_dr);
                double cse_p11_t1_g27946 = (cse_p565_t2_invr67235 * cse_p565_t2_invr67235 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g27946 * dx * dy) + (cse_p11_t1_g27947 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double h_x1_z1 = ((cse_p11_t1_g27946 * dx * dz) + (cse_p11_t1_g27947 * dz * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g27951 = (cse_p11_t1_g27947 * dz);
                double h_x1_x2 = ((cse_p11_t1_g27946 * dx * (-(dx))) + (cse_p11_t1_g27947 * (cse_p11_t41_g27893 + (-(cse_p11_t1_g27948)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g27949 = (cse_p11_t1_g27946 * (-(dx)));
                double h_x1_y2 = ((cse_p11_t1_g27946 * dx * (-(dy))) + (cse_p11_t1_g27947 * dx * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g27950 = (cse_p11_t1_g27946 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g27946 * dx * (-(dz))) + (cse_p11_t1_g27951 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g27952 = (cse_p11_t1_g27946 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g27946 * dy * dy) + (cse_p11_t1_g27947 * (cse_p11_t1_g27948 + (-(cse_p11_t42_g27894)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g27953 = (cse_p11_t1_g27947 * (cse_p11_t1_g27948 + (-(cse_p11_t42_g27894))));
                double h_y1_z1 = ((cse_p11_t1_g27946 * dy * dz) + (cse_p11_t1_g27951 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g27947 * dx * dy) + (cse_p11_t1_g27949 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g27947 * (cse_p11_t42_g27894 + (-(cse_p11_t1_g27948)))) + (cse_p11_t1_g27950 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g27951 * dy) + (cse_p11_t1_g27952 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g27946 * dz * dz) + (cse_p11_t1_g27947 * (cse_p11_t1_g27948 + (-(cse_p11_t43_g27895)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g27954 = (cse_p11_t1_g27947 * (cse_p11_t1_g27948 + (-(cse_p11_t43_g27895))));
                double h_z1_x2 = ((cse_p11_t1_g27949 * dz) + (cse_p11_t1_g27951 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g27950 * dz) + (cse_p11_t1_g27951 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g27947 * (cse_p11_t43_g27895 + (-(cse_p11_t1_g27948)))) + (cse_p11_t1_g27952 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g27947 * (cse_p11_t1_g27948 + (-(cse_p11_t41_g27893)))) + (cse_p11_t1_g27949 * (-(dx))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p11_t1_g27955 = (cse_p11_t1_g27947 * (cse_p11_t1_g27948 + (-(cse_p11_t41_g27893))));
                double h_x2_y2 = ((cse_p11_t1_g27947 * dy * (-(dx))) + (cse_p11_t1_g27949 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double h_x2_z2 = ((cse_p11_t1_g27949 * (-(dz))) + (cse_p11_t1_g27951 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g27953 + (cse_p11_t1_g27950 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g27950 * (-(dz))) + (cse_p11_t1_g27951 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g27954 + (cse_p11_t1_g27952 * (-(dz))));
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
