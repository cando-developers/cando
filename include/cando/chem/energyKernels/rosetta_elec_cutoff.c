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
      double cse_p1_t1_sqrt388386 = sqrt(r2);
      double cse_p1_t2_invsqrt388387 = (1.0 / (cse_p1_t1_sqrt388386));
      double rr = cse_p1_t1_sqrt388386;
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
          double cse_p528_t1_g395417 = (cse_p1_t2_invsqrt388387 * dx);
          double cse_p528_t2_g395405 = (cse_p1_t2_invsqrt388387 * dy);
          double cse_p528_t3_g395406 = (cse_p1_t2_invsqrt388387 * dz);
          double cse_p528_t4_g395407 = (rr * rr);
          double cse_p12_t1_g388427 = (term.aa_low * cse_p528_t4_g395407);
          double cse_p12_t4_g388419 = (3.00000000000000000e+0 * cse_p12_t1_g388427);
          double cse_p11_t1_g388411 = (term.bb_low * rr);
          double cse_p11_t12_g388401 = (term.cc_low + cse_p11_t1_g388411 + cse_p11_t1_g388411 + cse_p12_t4_g388419);
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p528_t4_g395407 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p12_t1_g388431 = (cse_p11_t12_g388401 * params.elec_weight);
          double g_x1 = (cse_p12_t1_g388431 * cse_p528_t1_g395417);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p12_t1_g388431 * cse_p528_t2_g395405);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p12_t1_g388431 * cse_p528_t3_g395406);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double cse_p63_t1_g389116 = (-2.00000000000000000e+0 * cse_p11_t1_g388411);
          double cse_p63_t1_g389118 = (-3.00000000000000000e+0 * cse_p12_t1_g388427);
          double g_x2 = (params.elec_weight * ((cse_p528_t1_g395417 * cse_p63_t1_g389116) + (cse_p528_t1_g395417 * cse_p63_t1_g389118) + (-((term.cc_low * cse_p528_t1_g395417)))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (params.elec_weight * ((cse_p528_t2_g395405 * cse_p63_t1_g389116) + (cse_p528_t2_g395405 * cse_p63_t1_g389118) + (-((term.cc_low * cse_p528_t2_g395405)))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (params.elec_weight * ((cse_p528_t3_g395406 * cse_p63_t1_g389116) + (cse_p528_t3_g395406 * cse_p63_t1_g389118) + (-((term.cc_low * cse_p528_t3_g395406)))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (params.elec_weight * (term.cc_low + cse_p12_t4_g388419 + (2.00000000000000000e+0 * cse_p11_t1_g388411)));
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t10_g388425 = (-(params.eps_solvent));
            double cse_p11_t13_g388402 = (cse_p12_t10_g388425 + params.eps_core);
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p12_t1_g388416 = (-1.25000000000000000e-1 * rx);
            double cse_p12_t2_g388417 = (1.25000000000000000e-1 * rx);
            double cse_p12_t3_g388418 = (5.00000000000000000e-1 * rx * rx);
            double cse_p12_t11_g388426 = (-(rx));
            double cse_p11_t8_g388397 = (-2.50000000000000000e-1 + cse_p12_t1_g388416 + cse_p12_t1_g388416);
            double cse_p11_t9_g388398 = (2.50000000000000000e-1 + cse_p12_t2_g388417 + cse_p12_t2_g388417);
            double cse_p11_t10_g388399 = (1.00000000000000000e+0 + cse_p12_t3_g388418 + rx);
            double cse_p11_t16_g388405 = exp(cse_p12_t11_g388426);
            double ex = cse_p11_t16_g388405;
            double gx = (cse_p11_t10_g388399 * ex);
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p12_t11_g388426 + cse_p12_t3_g388418));
            double ep = (params.eps_solvent + (cse_p11_t13_g388402 * gx));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t13_g388402 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t13_g388402 * gx_pp);
            double cse_p1_t3_invr388388 = (1.0 / (ep));
            double cse_p1_t4_invr2388389 = (cse_p1_t3_invr388388 * cse_p1_t3_invr388388);
            double ee = (cse_p1_t3_invr388388 * term.kqq);
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p528_t1_g395417 = (cse_p1_t2_invsqrt388387 * dx);
            double cse_p11_t1_g388412 = (cse_p11_t10_g388399 * cse_p11_t16_g388405);
            double cse_p11_t1_g388413 = (cse_p11_t13_g388402 * cse_p1_t4_invr2388389);
            double cse_p11_t1_g388414 = (params.elec_weight * term.kqq);
            double cse_p12_t1_g388428 = (-2.50000000000000000e-1 * cse_p11_t1_g388412);
            double cse_p12_t1_g388433 = (cse_p11_t9_g388398 * ex);
            double cse_p63_t1_g389114 = (cse_p11_t1_g388413 * cse_p11_t1_g388414);
            double g_x1 = (-((cse_p63_t1_g389114 * ((cse_p12_t1_g388428 * cse_p528_t1_g395417) + (cse_p12_t1_g388433 * cse_p528_t1_g395417)))));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p11_t1_g388407 = (cse_p1_t2_invsqrt388387 * dy);
            double g_y1 = (-((cse_p63_t1_g389114 * ((cse_p11_t1_g388407 * cse_p12_t1_g388428) + (cse_p11_t1_g388407 * cse_p12_t1_g388433)))));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p11_t1_g388408 = (cse_p1_t2_invsqrt388387 * dz);
            double g_z1 = (-((cse_p63_t1_g389114 * ((cse_p11_t1_g388408 * cse_p12_t1_g388428) + (cse_p11_t1_g388408 * cse_p12_t1_g388433)))));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double cse_p12_t1_g388429 = (2.50000000000000000e-1 * cse_p11_t1_g388412);
            double cse_p12_t1_g388432 = (cse_p11_t8_g388397 * ex);
            double g_x2 = (-((cse_p63_t1_g389114 * ((cse_p12_t1_g388429 * cse_p528_t1_g395417) + (cse_p12_t1_g388432 * cse_p528_t1_g395417)))));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (-((cse_p63_t1_g389114 * ((cse_p11_t1_g388407 * cse_p12_t1_g388429) + (cse_p11_t1_g388407 * cse_p12_t1_g388432)))));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (-((cse_p63_t1_g389114 * ((cse_p11_t1_g388408 * cse_p12_t1_g388429) + (cse_p11_t1_g388408 * cse_p12_t1_g388432)))));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (params.elec_weight * ep_p * (1.0 / ((ep * ep))) * (-(term.kqq)));
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p11_t1_g388409 = (rr * rr);
              double cse_p11_t1_g388415 = (term.aa_high * cse_p11_t1_g388409);
              double cse_p12_t5_g388420 = (3.00000000000000000e+0 * cse_p11_t1_g388415);
              double cse_p11_t1_g388410 = (term.bb_high * rr);
              double cse_p11_t11_g388400 = (term.cc_high + cse_p11_t1_g388410 + cse_p11_t1_g388410 + cse_p12_t5_g388420);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p11_t1_g388409 * (term.bb_high + (term.aa_high * rr)))));
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p528_t1_g395417 = (cse_p1_t2_invsqrt388387 * dx);
              double cse_p12_t1_g388430 = (cse_p11_t11_g388400 * params.elec_weight);
              double g_x1 = (cse_p12_t1_g388430 * cse_p528_t1_g395417);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double cse_p11_t1_g388407 = (cse_p1_t2_invsqrt388387 * dy);
              double g_y1 = (cse_p11_t1_g388407 * cse_p12_t1_g388430);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double cse_p11_t1_g388408 = (cse_p1_t2_invsqrt388387 * dz);
              double g_z1 = (cse_p11_t1_g388408 * cse_p12_t1_g388430);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p63_t1_g389115 = (-2.00000000000000000e+0 * cse_p11_t1_g388410);
              double cse_p63_t1_g389117 = (-3.00000000000000000e+0 * cse_p11_t1_g388415);
              double g_x2 = (params.elec_weight * ((cse_p528_t1_g395417 * cse_p63_t1_g389115) + (cse_p528_t1_g395417 * cse_p63_t1_g389117) + (-((term.cc_high * cse_p528_t1_g395417)))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (params.elec_weight * ((cse_p11_t1_g388407 * cse_p63_t1_g389115) + (cse_p11_t1_g388407 * cse_p63_t1_g389117) + (-((term.cc_high * cse_p11_t1_g388407)))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (params.elec_weight * ((cse_p11_t1_g388408 * cse_p63_t1_g389115) + (cse_p11_t1_g388408 * cse_p63_t1_g389117) + (-((term.cc_high * cse_p11_t1_g388408)))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double de_dr = (params.elec_weight * (term.cc_high + cse_p12_t5_g388420 + (2.00000000000000000e+0 * cse_p11_t1_g388410)));
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
          double cse_p528_t1_g400931 = (rr * rr);
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p528_t1_g400931 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t2_g395434 = (-(params.eps_solvent));
            double cse_p11_t2_g395420 = (cse_p12_t2_g395434 + params.eps_core);
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p11_t1_g395419 = (5.00000000000000000e-1 * rx * rx);
            double cse_p11_t4_g395422 = (-(rx));
            double ex = exp(cse_p11_t4_g395422);
            double gx = (ex * (1.00000000000000000e+0 + cse_p11_t1_g395419 + rx));
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p11_t1_g395419 + cse_p11_t4_g395422));
            double ep = (params.eps_solvent + (cse_p11_t2_g395420 * gx));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t2_g395420 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t2_g395420 * gx_pp);
            double ee = (term.kqq * (1.0 / (ep)));
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p528_t1_g400931 = (rr * rr);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p528_t1_g400931 * (term.bb_high + (term.aa_high * rr)))));
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
      double cse_p11_t41_g400977 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t42_g400978 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t43_g400979 = (dz * dz);
      double r2 = (cse_p11_t41_g400977 + cse_p11_t42_g400978 + cse_p11_t43_g400979);
      double cse_p1_t1_sqrt400933 = sqrt(r2);
      double cse_p1_t2_invsqrt400934 = (1.0 / (cse_p1_t1_sqrt400933));
      double rr = cse_p1_t1_sqrt400933;
      if ((rr <= params.rmin)) {
        /* !COND-1 */
        double cse_p60_t40_g404847 = (rr * rr);
        double cse_p60_t68_g404875 = (-(cse_p11_t41_g400977));
        double cse_p60_t69_g404876 = (-(cse_p11_t42_g400978));
        double cse_p60_t70_g404877 = (-(cse_p11_t43_g400979));
        double cse_p60_t72_g404879 = (-(dx));
        double cse_p60_t73_g404880 = (-(dy));
        double cse_p60_t74_g404881 = (-(dz));
        double cse_p565_t1_invr440318 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr440319 = (1.0 / (rr));
        double cse_p60_t77_g404884 = (rr * rr * rr);
        double cse_p15_t72_g401413 = (-(cse_p60_t40_g404847));
        double cse_p14_t76_g401317 = (cse_p11_t41_g400977 + cse_p15_t72_g401413);
        double cse_p14_t78_g401319 = (cse_p11_t42_g400978 + cse_p15_t72_g401413);
        double cse_p14_t80_g401321 = (cse_p11_t43_g400979 + cse_p15_t72_g401413);
        double cse_p13_t101_g401228 = (cse_p60_t40_g404847 + cse_p60_t68_g404875);
        double cse_p13_t102_g401229 = (cse_p60_t40_g404847 + cse_p60_t69_g404876);
        double cse_p13_t103_g401230 = (cse_p60_t40_g404847 + cse_p60_t70_g404877);
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
        double cse_p60_t39_g404846 = (cse_p565_t1_invr440318 * de_dr);
        double cse_p13_t52_g401179 = (cse_p14_t76_g401317 * cse_p60_t39_g404846);
        double cse_p13_t53_g401180 = (cse_p14_t78_g401319 * cse_p60_t39_g404846);
        double cse_p13_t54_g401181 = (cse_p14_t80_g401321 * cse_p60_t39_g404846);
        double cse_p13_t55_g401182 = (cse_p60_t39_g404846 * dx * dy);
        double cse_p13_t56_g401183 = (cse_p60_t39_g404846 * cse_p60_t72_g404879 * dy);
        double cse_p12_t20_g401060 = (cse_p13_t101_g401228 * cse_p60_t39_g404846);
        double cse_p12_t21_g401061 = (cse_p13_t102_g401229 * cse_p60_t39_g404846);
        double cse_p12_t22_g401062 = (cse_p13_t103_g401230 * cse_p60_t39_g404846);
        double cse_p12_t25_g401065 = (cse_p60_t39_g404846 * dz);
        double cse_p13_t42_g401169 = (cse_p12_t25_g401065 * cse_p60_t72_g404879);
        double cse_p13_t43_g401170 = (cse_p12_t25_g401065 * cse_p60_t73_g404880);
        double cse_p13_t44_g401171 = (cse_p12_t25_g401065 * dx);
        double cse_p13_t45_g401172 = (cse_p12_t25_g401065 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p60_t38_g404845 = (cse_p565_t2_invr440319 * cse_p565_t2_invr440319 * d2e_dr2);
        double cse_p13_t46_g401173 = (cse_p60_t38_g404845 * dx * dx);
        double cse_p13_t47_g401174 = (cse_p60_t38_g404845 * dx * dy);
        double cse_p13_t48_g401175 = (cse_p60_t38_g404845 * dx * dz);
        double cse_p13_t49_g401176 = (cse_p60_t38_g404845 * dy * dy);
        double cse_p13_t50_g401177 = (cse_p60_t38_g404845 * dy * dz);
        double cse_p13_t51_g401178 = (cse_p60_t38_g404845 * dz * dz);
        double cse_p12_t8_g401048 = (cse_p60_t38_g404845 * cse_p60_t72_g404879);
        double cse_p13_t57_g401184 = (cse_p12_t8_g401048 * cse_p60_t72_g404879);
        double cse_p13_t58_g401185 = (cse_p12_t8_g401048 * cse_p60_t73_g404880);
        double cse_p13_t59_g401186 = (cse_p12_t8_g401048 * cse_p60_t74_g404881);
        double cse_p13_t60_g401187 = (cse_p12_t8_g401048 * dx);
        double cse_p13_t61_g401188 = (cse_p12_t8_g401048 * dy);
        double cse_p13_t62_g401189 = (cse_p12_t8_g401048 * dz);
        double cse_p12_t9_g401049 = (cse_p60_t38_g404845 * cse_p60_t73_g404880);
        double cse_p13_t63_g401190 = (cse_p12_t9_g401049 * cse_p60_t73_g404880);
        double cse_p13_t64_g401191 = (cse_p12_t9_g401049 * cse_p60_t74_g404881);
        double cse_p13_t65_g401192 = (cse_p12_t9_g401049 * dx);
        double cse_p13_t66_g401193 = (cse_p12_t9_g401049 * dy);
        double cse_p13_t67_g401194 = (cse_p12_t9_g401049 * dz);
        double cse_p12_t10_g401050 = (cse_p60_t38_g404845 * cse_p60_t74_g404881);
        double cse_p13_t38_g401165 = (cse_p12_t10_g401050 * cse_p60_t74_g404881);
        double cse_p13_t39_g401166 = (cse_p12_t10_g401050 * dx);
        double cse_p13_t40_g401167 = (cse_p12_t10_g401050 * dy);
        double cse_p13_t41_g401168 = (cse_p12_t10_g401050 * dz);
        double cse_p12_t52_g401092 = (cse_p13_t47_g401174 + cse_p13_t56_g401183);
        double cse_p12_t53_g401093 = (cse_p13_t42_g401169 + cse_p13_t48_g401175);
        double cse_p12_t54_g401094 = (cse_p13_t43_g401170 + cse_p13_t50_g401177);
        double cse_p12_t70_g401110 = (cse_p12_t21_g401061 + cse_p13_t49_g401176);
        double cse_p12_t72_g401112 = (cse_p12_t22_g401062 + cse_p13_t51_g401178);
        double cse_p12_t74_g401114 = (cse_p12_t20_g401060 + cse_p13_t46_g401173);
        double cse_p12_t49_g401089 = (cse_p13_t52_g401179 + cse_p13_t60_g401187);
        double cse_p12_t55_g401095 = (cse_p13_t55_g401182 + cse_p13_t61_g401188);
        double cse_p12_t57_g401097 = (cse_p13_t56_g401183 + cse_p13_t58_g401185);
        double cse_p12_t58_g401098 = (cse_p13_t42_g401169 + cse_p13_t59_g401186);
        double cse_p12_t59_g401099 = (cse_p13_t44_g401171 + cse_p13_t62_g401189);
        double cse_p12_t75_g401115 = (cse_p12_t20_g401060 + cse_p13_t57_g401184);
        double cse_p12_t50_g401090 = (cse_p13_t53_g401180 + cse_p13_t66_g401193);
        double cse_p12_t56_g401096 = (cse_p13_t55_g401182 + cse_p13_t65_g401192);
        double cse_p12_t60_g401100 = (cse_p13_t43_g401170 + cse_p13_t64_g401191);
        double cse_p12_t61_g401101 = (cse_p13_t45_g401172 + cse_p13_t67_g401194);
        double cse_p12_t71_g401111 = (cse_p12_t21_g401061 + cse_p13_t63_g401190);
        double cse_p12_t51_g401091 = (cse_p13_t41_g401168 + cse_p13_t54_g401181);
        double cse_p12_t62_g401102 = (cse_p13_t39_g401166 + cse_p13_t44_g401171);
        double cse_p12_t63_g401103 = (cse_p13_t40_g401167 + cse_p13_t45_g401172);
        double cse_p12_t73_g401113 = (cse_p12_t22_g401062 + cse_p13_t38_g401165);
        double h_x1_x1 = cse_p12_t74_g401114;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t52_g401092;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t53_g401093;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t49_g401089;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t56_g401096;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t62_g401102;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t70_g401110;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t54_g401094;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t55_g401095;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t50_g401090;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t63_g401103;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t72_g401112;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t59_g401099;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t61_g401101;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t51_g401091;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t75_g401115;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t57_g401097;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t58_g401098;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t71_g401111;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t60_g401100;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t73_g401113;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= params.rlow)) {
          /* !COND-2 */
          double cse_p12_t7_g401047 = (term.bb_low * rr);
          double cse_p11_t1_g401032 = (rr * rr);
          double cse_p60_t56_g404863 = (cse_p11_t41_g400977 + (-(cse_p11_t1_g401032)));
          double cse_p60_t57_g404864 = (cse_p11_t42_g400978 + (-(cse_p11_t1_g401032)));
          double cse_p60_t58_g404865 = (cse_p11_t43_g400979 + (-(cse_p11_t1_g401032)));
          double cse_p60_t59_g404866 = (cse_p11_t1_g401032 + (-(cse_p11_t41_g400977)));
          double cse_p60_t60_g404867 = (cse_p11_t1_g401032 + (-(cse_p11_t42_g400978)));
          double cse_p60_t61_g404868 = (cse_p11_t1_g401032 + (-(cse_p11_t43_g400979)));
          double cse_p60_t71_g404878 = (-(cse_p11_t1_g401032));
          double cse_p12_t5_g401045 = (3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g401032);
          double cse_p11_t77_g401013 = (term.cc_low + cse_p12_t5_g401045 + cse_p12_t7_g401047 + cse_p12_t7_g401047);
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p11_t1_g401032 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (cse_p11_t77_g401013 * cse_p1_t2_invsqrt400934 * dx * params.elec_weight);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t77_g401013 * cse_p1_t2_invsqrt400934 * dy * params.elec_weight);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t77_g401013 * cse_p1_t2_invsqrt400934 * dz * params.elec_weight);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (params.elec_weight * ((cse_p1_t2_invsqrt400934 * dx * ((-2.00000000000000000e+0 * term.bb_low * rr) + (-3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g401032))) + (-((term.cc_low * cse_p1_t2_invsqrt400934 * dx)))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (params.elec_weight * ((cse_p1_t2_invsqrt400934 * dy * ((-2.00000000000000000e+0 * term.bb_low * rr) + (-3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g401032))) + (-((term.cc_low * cse_p1_t2_invsqrt400934 * dy)))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (params.elec_weight * ((cse_p1_t2_invsqrt400934 * dz * ((-2.00000000000000000e+0 * term.bb_low * rr) + (-3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g401032))) + (-((term.cc_low * cse_p1_t2_invsqrt400934 * dz)))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (params.elec_weight * (term.cc_low + cse_p12_t5_g401045 + (2.00000000000000000e+0 * term.bb_low * rr)));
          double d2e_dr2 = (params.elec_weight * ((2.00000000000000000e+0 * term.bb_low) + (6.00000000000000000e+0 * term.aa_low * rr)));
          double cse_p565_t2_invr440319 = (1.0 / (rr));
          double cse_p11_t1_g401030 = (cse_p565_t2_invr440319 * cse_p565_t2_invr440319 * d2e_dr2);
          double cse_p60_t1_g404808 = (cse_p11_t1_g401030 * (-(dx)));
          double cse_p60_t2_g404809 = (cse_p11_t1_g401030 * (-(dy)));
          double cse_p60_t3_g404810 = (cse_p11_t1_g401030 * (-(dz)));
          double cse_p60_t4_g404811 = (cse_p11_t1_g401030 * dx * dx);
          double cse_p60_t5_g404812 = (cse_p11_t1_g401030 * dx * dy);
          double cse_p60_t6_g404813 = (cse_p11_t1_g401030 * dx * dz);
          double cse_p60_t7_g404814 = (cse_p11_t1_g401030 * dy * dy);
          double cse_p60_t8_g404815 = (cse_p11_t1_g401030 * dy * dz);
          double cse_p60_t9_g404816 = (cse_p11_t1_g401030 * dz * dz);
          double cse_p14_t42_g401283 = (cse_p60_t1_g404808 * (-(dx)));
          double cse_p14_t43_g401284 = (cse_p60_t1_g404808 * (-(dy)));
          double cse_p14_t44_g401285 = (cse_p60_t1_g404808 * (-(dz)));
          double cse_p14_t45_g401286 = (cse_p60_t1_g404808 * dx);
          double cse_p14_t46_g401287 = (cse_p60_t1_g404808 * dy);
          double cse_p14_t47_g401288 = (cse_p60_t1_g404808 * dz);
          double cse_p14_t48_g401289 = (cse_p60_t2_g404809 * (-(dy)));
          double cse_p14_t49_g401290 = (cse_p60_t2_g404809 * (-(dz)));
          double cse_p14_t50_g401291 = (cse_p60_t2_g404809 * dx);
          double cse_p14_t51_g401292 = (cse_p60_t2_g404809 * dy);
          double cse_p14_t52_g401293 = (cse_p60_t2_g404809 * dz);
          double cse_p14_t53_g401294 = (cse_p60_t3_g404810 * (-(dz)));
          double cse_p14_t54_g401295 = (cse_p60_t3_g404810 * dx);
          double cse_p14_t55_g401296 = (cse_p60_t3_g404810 * dy);
          double cse_p14_t56_g401297 = (cse_p60_t3_g404810 * dz);
          double cse_p565_t1_invr440318 = (1.0 / ((rr * rr * rr)));
          double cse_p11_t1_g401031 = (cse_p565_t1_invr440318 * de_dr);
          double cse_p60_t10_g404817 = (cse_p11_t1_g401031 * (cse_p11_t41_g400977 + cse_p60_t71_g404878));
          double cse_p60_t11_g404818 = (cse_p11_t1_g401031 * (cse_p11_t42_g400978 + cse_p60_t71_g404878));
          double cse_p60_t12_g404819 = (cse_p11_t1_g401031 * (cse_p11_t43_g400979 + cse_p60_t71_g404878));
          double cse_p60_t13_g404820 = (cse_p11_t1_g401031 * cse_p60_t59_g404866);
          double cse_p60_t14_g404821 = (cse_p11_t1_g401031 * cse_p60_t60_g404867);
          double cse_p60_t15_g404822 = (cse_p11_t1_g401031 * cse_p60_t61_g404868);
          double cse_p60_t16_g404823 = (cse_p11_t1_g401031 * dx * dy);
          double cse_p60_t17_g404824 = (cse_p11_t1_g401031 * dy * (-(dx)));
          double cse_p60_t18_g404825 = (cse_p11_t1_g401031 * dz);
          double cse_p60_t44_g404851 = (cse_p60_t17_g404824 + cse_p60_t5_g404812);
          double cse_p14_t38_g401279 = (cse_p60_t18_g404825 * (-(dx)));
          double cse_p14_t39_g401280 = (cse_p60_t18_g404825 * (-(dy)));
          double cse_p14_t40_g401281 = (cse_p60_t18_g404825 * dx);
          double cse_p14_t41_g401282 = (cse_p60_t18_g404825 * dy);
          double cse_p13_t99_g401226 = (cse_p60_t13_g404820 + cse_p60_t4_g404811);
          double h_x1_x1 = cse_p13_t99_g401226;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t44_g404851;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t75_g401202 = (cse_p14_t38_g401279 + cse_p60_t6_g404813);
          double cse_p13_t76_g401203 = (cse_p14_t39_g401280 + cse_p60_t8_g404815);
          double h_x1_z1 = cse_p13_t75_g401202;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t71_g401198 = (cse_p14_t45_g401286 + cse_p60_t10_g404817);
          double cse_p13_t77_g401204 = (cse_p14_t46_g401287 + cse_p60_t16_g404823);
          double cse_p13_t79_g401206 = (cse_p14_t43_g401284 + cse_p60_t17_g404824);
          double cse_p13_t80_g401207 = (cse_p14_t38_g401279 + cse_p14_t44_g401285);
          double cse_p13_t81_g401208 = (cse_p14_t40_g401281 + cse_p14_t47_g401288);
          double cse_p13_t100_g401227 = (cse_p14_t42_g401283 + cse_p60_t13_g404820);
          double h_x1_x2 = cse_p13_t71_g401198;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t72_g401199 = (cse_p14_t51_g401292 + cse_p60_t11_g404818);
          double cse_p13_t78_g401205 = (cse_p14_t50_g401291 + cse_p60_t16_g404823);
          double cse_p13_t82_g401209 = (cse_p14_t39_g401280 + cse_p14_t49_g401290);
          double cse_p13_t83_g401210 = (cse_p14_t41_g401282 + cse_p14_t52_g401293);
          double h_x1_y2 = cse_p13_t78_g401205;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t73_g401200 = (cse_p14_t56_g401297 + cse_p60_t12_g404819);
          double cse_p13_t84_g401211 = (cse_p14_t40_g401281 + cse_p14_t54_g401295);
          double cse_p13_t85_g401212 = (cse_p14_t41_g401282 + cse_p14_t55_g401296);
          double h_x1_z2 = cse_p13_t84_g401211;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t95_g401222 = (cse_p60_t14_g404821 + cse_p60_t7_g404814);
          double cse_p13_t96_g401223 = (cse_p14_t48_g401289 + cse_p60_t14_g404821);
          double h_y1_y1 = cse_p13_t95_g401222;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t76_g401203;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t77_g401204;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t72_g401199;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t85_g401212;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t97_g401224 = (cse_p60_t15_g404822 + cse_p60_t9_g404816);
          double cse_p13_t98_g401225 = (cse_p14_t53_g401294 + cse_p60_t15_g404822);
          double h_z1_z1 = cse_p13_t97_g401224;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t81_g401208;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t83_g401210;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t73_g401200;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t100_g401227;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t79_g401206;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t80_g401207;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t96_g401223;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t82_g401209;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t98_g401225;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t83_g401123 = (-(params.eps_solvent));
            double cse_p11_t78_g401014 = (cse_p12_t83_g401123 + params.eps_core);
            double cse_p11_t87_g401023 = (-(term.kqq));
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p12_t1_g401041 = (-1.25000000000000000e-1 * rx);
            double cse_p12_t2_g401042 = (1.25000000000000000e-1 * rx);
            double cse_p12_t3_g401043 = (5.00000000000000000e-1 * rx * rx);
            double cse_p12_t84_g401124 = (-(rx));
            double cse_p11_t73_g401009 = (-2.50000000000000000e-1 + cse_p12_t1_g401041 + cse_p12_t1_g401041);
            double cse_p11_t74_g401010 = (2.50000000000000000e-1 + cse_p12_t2_g401042 + cse_p12_t2_g401042);
            double cse_p11_t75_g401011 = (1.00000000000000000e+0 + cse_p12_t3_g401043 + rx);
            double cse_p11_t89_g401025 = exp(cse_p12_t84_g401124);
            double ex = cse_p11_t89_g401025;
            double gx = (cse_p11_t75_g401011 * ex);
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p12_t3_g401043 + cse_p12_t84_g401124));
            double ep = (params.eps_solvent + (cse_p11_t78_g401014 * gx));
            double cse_p12_t47_g401087 = (ep * ep);
            double cse_p11_t90_g401026 = (1.0 / (cse_p12_t47_g401087));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t78_g401014 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t78_g401014 * gx_pp);
            double cse_p1_t3_invr400935 = (1.0 / (ep));
            double cse_p1_t4_invr2400936 = (cse_p1_t3_invr400935 * cse_p1_t3_invr400935);
            double ee = (cse_p1_t3_invr400935 * term.kqq);
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p11_t1_g401040 = (cse_p11_t75_g401011 * cse_p11_t89_g401025 * cse_p1_t2_invsqrt400934);
            double g_x1 = (-((cse_p11_t78_g401014 * cse_p1_t4_invr2400936 * params.elec_weight * term.kqq * ((-2.50000000000000000e-1 * cse_p11_t1_g401040 * dx) + (cse_p11_t74_g401010 * cse_p1_t2_invsqrt400934 * dx * ex)))));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double g_y1 = (-((cse_p11_t78_g401014 * cse_p1_t4_invr2400936 * params.elec_weight * term.kqq * ((-2.50000000000000000e-1 * cse_p11_t1_g401040 * dy) + (cse_p11_t74_g401010 * cse_p1_t2_invsqrt400934 * dy * ex)))));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double g_z1 = (-((cse_p11_t78_g401014 * cse_p1_t4_invr2400936 * params.elec_weight * term.kqq * ((-2.50000000000000000e-1 * cse_p11_t1_g401040 * dz) + (cse_p11_t74_g401010 * cse_p1_t2_invsqrt400934 * dz * ex)))));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (-((cse_p11_t78_g401014 * cse_p1_t4_invr2400936 * params.elec_weight * term.kqq * ((2.50000000000000000e-1 * cse_p11_t1_g401040 * dx) + (cse_p11_t73_g401009 * cse_p1_t2_invsqrt400934 * dx * ex)))));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (-((cse_p11_t78_g401014 * cse_p1_t4_invr2400936 * params.elec_weight * term.kqq * ((2.50000000000000000e-1 * cse_p11_t1_g401040 * dy) + (cse_p11_t73_g401009 * cse_p1_t2_invsqrt400934 * dy * ex)))));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (-((cse_p11_t78_g401014 * cse_p1_t4_invr2400936 * params.elec_weight * term.kqq * ((2.50000000000000000e-1 * cse_p11_t1_g401040 * dz) + (cse_p11_t73_g401009 * cse_p1_t2_invsqrt400934 * dz * ex)))));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (cse_p11_t87_g401023 * cse_p11_t90_g401026 * params.elec_weight * ep_p);
            double d2e_dr2 = (cse_p11_t87_g401023 * params.elec_weight * ((-2.00000000000000000e+0 * ep_p * ep_p * (1.0 / ((ep * ep * ep)))) + (cse_p11_t90_g401026 * ep_pp)));
            double cse_p565_t2_invr440319 = (1.0 / (rr));
            double cse_p11_t1_g401030 = (cse_p565_t2_invr440319 * cse_p565_t2_invr440319 * d2e_dr2);
            double cse_p565_t1_invr440318 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g401031 = (cse_p565_t1_invr440318 * de_dr);
            double cse_p11_t1_g401032 = (rr * rr);
            double cse_p11_t1_g401039 = (cse_p11_t1_g401031 * (cse_p11_t1_g401032 + (-(cse_p11_t41_g400977))));
            double cse_p60_t66_g404873 = (cse_p11_t1_g401039 + (cse_p11_t1_g401030 * dx * dx));
            double h_x1_x1 = cse_p60_t66_g404873;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = ((cse_p11_t1_g401030 * dx * dy) + (cse_p11_t1_g401031 * dy * (-(dx))));
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p11_t1_g401035 = (cse_p11_t1_g401031 * dz);
            double cse_p60_t30_g404837 = (cse_p11_t1_g401035 * (-(dx)));
            double cse_p60_t31_g404838 = (cse_p11_t1_g401035 * (-(dy)));
            double cse_p60_t32_g404839 = (cse_p11_t1_g401035 * dx);
            double cse_p60_t33_g404840 = (cse_p11_t1_g401035 * dy);
            double cse_p60_t45_g404852 = (cse_p60_t30_g404837 + (cse_p11_t1_g401030 * dx * dz));
            double cse_p60_t46_g404853 = (cse_p60_t31_g404838 + (cse_p11_t1_g401030 * dy * dz));
            double h_x1_z1 = cse_p60_t45_g404852;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p11_t1_g401033 = (cse_p11_t1_g401030 * (-(dx)));
            double cse_p60_t19_g404826 = (cse_p11_t1_g401033 * (-(dx)));
            double cse_p60_t20_g404827 = (cse_p11_t1_g401033 * (-(dy)));
            double cse_p60_t21_g404828 = (cse_p11_t1_g401033 * (-(dz)));
            double cse_p60_t22_g404829 = (cse_p11_t1_g401033 * dx);
            double cse_p60_t23_g404830 = (cse_p11_t1_g401033 * dy);
            double cse_p60_t24_g404831 = (cse_p11_t1_g401033 * dz);
            double cse_p60_t41_g404848 = (cse_p60_t22_g404829 + (cse_p11_t1_g401031 * (cse_p11_t41_g400977 + (-(cse_p11_t1_g401032)))));
            double cse_p60_t47_g404854 = (cse_p60_t23_g404830 + (cse_p11_t1_g401031 * dx * dy));
            double cse_p60_t49_g404856 = (cse_p60_t20_g404827 + (cse_p11_t1_g401031 * dy * (-(dx))));
            double cse_p60_t50_g404857 = (cse_p60_t21_g404828 + cse_p60_t30_g404837);
            double cse_p60_t51_g404858 = (cse_p60_t24_g404831 + cse_p60_t32_g404839);
            double cse_p60_t67_g404874 = (cse_p11_t1_g401039 + cse_p60_t19_g404826);
            double h_x1_x2 = cse_p60_t41_g404848;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p11_t1_g401034 = (cse_p11_t1_g401030 * (-(dy)));
            double cse_p60_t25_g404832 = (cse_p11_t1_g401034 * (-(dy)));
            double cse_p60_t26_g404833 = (cse_p11_t1_g401034 * (-(dz)));
            double cse_p60_t27_g404834 = (cse_p11_t1_g401034 * dx);
            double cse_p60_t28_g404835 = (cse_p11_t1_g401034 * dy);
            double cse_p60_t29_g404836 = (cse_p11_t1_g401034 * dz);
            double cse_p60_t42_g404849 = (cse_p60_t28_g404835 + (cse_p11_t1_g401031 * (cse_p11_t42_g400978 + (-(cse_p11_t1_g401032)))));
            double cse_p60_t48_g404855 = (cse_p60_t27_g404834 + (cse_p11_t1_g401031 * dx * dy));
            double cse_p60_t52_g404859 = (cse_p60_t26_g404833 + cse_p60_t31_g404838);
            double cse_p60_t53_g404860 = (cse_p60_t29_g404836 + cse_p60_t33_g404840);
            double h_x1_y2 = cse_p60_t48_g404855;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p11_t1_g401036 = (cse_p11_t1_g401030 * (-(dz)));
            double cse_p60_t34_g404841 = (cse_p11_t1_g401036 * (-(dz)));
            double cse_p60_t35_g404842 = (cse_p11_t1_g401036 * dx);
            double cse_p60_t36_g404843 = (cse_p11_t1_g401036 * dy);
            double cse_p60_t37_g404844 = (cse_p11_t1_g401036 * dz);
            double cse_p60_t43_g404850 = (cse_p60_t37_g404844 + (cse_p11_t1_g401031 * (cse_p11_t43_g400979 + (-(cse_p11_t1_g401032)))));
            double cse_p60_t54_g404861 = (cse_p60_t32_g404839 + cse_p60_t35_g404842);
            double cse_p60_t55_g404862 = (cse_p60_t33_g404840 + cse_p60_t36_g404843);
            double h_x1_z2 = cse_p60_t54_g404861;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p11_t1_g401037 = (cse_p11_t1_g401031 * (cse_p11_t1_g401032 + (-(cse_p11_t42_g400978))));
            double cse_p60_t62_g404869 = (cse_p11_t1_g401037 + (cse_p11_t1_g401030 * dy * dy));
            double cse_p60_t63_g404870 = (cse_p11_t1_g401037 + cse_p60_t25_g404832);
            double h_y1_y1 = cse_p60_t62_g404869;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t46_g404853;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t47_g404854;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t42_g404849;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t55_g404862;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p11_t1_g401038 = (cse_p11_t1_g401031 * (cse_p11_t1_g401032 + (-(cse_p11_t43_g400979))));
            double cse_p60_t64_g404871 = (cse_p11_t1_g401038 + (cse_p11_t1_g401030 * dz * dz));
            double cse_p60_t65_g404872 = (cse_p11_t1_g401038 + cse_p60_t34_g404841);
            double h_z1_z1 = cse_p60_t64_g404871;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t51_g404858;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t53_g404860;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t43_g404850;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t67_g404874;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t49_g404856;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t50_g404857;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t63_g404870;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t52_g404859;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t65_g404872;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p12_t6_g401046 = (term.bb_high * rr);
              double cse_p11_t1_g401032 = (rr * rr);
              double cse_p12_t4_g401044 = (3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g401032);
              double cse_p11_t76_g401012 = (term.cc_high + cse_p12_t4_g401044 + cse_p12_t6_g401046 + cse_p12_t6_g401046);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p11_t1_g401032 * (term.bb_high + (term.aa_high * rr)))));
              energy_added += energy;
              *energy_accumulate += energy;
              double g_x1 = (cse_p11_t76_g401012 * cse_p1_t2_invsqrt400934 * dx * params.elec_weight);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t76_g401012 * cse_p1_t2_invsqrt400934 * dy * params.elec_weight);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t76_g401012 * cse_p1_t2_invsqrt400934 * dz * params.elec_weight);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = (params.elec_weight * ((cse_p1_t2_invsqrt400934 * dx * ((-2.00000000000000000e+0 * term.bb_high * rr) + (-3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g401032))) + (-((term.cc_high * cse_p1_t2_invsqrt400934 * dx)))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (params.elec_weight * ((cse_p1_t2_invsqrt400934 * dy * ((-2.00000000000000000e+0 * term.bb_high * rr) + (-3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g401032))) + (-((term.cc_high * cse_p1_t2_invsqrt400934 * dy)))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (params.elec_weight * ((cse_p1_t2_invsqrt400934 * dz * ((-2.00000000000000000e+0 * term.bb_high * rr) + (-3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g401032))) + (-((term.cc_high * cse_p1_t2_invsqrt400934 * dz)))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double de_dr = (params.elec_weight * (term.cc_high + cse_p12_t4_g401044 + (2.00000000000000000e+0 * term.bb_high * rr)));
              double d2e_dr2 = (params.elec_weight * ((2.00000000000000000e+0 * term.bb_high) + (6.00000000000000000e+0 * term.aa_high * rr)));
              double cse_p565_t2_invr440319 = (1.0 / (rr));
              double cse_p11_t1_g401030 = (cse_p565_t2_invr440319 * cse_p565_t2_invr440319 * d2e_dr2);
              double cse_p565_t1_invr440318 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g401031 = (cse_p565_t1_invr440318 * de_dr);
              double cse_p11_t1_g401039 = (cse_p11_t1_g401031 * (cse_p11_t1_g401032 + (-(cse_p11_t41_g400977))));
              double h_x1_x1 = (cse_p11_t1_g401039 + (cse_p11_t1_g401030 * dx * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double h_x1_y1 = ((cse_p11_t1_g401030 * dx * dy) + (cse_p11_t1_g401031 * dy * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g401035 = (cse_p11_t1_g401031 * dz);
              double h_x1_z1 = ((cse_p11_t1_g401030 * dx * dz) + (cse_p11_t1_g401035 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g401033 = (cse_p11_t1_g401030 * (-(dx)));
              double h_x1_x2 = ((cse_p11_t1_g401031 * (cse_p11_t41_g400977 + (-(cse_p11_t1_g401032)))) + (cse_p11_t1_g401033 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g401034 = (cse_p11_t1_g401030 * (-(dy)));
              double h_x1_y2 = ((cse_p11_t1_g401031 * dx * dy) + (cse_p11_t1_g401034 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g401036 = (cse_p11_t1_g401030 * (-(dz)));
              double h_x1_z2 = ((cse_p11_t1_g401035 * dx) + (cse_p11_t1_g401036 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g401037 = (cse_p11_t1_g401031 * (cse_p11_t1_g401032 + (-(cse_p11_t42_g400978))));
              double h_y1_y1 = (cse_p11_t1_g401037 + (cse_p11_t1_g401030 * dy * dy));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = ((cse_p11_t1_g401030 * dy * dz) + (cse_p11_t1_g401035 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = ((cse_p11_t1_g401031 * dx * dy) + (cse_p11_t1_g401033 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = ((cse_p11_t1_g401031 * (cse_p11_t42_g400978 + (-(cse_p11_t1_g401032)))) + (cse_p11_t1_g401034 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = ((cse_p11_t1_g401035 * dy) + (cse_p11_t1_g401036 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g401038 = (cse_p11_t1_g401031 * (cse_p11_t1_g401032 + (-(cse_p11_t43_g400979))));
              double h_z1_z1 = (cse_p11_t1_g401038 + (cse_p11_t1_g401030 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = ((cse_p11_t1_g401033 * dz) + (cse_p11_t1_g401035 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = ((cse_p11_t1_g401034 * dz) + (cse_p11_t1_g401035 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = ((cse_p11_t1_g401031 * (cse_p11_t43_g400979 + (-(cse_p11_t1_g401032)))) + (cse_p11_t1_g401036 * dz));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p11_t1_g401039 + (cse_p11_t1_g401033 * (-(dx))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = ((cse_p11_t1_g401031 * dy * (-(dx))) + (cse_p11_t1_g401033 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = ((cse_p11_t1_g401033 * (-(dz))) + (cse_p11_t1_g401035 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = (cse_p11_t1_g401037 + (cse_p11_t1_g401034 * (-(dy))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = ((cse_p11_t1_g401034 * (-(dz))) + (cse_p11_t1_g401035 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = (cse_p11_t1_g401038 + (cse_p11_t1_g401036 * (-(dz))));
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
                double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t41_g400977))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t2_invr440319 = (1.0 / (rr));
                double cse_p565_t1_invr440318 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g401032 = (rr * rr);
                double cse_p11_t1_g401031 = (cse_p565_t1_invr440318 * de_dr);
                double cse_p11_t1_g401030 = (cse_p565_t2_invr440319 * cse_p565_t2_invr440319 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g401030 * dx * dy) + (cse_p11_t1_g401031 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double h_x1_z1 = ((cse_p11_t1_g401030 * dx * dz) + (cse_p11_t1_g401031 * dz * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g401035 = (cse_p11_t1_g401031 * dz);
                double h_x1_x2 = ((cse_p11_t1_g401030 * dx * (-(dx))) + (cse_p11_t1_g401031 * (cse_p11_t41_g400977 + (-(cse_p11_t1_g401032)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g401033 = (cse_p11_t1_g401030 * (-(dx)));
                double h_x1_y2 = ((cse_p11_t1_g401030 * dx * (-(dy))) + (cse_p11_t1_g401031 * dx * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g401034 = (cse_p11_t1_g401030 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g401030 * dx * (-(dz))) + (cse_p11_t1_g401035 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g401036 = (cse_p11_t1_g401030 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g401030 * dy * dy) + (cse_p11_t1_g401031 * (cse_p11_t1_g401032 + (-(cse_p11_t42_g400978)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g401037 = (cse_p11_t1_g401031 * (cse_p11_t1_g401032 + (-(cse_p11_t42_g400978))));
                double h_y1_z1 = ((cse_p11_t1_g401030 * dy * dz) + (cse_p11_t1_g401035 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g401031 * dx * dy) + (cse_p11_t1_g401033 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g401031 * (cse_p11_t42_g400978 + (-(cse_p11_t1_g401032)))) + (cse_p11_t1_g401034 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g401035 * dy) + (cse_p11_t1_g401036 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g401030 * dz * dz) + (cse_p11_t1_g401031 * (cse_p11_t1_g401032 + (-(cse_p11_t43_g400979)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g401038 = (cse_p11_t1_g401031 * (cse_p11_t1_g401032 + (-(cse_p11_t43_g400979))));
                double h_z1_x2 = ((cse_p11_t1_g401033 * dz) + (cse_p11_t1_g401035 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g401034 * dz) + (cse_p11_t1_g401035 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g401031 * (cse_p11_t43_g400979 + (-(cse_p11_t1_g401032)))) + (cse_p11_t1_g401036 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g401031 * (cse_p11_t1_g401032 + (-(cse_p11_t41_g400977)))) + (cse_p11_t1_g401033 * (-(dx))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p11_t1_g401039 = (cse_p11_t1_g401031 * (cse_p11_t1_g401032 + (-(cse_p11_t41_g400977))));
                double h_x2_y2 = ((cse_p11_t1_g401031 * dy * (-(dx))) + (cse_p11_t1_g401033 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double h_x2_z2 = ((cse_p11_t1_g401033 * (-(dz))) + (cse_p11_t1_g401035 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g401037 + (cse_p11_t1_g401034 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g401034 * (-(dz))) + (cse_p11_t1_g401035 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g401038 + (cse_p11_t1_g401036 * (-(dz))));
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
