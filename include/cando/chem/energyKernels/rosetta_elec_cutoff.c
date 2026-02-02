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
      double cse_p1_t1_sqrt280 = sqrt(r2);
      double cse_p1_t2_invsqrt281 = (1.0 / (cse_p1_t1_sqrt280));
      double rr = cse_p1_t1_sqrt280;
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
          double cse_p528_t1_g7343 = (cse_p1_t2_invsqrt281 * dx);
          double cse_p528_t2_g7331 = (cse_p1_t2_invsqrt281 * dy);
          double cse_p528_t3_g7332 = (cse_p1_t2_invsqrt281 * dz);
          double cse_p528_t4_g7333 = (rr * rr);
          double cse_p12_t1_g353 = (term.aa_low * cse_p528_t4_g7333);
          double cse_p12_t4_g345 = (3.00000000000000000e+0 * cse_p12_t1_g353);
          double cse_p11_t1_g331 = (term.bb_low * rr);
          double cse_p11_t12_g305 = (term.cc_low + cse_p11_t1_g331 + cse_p11_t1_g331 + cse_p12_t4_g345);
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p528_t4_g7333 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p12_t1_g357 = (cse_p11_t12_g305 * params.elec_weight);
          double g_x1 = (cse_p12_t1_g357 * cse_p528_t1_g7343);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p12_t1_g357 * cse_p528_t2_g7331);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p12_t1_g357 * cse_p528_t3_g7332);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double cse_p63_t1_g1042 = (-2.00000000000000000e+0 * cse_p11_t1_g331);
          double cse_p63_t1_g1044 = (-3.00000000000000000e+0 * cse_p12_t1_g353);
          double g_x2 = (params.elec_weight * ((cse_p528_t1_g7343 * cse_p63_t1_g1042) + (cse_p528_t1_g7343 * cse_p63_t1_g1044) + (-((term.cc_low * cse_p528_t1_g7343)))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (params.elec_weight * ((cse_p528_t2_g7331 * cse_p63_t1_g1042) + (cse_p528_t2_g7331 * cse_p63_t1_g1044) + (-((term.cc_low * cse_p528_t2_g7331)))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (params.elec_weight * ((cse_p528_t3_g7332 * cse_p63_t1_g1042) + (cse_p528_t3_g7332 * cse_p63_t1_g1044) + (-((term.cc_low * cse_p528_t3_g7332)))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (params.elec_weight * (term.cc_low + cse_p12_t4_g345 + (2.00000000000000000e+0 * cse_p11_t1_g331)));
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t10_g351 = (-(params.eps_solvent));
            double cse_p11_t13_g306 = (cse_p12_t10_g351 + params.eps_core);
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p12_t1_g342 = (-1.25000000000000000e-1 * rx);
            double cse_p12_t2_g343 = (1.25000000000000000e-1 * rx);
            double cse_p12_t3_g344 = (5.00000000000000000e-1 * rx * rx);
            double cse_p12_t11_g352 = (-(rx));
            double cse_p11_t8_g301 = (-2.50000000000000000e-1 + cse_p12_t1_g342 + cse_p12_t1_g342);
            double cse_p11_t9_g302 = (2.50000000000000000e-1 + cse_p12_t2_g343 + cse_p12_t2_g343);
            double cse_p11_t10_g303 = (1.00000000000000000e+0 + cse_p12_t3_g344 + rx);
            double cse_p11_t16_g309 = exp(cse_p12_t11_g352);
            double ex = cse_p11_t16_g309;
            double gx = (cse_p11_t10_g303 * ex);
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p12_t11_g352 + cse_p12_t3_g344));
            double ep = (params.eps_solvent + (cse_p11_t13_g306 * gx));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t13_g306 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t13_g306 * gx_pp);
            double cse_p1_t3_invr282 = (1.0 / (ep));
            double cse_p1_t4_invr2283 = (cse_p1_t3_invr282 * cse_p1_t3_invr282);
            double ee = (cse_p1_t3_invr282 * term.kqq);
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p528_t1_g7343 = (cse_p1_t2_invsqrt281 * dx);
            double cse_p11_t1_g332 = (cse_p11_t10_g303 * cse_p11_t16_g309);
            double cse_p11_t1_g333 = (cse_p11_t13_g306 * cse_p1_t4_invr2283);
            double cse_p11_t1_g334 = (params.elec_weight * term.kqq);
            double cse_p12_t1_g354 = (-2.50000000000000000e-1 * cse_p11_t1_g332);
            double cse_p12_t1_g359 = (cse_p11_t9_g302 * ex);
            double cse_p63_t1_g1040 = (cse_p11_t1_g333 * cse_p11_t1_g334);
            double g_x1 = (-((cse_p63_t1_g1040 * ((cse_p12_t1_g354 * cse_p528_t1_g7343) + (cse_p12_t1_g359 * cse_p528_t1_g7343)))));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p11_t1_g327 = (cse_p1_t2_invsqrt281 * dy);
            double g_y1 = (-((cse_p63_t1_g1040 * ((cse_p11_t1_g327 * cse_p12_t1_g354) + (cse_p11_t1_g327 * cse_p12_t1_g359)))));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p11_t1_g328 = (cse_p1_t2_invsqrt281 * dz);
            double g_z1 = (-((cse_p63_t1_g1040 * ((cse_p11_t1_g328 * cse_p12_t1_g354) + (cse_p11_t1_g328 * cse_p12_t1_g359)))));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double cse_p12_t1_g355 = (2.50000000000000000e-1 * cse_p11_t1_g332);
            double cse_p12_t1_g358 = (cse_p11_t8_g301 * ex);
            double g_x2 = (-((cse_p63_t1_g1040 * ((cse_p12_t1_g355 * cse_p528_t1_g7343) + (cse_p12_t1_g358 * cse_p528_t1_g7343)))));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (-((cse_p63_t1_g1040 * ((cse_p11_t1_g327 * cse_p12_t1_g355) + (cse_p11_t1_g327 * cse_p12_t1_g358)))));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (-((cse_p63_t1_g1040 * ((cse_p11_t1_g328 * cse_p12_t1_g355) + (cse_p11_t1_g328 * cse_p12_t1_g358)))));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (params.elec_weight * ep_p * (1.0 / ((ep * ep))) * (-(term.kqq)));
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p11_t1_g329 = (rr * rr);
              double cse_p11_t1_g335 = (term.aa_high * cse_p11_t1_g329);
              double cse_p12_t5_g346 = (3.00000000000000000e+0 * cse_p11_t1_g335);
              double cse_p11_t1_g330 = (term.bb_high * rr);
              double cse_p11_t11_g304 = (term.cc_high + cse_p11_t1_g330 + cse_p11_t1_g330 + cse_p12_t5_g346);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p11_t1_g329 * (term.bb_high + (term.aa_high * rr)))));
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p528_t1_g7343 = (cse_p1_t2_invsqrt281 * dx);
              double cse_p12_t1_g356 = (cse_p11_t11_g304 * params.elec_weight);
              double g_x1 = (cse_p12_t1_g356 * cse_p528_t1_g7343);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double cse_p11_t1_g327 = (cse_p1_t2_invsqrt281 * dy);
              double g_y1 = (cse_p11_t1_g327 * cse_p12_t1_g356);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double cse_p11_t1_g328 = (cse_p1_t2_invsqrt281 * dz);
              double g_z1 = (cse_p11_t1_g328 * cse_p12_t1_g356);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p63_t1_g1041 = (-2.00000000000000000e+0 * cse_p11_t1_g330);
              double cse_p63_t1_g1043 = (-3.00000000000000000e+0 * cse_p11_t1_g335);
              double g_x2 = (params.elec_weight * ((cse_p528_t1_g7343 * cse_p63_t1_g1041) + (cse_p528_t1_g7343 * cse_p63_t1_g1043) + (-((term.cc_high * cse_p528_t1_g7343)))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (params.elec_weight * ((cse_p11_t1_g327 * cse_p63_t1_g1041) + (cse_p11_t1_g327 * cse_p63_t1_g1043) + (-((term.cc_high * cse_p11_t1_g327)))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (params.elec_weight * ((cse_p11_t1_g328 * cse_p63_t1_g1041) + (cse_p11_t1_g328 * cse_p63_t1_g1043) + (-((term.cc_high * cse_p11_t1_g328)))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double de_dr = (params.elec_weight * (term.cc_high + cse_p12_t5_g346 + (2.00000000000000000e+0 * cse_p11_t1_g330)));
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
          double cse_p528_t1_g12857 = (rr * rr);
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p528_t1_g12857 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t2_g7360 = (-(params.eps_solvent));
            double cse_p11_t2_g7346 = (cse_p12_t2_g7360 + params.eps_core);
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p11_t1_g7345 = (5.00000000000000000e-1 * rx * rx);
            double cse_p11_t4_g7348 = (-(rx));
            double ex = exp(cse_p11_t4_g7348);
            double gx = (ex * (1.00000000000000000e+0 + cse_p11_t1_g7345 + rx));
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p11_t1_g7345 + cse_p11_t4_g7348));
            double ep = (params.eps_solvent + (cse_p11_t2_g7346 * gx));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t2_g7346 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t2_g7346 * gx_pp);
            double ee = (term.kqq * (1.0 / (ep)));
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p528_t1_g12857 = (rr * rr);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p528_t1_g12857 * (term.bb_high + (term.aa_high * rr)))));
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
      double cse_p11_t41_g12909 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t42_g12910 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t43_g12911 = (dz * dz);
      double r2 = (cse_p11_t41_g12909 + cse_p11_t42_g12910 + cse_p11_t43_g12911);
      double cse_p1_t1_sqrt12859 = sqrt(r2);
      double cse_p1_t2_invsqrt12860 = (1.0 / (cse_p1_t1_sqrt12859));
      double rr = cse_p1_t1_sqrt12859;
      if ((rr <= params.rmin)) {
        /* !COND-1 */
        double cse_p60_t40_g16779 = (rr * rr);
        double cse_p60_t68_g16807 = (-(cse_p11_t41_g12909));
        double cse_p60_t69_g16808 = (-(cse_p11_t42_g12910));
        double cse_p60_t70_g16809 = (-(cse_p11_t43_g12911));
        double cse_p60_t72_g16811 = (-(dx));
        double cse_p60_t73_g16812 = (-(dy));
        double cse_p60_t74_g16813 = (-(dz));
        double cse_p565_t1_invr52250 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr52251 = (1.0 / (rr));
        double cse_p60_t77_g16816 = (rr * rr * rr);
        double cse_p15_t72_g13345 = (-(cse_p60_t40_g16779));
        double cse_p14_t76_g13249 = (cse_p11_t41_g12909 + cse_p15_t72_g13345);
        double cse_p14_t78_g13251 = (cse_p11_t42_g12910 + cse_p15_t72_g13345);
        double cse_p14_t80_g13253 = (cse_p11_t43_g12911 + cse_p15_t72_g13345);
        double cse_p13_t101_g13160 = (cse_p60_t40_g16779 + cse_p60_t68_g16807);
        double cse_p13_t102_g13161 = (cse_p60_t40_g16779 + cse_p60_t69_g16808);
        double cse_p13_t103_g13162 = (cse_p60_t40_g16779 + cse_p60_t70_g16809);
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
        double cse_p60_t39_g16778 = (cse_p565_t1_invr52250 * de_dr);
        double cse_p13_t52_g13111 = (cse_p14_t76_g13249 * cse_p60_t39_g16778);
        double cse_p13_t53_g13112 = (cse_p14_t78_g13251 * cse_p60_t39_g16778);
        double cse_p13_t54_g13113 = (cse_p14_t80_g13253 * cse_p60_t39_g16778);
        double cse_p13_t55_g13114 = (cse_p60_t39_g16778 * dx * dy);
        double cse_p13_t56_g13115 = (cse_p60_t39_g16778 * cse_p60_t72_g16811 * dy);
        double cse_p12_t20_g12992 = (cse_p13_t101_g13160 * cse_p60_t39_g16778);
        double cse_p12_t21_g12993 = (cse_p13_t102_g13161 * cse_p60_t39_g16778);
        double cse_p12_t22_g12994 = (cse_p13_t103_g13162 * cse_p60_t39_g16778);
        double cse_p12_t25_g12997 = (cse_p60_t39_g16778 * dz);
        double cse_p13_t42_g13101 = (cse_p12_t25_g12997 * cse_p60_t72_g16811);
        double cse_p13_t43_g13102 = (cse_p12_t25_g12997 * cse_p60_t73_g16812);
        double cse_p13_t44_g13103 = (cse_p12_t25_g12997 * dx);
        double cse_p13_t45_g13104 = (cse_p12_t25_g12997 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p60_t38_g16777 = (cse_p565_t2_invr52251 * cse_p565_t2_invr52251 * d2e_dr2);
        double cse_p13_t46_g13105 = (cse_p60_t38_g16777 * dx * dx);
        double cse_p13_t47_g13106 = (cse_p60_t38_g16777 * dx * dy);
        double cse_p13_t48_g13107 = (cse_p60_t38_g16777 * dx * dz);
        double cse_p13_t49_g13108 = (cse_p60_t38_g16777 * dy * dy);
        double cse_p13_t50_g13109 = (cse_p60_t38_g16777 * dy * dz);
        double cse_p13_t51_g13110 = (cse_p60_t38_g16777 * dz * dz);
        double cse_p12_t8_g12980 = (cse_p60_t38_g16777 * cse_p60_t72_g16811);
        double cse_p13_t57_g13116 = (cse_p12_t8_g12980 * cse_p60_t72_g16811);
        double cse_p13_t58_g13117 = (cse_p12_t8_g12980 * cse_p60_t73_g16812);
        double cse_p13_t59_g13118 = (cse_p12_t8_g12980 * cse_p60_t74_g16813);
        double cse_p13_t60_g13119 = (cse_p12_t8_g12980 * dx);
        double cse_p13_t61_g13120 = (cse_p12_t8_g12980 * dy);
        double cse_p13_t62_g13121 = (cse_p12_t8_g12980 * dz);
        double cse_p12_t9_g12981 = (cse_p60_t38_g16777 * cse_p60_t73_g16812);
        double cse_p13_t63_g13122 = (cse_p12_t9_g12981 * cse_p60_t73_g16812);
        double cse_p13_t64_g13123 = (cse_p12_t9_g12981 * cse_p60_t74_g16813);
        double cse_p13_t65_g13124 = (cse_p12_t9_g12981 * dx);
        double cse_p13_t66_g13125 = (cse_p12_t9_g12981 * dy);
        double cse_p13_t67_g13126 = (cse_p12_t9_g12981 * dz);
        double cse_p12_t10_g12982 = (cse_p60_t38_g16777 * cse_p60_t74_g16813);
        double cse_p13_t38_g13097 = (cse_p12_t10_g12982 * cse_p60_t74_g16813);
        double cse_p13_t39_g13098 = (cse_p12_t10_g12982 * dx);
        double cse_p13_t40_g13099 = (cse_p12_t10_g12982 * dy);
        double cse_p13_t41_g13100 = (cse_p12_t10_g12982 * dz);
        double cse_p12_t52_g13024 = (cse_p13_t47_g13106 + cse_p13_t56_g13115);
        double cse_p12_t53_g13025 = (cse_p13_t42_g13101 + cse_p13_t48_g13107);
        double cse_p12_t54_g13026 = (cse_p13_t43_g13102 + cse_p13_t50_g13109);
        double cse_p12_t70_g13042 = (cse_p12_t21_g12993 + cse_p13_t49_g13108);
        double cse_p12_t72_g13044 = (cse_p12_t22_g12994 + cse_p13_t51_g13110);
        double cse_p12_t74_g13046 = (cse_p12_t20_g12992 + cse_p13_t46_g13105);
        double cse_p12_t49_g13021 = (cse_p13_t52_g13111 + cse_p13_t60_g13119);
        double cse_p12_t55_g13027 = (cse_p13_t55_g13114 + cse_p13_t61_g13120);
        double cse_p12_t57_g13029 = (cse_p13_t56_g13115 + cse_p13_t58_g13117);
        double cse_p12_t58_g13030 = (cse_p13_t42_g13101 + cse_p13_t59_g13118);
        double cse_p12_t59_g13031 = (cse_p13_t44_g13103 + cse_p13_t62_g13121);
        double cse_p12_t75_g13047 = (cse_p12_t20_g12992 + cse_p13_t57_g13116);
        double cse_p12_t50_g13022 = (cse_p13_t53_g13112 + cse_p13_t66_g13125);
        double cse_p12_t56_g13028 = (cse_p13_t55_g13114 + cse_p13_t65_g13124);
        double cse_p12_t60_g13032 = (cse_p13_t43_g13102 + cse_p13_t64_g13123);
        double cse_p12_t61_g13033 = (cse_p13_t45_g13104 + cse_p13_t67_g13126);
        double cse_p12_t71_g13043 = (cse_p12_t21_g12993 + cse_p13_t63_g13122);
        double cse_p12_t51_g13023 = (cse_p13_t41_g13100 + cse_p13_t54_g13113);
        double cse_p12_t62_g13034 = (cse_p13_t39_g13098 + cse_p13_t44_g13103);
        double cse_p12_t63_g13035 = (cse_p13_t40_g13099 + cse_p13_t45_g13104);
        double cse_p12_t73_g13045 = (cse_p12_t22_g12994 + cse_p13_t38_g13097);
        double h_x1_x1 = cse_p12_t74_g13046;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t52_g13024;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t53_g13025;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t49_g13021;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t56_g13028;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t62_g13034;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t70_g13042;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t54_g13026;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t55_g13027;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t50_g13022;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t63_g13035;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t72_g13044;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t59_g13031;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t61_g13033;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t51_g13023;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t75_g13047;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t57_g13029;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t58_g13030;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t71_g13043;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t60_g13032;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t73_g13045;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= params.rlow)) {
          /* !COND-2 */
          double cse_p12_t7_g12979 = (term.bb_low * rr);
          double cse_p11_t1_g12964 = (rr * rr);
          double cse_p60_t56_g16795 = (cse_p11_t41_g12909 + (-(cse_p11_t1_g12964)));
          double cse_p60_t57_g16796 = (cse_p11_t42_g12910 + (-(cse_p11_t1_g12964)));
          double cse_p60_t58_g16797 = (cse_p11_t43_g12911 + (-(cse_p11_t1_g12964)));
          double cse_p60_t59_g16798 = (cse_p11_t1_g12964 + (-(cse_p11_t41_g12909)));
          double cse_p60_t60_g16799 = (cse_p11_t1_g12964 + (-(cse_p11_t42_g12910)));
          double cse_p60_t61_g16800 = (cse_p11_t1_g12964 + (-(cse_p11_t43_g12911)));
          double cse_p60_t71_g16810 = (-(cse_p11_t1_g12964));
          double cse_p12_t5_g12977 = (3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g12964);
          double cse_p11_t77_g12945 = (term.cc_low + cse_p12_t5_g12977 + cse_p12_t7_g12979 + cse_p12_t7_g12979);
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p11_t1_g12964 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (cse_p11_t77_g12945 * cse_p1_t2_invsqrt12860 * dx * params.elec_weight);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t77_g12945 * cse_p1_t2_invsqrt12860 * dy * params.elec_weight);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t77_g12945 * cse_p1_t2_invsqrt12860 * dz * params.elec_weight);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (params.elec_weight * ((cse_p1_t2_invsqrt12860 * dx * ((-2.00000000000000000e+0 * term.bb_low * rr) + (-3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g12964))) + (-((term.cc_low * cse_p1_t2_invsqrt12860 * dx)))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (params.elec_weight * ((cse_p1_t2_invsqrt12860 * dy * ((-2.00000000000000000e+0 * term.bb_low * rr) + (-3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g12964))) + (-((term.cc_low * cse_p1_t2_invsqrt12860 * dy)))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (params.elec_weight * ((cse_p1_t2_invsqrt12860 * dz * ((-2.00000000000000000e+0 * term.bb_low * rr) + (-3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g12964))) + (-((term.cc_low * cse_p1_t2_invsqrt12860 * dz)))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (params.elec_weight * (term.cc_low + cse_p12_t5_g12977 + (2.00000000000000000e+0 * term.bb_low * rr)));
          double d2e_dr2 = (params.elec_weight * ((2.00000000000000000e+0 * term.bb_low) + (6.00000000000000000e+0 * term.aa_low * rr)));
          double cse_p565_t2_invr52251 = (1.0 / (rr));
          double cse_p11_t1_g12962 = (cse_p565_t2_invr52251 * cse_p565_t2_invr52251 * d2e_dr2);
          double cse_p60_t1_g16740 = (cse_p11_t1_g12962 * (-(dx)));
          double cse_p60_t2_g16741 = (cse_p11_t1_g12962 * (-(dy)));
          double cse_p60_t3_g16742 = (cse_p11_t1_g12962 * (-(dz)));
          double cse_p60_t4_g16743 = (cse_p11_t1_g12962 * dx * dx);
          double cse_p60_t5_g16744 = (cse_p11_t1_g12962 * dx * dy);
          double cse_p60_t6_g16745 = (cse_p11_t1_g12962 * dx * dz);
          double cse_p60_t7_g16746 = (cse_p11_t1_g12962 * dy * dy);
          double cse_p60_t8_g16747 = (cse_p11_t1_g12962 * dy * dz);
          double cse_p60_t9_g16748 = (cse_p11_t1_g12962 * dz * dz);
          double cse_p14_t42_g13215 = (cse_p60_t1_g16740 * (-(dx)));
          double cse_p14_t43_g13216 = (cse_p60_t1_g16740 * (-(dy)));
          double cse_p14_t44_g13217 = (cse_p60_t1_g16740 * (-(dz)));
          double cse_p14_t45_g13218 = (cse_p60_t1_g16740 * dx);
          double cse_p14_t46_g13219 = (cse_p60_t1_g16740 * dy);
          double cse_p14_t47_g13220 = (cse_p60_t1_g16740 * dz);
          double cse_p14_t48_g13221 = (cse_p60_t2_g16741 * (-(dy)));
          double cse_p14_t49_g13222 = (cse_p60_t2_g16741 * (-(dz)));
          double cse_p14_t50_g13223 = (cse_p60_t2_g16741 * dx);
          double cse_p14_t51_g13224 = (cse_p60_t2_g16741 * dy);
          double cse_p14_t52_g13225 = (cse_p60_t2_g16741 * dz);
          double cse_p14_t53_g13226 = (cse_p60_t3_g16742 * (-(dz)));
          double cse_p14_t54_g13227 = (cse_p60_t3_g16742 * dx);
          double cse_p14_t55_g13228 = (cse_p60_t3_g16742 * dy);
          double cse_p14_t56_g13229 = (cse_p60_t3_g16742 * dz);
          double cse_p565_t1_invr52250 = (1.0 / ((rr * rr * rr)));
          double cse_p11_t1_g12963 = (cse_p565_t1_invr52250 * de_dr);
          double cse_p60_t10_g16749 = (cse_p11_t1_g12963 * (cse_p11_t41_g12909 + cse_p60_t71_g16810));
          double cse_p60_t11_g16750 = (cse_p11_t1_g12963 * (cse_p11_t42_g12910 + cse_p60_t71_g16810));
          double cse_p60_t12_g16751 = (cse_p11_t1_g12963 * (cse_p11_t43_g12911 + cse_p60_t71_g16810));
          double cse_p60_t13_g16752 = (cse_p11_t1_g12963 * cse_p60_t59_g16798);
          double cse_p60_t14_g16753 = (cse_p11_t1_g12963 * cse_p60_t60_g16799);
          double cse_p60_t15_g16754 = (cse_p11_t1_g12963 * cse_p60_t61_g16800);
          double cse_p60_t16_g16755 = (cse_p11_t1_g12963 * dx * dy);
          double cse_p60_t17_g16756 = (cse_p11_t1_g12963 * dy * (-(dx)));
          double cse_p60_t18_g16757 = (cse_p11_t1_g12963 * dz);
          double cse_p60_t44_g16783 = (cse_p60_t17_g16756 + cse_p60_t5_g16744);
          double cse_p14_t38_g13211 = (cse_p60_t18_g16757 * (-(dx)));
          double cse_p14_t39_g13212 = (cse_p60_t18_g16757 * (-(dy)));
          double cse_p14_t40_g13213 = (cse_p60_t18_g16757 * dx);
          double cse_p14_t41_g13214 = (cse_p60_t18_g16757 * dy);
          double cse_p13_t99_g13158 = (cse_p60_t13_g16752 + cse_p60_t4_g16743);
          double h_x1_x1 = cse_p13_t99_g13158;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t44_g16783;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t75_g13134 = (cse_p14_t38_g13211 + cse_p60_t6_g16745);
          double cse_p13_t76_g13135 = (cse_p14_t39_g13212 + cse_p60_t8_g16747);
          double h_x1_z1 = cse_p13_t75_g13134;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t71_g13130 = (cse_p14_t45_g13218 + cse_p60_t10_g16749);
          double cse_p13_t77_g13136 = (cse_p14_t46_g13219 + cse_p60_t16_g16755);
          double cse_p13_t79_g13138 = (cse_p14_t43_g13216 + cse_p60_t17_g16756);
          double cse_p13_t80_g13139 = (cse_p14_t38_g13211 + cse_p14_t44_g13217);
          double cse_p13_t81_g13140 = (cse_p14_t40_g13213 + cse_p14_t47_g13220);
          double cse_p13_t100_g13159 = (cse_p14_t42_g13215 + cse_p60_t13_g16752);
          double h_x1_x2 = cse_p13_t71_g13130;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t72_g13131 = (cse_p14_t51_g13224 + cse_p60_t11_g16750);
          double cse_p13_t78_g13137 = (cse_p14_t50_g13223 + cse_p60_t16_g16755);
          double cse_p13_t82_g13141 = (cse_p14_t39_g13212 + cse_p14_t49_g13222);
          double cse_p13_t83_g13142 = (cse_p14_t41_g13214 + cse_p14_t52_g13225);
          double h_x1_y2 = cse_p13_t78_g13137;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t73_g13132 = (cse_p14_t56_g13229 + cse_p60_t12_g16751);
          double cse_p13_t84_g13143 = (cse_p14_t40_g13213 + cse_p14_t54_g13227);
          double cse_p13_t85_g13144 = (cse_p14_t41_g13214 + cse_p14_t55_g13228);
          double h_x1_z2 = cse_p13_t84_g13143;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t95_g13154 = (cse_p60_t14_g16753 + cse_p60_t7_g16746);
          double cse_p13_t96_g13155 = (cse_p14_t48_g13221 + cse_p60_t14_g16753);
          double h_y1_y1 = cse_p13_t95_g13154;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t76_g13135;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t77_g13136;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t72_g13131;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t85_g13144;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t97_g13156 = (cse_p60_t15_g16754 + cse_p60_t9_g16748);
          double cse_p13_t98_g13157 = (cse_p14_t53_g13226 + cse_p60_t15_g16754);
          double h_z1_z1 = cse_p13_t97_g13156;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t81_g13140;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t83_g13142;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t73_g13132;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t100_g13159;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t79_g13138;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t80_g13139;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t96_g13155;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t82_g13141;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t98_g13157;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t83_g13055 = (-(params.eps_solvent));
            double cse_p11_t78_g12946 = (cse_p12_t83_g13055 + params.eps_core);
            double cse_p11_t87_g12955 = (-(term.kqq));
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p12_t1_g12973 = (-1.25000000000000000e-1 * rx);
            double cse_p12_t2_g12974 = (1.25000000000000000e-1 * rx);
            double cse_p12_t3_g12975 = (5.00000000000000000e-1 * rx * rx);
            double cse_p12_t84_g13056 = (-(rx));
            double cse_p11_t73_g12941 = (-2.50000000000000000e-1 + cse_p12_t1_g12973 + cse_p12_t1_g12973);
            double cse_p11_t74_g12942 = (2.50000000000000000e-1 + cse_p12_t2_g12974 + cse_p12_t2_g12974);
            double cse_p11_t75_g12943 = (1.00000000000000000e+0 + cse_p12_t3_g12975 + rx);
            double cse_p11_t89_g12957 = exp(cse_p12_t84_g13056);
            double ex = cse_p11_t89_g12957;
            double gx = (cse_p11_t75_g12943 * ex);
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p12_t3_g12975 + cse_p12_t84_g13056));
            double ep = (params.eps_solvent + (cse_p11_t78_g12946 * gx));
            double cse_p12_t47_g13019 = (ep * ep);
            double cse_p11_t90_g12958 = (1.0 / (cse_p12_t47_g13019));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t78_g12946 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t78_g12946 * gx_pp);
            double cse_p1_t3_invr12861 = (1.0 / (ep));
            double cse_p1_t4_invr212862 = (cse_p1_t3_invr12861 * cse_p1_t3_invr12861);
            double ee = (cse_p1_t3_invr12861 * term.kqq);
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p11_t1_g12972 = (cse_p11_t75_g12943 * cse_p11_t89_g12957 * cse_p1_t2_invsqrt12860);
            double g_x1 = (-((cse_p11_t78_g12946 * cse_p1_t4_invr212862 * params.elec_weight * term.kqq * ((-2.50000000000000000e-1 * cse_p11_t1_g12972 * dx) + (cse_p11_t74_g12942 * cse_p1_t2_invsqrt12860 * dx * ex)))));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double g_y1 = (-((cse_p11_t78_g12946 * cse_p1_t4_invr212862 * params.elec_weight * term.kqq * ((-2.50000000000000000e-1 * cse_p11_t1_g12972 * dy) + (cse_p11_t74_g12942 * cse_p1_t2_invsqrt12860 * dy * ex)))));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double g_z1 = (-((cse_p11_t78_g12946 * cse_p1_t4_invr212862 * params.elec_weight * term.kqq * ((-2.50000000000000000e-1 * cse_p11_t1_g12972 * dz) + (cse_p11_t74_g12942 * cse_p1_t2_invsqrt12860 * dz * ex)))));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (-((cse_p11_t78_g12946 * cse_p1_t4_invr212862 * params.elec_weight * term.kqq * ((2.50000000000000000e-1 * cse_p11_t1_g12972 * dx) + (cse_p11_t73_g12941 * cse_p1_t2_invsqrt12860 * dx * ex)))));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (-((cse_p11_t78_g12946 * cse_p1_t4_invr212862 * params.elec_weight * term.kqq * ((2.50000000000000000e-1 * cse_p11_t1_g12972 * dy) + (cse_p11_t73_g12941 * cse_p1_t2_invsqrt12860 * dy * ex)))));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (-((cse_p11_t78_g12946 * cse_p1_t4_invr212862 * params.elec_weight * term.kqq * ((2.50000000000000000e-1 * cse_p11_t1_g12972 * dz) + (cse_p11_t73_g12941 * cse_p1_t2_invsqrt12860 * dz * ex)))));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (cse_p11_t87_g12955 * cse_p11_t90_g12958 * params.elec_weight * ep_p);
            double d2e_dr2 = (cse_p11_t87_g12955 * params.elec_weight * ((-2.00000000000000000e+0 * ep_p * ep_p * (1.0 / ((ep * ep * ep)))) + (cse_p11_t90_g12958 * ep_pp)));
            double cse_p565_t2_invr52251 = (1.0 / (rr));
            double cse_p11_t1_g12962 = (cse_p565_t2_invr52251 * cse_p565_t2_invr52251 * d2e_dr2);
            double cse_p565_t1_invr52250 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g12963 = (cse_p565_t1_invr52250 * de_dr);
            double cse_p11_t1_g12964 = (rr * rr);
            double cse_p11_t1_g12971 = (cse_p11_t1_g12963 * (cse_p11_t1_g12964 + (-(cse_p11_t41_g12909))));
            double cse_p60_t66_g16805 = (cse_p11_t1_g12971 + (cse_p11_t1_g12962 * dx * dx));
            double h_x1_x1 = cse_p60_t66_g16805;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = ((cse_p11_t1_g12962 * dx * dy) + (cse_p11_t1_g12963 * dy * (-(dx))));
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p11_t1_g12967 = (cse_p11_t1_g12963 * dz);
            double cse_p60_t30_g16769 = (cse_p11_t1_g12967 * (-(dx)));
            double cse_p60_t31_g16770 = (cse_p11_t1_g12967 * (-(dy)));
            double cse_p60_t32_g16771 = (cse_p11_t1_g12967 * dx);
            double cse_p60_t33_g16772 = (cse_p11_t1_g12967 * dy);
            double cse_p60_t45_g16784 = (cse_p60_t30_g16769 + (cse_p11_t1_g12962 * dx * dz));
            double cse_p60_t46_g16785 = (cse_p60_t31_g16770 + (cse_p11_t1_g12962 * dy * dz));
            double h_x1_z1 = cse_p60_t45_g16784;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p11_t1_g12965 = (cse_p11_t1_g12962 * (-(dx)));
            double cse_p60_t19_g16758 = (cse_p11_t1_g12965 * (-(dx)));
            double cse_p60_t20_g16759 = (cse_p11_t1_g12965 * (-(dy)));
            double cse_p60_t21_g16760 = (cse_p11_t1_g12965 * (-(dz)));
            double cse_p60_t22_g16761 = (cse_p11_t1_g12965 * dx);
            double cse_p60_t23_g16762 = (cse_p11_t1_g12965 * dy);
            double cse_p60_t24_g16763 = (cse_p11_t1_g12965 * dz);
            double cse_p60_t41_g16780 = (cse_p60_t22_g16761 + (cse_p11_t1_g12963 * (cse_p11_t41_g12909 + (-(cse_p11_t1_g12964)))));
            double cse_p60_t47_g16786 = (cse_p60_t23_g16762 + (cse_p11_t1_g12963 * dx * dy));
            double cse_p60_t49_g16788 = (cse_p60_t20_g16759 + (cse_p11_t1_g12963 * dy * (-(dx))));
            double cse_p60_t50_g16789 = (cse_p60_t21_g16760 + cse_p60_t30_g16769);
            double cse_p60_t51_g16790 = (cse_p60_t24_g16763 + cse_p60_t32_g16771);
            double cse_p60_t67_g16806 = (cse_p11_t1_g12971 + cse_p60_t19_g16758);
            double h_x1_x2 = cse_p60_t41_g16780;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p11_t1_g12966 = (cse_p11_t1_g12962 * (-(dy)));
            double cse_p60_t25_g16764 = (cse_p11_t1_g12966 * (-(dy)));
            double cse_p60_t26_g16765 = (cse_p11_t1_g12966 * (-(dz)));
            double cse_p60_t27_g16766 = (cse_p11_t1_g12966 * dx);
            double cse_p60_t28_g16767 = (cse_p11_t1_g12966 * dy);
            double cse_p60_t29_g16768 = (cse_p11_t1_g12966 * dz);
            double cse_p60_t42_g16781 = (cse_p60_t28_g16767 + (cse_p11_t1_g12963 * (cse_p11_t42_g12910 + (-(cse_p11_t1_g12964)))));
            double cse_p60_t48_g16787 = (cse_p60_t27_g16766 + (cse_p11_t1_g12963 * dx * dy));
            double cse_p60_t52_g16791 = (cse_p60_t26_g16765 + cse_p60_t31_g16770);
            double cse_p60_t53_g16792 = (cse_p60_t29_g16768 + cse_p60_t33_g16772);
            double h_x1_y2 = cse_p60_t48_g16787;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p11_t1_g12968 = (cse_p11_t1_g12962 * (-(dz)));
            double cse_p60_t34_g16773 = (cse_p11_t1_g12968 * (-(dz)));
            double cse_p60_t35_g16774 = (cse_p11_t1_g12968 * dx);
            double cse_p60_t36_g16775 = (cse_p11_t1_g12968 * dy);
            double cse_p60_t37_g16776 = (cse_p11_t1_g12968 * dz);
            double cse_p60_t43_g16782 = (cse_p60_t37_g16776 + (cse_p11_t1_g12963 * (cse_p11_t43_g12911 + (-(cse_p11_t1_g12964)))));
            double cse_p60_t54_g16793 = (cse_p60_t32_g16771 + cse_p60_t35_g16774);
            double cse_p60_t55_g16794 = (cse_p60_t33_g16772 + cse_p60_t36_g16775);
            double h_x1_z2 = cse_p60_t54_g16793;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p11_t1_g12969 = (cse_p11_t1_g12963 * (cse_p11_t1_g12964 + (-(cse_p11_t42_g12910))));
            double cse_p60_t62_g16801 = (cse_p11_t1_g12969 + (cse_p11_t1_g12962 * dy * dy));
            double cse_p60_t63_g16802 = (cse_p11_t1_g12969 + cse_p60_t25_g16764);
            double h_y1_y1 = cse_p60_t62_g16801;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t46_g16785;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t47_g16786;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t42_g16781;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t55_g16794;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p11_t1_g12970 = (cse_p11_t1_g12963 * (cse_p11_t1_g12964 + (-(cse_p11_t43_g12911))));
            double cse_p60_t64_g16803 = (cse_p11_t1_g12970 + (cse_p11_t1_g12962 * dz * dz));
            double cse_p60_t65_g16804 = (cse_p11_t1_g12970 + cse_p60_t34_g16773);
            double h_z1_z1 = cse_p60_t64_g16803;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t51_g16790;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t53_g16792;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t43_g16782;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t67_g16806;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t49_g16788;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t50_g16789;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t63_g16802;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t52_g16791;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t65_g16804;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p12_t6_g12978 = (term.bb_high * rr);
              double cse_p11_t1_g12964 = (rr * rr);
              double cse_p12_t4_g12976 = (3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g12964);
              double cse_p11_t76_g12944 = (term.cc_high + cse_p12_t4_g12976 + cse_p12_t6_g12978 + cse_p12_t6_g12978);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p11_t1_g12964 * (term.bb_high + (term.aa_high * rr)))));
              energy_added += energy;
              *energy_accumulate += energy;
              double g_x1 = (cse_p11_t76_g12944 * cse_p1_t2_invsqrt12860 * dx * params.elec_weight);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t76_g12944 * cse_p1_t2_invsqrt12860 * dy * params.elec_weight);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t76_g12944 * cse_p1_t2_invsqrt12860 * dz * params.elec_weight);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = (params.elec_weight * ((cse_p1_t2_invsqrt12860 * dx * ((-2.00000000000000000e+0 * term.bb_high * rr) + (-3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g12964))) + (-((term.cc_high * cse_p1_t2_invsqrt12860 * dx)))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (params.elec_weight * ((cse_p1_t2_invsqrt12860 * dy * ((-2.00000000000000000e+0 * term.bb_high * rr) + (-3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g12964))) + (-((term.cc_high * cse_p1_t2_invsqrt12860 * dy)))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (params.elec_weight * ((cse_p1_t2_invsqrt12860 * dz * ((-2.00000000000000000e+0 * term.bb_high * rr) + (-3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g12964))) + (-((term.cc_high * cse_p1_t2_invsqrt12860 * dz)))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double de_dr = (params.elec_weight * (term.cc_high + cse_p12_t4_g12976 + (2.00000000000000000e+0 * term.bb_high * rr)));
              double d2e_dr2 = (params.elec_weight * ((2.00000000000000000e+0 * term.bb_high) + (6.00000000000000000e+0 * term.aa_high * rr)));
              double cse_p565_t2_invr52251 = (1.0 / (rr));
              double cse_p11_t1_g12962 = (cse_p565_t2_invr52251 * cse_p565_t2_invr52251 * d2e_dr2);
              double cse_p565_t1_invr52250 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g12963 = (cse_p565_t1_invr52250 * de_dr);
              double cse_p11_t1_g12971 = (cse_p11_t1_g12963 * (cse_p11_t1_g12964 + (-(cse_p11_t41_g12909))));
              double h_x1_x1 = (cse_p11_t1_g12971 + (cse_p11_t1_g12962 * dx * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double h_x1_y1 = ((cse_p11_t1_g12962 * dx * dy) + (cse_p11_t1_g12963 * dy * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g12967 = (cse_p11_t1_g12963 * dz);
              double h_x1_z1 = ((cse_p11_t1_g12962 * dx * dz) + (cse_p11_t1_g12967 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g12965 = (cse_p11_t1_g12962 * (-(dx)));
              double h_x1_x2 = ((cse_p11_t1_g12963 * (cse_p11_t41_g12909 + (-(cse_p11_t1_g12964)))) + (cse_p11_t1_g12965 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g12966 = (cse_p11_t1_g12962 * (-(dy)));
              double h_x1_y2 = ((cse_p11_t1_g12963 * dx * dy) + (cse_p11_t1_g12966 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g12968 = (cse_p11_t1_g12962 * (-(dz)));
              double h_x1_z2 = ((cse_p11_t1_g12967 * dx) + (cse_p11_t1_g12968 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g12969 = (cse_p11_t1_g12963 * (cse_p11_t1_g12964 + (-(cse_p11_t42_g12910))));
              double h_y1_y1 = (cse_p11_t1_g12969 + (cse_p11_t1_g12962 * dy * dy));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = ((cse_p11_t1_g12962 * dy * dz) + (cse_p11_t1_g12967 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = ((cse_p11_t1_g12963 * dx * dy) + (cse_p11_t1_g12965 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = ((cse_p11_t1_g12963 * (cse_p11_t42_g12910 + (-(cse_p11_t1_g12964)))) + (cse_p11_t1_g12966 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = ((cse_p11_t1_g12967 * dy) + (cse_p11_t1_g12968 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g12970 = (cse_p11_t1_g12963 * (cse_p11_t1_g12964 + (-(cse_p11_t43_g12911))));
              double h_z1_z1 = (cse_p11_t1_g12970 + (cse_p11_t1_g12962 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = ((cse_p11_t1_g12965 * dz) + (cse_p11_t1_g12967 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = ((cse_p11_t1_g12966 * dz) + (cse_p11_t1_g12967 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = ((cse_p11_t1_g12963 * (cse_p11_t43_g12911 + (-(cse_p11_t1_g12964)))) + (cse_p11_t1_g12968 * dz));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p11_t1_g12971 + (cse_p11_t1_g12965 * (-(dx))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = ((cse_p11_t1_g12963 * dy * (-(dx))) + (cse_p11_t1_g12965 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = ((cse_p11_t1_g12965 * (-(dz))) + (cse_p11_t1_g12967 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = (cse_p11_t1_g12969 + (cse_p11_t1_g12966 * (-(dy))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = ((cse_p11_t1_g12966 * (-(dz))) + (cse_p11_t1_g12967 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = (cse_p11_t1_g12970 + (cse_p11_t1_g12968 * (-(dz))));
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
                double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t41_g12909))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t2_invr52251 = (1.0 / (rr));
                double cse_p565_t1_invr52250 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g12964 = (rr * rr);
                double cse_p11_t1_g12963 = (cse_p565_t1_invr52250 * de_dr);
                double cse_p11_t1_g12962 = (cse_p565_t2_invr52251 * cse_p565_t2_invr52251 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g12962 * dx * dy) + (cse_p11_t1_g12963 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double h_x1_z1 = ((cse_p11_t1_g12962 * dx * dz) + (cse_p11_t1_g12963 * dz * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g12967 = (cse_p11_t1_g12963 * dz);
                double h_x1_x2 = ((cse_p11_t1_g12962 * dx * (-(dx))) + (cse_p11_t1_g12963 * (cse_p11_t41_g12909 + (-(cse_p11_t1_g12964)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g12965 = (cse_p11_t1_g12962 * (-(dx)));
                double h_x1_y2 = ((cse_p11_t1_g12962 * dx * (-(dy))) + (cse_p11_t1_g12963 * dx * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g12966 = (cse_p11_t1_g12962 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g12962 * dx * (-(dz))) + (cse_p11_t1_g12967 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g12968 = (cse_p11_t1_g12962 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g12962 * dy * dy) + (cse_p11_t1_g12963 * (cse_p11_t1_g12964 + (-(cse_p11_t42_g12910)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g12969 = (cse_p11_t1_g12963 * (cse_p11_t1_g12964 + (-(cse_p11_t42_g12910))));
                double h_y1_z1 = ((cse_p11_t1_g12962 * dy * dz) + (cse_p11_t1_g12967 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g12963 * dx * dy) + (cse_p11_t1_g12965 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g12963 * (cse_p11_t42_g12910 + (-(cse_p11_t1_g12964)))) + (cse_p11_t1_g12966 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g12967 * dy) + (cse_p11_t1_g12968 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g12962 * dz * dz) + (cse_p11_t1_g12963 * (cse_p11_t1_g12964 + (-(cse_p11_t43_g12911)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g12970 = (cse_p11_t1_g12963 * (cse_p11_t1_g12964 + (-(cse_p11_t43_g12911))));
                double h_z1_x2 = ((cse_p11_t1_g12965 * dz) + (cse_p11_t1_g12967 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g12966 * dz) + (cse_p11_t1_g12967 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g12963 * (cse_p11_t43_g12911 + (-(cse_p11_t1_g12964)))) + (cse_p11_t1_g12968 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g12963 * (cse_p11_t1_g12964 + (-(cse_p11_t41_g12909)))) + (cse_p11_t1_g12965 * (-(dx))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p11_t1_g12971 = (cse_p11_t1_g12963 * (cse_p11_t1_g12964 + (-(cse_p11_t41_g12909))));
                double h_x2_y2 = ((cse_p11_t1_g12963 * dy * (-(dx))) + (cse_p11_t1_g12965 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double h_x2_z2 = ((cse_p11_t1_g12965 * (-(dz))) + (cse_p11_t1_g12967 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g12969 + (cse_p11_t1_g12966 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g12966 * (-(dz))) + (cse_p11_t1_g12967 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g12970 + (cse_p11_t1_g12968 * (-(dz))));
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
