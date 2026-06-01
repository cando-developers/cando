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
      double cse_p1_t1_sqrt1 = sqrt(r2);
      double cse_p1_t2_invsqrt2 = (1.0 / (cse_p1_t1_sqrt1));
      double rr = cse_p1_t1_sqrt1;
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
          double cse_p528_t1_g7032 = (cse_p1_t2_invsqrt2 * dx);
          double cse_p528_t2_g7020 = (cse_p1_t2_invsqrt2 * dy);
          double cse_p528_t3_g7021 = (cse_p1_t2_invsqrt2 * dz);
          double cse_p528_t4_g7022 = (rr * rr);
          double cse_p12_t1_g42 = (term.aa_low * cse_p528_t4_g7022);
          double cse_p12_t4_g34 = (3.00000000000000000e+0 * cse_p12_t1_g42);
          double cse_p11_t1_g26 = (term.bb_low * rr);
          double cse_p11_t12_g16 = (term.cc_low + cse_p11_t1_g26 + cse_p11_t1_g26 + cse_p12_t4_g34);
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p528_t4_g7022 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p12_t1_g46 = (cse_p11_t12_g16 * params.elec_weight);
          double g_x1 = (cse_p12_t1_g46 * cse_p528_t1_g7032);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p12_t1_g46 * cse_p528_t2_g7020);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p12_t1_g46 * cse_p528_t3_g7021);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double cse_p63_t1_g731 = (-2.00000000000000000e+0 * cse_p11_t1_g26);
          double cse_p63_t1_g733 = (-3.00000000000000000e+0 * cse_p12_t1_g42);
          double g_x2 = (params.elec_weight * ((cse_p528_t1_g7032 * cse_p63_t1_g731) + (cse_p528_t1_g7032 * cse_p63_t1_g733) + (-((term.cc_low * cse_p528_t1_g7032)))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (params.elec_weight * ((cse_p528_t2_g7020 * cse_p63_t1_g731) + (cse_p528_t2_g7020 * cse_p63_t1_g733) + (-((term.cc_low * cse_p528_t2_g7020)))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (params.elec_weight * ((cse_p528_t3_g7021 * cse_p63_t1_g731) + (cse_p528_t3_g7021 * cse_p63_t1_g733) + (-((term.cc_low * cse_p528_t3_g7021)))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (params.elec_weight * (term.cc_low + cse_p12_t4_g34 + (2.00000000000000000e+0 * cse_p11_t1_g26)));
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t10_g40 = (-(params.eps_solvent));
            double cse_p11_t13_g17 = (cse_p12_t10_g40 + params.eps_core);
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p12_t1_g31 = (-1.25000000000000000e-1 * rx);
            double cse_p12_t2_g32 = (1.25000000000000000e-1 * rx);
            double cse_p12_t3_g33 = (5.00000000000000000e-1 * rx * rx);
            double cse_p12_t11_g41 = (-(rx));
            double cse_p11_t8_g12 = (-2.50000000000000000e-1 + cse_p12_t1_g31 + cse_p12_t1_g31);
            double cse_p11_t9_g13 = (2.50000000000000000e-1 + cse_p12_t2_g32 + cse_p12_t2_g32);
            double cse_p11_t10_g14 = (1.00000000000000000e+0 + cse_p12_t3_g33 + rx);
            double cse_p11_t16_g20 = exp(cse_p12_t11_g41);
            double ex = cse_p11_t16_g20;
            double gx = (cse_p11_t10_g14 * ex);
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p12_t11_g41 + cse_p12_t3_g33));
            double ep = (params.eps_solvent + (cse_p11_t13_g17 * gx));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t13_g17 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t13_g17 * gx_pp);
            double cse_p1_t3_invr3 = (1.0 / (ep));
            double cse_p1_t4_invr24 = (cse_p1_t3_invr3 * cse_p1_t3_invr3);
            double ee = (cse_p1_t3_invr3 * term.kqq);
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p528_t1_g7032 = (cse_p1_t2_invsqrt2 * dx);
            double cse_p11_t1_g27 = (cse_p11_t10_g14 * cse_p11_t16_g20);
            double cse_p11_t1_g28 = (cse_p11_t13_g17 * cse_p1_t4_invr24);
            double cse_p11_t1_g29 = (params.elec_weight * term.kqq);
            double cse_p12_t1_g43 = (-2.50000000000000000e-1 * cse_p11_t1_g27);
            double cse_p12_t1_g48 = (cse_p11_t9_g13 * ex);
            double cse_p63_t1_g729 = (cse_p11_t1_g28 * cse_p11_t1_g29);
            double g_x1 = (-((cse_p63_t1_g729 * ((cse_p12_t1_g43 * cse_p528_t1_g7032) + (cse_p12_t1_g48 * cse_p528_t1_g7032)))));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p11_t1_g22 = (cse_p1_t2_invsqrt2 * dy);
            double g_y1 = (-((cse_p63_t1_g729 * ((cse_p11_t1_g22 * cse_p12_t1_g43) + (cse_p11_t1_g22 * cse_p12_t1_g48)))));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p11_t1_g23 = (cse_p1_t2_invsqrt2 * dz);
            double g_z1 = (-((cse_p63_t1_g729 * ((cse_p11_t1_g23 * cse_p12_t1_g43) + (cse_p11_t1_g23 * cse_p12_t1_g48)))));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double cse_p12_t1_g44 = (2.50000000000000000e-1 * cse_p11_t1_g27);
            double cse_p12_t1_g47 = (cse_p11_t8_g12 * ex);
            double g_x2 = (-((cse_p63_t1_g729 * ((cse_p12_t1_g44 * cse_p528_t1_g7032) + (cse_p12_t1_g47 * cse_p528_t1_g7032)))));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (-((cse_p63_t1_g729 * ((cse_p11_t1_g22 * cse_p12_t1_g44) + (cse_p11_t1_g22 * cse_p12_t1_g47)))));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (-((cse_p63_t1_g729 * ((cse_p11_t1_g23 * cse_p12_t1_g44) + (cse_p11_t1_g23 * cse_p12_t1_g47)))));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (params.elec_weight * ep_p * (1.0 / ((ep * ep))) * (-(term.kqq)));
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p11_t1_g24 = (rr * rr);
              double cse_p11_t1_g30 = (term.aa_high * cse_p11_t1_g24);
              double cse_p12_t5_g35 = (3.00000000000000000e+0 * cse_p11_t1_g30);
              double cse_p11_t1_g25 = (term.bb_high * rr);
              double cse_p11_t11_g15 = (term.cc_high + cse_p11_t1_g25 + cse_p11_t1_g25 + cse_p12_t5_g35);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p11_t1_g24 * (term.bb_high + (term.aa_high * rr)))));
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p528_t1_g7032 = (cse_p1_t2_invsqrt2 * dx);
              double cse_p12_t1_g45 = (cse_p11_t11_g15 * params.elec_weight);
              double g_x1 = (cse_p12_t1_g45 * cse_p528_t1_g7032);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double cse_p11_t1_g22 = (cse_p1_t2_invsqrt2 * dy);
              double g_y1 = (cse_p11_t1_g22 * cse_p12_t1_g45);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double cse_p11_t1_g23 = (cse_p1_t2_invsqrt2 * dz);
              double g_z1 = (cse_p11_t1_g23 * cse_p12_t1_g45);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p63_t1_g730 = (-2.00000000000000000e+0 * cse_p11_t1_g25);
              double cse_p63_t1_g732 = (-3.00000000000000000e+0 * cse_p11_t1_g30);
              double g_x2 = (params.elec_weight * ((cse_p528_t1_g7032 * cse_p63_t1_g730) + (cse_p528_t1_g7032 * cse_p63_t1_g732) + (-((term.cc_high * cse_p528_t1_g7032)))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (params.elec_weight * ((cse_p11_t1_g22 * cse_p63_t1_g730) + (cse_p11_t1_g22 * cse_p63_t1_g732) + (-((term.cc_high * cse_p11_t1_g22)))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (params.elec_weight * ((cse_p11_t1_g23 * cse_p63_t1_g730) + (cse_p11_t1_g23 * cse_p63_t1_g732) + (-((term.cc_high * cse_p11_t1_g23)))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double de_dr = (params.elec_weight * (term.cc_high + cse_p12_t5_g35 + (2.00000000000000000e+0 * cse_p11_t1_g25)));
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
      double cse_p11_t41_g45 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t42_g46 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t43_g47 = (dz * dz);
      double r2 = (cse_p11_t41_g45 + cse_p11_t42_g46 + cse_p11_t43_g47);
      double cse_p1_t1_sqrt1 = sqrt(r2);
      double cse_p1_t2_invsqrt2 = (1.0 / (cse_p1_t1_sqrt1));
      double rr = cse_p1_t1_sqrt1;
      if ((rr <= params.rmin)) {
        /* !COND-1 */
        double cse_p60_t40_g3915 = (rr * rr);
        double cse_p60_t68_g3943 = (-(cse_p11_t41_g45));
        double cse_p60_t69_g3944 = (-(cse_p11_t42_g46));
        double cse_p60_t70_g3945 = (-(cse_p11_t43_g47));
        double cse_p60_t72_g3947 = (-(dx));
        double cse_p60_t73_g3948 = (-(dy));
        double cse_p60_t74_g3949 = (-(dz));
        double cse_p565_t1_invr39386 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr39387 = (1.0 / (rr));
        double cse_p60_t77_g3952 = (rr * rr * rr);
        double cse_p15_t72_g481 = (-(cse_p60_t40_g3915));
        double cse_p14_t76_g385 = (cse_p11_t41_g45 + cse_p15_t72_g481);
        double cse_p14_t78_g387 = (cse_p11_t42_g46 + cse_p15_t72_g481);
        double cse_p14_t80_g389 = (cse_p11_t43_g47 + cse_p15_t72_g481);
        double cse_p13_t101_g296 = (cse_p60_t40_g3915 + cse_p60_t68_g3943);
        double cse_p13_t102_g297 = (cse_p60_t40_g3915 + cse_p60_t69_g3944);
        double cse_p13_t103_g298 = (cse_p60_t40_g3915 + cse_p60_t70_g3945);
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
        double cse_p60_t39_g3914 = (cse_p565_t1_invr39386 * de_dr);
        double cse_p13_t63_g258 = (cse_p14_t76_g385 * cse_p60_t39_g3914);
        double cse_p13_t64_g259 = (cse_p14_t78_g387 * cse_p60_t39_g3914);
        double cse_p13_t65_g260 = (cse_p14_t80_g389 * cse_p60_t39_g3914);
        double cse_p13_t66_g261 = (cse_p60_t39_g3914 * dx * dy);
        double cse_p13_t67_g262 = (cse_p60_t39_g3914 * cse_p60_t72_g3947 * dy);
        double cse_p12_t39_g147 = (cse_p13_t101_g296 * cse_p60_t39_g3914);
        double cse_p12_t40_g148 = (cse_p13_t102_g297 * cse_p60_t39_g3914);
        double cse_p12_t41_g149 = (cse_p13_t103_g298 * cse_p60_t39_g3914);
        double cse_p12_t44_g152 = (cse_p60_t39_g3914 * dz);
        double cse_p13_t53_g248 = (cse_p12_t44_g152 * cse_p60_t72_g3947);
        double cse_p13_t54_g249 = (cse_p12_t44_g152 * cse_p60_t73_g3948);
        double cse_p13_t55_g250 = (cse_p12_t44_g152 * dx);
        double cse_p13_t56_g251 = (cse_p12_t44_g152 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p60_t38_g3913 = (cse_p565_t2_invr39387 * cse_p565_t2_invr39387 * d2e_dr2);
        double cse_p13_t57_g252 = (cse_p60_t38_g3913 * dx * dx);
        double cse_p13_t58_g253 = (cse_p60_t38_g3913 * dx * dy);
        double cse_p13_t59_g254 = (cse_p60_t38_g3913 * dx * dz);
        double cse_p13_t60_g255 = (cse_p60_t38_g3913 * dy * dy);
        double cse_p13_t61_g256 = (cse_p60_t38_g3913 * dy * dz);
        double cse_p13_t62_g257 = (cse_p60_t38_g3913 * dz * dz);
        double cse_p12_t27_g135 = (cse_p60_t38_g3913 * cse_p60_t72_g3947);
        double cse_p13_t38_g233 = (cse_p12_t27_g135 * cse_p60_t72_g3947);
        double cse_p13_t39_g234 = (cse_p12_t27_g135 * cse_p60_t73_g3948);
        double cse_p13_t40_g235 = (cse_p12_t27_g135 * cse_p60_t74_g3949);
        double cse_p13_t41_g236 = (cse_p12_t27_g135 * dx);
        double cse_p13_t42_g237 = (cse_p12_t27_g135 * dy);
        double cse_p13_t43_g238 = (cse_p12_t27_g135 * dz);
        double cse_p12_t28_g136 = (cse_p60_t38_g3913 * cse_p60_t73_g3948);
        double cse_p13_t44_g239 = (cse_p12_t28_g136 * cse_p60_t73_g3948);
        double cse_p13_t45_g240 = (cse_p12_t28_g136 * cse_p60_t74_g3949);
        double cse_p13_t46_g241 = (cse_p12_t28_g136 * dx);
        double cse_p13_t47_g242 = (cse_p12_t28_g136 * dy);
        double cse_p13_t48_g243 = (cse_p12_t28_g136 * dz);
        double cse_p12_t29_g137 = (cse_p60_t38_g3913 * cse_p60_t74_g3949);
        double cse_p13_t49_g244 = (cse_p12_t29_g137 * cse_p60_t74_g3949);
        double cse_p13_t50_g245 = (cse_p12_t29_g137 * dx);
        double cse_p13_t51_g246 = (cse_p12_t29_g137 * dy);
        double cse_p13_t52_g247 = (cse_p12_t29_g137 * dz);
        double cse_p12_t58_g166 = (cse_p13_t53_g248 + cse_p13_t59_g254);
        double cse_p12_t59_g167 = (cse_p13_t54_g249 + cse_p13_t61_g256);
        double cse_p12_t63_g171 = (cse_p13_t58_g253 + cse_p13_t67_g262);
        double cse_p12_t70_g178 = (cse_p12_t40_g148 + cse_p13_t60_g255);
        double cse_p12_t72_g180 = (cse_p12_t41_g149 + cse_p13_t62_g257);
        double cse_p12_t74_g182 = (cse_p12_t39_g147 + cse_p13_t57_g252);
        double cse_p12_t49_g157 = (cse_p13_t39_g234 + cse_p13_t67_g262);
        double cse_p12_t50_g158 = (cse_p13_t40_g235 + cse_p13_t53_g248);
        double cse_p12_t51_g159 = (cse_p13_t41_g236 + cse_p13_t63_g258);
        double cse_p12_t52_g160 = (cse_p13_t42_g237 + cse_p13_t66_g261);
        double cse_p12_t53_g161 = (cse_p13_t43_g238 + cse_p13_t55_g250);
        double cse_p12_t75_g183 = (cse_p12_t39_g147 + cse_p13_t38_g233);
        double cse_p12_t54_g162 = (cse_p13_t45_g240 + cse_p13_t54_g249);
        double cse_p12_t55_g163 = (cse_p13_t46_g241 + cse_p13_t66_g261);
        double cse_p12_t56_g164 = (cse_p13_t47_g242 + cse_p13_t64_g259);
        double cse_p12_t57_g165 = (cse_p13_t48_g243 + cse_p13_t56_g251);
        double cse_p12_t71_g179 = (cse_p12_t40_g148 + cse_p13_t44_g239);
        double cse_p12_t60_g168 = (cse_p13_t50_g245 + cse_p13_t55_g250);
        double cse_p12_t61_g169 = (cse_p13_t51_g246 + cse_p13_t56_g251);
        double cse_p12_t62_g170 = (cse_p13_t52_g247 + cse_p13_t65_g260);
        double cse_p12_t73_g181 = (cse_p12_t41_g149 + cse_p13_t49_g244);
        double h_x1_x1 = cse_p12_t74_g182;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t63_g171;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t58_g166;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t51_g159;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t55_g163;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t60_g168;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t70_g178;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t59_g167;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t52_g160;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t56_g164;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t61_g169;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t72_g180;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t53_g161;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t57_g165;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t62_g170;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t75_g183;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t49_g157;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t50_g158;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t71_g179;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t54_g162;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t73_g181;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= params.rlow)) {
          /* !COND-2 */
          double cse_p12_t7_g115 = (term.bb_low * rr);
          double cse_p11_t1_g100 = (rr * rr);
          double cse_p60_t56_g3931 = (cse_p11_t41_g45 + (-(cse_p11_t1_g100)));
          double cse_p60_t57_g3932 = (cse_p11_t42_g46 + (-(cse_p11_t1_g100)));
          double cse_p60_t58_g3933 = (cse_p11_t43_g47 + (-(cse_p11_t1_g100)));
          double cse_p60_t59_g3934 = (cse_p11_t1_g100 + (-(cse_p11_t41_g45)));
          double cse_p60_t60_g3935 = (cse_p11_t1_g100 + (-(cse_p11_t42_g46)));
          double cse_p60_t61_g3936 = (cse_p11_t1_g100 + (-(cse_p11_t43_g47)));
          double cse_p60_t71_g3946 = (-(cse_p11_t1_g100));
          double cse_p12_t5_g113 = (3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g100);
          double cse_p11_t77_g81 = (term.cc_low + cse_p12_t5_g113 + cse_p12_t7_g115 + cse_p12_t7_g115);
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p11_t1_g100 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (cse_p11_t77_g81 * cse_p1_t2_invsqrt2 * dx * params.elec_weight);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t77_g81 * cse_p1_t2_invsqrt2 * dy * params.elec_weight);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t77_g81 * cse_p1_t2_invsqrt2 * dz * params.elec_weight);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (params.elec_weight * ((cse_p1_t2_invsqrt2 * dx * ((-2.00000000000000000e+0 * term.bb_low * rr) + (-3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g100))) + (-((term.cc_low * cse_p1_t2_invsqrt2 * dx)))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (params.elec_weight * ((cse_p1_t2_invsqrt2 * dy * ((-2.00000000000000000e+0 * term.bb_low * rr) + (-3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g100))) + (-((term.cc_low * cse_p1_t2_invsqrt2 * dy)))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (params.elec_weight * ((cse_p1_t2_invsqrt2 * dz * ((-2.00000000000000000e+0 * term.bb_low * rr) + (-3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g100))) + (-((term.cc_low * cse_p1_t2_invsqrt2 * dz)))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (params.elec_weight * (term.cc_low + cse_p12_t5_g113 + (2.00000000000000000e+0 * term.bb_low * rr)));
          double d2e_dr2 = (params.elec_weight * ((2.00000000000000000e+0 * term.bb_low) + (6.00000000000000000e+0 * term.aa_low * rr)));
          double cse_p565_t2_invr39387 = (1.0 / (rr));
          double cse_p11_t1_g98 = (cse_p565_t2_invr39387 * cse_p565_t2_invr39387 * d2e_dr2);
          double cse_p60_t20_g3895 = (cse_p11_t1_g98 * (-(dx)));
          double cse_p60_t21_g3896 = (cse_p11_t1_g98 * (-(dy)));
          double cse_p60_t22_g3897 = (cse_p11_t1_g98 * (-(dz)));
          double cse_p60_t23_g3898 = (cse_p11_t1_g98 * dx * dx);
          double cse_p60_t24_g3899 = (cse_p11_t1_g98 * dx * dy);
          double cse_p60_t25_g3900 = (cse_p11_t1_g98 * dx * dz);
          double cse_p60_t26_g3901 = (cse_p11_t1_g98 * dy * dy);
          double cse_p60_t27_g3902 = (cse_p11_t1_g98 * dy * dz);
          double cse_p60_t28_g3903 = (cse_p11_t1_g98 * dz * dz);
          double cse_p14_t38_g347 = (cse_p60_t20_g3895 * (-(dx)));
          double cse_p14_t39_g348 = (cse_p60_t20_g3895 * (-(dy)));
          double cse_p14_t40_g349 = (cse_p60_t20_g3895 * (-(dz)));
          double cse_p14_t41_g350 = (cse_p60_t20_g3895 * dx);
          double cse_p14_t42_g351 = (cse_p60_t20_g3895 * dy);
          double cse_p14_t43_g352 = (cse_p60_t20_g3895 * dz);
          double cse_p14_t44_g353 = (cse_p60_t21_g3896 * (-(dy)));
          double cse_p14_t45_g354 = (cse_p60_t21_g3896 * (-(dz)));
          double cse_p14_t46_g355 = (cse_p60_t21_g3896 * dx);
          double cse_p14_t47_g356 = (cse_p60_t21_g3896 * dy);
          double cse_p14_t48_g357 = (cse_p60_t21_g3896 * dz);
          double cse_p14_t49_g358 = (cse_p60_t22_g3897 * (-(dz)));
          double cse_p14_t50_g359 = (cse_p60_t22_g3897 * dx);
          double cse_p14_t51_g360 = (cse_p60_t22_g3897 * dy);
          double cse_p14_t52_g361 = (cse_p60_t22_g3897 * dz);
          double cse_p565_t1_invr39386 = (1.0 / ((rr * rr * rr)));
          double cse_p11_t1_g99 = (cse_p565_t1_invr39386 * de_dr);
          double cse_p60_t29_g3904 = (cse_p11_t1_g99 * (cse_p11_t41_g45 + cse_p60_t71_g3946));
          double cse_p60_t30_g3905 = (cse_p11_t1_g99 * (cse_p11_t42_g46 + cse_p60_t71_g3946));
          double cse_p60_t31_g3906 = (cse_p11_t1_g99 * (cse_p11_t43_g47 + cse_p60_t71_g3946));
          double cse_p60_t32_g3907 = (cse_p11_t1_g99 * cse_p60_t59_g3934);
          double cse_p60_t33_g3908 = (cse_p11_t1_g99 * cse_p60_t60_g3935);
          double cse_p60_t34_g3909 = (cse_p11_t1_g99 * cse_p60_t61_g3936);
          double cse_p60_t35_g3910 = (cse_p11_t1_g99 * dx * dy);
          double cse_p60_t36_g3911 = (cse_p11_t1_g99 * dy * (-(dx)));
          double cse_p60_t37_g3912 = (cse_p11_t1_g99 * dz);
          double cse_p60_t55_g3930 = (cse_p60_t24_g3899 + cse_p60_t36_g3911);
          double cse_p14_t53_g362 = (cse_p60_t37_g3912 * (-(dx)));
          double cse_p14_t54_g363 = (cse_p60_t37_g3912 * (-(dy)));
          double cse_p14_t55_g364 = (cse_p60_t37_g3912 * dx);
          double cse_p14_t56_g365 = (cse_p60_t37_g3912 * dy);
          double cse_p13_t99_g294 = (cse_p60_t23_g3898 + cse_p60_t32_g3907);
          double h_x1_x1 = cse_p13_t99_g294;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t55_g3930;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t80_g275 = (cse_p14_t53_g362 + cse_p60_t25_g3900);
          double cse_p13_t81_g276 = (cse_p14_t54_g363 + cse_p60_t27_g3902);
          double h_x1_z1 = cse_p13_t80_g275;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t71_g266 = (cse_p14_t39_g348 + cse_p60_t36_g3911);
          double cse_p13_t72_g267 = (cse_p14_t40_g349 + cse_p14_t53_g362);
          double cse_p13_t73_g268 = (cse_p14_t41_g350 + cse_p60_t29_g3904);
          double cse_p13_t74_g269 = (cse_p14_t42_g351 + cse_p60_t35_g3910);
          double cse_p13_t75_g270 = (cse_p14_t43_g352 + cse_p14_t55_g364);
          double cse_p13_t100_g295 = (cse_p14_t38_g347 + cse_p60_t32_g3907);
          double h_x1_x2 = cse_p13_t73_g268;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t76_g271 = (cse_p14_t45_g354 + cse_p14_t54_g363);
          double cse_p13_t77_g272 = (cse_p14_t46_g355 + cse_p60_t35_g3910);
          double cse_p13_t78_g273 = (cse_p14_t47_g356 + cse_p60_t30_g3905);
          double cse_p13_t79_g274 = (cse_p14_t48_g357 + cse_p14_t56_g365);
          double h_x1_y2 = cse_p13_t77_g272;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t82_g277 = (cse_p14_t50_g359 + cse_p14_t55_g364);
          double cse_p13_t83_g278 = (cse_p14_t51_g360 + cse_p14_t56_g365);
          double cse_p13_t84_g279 = (cse_p14_t52_g361 + cse_p60_t31_g3906);
          double h_x1_z2 = cse_p13_t82_g277;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t95_g290 = (cse_p60_t26_g3901 + cse_p60_t33_g3908);
          double cse_p13_t96_g291 = (cse_p14_t44_g353 + cse_p60_t33_g3908);
          double h_y1_y1 = cse_p13_t95_g290;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t81_g276;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t74_g269;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t78_g273;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t83_g278;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t97_g292 = (cse_p60_t28_g3903 + cse_p60_t34_g3909);
          double cse_p13_t98_g293 = (cse_p14_t49_g358 + cse_p60_t34_g3909);
          double h_z1_z1 = cse_p13_t97_g292;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t75_g270;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t79_g274;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t84_g279;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t100_g295;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t71_g266;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t72_g267;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t96_g291;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t76_g271;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t98_g293;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t83_g191 = (-(params.eps_solvent));
            double cse_p11_t78_g82 = (cse_p12_t83_g191 + params.eps_core);
            double cse_p11_t87_g91 = (-(term.kqq));
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p12_t1_g109 = (-1.25000000000000000e-1 * rx);
            double cse_p12_t2_g110 = (1.25000000000000000e-1 * rx);
            double cse_p12_t3_g111 = (5.00000000000000000e-1 * rx * rx);
            double cse_p12_t84_g192 = (-(rx));
            double cse_p11_t73_g77 = (-2.50000000000000000e-1 + cse_p12_t1_g109 + cse_p12_t1_g109);
            double cse_p11_t74_g78 = (2.50000000000000000e-1 + cse_p12_t2_g110 + cse_p12_t2_g110);
            double cse_p11_t75_g79 = (1.00000000000000000e+0 + cse_p12_t3_g111 + rx);
            double cse_p11_t89_g93 = exp(cse_p12_t84_g192);
            double ex = cse_p11_t89_g93;
            double gx = (cse_p11_t75_g79 * ex);
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p12_t3_g111 + cse_p12_t84_g192));
            double ep = (params.eps_solvent + (cse_p11_t78_g82 * gx));
            double cse_p12_t47_g155 = (ep * ep);
            double cse_p11_t90_g94 = (1.0 / (cse_p12_t47_g155));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t78_g82 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t78_g82 * gx_pp);
            double cse_p1_t3_invr3 = (1.0 / (ep));
            double cse_p1_t4_invr24 = (cse_p1_t3_invr3 * cse_p1_t3_invr3);
            double ee = (cse_p1_t3_invr3 * term.kqq);
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p11_t1_g108 = (cse_p11_t75_g79 * cse_p11_t89_g93 * cse_p1_t2_invsqrt2);
            double g_x1 = (-((cse_p11_t78_g82 * cse_p1_t4_invr24 * params.elec_weight * term.kqq * ((-2.50000000000000000e-1 * cse_p11_t1_g108 * dx) + (cse_p11_t74_g78 * cse_p1_t2_invsqrt2 * dx * ex)))));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double g_y1 = (-((cse_p11_t78_g82 * cse_p1_t4_invr24 * params.elec_weight * term.kqq * ((-2.50000000000000000e-1 * cse_p11_t1_g108 * dy) + (cse_p11_t74_g78 * cse_p1_t2_invsqrt2 * dy * ex)))));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double g_z1 = (-((cse_p11_t78_g82 * cse_p1_t4_invr24 * params.elec_weight * term.kqq * ((-2.50000000000000000e-1 * cse_p11_t1_g108 * dz) + (cse_p11_t74_g78 * cse_p1_t2_invsqrt2 * dz * ex)))));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (-((cse_p11_t78_g82 * cse_p1_t4_invr24 * params.elec_weight * term.kqq * ((2.50000000000000000e-1 * cse_p11_t1_g108 * dx) + (cse_p11_t73_g77 * cse_p1_t2_invsqrt2 * dx * ex)))));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (-((cse_p11_t78_g82 * cse_p1_t4_invr24 * params.elec_weight * term.kqq * ((2.50000000000000000e-1 * cse_p11_t1_g108 * dy) + (cse_p11_t73_g77 * cse_p1_t2_invsqrt2 * dy * ex)))));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (-((cse_p11_t78_g82 * cse_p1_t4_invr24 * params.elec_weight * term.kqq * ((2.50000000000000000e-1 * cse_p11_t1_g108 * dz) + (cse_p11_t73_g77 * cse_p1_t2_invsqrt2 * dz * ex)))));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (cse_p11_t87_g91 * cse_p11_t90_g94 * params.elec_weight * ep_p);
            double d2e_dr2 = (cse_p11_t87_g91 * params.elec_weight * ((-2.00000000000000000e+0 * ep_p * ep_p * (1.0 / ((ep * ep * ep)))) + (cse_p11_t90_g94 * ep_pp)));
            double cse_p565_t2_invr39387 = (1.0 / (rr));
            double cse_p11_t1_g98 = (cse_p565_t2_invr39387 * cse_p565_t2_invr39387 * d2e_dr2);
            double cse_p565_t1_invr39386 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g99 = (cse_p565_t1_invr39386 * de_dr);
            double cse_p11_t1_g100 = (rr * rr);
            double cse_p11_t1_g107 = (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t41_g45))));
            double cse_p60_t66_g3941 = (cse_p11_t1_g107 + (cse_p11_t1_g98 * dx * dx));
            double h_x1_x1 = cse_p60_t66_g3941;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = ((cse_p11_t1_g98 * dx * dy) + (cse_p11_t1_g99 * dy * (-(dx))));
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p11_t1_g103 = (cse_p11_t1_g99 * dz);
            double cse_p60_t12_g3887 = (cse_p11_t1_g103 * (-(dx)));
            double cse_p60_t13_g3888 = (cse_p11_t1_g103 * (-(dy)));
            double cse_p60_t14_g3889 = (cse_p11_t1_g103 * dx);
            double cse_p60_t15_g3890 = (cse_p11_t1_g103 * dy);
            double cse_p60_t50_g3925 = (cse_p60_t12_g3887 + (cse_p11_t1_g98 * dx * dz));
            double cse_p60_t51_g3926 = (cse_p60_t13_g3888 + (cse_p11_t1_g98 * dy * dz));
            double h_x1_z1 = cse_p60_t50_g3925;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p11_t1_g101 = (cse_p11_t1_g98 * (-(dx)));
            double cse_p60_t1_g3876 = (cse_p11_t1_g101 * (-(dx)));
            double cse_p60_t2_g3877 = (cse_p11_t1_g101 * (-(dy)));
            double cse_p60_t3_g3878 = (cse_p11_t1_g101 * (-(dz)));
            double cse_p60_t4_g3879 = (cse_p11_t1_g101 * dx);
            double cse_p60_t5_g3880 = (cse_p11_t1_g101 * dy);
            double cse_p60_t6_g3881 = (cse_p11_t1_g101 * dz);
            double cse_p60_t41_g3916 = (cse_p60_t2_g3877 + (cse_p11_t1_g99 * dy * (-(dx))));
            double cse_p60_t42_g3917 = (cse_p60_t12_g3887 + cse_p60_t3_g3878);
            double cse_p60_t43_g3918 = (cse_p60_t4_g3879 + (cse_p11_t1_g99 * (cse_p11_t41_g45 + (-(cse_p11_t1_g100)))));
            double cse_p60_t44_g3919 = (cse_p60_t5_g3880 + (cse_p11_t1_g99 * dx * dy));
            double cse_p60_t45_g3920 = (cse_p60_t14_g3889 + cse_p60_t6_g3881);
            double cse_p60_t67_g3942 = (cse_p11_t1_g107 + cse_p60_t1_g3876);
            double h_x1_x2 = cse_p60_t43_g3918;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p11_t1_g102 = (cse_p11_t1_g98 * (-(dy)));
            double cse_p60_t7_g3882 = (cse_p11_t1_g102 * (-(dy)));
            double cse_p60_t8_g3883 = (cse_p11_t1_g102 * (-(dz)));
            double cse_p60_t9_g3884 = (cse_p11_t1_g102 * dx);
            double cse_p60_t10_g3885 = (cse_p11_t1_g102 * dy);
            double cse_p60_t11_g3886 = (cse_p11_t1_g102 * dz);
            double cse_p60_t46_g3921 = (cse_p60_t13_g3888 + cse_p60_t8_g3883);
            double cse_p60_t47_g3922 = (cse_p60_t9_g3884 + (cse_p11_t1_g99 * dx * dy));
            double cse_p60_t48_g3923 = (cse_p60_t10_g3885 + (cse_p11_t1_g99 * (cse_p11_t42_g46 + (-(cse_p11_t1_g100)))));
            double cse_p60_t49_g3924 = (cse_p60_t11_g3886 + cse_p60_t15_g3890);
            double h_x1_y2 = cse_p60_t47_g3922;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p11_t1_g104 = (cse_p11_t1_g98 * (-(dz)));
            double cse_p60_t16_g3891 = (cse_p11_t1_g104 * (-(dz)));
            double cse_p60_t17_g3892 = (cse_p11_t1_g104 * dx);
            double cse_p60_t18_g3893 = (cse_p11_t1_g104 * dy);
            double cse_p60_t19_g3894 = (cse_p11_t1_g104 * dz);
            double cse_p60_t52_g3927 = (cse_p60_t14_g3889 + cse_p60_t17_g3892);
            double cse_p60_t53_g3928 = (cse_p60_t15_g3890 + cse_p60_t18_g3893);
            double cse_p60_t54_g3929 = (cse_p60_t19_g3894 + (cse_p11_t1_g99 * (cse_p11_t43_g47 + (-(cse_p11_t1_g100)))));
            double h_x1_z2 = cse_p60_t52_g3927;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p11_t1_g105 = (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t42_g46))));
            double cse_p60_t62_g3937 = (cse_p11_t1_g105 + (cse_p11_t1_g98 * dy * dy));
            double cse_p60_t63_g3938 = (cse_p11_t1_g105 + cse_p60_t7_g3882);
            double h_y1_y1 = cse_p60_t62_g3937;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t51_g3926;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t44_g3919;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t48_g3923;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t53_g3928;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p11_t1_g106 = (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t43_g47))));
            double cse_p60_t64_g3939 = (cse_p11_t1_g106 + (cse_p11_t1_g98 * dz * dz));
            double cse_p60_t65_g3940 = (cse_p11_t1_g106 + cse_p60_t16_g3891);
            double h_z1_z1 = cse_p60_t64_g3939;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t45_g3920;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t49_g3924;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t54_g3929;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t67_g3942;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t41_g3916;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t42_g3917;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t63_g3938;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t46_g3921;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t65_g3940;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p12_t6_g114 = (term.bb_high * rr);
              double cse_p11_t1_g100 = (rr * rr);
              double cse_p12_t4_g112 = (3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g100);
              double cse_p11_t76_g80 = (term.cc_high + cse_p12_t4_g112 + cse_p12_t6_g114 + cse_p12_t6_g114);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p11_t1_g100 * (term.bb_high + (term.aa_high * rr)))));
              energy_added += energy;
              *energy_accumulate += energy;
              double g_x1 = (cse_p11_t76_g80 * cse_p1_t2_invsqrt2 * dx * params.elec_weight);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t76_g80 * cse_p1_t2_invsqrt2 * dy * params.elec_weight);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t76_g80 * cse_p1_t2_invsqrt2 * dz * params.elec_weight);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = (params.elec_weight * ((cse_p1_t2_invsqrt2 * dx * ((-2.00000000000000000e+0 * term.bb_high * rr) + (-3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g100))) + (-((term.cc_high * cse_p1_t2_invsqrt2 * dx)))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (params.elec_weight * ((cse_p1_t2_invsqrt2 * dy * ((-2.00000000000000000e+0 * term.bb_high * rr) + (-3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g100))) + (-((term.cc_high * cse_p1_t2_invsqrt2 * dy)))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (params.elec_weight * ((cse_p1_t2_invsqrt2 * dz * ((-2.00000000000000000e+0 * term.bb_high * rr) + (-3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g100))) + (-((term.cc_high * cse_p1_t2_invsqrt2 * dz)))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double de_dr = (params.elec_weight * (term.cc_high + cse_p12_t4_g112 + (2.00000000000000000e+0 * term.bb_high * rr)));
              double d2e_dr2 = (params.elec_weight * ((2.00000000000000000e+0 * term.bb_high) + (6.00000000000000000e+0 * term.aa_high * rr)));
              double cse_p565_t2_invr39387 = (1.0 / (rr));
              double cse_p11_t1_g98 = (cse_p565_t2_invr39387 * cse_p565_t2_invr39387 * d2e_dr2);
              double cse_p565_t1_invr39386 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g99 = (cse_p565_t1_invr39386 * de_dr);
              double cse_p11_t1_g107 = (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t41_g45))));
              double h_x1_x1 = (cse_p11_t1_g107 + (cse_p11_t1_g98 * dx * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double h_x1_y1 = ((cse_p11_t1_g98 * dx * dy) + (cse_p11_t1_g99 * dy * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g103 = (cse_p11_t1_g99 * dz);
              double h_x1_z1 = ((cse_p11_t1_g103 * (-(dx))) + (cse_p11_t1_g98 * dx * dz));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g101 = (cse_p11_t1_g98 * (-(dx)));
              double h_x1_x2 = ((cse_p11_t1_g101 * dx) + (cse_p11_t1_g99 * (cse_p11_t41_g45 + (-(cse_p11_t1_g100)))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g102 = (cse_p11_t1_g98 * (-(dy)));
              double h_x1_y2 = ((cse_p11_t1_g102 * dx) + (cse_p11_t1_g99 * dx * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g104 = (cse_p11_t1_g98 * (-(dz)));
              double h_x1_z2 = ((cse_p11_t1_g103 * dx) + (cse_p11_t1_g104 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g105 = (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t42_g46))));
              double h_y1_y1 = (cse_p11_t1_g105 + (cse_p11_t1_g98 * dy * dy));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = ((cse_p11_t1_g103 * (-(dy))) + (cse_p11_t1_g98 * dy * dz));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = ((cse_p11_t1_g101 * dy) + (cse_p11_t1_g99 * dx * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = ((cse_p11_t1_g102 * dy) + (cse_p11_t1_g99 * (cse_p11_t42_g46 + (-(cse_p11_t1_g100)))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = ((cse_p11_t1_g103 * dy) + (cse_p11_t1_g104 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g106 = (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t43_g47))));
              double h_z1_z1 = (cse_p11_t1_g106 + (cse_p11_t1_g98 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = ((cse_p11_t1_g101 * dz) + (cse_p11_t1_g103 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = ((cse_p11_t1_g102 * dz) + (cse_p11_t1_g103 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = ((cse_p11_t1_g104 * dz) + (cse_p11_t1_g99 * (cse_p11_t43_g47 + (-(cse_p11_t1_g100)))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p11_t1_g107 + (cse_p11_t1_g101 * (-(dx))));
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
                double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t41_g45))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t2_invr39387 = (1.0 / (rr));
                double cse_p565_t1_invr39386 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g100 = (rr * rr);
                double cse_p11_t1_g99 = (cse_p565_t1_invr39386 * de_dr);
                double cse_p11_t1_g98 = (cse_p565_t2_invr39387 * cse_p565_t2_invr39387 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g98 * dx * dy) + (cse_p11_t1_g99 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double h_x1_z1 = ((cse_p11_t1_g98 * dx * dz) + (cse_p11_t1_g99 * dz * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g103 = (cse_p11_t1_g99 * dz);
                double h_x1_x2 = ((cse_p11_t1_g98 * dx * (-(dx))) + (cse_p11_t1_g99 * (cse_p11_t41_g45 + (-(cse_p11_t1_g100)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g101 = (cse_p11_t1_g98 * (-(dx)));
                double h_x1_y2 = ((cse_p11_t1_g98 * dx * (-(dy))) + (cse_p11_t1_g99 * dx * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g102 = (cse_p11_t1_g98 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g103 * dx) + (cse_p11_t1_g98 * dx * (-(dz))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g104 = (cse_p11_t1_g98 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g98 * dy * dy) + (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t42_g46)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g105 = (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t42_g46))));
                double h_y1_z1 = ((cse_p11_t1_g103 * (-(dy))) + (cse_p11_t1_g98 * dy * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g101 * dy) + (cse_p11_t1_g99 * dx * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g102 * dy) + (cse_p11_t1_g99 * (cse_p11_t42_g46 + (-(cse_p11_t1_g100)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g103 * dy) + (cse_p11_t1_g104 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g98 * dz * dz) + (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t43_g47)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g106 = (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t43_g47))));
                double h_z1_x2 = ((cse_p11_t1_g101 * dz) + (cse_p11_t1_g103 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g102 * dz) + (cse_p11_t1_g103 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g104 * dz) + (cse_p11_t1_g99 * (cse_p11_t43_g47 + (-(cse_p11_t1_g100)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g101 * (-(dx))) + (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t41_g45)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p11_t1_g107 = (cse_p11_t1_g99 * (cse_p11_t1_g100 + (-(cse_p11_t41_g45))));
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
