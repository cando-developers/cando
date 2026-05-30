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
      double cse_p1_t1_sqrt201642 = sqrt(r2);
      double cse_p1_t2_invsqrt201643 = (1.0 / (cse_p1_t1_sqrt201642));
      double rr = cse_p1_t1_sqrt201642;
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
          double cse_p528_t1_g208673 = (cse_p1_t2_invsqrt201643 * dx);
          double cse_p528_t2_g208661 = (cse_p1_t2_invsqrt201643 * dy);
          double cse_p528_t3_g208662 = (cse_p1_t2_invsqrt201643 * dz);
          double cse_p528_t4_g208663 = (rr * rr);
          double cse_p12_t1_g201683 = (term.aa_low * cse_p528_t4_g208663);
          double cse_p12_t4_g201675 = (3.00000000000000000e+0 * cse_p12_t1_g201683);
          double cse_p11_t1_g201667 = (term.bb_low * rr);
          double cse_p11_t12_g201657 = (term.cc_low + cse_p11_t1_g201667 + cse_p11_t1_g201667 + cse_p12_t4_g201675);
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p528_t4_g208663 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p12_t1_g201687 = (cse_p11_t12_g201657 * params.elec_weight);
          double g_x1 = (cse_p12_t1_g201687 * cse_p528_t1_g208673);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p12_t1_g201687 * cse_p528_t2_g208661);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p12_t1_g201687 * cse_p528_t3_g208662);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double cse_p63_t1_g202372 = (-2.00000000000000000e+0 * cse_p11_t1_g201667);
          double cse_p63_t1_g202374 = (-3.00000000000000000e+0 * cse_p12_t1_g201683);
          double g_x2 = (params.elec_weight * ((cse_p528_t1_g208673 * cse_p63_t1_g202372) + (cse_p528_t1_g208673 * cse_p63_t1_g202374) + (-((term.cc_low * cse_p528_t1_g208673)))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (params.elec_weight * ((cse_p528_t2_g208661 * cse_p63_t1_g202372) + (cse_p528_t2_g208661 * cse_p63_t1_g202374) + (-((term.cc_low * cse_p528_t2_g208661)))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (params.elec_weight * ((cse_p528_t3_g208662 * cse_p63_t1_g202372) + (cse_p528_t3_g208662 * cse_p63_t1_g202374) + (-((term.cc_low * cse_p528_t3_g208662)))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (params.elec_weight * (term.cc_low + cse_p12_t4_g201675 + (2.00000000000000000e+0 * cse_p11_t1_g201667)));
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t10_g201681 = (-(params.eps_solvent));
            double cse_p11_t13_g201658 = (cse_p12_t10_g201681 + params.eps_core);
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p12_t1_g201672 = (-1.25000000000000000e-1 * rx);
            double cse_p12_t2_g201673 = (1.25000000000000000e-1 * rx);
            double cse_p12_t3_g201674 = (5.00000000000000000e-1 * rx * rx);
            double cse_p12_t11_g201682 = (-(rx));
            double cse_p11_t8_g201653 = (-2.50000000000000000e-1 + cse_p12_t1_g201672 + cse_p12_t1_g201672);
            double cse_p11_t9_g201654 = (2.50000000000000000e-1 + cse_p12_t2_g201673 + cse_p12_t2_g201673);
            double cse_p11_t10_g201655 = (1.00000000000000000e+0 + cse_p12_t3_g201674 + rx);
            double cse_p11_t16_g201661 = exp(cse_p12_t11_g201682);
            double ex = cse_p11_t16_g201661;
            double gx = (cse_p11_t10_g201655 * ex);
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p12_t11_g201682 + cse_p12_t3_g201674));
            double ep = (params.eps_solvent + (cse_p11_t13_g201658 * gx));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t13_g201658 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t13_g201658 * gx_pp);
            double cse_p1_t3_invr201644 = (1.0 / (ep));
            double cse_p1_t4_invr2201645 = (cse_p1_t3_invr201644 * cse_p1_t3_invr201644);
            double ee = (cse_p1_t3_invr201644 * term.kqq);
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p528_t1_g208673 = (cse_p1_t2_invsqrt201643 * dx);
            double cse_p11_t1_g201668 = (cse_p11_t10_g201655 * cse_p11_t16_g201661);
            double cse_p11_t1_g201669 = (cse_p11_t13_g201658 * cse_p1_t4_invr2201645);
            double cse_p11_t1_g201670 = (params.elec_weight * term.kqq);
            double cse_p12_t1_g201684 = (-2.50000000000000000e-1 * cse_p11_t1_g201668);
            double cse_p12_t1_g201689 = (cse_p11_t9_g201654 * ex);
            double cse_p63_t1_g202370 = (cse_p11_t1_g201669 * cse_p11_t1_g201670);
            double g_x1 = (-((cse_p63_t1_g202370 * ((cse_p12_t1_g201684 * cse_p528_t1_g208673) + (cse_p12_t1_g201689 * cse_p528_t1_g208673)))));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p11_t1_g201663 = (cse_p1_t2_invsqrt201643 * dy);
            double g_y1 = (-((cse_p63_t1_g202370 * ((cse_p11_t1_g201663 * cse_p12_t1_g201684) + (cse_p11_t1_g201663 * cse_p12_t1_g201689)))));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p11_t1_g201664 = (cse_p1_t2_invsqrt201643 * dz);
            double g_z1 = (-((cse_p63_t1_g202370 * ((cse_p11_t1_g201664 * cse_p12_t1_g201684) + (cse_p11_t1_g201664 * cse_p12_t1_g201689)))));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double cse_p12_t1_g201685 = (2.50000000000000000e-1 * cse_p11_t1_g201668);
            double cse_p12_t1_g201688 = (cse_p11_t8_g201653 * ex);
            double g_x2 = (-((cse_p63_t1_g202370 * ((cse_p12_t1_g201685 * cse_p528_t1_g208673) + (cse_p12_t1_g201688 * cse_p528_t1_g208673)))));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (-((cse_p63_t1_g202370 * ((cse_p11_t1_g201663 * cse_p12_t1_g201685) + (cse_p11_t1_g201663 * cse_p12_t1_g201688)))));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (-((cse_p63_t1_g202370 * ((cse_p11_t1_g201664 * cse_p12_t1_g201685) + (cse_p11_t1_g201664 * cse_p12_t1_g201688)))));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (params.elec_weight * ep_p * (1.0 / ((ep * ep))) * (-(term.kqq)));
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p11_t1_g201665 = (rr * rr);
              double cse_p11_t1_g201671 = (term.aa_high * cse_p11_t1_g201665);
              double cse_p12_t5_g201676 = (3.00000000000000000e+0 * cse_p11_t1_g201671);
              double cse_p11_t1_g201666 = (term.bb_high * rr);
              double cse_p11_t11_g201656 = (term.cc_high + cse_p11_t1_g201666 + cse_p11_t1_g201666 + cse_p12_t5_g201676);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p11_t1_g201665 * (term.bb_high + (term.aa_high * rr)))));
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p528_t1_g208673 = (cse_p1_t2_invsqrt201643 * dx);
              double cse_p12_t1_g201686 = (cse_p11_t11_g201656 * params.elec_weight);
              double g_x1 = (cse_p12_t1_g201686 * cse_p528_t1_g208673);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double cse_p11_t1_g201663 = (cse_p1_t2_invsqrt201643 * dy);
              double g_y1 = (cse_p11_t1_g201663 * cse_p12_t1_g201686);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double cse_p11_t1_g201664 = (cse_p1_t2_invsqrt201643 * dz);
              double g_z1 = (cse_p11_t1_g201664 * cse_p12_t1_g201686);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p63_t1_g202371 = (-2.00000000000000000e+0 * cse_p11_t1_g201666);
              double cse_p63_t1_g202373 = (-3.00000000000000000e+0 * cse_p11_t1_g201671);
              double g_x2 = (params.elec_weight * ((cse_p528_t1_g208673 * cse_p63_t1_g202371) + (cse_p528_t1_g208673 * cse_p63_t1_g202373) + (-((term.cc_high * cse_p528_t1_g208673)))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (params.elec_weight * ((cse_p11_t1_g201663 * cse_p63_t1_g202371) + (cse_p11_t1_g201663 * cse_p63_t1_g202373) + (-((term.cc_high * cse_p11_t1_g201663)))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (params.elec_weight * ((cse_p11_t1_g201664 * cse_p63_t1_g202371) + (cse_p11_t1_g201664 * cse_p63_t1_g202373) + (-((term.cc_high * cse_p11_t1_g201664)))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double de_dr = (params.elec_weight * (term.cc_high + cse_p12_t5_g201676 + (2.00000000000000000e+0 * cse_p11_t1_g201666)));
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
          double cse_p528_t1_g214187 = (rr * rr);
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p528_t1_g214187 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t2_g208690 = (-(params.eps_solvent));
            double cse_p11_t2_g208676 = (cse_p12_t2_g208690 + params.eps_core);
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p11_t1_g208675 = (5.00000000000000000e-1 * rx * rx);
            double cse_p11_t4_g208678 = (-(rx));
            double ex = exp(cse_p11_t4_g208678);
            double gx = (ex * (1.00000000000000000e+0 + cse_p11_t1_g208675 + rx));
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p11_t1_g208675 + cse_p11_t4_g208678));
            double ep = (params.eps_solvent + (cse_p11_t2_g208676 * gx));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t2_g208676 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t2_g208676 * gx_pp);
            double ee = (term.kqq * (1.0 / (ep)));
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p528_t1_g214187 = (rr * rr);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p528_t1_g214187 * (term.bb_high + (term.aa_high * rr)))));
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
      double cse_p11_t41_g214233 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t42_g214234 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t43_g214235 = (dz * dz);
      double r2 = (cse_p11_t41_g214233 + cse_p11_t42_g214234 + cse_p11_t43_g214235);
      double cse_p1_t1_sqrt214189 = sqrt(r2);
      double cse_p1_t2_invsqrt214190 = (1.0 / (cse_p1_t1_sqrt214189));
      double rr = cse_p1_t1_sqrt214189;
      if ((rr <= params.rmin)) {
        /* !COND-1 */
        double cse_p60_t40_g218103 = (rr * rr);
        double cse_p60_t68_g218131 = (-(cse_p11_t41_g214233));
        double cse_p60_t69_g218132 = (-(cse_p11_t42_g214234));
        double cse_p60_t70_g218133 = (-(cse_p11_t43_g214235));
        double cse_p60_t72_g218135 = (-(dx));
        double cse_p60_t73_g218136 = (-(dy));
        double cse_p60_t74_g218137 = (-(dz));
        double cse_p565_t1_invr253574 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr253575 = (1.0 / (rr));
        double cse_p60_t77_g218140 = (rr * rr * rr);
        double cse_p15_t72_g214669 = (-(cse_p60_t40_g218103));
        double cse_p14_t76_g214573 = (cse_p11_t41_g214233 + cse_p15_t72_g214669);
        double cse_p14_t78_g214575 = (cse_p11_t42_g214234 + cse_p15_t72_g214669);
        double cse_p14_t80_g214577 = (cse_p11_t43_g214235 + cse_p15_t72_g214669);
        double cse_p13_t101_g214484 = (cse_p60_t40_g218103 + cse_p60_t68_g218131);
        double cse_p13_t102_g214485 = (cse_p60_t40_g218103 + cse_p60_t69_g218132);
        double cse_p13_t103_g214486 = (cse_p60_t40_g218103 + cse_p60_t70_g218133);
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
        double cse_p60_t39_g218102 = (cse_p565_t1_invr253574 * de_dr);
        double cse_p13_t52_g214435 = (cse_p14_t76_g214573 * cse_p60_t39_g218102);
        double cse_p13_t53_g214436 = (cse_p14_t78_g214575 * cse_p60_t39_g218102);
        double cse_p13_t54_g214437 = (cse_p14_t80_g214577 * cse_p60_t39_g218102);
        double cse_p13_t55_g214438 = (cse_p60_t39_g218102 * dx * dy);
        double cse_p13_t56_g214439 = (cse_p60_t39_g218102 * cse_p60_t72_g218135 * dy);
        double cse_p12_t20_g214316 = (cse_p13_t101_g214484 * cse_p60_t39_g218102);
        double cse_p12_t21_g214317 = (cse_p13_t102_g214485 * cse_p60_t39_g218102);
        double cse_p12_t22_g214318 = (cse_p13_t103_g214486 * cse_p60_t39_g218102);
        double cse_p12_t25_g214321 = (cse_p60_t39_g218102 * dz);
        double cse_p13_t42_g214425 = (cse_p12_t25_g214321 * cse_p60_t72_g218135);
        double cse_p13_t43_g214426 = (cse_p12_t25_g214321 * cse_p60_t73_g218136);
        double cse_p13_t44_g214427 = (cse_p12_t25_g214321 * dx);
        double cse_p13_t45_g214428 = (cse_p12_t25_g214321 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p60_t38_g218101 = (cse_p565_t2_invr253575 * cse_p565_t2_invr253575 * d2e_dr2);
        double cse_p13_t46_g214429 = (cse_p60_t38_g218101 * dx * dx);
        double cse_p13_t47_g214430 = (cse_p60_t38_g218101 * dx * dy);
        double cse_p13_t48_g214431 = (cse_p60_t38_g218101 * dx * dz);
        double cse_p13_t49_g214432 = (cse_p60_t38_g218101 * dy * dy);
        double cse_p13_t50_g214433 = (cse_p60_t38_g218101 * dy * dz);
        double cse_p13_t51_g214434 = (cse_p60_t38_g218101 * dz * dz);
        double cse_p12_t8_g214304 = (cse_p60_t38_g218101 * cse_p60_t72_g218135);
        double cse_p13_t57_g214440 = (cse_p12_t8_g214304 * cse_p60_t72_g218135);
        double cse_p13_t58_g214441 = (cse_p12_t8_g214304 * cse_p60_t73_g218136);
        double cse_p13_t59_g214442 = (cse_p12_t8_g214304 * cse_p60_t74_g218137);
        double cse_p13_t60_g214443 = (cse_p12_t8_g214304 * dx);
        double cse_p13_t61_g214444 = (cse_p12_t8_g214304 * dy);
        double cse_p13_t62_g214445 = (cse_p12_t8_g214304 * dz);
        double cse_p12_t9_g214305 = (cse_p60_t38_g218101 * cse_p60_t73_g218136);
        double cse_p13_t63_g214446 = (cse_p12_t9_g214305 * cse_p60_t73_g218136);
        double cse_p13_t64_g214447 = (cse_p12_t9_g214305 * cse_p60_t74_g218137);
        double cse_p13_t65_g214448 = (cse_p12_t9_g214305 * dx);
        double cse_p13_t66_g214449 = (cse_p12_t9_g214305 * dy);
        double cse_p13_t67_g214450 = (cse_p12_t9_g214305 * dz);
        double cse_p12_t10_g214306 = (cse_p60_t38_g218101 * cse_p60_t74_g218137);
        double cse_p13_t38_g214421 = (cse_p12_t10_g214306 * cse_p60_t74_g218137);
        double cse_p13_t39_g214422 = (cse_p12_t10_g214306 * dx);
        double cse_p13_t40_g214423 = (cse_p12_t10_g214306 * dy);
        double cse_p13_t41_g214424 = (cse_p12_t10_g214306 * dz);
        double cse_p12_t52_g214348 = (cse_p13_t47_g214430 + cse_p13_t56_g214439);
        double cse_p12_t53_g214349 = (cse_p13_t42_g214425 + cse_p13_t48_g214431);
        double cse_p12_t54_g214350 = (cse_p13_t43_g214426 + cse_p13_t50_g214433);
        double cse_p12_t70_g214366 = (cse_p12_t21_g214317 + cse_p13_t49_g214432);
        double cse_p12_t72_g214368 = (cse_p12_t22_g214318 + cse_p13_t51_g214434);
        double cse_p12_t74_g214370 = (cse_p12_t20_g214316 + cse_p13_t46_g214429);
        double cse_p12_t49_g214345 = (cse_p13_t52_g214435 + cse_p13_t60_g214443);
        double cse_p12_t55_g214351 = (cse_p13_t55_g214438 + cse_p13_t61_g214444);
        double cse_p12_t57_g214353 = (cse_p13_t56_g214439 + cse_p13_t58_g214441);
        double cse_p12_t58_g214354 = (cse_p13_t42_g214425 + cse_p13_t59_g214442);
        double cse_p12_t59_g214355 = (cse_p13_t44_g214427 + cse_p13_t62_g214445);
        double cse_p12_t75_g214371 = (cse_p12_t20_g214316 + cse_p13_t57_g214440);
        double cse_p12_t50_g214346 = (cse_p13_t53_g214436 + cse_p13_t66_g214449);
        double cse_p12_t56_g214352 = (cse_p13_t55_g214438 + cse_p13_t65_g214448);
        double cse_p12_t60_g214356 = (cse_p13_t43_g214426 + cse_p13_t64_g214447);
        double cse_p12_t61_g214357 = (cse_p13_t45_g214428 + cse_p13_t67_g214450);
        double cse_p12_t71_g214367 = (cse_p12_t21_g214317 + cse_p13_t63_g214446);
        double cse_p12_t51_g214347 = (cse_p13_t41_g214424 + cse_p13_t54_g214437);
        double cse_p12_t62_g214358 = (cse_p13_t39_g214422 + cse_p13_t44_g214427);
        double cse_p12_t63_g214359 = (cse_p13_t40_g214423 + cse_p13_t45_g214428);
        double cse_p12_t73_g214369 = (cse_p12_t22_g214318 + cse_p13_t38_g214421);
        double h_x1_x1 = cse_p12_t74_g214370;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t52_g214348;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t53_g214349;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t49_g214345;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t56_g214352;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t62_g214358;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t70_g214366;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t54_g214350;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t55_g214351;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t50_g214346;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t63_g214359;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t72_g214368;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t59_g214355;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t61_g214357;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t51_g214347;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t75_g214371;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t57_g214353;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t58_g214354;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t71_g214367;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t60_g214356;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t73_g214369;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= params.rlow)) {
          /* !COND-2 */
          double cse_p12_t7_g214303 = (term.bb_low * rr);
          double cse_p11_t1_g214288 = (rr * rr);
          double cse_p60_t56_g218119 = (cse_p11_t41_g214233 + (-(cse_p11_t1_g214288)));
          double cse_p60_t57_g218120 = (cse_p11_t42_g214234 + (-(cse_p11_t1_g214288)));
          double cse_p60_t58_g218121 = (cse_p11_t43_g214235 + (-(cse_p11_t1_g214288)));
          double cse_p60_t59_g218122 = (cse_p11_t1_g214288 + (-(cse_p11_t41_g214233)));
          double cse_p60_t60_g218123 = (cse_p11_t1_g214288 + (-(cse_p11_t42_g214234)));
          double cse_p60_t61_g218124 = (cse_p11_t1_g214288 + (-(cse_p11_t43_g214235)));
          double cse_p60_t71_g218134 = (-(cse_p11_t1_g214288));
          double cse_p12_t5_g214301 = (3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g214288);
          double cse_p11_t77_g214269 = (term.cc_low + cse_p12_t5_g214301 + cse_p12_t7_g214303 + cse_p12_t7_g214303);
          double energy = (params.elec_weight * (term.dd_low + (term.cc_low * rr) + (cse_p11_t1_g214288 * (term.bb_low + (term.aa_low * rr)))));
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (cse_p11_t77_g214269 * cse_p1_t2_invsqrt214190 * dx * params.elec_weight);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t77_g214269 * cse_p1_t2_invsqrt214190 * dy * params.elec_weight);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t77_g214269 * cse_p1_t2_invsqrt214190 * dz * params.elec_weight);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (params.elec_weight * ((cse_p1_t2_invsqrt214190 * dx * ((-2.00000000000000000e+0 * term.bb_low * rr) + (-3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g214288))) + (-((term.cc_low * cse_p1_t2_invsqrt214190 * dx)))));
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (params.elec_weight * ((cse_p1_t2_invsqrt214190 * dy * ((-2.00000000000000000e+0 * term.bb_low * rr) + (-3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g214288))) + (-((term.cc_low * cse_p1_t2_invsqrt214190 * dy)))));
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (params.elec_weight * ((cse_p1_t2_invsqrt214190 * dz * ((-2.00000000000000000e+0 * term.bb_low * rr) + (-3.00000000000000000e+0 * term.aa_low * cse_p11_t1_g214288))) + (-((term.cc_low * cse_p1_t2_invsqrt214190 * dz)))));
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double de_dr = (params.elec_weight * (term.cc_low + cse_p12_t5_g214301 + (2.00000000000000000e+0 * term.bb_low * rr)));
          double d2e_dr2 = (params.elec_weight * ((2.00000000000000000e+0 * term.bb_low) + (6.00000000000000000e+0 * term.aa_low * rr)));
          double cse_p565_t2_invr253575 = (1.0 / (rr));
          double cse_p11_t1_g214286 = (cse_p565_t2_invr253575 * cse_p565_t2_invr253575 * d2e_dr2);
          double cse_p60_t1_g218064 = (cse_p11_t1_g214286 * (-(dx)));
          double cse_p60_t2_g218065 = (cse_p11_t1_g214286 * (-(dy)));
          double cse_p60_t3_g218066 = (cse_p11_t1_g214286 * (-(dz)));
          double cse_p60_t4_g218067 = (cse_p11_t1_g214286 * dx * dx);
          double cse_p60_t5_g218068 = (cse_p11_t1_g214286 * dx * dy);
          double cse_p60_t6_g218069 = (cse_p11_t1_g214286 * dx * dz);
          double cse_p60_t7_g218070 = (cse_p11_t1_g214286 * dy * dy);
          double cse_p60_t8_g218071 = (cse_p11_t1_g214286 * dy * dz);
          double cse_p60_t9_g218072 = (cse_p11_t1_g214286 * dz * dz);
          double cse_p14_t42_g214539 = (cse_p60_t1_g218064 * (-(dx)));
          double cse_p14_t43_g214540 = (cse_p60_t1_g218064 * (-(dy)));
          double cse_p14_t44_g214541 = (cse_p60_t1_g218064 * (-(dz)));
          double cse_p14_t45_g214542 = (cse_p60_t1_g218064 * dx);
          double cse_p14_t46_g214543 = (cse_p60_t1_g218064 * dy);
          double cse_p14_t47_g214544 = (cse_p60_t1_g218064 * dz);
          double cse_p14_t48_g214545 = (cse_p60_t2_g218065 * (-(dy)));
          double cse_p14_t49_g214546 = (cse_p60_t2_g218065 * (-(dz)));
          double cse_p14_t50_g214547 = (cse_p60_t2_g218065 * dx);
          double cse_p14_t51_g214548 = (cse_p60_t2_g218065 * dy);
          double cse_p14_t52_g214549 = (cse_p60_t2_g218065 * dz);
          double cse_p14_t53_g214550 = (cse_p60_t3_g218066 * (-(dz)));
          double cse_p14_t54_g214551 = (cse_p60_t3_g218066 * dx);
          double cse_p14_t55_g214552 = (cse_p60_t3_g218066 * dy);
          double cse_p14_t56_g214553 = (cse_p60_t3_g218066 * dz);
          double cse_p565_t1_invr253574 = (1.0 / ((rr * rr * rr)));
          double cse_p11_t1_g214287 = (cse_p565_t1_invr253574 * de_dr);
          double cse_p60_t10_g218073 = (cse_p11_t1_g214287 * (cse_p11_t41_g214233 + cse_p60_t71_g218134));
          double cse_p60_t11_g218074 = (cse_p11_t1_g214287 * (cse_p11_t42_g214234 + cse_p60_t71_g218134));
          double cse_p60_t12_g218075 = (cse_p11_t1_g214287 * (cse_p11_t43_g214235 + cse_p60_t71_g218134));
          double cse_p60_t13_g218076 = (cse_p11_t1_g214287 * cse_p60_t59_g218122);
          double cse_p60_t14_g218077 = (cse_p11_t1_g214287 * cse_p60_t60_g218123);
          double cse_p60_t15_g218078 = (cse_p11_t1_g214287 * cse_p60_t61_g218124);
          double cse_p60_t16_g218079 = (cse_p11_t1_g214287 * dx * dy);
          double cse_p60_t17_g218080 = (cse_p11_t1_g214287 * dy * (-(dx)));
          double cse_p60_t18_g218081 = (cse_p11_t1_g214287 * dz);
          double cse_p60_t44_g218107 = (cse_p60_t17_g218080 + cse_p60_t5_g218068);
          double cse_p14_t38_g214535 = (cse_p60_t18_g218081 * (-(dx)));
          double cse_p14_t39_g214536 = (cse_p60_t18_g218081 * (-(dy)));
          double cse_p14_t40_g214537 = (cse_p60_t18_g218081 * dx);
          double cse_p14_t41_g214538 = (cse_p60_t18_g218081 * dy);
          double cse_p13_t99_g214482 = (cse_p60_t13_g218076 + cse_p60_t4_g218067);
          double h_x1_x1 = cse_p13_t99_g214482;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t44_g218107;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t75_g214458 = (cse_p14_t38_g214535 + cse_p60_t6_g218069);
          double cse_p13_t76_g214459 = (cse_p14_t39_g214536 + cse_p60_t8_g218071);
          double h_x1_z1 = cse_p13_t75_g214458;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t71_g214454 = (cse_p14_t45_g214542 + cse_p60_t10_g218073);
          double cse_p13_t77_g214460 = (cse_p14_t46_g214543 + cse_p60_t16_g218079);
          double cse_p13_t79_g214462 = (cse_p14_t43_g214540 + cse_p60_t17_g218080);
          double cse_p13_t80_g214463 = (cse_p14_t38_g214535 + cse_p14_t44_g214541);
          double cse_p13_t81_g214464 = (cse_p14_t40_g214537 + cse_p14_t47_g214544);
          double cse_p13_t100_g214483 = (cse_p14_t42_g214539 + cse_p60_t13_g218076);
          double h_x1_x2 = cse_p13_t71_g214454;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t72_g214455 = (cse_p14_t51_g214548 + cse_p60_t11_g218074);
          double cse_p13_t78_g214461 = (cse_p14_t50_g214547 + cse_p60_t16_g218079);
          double cse_p13_t82_g214465 = (cse_p14_t39_g214536 + cse_p14_t49_g214546);
          double cse_p13_t83_g214466 = (cse_p14_t41_g214538 + cse_p14_t52_g214549);
          double h_x1_y2 = cse_p13_t78_g214461;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t73_g214456 = (cse_p14_t56_g214553 + cse_p60_t12_g218075);
          double cse_p13_t84_g214467 = (cse_p14_t40_g214537 + cse_p14_t54_g214551);
          double cse_p13_t85_g214468 = (cse_p14_t41_g214538 + cse_p14_t55_g214552);
          double h_x1_z2 = cse_p13_t84_g214467;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t95_g214478 = (cse_p60_t14_g218077 + cse_p60_t7_g218070);
          double cse_p13_t96_g214479 = (cse_p14_t48_g214545 + cse_p60_t14_g218077);
          double h_y1_y1 = cse_p13_t95_g214478;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t76_g214459;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t77_g214460;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t72_g214455;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t85_g214468;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t97_g214480 = (cse_p60_t15_g218078 + cse_p60_t9_g218072);
          double cse_p13_t98_g214481 = (cse_p14_t53_g214550 + cse_p60_t15_g218078);
          double h_z1_z1 = cse_p13_t97_g214480;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t81_g214464;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t83_g214466;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t73_g214456;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t100_g214483;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t79_g214462;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t80_g214463;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t96_g214479;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t82_g214465;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t98_g214481;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rhi)) {
            /* !COND-3 */
            double cse_p12_t83_g214379 = (-(params.eps_solvent));
            double cse_p11_t78_g214270 = (cse_p12_t83_g214379 + params.eps_core);
            double cse_p11_t87_g214279 = (-(term.kqq));
            double rx = (2.50000000000000000e-1 * rr);
            double cse_p12_t1_g214297 = (-1.25000000000000000e-1 * rx);
            double cse_p12_t2_g214298 = (1.25000000000000000e-1 * rx);
            double cse_p12_t3_g214299 = (5.00000000000000000e-1 * rx * rx);
            double cse_p12_t84_g214380 = (-(rx));
            double cse_p11_t73_g214265 = (-2.50000000000000000e-1 + cse_p12_t1_g214297 + cse_p12_t1_g214297);
            double cse_p11_t74_g214266 = (2.50000000000000000e-1 + cse_p12_t2_g214298 + cse_p12_t2_g214298);
            double cse_p11_t75_g214267 = (1.00000000000000000e+0 + cse_p12_t3_g214299 + rx);
            double cse_p11_t89_g214281 = exp(cse_p12_t84_g214380);
            double ex = cse_p11_t89_g214281;
            double gx = (cse_p11_t75_g214267 * ex);
            double gx_p = (-5.00000000000000000e-1 * ex * rx * rx);
            double gx_pp = (ex * (cse_p12_t3_g214299 + cse_p12_t84_g214380));
            double ep = (params.eps_solvent + (cse_p11_t78_g214270 * gx));
            double cse_p12_t47_g214343 = (ep * ep);
            double cse_p11_t90_g214282 = (1.0 / (cse_p12_t47_g214343));
            double ep_p = (2.50000000000000000e-1 * cse_p11_t78_g214270 * gx_p);
            double ep_pp = (6.25000000000000000e-2 * cse_p11_t78_g214270 * gx_pp);
            double cse_p1_t3_invr214191 = (1.0 / (ep));
            double cse_p1_t4_invr2214192 = (cse_p1_t3_invr214191 * cse_p1_t3_invr214191);
            double ee = (cse_p1_t3_invr214191 * term.kqq);
            double energy = (ee * params.elec_weight);
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p11_t1_g214296 = (cse_p11_t75_g214267 * cse_p11_t89_g214281 * cse_p1_t2_invsqrt214190);
            double g_x1 = (-((cse_p11_t78_g214270 * cse_p1_t4_invr2214192 * params.elec_weight * term.kqq * ((-2.50000000000000000e-1 * cse_p11_t1_g214296 * dx) + (cse_p11_t74_g214266 * cse_p1_t2_invsqrt214190 * dx * ex)))));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double g_y1 = (-((cse_p11_t78_g214270 * cse_p1_t4_invr2214192 * params.elec_weight * term.kqq * ((-2.50000000000000000e-1 * cse_p11_t1_g214296 * dy) + (cse_p11_t74_g214266 * cse_p1_t2_invsqrt214190 * dy * ex)))));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double g_z1 = (-((cse_p11_t78_g214270 * cse_p1_t4_invr2214192 * params.elec_weight * term.kqq * ((-2.50000000000000000e-1 * cse_p11_t1_g214296 * dz) + (cse_p11_t74_g214266 * cse_p1_t2_invsqrt214190 * dz * ex)))));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (-((cse_p11_t78_g214270 * cse_p1_t4_invr2214192 * params.elec_weight * term.kqq * ((2.50000000000000000e-1 * cse_p11_t1_g214296 * dx) + (cse_p11_t73_g214265 * cse_p1_t2_invsqrt214190 * dx * ex)))));
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (-((cse_p11_t78_g214270 * cse_p1_t4_invr2214192 * params.elec_weight * term.kqq * ((2.50000000000000000e-1 * cse_p11_t1_g214296 * dy) + (cse_p11_t73_g214265 * cse_p1_t2_invsqrt214190 * dy * ex)))));
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (-((cse_p11_t78_g214270 * cse_p1_t4_invr2214192 * params.elec_weight * term.kqq * ((2.50000000000000000e-1 * cse_p11_t1_g214296 * dz) + (cse_p11_t73_g214265 * cse_p1_t2_invsqrt214190 * dz * ex)))));
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double de_dr = (cse_p11_t87_g214279 * cse_p11_t90_g214282 * params.elec_weight * ep_p);
            double d2e_dr2 = (cse_p11_t87_g214279 * params.elec_weight * ((-2.00000000000000000e+0 * ep_p * ep_p * (1.0 / ((ep * ep * ep)))) + (cse_p11_t90_g214282 * ep_pp)));
            double cse_p565_t2_invr253575 = (1.0 / (rr));
            double cse_p11_t1_g214286 = (cse_p565_t2_invr253575 * cse_p565_t2_invr253575 * d2e_dr2);
            double cse_p565_t1_invr253574 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g214287 = (cse_p565_t1_invr253574 * de_dr);
            double cse_p11_t1_g214288 = (rr * rr);
            double cse_p11_t1_g214295 = (cse_p11_t1_g214287 * (cse_p11_t1_g214288 + (-(cse_p11_t41_g214233))));
            double cse_p60_t66_g218129 = (cse_p11_t1_g214295 + (cse_p11_t1_g214286 * dx * dx));
            double h_x1_x1 = cse_p60_t66_g218129;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = ((cse_p11_t1_g214286 * dx * dy) + (cse_p11_t1_g214287 * dy * (-(dx))));
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p11_t1_g214291 = (cse_p11_t1_g214287 * dz);
            double cse_p60_t30_g218093 = (cse_p11_t1_g214291 * (-(dx)));
            double cse_p60_t31_g218094 = (cse_p11_t1_g214291 * (-(dy)));
            double cse_p60_t32_g218095 = (cse_p11_t1_g214291 * dx);
            double cse_p60_t33_g218096 = (cse_p11_t1_g214291 * dy);
            double cse_p60_t45_g218108 = (cse_p60_t30_g218093 + (cse_p11_t1_g214286 * dx * dz));
            double cse_p60_t46_g218109 = (cse_p60_t31_g218094 + (cse_p11_t1_g214286 * dy * dz));
            double h_x1_z1 = cse_p60_t45_g218108;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p11_t1_g214289 = (cse_p11_t1_g214286 * (-(dx)));
            double cse_p60_t19_g218082 = (cse_p11_t1_g214289 * (-(dx)));
            double cse_p60_t20_g218083 = (cse_p11_t1_g214289 * (-(dy)));
            double cse_p60_t21_g218084 = (cse_p11_t1_g214289 * (-(dz)));
            double cse_p60_t22_g218085 = (cse_p11_t1_g214289 * dx);
            double cse_p60_t23_g218086 = (cse_p11_t1_g214289 * dy);
            double cse_p60_t24_g218087 = (cse_p11_t1_g214289 * dz);
            double cse_p60_t41_g218104 = (cse_p60_t22_g218085 + (cse_p11_t1_g214287 * (cse_p11_t41_g214233 + (-(cse_p11_t1_g214288)))));
            double cse_p60_t47_g218110 = (cse_p60_t23_g218086 + (cse_p11_t1_g214287 * dx * dy));
            double cse_p60_t49_g218112 = (cse_p60_t20_g218083 + (cse_p11_t1_g214287 * dy * (-(dx))));
            double cse_p60_t50_g218113 = (cse_p60_t21_g218084 + cse_p60_t30_g218093);
            double cse_p60_t51_g218114 = (cse_p60_t24_g218087 + cse_p60_t32_g218095);
            double cse_p60_t67_g218130 = (cse_p11_t1_g214295 + cse_p60_t19_g218082);
            double h_x1_x2 = cse_p60_t41_g218104;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p11_t1_g214290 = (cse_p11_t1_g214286 * (-(dy)));
            double cse_p60_t25_g218088 = (cse_p11_t1_g214290 * (-(dy)));
            double cse_p60_t26_g218089 = (cse_p11_t1_g214290 * (-(dz)));
            double cse_p60_t27_g218090 = (cse_p11_t1_g214290 * dx);
            double cse_p60_t28_g218091 = (cse_p11_t1_g214290 * dy);
            double cse_p60_t29_g218092 = (cse_p11_t1_g214290 * dz);
            double cse_p60_t42_g218105 = (cse_p60_t28_g218091 + (cse_p11_t1_g214287 * (cse_p11_t42_g214234 + (-(cse_p11_t1_g214288)))));
            double cse_p60_t48_g218111 = (cse_p60_t27_g218090 + (cse_p11_t1_g214287 * dx * dy));
            double cse_p60_t52_g218115 = (cse_p60_t26_g218089 + cse_p60_t31_g218094);
            double cse_p60_t53_g218116 = (cse_p60_t29_g218092 + cse_p60_t33_g218096);
            double h_x1_y2 = cse_p60_t48_g218111;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p11_t1_g214292 = (cse_p11_t1_g214286 * (-(dz)));
            double cse_p60_t34_g218097 = (cse_p11_t1_g214292 * (-(dz)));
            double cse_p60_t35_g218098 = (cse_p11_t1_g214292 * dx);
            double cse_p60_t36_g218099 = (cse_p11_t1_g214292 * dy);
            double cse_p60_t37_g218100 = (cse_p11_t1_g214292 * dz);
            double cse_p60_t43_g218106 = (cse_p60_t37_g218100 + (cse_p11_t1_g214287 * (cse_p11_t43_g214235 + (-(cse_p11_t1_g214288)))));
            double cse_p60_t54_g218117 = (cse_p60_t32_g218095 + cse_p60_t35_g218098);
            double cse_p60_t55_g218118 = (cse_p60_t33_g218096 + cse_p60_t36_g218099);
            double h_x1_z2 = cse_p60_t54_g218117;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p11_t1_g214293 = (cse_p11_t1_g214287 * (cse_p11_t1_g214288 + (-(cse_p11_t42_g214234))));
            double cse_p60_t62_g218125 = (cse_p11_t1_g214293 + (cse_p11_t1_g214286 * dy * dy));
            double cse_p60_t63_g218126 = (cse_p11_t1_g214293 + cse_p60_t25_g218088);
            double h_y1_y1 = cse_p60_t62_g218125;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t46_g218109;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t47_g218110;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t42_g218105;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t55_g218118;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p11_t1_g214294 = (cse_p11_t1_g214287 * (cse_p11_t1_g214288 + (-(cse_p11_t43_g214235))));
            double cse_p60_t64_g218127 = (cse_p11_t1_g214294 + (cse_p11_t1_g214286 * dz * dz));
            double cse_p60_t65_g218128 = (cse_p11_t1_g214294 + cse_p60_t34_g218097);
            double h_z1_z1 = cse_p60_t64_g218127;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t51_g218114;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t53_g218116;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t43_g218106;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t67_g218130;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t49_g218112;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t50_g218113;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t63_g218126;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t52_g218115;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t65_g218128;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p12_t6_g214302 = (term.bb_high * rr);
              double cse_p11_t1_g214288 = (rr * rr);
              double cse_p12_t4_g214300 = (3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g214288);
              double cse_p11_t76_g214268 = (term.cc_high + cse_p12_t4_g214300 + cse_p12_t6_g214302 + cse_p12_t6_g214302);
              double energy = (params.elec_weight * (term.dd_high + (term.cc_high * rr) + (cse_p11_t1_g214288 * (term.bb_high + (term.aa_high * rr)))));
              energy_added += energy;
              *energy_accumulate += energy;
              double g_x1 = (cse_p11_t76_g214268 * cse_p1_t2_invsqrt214190 * dx * params.elec_weight);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t76_g214268 * cse_p1_t2_invsqrt214190 * dy * params.elec_weight);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t76_g214268 * cse_p1_t2_invsqrt214190 * dz * params.elec_weight);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = (params.elec_weight * ((cse_p1_t2_invsqrt214190 * dx * ((-2.00000000000000000e+0 * term.bb_high * rr) + (-3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g214288))) + (-((term.cc_high * cse_p1_t2_invsqrt214190 * dx)))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = (params.elec_weight * ((cse_p1_t2_invsqrt214190 * dy * ((-2.00000000000000000e+0 * term.bb_high * rr) + (-3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g214288))) + (-((term.cc_high * cse_p1_t2_invsqrt214190 * dy)))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = (params.elec_weight * ((cse_p1_t2_invsqrt214190 * dz * ((-2.00000000000000000e+0 * term.bb_high * rr) + (-3.00000000000000000e+0 * term.aa_high * cse_p11_t1_g214288))) + (-((term.cc_high * cse_p1_t2_invsqrt214190 * dz)))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double de_dr = (params.elec_weight * (term.cc_high + cse_p12_t4_g214300 + (2.00000000000000000e+0 * term.bb_high * rr)));
              double d2e_dr2 = (params.elec_weight * ((2.00000000000000000e+0 * term.bb_high) + (6.00000000000000000e+0 * term.aa_high * rr)));
              double cse_p565_t2_invr253575 = (1.0 / (rr));
              double cse_p11_t1_g214286 = (cse_p565_t2_invr253575 * cse_p565_t2_invr253575 * d2e_dr2);
              double cse_p565_t1_invr253574 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g214287 = (cse_p565_t1_invr253574 * de_dr);
              double cse_p11_t1_g214295 = (cse_p11_t1_g214287 * (cse_p11_t1_g214288 + (-(cse_p11_t41_g214233))));
              double h_x1_x1 = (cse_p11_t1_g214295 + (cse_p11_t1_g214286 * dx * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double h_x1_y1 = ((cse_p11_t1_g214286 * dx * dy) + (cse_p11_t1_g214287 * dy * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g214291 = (cse_p11_t1_g214287 * dz);
              double h_x1_z1 = ((cse_p11_t1_g214286 * dx * dz) + (cse_p11_t1_g214291 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g214289 = (cse_p11_t1_g214286 * (-(dx)));
              double h_x1_x2 = ((cse_p11_t1_g214287 * (cse_p11_t41_g214233 + (-(cse_p11_t1_g214288)))) + (cse_p11_t1_g214289 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g214290 = (cse_p11_t1_g214286 * (-(dy)));
              double h_x1_y2 = ((cse_p11_t1_g214287 * dx * dy) + (cse_p11_t1_g214290 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g214292 = (cse_p11_t1_g214286 * (-(dz)));
              double h_x1_z2 = ((cse_p11_t1_g214291 * dx) + (cse_p11_t1_g214292 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g214293 = (cse_p11_t1_g214287 * (cse_p11_t1_g214288 + (-(cse_p11_t42_g214234))));
              double h_y1_y1 = (cse_p11_t1_g214293 + (cse_p11_t1_g214286 * dy * dy));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = ((cse_p11_t1_g214286 * dy * dz) + (cse_p11_t1_g214291 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = ((cse_p11_t1_g214287 * dx * dy) + (cse_p11_t1_g214289 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = ((cse_p11_t1_g214287 * (cse_p11_t42_g214234 + (-(cse_p11_t1_g214288)))) + (cse_p11_t1_g214290 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = ((cse_p11_t1_g214291 * dy) + (cse_p11_t1_g214292 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g214294 = (cse_p11_t1_g214287 * (cse_p11_t1_g214288 + (-(cse_p11_t43_g214235))));
              double h_z1_z1 = (cse_p11_t1_g214294 + (cse_p11_t1_g214286 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = ((cse_p11_t1_g214289 * dz) + (cse_p11_t1_g214291 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = ((cse_p11_t1_g214290 * dz) + (cse_p11_t1_g214291 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = ((cse_p11_t1_g214287 * (cse_p11_t43_g214235 + (-(cse_p11_t1_g214288)))) + (cse_p11_t1_g214292 * dz));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p11_t1_g214295 + (cse_p11_t1_g214289 * (-(dx))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = ((cse_p11_t1_g214287 * dy * (-(dx))) + (cse_p11_t1_g214289 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = ((cse_p11_t1_g214289 * (-(dz))) + (cse_p11_t1_g214291 * (-(dx))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = (cse_p11_t1_g214293 + (cse_p11_t1_g214290 * (-(dy))));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = ((cse_p11_t1_g214290 * (-(dz))) + (cse_p11_t1_g214291 * (-(dy))));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = (cse_p11_t1_g214294 + (cse_p11_t1_g214292 * (-(dz))));
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
                double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t41_g214233))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t2_invr253575 = (1.0 / (rr));
                double cse_p565_t1_invr253574 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g214288 = (rr * rr);
                double cse_p11_t1_g214287 = (cse_p565_t1_invr253574 * de_dr);
                double cse_p11_t1_g214286 = (cse_p565_t2_invr253575 * cse_p565_t2_invr253575 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g214286 * dx * dy) + (cse_p11_t1_g214287 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double h_x1_z1 = ((cse_p11_t1_g214286 * dx * dz) + (cse_p11_t1_g214287 * dz * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g214291 = (cse_p11_t1_g214287 * dz);
                double h_x1_x2 = ((cse_p11_t1_g214286 * dx * (-(dx))) + (cse_p11_t1_g214287 * (cse_p11_t41_g214233 + (-(cse_p11_t1_g214288)))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g214289 = (cse_p11_t1_g214286 * (-(dx)));
                double h_x1_y2 = ((cse_p11_t1_g214286 * dx * (-(dy))) + (cse_p11_t1_g214287 * dx * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g214290 = (cse_p11_t1_g214286 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g214286 * dx * (-(dz))) + (cse_p11_t1_g214291 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g214292 = (cse_p11_t1_g214286 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g214286 * dy * dy) + (cse_p11_t1_g214287 * (cse_p11_t1_g214288 + (-(cse_p11_t42_g214234)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g214293 = (cse_p11_t1_g214287 * (cse_p11_t1_g214288 + (-(cse_p11_t42_g214234))));
                double h_y1_z1 = ((cse_p11_t1_g214286 * dy * dz) + (cse_p11_t1_g214291 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g214287 * dx * dy) + (cse_p11_t1_g214289 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g214287 * (cse_p11_t42_g214234 + (-(cse_p11_t1_g214288)))) + (cse_p11_t1_g214290 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g214291 * dy) + (cse_p11_t1_g214292 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g214286 * dz * dz) + (cse_p11_t1_g214287 * (cse_p11_t1_g214288 + (-(cse_p11_t43_g214235)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g214294 = (cse_p11_t1_g214287 * (cse_p11_t1_g214288 + (-(cse_p11_t43_g214235))));
                double h_z1_x2 = ((cse_p11_t1_g214289 * dz) + (cse_p11_t1_g214291 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g214290 * dz) + (cse_p11_t1_g214291 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g214287 * (cse_p11_t43_g214235 + (-(cse_p11_t1_g214288)))) + (cse_p11_t1_g214292 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g214287 * (cse_p11_t1_g214288 + (-(cse_p11_t41_g214233)))) + (cse_p11_t1_g214289 * (-(dx))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p11_t1_g214295 = (cse_p11_t1_g214287 * (cse_p11_t1_g214288 + (-(cse_p11_t41_g214233))));
                double h_x2_y2 = ((cse_p11_t1_g214287 * dy * (-(dx))) + (cse_p11_t1_g214289 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double h_x2_z2 = ((cse_p11_t1_g214289 * (-(dz))) + (cse_p11_t1_g214291 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g214293 + (cse_p11_t1_g214290 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g214290 * (-(dz))) + (cse_p11_t1_g214291 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g214294 + (cse_p11_t1_g214292 * (-(dz))));
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
