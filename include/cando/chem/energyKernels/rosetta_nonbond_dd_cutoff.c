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
      double cse_p1_t1_sqrt253581 = sqrt(r2);
      double cse_p1_t2_invsqrt253582 = (1.0 / (cse_p1_t1_sqrt253581));
      double rr = cse_p1_t1_sqrt253581;
      if ((rr <= (6.00000023841857900e-1 * term.sigma))) {
        /* !COND-1 */
        double cse_p528_t7_g267122 = (cse_p1_t2_invsqrt253582 * dx);
        double cse_p528_t8_g267123 = (cse_p1_t2_invsqrt253582 * dy);
        double cse_p528_t9_g267124 = (cse_p1_t2_invsqrt253582 * dz);
        double cse_p528_t16_g267131 = (-(term.epsilon));
        double cse_p12_t1_g253649 = (term.mij * params.rep_weight);
        double cse_p11_t9_g253593 = (cse_p12_t1_g253649 * cse_p528_t7_g267122);
        double cse_p11_t10_g253594 = (cse_p12_t1_g253649 * cse_p528_t8_g267123);
        double cse_p11_t11_g253595 = (cse_p12_t1_g253649 * cse_p528_t9_g267124);
        double fa_rep = (term.bij + (term.mij * rr));
        double cse_p528_t11_g267126 = (fa_rep * params.rep_weight);
        double fa_atr = cse_p528_t16_g267131;
        double cse_p13_t20_g253675 = (cse_p528_t11_g267126 + fa_atr);
        double de_dr = cse_p12_t1_g253649;
        double energy = cse_p13_t20_g253675;
        energy_added += energy;
        *energy_accumulate += energy;
        double g_x1 = cse_p11_t9_g253593;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = cse_p11_t10_g253594;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = cse_p11_t11_g253595;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = (-(cse_p11_t9_g253593));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (-(cse_p11_t10_g253594));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (-(cse_p11_t11_g253595));
        KernelGradientAcc(term.i3x2, 2, g_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.sigma)) {
          /* !COND-2 */
          double cse_p529_t1_invr267143 = (1.0 / (rr));
          double cse_p479_t11_g265801 = (cse_p529_t1_invr267143 * cse_p529_t1_invr267143);
          double cse_p479_t12_g265802 = (cse_p529_t1_invr267143 * term.sigma);
          double sr = cse_p479_t12_g265802;
          double cse_p528_t1_g267142 = (cse_p1_t2_invsqrt253582 * term.sigma * sr);
          double cse_p528_t12_g267127 = (sr * sr);
          double cse_p11_t1_g253613 = (cse_p479_t11_g265801 * cse_p528_t1_g267142);
          double cse_p528_t4_g267119 = (cse_p11_t1_g253613 * dx);
          double cse_p528_t5_g267120 = (cse_p11_t1_g253613 * dy);
          double cse_p528_t6_g267121 = (cse_p11_t1_g253613 * dz);
          double cse_p12_t19_g253641 = (-(cse_p528_t4_g267119));
          double cse_p12_t20_g253642 = (-(cse_p528_t5_g267120));
          double cse_p12_t21_g253643 = (-(cse_p528_t6_g267121));
          double sr2 = cse_p528_t12_g267127;
          double cse_p528_t13_g267128 = (sr2 * sr2);
          double cse_p13_t4_g253659 = (6.00000000000000000e+0 * cse_p528_t13_g267128 * term.epsilon);
          double cse_p13_t8_g253663 = (cse_p528_t13_g267128 * sr2);
          double sr6 = cse_p13_t8_g253663;
          double cse_p528_t2_g267117 = (-2.00000000000000000e+0 * sr6);
          double cse_p528_t3_g267118 = (1.20000000000000000e+1 * sr6);
          double cse_p528_t14_g267129 = (sr6 * sr6);
          double cse_p528_t15_g267130 = (-2.00000000000000000e+0 + sr6 + sr6);
          double sr12 = cse_p528_t14_g267129;
          double cse_p528_t1_g267116 = (-1.20000000000000000e+1 * sr12);
          double cse_p14_t15_g253702 = (cse_p528_t2_g267117 + sr12);
          double cse_p13_t1_g253687 = (cse_p14_t15_g253702 * term.epsilon);
          double cse_p13_t19_g253674 = (cse_p528_t1_g267116 + cse_p528_t3_g267118);
          double fa_rep_epsilon = cse_p13_t1_g253687;
          double fa_rep = (term.epsilon + fa_rep_epsilon);
          double fa_atr = (-(term.epsilon));
          double de_dr = (cse_p13_t19_g253674 * cse_p529_t1_invr267143 * term.epsilon * params.rep_weight);
          double cse_p12_t1_g253653 = (fa_rep * params.rep_weight);
          double energy = (cse_p12_t1_g253653 + fa_atr);
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p11_t1_g253616 = (cse_p13_t4_g253659 * cse_p528_t15_g267130 * params.rep_weight);
          double g_x1 = (cse_p11_t1_g253616 * cse_p12_t19_g253641);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t1_g253616 * cse_p12_t20_g253642);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t1_g253616 * cse_p12_t21_g253643);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p11_t1_g253616 * cse_p528_t4_g267119);
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p11_t1_g253616 * cse_p528_t5_g267120);
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p11_t1_g253616 * cse_p528_t6_g267121);
          KernelGradientAcc(term.i3x2, 2, g_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rswitch)) {
            /* !COND-3 */
            double cse_p529_t1_invr267143 = (1.0 / (rr));
            double cse_p1_t4_invr2253584 = (cse_p529_t1_invr267143 * cse_p529_t1_invr267143);
            double sr = (cse_p529_t1_invr267143 * term.sigma);
            double sr2 = (sr * sr);
            double cse_p12_t1_g253650 = (sr2 * sr2);
            double sr6 = (cse_p12_t1_g253650 * sr2);
            double sr12 = (sr6 * sr6);
            double cse_p12_t1_g253648 = (-2.00000000000000000e+0 * sr6);
            double cse_p13_t1_g253687 = (term.epsilon * (cse_p12_t1_g253648 + sr12));
            double fa_atr = cse_p13_t1_g253687;
            double cse_p12_t1_g253651 = (-1.20000000000000000e+1 * sr12);
            double cse_p12_t1_g253652 = (1.20000000000000000e+1 * sr6);
            double de_dr = (cse_p529_t1_invr267143 * term.epsilon * (cse_p12_t1_g253651 + cse_p12_t1_g253652));
            double energy = fa_atr;
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p528_t1_g267142 = (cse_p1_t2_invsqrt253582 * term.sigma * sr);
            double cse_p11_t1_g253613 = (cse_p1_t4_invr2253584 * cse_p528_t1_g267142);
            double cse_p11_t1_g253614 = (6.00000000000000000e+0 * cse_p12_t1_g253650 * term.epsilon);
            double cse_p11_t1_g253615 = (cse_p11_t1_g253614 * (-2.00000000000000000e+0 + sr6 + sr6));
            double cse_p11_t1_g253619 = (cse_p11_t1_g253613 * dx);
            double g_x1 = (cse_p11_t1_g253615 * (-(cse_p11_t1_g253619)));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p11_t1_g253620 = (cse_p11_t1_g253613 * dy);
            double g_y1 = (cse_p11_t1_g253615 * (-(cse_p11_t1_g253620)));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p11_t1_g253621 = (cse_p11_t1_g253613 * dz);
            double g_z1 = (cse_p11_t1_g253615 * (-(cse_p11_t1_g253621)));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (cse_p11_t1_g253615 * cse_p11_t1_g253619);
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (cse_p11_t1_g253615 * cse_p11_t1_g253620);
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (cse_p11_t1_g253615 * cse_p11_t1_g253621);
            KernelGradientAcc(term.i3x2, 2, g_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p11_t1_g253618 = (term.aa * rr * rr);
              double cse_p12_t4_g253626 = (3.00000000000000000e+0 * cse_p11_t1_g253618);
              double cse_p11_t1_g253617 = (term.bb * rr);
              double cse_p11_t21_g253605 = (term.cc + cse_p11_t1_g253617 + cse_p11_t1_g253617 + cse_p12_t4_g253626);
              double fa_atr = (term.dd + (term.cc * rr) + (rr * rr * (term.bb + (term.aa * rr))));
              double de_dr = (term.cc + cse_p12_t4_g253626 + (2.00000000000000000e+0 * cse_p11_t1_g253617));
              double energy = fa_atr;
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p11_t1_g253622 = (cse_p1_t2_invsqrt253582 * dx);
              double g_x1 = (cse_p11_t1_g253622 * cse_p11_t21_g253605);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double cse_p12_t1_g253646 = (cse_p1_t2_invsqrt253582 * dy);
              double g_y1 = (cse_p11_t21_g253605 * cse_p12_t1_g253646);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double cse_p12_t1_g253647 = (cse_p1_t2_invsqrt253582 * dz);
              double g_z1 = (cse_p11_t21_g253605 * cse_p12_t1_g253647);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p63_t1_g254978 = (-2.00000000000000000e+0 * cse_p11_t1_g253617);
              double cse_p63_t1_g254979 = (-3.00000000000000000e+0 * cse_p11_t1_g253618);
              double g_x2 = ((cse_p11_t1_g253622 * cse_p63_t1_g254978) + (cse_p11_t1_g253622 * cse_p63_t1_g254979) + (-((term.cc * cse_p11_t1_g253622))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = ((cse_p12_t1_g253646 * cse_p63_t1_g254978) + (cse_p12_t1_g253646 * cse_p63_t1_g254979) + (-((term.cc * cse_p12_t1_g253646))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = ((cse_p12_t1_g253647 * cse_p63_t1_g254978) + (cse_p12_t1_g253647 * cse_p63_t1_g254979) + (-((term.cc * cse_p12_t1_g253647))));
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
      double cse_p11_t49_g276277 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t50_g276278 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t51_g276279 = (dz * dz);
      double r2 = (cse_p11_t49_g276277 + cse_p11_t50_g276278 + cse_p11_t51_g276279);
      double cse_p1_t1_sqrt276225 = sqrt(r2);
      double cse_p1_t2_invsqrt276226 = (1.0 / (cse_p1_t1_sqrt276225));
      double rr = cse_p1_t1_sqrt276225;
      if ((rr <= (6.00000023841857900e-1 * term.sigma))) {
        /* !COND-1 */
        double cse_p564_t46_g316832 = (rr * rr);
        double cse_p564_t58_g316844 = (-(cse_p11_t49_g276277));
        double cse_p564_t59_g316845 = (-(cse_p11_t50_g276278));
        double cse_p564_t60_g316846 = (-(cse_p11_t51_g276279));
        double cse_p564_t61_g316847 = (-(dx));
        double cse_p564_t62_g316848 = (-(dy));
        double cse_p564_t63_g316849 = (-(dz));
        double cse_p564_t64_g316850 = (-(term.epsilon));
        double cse_p565_t1_invr316865 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr316866 = (1.0 / (rr));
        double cse_p564_t68_g316854 = (rr * rr * rr);
        double cse_p15_t89_g276785 = (-(cse_p564_t46_g316832));
        double cse_p14_t90_g276667 = (cse_p11_t49_g276277 + cse_p15_t89_g276785);
        double cse_p14_t92_g276669 = (cse_p11_t50_g276278 + cse_p15_t89_g276785);
        double cse_p14_t94_g276671 = (cse_p11_t51_g276279 + cse_p15_t89_g276785);
        double cse_p13_t116_g276560 = (cse_p564_t46_g316832 + cse_p564_t58_g316844);
        double cse_p13_t117_g276561 = (cse_p564_t46_g316832 + cse_p564_t59_g316845);
        double cse_p13_t118_g276562 = (cse_p564_t46_g316832 + cse_p564_t60_g316846);
        double cse_p11_t11_g276239 = (cse_p1_t2_invsqrt276226 * dx * term.mij * params.rep_weight);
        double cse_p11_t12_g276240 = (cse_p1_t2_invsqrt276226 * dy * term.mij * params.rep_weight);
        double cse_p11_t13_g276241 = (cse_p1_t2_invsqrt276226 * dz * term.mij * params.rep_weight);
        double fa_rep = (term.bij + (term.mij * rr));
        double cse_p564_t45_g316831 = (fa_rep * params.rep_weight);
        double fa_atr = cse_p564_t64_g316850;
        double cse_p60_t83_g281100 = (cse_p564_t45_g316831 + fa_atr);
        double de_dr = (term.mij * params.rep_weight);
        double cse_p515_t36_g313000 = (cse_p565_t1_invr316865 * de_dr);
        double cse_p13_t68_g276512 = (cse_p14_t90_g276667 * cse_p515_t36_g313000);
        double cse_p13_t69_g276513 = (cse_p14_t92_g276669 * cse_p515_t36_g313000);
        double cse_p13_t70_g276514 = (cse_p14_t94_g276671 * cse_p515_t36_g313000);
        double cse_p13_t71_g276515 = (cse_p515_t36_g313000 * dx * dy);
        double cse_p13_t72_g276516 = (cse_p515_t36_g313000 * cse_p564_t61_g316847 * dy);
        double cse_p12_t23_g276366 = (cse_p13_t116_g276560 * cse_p515_t36_g313000);
        double cse_p12_t24_g276367 = (cse_p13_t117_g276561 * cse_p515_t36_g313000);
        double cse_p12_t25_g276368 = (cse_p13_t118_g276562 * cse_p515_t36_g313000);
        double cse_p12_t28_g276371 = (cse_p515_t36_g313000 * dz);
        double cse_p13_t58_g276502 = (cse_p12_t28_g276371 * cse_p564_t61_g316847);
        double cse_p13_t59_g276503 = (cse_p12_t28_g276371 * cse_p564_t62_g316848);
        double cse_p13_t60_g276504 = (cse_p12_t28_g276371 * dx);
        double cse_p13_t61_g276505 = (cse_p12_t28_g276371 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p515_t37_g313001 = (cse_p565_t2_invr316866 * cse_p565_t2_invr316866 * d2e_dr2);
        double cse_p13_t62_g276506 = (cse_p515_t37_g313001 * dx * dx);
        double cse_p13_t63_g276507 = (cse_p515_t37_g313001 * dx * dy);
        double cse_p13_t64_g276508 = (cse_p515_t37_g313001 * dx * dz);
        double cse_p13_t65_g276509 = (cse_p515_t37_g313001 * dy * dy);
        double cse_p13_t66_g276510 = (cse_p515_t37_g313001 * dy * dz);
        double cse_p13_t67_g276511 = (cse_p515_t37_g313001 * dz * dz);
        double cse_p12_t11_g276354 = (cse_p515_t37_g313001 * cse_p564_t61_g316847);
        double cse_p13_t43_g276487 = (cse_p12_t11_g276354 * cse_p564_t61_g316847);
        double cse_p13_t44_g276488 = (cse_p12_t11_g276354 * cse_p564_t62_g316848);
        double cse_p13_t45_g276489 = (cse_p12_t11_g276354 * cse_p564_t63_g316849);
        double cse_p13_t46_g276490 = (cse_p12_t11_g276354 * dx);
        double cse_p13_t47_g276491 = (cse_p12_t11_g276354 * dy);
        double cse_p13_t48_g276492 = (cse_p12_t11_g276354 * dz);
        double cse_p12_t12_g276355 = (cse_p515_t37_g313001 * cse_p564_t62_g316848);
        double cse_p13_t49_g276493 = (cse_p12_t12_g276355 * cse_p564_t62_g316848);
        double cse_p13_t50_g276494 = (cse_p12_t12_g276355 * cse_p564_t63_g316849);
        double cse_p13_t51_g276495 = (cse_p12_t12_g276355 * dx);
        double cse_p13_t52_g276496 = (cse_p12_t12_g276355 * dy);
        double cse_p13_t53_g276497 = (cse_p12_t12_g276355 * dz);
        double cse_p12_t13_g276356 = (cse_p515_t37_g313001 * cse_p564_t63_g316849);
        double cse_p13_t54_g276498 = (cse_p12_t13_g276356 * cse_p564_t63_g316849);
        double cse_p13_t55_g276499 = (cse_p12_t13_g276356 * dx);
        double cse_p13_t56_g276500 = (cse_p12_t13_g276356 * dy);
        double cse_p13_t57_g276501 = (cse_p12_t13_g276356 * dz);
        double cse_p12_t63_g276406 = (cse_p13_t63_g276507 + cse_p13_t72_g276516);
        double cse_p12_t64_g276407 = (cse_p13_t58_g276502 + cse_p13_t64_g276508);
        double cse_p12_t65_g276408 = (cse_p13_t59_g276503 + cse_p13_t66_g276510);
        double cse_p12_t82_g276425 = (cse_p12_t24_g276367 + cse_p13_t65_g276509);
        double cse_p12_t84_g276427 = (cse_p12_t25_g276368 + cse_p13_t67_g276511);
        double cse_p12_t86_g276429 = (cse_p12_t23_g276366 + cse_p13_t62_g276506);
        double cse_p12_t58_g276401 = (cse_p13_t46_g276490 + cse_p13_t68_g276512);
        double cse_p12_t66_g276409 = (cse_p13_t47_g276491 + cse_p13_t71_g276515);
        double cse_p12_t68_g276411 = (cse_p13_t44_g276488 + cse_p13_t72_g276516);
        double cse_p12_t69_g276412 = (cse_p13_t45_g276489 + cse_p13_t58_g276502);
        double cse_p12_t70_g276413 = (cse_p13_t48_g276492 + cse_p13_t60_g276504);
        double cse_p12_t87_g276430 = (cse_p12_t23_g276366 + cse_p13_t43_g276487);
        double cse_p12_t59_g276402 = (cse_p13_t52_g276496 + cse_p13_t69_g276513);
        double cse_p12_t67_g276410 = (cse_p13_t51_g276495 + cse_p13_t71_g276515);
        double cse_p12_t71_g276414 = (cse_p13_t50_g276494 + cse_p13_t59_g276503);
        double cse_p12_t72_g276415 = (cse_p13_t53_g276497 + cse_p13_t61_g276505);
        double cse_p12_t83_g276426 = (cse_p12_t24_g276367 + cse_p13_t49_g276493);
        double cse_p12_t60_g276403 = (cse_p13_t57_g276501 + cse_p13_t70_g276514);
        double cse_p12_t73_g276416 = (cse_p13_t55_g276499 + cse_p13_t60_g276504);
        double cse_p12_t74_g276417 = (cse_p13_t56_g276500 + cse_p13_t61_g276505);
        double cse_p12_t85_g276428 = (cse_p12_t25_g276368 + cse_p13_t54_g276498);
        double energy = cse_p60_t83_g281100;
        energy_added += energy;
        *energy_accumulate += energy;
        double g_x1 = cse_p11_t11_g276239;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = cse_p11_t12_g276240;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = cse_p11_t13_g276241;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = (-(cse_p11_t11_g276239));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (-(cse_p11_t12_g276240));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (-(cse_p11_t13_g276241));
        KernelGradientAcc(term.i3x2, 2, g_z2);
        double h_x1_x1 = cse_p12_t86_g276429;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t63_g276406;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t64_g276407;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t58_g276401;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t67_g276410;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t73_g276416;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t82_g276425;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t65_g276408;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t66_g276409;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t59_g276402;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t74_g276417;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t84_g276427;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t70_g276413;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t72_g276415;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t60_g276403;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t87_g276430;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t68_g276411;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t69_g276412;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t83_g276426;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t71_g276414;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t85_g276428;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.sigma)) {
          /* !COND-2 */
          double cse_p11_t1_g276336 = (rr * rr);
          double cse_p564_t51_g316837 = (cse_p11_t1_g276336 + (-(cse_p11_t49_g276277)));
          double cse_p564_t52_g316838 = (cse_p11_t1_g276336 + (-(cse_p11_t50_g276278)));
          double cse_p564_t53_g316839 = (cse_p11_t1_g276336 + (-(cse_p11_t51_g276279)));
          double cse_p564_t54_g316840 = (cse_p11_t49_g276277 + (-(cse_p11_t1_g276336)));
          double cse_p564_t55_g316841 = (cse_p11_t50_g276278 + (-(cse_p11_t1_g276336)));
          double cse_p564_t56_g316842 = (cse_p11_t51_g276279 + (-(cse_p11_t1_g276336)));
          double cse_p564_t57_g316843 = (-(cse_p11_t1_g276336));
          double cse_p565_t3_invr316867 = (1.0 / (cse_p11_t1_g276336));
          double cse_p565_t2_invr316866 = (1.0 / (rr));
          double cse_p564_t38_g316824 = (cse_p565_t2_invr316866 * cse_p565_t2_invr316866);
          double cse_p564_t39_g316825 = (cse_p565_t2_invr316866 * term.sigma);
          double sr = cse_p564_t39_g316825;
          double cse_p564_t47_g316833 = (sr * sr);
          double cse_p67_t1_g281218 = (cse_p1_t2_invsqrt276226 * cse_p564_t38_g316824 * term.sigma * sr);
          double cse_p12_t8_g276351 = (cse_p67_t1_g281218 * dx);
          double cse_p12_t9_g276352 = (cse_p67_t1_g281218 * dy);
          double cse_p12_t10_g276353 = (cse_p67_t1_g281218 * dz);
          double cse_p11_t91_g276319 = (-(cse_p12_t8_g276351));
          double cse_p11_t92_g276320 = (-(cse_p12_t9_g276352));
          double cse_p11_t93_g276321 = (-(cse_p12_t10_g276353));
          double sr2 = cse_p564_t47_g316833;
          double cse_p564_t48_g316834 = (sr2 * sr2 * sr2);
          double sr6 = cse_p564_t48_g316834;
          double cse_p564_t2_g316788 = (-2.00000000000000000e+0 * sr6);
          double cse_p564_t3_g316789 = (-4.20000000000000040e+1 * sr6);
          double cse_p564_t4_g316790 = (1.20000000000000000e+1 * sr6);
          double cse_p564_t49_g316835 = (sr6 * sr6);
          double cse_p564_t50_g316836 = (-2.00000000000000000e+0 + sr6 + sr6);
          double sr12 = cse_p564_t49_g316835;
          double cse_p564_t1_g316787 = (-1.20000000000000000e+1 * sr12);
          double cse_p564_t5_g316791 = (1.56000000000000000e+2 * sr12);
          double cse_p67_t78_g281199 = (cse_p564_t2_g316788 + sr12);
          double cse_p60_t47_g281064 = (cse_p67_t78_g281199 * term.epsilon);
          double cse_p60_t56_g281073 = (cse_p564_t1_g316787 + cse_p564_t4_g316790);
          double cse_p60_t57_g281074 = (cse_p564_t3_g316789 + cse_p564_t5_g316791);
          double fa_rep_epsilon = cse_p60_t47_g281064;
          double fa_rep = (term.epsilon + fa_rep_epsilon);
          double fa_atr = (-(term.epsilon));
          double de_dr = (cse_p565_t2_invr316866 * cse_p60_t56_g281073 * term.epsilon * params.rep_weight);
          double d2e_dr2 = (cse_p565_t3_invr316867 * cse_p60_t57_g281074 * term.epsilon * params.rep_weight);
          double cse_p564_t37_g316823 = (cse_p565_t2_invr316866 * cse_p565_t2_invr316866 * d2e_dr2);
          double energy = (fa_atr + (fa_rep * params.rep_weight));
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p67_t1_g281215 = (6.00000000000000000e+0 * cse_p564_t50_g316836 * term.epsilon * params.rep_weight * sr2 * sr2);
          double g_x1 = (cse_p11_t91_g276319 * cse_p67_t1_g281215);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t92_g276320 * cse_p67_t1_g281215);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t93_g276321 * cse_p67_t1_g281215);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p67_t1_g281215 * cse_p67_t1_g281218 * dx);
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p67_t1_g281215 * cse_p67_t1_g281218 * dy);
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p67_t1_g281215 * cse_p67_t1_g281218 * dz);
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double cse_p68_t8_g281229 = (cse_p564_t37_g316823 * (-(dx)));
          double cse_p68_t9_g281230 = (cse_p564_t37_g316823 * (-(dy)));
          double cse_p68_t10_g281231 = (cse_p564_t37_g316823 * (-(dz)));
          double cse_p68_t11_g281232 = (cse_p564_t37_g316823 * dx);
          double cse_p68_t12_g281233 = (cse_p564_t37_g316823 * dy);
          double cse_p68_t13_g281234 = (cse_p564_t37_g316823 * dz * dz);
          double cse_p68_t45_g281266 = (cse_p68_t11_g281232 * dx);
          double cse_p68_t46_g281267 = (cse_p68_t11_g281232 * dy);
          double cse_p68_t47_g281268 = (cse_p68_t11_g281232 * dz);
          double cse_p68_t50_g281271 = (cse_p68_t12_g281233 * dy);
          double cse_p68_t51_g281272 = (cse_p68_t12_g281233 * dz);
          double cse_p14_t47_g276624 = (cse_p68_t8_g281229 * (-(dx)));
          double cse_p14_t48_g276625 = (cse_p68_t8_g281229 * (-(dy)));
          double cse_p14_t49_g276626 = (cse_p68_t8_g281229 * (-(dz)));
          double cse_p14_t50_g276627 = (cse_p68_t8_g281229 * dx);
          double cse_p14_t51_g276628 = (cse_p68_t8_g281229 * dy);
          double cse_p14_t52_g276629 = (cse_p68_t8_g281229 * dz);
          double cse_p14_t53_g276630 = (cse_p68_t9_g281230 * (-(dy)));
          double cse_p14_t54_g276631 = (cse_p68_t9_g281230 * (-(dz)));
          double cse_p14_t55_g276632 = (cse_p68_t9_g281230 * dx);
          double cse_p14_t56_g276633 = (cse_p68_t9_g281230 * dy);
          double cse_p14_t57_g276634 = (cse_p68_t9_g281230 * dz);
          double cse_p14_t58_g276635 = (cse_p68_t10_g281231 * (-(dz)));
          double cse_p14_t59_g276636 = (cse_p68_t10_g281231 * dx);
          double cse_p14_t60_g276637 = (cse_p68_t10_g281231 * dy);
          double cse_p14_t61_g276638 = (cse_p68_t10_g281231 * dz);
          double cse_p565_t1_invr316865 = (1.0 / ((rr * rr * rr)));
          double cse_p564_t36_g316822 = (cse_p565_t1_invr316865 * de_dr);
          double cse_p68_t20_g281241 = (cse_p564_t36_g316822 * dy);
          double cse_p68_t21_g281242 = (cse_p564_t36_g316822 * dz);
          double cse_p68_t48_g281269 = (cse_p68_t20_g281241 * (-(dx)));
          double cse_p68_t49_g281270 = (cse_p68_t20_g281241 * dx);
          double cse_p60_t15_g281032 = (cse_p564_t36_g316822 * (cse_p11_t49_g276277 + cse_p564_t57_g316843));
          double cse_p60_t16_g281033 = (cse_p564_t36_g316822 * (cse_p11_t50_g276278 + cse_p564_t57_g316843));
          double cse_p60_t17_g281034 = (cse_p564_t36_g316822 * (cse_p11_t51_g276279 + cse_p564_t57_g316843));
          double cse_p60_t18_g281035 = (cse_p564_t36_g316822 * cse_p564_t51_g316837);
          double cse_p60_t19_g281036 = (cse_p564_t36_g316822 * cse_p564_t52_g316838);
          double cse_p60_t20_g281037 = (cse_p564_t36_g316822 * cse_p564_t53_g316839);
          double cse_p60_t58_g281075 = (cse_p68_t46_g281267 + cse_p68_t48_g281269);
          double cse_p14_t43_g276620 = (cse_p68_t21_g281242 * (-(dx)));
          double cse_p14_t44_g276621 = (cse_p68_t21_g281242 * (-(dy)));
          double cse_p14_t45_g276622 = (cse_p68_t21_g281242 * dx);
          double cse_p14_t46_g276623 = (cse_p68_t21_g281242 * dy);
          double cse_p13_t114_g276558 = (cse_p60_t18_g281035 + cse_p68_t45_g281266);
          double h_x1_x1 = cse_p13_t114_g276558;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t58_g281075;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t89_g276533 = (cse_p14_t43_g276620 + cse_p68_t47_g281268);
          double cse_p13_t90_g276534 = (cse_p14_t44_g276621 + cse_p68_t51_g281272);
          double h_x1_z1 = cse_p13_t89_g276533;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t83_g276527 = (cse_p14_t50_g276627 + cse_p60_t15_g281032);
          double cse_p13_t91_g276535 = (cse_p14_t51_g276628 + cse_p68_t49_g281270);
          double cse_p13_t93_g276537 = (cse_p14_t48_g276625 + cse_p68_t48_g281269);
          double cse_p13_t94_g276538 = (cse_p14_t43_g276620 + cse_p14_t49_g276626);
          double cse_p13_t95_g276539 = (cse_p14_t45_g276622 + cse_p14_t52_g276629);
          double cse_p13_t115_g276559 = (cse_p14_t47_g276624 + cse_p60_t18_g281035);
          double h_x1_x2 = cse_p13_t83_g276527;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t84_g276528 = (cse_p14_t56_g276633 + cse_p60_t16_g281033);
          double cse_p13_t92_g276536 = (cse_p14_t55_g276632 + cse_p68_t49_g281270);
          double cse_p13_t96_g276540 = (cse_p14_t44_g276621 + cse_p14_t54_g276631);
          double cse_p13_t97_g276541 = (cse_p14_t46_g276623 + cse_p14_t57_g276634);
          double h_x1_y2 = cse_p13_t92_g276536;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t85_g276529 = (cse_p14_t61_g276638 + cse_p60_t17_g281034);
          double cse_p13_t98_g276542 = (cse_p14_t45_g276622 + cse_p14_t59_g276636);
          double cse_p13_t99_g276543 = (cse_p14_t46_g276623 + cse_p14_t60_g276637);
          double h_x1_z2 = cse_p13_t98_g276542;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t110_g276554 = (cse_p60_t19_g281036 + cse_p68_t50_g281271);
          double cse_p13_t111_g276555 = (cse_p14_t53_g276630 + cse_p60_t19_g281036);
          double h_y1_y1 = cse_p13_t110_g276554;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t90_g276534;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t91_g276535;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t84_g276528;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t99_g276543;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t112_g276556 = (cse_p60_t20_g281037 + cse_p68_t13_g281234);
          double cse_p13_t113_g276557 = (cse_p14_t58_g276635 + cse_p60_t20_g281037);
          double h_z1_z1 = cse_p13_t112_g276556;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t95_g276539;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t97_g276541;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t85_g276529;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t115_g276559;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t93_g276537;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t94_g276538;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t111_g276555;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t96_g276540;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t113_g276557;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rswitch)) {
            /* !COND-3 */
            double cse_p565_t2_invr316866 = (1.0 / (rr));
            double cse_p1_t4_invr2276228 = (cse_p565_t2_invr316866 * cse_p565_t2_invr316866);
            double sr = (cse_p565_t2_invr316866 * term.sigma);
            double sr2 = (sr * sr);
            double sr6 = (sr2 * sr2 * sr2);
            double sr12 = (sr6 * sr6);
            double fa_atr = (term.epsilon * (sr12 + (-2.00000000000000000e+0 * sr6)));
            double de_dr = (cse_p565_t2_invr316866 * term.epsilon * ((-1.20000000000000000e+1 * sr12) + (1.20000000000000000e+1 * sr6)));
            double cse_p11_t1_g276336 = (rr * rr);
            double cse_p565_t3_invr316867 = (1.0 / (cse_p11_t1_g276336));
            double d2e_dr2 = (cse_p565_t3_invr316867 * term.epsilon * ((-4.20000000000000040e+1 * sr6) + (1.56000000000000000e+2 * sr12)));
            double energy = fa_atr;
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p67_t1_g281212 = (cse_p1_t2_invsqrt276226 * cse_p1_t4_invr2276228 * term.epsilon * term.sigma * sr2 * sr2 * sr * (-2.00000000000000000e+0 + sr6 + sr6));
            double g_x1 = (-6.00000000000000000e+0 * cse_p67_t1_g281212 * dx);
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double g_y1 = (-6.00000000000000000e+0 * cse_p67_t1_g281212 * dy);
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double g_z1 = (-6.00000000000000000e+0 * cse_p67_t1_g281212 * dz);
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (6.00000000000000000e+0 * cse_p67_t1_g281212 * dx);
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (6.00000000000000000e+0 * cse_p67_t1_g281212 * dy);
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (6.00000000000000000e+0 * cse_p67_t1_g281212 * dz);
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double cse_p11_t1_g276334 = (cse_p565_t2_invr316866 * cse_p565_t2_invr316866 * d2e_dr2);
            double cse_p564_t6_g316792 = (cse_p11_t1_g276334 * (-(dx)));
            double cse_p564_t7_g316793 = (cse_p11_t1_g276334 * (-(dy)));
            double cse_p564_t8_g316794 = (cse_p11_t1_g276334 * (-(dz)));
            double cse_p564_t9_g316795 = (cse_p11_t1_g276334 * dx);
            double cse_p564_t10_g316796 = (cse_p11_t1_g276334 * dy);
            double cse_p564_t11_g316797 = (cse_p11_t1_g276334 * dz * dz);
            double cse_p70_t45_g281460 = (cse_p564_t9_g316795 * dx);
            double cse_p70_t46_g281461 = (cse_p564_t9_g316795 * dy);
            double cse_p565_t1_invr316865 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g276335 = (cse_p565_t1_invr316865 * de_dr);
            double cse_p564_t1_g316864 = (cse_p11_t1_g276335 * (cse_p11_t1_g276336 + (-(cse_p11_t49_g276277))));
            double cse_p564_t13_g316799 = (cse_p11_t1_g276335 * (cse_p11_t1_g276336 + (-(cse_p11_t50_g276278))));
            double cse_p564_t14_g316800 = (cse_p11_t1_g276335 * (cse_p11_t1_g276336 + (-(cse_p11_t51_g276279))));
            double cse_p564_t15_g316801 = (cse_p11_t1_g276335 * (cse_p11_t49_g276277 + (-(cse_p11_t1_g276336))));
            double cse_p564_t16_g316802 = (cse_p11_t1_g276335 * (cse_p11_t50_g276278 + (-(cse_p11_t1_g276336))));
            double cse_p564_t17_g316803 = (cse_p11_t1_g276335 * (cse_p11_t51_g276279 + (-(cse_p11_t1_g276336))));
            double cse_p564_t18_g316804 = (cse_p11_t1_g276335 * dy);
            double cse_p564_t19_g316805 = (cse_p11_t1_g276335 * dz);
            double cse_p69_t45_g281364 = (cse_p564_t9_g316795 * dz);
            double cse_p69_t74_g281393 = (cse_p564_t1_g316864 + cse_p70_t45_g281460);
            double cse_p69_t46_g281365 = (cse_p564_t18_g316804 * (-(dx)));
            double cse_p69_t47_g281366 = (cse_p564_t18_g316804 * dx);
            double cse_p69_t73_g281392 = (cse_p69_t46_g281365 + cse_p70_t46_g281461);
            double h_x1_x1 = cse_p69_t74_g281393;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = cse_p69_t73_g281392;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p69_t31_g281350 = (cse_p564_t19_g316805 * (-(dx)));
            double cse_p69_t32_g281351 = (cse_p564_t19_g316805 * (-(dy)));
            double cse_p69_t33_g281352 = (cse_p564_t19_g316805 * dx);
            double cse_p69_t34_g281353 = (cse_p564_t19_g316805 * dy);
            double cse_p60_t59_g281076 = (cse_p69_t31_g281350 + cse_p69_t45_g281364);
            double cse_p69_t48_g281367 = (cse_p564_t10_g316796 * dy);
            double cse_p69_t49_g281368 = (cse_p564_t10_g316796 * dz);
            double cse_p60_t60_g281077 = (cse_p69_t32_g281351 + cse_p69_t49_g281368);
            double h_x1_z1 = cse_p60_t59_g281076;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p69_t20_g281339 = (cse_p564_t6_g316792 * (-(dx)));
            double cse_p69_t21_g281340 = (cse_p564_t6_g316792 * (-(dy)));
            double cse_p69_t22_g281341 = (cse_p564_t6_g316792 * (-(dz)));
            double cse_p69_t23_g281342 = (cse_p564_t6_g316792 * dx);
            double cse_p69_t24_g281343 = (cse_p564_t6_g316792 * dy);
            double cse_p69_t25_g281344 = (cse_p564_t6_g316792 * dz);
            double cse_p60_t53_g281070 = (cse_p564_t15_g316801 + cse_p69_t23_g281342);
            double cse_p60_t61_g281078 = (cse_p69_t24_g281343 + cse_p69_t47_g281366);
            double cse_p60_t63_g281080 = (cse_p69_t21_g281340 + cse_p69_t46_g281365);
            double cse_p60_t64_g281081 = (cse_p69_t22_g281341 + cse_p69_t31_g281350);
            double cse_p60_t65_g281082 = (cse_p69_t25_g281344 + cse_p69_t33_g281352);
            double cse_p60_t82_g281099 = (cse_p564_t1_g316864 + cse_p69_t20_g281339);
            double h_x1_x2 = cse_p60_t53_g281070;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p69_t26_g281345 = (cse_p564_t7_g316793 * (-(dy)));
            double cse_p69_t27_g281346 = (cse_p564_t7_g316793 * (-(dz)));
            double cse_p69_t28_g281347 = (cse_p564_t7_g316793 * dx);
            double cse_p69_t29_g281348 = (cse_p564_t7_g316793 * dy);
            double cse_p69_t30_g281349 = (cse_p564_t7_g316793 * dz);
            double cse_p60_t54_g281071 = (cse_p564_t16_g316802 + cse_p69_t29_g281348);
            double cse_p60_t62_g281079 = (cse_p69_t28_g281347 + cse_p69_t47_g281366);
            double cse_p60_t66_g281083 = (cse_p69_t27_g281346 + cse_p69_t32_g281351);
            double cse_p60_t67_g281084 = (cse_p69_t30_g281349 + cse_p69_t34_g281353);
            double h_x1_y2 = cse_p60_t62_g281079;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p69_t35_g281354 = (cse_p564_t8_g316794 * (-(dz)));
            double cse_p69_t36_g281355 = (cse_p564_t8_g316794 * dx);
            double cse_p69_t37_g281356 = (cse_p564_t8_g316794 * dy);
            double cse_p69_t38_g281357 = (cse_p564_t8_g316794 * dz);
            double cse_p60_t55_g281072 = (cse_p564_t17_g316803 + cse_p69_t38_g281357);
            double cse_p60_t68_g281085 = (cse_p69_t33_g281352 + cse_p69_t36_g281355);
            double cse_p60_t69_g281086 = (cse_p69_t34_g281353 + cse_p69_t37_g281356);
            double h_x1_z2 = cse_p60_t68_g281085;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p67_t70_g281191 = (cse_p564_t13_g316799 + cse_p69_t48_g281367);
            double cse_p60_t78_g281095 = (cse_p564_t13_g316799 + cse_p69_t26_g281345);
            double h_y1_y1 = cse_p67_t70_g281191;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t60_g281077;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t61_g281078;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t54_g281071;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t69_g281086;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p67_t72_g281193 = (cse_p564_t11_g316797 + cse_p564_t14_g316800);
            double cse_p60_t80_g281097 = (cse_p564_t14_g316800 + cse_p69_t35_g281354);
            double h_z1_z1 = cse_p67_t72_g281193;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t65_g281082;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t67_g281084;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t55_g281072;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t82_g281099;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t63_g281080;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t64_g281081;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t78_g281095;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t66_g281083;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t80_g281097;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p68_t2_g281223 = (-2.00000000000000000e+0 * term.bb * rr);
              double cse_p12_t7_g276350 = (term.bb * rr);
              double cse_p11_t1_g276336 = (rr * rr);
              double cse_p68_t4_g281225 = (-3.00000000000000000e+0 * term.aa * cse_p11_t1_g276336);
              double cse_p67_t56_g281177 = (cse_p68_t2_g281223 + cse_p68_t4_g281225);
              double cse_p12_t6_g276349 = (3.00000000000000000e+0 * term.aa * cse_p11_t1_g276336);
              double cse_p11_t88_g276316 = (term.cc + cse_p12_t6_g276349 + cse_p12_t7_g276350 + cse_p12_t7_g276350);
              double fa_atr = (term.dd + (term.cc * rr) + (cse_p11_t1_g276336 * (term.bb + (term.aa * rr))));
              double de_dr = (term.cc + cse_p12_t6_g276349 + (2.00000000000000000e+0 * term.bb * rr));
              double d2e_dr2 = ((2.00000000000000000e+0 * term.bb) + (6.00000000000000000e+0 * term.aa * rr));
              double energy = fa_atr;
              energy_added += energy;
              *energy_accumulate += energy;
              double g_x1 = (cse_p11_t88_g276316 * cse_p1_t2_invsqrt276226 * dx);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t88_g276316 * cse_p1_t2_invsqrt276226 * dy);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t88_g276316 * cse_p1_t2_invsqrt276226 * dz);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = ((cse_p1_t2_invsqrt276226 * cse_p67_t56_g281177 * dx) + (-((term.cc * cse_p1_t2_invsqrt276226 * dx))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = ((cse_p1_t2_invsqrt276226 * cse_p67_t56_g281177 * dy) + (-((term.cc * cse_p1_t2_invsqrt276226 * dy))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = ((cse_p1_t2_invsqrt276226 * cse_p67_t56_g281177 * dz) + (-((term.cc * cse_p1_t2_invsqrt276226 * dz))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double cse_p565_t2_invr316866 = (1.0 / (rr));
              double cse_p11_t1_g276334 = (cse_p565_t2_invr316866 * cse_p565_t2_invr316866 * d2e_dr2);
              double cse_p565_t1_invr316865 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g276335 = (cse_p565_t1_invr316865 * de_dr);
              double cse_p564_t1_g316864 = (cse_p11_t1_g276335 * (cse_p11_t1_g276336 + (-(cse_p11_t49_g276277))));
              double cse_p67_t1_g281213 = (cse_p11_t1_g276334 * dx);
              double cse_p564_t40_g316826 = (cse_p67_t1_g281213 * dz);
              double h_x1_x1 = (cse_p564_t1_g316864 + (cse_p67_t1_g281213 * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double cse_p67_t1_g281214 = (cse_p11_t1_g276335 * dy);
              double cse_p564_t41_g316827 = (cse_p67_t1_g281214 * (-(dx)));
              double cse_p564_t42_g316828 = (cse_p67_t1_g281214 * dx);
              double h_x1_y1 = (cse_p564_t41_g316827 + (cse_p67_t1_g281213 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g276339 = (cse_p11_t1_g276335 * dz);
              double cse_p564_t29_g316815 = (cse_p11_t1_g276339 * (-(dx)));
              double cse_p564_t30_g316816 = (cse_p11_t1_g276339 * (-(dy)));
              double cse_p564_t31_g316817 = (cse_p11_t1_g276339 * dx);
              double cse_p564_t32_g316818 = (cse_p11_t1_g276339 * dy);
              double cse_p67_t1_g281216 = (cse_p11_t1_g276334 * dy);
              double cse_p564_t43_g316829 = (cse_p67_t1_g281216 * dy);
              double cse_p564_t44_g316830 = (cse_p67_t1_g281216 * dz);
              double cse_p116_t57_g285568 = (cse_p564_t30_g316816 + cse_p564_t44_g316830);
              double h_x1_z1 = (cse_p564_t29_g316815 + cse_p564_t40_g316826);
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g276337 = (cse_p11_t1_g276334 * (-(dx)));
              double cse_p564_t20_g316806 = (cse_p11_t1_g276337 * (-(dx)));
              double cse_p564_t21_g316807 = (cse_p11_t1_g276337 * (-(dy)));
              double cse_p564_t22_g316808 = (cse_p11_t1_g276337 * (-(dz)));
              double cse_p564_t23_g316809 = (cse_p11_t1_g276337 * dy);
              double cse_p564_t24_g316810 = (cse_p11_t1_g276337 * dz);
              double cse_p116_t52_g285563 = (cse_p564_t22_g316808 + cse_p564_t29_g316815);
              double cse_p116_t53_g285564 = (cse_p564_t23_g316809 + cse_p564_t42_g316828);
              double cse_p116_t54_g285565 = (cse_p564_t24_g316810 + cse_p564_t31_g316817);
              double cse_p67_t60_g281181 = (cse_p564_t21_g316807 + cse_p564_t41_g316827);
              double h_x1_x2 = ((cse_p11_t1_g276335 * (cse_p11_t49_g276277 + (-(cse_p11_t1_g276336)))) + (cse_p11_t1_g276337 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g276338 = (cse_p11_t1_g276334 * (-(dy)));
              double cse_p564_t25_g316811 = (cse_p11_t1_g276338 * (-(dy)));
              double cse_p564_t26_g316812 = (cse_p11_t1_g276338 * (-(dz)));
              double cse_p564_t27_g316813 = (cse_p11_t1_g276338 * dy);
              double cse_p564_t28_g316814 = (cse_p11_t1_g276338 * dz);
              double cse_p116_t50_g285561 = (cse_p564_t27_g316813 + (cse_p11_t1_g276335 * (cse_p11_t50_g276278 + (-(cse_p11_t1_g276336)))));
              double cse_p116_t55_g285566 = (cse_p564_t26_g316812 + cse_p564_t30_g316816);
              double cse_p116_t56_g285567 = (cse_p564_t28_g316814 + cse_p564_t32_g316818);
              double h_x1_y2 = (cse_p564_t42_g316828 + (cse_p11_t1_g276338 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g276340 = (cse_p11_t1_g276334 * (-(dz)));
              double cse_p564_t33_g316819 = (cse_p11_t1_g276340 * (-(dz)));
              double cse_p564_t34_g316820 = (cse_p11_t1_g276340 * dy);
              double cse_p564_t35_g316821 = (cse_p11_t1_g276340 * dz);
              double cse_p116_t51_g285562 = (cse_p564_t35_g316821 + (cse_p11_t1_g276335 * (cse_p11_t51_g276279 + (-(cse_p11_t1_g276336)))));
              double cse_p116_t58_g285569 = (cse_p564_t32_g316818 + cse_p564_t34_g316820);
              double h_x1_z2 = (cse_p564_t31_g316817 + (cse_p11_t1_g276340 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g276341 = (cse_p11_t1_g276335 * (cse_p11_t1_g276336 + (-(cse_p11_t50_g276278))));
              double cse_p116_t63_g285574 = (cse_p11_t1_g276341 + cse_p564_t25_g316811);
              double h_y1_y1 = (cse_p11_t1_g276341 + cse_p564_t43_g316829);
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = cse_p116_t57_g285568;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = cse_p116_t53_g285564;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = cse_p116_t50_g285561;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = cse_p116_t58_g285569;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g276342 = (cse_p11_t1_g276335 * (cse_p11_t1_g276336 + (-(cse_p11_t51_g276279))));
              double cse_p116_t64_g285575 = (cse_p11_t1_g276342 + cse_p564_t33_g316819);
              double h_z1_z1 = (cse_p11_t1_g276342 + (cse_p11_t1_g276334 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = cse_p116_t54_g285565;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = cse_p116_t56_g285567;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = cse_p116_t51_g285562;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p564_t1_g316864 + cse_p564_t20_g316806);
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = cse_p67_t60_g281181;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = cse_p116_t52_g285563;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = cse_p116_t63_g285574;
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = cse_p116_t55_g285566;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = cse_p116_t64_g285575;
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
                double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t49_g276277))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t2_invr316866 = (1.0 / (rr));
                double cse_p565_t1_invr316865 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g276336 = (rr * rr);
                double cse_p11_t1_g276335 = (cse_p565_t1_invr316865 * de_dr);
                double cse_p11_t1_g276334 = (cse_p565_t2_invr316866 * cse_p565_t2_invr316866 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g276334 * dx * dy) + (cse_p11_t1_g276335 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double cse_p67_t1_g281216 = (cse_p11_t1_g276334 * dy);
                double cse_p67_t1_g281214 = (cse_p11_t1_g276335 * dy);
                double cse_p67_t1_g281213 = (cse_p11_t1_g276334 * dx);
                double h_x1_z1 = ((cse_p11_t1_g276335 * dz * (-(dx))) + (cse_p67_t1_g281213 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g276339 = (cse_p11_t1_g276335 * dz);
                double h_x1_x2 = ((cse_p11_t1_g276335 * (cse_p11_t49_g276277 + (-(cse_p11_t1_g276336)))) + (cse_p67_t1_g281213 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g276337 = (cse_p11_t1_g276334 * (-(dx)));
                double h_x1_y2 = ((cse_p67_t1_g281213 * (-(dy))) + (cse_p67_t1_g281214 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g276338 = (cse_p11_t1_g276334 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g276339 * dx) + (cse_p67_t1_g281213 * (-(dz))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g276340 = (cse_p11_t1_g276334 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g276335 * (cse_p11_t1_g276336 + (-(cse_p11_t50_g276278)))) + (cse_p67_t1_g281216 * dy));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g276341 = (cse_p11_t1_g276335 * (cse_p11_t1_g276336 + (-(cse_p11_t50_g276278))));
                double h_y1_z1 = ((cse_p11_t1_g276339 * (-(dy))) + (cse_p67_t1_g281216 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g276337 * dy) + (cse_p67_t1_g281214 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g276335 * (cse_p11_t50_g276278 + (-(cse_p11_t1_g276336)))) + (cse_p11_t1_g276338 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g276339 * dy) + (cse_p11_t1_g276340 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g276334 * dz * dz) + (cse_p11_t1_g276335 * (cse_p11_t1_g276336 + (-(cse_p11_t51_g276279)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g276342 = (cse_p11_t1_g276335 * (cse_p11_t1_g276336 + (-(cse_p11_t51_g276279))));
                double h_z1_x2 = ((cse_p11_t1_g276337 * dz) + (cse_p11_t1_g276339 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g276338 * dz) + (cse_p11_t1_g276339 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g276335 * (cse_p11_t51_g276279 + (-(cse_p11_t1_g276336)))) + (cse_p11_t1_g276340 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g276335 * (cse_p11_t1_g276336 + (-(cse_p11_t49_g276277)))) + (cse_p11_t1_g276337 * (-(dx))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p564_t1_g316864 = (cse_p11_t1_g276335 * (cse_p11_t1_g276336 + (-(cse_p11_t49_g276277))));
                double h_x2_y2 = ((cse_p11_t1_g276337 * (-(dy))) + (cse_p67_t1_g281214 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double cse_p67_t1_g281217 = (cse_p67_t1_g281214 * (-(dx)));
                double h_x2_z2 = ((cse_p11_t1_g276337 * (-(dz))) + (cse_p11_t1_g276339 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g276341 + (cse_p11_t1_g276338 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g276338 * (-(dz))) + (cse_p11_t1_g276339 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g276342 + (cse_p11_t1_g276340 * (-(dz))));
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
