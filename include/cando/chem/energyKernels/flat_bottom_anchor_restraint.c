#include "flat_bottom_anchor_restraint.h"

template <typename HESSIAN>
struct Flat_Bottom_Anchor_Restraint {
  static constexpr size_t PositionSize = 3;
  static std::string description() { return "mathkernel-flat_bottom_anchor_restraint"; };
double energy(const flat_bottom_anchor_restraint_term& term, double* position, double* energy_accumulate) {
  double energy_added = 0.00000000000000000e+0;
  {
    /* !BASE */
    DOUBLE x1 = position[term.i3x1 + 0];
    DOUBLE y1 = position[term.i3x1 + 1];
    DOUBLE z1 = position[term.i3x1 + 2];
    {
      /* !BASE */
      double dx = (x1 + (-(term.xa)));
      double dy = (y1 + (-(term.ya)));
      double dz = (z1 + (-(term.za)));
      double r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      double rr = sqrt(r2);
      if ((rr <= term.r0)) {
        /* !COND-1 */
        double energy = 0.00000000000000000e+0;
        energy_added += energy;
        *energy_accumulate += energy;
      } else {
        /* !COND-ELSE-1 */
        {
          /* !COND-ELSE-2 */
          double dr = (rr + (-(term.r0)));
          double energy = (dr * dr * term.ka);
          energy_added += energy;
          *energy_accumulate += energy;
        }
      }
    }
  }
  return energy_added;
}
void energy_fd(const flat_bottom_anchor_restraint_term& term, double* position, double* energy_accumulate)
{
  energy(term, position, energy_accumulate);
}

double gradient(const flat_bottom_anchor_restraint_term& term, double* position, double* energy_accumulate, double* force) {
  double energy_added = 0.00000000000000000e+0;
  {
    /* !BASE */
    DOUBLE x1 = position[term.i3x1 + 0];
    DOUBLE y1 = position[term.i3x1 + 1];
    DOUBLE z1 = position[term.i3x1 + 2];
    {
      /* !BASE */
      double dx = (x1 + (-(term.xa)));
      double dy = (y1 + (-(term.ya)));
      double dz = (z1 + (-(term.za)));
      double r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      double cse_p1_t1_sqrt186603 = sqrt(r2);
      double cse_p1_t2_invsqrt186604 = (1.0 / (cse_p1_t1_sqrt186603));
      double rr = cse_p1_t1_sqrt186603;
      if ((rr <= term.r0)) {
        /* !COND-1 */
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
      } else {
        /* !COND-ELSE-1 */
        {
          /* !COND-ELSE-2 */
          double dr = (rr + (-(term.r0)));
          double de_dr = (2.00000000000000000e+0 * dr * term.ka);
          double cse_p11_t1_g186605 = (dr * term.ka);
          double energy = (cse_p11_t1_g186605 * dr);
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (2.00000000000000000e+0 * cse_p11_t1_g186605 * cse_p1_t2_invsqrt186604 * dx);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double cse_p528_t1_g191613 = (2.00000000000000000e+0 * cse_p11_t1_g186605 * cse_p1_t2_invsqrt186604);
          double g_y1 = (cse_p528_t1_g191613 * dy);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p528_t1_g191613 * dz);
          KernelGradientAcc(term.i3x1, 2, g_z1);
        }
      }
    }
  }
  return energy_added;
}
void gradient_fd(const flat_bottom_anchor_restraint_term& term, double* position, double* energy_accumulate, double* force)
{
  constexpr size_t PositionSize = 3;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  double e0 = 0.0;
  energy(term, position, &e0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[term.i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x1 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x1, 0, d);
  }
  {
    double saved = position[term.i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x1 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x1, 1, d);
  }
  {
    double saved = position[term.i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x1 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x1, 2, d);
  }
}

double hessian(const flat_bottom_anchor_restraint_term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double energy_added = 0.00000000000000000e+0;
  {
    /* !BASE */
    DOUBLE x1 = position[term.i3x1 + 0];
    DOUBLE y1 = position[term.i3x1 + 1];
    DOUBLE z1 = position[term.i3x1 + 2];
    {
      /* !BASE */
      double dx = (x1 + (-(term.xa)));
      double cse_p11_t13_g191629 = (dx * dx);
      double dy = (y1 + (-(term.ya)));
      double cse_p11_t16_g191632 = (dy * dy);
      double dz = (z1 + (-(term.za)));
      double cse_p11_t18_g191634 = (dz * dz);
      double r2 = (cse_p11_t13_g191629 + cse_p11_t16_g191632 + cse_p11_t18_g191634);
      double cse_p1_t1_sqrt191615 = sqrt(r2);
      double cse_p1_t2_invsqrt191616 = (1.0 / (cse_p1_t1_sqrt191615));
      double rr = cse_p1_t1_sqrt191615;
      if ((rr <= term.r0)) {
        /* !COND-1 */
        double cse_p60_t6_g193011 = (dx * dy);
        double cse_p60_t7_g193012 = (dx * dz);
        double cse_p60_t8_g193013 = (dy * dz);
        double cse_p60_t9_g193014 = (rr * rr);
        double cse_p60_t10_g193015 = (-(cse_p11_t13_g191629));
        double cse_p60_t11_g193016 = (-(cse_p60_t6_g193011));
        double cse_p60_t12_g193017 = (-(cse_p60_t7_g193012));
        double cse_p60_t13_g193018 = (-(cse_p11_t16_g191632));
        double cse_p60_t14_g193019 = (-(cse_p60_t8_g193013));
        double cse_p60_t15_g193020 = (-(cse_p11_t18_g191634));
        double cse_p538_t1_invr201638 = (1.0 / ((rr * rr * rr)));
        double cse_p538_t2_invr201639 = (1.0 / (rr));
        double cse_p60_t18_g193023 = (rr * rr * rr);
        double cse_p14_t7_g191725 = (cse_p60_t13_g193018 + cse_p60_t9_g193014);
        double cse_p14_t8_g191726 = (cse_p60_t15_g193020 + cse_p60_t9_g193014);
        double cse_p13_t15_g191707 = (cse_p60_t10_g193015 + cse_p60_t9_g193014);
        double de_dr = 0.00000000000000000e+0;
        double cse_p60_t5_g193010 = (cse_p538_t1_invr201638 * de_dr);
        double cse_p13_t3_g191695 = (cse_p14_t7_g191725 * cse_p60_t5_g193010);
        double cse_p13_t4_g191696 = (cse_p14_t8_g191726 * cse_p60_t5_g193010);
        double cse_p13_t5_g191697 = (cse_p60_t12_g193017 * cse_p60_t5_g193010);
        double cse_p13_t6_g191698 = (cse_p60_t14_g193019 * cse_p60_t5_g193010);
        double cse_p12_t5_g191662 = (cse_p13_t15_g191707 * cse_p60_t5_g193010);
        double cse_p12_t8_g191665 = (cse_p60_t11_g193016 * cse_p60_t5_g193010);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p60_t1_g193033 = (cse_p538_t2_invr201639 * cse_p538_t2_invr201639 * d2e_dr2);
        double cse_p60_t1_g193006 = (cse_p60_t1_g193033 * dx * dz);
        double cse_p60_t2_g193007 = (cse_p60_t1_g193033 * dy);
        double cse_p60_t3_g193008 = (cse_p60_t1_g193033 * dz * dz);
        double cse_p12_t1_g191658 = (cse_p60_t1_g193033 * dx * dx);
        double cse_p13_t7_g191699 = (cse_p60_t2_g193007 * dy);
        double cse_p13_t8_g191700 = (cse_p60_t2_g193007 * dz);
        double cse_p12_t21_g191678 = (cse_p13_t5_g191697 + cse_p60_t1_g193006);
        double cse_p12_t22_g191679 = (cse_p13_t4_g191696 + cse_p60_t3_g193008);
        double cse_p12_t11_g191668 = (cse_p60_t2_g193007 * dx);
        double cse_p12_t20_g191677 = (cse_p13_t3_g191695 + cse_p13_t7_g191699);
        double cse_p12_t23_g191680 = (cse_p13_t6_g191698 + cse_p13_t8_g191700);
        double cse_p11_t20_g191636 = (cse_p12_t1_g191658 + cse_p12_t5_g191662);
        double cse_p11_t21_g191637 = (cse_p12_t11_g191668 + cse_p12_t8_g191665);
        double energy = 0.00000000000000000e+0;
        energy_added += energy;
        *energy_accumulate += energy;
        double g_x1 = 0.00000000000000000e+0;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = 0.00000000000000000e+0;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = 0.00000000000000000e+0;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double h_x1_x1 = cse_p11_t20_g191636;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p11_t21_g191637;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t21_g191678;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_y1_y1 = cse_p12_t20_g191677;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t23_g191680;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_z1_z1 = cse_p12_t22_g191679;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
      } else {
        /* !COND-ELSE-1 */
        {
          /* !COND-ELSE-2 */
          double dr = (rr + (-(term.r0)));
          double de_dr = (2.00000000000000000e+0 * dr * term.ka);
          double d2e_dr2 = (2.00000000000000000e+0 * term.ka);
          double energy = (dr * dr * term.ka);
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (2.00000000000000000e+0 * cse_p1_t2_invsqrt191616 * dr * dx * term.ka);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (2.00000000000000000e+0 * cse_p1_t2_invsqrt191616 * dr * dy * term.ka);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (2.00000000000000000e+0 * cse_p1_t2_invsqrt191616 * dr * dz * term.ka);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t13_g191629))) * (1.0 / ((rr * rr * rr)))));
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double cse_p538_t2_invr201639 = (1.0 / (rr));
          double cse_p538_t1_invr201638 = (1.0 / ((rr * rr * rr)));
          double cse_p60_t1_g193033 = (cse_p538_t2_invr201639 * cse_p538_t2_invr201639 * d2e_dr2);
          double cse_p11_t1_g191656 = (rr * rr);
          double cse_p11_t1_g191655 = (cse_p538_t1_invr201638 * de_dr);
          double h_x1_y1 = ((cse_p60_t1_g193033 * dx * dy) + (-((cse_p11_t1_g191655 * dx * dy))));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p11_t1_g191657 = (cse_p60_t1_g193033 * dy);
          double h_x1_z1 = ((cse_p60_t1_g193033 * dx * dz) + (-((cse_p11_t1_g191655 * dx * dz))));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double h_y1_y1 = ((cse_p11_t1_g191655 * (cse_p11_t1_g191656 + (-(cse_p11_t16_g191632)))) + (cse_p11_t1_g191657 * dy));
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = ((cse_p11_t1_g191657 * dz) + (-((cse_p11_t1_g191655 * dy * dz))));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_z1_z1 = ((cse_p11_t1_g191655 * (cse_p11_t1_g191656 + (-(cse_p11_t18_g191634)))) + (cse_p60_t1_g193033 * dz * dz));
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        }
      }
    }
  }
  return energy_added;
}
void hessian_fd(const flat_bottom_anchor_restraint_term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 3;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  const double invh2 = 1.0/((h*h));
  double e0 = 0.0;
  energy(term, position, &e0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[term.i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x1 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x1, 0, d);
  }
  {
    double saved = position[term.i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x1 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x1, 1, d);
  }
  {
    double saved = position[term.i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x1 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x1, 2, d);
  }
  {
    double saved = position[term.i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x1 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x1 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 0, term.i3x1, 0, hval);
  }
  {
    double saved = position[term.i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x1 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x1 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 1, term.i3x1, 1, hval);
  }
  {
    double saved = position[term.i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x1 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x1 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x1 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 2, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x1 + 1];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x1 + 1] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x1 + 1] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
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
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x1 + 2] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
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
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x1 + 2] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x1 + 2] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 2, term.i3x1, 1, hval);
  }
}

};
