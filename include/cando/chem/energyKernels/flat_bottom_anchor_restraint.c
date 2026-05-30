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
      double cse_p1_t1_sqrt373347 = sqrt(r2);
      double cse_p1_t2_invsqrt373348 = (1.0 / (cse_p1_t1_sqrt373347));
      double rr = cse_p1_t1_sqrt373347;
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
          double cse_p11_t1_g373349 = (dr * term.ka);
          double energy = (cse_p11_t1_g373349 * dr);
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (2.00000000000000000e+0 * cse_p11_t1_g373349 * cse_p1_t2_invsqrt373348 * dx);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double cse_p528_t1_g378357 = (2.00000000000000000e+0 * cse_p11_t1_g373349 * cse_p1_t2_invsqrt373348);
          double g_y1 = (cse_p528_t1_g378357 * dy);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p528_t1_g378357 * dz);
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
      double cse_p11_t13_g378373 = (dx * dx);
      double dy = (y1 + (-(term.ya)));
      double cse_p11_t16_g378376 = (dy * dy);
      double dz = (z1 + (-(term.za)));
      double cse_p11_t18_g378378 = (dz * dz);
      double r2 = (cse_p11_t13_g378373 + cse_p11_t16_g378376 + cse_p11_t18_g378378);
      double cse_p1_t1_sqrt378359 = sqrt(r2);
      double cse_p1_t2_invsqrt378360 = (1.0 / (cse_p1_t1_sqrt378359));
      double rr = cse_p1_t1_sqrt378359;
      if ((rr <= term.r0)) {
        /* !COND-1 */
        double cse_p60_t6_g379755 = (dx * dy);
        double cse_p60_t7_g379756 = (dx * dz);
        double cse_p60_t8_g379757 = (dy * dz);
        double cse_p60_t9_g379758 = (rr * rr);
        double cse_p60_t10_g379759 = (-(cse_p11_t13_g378373));
        double cse_p60_t11_g379760 = (-(cse_p60_t6_g379755));
        double cse_p60_t12_g379761 = (-(cse_p60_t7_g379756));
        double cse_p60_t13_g379762 = (-(cse_p11_t16_g378376));
        double cse_p60_t14_g379763 = (-(cse_p60_t8_g379757));
        double cse_p60_t15_g379764 = (-(cse_p11_t18_g378378));
        double cse_p538_t1_invr388382 = (1.0 / ((rr * rr * rr)));
        double cse_p538_t2_invr388383 = (1.0 / (rr));
        double cse_p60_t18_g379767 = (rr * rr * rr);
        double cse_p14_t7_g378469 = (cse_p60_t13_g379762 + cse_p60_t9_g379758);
        double cse_p14_t8_g378470 = (cse_p60_t15_g379764 + cse_p60_t9_g379758);
        double cse_p13_t15_g378451 = (cse_p60_t10_g379759 + cse_p60_t9_g379758);
        double de_dr = 0.00000000000000000e+0;
        double cse_p60_t5_g379754 = (cse_p538_t1_invr388382 * de_dr);
        double cse_p13_t3_g378439 = (cse_p14_t7_g378469 * cse_p60_t5_g379754);
        double cse_p13_t4_g378440 = (cse_p14_t8_g378470 * cse_p60_t5_g379754);
        double cse_p13_t5_g378441 = (cse_p60_t12_g379761 * cse_p60_t5_g379754);
        double cse_p13_t6_g378442 = (cse_p60_t14_g379763 * cse_p60_t5_g379754);
        double cse_p12_t5_g378406 = (cse_p13_t15_g378451 * cse_p60_t5_g379754);
        double cse_p12_t8_g378409 = (cse_p60_t11_g379760 * cse_p60_t5_g379754);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p60_t1_g379777 = (cse_p538_t2_invr388383 * cse_p538_t2_invr388383 * d2e_dr2);
        double cse_p60_t1_g379750 = (cse_p60_t1_g379777 * dx * dz);
        double cse_p60_t2_g379751 = (cse_p60_t1_g379777 * dy);
        double cse_p60_t3_g379752 = (cse_p60_t1_g379777 * dz * dz);
        double cse_p12_t1_g378402 = (cse_p60_t1_g379777 * dx * dx);
        double cse_p13_t7_g378443 = (cse_p60_t2_g379751 * dy);
        double cse_p13_t8_g378444 = (cse_p60_t2_g379751 * dz);
        double cse_p12_t21_g378422 = (cse_p13_t5_g378441 + cse_p60_t1_g379750);
        double cse_p12_t22_g378423 = (cse_p13_t4_g378440 + cse_p60_t3_g379752);
        double cse_p12_t11_g378412 = (cse_p60_t2_g379751 * dx);
        double cse_p12_t20_g378421 = (cse_p13_t3_g378439 + cse_p13_t7_g378443);
        double cse_p12_t23_g378424 = (cse_p13_t6_g378442 + cse_p13_t8_g378444);
        double cse_p11_t20_g378380 = (cse_p12_t1_g378402 + cse_p12_t5_g378406);
        double cse_p11_t21_g378381 = (cse_p12_t11_g378412 + cse_p12_t8_g378409);
        double energy = 0.00000000000000000e+0;
        energy_added += energy;
        *energy_accumulate += energy;
        double g_x1 = 0.00000000000000000e+0;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = 0.00000000000000000e+0;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = 0.00000000000000000e+0;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double h_x1_x1 = cse_p11_t20_g378380;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p11_t21_g378381;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t21_g378422;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_y1_y1 = cse_p12_t20_g378421;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t23_g378424;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_z1_z1 = cse_p12_t22_g378423;
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
          double g_x1 = (2.00000000000000000e+0 * cse_p1_t2_invsqrt378360 * dr * dx * term.ka);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (2.00000000000000000e+0 * cse_p1_t2_invsqrt378360 * dr * dy * term.ka);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (2.00000000000000000e+0 * cse_p1_t2_invsqrt378360 * dr * dz * term.ka);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t13_g378373))) * (1.0 / ((rr * rr * rr)))));
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double cse_p538_t2_invr388383 = (1.0 / (rr));
          double cse_p538_t1_invr388382 = (1.0 / ((rr * rr * rr)));
          double cse_p60_t1_g379777 = (cse_p538_t2_invr388383 * cse_p538_t2_invr388383 * d2e_dr2);
          double cse_p11_t1_g378400 = (rr * rr);
          double cse_p11_t1_g378399 = (cse_p538_t1_invr388382 * de_dr);
          double h_x1_y1 = ((cse_p60_t1_g379777 * dx * dy) + (-((cse_p11_t1_g378399 * dx * dy))));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p11_t1_g378401 = (cse_p60_t1_g379777 * dy);
          double h_x1_z1 = ((cse_p60_t1_g379777 * dx * dz) + (-((cse_p11_t1_g378399 * dx * dz))));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double h_y1_y1 = ((cse_p11_t1_g378399 * (cse_p11_t1_g378400 + (-(cse_p11_t16_g378376)))) + (cse_p11_t1_g378401 * dy));
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = ((cse_p11_t1_g378401 * dz) + (-((cse_p11_t1_g378399 * dy * dz))));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_z1_z1 = ((cse_p11_t1_g378399 * (cse_p11_t1_g378400 + (-(cse_p11_t18_g378378)))) + (cse_p60_t1_g379777 * dz * dz));
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
