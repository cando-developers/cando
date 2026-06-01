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
      double rr = sqrt(r2);
      if ((rr <= term.r0)) {
        /* !COND-1 */
        double de_dr = 0.00000000000000000e+0;
        double energy = 0.00000000000000000e+0;
        energy_added += energy;
        *energy_accumulate += energy;
        double cse_p538_t1_invr5568 = (1.0 / (rr));
        double g_x1 = (cse_p538_t1_invr5568 * de_dr * dx);
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = (cse_p538_t1_invr5568 * de_dr * dy);
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = (cse_p538_t1_invr5568 * de_dr * dz);
        KernelGradientAcc(term.i3x1, 2, g_z1);
      } else {
        /* !COND-ELSE-1 */
        {
          /* !COND-ELSE-2 */
          double dr = (rr + (-(term.r0)));
          double de_dr = (2.00000000000000000e+0 * dr * term.ka);
          double energy = (dr * dr * term.ka);
          energy_added += energy;
          *energy_accumulate += energy;
          double g_x1 = (de_dr * dx * (1.0 / (rr)));
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double cse_p538_t1_invr5568 = (1.0 / (rr));
          double g_y1 = (cse_p538_t1_invr5568 * de_dr * dy);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p538_t1_invr5568 * de_dr * dz);
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
      double cse_p11_t15_g25 = (dx * dx);
      double dy = (y1 + (-(term.ya)));
      double cse_p11_t18_g28 = (dy * dy);
      double dz = (z1 + (-(term.za)));
      double cse_p11_t20_g30 = (dz * dz);
      double r2 = (cse_p11_t15_g25 + cse_p11_t18_g28 + cse_p11_t20_g30);
      double rr = sqrt(r2);
      if ((rr <= term.r0)) {
        /* !COND-1 */
        double cse_p60_t6_g1407 = (dx * dy);
        double cse_p60_t7_g1408 = (dx * dz);
        double cse_p60_t8_g1409 = (dy * dz);
        double cse_p60_t9_g1410 = (rr * rr);
        double cse_p60_t10_g1411 = (-(cse_p11_t15_g25));
        double cse_p60_t11_g1412 = (-(cse_p60_t6_g1407));
        double cse_p60_t12_g1413 = (-(cse_p60_t7_g1408));
        double cse_p60_t13_g1414 = (-(cse_p11_t18_g28));
        double cse_p60_t14_g1415 = (-(cse_p60_t8_g1409));
        double cse_p60_t15_g1416 = (-(cse_p11_t20_g30));
        double cse_p538_t1_invr11384 = (1.0 / ((rr * rr * rr)));
        double cse_p538_t2_invr11385 = (1.0 / (rr));
        double cse_p60_t18_g1419 = (rr * rr * rr);
        double cse_p14_t7_g121 = (cse_p60_t13_g1414 + cse_p60_t9_g1410);
        double cse_p14_t8_g122 = (cse_p60_t15_g1416 + cse_p60_t9_g1410);
        double cse_p13_t15_g103 = (cse_p60_t10_g1411 + cse_p60_t9_g1410);
        double de_dr = 0.00000000000000000e+0;
        double cse_p60_t5_g1406 = (cse_p538_t1_invr11384 * de_dr);
        double cse_p13_t3_g91 = (cse_p14_t7_g121 * cse_p60_t5_g1406);
        double cse_p13_t4_g92 = (cse_p14_t8_g122 * cse_p60_t5_g1406);
        double cse_p13_t5_g93 = (cse_p60_t12_g1413 * cse_p60_t5_g1406);
        double cse_p13_t6_g94 = (cse_p60_t14_g1415 * cse_p60_t5_g1406);
        double cse_p12_t5_g58 = (cse_p13_t15_g103 * cse_p60_t5_g1406);
        double cse_p12_t8_g61 = (cse_p60_t11_g1412 * cse_p60_t5_g1406);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p60_t1_g1429 = (cse_p538_t2_invr11385 * cse_p538_t2_invr11385 * d2e_dr2);
        double cse_p60_t1_g1402 = (cse_p60_t1_g1429 * dx * dz);
        double cse_p60_t2_g1403 = (cse_p60_t1_g1429 * dy);
        double cse_p60_t3_g1404 = (cse_p60_t1_g1429 * dz * dz);
        double cse_p12_t1_g54 = (cse_p60_t1_g1429 * dx * dx);
        double cse_p13_t7_g95 = (cse_p60_t2_g1403 * dy);
        double cse_p13_t8_g96 = (cse_p60_t2_g1403 * dz);
        double cse_p12_t21_g74 = (cse_p13_t5_g93 + cse_p60_t1_g1402);
        double cse_p12_t22_g75 = (cse_p13_t4_g92 + cse_p60_t3_g1404);
        double cse_p12_t11_g64 = (cse_p60_t2_g1403 * dx);
        double cse_p12_t20_g73 = (cse_p13_t3_g91 + cse_p13_t7_g95);
        double cse_p12_t23_g76 = (cse_p13_t6_g94 + cse_p13_t8_g96);
        double cse_p11_t22_g32 = (cse_p12_t1_g54 + cse_p12_t5_g58);
        double cse_p11_t23_g33 = (cse_p12_t11_g64 + cse_p12_t8_g61);
        double energy = 0.00000000000000000e+0;
        energy_added += energy;
        *energy_accumulate += energy;
        double cse_p11_t1_g11 = (cse_p538_t2_invr11385 * de_dr * dy);
        double cse_p11_t2_g12 = (cse_p538_t2_invr11385 * de_dr * dz);
        double g_x1 = (cse_p538_t2_invr11385 * de_dr * dx);
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = cse_p11_t1_g11;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = cse_p11_t2_g12;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double h_x1_x1 = cse_p11_t22_g32;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p11_t23_g33;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t21_g74;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_y1_y1 = cse_p12_t20_g73;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t23_g76;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_z1_z1 = cse_p12_t22_g75;
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
          double g_x1 = (de_dr * dx * (1.0 / (rr)));
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double cse_p538_t2_invr11385 = (1.0 / (rr));
          double g_y1 = (cse_p538_t2_invr11385 * de_dr * dy);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p538_t2_invr11385 * de_dr * dz);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double h_x1_x1 = ((cse_p538_t2_invr11385 * cse_p538_t2_invr11385 * d2e_dr2 * dx * dx) + (de_dr * ((rr * rr) + (-(cse_p11_t15_g25))) * (1.0 / ((rr * rr * rr)))));
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double cse_p538_t1_invr11384 = (1.0 / ((rr * rr * rr)));
          double cse_p60_t1_g1429 = (cse_p538_t2_invr11385 * cse_p538_t2_invr11385 * d2e_dr2);
          double cse_p11_t1_g52 = (rr * rr);
          double cse_p11_t1_g51 = (cse_p538_t1_invr11384 * de_dr);
          double h_x1_y1 = ((cse_p60_t1_g1429 * dx * dy) + (-((cse_p11_t1_g51 * dx * dy))));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p11_t1_g53 = (cse_p60_t1_g1429 * dy);
          double h_x1_z1 = ((cse_p60_t1_g1429 * dx * dz) + (-((cse_p11_t1_g51 * dx * dz))));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double h_y1_y1 = ((cse_p11_t1_g51 * (cse_p11_t1_g52 + (-(cse_p11_t18_g28)))) + (cse_p11_t1_g53 * dy));
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = ((cse_p11_t1_g53 * dz) + (-((cse_p11_t1_g51 * dy * dz))));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_z1_z1 = ((cse_p11_t1_g51 * (cse_p11_t1_g52 + (-(cse_p11_t20_g30)))) + (cse_p60_t1_g1429 * dz * dz));
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
