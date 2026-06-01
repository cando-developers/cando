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
      double cse_p1_t1_sqrt1 = sqrt(r2);
      double cse_p1_t2_invsqrt2 = (1.0 / (cse_p1_t1_sqrt1));
      double rr = cse_p1_t1_sqrt1;
      if ((rr <= (6.00000023841857900e-1 * term.sigma))) {
        /* !COND-1 */
        double cse_p528_t7_g13542 = (cse_p1_t2_invsqrt2 * dx);
        double cse_p528_t8_g13543 = (cse_p1_t2_invsqrt2 * dy);
        double cse_p528_t9_g13544 = (cse_p1_t2_invsqrt2 * dz);
        double cse_p528_t16_g13551 = (-(term.epsilon));
        double cse_p12_t1_g69 = (term.mij * params.rep_weight);
        double cse_p11_t9_g13 = (cse_p12_t1_g69 * cse_p528_t7_g13542);
        double cse_p11_t10_g14 = (cse_p12_t1_g69 * cse_p528_t8_g13543);
        double cse_p11_t11_g15 = (cse_p12_t1_g69 * cse_p528_t9_g13544);
        double fa_rep = (term.bij + (term.mij * rr));
        double cse_p528_t11_g13546 = (fa_rep * params.rep_weight);
        double fa_atr = cse_p528_t16_g13551;
        double cse_p13_t20_g95 = (cse_p528_t11_g13546 + fa_atr);
        double de_dr = cse_p12_t1_g69;
        double energy = cse_p13_t20_g95;
        energy_added += energy;
        *energy_accumulate += energy;
        double g_x1 = cse_p11_t9_g13;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = cse_p11_t10_g14;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = cse_p11_t11_g15;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = (-(cse_p11_t9_g13));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (-(cse_p11_t10_g14));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (-(cse_p11_t11_g15));
        KernelGradientAcc(term.i3x2, 2, g_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.sigma)) {
          /* !COND-2 */
          double cse_p529_t1_invr13563 = (1.0 / (rr));
          double cse_p479_t11_g12221 = (cse_p529_t1_invr13563 * cse_p529_t1_invr13563);
          double cse_p479_t12_g12222 = (cse_p529_t1_invr13563 * term.sigma);
          double sr = cse_p479_t12_g12222;
          double cse_p528_t1_g13562 = (cse_p1_t2_invsqrt2 * term.sigma * sr);
          double cse_p528_t12_g13547 = (sr * sr);
          double cse_p11_t1_g33 = (cse_p479_t11_g12221 * cse_p528_t1_g13562);
          double cse_p528_t4_g13539 = (cse_p11_t1_g33 * dx);
          double cse_p528_t5_g13540 = (cse_p11_t1_g33 * dy);
          double cse_p528_t6_g13541 = (cse_p11_t1_g33 * dz);
          double cse_p12_t19_g61 = (-(cse_p528_t4_g13539));
          double cse_p12_t20_g62 = (-(cse_p528_t5_g13540));
          double cse_p12_t21_g63 = (-(cse_p528_t6_g13541));
          double sr2 = cse_p528_t12_g13547;
          double cse_p528_t13_g13548 = (sr2 * sr2);
          double cse_p13_t4_g79 = (6.00000000000000000e+0 * cse_p528_t13_g13548 * term.epsilon);
          double cse_p13_t8_g83 = (cse_p528_t13_g13548 * sr2);
          double sr6 = cse_p13_t8_g83;
          double cse_p528_t2_g13537 = (-2.00000000000000000e+0 * sr6);
          double cse_p528_t3_g13538 = (1.20000000000000000e+1 * sr6);
          double cse_p528_t14_g13549 = (sr6 * sr6);
          double cse_p528_t15_g13550 = (-2.00000000000000000e+0 + sr6 + sr6);
          double sr12 = cse_p528_t14_g13549;
          double cse_p528_t1_g13536 = (-1.20000000000000000e+1 * sr12);
          double cse_p14_t15_g122 = (cse_p528_t2_g13537 + sr12);
          double cse_p13_t1_g107 = (cse_p14_t15_g122 * term.epsilon);
          double cse_p13_t19_g94 = (cse_p528_t1_g13536 + cse_p528_t3_g13538);
          double fa_rep_epsilon = cse_p13_t1_g107;
          double fa_rep = (term.epsilon + fa_rep_epsilon);
          double fa_atr = (-(term.epsilon));
          double de_dr = (cse_p13_t19_g94 * cse_p529_t1_invr13563 * term.epsilon * params.rep_weight);
          double cse_p12_t1_g73 = (fa_rep * params.rep_weight);
          double energy = (cse_p12_t1_g73 + fa_atr);
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p11_t1_g36 = (cse_p13_t4_g79 * cse_p528_t15_g13550 * params.rep_weight);
          double g_x1 = (cse_p11_t1_g36 * cse_p12_t19_g61);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t1_g36 * cse_p12_t20_g62);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t1_g36 * cse_p12_t21_g63);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p11_t1_g36 * cse_p528_t4_g13539);
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p11_t1_g36 * cse_p528_t5_g13540);
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p11_t1_g36 * cse_p528_t6_g13541);
          KernelGradientAcc(term.i3x2, 2, g_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rswitch)) {
            /* !COND-3 */
            double cse_p529_t1_invr13563 = (1.0 / (rr));
            double cse_p1_t4_invr24 = (cse_p529_t1_invr13563 * cse_p529_t1_invr13563);
            double sr = (cse_p529_t1_invr13563 * term.sigma);
            double sr2 = (sr * sr);
            double cse_p12_t1_g70 = (sr2 * sr2);
            double sr6 = (cse_p12_t1_g70 * sr2);
            double sr12 = (sr6 * sr6);
            double cse_p12_t1_g68 = (-2.00000000000000000e+0 * sr6);
            double cse_p13_t1_g107 = (term.epsilon * (cse_p12_t1_g68 + sr12));
            double fa_atr = cse_p13_t1_g107;
            double cse_p12_t1_g71 = (-1.20000000000000000e+1 * sr12);
            double cse_p12_t1_g72 = (1.20000000000000000e+1 * sr6);
            double de_dr = (cse_p529_t1_invr13563 * term.epsilon * (cse_p12_t1_g71 + cse_p12_t1_g72));
            double energy = fa_atr;
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p528_t1_g13562 = (cse_p1_t2_invsqrt2 * term.sigma * sr);
            double cse_p11_t1_g33 = (cse_p1_t4_invr24 * cse_p528_t1_g13562);
            double cse_p11_t1_g34 = (6.00000000000000000e+0 * cse_p12_t1_g70 * term.epsilon);
            double cse_p11_t1_g35 = (cse_p11_t1_g34 * (-2.00000000000000000e+0 + sr6 + sr6));
            double cse_p11_t1_g39 = (cse_p11_t1_g33 * dx);
            double g_x1 = (cse_p11_t1_g35 * (-(cse_p11_t1_g39)));
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double cse_p11_t1_g40 = (cse_p11_t1_g33 * dy);
            double g_y1 = (cse_p11_t1_g35 * (-(cse_p11_t1_g40)));
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double cse_p11_t1_g41 = (cse_p11_t1_g33 * dz);
            double g_z1 = (cse_p11_t1_g35 * (-(cse_p11_t1_g41)));
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (cse_p11_t1_g35 * cse_p11_t1_g39);
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (cse_p11_t1_g35 * cse_p11_t1_g40);
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (cse_p11_t1_g35 * cse_p11_t1_g41);
            KernelGradientAcc(term.i3x2, 2, g_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p11_t1_g38 = (term.aa * rr * rr);
              double cse_p12_t4_g46 = (3.00000000000000000e+0 * cse_p11_t1_g38);
              double cse_p11_t1_g37 = (term.bb * rr);
              double cse_p11_t21_g25 = (term.cc + cse_p11_t1_g37 + cse_p11_t1_g37 + cse_p12_t4_g46);
              double fa_atr = (term.dd + (term.cc * rr) + (rr * rr * (term.bb + (term.aa * rr))));
              double de_dr = (term.cc + cse_p12_t4_g46 + (2.00000000000000000e+0 * cse_p11_t1_g37));
              double energy = fa_atr;
              energy_added += energy;
              *energy_accumulate += energy;
              double cse_p11_t1_g42 = (cse_p1_t2_invsqrt2 * dx);
              double g_x1 = (cse_p11_t1_g42 * cse_p11_t21_g25);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double cse_p12_t1_g66 = (cse_p1_t2_invsqrt2 * dy);
              double g_y1 = (cse_p11_t21_g25 * cse_p12_t1_g66);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double cse_p12_t1_g67 = (cse_p1_t2_invsqrt2 * dz);
              double g_z1 = (cse_p11_t21_g25 * cse_p12_t1_g67);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double cse_p63_t1_g1398 = (-2.00000000000000000e+0 * cse_p11_t1_g37);
              double cse_p63_t1_g1399 = (-3.00000000000000000e+0 * cse_p11_t1_g38);
              double g_x2 = ((cse_p11_t1_g42 * cse_p63_t1_g1398) + (cse_p11_t1_g42 * cse_p63_t1_g1399) + (-((term.cc * cse_p11_t1_g42))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = ((cse_p12_t1_g66 * cse_p63_t1_g1398) + (cse_p12_t1_g66 * cse_p63_t1_g1399) + (-((term.cc * cse_p12_t1_g66))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = ((cse_p12_t1_g67 * cse_p63_t1_g1398) + (cse_p12_t1_g67 * cse_p63_t1_g1399) + (-((term.cc * cse_p12_t1_g67))));
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
      double cse_p11_t49_g53 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t50_g54 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t51_g55 = (dz * dz);
      double r2 = (cse_p11_t49_g53 + cse_p11_t50_g54 + cse_p11_t51_g55);
      double cse_p1_t1_sqrt1 = sqrt(r2);
      double cse_p1_t2_invsqrt2 = (1.0 / (cse_p1_t1_sqrt1));
      double rr = cse_p1_t1_sqrt1;
      if ((rr <= (6.00000023841857900e-1 * term.sigma))) {
        /* !COND-1 */
        double cse_p564_t46_g40608 = (rr * rr);
        double cse_p564_t58_g40620 = (-(cse_p11_t49_g53));
        double cse_p564_t59_g40621 = (-(cse_p11_t50_g54));
        double cse_p564_t60_g40622 = (-(cse_p11_t51_g55));
        double cse_p564_t61_g40623 = (-(dx));
        double cse_p564_t62_g40624 = (-(dy));
        double cse_p564_t63_g40625 = (-(dz));
        double cse_p564_t64_g40626 = (-(term.epsilon));
        double cse_p565_t1_invr40641 = (1.0 / ((rr * rr * rr)));
        double cse_p565_t2_invr40642 = (1.0 / (rr));
        double cse_p564_t68_g40630 = (rr * rr * rr);
        double cse_p15_t89_g561 = (-(cse_p564_t46_g40608));
        double cse_p14_t90_g443 = (cse_p11_t49_g53 + cse_p15_t89_g561);
        double cse_p14_t92_g445 = (cse_p11_t50_g54 + cse_p15_t89_g561);
        double cse_p14_t94_g447 = (cse_p11_t51_g55 + cse_p15_t89_g561);
        double cse_p13_t116_g336 = (cse_p564_t46_g40608 + cse_p564_t58_g40620);
        double cse_p13_t117_g337 = (cse_p564_t46_g40608 + cse_p564_t59_g40621);
        double cse_p13_t118_g338 = (cse_p564_t46_g40608 + cse_p564_t60_g40622);
        double cse_p11_t11_g15 = (cse_p1_t2_invsqrt2 * dx * term.mij * params.rep_weight);
        double cse_p11_t12_g16 = (cse_p1_t2_invsqrt2 * dy * term.mij * params.rep_weight);
        double cse_p11_t13_g17 = (cse_p1_t2_invsqrt2 * dz * term.mij * params.rep_weight);
        double fa_rep = (term.bij + (term.mij * rr));
        double cse_p564_t45_g40607 = (fa_rep * params.rep_weight);
        double fa_atr = cse_p564_t64_g40626;
        double cse_p60_t83_g4876 = (cse_p564_t45_g40607 + fa_atr);
        double de_dr = (term.mij * params.rep_weight);
        double cse_p515_t36_g36776 = (cse_p565_t1_invr40641 * de_dr);
        double cse_p13_t68_g288 = (cse_p14_t90_g443 * cse_p515_t36_g36776);
        double cse_p13_t69_g289 = (cse_p14_t92_g445 * cse_p515_t36_g36776);
        double cse_p13_t70_g290 = (cse_p14_t94_g447 * cse_p515_t36_g36776);
        double cse_p13_t71_g291 = (cse_p515_t36_g36776 * dx * dy);
        double cse_p13_t72_g292 = (cse_p515_t36_g36776 * cse_p564_t61_g40623 * dy);
        double cse_p12_t23_g142 = (cse_p13_t116_g336 * cse_p515_t36_g36776);
        double cse_p12_t24_g143 = (cse_p13_t117_g337 * cse_p515_t36_g36776);
        double cse_p12_t25_g144 = (cse_p13_t118_g338 * cse_p515_t36_g36776);
        double cse_p12_t28_g147 = (cse_p515_t36_g36776 * dz);
        double cse_p13_t58_g278 = (cse_p12_t28_g147 * cse_p564_t61_g40623);
        double cse_p13_t59_g279 = (cse_p12_t28_g147 * cse_p564_t62_g40624);
        double cse_p13_t60_g280 = (cse_p12_t28_g147 * dx);
        double cse_p13_t61_g281 = (cse_p12_t28_g147 * dy);
        double d2e_dr2 = 0.00000000000000000e+0;
        double cse_p515_t37_g36777 = (cse_p565_t2_invr40642 * cse_p565_t2_invr40642 * d2e_dr2);
        double cse_p13_t62_g282 = (cse_p515_t37_g36777 * dx * dx);
        double cse_p13_t63_g283 = (cse_p515_t37_g36777 * dx * dy);
        double cse_p13_t64_g284 = (cse_p515_t37_g36777 * dx * dz);
        double cse_p13_t65_g285 = (cse_p515_t37_g36777 * dy * dy);
        double cse_p13_t66_g286 = (cse_p515_t37_g36777 * dy * dz);
        double cse_p13_t67_g287 = (cse_p515_t37_g36777 * dz * dz);
        double cse_p12_t11_g130 = (cse_p515_t37_g36777 * cse_p564_t61_g40623);
        double cse_p13_t43_g263 = (cse_p12_t11_g130 * cse_p564_t61_g40623);
        double cse_p13_t44_g264 = (cse_p12_t11_g130 * cse_p564_t62_g40624);
        double cse_p13_t45_g265 = (cse_p12_t11_g130 * cse_p564_t63_g40625);
        double cse_p13_t46_g266 = (cse_p12_t11_g130 * dx);
        double cse_p13_t47_g267 = (cse_p12_t11_g130 * dy);
        double cse_p13_t48_g268 = (cse_p12_t11_g130 * dz);
        double cse_p12_t12_g131 = (cse_p515_t37_g36777 * cse_p564_t62_g40624);
        double cse_p13_t49_g269 = (cse_p12_t12_g131 * cse_p564_t62_g40624);
        double cse_p13_t50_g270 = (cse_p12_t12_g131 * cse_p564_t63_g40625);
        double cse_p13_t51_g271 = (cse_p12_t12_g131 * dx);
        double cse_p13_t52_g272 = (cse_p12_t12_g131 * dy);
        double cse_p13_t53_g273 = (cse_p12_t12_g131 * dz);
        double cse_p12_t13_g132 = (cse_p515_t37_g36777 * cse_p564_t63_g40625);
        double cse_p13_t54_g274 = (cse_p12_t13_g132 * cse_p564_t63_g40625);
        double cse_p13_t55_g275 = (cse_p12_t13_g132 * dx);
        double cse_p13_t56_g276 = (cse_p12_t13_g132 * dy);
        double cse_p13_t57_g277 = (cse_p12_t13_g132 * dz);
        double cse_p12_t63_g182 = (cse_p13_t63_g283 + cse_p13_t72_g292);
        double cse_p12_t64_g183 = (cse_p13_t58_g278 + cse_p13_t64_g284);
        double cse_p12_t65_g184 = (cse_p13_t59_g279 + cse_p13_t66_g286);
        double cse_p12_t82_g201 = (cse_p12_t24_g143 + cse_p13_t65_g285);
        double cse_p12_t84_g203 = (cse_p12_t25_g144 + cse_p13_t67_g287);
        double cse_p12_t86_g205 = (cse_p12_t23_g142 + cse_p13_t62_g282);
        double cse_p12_t58_g177 = (cse_p13_t46_g266 + cse_p13_t68_g288);
        double cse_p12_t66_g185 = (cse_p13_t47_g267 + cse_p13_t71_g291);
        double cse_p12_t68_g187 = (cse_p13_t44_g264 + cse_p13_t72_g292);
        double cse_p12_t69_g188 = (cse_p13_t45_g265 + cse_p13_t58_g278);
        double cse_p12_t70_g189 = (cse_p13_t48_g268 + cse_p13_t60_g280);
        double cse_p12_t87_g206 = (cse_p12_t23_g142 + cse_p13_t43_g263);
        double cse_p12_t59_g178 = (cse_p13_t52_g272 + cse_p13_t69_g289);
        double cse_p12_t67_g186 = (cse_p13_t51_g271 + cse_p13_t71_g291);
        double cse_p12_t71_g190 = (cse_p13_t50_g270 + cse_p13_t59_g279);
        double cse_p12_t72_g191 = (cse_p13_t53_g273 + cse_p13_t61_g281);
        double cse_p12_t83_g202 = (cse_p12_t24_g143 + cse_p13_t49_g269);
        double cse_p12_t60_g179 = (cse_p13_t57_g277 + cse_p13_t70_g290);
        double cse_p12_t73_g192 = (cse_p13_t55_g275 + cse_p13_t60_g280);
        double cse_p12_t74_g193 = (cse_p13_t56_g276 + cse_p13_t61_g281);
        double cse_p12_t85_g204 = (cse_p12_t25_g144 + cse_p13_t54_g274);
        double energy = cse_p60_t83_g4876;
        energy_added += energy;
        *energy_accumulate += energy;
        double g_x1 = cse_p11_t11_g15;
        KernelGradientAcc(term.i3x1, 0, g_x1);
        double g_y1 = cse_p11_t12_g16;
        KernelGradientAcc(term.i3x1, 1, g_y1);
        double g_z1 = cse_p11_t13_g17;
        KernelGradientAcc(term.i3x1, 2, g_z1);
        double g_x2 = (-(cse_p11_t11_g15));
        KernelGradientAcc(term.i3x2, 0, g_x2);
        double g_y2 = (-(cse_p11_t12_g16));
        KernelGradientAcc(term.i3x2, 1, g_y2);
        double g_z2 = (-(cse_p11_t13_g17));
        KernelGradientAcc(term.i3x2, 2, g_z2);
        double h_x1_x1 = cse_p12_t86_g205;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
        double h_x1_y1 = cse_p12_t63_g182;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
        double h_x1_z1 = cse_p12_t64_g183;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
        double h_x1_x2 = cse_p12_t58_g177;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
        double h_x1_y2 = cse_p12_t67_g186;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
        double h_x1_z2 = cse_p12_t73_g192;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
        double h_y1_y1 = cse_p12_t82_g201;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
        double h_y1_z1 = cse_p12_t65_g184;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
        double h_y1_x2 = cse_p12_t66_g185;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
        double h_y1_y2 = cse_p12_t59_g178;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
        double h_y1_z2 = cse_p12_t74_g193;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
        double h_z1_z1 = cse_p12_t84_g203;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
        double h_z1_x2 = cse_p12_t70_g189;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
        double h_z1_y2 = cse_p12_t72_g191;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
        double h_z1_z2 = cse_p12_t60_g179;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
        double h_x2_x2 = cse_p12_t87_g206;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
        double h_x2_y2 = cse_p12_t68_g187;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
        double h_x2_z2 = cse_p12_t69_g188;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
        double h_y2_y2 = cse_p12_t83_g202;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
        double h_y2_z2 = cse_p12_t71_g190;
        KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
        double h_z2_z2 = cse_p12_t85_g204;
        KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      } else {
        /* !COND-ELSE-1 */
        if ((rr <= term.sigma)) {
          /* !COND-2 */
          double cse_p11_t1_g112 = (rr * rr);
          double cse_p564_t51_g40613 = (cse_p11_t1_g112 + (-(cse_p11_t49_g53)));
          double cse_p564_t52_g40614 = (cse_p11_t1_g112 + (-(cse_p11_t50_g54)));
          double cse_p564_t53_g40615 = (cse_p11_t1_g112 + (-(cse_p11_t51_g55)));
          double cse_p564_t54_g40616 = (cse_p11_t49_g53 + (-(cse_p11_t1_g112)));
          double cse_p564_t55_g40617 = (cse_p11_t50_g54 + (-(cse_p11_t1_g112)));
          double cse_p564_t56_g40618 = (cse_p11_t51_g55 + (-(cse_p11_t1_g112)));
          double cse_p564_t57_g40619 = (-(cse_p11_t1_g112));
          double cse_p565_t3_invr40643 = (1.0 / (cse_p11_t1_g112));
          double cse_p565_t2_invr40642 = (1.0 / (rr));
          double cse_p564_t38_g40600 = (cse_p565_t2_invr40642 * cse_p565_t2_invr40642);
          double cse_p564_t39_g40601 = (cse_p565_t2_invr40642 * term.sigma);
          double sr = cse_p564_t39_g40601;
          double cse_p564_t47_g40609 = (sr * sr);
          double cse_p67_t1_g4994 = (cse_p1_t2_invsqrt2 * cse_p564_t38_g40600 * term.sigma * sr);
          double cse_p12_t8_g127 = (cse_p67_t1_g4994 * dx);
          double cse_p12_t9_g128 = (cse_p67_t1_g4994 * dy);
          double cse_p12_t10_g129 = (cse_p67_t1_g4994 * dz);
          double cse_p11_t91_g95 = (-(cse_p12_t8_g127));
          double cse_p11_t92_g96 = (-(cse_p12_t9_g128));
          double cse_p11_t93_g97 = (-(cse_p12_t10_g129));
          double sr2 = cse_p564_t47_g40609;
          double cse_p564_t48_g40610 = (sr2 * sr2 * sr2);
          double sr6 = cse_p564_t48_g40610;
          double cse_p564_t2_g40564 = (-2.00000000000000000e+0 * sr6);
          double cse_p564_t3_g40565 = (-4.20000000000000040e+1 * sr6);
          double cse_p564_t4_g40566 = (1.20000000000000000e+1 * sr6);
          double cse_p564_t49_g40611 = (sr6 * sr6);
          double cse_p564_t50_g40612 = (-2.00000000000000000e+0 + sr6 + sr6);
          double sr12 = cse_p564_t49_g40611;
          double cse_p564_t1_g40563 = (-1.20000000000000000e+1 * sr12);
          double cse_p564_t5_g40567 = (1.56000000000000000e+2 * sr12);
          double cse_p67_t78_g4975 = (cse_p564_t2_g40564 + sr12);
          double cse_p60_t47_g4840 = (cse_p67_t78_g4975 * term.epsilon);
          double cse_p60_t56_g4849 = (cse_p564_t1_g40563 + cse_p564_t4_g40566);
          double cse_p60_t57_g4850 = (cse_p564_t3_g40565 + cse_p564_t5_g40567);
          double fa_rep_epsilon = cse_p60_t47_g4840;
          double fa_rep = (term.epsilon + fa_rep_epsilon);
          double fa_atr = (-(term.epsilon));
          double de_dr = (cse_p565_t2_invr40642 * cse_p60_t56_g4849 * term.epsilon * params.rep_weight);
          double d2e_dr2 = (cse_p565_t3_invr40643 * cse_p60_t57_g4850 * term.epsilon * params.rep_weight);
          double cse_p564_t37_g40599 = (cse_p565_t2_invr40642 * cse_p565_t2_invr40642 * d2e_dr2);
          double energy = (fa_atr + (fa_rep * params.rep_weight));
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p67_t1_g4991 = (6.00000000000000000e+0 * cse_p564_t50_g40612 * term.epsilon * params.rep_weight * sr2 * sr2);
          double g_x1 = (cse_p11_t91_g95 * cse_p67_t1_g4991);
          KernelGradientAcc(term.i3x1, 0, g_x1);
          double g_y1 = (cse_p11_t92_g96 * cse_p67_t1_g4991);
          KernelGradientAcc(term.i3x1, 1, g_y1);
          double g_z1 = (cse_p11_t93_g97 * cse_p67_t1_g4991);
          KernelGradientAcc(term.i3x1, 2, g_z1);
          double g_x2 = (cse_p67_t1_g4991 * cse_p67_t1_g4994 * dx);
          KernelGradientAcc(term.i3x2, 0, g_x2);
          double g_y2 = (cse_p67_t1_g4991 * cse_p67_t1_g4994 * dy);
          KernelGradientAcc(term.i3x2, 1, g_y2);
          double g_z2 = (cse_p67_t1_g4991 * cse_p67_t1_g4994 * dz);
          KernelGradientAcc(term.i3x2, 2, g_z2);
          double cse_p68_t8_g5005 = (cse_p564_t37_g40599 * (-(dx)));
          double cse_p68_t9_g5006 = (cse_p564_t37_g40599 * (-(dy)));
          double cse_p68_t10_g5007 = (cse_p564_t37_g40599 * (-(dz)));
          double cse_p68_t11_g5008 = (cse_p564_t37_g40599 * dx);
          double cse_p68_t12_g5009 = (cse_p564_t37_g40599 * dy);
          double cse_p68_t13_g5010 = (cse_p564_t37_g40599 * dz * dz);
          double cse_p68_t45_g5042 = (cse_p68_t11_g5008 * dx);
          double cse_p68_t46_g5043 = (cse_p68_t11_g5008 * dy);
          double cse_p68_t47_g5044 = (cse_p68_t11_g5008 * dz);
          double cse_p68_t50_g5047 = (cse_p68_t12_g5009 * dy);
          double cse_p68_t51_g5048 = (cse_p68_t12_g5009 * dz);
          double cse_p14_t47_g400 = (cse_p68_t8_g5005 * (-(dx)));
          double cse_p14_t48_g401 = (cse_p68_t8_g5005 * (-(dy)));
          double cse_p14_t49_g402 = (cse_p68_t8_g5005 * (-(dz)));
          double cse_p14_t50_g403 = (cse_p68_t8_g5005 * dx);
          double cse_p14_t51_g404 = (cse_p68_t8_g5005 * dy);
          double cse_p14_t52_g405 = (cse_p68_t8_g5005 * dz);
          double cse_p14_t53_g406 = (cse_p68_t9_g5006 * (-(dy)));
          double cse_p14_t54_g407 = (cse_p68_t9_g5006 * (-(dz)));
          double cse_p14_t55_g408 = (cse_p68_t9_g5006 * dx);
          double cse_p14_t56_g409 = (cse_p68_t9_g5006 * dy);
          double cse_p14_t57_g410 = (cse_p68_t9_g5006 * dz);
          double cse_p14_t58_g411 = (cse_p68_t10_g5007 * (-(dz)));
          double cse_p14_t59_g412 = (cse_p68_t10_g5007 * dx);
          double cse_p14_t60_g413 = (cse_p68_t10_g5007 * dy);
          double cse_p14_t61_g414 = (cse_p68_t10_g5007 * dz);
          double cse_p565_t1_invr40641 = (1.0 / ((rr * rr * rr)));
          double cse_p564_t36_g40598 = (cse_p565_t1_invr40641 * de_dr);
          double cse_p68_t20_g5017 = (cse_p564_t36_g40598 * dy);
          double cse_p68_t21_g5018 = (cse_p564_t36_g40598 * dz);
          double cse_p68_t48_g5045 = (cse_p68_t20_g5017 * (-(dx)));
          double cse_p68_t49_g5046 = (cse_p68_t20_g5017 * dx);
          double cse_p60_t15_g4808 = (cse_p564_t36_g40598 * (cse_p11_t49_g53 + cse_p564_t57_g40619));
          double cse_p60_t16_g4809 = (cse_p564_t36_g40598 * (cse_p11_t50_g54 + cse_p564_t57_g40619));
          double cse_p60_t17_g4810 = (cse_p564_t36_g40598 * (cse_p11_t51_g55 + cse_p564_t57_g40619));
          double cse_p60_t18_g4811 = (cse_p564_t36_g40598 * cse_p564_t51_g40613);
          double cse_p60_t19_g4812 = (cse_p564_t36_g40598 * cse_p564_t52_g40614);
          double cse_p60_t20_g4813 = (cse_p564_t36_g40598 * cse_p564_t53_g40615);
          double cse_p60_t58_g4851 = (cse_p68_t46_g5043 + cse_p68_t48_g5045);
          double cse_p14_t43_g396 = (cse_p68_t21_g5018 * (-(dx)));
          double cse_p14_t44_g397 = (cse_p68_t21_g5018 * (-(dy)));
          double cse_p14_t45_g398 = (cse_p68_t21_g5018 * dx);
          double cse_p14_t46_g399 = (cse_p68_t21_g5018 * dy);
          double cse_p13_t114_g334 = (cse_p60_t18_g4811 + cse_p68_t45_g5042);
          double h_x1_x1 = cse_p13_t114_g334;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p60_t58_g4851;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
          double cse_p13_t89_g309 = (cse_p14_t43_g396 + cse_p68_t47_g5044);
          double cse_p13_t90_g310 = (cse_p14_t44_g397 + cse_p68_t51_g5048);
          double h_x1_z1 = cse_p13_t89_g309;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
          double cse_p13_t83_g303 = (cse_p14_t50_g403 + cse_p60_t15_g4808);
          double cse_p13_t91_g311 = (cse_p14_t51_g404 + cse_p68_t49_g5046);
          double cse_p13_t93_g313 = (cse_p14_t48_g401 + cse_p68_t48_g5045);
          double cse_p13_t94_g314 = (cse_p14_t43_g396 + cse_p14_t49_g402);
          double cse_p13_t95_g315 = (cse_p14_t45_g398 + cse_p14_t52_g405);
          double cse_p13_t115_g335 = (cse_p14_t47_g400 + cse_p60_t18_g4811);
          double h_x1_x2 = cse_p13_t83_g303;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
          double cse_p13_t84_g304 = (cse_p14_t56_g409 + cse_p60_t16_g4809);
          double cse_p13_t92_g312 = (cse_p14_t55_g408 + cse_p68_t49_g5046);
          double cse_p13_t96_g316 = (cse_p14_t44_g397 + cse_p14_t54_g407);
          double cse_p13_t97_g317 = (cse_p14_t46_g399 + cse_p14_t57_g410);
          double h_x1_y2 = cse_p13_t92_g312;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
          double cse_p13_t85_g305 = (cse_p14_t61_g414 + cse_p60_t17_g4810);
          double cse_p13_t98_g318 = (cse_p14_t45_g398 + cse_p14_t59_g412);
          double cse_p13_t99_g319 = (cse_p14_t46_g399 + cse_p14_t60_g413);
          double h_x1_z2 = cse_p13_t98_g318;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
          double cse_p13_t110_g330 = (cse_p60_t19_g4812 + cse_p68_t50_g5047);
          double cse_p13_t111_g331 = (cse_p14_t53_g406 + cse_p60_t19_g4812);
          double h_y1_y1 = cse_p13_t110_g330;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p13_t90_g310;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p13_t91_g311;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p13_t84_g304;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p13_t99_g319;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
          double cse_p13_t112_g332 = (cse_p60_t20_g4813 + cse_p68_t13_g5010);
          double cse_p13_t113_g333 = (cse_p14_t58_g411 + cse_p60_t20_g4813);
          double h_z1_z1 = cse_p13_t112_g332;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p13_t95_g315;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p13_t97_g317;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p13_t85_g305;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p13_t115_g335;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p13_t93_g313;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p13_t94_g314;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p13_t111_g331;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p13_t96_g316;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p13_t113_g333;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
        } else {
          /* !COND-ELSE-2 */
          if ((rr <= params.rswitch)) {
            /* !COND-3 */
            double cse_p565_t2_invr40642 = (1.0 / (rr));
            double cse_p1_t4_invr24 = (cse_p565_t2_invr40642 * cse_p565_t2_invr40642);
            double sr = (cse_p565_t2_invr40642 * term.sigma);
            double sr2 = (sr * sr);
            double sr6 = (sr2 * sr2 * sr2);
            double sr12 = (sr6 * sr6);
            double fa_atr = (term.epsilon * (sr12 + (-2.00000000000000000e+0 * sr6)));
            double de_dr = (cse_p565_t2_invr40642 * term.epsilon * ((-1.20000000000000000e+1 * sr12) + (1.20000000000000000e+1 * sr6)));
            double cse_p11_t1_g112 = (rr * rr);
            double cse_p565_t3_invr40643 = (1.0 / (cse_p11_t1_g112));
            double d2e_dr2 = (cse_p565_t3_invr40643 * term.epsilon * ((-4.20000000000000040e+1 * sr6) + (1.56000000000000000e+2 * sr12)));
            double energy = fa_atr;
            energy_added += energy;
            *energy_accumulate += energy;
            double cse_p67_t1_g4988 = (cse_p1_t2_invsqrt2 * cse_p1_t4_invr24 * term.epsilon * term.sigma * sr2 * sr2 * sr * (-2.00000000000000000e+0 + sr6 + sr6));
            double g_x1 = (-6.00000000000000000e+0 * cse_p67_t1_g4988 * dx);
            KernelGradientAcc(term.i3x1, 0, g_x1);
            double g_y1 = (-6.00000000000000000e+0 * cse_p67_t1_g4988 * dy);
            KernelGradientAcc(term.i3x1, 1, g_y1);
            double g_z1 = (-6.00000000000000000e+0 * cse_p67_t1_g4988 * dz);
            KernelGradientAcc(term.i3x1, 2, g_z1);
            double g_x2 = (6.00000000000000000e+0 * cse_p67_t1_g4988 * dx);
            KernelGradientAcc(term.i3x2, 0, g_x2);
            double g_y2 = (6.00000000000000000e+0 * cse_p67_t1_g4988 * dy);
            KernelGradientAcc(term.i3x2, 1, g_y2);
            double g_z2 = (6.00000000000000000e+0 * cse_p67_t1_g4988 * dz);
            KernelGradientAcc(term.i3x2, 2, g_z2);
            double cse_p11_t1_g110 = (cse_p565_t2_invr40642 * cse_p565_t2_invr40642 * d2e_dr2);
            double cse_p564_t6_g40568 = (cse_p11_t1_g110 * (-(dx)));
            double cse_p564_t7_g40569 = (cse_p11_t1_g110 * (-(dy)));
            double cse_p564_t8_g40570 = (cse_p11_t1_g110 * (-(dz)));
            double cse_p564_t9_g40571 = (cse_p11_t1_g110 * dx);
            double cse_p564_t10_g40572 = (cse_p11_t1_g110 * dy);
            double cse_p564_t11_g40573 = (cse_p11_t1_g110 * dz * dz);
            double cse_p70_t45_g5236 = (cse_p564_t9_g40571 * dx);
            double cse_p70_t46_g5237 = (cse_p564_t9_g40571 * dy);
            double cse_p565_t1_invr40641 = (1.0 / ((rr * rr * rr)));
            double cse_p11_t1_g111 = (cse_p565_t1_invr40641 * de_dr);
            double cse_p564_t1_g40640 = (cse_p11_t1_g111 * (cse_p11_t1_g112 + (-(cse_p11_t49_g53))));
            double cse_p564_t13_g40575 = (cse_p11_t1_g111 * (cse_p11_t1_g112 + (-(cse_p11_t50_g54))));
            double cse_p564_t14_g40576 = (cse_p11_t1_g111 * (cse_p11_t1_g112 + (-(cse_p11_t51_g55))));
            double cse_p564_t15_g40577 = (cse_p11_t1_g111 * (cse_p11_t49_g53 + (-(cse_p11_t1_g112))));
            double cse_p564_t16_g40578 = (cse_p11_t1_g111 * (cse_p11_t50_g54 + (-(cse_p11_t1_g112))));
            double cse_p564_t17_g40579 = (cse_p11_t1_g111 * (cse_p11_t51_g55 + (-(cse_p11_t1_g112))));
            double cse_p564_t18_g40580 = (cse_p11_t1_g111 * dy);
            double cse_p564_t19_g40581 = (cse_p11_t1_g111 * dz);
            double cse_p69_t45_g5140 = (cse_p564_t9_g40571 * dz);
            double cse_p69_t74_g5169 = (cse_p564_t1_g40640 + cse_p70_t45_g5236);
            double cse_p69_t46_g5141 = (cse_p564_t18_g40580 * (-(dx)));
            double cse_p69_t47_g5142 = (cse_p564_t18_g40580 * dx);
            double cse_p69_t73_g5168 = (cse_p69_t46_g5141 + cse_p70_t46_g5237);
            double h_x1_x1 = cse_p69_t74_g5169;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
            double h_x1_y1 = cse_p69_t73_g5168;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
            double cse_p69_t31_g5126 = (cse_p564_t19_g40581 * (-(dx)));
            double cse_p69_t32_g5127 = (cse_p564_t19_g40581 * (-(dy)));
            double cse_p69_t33_g5128 = (cse_p564_t19_g40581 * dx);
            double cse_p69_t34_g5129 = (cse_p564_t19_g40581 * dy);
            double cse_p60_t59_g4852 = (cse_p69_t31_g5126 + cse_p69_t45_g5140);
            double cse_p69_t48_g5143 = (cse_p564_t10_g40572 * dy);
            double cse_p69_t49_g5144 = (cse_p564_t10_g40572 * dz);
            double cse_p60_t60_g4853 = (cse_p69_t32_g5127 + cse_p69_t49_g5144);
            double h_x1_z1 = cse_p60_t59_g4852;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
            double cse_p69_t20_g5115 = (cse_p564_t6_g40568 * (-(dx)));
            double cse_p69_t21_g5116 = (cse_p564_t6_g40568 * (-(dy)));
            double cse_p69_t22_g5117 = (cse_p564_t6_g40568 * (-(dz)));
            double cse_p69_t23_g5118 = (cse_p564_t6_g40568 * dx);
            double cse_p69_t24_g5119 = (cse_p564_t6_g40568 * dy);
            double cse_p69_t25_g5120 = (cse_p564_t6_g40568 * dz);
            double cse_p60_t53_g4846 = (cse_p564_t15_g40577 + cse_p69_t23_g5118);
            double cse_p60_t61_g4854 = (cse_p69_t24_g5119 + cse_p69_t47_g5142);
            double cse_p60_t63_g4856 = (cse_p69_t21_g5116 + cse_p69_t46_g5141);
            double cse_p60_t64_g4857 = (cse_p69_t22_g5117 + cse_p69_t31_g5126);
            double cse_p60_t65_g4858 = (cse_p69_t25_g5120 + cse_p69_t33_g5128);
            double cse_p60_t82_g4875 = (cse_p564_t1_g40640 + cse_p69_t20_g5115);
            double h_x1_x2 = cse_p60_t53_g4846;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
            double cse_p69_t26_g5121 = (cse_p564_t7_g40569 * (-(dy)));
            double cse_p69_t27_g5122 = (cse_p564_t7_g40569 * (-(dz)));
            double cse_p69_t28_g5123 = (cse_p564_t7_g40569 * dx);
            double cse_p69_t29_g5124 = (cse_p564_t7_g40569 * dy);
            double cse_p69_t30_g5125 = (cse_p564_t7_g40569 * dz);
            double cse_p60_t54_g4847 = (cse_p564_t16_g40578 + cse_p69_t29_g5124);
            double cse_p60_t62_g4855 = (cse_p69_t28_g5123 + cse_p69_t47_g5142);
            double cse_p60_t66_g4859 = (cse_p69_t27_g5122 + cse_p69_t32_g5127);
            double cse_p60_t67_g4860 = (cse_p69_t30_g5125 + cse_p69_t34_g5129);
            double h_x1_y2 = cse_p60_t62_g4855;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
            double cse_p69_t35_g5130 = (cse_p564_t8_g40570 * (-(dz)));
            double cse_p69_t36_g5131 = (cse_p564_t8_g40570 * dx);
            double cse_p69_t37_g5132 = (cse_p564_t8_g40570 * dy);
            double cse_p69_t38_g5133 = (cse_p564_t8_g40570 * dz);
            double cse_p60_t55_g4848 = (cse_p564_t17_g40579 + cse_p69_t38_g5133);
            double cse_p60_t68_g4861 = (cse_p69_t33_g5128 + cse_p69_t36_g5131);
            double cse_p60_t69_g4862 = (cse_p69_t34_g5129 + cse_p69_t37_g5132);
            double h_x1_z2 = cse_p60_t68_g4861;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
            double cse_p67_t70_g4967 = (cse_p564_t13_g40575 + cse_p69_t48_g5143);
            double cse_p60_t78_g4871 = (cse_p564_t13_g40575 + cse_p69_t26_g5121);
            double h_y1_y1 = cse_p67_t70_g4967;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
            double h_y1_z1 = cse_p60_t60_g4853;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
            double h_y1_x2 = cse_p60_t61_g4854;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
            double h_y1_y2 = cse_p60_t54_g4847;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
            double h_y1_z2 = cse_p60_t69_g4862;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
            double cse_p67_t72_g4969 = (cse_p564_t11_g40573 + cse_p564_t14_g40576);
            double cse_p60_t80_g4873 = (cse_p564_t14_g40576 + cse_p69_t35_g5130);
            double h_z1_z1 = cse_p67_t72_g4969;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
            double h_z1_x2 = cse_p60_t65_g4858;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
            double h_z1_y2 = cse_p60_t67_g4860;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
            double h_z1_z2 = cse_p60_t55_g4848;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
            double h_x2_x2 = cse_p60_t82_g4875;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
            double h_x2_y2 = cse_p60_t63_g4856;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
            double h_x2_z2 = cse_p60_t64_g4857;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
            double h_y2_y2 = cse_p60_t78_g4871;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
            double h_y2_z2 = cse_p60_t66_g4859;
            KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
            double h_z2_z2 = cse_p60_t80_g4873;
            KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
          } else {
            /* !COND-ELSE-3 */
            if ((rr <= params.rcut)) {
              /* !COND-4 */
              double cse_p68_t2_g4999 = (-2.00000000000000000e+0 * term.bb * rr);
              double cse_p12_t7_g126 = (term.bb * rr);
              double cse_p11_t1_g112 = (rr * rr);
              double cse_p68_t4_g5001 = (-3.00000000000000000e+0 * term.aa * cse_p11_t1_g112);
              double cse_p67_t56_g4953 = (cse_p68_t2_g4999 + cse_p68_t4_g5001);
              double cse_p12_t6_g125 = (3.00000000000000000e+0 * term.aa * cse_p11_t1_g112);
              double cse_p11_t88_g92 = (term.cc + cse_p12_t6_g125 + cse_p12_t7_g126 + cse_p12_t7_g126);
              double fa_atr = (term.dd + (term.cc * rr) + (cse_p11_t1_g112 * (term.bb + (term.aa * rr))));
              double de_dr = (term.cc + cse_p12_t6_g125 + (2.00000000000000000e+0 * term.bb * rr));
              double d2e_dr2 = ((2.00000000000000000e+0 * term.bb) + (6.00000000000000000e+0 * term.aa * rr));
              double energy = fa_atr;
              energy_added += energy;
              *energy_accumulate += energy;
              double g_x1 = (cse_p11_t88_g92 * cse_p1_t2_invsqrt2 * dx);
              KernelGradientAcc(term.i3x1, 0, g_x1);
              double g_y1 = (cse_p11_t88_g92 * cse_p1_t2_invsqrt2 * dy);
              KernelGradientAcc(term.i3x1, 1, g_y1);
              double g_z1 = (cse_p11_t88_g92 * cse_p1_t2_invsqrt2 * dz);
              KernelGradientAcc(term.i3x1, 2, g_z1);
              double g_x2 = ((cse_p1_t2_invsqrt2 * cse_p67_t56_g4953 * dx) + (-((term.cc * cse_p1_t2_invsqrt2 * dx))));
              KernelGradientAcc(term.i3x2, 0, g_x2);
              double g_y2 = ((cse_p1_t2_invsqrt2 * cse_p67_t56_g4953 * dy) + (-((term.cc * cse_p1_t2_invsqrt2 * dy))));
              KernelGradientAcc(term.i3x2, 1, g_y2);
              double g_z2 = ((cse_p1_t2_invsqrt2 * cse_p67_t56_g4953 * dz) + (-((term.cc * cse_p1_t2_invsqrt2 * dz))));
              KernelGradientAcc(term.i3x2, 2, g_z2);
              double cse_p565_t2_invr40642 = (1.0 / (rr));
              double cse_p11_t1_g110 = (cse_p565_t2_invr40642 * cse_p565_t2_invr40642 * d2e_dr2);
              double cse_p565_t1_invr40641 = (1.0 / ((rr * rr * rr)));
              double cse_p11_t1_g111 = (cse_p565_t1_invr40641 * de_dr);
              double cse_p564_t1_g40640 = (cse_p11_t1_g111 * (cse_p11_t1_g112 + (-(cse_p11_t49_g53))));
              double cse_p67_t1_g4989 = (cse_p11_t1_g110 * dx);
              double cse_p564_t40_g40602 = (cse_p67_t1_g4989 * dz);
              double h_x1_x1 = (cse_p564_t1_g40640 + (cse_p67_t1_g4989 * dx));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
              double cse_p67_t1_g4990 = (cse_p11_t1_g111 * dy);
              double cse_p564_t41_g40603 = (cse_p67_t1_g4990 * (-(dx)));
              double cse_p564_t42_g40604 = (cse_p67_t1_g4990 * dx);
              double h_x1_y1 = (cse_p564_t41_g40603 + (cse_p67_t1_g4989 * dy));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
              double cse_p11_t1_g115 = (cse_p11_t1_g111 * dz);
              double cse_p564_t29_g40591 = (cse_p11_t1_g115 * (-(dx)));
              double cse_p564_t30_g40592 = (cse_p11_t1_g115 * (-(dy)));
              double cse_p564_t31_g40593 = (cse_p11_t1_g115 * dx);
              double cse_p564_t32_g40594 = (cse_p11_t1_g115 * dy);
              double cse_p67_t1_g4992 = (cse_p11_t1_g110 * dy);
              double cse_p564_t43_g40605 = (cse_p67_t1_g4992 * dy);
              double cse_p564_t44_g40606 = (cse_p67_t1_g4992 * dz);
              double cse_p116_t57_g9344 = (cse_p564_t30_g40592 + cse_p564_t44_g40606);
              double h_x1_z1 = (cse_p564_t29_g40591 + cse_p564_t40_g40602);
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
              double cse_p11_t1_g113 = (cse_p11_t1_g110 * (-(dx)));
              double cse_p564_t20_g40582 = (cse_p11_t1_g113 * (-(dx)));
              double cse_p564_t21_g40583 = (cse_p11_t1_g113 * (-(dy)));
              double cse_p564_t22_g40584 = (cse_p11_t1_g113 * (-(dz)));
              double cse_p564_t23_g40585 = (cse_p11_t1_g113 * dy);
              double cse_p564_t24_g40586 = (cse_p11_t1_g113 * dz);
              double cse_p116_t52_g9339 = (cse_p564_t22_g40584 + cse_p564_t29_g40591);
              double cse_p116_t53_g9340 = (cse_p564_t23_g40585 + cse_p564_t42_g40604);
              double cse_p116_t54_g9341 = (cse_p564_t24_g40586 + cse_p564_t31_g40593);
              double cse_p67_t60_g4957 = (cse_p564_t21_g40583 + cse_p564_t41_g40603);
              double h_x1_x2 = ((cse_p11_t1_g111 * (cse_p11_t49_g53 + (-(cse_p11_t1_g112)))) + (cse_p11_t1_g113 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
              double cse_p11_t1_g114 = (cse_p11_t1_g110 * (-(dy)));
              double cse_p564_t25_g40587 = (cse_p11_t1_g114 * (-(dy)));
              double cse_p564_t26_g40588 = (cse_p11_t1_g114 * (-(dz)));
              double cse_p564_t27_g40589 = (cse_p11_t1_g114 * dy);
              double cse_p564_t28_g40590 = (cse_p11_t1_g114 * dz);
              double cse_p116_t50_g9337 = (cse_p564_t27_g40589 + (cse_p11_t1_g111 * (cse_p11_t50_g54 + (-(cse_p11_t1_g112)))));
              double cse_p116_t55_g9342 = (cse_p564_t26_g40588 + cse_p564_t30_g40592);
              double cse_p116_t56_g9343 = (cse_p564_t28_g40590 + cse_p564_t32_g40594);
              double h_x1_y2 = (cse_p564_t42_g40604 + (cse_p11_t1_g114 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
              double cse_p11_t1_g116 = (cse_p11_t1_g110 * (-(dz)));
              double cse_p564_t33_g40595 = (cse_p11_t1_g116 * (-(dz)));
              double cse_p564_t34_g40596 = (cse_p11_t1_g116 * dy);
              double cse_p564_t35_g40597 = (cse_p11_t1_g116 * dz);
              double cse_p116_t51_g9338 = (cse_p564_t35_g40597 + (cse_p11_t1_g111 * (cse_p11_t51_g55 + (-(cse_p11_t1_g112)))));
              double cse_p116_t58_g9345 = (cse_p564_t32_g40594 + cse_p564_t34_g40596);
              double h_x1_z2 = (cse_p564_t31_g40593 + (cse_p11_t1_g116 * dx));
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
              double cse_p11_t1_g117 = (cse_p11_t1_g111 * (cse_p11_t1_g112 + (-(cse_p11_t50_g54))));
              double cse_p116_t63_g9350 = (cse_p11_t1_g117 + cse_p564_t25_g40587);
              double h_y1_y1 = (cse_p11_t1_g117 + cse_p564_t43_g40605);
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
              double h_y1_z1 = cse_p116_t57_g9344;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
              double h_y1_x2 = cse_p116_t53_g9340;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
              double h_y1_y2 = cse_p116_t50_g9337;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
              double h_y1_z2 = cse_p116_t58_g9345;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
              double cse_p11_t1_g118 = (cse_p11_t1_g111 * (cse_p11_t1_g112 + (-(cse_p11_t51_g55))));
              double cse_p116_t64_g9351 = (cse_p11_t1_g118 + cse_p564_t33_g40595);
              double h_z1_z1 = (cse_p11_t1_g118 + (cse_p11_t1_g110 * dz * dz));
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
              double h_z1_x2 = cse_p116_t54_g9341;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
              double h_z1_y2 = cse_p116_t56_g9343;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
              double h_z1_z2 = cse_p116_t51_g9338;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
              double h_x2_x2 = (cse_p564_t1_g40640 + cse_p564_t20_g40582);
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
              double h_x2_y2 = cse_p67_t60_g4957;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
              double h_x2_z2 = cse_p116_t52_g9339;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
              double h_y2_y2 = cse_p116_t63_g9350;
              KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
              double h_y2_z2 = cse_p116_t55_g9342;
              KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
              double h_z2_z2 = cse_p116_t64_g9351;
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
                double h_x1_x1 = ((d2e_dr2 * dx * dx * (1.0 / (rr)) * (1.0 / (rr))) + (de_dr * ((rr * rr) + (-(cse_p11_t49_g53))) * (1.0 / ((rr * rr * rr)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
                double cse_p565_t2_invr40642 = (1.0 / (rr));
                double cse_p565_t1_invr40641 = (1.0 / ((rr * rr * rr)));
                double cse_p11_t1_g112 = (rr * rr);
                double cse_p11_t1_g111 = (cse_p565_t1_invr40641 * de_dr);
                double cse_p11_t1_g110 = (cse_p565_t2_invr40642 * cse_p565_t2_invr40642 * d2e_dr2);
                double h_x1_y1 = ((cse_p11_t1_g110 * dx * dy) + (cse_p11_t1_g111 * dy * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
                double cse_p67_t1_g4992 = (cse_p11_t1_g110 * dy);
                double cse_p67_t1_g4990 = (cse_p11_t1_g111 * dy);
                double cse_p67_t1_g4989 = (cse_p11_t1_g110 * dx);
                double h_x1_z1 = ((cse_p11_t1_g111 * dz * (-(dx))) + (cse_p67_t1_g4989 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
                double cse_p11_t1_g115 = (cse_p11_t1_g111 * dz);
                double h_x1_x2 = ((cse_p11_t1_g111 * (cse_p11_t49_g53 + (-(cse_p11_t1_g112)))) + (cse_p67_t1_g4989 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
                double cse_p11_t1_g113 = (cse_p11_t1_g110 * (-(dx)));
                double h_x1_y2 = ((cse_p67_t1_g4989 * (-(dy))) + (cse_p67_t1_g4990 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
                double cse_p11_t1_g114 = (cse_p11_t1_g110 * (-(dy)));
                double h_x1_z2 = ((cse_p11_t1_g115 * dx) + (cse_p67_t1_g4989 * (-(dz))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
                double cse_p11_t1_g116 = (cse_p11_t1_g110 * (-(dz)));
                double h_y1_y1 = ((cse_p11_t1_g111 * (cse_p11_t1_g112 + (-(cse_p11_t50_g54)))) + (cse_p67_t1_g4992 * dy));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
                double cse_p11_t1_g117 = (cse_p11_t1_g111 * (cse_p11_t1_g112 + (-(cse_p11_t50_g54))));
                double h_y1_z1 = ((cse_p11_t1_g115 * (-(dy))) + (cse_p67_t1_g4992 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
                double h_y1_x2 = ((cse_p11_t1_g113 * dy) + (cse_p67_t1_g4990 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
                double h_y1_y2 = ((cse_p11_t1_g111 * (cse_p11_t50_g54 + (-(cse_p11_t1_g112)))) + (cse_p11_t1_g114 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
                double h_y1_z2 = ((cse_p11_t1_g115 * dy) + (cse_p11_t1_g116 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
                double h_z1_z1 = ((cse_p11_t1_g110 * dz * dz) + (cse_p11_t1_g111 * (cse_p11_t1_g112 + (-(cse_p11_t51_g55)))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
                double cse_p11_t1_g118 = (cse_p11_t1_g111 * (cse_p11_t1_g112 + (-(cse_p11_t51_g55))));
                double h_z1_x2 = ((cse_p11_t1_g113 * dz) + (cse_p11_t1_g115 * dx));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
                double h_z1_y2 = ((cse_p11_t1_g114 * dz) + (cse_p11_t1_g115 * dy));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
                double h_z1_z2 = ((cse_p11_t1_g111 * (cse_p11_t51_g55 + (-(cse_p11_t1_g112)))) + (cse_p11_t1_g116 * dz));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
                double h_x2_x2 = ((cse_p11_t1_g111 * (cse_p11_t1_g112 + (-(cse_p11_t49_g53)))) + (cse_p11_t1_g113 * (-(dx))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
                double cse_p564_t1_g40640 = (cse_p11_t1_g111 * (cse_p11_t1_g112 + (-(cse_p11_t49_g53))));
                double h_x2_y2 = ((cse_p11_t1_g113 * (-(dy))) + (cse_p67_t1_g4990 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
                double cse_p67_t1_g4993 = (cse_p67_t1_g4990 * (-(dx)));
                double h_x2_z2 = ((cse_p11_t1_g113 * (-(dz))) + (cse_p11_t1_g115 * (-(dx))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
                double h_y2_y2 = (cse_p11_t1_g117 + (cse_p11_t1_g114 * (-(dy))));
                KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
                double h_y2_z2 = ((cse_p11_t1_g114 * (-(dz))) + (cse_p11_t1_g115 * (-(dy))));
                KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
                double h_z2_z2 = (cse_p11_t1_g118 + (cse_p11_t1_g116 * (-(dz))));
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
