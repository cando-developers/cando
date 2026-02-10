template <typename HESSIAN>
struct Nonbond_Dd_Cutoff {
  static constexpr size_t PositionSize = 6;
  static std::string description() { return "mathkernel-nonbond_dd_cutoff"; };
double gradient(double a, double b, double qq, double invdd, double r_switch, double r_switch2, double r_cut, double r_cut2, double inv_range, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double energy_added = 0.0000000000000000    ;
  {
    /* !BASE */
    DOUBLE x1 = position[i3x1 + 0];
    DOUBLE y1 = position[i3x1 + 1];
    DOUBLE z1 = position[i3x1 + 2];
    DOUBLE x2 = position[i3x2 + 0];
    DOUBLE y2 = position[i3x2 + 1];
    DOUBLE z2 = position[i3x2 + 2];
    {
      /* !BASE */
      double dx = (x1 + (-(x2)));
      double dy = (y1 + (-(y2)));
      double dz = (z1 + (-(z2)));
      double r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      if ((r2 < r_cut2)) {
        /* !CUT */
        double cse_p1_t1_invsqrt28383 = (1.0 / sqrt(r2));
        double cse_p12_t21_g28457 = (cse_p1_t1_invsqrt28383 * cse_p1_t1_invsqrt28383);
        double cse_p13_t16_g28500 = (cse_p12_t21_g28457 * cse_p1_t1_invsqrt28383);
        double r = (cse_p1_t1_invsqrt28383 * r2);
        double cse_p14_t1_g28535 = (cse_p12_t21_g28457 * cse_p12_t21_g28457);
        double invr6 = (cse_p12_t21_g28457 * cse_p14_t1_g28535);
        double cse_p12_t1_g28476 = (a * invr6);
        double e_lj = ((cse_p12_t1_g28476 * invr6) + (-((b * invr6))));
        double cse_p12_t1_g28482 = (invdd * qq);
        double e_coul = (cse_p12_t1_g28482 * cse_p12_t21_g28457);
        double e_base = (e_coul + e_lj);
        double de_base_dr = ((-2.0000000000000000     * cse_p12_t1_g28482 * cse_p13_t16_g28500) + (cse_p1_t1_invsqrt28383 * invr6 * ((-12.000000000000000     * cse_p12_t1_g28476) + (6.0000000000000000     * b))));
        if ((r2 < r_switch2)) {
          /* !SWITCH */
          double cse_p519_t1_g34992 = (-3.0000000000000000     * b);
          double cse_p519_t2_g34981 = (3.0000000000000000     * cse_p12_t1_g28476);
          double cse_p14_t16_g28534 = (cse_p519_t1_g34992 + cse_p519_t2_g34981 + cse_p519_t2_g34981);
          double cse_p13_t1_g28518 = (cse_p14_t16_g28534 * cse_p14_t1_g28535);
          double cse_p519_t3_g34982 = (cse_p12_t1_g28482 + cse_p13_t1_g28518);
          double energy = e_base;
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p12_t1_g28478 = (cse_p13_t16_g28500 * cse_p1_t1_invsqrt28383 * cse_p519_t3_g34982);
          double cse_p14_t1_g28537 = (-2.0000000000000000     * cse_p12_t1_g28478);
          double g_x1 = (cse_p14_t1_g28537 * dx);
          KernelGradientAcc(i3x1, 0, g_x1);
          double g_y1 = (cse_p14_t1_g28537 * dy);
          KernelGradientAcc(i3x1, 1, g_y1);
          double g_z1 = (cse_p14_t1_g28537 * dz);
          KernelGradientAcc(i3x1, 2, g_z1);
          double cse_p14_t1_g28540 = (2.0000000000000000     * cse_p12_t1_g28478);
          double g_x2 = (cse_p14_t1_g28540 * dx);
          KernelGradientAcc(i3x2, 0, g_x2);
          double g_y2 = (cse_p14_t1_g28540 * dy);
          KernelGradientAcc(i3x2, 1, g_y2);
          double g_z2 = (cse_p14_t1_g28540 * dz);
          KernelGradientAcc(i3x2, 2, g_z2);
          double de_dr = de_base_dr;
        } else {
          /* !SKIN */
          double cse_p11_t1_g28429 = (cse_p13_t16_g28500 * r2);
          double cse_p11_t1_g28431 = (-(cse_p11_t1_g28429));
          double cse_p14_t3_g28521 = (cse_p11_t1_g28431 * dx);
          double cse_p14_t4_g28522 = (cse_p11_t1_g28431 * dy);
          double cse_p14_t5_g28523 = (cse_p11_t1_g28431 * dz);
          double cse_p11_t1_g28432 = (-2.0000000000000000     * cse_p1_t1_invsqrt28383);
          double cse_p14_t6_g28524 = (cse_p11_t1_g28432 * dx);
          double cse_p14_t7_g28525 = (cse_p11_t1_g28432 * dy);
          double cse_p14_t8_g28526 = (cse_p11_t1_g28432 * dz);
          double cse_p14_t9_g28527 = (cse_p11_t1_g28429 * dx);
          double cse_p14_t10_g28528 = (cse_p11_t1_g28429 * dy);
          double cse_p14_t11_g28529 = (cse_p11_t1_g28429 * dz);
          double cse_p13_t20_g28504 = (cse_p14_t6_g28524 + cse_p14_t9_g28527);
          double cse_p13_t21_g28505 = (cse_p14_t10_g28528 + cse_p14_t7_g28525);
          double cse_p13_t22_g28506 = (cse_p14_t11_g28529 + cse_p14_t8_g28526);
          double cse_p11_t1_g28435 = (2.0000000000000000     * cse_p1_t1_invsqrt28383);
          double cse_p14_t12_g28530 = (cse_p11_t1_g28435 * dx);
          double cse_p14_t13_g28531 = (cse_p11_t1_g28435 * dy);
          double cse_p14_t14_g28532 = (cse_p11_t1_g28435 * dz);
          double cse_p13_t17_g28501 = (cse_p14_t12_g28530 + cse_p14_t3_g28521);
          double cse_p13_t18_g28502 = (cse_p14_t13_g28531 + cse_p14_t4_g28522);
          double cse_p13_t19_g28503 = (cse_p14_t14_g28532 + cse_p14_t5_g28523);
          double drs = (r + (-(r_switch)));
          double t1 = (drs * inv_range);
          double cse_p11_t1_g28430 = (inv_range * t1);
          double cse_p11_t1_g28434 = (2.0000000000000000     * cse_p11_t1_g28430);
          double cse_p11_t1_g28436 = (cse_p11_t1_g28434 * t1);
          double cse_p12_t15_g28451 = (cse_p11_t1_g28436 * cse_p13_t17_g28501);
          double cse_p12_t16_g28452 = (cse_p11_t1_g28436 * cse_p13_t18_g28502);
          double cse_p12_t17_g28453 = (cse_p11_t1_g28436 * cse_p13_t19_g28503);
          double cse_p12_t18_g28454 = (cse_p11_t1_g28436 * cse_p13_t20_g28504);
          double cse_p12_t19_g28455 = (cse_p11_t1_g28436 * cse_p13_t21_g28505);
          double cse_p12_t20_g28456 = (cse_p11_t1_g28436 * cse_p13_t22_g28506);
          double t2 = (t1 * t1);
          double cse_p12_t1_g28481 = (inv_range * t2);
          double cse_p12_t25_g28461 = (cse_p12_t1_g28481 * cse_p13_t17_g28501);
          double cse_p12_t26_g28462 = (cse_p12_t1_g28481 * cse_p13_t18_g28502);
          double cse_p12_t27_g28463 = (cse_p12_t1_g28481 * cse_p13_t19_g28503);
          double cse_p12_t28_g28464 = (cse_p12_t1_g28481 * cse_p13_t20_g28504);
          double cse_p12_t29_g28465 = (cse_p12_t1_g28481 * cse_p13_t21_g28505);
          double cse_p12_t30_g28466 = (cse_p12_t1_g28481 * cse_p13_t22_g28506);
          double cse_p11_t32_g28416 = (cse_p12_t15_g28451 + cse_p12_t25_g28461);
          double cse_p11_t33_g28417 = (cse_p12_t16_g28452 + cse_p12_t26_g28462);
          double cse_p11_t34_g28418 = (cse_p12_t17_g28453 + cse_p12_t27_g28463);
          double cse_p11_t35_g28419 = (cse_p12_t18_g28454 + cse_p12_t28_g28464);
          double cse_p11_t36_g28420 = (cse_p12_t19_g28455 + cse_p12_t29_g28465);
          double cse_p11_t37_g28421 = (cse_p12_t20_g28456 + cse_p12_t30_g28466);
          double t3 = (t1 * t2);
          double t4 = (t2 * t2);
          double t5 = (t2 * t3);
          double s = (1.0000000000000000     + (-10.000000000000000     * t3) + (-6.0000000000000000     * t5) + (15.000000000000000     * t4));
          double ds_dt = ((-30.000000000000000     * t2) + (-30.000000000000000     * t4) + (60.000000000000000     * t3));
          double ds_dr = (ds_dt * inv_range);
          double energy = (e_base * s);
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p12_t1_g28479 = (3.0000000000000000     * cse_p12_t1_g28476);
          double cse_p519_t1_g34992 = (-3.0000000000000000     * b);
          double cse_p13_t1_g28518 = (cse_p14_t1_g28535 * (cse_p12_t1_g28479 + cse_p12_t1_g28479 + cse_p519_t1_g34992));
          double cse_p11_t1_g28428 = (cse_p13_t16_g28500 * cse_p1_t1_invsqrt28383 * s * (cse_p12_t1_g28482 + cse_p13_t1_g28518));
          double cse_p12_t1_g28477 = (60.000000000000000     * cse_p11_t1_g28430 * t2);
          double cse_p12_t1_g28480 = (cse_p11_t1_g28434 * t3);
          double cse_p14_t1_g28536 = (-2.0000000000000000     * cse_p11_t1_g28428);
          double g_x1 = ((cse_p14_t1_g28536 * dx) + (e_base * ((-10.000000000000000     * cse_p11_t32_g28416) + (-6.0000000000000000     * ((cse_p11_t32_g28416 * t2) + (cse_p12_t1_g28480 * cse_p13_t17_g28501))) + (cse_p12_t1_g28477 * cse_p13_t17_g28501))));
          KernelGradientAcc(i3x1, 0, g_x1);
          double g_y1 = ((cse_p14_t1_g28536 * dy) + (e_base * ((-10.000000000000000     * cse_p11_t33_g28417) + (-6.0000000000000000     * ((cse_p11_t33_g28417 * t2) + (cse_p12_t1_g28480 * cse_p13_t18_g28502))) + (cse_p12_t1_g28477 * cse_p13_t18_g28502))));
          KernelGradientAcc(i3x1, 1, g_y1);
          double g_z1 = ((cse_p14_t1_g28536 * dz) + (e_base * ((-10.000000000000000     * cse_p11_t34_g28418) + (-6.0000000000000000     * ((cse_p11_t34_g28418 * t2) + (cse_p12_t1_g28480 * cse_p13_t19_g28503))) + (cse_p12_t1_g28477 * cse_p13_t19_g28503))));
          KernelGradientAcc(i3x1, 2, g_z1);
          double cse_p14_t1_g28539 = (2.0000000000000000     * cse_p11_t1_g28428);
          double g_x2 = ((cse_p14_t1_g28539 * dx) + (e_base * ((-10.000000000000000     * cse_p11_t35_g28419) + (-6.0000000000000000     * ((cse_p11_t35_g28419 * t2) + (cse_p12_t1_g28480 * cse_p13_t20_g28504))) + (cse_p12_t1_g28477 * cse_p13_t20_g28504))));
          KernelGradientAcc(i3x2, 0, g_x2);
          double g_y2 = ((cse_p14_t1_g28539 * dy) + (e_base * ((-10.000000000000000     * cse_p11_t36_g28420) + (-6.0000000000000000     * ((cse_p11_t36_g28420 * t2) + (cse_p12_t1_g28480 * cse_p13_t21_g28505))) + (cse_p12_t1_g28477 * cse_p13_t21_g28505))));
          KernelGradientAcc(i3x2, 1, g_y2);
          double g_z2 = ((cse_p14_t1_g28539 * dz) + (e_base * ((-10.000000000000000     * cse_p11_t37_g28421) + (-6.0000000000000000     * ((cse_p11_t37_g28421 * t2) + (cse_p12_t1_g28480 * cse_p13_t22_g28506))) + (cse_p12_t1_g28477 * cse_p13_t22_g28506))));
          KernelGradientAcc(i3x2, 2, g_z2);
          double de_dr = ((de_base_dr * s) + (ds_dr * e_base));
        }
      }
    }
  }
  return energy_added;
}
void gradient_fd(double a, double b, double qq, double invdd, double r_switch, double r_switch2, double r_cut, double r_cut2, double inv_range, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 6;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  double e0 = 0.0;
  energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
}

double energy(double a, double b, double qq, double invdd, double r_switch, double r_switch2, double r_cut, double r_cut2, double inv_range, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double energy_added = 0.0000000000000000    ;
  {
    /* !BASE */
    DOUBLE x1 = position[i3x1 + 0];
    DOUBLE y1 = position[i3x1 + 1];
    DOUBLE z1 = position[i3x1 + 2];
    DOUBLE x2 = position[i3x2 + 0];
    DOUBLE y2 = position[i3x2 + 1];
    DOUBLE z2 = position[i3x2 + 2];
    {
      /* !BASE */
      double dx = (x1 + (-(x2)));
      double dy = (y1 + (-(y2)));
      double dz = (z1 + (-(z2)));
      double r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      if ((r2 < r_cut2)) {
        /* !CUT */
        double invr = (1.0 / sqrt(r2));
        double r = (invr * r2);
        double invr2 = (invr * invr);
        double invr3 = (invr2 * invr);
        double cse_p11_t1_g34993 = (invr2 * invr2);
        double invr6 = (cse_p11_t1_g34993 * invr2);
        double e_lj = ((a * invr6 * invr6) + (-((b * invr6))));
        double e_coul = (invdd * invr2 * qq);
        double e_base = (e_coul + e_lj);
        if ((r2 < r_switch2)) {
          /* !SWITCH */
          double energy = e_base;
          energy_added += energy;
          *energy_accumulate += energy;
        } else {
          /* !SKIN */
          double drs = (r + (-(r_switch)));
          double t1 = (drs * inv_range);
          double t2 = (t1 * t1);
          double t3 = (t1 * t2);
          double t4 = (t2 * t2);
          double t5 = (t2 * t3);
          double s = (1.0000000000000000     + (-10.000000000000000     * t3) + (-6.0000000000000000     * t5) + (15.000000000000000     * t4));
          double energy = (e_base * s);
          energy_added += energy;
          *energy_accumulate += energy;
        }
      }
    }
  }
  return energy_added;
}
void energy_fd(double a, double b, double qq, double invdd, double r_switch, double r_switch2, double r_cut, double r_cut2, double inv_range, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, energy_accumulate, force, hessian, dvec, hdvec);
}

double hessian(double a, double b, double qq, double invdd, double r_switch, double r_switch2, double r_cut, double r_cut2, double inv_range, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double energy_added = 0.0000000000000000    ;
  {
    /* !BASE */
    DOUBLE x1 = position[i3x1 + 0];
    DOUBLE y1 = position[i3x1 + 1];
    DOUBLE z1 = position[i3x1 + 2];
    DOUBLE x2 = position[i3x2 + 0];
    DOUBLE y2 = position[i3x2 + 1];
    DOUBLE z2 = position[i3x2 + 2];
    {
      /* !BASE */
      double dx = (x1 + (-(x2)));
      double cse_p11_t54_g35049 = (dx * dx);
      double dy = (y1 + (-(y2)));
      double cse_p11_t55_g35050 = (dy * dy);
      double dz = (z1 + (-(z2)));
      double cse_p11_t56_g35051 = (dz * dz);
      double r2 = (cse_p11_t54_g35049 + cse_p11_t55_g35050 + cse_p11_t56_g35051);
      if ((r2 < r_cut2)) {
        /* !CUT */
        double cse_p1_t1_invsqrt34994 = (1.0 / sqrt(r2));
        double cse_p12_t57_g35178 = (cse_p1_t1_invsqrt34994 * cse_p1_t1_invsqrt34994);
        double cse_p13_t60_g35295 = (cse_p12_t57_g35178 * cse_p1_t1_invsqrt34994);
        double r = (cse_p1_t1_invsqrt34994 * r2);
        double cse_p15_t1_g35420 = (cse_p12_t57_g35178 * cse_p12_t57_g35178);
        double invr6 = (cse_p12_t57_g35178 * cse_p15_t1_g35420);
        double cse_p13_t1_g35328 = (a * invr6);
        double e_lj = ((cse_p13_t1_g35328 * invr6) + (-((b * invr6))));
        double cse_p65_t1_g36441 = (invdd * qq);
        double e_coul = (cse_p12_t57_g35178 * cse_p65_t1_g36441);
        double e_base = (e_coul + e_lj);
        double de_base_dr = ((-2.0000000000000000     * cse_p13_t60_g35295 * cse_p65_t1_g36441) + (cse_p1_t1_invsqrt34994 * invr6 * ((-12.000000000000000     * cse_p13_t1_g35328) + (6.0000000000000000     * b))));
        double d2e_base_dr2 = ((6.0000000000000000     * cse_p15_t1_g35420 * cse_p65_t1_g36441) + (cse_p12_t57_g35178 * invr6 * ((-42.000000000000000     * b) + (156.00000000000000     * cse_p13_t1_g35328))));
        if ((r2 < r_switch2)) {
          /* !SWITCH */
          double cse_p522_t1_g44083 = (-3.0000000000000000     * b);
          double cse_p522_t2_g44084 = (3.0000000000000000     * cse_p13_t1_g35328);
          double cse_p522_t8_g44090 = (r * r);
          double cse_p522_t9_g44091 = (-(cse_p11_t54_g35049));
          double cse_p522_t10_g44092 = (-(cse_p11_t55_g35050));
          double cse_p522_t11_g44093 = (-(cse_p11_t56_g35051));
          double cse_p522_t12_g44094 = (-(dx));
          double cse_p522_t13_g44095 = (-(dy));
          double cse_p522_t14_g44096 = (-(dz));
          double cse_p523_t1_invr44100 = (1.0 / ((r * r * r)));
          double cse_p523_t2_invr44101 = (1.0 / (r));
          double cse_p522_t17_g44099 = (r * r * r);
          double cse_p16_t16_g35436 = (cse_p522_t1_g44083 + cse_p522_t2_g44084 + cse_p522_t2_g44084);
          double cse_p15_t12_g35396 = (cse_p15_t1_g35420 * cse_p16_t16_g35436);
          double cse_p15_t25_g35409 = (-(cse_p522_t8_g44090));
          double cse_p14_t38_g35367 = (cse_p15_t12_g35396 + cse_p65_t1_g36441);
          double cse_p14_t32_g35361 = (cse_p11_t54_g35049 + cse_p15_t25_g35409);
          double cse_p14_t33_g35362 = (cse_p11_t55_g35050 + cse_p15_t25_g35409);
          double cse_p14_t34_g35363 = (cse_p11_t56_g35051 + cse_p15_t25_g35409);
          double cse_p14_t39_g35368 = (cse_p522_t8_g44090 + cse_p522_t9_g44091);
          double cse_p14_t40_g35369 = (cse_p522_t10_g44092 + cse_p522_t8_g44090);
          double cse_p14_t41_g35370 = (cse_p522_t11_g44093 + cse_p522_t8_g44090);
          double energy = e_base;
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p12_t1_g35235 = (cse_p13_t60_g35295 * cse_p14_t38_g35367 * cse_p1_t1_invsqrt34994);
          double cse_p65_t1_g36445 = (-2.0000000000000000     * cse_p12_t1_g35235);
          double g_x1 = (cse_p65_t1_g36445 * dx);
          KernelGradientAcc(i3x1, 0, g_x1);
          double g_y1 = (cse_p65_t1_g36445 * dy);
          KernelGradientAcc(i3x1, 1, g_y1);
          double g_z1 = (cse_p65_t1_g36445 * dz);
          KernelGradientAcc(i3x1, 2, g_z1);
          double cse_p65_t1_g36448 = (2.0000000000000000     * cse_p12_t1_g35235);
          double g_x2 = (cse_p65_t1_g36448 * dx);
          KernelGradientAcc(i3x2, 0, g_x2);
          double g_y2 = (cse_p65_t1_g36448 * dy);
          KernelGradientAcc(i3x2, 1, g_y2);
          double g_z2 = (cse_p65_t1_g36448 * dz);
          KernelGradientAcc(i3x2, 2, g_z2);
          double de_dr = de_base_dr;
          double cse_p11_t1_g35113 = (cse_p523_t1_invr44100 * de_base_dr);
          double cse_p522_t6_g44088 = (cse_p11_t1_g35113 * dy);
          double cse_p522_t7_g44089 = (cse_p11_t1_g35113 * dz);
          double cse_p15_t1_g35419 = (cse_p522_t12_g44094 * cse_p522_t6_g44088);
          double cse_p15_t17_g35401 = (cse_p522_t6_g44088 * dx);
          double cse_p13_t24_g35259 = (cse_p11_t1_g35113 * cse_p14_t32_g35361);
          double cse_p13_t25_g35260 = (cse_p11_t1_g35113 * cse_p14_t33_g35362);
          double cse_p13_t26_g35261 = (cse_p11_t1_g35113 * cse_p14_t34_g35363);
          double cse_p13_t30_g35265 = (cse_p11_t1_g35113 * cse_p14_t39_g35368);
          double cse_p13_t31_g35266 = (cse_p11_t1_g35113 * cse_p14_t40_g35369);
          double cse_p13_t32_g35267 = (cse_p11_t1_g35113 * cse_p14_t41_g35370);
          double cse_p13_t55_g35290 = (cse_p522_t12_g44094 * cse_p522_t7_g44089);
          double cse_p13_t56_g35291 = (cse_p522_t13_g44095 * cse_p522_t7_g44089);
          double cse_p13_t57_g35292 = (cse_p522_t7_g44089 * dx);
          double cse_p13_t58_g35293 = (cse_p522_t7_g44089 * dy);
          double d2e_dr2 = d2e_base_dr2;
          double cse_p11_t1_g35112 = (cse_p523_t2_invr44101 * cse_p523_t2_invr44101 * d2e_base_dr2);
          double cse_p522_t3_g44085 = (cse_p11_t1_g35112 * dx);
          double cse_p522_t4_g44086 = (cse_p11_t1_g35112 * dy);
          double cse_p522_t5_g44087 = (cse_p11_t1_g35112 * dz * dz);
          double cse_p66_t8_g36456 = (cse_p522_t4_g44086 * dy);
          double cse_p66_t9_g36457 = (cse_p522_t4_g44086 * dz);
          double cse_p60_t3_g36394 = (cse_p11_t1_g35112 * cse_p522_t12_g44094);
          double cse_p60_t4_g36395 = (cse_p11_t1_g35112 * cse_p522_t13_g44095);
          double cse_p60_t5_g36396 = (cse_p11_t1_g35112 * cse_p522_t14_g44096);
          double cse_p16_t13_g35433 = (cse_p522_t3_g44085 * dy);
          double cse_p15_t13_g35397 = (cse_p522_t3_g44085 * dx);
          double cse_p15_t15_g35399 = (cse_p522_t3_g44085 * dz);
          double cse_p15_t20_g35404 = (cse_p15_t1_g35419 + cse_p16_t13_g35433);
          double cse_p13_t44_g35279 = (cse_p522_t12_g44094 * cse_p60_t3_g36394);
          double cse_p13_t45_g35280 = (cse_p522_t13_g44095 * cse_p60_t3_g36394);
          double cse_p13_t46_g35281 = (cse_p522_t14_g44096 * cse_p60_t3_g36394);
          double cse_p13_t47_g35282 = (cse_p60_t3_g36394 * dx);
          double cse_p13_t48_g35283 = (cse_p60_t3_g36394 * dy);
          double cse_p13_t49_g35284 = (cse_p60_t3_g36394 * dz);
          double cse_p13_t50_g35285 = (cse_p522_t13_g44095 * cse_p60_t4_g36395);
          double cse_p13_t51_g35286 = (cse_p522_t14_g44096 * cse_p60_t4_g36395);
          double cse_p13_t52_g35287 = (cse_p60_t4_g36395 * dx);
          double cse_p13_t53_g35288 = (cse_p60_t4_g36395 * dy);
          double cse_p13_t54_g35289 = (cse_p60_t4_g36395 * dz);
          double cse_p13_t40_g35275 = (cse_p522_t14_g44096 * cse_p60_t5_g36396);
          double cse_p13_t41_g35276 = (cse_p60_t5_g36396 * dx);
          double cse_p13_t42_g35277 = (cse_p60_t5_g36396 * dy);
          double cse_p13_t43_g35278 = (cse_p60_t5_g36396 * dz);
          double cse_p12_t79_g35200 = (cse_p13_t32_g35267 + cse_p13_t40_g35275);
          double cse_p12_t80_g35201 = (cse_p13_t41_g35276 + cse_p13_t57_g35292);
          double cse_p12_t81_g35202 = (cse_p13_t30_g35265 + cse_p15_t13_g35397);
          double cse_p12_t83_g35204 = (cse_p13_t55_g35290 + cse_p15_t15_g35399);
          double cse_p12_t84_g35205 = (cse_p13_t42_g35277 + cse_p13_t58_g35293);
          double cse_p12_t85_g35206 = (cse_p13_t31_g35266 + cse_p66_t8_g36456);
          double cse_p12_t86_g35207 = (cse_p13_t56_g35291 + cse_p66_t9_g36457);
          double cse_p12_t87_g35208 = (cse_p13_t26_g35261 + cse_p13_t43_g35278);
          double cse_p12_t88_g35209 = (cse_p13_t32_g35267 + cse_p522_t5_g44087);
          double cse_p12_t68_g35189 = (cse_p13_t24_g35259 + cse_p13_t47_g35282);
          double cse_p12_t70_g35191 = (cse_p13_t30_g35265 + cse_p13_t44_g35279);
          double cse_p12_t89_g35210 = (cse_p13_t48_g35283 + cse_p15_t17_g35401);
          double cse_p12_t91_g35212 = (cse_p13_t45_g35280 + cse_p15_t1_g35419);
          double cse_p12_t92_g35213 = (cse_p13_t46_g35281 + cse_p13_t55_g35290);
          double cse_p12_t93_g35214 = (cse_p13_t49_g35284 + cse_p13_t57_g35292);
          double cse_p12_t69_g35190 = (cse_p13_t25_g35260 + cse_p13_t53_g35288);
          double cse_p12_t71_g35192 = (cse_p13_t31_g35266 + cse_p13_t50_g35285);
          double cse_p12_t90_g35211 = (cse_p13_t52_g35287 + cse_p15_t17_g35401);
          double cse_p12_t94_g35215 = (cse_p13_t51_g35286 + cse_p13_t56_g35291);
          double cse_p12_t95_g35216 = (cse_p13_t54_g35289 + cse_p13_t58_g35293);
          double h_x1_x1 = cse_p12_t81_g35202;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 0, h_x1_x1);
          double h_x1_y1 = cse_p15_t20_g35404;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 1, h_x1_y1);
          double h_x1_z1 = cse_p12_t83_g35204;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 2, h_x1_z1);
          double h_x1_x2 = cse_p12_t68_g35189;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 0, h_x1_x2);
          double h_x1_y2 = cse_p12_t90_g35211;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 1, h_x1_y2);
          double h_x1_z2 = cse_p12_t80_g35201;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 2, h_x1_z2);
          double h_y1_y1 = cse_p12_t85_g35206;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 1, h_y1_y1);
          double h_y1_z1 = cse_p12_t86_g35207;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 2, h_y1_z1);
          double h_y1_x2 = cse_p12_t89_g35210;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 0, h_y1_x2);
          double h_y1_y2 = cse_p12_t69_g35190;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 1, h_y1_y2);
          double h_y1_z2 = cse_p12_t84_g35205;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 2, h_y1_z2);
          double h_z1_z1 = cse_p12_t88_g35209;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x1, 2, h_z1_z1);
          double h_z1_x2 = cse_p12_t93_g35214;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 0, h_z1_x2);
          double h_z1_y2 = cse_p12_t95_g35216;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 1, h_z1_y2);
          double h_z1_z2 = cse_p12_t87_g35208;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 2, h_z1_z2);
          double h_x2_x2 = cse_p12_t70_g35191;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 0, h_x2_x2);
          double h_x2_y2 = cse_p12_t91_g35212;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 1, h_x2_y2);
          double h_x2_z2 = cse_p12_t92_g35213;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 2, h_x2_z2);
          double h_y2_y2 = cse_p12_t71_g35192;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 1, h_y2_y2);
          double h_y2_z2 = cse_p12_t94_g35215;
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 2, h_y2_z2);
          double h_z2_z2 = cse_p12_t79_g35200;
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x2, 2, h_z2_z2);
        } else {
          /* !SKIN */
          double cse_p65_t1_g36443 = (-2.0000000000000000     * cse_p1_t1_invsqrt34994);
          double cse_p14_t4_g35333 = (cse_p65_t1_g36443 * dx);
          double cse_p14_t5_g35334 = (cse_p65_t1_g36443 * dy);
          double cse_p14_t6_g35335 = (cse_p65_t1_g36443 * dz);
          double cse_p65_t1_g36446 = (2.0000000000000000     * cse_p1_t1_invsqrt34994);
          double cse_p14_t11_g35340 = (cse_p65_t1_g36446 * dx);
          double cse_p14_t12_g35341 = (cse_p65_t1_g36446 * dy);
          double cse_p14_t13_g35342 = (cse_p65_t1_g36446 * dz);
          double cse_p11_t1_g35116 = (cse_p13_t60_g35295 * r2);
          double cse_p65_t9_g36427 = (cse_p11_t1_g35116 * dx);
          double cse_p65_t10_g36428 = (cse_p11_t1_g35116 * dy);
          double cse_p65_t11_g36429 = (cse_p11_t1_g35116 * dz);
          double cse_p14_t1_g35330 = (-(cse_p65_t9_g36427));
          double cse_p14_t2_g35331 = (-(cse_p65_t10_g36428));
          double cse_p14_t3_g35332 = (-(cse_p65_t11_g36429));
          double cse_p13_t64_g35299 = (cse_p14_t11_g35340 + cse_p14_t1_g35330);
          double cse_p13_t65_g35300 = (cse_p14_t12_g35341 + cse_p14_t2_g35331);
          double cse_p13_t66_g35301 = (cse_p14_t13_g35342 + cse_p14_t3_g35332);
          double cse_p13_t67_g35302 = (cse_p14_t4_g35333 + cse_p65_t9_g36427);
          double cse_p13_t68_g35303 = (cse_p14_t5_g35334 + cse_p65_t10_g36428);
          double cse_p13_t69_g35304 = (cse_p14_t6_g35335 + cse_p65_t11_g36429);
          double drs = (r + (-(r_switch)));
          double t1 = (drs * inv_range);
          double cse_p11_t1_g35119 = (inv_range * t1);
          double cse_p13_t1_g35327 = (2.0000000000000000     * cse_p11_t1_g35119);
          double cse_p15_t1_g35416 = (cse_p13_t1_g35327 * t1);
          double cse_p12_t11_g35132 = (cse_p13_t64_g35299 * cse_p15_t1_g35416);
          double cse_p12_t12_g35133 = (cse_p13_t65_g35300 * cse_p15_t1_g35416);
          double cse_p12_t13_g35134 = (cse_p13_t66_g35301 * cse_p15_t1_g35416);
          double cse_p12_t14_g35135 = (cse_p13_t67_g35302 * cse_p15_t1_g35416);
          double cse_p12_t15_g35136 = (cse_p13_t68_g35303 * cse_p15_t1_g35416);
          double cse_p12_t16_g35137 = (cse_p13_t69_g35304 * cse_p15_t1_g35416);
          double t2 = (t1 * t1);
          double cse_p15_t1_g35418 = (inv_range * t2);
          double cse_p12_t61_g35182 = (cse_p13_t64_g35299 * cse_p15_t1_g35418);
          double cse_p12_t62_g35183 = (cse_p13_t65_g35300 * cse_p15_t1_g35418);
          double cse_p12_t63_g35184 = (cse_p13_t66_g35301 * cse_p15_t1_g35418);
          double cse_p12_t64_g35185 = (cse_p13_t67_g35302 * cse_p15_t1_g35418);
          double cse_p12_t65_g35186 = (cse_p13_t68_g35303 * cse_p15_t1_g35418);
          double cse_p12_t66_g35187 = (cse_p13_t69_g35304 * cse_p15_t1_g35418);
          double cse_p11_t69_g35064 = (cse_p12_t11_g35132 + cse_p12_t61_g35182);
          double cse_p11_t70_g35065 = (cse_p12_t12_g35133 + cse_p12_t62_g35183);
          double cse_p11_t71_g35066 = (cse_p12_t13_g35134 + cse_p12_t63_g35184);
          double cse_p11_t72_g35067 = (cse_p12_t14_g35135 + cse_p12_t64_g35185);
          double cse_p11_t73_g35068 = (cse_p12_t15_g35136 + cse_p12_t65_g35186);
          double cse_p11_t74_g35069 = (cse_p12_t16_g35137 + cse_p12_t66_g35187);
          double t3 = (t1 * t2);
          double t4 = (t2 * t2);
          double t5 = (t2 * t3);
          double s = (1.0000000000000000     + (-10.000000000000000     * t3) + (-6.0000000000000000     * t5) + (15.000000000000000     * t4));
          double ds_dt = ((-30.000000000000000     * t2) + (-30.000000000000000     * t4) + (60.000000000000000     * t3));
          double d2s_dt2 = ((-120.00000000000000     * t3) + (-60.000000000000000     * t1) + (180.00000000000000     * t2));
          double ds_dr = (ds_dt * inv_range);
          double d2s_dr2 = (d2s_dt2 * inv_range * inv_range);
          double energy = (e_base * s);
          energy_added += energy;
          *energy_accumulate += energy;
          double cse_p14_t1_g35382 = (3.0000000000000000     * cse_p13_t1_g35328);
          double cse_p13_t1_g35329 = (cse_p15_t1_g35420 * (cse_p14_t1_g35382 + cse_p14_t1_g35382 + (-3.0000000000000000     * b)));
          double cse_p11_t1_g35115 = (cse_p13_t60_g35295 * cse_p1_t1_invsqrt34994 * s * (cse_p13_t1_g35329 + cse_p65_t1_g36441));
          double cse_p14_t1_g35381 = (60.000000000000000     * cse_p11_t1_g35119 * t2);
          double cse_p15_t1_g35417 = (cse_p13_t1_g35327 * t3);
          double cse_p65_t1_g36444 = (-2.0000000000000000     * cse_p11_t1_g35115);
          double g_x1 = ((cse_p65_t1_g36444 * dx) + (e_base * ((-10.000000000000000     * cse_p11_t69_g35064) + (-6.0000000000000000     * ((cse_p11_t69_g35064 * t2) + (cse_p13_t64_g35299 * cse_p15_t1_g35417))) + (cse_p13_t64_g35299 * cse_p14_t1_g35381))));
          KernelGradientAcc(i3x1, 0, g_x1);
          double g_y1 = ((cse_p65_t1_g36444 * dy) + (e_base * ((-10.000000000000000     * cse_p11_t70_g35065) + (-6.0000000000000000     * ((cse_p11_t70_g35065 * t2) + (cse_p13_t65_g35300 * cse_p15_t1_g35417))) + (cse_p13_t65_g35300 * cse_p14_t1_g35381))));
          KernelGradientAcc(i3x1, 1, g_y1);
          double g_z1 = ((cse_p65_t1_g36444 * dz) + (e_base * ((-10.000000000000000     * cse_p11_t71_g35066) + (-6.0000000000000000     * ((cse_p11_t71_g35066 * t2) + (cse_p13_t66_g35301 * cse_p15_t1_g35417))) + (cse_p13_t66_g35301 * cse_p14_t1_g35381))));
          KernelGradientAcc(i3x1, 2, g_z1);
          double cse_p65_t1_g36447 = (2.0000000000000000     * cse_p11_t1_g35115);
          double g_x2 = ((cse_p65_t1_g36447 * dx) + (e_base * ((-10.000000000000000     * cse_p11_t72_g35067) + (-6.0000000000000000     * ((cse_p11_t72_g35067 * t2) + (cse_p13_t67_g35302 * cse_p15_t1_g35417))) + (cse_p13_t67_g35302 * cse_p14_t1_g35381))));
          KernelGradientAcc(i3x2, 0, g_x2);
          double g_y2 = ((cse_p65_t1_g36447 * dy) + (e_base * ((-10.000000000000000     * cse_p11_t73_g35068) + (-6.0000000000000000     * ((cse_p11_t73_g35068 * t2) + (cse_p13_t68_g35303 * cse_p15_t1_g35417))) + (cse_p13_t68_g35303 * cse_p14_t1_g35381))));
          KernelGradientAcc(i3x2, 1, g_y2);
          double g_z2 = ((cse_p65_t1_g36447 * dz) + (e_base * ((-10.000000000000000     * cse_p11_t74_g35069) + (-6.0000000000000000     * ((cse_p11_t74_g35069 * t2) + (cse_p13_t69_g35304 * cse_p15_t1_g35417))) + (cse_p13_t69_g35304 * cse_p14_t1_g35381))));
          KernelGradientAcc(i3x2, 2, g_z2);
          double de_dr = ((de_base_dr * s) + (ds_dr * e_base));
          double d2e_dr2 = ((2.0000000000000000     * de_base_dr * ds_dr) + (d2e_base_dr2 * s) + (d2s_dr2 * e_base));
          double cse_p523_t2_invr44101 = (1.0 / (r));
          double cse_p11_t1_g35112 = (cse_p523_t2_invr44101 * cse_p523_t2_invr44101 * d2e_dr2);
          double cse_p523_t1_invr44100 = (1.0 / ((r * r * r)));
          double cse_p11_t1_g35113 = (cse_p523_t1_invr44100 * de_dr);
          double cse_p11_t1_g35117 = (r * r);
          double cse_p14_t35_g35364 = (cse_p11_t1_g35117 + (-(cse_p11_t54_g35049)));
          double cse_p14_t36_g35365 = (cse_p11_t1_g35117 + (-(cse_p11_t55_g35050)));
          double cse_p14_t37_g35366 = (cse_p11_t1_g35117 + (-(cse_p11_t56_g35051)));
          double cse_p13_t27_g35262 = (cse_p11_t1_g35113 * cse_p14_t35_g35364);
          double cse_p13_t28_g35263 = (cse_p11_t1_g35113 * cse_p14_t36_g35365);
          double cse_p13_t29_g35264 = (cse_p11_t1_g35113 * cse_p14_t37_g35366);
          double cse_p13_t84_g35319 = (-(cse_p11_t1_g35117));
          double cse_p14_t1_g35383 = (cse_p11_t1_g35112 * dx);
          double h_x1_x1 = (cse_p13_t27_g35262 + (cse_p14_t1_g35383 * dx));
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 0, h_x1_x1);
          double cse_p14_t1_g35384 = (cse_p11_t1_g35113 * dy);
          double cse_p16_t12_g35432 = (cse_p14_t1_g35384 * dx);
          double cse_p15_t1_g35419 = (cse_p14_t1_g35384 * (-(dx)));
          double h_x1_y1 = (cse_p15_t1_g35419 + (cse_p14_t1_g35383 * dy));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 1, h_x1_y1);
          double cse_p11_t1_g35121 = (cse_p11_t1_g35113 * dz);
          double cse_p13_t36_g35271 = (cse_p11_t1_g35121 * (-(dx)));
          double cse_p13_t37_g35272 = (cse_p11_t1_g35121 * (-(dy)));
          double cse_p13_t38_g35273 = (cse_p11_t1_g35121 * dx);
          double cse_p13_t39_g35274 = (cse_p11_t1_g35121 * dy);
          double h_x1_z1 = (cse_p13_t36_g35271 + (cse_p14_t1_g35383 * dz));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 2, h_x1_z1);
          double cse_p11_t1_g35118 = (cse_p11_t1_g35112 * (-(dx)));
          double h_x1_x2 = ((cse_p11_t1_g35113 * (cse_p11_t54_g35049 + cse_p13_t84_g35319)) + (cse_p11_t1_g35118 * dx));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 0, h_x1_x2);
          double cse_p11_t1_g35120 = (cse_p11_t1_g35112 * (-(dy)));
          double h_x1_y2 = (cse_p16_t12_g35432 + (cse_p11_t1_g35120 * dx));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 1, h_x1_y2);
          double cse_p12_t1_g35234 = (cse_p11_t1_g35112 * (-(dz)));
          double h_x1_z2 = (cse_p13_t38_g35273 + (cse_p12_t1_g35234 * dx));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 2, h_x1_z2);
          double cse_p65_t1_g36442 = (cse_p11_t1_g35112 * dy);
          double h_y1_y1 = (cse_p13_t28_g35263 + (cse_p65_t1_g36442 * dy));
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 1, h_y1_y1);
          double h_y1_z1 = (cse_p13_t37_g35272 + (cse_p65_t1_g36442 * dz));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 2, h_y1_z1);
          double h_y1_x2 = (cse_p16_t12_g35432 + (cse_p11_t1_g35118 * dy));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 0, h_y1_x2);
          double h_y1_y2 = ((cse_p11_t1_g35113 * (cse_p11_t55_g35050 + cse_p13_t84_g35319)) + (cse_p11_t1_g35120 * dy));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 1, h_y1_y2);
          double h_y1_z2 = (cse_p13_t39_g35274 + (cse_p12_t1_g35234 * dy));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 2, h_y1_z2);
          double h_z1_z1 = (cse_p13_t29_g35264 + (cse_p11_t1_g35112 * dz * dz));
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x1, 2, h_z1_z1);
          double h_z1_x2 = (cse_p13_t38_g35273 + (cse_p11_t1_g35118 * dz));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 0, h_z1_x2);
          double h_z1_y2 = (cse_p13_t39_g35274 + (cse_p11_t1_g35120 * dz));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 1, h_z1_y2);
          double h_z1_z2 = ((cse_p11_t1_g35113 * (cse_p11_t56_g35051 + cse_p13_t84_g35319)) + (cse_p12_t1_g35234 * dz));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 2, h_z1_z2);
          double h_x2_x2 = (cse_p13_t27_g35262 + (cse_p11_t1_g35118 * (-(dx))));
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 0, h_x2_x2);
          double h_x2_y2 = (cse_p15_t1_g35419 + (cse_p11_t1_g35118 * (-(dy))));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 1, h_x2_y2);
          double h_x2_z2 = (cse_p13_t36_g35271 + (cse_p11_t1_g35118 * (-(dz))));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 2, h_x2_z2);
          double h_y2_y2 = (cse_p13_t28_g35263 + (cse_p11_t1_g35120 * (-(dy))));
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 1, h_y2_y2);
          double h_y2_z2 = (cse_p13_t37_g35272 + (cse_p11_t1_g35120 * (-(dz))));
          KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 2, h_y2_z2);
          double h_z2_z2 = (cse_p13_t29_g35264 + (cse_p12_t1_g35234 * (-(dz))));
          KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x2, 2, h_z2_z2);
        }
      }
    }
  }
  return energy_added;
}
void hessian_fd(double a, double b, double qq, double invdd, double r_switch, double r_switch2, double r_cut, double r_cut2, double inv_range, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 6;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  const double invh2 = 1.0/((h*h));
  double e0 = 0.0;
  energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 0, i3x1, 0, hval);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 1, i3x1, 1, hval);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 2, i3x1, 2, hval);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 0, i3x2, 0, hval);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 1, i3x2, 1, hval);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x2, 2, hval);
  }
  {
    double saved_i = position[i3x1 + 1];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x1 + 1] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 1, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x1 + 2];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x1 + 2] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 2, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x1 + 2];
    double saved_j = position[i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x1 + 2] = saved_i + h; position[i3x1 + 1] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i; position[i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 2, i3x1, 1, hval);
  }
  {
    double saved_i = position[i3x2 + 0];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 0] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 0, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x2 + 0];
    double saved_j = position[i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 0] = saved_i + h; position[i3x1 + 1] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i; position[i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 0, i3x1, 1, hval);
  }
  {
    double saved_i = position[i3x2 + 0];
    double saved_j = position[i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 0] = saved_i + h; position[i3x1 + 2] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i; position[i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 0, i3x1, 2, hval);
  }
  {
    double saved_i = position[i3x2 + 1];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 1] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 1, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x2 + 1];
    double saved_j = position[i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 1] = saved_i + h; position[i3x1 + 1] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i; position[i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 1, i3x1, 1, hval);
  }
  {
    double saved_i = position[i3x2 + 1];
    double saved_j = position[i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 1] = saved_i + h; position[i3x1 + 2] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i; position[i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 1, i3x1, 2, hval);
  }
  {
    double saved_i = position[i3x2 + 1];
    double saved_j = position[i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 1] = saved_i + h; position[i3x2 + 0] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i; position[i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 1, i3x2, 0, hval);
  }
  {
    double saved_i = position[i3x2 + 2];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 2] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x2 + 2];
    double saved_j = position[i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 2] = saved_i + h; position[i3x1 + 1] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i; position[i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x1, 1, hval);
  }
  {
    double saved_i = position[i3x2 + 2];
    double saved_j = position[i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 2] = saved_i + h; position[i3x1 + 2] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i; position[i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x1, 2, hval);
  }
  {
    double saved_i = position[i3x2 + 2];
    double saved_j = position[i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 2] = saved_i + h; position[i3x2 + 0] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i; position[i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x2, 0, hval);
  }
  {
    double saved_i = position[i3x2 + 2];
    double saved_j = position[i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x2 + 2] = saved_i + h; position[i3x2 + 1] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i; position[i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x2, 1, hval);
  }
}

};
