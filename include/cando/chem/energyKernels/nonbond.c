template <typename HESSIAN>
struct Nonbond {
  static constexpr size_t PositionSize = 6;
  static std::string description() { return "mathkernel-nonbond"; };
double gradient(double a, double b, double qq, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
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
      double cse_p11_t1_g28260 = (-3.0000000000000000     * b);
      double dx = (x1 + (-(x2)));
      double dy = (y1 + (-(y2)));
      double dz = (z1 + (-(z2)));
      double r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      double cse_p1_t1_sqrt28256 = sqrt(r2);
      double cse_p61_t1_invr28272 = (1.0 / (cse_p1_t1_sqrt28256));
      double cse_p113_t1_g28273 = (cse_p61_t1_invr28272 * cse_p61_t1_invr28272);
      double cse_p12_t1_g28271 = (cse_p113_t1_g28273 * cse_p113_t1_g28273);
      double invr6 = (cse_p113_t1_g28273 * cse_p12_t1_g28271);
      double cse_p11_t1_g28270 = (a * invr6);
      double cse_p11_t2_g28261 = (6.0000000000000000     * cse_p11_t1_g28270);
      double cse_p11_t3_g28262 = (cse_p11_t1_g28260 + cse_p11_t2_g28261);
      double e_lj = ((cse_p11_t1_g28270 * invr6) + (-((b * invr6))));
      double e_coul = (cse_p61_t1_invr28272 * qq);
      double energy = (e_coul + e_lj);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t1_g28263 = (cse_p113_t1_g28273 * cse_p61_t1_invr28272);
      double cse_p11_t1_g28264 = (cse_p11_t1_g28263 * cse_p11_t3_g28262 * cse_p12_t1_g28271 * cse_p61_t1_invr28272);
      double cse_p11_t1_g28265 = (cse_p11_t1_g28263 * qq);
      double cse_p11_t1_g28266 = (-(cse_p11_t1_g28265));
      double cse_p11_t1_g28267 = (-2.0000000000000000     * cse_p11_t1_g28264);
      double g_x1 = ((cse_p11_t1_g28266 * dx) + (cse_p11_t1_g28267 * dx));
      KernelGradientAcc(i3x1, 0, g_x1);
      double g_y1 = ((cse_p11_t1_g28266 * dy) + (cse_p11_t1_g28267 * dy));
      KernelGradientAcc(i3x1, 1, g_y1);
      double g_z1 = ((cse_p11_t1_g28266 * dz) + (cse_p11_t1_g28267 * dz));
      KernelGradientAcc(i3x1, 2, g_z1);
      double cse_p11_t1_g28269 = (2.0000000000000000     * cse_p11_t1_g28264);
      double g_x2 = ((cse_p11_t1_g28265 * dx) + (cse_p11_t1_g28269 * dx));
      KernelGradientAcc(i3x2, 0, g_x2);
      double g_y2 = ((cse_p11_t1_g28265 * dy) + (cse_p11_t1_g28269 * dy));
      KernelGradientAcc(i3x2, 1, g_y2);
      double g_z2 = ((cse_p11_t1_g28265 * dz) + (cse_p11_t1_g28269 * dz));
      KernelGradientAcc(i3x2, 2, g_z2);
    }
  }
  return energy_added;
}
void gradient_fd(double a, double b, double qq, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 6;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  double e0 = 0.0;
  energy(a, b, qq, i3x1, i3x2, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
}

double energy(double a, double b, double qq, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
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
      double r = sqrt(r2);
      double invr = (1.0 / (r));
      double invr2 = (invr * invr);
      double invr6 = (invr2 * invr2 * invr2);
      double e_lj = ((a * invr6 * invr6) + (-((b * invr6))));
      double e_coul = (invr * qq);
      double energy = (e_coul + e_lj);
      energy_added += energy;
      *energy_accumulate += energy;
    }
  }
  return energy_added;
}
void energy_fd(double a, double b, double qq, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  energy(a, b, qq, i3x1, i3x2, position, energy_accumulate, force, hessian, dvec, hdvec);
}

double hessian(double a, double b, double qq, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
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
      double cse_p11_t5_g28282 = (-3.0000000000000000     * b);
      double dx = (x1 + (-(x2)));
      double cse_p11_t11_g28288 = (dx * dx);
      double cse_p11_t28_g28305 = (-(cse_p11_t11_g28288));
      double cse_p11_t33_g28310 = (-(dx));
      double dy = (y1 + (-(y2)));
      double cse_p11_t15_g28292 = (dy * dy);
      double cse_p11_t29_g28306 = (-(cse_p11_t15_g28292));
      double cse_p11_t34_g28311 = (-(dy));
      double dz = (z1 + (-(z2)));
      double cse_p11_t19_g28296 = (dz * dz);
      double cse_p11_t30_g28307 = (-(cse_p11_t19_g28296));
      double cse_p11_t35_g28312 = (-(dz));
      double r2 = (cse_p11_t11_g28288 + cse_p11_t15_g28292 + cse_p11_t19_g28296);
      double cse_p1_t1_sqrt28274 = sqrt(r2);
      double cse_p62_t1_g28337 = (cse_p1_t1_sqrt28274 * cse_p1_t1_sqrt28274);
      double cse_p61_t1_invr28330 = (1.0 / (cse_p1_t1_sqrt28274));
      double cse_p114_t1_g28350 = (cse_p61_t1_invr28330 * cse_p61_t1_invr28330);
      double cse_p61_t2_invr28331 = (1.0 / (((cse_p1_t1_sqrt28274) * (cse_p1_t1_sqrt28274))));
      double cse_p11_t1_g28329 = (cse_p61_t2_invr28331 * qq);
      double cse_p61_t5_invr28334 = (1.0 / ((cse_p1_t1_sqrt28274 * cse_p1_t1_sqrt28274 * cse_p1_t1_sqrt28274)));
      double cse_p61_t4_invr28333 = (1.0 / (pow(cse_p1_t1_sqrt28274, 7)));
      double cse_p218_t1_g28352 = (6.0000000000000000     * b * cse_p61_t4_invr28333);
      double cse_p61_t3_invr28332 = (1.0 / (pow(cse_p1_t1_sqrt28274, 13)));
      double cse_p166_t1_g28351 = (-12.000000000000000     * a * cse_p61_t3_invr28332);
      double cse_p11_t1_g28325 = (cse_p61_t5_invr28334 * (cse_p166_t1_g28351 + cse_p218_t1_g28352 + (-(cse_p11_t1_g28329))));
      double cse_p11_t1_g28278 = (cse_p11_t1_g28325 * (cse_p11_t28_g28305 + cse_p62_t1_g28337));
      double cse_p11_t2_g28279 = (cse_p11_t1_g28325 * (cse_p11_t29_g28306 + cse_p62_t1_g28337));
      double cse_p11_t3_g28280 = (cse_p11_t1_g28325 * (cse_p11_t30_g28307 + cse_p62_t1_g28337));
      double cse_p61_t6_invr28335 = (1.0 / (pow(cse_p1_t1_sqrt28274, 8)));
      double cse_p11_t6_g28283 = (-42.000000000000000     * b * cse_p61_t6_invr28335);
      double cse_p61_t7_invr28336 = (1.0 / (pow(cse_p1_t1_sqrt28274, 14)));
      double cse_p11_t7_g28284 = (156.00000000000000     * a * cse_p61_t7_invr28336);
      double cse_p11_t8_g28285 = (2.0000000000000000     * cse_p61_t5_invr28334 * qq);
      double cse_p11_t12_g28289 = (cse_p11_t1_g28325 * dx * dy);
      double cse_p63_t1_g28343 = (cse_p11_t1_g28325 * dz);
      double cse_p11_t13_g28290 = (cse_p63_t1_g28343 * dx);
      double cse_p11_t14_g28291 = (cse_p11_t1_g28325 * cse_p11_t33_g28310 * dy);
      double cse_p11_t16_g28293 = (cse_p63_t1_g28343 * dy);
      double cse_p11_t17_g28294 = (cse_p11_t33_g28310 * cse_p63_t1_g28343);
      double cse_p11_t18_g28295 = (cse_p11_t34_g28311 * cse_p63_t1_g28343);
      double cse_p11_t22_g28299 = (cse_p166_t1_g28351 + cse_p218_t1_g28352 + (-(cse_p11_t1_g28329)));
      double cse_p11_t24_g28301 = (cse_p11_t6_g28283 + cse_p11_t7_g28284 + cse_p11_t8_g28285);
      double cse_p11_t25_g28302 = (cse_p11_t28_g28305 + cse_p62_t1_g28337);
      double cse_p11_t26_g28303 = (cse_p11_t29_g28306 + cse_p62_t1_g28337);
      double cse_p11_t27_g28304 = (cse_p11_t30_g28307 + cse_p62_t1_g28337);
      double cse_p11_t31_g28308 = (-(cse_p11_t1_g28329));
      double cse_p11_t32_g28309 = (-(cse_p62_t1_g28337));
      double cse_p11_t42_g28319 = pow(cse_p1_t1_sqrt28274, 13);
      double cse_p11_t43_g28320 = pow(cse_p1_t1_sqrt28274, 14);
      double cse_p11_t44_g28321 = ((cse_p1_t1_sqrt28274) * (cse_p1_t1_sqrt28274));
      double cse_p11_t45_g28322 = (cse_p1_t1_sqrt28274 * cse_p1_t1_sqrt28274 * cse_p1_t1_sqrt28274);
      double cse_p11_t46_g28323 = pow(cse_p1_t1_sqrt28274, 7);
      double cse_p11_t47_g28324 = pow(cse_p1_t1_sqrt28274, 8);
      double invr6 = (cse_p114_t1_g28350 * cse_p114_t1_g28350 * cse_p114_t1_g28350);
      double cse_p11_t9_g28286 = (6.0000000000000000     * a * invr6);
      double cse_p11_t23_g28300 = (cse_p11_t5_g28282 + cse_p11_t9_g28286);
      double e_lj = ((a * invr6 * invr6) + (-((b * invr6))));
      double e_coul = (cse_p61_t1_invr28330 * qq);
      double energy = (e_coul + e_lj);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t1_g28327 = (cse_p114_t1_g28350 * cse_p61_t1_invr28330);
      double cse_p63_t1_g28348 = (cse_p11_t1_g28327 * qq);
      double cse_p63_t1_g28338 = (cse_p63_t1_g28348 * dx);
      double cse_p63_t2_g28339 = (cse_p63_t1_g28348 * dy);
      double cse_p63_t3_g28340 = (cse_p63_t1_g28348 * dz);
      double cse_p11_t1_g28328 = (cse_p114_t1_g28350 * cse_p114_t1_g28350 * cse_p11_t1_g28327 * cse_p11_t23_g28300 * cse_p61_t1_invr28330);
      double cse_p63_t1_g28345 = (-2.0000000000000000     * cse_p11_t1_g28328);
      double g_x1 = ((cse_p63_t1_g28345 * dx) + (-(cse_p63_t1_g28338)));
      KernelGradientAcc(i3x1, 0, g_x1);
      double g_y1 = ((cse_p63_t1_g28345 * dy) + (-(cse_p63_t2_g28339)));
      KernelGradientAcc(i3x1, 1, g_y1);
      double g_z1 = ((cse_p63_t1_g28345 * dz) + (-(cse_p63_t3_g28340)));
      KernelGradientAcc(i3x1, 2, g_z1);
      double cse_p63_t1_g28346 = (2.0000000000000000     * cse_p11_t1_g28328);
      double g_x2 = (cse_p63_t1_g28338 + (cse_p63_t1_g28346 * dx));
      KernelGradientAcc(i3x2, 0, g_x2);
      double g_y2 = (cse_p63_t2_g28339 + (cse_p63_t1_g28346 * dy));
      KernelGradientAcc(i3x2, 1, g_y2);
      double g_z2 = (cse_p63_t3_g28340 + (cse_p63_t1_g28346 * dz));
      KernelGradientAcc(i3x2, 2, g_z2);
      double cse_p11_t1_g28326 = (cse_p114_t1_g28350 * cse_p11_t24_g28301);
      double cse_p63_t1_g28347 = (cse_p11_t1_g28326 * dx);
      double h_x1_x1 = (cse_p11_t1_g28278 + (cse_p63_t1_g28347 * dx));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 0, h_x1_x1);
      double h_x1_y1 = (cse_p11_t14_g28291 + (cse_p63_t1_g28347 * dy));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 1, h_x1_y1);
      double h_x1_z1 = (cse_p11_t17_g28294 + (cse_p63_t1_g28347 * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 2, h_x1_z1);
      double cse_p63_t1_g28341 = (cse_p11_t1_g28326 * cse_p11_t33_g28310);
      double cse_p63_t1_g28349 = (cse_p11_t22_g28299 * cse_p61_t5_invr28334);
      double h_x1_x2 = ((cse_p63_t1_g28341 * dx) + (cse_p63_t1_g28349 * (cse_p11_t11_g28288 + cse_p11_t32_g28309)));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 0, h_x1_x2);
      double cse_p63_t1_g28342 = (cse_p11_t1_g28326 * cse_p11_t34_g28311);
      double h_x1_y2 = (cse_p11_t12_g28289 + (cse_p63_t1_g28342 * dx));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 1, h_x1_y2);
      double cse_p63_t1_g28344 = (cse_p11_t1_g28326 * cse_p11_t35_g28312);
      double h_x1_z2 = (cse_p11_t13_g28290 + (cse_p63_t1_g28344 * dx));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 2, h_x1_z2);
      double h_y1_y1 = (cse_p11_t2_g28279 + (cse_p11_t1_g28326 * dy * dy));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 1, h_y1_y1);
      double h_y1_z1 = (cse_p11_t18_g28295 + (cse_p11_t1_g28326 * dy * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 2, h_y1_z1);
      double h_y1_x2 = (cse_p11_t12_g28289 + (cse_p63_t1_g28341 * dy));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 0, h_y1_x2);
      double h_y1_y2 = ((cse_p63_t1_g28342 * dy) + (cse_p63_t1_g28349 * (cse_p11_t15_g28292 + cse_p11_t32_g28309)));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 1, h_y1_y2);
      double h_y1_z2 = (cse_p11_t16_g28293 + (cse_p63_t1_g28344 * dy));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 2, h_y1_z2);
      double h_z1_z1 = (cse_p11_t3_g28280 + (cse_p11_t1_g28326 * dz * dz));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x1, 2, h_z1_z1);
      double h_z1_x2 = (cse_p11_t13_g28290 + (cse_p63_t1_g28341 * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 0, h_z1_x2);
      double h_z1_y2 = (cse_p11_t16_g28293 + (cse_p63_t1_g28342 * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 1, h_z1_y2);
      double h_z1_z2 = ((cse_p63_t1_g28344 * dz) + (cse_p63_t1_g28349 * (cse_p11_t19_g28296 + cse_p11_t32_g28309)));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 2, h_z1_z2);
      double h_x2_x2 = (cse_p11_t1_g28278 + (cse_p11_t33_g28310 * cse_p63_t1_g28341));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 0, h_x2_x2);
      double h_x2_y2 = (cse_p11_t14_g28291 + (cse_p11_t34_g28311 * cse_p63_t1_g28341));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 1, h_x2_y2);
      double h_x2_z2 = (cse_p11_t17_g28294 + (cse_p11_t35_g28312 * cse_p63_t1_g28341));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 2, h_x2_z2);
      double h_y2_y2 = (cse_p11_t2_g28279 + (cse_p11_t34_g28311 * cse_p63_t1_g28342));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 1, h_y2_y2);
      double h_y2_z2 = (cse_p11_t18_g28295 + (cse_p11_t35_g28312 * cse_p63_t1_g28342));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 2, h_y2_z2);
      double h_z2_z2 = (cse_p11_t3_g28280 + (cse_p11_t35_g28312 * cse_p63_t1_g28344));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x2, 2, h_z2_z2);
    }
  }
  return energy_added;
}
void hessian_fd(double a, double b, double qq, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 6;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  const double invh2 = 1.0/((h*h));
  double e0 = 0.0;
  energy(a, b, qq, i3x1, i3x2, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 0, i3x1, 0, hval);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 1, i3x1, 1, hval);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 2, i3x1, 2, hval);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 0, i3x2, 0, hval);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 1, i3x2, 1, hval);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
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
    energy(a, b, qq, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(a, b, qq, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(a, b, qq, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(a, b, qq, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(a, b, qq, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(a, b, qq, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(a, b, qq, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(a, b, qq, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(a, b, qq, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(a, b, qq, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(a, b, qq, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(a, b, qq, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(a, b, qq, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(a, b, qq, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(a, b, qq, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(a, b, qq, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i; position[i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x2, 1, hval);
  }
}

};
