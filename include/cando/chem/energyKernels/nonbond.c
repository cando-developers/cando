template <typename HESSIAN>
struct Nonbond {
  static constexpr size_t PositionSize = 6;
  static std::string description() { return "mathkernel-nonbond"; };
double gradient(double a, double b, double qq, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double cse_p113_t1_g35452;
  double cse_p11_t1_g35439;
  double cse_p11_t1_g35442;
  double cse_p11_t1_g35443;
  double cse_p11_t1_g35444;
  double cse_p11_t1_g35445;
  double cse_p11_t1_g35446;
  double cse_p11_t1_g35448;
  double cse_p11_t1_g35449;
  double cse_p11_t2_g35440;
  double cse_p11_t3_g35441;
  double cse_p12_t1_g35450;
  double cse_p1_t1_sqrt35435;
  double cse_p61_t1_invr35451;
  double dx;
  double dy;
  double dz;
  double energy;
  double energy_added;
  double e_coul;
  double e_lj;
  double g_x1;
  double g_x2;
  double g_y1;
  double g_y2;
  double g_z1;
  double g_z2;
  double invr6;
  double r2;
  energy_added = 0.0000000000000000    ;
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
      cse_p11_t1_g35439 = (-3.0000000000000000     * b);
      dx = (x1 + (-(x2)));
      dy = (y1 + (-(y2)));
      dz = (z1 + (-(z2)));
      r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      cse_p1_t1_sqrt35435 = sqrt(r2);
      cse_p61_t1_invr35451 = (1.0 / (cse_p1_t1_sqrt35435));
      cse_p113_t1_g35452 = (cse_p61_t1_invr35451 * cse_p61_t1_invr35451);
      cse_p12_t1_g35450 = (cse_p113_t1_g35452 * cse_p113_t1_g35452);
      invr6 = (cse_p113_t1_g35452 * cse_p12_t1_g35450);
      cse_p11_t1_g35449 = (a * invr6);
      cse_p11_t2_g35440 = (6.0000000000000000     * cse_p11_t1_g35449);
      cse_p11_t3_g35441 = (cse_p11_t1_g35439 + cse_p11_t2_g35440);
      e_lj = ((cse_p11_t1_g35449 * invr6) + (-((b * invr6))));
      e_coul = (cse_p61_t1_invr35451 * qq);
      energy = (e_coul + e_lj);
      energy_added += energy;
      *energy_accumulate += energy;
      cse_p11_t1_g35442 = (cse_p113_t1_g35452 * cse_p61_t1_invr35451);
      cse_p11_t1_g35443 = (cse_p11_t1_g35442 * cse_p11_t3_g35441 * cse_p12_t1_g35450 * cse_p61_t1_invr35451);
      cse_p11_t1_g35444 = (cse_p11_t1_g35442 * qq);
      cse_p11_t1_g35445 = (-(cse_p11_t1_g35444));
      cse_p11_t1_g35446 = (-2.0000000000000000     * cse_p11_t1_g35443);
      g_x1 = ((cse_p11_t1_g35445 * dx) + (cse_p11_t1_g35446 * dx));
      KernelGradientAcc(i3x1, 0, g_x1);
      g_y1 = ((cse_p11_t1_g35445 * dy) + (cse_p11_t1_g35446 * dy));
      KernelGradientAcc(i3x1, 1, g_y1);
      g_z1 = ((cse_p11_t1_g35445 * dz) + (cse_p11_t1_g35446 * dz));
      KernelGradientAcc(i3x1, 2, g_z1);
      cse_p11_t1_g35448 = (2.0000000000000000     * cse_p11_t1_g35443);
      g_x2 = ((cse_p11_t1_g35444 * dx) + (cse_p11_t1_g35448 * dx));
      KernelGradientAcc(i3x2, 0, g_x2);
      g_y2 = ((cse_p11_t1_g35444 * dy) + (cse_p11_t1_g35448 * dy));
      KernelGradientAcc(i3x2, 1, g_y2);
      g_z2 = ((cse_p11_t1_g35444 * dz) + (cse_p11_t1_g35448 * dz));
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
  double dx;
  double dy;
  double dz;
  double energy;
  double energy_added;
  double e_coul;
  double e_lj;
  double invr;
  double invr2;
  double invr6;
  double r;
  double r2;
  energy_added = 0.0000000000000000    ;
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
      dx = (x1 + (-(x2)));
      dy = (y1 + (-(y2)));
      dz = (z1 + (-(z2)));
      r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      r = sqrt(r2);
      invr = (1.0 / (r));
      invr2 = (invr * invr);
      invr6 = (invr2 * invr2 * invr2);
      e_lj = ((a * invr6 * invr6) + (-((b * invr6))));
      e_coul = (invr * qq);
      energy = (e_coul + e_lj);
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
  double cse_p114_t1_g35529;
  double cse_p11_t11_g35467;
  double cse_p11_t12_g35468;
  double cse_p11_t13_g35469;
  double cse_p11_t14_g35470;
  double cse_p11_t15_g35471;
  double cse_p11_t16_g35472;
  double cse_p11_t17_g35473;
  double cse_p11_t18_g35474;
  double cse_p11_t19_g35475;
  double cse_p11_t1_g35457;
  double cse_p11_t1_g35504;
  double cse_p11_t1_g35505;
  double cse_p11_t1_g35506;
  double cse_p11_t1_g35507;
  double cse_p11_t1_g35508;
  double cse_p11_t22_g35478;
  double cse_p11_t23_g35479;
  double cse_p11_t24_g35480;
  double cse_p11_t25_g35481;
  double cse_p11_t26_g35482;
  double cse_p11_t27_g35483;
  double cse_p11_t28_g35484;
  double cse_p11_t29_g35485;
  double cse_p11_t2_g35458;
  double cse_p11_t30_g35486;
  double cse_p11_t31_g35487;
  double cse_p11_t32_g35488;
  double cse_p11_t33_g35489;
  double cse_p11_t34_g35490;
  double cse_p11_t35_g35491;
  double cse_p11_t3_g35459;
  double cse_p11_t42_g35498;
  double cse_p11_t43_g35499;
  double cse_p11_t44_g35500;
  double cse_p11_t45_g35501;
  double cse_p11_t46_g35502;
  double cse_p11_t47_g35503;
  double cse_p11_t5_g35461;
  double cse_p11_t6_g35462;
  double cse_p11_t7_g35463;
  double cse_p11_t8_g35464;
  double cse_p11_t9_g35465;
  double cse_p166_t1_g35530;
  double cse_p1_t1_sqrt35453;
  double cse_p218_t1_g35531;
  double cse_p61_t1_invr35509;
  double cse_p61_t2_invr35510;
  double cse_p61_t3_invr35511;
  double cse_p61_t4_invr35512;
  double cse_p61_t5_invr35513;
  double cse_p61_t6_invr35514;
  double cse_p61_t7_invr35515;
  double cse_p62_t1_g35516;
  double cse_p63_t1_g35517;
  double cse_p63_t1_g35520;
  double cse_p63_t1_g35521;
  double cse_p63_t1_g35522;
  double cse_p63_t1_g35523;
  double cse_p63_t1_g35524;
  double cse_p63_t1_g35525;
  double cse_p63_t1_g35526;
  double cse_p63_t1_g35527;
  double cse_p63_t1_g35528;
  double cse_p63_t2_g35518;
  double cse_p63_t3_g35519;
  double dx;
  double dy;
  double dz;
  double energy;
  double energy_added;
  double e_coul;
  double e_lj;
  double g_x1;
  double g_x2;
  double g_y1;
  double g_y2;
  double g_z1;
  double g_z2;
  double h_x1_x1;
  double h_x1_x2;
  double h_x1_y1;
  double h_x1_y2;
  double h_x1_z1;
  double h_x1_z2;
  double h_x2_x2;
  double h_x2_y2;
  double h_x2_z2;
  double h_y1_x2;
  double h_y1_y1;
  double h_y1_y2;
  double h_y1_z1;
  double h_y1_z2;
  double h_y2_y2;
  double h_y2_z2;
  double h_z1_x2;
  double h_z1_y2;
  double h_z1_z1;
  double h_z1_z2;
  double h_z2_z2;
  double invr6;
  double r2;
  energy_added = 0.0000000000000000    ;
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
      cse_p11_t5_g35461 = (-3.0000000000000000     * b);
      dx = (x1 + (-(x2)));
      cse_p11_t11_g35467 = (dx * dx);
      cse_p11_t28_g35484 = (-(cse_p11_t11_g35467));
      cse_p11_t33_g35489 = (-(dx));
      dy = (y1 + (-(y2)));
      cse_p11_t15_g35471 = (dy * dy);
      cse_p11_t29_g35485 = (-(cse_p11_t15_g35471));
      cse_p11_t34_g35490 = (-(dy));
      dz = (z1 + (-(z2)));
      cse_p11_t19_g35475 = (dz * dz);
      cse_p11_t30_g35486 = (-(cse_p11_t19_g35475));
      cse_p11_t35_g35491 = (-(dz));
      r2 = (cse_p11_t11_g35467 + cse_p11_t15_g35471 + cse_p11_t19_g35475);
      cse_p1_t1_sqrt35453 = sqrt(r2);
      cse_p62_t1_g35516 = (cse_p1_t1_sqrt35453 * cse_p1_t1_sqrt35453);
      cse_p61_t1_invr35509 = (1.0 / (cse_p1_t1_sqrt35453));
      cse_p114_t1_g35529 = (cse_p61_t1_invr35509 * cse_p61_t1_invr35509);
      cse_p61_t2_invr35510 = (1.0 / (((cse_p1_t1_sqrt35453) * (cse_p1_t1_sqrt35453))));
      cse_p11_t1_g35508 = (cse_p61_t2_invr35510 * qq);
      cse_p61_t5_invr35513 = (1.0 / ((cse_p1_t1_sqrt35453 * cse_p1_t1_sqrt35453 * cse_p1_t1_sqrt35453)));
      cse_p61_t4_invr35512 = (1.0 / (pow(cse_p1_t1_sqrt35453, 7)));
      cse_p218_t1_g35531 = (6.0000000000000000     * b * cse_p61_t4_invr35512);
      cse_p61_t3_invr35511 = (1.0 / (pow(cse_p1_t1_sqrt35453, 13)));
      cse_p166_t1_g35530 = (-12.000000000000000     * a * cse_p61_t3_invr35511);
      cse_p11_t1_g35504 = (cse_p61_t5_invr35513 * (cse_p166_t1_g35530 + cse_p218_t1_g35531 + (-(cse_p11_t1_g35508))));
      cse_p11_t1_g35457 = (cse_p11_t1_g35504 * (cse_p11_t28_g35484 + cse_p62_t1_g35516));
      cse_p11_t2_g35458 = (cse_p11_t1_g35504 * (cse_p11_t29_g35485 + cse_p62_t1_g35516));
      cse_p11_t3_g35459 = (cse_p11_t1_g35504 * (cse_p11_t30_g35486 + cse_p62_t1_g35516));
      cse_p61_t6_invr35514 = (1.0 / (pow(cse_p1_t1_sqrt35453, 8)));
      cse_p11_t6_g35462 = (-42.000000000000000     * b * cse_p61_t6_invr35514);
      cse_p61_t7_invr35515 = (1.0 / (pow(cse_p1_t1_sqrt35453, 14)));
      cse_p11_t7_g35463 = (156.00000000000000     * a * cse_p61_t7_invr35515);
      cse_p11_t8_g35464 = (2.0000000000000000     * cse_p61_t5_invr35513 * qq);
      cse_p11_t12_g35468 = (cse_p11_t1_g35504 * dx * dy);
      cse_p63_t1_g35522 = (cse_p11_t1_g35504 * dz);
      cse_p11_t13_g35469 = (cse_p63_t1_g35522 * dx);
      cse_p11_t14_g35470 = (cse_p11_t1_g35504 * cse_p11_t33_g35489 * dy);
      cse_p11_t16_g35472 = (cse_p63_t1_g35522 * dy);
      cse_p11_t17_g35473 = (cse_p11_t33_g35489 * cse_p63_t1_g35522);
      cse_p11_t18_g35474 = (cse_p11_t34_g35490 * cse_p63_t1_g35522);
      cse_p11_t22_g35478 = (cse_p166_t1_g35530 + cse_p218_t1_g35531 + (-(cse_p11_t1_g35508)));
      cse_p11_t24_g35480 = (cse_p11_t6_g35462 + cse_p11_t7_g35463 + cse_p11_t8_g35464);
      cse_p11_t25_g35481 = (cse_p11_t28_g35484 + cse_p62_t1_g35516);
      cse_p11_t26_g35482 = (cse_p11_t29_g35485 + cse_p62_t1_g35516);
      cse_p11_t27_g35483 = (cse_p11_t30_g35486 + cse_p62_t1_g35516);
      cse_p11_t31_g35487 = (-(cse_p11_t1_g35508));
      cse_p11_t32_g35488 = (-(cse_p62_t1_g35516));
      cse_p11_t42_g35498 = pow(cse_p1_t1_sqrt35453, 13);
      cse_p11_t43_g35499 = pow(cse_p1_t1_sqrt35453, 14);
      cse_p11_t44_g35500 = ((cse_p1_t1_sqrt35453) * (cse_p1_t1_sqrt35453));
      cse_p11_t45_g35501 = (cse_p1_t1_sqrt35453 * cse_p1_t1_sqrt35453 * cse_p1_t1_sqrt35453);
      cse_p11_t46_g35502 = pow(cse_p1_t1_sqrt35453, 7);
      cse_p11_t47_g35503 = pow(cse_p1_t1_sqrt35453, 8);
      invr6 = (cse_p114_t1_g35529 * cse_p114_t1_g35529 * cse_p114_t1_g35529);
      cse_p11_t9_g35465 = (6.0000000000000000     * a * invr6);
      cse_p11_t23_g35479 = (cse_p11_t5_g35461 + cse_p11_t9_g35465);
      e_lj = ((a * invr6 * invr6) + (-((b * invr6))));
      e_coul = (cse_p61_t1_invr35509 * qq);
      energy = (e_coul + e_lj);
      energy_added += energy;
      *energy_accumulate += energy;
      cse_p11_t1_g35506 = (cse_p114_t1_g35529 * cse_p61_t1_invr35509);
      cse_p63_t1_g35527 = (cse_p11_t1_g35506 * qq);
      cse_p63_t1_g35517 = (cse_p63_t1_g35527 * dx);
      cse_p63_t2_g35518 = (cse_p63_t1_g35527 * dy);
      cse_p63_t3_g35519 = (cse_p63_t1_g35527 * dz);
      cse_p11_t1_g35507 = (cse_p114_t1_g35529 * cse_p114_t1_g35529 * cse_p11_t1_g35506 * cse_p11_t23_g35479 * cse_p61_t1_invr35509);
      cse_p63_t1_g35524 = (-2.0000000000000000     * cse_p11_t1_g35507);
      g_x1 = ((cse_p63_t1_g35524 * dx) + (-(cse_p63_t1_g35517)));
      KernelGradientAcc(i3x1, 0, g_x1);
      g_y1 = ((cse_p63_t1_g35524 * dy) + (-(cse_p63_t2_g35518)));
      KernelGradientAcc(i3x1, 1, g_y1);
      g_z1 = ((cse_p63_t1_g35524 * dz) + (-(cse_p63_t3_g35519)));
      KernelGradientAcc(i3x1, 2, g_z1);
      cse_p63_t1_g35525 = (2.0000000000000000     * cse_p11_t1_g35507);
      g_x2 = (cse_p63_t1_g35517 + (cse_p63_t1_g35525 * dx));
      KernelGradientAcc(i3x2, 0, g_x2);
      g_y2 = (cse_p63_t2_g35518 + (cse_p63_t1_g35525 * dy));
      KernelGradientAcc(i3x2, 1, g_y2);
      g_z2 = (cse_p63_t3_g35519 + (cse_p63_t1_g35525 * dz));
      KernelGradientAcc(i3x2, 2, g_z2);
      cse_p11_t1_g35505 = (cse_p114_t1_g35529 * cse_p11_t24_g35480);
      cse_p63_t1_g35526 = (cse_p11_t1_g35505 * dx);
      h_x1_x1 = (cse_p11_t1_g35457 + (cse_p63_t1_g35526 * dx));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 0, h_x1_x1);
      h_x1_y1 = (cse_p11_t14_g35470 + (cse_p63_t1_g35526 * dy));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 1, h_x1_y1);
      h_x1_z1 = (cse_p11_t17_g35473 + (cse_p63_t1_g35526 * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 2, h_x1_z1);
      cse_p63_t1_g35520 = (cse_p11_t1_g35505 * cse_p11_t33_g35489);
      cse_p63_t1_g35528 = (cse_p11_t22_g35478 * cse_p61_t5_invr35513);
      h_x1_x2 = ((cse_p63_t1_g35520 * dx) + (cse_p63_t1_g35528 * (cse_p11_t11_g35467 + cse_p11_t32_g35488)));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 0, h_x1_x2);
      cse_p63_t1_g35521 = (cse_p11_t1_g35505 * cse_p11_t34_g35490);
      h_x1_y2 = (cse_p11_t12_g35468 + (cse_p63_t1_g35521 * dx));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 1, h_x1_y2);
      cse_p63_t1_g35523 = (cse_p11_t1_g35505 * cse_p11_t35_g35491);
      h_x1_z2 = (cse_p11_t13_g35469 + (cse_p63_t1_g35523 * dx));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 2, h_x1_z2);
      h_y1_y1 = (cse_p11_t2_g35458 + (cse_p11_t1_g35505 * dy * dy));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 1, h_y1_y1);
      h_y1_z1 = (cse_p11_t18_g35474 + (cse_p11_t1_g35505 * dy * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 2, h_y1_z1);
      h_y1_x2 = (cse_p11_t12_g35468 + (cse_p63_t1_g35520 * dy));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 0, h_y1_x2);
      h_y1_y2 = ((cse_p63_t1_g35521 * dy) + (cse_p63_t1_g35528 * (cse_p11_t15_g35471 + cse_p11_t32_g35488)));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 1, h_y1_y2);
      h_y1_z2 = (cse_p11_t16_g35472 + (cse_p63_t1_g35523 * dy));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 2, h_y1_z2);
      h_z1_z1 = (cse_p11_t3_g35459 + (cse_p11_t1_g35505 * dz * dz));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x1, 2, h_z1_z1);
      h_z1_x2 = (cse_p11_t13_g35469 + (cse_p63_t1_g35520 * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 0, h_z1_x2);
      h_z1_y2 = (cse_p11_t16_g35472 + (cse_p63_t1_g35521 * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 1, h_z1_y2);
      h_z1_z2 = ((cse_p63_t1_g35523 * dz) + (cse_p63_t1_g35528 * (cse_p11_t19_g35475 + cse_p11_t32_g35488)));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 2, h_z1_z2);
      h_x2_x2 = (cse_p11_t1_g35457 + (cse_p11_t33_g35489 * cse_p63_t1_g35520));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 0, h_x2_x2);
      h_x2_y2 = (cse_p11_t14_g35470 + (cse_p11_t34_g35490 * cse_p63_t1_g35520));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 1, h_x2_y2);
      h_x2_z2 = (cse_p11_t17_g35473 + (cse_p11_t35_g35491 * cse_p63_t1_g35520));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 2, h_x2_z2);
      h_y2_y2 = (cse_p11_t2_g35458 + (cse_p11_t34_g35490 * cse_p63_t1_g35521));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 1, h_y2_y2);
      h_y2_z2 = (cse_p11_t18_g35474 + (cse_p11_t35_g35491 * cse_p63_t1_g35521));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 2, h_y2_z2);
      h_z2_z2 = (cse_p11_t3_g35459 + (cse_p11_t35_g35491 * cse_p63_t1_g35523));
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
