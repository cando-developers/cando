template <typename HESSIAN>
struct Chiral_Restraint {
  static constexpr size_t PositionSize = 12;
  static std::string description() { return "mathkernel-chiral_restraint"; };
double energy(double k, double co, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double cse_p11_t1_g51299;
  double cse_p11_t2_g51300;
  double cse_p11_t3_g51301;
  double cx;
  double cy;
  double cz;
  double denom;
  double dx13;
  double dx23;
  double dx43;
  double dy13;
  double dy23;
  double dy43;
  double dz13;
  double dz23;
  double dz43;
  double energy;
  double energy_added;
  double q;
  double r13;
  double r13_2;
  double r23;
  double r23_2;
  double r43;
  double r43_2;
  double v;
  energy_added = 0.0000000000000000    ;
  {
    /* !BASE */
    DOUBLE x1 = position[i3x1 + 0];
    DOUBLE y1 = position[i3x1 + 1];
    DOUBLE z1 = position[i3x1 + 2];
    DOUBLE x2 = position[i3x2 + 0];
    DOUBLE y2 = position[i3x2 + 1];
    DOUBLE z2 = position[i3x2 + 2];
    DOUBLE x3 = position[i3x3 + 0];
    DOUBLE y3 = position[i3x3 + 1];
    DOUBLE z3 = position[i3x3 + 2];
    DOUBLE x4 = position[i3x4 + 0];
    DOUBLE y4 = position[i3x4 + 1];
    DOUBLE z4 = position[i3x4 + 2];
    {
      /* !BASE */
      cse_p11_t1_g51299 = (-(x3));
      cse_p11_t2_g51300 = (-(y3));
      cse_p11_t3_g51301 = (-(z3));
      dx13 = (cse_p11_t1_g51299 + x1);
      dy13 = (cse_p11_t2_g51300 + y1);
      dz13 = (cse_p11_t3_g51301 + z1);
      dx23 = (cse_p11_t1_g51299 + x2);
      dy23 = (cse_p11_t2_g51300 + y2);
      dz23 = (cse_p11_t3_g51301 + z2);
      dx43 = (cse_p11_t1_g51299 + x4);
      dy43 = (cse_p11_t2_g51300 + y4);
      dz43 = (cse_p11_t3_g51301 + z4);
      r13_2 = ((dx13 * dx13) + (dy13 * dy13) + (dz13 * dz13));
      r23_2 = ((dx23 * dx23) + (dy23 * dy23) + (dz23 * dz23));
      r43_2 = ((dx43 * dx43) + (dy43 * dy43) + (dz43 * dz43));
      r13 = sqrt(r13_2);
      r23 = sqrt(r23_2);
      r43 = sqrt(r43_2);
      cx = ((dy13 * dz23) + (-((dy23 * dz13))));
      cy = ((dx23 * dz13) + (-((dx13 * dz23))));
      cz = ((dx13 * dy23) + (-((dx23 * dy13))));
      v = ((cx * dx43) + (cy * dy43) + (cz * dz43));
      denom = (r13 * r23 * r43);
      q = (v * (1.0 / (denom)));
      energy = (k * ((co + q) * (co + q) * (co + q)));
      energy_added += energy;
      *energy_accumulate += energy;
    }
  }
  return energy_added;
}
void energy_fd(double k, double co, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  energy(k, co, i3x1, i3x2, i3x3, i3x4, position, energy_accumulate, force, hessian, dvec, hdvec);
}

double gradient(double k, double co, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double cse_p113_t1_g51330;
  double cse_p11_t10_g51319;
  double cse_p11_t11_g51320;
  double cse_p11_t1_g51310;
  double cse_p11_t1_g51321;
  double cse_p11_t1_g51322;
  double cse_p11_t1_g51324;
  double cse_p11_t1_g51325;
  double cse_p11_t1_g51326;
  double cse_p11_t1_g51327;
  double cse_p11_t1_g51328;
  double cse_p11_t1_g51329;
  double cse_p11_t2_g51311;
  double cse_p11_t3_g51312;
  double cse_p11_t4_g51313;
  double cse_p11_t5_g51314;
  double cse_p11_t6_g51315;
  double cse_p11_t7_g51316;
  double cse_p11_t8_g51317;
  double cse_p11_t9_g51318;
  double cse_p1_t1_sqrt51302;
  double cse_p1_t2_invsqrt51303;
  double cse_p1_t3_sqrt51304;
  double cse_p1_t4_invsqrt51305;
  double cse_p1_t5_sqrt51306;
  double cse_p1_t6_invsqrt51307;
  double cse_p1_t7_invr51308;
  double cse_p1_t8_invr251309;
  double cx;
  double cy;
  double cz;
  double denom;
  double dx13;
  double dx23;
  double dx43;
  double dy13;
  double dy23;
  double dy43;
  double dz13;
  double dz23;
  double dz43;
  double energy;
  double energy_added;
  double g_x1;
  double g_x2;
  double g_x3;
  double g_x4;
  double g_y1;
  double g_y2;
  double g_y3;
  double g_y4;
  double g_z1;
  double g_z2;
  double g_z3;
  double g_z4;
  double q;
  double r13_2;
  double r23_2;
  double r43_2;
  double v;
  energy_added = 0.0000000000000000    ;
  {
    /* !BASE */
    DOUBLE x1 = position[i3x1 + 0];
    DOUBLE y1 = position[i3x1 + 1];
    DOUBLE z1 = position[i3x1 + 2];
    DOUBLE x2 = position[i3x2 + 0];
    DOUBLE y2 = position[i3x2 + 1];
    DOUBLE z2 = position[i3x2 + 2];
    DOUBLE x3 = position[i3x3 + 0];
    DOUBLE y3 = position[i3x3 + 1];
    DOUBLE z3 = position[i3x3 + 2];
    DOUBLE x4 = position[i3x4 + 0];
    DOUBLE y4 = position[i3x4 + 1];
    DOUBLE z4 = position[i3x4 + 2];
    {
      /* !BASE */
      cse_p11_t8_g51317 = (-(x3));
      cse_p11_t9_g51318 = (-(y3));
      cse_p11_t10_g51319 = (-(z3));
      dx13 = (cse_p11_t8_g51317 + x1);
      cse_p11_t2_g51311 = (-(dx13));
      dy13 = (cse_p11_t9_g51318 + y1);
      cse_p11_t4_g51313 = (-(dy13));
      dz13 = (cse_p11_t10_g51319 + z1);
      cse_p11_t6_g51315 = (-(dz13));
      dx23 = (cse_p11_t8_g51317 + x2);
      cse_p11_t3_g51312 = (-(dx23));
      dy23 = (cse_p11_t9_g51318 + y2);
      cse_p11_t5_g51314 = (-(dy23));
      dz23 = (cse_p11_t10_g51319 + z2);
      cse_p11_t7_g51316 = (-(dz23));
      dx43 = (cse_p11_t8_g51317 + x4);
      dy43 = (cse_p11_t9_g51318 + y4);
      dz43 = (cse_p11_t10_g51319 + z4);
      r13_2 = ((dx13 * dx13) + (dy13 * dy13) + (dz13 * dz13));
      r23_2 = ((dx23 * dx23) + (dy23 * dy23) + (dz23 * dz23));
      r43_2 = ((dx43 * dx43) + (dy43 * dy43) + (dz43 * dz43));
      cse_p1_t1_sqrt51302 = sqrt(r13_2);
      cse_p1_t2_invsqrt51303 = (1.0 / (cse_p1_t1_sqrt51302));
      cse_p1_t3_sqrt51304 = sqrt(r23_2);
      cse_p1_t4_invsqrt51305 = (1.0 / (cse_p1_t3_sqrt51304));
      cse_p1_t5_sqrt51306 = sqrt(r43_2);
      cse_p1_t6_invsqrt51307 = (1.0 / (cse_p1_t5_sqrt51306));
      cx = ((dy13 * dz23) + (-((dy23 * dz13))));
      cy = ((dx23 * dz13) + (-((dx13 * dz23))));
      cz = ((dx13 * dy23) + (-((dx23 * dy13))));
      v = ((cx * dx43) + (cy * dy43) + (cz * dz43));
      cse_p11_t1_g51329 = (cse_p1_t1_sqrt51302 * cse_p1_t3_sqrt51304);
      denom = (cse_p11_t1_g51329 * cse_p1_t5_sqrt51306);
      cse_p1_t7_invr51308 = (1.0 / (denom));
      cse_p1_t8_invr251309 = (cse_p1_t7_invr51308 * cse_p1_t7_invr51308);
      cse_p113_t1_g51330 = (cse_p1_t8_invr251309 * v);
      q = (cse_p1_t7_invr51308 * v);
      cse_p11_t1_g51310 = (co + q);
      cse_p11_t11_g51320 = ((cse_p11_t1_g51310) * (cse_p11_t1_g51310));
      energy = (k * (cse_p11_t1_g51310 * cse_p11_t1_g51310 * cse_p11_t1_g51310));
      energy_added += energy;
      *energy_accumulate += energy;
      cse_p11_t1_g51321 = (3.0000000000000000     * cse_p11_t11_g51320 * k);
      cse_p11_t1_g51322 = (cse_p1_t2_invsqrt51303 * cse_p1_t3_sqrt51304 * cse_p1_t5_sqrt51306);
      cse_p11_t1_g51326 = (cse_p113_t1_g51330 * cse_p11_t1_g51322);
      g_x1 = (cse_p11_t1_g51321 * ((cse_p1_t7_invr51308 * ((cse_p11_t7_g51316 * dy43) + (dy23 * dz43))) + (-((cse_p11_t1_g51326 * dx13)))));
      KernelGradientAcc(i3x1, 0, g_x1);
      g_y1 = (cse_p11_t1_g51321 * ((cse_p1_t7_invr51308 * ((cse_p11_t3_g51312 * dz43) + (dx43 * dz23))) + (-((cse_p11_t1_g51326 * dy13)))));
      KernelGradientAcc(i3x1, 1, g_y1);
      g_z1 = (cse_p11_t1_g51321 * ((cse_p1_t7_invr51308 * ((cse_p11_t5_g51314 * dx43) + (dx23 * dy43))) + (-((cse_p11_t1_g51326 * dz13)))));
      KernelGradientAcc(i3x1, 2, g_z1);
      cse_p11_t1_g51324 = (cse_p1_t1_sqrt51302 * cse_p1_t4_invsqrt51305 * cse_p1_t5_sqrt51306);
      cse_p11_t1_g51327 = (cse_p113_t1_g51330 * cse_p11_t1_g51324);
      g_x2 = (cse_p11_t1_g51321 * ((cse_p1_t7_invr51308 * ((cse_p11_t4_g51313 * dz43) + (dy43 * dz13))) + (-((cse_p11_t1_g51327 * dx23)))));
      KernelGradientAcc(i3x2, 0, g_x2);
      g_y2 = (cse_p11_t1_g51321 * ((cse_p1_t7_invr51308 * ((cse_p11_t6_g51315 * dx43) + (dx13 * dz43))) + (-((cse_p11_t1_g51327 * dy23)))));
      KernelGradientAcc(i3x2, 1, g_y2);
      g_z2 = (cse_p11_t1_g51321 * ((cse_p1_t7_invr51308 * ((cse_p11_t2_g51311 * dy43) + (dx43 * dy13))) + (-((cse_p11_t1_g51327 * dz23)))));
      KernelGradientAcc(i3x2, 2, g_z2);
      cse_p11_t1_g51325 = (cse_p11_t1_g51329 * cse_p1_t6_invsqrt51307);
      g_x3 = (cse_p11_t1_g51321 * ((cse_p1_t7_invr51308 * ((dy43 * (cse_p11_t6_g51315 + dz23)) + (dz43 * (cse_p11_t5_g51314 + dy13)) + (-(cx)))) + (-((cse_p113_t1_g51330 * ((-((cse_p11_t1_g51322 * dx13))) + (-((cse_p11_t1_g51324 * dx23))) + (-((cse_p11_t1_g51325 * dx43)))))))));
      KernelGradientAcc(i3x3, 0, g_x3);
      g_y3 = (cse_p11_t1_g51321 * ((cse_p1_t7_invr51308 * ((dx43 * (cse_p11_t7_g51316 + dz13)) + (dz43 * (cse_p11_t2_g51311 + dx23)) + (-(cy)))) + (-((cse_p113_t1_g51330 * ((-((cse_p11_t1_g51322 * dy13))) + (-((cse_p11_t1_g51324 * dy23))) + (-((cse_p11_t1_g51325 * dy43)))))))));
      KernelGradientAcc(i3x3, 1, g_y3);
      g_z3 = (cse_p11_t1_g51321 * ((cse_p1_t7_invr51308 * ((dx43 * (cse_p11_t4_g51313 + dy23)) + (dy43 * (cse_p11_t3_g51312 + dx13)) + (-(cz)))) + (-((cse_p113_t1_g51330 * ((-((cse_p11_t1_g51322 * dz13))) + (-((cse_p11_t1_g51324 * dz23))) + (-((cse_p11_t1_g51325 * dz43)))))))));
      KernelGradientAcc(i3x3, 2, g_z3);
      cse_p11_t1_g51328 = (cse_p113_t1_g51330 * cse_p11_t1_g51325);
      g_x4 = (cse_p11_t1_g51321 * ((cse_p1_t7_invr51308 * cx) + (-((cse_p11_t1_g51328 * dx43)))));
      KernelGradientAcc(i3x4, 0, g_x4);
      g_y4 = (cse_p11_t1_g51321 * ((cse_p1_t7_invr51308 * cy) + (-((cse_p11_t1_g51328 * dy43)))));
      KernelGradientAcc(i3x4, 1, g_y4);
      g_z4 = (cse_p11_t1_g51321 * ((cse_p1_t7_invr51308 * cz) + (-((cse_p11_t1_g51328 * dz43)))));
      KernelGradientAcc(i3x4, 2, g_z4);
    }
  }
  return energy_added;
}
void gradient_fd(double k, double co, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 12;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  double e0 = 0.0;
  energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
  {
    double saved = position[i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 0] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 0, d);
  }
  {
    double saved = position[i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 1] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 1, d);
  }
  {
    double saved = position[i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 2] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 2, d);
  }
  {
    double saved = position[i3x4 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 0] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 0, d);
  }
  {
    double saved = position[i3x4 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 1] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 1, d);
  }
  {
    double saved = position[i3x4 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 2] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 2, d);
  }
}

double hessian(double k, double co, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double cse_p113_t1_g51359;
  double cse_p11_t10_g51348;
  double cse_p11_t11_g51349;
  double cse_p11_t1_g51339;
  double cse_p11_t1_g51350;
  double cse_p11_t1_g51351;
  double cse_p11_t1_g51353;
  double cse_p11_t1_g51354;
  double cse_p11_t1_g51355;
  double cse_p11_t1_g51356;
  double cse_p11_t1_g51357;
  double cse_p11_t1_g51358;
  double cse_p11_t2_g51340;
  double cse_p11_t3_g51341;
  double cse_p11_t4_g51342;
  double cse_p11_t5_g51343;
  double cse_p11_t6_g51344;
  double cse_p11_t7_g51345;
  double cse_p11_t8_g51346;
  double cse_p11_t9_g51347;
  double cse_p1_t1_sqrt51331;
  double cse_p1_t2_invsqrt51332;
  double cse_p1_t3_sqrt51333;
  double cse_p1_t4_invsqrt51334;
  double cse_p1_t5_sqrt51335;
  double cse_p1_t6_invsqrt51336;
  double cse_p1_t7_invr51337;
  double cse_p1_t8_invr251338;
  double cx;
  double cy;
  double cz;
  double denom;
  double dx13;
  double dx23;
  double dx43;
  double dy13;
  double dy23;
  double dy43;
  double dz13;
  double dz23;
  double dz43;
  double energy;
  double energy_added;
  double g_x1;
  double g_x2;
  double g_x3;
  double g_x4;
  double g_y1;
  double g_y2;
  double g_y3;
  double g_y4;
  double g_z1;
  double g_z2;
  double g_z3;
  double g_z4;
  double q;
  double r13_2;
  double r23_2;
  double r43_2;
  double v;
  energy_added = 0.0000000000000000    ;
  {
    /* !BASE */
    DOUBLE x1 = position[i3x1 + 0];
    DOUBLE y1 = position[i3x1 + 1];
    DOUBLE z1 = position[i3x1 + 2];
    DOUBLE x2 = position[i3x2 + 0];
    DOUBLE y2 = position[i3x2 + 1];
    DOUBLE z2 = position[i3x2 + 2];
    DOUBLE x3 = position[i3x3 + 0];
    DOUBLE y3 = position[i3x3 + 1];
    DOUBLE z3 = position[i3x3 + 2];
    DOUBLE x4 = position[i3x4 + 0];
    DOUBLE y4 = position[i3x4 + 1];
    DOUBLE z4 = position[i3x4 + 2];
    {
      /* !BASE */
      cse_p11_t8_g51346 = (-(x3));
      cse_p11_t9_g51347 = (-(y3));
      cse_p11_t10_g51348 = (-(z3));
      dx13 = (cse_p11_t8_g51346 + x1);
      cse_p11_t2_g51340 = (-(dx13));
      dy13 = (cse_p11_t9_g51347 + y1);
      cse_p11_t4_g51342 = (-(dy13));
      dz13 = (cse_p11_t10_g51348 + z1);
      cse_p11_t6_g51344 = (-(dz13));
      dx23 = (cse_p11_t8_g51346 + x2);
      cse_p11_t3_g51341 = (-(dx23));
      dy23 = (cse_p11_t9_g51347 + y2);
      cse_p11_t5_g51343 = (-(dy23));
      dz23 = (cse_p11_t10_g51348 + z2);
      cse_p11_t7_g51345 = (-(dz23));
      dx43 = (cse_p11_t8_g51346 + x4);
      dy43 = (cse_p11_t9_g51347 + y4);
      dz43 = (cse_p11_t10_g51348 + z4);
      r13_2 = ((dx13 * dx13) + (dy13 * dy13) + (dz13 * dz13));
      r23_2 = ((dx23 * dx23) + (dy23 * dy23) + (dz23 * dz23));
      r43_2 = ((dx43 * dx43) + (dy43 * dy43) + (dz43 * dz43));
      cse_p1_t1_sqrt51331 = sqrt(r13_2);
      cse_p1_t2_invsqrt51332 = (1.0 / (cse_p1_t1_sqrt51331));
      cse_p1_t3_sqrt51333 = sqrt(r23_2);
      cse_p1_t4_invsqrt51334 = (1.0 / (cse_p1_t3_sqrt51333));
      cse_p1_t5_sqrt51335 = sqrt(r43_2);
      cse_p1_t6_invsqrt51336 = (1.0 / (cse_p1_t5_sqrt51335));
      cx = ((dy13 * dz23) + (-((dy23 * dz13))));
      cy = ((dx23 * dz13) + (-((dx13 * dz23))));
      cz = ((dx13 * dy23) + (-((dx23 * dy13))));
      v = ((cx * dx43) + (cy * dy43) + (cz * dz43));
      cse_p11_t1_g51358 = (cse_p1_t1_sqrt51331 * cse_p1_t3_sqrt51333);
      denom = (cse_p11_t1_g51358 * cse_p1_t5_sqrt51335);
      cse_p1_t7_invr51337 = (1.0 / (denom));
      cse_p1_t8_invr251338 = (cse_p1_t7_invr51337 * cse_p1_t7_invr51337);
      cse_p113_t1_g51359 = (cse_p1_t8_invr251338 * v);
      q = (cse_p1_t7_invr51337 * v);
      cse_p11_t1_g51339 = (co + q);
      cse_p11_t11_g51349 = ((cse_p11_t1_g51339) * (cse_p11_t1_g51339));
      energy = (k * (cse_p11_t1_g51339 * cse_p11_t1_g51339 * cse_p11_t1_g51339));
      energy_added += energy;
      *energy_accumulate += energy;
      cse_p11_t1_g51350 = (3.0000000000000000     * cse_p11_t11_g51349 * k);
      cse_p11_t1_g51351 = (cse_p1_t2_invsqrt51332 * cse_p1_t3_sqrt51333 * cse_p1_t5_sqrt51335);
      cse_p11_t1_g51355 = (cse_p113_t1_g51359 * cse_p11_t1_g51351);
      g_x1 = (cse_p11_t1_g51350 * ((cse_p1_t7_invr51337 * ((cse_p11_t7_g51345 * dy43) + (dy23 * dz43))) + (-((cse_p11_t1_g51355 * dx13)))));
      KernelGradientAcc(i3x1, 0, g_x1);
      g_y1 = (cse_p11_t1_g51350 * ((cse_p1_t7_invr51337 * ((cse_p11_t3_g51341 * dz43) + (dx43 * dz23))) + (-((cse_p11_t1_g51355 * dy13)))));
      KernelGradientAcc(i3x1, 1, g_y1);
      g_z1 = (cse_p11_t1_g51350 * ((cse_p1_t7_invr51337 * ((cse_p11_t5_g51343 * dx43) + (dx23 * dy43))) + (-((cse_p11_t1_g51355 * dz13)))));
      KernelGradientAcc(i3x1, 2, g_z1);
      cse_p11_t1_g51353 = (cse_p1_t1_sqrt51331 * cse_p1_t4_invsqrt51334 * cse_p1_t5_sqrt51335);
      cse_p11_t1_g51356 = (cse_p113_t1_g51359 * cse_p11_t1_g51353);
      g_x2 = (cse_p11_t1_g51350 * ((cse_p1_t7_invr51337 * ((cse_p11_t4_g51342 * dz43) + (dy43 * dz13))) + (-((cse_p11_t1_g51356 * dx23)))));
      KernelGradientAcc(i3x2, 0, g_x2);
      g_y2 = (cse_p11_t1_g51350 * ((cse_p1_t7_invr51337 * ((cse_p11_t6_g51344 * dx43) + (dx13 * dz43))) + (-((cse_p11_t1_g51356 * dy23)))));
      KernelGradientAcc(i3x2, 1, g_y2);
      g_z2 = (cse_p11_t1_g51350 * ((cse_p1_t7_invr51337 * ((cse_p11_t2_g51340 * dy43) + (dx43 * dy13))) + (-((cse_p11_t1_g51356 * dz23)))));
      KernelGradientAcc(i3x2, 2, g_z2);
      cse_p11_t1_g51354 = (cse_p11_t1_g51358 * cse_p1_t6_invsqrt51336);
      g_x3 = (cse_p11_t1_g51350 * ((cse_p1_t7_invr51337 * ((dy43 * (cse_p11_t6_g51344 + dz23)) + (dz43 * (cse_p11_t5_g51343 + dy13)) + (-(cx)))) + (-((cse_p113_t1_g51359 * ((-((cse_p11_t1_g51351 * dx13))) + (-((cse_p11_t1_g51353 * dx23))) + (-((cse_p11_t1_g51354 * dx43)))))))));
      KernelGradientAcc(i3x3, 0, g_x3);
      g_y3 = (cse_p11_t1_g51350 * ((cse_p1_t7_invr51337 * ((dx43 * (cse_p11_t7_g51345 + dz13)) + (dz43 * (cse_p11_t2_g51340 + dx23)) + (-(cy)))) + (-((cse_p113_t1_g51359 * ((-((cse_p11_t1_g51351 * dy13))) + (-((cse_p11_t1_g51353 * dy23))) + (-((cse_p11_t1_g51354 * dy43)))))))));
      KernelGradientAcc(i3x3, 1, g_y3);
      g_z3 = (cse_p11_t1_g51350 * ((cse_p1_t7_invr51337 * ((dx43 * (cse_p11_t4_g51342 + dy23)) + (dy43 * (cse_p11_t3_g51341 + dx13)) + (-(cz)))) + (-((cse_p113_t1_g51359 * ((-((cse_p11_t1_g51351 * dz13))) + (-((cse_p11_t1_g51353 * dz23))) + (-((cse_p11_t1_g51354 * dz43)))))))));
      KernelGradientAcc(i3x3, 2, g_z3);
      cse_p11_t1_g51357 = (cse_p113_t1_g51359 * cse_p11_t1_g51354);
      g_x4 = (cse_p11_t1_g51350 * ((cse_p1_t7_invr51337 * cx) + (-((cse_p11_t1_g51357 * dx43)))));
      KernelGradientAcc(i3x4, 0, g_x4);
      g_y4 = (cse_p11_t1_g51350 * ((cse_p1_t7_invr51337 * cy) + (-((cse_p11_t1_g51357 * dy43)))));
      KernelGradientAcc(i3x4, 1, g_y4);
      g_z4 = (cse_p11_t1_g51350 * ((cse_p1_t7_invr51337 * cz) + (-((cse_p11_t1_g51357 * dz43)))));
      KernelGradientAcc(i3x4, 2, g_z4);
    }
  }
  return energy_added;
}
void hessian_fd(double k, double co, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 12;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  const double invh2 = 1.0/((h*h));
  double e0 = 0.0;
  energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
  {
    double saved = position[i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 0] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 0, d);
  }
  {
    double saved = position[i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 1] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 1, d);
  }
  {
    double saved = position[i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 2] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 2, d);
  }
  {
    double saved = position[i3x4 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 0] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 0, d);
  }
  {
    double saved = position[i3x4 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 1] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 1, d);
  }
  {
    double saved = position[i3x4 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 2] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 2, d);
  }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 0, i3x1, 0, hval);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 1, i3x1, 1, hval);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 2, i3x1, 2, hval);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 0, i3x2, 0, hval);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 1, i3x2, 1, hval);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x2, 2, hval);
  }
  {
    double saved = position[i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 0] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 0, i3x3, 0, hval);
  }
  {
    double saved = position[i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 1] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x3, 1, hval);
  }
  {
    double saved = position[i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 2] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x3, 2, hval);
  }
  {
    double saved = position[i3x4 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 0] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 0, i3x4, 0, hval);
  }
  {
    double saved = position[i3x4 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 1] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 1, i3x4, 1, hval);
  }
  {
    double saved = position[i3x4 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 2] = saved + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 2, i3x4, 2, hval);
  }
  {
    double saved_i = position[i3x1 + 1];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x1 + 1] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i; position[i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x2, 1, hval);
  }
  {
    double saved_i = position[i3x3 + 0];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 0] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 0, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x3 + 0];
    double saved_j = position[i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 0] = saved_i + h; position[i3x1 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i; position[i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 0, i3x1, 1, hval);
  }
  {
    double saved_i = position[i3x3 + 0];
    double saved_j = position[i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 0] = saved_i + h; position[i3x1 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i; position[i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 0, i3x1, 2, hval);
  }
  {
    double saved_i = position[i3x3 + 0];
    double saved_j = position[i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 0] = saved_i + h; position[i3x2 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i; position[i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 0, i3x2, 0, hval);
  }
  {
    double saved_i = position[i3x3 + 0];
    double saved_j = position[i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 0] = saved_i + h; position[i3x2 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i; position[i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 0, i3x2, 1, hval);
  }
  {
    double saved_i = position[i3x3 + 0];
    double saved_j = position[i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 0] = saved_i + h; position[i3x2 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i; position[i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 0, i3x2, 2, hval);
  }
  {
    double saved_i = position[i3x3 + 1];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 1] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x3 + 1];
    double saved_j = position[i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 1] = saved_i + h; position[i3x1 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i; position[i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x1, 1, hval);
  }
  {
    double saved_i = position[i3x3 + 1];
    double saved_j = position[i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 1] = saved_i + h; position[i3x1 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i; position[i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x1, 2, hval);
  }
  {
    double saved_i = position[i3x3 + 1];
    double saved_j = position[i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 1] = saved_i + h; position[i3x2 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i; position[i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x2, 0, hval);
  }
  {
    double saved_i = position[i3x3 + 1];
    double saved_j = position[i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 1] = saved_i + h; position[i3x2 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i; position[i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x2, 1, hval);
  }
  {
    double saved_i = position[i3x3 + 1];
    double saved_j = position[i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 1] = saved_i + h; position[i3x2 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i; position[i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x2, 2, hval);
  }
  {
    double saved_i = position[i3x3 + 1];
    double saved_j = position[i3x3 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 1] = saved_i + h; position[i3x3 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i; position[i3x3 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x3, 0, hval);
  }
  {
    double saved_i = position[i3x3 + 2];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 2] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x3 + 2];
    double saved_j = position[i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 2] = saved_i + h; position[i3x1 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i; position[i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x1, 1, hval);
  }
  {
    double saved_i = position[i3x3 + 2];
    double saved_j = position[i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 2] = saved_i + h; position[i3x1 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i; position[i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x1, 2, hval);
  }
  {
    double saved_i = position[i3x3 + 2];
    double saved_j = position[i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 2] = saved_i + h; position[i3x2 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i; position[i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x2, 0, hval);
  }
  {
    double saved_i = position[i3x3 + 2];
    double saved_j = position[i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 2] = saved_i + h; position[i3x2 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i; position[i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x2, 1, hval);
  }
  {
    double saved_i = position[i3x3 + 2];
    double saved_j = position[i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 2] = saved_i + h; position[i3x2 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i; position[i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x2, 2, hval);
  }
  {
    double saved_i = position[i3x3 + 2];
    double saved_j = position[i3x3 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 2] = saved_i + h; position[i3x3 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i; position[i3x3 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x3, 0, hval);
  }
  {
    double saved_i = position[i3x3 + 2];
    double saved_j = position[i3x3 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x3 + 2] = saved_i + h; position[i3x3 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x3 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i; position[i3x3 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x3, 1, hval);
  }
  {
    double saved_i = position[i3x4 + 0];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 0] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 0, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x4 + 0];
    double saved_j = position[i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 0] = saved_i + h; position[i3x1 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i; position[i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 0, i3x1, 1, hval);
  }
  {
    double saved_i = position[i3x4 + 0];
    double saved_j = position[i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 0] = saved_i + h; position[i3x1 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i; position[i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 0, i3x1, 2, hval);
  }
  {
    double saved_i = position[i3x4 + 0];
    double saved_j = position[i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 0] = saved_i + h; position[i3x2 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i; position[i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 0, i3x2, 0, hval);
  }
  {
    double saved_i = position[i3x4 + 0];
    double saved_j = position[i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 0] = saved_i + h; position[i3x2 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i; position[i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 0, i3x2, 1, hval);
  }
  {
    double saved_i = position[i3x4 + 0];
    double saved_j = position[i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 0] = saved_i + h; position[i3x2 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i; position[i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 0, i3x2, 2, hval);
  }
  {
    double saved_i = position[i3x4 + 0];
    double saved_j = position[i3x3 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 0] = saved_i + h; position[i3x3 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i; position[i3x3 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 0, i3x3, 0, hval);
  }
  {
    double saved_i = position[i3x4 + 0];
    double saved_j = position[i3x3 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 0] = saved_i + h; position[i3x3 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x3 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i; position[i3x3 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 0, i3x3, 1, hval);
  }
  {
    double saved_i = position[i3x4 + 0];
    double saved_j = position[i3x3 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 0] = saved_i + h; position[i3x3 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x3 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i; position[i3x3 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 0, i3x3, 2, hval);
  }
  {
    double saved_i = position[i3x4 + 1];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 1] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 1, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x4 + 1];
    double saved_j = position[i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 1] = saved_i + h; position[i3x1 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i; position[i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 1, i3x1, 1, hval);
  }
  {
    double saved_i = position[i3x4 + 1];
    double saved_j = position[i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 1] = saved_i + h; position[i3x1 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i; position[i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 1, i3x1, 2, hval);
  }
  {
    double saved_i = position[i3x4 + 1];
    double saved_j = position[i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 1] = saved_i + h; position[i3x2 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i; position[i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 1, i3x2, 0, hval);
  }
  {
    double saved_i = position[i3x4 + 1];
    double saved_j = position[i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 1] = saved_i + h; position[i3x2 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i; position[i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 1, i3x2, 1, hval);
  }
  {
    double saved_i = position[i3x4 + 1];
    double saved_j = position[i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 1] = saved_i + h; position[i3x2 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i; position[i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 1, i3x2, 2, hval);
  }
  {
    double saved_i = position[i3x4 + 1];
    double saved_j = position[i3x3 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 1] = saved_i + h; position[i3x3 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i; position[i3x3 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 1, i3x3, 0, hval);
  }
  {
    double saved_i = position[i3x4 + 1];
    double saved_j = position[i3x3 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 1] = saved_i + h; position[i3x3 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x3 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i; position[i3x3 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 1, i3x3, 1, hval);
  }
  {
    double saved_i = position[i3x4 + 1];
    double saved_j = position[i3x3 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 1] = saved_i + h; position[i3x3 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x3 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i; position[i3x3 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 1, i3x3, 2, hval);
  }
  {
    double saved_i = position[i3x4 + 1];
    double saved_j = position[i3x4 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 1] = saved_i + h; position[i3x4 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x4 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i; position[i3x4 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 1, i3x4, 0, hval);
  }
  {
    double saved_i = position[i3x4 + 2];
    double saved_j = position[i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 2] = saved_i + h; position[i3x1 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i; position[i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 2, i3x1, 0, hval);
  }
  {
    double saved_i = position[i3x4 + 2];
    double saved_j = position[i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 2] = saved_i + h; position[i3x1 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i; position[i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 2, i3x1, 1, hval);
  }
  {
    double saved_i = position[i3x4 + 2];
    double saved_j = position[i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 2] = saved_i + h; position[i3x1 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i; position[i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 2, i3x1, 2, hval);
  }
  {
    double saved_i = position[i3x4 + 2];
    double saved_j = position[i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 2] = saved_i + h; position[i3x2 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i; position[i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 2, i3x2, 0, hval);
  }
  {
    double saved_i = position[i3x4 + 2];
    double saved_j = position[i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 2] = saved_i + h; position[i3x2 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i; position[i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 2, i3x2, 1, hval);
  }
  {
    double saved_i = position[i3x4 + 2];
    double saved_j = position[i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 2] = saved_i + h; position[i3x2 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i; position[i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 2, i3x2, 2, hval);
  }
  {
    double saved_i = position[i3x4 + 2];
    double saved_j = position[i3x3 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 2] = saved_i + h; position[i3x3 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i; position[i3x3 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 2, i3x3, 0, hval);
  }
  {
    double saved_i = position[i3x4 + 2];
    double saved_j = position[i3x3 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 2] = saved_i + h; position[i3x3 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x3 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i; position[i3x3 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 2, i3x3, 1, hval);
  }
  {
    double saved_i = position[i3x4 + 2];
    double saved_j = position[i3x3 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 2] = saved_i + h; position[i3x3 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x3 + 2] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i; position[i3x3 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 2, i3x3, 2, hval);
  }
  {
    double saved_i = position[i3x4 + 2];
    double saved_j = position[i3x4 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 2] = saved_i + h; position[i3x4 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x4 + 0] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i; position[i3x4 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 2, i3x4, 0, hval);
  }
  {
    double saved_i = position[i3x4 + 2];
    double saved_j = position[i3x4 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[i3x4 + 2] = saved_i + h; position[i3x4 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x4 + 1] = saved_j + h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_j - h;
    energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i; position[i3x4 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 2, i3x4, 1, hval);
  }
}

};
