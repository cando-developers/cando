template <typename HESSIAN>
struct Stretch {
  static constexpr size_t PositionSize = 6;
  static std::string description() { return "mathkernel-stretch"; };
void hessian(double kb, double r0, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  /* !BASE */
  double cse_p114_t1_g340;
  double cse_p11_t10_g246;
  double cse_p11_t11_g247;
  double cse_p11_t12_g248;
  double cse_p11_t13_g249;
  double cse_p11_t14_g250;
  double cse_p11_t15_g251;
  double cse_p11_t16_g252;
  double cse_p11_t17_g253;
  double cse_p11_t18_g254;
  double cse_p11_t19_g255;
  double cse_p11_t1_g237;
  double cse_p11_t1_g296;
  double cse_p11_t1_g301;
  double cse_p11_t1_g302;
  double cse_p11_t1_g303;
  double cse_p11_t1_g304;
  double cse_p11_t1_g305;
  double cse_p11_t1_g306;
  double cse_p11_t1_g307;
  double cse_p11_t1_g308;
  double cse_p11_t1_g309;
  double cse_p11_t20_g256;
  double cse_p11_t21_g257;
  double cse_p11_t23_g259;
  double cse_p11_t24_g260;
  double cse_p11_t25_g261;
  double cse_p11_t26_g262;
  double cse_p11_t27_g263;
  double cse_p11_t28_g264;
  double cse_p11_t29_g265;
  double cse_p11_t2_g238;
  double cse_p11_t30_g266;
  double cse_p11_t31_g267;
  double cse_p11_t34_g270;
  double cse_p11_t38_g274;
  double cse_p11_t39_g275;
  double cse_p11_t3_g239;
  double cse_p11_t40_g276;
  double cse_p11_t41_g277;
  double cse_p11_t42_g278;
  double cse_p11_t43_g279;
  double cse_p11_t44_g280;
  double cse_p11_t45_g281;
  double cse_p11_t46_g282;
  double cse_p11_t47_g283;
  double cse_p11_t48_g284;
  double cse_p11_t49_g285;
  double cse_p11_t4_g240;
  double cse_p11_t50_g286;
  double cse_p11_t51_g287;
  double cse_p11_t5_g241;
  double cse_p11_t6_g242;
  double cse_p11_t7_g243;
  double cse_p11_t8_g244;
  double cse_p11_t9_g245;
  double cse_p12_t1_g314;
  double cse_p12_t1_g315;
  double cse_p12_t1_g316;
  double cse_p12_t1_g317;
  double cse_p12_t1_g318;
  double cse_p12_t1_g319;
  double cse_p12_t1_g320;
  double cse_p12_t1_g321;
  double cse_p12_t1_g322;
  double cse_p12_t1_g323;
  double cse_p13_t1_g324;
  double cse_p13_t1_g325;
  double cse_p13_t1_g326;
  double cse_p166_t1_g345;
  double cse_p1_t2_invsqrt226;
  double cse_p529_t1_sqrt378;
  double cse_p529_t2_invr379;
  double cse_p61_t3_invr2329;
  double cse_p62_t1_g330;
  double cse_p63_t1_g334;
  double cse_p63_t1_g335;
  double cse_p63_t1_g336;
  double cse_p63_t1_g337;
  double cse_p63_t2_g332;
  double cse_p63_t3_g333;
  double dr;
  double dx;
  double dy;
  double dz;
  double energy;
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
  double r2;
  DOUBLE x1 = position[i3x1 + 0];
  DOUBLE y1 = position[i3x1 + 1];
  DOUBLE z1 = position[i3x1 + 2];
  DOUBLE x2 = position[i3x2 + 0];
  DOUBLE y2 = position[i3x2 + 1];
  DOUBLE z2 = position[i3x2 + 2];
  {
    /* !BASE */
    dx = (x2 + (-(x1)));
    dy = (y2 + (-(y1)));
    dz = (z2 + (-(z1)));
    cse_p11_t1_g306 = (dy * dy);
    cse_p11_t1_g307 = (dz * dz);
    cse_p11_t1_g308 = (dx * dx);
    r2 = (cse_p11_t1_g306 + cse_p11_t1_g307 + cse_p11_t1_g308);
    cse_p529_t1_sqrt378 = sqrt(r2);
    cse_p529_t2_invr379 = (1.0 / (r2));
    cse_p61_t3_invr2329 = (cse_p529_t2_invr379 * cse_p529_t2_invr379);
    cse_p11_t1_g296 = (cse_p529_t1_sqrt378 * cse_p529_t2_invr379);
    cse_p63_t1_g336 = (cse_p11_t1_g296 * kb);
    cse_p11_t1_g302 = (cse_p529_t2_invr379 * cse_p63_t1_g336);
    cse_p11_t1_g305 = (cse_p11_t1_g302 * dx);
    cse_p12_t1_g316 = (cse_p11_t1_g305 * cse_p529_t1_sqrt378);
    cse_p11_t7_g243 = (-2.0000000000000000     * cse_p12_t1_g316 * dy);
    cse_p11_t8_g244 = (-2.0000000000000000     * cse_p12_t1_g316 * dz);
    cse_p11_t1_g309 = (dy * dz);
    cse_p12_t1_g317 = (cse_p11_t1_g302 * cse_p11_t1_g309);
    cse_p13_t1_g326 = (cse_p12_t1_g317 * cse_p529_t1_sqrt378);
    cse_p11_t9_g245 = (2.0000000000000000     * (-(cse_p13_t1_g326)));
    cse_p11_t26_g262 = (cse_p11_t1_g296 * cse_p11_t1_g308);
    cse_p11_t28_g264 = (cse_p12_t1_g316 * dy);
    cse_p11_t30_g266 = (cse_p12_t1_g316 * dz);
    cse_p11_t31_g267 = (cse_p11_t1_g296 * cse_p11_t1_g306);
    cse_p11_t34_g270 = (cse_p11_t1_g296 * cse_p11_t1_g307);
    cse_p11_t48_g284 = (-(cse_p11_t28_g264));
    cse_p114_t1_g340 = (2.0000000000000000     * cse_p11_t48_g284);
    cse_p11_t49_g285 = (-(cse_p11_t30_g266));
    cse_p63_t2_g332 = (2.0000000000000000     * cse_p11_t49_g285);
    cse_p11_t50_g286 = (-(cse_p13_t1_g326));
    cse_p63_t3_g333 = (2.0000000000000000     * cse_p11_t50_g286);
    cse_p1_t2_invsqrt226 = (1.0 / (cse_p529_t1_sqrt378));
    dr = (cse_p529_t1_sqrt378 + (-(r0)));
    cse_p62_t1_g330 = (cse_p12_t1_g317 * dr);
    cse_p13_t1_g324 = (dr * kb);
    cse_p11_t1_g301 = (cse_p13_t1_g324 * cse_p529_t1_sqrt378 * cse_p61_t3_invr2329);
    cse_p11_t1_g304 = (-2.0000000000000000     * cse_p11_t1_g301);
    cse_p11_t1_g237 = (cse_p11_t1_g304 * cse_p11_t1_g308);
    cse_p12_t1_g322 = (cse_p11_t1_g304 * dx);
    cse_p11_t2_g238 = (cse_p12_t1_g322 * dy);
    cse_p11_t3_g239 = (cse_p12_t1_g322 * dz);
    cse_p11_t4_g240 = (cse_p11_t1_g304 * cse_p11_t1_g306);
    cse_p11_t5_g241 = (cse_p11_t1_g304 * cse_p11_t1_g309);
    cse_p11_t6_g242 = (cse_p11_t1_g304 * cse_p11_t1_g307);
    cse_p12_t1_g315 = (2.0000000000000000     * cse_p11_t1_g301);
    cse_p11_t10_g246 = (cse_p11_t1_g308 * cse_p12_t1_g315);
    cse_p12_t1_g323 = (cse_p12_t1_g315 * dx);
    cse_p11_t11_g247 = (cse_p12_t1_g323 * dy);
    cse_p11_t12_g248 = (cse_p12_t1_g323 * dz);
    cse_p11_t13_g249 = (cse_p11_t1_g306 * cse_p12_t1_g315);
    cse_p11_t14_g250 = (cse_p11_t1_g309 * cse_p12_t1_g315);
    cse_p11_t15_g251 = (cse_p11_t1_g307 * cse_p12_t1_g315);
    cse_p11_t1_g303 = (cse_p529_t2_invr379 * dr);
    cse_p11_t16_g252 = (cse_p11_t1_g303 * cse_p11_t1_g308);
    cse_p12_t1_g321 = (cse_p11_t1_g305 * dr);
    cse_p11_t17_g253 = (cse_p12_t1_g321 * dx);
    cse_p11_t18_g254 = (cse_p12_t1_g321 * dy);
    cse_p11_t19_g255 = (cse_p12_t1_g321 * dz);
    cse_p11_t20_g256 = (cse_p11_t1_g303 * cse_p11_t1_g306);
    cse_p11_t21_g257 = (cse_p11_t1_g302 * cse_p11_t1_g306 * dr);
    cse_p11_t23_g259 = (cse_p11_t1_g303 * cse_p11_t1_g307);
    cse_p11_t24_g260 = (cse_p11_t1_g302 * cse_p11_t1_g307 * dr);
    cse_p11_t25_g261 = (cse_p11_t1_g296 * cse_p13_t1_g324);
    cse_p12_t1_g314 = (cse_p11_t1_g305 * (cse_p529_t1_sqrt378 + cse_p529_t1_sqrt378 + dr + dr));
    cse_p11_t27_g263 = (cse_p12_t1_g314 * dy);
    cse_p11_t29_g265 = (cse_p12_t1_g314 * dz);
    cse_p13_t1_g325 = (cse_p12_t1_g317 * (cse_p529_t1_sqrt378 + cse_p529_t1_sqrt378 + dr + dr));
    cse_p12_t1_g320 = (cse_p63_t1_g336 * (cse_p11_t16_g252 + cse_p11_t16_g252 + cse_p11_t26_g262 + cse_p11_t26_g262 + dr + dr));
    cse_p12_t1_g318 = (cse_p63_t1_g336 * (cse_p11_t20_g256 + cse_p11_t20_g256 + cse_p11_t31_g267 + cse_p11_t31_g267 + dr + dr));
    cse_p12_t1_g319 = (cse_p63_t1_g336 * (cse_p11_t23_g259 + cse_p11_t23_g259 + cse_p11_t34_g270 + cse_p11_t34_g270 + dr + dr));
    cse_p11_t38_g274 = (cse_p11_t16_g252 + cse_p11_t16_g252 + cse_p11_t26_g262 + cse_p11_t26_g262 + dr + dr);
    cse_p11_t39_g275 = (cse_p11_t20_g256 + cse_p11_t20_g256 + cse_p11_t31_g267 + cse_p11_t31_g267 + dr + dr);
    cse_p11_t40_g276 = (cse_p11_t23_g259 + cse_p11_t23_g259 + cse_p11_t34_g270 + cse_p11_t34_g270 + dr + dr);
    cse_p11_t41_g277 = (cse_p529_t1_sqrt378 + cse_p529_t1_sqrt378 + dr + dr);
    cse_p11_t42_g278 = (-(cse_p11_t17_g253));
    cse_p11_t43_g279 = (-(cse_p11_t18_g254));
    cse_p11_t44_g280 = (-(cse_p11_t19_g255));
    cse_p11_t45_g281 = (-(cse_p11_t21_g257));
    cse_p11_t46_g282 = (-(cse_p62_t1_g330));
    cse_p11_t47_g283 = (-(cse_p11_t24_g260));
    cse_p11_t51_g287 = (-(cse_p11_t25_g261));
    energy = (cse_p13_t1_g324 * dr);
    *energy_accumulate += energy;
    cse_p63_t1_g334 = (cse_p13_t1_g324 * cse_p1_t2_invsqrt226);
    cse_p166_t1_g345 = (2.0000000000000000     * cse_p63_t1_g334);
    cse_p63_t1_g337 = (-2.0000000000000000     * cse_p63_t1_g334);
    g_x1 = (cse_p63_t1_g337 * dx);
    KernelGradientAcc(i3x1, 0, g_x1);
    g_y1 = (cse_p63_t1_g337 * dy);
    KernelGradientAcc(i3x1, 1, g_y1);
    g_z1 = (cse_p63_t1_g337 * dz);
    KernelGradientAcc(i3x1, 2, g_z1);
    g_x2 = (cse_p166_t1_g345 * dx);
    KernelGradientAcc(i3x2, 0, g_x2);
    g_y2 = (cse_p166_t1_g345 * dy);
    KernelGradientAcc(i3x2, 1, g_y2);
    g_z2 = (cse_p166_t1_g345 * dz);
    KernelGradientAcc(i3x2, 2, g_z2);
    h_x1_x1 = (cse_p11_t1_g237 + cse_p11_t1_g237 + cse_p12_t1_g320);
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 0, h_x1_x1);
    h_x1_y1 = (cse_p11_t27_g263 + cse_p11_t2_g238 + cse_p11_t2_g238);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 1, h_x1_y1);
    h_x1_z1 = (cse_p11_t29_g265 + cse_p11_t3_g239 + cse_p11_t3_g239);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 2, h_x1_z1);
    cse_p63_t1_g335 = (-2.0000000000000000     * cse_p529_t1_sqrt378 * cse_p529_t1_sqrt378 * cse_p529_t2_invr379 * cse_p529_t2_invr379 * kb);
    h_x1_x2 = (cse_p11_t10_g246 + cse_p11_t10_g246 + cse_p11_t42_g278 + cse_p11_t42_g278 + cse_p11_t51_g287 + cse_p11_t51_g287 + (cse_p11_t1_g308 * cse_p63_t1_g335));
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 0, h_x1_x2);
    h_x1_y2 = (cse_p114_t1_g340 + cse_p11_t11_g247 + cse_p11_t11_g247 + cse_p11_t43_g279 + cse_p11_t43_g279);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 1, h_x1_y2);
    h_x1_z2 = (cse_p11_t12_g248 + cse_p11_t12_g248 + cse_p11_t44_g280 + cse_p11_t44_g280 + cse_p63_t2_g332);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 2, h_x1_z2);
    h_y1_y1 = (cse_p11_t4_g240 + cse_p11_t4_g240 + cse_p12_t1_g318);
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 1, h_y1_y1);
    h_y1_z1 = (cse_p11_t5_g241 + cse_p11_t5_g241 + cse_p13_t1_g325);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 2, h_y1_z1);
    h_y1_x2 = h_x1_y2;
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 0, h_y1_x2);
    h_y1_y2 = (cse_p11_t13_g249 + cse_p11_t13_g249 + cse_p11_t45_g281 + cse_p11_t45_g281 + cse_p11_t51_g287 + cse_p11_t51_g287 + (cse_p11_t1_g306 * cse_p63_t1_g335));
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 1, h_y1_y2);
    h_y1_z2 = (cse_p11_t14_g250 + cse_p11_t14_g250 + cse_p11_t46_g282 + cse_p11_t46_g282 + cse_p63_t3_g333);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 2, h_y1_z2);
    h_z1_z1 = (cse_p11_t6_g242 + cse_p11_t6_g242 + cse_p12_t1_g319);
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x1, 2, h_z1_z1);
    h_z1_x2 = h_x1_z2;
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 0, h_z1_x2);
    h_z1_y2 = h_y1_z2;
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 1, h_z1_y2);
    h_z1_z2 = (cse_p11_t15_g251 + cse_p11_t15_g251 + cse_p11_t47_g283 + cse_p11_t47_g283 + cse_p11_t51_g287 + cse_p11_t51_g287 + (cse_p11_t1_g307 * cse_p63_t1_g335));
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 2, h_z1_z2);
    h_x2_x2 = h_x1_x1;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 0, h_x2_x2);
    h_x2_y2 = h_x1_y1;
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 1, h_x2_y2);
    h_x2_z2 = h_x1_z1;
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 2, h_x2_z2);
    h_y2_y2 = h_y1_y1;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 1, h_y2_y2);
    h_y2_z2 = h_y1_z1;
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 2, h_y2_z2);
    h_z2_z2 = h_z1_z1;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x2, 2, h_z2_z2);
  }
}
void hessian_fd(double kb, double r0, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 6;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  const double invh2 = 1.0/((h*h));
  double e0 = 0.0;
  energy(kb, r0, i3x1, i3x2, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 0, i3x1, 0, hval);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 1, i3x1, 1, hval);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 2, i3x1, 2, hval);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 0, i3x2, 0, hval);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 1, i3x2, 1, hval);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
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
    energy(kb, r0, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kb, r0, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(kb, r0, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kb, r0, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(kb, r0, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kb, r0, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(kb, r0, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kb, r0, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(kb, r0, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kb, r0, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(kb, r0, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(kb, r0, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(kb, r0, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kb, r0, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(kb, r0, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kb, r0, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(kb, r0, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(kb, r0, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(kb, r0, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(kb, r0, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(kb, r0, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kb, r0, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(kb, r0, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kb, r0, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(kb, r0, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(kb, r0, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(kb, r0, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(kb, r0, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
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
    energy(kb, r0, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(kb, r0, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kb, r0, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i; position[i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x2, 1, hval);
  }
}

void gradient(double kb, double r0, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  /* !BASE */
  double cse_p11_t1_g383;
  double cse_p11_t1_g384;
  double cse_p11_t1_g385;
  double cse_p1_t1_sqrt381;
  double cse_p1_t2_invsqrt382;
  double cse_p62_t1_g386;
  double dr;
  double dx;
  double dy;
  double dz;
  double energy;
  double g_x1;
  double g_x2;
  double g_y1;
  double g_y2;
  double g_z1;
  double g_z2;
  double r2;
  DOUBLE x1 = position[i3x1 + 0];
  DOUBLE y1 = position[i3x1 + 1];
  DOUBLE z1 = position[i3x1 + 2];
  DOUBLE x2 = position[i3x2 + 0];
  DOUBLE y2 = position[i3x2 + 1];
  DOUBLE z2 = position[i3x2 + 2];
  {
    /* !BASE */
    dx = (x2 + (-(x1)));
    dy = (y2 + (-(y1)));
    dz = (z2 + (-(z1)));
    r2 = ((dx * dx) + (dy * dy) + (dz * dz));
    cse_p1_t1_sqrt381 = sqrt(r2);
    cse_p1_t2_invsqrt382 = (1.0 / (cse_p1_t1_sqrt381));
    dr = (cse_p1_t1_sqrt381 + (-(r0)));
    cse_p11_t1_g385 = (dr * kb);
    energy = (cse_p11_t1_g385 * dr);
    *energy_accumulate += energy;
    cse_p11_t1_g383 = (cse_p11_t1_g385 * cse_p1_t2_invsqrt382);
    cse_p62_t1_g386 = (-2.0000000000000000     * cse_p11_t1_g383);
    g_x1 = (cse_p62_t1_g386 * dx);
    KernelGradientAcc(i3x1, 0, g_x1);
    g_y1 = (cse_p62_t1_g386 * dy);
    KernelGradientAcc(i3x1, 1, g_y1);
    g_z1 = (cse_p62_t1_g386 * dz);
    KernelGradientAcc(i3x1, 2, g_z1);
    cse_p11_t1_g384 = (2.0000000000000000     * cse_p11_t1_g383);
    g_x2 = (cse_p11_t1_g384 * dx);
    KernelGradientAcc(i3x2, 0, g_x2);
    g_y2 = (cse_p11_t1_g384 * dy);
    KernelGradientAcc(i3x2, 1, g_y2);
    g_z2 = (cse_p11_t1_g384 * dz);
    KernelGradientAcc(i3x2, 2, g_z2);
  }
}
void gradient_fd(double kb, double r0, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 6;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  double e0 = 0.0;
  energy(kb, r0, i3x1, i3x2, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(kb, r0, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(kb, r0, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
}

void energy(double kb, double r0, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  /* !BASE */
  double dr;
  double dx;
  double dy;
  double dz;
  double energy;
  double r;
  double r2;
  DOUBLE x1 = position[i3x1 + 0];
  DOUBLE y1 = position[i3x1 + 1];
  DOUBLE z1 = position[i3x1 + 2];
  DOUBLE x2 = position[i3x2 + 0];
  DOUBLE y2 = position[i3x2 + 1];
  DOUBLE z2 = position[i3x2 + 2];
  {
    /* !BASE */
    dx = (x2 + (-(x1)));
    dy = (y2 + (-(y1)));
    dz = (z2 + (-(z1)));
    r2 = ((dx * dx) + (dy * dy) + (dz * dz));
    r = sqrt(r2);
    dr = (r + (-(r0)));
    energy = (dr * dr * kb);
    *energy_accumulate += energy;
  }
}
void energy_fd(double kb, double r0, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  energy(kb, r0, i3x1, i3x2, position, energy_accumulate, force, hessian, dvec, hdvec);
}

};
