template <typename HESSIAN>
struct Stretch {
  static constexpr size_t PositionSize = 6;
  static std::string description() { return "mathkernel-stretch"; };
double hessian(double kb, double r0, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
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
      double dx = (x2 + (-(x1)));
      double dy = (y2 + (-(y1)));
      double dz = (z2 + (-(z1)));
      double cse_p11_t1_g27495 = (dy * dy);
      double cse_p11_t1_g27496 = (dz * dz);
      double cse_p11_t1_g27497 = (dx * dx);
      double r2 = (cse_p11_t1_g27495 + cse_p11_t1_g27496 + cse_p11_t1_g27497);
      double cse_p529_t1_sqrt27563 = sqrt(r2);
      double cse_p529_t2_invr27564 = (1.0 / (r2));
      double cse_p61_t3_invr227514 = (cse_p529_t2_invr27564 * cse_p529_t2_invr27564);
      double cse_p11_t1_g27489 = (cse_p529_t1_sqrt27563 * cse_p529_t2_invr27564);
      double cse_p63_t1_g27521 = (cse_p11_t1_g27489 * kb);
      double cse_p11_t1_g27491 = (cse_p529_t2_invr27564 * cse_p63_t1_g27521);
      double cse_p11_t1_g27494 = (cse_p11_t1_g27491 * dx);
      double cse_p12_t1_g27501 = (cse_p11_t1_g27494 * cse_p529_t1_sqrt27563);
      double cse_p11_t7_g27442 = (-2.0000000000000000     * cse_p12_t1_g27501 * dy);
      double cse_p11_t8_g27443 = (-2.0000000000000000     * cse_p12_t1_g27501 * dz);
      double cse_p11_t1_g27498 = (dy * dz);
      double cse_p12_t1_g27502 = (cse_p11_t1_g27491 * cse_p11_t1_g27498);
      double cse_p13_t1_g27511 = (cse_p12_t1_g27502 * cse_p529_t1_sqrt27563);
      double cse_p11_t9_g27444 = (2.0000000000000000     * (-(cse_p13_t1_g27511)));
      double cse_p11_t26_g27461 = (cse_p11_t1_g27489 * cse_p11_t1_g27497);
      double cse_p11_t28_g27463 = (cse_p12_t1_g27501 * dy);
      double cse_p11_t30_g27465 = (cse_p12_t1_g27501 * dz);
      double cse_p11_t31_g27466 = (cse_p11_t1_g27489 * cse_p11_t1_g27495);
      double cse_p11_t34_g27469 = (cse_p11_t1_g27489 * cse_p11_t1_g27496);
      double cse_p11_t48_g27483 = (-(cse_p11_t28_g27463));
      double cse_p114_t1_g27525 = (2.0000000000000000     * cse_p11_t48_g27483);
      double cse_p11_t49_g27484 = (-(cse_p11_t30_g27465));
      double cse_p63_t2_g27517 = (2.0000000000000000     * cse_p11_t49_g27484);
      double cse_p11_t50_g27485 = (-(cse_p13_t1_g27511));
      double cse_p63_t3_g27518 = (2.0000000000000000     * cse_p11_t50_g27485);
      double cse_p1_t2_invsqrt27435 = (1.0 / (cse_p529_t1_sqrt27563));
      double dr = (cse_p529_t1_sqrt27563 + (-(r0)));
      double cse_p62_t1_g27515 = (cse_p12_t1_g27502 * dr);
      double cse_p13_t1_g27509 = (dr * kb);
      double cse_p11_t1_g27490 = (cse_p13_t1_g27509 * cse_p529_t1_sqrt27563 * cse_p61_t3_invr227514);
      double cse_p11_t1_g27493 = (-2.0000000000000000     * cse_p11_t1_g27490);
      double cse_p11_t1_g27436 = (cse_p11_t1_g27493 * cse_p11_t1_g27497);
      double cse_p12_t1_g27507 = (cse_p11_t1_g27493 * dx);
      double cse_p11_t2_g27437 = (cse_p12_t1_g27507 * dy);
      double cse_p11_t3_g27438 = (cse_p12_t1_g27507 * dz);
      double cse_p11_t4_g27439 = (cse_p11_t1_g27493 * cse_p11_t1_g27495);
      double cse_p11_t5_g27440 = (cse_p11_t1_g27493 * cse_p11_t1_g27498);
      double cse_p11_t6_g27441 = (cse_p11_t1_g27493 * cse_p11_t1_g27496);
      double cse_p12_t1_g27500 = (2.0000000000000000     * cse_p11_t1_g27490);
      double cse_p11_t10_g27445 = (cse_p11_t1_g27497 * cse_p12_t1_g27500);
      double cse_p12_t1_g27508 = (cse_p12_t1_g27500 * dx);
      double cse_p11_t11_g27446 = (cse_p12_t1_g27508 * dy);
      double cse_p11_t12_g27447 = (cse_p12_t1_g27508 * dz);
      double cse_p11_t13_g27448 = (cse_p11_t1_g27495 * cse_p12_t1_g27500);
      double cse_p11_t14_g27449 = (cse_p11_t1_g27498 * cse_p12_t1_g27500);
      double cse_p11_t15_g27450 = (cse_p11_t1_g27496 * cse_p12_t1_g27500);
      double cse_p11_t1_g27492 = (cse_p529_t2_invr27564 * dr);
      double cse_p11_t16_g27451 = (cse_p11_t1_g27492 * cse_p11_t1_g27497);
      double cse_p12_t1_g27506 = (cse_p11_t1_g27494 * dr);
      double cse_p11_t17_g27452 = (cse_p12_t1_g27506 * dx);
      double cse_p11_t18_g27453 = (cse_p12_t1_g27506 * dy);
      double cse_p11_t19_g27454 = (cse_p12_t1_g27506 * dz);
      double cse_p11_t20_g27455 = (cse_p11_t1_g27492 * cse_p11_t1_g27495);
      double cse_p11_t21_g27456 = (cse_p11_t1_g27491 * cse_p11_t1_g27495 * dr);
      double cse_p11_t23_g27458 = (cse_p11_t1_g27492 * cse_p11_t1_g27496);
      double cse_p11_t24_g27459 = (cse_p11_t1_g27491 * cse_p11_t1_g27496 * dr);
      double cse_p11_t25_g27460 = (cse_p11_t1_g27489 * cse_p13_t1_g27509);
      double cse_p12_t1_g27499 = (cse_p11_t1_g27494 * (cse_p529_t1_sqrt27563 + cse_p529_t1_sqrt27563 + dr + dr));
      double cse_p11_t27_g27462 = (cse_p12_t1_g27499 * dy);
      double cse_p11_t29_g27464 = (cse_p12_t1_g27499 * dz);
      double cse_p13_t1_g27510 = (cse_p12_t1_g27502 * (cse_p529_t1_sqrt27563 + cse_p529_t1_sqrt27563 + dr + dr));
      double cse_p12_t1_g27505 = (cse_p63_t1_g27521 * (cse_p11_t16_g27451 + cse_p11_t16_g27451 + cse_p11_t26_g27461 + cse_p11_t26_g27461 + dr + dr));
      double cse_p12_t1_g27503 = (cse_p63_t1_g27521 * (cse_p11_t20_g27455 + cse_p11_t20_g27455 + cse_p11_t31_g27466 + cse_p11_t31_g27466 + dr + dr));
      double cse_p12_t1_g27504 = (cse_p63_t1_g27521 * (cse_p11_t23_g27458 + cse_p11_t23_g27458 + cse_p11_t34_g27469 + cse_p11_t34_g27469 + dr + dr));
      double cse_p11_t38_g27473 = (cse_p11_t16_g27451 + cse_p11_t16_g27451 + cse_p11_t26_g27461 + cse_p11_t26_g27461 + dr + dr);
      double cse_p11_t39_g27474 = (cse_p11_t20_g27455 + cse_p11_t20_g27455 + cse_p11_t31_g27466 + cse_p11_t31_g27466 + dr + dr);
      double cse_p11_t40_g27475 = (cse_p11_t23_g27458 + cse_p11_t23_g27458 + cse_p11_t34_g27469 + cse_p11_t34_g27469 + dr + dr);
      double cse_p11_t41_g27476 = (cse_p529_t1_sqrt27563 + cse_p529_t1_sqrt27563 + dr + dr);
      double cse_p11_t42_g27477 = (-(cse_p11_t17_g27452));
      double cse_p11_t43_g27478 = (-(cse_p11_t18_g27453));
      double cse_p11_t44_g27479 = (-(cse_p11_t19_g27454));
      double cse_p11_t45_g27480 = (-(cse_p11_t21_g27456));
      double cse_p11_t46_g27481 = (-(cse_p62_t1_g27515));
      double cse_p11_t47_g27482 = (-(cse_p11_t24_g27459));
      double cse_p11_t51_g27486 = (-(cse_p11_t25_g27460));
      double energy = (cse_p13_t1_g27509 * dr);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p63_t1_g27519 = (cse_p13_t1_g27509 * cse_p1_t2_invsqrt27435);
      double cse_p166_t1_g27530 = (2.0000000000000000     * cse_p63_t1_g27519);
      double cse_p63_t1_g27522 = (-2.0000000000000000     * cse_p63_t1_g27519);
      double g_x1 = (cse_p63_t1_g27522 * dx);
      KernelGradientAcc(i3x1, 0, g_x1);
      double g_y1 = (cse_p63_t1_g27522 * dy);
      KernelGradientAcc(i3x1, 1, g_y1);
      double g_z1 = (cse_p63_t1_g27522 * dz);
      KernelGradientAcc(i3x1, 2, g_z1);
      double g_x2 = (cse_p166_t1_g27530 * dx);
      KernelGradientAcc(i3x2, 0, g_x2);
      double g_y2 = (cse_p166_t1_g27530 * dy);
      KernelGradientAcc(i3x2, 1, g_y2);
      double g_z2 = (cse_p166_t1_g27530 * dz);
      KernelGradientAcc(i3x2, 2, g_z2);
      double h_x1_x1 = (cse_p11_t1_g27436 + cse_p11_t1_g27436 + cse_p12_t1_g27505);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 0, h_x1_x1);
      double h_x1_y1 = (cse_p11_t27_g27462 + cse_p11_t2_g27437 + cse_p11_t2_g27437);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 1, h_x1_y1);
      double h_x1_z1 = (cse_p11_t29_g27464 + cse_p11_t3_g27438 + cse_p11_t3_g27438);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 2, h_x1_z1);
      double cse_p63_t1_g27520 = (-2.0000000000000000     * cse_p529_t1_sqrt27563 * cse_p529_t1_sqrt27563 * cse_p529_t2_invr27564 * cse_p529_t2_invr27564 * kb);
      double h_x1_x2 = (cse_p11_t10_g27445 + cse_p11_t10_g27445 + cse_p11_t42_g27477 + cse_p11_t42_g27477 + cse_p11_t51_g27486 + cse_p11_t51_g27486 + (cse_p11_t1_g27497 * cse_p63_t1_g27520));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 0, h_x1_x2);
      double h_x1_y2 = (cse_p114_t1_g27525 + cse_p11_t11_g27446 + cse_p11_t11_g27446 + cse_p11_t43_g27478 + cse_p11_t43_g27478);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 1, h_x1_y2);
      double h_x1_z2 = (cse_p11_t12_g27447 + cse_p11_t12_g27447 + cse_p11_t44_g27479 + cse_p11_t44_g27479 + cse_p63_t2_g27517);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 2, h_x1_z2);
      double h_y1_y1 = (cse_p11_t4_g27439 + cse_p11_t4_g27439 + cse_p12_t1_g27503);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 1, h_y1_y1);
      double h_y1_z1 = (cse_p11_t5_g27440 + cse_p11_t5_g27440 + cse_p13_t1_g27510);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 2, h_y1_z1);
      double h_y1_x2 = h_x1_y2;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 0, h_y1_x2);
      double h_y1_y2 = (cse_p11_t13_g27448 + cse_p11_t13_g27448 + cse_p11_t45_g27480 + cse_p11_t45_g27480 + cse_p11_t51_g27486 + cse_p11_t51_g27486 + (cse_p11_t1_g27495 * cse_p63_t1_g27520));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 1, h_y1_y2);
      double h_y1_z2 = (cse_p11_t14_g27449 + cse_p11_t14_g27449 + cse_p11_t46_g27481 + cse_p11_t46_g27481 + cse_p63_t3_g27518);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 2, h_y1_z2);
      double h_z1_z1 = (cse_p11_t6_g27441 + cse_p11_t6_g27441 + cse_p12_t1_g27504);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x1, 2, h_z1_z1);
      double h_z1_x2 = h_x1_z2;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 0, h_z1_x2);
      double h_z1_y2 = h_y1_z2;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 1, h_z1_y2);
      double h_z1_z2 = (cse_p11_t15_g27450 + cse_p11_t15_g27450 + cse_p11_t47_g27482 + cse_p11_t47_g27482 + cse_p11_t51_g27486 + cse_p11_t51_g27486 + (cse_p11_t1_g27496 * cse_p63_t1_g27520));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 2, h_z1_z2);
      double h_x2_x2 = h_x1_x1;
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 0, h_x2_x2);
      double h_x2_y2 = h_x1_y1;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 1, h_x2_y2);
      double h_x2_z2 = h_x1_z1;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 2, h_x2_z2);
      double h_y2_y2 = h_y1_y1;
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 1, h_y2_y2);
      double h_y2_z2 = h_y1_z1;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 2, h_y2_z2);
      double h_z2_z2 = h_z1_z1;
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x2, 2, h_z2_z2);
    }
  }
  return energy_added;
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

double gradient(double kb, double r0, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
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
      double dx = (x2 + (-(x1)));
      double dy = (y2 + (-(y1)));
      double dz = (z2 + (-(z1)));
      double r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      double cse_p1_t1_sqrt27566 = sqrt(r2);
      double cse_p1_t2_invsqrt27567 = (1.0 / (cse_p1_t1_sqrt27566));
      double dr = (cse_p1_t1_sqrt27566 + (-(r0)));
      double cse_p11_t1_g27570 = (dr * kb);
      double energy = (cse_p11_t1_g27570 * dr);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t1_g27568 = (cse_p11_t1_g27570 * cse_p1_t2_invsqrt27567);
      double cse_p62_t1_g27571 = (-2.0000000000000000     * cse_p11_t1_g27568);
      double g_x1 = (cse_p62_t1_g27571 * dx);
      KernelGradientAcc(i3x1, 0, g_x1);
      double g_y1 = (cse_p62_t1_g27571 * dy);
      KernelGradientAcc(i3x1, 1, g_y1);
      double g_z1 = (cse_p62_t1_g27571 * dz);
      KernelGradientAcc(i3x1, 2, g_z1);
      double cse_p11_t1_g27569 = (2.0000000000000000     * cse_p11_t1_g27568);
      double g_x2 = (cse_p11_t1_g27569 * dx);
      KernelGradientAcc(i3x2, 0, g_x2);
      double g_y2 = (cse_p11_t1_g27569 * dy);
      KernelGradientAcc(i3x2, 1, g_y2);
      double g_z2 = (cse_p11_t1_g27569 * dz);
      KernelGradientAcc(i3x2, 2, g_z2);
    }
  }
  return energy_added;
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

double energy(double kb, double r0, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
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
      double dx = (x2 + (-(x1)));
      double dy = (y2 + (-(y1)));
      double dz = (z2 + (-(z1)));
      double r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      double r = sqrt(r2);
      double dr = (r + (-(r0)));
      double energy = (dr * dr * kb);
      energy_added += energy;
      *energy_accumulate += energy;
    }
  }
  return energy_added;
}
void energy_fd(double kb, double r0, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  energy(kb, r0, i3x1, i3x2, position, energy_accumulate, force, hessian, dvec, hdvec);
}

};
