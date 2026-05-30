#include "stretch.h"

template <typename HESSIAN>
struct Stretch {
  static constexpr size_t PositionSize = 6;
  static std::string description() { return "mathkernel-stretch"; };
double hessian(const TermStretch& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
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
      double dx = (x2 + (-(x1)));
      double dy = (y2 + (-(y1)));
      double dz = (z2 + (-(z1)));
      double cse_p11_t1_g558439 = (dy * dy);
      double cse_p11_t1_g558440 = (dz * dz);
      double cse_p11_t1_g558441 = (dx * dx);
      double r2 = (cse_p11_t1_g558439 + cse_p11_t1_g558440 + cse_p11_t1_g558441);
      double cse_p529_t1_sqrt558507 = sqrt(r2);
      double cse_p529_t2_invr558508 = (1.0 / (r2));
      double cse_p61_t3_invr2558458 = (cse_p529_t2_invr558508 * cse_p529_t2_invr558508);
      double cse_p11_t1_g558433 = (cse_p529_t1_sqrt558507 * cse_p529_t2_invr558508);
      double cse_p63_t1_g558465 = (cse_p11_t1_g558433 * term.kb);
      double cse_p11_t1_g558435 = (cse_p529_t2_invr558508 * cse_p63_t1_g558465);
      double cse_p11_t1_g558438 = (cse_p11_t1_g558435 * dx);
      double cse_p12_t1_g558445 = (cse_p11_t1_g558438 * cse_p529_t1_sqrt558507);
      double cse_p11_t7_g558386 = (-2.00000000000000000e+0 * cse_p12_t1_g558445 * dy);
      double cse_p11_t8_g558387 = (-2.00000000000000000e+0 * cse_p12_t1_g558445 * dz);
      double cse_p11_t1_g558442 = (dy * dz);
      double cse_p12_t1_g558446 = (cse_p11_t1_g558435 * cse_p11_t1_g558442);
      double cse_p13_t1_g558455 = (cse_p12_t1_g558446 * cse_p529_t1_sqrt558507);
      double cse_p11_t9_g558388 = (2.00000000000000000e+0 * (-(cse_p13_t1_g558455)));
      double cse_p11_t26_g558405 = (cse_p11_t1_g558433 * cse_p11_t1_g558441);
      double cse_p11_t28_g558407 = (cse_p12_t1_g558445 * dy);
      double cse_p11_t30_g558409 = (cse_p12_t1_g558445 * dz);
      double cse_p11_t31_g558410 = (cse_p11_t1_g558433 * cse_p11_t1_g558439);
      double cse_p11_t34_g558413 = (cse_p11_t1_g558433 * cse_p11_t1_g558440);
      double cse_p11_t48_g558427 = (-(cse_p11_t28_g558407));
      double cse_p114_t1_g558469 = (2.00000000000000000e+0 * cse_p11_t48_g558427);
      double cse_p11_t49_g558428 = (-(cse_p11_t30_g558409));
      double cse_p63_t2_g558461 = (2.00000000000000000e+0 * cse_p11_t49_g558428);
      double cse_p11_t50_g558429 = (-(cse_p13_t1_g558455));
      double cse_p63_t3_g558462 = (2.00000000000000000e+0 * cse_p11_t50_g558429);
      double cse_p1_t2_invsqrt558379 = (1.0 / (cse_p529_t1_sqrt558507));
      double r = cse_p529_t1_sqrt558507;
      double dr = (r + (-(term.r0)));
      double cse_p62_t1_g558459 = (cse_p12_t1_g558446 * dr);
      double cse_p13_t1_g558453 = (dr * term.kb);
      double cse_p11_t1_g558434 = (cse_p13_t1_g558453 * cse_p529_t1_sqrt558507 * cse_p61_t3_invr2558458);
      double cse_p11_t1_g558437 = (-2.00000000000000000e+0 * cse_p11_t1_g558434);
      double cse_p11_t1_g558380 = (cse_p11_t1_g558437 * cse_p11_t1_g558441);
      double cse_p12_t1_g558451 = (cse_p11_t1_g558437 * dx);
      double cse_p11_t2_g558381 = (cse_p12_t1_g558451 * dy);
      double cse_p11_t3_g558382 = (cse_p12_t1_g558451 * dz);
      double cse_p11_t4_g558383 = (cse_p11_t1_g558437 * cse_p11_t1_g558439);
      double cse_p11_t5_g558384 = (cse_p11_t1_g558437 * cse_p11_t1_g558442);
      double cse_p11_t6_g558385 = (cse_p11_t1_g558437 * cse_p11_t1_g558440);
      double cse_p12_t1_g558444 = (2.00000000000000000e+0 * cse_p11_t1_g558434);
      double cse_p11_t10_g558389 = (cse_p11_t1_g558441 * cse_p12_t1_g558444);
      double cse_p12_t1_g558452 = (cse_p12_t1_g558444 * dx);
      double cse_p11_t11_g558390 = (cse_p12_t1_g558452 * dy);
      double cse_p11_t12_g558391 = (cse_p12_t1_g558452 * dz);
      double cse_p11_t13_g558392 = (cse_p11_t1_g558439 * cse_p12_t1_g558444);
      double cse_p11_t14_g558393 = (cse_p11_t1_g558442 * cse_p12_t1_g558444);
      double cse_p11_t15_g558394 = (cse_p11_t1_g558440 * cse_p12_t1_g558444);
      double cse_p11_t1_g558436 = (cse_p529_t2_invr558508 * dr);
      double cse_p11_t16_g558395 = (cse_p11_t1_g558436 * cse_p11_t1_g558441);
      double cse_p12_t1_g558450 = (cse_p11_t1_g558438 * dr);
      double cse_p11_t17_g558396 = (cse_p12_t1_g558450 * dx);
      double cse_p11_t18_g558397 = (cse_p12_t1_g558450 * dy);
      double cse_p11_t19_g558398 = (cse_p12_t1_g558450 * dz);
      double cse_p11_t20_g558399 = (cse_p11_t1_g558436 * cse_p11_t1_g558439);
      double cse_p11_t21_g558400 = (cse_p11_t1_g558435 * cse_p11_t1_g558439 * dr);
      double cse_p11_t23_g558402 = (cse_p11_t1_g558436 * cse_p11_t1_g558440);
      double cse_p11_t24_g558403 = (cse_p11_t1_g558435 * cse_p11_t1_g558440 * dr);
      double cse_p11_t25_g558404 = (cse_p11_t1_g558433 * cse_p13_t1_g558453);
      double cse_p12_t1_g558443 = (cse_p11_t1_g558438 * (cse_p529_t1_sqrt558507 + cse_p529_t1_sqrt558507 + dr + dr));
      double cse_p11_t27_g558406 = (cse_p12_t1_g558443 * dy);
      double cse_p11_t29_g558408 = (cse_p12_t1_g558443 * dz);
      double cse_p13_t1_g558454 = (cse_p12_t1_g558446 * (cse_p529_t1_sqrt558507 + cse_p529_t1_sqrt558507 + dr + dr));
      double cse_p12_t1_g558449 = (cse_p63_t1_g558465 * (cse_p11_t16_g558395 + cse_p11_t16_g558395 + cse_p11_t26_g558405 + cse_p11_t26_g558405 + dr + dr));
      double cse_p12_t1_g558447 = (cse_p63_t1_g558465 * (cse_p11_t20_g558399 + cse_p11_t20_g558399 + cse_p11_t31_g558410 + cse_p11_t31_g558410 + dr + dr));
      double cse_p12_t1_g558448 = (cse_p63_t1_g558465 * (cse_p11_t23_g558402 + cse_p11_t23_g558402 + cse_p11_t34_g558413 + cse_p11_t34_g558413 + dr + dr));
      double cse_p11_t38_g558417 = (cse_p11_t16_g558395 + cse_p11_t16_g558395 + cse_p11_t26_g558405 + cse_p11_t26_g558405 + dr + dr);
      double cse_p11_t39_g558418 = (cse_p11_t20_g558399 + cse_p11_t20_g558399 + cse_p11_t31_g558410 + cse_p11_t31_g558410 + dr + dr);
      double cse_p11_t40_g558419 = (cse_p11_t23_g558402 + cse_p11_t23_g558402 + cse_p11_t34_g558413 + cse_p11_t34_g558413 + dr + dr);
      double cse_p11_t41_g558420 = (cse_p529_t1_sqrt558507 + cse_p529_t1_sqrt558507 + dr + dr);
      double cse_p11_t42_g558421 = (-(cse_p11_t17_g558396));
      double cse_p11_t43_g558422 = (-(cse_p11_t18_g558397));
      double cse_p11_t44_g558423 = (-(cse_p11_t19_g558398));
      double cse_p11_t45_g558424 = (-(cse_p11_t21_g558400));
      double cse_p11_t46_g558425 = (-(cse_p62_t1_g558459));
      double cse_p11_t47_g558426 = (-(cse_p11_t24_g558403));
      double cse_p11_t51_g558430 = (-(cse_p11_t25_g558404));
      double energy = (cse_p13_t1_g558453 * dr);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p63_t1_g558463 = (cse_p13_t1_g558453 * cse_p1_t2_invsqrt558379);
      double cse_p166_t1_g558474 = (2.00000000000000000e+0 * cse_p63_t1_g558463);
      double cse_p63_t1_g558466 = (-2.00000000000000000e+0 * cse_p63_t1_g558463);
      double g_x1 = (cse_p63_t1_g558466 * dx);
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = (cse_p63_t1_g558466 * dy);
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = (cse_p63_t1_g558466 * dz);
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double g_x2 = (cse_p166_t1_g558474 * dx);
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (cse_p166_t1_g558474 * dy);
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (cse_p166_t1_g558474 * dz);
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double h_x1_x1 = (cse_p11_t1_g558380 + cse_p11_t1_g558380 + cse_p12_t1_g558449);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
      double h_x1_y1 = (cse_p11_t27_g558406 + cse_p11_t2_g558381 + cse_p11_t2_g558381);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
      double h_x1_z1 = (cse_p11_t29_g558408 + cse_p11_t3_g558382 + cse_p11_t3_g558382);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
      double cse_p63_t1_g558464 = (-2.00000000000000000e+0 * cse_p529_t1_sqrt558507 * cse_p529_t1_sqrt558507 * cse_p529_t2_invr558508 * cse_p529_t2_invr558508 * term.kb);
      double h_x1_x2 = (cse_p11_t10_g558389 + cse_p11_t10_g558389 + cse_p11_t42_g558421 + cse_p11_t42_g558421 + cse_p11_t51_g558430 + cse_p11_t51_g558430 + (cse_p11_t1_g558441 * cse_p63_t1_g558464));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
      double h_x1_y2 = (cse_p114_t1_g558469 + cse_p11_t11_g558390 + cse_p11_t11_g558390 + cse_p11_t43_g558422 + cse_p11_t43_g558422);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
      double h_x1_z2 = (cse_p11_t12_g558391 + cse_p11_t12_g558391 + cse_p11_t44_g558423 + cse_p11_t44_g558423 + cse_p63_t2_g558461);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
      double h_y1_y1 = (cse_p11_t4_g558383 + cse_p11_t4_g558383 + cse_p12_t1_g558447);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
      double h_y1_z1 = (cse_p11_t5_g558384 + cse_p11_t5_g558384 + cse_p13_t1_g558454);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
      double h_y1_x2 = h_x1_y2;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
      double h_y1_y2 = (cse_p11_t13_g558392 + cse_p11_t13_g558392 + cse_p11_t45_g558424 + cse_p11_t45_g558424 + cse_p11_t51_g558430 + cse_p11_t51_g558430 + (cse_p11_t1_g558439 * cse_p63_t1_g558464));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
      double h_y1_z2 = (cse_p11_t14_g558393 + cse_p11_t14_g558393 + cse_p11_t46_g558425 + cse_p11_t46_g558425 + cse_p63_t3_g558462);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
      double h_z1_z1 = (cse_p11_t6_g558385 + cse_p11_t6_g558385 + cse_p12_t1_g558448);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
      double h_z1_x2 = h_x1_z2;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
      double h_z1_y2 = h_y1_z2;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
      double h_z1_z2 = (cse_p11_t15_g558394 + cse_p11_t15_g558394 + cse_p11_t47_g558426 + cse_p11_t47_g558426 + cse_p11_t51_g558430 + cse_p11_t51_g558430 + (cse_p11_t1_g558440 * cse_p63_t1_g558464));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
      double h_x2_x2 = h_x1_x1;
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
      double h_x2_y2 = h_x1_y1;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
      double h_x2_z2 = h_x1_z1;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
      double h_y2_y2 = h_y1_y1;
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
      double h_y2_z2 = h_y1_z1;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
      double h_z2_z2 = h_z1_z1;
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
    }
  }
  return energy_added;
}
void hessian_fd(const TermStretch& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 6;
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
    double saved = position[term.i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x2 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x2, 0, d);
  }
  {
    double saved = position[term.i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x2 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x2, 1, d);
  }
  {
    double saved = position[term.i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x2 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x2, 2, d);
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
    double saved = position[term.i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x2 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x2 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 0, term.i3x2, 0, hval);
  }
  {
    double saved = position[term.i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x2 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x2 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 1, term.i3x2, 1, hval);
  }
  {
    double saved = position[term.i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x2 + 2] = saved - h;
    energy(term, position, &e_minus);
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
  {
    double saved_i = position[term.i3x2 + 0];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x2 + 0] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 0] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
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
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 0] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
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
    energy(term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 0] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_mm);
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
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 1] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
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
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 1] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
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
    energy(term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 1] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_mm);
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
    energy(term, position, &e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 1] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_mm);
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
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
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
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
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
    energy(term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_mm);
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
    energy(term, position, &e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_mm);
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
    energy(term, position, &e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x2 + 2] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 2, term.i3x2, 1, hval);
  }
}

double gradient(const TermStretch& term, double* position, double* energy_accumulate, double* force) {
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
      double dx = (x2 + (-(x1)));
      double dy = (y2 + (-(y1)));
      double dz = (z2 + (-(z1)));
      double r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      double cse_p1_t1_sqrt558510 = sqrt(r2);
      double cse_p1_t2_invsqrt558511 = (1.0 / (cse_p1_t1_sqrt558510));
      double r = cse_p1_t1_sqrt558510;
      double dr = (r + (-(term.r0)));
      double cse_p11_t1_g558514 = (dr * term.kb);
      double energy = (cse_p11_t1_g558514 * dr);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t1_g558512 = (cse_p11_t1_g558514 * cse_p1_t2_invsqrt558511);
      double cse_p62_t1_g558515 = (-2.00000000000000000e+0 * cse_p11_t1_g558512);
      double g_x1 = (cse_p62_t1_g558515 * dx);
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = (cse_p62_t1_g558515 * dy);
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = (cse_p62_t1_g558515 * dz);
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p11_t1_g558513 = (2.00000000000000000e+0 * cse_p11_t1_g558512);
      double g_x2 = (cse_p11_t1_g558513 * dx);
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (cse_p11_t1_g558513 * dy);
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (cse_p11_t1_g558513 * dz);
      KernelGradientAcc(term.i3x2, 2, g_z2);
    }
  }
  return energy_added;
}
void gradient_fd(const TermStretch& term, double* position, double* energy_accumulate, double* force)
{
  constexpr size_t PositionSize = 6;
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
  {
    double saved = position[term.i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x2 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x2, 0, d);
  }
  {
    double saved = position[term.i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x2 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x2, 1, d);
  }
  {
    double saved = position[term.i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x2 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x2 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x2, 2, d);
  }
}

double energy(const TermStretch& term, double* position, double* energy_accumulate) {
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
      double dx = (x2 + (-(x1)));
      double dy = (y2 + (-(y1)));
      double dz = (z2 + (-(z1)));
      double r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      double r = sqrt(r2);
      double dr = (r + (-(term.r0)));
      double energy = (dr * dr * term.kb);
      energy_added += energy;
      *energy_accumulate += energy;
    }
  }
  return energy_added;
}
void energy_fd(const TermStretch& term, double* position, double* energy_accumulate)
{
  energy(term, position, energy_accumulate);
}

};
