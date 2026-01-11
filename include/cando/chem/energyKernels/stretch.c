template <typename HESSIAN>
struct Stretch {
  static constexpr size_t PositionSize = 6;
  static std::string description() { return "mathkernel-stretch"; };
double hessian(double kb, double r0, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double cse_p114_t1_g34704;
  double cse_p11_t10_g34610;
  double cse_p11_t11_g34611;
  double cse_p11_t12_g34612;
  double cse_p11_t13_g34613;
  double cse_p11_t14_g34614;
  double cse_p11_t15_g34615;
  double cse_p11_t16_g34616;
  double cse_p11_t17_g34617;
  double cse_p11_t18_g34618;
  double cse_p11_t19_g34619;
  double cse_p11_t1_g34601;
  double cse_p11_t1_g34660;
  double cse_p11_t1_g34665;
  double cse_p11_t1_g34666;
  double cse_p11_t1_g34667;
  double cse_p11_t1_g34668;
  double cse_p11_t1_g34669;
  double cse_p11_t1_g34670;
  double cse_p11_t1_g34671;
  double cse_p11_t1_g34672;
  double cse_p11_t1_g34673;
  double cse_p11_t20_g34620;
  double cse_p11_t21_g34621;
  double cse_p11_t23_g34623;
  double cse_p11_t24_g34624;
  double cse_p11_t25_g34625;
  double cse_p11_t26_g34626;
  double cse_p11_t27_g34627;
  double cse_p11_t28_g34628;
  double cse_p11_t29_g34629;
  double cse_p11_t2_g34602;
  double cse_p11_t30_g34630;
  double cse_p11_t31_g34631;
  double cse_p11_t34_g34634;
  double cse_p11_t38_g34638;
  double cse_p11_t39_g34639;
  double cse_p11_t3_g34603;
  double cse_p11_t40_g34640;
  double cse_p11_t41_g34641;
  double cse_p11_t42_g34642;
  double cse_p11_t43_g34643;
  double cse_p11_t44_g34644;
  double cse_p11_t45_g34645;
  double cse_p11_t46_g34646;
  double cse_p11_t47_g34647;
  double cse_p11_t48_g34648;
  double cse_p11_t49_g34649;
  double cse_p11_t4_g34604;
  double cse_p11_t50_g34650;
  double cse_p11_t51_g34651;
  double cse_p11_t5_g34605;
  double cse_p11_t6_g34606;
  double cse_p11_t7_g34607;
  double cse_p11_t8_g34608;
  double cse_p11_t9_g34609;
  double cse_p12_t1_g34678;
  double cse_p12_t1_g34679;
  double cse_p12_t1_g34680;
  double cse_p12_t1_g34681;
  double cse_p12_t1_g34682;
  double cse_p12_t1_g34683;
  double cse_p12_t1_g34684;
  double cse_p12_t1_g34685;
  double cse_p12_t1_g34686;
  double cse_p12_t1_g34687;
  double cse_p13_t1_g34688;
  double cse_p13_t1_g34689;
  double cse_p13_t1_g34690;
  double cse_p166_t1_g34709;
  double cse_p1_t2_invsqrt34590;
  double cse_p529_t1_sqrt34742;
  double cse_p529_t2_invr34743;
  double cse_p61_t3_invr234693;
  double cse_p62_t1_g34694;
  double cse_p63_t1_g34698;
  double cse_p63_t1_g34699;
  double cse_p63_t1_g34700;
  double cse_p63_t1_g34701;
  double cse_p63_t2_g34696;
  double cse_p63_t3_g34697;
  double dr;
  double dx;
  double dy;
  double dz;
  double energy;
  double energy_added;
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
      dx = (x2 + (-(x1)));
      dy = (y2 + (-(y1)));
      dz = (z2 + (-(z1)));
      cse_p11_t1_g34670 = (dy * dy);
      cse_p11_t1_g34671 = (dz * dz);
      cse_p11_t1_g34672 = (dx * dx);
      r2 = (cse_p11_t1_g34670 + cse_p11_t1_g34671 + cse_p11_t1_g34672);
      cse_p529_t1_sqrt34742 = sqrt(r2);
      cse_p529_t2_invr34743 = (1.0 / (r2));
      cse_p61_t3_invr234693 = (cse_p529_t2_invr34743 * cse_p529_t2_invr34743);
      cse_p11_t1_g34660 = (cse_p529_t1_sqrt34742 * cse_p529_t2_invr34743);
      cse_p63_t1_g34700 = (cse_p11_t1_g34660 * kb);
      cse_p11_t1_g34666 = (cse_p529_t2_invr34743 * cse_p63_t1_g34700);
      cse_p11_t1_g34669 = (cse_p11_t1_g34666 * dx);
      cse_p12_t1_g34680 = (cse_p11_t1_g34669 * cse_p529_t1_sqrt34742);
      cse_p11_t7_g34607 = (-2.0000000000000000     * cse_p12_t1_g34680 * dy);
      cse_p11_t8_g34608 = (-2.0000000000000000     * cse_p12_t1_g34680 * dz);
      cse_p11_t1_g34673 = (dy * dz);
      cse_p12_t1_g34681 = (cse_p11_t1_g34666 * cse_p11_t1_g34673);
      cse_p13_t1_g34690 = (cse_p12_t1_g34681 * cse_p529_t1_sqrt34742);
      cse_p11_t9_g34609 = (2.0000000000000000     * (-(cse_p13_t1_g34690)));
      cse_p11_t26_g34626 = (cse_p11_t1_g34660 * cse_p11_t1_g34672);
      cse_p11_t28_g34628 = (cse_p12_t1_g34680 * dy);
      cse_p11_t30_g34630 = (cse_p12_t1_g34680 * dz);
      cse_p11_t31_g34631 = (cse_p11_t1_g34660 * cse_p11_t1_g34670);
      cse_p11_t34_g34634 = (cse_p11_t1_g34660 * cse_p11_t1_g34671);
      cse_p11_t48_g34648 = (-(cse_p11_t28_g34628));
      cse_p114_t1_g34704 = (2.0000000000000000     * cse_p11_t48_g34648);
      cse_p11_t49_g34649 = (-(cse_p11_t30_g34630));
      cse_p63_t2_g34696 = (2.0000000000000000     * cse_p11_t49_g34649);
      cse_p11_t50_g34650 = (-(cse_p13_t1_g34690));
      cse_p63_t3_g34697 = (2.0000000000000000     * cse_p11_t50_g34650);
      cse_p1_t2_invsqrt34590 = (1.0 / (cse_p529_t1_sqrt34742));
      dr = (cse_p529_t1_sqrt34742 + (-(r0)));
      cse_p62_t1_g34694 = (cse_p12_t1_g34681 * dr);
      cse_p13_t1_g34688 = (dr * kb);
      cse_p11_t1_g34665 = (cse_p13_t1_g34688 * cse_p529_t1_sqrt34742 * cse_p61_t3_invr234693);
      cse_p11_t1_g34668 = (-2.0000000000000000     * cse_p11_t1_g34665);
      cse_p11_t1_g34601 = (cse_p11_t1_g34668 * cse_p11_t1_g34672);
      cse_p12_t1_g34686 = (cse_p11_t1_g34668 * dx);
      cse_p11_t2_g34602 = (cse_p12_t1_g34686 * dy);
      cse_p11_t3_g34603 = (cse_p12_t1_g34686 * dz);
      cse_p11_t4_g34604 = (cse_p11_t1_g34668 * cse_p11_t1_g34670);
      cse_p11_t5_g34605 = (cse_p11_t1_g34668 * cse_p11_t1_g34673);
      cse_p11_t6_g34606 = (cse_p11_t1_g34668 * cse_p11_t1_g34671);
      cse_p12_t1_g34679 = (2.0000000000000000     * cse_p11_t1_g34665);
      cse_p11_t10_g34610 = (cse_p11_t1_g34672 * cse_p12_t1_g34679);
      cse_p12_t1_g34687 = (cse_p12_t1_g34679 * dx);
      cse_p11_t11_g34611 = (cse_p12_t1_g34687 * dy);
      cse_p11_t12_g34612 = (cse_p12_t1_g34687 * dz);
      cse_p11_t13_g34613 = (cse_p11_t1_g34670 * cse_p12_t1_g34679);
      cse_p11_t14_g34614 = (cse_p11_t1_g34673 * cse_p12_t1_g34679);
      cse_p11_t15_g34615 = (cse_p11_t1_g34671 * cse_p12_t1_g34679);
      cse_p11_t1_g34667 = (cse_p529_t2_invr34743 * dr);
      cse_p11_t16_g34616 = (cse_p11_t1_g34667 * cse_p11_t1_g34672);
      cse_p12_t1_g34685 = (cse_p11_t1_g34669 * dr);
      cse_p11_t17_g34617 = (cse_p12_t1_g34685 * dx);
      cse_p11_t18_g34618 = (cse_p12_t1_g34685 * dy);
      cse_p11_t19_g34619 = (cse_p12_t1_g34685 * dz);
      cse_p11_t20_g34620 = (cse_p11_t1_g34667 * cse_p11_t1_g34670);
      cse_p11_t21_g34621 = (cse_p11_t1_g34666 * cse_p11_t1_g34670 * dr);
      cse_p11_t23_g34623 = (cse_p11_t1_g34667 * cse_p11_t1_g34671);
      cse_p11_t24_g34624 = (cse_p11_t1_g34666 * cse_p11_t1_g34671 * dr);
      cse_p11_t25_g34625 = (cse_p11_t1_g34660 * cse_p13_t1_g34688);
      cse_p12_t1_g34678 = (cse_p11_t1_g34669 * (cse_p529_t1_sqrt34742 + cse_p529_t1_sqrt34742 + dr + dr));
      cse_p11_t27_g34627 = (cse_p12_t1_g34678 * dy);
      cse_p11_t29_g34629 = (cse_p12_t1_g34678 * dz);
      cse_p13_t1_g34689 = (cse_p12_t1_g34681 * (cse_p529_t1_sqrt34742 + cse_p529_t1_sqrt34742 + dr + dr));
      cse_p12_t1_g34684 = (cse_p63_t1_g34700 * (cse_p11_t16_g34616 + cse_p11_t16_g34616 + cse_p11_t26_g34626 + cse_p11_t26_g34626 + dr + dr));
      cse_p12_t1_g34682 = (cse_p63_t1_g34700 * (cse_p11_t20_g34620 + cse_p11_t20_g34620 + cse_p11_t31_g34631 + cse_p11_t31_g34631 + dr + dr));
      cse_p12_t1_g34683 = (cse_p63_t1_g34700 * (cse_p11_t23_g34623 + cse_p11_t23_g34623 + cse_p11_t34_g34634 + cse_p11_t34_g34634 + dr + dr));
      cse_p11_t38_g34638 = (cse_p11_t16_g34616 + cse_p11_t16_g34616 + cse_p11_t26_g34626 + cse_p11_t26_g34626 + dr + dr);
      cse_p11_t39_g34639 = (cse_p11_t20_g34620 + cse_p11_t20_g34620 + cse_p11_t31_g34631 + cse_p11_t31_g34631 + dr + dr);
      cse_p11_t40_g34640 = (cse_p11_t23_g34623 + cse_p11_t23_g34623 + cse_p11_t34_g34634 + cse_p11_t34_g34634 + dr + dr);
      cse_p11_t41_g34641 = (cse_p529_t1_sqrt34742 + cse_p529_t1_sqrt34742 + dr + dr);
      cse_p11_t42_g34642 = (-(cse_p11_t17_g34617));
      cse_p11_t43_g34643 = (-(cse_p11_t18_g34618));
      cse_p11_t44_g34644 = (-(cse_p11_t19_g34619));
      cse_p11_t45_g34645 = (-(cse_p11_t21_g34621));
      cse_p11_t46_g34646 = (-(cse_p62_t1_g34694));
      cse_p11_t47_g34647 = (-(cse_p11_t24_g34624));
      cse_p11_t51_g34651 = (-(cse_p11_t25_g34625));
      energy = (cse_p13_t1_g34688 * dr);
      energy_added += energy;
      *energy_accumulate += energy;
      cse_p63_t1_g34698 = (cse_p13_t1_g34688 * cse_p1_t2_invsqrt34590);
      cse_p166_t1_g34709 = (2.0000000000000000     * cse_p63_t1_g34698);
      cse_p63_t1_g34701 = (-2.0000000000000000     * cse_p63_t1_g34698);
      g_x1 = (cse_p63_t1_g34701 * dx);
      KernelGradientAcc(i3x1, 0, g_x1);
      g_y1 = (cse_p63_t1_g34701 * dy);
      KernelGradientAcc(i3x1, 1, g_y1);
      g_z1 = (cse_p63_t1_g34701 * dz);
      KernelGradientAcc(i3x1, 2, g_z1);
      g_x2 = (cse_p166_t1_g34709 * dx);
      KernelGradientAcc(i3x2, 0, g_x2);
      g_y2 = (cse_p166_t1_g34709 * dy);
      KernelGradientAcc(i3x2, 1, g_y2);
      g_z2 = (cse_p166_t1_g34709 * dz);
      KernelGradientAcc(i3x2, 2, g_z2);
      h_x1_x1 = (cse_p11_t1_g34601 + cse_p11_t1_g34601 + cse_p12_t1_g34684);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 0, h_x1_x1);
      h_x1_y1 = (cse_p11_t27_g34627 + cse_p11_t2_g34602 + cse_p11_t2_g34602);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 1, h_x1_y1);
      h_x1_z1 = (cse_p11_t29_g34629 + cse_p11_t3_g34603 + cse_p11_t3_g34603);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 2, h_x1_z1);
      cse_p63_t1_g34699 = (-2.0000000000000000     * cse_p529_t1_sqrt34742 * cse_p529_t1_sqrt34742 * cse_p529_t2_invr34743 * cse_p529_t2_invr34743 * kb);
      h_x1_x2 = (cse_p11_t10_g34610 + cse_p11_t10_g34610 + cse_p11_t42_g34642 + cse_p11_t42_g34642 + cse_p11_t51_g34651 + cse_p11_t51_g34651 + (cse_p11_t1_g34672 * cse_p63_t1_g34699));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 0, h_x1_x2);
      h_x1_y2 = (cse_p114_t1_g34704 + cse_p11_t11_g34611 + cse_p11_t11_g34611 + cse_p11_t43_g34643 + cse_p11_t43_g34643);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 1, h_x1_y2);
      h_x1_z2 = (cse_p11_t12_g34612 + cse_p11_t12_g34612 + cse_p11_t44_g34644 + cse_p11_t44_g34644 + cse_p63_t2_g34696);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 2, h_x1_z2);
      h_y1_y1 = (cse_p11_t4_g34604 + cse_p11_t4_g34604 + cse_p12_t1_g34682);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 1, h_y1_y1);
      h_y1_z1 = (cse_p11_t5_g34605 + cse_p11_t5_g34605 + cse_p13_t1_g34689);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 2, h_y1_z1);
      h_y1_x2 = h_x1_y2;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 0, h_y1_x2);
      h_y1_y2 = (cse_p11_t13_g34613 + cse_p11_t13_g34613 + cse_p11_t45_g34645 + cse_p11_t45_g34645 + cse_p11_t51_g34651 + cse_p11_t51_g34651 + (cse_p11_t1_g34670 * cse_p63_t1_g34699));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 1, h_y1_y2);
      h_y1_z2 = (cse_p11_t14_g34614 + cse_p11_t14_g34614 + cse_p11_t46_g34646 + cse_p11_t46_g34646 + cse_p63_t3_g34697);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 2, h_y1_z2);
      h_z1_z1 = (cse_p11_t6_g34606 + cse_p11_t6_g34606 + cse_p12_t1_g34683);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x1, 2, h_z1_z1);
      h_z1_x2 = h_x1_z2;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 0, h_z1_x2);
      h_z1_y2 = h_y1_z2;
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 1, h_z1_y2);
      h_z1_z2 = (cse_p11_t15_g34615 + cse_p11_t15_g34615 + cse_p11_t47_g34647 + cse_p11_t47_g34647 + cse_p11_t51_g34651 + cse_p11_t51_g34651 + (cse_p11_t1_g34671 * cse_p63_t1_g34699));
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
  double cse_p11_t1_g34747;
  double cse_p11_t1_g34748;
  double cse_p11_t1_g34749;
  double cse_p1_t1_sqrt34745;
  double cse_p1_t2_invsqrt34746;
  double cse_p62_t1_g34750;
  double dr;
  double dx;
  double dy;
  double dz;
  double energy;
  double energy_added;
  double g_x1;
  double g_x2;
  double g_y1;
  double g_y2;
  double g_z1;
  double g_z2;
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
      dx = (x2 + (-(x1)));
      dy = (y2 + (-(y1)));
      dz = (z2 + (-(z1)));
      r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      cse_p1_t1_sqrt34745 = sqrt(r2);
      cse_p1_t2_invsqrt34746 = (1.0 / (cse_p1_t1_sqrt34745));
      dr = (cse_p1_t1_sqrt34745 + (-(r0)));
      cse_p11_t1_g34749 = (dr * kb);
      energy = (cse_p11_t1_g34749 * dr);
      energy_added += energy;
      *energy_accumulate += energy;
      cse_p11_t1_g34747 = (cse_p11_t1_g34749 * cse_p1_t2_invsqrt34746);
      cse_p62_t1_g34750 = (-2.0000000000000000     * cse_p11_t1_g34747);
      g_x1 = (cse_p62_t1_g34750 * dx);
      KernelGradientAcc(i3x1, 0, g_x1);
      g_y1 = (cse_p62_t1_g34750 * dy);
      KernelGradientAcc(i3x1, 1, g_y1);
      g_z1 = (cse_p62_t1_g34750 * dz);
      KernelGradientAcc(i3x1, 2, g_z1);
      cse_p11_t1_g34748 = (2.0000000000000000     * cse_p11_t1_g34747);
      g_x2 = (cse_p11_t1_g34748 * dx);
      KernelGradientAcc(i3x2, 0, g_x2);
      g_y2 = (cse_p11_t1_g34748 * dy);
      KernelGradientAcc(i3x2, 1, g_y2);
      g_z2 = (cse_p11_t1_g34748 * dz);
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
  double dr;
  double dx;
  double dy;
  double dz;
  double energy;
  double energy_added;
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
      dx = (x2 + (-(x1)));
      dy = (y2 + (-(y1)));
      dz = (z2 + (-(z1)));
      r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      r = sqrt(r2);
      dr = (r + (-(r0)));
      energy = (dr * dr * kb);
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
