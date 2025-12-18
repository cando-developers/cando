void dihedral_gradient(double v, double n, double sinphase, double cosphase, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec) {
  constexpr size_t PositionSize = 12;
  /* !BASE */
  double c1x;
  double c1y;
  double c1z;
  double c1_sq;
  double c2x;
  double c2y;
  double c2z;
  double c2_sq;
  double cos_angle;
  double cos_nphi;
  double cse_p1_t10_g630;
  double cse_p1_t11_g631;
  double cse_p1_t12_g632;
  double cse_p1_t13_g633;
  double cse_p1_t14_g634;
  double cse_p1_t15_g635;
  double cse_p1_t16_g636;
  double cse_p1_t17_g637;
  double cse_p1_t18_g638;
  double cse_p1_t19_g639;
  double cse_p1_t1_g621;
  double cse_p1_t1_g653;
  double cse_p1_t1_g654;
  double cse_p1_t1_g655;
  double cse_p1_t1_g656;
  double cse_p1_t1_g657;
  double cse_p1_t20_g640;
  double cse_p1_t21_g641;
  double cse_p1_t22_g642;
  double cse_p1_t23_g643;
  double cse_p1_t24_g644;
  double cse_p1_t25_g645;
  double cse_p1_t26_g646;
  double cse_p1_t27_g647;
  double cse_p1_t28_g648;
  double cse_p1_t29_g649;
  double cse_p1_t2_g622;
  double cse_p1_t30_g650;
  double cse_p1_t31_g651;
  double cse_p1_t32_g652;
  double cse_p1_t3_g623;
  double cse_p1_t6_g626;
  double cse_p1_t7_g627;
  double cse_p1_t8_g628;
  double cse_p1_t9_g629;
  double dot12;
  double dot23;
  double energy;
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
  double nphi;
  double phi;
  double sin_angle;
  double sin_nphi;
  double t1;
  double t2;
  double v1x;
  double v1y;
  double v1z;
  double v2x;
  double v2y;
  double v2z;
  double v2_len;
  double v2_sq;
  double v3x;
  double v3y;
  double v3z;
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
    v1x = (x2 + (-(x1)));
    cse_p1_t20_g640 = (-(v1x));
    v1y = (y2 + (-(y1)));
    cse_p1_t21_g641 = (-(v1y));
    v1z = (z2 + (-(z1)));
    cse_p1_t22_g642 = (-(v1z));
    v2x = (x3 + (-(x2)));
    cse_p1_t8_g628 = (v1z * v2x);
    cse_p1_t23_g643 = (-(v2x));
    v2y = (y3 + (-(y2)));
    cse_p1_t6_g626 = (v1x * v2y);
    cse_p1_t24_g644 = (-(v2y));
    v2z = (z3 + (-(z2)));
    cse_p1_t7_g627 = (v1y * v2z);
    cse_p1_t25_g645 = (-(v2z));
    v3x = (x4 + (-(x3)));
    cse_p1_t11_g631 = (cse_p1_t23_g643 + (-(v3x)));
    cse_p1_t15_g635 = (v2x + v3x);
    cse_p1_t26_g646 = (-(v3x));
    v3y = (y4 + (-(y3)));
    cse_p1_t12_g632 = (cse_p1_t24_g644 + (-(v3y)));
    cse_p1_t16_g636 = (v2y + v3y);
    cse_p1_t27_g647 = (-(v3y));
    v3z = (z4 + (-(z3)));
    cse_p1_t13_g633 = (cse_p1_t25_g645 + (-(v3z)));
    cse_p1_t17_g637 = (v2z + v3z);
    cse_p1_t28_g648 = (-(v3z));
    c1x = (cse_p1_t7_g627 + (-((v1z * v2y))));
    c1y = (cse_p1_t8_g628 + (-((v1x * v2z))));
    c1z = (cse_p1_t6_g626 + (-((v1y * v2x))));
    c2x = ((v2y * v3z) + (-((v2z * v3y))));
    cse_p1_t1_g621 = (c2x * v1x);
    c2y = ((v2z * v3x) + (-((v2x * v3z))));
    cse_p1_t2_g622 = (c2y * v1y);
    c2z = ((v2x * v3y) + (-((v2y * v3x))));
    cse_p1_t3_g623 = (c2z * v1z);
    cse_p1_t9_g629 = (cse_p1_t1_g621 + cse_p1_t2_g622 + cse_p1_t3_g623);
    c1_sq = ((c1x * c1x) + (c1y * c1y) + (c1z * c1z));
    c2_sq = ((c2x * c2x) + (c2y * c2y) + (c2z * c2z));
    v2_sq = ((v2x * v2x) + (v2y * v2y) + (v2z * v2z));
    cse_p1_t32_g652 = pow(v2_sq, -0.50000000000000000    );
    v2_len = sqrt(v2_sq);
    dot12 = ((v1x * v2x) + (v1y * v2y) + (v1z * v2z));
    dot23 = ((v2x * v3x) + (v2y * v3y) + (v2z * v3z));
    t1 = (cse_p1_t9_g629 * v2_len);
    cse_p1_t19_g639 = (-(t1));
    cse_p1_t30_g650 = ((t1) * (t1));
    t2 = ((c1x * c2x) + (c1y * c2y) + (c1z * c2z));
    cse_p1_t14_g634 = (cse_p1_t30_g650 + ((t2) * (t2)));
    cse_p1_t29_g649 = pow(cse_p1_t14_g634, -1);
    cse_p1_t31_g651 = ((t2) * (t2));
    phi = atan2(t1, t2);
    nphi = (n * phi);
    cse_p1_t1_g656 = (cosphase * (-(sin(nphi))));
    cse_p1_t1_g657 = (sinphase * cos(nphi));
    cse_p1_t10_g630 = (cse_p1_t1_g656 + cse_p1_t1_g657);
    cse_p1_t18_g638 = (-(sin(nphi)));
    sin_nphi = sin(nphi);
    cos_nphi = cos(nphi);
    cos_angle = ((cosphase * cos_nphi) + (sinphase * sin_nphi));
    sin_angle = ((cosphase * sin_nphi) + (-((cos_nphi * sinphase))));
    energy = (v * (1.0000000000000000     + cos_angle));
    *energy_accumulate += energy;
    cse_p1_t1_g653 = (cse_p1_t10_g630 * cse_p1_t29_g649 * n * v);
    cse_p1_t1_g655 = (t2 * v2_len);
    g_x1 = (cse_p1_t1_g653 * ((cse_p1_t19_g639 * ((c2y * v2z) + (c2z * cse_p1_t24_g644))) + (cse_p1_t1_g655 * (-(c2x)))));
    KernelGradientAcc(i3x1, 0, g_x1);
    g_y1 = (cse_p1_t1_g653 * ((cse_p1_t19_g639 * ((c2x * cse_p1_t25_g645) + (c2z * v2x))) + (cse_p1_t1_g655 * (-(c2y)))));
    KernelGradientAcc(i3x1, 1, g_y1);
    g_z1 = (cse_p1_t1_g653 * ((cse_p1_t19_g639 * ((c2x * v2y) + (c2y * cse_p1_t23_g643))) + (cse_p1_t1_g655 * (-(c2z)))));
    KernelGradientAcc(i3x1, 2, g_z1);
    cse_p1_t1_g654 = (cse_p1_t32_g652 * cse_p1_t9_g629);
    g_x2 = (cse_p1_t1_g653 * ((cse_p1_t19_g639 * ((c1y * v3z) + (c1z * cse_p1_t27_g647) + (c2y * (cse_p1_t22_g642 + cse_p1_t25_g645)) + (c2z * (v1y + v2y)))) + (t2 * ((cse_p1_t1_g654 * cse_p1_t23_g643) + (v2_len * (c2x + (cse_p1_t27_g647 * v1z) + (v1y * v3z)))))));
    KernelGradientAcc(i3x2, 0, g_x2);
    g_y2 = (cse_p1_t1_g653 * ((cse_p1_t19_g639 * ((c1x * cse_p1_t28_g648) + (c1z * v3x) + (c2x * (v1z + v2z)) + (c2z * (cse_p1_t20_g640 + cse_p1_t23_g643)))) + (t2 * ((cse_p1_t1_g654 * cse_p1_t24_g644) + (v2_len * (c2y + (cse_p1_t28_g648 * v1x) + (v1z * v3x)))))));
    KernelGradientAcc(i3x2, 1, g_y2);
    g_z2 = (cse_p1_t1_g653 * ((cse_p1_t19_g639 * ((c1x * v3y) + (c1y * cse_p1_t26_g646) + (c2x * (cse_p1_t21_g641 + cse_p1_t24_g644)) + (c2y * (v1x + v2x)))) + (t2 * ((cse_p1_t1_g654 * cse_p1_t25_g645) + (v2_len * (c2z + (cse_p1_t26_g646 * v1y) + (v1x * v3y)))))));
    KernelGradientAcc(i3x2, 2, g_z2);
    g_x3 = (cse_p1_t1_g653 * ((cse_p1_t19_g639 * ((c1y * cse_p1_t13_g633) + (c1z * cse_p1_t16_g636) + (c2y * v1z) + (c2z * cse_p1_t21_g641))) + (t2 * ((cse_p1_t1_g654 * v2x) + (v2_len * ((cse_p1_t13_g633 * v1y) + (cse_p1_t16_g636 * v1z)))))));
    KernelGradientAcc(i3x3, 0, g_x3);
    g_y3 = (cse_p1_t1_g653 * ((cse_p1_t19_g639 * ((c1x * cse_p1_t17_g637) + (c1z * cse_p1_t11_g631) + (c2x * cse_p1_t22_g642) + (c2z * v1x))) + (t2 * ((cse_p1_t1_g654 * v2y) + (v2_len * ((cse_p1_t11_g631 * v1z) + (cse_p1_t17_g637 * v1x)))))));
    KernelGradientAcc(i3x3, 1, g_y3);
    g_z3 = (cse_p1_t1_g653 * ((cse_p1_t19_g639 * ((c1x * cse_p1_t12_g632) + (c1y * cse_p1_t15_g635) + (c2x * v1y) + (c2y * cse_p1_t20_g640))) + (t2 * ((cse_p1_t1_g654 * v2z) + (v2_len * ((cse_p1_t12_g632 * v1x) + (cse_p1_t15_g635 * v1y)))))));
    KernelGradientAcc(i3x3, 2, g_z3);
    g_x4 = (cse_p1_t1_g653 * ((cse_p1_t19_g639 * ((c1y * v2z) + (c1z * cse_p1_t24_g644))) + (cse_p1_t1_g655 * (cse_p1_t7_g627 + (cse_p1_t24_g644 * v1z)))));
    KernelGradientAcc(i3x4, 0, g_x4);
    g_y4 = (cse_p1_t1_g653 * ((cse_p1_t19_g639 * ((c1x * cse_p1_t25_g645) + (c1z * v2x))) + (cse_p1_t1_g655 * (cse_p1_t8_g628 + (cse_p1_t25_g645 * v1x)))));
    KernelGradientAcc(i3x4, 1, g_y4);
    g_z4 = (cse_p1_t1_g653 * ((cse_p1_t19_g639 * ((c1x * v2y) + (c1y * cse_p1_t23_g643))) + (cse_p1_t1_g655 * (cse_p1_t6_g626 + (cse_p1_t23_g643 * v1y)))));
    KernelGradientAcc(i3x4, 2, g_z4);
  }
}
void dihedral_gradient_fd(double v, double n, double sinphase, double cosphase, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 12;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  double e0 = 0.0;
  dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
  {
    double saved = position[i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 0] = saved + h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved - h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 0, d);
  }
  {
    double saved = position[i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 1] = saved + h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved - h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 1, d);
  }
  {
    double saved = position[i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 2] = saved + h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved - h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 2, d);
  }
  {
    double saved = position[i3x4 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 0] = saved + h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved - h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 0, d);
  }
  {
    double saved = position[i3x4 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 1] = saved + h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved - h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 1, d);
  }
  {
    double saved = position[i3x4 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 2] = saved + h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved - h;
    dihedral_energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 2, d);
  }
}

