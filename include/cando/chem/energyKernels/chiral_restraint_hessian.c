void chiral_restraint_hessian(double k, double co, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec) {
  constexpr size_t PositionSize = 12;
  /* !BASE */
  double cse_p101_t1_g53016;
  double cse_p1_t10_g52997;
  double cse_p1_t11_g52998;
  double cse_p1_t12_g52999;
  double cse_p1_t13_g53000;
  double cse_p1_t14_g53001;
  double cse_p1_t15_g53002;
  double cse_p1_t16_g53003;
  double cse_p1_t1_g52988;
  double cse_p1_t1_g53004;
  double cse_p1_t1_g53005;
  double cse_p1_t1_g53006;
  double cse_p1_t1_g53007;
  double cse_p1_t1_g53008;
  double cse_p1_t1_g53010;
  double cse_p1_t1_g53011;
  double cse_p1_t1_g53012;
  double cse_p1_t1_g53013;
  double cse_p1_t2_g52989;
  double cse_p1_t3_g52990;
  double cse_p1_t4_g52991;
  double cse_p1_t5_g52992;
  double cse_p1_t6_g52993;
  double cse_p1_t7_g52994;
  double cse_p1_t8_g52995;
  double cse_p1_t9_g52996;
  double cse_p2_t1_g53014;
  double cse_p2_t1_g53015;
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
  double r13;
  double r13_2;
  double r23;
  double r23_2;
  double r43;
  double r43_2;
  double v;
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
    cse_p1_t8_g52995 = (-(x3));
    cse_p1_t9_g52996 = (-(y3));
    cse_p1_t10_g52997 = (-(z3));
    dx13 = (cse_p1_t8_g52995 + x1);
    cse_p1_t2_g52989 = (-(dx13));
    dy13 = (cse_p1_t9_g52996 + y1);
    cse_p1_t4_g52991 = (-(dy13));
    dz13 = (cse_p1_t10_g52997 + z1);
    cse_p1_t6_g52993 = (-(dz13));
    dx23 = (cse_p1_t8_g52995 + x2);
    cse_p1_t3_g52990 = (-(dx23));
    dy23 = (cse_p1_t9_g52996 + y2);
    cse_p1_t5_g52992 = (-(dy23));
    dz23 = (cse_p1_t10_g52997 + z2);
    cse_p1_t7_g52994 = (-(dz23));
    dx43 = (cse_p1_t8_g52995 + x4);
    dy43 = (cse_p1_t9_g52996 + y4);
    dz43 = (cse_p1_t10_g52997 + z4);
    r13_2 = ((dx13 * dx13) + (dy13 * dy13) + (dz13 * dz13));
    cse_p1_t14_g53001 = pow(r13_2, -0.50000000000000000    );
    r23_2 = ((dx23 * dx23) + (dy23 * dy23) + (dz23 * dz23));
    cse_p1_t15_g53002 = pow(r23_2, -0.50000000000000000    );
    r43_2 = ((dx43 * dx43) + (dy43 * dy43) + (dz43 * dz43));
    cse_p1_t16_g53003 = pow(r43_2, -0.50000000000000000    );
    r13 = sqrt(r13_2);
    r23 = sqrt(r23_2);
    r43 = sqrt(r43_2);
    cx = ((dy13 * dz23) + (-((dy23 * dz13))));
    cy = ((dx23 * dz13) + (-((dx13 * dz23))));
    cz = ((dx13 * dy23) + (-((dx23 * dy13))));
    v = ((cx * dx43) + (cy * dy43) + (cz * dz43));
    cse_p2_t1_g53015 = (r23 * r43);
    denom = (cse_p2_t1_g53015 * r13);
    cse_p1_t12_g52999 = pow(denom, -1);
    cse_p1_t13_g53000 = pow(denom, -2);
    cse_p101_t1_g53016 = (cse_p1_t13_g53000 * v);
    q = (cse_p1_t12_g52999 * v);
    cse_p1_t1_g52988 = (co + q);
    cse_p1_t11_g52998 = ((cse_p1_t1_g52988) * (cse_p1_t1_g52988));
    energy = (k * (cse_p1_t1_g52988 * cse_p1_t1_g52988 * cse_p1_t1_g52988));
    *energy_accumulate += energy;
    cse_p1_t1_g53004 = (3.0000000000000000     * cse_p1_t11_g52998 * k);
    cse_p1_t1_g53006 = (cse_p1_t14_g53001 * cse_p2_t1_g53015);
    cse_p1_t1_g53012 = (cse_p101_t1_g53016 * cse_p1_t1_g53006);
    g_x1 = (cse_p1_t1_g53004 * ((cse_p1_t12_g52999 * ((cse_p1_t7_g52994 * dy43) + (dy23 * dz43))) + (-((cse_p1_t1_g53012 * dx13)))));
    KernelGradientAcc(i3x1, 0, g_x1);
    g_y1 = (cse_p1_t1_g53004 * ((cse_p1_t12_g52999 * ((cse_p1_t3_g52990 * dz43) + (dx43 * dz23))) + (-((cse_p1_t1_g53012 * dy13)))));
    KernelGradientAcc(i3x1, 1, g_y1);
    g_z1 = (cse_p1_t1_g53004 * ((cse_p1_t12_g52999 * ((cse_p1_t5_g52992 * dx43) + (dx23 * dy43))) + (-((cse_p1_t1_g53012 * dz13)))));
    KernelGradientAcc(i3x1, 2, g_z1);
    cse_p1_t1_g53005 = (cse_p101_t1_g53016 * r13);
    cse_p1_t1_g53007 = (cse_p1_t15_g53002 * r43);
    cse_p1_t1_g53010 = (cse_p1_t1_g53005 * cse_p1_t1_g53007);
    g_x2 = (cse_p1_t1_g53004 * ((cse_p1_t12_g52999 * ((cse_p1_t4_g52991 * dz43) + (dy43 * dz13))) + (-((cse_p1_t1_g53010 * dx23)))));
    KernelGradientAcc(i3x2, 0, g_x2);
    g_y2 = (cse_p1_t1_g53004 * ((cse_p1_t12_g52999 * ((cse_p1_t6_g52993 * dx43) + (dx13 * dz43))) + (-((cse_p1_t1_g53010 * dy23)))));
    KernelGradientAcc(i3x2, 1, g_y2);
    g_z2 = (cse_p1_t1_g53004 * ((cse_p1_t12_g52999 * ((cse_p1_t2_g52989 * dy43) + (dx43 * dy13))) + (-((cse_p1_t1_g53010 * dz23)))));
    KernelGradientAcc(i3x2, 2, g_z2);
    cse_p1_t1_g53008 = (cse_p1_t16_g53003 * r23);
    cse_p1_t1_g53013 = (cse_p1_t1_g53007 * r13);
    cse_p2_t1_g53014 = (cse_p1_t1_g53008 * r13);
    g_x3 = (cse_p1_t1_g53004 * ((cse_p1_t12_g52999 * ((dy43 * (cse_p1_t6_g52993 + dz23)) + (dz43 * (cse_p1_t5_g52992 + dy13)) + (-(cx)))) + (-((cse_p101_t1_g53016 * ((-((cse_p1_t1_g53006 * dx13))) + (-((cse_p1_t1_g53013 * dx23))) + (-((cse_p2_t1_g53014 * dx43)))))))));
    KernelGradientAcc(i3x3, 0, g_x3);
    g_y3 = (cse_p1_t1_g53004 * ((cse_p1_t12_g52999 * ((dx43 * (cse_p1_t7_g52994 + dz13)) + (dz43 * (cse_p1_t2_g52989 + dx23)) + (-(cy)))) + (-((cse_p101_t1_g53016 * ((-((cse_p1_t1_g53006 * dy13))) + (-((cse_p1_t1_g53013 * dy23))) + (-((cse_p2_t1_g53014 * dy43)))))))));
    KernelGradientAcc(i3x3, 1, g_y3);
    g_z3 = (cse_p1_t1_g53004 * ((cse_p1_t12_g52999 * ((dx43 * (cse_p1_t4_g52991 + dy23)) + (dy43 * (cse_p1_t3_g52990 + dx13)) + (-(cz)))) + (-((cse_p101_t1_g53016 * ((-((cse_p1_t1_g53006 * dz13))) + (-((cse_p1_t1_g53013 * dz23))) + (-((cse_p2_t1_g53014 * dz43)))))))));
    KernelGradientAcc(i3x3, 2, g_z3);
    cse_p1_t1_g53011 = (cse_p1_t1_g53005 * cse_p1_t1_g53008);
    g_x4 = (cse_p1_t1_g53004 * ((cse_p1_t12_g52999 * cx) + (-((cse_p1_t1_g53011 * dx43)))));
    KernelGradientAcc(i3x4, 0, g_x4);
    g_y4 = (cse_p1_t1_g53004 * ((cse_p1_t12_g52999 * cy) + (-((cse_p1_t1_g53011 * dy43)))));
    KernelGradientAcc(i3x4, 1, g_y4);
    g_z4 = (cse_p1_t1_g53004 * ((cse_p1_t12_g52999 * cz) + (-((cse_p1_t1_g53011 * dz43)))));
    KernelGradientAcc(i3x4, 2, g_z4);
  }
}
void chiral_restraint_hessian_fd(double k, double co, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 12;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  const double invh2 = 1.0/((h*h));
  double e0 = 0.0;
  chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
  {
    double saved = position[i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 0] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 0, d);
  }
  {
    double saved = position[i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 1] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 1, d);
  }
  {
    double saved = position[i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 2] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 2, d);
  }
  {
    double saved = position[i3x4 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 0] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 0, d);
  }
  {
    double saved = position[i3x4 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 1] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 1, d);
  }
  {
    double saved = position[i3x4 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 2] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 2, d);
  }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 0, i3x1, 0, hval);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 1, i3x1, 1, hval);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 2, i3x1, 2, hval);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 0, i3x2, 0, hval);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 1, i3x2, 1, hval);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x2, 2, hval);
  }
  {
    double saved = position[i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 0] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 0, i3x3, 0, hval);
  }
  {
    double saved = position[i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 1] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x3, 1, hval);
  }
  {
    double saved = position[i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 2] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x3, 2, hval);
  }
  {
    double saved = position[i3x4 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 0] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 0, i3x4, 0, hval);
  }
  {
    double saved = position[i3x4 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 1] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 1, i3x4, 1, hval);
  }
  {
    double saved = position[i3x4 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 2] = saved + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x3 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x3 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x3 + 2] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x3 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x3 + 2] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x4 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x3 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x3 + 2] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x4 + 0] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x4 + 1] = saved_j + h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_j - h;
    chiral_restraint_energy(k, co, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i; position[i3x4 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 2, i3x4, 1, hval);
  }
}

