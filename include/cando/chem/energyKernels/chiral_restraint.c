template <typename HESSIAN>
struct Chiral_Restraint {
  static constexpr size_t PositionSize = 12;
  static std::string description() { return "mathkernel-chiral_restraint"; };
void energy(double k, double co, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  /* !BASE */
  double cse_p11_t1_g16935;
  double cse_p11_t2_g16936;
  double cse_p11_t3_g16937;
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
    cse_p11_t1_g16935 = (-(x3));
    cse_p11_t2_g16936 = (-(y3));
    cse_p11_t3_g16937 = (-(z3));
    dx13 = (cse_p11_t1_g16935 + x1);
    dy13 = (cse_p11_t2_g16936 + y1);
    dz13 = (cse_p11_t3_g16937 + z1);
    dx23 = (cse_p11_t1_g16935 + x2);
    dy23 = (cse_p11_t2_g16936 + y2);
    dz23 = (cse_p11_t3_g16937 + z2);
    dx43 = (cse_p11_t1_g16935 + x4);
    dy43 = (cse_p11_t2_g16936 + y4);
    dz43 = (cse_p11_t3_g16937 + z4);
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
    *energy_accumulate += energy;
  }
}
void energy_fd(double k, double co, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  energy(k, co, i3x1, i3x2, i3x3, i3x4, position, energy_accumulate, force, hessian, dvec, hdvec);
}

void gradient(double k, double co, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  /* !BASE */
  double cse_p113_t1_g16966;
  double cse_p11_t10_g16955;
  double cse_p11_t11_g16956;
  double cse_p11_t1_g16946;
  double cse_p11_t1_g16957;
  double cse_p11_t1_g16958;
  double cse_p11_t1_g16960;
  double cse_p11_t1_g16961;
  double cse_p11_t1_g16962;
  double cse_p11_t1_g16963;
  double cse_p11_t1_g16964;
  double cse_p11_t1_g16965;
  double cse_p11_t2_g16947;
  double cse_p11_t3_g16948;
  double cse_p11_t4_g16949;
  double cse_p11_t5_g16950;
  double cse_p11_t6_g16951;
  double cse_p11_t7_g16952;
  double cse_p11_t8_g16953;
  double cse_p11_t9_g16954;
  double cse_p1_t1_sqrt16938;
  double cse_p1_t2_invsqrt16939;
  double cse_p1_t3_sqrt16940;
  double cse_p1_t4_invsqrt16941;
  double cse_p1_t5_sqrt16942;
  double cse_p1_t6_invsqrt16943;
  double cse_p1_t7_invr16944;
  double cse_p1_t8_invr216945;
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
  double r13_2;
  double r23_2;
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
    cse_p11_t8_g16953 = (-(x3));
    cse_p11_t9_g16954 = (-(y3));
    cse_p11_t10_g16955 = (-(z3));
    dx13 = (cse_p11_t8_g16953 + x1);
    cse_p11_t2_g16947 = (-(dx13));
    dy13 = (cse_p11_t9_g16954 + y1);
    cse_p11_t4_g16949 = (-(dy13));
    dz13 = (cse_p11_t10_g16955 + z1);
    cse_p11_t6_g16951 = (-(dz13));
    dx23 = (cse_p11_t8_g16953 + x2);
    cse_p11_t3_g16948 = (-(dx23));
    dy23 = (cse_p11_t9_g16954 + y2);
    cse_p11_t5_g16950 = (-(dy23));
    dz23 = (cse_p11_t10_g16955 + z2);
    cse_p11_t7_g16952 = (-(dz23));
    dx43 = (cse_p11_t8_g16953 + x4);
    dy43 = (cse_p11_t9_g16954 + y4);
    dz43 = (cse_p11_t10_g16955 + z4);
    r13_2 = ((dx13 * dx13) + (dy13 * dy13) + (dz13 * dz13));
    r23_2 = ((dx23 * dx23) + (dy23 * dy23) + (dz23 * dz23));
    r43_2 = ((dx43 * dx43) + (dy43 * dy43) + (dz43 * dz43));
    cse_p1_t1_sqrt16938 = sqrt(r13_2);
    cse_p1_t2_invsqrt16939 = (1.0 / (cse_p1_t1_sqrt16938));
    cse_p1_t3_sqrt16940 = sqrt(r23_2);
    cse_p1_t4_invsqrt16941 = (1.0 / (cse_p1_t3_sqrt16940));
    cse_p1_t5_sqrt16942 = sqrt(r43_2);
    cse_p1_t6_invsqrt16943 = (1.0 / (cse_p1_t5_sqrt16942));
    cx = ((dy13 * dz23) + (-((dy23 * dz13))));
    cy = ((dx23 * dz13) + (-((dx13 * dz23))));
    cz = ((dx13 * dy23) + (-((dx23 * dy13))));
    v = ((cx * dx43) + (cy * dy43) + (cz * dz43));
    cse_p11_t1_g16965 = (cse_p1_t1_sqrt16938 * cse_p1_t3_sqrt16940);
    denom = (cse_p11_t1_g16965 * cse_p1_t5_sqrt16942);
    cse_p1_t7_invr16944 = (1.0 / (denom));
    cse_p1_t8_invr216945 = (cse_p1_t7_invr16944 * cse_p1_t7_invr16944);
    cse_p113_t1_g16966 = (cse_p1_t8_invr216945 * v);
    q = (cse_p1_t7_invr16944 * v);
    cse_p11_t1_g16946 = (co + q);
    cse_p11_t11_g16956 = ((cse_p11_t1_g16946) * (cse_p11_t1_g16946));
    energy = (k * (cse_p11_t1_g16946 * cse_p11_t1_g16946 * cse_p11_t1_g16946));
    *energy_accumulate += energy;
    cse_p11_t1_g16957 = (3.0000000000000000     * cse_p11_t11_g16956 * k);
    cse_p11_t1_g16958 = (cse_p1_t2_invsqrt16939 * cse_p1_t3_sqrt16940 * cse_p1_t5_sqrt16942);
    cse_p11_t1_g16962 = (cse_p113_t1_g16966 * cse_p11_t1_g16958);
    g_x1 = (cse_p11_t1_g16957 * ((cse_p1_t7_invr16944 * ((cse_p11_t7_g16952 * dy43) + (dy23 * dz43))) + (-((cse_p11_t1_g16962 * dx13)))));
    KernelGradientAcc(i3x1, 0, g_x1);
    g_y1 = (cse_p11_t1_g16957 * ((cse_p1_t7_invr16944 * ((cse_p11_t3_g16948 * dz43) + (dx43 * dz23))) + (-((cse_p11_t1_g16962 * dy13)))));
    KernelGradientAcc(i3x1, 1, g_y1);
    g_z1 = (cse_p11_t1_g16957 * ((cse_p1_t7_invr16944 * ((cse_p11_t5_g16950 * dx43) + (dx23 * dy43))) + (-((cse_p11_t1_g16962 * dz13)))));
    KernelGradientAcc(i3x1, 2, g_z1);
    cse_p11_t1_g16960 = (cse_p1_t1_sqrt16938 * cse_p1_t4_invsqrt16941 * cse_p1_t5_sqrt16942);
    cse_p11_t1_g16963 = (cse_p113_t1_g16966 * cse_p11_t1_g16960);
    g_x2 = (cse_p11_t1_g16957 * ((cse_p1_t7_invr16944 * ((cse_p11_t4_g16949 * dz43) + (dy43 * dz13))) + (-((cse_p11_t1_g16963 * dx23)))));
    KernelGradientAcc(i3x2, 0, g_x2);
    g_y2 = (cse_p11_t1_g16957 * ((cse_p1_t7_invr16944 * ((cse_p11_t6_g16951 * dx43) + (dx13 * dz43))) + (-((cse_p11_t1_g16963 * dy23)))));
    KernelGradientAcc(i3x2, 1, g_y2);
    g_z2 = (cse_p11_t1_g16957 * ((cse_p1_t7_invr16944 * ((cse_p11_t2_g16947 * dy43) + (dx43 * dy13))) + (-((cse_p11_t1_g16963 * dz23)))));
    KernelGradientAcc(i3x2, 2, g_z2);
    cse_p11_t1_g16961 = (cse_p11_t1_g16965 * cse_p1_t6_invsqrt16943);
    g_x3 = (cse_p11_t1_g16957 * ((cse_p1_t7_invr16944 * ((dy43 * (cse_p11_t6_g16951 + dz23)) + (dz43 * (cse_p11_t5_g16950 + dy13)) + (-(cx)))) + (-((cse_p113_t1_g16966 * ((-((cse_p11_t1_g16958 * dx13))) + (-((cse_p11_t1_g16960 * dx23))) + (-((cse_p11_t1_g16961 * dx43)))))))));
    KernelGradientAcc(i3x3, 0, g_x3);
    g_y3 = (cse_p11_t1_g16957 * ((cse_p1_t7_invr16944 * ((dx43 * (cse_p11_t7_g16952 + dz13)) + (dz43 * (cse_p11_t2_g16947 + dx23)) + (-(cy)))) + (-((cse_p113_t1_g16966 * ((-((cse_p11_t1_g16958 * dy13))) + (-((cse_p11_t1_g16960 * dy23))) + (-((cse_p11_t1_g16961 * dy43)))))))));
    KernelGradientAcc(i3x3, 1, g_y3);
    g_z3 = (cse_p11_t1_g16957 * ((cse_p1_t7_invr16944 * ((dx43 * (cse_p11_t4_g16949 + dy23)) + (dy43 * (cse_p11_t3_g16948 + dx13)) + (-(cz)))) + (-((cse_p113_t1_g16966 * ((-((cse_p11_t1_g16958 * dz13))) + (-((cse_p11_t1_g16960 * dz23))) + (-((cse_p11_t1_g16961 * dz43)))))))));
    KernelGradientAcc(i3x3, 2, g_z3);
    cse_p11_t1_g16964 = (cse_p113_t1_g16966 * cse_p11_t1_g16961);
    g_x4 = (cse_p11_t1_g16957 * ((cse_p1_t7_invr16944 * cx) + (-((cse_p11_t1_g16964 * dx43)))));
    KernelGradientAcc(i3x4, 0, g_x4);
    g_y4 = (cse_p11_t1_g16957 * ((cse_p1_t7_invr16944 * cy) + (-((cse_p11_t1_g16964 * dy43)))));
    KernelGradientAcc(i3x4, 1, g_y4);
    g_z4 = (cse_p11_t1_g16957 * ((cse_p1_t7_invr16944 * cz) + (-((cse_p11_t1_g16964 * dz43)))));
    KernelGradientAcc(i3x4, 2, g_z4);
  }
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

void hessian(double k, double co, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  /* !BASE */
  double cse_p113_t1_g16995;
  double cse_p11_t10_g16984;
  double cse_p11_t11_g16985;
  double cse_p11_t1_g16975;
  double cse_p11_t1_g16986;
  double cse_p11_t1_g16987;
  double cse_p11_t1_g16989;
  double cse_p11_t1_g16990;
  double cse_p11_t1_g16991;
  double cse_p11_t1_g16992;
  double cse_p11_t1_g16993;
  double cse_p11_t1_g16994;
  double cse_p11_t2_g16976;
  double cse_p11_t3_g16977;
  double cse_p11_t4_g16978;
  double cse_p11_t5_g16979;
  double cse_p11_t6_g16980;
  double cse_p11_t7_g16981;
  double cse_p11_t8_g16982;
  double cse_p11_t9_g16983;
  double cse_p1_t1_sqrt16967;
  double cse_p1_t2_invsqrt16968;
  double cse_p1_t3_sqrt16969;
  double cse_p1_t4_invsqrt16970;
  double cse_p1_t5_sqrt16971;
  double cse_p1_t6_invsqrt16972;
  double cse_p1_t7_invr16973;
  double cse_p1_t8_invr216974;
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
  double r13_2;
  double r23_2;
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
    cse_p11_t8_g16982 = (-(x3));
    cse_p11_t9_g16983 = (-(y3));
    cse_p11_t10_g16984 = (-(z3));
    dx13 = (cse_p11_t8_g16982 + x1);
    cse_p11_t2_g16976 = (-(dx13));
    dy13 = (cse_p11_t9_g16983 + y1);
    cse_p11_t4_g16978 = (-(dy13));
    dz13 = (cse_p11_t10_g16984 + z1);
    cse_p11_t6_g16980 = (-(dz13));
    dx23 = (cse_p11_t8_g16982 + x2);
    cse_p11_t3_g16977 = (-(dx23));
    dy23 = (cse_p11_t9_g16983 + y2);
    cse_p11_t5_g16979 = (-(dy23));
    dz23 = (cse_p11_t10_g16984 + z2);
    cse_p11_t7_g16981 = (-(dz23));
    dx43 = (cse_p11_t8_g16982 + x4);
    dy43 = (cse_p11_t9_g16983 + y4);
    dz43 = (cse_p11_t10_g16984 + z4);
    r13_2 = ((dx13 * dx13) + (dy13 * dy13) + (dz13 * dz13));
    r23_2 = ((dx23 * dx23) + (dy23 * dy23) + (dz23 * dz23));
    r43_2 = ((dx43 * dx43) + (dy43 * dy43) + (dz43 * dz43));
    cse_p1_t1_sqrt16967 = sqrt(r13_2);
    cse_p1_t2_invsqrt16968 = (1.0 / (cse_p1_t1_sqrt16967));
    cse_p1_t3_sqrt16969 = sqrt(r23_2);
    cse_p1_t4_invsqrt16970 = (1.0 / (cse_p1_t3_sqrt16969));
    cse_p1_t5_sqrt16971 = sqrt(r43_2);
    cse_p1_t6_invsqrt16972 = (1.0 / (cse_p1_t5_sqrt16971));
    cx = ((dy13 * dz23) + (-((dy23 * dz13))));
    cy = ((dx23 * dz13) + (-((dx13 * dz23))));
    cz = ((dx13 * dy23) + (-((dx23 * dy13))));
    v = ((cx * dx43) + (cy * dy43) + (cz * dz43));
    cse_p11_t1_g16994 = (cse_p1_t1_sqrt16967 * cse_p1_t3_sqrt16969);
    denom = (cse_p11_t1_g16994 * cse_p1_t5_sqrt16971);
    cse_p1_t7_invr16973 = (1.0 / (denom));
    cse_p1_t8_invr216974 = (cse_p1_t7_invr16973 * cse_p1_t7_invr16973);
    cse_p113_t1_g16995 = (cse_p1_t8_invr216974 * v);
    q = (cse_p1_t7_invr16973 * v);
    cse_p11_t1_g16975 = (co + q);
    cse_p11_t11_g16985 = ((cse_p11_t1_g16975) * (cse_p11_t1_g16975));
    energy = (k * (cse_p11_t1_g16975 * cse_p11_t1_g16975 * cse_p11_t1_g16975));
    *energy_accumulate += energy;
    cse_p11_t1_g16986 = (3.0000000000000000     * cse_p11_t11_g16985 * k);
    cse_p11_t1_g16987 = (cse_p1_t2_invsqrt16968 * cse_p1_t3_sqrt16969 * cse_p1_t5_sqrt16971);
    cse_p11_t1_g16991 = (cse_p113_t1_g16995 * cse_p11_t1_g16987);
    g_x1 = (cse_p11_t1_g16986 * ((cse_p1_t7_invr16973 * ((cse_p11_t7_g16981 * dy43) + (dy23 * dz43))) + (-((cse_p11_t1_g16991 * dx13)))));
    KernelGradientAcc(i3x1, 0, g_x1);
    g_y1 = (cse_p11_t1_g16986 * ((cse_p1_t7_invr16973 * ((cse_p11_t3_g16977 * dz43) + (dx43 * dz23))) + (-((cse_p11_t1_g16991 * dy13)))));
    KernelGradientAcc(i3x1, 1, g_y1);
    g_z1 = (cse_p11_t1_g16986 * ((cse_p1_t7_invr16973 * ((cse_p11_t5_g16979 * dx43) + (dx23 * dy43))) + (-((cse_p11_t1_g16991 * dz13)))));
    KernelGradientAcc(i3x1, 2, g_z1);
    cse_p11_t1_g16989 = (cse_p1_t1_sqrt16967 * cse_p1_t4_invsqrt16970 * cse_p1_t5_sqrt16971);
    cse_p11_t1_g16992 = (cse_p113_t1_g16995 * cse_p11_t1_g16989);
    g_x2 = (cse_p11_t1_g16986 * ((cse_p1_t7_invr16973 * ((cse_p11_t4_g16978 * dz43) + (dy43 * dz13))) + (-((cse_p11_t1_g16992 * dx23)))));
    KernelGradientAcc(i3x2, 0, g_x2);
    g_y2 = (cse_p11_t1_g16986 * ((cse_p1_t7_invr16973 * ((cse_p11_t6_g16980 * dx43) + (dx13 * dz43))) + (-((cse_p11_t1_g16992 * dy23)))));
    KernelGradientAcc(i3x2, 1, g_y2);
    g_z2 = (cse_p11_t1_g16986 * ((cse_p1_t7_invr16973 * ((cse_p11_t2_g16976 * dy43) + (dx43 * dy13))) + (-((cse_p11_t1_g16992 * dz23)))));
    KernelGradientAcc(i3x2, 2, g_z2);
    cse_p11_t1_g16990 = (cse_p11_t1_g16994 * cse_p1_t6_invsqrt16972);
    g_x3 = (cse_p11_t1_g16986 * ((cse_p1_t7_invr16973 * ((dy43 * (cse_p11_t6_g16980 + dz23)) + (dz43 * (cse_p11_t5_g16979 + dy13)) + (-(cx)))) + (-((cse_p113_t1_g16995 * ((-((cse_p11_t1_g16987 * dx13))) + (-((cse_p11_t1_g16989 * dx23))) + (-((cse_p11_t1_g16990 * dx43)))))))));
    KernelGradientAcc(i3x3, 0, g_x3);
    g_y3 = (cse_p11_t1_g16986 * ((cse_p1_t7_invr16973 * ((dx43 * (cse_p11_t7_g16981 + dz13)) + (dz43 * (cse_p11_t2_g16976 + dx23)) + (-(cy)))) + (-((cse_p113_t1_g16995 * ((-((cse_p11_t1_g16987 * dy13))) + (-((cse_p11_t1_g16989 * dy23))) + (-((cse_p11_t1_g16990 * dy43)))))))));
    KernelGradientAcc(i3x3, 1, g_y3);
    g_z3 = (cse_p11_t1_g16986 * ((cse_p1_t7_invr16973 * ((dx43 * (cse_p11_t4_g16978 + dy23)) + (dy43 * (cse_p11_t3_g16977 + dx13)) + (-(cz)))) + (-((cse_p113_t1_g16995 * ((-((cse_p11_t1_g16987 * dz13))) + (-((cse_p11_t1_g16989 * dz23))) + (-((cse_p11_t1_g16990 * dz43)))))))));
    KernelGradientAcc(i3x3, 2, g_z3);
    cse_p11_t1_g16993 = (cse_p113_t1_g16995 * cse_p11_t1_g16990);
    g_x4 = (cse_p11_t1_g16986 * ((cse_p1_t7_invr16973 * cx) + (-((cse_p11_t1_g16993 * dx43)))));
    KernelGradientAcc(i3x4, 0, g_x4);
    g_y4 = (cse_p11_t1_g16986 * ((cse_p1_t7_invr16973 * cy) + (-((cse_p11_t1_g16993 * dy43)))));
    KernelGradientAcc(i3x4, 1, g_y4);
    g_z4 = (cse_p11_t1_g16986 * ((cse_p1_t7_invr16973 * cz) + (-((cse_p11_t1_g16993 * dz43)))));
    KernelGradientAcc(i3x4, 2, g_z4);
  }
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
