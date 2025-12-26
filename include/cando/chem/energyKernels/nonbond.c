template <typename HESSIAN>
struct Nonbond {
  static constexpr size_t PositionSize = 6;
  static std::string description() { return "mathkernel-nonbond"; };
void gradient(double a, double b, double qq, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  /* !BASE */
  double cse_p113_t1_g1088;
  double cse_p11_t1_g1075;
  double cse_p11_t1_g1078;
  double cse_p11_t1_g1079;
  double cse_p11_t1_g1080;
  double cse_p11_t1_g1081;
  double cse_p11_t1_g1082;
  double cse_p11_t1_g1084;
  double cse_p11_t1_g1085;
  double cse_p11_t2_g1076;
  double cse_p11_t3_g1077;
  double cse_p12_t1_g1086;
  double cse_p1_t1_sqrt1071;
  double cse_p61_t1_invr1087;
  double dx;
  double dy;
  double dz;
  double energy;
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
  DOUBLE x1 = position[i3x1 + 0];
  DOUBLE y1 = position[i3x1 + 1];
  DOUBLE z1 = position[i3x1 + 2];
  DOUBLE x2 = position[i3x2 + 0];
  DOUBLE y2 = position[i3x2 + 1];
  DOUBLE z2 = position[i3x2 + 2];
  {
    /* !BASE */
    cse_p11_t1_g1075 = (-3.0000000000000000     * b);
    dx = (x1 + (-(x2)));
    dy = (y1 + (-(y2)));
    dz = (z1 + (-(z2)));
    r2 = ((dx * dx) + (dy * dy) + (dz * dz));
    cse_p1_t1_sqrt1071 = sqrt(r2);
    cse_p61_t1_invr1087 = (1.0 / (cse_p1_t1_sqrt1071));
    cse_p113_t1_g1088 = (cse_p61_t1_invr1087 * cse_p61_t1_invr1087);
    cse_p12_t1_g1086 = (cse_p113_t1_g1088 * cse_p113_t1_g1088);
    invr6 = (cse_p113_t1_g1088 * cse_p12_t1_g1086);
    cse_p11_t1_g1085 = (a * invr6);
    cse_p11_t2_g1076 = (6.0000000000000000     * cse_p11_t1_g1085);
    cse_p11_t3_g1077 = (cse_p11_t1_g1075 + cse_p11_t2_g1076);
    e_lj = ((cse_p11_t1_g1085 * invr6) + (-((b * invr6))));
    e_coul = (cse_p61_t1_invr1087 * qq);
    energy = (e_coul + e_lj);
    *energy_accumulate += energy;
    cse_p11_t1_g1078 = (cse_p113_t1_g1088 * cse_p61_t1_invr1087);
    cse_p11_t1_g1079 = (cse_p11_t1_g1078 * cse_p11_t3_g1077 * cse_p12_t1_g1086 * cse_p61_t1_invr1087);
    cse_p11_t1_g1080 = (cse_p11_t1_g1078 * qq);
    cse_p11_t1_g1081 = (-(cse_p11_t1_g1080));
    cse_p11_t1_g1082 = (-2.0000000000000000     * cse_p11_t1_g1079);
    g_x1 = ((cse_p11_t1_g1081 * dx) + (cse_p11_t1_g1082 * dx));
    KernelGradientAcc(i3x1, 0, g_x1);
    g_y1 = ((cse_p11_t1_g1081 * dy) + (cse_p11_t1_g1082 * dy));
    KernelGradientAcc(i3x1, 1, g_y1);
    g_z1 = ((cse_p11_t1_g1081 * dz) + (cse_p11_t1_g1082 * dz));
    KernelGradientAcc(i3x1, 2, g_z1);
    cse_p11_t1_g1084 = (2.0000000000000000     * cse_p11_t1_g1079);
    g_x2 = ((cse_p11_t1_g1080 * dx) + (cse_p11_t1_g1084 * dx));
    KernelGradientAcc(i3x2, 0, g_x2);
    g_y2 = ((cse_p11_t1_g1080 * dy) + (cse_p11_t1_g1084 * dy));
    KernelGradientAcc(i3x2, 1, g_y2);
    g_z2 = ((cse_p11_t1_g1080 * dz) + (cse_p11_t1_g1084 * dz));
    KernelGradientAcc(i3x2, 2, g_z2);
  }
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

void energy(double a, double b, double qq, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  /* !BASE */
  double dx;
  double dy;
  double dz;
  double energy;
  double e_coul;
  double e_lj;
  double invr;
  double invr2;
  double invr6;
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
    *energy_accumulate += energy;
  }
}
void energy_fd(double a, double b, double qq, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  energy(a, b, qq, i3x1, i3x2, position, energy_accumulate, force, hessian, dvec, hdvec);
}

void hessian(double a, double b, double qq, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  /* !BASE */
  double cse_p114_t1_g1165;
  double cse_p11_t11_g1103;
  double cse_p11_t12_g1104;
  double cse_p11_t13_g1105;
  double cse_p11_t14_g1106;
  double cse_p11_t15_g1107;
  double cse_p11_t16_g1108;
  double cse_p11_t17_g1109;
  double cse_p11_t18_g1110;
  double cse_p11_t19_g1111;
  double cse_p11_t1_g1093;
  double cse_p11_t1_g1140;
  double cse_p11_t1_g1141;
  double cse_p11_t1_g1142;
  double cse_p11_t1_g1143;
  double cse_p11_t1_g1144;
  double cse_p11_t22_g1114;
  double cse_p11_t23_g1115;
  double cse_p11_t24_g1116;
  double cse_p11_t25_g1117;
  double cse_p11_t26_g1118;
  double cse_p11_t27_g1119;
  double cse_p11_t28_g1120;
  double cse_p11_t29_g1121;
  double cse_p11_t2_g1094;
  double cse_p11_t30_g1122;
  double cse_p11_t31_g1123;
  double cse_p11_t32_g1124;
  double cse_p11_t33_g1125;
  double cse_p11_t34_g1126;
  double cse_p11_t35_g1127;
  double cse_p11_t3_g1095;
  double cse_p11_t42_g1134;
  double cse_p11_t43_g1135;
  double cse_p11_t44_g1136;
  double cse_p11_t45_g1137;
  double cse_p11_t46_g1138;
  double cse_p11_t47_g1139;
  double cse_p11_t5_g1097;
  double cse_p11_t6_g1098;
  double cse_p11_t7_g1099;
  double cse_p11_t8_g1100;
  double cse_p11_t9_g1101;
  double cse_p166_t1_g1166;
  double cse_p1_t1_sqrt1089;
  double cse_p218_t1_g1167;
  double cse_p61_t1_invr1145;
  double cse_p61_t2_invr1146;
  double cse_p61_t3_invr1147;
  double cse_p61_t4_invr1148;
  double cse_p61_t5_invr1149;
  double cse_p61_t6_invr1150;
  double cse_p61_t7_invr1151;
  double cse_p62_t1_g1152;
  double cse_p63_t1_g1153;
  double cse_p63_t1_g1156;
  double cse_p63_t1_g1157;
  double cse_p63_t1_g1158;
  double cse_p63_t1_g1159;
  double cse_p63_t1_g1160;
  double cse_p63_t1_g1161;
  double cse_p63_t1_g1162;
  double cse_p63_t1_g1163;
  double cse_p63_t1_g1164;
  double cse_p63_t2_g1154;
  double cse_p63_t3_g1155;
  double dx;
  double dy;
  double dz;
  double energy;
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
  DOUBLE x1 = position[i3x1 + 0];
  DOUBLE y1 = position[i3x1 + 1];
  DOUBLE z1 = position[i3x1 + 2];
  DOUBLE x2 = position[i3x2 + 0];
  DOUBLE y2 = position[i3x2 + 1];
  DOUBLE z2 = position[i3x2 + 2];
  {
    /* !BASE */
    cse_p11_t5_g1097 = (-3.0000000000000000     * b);
    dx = (x1 + (-(x2)));
    cse_p11_t11_g1103 = (dx * dx);
    cse_p11_t28_g1120 = (-(cse_p11_t11_g1103));
    cse_p11_t33_g1125 = (-(dx));
    dy = (y1 + (-(y2)));
    cse_p11_t15_g1107 = (dy * dy);
    cse_p11_t29_g1121 = (-(cse_p11_t15_g1107));
    cse_p11_t34_g1126 = (-(dy));
    dz = (z1 + (-(z2)));
    cse_p11_t19_g1111 = (dz * dz);
    cse_p11_t30_g1122 = (-(cse_p11_t19_g1111));
    cse_p11_t35_g1127 = (-(dz));
    r2 = (cse_p11_t11_g1103 + cse_p11_t15_g1107 + cse_p11_t19_g1111);
    cse_p1_t1_sqrt1089 = sqrt(r2);
    cse_p62_t1_g1152 = (cse_p1_t1_sqrt1089 * cse_p1_t1_sqrt1089);
    cse_p61_t1_invr1145 = (1.0 / (cse_p1_t1_sqrt1089));
    cse_p114_t1_g1165 = (cse_p61_t1_invr1145 * cse_p61_t1_invr1145);
    cse_p61_t2_invr1146 = (1.0 / (((cse_p1_t1_sqrt1089) * (cse_p1_t1_sqrt1089))));
    cse_p11_t1_g1144 = (cse_p61_t2_invr1146 * qq);
    cse_p61_t5_invr1149 = (1.0 / ((cse_p1_t1_sqrt1089 * cse_p1_t1_sqrt1089 * cse_p1_t1_sqrt1089)));
    cse_p61_t4_invr1148 = (1.0 / (pow(cse_p1_t1_sqrt1089, 7)));
    cse_p218_t1_g1167 = (6.0000000000000000     * b * cse_p61_t4_invr1148);
    cse_p61_t3_invr1147 = (1.0 / (pow(cse_p1_t1_sqrt1089, 13)));
    cse_p166_t1_g1166 = (-12.000000000000000     * a * cse_p61_t3_invr1147);
    cse_p11_t1_g1140 = (cse_p61_t5_invr1149 * (cse_p166_t1_g1166 + cse_p218_t1_g1167 + (-(cse_p11_t1_g1144))));
    cse_p11_t1_g1093 = (cse_p11_t1_g1140 * (cse_p11_t28_g1120 + cse_p62_t1_g1152));
    cse_p11_t2_g1094 = (cse_p11_t1_g1140 * (cse_p11_t29_g1121 + cse_p62_t1_g1152));
    cse_p11_t3_g1095 = (cse_p11_t1_g1140 * (cse_p11_t30_g1122 + cse_p62_t1_g1152));
    cse_p61_t6_invr1150 = (1.0 / (pow(cse_p1_t1_sqrt1089, 8)));
    cse_p11_t6_g1098 = (-42.000000000000000     * b * cse_p61_t6_invr1150);
    cse_p61_t7_invr1151 = (1.0 / (pow(cse_p1_t1_sqrt1089, 14)));
    cse_p11_t7_g1099 = (156.00000000000000     * a * cse_p61_t7_invr1151);
    cse_p11_t8_g1100 = (2.0000000000000000     * cse_p61_t5_invr1149 * qq);
    cse_p11_t12_g1104 = (cse_p11_t1_g1140 * dx * dy);
    cse_p63_t1_g1158 = (cse_p11_t1_g1140 * dz);
    cse_p11_t13_g1105 = (cse_p63_t1_g1158 * dx);
    cse_p11_t14_g1106 = (cse_p11_t1_g1140 * cse_p11_t33_g1125 * dy);
    cse_p11_t16_g1108 = (cse_p63_t1_g1158 * dy);
    cse_p11_t17_g1109 = (cse_p11_t33_g1125 * cse_p63_t1_g1158);
    cse_p11_t18_g1110 = (cse_p11_t34_g1126 * cse_p63_t1_g1158);
    cse_p11_t22_g1114 = (cse_p166_t1_g1166 + cse_p218_t1_g1167 + (-(cse_p11_t1_g1144)));
    cse_p11_t24_g1116 = (cse_p11_t6_g1098 + cse_p11_t7_g1099 + cse_p11_t8_g1100);
    cse_p11_t25_g1117 = (cse_p11_t28_g1120 + cse_p62_t1_g1152);
    cse_p11_t26_g1118 = (cse_p11_t29_g1121 + cse_p62_t1_g1152);
    cse_p11_t27_g1119 = (cse_p11_t30_g1122 + cse_p62_t1_g1152);
    cse_p11_t31_g1123 = (-(cse_p11_t1_g1144));
    cse_p11_t32_g1124 = (-(cse_p62_t1_g1152));
    cse_p11_t42_g1134 = pow(cse_p1_t1_sqrt1089, 13);
    cse_p11_t43_g1135 = pow(cse_p1_t1_sqrt1089, 14);
    cse_p11_t44_g1136 = ((cse_p1_t1_sqrt1089) * (cse_p1_t1_sqrt1089));
    cse_p11_t45_g1137 = (cse_p1_t1_sqrt1089 * cse_p1_t1_sqrt1089 * cse_p1_t1_sqrt1089);
    cse_p11_t46_g1138 = pow(cse_p1_t1_sqrt1089, 7);
    cse_p11_t47_g1139 = pow(cse_p1_t1_sqrt1089, 8);
    invr6 = (cse_p114_t1_g1165 * cse_p114_t1_g1165 * cse_p114_t1_g1165);
    cse_p11_t9_g1101 = (6.0000000000000000     * a * invr6);
    cse_p11_t23_g1115 = (cse_p11_t5_g1097 + cse_p11_t9_g1101);
    e_lj = ((a * invr6 * invr6) + (-((b * invr6))));
    e_coul = (cse_p61_t1_invr1145 * qq);
    energy = (e_coul + e_lj);
    *energy_accumulate += energy;
    cse_p11_t1_g1142 = (cse_p114_t1_g1165 * cse_p61_t1_invr1145);
    cse_p63_t1_g1163 = (cse_p11_t1_g1142 * qq);
    cse_p63_t1_g1153 = (cse_p63_t1_g1163 * dx);
    cse_p63_t2_g1154 = (cse_p63_t1_g1163 * dy);
    cse_p63_t3_g1155 = (cse_p63_t1_g1163 * dz);
    cse_p11_t1_g1143 = (cse_p114_t1_g1165 * cse_p114_t1_g1165 * cse_p11_t1_g1142 * cse_p11_t23_g1115 * cse_p61_t1_invr1145);
    cse_p63_t1_g1160 = (-2.0000000000000000     * cse_p11_t1_g1143);
    g_x1 = ((cse_p63_t1_g1160 * dx) + (-(cse_p63_t1_g1153)));
    KernelGradientAcc(i3x1, 0, g_x1);
    g_y1 = ((cse_p63_t1_g1160 * dy) + (-(cse_p63_t2_g1154)));
    KernelGradientAcc(i3x1, 1, g_y1);
    g_z1 = ((cse_p63_t1_g1160 * dz) + (-(cse_p63_t3_g1155)));
    KernelGradientAcc(i3x1, 2, g_z1);
    cse_p63_t1_g1161 = (2.0000000000000000     * cse_p11_t1_g1143);
    g_x2 = (cse_p63_t1_g1153 + (cse_p63_t1_g1161 * dx));
    KernelGradientAcc(i3x2, 0, g_x2);
    g_y2 = (cse_p63_t2_g1154 + (cse_p63_t1_g1161 * dy));
    KernelGradientAcc(i3x2, 1, g_y2);
    g_z2 = (cse_p63_t3_g1155 + (cse_p63_t1_g1161 * dz));
    KernelGradientAcc(i3x2, 2, g_z2);
    cse_p11_t1_g1141 = (cse_p114_t1_g1165 * cse_p11_t24_g1116);
    cse_p63_t1_g1162 = (cse_p11_t1_g1141 * dx);
    h_x1_x1 = (cse_p11_t1_g1093 + (cse_p63_t1_g1162 * dx));
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 0, h_x1_x1);
    h_x1_y1 = (cse_p11_t14_g1106 + (cse_p63_t1_g1162 * dy));
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 1, h_x1_y1);
    h_x1_z1 = (cse_p11_t17_g1109 + (cse_p63_t1_g1162 * dz));
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 2, h_x1_z1);
    cse_p63_t1_g1156 = (cse_p11_t1_g1141 * cse_p11_t33_g1125);
    cse_p63_t1_g1164 = (cse_p11_t22_g1114 * cse_p61_t5_invr1149);
    h_x1_x2 = ((cse_p63_t1_g1156 * dx) + (cse_p63_t1_g1164 * (cse_p11_t11_g1103 + cse_p11_t32_g1124)));
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 0, h_x1_x2);
    cse_p63_t1_g1157 = (cse_p11_t1_g1141 * cse_p11_t34_g1126);
    h_x1_y2 = (cse_p11_t12_g1104 + (cse_p63_t1_g1157 * dx));
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 1, h_x1_y2);
    cse_p63_t1_g1159 = (cse_p11_t1_g1141 * cse_p11_t35_g1127);
    h_x1_z2 = (cse_p11_t13_g1105 + (cse_p63_t1_g1159 * dx));
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 2, h_x1_z2);
    h_y1_y1 = (cse_p11_t2_g1094 + (cse_p11_t1_g1141 * dy * dy));
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 1, h_y1_y1);
    h_y1_z1 = (cse_p11_t18_g1110 + (cse_p11_t1_g1141 * dy * dz));
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 2, h_y1_z1);
    h_y1_x2 = (cse_p11_t12_g1104 + (cse_p63_t1_g1156 * dy));
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 0, h_y1_x2);
    h_y1_y2 = ((cse_p63_t1_g1157 * dy) + (cse_p63_t1_g1164 * (cse_p11_t15_g1107 + cse_p11_t32_g1124)));
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 1, h_y1_y2);
    h_y1_z2 = (cse_p11_t16_g1108 + (cse_p63_t1_g1159 * dy));
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 2, h_y1_z2);
    h_z1_z1 = (cse_p11_t3_g1095 + (cse_p11_t1_g1141 * dz * dz));
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x1, 2, h_z1_z1);
    h_z1_x2 = (cse_p11_t13_g1105 + (cse_p63_t1_g1156 * dz));
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 0, h_z1_x2);
    h_z1_y2 = (cse_p11_t16_g1108 + (cse_p63_t1_g1157 * dz));
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 1, h_z1_y2);
    h_z1_z2 = ((cse_p63_t1_g1159 * dz) + (cse_p63_t1_g1164 * (cse_p11_t19_g1111 + cse_p11_t32_g1124)));
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 2, h_z1_z2);
    h_x2_x2 = (cse_p11_t1_g1093 + (cse_p11_t33_g1125 * cse_p63_t1_g1156));
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 0, h_x2_x2);
    h_x2_y2 = (cse_p11_t14_g1106 + (cse_p11_t34_g1126 * cse_p63_t1_g1156));
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 1, h_x2_y2);
    h_x2_z2 = (cse_p11_t17_g1109 + (cse_p11_t35_g1127 * cse_p63_t1_g1156));
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 2, h_x2_z2);
    h_y2_y2 = (cse_p11_t2_g1094 + (cse_p11_t34_g1126 * cse_p63_t1_g1157));
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 1, h_y2_y2);
    h_y2_z2 = (cse_p11_t18_g1110 + (cse_p11_t35_g1127 * cse_p63_t1_g1157));
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 2, h_y2_z2);
    h_z2_z2 = (cse_p11_t3_g1095 + (cse_p11_t35_g1127 * cse_p63_t1_g1159));
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x2, 2, h_z2_z2);
  }
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
