template <typename HESSIAN>
struct Dihedral_Restraint_Fast {
  static constexpr size_t PositionSize = 12;
  static std::string description() { return "mathkernel-dihedral_restraint_fast"; };
double energy(double kdh, double phi0, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double energy_added = 0.0000000000000000    ;
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
      double dx12 = (x2 + (-(x1)));
      double dy12 = (y2 + (-(y1)));
      double dz12 = (z2 + (-(z1)));
      double dx23 = (x3 + (-(x2)));
      double dy23 = (y3 + (-(y2)));
      double dz23 = (z3 + (-(z2)));
      double dx34 = (x4 + (-(x3)));
      double dy34 = (y4 + (-(y3)));
      double dz34 = (z4 + (-(z3)));
      double n1x = ((dy12 * dz23) + (-((dy23 * dz12))));
      double n1y = ((dx23 * dz12) + (-((dx12 * dz23))));
      double n1z = ((dx12 * dy23) + (-((dx23 * dy12))));
      double n2x = ((dy23 * dz34) + (-((dy34 * dz23))));
      double n2y = ((dx34 * dz23) + (-((dx23 * dz34))));
      double n2z = ((dx23 * dy34) + (-((dx34 * dy23))));
      double n1_2 = ((n1x * n1x) + (n1y * n1y) + (n1z * n1z));
      double n2_2 = ((n2x * n2x) + (n2y * n2y) + (n2z * n2z));
      double inv_n1n2 = (1.0 / sqrt((n1_2 * n2_2)));
      double dot_n1n2 = ((n1x * n2x) + (n1y * n2y) + (n1z * n2z));
      double cosphi = (dot_n1n2 * inv_n1n2);
      double phi = acos(cosphi);
      double deltaphi = (phi + (-(phi0)));
      double energy = (deltaphi * deltaphi * kdh);
      energy_added += energy;
      *energy_accumulate += energy;
    }
  }
  return energy_added;
}
void energy_fd(double kdh, double phi0, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, energy_accumulate, force, hessian, dvec, hdvec);
}

double gradient(double kdh, double phi0, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double energy_added = 0.0000000000000000    ;
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
      double dx12 = (x2 + (-(x1)));
      double dy12 = (y2 + (-(y1)));
      double dz12 = (z2 + (-(z1)));
      double dx23 = (x3 + (-(x2)));
      double cse_p12_t37_g44204 = ((-(dx12)) + (-(dx23)));
      double cse_p12_t43_g44210 = (dx12 + dx23);
      double dy23 = (y3 + (-(y2)));
      double cse_p12_t39_g44206 = ((-(dy12)) + (-(dy23)));
      double cse_p12_t45_g44212 = (dy12 + dy23);
      double dz23 = (z3 + (-(z2)));
      double cse_p12_t41_g44208 = ((-(dz12)) + (-(dz23)));
      double cse_p12_t47_g44214 = (dz12 + dz23);
      double dx34 = (x4 + (-(x3)));
      double cse_p12_t38_g44205 = ((-(dx23)) + (-(dx34)));
      double cse_p12_t44_g44211 = (dx23 + dx34);
      double dy34 = (y4 + (-(y3)));
      double cse_p12_t40_g44207 = ((-(dy23)) + (-(dy34)));
      double cse_p12_t46_g44213 = (dy23 + dy34);
      double dz34 = (z4 + (-(z3)));
      double cse_p12_t42_g44209 = ((-(dz23)) + (-(dz34)));
      double cse_p12_t48_g44215 = (dz23 + dz34);
      double n1x = ((dy12 * dz23) + (-((dy23 * dz12))));
      double cse_p12_t5_g44172 = (dy12 * n1x);
      double cse_p12_t6_g44173 = (dy23 * n1x);
      double cse_p12_t13_g44180 = (cse_p12_t39_g44206 * n1x);
      double cse_p12_t14_g44181 = (cse_p12_t47_g44214 * n1x);
      double cse_p12_t15_g44182 = (n1x * (-(dz12)));
      double cse_p12_t16_g44183 = (n1x * (-(dz23)));
      double n1y = ((dx23 * dz12) + (-((dx12 * dz23))));
      double cse_p12_t9_g44176 = (dz12 * n1y);
      double cse_p12_t10_g44177 = (dz23 * n1y);
      double cse_p12_t17_g44184 = (cse_p12_t41_g44208 * n1y);
      double cse_p12_t18_g44185 = (cse_p12_t43_g44210 * n1y);
      double cse_p12_t19_g44186 = (n1y * (-(dx12)));
      double cse_p12_t20_g44187 = (n1y * (-(dx23)));
      double n1z = ((dx12 * dy23) + (-((dx23 * dy12))));
      double cse_p12_t1_g44168 = (dx12 * n1z);
      double cse_p12_t2_g44169 = (dx23 * n1z);
      double cse_p12_t21_g44188 = (cse_p12_t37_g44204 * n1z);
      double cse_p12_t22_g44189 = (cse_p12_t45_g44212 * n1z);
      double cse_p12_t23_g44190 = (n1z * (-(dy12)));
      double cse_p12_t24_g44191 = (n1z * (-(dy23)));
      double n2x = ((dy23 * dz34) + (-((dy34 * dz23))));
      double cse_p12_t7_g44174 = (dy23 * n2x);
      double cse_p12_t8_g44175 = (dy34 * n2x);
      double cse_p12_t25_g44192 = (cse_p12_t40_g44207 * n2x);
      double cse_p12_t26_g44193 = (cse_p12_t48_g44215 * n2x);
      double cse_p12_t27_g44194 = (n2x * (-(dz23)));
      double cse_p12_t28_g44195 = (n2x * (-(dz34)));
      double n2y = ((dx34 * dz23) + (-((dx23 * dz34))));
      double cse_p12_t11_g44178 = (dz23 * n2y);
      double cse_p12_t12_g44179 = (dz34 * n2y);
      double cse_p12_t29_g44196 = (cse_p12_t42_g44209 * n2y);
      double cse_p12_t30_g44197 = (cse_p12_t44_g44211 * n2y);
      double cse_p12_t31_g44198 = (n2y * (-(dx23)));
      double cse_p12_t32_g44199 = (n2y * (-(dx34)));
      double n2z = ((dx23 * dy34) + (-((dx34 * dy23))));
      double cse_p12_t3_g44170 = (dx23 * n2z);
      double cse_p12_t4_g44171 = (dx34 * n2z);
      double cse_p12_t33_g44200 = (cse_p12_t38_g44205 * n2z);
      double cse_p12_t34_g44201 = (cse_p12_t46_g44213 * n2z);
      double cse_p12_t35_g44202 = (n2z * (-(dy23)));
      double cse_p12_t36_g44203 = (n2z * (-(dy34)));
      double n1_2 = ((n1x * n1x) + (n1y * n1y) + (n1z * n1z));
      double n2_2 = ((n2x * n2x) + (n2y * n2y) + (n2z * n2z));
      double cse_p1_t1_invsqrt44165 = (1.0 / sqrt((n1_2 * n2_2)));
      double cse_p1_t2_invr44166 = (cse_p1_t1_invsqrt44165 * cse_p1_t1_invsqrt44165);
      double dot_n1n2 = ((n1x * n2x) + (n1y * n2y) + (n1z * n2z));
      double cosphi = (cse_p1_t1_invsqrt44165 * dot_n1n2);
      double phi = acos(cosphi);
      double deltaphi = (phi + (-(phi0)));
      double cse_p12_t1_g44217 = (deltaphi * kdh);
      double energy = (cse_p12_t1_g44217 * deltaphi);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p1_t3_invsqrt44167 = (1.0 / sqrt((1.0000000000000000     + (-(((cosphi) * (cosphi)))))));
      double cse_p12_t1_g44216 = (-0.50000000000000000     * cse_p1_t1_invsqrt44165 * cse_p1_t2_invr44166 * dot_n1n2);
      double cse_p12_t1_g44219 = (cse_p12_t1_g44216 * n2_2);
      double cse_p21_t1_g44220 = (cse_p12_t1_g44217 * cse_p1_t3_invsqrt44167);
      double g_x1 = ((-((cse_p21_t1_g44220 * ((cse_p12_t1_g44219 * (cse_p12_t10_g44177 + cse_p12_t10_g44177 + cse_p12_t24_g44191 + cse_p12_t24_g44191)) + (cse_p1_t1_invsqrt44165 * (cse_p12_t11_g44178 + cse_p12_t35_g44202)))))) + (-((cse_p21_t1_g44220 * ((cse_p12_t1_g44219 * (cse_p12_t10_g44177 + cse_p12_t10_g44177 + cse_p12_t24_g44191 + cse_p12_t24_g44191)) + (cse_p1_t1_invsqrt44165 * (cse_p12_t11_g44178 + cse_p12_t35_g44202)))))));
      KernelGradientAcc(i3x1, 0, g_x1);
      double g_y1 = ((-((cse_p21_t1_g44220 * ((cse_p12_t1_g44219 * (cse_p12_t16_g44183 + cse_p12_t16_g44183 + cse_p12_t2_g44169 + cse_p12_t2_g44169)) + (cse_p1_t1_invsqrt44165 * (cse_p12_t27_g44194 + cse_p12_t3_g44170)))))) + (-((cse_p21_t1_g44220 * ((cse_p12_t1_g44219 * (cse_p12_t16_g44183 + cse_p12_t16_g44183 + cse_p12_t2_g44169 + cse_p12_t2_g44169)) + (cse_p1_t1_invsqrt44165 * (cse_p12_t27_g44194 + cse_p12_t3_g44170)))))));
      KernelGradientAcc(i3x1, 1, g_y1);
      double g_z1 = ((-((cse_p21_t1_g44220 * ((cse_p12_t1_g44219 * (cse_p12_t20_g44187 + cse_p12_t20_g44187 + cse_p12_t6_g44173 + cse_p12_t6_g44173)) + (cse_p1_t1_invsqrt44165 * (cse_p12_t31_g44198 + cse_p12_t7_g44174)))))) + (-((cse_p21_t1_g44220 * ((cse_p12_t1_g44219 * (cse_p12_t20_g44187 + cse_p12_t20_g44187 + cse_p12_t6_g44173 + cse_p12_t6_g44173)) + (cse_p1_t1_invsqrt44165 * (cse_p12_t31_g44198 + cse_p12_t7_g44174)))))));
      KernelGradientAcc(i3x1, 2, g_z1);
      double g_x2 = ((-((cse_p21_t1_g44220 * ((cse_p12_t1_g44216 * ((n1_2 * (cse_p12_t12_g44179 + cse_p12_t12_g44179 + cse_p12_t36_g44203 + cse_p12_t36_g44203)) + (n2_2 * (cse_p12_t17_g44184 + cse_p12_t17_g44184 + cse_p12_t22_g44189 + cse_p12_t22_g44189)))) + (cse_p1_t1_invsqrt44165 * ((cse_p12_t41_g44208 * n2y) + (cse_p12_t45_g44212 * n2z) + (dz34 * n1y) + (n1z * (-(dy34))))))))) + (-((cse_p21_t1_g44220 * ((cse_p12_t1_g44216 * ((n1_2 * (cse_p12_t12_g44179 + cse_p12_t12_g44179 + cse_p12_t36_g44203 + cse_p12_t36_g44203)) + (n2_2 * (cse_p12_t17_g44184 + cse_p12_t17_g44184 + cse_p12_t22_g44189 + cse_p12_t22_g44189)))) + (cse_p1_t1_invsqrt44165 * ((cse_p12_t41_g44208 * n2y) + (cse_p12_t45_g44212 * n2z) + (dz34 * n1y) + (n1z * (-(dy34))))))))));
      KernelGradientAcc(i3x2, 0, g_x2);
      double g_y2 = ((-((cse_p21_t1_g44220 * ((cse_p12_t1_g44216 * ((n1_2 * (cse_p12_t28_g44195 + cse_p12_t28_g44195 + cse_p12_t4_g44171 + cse_p12_t4_g44171)) + (n2_2 * (cse_p12_t14_g44181 + cse_p12_t14_g44181 + cse_p12_t21_g44188 + cse_p12_t21_g44188)))) + (cse_p1_t1_invsqrt44165 * ((cse_p12_t37_g44204 * n2z) + (cse_p12_t47_g44214 * n2x) + (dx34 * n1z) + (n1x * (-(dz34))))))))) + (-((cse_p21_t1_g44220 * ((cse_p12_t1_g44216 * ((n1_2 * (cse_p12_t28_g44195 + cse_p12_t28_g44195 + cse_p12_t4_g44171 + cse_p12_t4_g44171)) + (n2_2 * (cse_p12_t14_g44181 + cse_p12_t14_g44181 + cse_p12_t21_g44188 + cse_p12_t21_g44188)))) + (cse_p1_t1_invsqrt44165 * ((cse_p12_t37_g44204 * n2z) + (cse_p12_t47_g44214 * n2x) + (dx34 * n1z) + (n1x * (-(dz34))))))))));
      KernelGradientAcc(i3x2, 1, g_y2);
      double g_z2 = ((-((cse_p21_t1_g44220 * ((cse_p12_t1_g44216 * ((n1_2 * (cse_p12_t32_g44199 + cse_p12_t32_g44199 + cse_p12_t8_g44175 + cse_p12_t8_g44175)) + (n2_2 * (cse_p12_t13_g44180 + cse_p12_t13_g44180 + cse_p12_t18_g44185 + cse_p12_t18_g44185)))) + (cse_p1_t1_invsqrt44165 * ((cse_p12_t39_g44206 * n2x) + (cse_p12_t43_g44210 * n2y) + (dy34 * n1x) + (n1y * (-(dx34))))))))) + (-((cse_p21_t1_g44220 * ((cse_p12_t1_g44216 * ((n1_2 * (cse_p12_t32_g44199 + cse_p12_t32_g44199 + cse_p12_t8_g44175 + cse_p12_t8_g44175)) + (n2_2 * (cse_p12_t13_g44180 + cse_p12_t13_g44180 + cse_p12_t18_g44185 + cse_p12_t18_g44185)))) + (cse_p1_t1_invsqrt44165 * ((cse_p12_t39_g44206 * n2x) + (cse_p12_t43_g44210 * n2y) + (dy34 * n1x) + (n1y * (-(dx34))))))))));
      KernelGradientAcc(i3x2, 2, g_z2);
      double g_x3 = ((-((cse_p21_t1_g44220 * ((cse_p12_t1_g44216 * ((n1_2 * (cse_p12_t29_g44196 + cse_p12_t29_g44196 + cse_p12_t34_g44201 + cse_p12_t34_g44201)) + (n2_2 * (cse_p12_t23_g44190 + cse_p12_t23_g44190 + cse_p12_t9_g44176 + cse_p12_t9_g44176)))) + (cse_p1_t1_invsqrt44165 * ((cse_p12_t42_g44209 * n1y) + (cse_p12_t46_g44213 * n1z) + (dz12 * n2y) + (n2z * (-(dy12))))))))) + (-((cse_p21_t1_g44220 * ((cse_p12_t1_g44216 * ((n1_2 * (cse_p12_t29_g44196 + cse_p12_t29_g44196 + cse_p12_t34_g44201 + cse_p12_t34_g44201)) + (n2_2 * (cse_p12_t23_g44190 + cse_p12_t23_g44190 + cse_p12_t9_g44176 + cse_p12_t9_g44176)))) + (cse_p1_t1_invsqrt44165 * ((cse_p12_t42_g44209 * n1y) + (cse_p12_t46_g44213 * n1z) + (dz12 * n2y) + (n2z * (-(dy12))))))))));
      KernelGradientAcc(i3x3, 0, g_x3);
      double g_y3 = ((-((cse_p21_t1_g44220 * ((cse_p12_t1_g44216 * ((n1_2 * (cse_p12_t26_g44193 + cse_p12_t26_g44193 + cse_p12_t33_g44200 + cse_p12_t33_g44200)) + (n2_2 * (cse_p12_t15_g44182 + cse_p12_t15_g44182 + cse_p12_t1_g44168 + cse_p12_t1_g44168)))) + (cse_p1_t1_invsqrt44165 * ((cse_p12_t38_g44205 * n1z) + (cse_p12_t48_g44215 * n1x) + (dx12 * n2z) + (n2x * (-(dz12))))))))) + (-((cse_p21_t1_g44220 * ((cse_p12_t1_g44216 * ((n1_2 * (cse_p12_t26_g44193 + cse_p12_t26_g44193 + cse_p12_t33_g44200 + cse_p12_t33_g44200)) + (n2_2 * (cse_p12_t15_g44182 + cse_p12_t15_g44182 + cse_p12_t1_g44168 + cse_p12_t1_g44168)))) + (cse_p1_t1_invsqrt44165 * ((cse_p12_t38_g44205 * n1z) + (cse_p12_t48_g44215 * n1x) + (dx12 * n2z) + (n2x * (-(dz12))))))))));
      KernelGradientAcc(i3x3, 1, g_y3);
      double g_z3 = ((-((cse_p21_t1_g44220 * ((cse_p12_t1_g44216 * ((n1_2 * (cse_p12_t25_g44192 + cse_p12_t25_g44192 + cse_p12_t30_g44197 + cse_p12_t30_g44197)) + (n2_2 * (cse_p12_t19_g44186 + cse_p12_t19_g44186 + cse_p12_t5_g44172 + cse_p12_t5_g44172)))) + (cse_p1_t1_invsqrt44165 * ((cse_p12_t40_g44207 * n1x) + (cse_p12_t44_g44211 * n1y) + (dy12 * n2x) + (n2y * (-(dx12))))))))) + (-((cse_p21_t1_g44220 * ((cse_p12_t1_g44216 * ((n1_2 * (cse_p12_t25_g44192 + cse_p12_t25_g44192 + cse_p12_t30_g44197 + cse_p12_t30_g44197)) + (n2_2 * (cse_p12_t19_g44186 + cse_p12_t19_g44186 + cse_p12_t5_g44172 + cse_p12_t5_g44172)))) + (cse_p1_t1_invsqrt44165 * ((cse_p12_t40_g44207 * n1x) + (cse_p12_t44_g44211 * n1y) + (dy12 * n2x) + (n2y * (-(dx12))))))))));
      KernelGradientAcc(i3x3, 2, g_z3);
      double cse_p12_t1_g44218 = (cse_p12_t1_g44216 * n1_2);
      double g_x4 = ((-((cse_p21_t1_g44220 * ((cse_p12_t1_g44218 * (cse_p12_t11_g44178 + cse_p12_t11_g44178 + cse_p12_t35_g44202 + cse_p12_t35_g44202)) + (cse_p1_t1_invsqrt44165 * (cse_p12_t10_g44177 + cse_p12_t24_g44191)))))) + (-((cse_p21_t1_g44220 * ((cse_p12_t1_g44218 * (cse_p12_t11_g44178 + cse_p12_t11_g44178 + cse_p12_t35_g44202 + cse_p12_t35_g44202)) + (cse_p1_t1_invsqrt44165 * (cse_p12_t10_g44177 + cse_p12_t24_g44191)))))));
      KernelGradientAcc(i3x4, 0, g_x4);
      double g_y4 = ((-((cse_p21_t1_g44220 * ((cse_p12_t1_g44218 * (cse_p12_t27_g44194 + cse_p12_t27_g44194 + cse_p12_t3_g44170 + cse_p12_t3_g44170)) + (cse_p1_t1_invsqrt44165 * (cse_p12_t16_g44183 + cse_p12_t2_g44169)))))) + (-((cse_p21_t1_g44220 * ((cse_p12_t1_g44218 * (cse_p12_t27_g44194 + cse_p12_t27_g44194 + cse_p12_t3_g44170 + cse_p12_t3_g44170)) + (cse_p1_t1_invsqrt44165 * (cse_p12_t16_g44183 + cse_p12_t2_g44169)))))));
      KernelGradientAcc(i3x4, 1, g_y4);
      double g_z4 = ((-((cse_p21_t1_g44220 * ((cse_p12_t1_g44218 * (cse_p12_t31_g44198 + cse_p12_t31_g44198 + cse_p12_t7_g44174 + cse_p12_t7_g44174)) + (cse_p1_t1_invsqrt44165 * (cse_p12_t20_g44187 + cse_p12_t6_g44173)))))) + (-((cse_p21_t1_g44220 * ((cse_p12_t1_g44218 * (cse_p12_t31_g44198 + cse_p12_t31_g44198 + cse_p12_t7_g44174 + cse_p12_t7_g44174)) + (cse_p1_t1_invsqrt44165 * (cse_p12_t20_g44187 + cse_p12_t6_g44173)))))));
      KernelGradientAcc(i3x4, 2, g_z4);
    }
  }
  return energy_added;
}
void gradient_fd(double kdh, double phi0, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 12;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  double e0 = 0.0;
  energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
  {
    double saved = position[i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 0] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 0, d);
  }
  {
    double saved = position[i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 1] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 1, d);
  }
  {
    double saved = position[i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 2] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 2, d);
  }
  {
    double saved = position[i3x4 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 0] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 0, d);
  }
  {
    double saved = position[i3x4 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 1] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 1, d);
  }
  {
    double saved = position[i3x4 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 2] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 2, d);
  }
}

double hessian(double kdh, double phi0, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double energy_added = 0.0000000000000000    ;
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
      double dx12 = (x2 + (-(x1)));
      double dy12 = (y2 + (-(y1)));
      double dz12 = (z2 + (-(z1)));
      double dx23 = (x3 + (-(x2)));
      double cse_p12_t157_g44383 = ((-(dx12)) + (-(dx23)));
      double cse_p12_t163_g44389 = (dx12 + dx23);
      double dy23 = (y3 + (-(y2)));
      double cse_p12_t159_g44385 = ((-(dy12)) + (-(dy23)));
      double cse_p12_t165_g44391 = (dy12 + dy23);
      double dz23 = (z3 + (-(z2)));
      double cse_p12_t161_g44387 = ((-(dz12)) + (-(dz23)));
      double cse_p12_t167_g44393 = (dz12 + dz23);
      double dx34 = (x4 + (-(x3)));
      double cse_p12_t158_g44384 = ((-(dx23)) + (-(dx34)));
      double cse_p12_t164_g44390 = (dx23 + dx34);
      double dy34 = (y4 + (-(y3)));
      double cse_p12_t160_g44386 = ((-(dy23)) + (-(dy34)));
      double cse_p12_t166_g44392 = (dy23 + dy34);
      double dz34 = (z4 + (-(z3)));
      double cse_p12_t162_g44388 = ((-(dz23)) + (-(dz34)));
      double cse_p12_t168_g44394 = (dz23 + dz34);
      double n1x = ((dy12 * dz23) + (-((dy23 * dz12))));
      double cse_p12_t31_g44257 = (dy12 * n1x);
      double cse_p12_t33_g44259 = (dy23 * n1x);
      double cse_p12_t35_g44261 = (dy34 * n1x);
      double cse_p13_t1_g44426 = (cse_p12_t159_g44385 * n1x);
      double cse_p12_t56_g44282 = (cse_p12_t160_g44386 * n1x);
      double cse_p12_t57_g44283 = (cse_p12_t167_g44393 * n1x);
      double cse_p12_t58_g44284 = (cse_p12_t168_g44394 * n1x);
      double cse_p12_t59_g44285 = (n1x * (-(dz12)));
      double cse_p13_t1_g44422 = (n1x * (-(dz23)));
      double cse_p12_t61_g44287 = (n1x * (-(dz34)));
      double n1y = ((dx23 * dz12) + (-((dx12 * dz23))));
      double cse_p12_t37_g44263 = (dz12 * n1y);
      double cse_p12_t39_g44265 = (dz23 * n1y);
      double cse_p12_t41_g44267 = (dz34 * n1y);
      double cse_p12_t62_g44288 = (cse_p12_t161_g44387 * n1y);
      double cse_p12_t63_g44289 = (cse_p12_t162_g44388 * n1y);
      double cse_p12_t64_g44290 = (cse_p12_t163_g44389 * n1y);
      double cse_p12_t65_g44291 = (cse_p12_t164_g44390 * n1y);
      double cse_p12_t66_g44292 = (n1y * (-(dx12)));
      double cse_p13_t1_g44418 = (n1y * (-(dx23)));
      double cse_p12_t68_g44294 = (n1y * (-(dx34)));
      double cse_p12_t135_g44361 = (cse_p12_t31_g44257 + cse_p12_t31_g44257 + cse_p12_t66_g44292 + cse_p12_t66_g44292);
      double cse_p12_t137_g44363 = (cse_p12_t33_g44259 + cse_p12_t33_g44259 + cse_p13_t1_g44418 + cse_p13_t1_g44418);
      double cse_p12_t138_g44364 = (cse_p12_t33_g44259 + cse_p13_t1_g44418);
      double cse_p12_t151_g44377 = (cse_p12_t64_g44290 + cse_p12_t64_g44290 + cse_p13_t1_g44426 + cse_p13_t1_g44426);
      double n1z = ((dx12 * dy23) + (-((dx23 * dy12))));
      double cse_p12_t25_g44251 = (dx12 * n1z);
      double cse_p12_t27_g44253 = (dx23 * n1z);
      double cse_p12_t29_g44255 = (dx34 * n1z);
      double cse_p13_t1_g44424 = (cse_p12_t157_g44383 * n1z);
      double cse_p12_t70_g44296 = (cse_p12_t158_g44384 * n1z);
      double cse_p12_t71_g44297 = (cse_p12_t165_g44391 * n1z);
      double cse_p12_t72_g44298 = (cse_p12_t166_g44392 * n1z);
      double cse_p12_t73_g44299 = (n1z * (-(dy12)));
      double cse_p13_t1_g44420 = (n1z * (-(dy23)));
      double cse_p12_t75_g44301 = (n1z * (-(dy34)));
      double cse_p12_t127_g44353 = (cse_p12_t25_g44251 + cse_p12_t25_g44251 + cse_p12_t59_g44285 + cse_p12_t59_g44285);
      double cse_p12_t129_g44355 = (cse_p12_t27_g44253 + cse_p12_t27_g44253 + cse_p13_t1_g44422 + cse_p13_t1_g44422);
      double cse_p12_t130_g44356 = (cse_p12_t27_g44253 + cse_p13_t1_g44422);
      double cse_p12_t143_g44369 = (cse_p12_t37_g44263 + cse_p12_t37_g44263 + cse_p12_t73_g44299 + cse_p12_t73_g44299);
      double cse_p12_t145_g44371 = (cse_p12_t39_g44265 + cse_p12_t39_g44265 + cse_p13_t1_g44420 + cse_p13_t1_g44420);
      double cse_p12_t146_g44372 = (cse_p12_t39_g44265 + cse_p13_t1_g44420);
      double cse_p12_t152_g44378 = (cse_p12_t57_g44283 + cse_p12_t57_g44283 + cse_p13_t1_g44424 + cse_p13_t1_g44424);
      double cse_p12_t153_g44379 = (cse_p12_t62_g44288 + cse_p12_t62_g44288 + cse_p12_t71_g44297 + cse_p12_t71_g44297);
      double n2x = ((dy23 * dz34) + (-((dy34 * dz23))));
      double cse_p12_t32_g44258 = (dy12 * n2x);
      double cse_p12_t34_g44260 = (dy23 * n2x);
      double cse_p12_t36_g44262 = (dy34 * n2x);
      double cse_p12_t82_g44308 = (cse_p12_t159_g44385 * n2x);
      double cse_p12_t83_g44309 = (cse_p12_t160_g44386 * n2x);
      double cse_p12_t84_g44310 = (cse_p12_t167_g44393 * n2x);
      double cse_p12_t85_g44311 = (cse_p12_t168_g44394 * n2x);
      double cse_p12_t86_g44312 = (n2x * (-(dz12)));
      double cse_p13_t1_g44423 = (n2x * (-(dz23)));
      double cse_p12_t88_g44314 = (n2x * (-(dz34)));
      double n2y = ((dx34 * dz23) + (-((dx23 * dz34))));
      double cse_p12_t38_g44264 = (dz12 * n2y);
      double cse_p12_t40_g44266 = (dz23 * n2y);
      double cse_p12_t42_g44268 = (dz34 * n2y);
      double cse_p12_t89_g44315 = (cse_p12_t161_g44387 * n2y);
      double cse_p12_t90_g44316 = (cse_p12_t162_g44388 * n2y);
      double cse_p12_t91_g44317 = (cse_p12_t163_g44389 * n2y);
      double cse_p12_t92_g44318 = (cse_p12_t164_g44390 * n2y);
      double cse_p12_t93_g44319 = (n2y * (-(dx12)));
      double cse_p13_t1_g44419 = (n2y * (-(dx23)));
      double cse_p12_t95_g44321 = (n2y * (-(dx34)));
      double cse_p12_t136_g44362 = (cse_p12_t32_g44258 + cse_p12_t56_g44282 + cse_p12_t65_g44291 + cse_p12_t93_g44319);
      double cse_p12_t139_g44365 = (cse_p12_t34_g44260 + cse_p12_t34_g44260 + cse_p13_t1_g44419 + cse_p13_t1_g44419);
      double cse_p12_t140_g44366 = (cse_p12_t34_g44260 + cse_p13_t1_g44419);
      double cse_p12_t141_g44367 = (cse_p12_t35_g44261 + cse_p12_t68_g44294 + cse_p12_t82_g44308 + cse_p12_t91_g44317);
      double cse_p12_t142_g44368 = (cse_p12_t36_g44262 + cse_p12_t36_g44262 + cse_p12_t95_g44321 + cse_p12_t95_g44321);
      double cse_p12_t154_g44380 = (cse_p12_t83_g44309 + cse_p12_t83_g44309 + cse_p12_t92_g44318 + cse_p12_t92_g44318);
      double n2z = ((dx23 * dy34) + (-((dx34 * dy23))));
      double cse_p12_t26_g44252 = (dx12 * n2z);
      double cse_p12_t28_g44254 = (dx23 * n2z);
      double cse_p12_t30_g44256 = (dx34 * n2z);
      double cse_p12_t96_g44322 = (cse_p12_t157_g44383 * n2z);
      double cse_p13_t1_g44425 = (cse_p12_t158_g44384 * n2z);
      double cse_p12_t98_g44324 = (cse_p12_t165_g44391 * n2z);
      double cse_p12_t99_g44325 = (cse_p12_t166_g44392 * n2z);
      double cse_p12_t100_g44326 = (n2z * (-(dy12)));
      double cse_p13_t1_g44421 = (n2z * (-(dy23)));
      double cse_p12_t102_g44328 = (n2z * (-(dy34)));
      double cse_p12_t128_g44354 = (cse_p12_t26_g44252 + cse_p12_t58_g44284 + cse_p12_t70_g44296 + cse_p12_t86_g44312);
      double cse_p12_t131_g44357 = (cse_p12_t28_g44254 + cse_p12_t28_g44254 + cse_p13_t1_g44423 + cse_p13_t1_g44423);
      double cse_p12_t132_g44358 = (cse_p12_t28_g44254 + cse_p13_t1_g44423);
      double cse_p12_t133_g44359 = (cse_p12_t29_g44255 + cse_p12_t61_g44287 + cse_p12_t84_g44310 + cse_p12_t96_g44322);
      double cse_p12_t134_g44360 = (cse_p12_t30_g44256 + cse_p12_t30_g44256 + cse_p12_t88_g44314 + cse_p12_t88_g44314);
      double cse_p12_t144_g44370 = (cse_p12_t100_g44326 + cse_p12_t38_g44264 + cse_p12_t63_g44289 + cse_p12_t72_g44298);
      double cse_p12_t147_g44373 = (cse_p12_t40_g44266 + cse_p12_t40_g44266 + cse_p13_t1_g44421 + cse_p13_t1_g44421);
      double cse_p12_t148_g44374 = (cse_p12_t40_g44266 + cse_p13_t1_g44421);
      double cse_p12_t149_g44375 = (cse_p12_t41_g44267 + cse_p12_t75_g44301 + cse_p12_t89_g44315 + cse_p12_t98_g44324);
      double cse_p12_t150_g44376 = (cse_p12_t102_g44328 + cse_p12_t102_g44328 + cse_p12_t42_g44268 + cse_p12_t42_g44268);
      double cse_p12_t155_g44381 = (cse_p12_t85_g44311 + cse_p12_t85_g44311 + cse_p13_t1_g44425 + cse_p13_t1_g44425);
      double cse_p12_t156_g44382 = (cse_p12_t90_g44316 + cse_p12_t90_g44316 + cse_p12_t99_g44325 + cse_p12_t99_g44325);
      double n1_2 = ((n1x * n1x) + (n1y * n1y) + (n1z * n1z));
      double cse_p12_t76_g44302 = (cse_p12_t134_g44360 * n1_2);
      double cse_p12_t77_g44303 = (cse_p12_t142_g44368 * n1_2);
      double cse_p12_t78_g44304 = (cse_p12_t150_g44376 * n1_2);
      double cse_p12_t79_g44305 = (cse_p12_t154_g44380 * n1_2);
      double cse_p12_t80_g44306 = (cse_p12_t155_g44381 * n1_2);
      double cse_p12_t81_g44307 = (cse_p12_t156_g44382 * n1_2);
      double n2_2 = ((n2x * n2x) + (n2y * n2y) + (n2z * n2z));
      double cse_p12_t103_g44329 = (cse_p12_t127_g44353 * n2_2);
      double cse_p12_t104_g44330 = (cse_p12_t135_g44361 * n2_2);
      double cse_p12_t105_g44331 = (cse_p12_t143_g44369 * n2_2);
      double cse_p12_t106_g44332 = (cse_p12_t151_g44377 * n2_2);
      double cse_p12_t107_g44333 = (cse_p12_t152_g44378 * n2_2);
      double cse_p12_t108_g44334 = (cse_p12_t153_g44379 * n2_2);
      double cse_p12_t121_g44347 = (cse_p12_t107_g44333 + cse_p12_t76_g44302);
      double cse_p12_t122_g44348 = (cse_p12_t106_g44332 + cse_p12_t77_g44303);
      double cse_p12_t123_g44349 = (cse_p12_t108_g44334 + cse_p12_t78_g44304);
      double cse_p12_t124_g44350 = (cse_p12_t104_g44330 + cse_p12_t79_g44305);
      double cse_p12_t125_g44351 = (cse_p12_t103_g44329 + cse_p12_t80_g44306);
      double cse_p12_t126_g44352 = (cse_p12_t105_g44331 + cse_p12_t81_g44307);
      double cse_p11_t1_invsqrt44224 = (1.0 / sqrt((n1_2 * n2_2)));
      double cse_p11_t2_invr44225 = (cse_p11_t1_invsqrt44224 * cse_p11_t1_invsqrt44224);
      double cse_p12_t43_g44269 = (cse_p11_t1_invsqrt44224 * cse_p12_t128_g44354);
      double cse_p12_t44_g44270 = (cse_p11_t1_invsqrt44224 * cse_p12_t130_g44356);
      double cse_p12_t45_g44271 = (cse_p11_t1_invsqrt44224 * cse_p12_t132_g44358);
      double cse_p12_t46_g44272 = (cse_p11_t1_invsqrt44224 * cse_p12_t133_g44359);
      double cse_p12_t47_g44273 = (cse_p11_t1_invsqrt44224 * cse_p12_t136_g44362);
      double cse_p12_t48_g44274 = (cse_p11_t1_invsqrt44224 * cse_p12_t138_g44364);
      double cse_p12_t49_g44275 = (cse_p11_t1_invsqrt44224 * cse_p12_t140_g44366);
      double cse_p12_t50_g44276 = (cse_p11_t1_invsqrt44224 * cse_p12_t141_g44367);
      double cse_p12_t51_g44277 = (cse_p11_t1_invsqrt44224 * cse_p12_t144_g44370);
      double cse_p12_t52_g44278 = (cse_p11_t1_invsqrt44224 * cse_p12_t146_g44372);
      double cse_p12_t53_g44279 = (cse_p11_t1_invsqrt44224 * cse_p12_t148_g44374);
      double cse_p12_t54_g44280 = (cse_p11_t1_invsqrt44224 * cse_p12_t149_g44375);
      double dot_n1n2 = ((n1x * n2x) + (n1y * n2y) + (n1z * n2z));
      double cse_p12_t1_g44408 = (-0.50000000000000000     * dot_n1n2);
      double cse_p12_t1_g44409 = (cse_p11_t1_invsqrt44224 * cse_p11_t2_invr44225 * cse_p12_t1_g44408);
      double cse_p12_t1_g44410 = (cse_p12_t121_g44347 * cse_p12_t1_g44409);
      double cse_p12_t1_g44411 = (cse_p12_t122_g44348 * cse_p12_t1_g44409);
      double cse_p12_t1_g44412 = (cse_p12_t123_g44349 * cse_p12_t1_g44409);
      double cse_p12_t1_g44413 = (cse_p12_t124_g44350 * cse_p12_t1_g44409);
      double cse_p12_t1_g44414 = (cse_p12_t125_g44351 * cse_p12_t1_g44409);
      double cse_p12_t1_g44415 = (cse_p12_t126_g44352 * cse_p12_t1_g44409);
      double cse_p14_t1_g44434 = (cse_p12_t1_g44409 * n1_2);
      double cse_p12_t7_g44233 = (cse_p12_t131_g44357 * cse_p14_t1_g44434);
      double cse_p12_t8_g44234 = (cse_p12_t139_g44365 * cse_p14_t1_g44434);
      double cse_p12_t9_g44235 = (cse_p12_t147_g44373 * cse_p14_t1_g44434);
      double cse_p14_t1_g44435 = (cse_p12_t1_g44409 * n2_2);
      double cse_p12_t10_g44236 = (cse_p12_t129_g44355 * cse_p14_t1_g44435);
      double cse_p12_t11_g44237 = (cse_p12_t137_g44363 * cse_p14_t1_g44435);
      double cse_p12_t12_g44238 = (cse_p12_t145_g44371 * cse_p14_t1_g44435);
      double cse_p12_t109_g44335 = (cse_p12_t1_g44410 + cse_p12_t46_g44272);
      double cse_p12_t110_g44336 = (cse_p12_t1_g44411 + cse_p12_t50_g44276);
      double cse_p12_t111_g44337 = (cse_p12_t1_g44412 + cse_p12_t54_g44280);
      double cse_p12_t112_g44338 = (cse_p12_t1_g44413 + cse_p12_t47_g44273);
      double cse_p12_t113_g44339 = (cse_p12_t1_g44414 + cse_p12_t43_g44269);
      double cse_p12_t114_g44340 = (cse_p12_t1_g44415 + cse_p12_t51_g44277);
      double cse_p12_t115_g44341 = (cse_p12_t44_g44270 + cse_p12_t7_g44233);
      double cse_p12_t116_g44342 = (cse_p12_t48_g44274 + cse_p12_t8_g44234);
      double cse_p12_t117_g44343 = (cse_p12_t52_g44278 + cse_p12_t9_g44235);
      double cse_p12_t118_g44344 = (cse_p12_t10_g44236 + cse_p12_t45_g44271);
      double cse_p12_t119_g44345 = (cse_p12_t11_g44237 + cse_p12_t49_g44275);
      double cse_p12_t120_g44346 = (cse_p12_t12_g44238 + cse_p12_t53_g44279);
      double cosphi = (cse_p11_t1_invsqrt44224 * dot_n1n2);
      double phi = acos(cosphi);
      double deltaphi = (phi + (-(phi0)));
      double cse_p13_t1_g44417 = (deltaphi * kdh);
      double energy = (cse_p13_t1_g44417 * deltaphi);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t3_invsqrt44226 = (1.0 / sqrt((1.0000000000000000     + (-(((cosphi) * (cosphi)))))));
      double cse_p12_t13_g44239 = (cse_p11_t3_invsqrt44226 * cse_p12_t109_g44335);
      double cse_p12_t14_g44240 = (cse_p11_t3_invsqrt44226 * cse_p12_t110_g44336);
      double cse_p12_t15_g44241 = (cse_p11_t3_invsqrt44226 * cse_p12_t111_g44337);
      double cse_p12_t16_g44242 = (cse_p11_t3_invsqrt44226 * cse_p12_t112_g44338);
      double cse_p12_t17_g44243 = (cse_p11_t3_invsqrt44226 * cse_p12_t113_g44339);
      double cse_p12_t18_g44244 = (cse_p11_t3_invsqrt44226 * cse_p12_t114_g44340);
      double cse_p12_t19_g44245 = (cse_p11_t3_invsqrt44226 * cse_p12_t115_g44341);
      double cse_p12_t20_g44246 = (cse_p11_t3_invsqrt44226 * cse_p12_t116_g44342);
      double cse_p12_t21_g44247 = (cse_p11_t3_invsqrt44226 * cse_p12_t117_g44343);
      double cse_p12_t22_g44248 = (cse_p11_t3_invsqrt44226 * cse_p12_t118_g44344);
      double cse_p12_t23_g44249 = (cse_p11_t3_invsqrt44226 * cse_p12_t119_g44345);
      double cse_p12_t24_g44250 = (cse_p11_t3_invsqrt44226 * cse_p12_t120_g44346);
      double cse_p12_t169_g44395 = (-(cse_p12_t13_g44239));
      double cse_p12_t170_g44396 = (-(cse_p12_t14_g44240));
      double cse_p12_t171_g44397 = (-(cse_p12_t15_g44241));
      double cse_p12_t172_g44398 = (-(cse_p12_t16_g44242));
      double cse_p12_t173_g44399 = (-(cse_p12_t17_g44243));
      double cse_p12_t174_g44400 = (-(cse_p12_t18_g44244));
      double cse_p12_t175_g44401 = (-(cse_p12_t19_g44245));
      double cse_p12_t176_g44402 = (-(cse_p12_t20_g44246));
      double cse_p12_t177_g44403 = (-(cse_p12_t21_g44247));
      double cse_p12_t178_g44404 = (-(cse_p12_t22_g44248));
      double cse_p12_t179_g44405 = (-(cse_p12_t23_g44249));
      double cse_p12_t180_g44406 = (-(cse_p12_t24_g44250));
      double cse_p21_t1_g44441 = (cse_p11_t3_invsqrt44226 * cse_p13_t1_g44417);
      double g_x1 = ((-((cse_p12_t120_g44346 * cse_p21_t1_g44441))) + (-((cse_p12_t120_g44346 * cse_p21_t1_g44441))));
      KernelGradientAcc(i3x1, 0, g_x1);
      double g_y1 = ((-((cse_p12_t118_g44344 * cse_p21_t1_g44441))) + (-((cse_p12_t118_g44344 * cse_p21_t1_g44441))));
      KernelGradientAcc(i3x1, 1, g_y1);
      double g_z1 = ((-((cse_p12_t119_g44345 * cse_p21_t1_g44441))) + (-((cse_p12_t119_g44345 * cse_p21_t1_g44441))));
      KernelGradientAcc(i3x1, 2, g_z1);
      double g_x2 = ((-((cse_p12_t111_g44337 * cse_p21_t1_g44441))) + (-((cse_p12_t111_g44337 * cse_p21_t1_g44441))));
      KernelGradientAcc(i3x2, 0, g_x2);
      double g_y2 = ((-((cse_p12_t109_g44335 * cse_p21_t1_g44441))) + (-((cse_p12_t109_g44335 * cse_p21_t1_g44441))));
      KernelGradientAcc(i3x2, 1, g_y2);
      double g_z2 = ((-((cse_p12_t110_g44336 * cse_p21_t1_g44441))) + (-((cse_p12_t110_g44336 * cse_p21_t1_g44441))));
      KernelGradientAcc(i3x2, 2, g_z2);
      double g_x3 = ((-((cse_p12_t114_g44340 * cse_p21_t1_g44441))) + (-((cse_p12_t114_g44340 * cse_p21_t1_g44441))));
      KernelGradientAcc(i3x3, 0, g_x3);
      double g_y3 = ((-((cse_p12_t113_g44339 * cse_p21_t1_g44441))) + (-((cse_p12_t113_g44339 * cse_p21_t1_g44441))));
      KernelGradientAcc(i3x3, 1, g_y3);
      double g_z3 = ((-((cse_p12_t112_g44338 * cse_p21_t1_g44441))) + (-((cse_p12_t112_g44338 * cse_p21_t1_g44441))));
      KernelGradientAcc(i3x3, 2, g_z3);
      double g_x4 = ((-((cse_p12_t117_g44343 * cse_p21_t1_g44441))) + (-((cse_p12_t117_g44343 * cse_p21_t1_g44441))));
      KernelGradientAcc(i3x4, 0, g_x4);
      double g_y4 = ((-((cse_p12_t115_g44341 * cse_p21_t1_g44441))) + (-((cse_p12_t115_g44341 * cse_p21_t1_g44441))));
      KernelGradientAcc(i3x4, 1, g_y4);
      double g_z4 = ((-((cse_p12_t116_g44342 * cse_p21_t1_g44441))) + (-((cse_p12_t116_g44342 * cse_p21_t1_g44441))));
      KernelGradientAcc(i3x4, 2, g_z4);
      double cse_p12_t1_g44407 = (2.0000000000000000     * kdh);
      double h_x1_x1 = (cse_p12_t180_g44406 * cse_p12_t180_g44406 * cse_p12_t1_g44407);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 0, h_x1_x1);
      double cse_p15_t1_g44440 = (cse_p12_t178_g44404 * cse_p12_t1_g44407);
      double h_x1_y1 = (cse_p12_t180_g44406 * cse_p15_t1_g44440);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 1, h_x1_y1);
      double h_x1_z1 = (cse_p12_t179_g44405 * cse_p12_t180_g44406 * cse_p12_t1_g44407);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 2, h_x1_z1);
      double cse_p14_t1_g44429 = (cse_p12_t171_g44397 * cse_p12_t1_g44407);
      double h_x1_x2 = (cse_p12_t180_g44406 * cse_p14_t1_g44429);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 0, h_x1_x2);
      double cse_p14_t1_g44427 = (cse_p12_t169_g44395 * cse_p12_t1_g44407);
      double h_x1_y2 = (cse_p12_t180_g44406 * cse_p14_t1_g44427);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 1, h_x1_y2);
      double cse_p14_t1_g44428 = (cse_p12_t170_g44396 * cse_p12_t1_g44407);
      double h_x1_z2 = (cse_p12_t180_g44406 * cse_p14_t1_g44428);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 2, h_x1_z2);
      double cse_p14_t1_g44432 = (cse_p12_t174_g44400 * cse_p12_t1_g44407);
      double h_x1_x3 = (cse_p12_t180_g44406 * cse_p14_t1_g44432);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x3, 0, h_x1_x3);
      double cse_p14_t1_g44431 = (cse_p12_t173_g44399 * cse_p12_t1_g44407);
      double h_x1_y3 = (cse_p12_t180_g44406 * cse_p14_t1_g44431);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x3, 1, h_x1_y3);
      double cse_p14_t1_g44430 = (cse_p12_t172_g44398 * cse_p12_t1_g44407);
      double h_x1_z3 = (cse_p12_t180_g44406 * cse_p14_t1_g44430);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x3, 2, h_x1_z3);
      double cse_p15_t1_g44439 = (cse_p12_t177_g44403 * cse_p12_t1_g44407);
      double h_x1_x4 = (cse_p12_t180_g44406 * cse_p15_t1_g44439);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x4, 0, h_x1_x4);
      double cse_p14_t1_g44433 = (cse_p12_t175_g44401 * cse_p12_t1_g44407);
      double h_x1_y4 = (cse_p12_t180_g44406 * cse_p14_t1_g44433);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x4, 1, h_x1_y4);
      double cse_p15_t1_g44438 = (cse_p12_t176_g44402 * cse_p12_t1_g44407);
      double h_x1_z4 = (cse_p12_t180_g44406 * cse_p15_t1_g44438);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x4, 2, h_x1_z4);
      double h_y1_y1 = (cse_p12_t178_g44404 * cse_p15_t1_g44440);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 1, h_y1_y1);
      double h_y1_z1 = (cse_p12_t179_g44405 * cse_p15_t1_g44440);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 2, h_y1_z1);
      double h_y1_x2 = (cse_p12_t178_g44404 * cse_p14_t1_g44429);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 0, h_y1_x2);
      double h_y1_y2 = (cse_p12_t178_g44404 * cse_p14_t1_g44427);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 1, h_y1_y2);
      double h_y1_z2 = (cse_p12_t178_g44404 * cse_p14_t1_g44428);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 2, h_y1_z2);
      double h_y1_x3 = (cse_p12_t178_g44404 * cse_p14_t1_g44432);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x3, 0, h_y1_x3);
      double h_y1_y3 = (cse_p12_t178_g44404 * cse_p14_t1_g44431);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x3, 1, h_y1_y3);
      double h_y1_z3 = (cse_p12_t178_g44404 * cse_p14_t1_g44430);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x3, 2, h_y1_z3);
      double h_y1_x4 = (cse_p12_t178_g44404 * cse_p15_t1_g44439);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x4, 0, h_y1_x4);
      double h_y1_y4 = (cse_p12_t178_g44404 * cse_p14_t1_g44433);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x4, 1, h_y1_y4);
      double h_y1_z4 = (cse_p12_t178_g44404 * cse_p15_t1_g44438);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x4, 2, h_y1_z4);
      double h_z1_z1 = (cse_p12_t179_g44405 * cse_p12_t179_g44405 * cse_p12_t1_g44407);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x1, 2, h_z1_z1);
      double h_z1_x2 = (cse_p12_t179_g44405 * cse_p14_t1_g44429);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 0, h_z1_x2);
      double h_z1_y2 = (cse_p12_t179_g44405 * cse_p14_t1_g44427);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 1, h_z1_y2);
      double h_z1_z2 = (cse_p12_t179_g44405 * cse_p14_t1_g44428);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 2, h_z1_z2);
      double h_z1_x3 = (cse_p12_t179_g44405 * cse_p14_t1_g44432);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x3, 0, h_z1_x3);
      double h_z1_y3 = (cse_p12_t179_g44405 * cse_p14_t1_g44431);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x3, 1, h_z1_y3);
      double h_z1_z3 = (cse_p12_t179_g44405 * cse_p14_t1_g44430);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x3, 2, h_z1_z3);
      double h_z1_x4 = (cse_p12_t179_g44405 * cse_p15_t1_g44439);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x4, 0, h_z1_x4);
      double h_z1_y4 = (cse_p12_t179_g44405 * cse_p14_t1_g44433);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x4, 1, h_z1_y4);
      double h_z1_z4 = (cse_p12_t179_g44405 * cse_p15_t1_g44438);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x4, 2, h_z1_z4);
      double h_x2_x2 = (cse_p12_t171_g44397 * cse_p14_t1_g44429);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 0, h_x2_x2);
      double h_x2_y2 = (cse_p12_t171_g44397 * cse_p14_t1_g44427);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 1, h_x2_y2);
      double h_x2_z2 = (cse_p12_t171_g44397 * cse_p14_t1_g44428);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 2, h_x2_z2);
      double h_x2_x3 = (cse_p12_t174_g44400 * cse_p14_t1_g44429);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x3, 0, h_x2_x3);
      double h_x2_y3 = (cse_p12_t173_g44399 * cse_p14_t1_g44429);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x3, 1, h_x2_y3);
      double h_x2_z3 = (cse_p12_t172_g44398 * cse_p14_t1_g44429);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x3, 2, h_x2_z3);
      double h_x2_x4 = (cse_p12_t177_g44403 * cse_p14_t1_g44429);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x4, 0, h_x2_x4);
      double h_x2_y4 = (cse_p12_t175_g44401 * cse_p14_t1_g44429);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x4, 1, h_x2_y4);
      double h_x2_z4 = (cse_p12_t176_g44402 * cse_p14_t1_g44429);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x4, 2, h_x2_z4);
      double h_y2_y2 = (cse_p12_t169_g44395 * cse_p14_t1_g44427);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 1, h_y2_y2);
      double h_y2_z2 = (cse_p12_t170_g44396 * cse_p14_t1_g44427);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 2, h_y2_z2);
      double h_y2_x3 = (cse_p12_t174_g44400 * cse_p14_t1_g44427);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x3, 0, h_y2_x3);
      double h_y2_y3 = (cse_p12_t173_g44399 * cse_p14_t1_g44427);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x3, 1, h_y2_y3);
      double h_y2_z3 = (cse_p12_t172_g44398 * cse_p14_t1_g44427);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x3, 2, h_y2_z3);
      double h_y2_x4 = (cse_p12_t177_g44403 * cse_p14_t1_g44427);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x4, 0, h_y2_x4);
      double h_y2_y4 = (cse_p12_t175_g44401 * cse_p14_t1_g44427);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x4, 1, h_y2_y4);
      double h_y2_z4 = (cse_p12_t176_g44402 * cse_p14_t1_g44427);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x4, 2, h_y2_z4);
      double h_z2_z2 = (cse_p12_t170_g44396 * cse_p14_t1_g44428);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x2, 2, h_z2_z2);
      double h_z2_x3 = (cse_p12_t174_g44400 * cse_p14_t1_g44428);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x3, 0, h_z2_x3);
      double h_z2_y3 = (cse_p12_t173_g44399 * cse_p14_t1_g44428);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x3, 1, h_z2_y3);
      double h_z2_z3 = (cse_p12_t172_g44398 * cse_p14_t1_g44428);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x3, 2, h_z2_z3);
      double h_z2_x4 = (cse_p12_t177_g44403 * cse_p14_t1_g44428);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x4, 0, h_z2_x4);
      double h_z2_y4 = (cse_p12_t175_g44401 * cse_p14_t1_g44428);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x4, 1, h_z2_y4);
      double h_z2_z4 = (cse_p12_t176_g44402 * cse_p14_t1_g44428);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x4, 2, h_z2_z4);
      double h_x3_x3 = (cse_p12_t174_g44400 * cse_p14_t1_g44432);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 0, i3x3, 0, h_x3_x3);
      double h_x3_y3 = (cse_p12_t174_g44400 * cse_p14_t1_g44431);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 0, i3x3, 1, h_x3_y3);
      double h_x3_z3 = (cse_p12_t174_g44400 * cse_p14_t1_g44430);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 0, i3x3, 2, h_x3_z3);
      double h_x3_x4 = (cse_p12_t177_g44403 * cse_p14_t1_g44432);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 0, i3x4, 0, h_x3_x4);
      double h_x3_y4 = (cse_p12_t175_g44401 * cse_p14_t1_g44432);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 0, i3x4, 1, h_x3_y4);
      double h_x3_z4 = (cse_p12_t176_g44402 * cse_p14_t1_g44432);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 0, i3x4, 2, h_x3_z4);
      double h_y3_y3 = (cse_p12_t173_g44399 * cse_p14_t1_g44431);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 1, i3x3, 1, h_y3_y3);
      double h_y3_z3 = (cse_p12_t173_g44399 * cse_p14_t1_g44430);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 1, i3x3, 2, h_y3_z3);
      double h_y3_x4 = (cse_p12_t177_g44403 * cse_p14_t1_g44431);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 1, i3x4, 0, h_y3_x4);
      double h_y3_y4 = (cse_p12_t175_g44401 * cse_p14_t1_g44431);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 1, i3x4, 1, h_y3_y4);
      double h_y3_z4 = (cse_p12_t176_g44402 * cse_p14_t1_g44431);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 1, i3x4, 2, h_y3_z4);
      double h_z3_z3 = (cse_p12_t172_g44398 * cse_p14_t1_g44430);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 2, i3x3, 2, h_z3_z3);
      double h_z3_x4 = (cse_p12_t177_g44403 * cse_p14_t1_g44430);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 2, i3x4, 0, h_z3_x4);
      double h_z3_y4 = (cse_p12_t175_g44401 * cse_p14_t1_g44430);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 2, i3x4, 1, h_z3_y4);
      double h_z3_z4 = (cse_p12_t176_g44402 * cse_p14_t1_g44430);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 2, i3x4, 2, h_z3_z4);
      double h_x4_x4 = (cse_p12_t177_g44403 * cse_p15_t1_g44439);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x4, 0, i3x4, 0, h_x4_x4);
      double h_x4_y4 = (cse_p12_t177_g44403 * cse_p14_t1_g44433);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x4, 0, i3x4, 1, h_x4_y4);
      double h_x4_z4 = (cse_p12_t177_g44403 * cse_p15_t1_g44438);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x4, 0, i3x4, 2, h_x4_z4);
      double h_y4_y4 = (cse_p12_t175_g44401 * cse_p14_t1_g44433);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x4, 1, i3x4, 1, h_y4_y4);
      double h_y4_z4 = (cse_p12_t176_g44402 * cse_p14_t1_g44433);
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x4, 1, i3x4, 2, h_y4_z4);
      double h_z4_z4 = (cse_p12_t176_g44402 * cse_p15_t1_g44438);
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x4, 2, i3x4, 2, h_z4_z4);
    }
  }
  return energy_added;
}
void hessian_fd(double kdh, double phi0, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 12;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  const double invh2 = 1.0/((h*h));
  double e0 = 0.0;
  energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
  {
    double saved = position[i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 0] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 0, d);
  }
  {
    double saved = position[i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 1] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 1, d);
  }
  {
    double saved = position[i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 2] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 2, d);
  }
  {
    double saved = position[i3x4 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 0] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 0, d);
  }
  {
    double saved = position[i3x4 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 1] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 1, d);
  }
  {
    double saved = position[i3x4 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 2] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 2, d);
  }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 0, i3x1, 0, hval);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 1, i3x1, 1, hval);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 2, i3x1, 2, hval);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 0, i3x2, 0, hval);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 1, i3x2, 1, hval);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x2, 2, hval);
  }
  {
    double saved = position[i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 0] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 0, i3x3, 0, hval);
  }
  {
    double saved = position[i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 1] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x3, 1, hval);
  }
  {
    double saved = position[i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 2] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x3, 2, hval);
  }
  {
    double saved = position[i3x4 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 0] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 0, i3x4, 0, hval);
  }
  {
    double saved = position[i3x4 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 1] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 1, i3x4, 1, hval);
  }
  {
    double saved = position[i3x4 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 2] = saved + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x3 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x3 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x3 + 2] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x3 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x3 + 2] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x4 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x3 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x3 + 2] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x4 + 0] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x4 + 1] = saved_j + h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_j - h;
    energy(kdh, phi0, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i; position[i3x4 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 2, i3x4, 1, hval);
  }
}

};
