#include "amber_nonbond14.h"

template <typename HESSIAN>
struct Amber_Nonbond14 {
  static constexpr size_t PositionSize = 6;
  static std::string description() { return "mathkernel-amber_nonbond14"; };
double gradient(const AmberNonbond14Term& term, double* position, double* energy_accumulate, double* force) {
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
      double cse_p11_t1_g340191 = (-3.00000000000000000e+0 * term.c);
      double dx = (x1 + (-(x2)));
      double dy = (y1 + (-(y2)));
      double dz = (z1 + (-(z2)));
      double r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      double cse_p1_t1_sqrt340187 = sqrt(r2);
      double cse_p1_t2_invsqrt340188 = (1.0 / (cse_p1_t1_sqrt340187));
      double r = cse_p1_t1_sqrt340187;
      double cse_p1_t3_invr340189 = (1.0 / (r));
      double cse_p1_t4_invr2340190 = (cse_p1_t3_invr340189 * cse_p1_t3_invr340189);
      double invr = cse_p1_t3_invr340189;
      double invr2 = (invr * invr);
      double cse_p12_t1_g340202 = (invr2 * invr2);
      double invr6 = (cse_p12_t1_g340202 * invr2);
      double cse_p11_t1_g340201 = (term.a * invr6);
      double cse_p11_t2_g340192 = (6.00000000000000000e+0 * cse_p11_t1_g340201);
      double cse_p11_t3_g340193 = (cse_p11_t1_g340191 + cse_p11_t2_g340192);
      double e_lj = ((cse_p11_t1_g340201 * invr6) + (-((term.c * invr6))));
      double e_coul = (invr * term.kqq);
      double energy = (e_coul + e_lj);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t1_g340194 = (cse_p1_t2_invsqrt340188 * cse_p1_t4_invr2340190);
      double cse_p11_t1_g340195 = (cse_p11_t1_g340194 * cse_p11_t3_g340193 * cse_p12_t1_g340202 * invr);
      double cse_p11_t1_g340196 = (cse_p11_t1_g340194 * term.kqq);
      double cse_p11_t1_g340197 = (-(cse_p11_t1_g340196));
      double cse_p11_t1_g340198 = (-2.00000000000000000e+0 * cse_p11_t1_g340195);
      double g_x1 = ((cse_p11_t1_g340197 * dx) + (cse_p11_t1_g340198 * dx));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = ((cse_p11_t1_g340197 * dy) + (cse_p11_t1_g340198 * dy));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = ((cse_p11_t1_g340197 * dz) + (cse_p11_t1_g340198 * dz));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p11_t1_g340200 = (2.00000000000000000e+0 * cse_p11_t1_g340195);
      double g_x2 = ((cse_p11_t1_g340196 * dx) + (cse_p11_t1_g340200 * dx));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = ((cse_p11_t1_g340196 * dy) + (cse_p11_t1_g340200 * dy));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = ((cse_p11_t1_g340196 * dz) + (cse_p11_t1_g340200 * dz));
      KernelGradientAcc(term.i3x2, 2, g_z2);
    }
  }
  return energy_added;
}
void gradient_fd(const AmberNonbond14Term& term, double* position, double* energy_accumulate, double* force)
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

double energy(const AmberNonbond14Term& term, double* position, double* energy_accumulate) {
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
      double dx = (x1 + (-(x2)));
      double dy = (y1 + (-(y2)));
      double dz = (z1 + (-(z2)));
      double r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      double r = sqrt(r2);
      double invr = (1.0 / (r));
      double invr2 = (invr * invr);
      double invr6 = (invr2 * invr2 * invr2);
      double e_lj = ((term.a * invr6 * invr6) + (-((term.c * invr6))));
      double e_coul = (invr * term.kqq);
      double energy = (e_coul + e_lj);
      energy_added += energy;
      *energy_accumulate += energy;
    }
  }
  return energy_added;
}
void energy_fd(const AmberNonbond14Term& term, double* position, double* energy_accumulate)
{
  energy(term, position, energy_accumulate);
}

double hessian(const AmberNonbond14Term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
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
      double cse_p11_t5_g340211 = (-3.00000000000000000e+0 * term.c);
      double dx = (x1 + (-(x2)));
      double cse_p11_t11_g340217 = (dx * dx);
      double cse_p11_t28_g340234 = (-(cse_p11_t11_g340217));
      double cse_p11_t33_g340239 = (-(dx));
      double dy = (y1 + (-(y2)));
      double cse_p11_t15_g340221 = (dy * dy);
      double cse_p11_t29_g340235 = (-(cse_p11_t15_g340221));
      double cse_p11_t34_g340240 = (-(dy));
      double dz = (z1 + (-(z2)));
      double cse_p11_t19_g340225 = (dz * dz);
      double cse_p11_t30_g340236 = (-(cse_p11_t19_g340225));
      double cse_p11_t35_g340241 = (-(dz));
      double r2 = (cse_p11_t11_g340217 + cse_p11_t15_g340221 + cse_p11_t19_g340225);
      double cse_p1_t1_sqrt340203 = sqrt(r2);
      double cse_p1_t2_invsqrt340204 = (1.0 / (cse_p1_t1_sqrt340203));
      double r = cse_p1_t1_sqrt340203;
      double cse_p62_t1_g340265 = (r * r);
      double cse_p61_t1_invr340259 = (1.0 / (((r) * (r))));
      double cse_p11_t1_g340258 = (cse_p61_t1_invr340259 * term.kqq);
      double cse_p61_t4_invr340262 = (1.0 / ((r * r * r)));
      double cse_p61_t3_invr340261 = (1.0 / (pow(r, 7)));
      double cse_p166_t1_g340279 = (6.00000000000000000e+0 * term.c * cse_p61_t3_invr340261);
      double cse_p61_t2_invr340260 = (1.0 / (pow(r, 13)));
      double cse_p114_t1_g340278 = (-1.20000000000000000e+1 * term.a * cse_p61_t2_invr340260);
      double cse_p11_t1_g340254 = (cse_p61_t4_invr340262 * (cse_p114_t1_g340278 + cse_p166_t1_g340279 + (-(cse_p11_t1_g340258))));
      double cse_p11_t1_g340207 = (cse_p11_t1_g340254 * (cse_p11_t28_g340234 + cse_p62_t1_g340265));
      double cse_p11_t2_g340208 = (cse_p11_t1_g340254 * (cse_p11_t29_g340235 + cse_p62_t1_g340265));
      double cse_p11_t3_g340209 = (cse_p11_t1_g340254 * (cse_p11_t30_g340236 + cse_p62_t1_g340265));
      double cse_p61_t5_invr340263 = (1.0 / (pow(r, 8)));
      double cse_p11_t6_g340212 = (-4.20000000000000040e+1 * term.c * cse_p61_t5_invr340263);
      double cse_p61_t6_invr340264 = (1.0 / (pow(r, 14)));
      double cse_p11_t7_g340213 = (1.56000000000000000e+2 * term.a * cse_p61_t6_invr340264);
      double cse_p11_t8_g340214 = (2.00000000000000000e+0 * cse_p61_t4_invr340262 * term.kqq);
      double cse_p11_t12_g340218 = (cse_p11_t1_g340254 * dx * dy);
      double cse_p63_t1_g340271 = (cse_p11_t1_g340254 * dz);
      double cse_p11_t13_g340219 = (cse_p63_t1_g340271 * dx);
      double cse_p11_t14_g340220 = (cse_p11_t1_g340254 * cse_p11_t33_g340239 * dy);
      double cse_p11_t16_g340222 = (cse_p63_t1_g340271 * dy);
      double cse_p11_t17_g340223 = (cse_p11_t33_g340239 * cse_p63_t1_g340271);
      double cse_p11_t18_g340224 = (cse_p11_t34_g340240 * cse_p63_t1_g340271);
      double cse_p11_t22_g340228 = (cse_p114_t1_g340278 + cse_p166_t1_g340279 + (-(cse_p11_t1_g340258)));
      double cse_p11_t24_g340230 = (cse_p11_t6_g340212 + cse_p11_t7_g340213 + cse_p11_t8_g340214);
      double cse_p11_t25_g340231 = (cse_p11_t28_g340234 + cse_p62_t1_g340265);
      double cse_p11_t26_g340232 = (cse_p11_t29_g340235 + cse_p62_t1_g340265);
      double cse_p11_t27_g340233 = (cse_p11_t30_g340236 + cse_p62_t1_g340265);
      double cse_p11_t31_g340237 = (-(cse_p11_t1_g340258));
      double cse_p11_t32_g340238 = (-(cse_p62_t1_g340265));
      double cse_p11_t42_g340248 = pow(r, 13);
      double cse_p11_t43_g340249 = pow(r, 14);
      double cse_p11_t44_g340250 = ((r) * (r));
      double cse_p11_t45_g340251 = (r * r * r);
      double cse_p11_t46_g340252 = pow(r, 7);
      double cse_p11_t47_g340253 = pow(r, 8);
      double cse_p1_t3_invr340205 = (1.0 / (r));
      double cse_p1_t4_invr2340206 = (cse_p1_t3_invr340205 * cse_p1_t3_invr340205);
      double invr = cse_p1_t3_invr340205;
      double invr2 = (invr * invr);
      double invr6 = (invr2 * invr2 * invr2);
      double cse_p11_t9_g340215 = (6.00000000000000000e+0 * term.a * invr6);
      double cse_p11_t23_g340229 = (cse_p11_t5_g340211 + cse_p11_t9_g340215);
      double e_lj = ((term.a * invr6 * invr6) + (-((term.c * invr6))));
      double e_coul = (invr * term.kqq);
      double energy = (e_coul + e_lj);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t1_g340256 = (cse_p1_t2_invsqrt340204 * cse_p1_t4_invr2340206);
      double cse_p63_t1_g340276 = (cse_p11_t1_g340256 * term.kqq);
      double cse_p63_t1_g340266 = (cse_p63_t1_g340276 * dx);
      double cse_p63_t2_g340267 = (cse_p63_t1_g340276 * dy);
      double cse_p63_t3_g340268 = (cse_p63_t1_g340276 * dz);
      double cse_p11_t1_g340257 = (cse_p11_t1_g340256 * cse_p11_t23_g340229 * invr2 * invr2 * invr);
      double cse_p63_t1_g340273 = (-2.00000000000000000e+0 * cse_p11_t1_g340257);
      double g_x1 = ((cse_p63_t1_g340273 * dx) + (-(cse_p63_t1_g340266)));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = ((cse_p63_t1_g340273 * dy) + (-(cse_p63_t2_g340267)));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = ((cse_p63_t1_g340273 * dz) + (-(cse_p63_t3_g340268)));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p63_t1_g340274 = (2.00000000000000000e+0 * cse_p11_t1_g340257);
      double g_x2 = (cse_p63_t1_g340266 + (cse_p63_t1_g340274 * dx));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (cse_p63_t2_g340267 + (cse_p63_t1_g340274 * dy));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (cse_p63_t3_g340268 + (cse_p63_t1_g340274 * dz));
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double cse_p11_t1_g340255 = (cse_p11_t24_g340230 * cse_p1_t3_invr340205 * cse_p1_t3_invr340205);
      double cse_p63_t1_g340275 = (cse_p11_t1_g340255 * dx);
      double h_x1_x1 = (cse_p11_t1_g340207 + (cse_p63_t1_g340275 * dx));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
      double h_x1_y1 = (cse_p11_t14_g340220 + (cse_p63_t1_g340275 * dy));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
      double h_x1_z1 = (cse_p11_t17_g340223 + (cse_p63_t1_g340275 * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
      double cse_p63_t1_g340269 = (cse_p11_t1_g340255 * cse_p11_t33_g340239);
      double cse_p63_t1_g340277 = (cse_p11_t22_g340228 * cse_p61_t4_invr340262);
      double h_x1_x2 = ((cse_p63_t1_g340269 * dx) + (cse_p63_t1_g340277 * (cse_p11_t11_g340217 + cse_p11_t32_g340238)));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
      double cse_p63_t1_g340270 = (cse_p11_t1_g340255 * cse_p11_t34_g340240);
      double h_x1_y2 = (cse_p11_t12_g340218 + (cse_p63_t1_g340270 * dx));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
      double cse_p63_t1_g340272 = (cse_p11_t1_g340255 * cse_p11_t35_g340241);
      double h_x1_z2 = (cse_p11_t13_g340219 + (cse_p63_t1_g340272 * dx));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
      double h_y1_y1 = (cse_p11_t2_g340208 + (cse_p11_t1_g340255 * dy * dy));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
      double h_y1_z1 = (cse_p11_t18_g340224 + (cse_p11_t1_g340255 * dy * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
      double h_y1_x2 = (cse_p11_t12_g340218 + (cse_p63_t1_g340269 * dy));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
      double h_y1_y2 = ((cse_p63_t1_g340270 * dy) + (cse_p63_t1_g340277 * (cse_p11_t15_g340221 + cse_p11_t32_g340238)));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
      double h_y1_z2 = (cse_p11_t16_g340222 + (cse_p63_t1_g340272 * dy));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
      double h_z1_z1 = (cse_p11_t3_g340209 + (cse_p11_t1_g340255 * dz * dz));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
      double h_z1_x2 = (cse_p11_t13_g340219 + (cse_p63_t1_g340269 * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
      double h_z1_y2 = (cse_p11_t16_g340222 + (cse_p63_t1_g340270 * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
      double h_z1_z2 = ((cse_p63_t1_g340272 * dz) + (cse_p63_t1_g340277 * (cse_p11_t19_g340225 + cse_p11_t32_g340238)));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
      double h_x2_x2 = (cse_p11_t1_g340207 + (cse_p11_t33_g340239 * cse_p63_t1_g340269));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
      double h_x2_y2 = (cse_p11_t14_g340220 + (cse_p11_t34_g340240 * cse_p63_t1_g340269));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
      double h_x2_z2 = (cse_p11_t17_g340223 + (cse_p11_t35_g340241 * cse_p63_t1_g340269));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
      double h_y2_y2 = (cse_p11_t2_g340208 + (cse_p11_t34_g340240 * cse_p63_t1_g340270));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
      double h_y2_z2 = (cse_p11_t18_g340224 + (cse_p11_t35_g340241 * cse_p63_t1_g340270));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
      double h_z2_z2 = (cse_p11_t3_g340209 + (cse_p11_t35_g340241 * cse_p63_t1_g340272));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
    }
  }
  return energy_added;
}
void hessian_fd(const AmberNonbond14Term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
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

};
