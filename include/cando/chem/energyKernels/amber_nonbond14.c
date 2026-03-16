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
      double cse_p11_t1_g185177 = (-3.00000000000000000e+0 * term.c);
      double dx = (x1 + (-(x2)));
      double dy = (y1 + (-(y2)));
      double dz = (z1 + (-(z2)));
      double r2 = ((dx * dx) + (dy * dy) + (dz * dz));
      double cse_p1_t1_sqrt185173 = sqrt(r2);
      double cse_p1_t2_invsqrt185174 = (1.0 / (cse_p1_t1_sqrt185173));
      double r = cse_p1_t1_sqrt185173;
      double cse_p1_t3_invr185175 = (1.0 / (r));
      double cse_p1_t4_invr2185176 = (cse_p1_t3_invr185175 * cse_p1_t3_invr185175);
      double invr = cse_p1_t3_invr185175;
      double invr2 = (invr * invr);
      double cse_p12_t1_g185188 = (invr2 * invr2);
      double invr6 = (cse_p12_t1_g185188 * invr2);
      double cse_p11_t1_g185187 = (term.a * invr6);
      double cse_p11_t2_g185178 = (6.00000000000000000e+0 * cse_p11_t1_g185187);
      double cse_p11_t3_g185179 = (cse_p11_t1_g185177 + cse_p11_t2_g185178);
      double e_lj = ((cse_p11_t1_g185187 * invr6) + (-((term.c * invr6))));
      double e_coul = (invr * term.kqq);
      double energy = (e_coul + e_lj);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t1_g185180 = (cse_p1_t2_invsqrt185174 * cse_p1_t4_invr2185176);
      double cse_p11_t1_g185181 = (cse_p11_t1_g185180 * cse_p11_t3_g185179 * cse_p12_t1_g185188 * invr);
      double cse_p11_t1_g185182 = (cse_p11_t1_g185180 * term.kqq);
      double cse_p11_t1_g185183 = (-(cse_p11_t1_g185182));
      double cse_p11_t1_g185184 = (-2.00000000000000000e+0 * cse_p11_t1_g185181);
      double g_x1 = ((cse_p11_t1_g185183 * dx) + (cse_p11_t1_g185184 * dx));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = ((cse_p11_t1_g185183 * dy) + (cse_p11_t1_g185184 * dy));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = ((cse_p11_t1_g185183 * dz) + (cse_p11_t1_g185184 * dz));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p11_t1_g185186 = (2.00000000000000000e+0 * cse_p11_t1_g185181);
      double g_x2 = ((cse_p11_t1_g185182 * dx) + (cse_p11_t1_g185186 * dx));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = ((cse_p11_t1_g185182 * dy) + (cse_p11_t1_g185186 * dy));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = ((cse_p11_t1_g185182 * dz) + (cse_p11_t1_g185186 * dz));
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
      double cse_p11_t5_g185197 = (-3.00000000000000000e+0 * term.c);
      double dx = (x1 + (-(x2)));
      double cse_p11_t11_g185203 = (dx * dx);
      double cse_p11_t28_g185220 = (-(cse_p11_t11_g185203));
      double cse_p11_t33_g185225 = (-(dx));
      double dy = (y1 + (-(y2)));
      double cse_p11_t15_g185207 = (dy * dy);
      double cse_p11_t29_g185221 = (-(cse_p11_t15_g185207));
      double cse_p11_t34_g185226 = (-(dy));
      double dz = (z1 + (-(z2)));
      double cse_p11_t19_g185211 = (dz * dz);
      double cse_p11_t30_g185222 = (-(cse_p11_t19_g185211));
      double cse_p11_t35_g185227 = (-(dz));
      double r2 = (cse_p11_t11_g185203 + cse_p11_t15_g185207 + cse_p11_t19_g185211);
      double cse_p1_t1_sqrt185189 = sqrt(r2);
      double cse_p1_t2_invsqrt185190 = (1.0 / (cse_p1_t1_sqrt185189));
      double r = cse_p1_t1_sqrt185189;
      double cse_p62_t1_g185251 = (r * r);
      double cse_p61_t1_invr185245 = (1.0 / (((r) * (r))));
      double cse_p11_t1_g185244 = (cse_p61_t1_invr185245 * term.kqq);
      double cse_p61_t4_invr185248 = (1.0 / ((r * r * r)));
      double cse_p61_t3_invr185247 = (1.0 / (pow(r, 7)));
      double cse_p166_t1_g185265 = (6.00000000000000000e+0 * term.c * cse_p61_t3_invr185247);
      double cse_p61_t2_invr185246 = (1.0 / (pow(r, 13)));
      double cse_p114_t1_g185264 = (-1.20000000000000000e+1 * term.a * cse_p61_t2_invr185246);
      double cse_p11_t1_g185240 = (cse_p61_t4_invr185248 * (cse_p114_t1_g185264 + cse_p166_t1_g185265 + (-(cse_p11_t1_g185244))));
      double cse_p11_t1_g185193 = (cse_p11_t1_g185240 * (cse_p11_t28_g185220 + cse_p62_t1_g185251));
      double cse_p11_t2_g185194 = (cse_p11_t1_g185240 * (cse_p11_t29_g185221 + cse_p62_t1_g185251));
      double cse_p11_t3_g185195 = (cse_p11_t1_g185240 * (cse_p11_t30_g185222 + cse_p62_t1_g185251));
      double cse_p61_t5_invr185249 = (1.0 / (pow(r, 8)));
      double cse_p11_t6_g185198 = (-4.20000000000000040e+1 * term.c * cse_p61_t5_invr185249);
      double cse_p61_t6_invr185250 = (1.0 / (pow(r, 14)));
      double cse_p11_t7_g185199 = (1.56000000000000000e+2 * term.a * cse_p61_t6_invr185250);
      double cse_p11_t8_g185200 = (2.00000000000000000e+0 * cse_p61_t4_invr185248 * term.kqq);
      double cse_p11_t12_g185204 = (cse_p11_t1_g185240 * dx * dy);
      double cse_p63_t1_g185257 = (cse_p11_t1_g185240 * dz);
      double cse_p11_t13_g185205 = (cse_p63_t1_g185257 * dx);
      double cse_p11_t14_g185206 = (cse_p11_t1_g185240 * cse_p11_t33_g185225 * dy);
      double cse_p11_t16_g185208 = (cse_p63_t1_g185257 * dy);
      double cse_p11_t17_g185209 = (cse_p11_t33_g185225 * cse_p63_t1_g185257);
      double cse_p11_t18_g185210 = (cse_p11_t34_g185226 * cse_p63_t1_g185257);
      double cse_p11_t22_g185214 = (cse_p114_t1_g185264 + cse_p166_t1_g185265 + (-(cse_p11_t1_g185244)));
      double cse_p11_t24_g185216 = (cse_p11_t6_g185198 + cse_p11_t7_g185199 + cse_p11_t8_g185200);
      double cse_p11_t25_g185217 = (cse_p11_t28_g185220 + cse_p62_t1_g185251);
      double cse_p11_t26_g185218 = (cse_p11_t29_g185221 + cse_p62_t1_g185251);
      double cse_p11_t27_g185219 = (cse_p11_t30_g185222 + cse_p62_t1_g185251);
      double cse_p11_t31_g185223 = (-(cse_p11_t1_g185244));
      double cse_p11_t32_g185224 = (-(cse_p62_t1_g185251));
      double cse_p11_t42_g185234 = pow(r, 13);
      double cse_p11_t43_g185235 = pow(r, 14);
      double cse_p11_t44_g185236 = ((r) * (r));
      double cse_p11_t45_g185237 = (r * r * r);
      double cse_p11_t46_g185238 = pow(r, 7);
      double cse_p11_t47_g185239 = pow(r, 8);
      double cse_p1_t3_invr185191 = (1.0 / (r));
      double cse_p1_t4_invr2185192 = (cse_p1_t3_invr185191 * cse_p1_t3_invr185191);
      double invr = cse_p1_t3_invr185191;
      double invr2 = (invr * invr);
      double invr6 = (invr2 * invr2 * invr2);
      double cse_p11_t9_g185201 = (6.00000000000000000e+0 * term.a * invr6);
      double cse_p11_t23_g185215 = (cse_p11_t5_g185197 + cse_p11_t9_g185201);
      double e_lj = ((term.a * invr6 * invr6) + (-((term.c * invr6))));
      double e_coul = (invr * term.kqq);
      double energy = (e_coul + e_lj);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t1_g185242 = (cse_p1_t2_invsqrt185190 * cse_p1_t4_invr2185192);
      double cse_p63_t1_g185262 = (cse_p11_t1_g185242 * term.kqq);
      double cse_p63_t1_g185252 = (cse_p63_t1_g185262 * dx);
      double cse_p63_t2_g185253 = (cse_p63_t1_g185262 * dy);
      double cse_p63_t3_g185254 = (cse_p63_t1_g185262 * dz);
      double cse_p11_t1_g185243 = (cse_p11_t1_g185242 * cse_p11_t23_g185215 * invr2 * invr2 * invr);
      double cse_p63_t1_g185259 = (-2.00000000000000000e+0 * cse_p11_t1_g185243);
      double g_x1 = ((cse_p63_t1_g185259 * dx) + (-(cse_p63_t1_g185252)));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = ((cse_p63_t1_g185259 * dy) + (-(cse_p63_t2_g185253)));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = ((cse_p63_t1_g185259 * dz) + (-(cse_p63_t3_g185254)));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p63_t1_g185260 = (2.00000000000000000e+0 * cse_p11_t1_g185243);
      double g_x2 = (cse_p63_t1_g185252 + (cse_p63_t1_g185260 * dx));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (cse_p63_t2_g185253 + (cse_p63_t1_g185260 * dy));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (cse_p63_t3_g185254 + (cse_p63_t1_g185260 * dz));
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double cse_p11_t1_g185241 = (cse_p11_t24_g185216 * cse_p1_t3_invr185191 * cse_p1_t3_invr185191);
      double cse_p63_t1_g185261 = (cse_p11_t1_g185241 * dx);
      double h_x1_x1 = (cse_p11_t1_g185193 + (cse_p63_t1_g185261 * dx));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
      double h_x1_y1 = (cse_p11_t14_g185206 + (cse_p63_t1_g185261 * dy));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
      double h_x1_z1 = (cse_p11_t17_g185209 + (cse_p63_t1_g185261 * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
      double cse_p63_t1_g185255 = (cse_p11_t1_g185241 * cse_p11_t33_g185225);
      double cse_p63_t1_g185263 = (cse_p11_t22_g185214 * cse_p61_t4_invr185248);
      double h_x1_x2 = ((cse_p63_t1_g185255 * dx) + (cse_p63_t1_g185263 * (cse_p11_t11_g185203 + cse_p11_t32_g185224)));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
      double cse_p63_t1_g185256 = (cse_p11_t1_g185241 * cse_p11_t34_g185226);
      double h_x1_y2 = (cse_p11_t12_g185204 + (cse_p63_t1_g185256 * dx));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
      double cse_p63_t1_g185258 = (cse_p11_t1_g185241 * cse_p11_t35_g185227);
      double h_x1_z2 = (cse_p11_t13_g185205 + (cse_p63_t1_g185258 * dx));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
      double h_y1_y1 = (cse_p11_t2_g185194 + (cse_p11_t1_g185241 * dy * dy));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
      double h_y1_z1 = (cse_p11_t18_g185210 + (cse_p11_t1_g185241 * dy * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
      double h_y1_x2 = (cse_p11_t12_g185204 + (cse_p63_t1_g185255 * dy));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
      double h_y1_y2 = ((cse_p63_t1_g185256 * dy) + (cse_p63_t1_g185263 * (cse_p11_t15_g185207 + cse_p11_t32_g185224)));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
      double h_y1_z2 = (cse_p11_t16_g185208 + (cse_p63_t1_g185258 * dy));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
      double h_z1_z1 = (cse_p11_t3_g185195 + (cse_p11_t1_g185241 * dz * dz));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
      double h_z1_x2 = (cse_p11_t13_g185205 + (cse_p63_t1_g185255 * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
      double h_z1_y2 = (cse_p11_t16_g185208 + (cse_p63_t1_g185256 * dz));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
      double h_z1_z2 = ((cse_p63_t1_g185258 * dz) + (cse_p63_t1_g185263 * (cse_p11_t19_g185211 + cse_p11_t32_g185224)));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
      double h_x2_x2 = (cse_p11_t1_g185193 + (cse_p11_t33_g185225 * cse_p63_t1_g185255));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
      double h_x2_y2 = (cse_p11_t14_g185206 + (cse_p11_t34_g185226 * cse_p63_t1_g185255));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
      double h_x2_z2 = (cse_p11_t17_g185209 + (cse_p11_t35_g185227 * cse_p63_t1_g185255));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
      double h_y2_y2 = (cse_p11_t2_g185194 + (cse_p11_t34_g185226 * cse_p63_t1_g185256));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
      double h_y2_z2 = (cse_p11_t18_g185210 + (cse_p11_t35_g185227 * cse_p63_t1_g185256));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
      double h_z2_z2 = (cse_p11_t3_g185195 + (cse_p11_t35_g185227 * cse_p63_t1_g185258));
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
