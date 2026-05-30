#include "dihedral_fast.h"

template <typename HESSIAN>
struct Dihedral_Fast {
  static constexpr size_t PositionSize = 12;
  static std::string description() { return "mathkernel-dihedral_fast"; };
double energy(const dihedral_term& term, double* position, double* energy_accumulate) {
  double energy_added = 0.00000000000000000e+0;
  {
    /* !BASE */
    DOUBLE x1 = position[term.i3x1 + 0];
    DOUBLE y1 = position[term.i3x1 + 1];
    DOUBLE z1 = position[term.i3x1 + 2];
    DOUBLE x2 = position[term.i3x2 + 0];
    DOUBLE y2 = position[term.i3x2 + 1];
    DOUBLE z2 = position[term.i3x2 + 2];
    DOUBLE x3 = position[term.i3x3 + 0];
    DOUBLE y3 = position[term.i3x3 + 1];
    DOUBLE z3 = position[term.i3x3 + 2];
    DOUBLE x4 = position[term.i3x4 + 0];
    DOUBLE y4 = position[term.i3x4 + 1];
    DOUBLE z4 = position[term.i3x4 + 2];
    {
      /* !BASE */
      double v1x = (x2 + (-(x1)));
      double v1y = (y2 + (-(y1)));
      double v1z = (z2 + (-(z1)));
      double v2x = (x3 + (-(x2)));
      double v2y = (y3 + (-(y2)));
      double v2z = (z3 + (-(z2)));
      double v3x = (x4 + (-(x3)));
      double v3y = (y4 + (-(y3)));
      double v3z = (z4 + (-(z3)));
      double c1x = ((v1y * v2z) + (-((v1z * v2y))));
      double c1y = ((v1z * v2x) + (-((v1x * v2z))));
      double c1z = ((v1x * v2y) + (-((v1y * v2x))));
      double c2x = ((v2y * v3z) + (-((v2z * v3y))));
      double c2y = ((v2z * v3x) + (-((v2x * v3z))));
      double c2z = ((v2x * v3y) + (-((v2y * v3x))));
      double v2_sq = ((v2x * v2x) + (v2y * v2y) + (v2z * v2z));
      double v2_len = sqrt(v2_sq);
      double t1 = (v2_len * ((c2x * v1x) + (c2y * v1y) + (c2z * v1z)));
      double t2 = ((c1x * c2x) + (c1y * c2y) + (c1z * c2z));
      double phi = atan2(t1, t2);
      double nphi = (term.n * phi);
      double sin_nphi = sin(nphi);
      double cos_nphi = cos(nphi);
      double cos_angle = ((term.cosphase * cos_nphi) + (term.sinphase * sin_nphi));
      double sin_angle = ((term.cosphase * sin_nphi) + (-((cos_nphi * term.sinphase))));
      double energy = (term.v * (1.00000000000000000e+0 + cos_angle));
      energy_added += energy;
      *energy_accumulate += energy;
    }
  }
  return energy_added;
}
void energy_fd(const dihedral_term& term, double* position, double* energy_accumulate)
{
  energy(term, position, energy_accumulate);
}

double gradient(const dihedral_term& term, double* position, double* energy_accumulate, double* force) {
  double energy_added = 0.00000000000000000e+0;
  {
    /* !BASE */
    DOUBLE x1 = position[term.i3x1 + 0];
    DOUBLE y1 = position[term.i3x1 + 1];
    DOUBLE z1 = position[term.i3x1 + 2];
    DOUBLE x2 = position[term.i3x2 + 0];
    DOUBLE y2 = position[term.i3x2 + 1];
    DOUBLE z2 = position[term.i3x2 + 2];
    DOUBLE x3 = position[term.i3x3 + 0];
    DOUBLE y3 = position[term.i3x3 + 1];
    DOUBLE z3 = position[term.i3x3 + 2];
    DOUBLE x4 = position[term.i3x4 + 0];
    DOUBLE y4 = position[term.i3x4 + 1];
    DOUBLE z4 = position[term.i3x4 + 2];
    {
      /* !BASE */
      double v1x = (x2 + (-(x1)));
      double v1y = (y2 + (-(y1)));
      double v1z = (z2 + (-(z1)));
      double v2x = (x3 + (-(x2)));
      double cse_p12_t6_g559213 = (v1z * v2x);
      double v2y = (y3 + (-(y2)));
      double cse_p12_t4_g559211 = (v1x * v2y);
      double v2z = (z3 + (-(z2)));
      double cse_p12_t5_g559212 = (v1y * v2z);
      double v3x = (x4 + (-(x3)));
      double cse_p12_t8_g559215 = ((-(v2x)) + (-(v3x)));
      double cse_p12_t11_g559218 = (v2x + v3x);
      double v3y = (y4 + (-(y3)));
      double cse_p12_t9_g559216 = ((-(v2y)) + (-(v3y)));
      double cse_p12_t12_g559219 = (v2y + v3y);
      double v3z = (z4 + (-(z3)));
      double cse_p12_t10_g559217 = ((-(v2z)) + (-(v3z)));
      double cse_p12_t13_g559220 = (v2z + v3z);
      double c1x = (cse_p12_t5_g559212 + (-((v1z * v2y))));
      double c1y = (cse_p12_t6_g559213 + (-((v1x * v2z))));
      double c1z = (cse_p12_t4_g559211 + (-((v1y * v2x))));
      double c2x = ((v2y * v3z) + (-((v2z * v3y))));
      double cse_p12_t1_g559208 = (c2x * v1x);
      double c2y = ((v2z * v3x) + (-((v2x * v3z))));
      double cse_p12_t2_g559209 = (c2y * v1y);
      double c2z = ((v2x * v3y) + (-((v2y * v3x))));
      double cse_p12_t3_g559210 = (c2z * v1z);
      double cse_p12_t7_g559214 = (cse_p12_t1_g559208 + cse_p12_t2_g559209 + cse_p12_t3_g559210);
      double v2_sq = ((v2x * v2x) + (v2y * v2y) + (v2z * v2z));
      double cse_p1_t1_sqrt559205 = sqrt(v2_sq);
      double cse_p1_t2_invsqrt559206 = (1.0 / (cse_p1_t1_sqrt559205));
      double v2_len = cse_p1_t1_sqrt559205;
      double t1 = (cse_p12_t7_g559214 * v2_len);
      double t2 = ((c1x * c2x) + (c1y * c2y) + (c1z * c2z));
      double phi = atan2(t1, t2);
      double nphi = (term.n * phi);
      double cse_p12_t14_g559221 = (-(sin(nphi)));
      double sin_nphi = sin(nphi);
      double cos_nphi = cos(nphi);
      double cos_angle = ((term.cosphase * cos_nphi) + (term.sinphase * sin_nphi));
      double sin_angle = ((term.cosphase * sin_nphi) + (-((cos_nphi * term.sinphase))));
      double energy = (term.v * (1.00000000000000000e+0 + cos_angle));
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p1_t3_invr559207 = (1.0 / ((((t1) * (t1)) + ((t2) * (t2)))));
      double cse_p12_t1_g559222 = (cse_p1_t3_invr559207 * (-(t1)));
      double cse_p12_t1_g559223 = (cse_p1_t3_invr559207 * t2);
      double cse_p12_t1_g559224 = (cos_nphi * term.n * term.sinphase);
      double cse_p12_t1_g559225 = (term.cosphase * cse_p12_t14_g559221 * term.n);
      double cse_p12_t1_g559226 = (cse_p12_t1_g559223 * v2_len);
      double g_x1 = (term.v * ((cse_p12_t1_g559224 * ((cse_p12_t1_g559222 * ((c2y * v2z) + (c2z * (-(v2y))))) + (cse_p12_t1_g559226 * (-(c2x))))) + (cse_p12_t1_g559225 * ((cse_p12_t1_g559222 * ((c2y * v2z) + (c2z * (-(v2y))))) + (cse_p12_t1_g559226 * (-(c2x)))))));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = (term.v * ((cse_p12_t1_g559224 * ((cse_p12_t1_g559222 * ((c2x * (-(v2z))) + (c2z * v2x))) + (cse_p12_t1_g559226 * (-(c2y))))) + (cse_p12_t1_g559225 * ((cse_p12_t1_g559222 * ((c2x * (-(v2z))) + (c2z * v2x))) + (cse_p12_t1_g559226 * (-(c2y)))))));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = (term.v * ((cse_p12_t1_g559224 * ((cse_p12_t1_g559222 * ((c2x * v2y) + (c2y * (-(v2x))))) + (cse_p12_t1_g559226 * (-(c2z))))) + (cse_p12_t1_g559225 * ((cse_p12_t1_g559222 * ((c2x * v2y) + (c2y * (-(v2x))))) + (cse_p12_t1_g559226 * (-(c2z)))))));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p12_t1_g559227 = (cse_p12_t7_g559214 * cse_p1_t2_invsqrt559206);
      double g_x2 = (term.v * ((cse_p12_t1_g559224 * ((cse_p12_t1_g559222 * ((c1y * v3z) + (c1z * (-(v3y))) + (c2y * ((-(v1z)) + (-(v2z)))) + (c2z * (v1y + v2y)))) + (cse_p12_t1_g559223 * ((cse_p12_t1_g559227 * (-(v2x))) + (v2_len * (c2x + (v1y * v3z) + (v1z * (-(v3y))))))))) + (cse_p12_t1_g559225 * ((cse_p12_t1_g559222 * ((c1y * v3z) + (c1z * (-(v3y))) + (c2y * ((-(v1z)) + (-(v2z)))) + (c2z * (v1y + v2y)))) + (cse_p12_t1_g559223 * ((cse_p12_t1_g559227 * (-(v2x))) + (v2_len * (c2x + (v1y * v3z) + (v1z * (-(v3y)))))))))));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (term.v * ((cse_p12_t1_g559224 * ((cse_p12_t1_g559222 * ((c1x * (-(v3z))) + (c1z * v3x) + (c2x * (v1z + v2z)) + (c2z * ((-(v1x)) + (-(v2x)))))) + (cse_p12_t1_g559223 * ((cse_p12_t1_g559227 * (-(v2y))) + (v2_len * (c2y + (v1x * (-(v3z))) + (v1z * v3x))))))) + (cse_p12_t1_g559225 * ((cse_p12_t1_g559222 * ((c1x * (-(v3z))) + (c1z * v3x) + (c2x * (v1z + v2z)) + (c2z * ((-(v1x)) + (-(v2x)))))) + (cse_p12_t1_g559223 * ((cse_p12_t1_g559227 * (-(v2y))) + (v2_len * (c2y + (v1x * (-(v3z))) + (v1z * v3x)))))))));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (term.v * ((cse_p12_t1_g559224 * ((cse_p12_t1_g559222 * ((c1x * v3y) + (c1y * (-(v3x))) + (c2x * ((-(v1y)) + (-(v2y)))) + (c2y * (v1x + v2x)))) + (cse_p12_t1_g559223 * ((cse_p12_t1_g559227 * (-(v2z))) + (v2_len * (c2z + (v1x * v3y) + (v1y * (-(v3x))))))))) + (cse_p12_t1_g559225 * ((cse_p12_t1_g559222 * ((c1x * v3y) + (c1y * (-(v3x))) + (c2x * ((-(v1y)) + (-(v2y)))) + (c2y * (v1x + v2x)))) + (cse_p12_t1_g559223 * ((cse_p12_t1_g559227 * (-(v2z))) + (v2_len * (c2z + (v1x * v3y) + (v1y * (-(v3x)))))))))));
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double g_x3 = (term.v * ((cse_p12_t1_g559224 * ((cse_p12_t1_g559222 * ((c1y * cse_p12_t10_g559217) + (c1z * cse_p12_t12_g559219) + (c2y * v1z) + (c2z * (-(v1y))))) + (cse_p12_t1_g559223 * ((cse_p12_t1_g559227 * v2x) + (v2_len * ((cse_p12_t10_g559217 * v1y) + (cse_p12_t12_g559219 * v1z))))))) + (cse_p12_t1_g559225 * ((cse_p12_t1_g559222 * ((c1y * cse_p12_t10_g559217) + (c1z * cse_p12_t12_g559219) + (c2y * v1z) + (c2z * (-(v1y))))) + (cse_p12_t1_g559223 * ((cse_p12_t1_g559227 * v2x) + (v2_len * ((cse_p12_t10_g559217 * v1y) + (cse_p12_t12_g559219 * v1z)))))))));
      KernelGradientAcc(term.i3x3, 0, g_x3);
      double g_y3 = (term.v * ((cse_p12_t1_g559224 * ((cse_p12_t1_g559222 * ((c1x * cse_p12_t13_g559220) + (c1z * cse_p12_t8_g559215) + (c2x * (-(v1z))) + (c2z * v1x))) + (cse_p12_t1_g559223 * ((cse_p12_t1_g559227 * v2y) + (v2_len * ((cse_p12_t13_g559220 * v1x) + (cse_p12_t8_g559215 * v1z))))))) + (cse_p12_t1_g559225 * ((cse_p12_t1_g559222 * ((c1x * cse_p12_t13_g559220) + (c1z * cse_p12_t8_g559215) + (c2x * (-(v1z))) + (c2z * v1x))) + (cse_p12_t1_g559223 * ((cse_p12_t1_g559227 * v2y) + (v2_len * ((cse_p12_t13_g559220 * v1x) + (cse_p12_t8_g559215 * v1z)))))))));
      KernelGradientAcc(term.i3x3, 1, g_y3);
      double g_z3 = (term.v * ((cse_p12_t1_g559224 * ((cse_p12_t1_g559222 * ((c1x * cse_p12_t9_g559216) + (c1y * cse_p12_t11_g559218) + (c2x * v1y) + (c2y * (-(v1x))))) + (cse_p12_t1_g559223 * ((cse_p12_t1_g559227 * v2z) + (v2_len * ((cse_p12_t11_g559218 * v1y) + (cse_p12_t9_g559216 * v1x))))))) + (cse_p12_t1_g559225 * ((cse_p12_t1_g559222 * ((c1x * cse_p12_t9_g559216) + (c1y * cse_p12_t11_g559218) + (c2x * v1y) + (c2y * (-(v1x))))) + (cse_p12_t1_g559223 * ((cse_p12_t1_g559227 * v2z) + (v2_len * ((cse_p12_t11_g559218 * v1y) + (cse_p12_t9_g559216 * v1x)))))))));
      KernelGradientAcc(term.i3x3, 2, g_z3);
      double g_x4 = (term.v * ((cse_p12_t1_g559224 * ((cse_p12_t1_g559222 * ((c1y * v2z) + (c1z * (-(v2y))))) + (cse_p12_t1_g559226 * (cse_p12_t5_g559212 + (v1z * (-(v2y))))))) + (cse_p12_t1_g559225 * ((cse_p12_t1_g559222 * ((c1y * v2z) + (c1z * (-(v2y))))) + (cse_p12_t1_g559226 * (cse_p12_t5_g559212 + (v1z * (-(v2y)))))))));
      KernelGradientAcc(term.i3x4, 0, g_x4);
      double g_y4 = (term.v * ((cse_p12_t1_g559224 * ((cse_p12_t1_g559222 * ((c1x * (-(v2z))) + (c1z * v2x))) + (cse_p12_t1_g559226 * (cse_p12_t6_g559213 + (v1x * (-(v2z))))))) + (cse_p12_t1_g559225 * ((cse_p12_t1_g559222 * ((c1x * (-(v2z))) + (c1z * v2x))) + (cse_p12_t1_g559226 * (cse_p12_t6_g559213 + (v1x * (-(v2z)))))))));
      KernelGradientAcc(term.i3x4, 1, g_y4);
      double g_z4 = (term.v * ((cse_p12_t1_g559224 * ((cse_p12_t1_g559222 * ((c1x * v2y) + (c1y * (-(v2x))))) + (cse_p12_t1_g559226 * (cse_p12_t4_g559211 + (v1y * (-(v2x))))))) + (cse_p12_t1_g559225 * ((cse_p12_t1_g559222 * ((c1x * v2y) + (c1y * (-(v2x))))) + (cse_p12_t1_g559226 * (cse_p12_t4_g559211 + (v1y * (-(v2x)))))))));
      KernelGradientAcc(term.i3x4, 2, g_z4);
    }
  }
  return energy_added;
}
void gradient_fd(const dihedral_term& term, double* position, double* energy_accumulate, double* force)
{
  constexpr size_t PositionSize = 12;
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
  {
    double saved = position[term.i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x3 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x3 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x3 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x3, 0, d);
  }
  {
    double saved = position[term.i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x3 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x3 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x3 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x3, 1, d);
  }
  {
    double saved = position[term.i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x3 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x3 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x3 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x3, 2, d);
  }
  {
    double saved = position[term.i3x4 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x4 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x4 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x4 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x4, 0, d);
  }
  {
    double saved = position[term.i3x4 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x4 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x4 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x4 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x4, 1, d);
  }
  {
    double saved = position[term.i3x4 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x4 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x4 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x4 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x4, 2, d);
  }
}

double hessian(const dihedral_term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double energy_added = 0.00000000000000000e+0;
  {
    /* !BASE */
    DOUBLE x1 = position[term.i3x1 + 0];
    DOUBLE y1 = position[term.i3x1 + 1];
    DOUBLE z1 = position[term.i3x1 + 2];
    DOUBLE x2 = position[term.i3x2 + 0];
    DOUBLE y2 = position[term.i3x2 + 1];
    DOUBLE z2 = position[term.i3x2 + 2];
    DOUBLE x3 = position[term.i3x3 + 0];
    DOUBLE y3 = position[term.i3x3 + 1];
    DOUBLE z3 = position[term.i3x3 + 2];
    DOUBLE x4 = position[term.i3x4 + 0];
    DOUBLE y4 = position[term.i3x4 + 1];
    DOUBLE z4 = position[term.i3x4 + 2];
    {
      /* !BASE */
      double v1x = (x2 + (-(x1)));
      double v1y = (y2 + (-(y1)));
      double v1z = (z2 + (-(z1)));
      double v2x = (x3 + (-(x2)));
      double cse_p12_t118_g559353 = (v1y * (-(v2x)));
      double cse_p12_t126_g559361 = (v1z * v2x);
      double cse_p12_t186_g559421 = ((-(v1x)) + (-(v2x)));
      double cse_p12_t195_g559430 = (v1x + v2x);
      double v2y = (y3 + (-(y2)));
      double cse_p12_t114_g559349 = (v1x * v2y);
      double cse_p12_t124_g559359 = (v1z * (-(v2y)));
      double cse_p12_t183_g559418 = (cse_p12_t114_g559349 + cse_p12_t118_g559353);
      double cse_p12_t187_g559422 = ((-(v1y)) + (-(v2y)));
      double cse_p12_t196_g559431 = (v1y + v2y);
      double v2z = (z3 + (-(z2)));
      double cse_p12_t112_g559347 = (v1x * (-(v2z)));
      double cse_p12_t120_g559355 = (v1y * v2z);
      double cse_p12_t182_g559417 = (cse_p12_t112_g559347 + cse_p12_t126_g559361);
      double cse_p12_t185_g559420 = (cse_p12_t120_g559355 + cse_p12_t124_g559359);
      double cse_p12_t188_g559423 = ((-(v1z)) + (-(v2z)));
      double cse_p12_t197_g559432 = (v1z + v2z);
      double v3x = (x4 + (-(x3)));
      double cse_p12_t117_g559352 = (v1y * (v2x + v3x));
      double cse_p12_t119_g559354 = (v1y * (-(v3x)));
      double cse_p12_t122_g559357 = (v1z * ((-(v2x)) + (-(v3x))));
      double cse_p12_t127_g559362 = (v1z * v3x);
      double cse_p12_t189_g559424 = ((-(v2x)) + (-(v3x)));
      double cse_p12_t198_g559433 = (v2x + v3x);
      double v3y = (y4 + (-(y3)));
      double cse_p12_t110_g559345 = (v1x * ((-(v2y)) + (-(v3y))));
      double cse_p12_t115_g559350 = (v1x * v3y);
      double cse_p12_t123_g559358 = (v1z * (v2y + v3y));
      double cse_p12_t125_g559360 = (v1z * (-(v3y)));
      double cse_p12_t180_g559415 = (cse_p12_t110_g559345 + cse_p12_t117_g559352);
      double cse_p12_t190_g559425 = ((-(v2y)) + (-(v3y)));
      double cse_p12_t199_g559434 = (v2y + v3y);
      double v3z = (z4 + (-(z3)));
      double cse_p12_t111_g559346 = (v1x * (v2z + v3z));
      double cse_p12_t113_g559348 = (v1x * (-(v3z)));
      double cse_p12_t116_g559351 = (v1y * ((-(v2z)) + (-(v3z))));
      double cse_p12_t121_g559356 = (v1y * v3z);
      double cse_p12_t181_g559416 = (cse_p12_t111_g559346 + cse_p12_t122_g559357);
      double cse_p12_t184_g559419 = (cse_p12_t116_g559351 + cse_p12_t123_g559358);
      double cse_p12_t191_g559426 = ((-(v2z)) + (-(v3z)));
      double cse_p12_t200_g559435 = (v2z + v3z);
      double c1x = (cse_p12_t120_g559355 + (-((v1z * v2y))));
      double cse_p12_t23_g559258 = (c1x * cse_p12_t190_g559425);
      double cse_p12_t24_g559259 = (c1x * cse_p12_t200_g559435);
      double cse_p12_t25_g559260 = (c1x * (-(v2z)));
      double cse_p12_t26_g559261 = (c1x * (-(v3z)));
      double cse_p12_t27_g559262 = (c1x * v2y);
      double cse_p12_t28_g559263 = (c1x * v3y);
      double c1y = (cse_p12_t126_g559361 + (-((v1x * v2z))));
      double cse_p12_t29_g559264 = (c1y * cse_p12_t191_g559426);
      double cse_p12_t30_g559265 = (c1y * cse_p12_t198_g559433);
      double cse_p12_t31_g559266 = (c1y * (-(v2x)));
      double cse_p12_t32_g559267 = (c1y * (-(v3x)));
      double cse_p12_t33_g559268 = (c1y * v2z);
      double cse_p12_t34_g559269 = (c1y * v3z);
      double cse_p12_t171_g559406 = (cse_p12_t27_g559262 + cse_p12_t31_g559266);
      double c1z = (cse_p12_t114_g559349 + (-((v1y * v2x))));
      double cse_p12_t35_g559270 = (c1z * cse_p12_t189_g559424);
      double cse_p12_t36_g559271 = (c1z * cse_p12_t199_g559434);
      double cse_p12_t37_g559272 = (c1z * (-(v2y)));
      double cse_p12_t38_g559273 = (c1z * (-(v3y)));
      double cse_p12_t39_g559274 = (c1z * v2x);
      double cse_p12_t40_g559275 = (c1z * v3x);
      double cse_p12_t169_g559404 = (cse_p12_t25_g559260 + cse_p12_t39_g559274);
      double cse_p12_t174_g559409 = (cse_p12_t33_g559268 + cse_p12_t37_g559272);
      double c2x = ((v2y * v3z) + (-((v2z * v3y))));
      double cse_p12_t41_g559276 = (c2x * cse_p12_t187_g559422);
      double cse_p12_t42_g559277 = (c2x * cse_p12_t197_g559432);
      double cse_p12_t43_g559278 = (c2x * (-(v1z)));
      double cse_p12_t44_g559279 = (c2x * (-(v2z)));
      double cse_p12_t45_g559280 = (c2x * v1x);
      double cse_p12_t46_g559281 = (c2x * v1y);
      double cse_p12_t47_g559282 = (c2x * v2y);
      double cse_p12_t192_g559427 = (c2x + cse_p12_t121_g559356 + cse_p12_t125_g559360);
      double c2y = ((v2z * v3x) + (-((v2x * v3z))));
      double cse_p12_t48_g559283 = (c2y * cse_p12_t188_g559423);
      double cse_p12_t49_g559284 = (c2y * cse_p12_t195_g559430);
      double cse_p12_t50_g559285 = (c2y * (-(v1x)));
      double cse_p12_t51_g559286 = (c2y * (-(v2x)));
      double cse_p12_t52_g559287 = (c2y * v1y);
      double cse_p12_t53_g559288 = (c2y * v1z);
      double cse_p12_t54_g559289 = (c2y * v2z);
      double cse_p12_t167_g559402 = (cse_p12_t23_g559258 + cse_p12_t30_g559265 + cse_p12_t46_g559281 + cse_p12_t50_g559285);
      double cse_p12_t172_g559407 = (cse_p12_t28_g559263 + cse_p12_t32_g559267 + cse_p12_t41_g559276 + cse_p12_t49_g559284);
      double cse_p12_t178_g559413 = (cse_p12_t47_g559282 + cse_p12_t51_g559286);
      double cse_p12_t193_g559428 = (c2y + cse_p12_t113_g559348 + cse_p12_t127_g559362);
      double c2z = ((v2x * v3y) + (-((v2y * v3x))));
      double cse_p12_t55_g559290 = (c2z * cse_p12_t186_g559421);
      double cse_p12_t56_g559291 = (c2z * cse_p12_t196_g559431);
      double cse_p12_t57_g559292 = (c2z * (-(v1y)));
      double cse_p12_t58_g559293 = (c2z * (-(v2y)));
      double cse_p12_t59_g559294 = (c2z * v1x);
      double cse_p12_t60_g559295 = (c2z * v1z);
      double cse_p12_t61_g559296 = (c2z * v2x);
      double cse_p12_t168_g559403 = (cse_p12_t24_g559259 + cse_p12_t35_g559270 + cse_p12_t43_g559278 + cse_p12_t59_g559294);
      double cse_p12_t170_g559405 = (cse_p12_t26_g559261 + cse_p12_t40_g559275 + cse_p12_t42_g559277 + cse_p12_t55_g559290);
      double cse_p12_t173_g559408 = (cse_p12_t29_g559264 + cse_p12_t36_g559271 + cse_p12_t53_g559288 + cse_p12_t57_g559292);
      double cse_p12_t175_g559410 = (cse_p12_t34_g559269 + cse_p12_t38_g559273 + cse_p12_t48_g559283 + cse_p12_t56_g559291);
      double cse_p12_t176_g559411 = (cse_p12_t44_g559279 + cse_p12_t61_g559296);
      double cse_p12_t177_g559412 = (cse_p12_t45_g559280 + cse_p12_t52_g559287 + cse_p12_t60_g559295);
      double cse_p12_t179_g559414 = (cse_p12_t54_g559289 + cse_p12_t58_g559293);
      double cse_p12_t194_g559429 = (c2z + cse_p12_t115_g559350 + cse_p12_t119_g559354);
      double v2_sq = ((v2x * v2x) + (v2y * v2y) + (v2z * v2z));
      double cse_p11_t1_sqrt559231 = sqrt(v2_sq);
      double cse_p11_t2_invsqrt559232 = (1.0 / (cse_p11_t1_sqrt559231));
      double cse_p12_t1_g559459 = (cse_p11_t2_invsqrt559232 * cse_p12_t177_g559412);
      double cse_p12_t62_g559297 = (cse_p12_t1_g559459 * (-(v2x)));
      double cse_p12_t63_g559298 = (cse_p12_t1_g559459 * (-(v2y)));
      double cse_p12_t64_g559299 = (cse_p12_t1_g559459 * (-(v2z)));
      double cse_p12_t72_g559307 = (cse_p12_t1_g559459 * v2x);
      double cse_p12_t73_g559308 = (cse_p12_t1_g559459 * v2y);
      double cse_p12_t74_g559309 = (cse_p12_t1_g559459 * v2z);
      double cse_p12_t201_g559436 = (-(cse_p12_t1_g559459));
      double cse_p11_t3_invr559233 = (cse_p11_t2_invsqrt559232 * cse_p11_t2_invsqrt559232);
      double cse_p12_t66_g559301 = (cse_p11_t3_invr559233 * cse_p12_t1_g559459 * v2x * v2x);
      double cse_p12_t67_g559302 = (cse_p11_t3_invr559233 * cse_p12_t1_g559459 * v2x * v2y);
      double cse_p12_t68_g559303 = (cse_p11_t3_invr559233 * cse_p12_t1_g559459 * v2x * v2z);
      double cse_p12_t69_g559304 = (cse_p11_t3_invr559233 * cse_p12_t1_g559459 * v2y * v2y);
      double cse_p12_t70_g559305 = (cse_p11_t3_invr559233 * cse_p12_t1_g559459 * v2y * v2z);
      double cse_p12_t71_g559306 = (cse_p11_t3_invr559233 * cse_p12_t1_g559459 * v2z * v2z);
      double v2_len = cse_p11_t1_sqrt559231;
      double cse_p12_t1_g559466 = (cse_p12_t180_g559415 * v2_len);
      double cse_p12_t1_g559467 = (cse_p12_t181_g559416 * v2_len);
      double cse_p12_t130_g559365 = (cse_p12_t182_g559417 * v2_len);
      double cse_p12_t131_g559366 = (cse_p12_t183_g559418 * v2_len);
      double cse_p12_t1_g559468 = (cse_p12_t184_g559419 * v2_len);
      double cse_p12_t133_g559368 = (cse_p12_t185_g559420 * v2_len);
      double cse_p12_t134_g559369 = (cse_p12_t192_g559427 * v2_len);
      double cse_p12_t135_g559370 = (cse_p12_t193_g559428 * v2_len);
      double cse_p12_t136_g559371 = (cse_p12_t194_g559429 * v2_len);
      double cse_p12_t137_g559372 = (v2_len * (-(c2z)));
      double cse_p12_t147_g559382 = (cse_p12_t134_g559369 + cse_p12_t62_g559297);
      double cse_p12_t148_g559383 = (cse_p12_t135_g559370 + cse_p12_t63_g559298);
      double cse_p12_t149_g559384 = (cse_p12_t136_g559371 + cse_p12_t64_g559299);
      double cse_p12_t150_g559385 = (cse_p12_t1_g559468 + cse_p12_t72_g559307);
      double cse_p12_t151_g559386 = (cse_p12_t1_g559467 + cse_p12_t73_g559308);
      double cse_p12_t152_g559387 = (cse_p12_t1_g559466 + cse_p12_t74_g559309);
      double cse_p12_t213_g559448 = (-(cse_p12_t130_g559365));
      double cse_p12_t214_g559449 = (-(cse_p12_t131_g559366));
      double cse_p12_t215_g559450 = (-(cse_p12_t133_g559368));
      double cse_p12_t216_g559451 = (-(cse_p12_t147_g559382));
      double cse_p12_t217_g559452 = (-(cse_p12_t148_g559383));
      double cse_p12_t218_g559453 = (-(cse_p12_t149_g559384));
      double cse_p12_t219_g559454 = (-(cse_p12_t150_g559385));
      double cse_p12_t220_g559455 = (-(cse_p12_t151_g559386));
      double cse_p12_t221_g559456 = (-(cse_p12_t152_g559387));
      double cse_p12_t222_g559457 = (-(cse_p12_t137_g559372));
      double t1 = (cse_p12_t177_g559412 * v2_len);
      double cse_p13_t1_g559475 = (2.00000000000000000e+0 * t1);
      double cse_p12_t1_g559236 = (cse_p12_t147_g559382 * cse_p13_t1_g559475);
      double cse_p12_t2_g559237 = (cse_p12_t148_g559383 * cse_p13_t1_g559475);
      double cse_p12_t3_g559238 = (cse_p12_t149_g559384 * cse_p13_t1_g559475);
      double cse_p12_t4_g559239 = (cse_p12_t150_g559385 * cse_p13_t1_g559475);
      double cse_p12_t5_g559240 = (cse_p12_t151_g559386 * cse_p13_t1_g559475);
      double cse_p12_t6_g559241 = (cse_p12_t152_g559387 * cse_p13_t1_g559475);
      double cse_p12_t7_g559242 = (cse_p12_t182_g559417 * cse_p13_t1_g559475 * v2_len);
      double cse_p12_t8_g559243 = (cse_p12_t183_g559418 * cse_p13_t1_g559475 * v2_len);
      double cse_p12_t9_g559244 = (cse_p12_t185_g559420 * cse_p13_t1_g559475 * v2_len);
      double cse_p12_t10_g559245 = (cse_p13_t1_g559475 * v2_len * (-(c2y)));
      double cse_p12_t11_g559246 = (cse_p13_t1_g559475 * v2_len * (-(c2z)));
      double t2 = ((c1x * c2x) + (c1y * c2y) + (c1z * c2z));
      double cse_p13_t1_g559476 = (2.00000000000000000e+0 * t2);
      double cse_p12_t12_g559247 = (cse_p12_t167_g559402 * cse_p13_t1_g559476);
      double cse_p12_t13_g559248 = (cse_p12_t168_g559403 * cse_p13_t1_g559476);
      double cse_p12_t14_g559249 = (cse_p12_t169_g559404 * cse_p13_t1_g559476);
      double cse_p12_t15_g559250 = (cse_p12_t170_g559405 * cse_p13_t1_g559476);
      double cse_p12_t16_g559251 = (cse_p12_t171_g559406 * cse_p13_t1_g559476);
      double cse_p12_t17_g559252 = (cse_p12_t172_g559407 * cse_p13_t1_g559476);
      double cse_p12_t18_g559253 = (cse_p12_t173_g559408 * cse_p13_t1_g559476);
      double cse_p12_t19_g559254 = (cse_p12_t174_g559409 * cse_p13_t1_g559476);
      double cse_p12_t20_g559255 = (cse_p12_t175_g559410 * cse_p13_t1_g559476);
      double cse_p12_t21_g559256 = (cse_p12_t176_g559411 * cse_p13_t1_g559476);
      double cse_p12_t22_g559257 = (cse_p12_t178_g559413 * cse_p13_t1_g559476);
      double cse_p12_t138_g559373 = (cse_p12_t1_g559236 + cse_p12_t20_g559255);
      double cse_p12_t139_g559374 = (cse_p12_t15_g559250 + cse_p12_t2_g559237);
      double cse_p12_t140_g559375 = (cse_p12_t17_g559252 + cse_p12_t3_g559238);
      double cse_p12_t141_g559376 = (cse_p12_t18_g559253 + cse_p12_t4_g559239);
      double cse_p12_t142_g559377 = (cse_p12_t13_g559248 + cse_p12_t5_g559240);
      double cse_p12_t143_g559378 = (cse_p12_t12_g559247 + cse_p12_t6_g559241);
      double cse_p12_t144_g559379 = (cse_p12_t14_g559249 + cse_p12_t7_g559242);
      double cse_p12_t145_g559380 = (cse_p12_t16_g559251 + cse_p12_t8_g559243);
      double cse_p12_t146_g559381 = (cse_p12_t19_g559254 + cse_p12_t9_g559244);
      double cse_p12_t165_g559400 = (cse_p12_t10_g559245 + cse_p12_t21_g559256);
      double cse_p12_t166_g559401 = (cse_p12_t11_g559246 + cse_p12_t22_g559257);
      double phi = atan2(t1, t2);
      double nphi = (term.n * phi);
      double cse_p12_t223_g559458 = (-(sin(nphi)));
      double sin_nphi = sin(nphi);
      double cos_nphi = cos(nphi);
      double cos_angle = ((term.cosphase * cos_nphi) + (term.sinphase * sin_nphi));
      double sin_angle = ((term.cosphase * sin_nphi) + (-((cos_nphi * term.sinphase))));
      double energy = (term.v * (1.00000000000000000e+0 + cos_angle));
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t4_invr559234 = (1.0 / ((((t1) * (t1)) + ((t2) * (t2)))));
      double cse_p13_t1_g559473 = (cse_p11_t4_invr559234 * (-(t1)));
      double cse_p12_t75_g559310 = (cse_p12_t167_g559402 * cse_p13_t1_g559473);
      double cse_p12_t76_g559311 = (cse_p12_t168_g559403 * cse_p13_t1_g559473);
      double cse_p12_t77_g559312 = (cse_p12_t169_g559404 * cse_p13_t1_g559473);
      double cse_p12_t78_g559313 = (cse_p12_t170_g559405 * cse_p13_t1_g559473);
      double cse_p12_t79_g559314 = (cse_p12_t171_g559406 * cse_p13_t1_g559473);
      double cse_p12_t80_g559315 = (cse_p12_t172_g559407 * cse_p13_t1_g559473);
      double cse_p12_t81_g559316 = (cse_p12_t173_g559408 * cse_p13_t1_g559473);
      double cse_p12_t82_g559317 = (cse_p12_t174_g559409 * cse_p13_t1_g559473);
      double cse_p12_t83_g559318 = (cse_p12_t175_g559410 * cse_p13_t1_g559473);
      double cse_p12_t84_g559319 = (cse_p12_t176_g559411 * cse_p13_t1_g559473);
      double cse_p12_t85_g559320 = (cse_p12_t178_g559413 * cse_p13_t1_g559473);
      double cse_p12_t86_g559321 = (cse_p12_t179_g559414 * cse_p13_t1_g559473);
      double cse_p12_t1_g559462 = (cse_p11_t4_invr559234 * t2);
      double cse_p12_t87_g559322 = (cse_p12_t147_g559382 * cse_p12_t1_g559462);
      double cse_p12_t88_g559323 = (cse_p12_t148_g559383 * cse_p12_t1_g559462);
      double cse_p12_t89_g559324 = (cse_p12_t149_g559384 * cse_p12_t1_g559462);
      double cse_p12_t90_g559325 = (cse_p12_t150_g559385 * cse_p12_t1_g559462);
      double cse_p12_t91_g559326 = (cse_p12_t151_g559386 * cse_p12_t1_g559462);
      double cse_p12_t92_g559327 = (cse_p12_t152_g559387 * cse_p12_t1_g559462);
      double cse_p13_t1_g559474 = (cse_p12_t1_g559462 * v2_len);
      double cse_p12_t93_g559328 = (cse_p12_t182_g559417 * cse_p13_t1_g559474);
      double cse_p12_t94_g559329 = (cse_p12_t183_g559418 * cse_p13_t1_g559474);
      double cse_p12_t95_g559330 = (cse_p12_t185_g559420 * cse_p13_t1_g559474);
      double cse_p12_t96_g559331 = (cse_p13_t1_g559474 * (-(c2x)));
      double cse_p12_t97_g559332 = (cse_p13_t1_g559474 * (-(c2y)));
      double cse_p12_t98_g559333 = (cse_p13_t1_g559474 * (-(c2z)));
      double cse_p12_t153_g559388 = (cse_p12_t75_g559310 + cse_p12_t92_g559327);
      double cse_p12_t154_g559389 = (cse_p12_t76_g559311 + cse_p12_t91_g559326);
      double cse_p12_t155_g559390 = (cse_p12_t77_g559312 + cse_p12_t93_g559328);
      double cse_p12_t156_g559391 = (cse_p12_t78_g559313 + cse_p12_t88_g559323);
      double cse_p12_t157_g559392 = (cse_p12_t79_g559314 + cse_p12_t94_g559329);
      double cse_p12_t158_g559393 = (cse_p12_t80_g559315 + cse_p12_t89_g559324);
      double cse_p12_t159_g559394 = (cse_p12_t81_g559316 + cse_p12_t90_g559325);
      double cse_p12_t160_g559395 = (cse_p12_t82_g559317 + cse_p12_t95_g559330);
      double cse_p12_t161_g559396 = (cse_p12_t83_g559318 + cse_p12_t87_g559322);
      double cse_p12_t162_g559397 = (cse_p12_t84_g559319 + cse_p12_t97_g559332);
      double cse_p12_t163_g559398 = (cse_p12_t85_g559320 + cse_p12_t98_g559333);
      double cse_p12_t164_g559399 = (cse_p12_t86_g559321 + cse_p12_t96_g559331);
      double cse_p11_t5_invr2559235 = (cse_p11_t4_invr559234 * cse_p11_t4_invr559234);
      double cse_p12_t99_g559334 = (cse_p11_t5_invr2559235 * cse_p12_t138_g559373);
      double cse_p12_t100_g559335 = (cse_p11_t5_invr2559235 * cse_p12_t139_g559374);
      double cse_p12_t101_g559336 = (cse_p11_t5_invr2559235 * cse_p12_t140_g559375);
      double cse_p12_t102_g559337 = (cse_p11_t5_invr2559235 * cse_p12_t141_g559376);
      double cse_p12_t103_g559338 = (cse_p11_t5_invr2559235 * cse_p12_t142_g559377);
      double cse_p12_t104_g559339 = (cse_p11_t5_invr2559235 * cse_p12_t143_g559378);
      double cse_p12_t105_g559340 = (cse_p11_t5_invr2559235 * cse_p12_t144_g559379);
      double cse_p12_t106_g559341 = (cse_p11_t5_invr2559235 * cse_p12_t145_g559380);
      double cse_p12_t107_g559342 = (cse_p11_t5_invr2559235 * cse_p12_t146_g559381);
      double cse_p12_t108_g559343 = (cse_p11_t5_invr2559235 * cse_p12_t165_g559400);
      double cse_p12_t109_g559344 = (cse_p11_t5_invr2559235 * cse_p12_t166_g559401);
      double cse_p12_t202_g559437 = (-(cse_p12_t99_g559334));
      double cse_p12_t203_g559438 = (-(cse_p12_t100_g559335));
      double cse_p12_t204_g559439 = (-(cse_p12_t101_g559336));
      double cse_p12_t205_g559440 = (-(cse_p12_t102_g559337));
      double cse_p12_t206_g559441 = (-(cse_p12_t103_g559338));
      double cse_p12_t207_g559442 = (-(cse_p12_t104_g559339));
      double cse_p12_t208_g559443 = (-(cse_p12_t105_g559340));
      double cse_p12_t209_g559444 = (-(cse_p12_t106_g559341));
      double cse_p12_t210_g559445 = (-(cse_p12_t107_g559342));
      double cse_p12_t211_g559446 = (-(cse_p12_t108_g559343));
      double cse_p12_t212_g559447 = (-(cse_p12_t109_g559344));
      double cse_p12_t1_g559465 = (t2 * v2_len);
      double cse_p13_t1_g559478 = (cse_p12_t1_g559465 * (-(c2x)));
      double cse_p13_t1_g559479 = (cos_nphi * term.n * term.sinphase);
      double cse_p14_t1_g559486 = (term.cosphase * cse_p12_t223_g559458 * term.n);
      double g_x1 = (term.v * ((cse_p13_t1_g559479 * (cse_p12_t86_g559321 + (cse_p11_t4_invr559234 * cse_p13_t1_g559478))) + (cse_p14_t1_g559486 * (cse_p12_t86_g559321 + (cse_p11_t4_invr559234 * cse_p13_t1_g559478)))));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double cse_p13_t1_g559480 = (cse_p12_t1_g559465 * (-(c2y)));
      double g_y1 = (term.v * ((cse_p13_t1_g559479 * (cse_p12_t84_g559319 + (cse_p11_t4_invr559234 * cse_p13_t1_g559480))) + (cse_p14_t1_g559486 * (cse_p12_t84_g559319 + (cse_p11_t4_invr559234 * cse_p13_t1_g559480)))));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double cse_p14_t1_g559482 = (cse_p12_t1_g559465 * (-(c2z)));
      double g_z1 = (term.v * ((cse_p13_t1_g559479 * (cse_p12_t85_g559320 + (cse_p11_t4_invr559234 * cse_p14_t1_g559482))) + (cse_p14_t1_g559486 * (cse_p12_t85_g559320 + (cse_p11_t4_invr559234 * cse_p14_t1_g559482)))));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double g_x2 = (term.v * ((cse_p13_t1_g559479 * (cse_p12_t83_g559318 + (cse_p11_t4_invr559234 * t2 * (cse_p12_t134_g559369 + (cse_p11_t2_invsqrt559232 * cse_p12_t177_g559412 * (-(v2x))))))) + (cse_p14_t1_g559486 * (cse_p12_t83_g559318 + (cse_p11_t4_invr559234 * t2 * (cse_p12_t134_g559369 + (cse_p11_t2_invsqrt559232 * cse_p12_t177_g559412 * (-(v2x)))))))));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (term.v * ((cse_p13_t1_g559479 * (cse_p12_t78_g559313 + (cse_p11_t4_invr559234 * t2 * (cse_p12_t135_g559370 + (cse_p11_t2_invsqrt559232 * cse_p12_t177_g559412 * (-(v2y))))))) + (cse_p14_t1_g559486 * (cse_p12_t78_g559313 + (cse_p11_t4_invr559234 * t2 * (cse_p12_t135_g559370 + (cse_p11_t2_invsqrt559232 * cse_p12_t177_g559412 * (-(v2y)))))))));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (term.v * ((cse_p13_t1_g559479 * (cse_p12_t80_g559315 + (cse_p11_t4_invr559234 * t2 * (cse_p12_t136_g559371 + (cse_p11_t2_invsqrt559232 * cse_p12_t177_g559412 * (-(v2z))))))) + (cse_p14_t1_g559486 * (cse_p12_t80_g559315 + (cse_p11_t4_invr559234 * t2 * (cse_p12_t136_g559371 + (cse_p11_t2_invsqrt559232 * cse_p12_t177_g559412 * (-(v2z)))))))));
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double g_x3 = (term.v * ((cse_p13_t1_g559479 * (cse_p12_t81_g559316 + (cse_p11_t4_invr559234 * t2 * (cse_p12_t1_g559468 + (cse_p11_t2_invsqrt559232 * cse_p12_t177_g559412 * v2x))))) + (cse_p14_t1_g559486 * (cse_p12_t81_g559316 + (cse_p11_t4_invr559234 * t2 * (cse_p12_t1_g559468 + (cse_p11_t2_invsqrt559232 * cse_p12_t177_g559412 * v2x)))))));
      KernelGradientAcc(term.i3x3, 0, g_x3);
      double g_y3 = (term.v * ((cse_p13_t1_g559479 * (cse_p12_t76_g559311 + (cse_p11_t4_invr559234 * t2 * (cse_p12_t1_g559467 + (cse_p11_t2_invsqrt559232 * cse_p12_t177_g559412 * v2y))))) + (cse_p14_t1_g559486 * (cse_p12_t76_g559311 + (cse_p11_t4_invr559234 * t2 * (cse_p12_t1_g559467 + (cse_p11_t2_invsqrt559232 * cse_p12_t177_g559412 * v2y)))))));
      KernelGradientAcc(term.i3x3, 1, g_y3);
      double g_z3 = (term.v * ((cse_p13_t1_g559479 * (cse_p12_t75_g559310 + (cse_p11_t4_invr559234 * t2 * (cse_p12_t1_g559466 + (cse_p11_t2_invsqrt559232 * cse_p12_t177_g559412 * v2z))))) + (cse_p14_t1_g559486 * (cse_p12_t75_g559310 + (cse_p11_t4_invr559234 * t2 * (cse_p12_t1_g559466 + (cse_p11_t2_invsqrt559232 * cse_p12_t177_g559412 * v2z)))))));
      KernelGradientAcc(term.i3x3, 2, g_z3);
      double cse_p12_t1_g559471 = (cse_p11_t4_invr559234 * cse_p12_t185_g559420 * cse_p12_t1_g559465);
      double g_x4 = (term.v * ((cse_p13_t1_g559479 * (cse_p12_t1_g559471 + cse_p12_t82_g559317)) + (cse_p14_t1_g559486 * (cse_p12_t1_g559471 + cse_p12_t82_g559317))));
      KernelGradientAcc(term.i3x4, 0, g_x4);
      double cse_p12_t1_g559469 = (cse_p11_t4_invr559234 * cse_p12_t182_g559417 * cse_p12_t1_g559465);
      double g_y4 = (term.v * ((cse_p13_t1_g559479 * (cse_p12_t1_g559469 + cse_p12_t77_g559312)) + (cse_p14_t1_g559486 * (cse_p12_t1_g559469 + cse_p12_t77_g559312))));
      KernelGradientAcc(term.i3x4, 1, g_y4);
      double cse_p12_t1_g559470 = (cse_p11_t4_invr559234 * cse_p12_t183_g559418 * cse_p12_t1_g559465);
      double g_z4 = (term.v * ((cse_p13_t1_g559479 * (cse_p12_t1_g559470 + cse_p12_t79_g559314)) + (cse_p14_t1_g559486 * (cse_p12_t1_g559470 + cse_p12_t79_g559314))));
      KernelGradientAcc(term.i3x4, 2, g_z4);
      double cse_p12_t1_g559460 = (term.n * (-(term.v)));
      double cse_p12_t1_g559463 = (cos_angle * cse_p12_t1_g559460 * term.n);
      double cse_p12_t1_g559464 = (cse_p12_t1_g559460 * sin_angle);
      double cse_p13_t1_g559472 = (cse_p11_t4_invr559234 * v2_len);
      double cse_p13_t1_g559481 = (cse_p13_t1_g559472 * (-(c2x)));
      double h_x1_x1 = ((cse_p12_t164_g559399 * cse_p12_t164_g559399 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p12_t179_g559414 * (cse_p13_t1_g559481 + (-((cse_p11_t4_invr559234 * v2_len * (-(c2x))))) + (-((cse_p11_t5_invr2559235 * ((cse_p12_t179_g559414 * cse_p13_t1_g559476) + (cse_p13_t1_g559475 * v2_len * (-(c2x)))) * (-(t1))))))) + (-((cse_p11_t5_invr2559235 * cse_p13_t1_g559478 * ((cse_p12_t179_g559414 * cse_p13_t1_g559476) + (cse_p13_t1_g559475 * v2_len * (-(c2x))))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
      double h_x1_y1 = ((cse_p12_t162_g559397 * cse_p12_t164_g559399 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p12_t176_g559411 * cse_p13_t1_g559481) + (cse_p12_t179_g559414 * cse_p12_t211_g559446 * (-(t1))) + (cse_p12_t211_g559446 * cse_p13_t1_g559478) + (-((cse_p11_t4_invr559234 * cse_p12_t179_g559414 * v2_len * (-(c2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
      double h_x1_z1 = ((cse_p12_t163_g559398 * cse_p12_t164_g559399 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t179_g559414 * cse_p12_t222_g559457) + (cse_p12_t178_g559413 * cse_p13_t1_g559481) + (cse_p12_t179_g559414 * cse_p12_t212_g559447 * (-(t1))) + (cse_p12_t212_g559447 * cse_p13_t1_g559478))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
      double cse_p13_t1_g559477 = (cse_p11_t2_invsqrt559232 * cse_p12_t1_g559462);
      double h_x1_x2 = ((cse_p12_t161_g559396 * cse_p12_t164_g559399 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t179_g559414 * cse_p12_t216_g559451) + (cse_p12_t175_g559410 * cse_p13_t1_g559481) + (cse_p12_t179_g559414 * cse_p12_t202_g559437 * (-(t1))) + (cse_p12_t202_g559437 * cse_p13_t1_g559478) + (cse_p13_t1_g559473 * (((-(v2y)) * (-(v3y))) + (v2z * v3z))) + (-((cse_p13_t1_g559477 * v2x * (-(c2x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
      double h_x1_y2 = ((cse_p12_t156_g559391 * cse_p12_t164_g559399 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t179_g559414 * cse_p12_t217_g559452) + (cse_p12_t170_g559405 * cse_p13_t1_g559481) + (cse_p12_t179_g559414 * cse_p12_t203_g559438 * (-(t1))) + (cse_p12_t203_g559438 * cse_p13_t1_g559478) + (cse_p13_t1_g559473 * (c2z + (v3x * (-(v2y))))) + (cse_p13_t1_g559474 * v3z) + (-((cse_p13_t1_g559477 * v2y * (-(c2x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
      double h_x1_z2 = ((cse_p12_t158_g559393 * cse_p12_t164_g559399 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t179_g559414 * cse_p12_t218_g559453) + (cse_p12_t172_g559407 * cse_p13_t1_g559481) + (cse_p12_t179_g559414 * cse_p12_t204_g559439 * (-(t1))) + (cse_p12_t204_g559439 * cse_p13_t1_g559478) + (cse_p13_t1_g559473 * ((v2z * (-(v3x))) + (-(c2y)))) + (cse_p13_t1_g559474 * (-(v3y))) + (-((cse_p13_t1_g559477 * v2z * (-(c2x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
      double h_x1_x3 = ((cse_p12_t159_g559394 * cse_p12_t164_g559399 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t179_g559414 * cse_p12_t219_g559454) + (cse_p12_t173_g559408 * cse_p13_t1_g559481) + (cse_p12_t179_g559414 * cse_p12_t205_g559440 * (-(t1))) + (cse_p12_t205_g559440 * cse_p13_t1_g559478) + (cse_p13_t1_g559473 * ((cse_p12_t191_g559426 * v2z) + (cse_p12_t199_g559434 * (-(v2y))))) + (cse_p13_t1_g559477 * v2x * (-(c2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 0, h_x1_x3);
      double h_x1_y3 = ((cse_p12_t154_g559389 * cse_p12_t164_g559399 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t179_g559414 * cse_p12_t220_g559455) + (cse_p12_t168_g559403 * cse_p13_t1_g559481) + (cse_p12_t179_g559414 * cse_p12_t206_g559441 * (-(t1))) + (cse_p12_t191_g559426 * cse_p13_t1_g559474) + (cse_p12_t206_g559441 * cse_p13_t1_g559478) + (cse_p13_t1_g559473 * ((cse_p12_t189_g559424 * (-(v2y))) + (-(c2z)))) + (cse_p13_t1_g559477 * v2y * (-(c2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 1, h_x1_y3);
      double cse_p14_t1_g559489 = (cse_p12_t207_g559442 * (-(t1)));
      double h_x1_z3 = ((cse_p12_t153_g559388 * cse_p12_t164_g559399 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t179_g559414 * cse_p12_t221_g559456) + (cse_p12_t167_g559402 * cse_p13_t1_g559481) + (cse_p12_t179_g559414 * cse_p14_t1_g559489) + (cse_p12_t199_g559434 * cse_p13_t1_g559474) + (cse_p12_t207_g559442 * cse_p13_t1_g559478) + (cse_p13_t1_g559473 * (c2y + (cse_p12_t198_g559433 * v2z))) + (cse_p13_t1_g559477 * v2z * (-(c2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 2, h_x1_z3);
      double cse_p14_t1_g559488 = (cse_p12_t210_g559445 * (-(t1)));
      double h_x1_x4 = ((cse_p12_t160_g559395 * cse_p12_t164_g559399 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t179_g559414 * cse_p12_t215_g559450) + (cse_p12_t174_g559409 * cse_p13_t1_g559481) + (cse_p12_t179_g559414 * cse_p14_t1_g559488) + (cse_p12_t210_g559445 * cse_p13_t1_g559478) + (cse_p13_t1_g559473 * (((-(v2y)) * (-(v2y))) + (v2z * v2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x4, 0, h_x1_x4);
      double cse_p14_t1_g559485 = (cse_p12_t208_g559443 * (-(t1)));
      double h_x1_y4 = ((cse_p12_t155_g559390 * cse_p12_t164_g559399 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t179_g559414 * cse_p12_t213_g559448) + (cse_p12_t169_g559404 * cse_p13_t1_g559481) + (cse_p12_t179_g559414 * cse_p14_t1_g559485) + (cse_p12_t208_g559443 * cse_p13_t1_g559478) + (cse_p13_t1_g559473 * v2x * (-(v2y))) + (cse_p13_t1_g559474 * v2z))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x4, 1, h_x1_y4);
      double cse_p14_t1_g559483 = (cse_p12_t209_g559444 * (-(t1)));
      double h_x1_z4 = ((cse_p12_t157_g559392 * cse_p12_t164_g559399 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t179_g559414 * cse_p12_t214_g559449) + (cse_p12_t171_g559406 * cse_p13_t1_g559481) + (cse_p12_t179_g559414 * cse_p14_t1_g559483) + (cse_p12_t209_g559444 * cse_p13_t1_g559478) + (cse_p13_t1_g559473 * v2z * (-(v2x))) + (cse_p13_t1_g559474 * (-(v2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x4, 2, h_x1_z4);
      double cse_p14_t1_g559484 = (cse_p13_t1_g559472 * (-(c2y)));
      double h_y1_y1 = ((cse_p12_t162_g559397 * cse_p12_t162_g559397 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p12_t176_g559411 * (cse_p14_t1_g559484 + (cse_p12_t211_g559446 * (-(t1))) + (-((cse_p11_t4_invr559234 * v2_len * (-(c2y))))))) + (cse_p12_t211_g559446 * cse_p13_t1_g559480))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
      double h_y1_z1 = ((cse_p12_t162_g559397 * cse_p12_t163_g559398 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t176_g559411 * cse_p12_t222_g559457) + (cse_p12_t176_g559411 * cse_p12_t212_g559447 * (-(t1))) + (cse_p12_t178_g559413 * cse_p14_t1_g559484) + (cse_p12_t212_g559447 * cse_p13_t1_g559480))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
      double h_y1_x2 = ((cse_p12_t161_g559396 * cse_p12_t162_g559397 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t176_g559411 * cse_p12_t216_g559451) + (cse_p12_t175_g559410 * cse_p14_t1_g559484) + (cse_p12_t176_g559411 * cse_p12_t202_g559437 * (-(t1))) + (cse_p12_t202_g559437 * cse_p13_t1_g559480) + (cse_p13_t1_g559473 * ((v2x * (-(v3y))) + (-(c2z)))) + (cse_p13_t1_g559474 * (-(v3z))) + (-((cse_p13_t1_g559477 * v2x * (-(c2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
      double h_y1_y2 = ((cse_p12_t156_g559391 * cse_p12_t162_g559397 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t176_g559411 * cse_p12_t217_g559452) + (cse_p12_t170_g559405 * cse_p14_t1_g559484) + (cse_p12_t176_g559411 * cse_p12_t203_g559438 * (-(t1))) + (cse_p12_t203_g559438 * cse_p13_t1_g559480) + (cse_p13_t1_g559473 * (((-(v2z)) * (-(v3z))) + (v2x * v3x))) + (-((cse_p13_t1_g559477 * v2y * (-(c2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
      double h_y1_z2 = ((cse_p12_t158_g559393 * cse_p12_t162_g559397 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t176_g559411 * cse_p12_t218_g559453) + (cse_p12_t172_g559407 * cse_p14_t1_g559484) + (cse_p12_t176_g559411 * cse_p12_t204_g559439 * (-(t1))) + (cse_p12_t204_g559439 * cse_p13_t1_g559480) + (cse_p13_t1_g559473 * (c2x + (v3y * (-(v2z))))) + (cse_p13_t1_g559474 * v3x) + (-((cse_p13_t1_g559477 * v2z * (-(c2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
      double h_y1_x3 = ((cse_p12_t159_g559394 * cse_p12_t162_g559397 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t176_g559411 * cse_p12_t219_g559454) + (cse_p12_t173_g559408 * cse_p14_t1_g559484) + (cse_p12_t176_g559411 * cse_p12_t205_g559440 * (-(t1))) + (cse_p12_t200_g559435 * cse_p13_t1_g559474) + (cse_p12_t205_g559440 * cse_p13_t1_g559480) + (cse_p13_t1_g559473 * (c2z + (cse_p12_t199_g559434 * v2x))) + (cse_p13_t1_g559477 * v2x * (-(c2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 0, h_y1_x3);
      double h_y1_y3 = ((cse_p12_t154_g559389 * cse_p12_t162_g559397 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t176_g559411 * cse_p12_t220_g559455) + (cse_p12_t168_g559403 * cse_p14_t1_g559484) + (cse_p12_t176_g559411 * cse_p12_t206_g559441 * (-(t1))) + (cse_p12_t206_g559441 * cse_p13_t1_g559480) + (cse_p13_t1_g559473 * ((cse_p12_t189_g559424 * v2x) + (cse_p12_t200_g559435 * (-(v2z))))) + (cse_p13_t1_g559477 * v2y * (-(c2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 1, h_y1_y3);
      double h_y1_z3 = ((cse_p12_t153_g559388 * cse_p12_t162_g559397 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t176_g559411 * cse_p12_t221_g559456) + (cse_p12_t167_g559402 * cse_p14_t1_g559484) + (cse_p12_t176_g559411 * cse_p14_t1_g559489) + (cse_p12_t189_g559424 * cse_p13_t1_g559474) + (cse_p12_t207_g559442 * cse_p13_t1_g559480) + (cse_p13_t1_g559473 * ((cse_p12_t190_g559425 * (-(v2z))) + (-(c2x)))) + (cse_p13_t1_g559477 * v2z * (-(c2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 2, h_y1_z3);
      double h_y1_x4 = ((cse_p12_t160_g559395 * cse_p12_t162_g559397 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t176_g559411 * cse_p12_t215_g559450) + (cse_p12_t174_g559409 * cse_p14_t1_g559484) + (cse_p12_t176_g559411 * cse_p14_t1_g559488) + (cse_p12_t210_g559445 * cse_p13_t1_g559480) + (cse_p13_t1_g559473 * v2x * (-(v2y))) + (cse_p13_t1_g559474 * (-(v2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x4, 0, h_y1_x4);
      double h_y1_y4 = ((cse_p12_t155_g559390 * cse_p12_t162_g559397 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t176_g559411 * cse_p12_t213_g559448) + (cse_p12_t169_g559404 * cse_p14_t1_g559484) + (cse_p12_t176_g559411 * cse_p14_t1_g559485) + (cse_p12_t208_g559443 * cse_p13_t1_g559480) + (cse_p13_t1_g559473 * (((-(v2z)) * (-(v2z))) + (v2x * v2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x4, 1, h_y1_y4);
      double h_y1_z4 = ((cse_p12_t157_g559392 * cse_p12_t162_g559397 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t176_g559411 * cse_p12_t214_g559449) + (cse_p12_t171_g559406 * cse_p14_t1_g559484) + (cse_p12_t176_g559411 * cse_p14_t1_g559483) + (cse_p12_t209_g559444 * cse_p13_t1_g559480) + (cse_p13_t1_g559473 * v2y * (-(v2z))) + (cse_p13_t1_g559474 * v2x))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x4, 2, h_y1_z4);
      double cse_p14_t1_g559487 = (cse_p13_t1_g559472 * (-(c2z)));
      double h_z1_z1 = ((cse_p12_t163_g559398 * cse_p12_t163_g559398 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p12_t178_g559413 * (cse_p14_t1_g559487 + (cse_p11_t4_invr559234 * cse_p12_t222_g559457) + (cse_p12_t212_g559447 * (-(t1))))) + (cse_p12_t212_g559447 * cse_p14_t1_g559482))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
      double h_z1_x2 = ((cse_p12_t161_g559396 * cse_p12_t163_g559398 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t178_g559413 * cse_p12_t216_g559451) + (cse_p12_t175_g559410 * cse_p14_t1_g559487) + (cse_p12_t178_g559413 * cse_p12_t202_g559437 * (-(t1))) + (cse_p12_t202_g559437 * cse_p14_t1_g559482) + (cse_p13_t1_g559473 * (c2y + (v3z * (-(v2x))))) + (cse_p13_t1_g559474 * v3y) + (-((cse_p13_t1_g559477 * v2x * (-(c2z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
      double h_z1_y2 = ((cse_p12_t156_g559391 * cse_p12_t163_g559398 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t178_g559413 * cse_p12_t217_g559452) + (cse_p12_t170_g559405 * cse_p14_t1_g559487) + (cse_p12_t178_g559413 * cse_p12_t203_g559438 * (-(t1))) + (cse_p12_t203_g559438 * cse_p14_t1_g559482) + (cse_p13_t1_g559473 * ((v2y * (-(v3z))) + (-(c2x)))) + (cse_p13_t1_g559474 * (-(v3x))) + (-((cse_p13_t1_g559477 * v2y * (-(c2z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
      double h_z1_z2 = ((cse_p12_t158_g559393 * cse_p12_t163_g559398 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t178_g559413 * cse_p12_t218_g559453) + (cse_p12_t172_g559407 * cse_p14_t1_g559487) + (cse_p12_t178_g559413 * cse_p12_t204_g559439 * (-(t1))) + (cse_p12_t204_g559439 * cse_p14_t1_g559482) + (cse_p13_t1_g559473 * (((-(v2x)) * (-(v3x))) + (v2y * v3y))) + (-((cse_p13_t1_g559477 * v2z * (-(c2z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
      double h_z1_x3 = ((cse_p12_t159_g559394 * cse_p12_t163_g559398 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t178_g559413 * cse_p12_t219_g559454) + (cse_p12_t173_g559408 * cse_p14_t1_g559487) + (cse_p12_t178_g559413 * cse_p12_t205_g559440 * (-(t1))) + (cse_p12_t190_g559425 * cse_p13_t1_g559474) + (cse_p12_t205_g559440 * cse_p14_t1_g559482) + (cse_p13_t1_g559473 * ((cse_p12_t191_g559426 * (-(v2x))) + (-(c2y)))) + (cse_p13_t1_g559477 * v2x * (-(c2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 0, h_z1_x3);
      double h_z1_y3 = ((cse_p12_t154_g559389 * cse_p12_t163_g559398 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t178_g559413 * cse_p12_t220_g559455) + (cse_p12_t168_g559403 * cse_p14_t1_g559487) + (cse_p12_t178_g559413 * cse_p12_t206_g559441 * (-(t1))) + (cse_p12_t198_g559433 * cse_p13_t1_g559474) + (cse_p12_t206_g559441 * cse_p14_t1_g559482) + (cse_p13_t1_g559473 * (c2x + (cse_p12_t200_g559435 * v2y))) + (cse_p13_t1_g559477 * v2y * (-(c2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 1, h_z1_y3);
      double h_z1_z3 = ((cse_p12_t153_g559388 * cse_p12_t163_g559398 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t178_g559413 * cse_p12_t221_g559456) + (cse_p12_t167_g559402 * cse_p14_t1_g559487) + (cse_p12_t178_g559413 * cse_p14_t1_g559489) + (cse_p12_t207_g559442 * cse_p14_t1_g559482) + (cse_p13_t1_g559473 * ((cse_p12_t190_g559425 * v2y) + (cse_p12_t198_g559433 * (-(v2x))))) + (cse_p13_t1_g559477 * v2z * (-(c2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 2, h_z1_z3);
      double h_z1_x4 = ((cse_p12_t160_g559395 * cse_p12_t163_g559398 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t178_g559413 * cse_p12_t215_g559450) + (cse_p12_t174_g559409 * cse_p14_t1_g559487) + (cse_p12_t178_g559413 * cse_p14_t1_g559488) + (cse_p12_t210_g559445 * cse_p14_t1_g559482) + (cse_p13_t1_g559473 * v2z * (-(v2x))) + (cse_p13_t1_g559474 * v2y))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x4, 0, h_z1_x4);
      double h_z1_y4 = ((cse_p12_t155_g559390 * cse_p12_t163_g559398 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t178_g559413 * cse_p12_t213_g559448) + (cse_p12_t169_g559404 * cse_p14_t1_g559487) + (cse_p12_t178_g559413 * cse_p14_t1_g559485) + (cse_p12_t208_g559443 * cse_p14_t1_g559482) + (cse_p13_t1_g559473 * v2y * (-(v2z))) + (cse_p13_t1_g559474 * (-(v2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x4, 1, h_z1_y4);
      double h_z1_z4 = ((cse_p12_t157_g559392 * cse_p12_t163_g559398 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t178_g559413 * cse_p12_t214_g559449) + (cse_p12_t171_g559406 * cse_p14_t1_g559487) + (cse_p12_t178_g559413 * cse_p14_t1_g559483) + (cse_p12_t209_g559444 * cse_p14_t1_g559482) + (cse_p13_t1_g559473 * (((-(v2x)) * (-(v2x))) + (v2y * v2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x4, 2, h_z1_z4);
      double cse_p14_t1_g559490 = (cse_p11_t2_invsqrt559232 * (-(v2x)));
      double h_x2_x2 = ((cse_p12_t161_g559396 * cse_p12_t161_g559396 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p12_t147_g559382 * cse_p12_t202_g559437 * t2) + (cse_p12_t175_g559410 * ((cse_p11_t4_invr559234 * cse_p12_t147_g559382) + (cse_p11_t4_invr559234 * cse_p12_t216_g559451) + (cse_p12_t202_g559437 * (-(t1))))) + (cse_p12_t1_g559462 * (cse_p12_t1_g559459 + (cse_p12_t192_g559427 * cse_p14_t1_g559490) + (cse_p12_t192_g559427 * cse_p14_t1_g559490) + (-(cse_p12_t66_g559301)))) + (cse_p13_t1_g559473 * ((cse_p12_t188_g559423 * v3z) + (cse_p12_t188_g559423 * v3z) + (cse_p12_t196_g559431 * (-(v3y))) + (cse_p12_t196_g559431 * (-(v3y))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
      double cse_p14_t1_g559491 = (cse_p11_t2_invsqrt559232 * (-(v2y)));
      double h_x2_y2 = ((cse_p12_t156_g559391 * cse_p12_t161_g559396 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t147_g559382 * cse_p12_t170_g559405) + (cse_p11_t4_invr559234 * cse_p12_t175_g559410 * cse_p12_t217_g559452) + (cse_p12_t147_g559382 * cse_p12_t203_g559438 * t2) + (cse_p12_t175_g559410 * cse_p12_t203_g559438 * (-(t1))) + (cse_p12_t1_g559462 * ((cse_p12_t192_g559427 * cse_p14_t1_g559491) + (cse_p12_t193_g559428 * cse_p14_t1_g559490) + (-(cse_p12_t67_g559302)))) + (cse_p13_t1_g559473 * ((cse_p12_t186_g559421 * (-(v3y))) + (cse_p12_t196_g559431 * v3x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
      double cse_p15_t1_g559492 = (cse_p11_t2_invsqrt559232 * (-(v2z)));
      double h_x2_z2 = ((cse_p12_t158_g559393 * cse_p12_t161_g559396 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t147_g559382 * cse_p12_t172_g559407) + (cse_p11_t4_invr559234 * cse_p12_t175_g559410 * cse_p12_t218_g559453) + (cse_p12_t147_g559382 * cse_p12_t204_g559439 * t2) + (cse_p12_t175_g559410 * cse_p12_t204_g559439 * (-(t1))) + (cse_p12_t1_g559462 * ((cse_p12_t192_g559427 * cse_p15_t1_g559492) + (cse_p12_t194_g559429 * cse_p14_t1_g559490) + (-(cse_p12_t68_g559303)))) + (cse_p13_t1_g559473 * ((cse_p12_t188_g559423 * (-(v3x))) + (cse_p12_t195_g559430 * v3z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
      double h_x2_x3 = ((cse_p12_t159_g559394 * cse_p12_t161_g559396 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t147_g559382 * cse_p12_t173_g559408) + (cse_p11_t4_invr559234 * cse_p12_t175_g559410 * cse_p12_t219_g559454) + (cse_p12_t147_g559382 * cse_p12_t205_g559440 * t2) + (cse_p12_t175_g559410 * cse_p12_t205_g559440 * (-(t1))) + (cse_p12_t1_g559462 * (cse_p12_t201_g559436 + cse_p12_t66_g559301 + (cse_p11_t2_invsqrt559232 * cse_p12_t192_g559427 * v2x) + (cse_p12_t184_g559419 * cse_p14_t1_g559490))) + (cse_p13_t1_g559473 * (((-(v1y)) * (-(v3y))) + (cse_p12_t188_g559423 * cse_p12_t191_g559426) + (cse_p12_t196_g559431 * cse_p12_t199_g559434) + (v1z * v3z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 0, h_x2_x3);
      double h_x2_y3 = ((cse_p12_t154_g559389 * cse_p12_t161_g559396 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t147_g559382 * cse_p12_t168_g559403) + (cse_p11_t4_invr559234 * cse_p12_t175_g559410 * cse_p12_t220_g559455) + (cse_p12_t147_g559382 * cse_p12_t206_g559441 * t2) + (cse_p12_t175_g559410 * cse_p12_t206_g559441 * (-(t1))) + (cse_p12_t1_g559462 * (cse_p12_t67_g559302 + (cse_p11_t2_invsqrt559232 * cse_p12_t192_g559427 * v2y) + (cse_p12_t181_g559416 * cse_p14_t1_g559490) + (v2_len * (v1z + v2z + v3z)))) + (cse_p13_t1_g559473 * (c1z + c2z + (cse_p12_t189_g559424 * cse_p12_t196_g559431) + (v1x * (-(v3y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 1, h_x2_y3);
      double h_x2_z3 = ((cse_p12_t153_g559388 * cse_p12_t161_g559396 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t147_g559382 * cse_p12_t167_g559402) + (cse_p11_t4_invr559234 * cse_p12_t175_g559410 * cse_p12_t221_g559456) + (cse_p12_t147_g559382 * cse_p12_t207_g559442 * t2) + (cse_p12_t175_g559410 * cse_p14_t1_g559489) + (cse_p12_t1_g559462 * (cse_p12_t68_g559303 + (cse_p11_t2_invsqrt559232 * cse_p12_t192_g559427 * v2z) + (cse_p12_t180_g559415 * cse_p14_t1_g559490) + (v2_len * ((-(v1y)) + (-(v2y)) + (-(v3y)))))) + (cse_p13_t1_g559473 * ((cse_p12_t188_g559423 * cse_p12_t198_g559433) + (v3z * (-(v1x))) + (-(c1y)) + (-(c2y)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 2, h_x2_z3);
      double h_x2_x4 = ((cse_p12_t160_g559395 * cse_p12_t161_g559396 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t147_g559382 * cse_p12_t174_g559409) + (cse_p11_t4_invr559234 * cse_p12_t175_g559410 * cse_p12_t215_g559450) + (cse_p12_t147_g559382 * cse_p12_t210_g559445 * t2) + (cse_p12_t175_g559410 * cse_p14_t1_g559488) + (cse_p13_t1_g559473 * ((cse_p12_t188_g559423 * v2z) + (cse_p12_t196_g559431 * (-(v2y))))) + (-((cse_p12_t185_g559420 * cse_p13_t1_g559477 * v2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x4, 0, h_x2_x4);
      double h_x2_y4 = ((cse_p12_t155_g559390 * cse_p12_t161_g559396 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t147_g559382 * cse_p12_t169_g559404) + (cse_p11_t4_invr559234 * cse_p12_t175_g559410 * cse_p12_t213_g559448) + (cse_p12_t147_g559382 * cse_p12_t208_g559443 * t2) + (cse_p12_t175_g559410 * cse_p14_t1_g559485) + (cse_p12_t1_g559462 * ((cse_p12_t182_g559417 * cse_p14_t1_g559490) + (cse_p12_t188_g559423 * v2_len))) + (cse_p13_t1_g559473 * ((cse_p12_t196_g559431 * v2x) + (-(c1z)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x4, 1, h_x2_y4);
      double h_x2_z4 = ((cse_p12_t157_g559392 * cse_p12_t161_g559396 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t147_g559382 * cse_p12_t171_g559406) + (cse_p11_t4_invr559234 * cse_p12_t175_g559410 * cse_p12_t214_g559449) + (cse_p12_t147_g559382 * cse_p12_t209_g559444 * t2) + (cse_p12_t175_g559410 * cse_p14_t1_g559483) + (cse_p12_t1_g559462 * ((cse_p12_t183_g559418 * cse_p14_t1_g559490) + (cse_p12_t196_g559431 * v2_len))) + (cse_p13_t1_g559473 * (c1y + (cse_p12_t188_g559423 * (-(v2x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x4, 2, h_x2_z4);
      double h_y2_y2 = ((cse_p12_t156_g559391 * cse_p12_t156_g559391 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p12_t148_g559383 * cse_p12_t203_g559438 * t2) + (cse_p12_t170_g559405 * ((cse_p11_t4_invr559234 * cse_p12_t148_g559383) + (cse_p11_t4_invr559234 * cse_p12_t217_g559452) + (cse_p12_t203_g559438 * (-(t1))))) + (cse_p12_t1_g559462 * (cse_p12_t1_g559459 + (cse_p12_t193_g559428 * cse_p14_t1_g559491) + (cse_p12_t193_g559428 * cse_p14_t1_g559491) + (-(cse_p12_t69_g559304)))) + (cse_p13_t1_g559473 * ((cse_p12_t186_g559421 * v3x) + (cse_p12_t186_g559421 * v3x) + (cse_p12_t197_g559432 * (-(v3z))) + (cse_p12_t197_g559432 * (-(v3z))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
      double h_y2_z2 = ((cse_p12_t156_g559391 * cse_p12_t158_g559393 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t148_g559383 * cse_p12_t172_g559407) + (cse_p11_t4_invr559234 * cse_p12_t170_g559405 * cse_p12_t218_g559453) + (cse_p12_t148_g559383 * cse_p12_t204_g559439 * t2) + (cse_p12_t170_g559405 * cse_p12_t204_g559439 * (-(t1))) + (cse_p12_t1_g559462 * ((cse_p12_t193_g559428 * cse_p15_t1_g559492) + (cse_p12_t194_g559429 * cse_p14_t1_g559491) + (-(cse_p12_t70_g559305)))) + (cse_p13_t1_g559473 * ((cse_p12_t187_g559422 * (-(v3z))) + (cse_p12_t197_g559432 * v3y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
      double h_y2_x3 = ((cse_p12_t156_g559391 * cse_p12_t159_g559394 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t148_g559383 * cse_p12_t173_g559408) + (cse_p11_t4_invr559234 * cse_p12_t170_g559405 * cse_p12_t219_g559454) + (cse_p12_t148_g559383 * cse_p12_t205_g559440 * t2) + (cse_p12_t170_g559405 * cse_p12_t205_g559440 * (-(t1))) + (cse_p12_t1_g559462 * (cse_p12_t67_g559302 + (cse_p11_t2_invsqrt559232 * cse_p12_t193_g559428 * v2x) + (cse_p12_t184_g559419 * cse_p14_t1_g559491) + (v2_len * ((-(v1z)) + (-(v2z)) + (-(v3z)))))) + (cse_p13_t1_g559473 * ((cse_p12_t186_g559421 * cse_p12_t199_g559434) + (v3x * (-(v1y))) + (-(c1z)) + (-(c2z)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 0, h_y2_x3);
      double h_y2_y3 = ((cse_p12_t154_g559389 * cse_p12_t156_g559391 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t148_g559383 * cse_p12_t168_g559403) + (cse_p11_t4_invr559234 * cse_p12_t170_g559405 * cse_p12_t220_g559455) + (cse_p12_t148_g559383 * cse_p12_t206_g559441 * t2) + (cse_p12_t170_g559405 * cse_p12_t206_g559441 * (-(t1))) + (cse_p12_t1_g559462 * (cse_p12_t201_g559436 + cse_p12_t69_g559304 + (cse_p11_t2_invsqrt559232 * cse_p12_t193_g559428 * v2y) + (cse_p12_t181_g559416 * cse_p14_t1_g559491))) + (cse_p13_t1_g559473 * (((-(v1z)) * (-(v3z))) + (cse_p12_t186_g559421 * cse_p12_t189_g559424) + (cse_p12_t197_g559432 * cse_p12_t200_g559435) + (v1x * v3x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 1, h_y2_y3);
      double h_y2_z3 = ((cse_p12_t153_g559388 * cse_p12_t156_g559391 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t148_g559383 * cse_p12_t167_g559402) + (cse_p11_t4_invr559234 * cse_p12_t170_g559405 * cse_p12_t221_g559456) + (cse_p12_t148_g559383 * cse_p12_t207_g559442 * t2) + (cse_p12_t170_g559405 * cse_p14_t1_g559489) + (cse_p12_t1_g559462 * (cse_p12_t70_g559305 + (cse_p11_t2_invsqrt559232 * cse_p12_t193_g559428 * v2z) + (cse_p12_t180_g559415 * cse_p14_t1_g559491) + (v2_len * (v1x + v2x + v3x)))) + (cse_p13_t1_g559473 * (c1x + c2x + (cse_p12_t190_g559425 * cse_p12_t197_g559432) + (v1y * (-(v3z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 2, h_y2_z3);
      double h_y2_x4 = ((cse_p12_t156_g559391 * cse_p12_t160_g559395 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t148_g559383 * cse_p12_t174_g559409) + (cse_p11_t4_invr559234 * cse_p12_t170_g559405 * cse_p12_t215_g559450) + (cse_p12_t148_g559383 * cse_p12_t210_g559445 * t2) + (cse_p12_t170_g559405 * cse_p14_t1_g559488) + (cse_p12_t1_g559462 * ((cse_p12_t185_g559420 * cse_p14_t1_g559491) + (cse_p12_t197_g559432 * v2_len))) + (cse_p13_t1_g559473 * (c1z + (cse_p12_t186_g559421 * (-(v2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x4, 0, h_y2_x4);
      double h_y2_y4 = ((cse_p12_t155_g559390 * cse_p12_t156_g559391 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t148_g559383 * cse_p12_t169_g559404) + (cse_p11_t4_invr559234 * cse_p12_t170_g559405 * cse_p12_t213_g559448) + (cse_p12_t148_g559383 * cse_p12_t208_g559443 * t2) + (cse_p12_t170_g559405 * cse_p14_t1_g559485) + (cse_p13_t1_g559473 * ((cse_p12_t186_g559421 * v2x) + (cse_p12_t197_g559432 * (-(v2z))))) + (-((cse_p12_t182_g559417 * cse_p13_t1_g559477 * v2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x4, 1, h_y2_y4);
      double h_y2_z4 = ((cse_p12_t156_g559391 * cse_p12_t157_g559392 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t148_g559383 * cse_p12_t171_g559406) + (cse_p11_t4_invr559234 * cse_p12_t170_g559405 * cse_p12_t214_g559449) + (cse_p12_t148_g559383 * cse_p12_t209_g559444 * t2) + (cse_p12_t170_g559405 * cse_p14_t1_g559483) + (cse_p12_t1_g559462 * ((cse_p12_t183_g559418 * cse_p14_t1_g559491) + (cse_p12_t186_g559421 * v2_len))) + (cse_p13_t1_g559473 * ((cse_p12_t197_g559432 * v2y) + (-(c1x)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x4, 2, h_y2_z4);
      double h_z2_z2 = ((cse_p12_t158_g559393 * cse_p12_t158_g559393 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p12_t149_g559384 * cse_p12_t204_g559439 * t2) + (cse_p12_t172_g559407 * ((cse_p11_t4_invr559234 * cse_p12_t149_g559384) + (cse_p11_t4_invr559234 * cse_p12_t218_g559453) + (cse_p12_t204_g559439 * (-(t1))))) + (cse_p12_t1_g559462 * (cse_p12_t1_g559459 + (cse_p12_t194_g559429 * cse_p15_t1_g559492) + (cse_p12_t194_g559429 * cse_p15_t1_g559492) + (-(cse_p12_t71_g559306)))) + (cse_p13_t1_g559473 * ((cse_p12_t187_g559422 * v3y) + (cse_p12_t187_g559422 * v3y) + (cse_p12_t195_g559430 * (-(v3x))) + (cse_p12_t195_g559430 * (-(v3x))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      double h_z2_x3 = ((cse_p12_t158_g559393 * cse_p12_t159_g559394 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t149_g559384 * cse_p12_t173_g559408) + (cse_p11_t4_invr559234 * cse_p12_t172_g559407 * cse_p12_t219_g559454) + (cse_p12_t149_g559384 * cse_p12_t205_g559440 * t2) + (cse_p12_t172_g559407 * cse_p12_t205_g559440 * (-(t1))) + (cse_p12_t1_g559462 * (cse_p12_t68_g559303 + (cse_p11_t2_invsqrt559232 * cse_p12_t194_g559429 * v2x) + (cse_p12_t184_g559419 * cse_p15_t1_g559492) + (v2_len * (v1y + v2y + v3y)))) + (cse_p13_t1_g559473 * (c1y + c2y + (cse_p12_t191_g559426 * cse_p12_t195_g559430) + (v1z * (-(v3x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 0, h_z2_x3);
      double h_z2_y3 = ((cse_p12_t154_g559389 * cse_p12_t158_g559393 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t149_g559384 * cse_p12_t168_g559403) + (cse_p11_t4_invr559234 * cse_p12_t172_g559407 * cse_p12_t220_g559455) + (cse_p12_t149_g559384 * cse_p12_t206_g559441 * t2) + (cse_p12_t172_g559407 * cse_p12_t206_g559441 * (-(t1))) + (cse_p12_t1_g559462 * (cse_p12_t70_g559305 + (cse_p11_t2_invsqrt559232 * cse_p12_t194_g559429 * v2y) + (cse_p12_t181_g559416 * cse_p15_t1_g559492) + (v2_len * ((-(v1x)) + (-(v2x)) + (-(v3x)))))) + (cse_p13_t1_g559473 * ((cse_p12_t187_g559422 * cse_p12_t200_g559435) + (v3y * (-(v1z))) + (-(c1x)) + (-(c2x)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 1, h_z2_y3);
      double h_z2_z3 = ((cse_p12_t153_g559388 * cse_p12_t158_g559393 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t149_g559384 * cse_p12_t167_g559402) + (cse_p11_t4_invr559234 * cse_p12_t172_g559407 * cse_p12_t221_g559456) + (cse_p12_t149_g559384 * cse_p12_t207_g559442 * t2) + (cse_p12_t172_g559407 * cse_p14_t1_g559489) + (cse_p12_t1_g559462 * (cse_p12_t201_g559436 + cse_p12_t71_g559306 + (cse_p11_t2_invsqrt559232 * cse_p12_t194_g559429 * v2z) + (cse_p12_t180_g559415 * cse_p15_t1_g559492))) + (cse_p13_t1_g559473 * (((-(v1x)) * (-(v3x))) + (cse_p12_t187_g559422 * cse_p12_t190_g559425) + (cse_p12_t195_g559430 * cse_p12_t198_g559433) + (v1y * v3y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 2, h_z2_z3);
      double h_z2_x4 = ((cse_p12_t158_g559393 * cse_p12_t160_g559395 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t149_g559384 * cse_p12_t174_g559409) + (cse_p11_t4_invr559234 * cse_p12_t172_g559407 * cse_p12_t215_g559450) + (cse_p12_t149_g559384 * cse_p12_t210_g559445 * t2) + (cse_p12_t172_g559407 * cse_p14_t1_g559488) + (cse_p12_t1_g559462 * ((cse_p12_t185_g559420 * cse_p15_t1_g559492) + (cse_p12_t187_g559422 * v2_len))) + (cse_p13_t1_g559473 * ((cse_p12_t195_g559430 * v2z) + (-(c1y)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x4, 0, h_z2_x4);
      double h_z2_y4 = ((cse_p12_t155_g559390 * cse_p12_t158_g559393 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t149_g559384 * cse_p12_t169_g559404) + (cse_p11_t4_invr559234 * cse_p12_t172_g559407 * cse_p12_t213_g559448) + (cse_p12_t149_g559384 * cse_p12_t208_g559443 * t2) + (cse_p12_t172_g559407 * cse_p14_t1_g559485) + (cse_p12_t1_g559462 * ((cse_p12_t182_g559417 * cse_p15_t1_g559492) + (cse_p12_t195_g559430 * v2_len))) + (cse_p13_t1_g559473 * (c1x + (cse_p12_t187_g559422 * (-(v2z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x4, 1, h_z2_y4);
      double h_z2_z4 = ((cse_p12_t157_g559392 * cse_p12_t158_g559393 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t149_g559384 * cse_p12_t171_g559406) + (cse_p11_t4_invr559234 * cse_p12_t172_g559407 * cse_p12_t214_g559449) + (cse_p12_t149_g559384 * cse_p12_t209_g559444 * t2) + (cse_p12_t172_g559407 * cse_p14_t1_g559483) + (cse_p13_t1_g559473 * ((cse_p12_t187_g559422 * v2y) + (cse_p12_t195_g559430 * (-(v2x))))) + (-((cse_p12_t183_g559418 * cse_p13_t1_g559477 * v2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x4, 2, h_z2_z4);
      double h_x3_x3 = ((cse_p12_t159_g559394 * cse_p12_t159_g559394 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p12_t150_g559385 * cse_p12_t205_g559440 * t2) + (cse_p12_t173_g559408 * ((cse_p11_t4_invr559234 * cse_p12_t150_g559385) + (cse_p11_t4_invr559234 * cse_p12_t219_g559454) + (cse_p12_t205_g559440 * (-(t1))))) + (cse_p12_t1_g559462 * (cse_p12_t1_g559459 + (cse_p11_t2_invsqrt559232 * cse_p12_t184_g559419 * v2x) + (cse_p11_t2_invsqrt559232 * cse_p12_t184_g559419 * v2x) + (-(cse_p12_t66_g559301)))) + (cse_p13_t1_g559473 * ((cse_p12_t191_g559426 * v1z) + (cse_p12_t191_g559426 * v1z) + (cse_p12_t199_g559434 * (-(v1y))) + (cse_p12_t199_g559434 * (-(v1y))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 0, h_x3_x3);
      double h_x3_y3 = ((cse_p12_t154_g559389 * cse_p12_t159_g559394 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t150_g559385 * cse_p12_t168_g559403) + (cse_p11_t4_invr559234 * cse_p12_t173_g559408 * cse_p12_t220_g559455) + (cse_p12_t150_g559385 * cse_p12_t206_g559441 * t2) + (cse_p12_t173_g559408 * cse_p12_t206_g559441 * (-(t1))) + (cse_p12_t1_g559462 * ((cse_p11_t2_invsqrt559232 * cse_p12_t181_g559416 * v2x) + (cse_p11_t2_invsqrt559232 * cse_p12_t184_g559419 * v2y) + (-(cse_p12_t67_g559302)))) + (cse_p13_t1_g559473 * ((cse_p12_t189_g559424 * (-(v1y))) + (cse_p12_t199_g559434 * v1x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 1, h_x3_y3);
      double h_x3_z3 = ((cse_p12_t153_g559388 * cse_p12_t159_g559394 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t150_g559385 * cse_p12_t167_g559402) + (cse_p11_t4_invr559234 * cse_p12_t173_g559408 * cse_p12_t221_g559456) + (cse_p12_t150_g559385 * cse_p12_t207_g559442 * t2) + (cse_p12_t173_g559408 * cse_p14_t1_g559489) + (cse_p12_t1_g559462 * ((cse_p11_t2_invsqrt559232 * cse_p12_t180_g559415 * v2x) + (cse_p11_t2_invsqrt559232 * cse_p12_t184_g559419 * v2z) + (-(cse_p12_t68_g559303)))) + (cse_p13_t1_g559473 * ((cse_p12_t191_g559426 * (-(v1x))) + (cse_p12_t198_g559433 * v1z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 2, h_x3_z3);
      double h_x3_x4 = ((cse_p12_t159_g559394 * cse_p12_t160_g559395 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t150_g559385 * cse_p12_t174_g559409) + (cse_p11_t4_invr559234 * cse_p12_t173_g559408 * cse_p12_t215_g559450) + (cse_p12_t150_g559385 * cse_p12_t210_g559445 * t2) + (cse_p12_t173_g559408 * cse_p14_t1_g559488) + (cse_p12_t185_g559420 * cse_p13_t1_g559477 * v2x) + (cse_p13_t1_g559473 * (((-(v1y)) * (-(v2y))) + (v1z * v2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x4, 0, h_x3_x4);
      double h_x3_y4 = ((cse_p12_t155_g559390 * cse_p12_t159_g559394 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t150_g559385 * cse_p12_t169_g559404) + (cse_p11_t4_invr559234 * cse_p12_t173_g559408 * cse_p12_t213_g559448) + (cse_p12_t150_g559385 * cse_p12_t208_g559443 * t2) + (cse_p12_t173_g559408 * cse_p14_t1_g559485) + (cse_p12_t1_g559462 * ((cse_p11_t2_invsqrt559232 * cse_p12_t182_g559417 * v2x) + (v1z * v2_len))) + (cse_p13_t1_g559473 * (c1z + (v2x * (-(v1y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x4, 1, h_x3_y4);
      double h_x3_z4 = ((cse_p12_t157_g559392 * cse_p12_t159_g559394 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t150_g559385 * cse_p12_t171_g559406) + (cse_p11_t4_invr559234 * cse_p12_t173_g559408 * cse_p12_t214_g559449) + (cse_p12_t150_g559385 * cse_p12_t209_g559444 * t2) + (cse_p12_t173_g559408 * cse_p14_t1_g559483) + (cse_p12_t1_g559462 * ((cse_p11_t2_invsqrt559232 * cse_p12_t183_g559418 * v2x) + (v2_len * (-(v1y))))) + (cse_p13_t1_g559473 * ((v1z * (-(v2x))) + (-(c1y)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x4, 2, h_x3_z4);
      double h_y3_y3 = ((cse_p12_t154_g559389 * cse_p12_t154_g559389 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p12_t151_g559386 * cse_p12_t206_g559441 * t2) + (cse_p12_t168_g559403 * ((cse_p11_t4_invr559234 * cse_p12_t151_g559386) + (cse_p11_t4_invr559234 * cse_p12_t220_g559455) + (cse_p12_t206_g559441 * (-(t1))))) + (cse_p12_t1_g559462 * (cse_p12_t1_g559459 + (cse_p11_t2_invsqrt559232 * cse_p12_t181_g559416 * v2y) + (cse_p11_t2_invsqrt559232 * cse_p12_t181_g559416 * v2y) + (-(cse_p12_t69_g559304)))) + (cse_p13_t1_g559473 * ((cse_p12_t189_g559424 * v1x) + (cse_p12_t189_g559424 * v1x) + (cse_p12_t200_g559435 * (-(v1z))) + (cse_p12_t200_g559435 * (-(v1z))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 1, h_y3_y3);
      double h_y3_z3 = ((cse_p12_t153_g559388 * cse_p12_t154_g559389 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t151_g559386 * cse_p12_t167_g559402) + (cse_p11_t4_invr559234 * cse_p12_t168_g559403 * cse_p12_t221_g559456) + (cse_p12_t151_g559386 * cse_p12_t207_g559442 * t2) + (cse_p12_t168_g559403 * cse_p14_t1_g559489) + (cse_p12_t1_g559462 * ((cse_p11_t2_invsqrt559232 * cse_p12_t180_g559415 * v2y) + (cse_p11_t2_invsqrt559232 * cse_p12_t181_g559416 * v2z) + (-(cse_p12_t70_g559305)))) + (cse_p13_t1_g559473 * ((cse_p12_t190_g559425 * (-(v1z))) + (cse_p12_t200_g559435 * v1y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 2, h_y3_z3);
      double h_y3_x4 = ((cse_p12_t154_g559389 * cse_p12_t160_g559395 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t151_g559386 * cse_p12_t174_g559409) + (cse_p11_t4_invr559234 * cse_p12_t168_g559403 * cse_p12_t215_g559450) + (cse_p12_t151_g559386 * cse_p12_t210_g559445 * t2) + (cse_p12_t168_g559403 * cse_p14_t1_g559488) + (cse_p12_t1_g559462 * ((cse_p11_t2_invsqrt559232 * cse_p12_t185_g559420 * v2y) + (v2_len * (-(v1z))))) + (cse_p13_t1_g559473 * ((v1x * (-(v2y))) + (-(c1z)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x4, 0, h_y3_x4);
      double h_y3_y4 = ((cse_p12_t154_g559389 * cse_p12_t155_g559390 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t151_g559386 * cse_p12_t169_g559404) + (cse_p11_t4_invr559234 * cse_p12_t168_g559403 * cse_p12_t213_g559448) + (cse_p12_t151_g559386 * cse_p12_t208_g559443 * t2) + (cse_p12_t168_g559403 * cse_p14_t1_g559485) + (cse_p12_t182_g559417 * cse_p13_t1_g559477 * v2y) + (cse_p13_t1_g559473 * (((-(v1z)) * (-(v2z))) + (v1x * v2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x4, 1, h_y3_y4);
      double h_y3_z4 = ((cse_p12_t154_g559389 * cse_p12_t157_g559392 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t151_g559386 * cse_p12_t171_g559406) + (cse_p11_t4_invr559234 * cse_p12_t168_g559403 * cse_p12_t214_g559449) + (cse_p12_t151_g559386 * cse_p12_t209_g559444 * t2) + (cse_p12_t168_g559403 * cse_p14_t1_g559483) + (cse_p12_t1_g559462 * ((cse_p11_t2_invsqrt559232 * cse_p12_t183_g559418 * v2y) + (v1x * v2_len))) + (cse_p13_t1_g559473 * (c1x + (v2y * (-(v1z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x4, 2, h_y3_z4);
      double h_z3_z3 = ((cse_p12_t153_g559388 * cse_p12_t153_g559388 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p12_t152_g559387 * cse_p12_t207_g559442 * t2) + (cse_p12_t167_g559402 * (cse_p14_t1_g559489 + (cse_p11_t4_invr559234 * cse_p12_t152_g559387) + (cse_p11_t4_invr559234 * cse_p12_t221_g559456))) + (cse_p12_t1_g559462 * (cse_p12_t1_g559459 + (cse_p11_t2_invsqrt559232 * cse_p12_t180_g559415 * v2z) + (cse_p11_t2_invsqrt559232 * cse_p12_t180_g559415 * v2z) + (-(cse_p12_t71_g559306)))) + (cse_p13_t1_g559473 * ((cse_p12_t190_g559425 * v1y) + (cse_p12_t190_g559425 * v1y) + (cse_p12_t198_g559433 * (-(v1x))) + (cse_p12_t198_g559433 * (-(v1x))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x3, 2, h_z3_z3);
      double h_z3_x4 = ((cse_p12_t153_g559388 * cse_p12_t160_g559395 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t152_g559387 * cse_p12_t174_g559409) + (cse_p11_t4_invr559234 * cse_p12_t167_g559402 * cse_p12_t215_g559450) + (cse_p12_t152_g559387 * cse_p12_t210_g559445 * t2) + (cse_p12_t167_g559402 * cse_p14_t1_g559488) + (cse_p12_t1_g559462 * ((cse_p11_t2_invsqrt559232 * cse_p12_t185_g559420 * v2z) + (v1y * v2_len))) + (cse_p13_t1_g559473 * (c1y + (v2z * (-(v1x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x4, 0, h_z3_x4);
      double h_z3_y4 = ((cse_p12_t153_g559388 * cse_p12_t155_g559390 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t152_g559387 * cse_p12_t169_g559404) + (cse_p11_t4_invr559234 * cse_p12_t167_g559402 * cse_p12_t213_g559448) + (cse_p12_t152_g559387 * cse_p12_t208_g559443 * t2) + (cse_p12_t167_g559402 * cse_p14_t1_g559485) + (cse_p12_t1_g559462 * ((cse_p11_t2_invsqrt559232 * cse_p12_t182_g559417 * v2z) + (v2_len * (-(v1x))))) + (cse_p13_t1_g559473 * ((v1y * (-(v2z))) + (-(c1x)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x4, 1, h_z3_y4);
      double h_z3_z4 = ((cse_p12_t153_g559388 * cse_p12_t157_g559392 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t152_g559387 * cse_p12_t171_g559406) + (cse_p11_t4_invr559234 * cse_p12_t167_g559402 * cse_p12_t214_g559449) + (cse_p12_t152_g559387 * cse_p12_t209_g559444 * t2) + (cse_p12_t167_g559402 * cse_p14_t1_g559483) + (cse_p12_t183_g559418 * cse_p13_t1_g559477 * v2z) + (cse_p13_t1_g559473 * (((-(v1x)) * (-(v2x))) + (v1y * v2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x4, 2, h_z3_z4);
      double h_x4_x4 = ((cse_p12_t160_g559395 * cse_p12_t160_g559395 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p12_t174_g559409 * (cse_p14_t1_g559488 + (cse_p11_t4_invr559234 * cse_p12_t215_g559450) + (cse_p12_t185_g559420 * cse_p13_t1_g559472))) + (cse_p12_t185_g559420 * cse_p12_t1_g559465 * cse_p12_t210_g559445))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 0, term.i3x4, 0, h_x4_x4);
      double h_x4_y4 = ((cse_p12_t155_g559390 * cse_p12_t160_g559395 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t174_g559409 * cse_p12_t213_g559448) + (cse_p12_t169_g559404 * cse_p12_t185_g559420 * cse_p13_t1_g559472) + (cse_p12_t174_g559409 * cse_p14_t1_g559485) + (cse_p12_t185_g559420 * cse_p12_t1_g559465 * cse_p12_t208_g559443))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 0, term.i3x4, 1, h_x4_y4);
      double h_x4_z4 = ((cse_p12_t157_g559392 * cse_p12_t160_g559395 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t174_g559409 * cse_p12_t214_g559449) + (cse_p12_t171_g559406 * cse_p12_t185_g559420 * cse_p13_t1_g559472) + (cse_p12_t174_g559409 * cse_p14_t1_g559483) + (cse_p12_t185_g559420 * cse_p12_t1_g559465 * cse_p12_t209_g559444))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 0, term.i3x4, 2, h_x4_z4);
      double h_y4_y4 = ((cse_p12_t155_g559390 * cse_p12_t155_g559390 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p12_t169_g559404 * (cse_p14_t1_g559485 + (cse_p11_t4_invr559234 * cse_p12_t213_g559448) + (cse_p12_t182_g559417 * cse_p13_t1_g559472))) + (cse_p12_t182_g559417 * cse_p12_t1_g559465 * cse_p12_t208_g559443))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 1, term.i3x4, 1, h_y4_y4);
      double h_y4_z4 = ((cse_p12_t155_g559390 * cse_p12_t157_g559392 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p11_t4_invr559234 * cse_p12_t169_g559404 * cse_p12_t214_g559449) + (cse_p12_t169_g559404 * cse_p14_t1_g559483) + (cse_p12_t171_g559406 * cse_p12_t182_g559417 * cse_p13_t1_g559472) + (cse_p12_t182_g559417 * cse_p12_t1_g559465 * cse_p12_t209_g559444))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 1, term.i3x4, 2, h_y4_z4);
      double h_z4_z4 = ((cse_p12_t157_g559392 * cse_p12_t157_g559392 * cse_p12_t1_g559463) + (cse_p12_t1_g559464 * ((cse_p12_t171_g559406 * (cse_p14_t1_g559483 + (cse_p11_t4_invr559234 * cse_p12_t214_g559449) + (cse_p12_t183_g559418 * cse_p13_t1_g559472))) + (cse_p12_t183_g559418 * cse_p12_t1_g559465 * cse_p12_t209_g559444))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 2, term.i3x4, 2, h_z4_z4);
    }
  }
  return energy_added;
}
void hessian_fd(const dihedral_term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 12;
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
    double saved = position[term.i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x3 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x3 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x3 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x3, 0, d);
  }
  {
    double saved = position[term.i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x3 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x3 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x3 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x3, 1, d);
  }
  {
    double saved = position[term.i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x3 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x3 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x3 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x3, 2, d);
  }
  {
    double saved = position[term.i3x4 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x4 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x4 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x4 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x4, 0, d);
  }
  {
    double saved = position[term.i3x4 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x4 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x4 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x4 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x4, 1, d);
  }
  {
    double saved = position[term.i3x4 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x4 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x4 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x4 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(term.i3x4, 2, d);
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
    double saved = position[term.i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x3 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x3 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x3 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x3, 0, hval);
  }
  {
    double saved = position[term.i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x3 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x3 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x3 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x3, 1, hval);
  }
  {
    double saved = position[term.i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x3 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x3 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x3 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x3, 2, hval);
  }
  {
    double saved = position[term.i3x4 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x4 + 0] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x4 + 0] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x4 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x4, 0, hval);
  }
  {
    double saved = position[term.i3x4 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x4 + 1] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x4 + 1] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x4 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x4, 1, hval);
  }
  {
    double saved = position[term.i3x4 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[term.i3x4 + 2] = saved + h;
    energy(term, position, &e_plus);
    position[term.i3x4 + 2] = saved - h;
    energy(term, position, &e_minus);
    position[term.i3x4 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x4, 2, hval);
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
  {
    double saved_i = position[term.i3x3 + 0];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 0] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 0] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 0] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 0];
    double saved_j = position[term.i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 0] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 0] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 0] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 0];
    double saved_j = position[term.i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 0] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 0] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 0] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x3 + 0];
    double saved_j = position[term.i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 0] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 0] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 0] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 0];
    double saved_j = position[term.i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 0] = saved_i + h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 0] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 0] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x2, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 0];
    double saved_j = position[term.i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 0] = saved_i + h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 0] = saved_i - h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 0] = saved_i; position[term.i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x2, 2, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x2, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x2, 2, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x3 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x3 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x3, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x2, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x2, 2, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x3 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x3 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x3, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x3 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x3 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x3 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x3 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x3, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 0];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 0] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 0] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 0] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 0];
    double saved_j = position[term.i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 0] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 0] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 0] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 0];
    double saved_j = position[term.i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 0] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 0] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 0] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x4 + 0];
    double saved_j = position[term.i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 0] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 0] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 0] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 0];
    double saved_j = position[term.i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 0] = saved_i + h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 0] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 0] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x2, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 0];
    double saved_j = position[term.i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 0] = saved_i + h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 0] = saved_i - h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 0] = saved_i; position[term.i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x2, 2, hval);
  }
  {
    double saved_i = position[term.i3x4 + 0];
    double saved_j = position[term.i3x3 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 0] = saved_i + h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 0] = saved_i - h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 0] = saved_i; position[term.i3x3 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x3, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 0];
    double saved_j = position[term.i3x3 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 0] = saved_i + h; position[term.i3x3 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 0] = saved_i - h; position[term.i3x3 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 0] = saved_i; position[term.i3x3 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x3, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 0];
    double saved_j = position[term.i3x3 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 0] = saved_i + h; position[term.i3x3 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 0] = saved_i - h; position[term.i3x3 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 0] = saved_i; position[term.i3x3 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 0, term.i3x3, 2, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x2, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x2, 2, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x3 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x3 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x3, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x3 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x3 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x3 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x3 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x3, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x3 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x3 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x3 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x3 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x3, 2, hval);
  }
  {
    double saved_i = position[term.i3x4 + 1];
    double saved_j = position[term.i3x4 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 1] = saved_i + h; position[term.i3x4 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x4 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 1] = saved_i - h; position[term.i3x4 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x4 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 1] = saved_i; position[term.i3x4 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 1, term.i3x4, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x1 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x1 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x1 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x2 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x2 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x2, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x2 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x2 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x2, 2, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x3 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x3 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x3, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x3 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x3 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x3 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x3 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x3, 1, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x3 + 2];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x3 + 2] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x3 + 2] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x3 + 2] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x3 + 2] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x3 + 2] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x3, 2, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x4 + 0];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x4 + 0] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x4 + 0] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x4 + 0] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x4 + 0] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x4 + 0] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x4, 0, hval);
  }
  {
    double saved_i = position[term.i3x4 + 2];
    double saved_j = position[term.i3x4 + 1];
    double e_pp = 0.0;
    double e_pm = 0.0;
    double e_mp = 0.0;
    double e_mm = 0.0;
    position[term.i3x4 + 2] = saved_i + h; position[term.i3x4 + 1] = saved_j + h;
    energy(term, position, &e_pp);
    position[term.i3x4 + 1] = saved_j - h;
    energy(term, position, &e_pm);
    position[term.i3x4 + 2] = saved_i - h; position[term.i3x4 + 1] = saved_j + h;
    energy(term, position, &e_mp);
    position[term.i3x4 + 1] = saved_j - h;
    energy(term, position, &e_mm);
    position[term.i3x4 + 2] = saved_i; position[term.i3x4 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x4, 2, term.i3x4, 1, hval);
  }
}

};
