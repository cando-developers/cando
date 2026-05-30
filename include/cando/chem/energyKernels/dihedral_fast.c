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
      double cse_p12_t6_g372198 = (v1z * v2x);
      double v2y = (y3 + (-(y2)));
      double cse_p12_t4_g372196 = (v1x * v2y);
      double v2z = (z3 + (-(z2)));
      double cse_p12_t5_g372197 = (v1y * v2z);
      double v3x = (x4 + (-(x3)));
      double cse_p12_t8_g372200 = ((-(v2x)) + (-(v3x)));
      double cse_p12_t11_g372203 = (v2x + v3x);
      double v3y = (y4 + (-(y3)));
      double cse_p12_t9_g372201 = ((-(v2y)) + (-(v3y)));
      double cse_p12_t12_g372204 = (v2y + v3y);
      double v3z = (z4 + (-(z3)));
      double cse_p12_t10_g372202 = ((-(v2z)) + (-(v3z)));
      double cse_p12_t13_g372205 = (v2z + v3z);
      double c1x = (cse_p12_t5_g372197 + (-((v1z * v2y))));
      double c1y = (cse_p12_t6_g372198 + (-((v1x * v2z))));
      double c1z = (cse_p12_t4_g372196 + (-((v1y * v2x))));
      double c2x = ((v2y * v3z) + (-((v2z * v3y))));
      double cse_p12_t1_g372193 = (c2x * v1x);
      double c2y = ((v2z * v3x) + (-((v2x * v3z))));
      double cse_p12_t2_g372194 = (c2y * v1y);
      double c2z = ((v2x * v3y) + (-((v2y * v3x))));
      double cse_p12_t3_g372195 = (c2z * v1z);
      double cse_p12_t7_g372199 = (cse_p12_t1_g372193 + cse_p12_t2_g372194 + cse_p12_t3_g372195);
      double v2_sq = ((v2x * v2x) + (v2y * v2y) + (v2z * v2z));
      double cse_p1_t1_sqrt372190 = sqrt(v2_sq);
      double cse_p1_t2_invsqrt372191 = (1.0 / (cse_p1_t1_sqrt372190));
      double v2_len = cse_p1_t1_sqrt372190;
      double t1 = (cse_p12_t7_g372199 * v2_len);
      double t2 = ((c1x * c2x) + (c1y * c2y) + (c1z * c2z));
      double phi = atan2(t1, t2);
      double nphi = (term.n * phi);
      double cse_p12_t14_g372206 = (-(sin(nphi)));
      double sin_nphi = sin(nphi);
      double cos_nphi = cos(nphi);
      double cos_angle = ((term.cosphase * cos_nphi) + (term.sinphase * sin_nphi));
      double sin_angle = ((term.cosphase * sin_nphi) + (-((cos_nphi * term.sinphase))));
      double energy = (term.v * (1.00000000000000000e+0 + cos_angle));
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p1_t3_invr372192 = (1.0 / ((((t1) * (t1)) + ((t2) * (t2)))));
      double cse_p12_t1_g372207 = (cse_p1_t3_invr372192 * (-(t1)));
      double cse_p12_t1_g372208 = (cse_p1_t3_invr372192 * t2);
      double cse_p12_t1_g372209 = (cos_nphi * term.n * term.sinphase);
      double cse_p12_t1_g372210 = (term.cosphase * cse_p12_t14_g372206 * term.n);
      double cse_p12_t1_g372211 = (cse_p12_t1_g372208 * v2_len);
      double g_x1 = (term.v * ((cse_p12_t1_g372209 * ((cse_p12_t1_g372207 * ((c2y * v2z) + (c2z * (-(v2y))))) + (cse_p12_t1_g372211 * (-(c2x))))) + (cse_p12_t1_g372210 * ((cse_p12_t1_g372207 * ((c2y * v2z) + (c2z * (-(v2y))))) + (cse_p12_t1_g372211 * (-(c2x)))))));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = (term.v * ((cse_p12_t1_g372209 * ((cse_p12_t1_g372207 * ((c2x * (-(v2z))) + (c2z * v2x))) + (cse_p12_t1_g372211 * (-(c2y))))) + (cse_p12_t1_g372210 * ((cse_p12_t1_g372207 * ((c2x * (-(v2z))) + (c2z * v2x))) + (cse_p12_t1_g372211 * (-(c2y)))))));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = (term.v * ((cse_p12_t1_g372209 * ((cse_p12_t1_g372207 * ((c2x * v2y) + (c2y * (-(v2x))))) + (cse_p12_t1_g372211 * (-(c2z))))) + (cse_p12_t1_g372210 * ((cse_p12_t1_g372207 * ((c2x * v2y) + (c2y * (-(v2x))))) + (cse_p12_t1_g372211 * (-(c2z)))))));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p12_t1_g372212 = (cse_p12_t7_g372199 * cse_p1_t2_invsqrt372191);
      double g_x2 = (term.v * ((cse_p12_t1_g372209 * ((cse_p12_t1_g372207 * ((c1y * v3z) + (c1z * (-(v3y))) + (c2y * ((-(v1z)) + (-(v2z)))) + (c2z * (v1y + v2y)))) + (cse_p12_t1_g372208 * ((cse_p12_t1_g372212 * (-(v2x))) + (v2_len * (c2x + (v1y * v3z) + (v1z * (-(v3y))))))))) + (cse_p12_t1_g372210 * ((cse_p12_t1_g372207 * ((c1y * v3z) + (c1z * (-(v3y))) + (c2y * ((-(v1z)) + (-(v2z)))) + (c2z * (v1y + v2y)))) + (cse_p12_t1_g372208 * ((cse_p12_t1_g372212 * (-(v2x))) + (v2_len * (c2x + (v1y * v3z) + (v1z * (-(v3y)))))))))));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (term.v * ((cse_p12_t1_g372209 * ((cse_p12_t1_g372207 * ((c1x * (-(v3z))) + (c1z * v3x) + (c2x * (v1z + v2z)) + (c2z * ((-(v1x)) + (-(v2x)))))) + (cse_p12_t1_g372208 * ((cse_p12_t1_g372212 * (-(v2y))) + (v2_len * (c2y + (v1x * (-(v3z))) + (v1z * v3x))))))) + (cse_p12_t1_g372210 * ((cse_p12_t1_g372207 * ((c1x * (-(v3z))) + (c1z * v3x) + (c2x * (v1z + v2z)) + (c2z * ((-(v1x)) + (-(v2x)))))) + (cse_p12_t1_g372208 * ((cse_p12_t1_g372212 * (-(v2y))) + (v2_len * (c2y + (v1x * (-(v3z))) + (v1z * v3x)))))))));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (term.v * ((cse_p12_t1_g372209 * ((cse_p12_t1_g372207 * ((c1x * v3y) + (c1y * (-(v3x))) + (c2x * ((-(v1y)) + (-(v2y)))) + (c2y * (v1x + v2x)))) + (cse_p12_t1_g372208 * ((cse_p12_t1_g372212 * (-(v2z))) + (v2_len * (c2z + (v1x * v3y) + (v1y * (-(v3x))))))))) + (cse_p12_t1_g372210 * ((cse_p12_t1_g372207 * ((c1x * v3y) + (c1y * (-(v3x))) + (c2x * ((-(v1y)) + (-(v2y)))) + (c2y * (v1x + v2x)))) + (cse_p12_t1_g372208 * ((cse_p12_t1_g372212 * (-(v2z))) + (v2_len * (c2z + (v1x * v3y) + (v1y * (-(v3x)))))))))));
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double g_x3 = (term.v * ((cse_p12_t1_g372209 * ((cse_p12_t1_g372207 * ((c1y * cse_p12_t10_g372202) + (c1z * cse_p12_t12_g372204) + (c2y * v1z) + (c2z * (-(v1y))))) + (cse_p12_t1_g372208 * ((cse_p12_t1_g372212 * v2x) + (v2_len * ((cse_p12_t10_g372202 * v1y) + (cse_p12_t12_g372204 * v1z))))))) + (cse_p12_t1_g372210 * ((cse_p12_t1_g372207 * ((c1y * cse_p12_t10_g372202) + (c1z * cse_p12_t12_g372204) + (c2y * v1z) + (c2z * (-(v1y))))) + (cse_p12_t1_g372208 * ((cse_p12_t1_g372212 * v2x) + (v2_len * ((cse_p12_t10_g372202 * v1y) + (cse_p12_t12_g372204 * v1z)))))))));
      KernelGradientAcc(term.i3x3, 0, g_x3);
      double g_y3 = (term.v * ((cse_p12_t1_g372209 * ((cse_p12_t1_g372207 * ((c1x * cse_p12_t13_g372205) + (c1z * cse_p12_t8_g372200) + (c2x * (-(v1z))) + (c2z * v1x))) + (cse_p12_t1_g372208 * ((cse_p12_t1_g372212 * v2y) + (v2_len * ((cse_p12_t13_g372205 * v1x) + (cse_p12_t8_g372200 * v1z))))))) + (cse_p12_t1_g372210 * ((cse_p12_t1_g372207 * ((c1x * cse_p12_t13_g372205) + (c1z * cse_p12_t8_g372200) + (c2x * (-(v1z))) + (c2z * v1x))) + (cse_p12_t1_g372208 * ((cse_p12_t1_g372212 * v2y) + (v2_len * ((cse_p12_t13_g372205 * v1x) + (cse_p12_t8_g372200 * v1z)))))))));
      KernelGradientAcc(term.i3x3, 1, g_y3);
      double g_z3 = (term.v * ((cse_p12_t1_g372209 * ((cse_p12_t1_g372207 * ((c1x * cse_p12_t9_g372201) + (c1y * cse_p12_t11_g372203) + (c2x * v1y) + (c2y * (-(v1x))))) + (cse_p12_t1_g372208 * ((cse_p12_t1_g372212 * v2z) + (v2_len * ((cse_p12_t11_g372203 * v1y) + (cse_p12_t9_g372201 * v1x))))))) + (cse_p12_t1_g372210 * ((cse_p12_t1_g372207 * ((c1x * cse_p12_t9_g372201) + (c1y * cse_p12_t11_g372203) + (c2x * v1y) + (c2y * (-(v1x))))) + (cse_p12_t1_g372208 * ((cse_p12_t1_g372212 * v2z) + (v2_len * ((cse_p12_t11_g372203 * v1y) + (cse_p12_t9_g372201 * v1x)))))))));
      KernelGradientAcc(term.i3x3, 2, g_z3);
      double g_x4 = (term.v * ((cse_p12_t1_g372209 * ((cse_p12_t1_g372207 * ((c1y * v2z) + (c1z * (-(v2y))))) + (cse_p12_t1_g372211 * (cse_p12_t5_g372197 + (v1z * (-(v2y))))))) + (cse_p12_t1_g372210 * ((cse_p12_t1_g372207 * ((c1y * v2z) + (c1z * (-(v2y))))) + (cse_p12_t1_g372211 * (cse_p12_t5_g372197 + (v1z * (-(v2y)))))))));
      KernelGradientAcc(term.i3x4, 0, g_x4);
      double g_y4 = (term.v * ((cse_p12_t1_g372209 * ((cse_p12_t1_g372207 * ((c1x * (-(v2z))) + (c1z * v2x))) + (cse_p12_t1_g372211 * (cse_p12_t6_g372198 + (v1x * (-(v2z))))))) + (cse_p12_t1_g372210 * ((cse_p12_t1_g372207 * ((c1x * (-(v2z))) + (c1z * v2x))) + (cse_p12_t1_g372211 * (cse_p12_t6_g372198 + (v1x * (-(v2z)))))))));
      KernelGradientAcc(term.i3x4, 1, g_y4);
      double g_z4 = (term.v * ((cse_p12_t1_g372209 * ((cse_p12_t1_g372207 * ((c1x * v2y) + (c1y * (-(v2x))))) + (cse_p12_t1_g372211 * (cse_p12_t4_g372196 + (v1y * (-(v2x))))))) + (cse_p12_t1_g372210 * ((cse_p12_t1_g372207 * ((c1x * v2y) + (c1y * (-(v2x))))) + (cse_p12_t1_g372211 * (cse_p12_t4_g372196 + (v1y * (-(v2x)))))))));
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
      double cse_p12_t118_g372338 = (v1y * (-(v2x)));
      double cse_p12_t126_g372346 = (v1z * v2x);
      double cse_p12_t186_g372406 = ((-(v1x)) + (-(v2x)));
      double cse_p12_t195_g372415 = (v1x + v2x);
      double v2y = (y3 + (-(y2)));
      double cse_p12_t114_g372334 = (v1x * v2y);
      double cse_p12_t124_g372344 = (v1z * (-(v2y)));
      double cse_p12_t183_g372403 = (cse_p12_t114_g372334 + cse_p12_t118_g372338);
      double cse_p12_t187_g372407 = ((-(v1y)) + (-(v2y)));
      double cse_p12_t196_g372416 = (v1y + v2y);
      double v2z = (z3 + (-(z2)));
      double cse_p12_t112_g372332 = (v1x * (-(v2z)));
      double cse_p12_t120_g372340 = (v1y * v2z);
      double cse_p12_t182_g372402 = (cse_p12_t112_g372332 + cse_p12_t126_g372346);
      double cse_p12_t185_g372405 = (cse_p12_t120_g372340 + cse_p12_t124_g372344);
      double cse_p12_t188_g372408 = ((-(v1z)) + (-(v2z)));
      double cse_p12_t197_g372417 = (v1z + v2z);
      double v3x = (x4 + (-(x3)));
      double cse_p12_t117_g372337 = (v1y * (v2x + v3x));
      double cse_p12_t119_g372339 = (v1y * (-(v3x)));
      double cse_p12_t122_g372342 = (v1z * ((-(v2x)) + (-(v3x))));
      double cse_p12_t127_g372347 = (v1z * v3x);
      double cse_p12_t189_g372409 = ((-(v2x)) + (-(v3x)));
      double cse_p12_t198_g372418 = (v2x + v3x);
      double v3y = (y4 + (-(y3)));
      double cse_p12_t110_g372330 = (v1x * ((-(v2y)) + (-(v3y))));
      double cse_p12_t115_g372335 = (v1x * v3y);
      double cse_p12_t123_g372343 = (v1z * (v2y + v3y));
      double cse_p12_t125_g372345 = (v1z * (-(v3y)));
      double cse_p12_t180_g372400 = (cse_p12_t110_g372330 + cse_p12_t117_g372337);
      double cse_p12_t190_g372410 = ((-(v2y)) + (-(v3y)));
      double cse_p12_t199_g372419 = (v2y + v3y);
      double v3z = (z4 + (-(z3)));
      double cse_p12_t111_g372331 = (v1x * (v2z + v3z));
      double cse_p12_t113_g372333 = (v1x * (-(v3z)));
      double cse_p12_t116_g372336 = (v1y * ((-(v2z)) + (-(v3z))));
      double cse_p12_t121_g372341 = (v1y * v3z);
      double cse_p12_t181_g372401 = (cse_p12_t111_g372331 + cse_p12_t122_g372342);
      double cse_p12_t184_g372404 = (cse_p12_t116_g372336 + cse_p12_t123_g372343);
      double cse_p12_t191_g372411 = ((-(v2z)) + (-(v3z)));
      double cse_p12_t200_g372420 = (v2z + v3z);
      double c1x = (cse_p12_t120_g372340 + (-((v1z * v2y))));
      double cse_p12_t23_g372243 = (c1x * cse_p12_t190_g372410);
      double cse_p12_t24_g372244 = (c1x * cse_p12_t200_g372420);
      double cse_p12_t25_g372245 = (c1x * (-(v2z)));
      double cse_p12_t26_g372246 = (c1x * (-(v3z)));
      double cse_p12_t27_g372247 = (c1x * v2y);
      double cse_p12_t28_g372248 = (c1x * v3y);
      double c1y = (cse_p12_t126_g372346 + (-((v1x * v2z))));
      double cse_p12_t29_g372249 = (c1y * cse_p12_t191_g372411);
      double cse_p12_t30_g372250 = (c1y * cse_p12_t198_g372418);
      double cse_p12_t31_g372251 = (c1y * (-(v2x)));
      double cse_p12_t32_g372252 = (c1y * (-(v3x)));
      double cse_p12_t33_g372253 = (c1y * v2z);
      double cse_p12_t34_g372254 = (c1y * v3z);
      double cse_p12_t171_g372391 = (cse_p12_t27_g372247 + cse_p12_t31_g372251);
      double c1z = (cse_p12_t114_g372334 + (-((v1y * v2x))));
      double cse_p12_t35_g372255 = (c1z * cse_p12_t189_g372409);
      double cse_p12_t36_g372256 = (c1z * cse_p12_t199_g372419);
      double cse_p12_t37_g372257 = (c1z * (-(v2y)));
      double cse_p12_t38_g372258 = (c1z * (-(v3y)));
      double cse_p12_t39_g372259 = (c1z * v2x);
      double cse_p12_t40_g372260 = (c1z * v3x);
      double cse_p12_t169_g372389 = (cse_p12_t25_g372245 + cse_p12_t39_g372259);
      double cse_p12_t174_g372394 = (cse_p12_t33_g372253 + cse_p12_t37_g372257);
      double c2x = ((v2y * v3z) + (-((v2z * v3y))));
      double cse_p12_t41_g372261 = (c2x * cse_p12_t187_g372407);
      double cse_p12_t42_g372262 = (c2x * cse_p12_t197_g372417);
      double cse_p12_t43_g372263 = (c2x * (-(v1z)));
      double cse_p12_t44_g372264 = (c2x * (-(v2z)));
      double cse_p12_t45_g372265 = (c2x * v1x);
      double cse_p12_t46_g372266 = (c2x * v1y);
      double cse_p12_t47_g372267 = (c2x * v2y);
      double cse_p12_t192_g372412 = (c2x + cse_p12_t121_g372341 + cse_p12_t125_g372345);
      double c2y = ((v2z * v3x) + (-((v2x * v3z))));
      double cse_p12_t48_g372268 = (c2y * cse_p12_t188_g372408);
      double cse_p12_t49_g372269 = (c2y * cse_p12_t195_g372415);
      double cse_p12_t50_g372270 = (c2y * (-(v1x)));
      double cse_p12_t51_g372271 = (c2y * (-(v2x)));
      double cse_p12_t52_g372272 = (c2y * v1y);
      double cse_p12_t53_g372273 = (c2y * v1z);
      double cse_p12_t54_g372274 = (c2y * v2z);
      double cse_p12_t167_g372387 = (cse_p12_t23_g372243 + cse_p12_t30_g372250 + cse_p12_t46_g372266 + cse_p12_t50_g372270);
      double cse_p12_t172_g372392 = (cse_p12_t28_g372248 + cse_p12_t32_g372252 + cse_p12_t41_g372261 + cse_p12_t49_g372269);
      double cse_p12_t178_g372398 = (cse_p12_t47_g372267 + cse_p12_t51_g372271);
      double cse_p12_t193_g372413 = (c2y + cse_p12_t113_g372333 + cse_p12_t127_g372347);
      double c2z = ((v2x * v3y) + (-((v2y * v3x))));
      double cse_p12_t55_g372275 = (c2z * cse_p12_t186_g372406);
      double cse_p12_t56_g372276 = (c2z * cse_p12_t196_g372416);
      double cse_p12_t57_g372277 = (c2z * (-(v1y)));
      double cse_p12_t58_g372278 = (c2z * (-(v2y)));
      double cse_p12_t59_g372279 = (c2z * v1x);
      double cse_p12_t60_g372280 = (c2z * v1z);
      double cse_p12_t61_g372281 = (c2z * v2x);
      double cse_p12_t168_g372388 = (cse_p12_t24_g372244 + cse_p12_t35_g372255 + cse_p12_t43_g372263 + cse_p12_t59_g372279);
      double cse_p12_t170_g372390 = (cse_p12_t26_g372246 + cse_p12_t40_g372260 + cse_p12_t42_g372262 + cse_p12_t55_g372275);
      double cse_p12_t173_g372393 = (cse_p12_t29_g372249 + cse_p12_t36_g372256 + cse_p12_t53_g372273 + cse_p12_t57_g372277);
      double cse_p12_t175_g372395 = (cse_p12_t34_g372254 + cse_p12_t38_g372258 + cse_p12_t48_g372268 + cse_p12_t56_g372276);
      double cse_p12_t176_g372396 = (cse_p12_t44_g372264 + cse_p12_t61_g372281);
      double cse_p12_t177_g372397 = (cse_p12_t45_g372265 + cse_p12_t52_g372272 + cse_p12_t60_g372280);
      double cse_p12_t179_g372399 = (cse_p12_t54_g372274 + cse_p12_t58_g372278);
      double cse_p12_t194_g372414 = (c2z + cse_p12_t115_g372335 + cse_p12_t119_g372339);
      double v2_sq = ((v2x * v2x) + (v2y * v2y) + (v2z * v2z));
      double cse_p11_t1_sqrt372216 = sqrt(v2_sq);
      double cse_p11_t2_invsqrt372217 = (1.0 / (cse_p11_t1_sqrt372216));
      double cse_p12_t1_g372444 = (cse_p11_t2_invsqrt372217 * cse_p12_t177_g372397);
      double cse_p12_t62_g372282 = (cse_p12_t1_g372444 * (-(v2x)));
      double cse_p12_t63_g372283 = (cse_p12_t1_g372444 * (-(v2y)));
      double cse_p12_t64_g372284 = (cse_p12_t1_g372444 * (-(v2z)));
      double cse_p12_t72_g372292 = (cse_p12_t1_g372444 * v2x);
      double cse_p12_t73_g372293 = (cse_p12_t1_g372444 * v2y);
      double cse_p12_t74_g372294 = (cse_p12_t1_g372444 * v2z);
      double cse_p12_t201_g372421 = (-(cse_p12_t1_g372444));
      double cse_p11_t3_invr372218 = (cse_p11_t2_invsqrt372217 * cse_p11_t2_invsqrt372217);
      double cse_p12_t66_g372286 = (cse_p11_t3_invr372218 * cse_p12_t1_g372444 * v2x * v2x);
      double cse_p12_t67_g372287 = (cse_p11_t3_invr372218 * cse_p12_t1_g372444 * v2x * v2y);
      double cse_p12_t68_g372288 = (cse_p11_t3_invr372218 * cse_p12_t1_g372444 * v2x * v2z);
      double cse_p12_t69_g372289 = (cse_p11_t3_invr372218 * cse_p12_t1_g372444 * v2y * v2y);
      double cse_p12_t70_g372290 = (cse_p11_t3_invr372218 * cse_p12_t1_g372444 * v2y * v2z);
      double cse_p12_t71_g372291 = (cse_p11_t3_invr372218 * cse_p12_t1_g372444 * v2z * v2z);
      double v2_len = cse_p11_t1_sqrt372216;
      double cse_p12_t1_g372451 = (cse_p12_t180_g372400 * v2_len);
      double cse_p12_t1_g372452 = (cse_p12_t181_g372401 * v2_len);
      double cse_p12_t130_g372350 = (cse_p12_t182_g372402 * v2_len);
      double cse_p12_t131_g372351 = (cse_p12_t183_g372403 * v2_len);
      double cse_p12_t1_g372453 = (cse_p12_t184_g372404 * v2_len);
      double cse_p12_t133_g372353 = (cse_p12_t185_g372405 * v2_len);
      double cse_p12_t134_g372354 = (cse_p12_t192_g372412 * v2_len);
      double cse_p12_t135_g372355 = (cse_p12_t193_g372413 * v2_len);
      double cse_p12_t136_g372356 = (cse_p12_t194_g372414 * v2_len);
      double cse_p12_t137_g372357 = (v2_len * (-(c2z)));
      double cse_p12_t147_g372367 = (cse_p12_t134_g372354 + cse_p12_t62_g372282);
      double cse_p12_t148_g372368 = (cse_p12_t135_g372355 + cse_p12_t63_g372283);
      double cse_p12_t149_g372369 = (cse_p12_t136_g372356 + cse_p12_t64_g372284);
      double cse_p12_t150_g372370 = (cse_p12_t1_g372453 + cse_p12_t72_g372292);
      double cse_p12_t151_g372371 = (cse_p12_t1_g372452 + cse_p12_t73_g372293);
      double cse_p12_t152_g372372 = (cse_p12_t1_g372451 + cse_p12_t74_g372294);
      double cse_p12_t213_g372433 = (-(cse_p12_t130_g372350));
      double cse_p12_t214_g372434 = (-(cse_p12_t131_g372351));
      double cse_p12_t215_g372435 = (-(cse_p12_t133_g372353));
      double cse_p12_t216_g372436 = (-(cse_p12_t147_g372367));
      double cse_p12_t217_g372437 = (-(cse_p12_t148_g372368));
      double cse_p12_t218_g372438 = (-(cse_p12_t149_g372369));
      double cse_p12_t219_g372439 = (-(cse_p12_t150_g372370));
      double cse_p12_t220_g372440 = (-(cse_p12_t151_g372371));
      double cse_p12_t221_g372441 = (-(cse_p12_t152_g372372));
      double cse_p12_t222_g372442 = (-(cse_p12_t137_g372357));
      double t1 = (cse_p12_t177_g372397 * v2_len);
      double cse_p13_t1_g372460 = (2.00000000000000000e+0 * t1);
      double cse_p12_t1_g372221 = (cse_p12_t147_g372367 * cse_p13_t1_g372460);
      double cse_p12_t2_g372222 = (cse_p12_t148_g372368 * cse_p13_t1_g372460);
      double cse_p12_t3_g372223 = (cse_p12_t149_g372369 * cse_p13_t1_g372460);
      double cse_p12_t4_g372224 = (cse_p12_t150_g372370 * cse_p13_t1_g372460);
      double cse_p12_t5_g372225 = (cse_p12_t151_g372371 * cse_p13_t1_g372460);
      double cse_p12_t6_g372226 = (cse_p12_t152_g372372 * cse_p13_t1_g372460);
      double cse_p12_t7_g372227 = (cse_p12_t182_g372402 * cse_p13_t1_g372460 * v2_len);
      double cse_p12_t8_g372228 = (cse_p12_t183_g372403 * cse_p13_t1_g372460 * v2_len);
      double cse_p12_t9_g372229 = (cse_p12_t185_g372405 * cse_p13_t1_g372460 * v2_len);
      double cse_p12_t10_g372230 = (cse_p13_t1_g372460 * v2_len * (-(c2y)));
      double cse_p12_t11_g372231 = (cse_p13_t1_g372460 * v2_len * (-(c2z)));
      double t2 = ((c1x * c2x) + (c1y * c2y) + (c1z * c2z));
      double cse_p13_t1_g372461 = (2.00000000000000000e+0 * t2);
      double cse_p12_t12_g372232 = (cse_p12_t167_g372387 * cse_p13_t1_g372461);
      double cse_p12_t13_g372233 = (cse_p12_t168_g372388 * cse_p13_t1_g372461);
      double cse_p12_t14_g372234 = (cse_p12_t169_g372389 * cse_p13_t1_g372461);
      double cse_p12_t15_g372235 = (cse_p12_t170_g372390 * cse_p13_t1_g372461);
      double cse_p12_t16_g372236 = (cse_p12_t171_g372391 * cse_p13_t1_g372461);
      double cse_p12_t17_g372237 = (cse_p12_t172_g372392 * cse_p13_t1_g372461);
      double cse_p12_t18_g372238 = (cse_p12_t173_g372393 * cse_p13_t1_g372461);
      double cse_p12_t19_g372239 = (cse_p12_t174_g372394 * cse_p13_t1_g372461);
      double cse_p12_t20_g372240 = (cse_p12_t175_g372395 * cse_p13_t1_g372461);
      double cse_p12_t21_g372241 = (cse_p12_t176_g372396 * cse_p13_t1_g372461);
      double cse_p12_t22_g372242 = (cse_p12_t178_g372398 * cse_p13_t1_g372461);
      double cse_p12_t138_g372358 = (cse_p12_t1_g372221 + cse_p12_t20_g372240);
      double cse_p12_t139_g372359 = (cse_p12_t15_g372235 + cse_p12_t2_g372222);
      double cse_p12_t140_g372360 = (cse_p12_t17_g372237 + cse_p12_t3_g372223);
      double cse_p12_t141_g372361 = (cse_p12_t18_g372238 + cse_p12_t4_g372224);
      double cse_p12_t142_g372362 = (cse_p12_t13_g372233 + cse_p12_t5_g372225);
      double cse_p12_t143_g372363 = (cse_p12_t12_g372232 + cse_p12_t6_g372226);
      double cse_p12_t144_g372364 = (cse_p12_t14_g372234 + cse_p12_t7_g372227);
      double cse_p12_t145_g372365 = (cse_p12_t16_g372236 + cse_p12_t8_g372228);
      double cse_p12_t146_g372366 = (cse_p12_t19_g372239 + cse_p12_t9_g372229);
      double cse_p12_t165_g372385 = (cse_p12_t10_g372230 + cse_p12_t21_g372241);
      double cse_p12_t166_g372386 = (cse_p12_t11_g372231 + cse_p12_t22_g372242);
      double phi = atan2(t1, t2);
      double nphi = (term.n * phi);
      double cse_p12_t223_g372443 = (-(sin(nphi)));
      double sin_nphi = sin(nphi);
      double cos_nphi = cos(nphi);
      double cos_angle = ((term.cosphase * cos_nphi) + (term.sinphase * sin_nphi));
      double sin_angle = ((term.cosphase * sin_nphi) + (-((cos_nphi * term.sinphase))));
      double energy = (term.v * (1.00000000000000000e+0 + cos_angle));
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t4_invr372219 = (1.0 / ((((t1) * (t1)) + ((t2) * (t2)))));
      double cse_p13_t1_g372458 = (cse_p11_t4_invr372219 * (-(t1)));
      double cse_p12_t75_g372295 = (cse_p12_t167_g372387 * cse_p13_t1_g372458);
      double cse_p12_t76_g372296 = (cse_p12_t168_g372388 * cse_p13_t1_g372458);
      double cse_p12_t77_g372297 = (cse_p12_t169_g372389 * cse_p13_t1_g372458);
      double cse_p12_t78_g372298 = (cse_p12_t170_g372390 * cse_p13_t1_g372458);
      double cse_p12_t79_g372299 = (cse_p12_t171_g372391 * cse_p13_t1_g372458);
      double cse_p12_t80_g372300 = (cse_p12_t172_g372392 * cse_p13_t1_g372458);
      double cse_p12_t81_g372301 = (cse_p12_t173_g372393 * cse_p13_t1_g372458);
      double cse_p12_t82_g372302 = (cse_p12_t174_g372394 * cse_p13_t1_g372458);
      double cse_p12_t83_g372303 = (cse_p12_t175_g372395 * cse_p13_t1_g372458);
      double cse_p12_t84_g372304 = (cse_p12_t176_g372396 * cse_p13_t1_g372458);
      double cse_p12_t85_g372305 = (cse_p12_t178_g372398 * cse_p13_t1_g372458);
      double cse_p12_t86_g372306 = (cse_p12_t179_g372399 * cse_p13_t1_g372458);
      double cse_p12_t1_g372447 = (cse_p11_t4_invr372219 * t2);
      double cse_p12_t87_g372307 = (cse_p12_t147_g372367 * cse_p12_t1_g372447);
      double cse_p12_t88_g372308 = (cse_p12_t148_g372368 * cse_p12_t1_g372447);
      double cse_p12_t89_g372309 = (cse_p12_t149_g372369 * cse_p12_t1_g372447);
      double cse_p12_t90_g372310 = (cse_p12_t150_g372370 * cse_p12_t1_g372447);
      double cse_p12_t91_g372311 = (cse_p12_t151_g372371 * cse_p12_t1_g372447);
      double cse_p12_t92_g372312 = (cse_p12_t152_g372372 * cse_p12_t1_g372447);
      double cse_p13_t1_g372459 = (cse_p12_t1_g372447 * v2_len);
      double cse_p12_t93_g372313 = (cse_p12_t182_g372402 * cse_p13_t1_g372459);
      double cse_p12_t94_g372314 = (cse_p12_t183_g372403 * cse_p13_t1_g372459);
      double cse_p12_t95_g372315 = (cse_p12_t185_g372405 * cse_p13_t1_g372459);
      double cse_p12_t96_g372316 = (cse_p13_t1_g372459 * (-(c2x)));
      double cse_p12_t97_g372317 = (cse_p13_t1_g372459 * (-(c2y)));
      double cse_p12_t98_g372318 = (cse_p13_t1_g372459 * (-(c2z)));
      double cse_p12_t153_g372373 = (cse_p12_t75_g372295 + cse_p12_t92_g372312);
      double cse_p12_t154_g372374 = (cse_p12_t76_g372296 + cse_p12_t91_g372311);
      double cse_p12_t155_g372375 = (cse_p12_t77_g372297 + cse_p12_t93_g372313);
      double cse_p12_t156_g372376 = (cse_p12_t78_g372298 + cse_p12_t88_g372308);
      double cse_p12_t157_g372377 = (cse_p12_t79_g372299 + cse_p12_t94_g372314);
      double cse_p12_t158_g372378 = (cse_p12_t80_g372300 + cse_p12_t89_g372309);
      double cse_p12_t159_g372379 = (cse_p12_t81_g372301 + cse_p12_t90_g372310);
      double cse_p12_t160_g372380 = (cse_p12_t82_g372302 + cse_p12_t95_g372315);
      double cse_p12_t161_g372381 = (cse_p12_t83_g372303 + cse_p12_t87_g372307);
      double cse_p12_t162_g372382 = (cse_p12_t84_g372304 + cse_p12_t97_g372317);
      double cse_p12_t163_g372383 = (cse_p12_t85_g372305 + cse_p12_t98_g372318);
      double cse_p12_t164_g372384 = (cse_p12_t86_g372306 + cse_p12_t96_g372316);
      double cse_p11_t5_invr2372220 = (cse_p11_t4_invr372219 * cse_p11_t4_invr372219);
      double cse_p12_t99_g372319 = (cse_p11_t5_invr2372220 * cse_p12_t138_g372358);
      double cse_p12_t100_g372320 = (cse_p11_t5_invr2372220 * cse_p12_t139_g372359);
      double cse_p12_t101_g372321 = (cse_p11_t5_invr2372220 * cse_p12_t140_g372360);
      double cse_p12_t102_g372322 = (cse_p11_t5_invr2372220 * cse_p12_t141_g372361);
      double cse_p12_t103_g372323 = (cse_p11_t5_invr2372220 * cse_p12_t142_g372362);
      double cse_p12_t104_g372324 = (cse_p11_t5_invr2372220 * cse_p12_t143_g372363);
      double cse_p12_t105_g372325 = (cse_p11_t5_invr2372220 * cse_p12_t144_g372364);
      double cse_p12_t106_g372326 = (cse_p11_t5_invr2372220 * cse_p12_t145_g372365);
      double cse_p12_t107_g372327 = (cse_p11_t5_invr2372220 * cse_p12_t146_g372366);
      double cse_p12_t108_g372328 = (cse_p11_t5_invr2372220 * cse_p12_t165_g372385);
      double cse_p12_t109_g372329 = (cse_p11_t5_invr2372220 * cse_p12_t166_g372386);
      double cse_p12_t202_g372422 = (-(cse_p12_t99_g372319));
      double cse_p12_t203_g372423 = (-(cse_p12_t100_g372320));
      double cse_p12_t204_g372424 = (-(cse_p12_t101_g372321));
      double cse_p12_t205_g372425 = (-(cse_p12_t102_g372322));
      double cse_p12_t206_g372426 = (-(cse_p12_t103_g372323));
      double cse_p12_t207_g372427 = (-(cse_p12_t104_g372324));
      double cse_p12_t208_g372428 = (-(cse_p12_t105_g372325));
      double cse_p12_t209_g372429 = (-(cse_p12_t106_g372326));
      double cse_p12_t210_g372430 = (-(cse_p12_t107_g372327));
      double cse_p12_t211_g372431 = (-(cse_p12_t108_g372328));
      double cse_p12_t212_g372432 = (-(cse_p12_t109_g372329));
      double cse_p12_t1_g372450 = (t2 * v2_len);
      double cse_p13_t1_g372463 = (cse_p12_t1_g372450 * (-(c2x)));
      double cse_p13_t1_g372464 = (cos_nphi * term.n * term.sinphase);
      double cse_p14_t1_g372471 = (term.cosphase * cse_p12_t223_g372443 * term.n);
      double g_x1 = (term.v * ((cse_p13_t1_g372464 * (cse_p12_t86_g372306 + (cse_p11_t4_invr372219 * cse_p13_t1_g372463))) + (cse_p14_t1_g372471 * (cse_p12_t86_g372306 + (cse_p11_t4_invr372219 * cse_p13_t1_g372463)))));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double cse_p13_t1_g372465 = (cse_p12_t1_g372450 * (-(c2y)));
      double g_y1 = (term.v * ((cse_p13_t1_g372464 * (cse_p12_t84_g372304 + (cse_p11_t4_invr372219 * cse_p13_t1_g372465))) + (cse_p14_t1_g372471 * (cse_p12_t84_g372304 + (cse_p11_t4_invr372219 * cse_p13_t1_g372465)))));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double cse_p14_t1_g372467 = (cse_p12_t1_g372450 * (-(c2z)));
      double g_z1 = (term.v * ((cse_p13_t1_g372464 * (cse_p12_t85_g372305 + (cse_p11_t4_invr372219 * cse_p14_t1_g372467))) + (cse_p14_t1_g372471 * (cse_p12_t85_g372305 + (cse_p11_t4_invr372219 * cse_p14_t1_g372467)))));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double g_x2 = (term.v * ((cse_p13_t1_g372464 * (cse_p12_t83_g372303 + (cse_p11_t4_invr372219 * t2 * (cse_p12_t134_g372354 + (cse_p11_t2_invsqrt372217 * cse_p12_t177_g372397 * (-(v2x))))))) + (cse_p14_t1_g372471 * (cse_p12_t83_g372303 + (cse_p11_t4_invr372219 * t2 * (cse_p12_t134_g372354 + (cse_p11_t2_invsqrt372217 * cse_p12_t177_g372397 * (-(v2x)))))))));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (term.v * ((cse_p13_t1_g372464 * (cse_p12_t78_g372298 + (cse_p11_t4_invr372219 * t2 * (cse_p12_t135_g372355 + (cse_p11_t2_invsqrt372217 * cse_p12_t177_g372397 * (-(v2y))))))) + (cse_p14_t1_g372471 * (cse_p12_t78_g372298 + (cse_p11_t4_invr372219 * t2 * (cse_p12_t135_g372355 + (cse_p11_t2_invsqrt372217 * cse_p12_t177_g372397 * (-(v2y)))))))));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (term.v * ((cse_p13_t1_g372464 * (cse_p12_t80_g372300 + (cse_p11_t4_invr372219 * t2 * (cse_p12_t136_g372356 + (cse_p11_t2_invsqrt372217 * cse_p12_t177_g372397 * (-(v2z))))))) + (cse_p14_t1_g372471 * (cse_p12_t80_g372300 + (cse_p11_t4_invr372219 * t2 * (cse_p12_t136_g372356 + (cse_p11_t2_invsqrt372217 * cse_p12_t177_g372397 * (-(v2z)))))))));
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double g_x3 = (term.v * ((cse_p13_t1_g372464 * (cse_p12_t81_g372301 + (cse_p11_t4_invr372219 * t2 * (cse_p12_t1_g372453 + (cse_p11_t2_invsqrt372217 * cse_p12_t177_g372397 * v2x))))) + (cse_p14_t1_g372471 * (cse_p12_t81_g372301 + (cse_p11_t4_invr372219 * t2 * (cse_p12_t1_g372453 + (cse_p11_t2_invsqrt372217 * cse_p12_t177_g372397 * v2x)))))));
      KernelGradientAcc(term.i3x3, 0, g_x3);
      double g_y3 = (term.v * ((cse_p13_t1_g372464 * (cse_p12_t76_g372296 + (cse_p11_t4_invr372219 * t2 * (cse_p12_t1_g372452 + (cse_p11_t2_invsqrt372217 * cse_p12_t177_g372397 * v2y))))) + (cse_p14_t1_g372471 * (cse_p12_t76_g372296 + (cse_p11_t4_invr372219 * t2 * (cse_p12_t1_g372452 + (cse_p11_t2_invsqrt372217 * cse_p12_t177_g372397 * v2y)))))));
      KernelGradientAcc(term.i3x3, 1, g_y3);
      double g_z3 = (term.v * ((cse_p13_t1_g372464 * (cse_p12_t75_g372295 + (cse_p11_t4_invr372219 * t2 * (cse_p12_t1_g372451 + (cse_p11_t2_invsqrt372217 * cse_p12_t177_g372397 * v2z))))) + (cse_p14_t1_g372471 * (cse_p12_t75_g372295 + (cse_p11_t4_invr372219 * t2 * (cse_p12_t1_g372451 + (cse_p11_t2_invsqrt372217 * cse_p12_t177_g372397 * v2z)))))));
      KernelGradientAcc(term.i3x3, 2, g_z3);
      double cse_p12_t1_g372456 = (cse_p11_t4_invr372219 * cse_p12_t185_g372405 * cse_p12_t1_g372450);
      double g_x4 = (term.v * ((cse_p13_t1_g372464 * (cse_p12_t1_g372456 + cse_p12_t82_g372302)) + (cse_p14_t1_g372471 * (cse_p12_t1_g372456 + cse_p12_t82_g372302))));
      KernelGradientAcc(term.i3x4, 0, g_x4);
      double cse_p12_t1_g372454 = (cse_p11_t4_invr372219 * cse_p12_t182_g372402 * cse_p12_t1_g372450);
      double g_y4 = (term.v * ((cse_p13_t1_g372464 * (cse_p12_t1_g372454 + cse_p12_t77_g372297)) + (cse_p14_t1_g372471 * (cse_p12_t1_g372454 + cse_p12_t77_g372297))));
      KernelGradientAcc(term.i3x4, 1, g_y4);
      double cse_p12_t1_g372455 = (cse_p11_t4_invr372219 * cse_p12_t183_g372403 * cse_p12_t1_g372450);
      double g_z4 = (term.v * ((cse_p13_t1_g372464 * (cse_p12_t1_g372455 + cse_p12_t79_g372299)) + (cse_p14_t1_g372471 * (cse_p12_t1_g372455 + cse_p12_t79_g372299))));
      KernelGradientAcc(term.i3x4, 2, g_z4);
      double cse_p12_t1_g372445 = (term.n * (-(term.v)));
      double cse_p12_t1_g372448 = (cos_angle * cse_p12_t1_g372445 * term.n);
      double cse_p12_t1_g372449 = (cse_p12_t1_g372445 * sin_angle);
      double cse_p13_t1_g372457 = (cse_p11_t4_invr372219 * v2_len);
      double cse_p13_t1_g372466 = (cse_p13_t1_g372457 * (-(c2x)));
      double h_x1_x1 = ((cse_p12_t164_g372384 * cse_p12_t164_g372384 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p12_t179_g372399 * (cse_p13_t1_g372466 + (-((cse_p11_t4_invr372219 * v2_len * (-(c2x))))) + (-((cse_p11_t5_invr2372220 * ((cse_p12_t179_g372399 * cse_p13_t1_g372461) + (cse_p13_t1_g372460 * v2_len * (-(c2x)))) * (-(t1))))))) + (-((cse_p11_t5_invr2372220 * cse_p13_t1_g372463 * ((cse_p12_t179_g372399 * cse_p13_t1_g372461) + (cse_p13_t1_g372460 * v2_len * (-(c2x))))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
      double h_x1_y1 = ((cse_p12_t162_g372382 * cse_p12_t164_g372384 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p12_t176_g372396 * cse_p13_t1_g372466) + (cse_p12_t179_g372399 * cse_p12_t211_g372431 * (-(t1))) + (cse_p12_t211_g372431 * cse_p13_t1_g372463) + (-((cse_p11_t4_invr372219 * cse_p12_t179_g372399 * v2_len * (-(c2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
      double h_x1_z1 = ((cse_p12_t163_g372383 * cse_p12_t164_g372384 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t179_g372399 * cse_p12_t222_g372442) + (cse_p12_t178_g372398 * cse_p13_t1_g372466) + (cse_p12_t179_g372399 * cse_p12_t212_g372432 * (-(t1))) + (cse_p12_t212_g372432 * cse_p13_t1_g372463))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
      double cse_p13_t1_g372462 = (cse_p11_t2_invsqrt372217 * cse_p12_t1_g372447);
      double h_x1_x2 = ((cse_p12_t161_g372381 * cse_p12_t164_g372384 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t179_g372399 * cse_p12_t216_g372436) + (cse_p12_t175_g372395 * cse_p13_t1_g372466) + (cse_p12_t179_g372399 * cse_p12_t202_g372422 * (-(t1))) + (cse_p12_t202_g372422 * cse_p13_t1_g372463) + (cse_p13_t1_g372458 * (((-(v2y)) * (-(v3y))) + (v2z * v3z))) + (-((cse_p13_t1_g372462 * v2x * (-(c2x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
      double h_x1_y2 = ((cse_p12_t156_g372376 * cse_p12_t164_g372384 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t179_g372399 * cse_p12_t217_g372437) + (cse_p12_t170_g372390 * cse_p13_t1_g372466) + (cse_p12_t179_g372399 * cse_p12_t203_g372423 * (-(t1))) + (cse_p12_t203_g372423 * cse_p13_t1_g372463) + (cse_p13_t1_g372458 * (c2z + (v3x * (-(v2y))))) + (cse_p13_t1_g372459 * v3z) + (-((cse_p13_t1_g372462 * v2y * (-(c2x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
      double h_x1_z2 = ((cse_p12_t158_g372378 * cse_p12_t164_g372384 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t179_g372399 * cse_p12_t218_g372438) + (cse_p12_t172_g372392 * cse_p13_t1_g372466) + (cse_p12_t179_g372399 * cse_p12_t204_g372424 * (-(t1))) + (cse_p12_t204_g372424 * cse_p13_t1_g372463) + (cse_p13_t1_g372458 * ((v2z * (-(v3x))) + (-(c2y)))) + (cse_p13_t1_g372459 * (-(v3y))) + (-((cse_p13_t1_g372462 * v2z * (-(c2x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
      double h_x1_x3 = ((cse_p12_t159_g372379 * cse_p12_t164_g372384 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t179_g372399 * cse_p12_t219_g372439) + (cse_p12_t173_g372393 * cse_p13_t1_g372466) + (cse_p12_t179_g372399 * cse_p12_t205_g372425 * (-(t1))) + (cse_p12_t205_g372425 * cse_p13_t1_g372463) + (cse_p13_t1_g372458 * ((cse_p12_t191_g372411 * v2z) + (cse_p12_t199_g372419 * (-(v2y))))) + (cse_p13_t1_g372462 * v2x * (-(c2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 0, h_x1_x3);
      double h_x1_y3 = ((cse_p12_t154_g372374 * cse_p12_t164_g372384 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t179_g372399 * cse_p12_t220_g372440) + (cse_p12_t168_g372388 * cse_p13_t1_g372466) + (cse_p12_t179_g372399 * cse_p12_t206_g372426 * (-(t1))) + (cse_p12_t191_g372411 * cse_p13_t1_g372459) + (cse_p12_t206_g372426 * cse_p13_t1_g372463) + (cse_p13_t1_g372458 * ((cse_p12_t189_g372409 * (-(v2y))) + (-(c2z)))) + (cse_p13_t1_g372462 * v2y * (-(c2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 1, h_x1_y3);
      double cse_p14_t1_g372474 = (cse_p12_t207_g372427 * (-(t1)));
      double h_x1_z3 = ((cse_p12_t153_g372373 * cse_p12_t164_g372384 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t179_g372399 * cse_p12_t221_g372441) + (cse_p12_t167_g372387 * cse_p13_t1_g372466) + (cse_p12_t179_g372399 * cse_p14_t1_g372474) + (cse_p12_t199_g372419 * cse_p13_t1_g372459) + (cse_p12_t207_g372427 * cse_p13_t1_g372463) + (cse_p13_t1_g372458 * (c2y + (cse_p12_t198_g372418 * v2z))) + (cse_p13_t1_g372462 * v2z * (-(c2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 2, h_x1_z3);
      double cse_p14_t1_g372473 = (cse_p12_t210_g372430 * (-(t1)));
      double h_x1_x4 = ((cse_p12_t160_g372380 * cse_p12_t164_g372384 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t179_g372399 * cse_p12_t215_g372435) + (cse_p12_t174_g372394 * cse_p13_t1_g372466) + (cse_p12_t179_g372399 * cse_p14_t1_g372473) + (cse_p12_t210_g372430 * cse_p13_t1_g372463) + (cse_p13_t1_g372458 * (((-(v2y)) * (-(v2y))) + (v2z * v2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x4, 0, h_x1_x4);
      double cse_p14_t1_g372470 = (cse_p12_t208_g372428 * (-(t1)));
      double h_x1_y4 = ((cse_p12_t155_g372375 * cse_p12_t164_g372384 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t179_g372399 * cse_p12_t213_g372433) + (cse_p12_t169_g372389 * cse_p13_t1_g372466) + (cse_p12_t179_g372399 * cse_p14_t1_g372470) + (cse_p12_t208_g372428 * cse_p13_t1_g372463) + (cse_p13_t1_g372458 * v2x * (-(v2y))) + (cse_p13_t1_g372459 * v2z))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x4, 1, h_x1_y4);
      double cse_p14_t1_g372468 = (cse_p12_t209_g372429 * (-(t1)));
      double h_x1_z4 = ((cse_p12_t157_g372377 * cse_p12_t164_g372384 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t179_g372399 * cse_p12_t214_g372434) + (cse_p12_t171_g372391 * cse_p13_t1_g372466) + (cse_p12_t179_g372399 * cse_p14_t1_g372468) + (cse_p12_t209_g372429 * cse_p13_t1_g372463) + (cse_p13_t1_g372458 * v2z * (-(v2x))) + (cse_p13_t1_g372459 * (-(v2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x4, 2, h_x1_z4);
      double cse_p14_t1_g372469 = (cse_p13_t1_g372457 * (-(c2y)));
      double h_y1_y1 = ((cse_p12_t162_g372382 * cse_p12_t162_g372382 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p12_t176_g372396 * (cse_p14_t1_g372469 + (cse_p12_t211_g372431 * (-(t1))) + (-((cse_p11_t4_invr372219 * v2_len * (-(c2y))))))) + (cse_p12_t211_g372431 * cse_p13_t1_g372465))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
      double h_y1_z1 = ((cse_p12_t162_g372382 * cse_p12_t163_g372383 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t176_g372396 * cse_p12_t222_g372442) + (cse_p12_t176_g372396 * cse_p12_t212_g372432 * (-(t1))) + (cse_p12_t178_g372398 * cse_p14_t1_g372469) + (cse_p12_t212_g372432 * cse_p13_t1_g372465))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
      double h_y1_x2 = ((cse_p12_t161_g372381 * cse_p12_t162_g372382 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t176_g372396 * cse_p12_t216_g372436) + (cse_p12_t175_g372395 * cse_p14_t1_g372469) + (cse_p12_t176_g372396 * cse_p12_t202_g372422 * (-(t1))) + (cse_p12_t202_g372422 * cse_p13_t1_g372465) + (cse_p13_t1_g372458 * ((v2x * (-(v3y))) + (-(c2z)))) + (cse_p13_t1_g372459 * (-(v3z))) + (-((cse_p13_t1_g372462 * v2x * (-(c2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
      double h_y1_y2 = ((cse_p12_t156_g372376 * cse_p12_t162_g372382 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t176_g372396 * cse_p12_t217_g372437) + (cse_p12_t170_g372390 * cse_p14_t1_g372469) + (cse_p12_t176_g372396 * cse_p12_t203_g372423 * (-(t1))) + (cse_p12_t203_g372423 * cse_p13_t1_g372465) + (cse_p13_t1_g372458 * (((-(v2z)) * (-(v3z))) + (v2x * v3x))) + (-((cse_p13_t1_g372462 * v2y * (-(c2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
      double h_y1_z2 = ((cse_p12_t158_g372378 * cse_p12_t162_g372382 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t176_g372396 * cse_p12_t218_g372438) + (cse_p12_t172_g372392 * cse_p14_t1_g372469) + (cse_p12_t176_g372396 * cse_p12_t204_g372424 * (-(t1))) + (cse_p12_t204_g372424 * cse_p13_t1_g372465) + (cse_p13_t1_g372458 * (c2x + (v3y * (-(v2z))))) + (cse_p13_t1_g372459 * v3x) + (-((cse_p13_t1_g372462 * v2z * (-(c2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
      double h_y1_x3 = ((cse_p12_t159_g372379 * cse_p12_t162_g372382 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t176_g372396 * cse_p12_t219_g372439) + (cse_p12_t173_g372393 * cse_p14_t1_g372469) + (cse_p12_t176_g372396 * cse_p12_t205_g372425 * (-(t1))) + (cse_p12_t200_g372420 * cse_p13_t1_g372459) + (cse_p12_t205_g372425 * cse_p13_t1_g372465) + (cse_p13_t1_g372458 * (c2z + (cse_p12_t199_g372419 * v2x))) + (cse_p13_t1_g372462 * v2x * (-(c2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 0, h_y1_x3);
      double h_y1_y3 = ((cse_p12_t154_g372374 * cse_p12_t162_g372382 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t176_g372396 * cse_p12_t220_g372440) + (cse_p12_t168_g372388 * cse_p14_t1_g372469) + (cse_p12_t176_g372396 * cse_p12_t206_g372426 * (-(t1))) + (cse_p12_t206_g372426 * cse_p13_t1_g372465) + (cse_p13_t1_g372458 * ((cse_p12_t189_g372409 * v2x) + (cse_p12_t200_g372420 * (-(v2z))))) + (cse_p13_t1_g372462 * v2y * (-(c2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 1, h_y1_y3);
      double h_y1_z3 = ((cse_p12_t153_g372373 * cse_p12_t162_g372382 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t176_g372396 * cse_p12_t221_g372441) + (cse_p12_t167_g372387 * cse_p14_t1_g372469) + (cse_p12_t176_g372396 * cse_p14_t1_g372474) + (cse_p12_t189_g372409 * cse_p13_t1_g372459) + (cse_p12_t207_g372427 * cse_p13_t1_g372465) + (cse_p13_t1_g372458 * ((cse_p12_t190_g372410 * (-(v2z))) + (-(c2x)))) + (cse_p13_t1_g372462 * v2z * (-(c2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 2, h_y1_z3);
      double h_y1_x4 = ((cse_p12_t160_g372380 * cse_p12_t162_g372382 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t176_g372396 * cse_p12_t215_g372435) + (cse_p12_t174_g372394 * cse_p14_t1_g372469) + (cse_p12_t176_g372396 * cse_p14_t1_g372473) + (cse_p12_t210_g372430 * cse_p13_t1_g372465) + (cse_p13_t1_g372458 * v2x * (-(v2y))) + (cse_p13_t1_g372459 * (-(v2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x4, 0, h_y1_x4);
      double h_y1_y4 = ((cse_p12_t155_g372375 * cse_p12_t162_g372382 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t176_g372396 * cse_p12_t213_g372433) + (cse_p12_t169_g372389 * cse_p14_t1_g372469) + (cse_p12_t176_g372396 * cse_p14_t1_g372470) + (cse_p12_t208_g372428 * cse_p13_t1_g372465) + (cse_p13_t1_g372458 * (((-(v2z)) * (-(v2z))) + (v2x * v2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x4, 1, h_y1_y4);
      double h_y1_z4 = ((cse_p12_t157_g372377 * cse_p12_t162_g372382 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t176_g372396 * cse_p12_t214_g372434) + (cse_p12_t171_g372391 * cse_p14_t1_g372469) + (cse_p12_t176_g372396 * cse_p14_t1_g372468) + (cse_p12_t209_g372429 * cse_p13_t1_g372465) + (cse_p13_t1_g372458 * v2y * (-(v2z))) + (cse_p13_t1_g372459 * v2x))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x4, 2, h_y1_z4);
      double cse_p14_t1_g372472 = (cse_p13_t1_g372457 * (-(c2z)));
      double h_z1_z1 = ((cse_p12_t163_g372383 * cse_p12_t163_g372383 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p12_t178_g372398 * (cse_p14_t1_g372472 + (cse_p11_t4_invr372219 * cse_p12_t222_g372442) + (cse_p12_t212_g372432 * (-(t1))))) + (cse_p12_t212_g372432 * cse_p14_t1_g372467))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
      double h_z1_x2 = ((cse_p12_t161_g372381 * cse_p12_t163_g372383 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t178_g372398 * cse_p12_t216_g372436) + (cse_p12_t175_g372395 * cse_p14_t1_g372472) + (cse_p12_t178_g372398 * cse_p12_t202_g372422 * (-(t1))) + (cse_p12_t202_g372422 * cse_p14_t1_g372467) + (cse_p13_t1_g372458 * (c2y + (v3z * (-(v2x))))) + (cse_p13_t1_g372459 * v3y) + (-((cse_p13_t1_g372462 * v2x * (-(c2z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
      double h_z1_y2 = ((cse_p12_t156_g372376 * cse_p12_t163_g372383 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t178_g372398 * cse_p12_t217_g372437) + (cse_p12_t170_g372390 * cse_p14_t1_g372472) + (cse_p12_t178_g372398 * cse_p12_t203_g372423 * (-(t1))) + (cse_p12_t203_g372423 * cse_p14_t1_g372467) + (cse_p13_t1_g372458 * ((v2y * (-(v3z))) + (-(c2x)))) + (cse_p13_t1_g372459 * (-(v3x))) + (-((cse_p13_t1_g372462 * v2y * (-(c2z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
      double h_z1_z2 = ((cse_p12_t158_g372378 * cse_p12_t163_g372383 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t178_g372398 * cse_p12_t218_g372438) + (cse_p12_t172_g372392 * cse_p14_t1_g372472) + (cse_p12_t178_g372398 * cse_p12_t204_g372424 * (-(t1))) + (cse_p12_t204_g372424 * cse_p14_t1_g372467) + (cse_p13_t1_g372458 * (((-(v2x)) * (-(v3x))) + (v2y * v3y))) + (-((cse_p13_t1_g372462 * v2z * (-(c2z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
      double h_z1_x3 = ((cse_p12_t159_g372379 * cse_p12_t163_g372383 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t178_g372398 * cse_p12_t219_g372439) + (cse_p12_t173_g372393 * cse_p14_t1_g372472) + (cse_p12_t178_g372398 * cse_p12_t205_g372425 * (-(t1))) + (cse_p12_t190_g372410 * cse_p13_t1_g372459) + (cse_p12_t205_g372425 * cse_p14_t1_g372467) + (cse_p13_t1_g372458 * ((cse_p12_t191_g372411 * (-(v2x))) + (-(c2y)))) + (cse_p13_t1_g372462 * v2x * (-(c2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 0, h_z1_x3);
      double h_z1_y3 = ((cse_p12_t154_g372374 * cse_p12_t163_g372383 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t178_g372398 * cse_p12_t220_g372440) + (cse_p12_t168_g372388 * cse_p14_t1_g372472) + (cse_p12_t178_g372398 * cse_p12_t206_g372426 * (-(t1))) + (cse_p12_t198_g372418 * cse_p13_t1_g372459) + (cse_p12_t206_g372426 * cse_p14_t1_g372467) + (cse_p13_t1_g372458 * (c2x + (cse_p12_t200_g372420 * v2y))) + (cse_p13_t1_g372462 * v2y * (-(c2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 1, h_z1_y3);
      double h_z1_z3 = ((cse_p12_t153_g372373 * cse_p12_t163_g372383 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t178_g372398 * cse_p12_t221_g372441) + (cse_p12_t167_g372387 * cse_p14_t1_g372472) + (cse_p12_t178_g372398 * cse_p14_t1_g372474) + (cse_p12_t207_g372427 * cse_p14_t1_g372467) + (cse_p13_t1_g372458 * ((cse_p12_t190_g372410 * v2y) + (cse_p12_t198_g372418 * (-(v2x))))) + (cse_p13_t1_g372462 * v2z * (-(c2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 2, h_z1_z3);
      double h_z1_x4 = ((cse_p12_t160_g372380 * cse_p12_t163_g372383 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t178_g372398 * cse_p12_t215_g372435) + (cse_p12_t174_g372394 * cse_p14_t1_g372472) + (cse_p12_t178_g372398 * cse_p14_t1_g372473) + (cse_p12_t210_g372430 * cse_p14_t1_g372467) + (cse_p13_t1_g372458 * v2z * (-(v2x))) + (cse_p13_t1_g372459 * v2y))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x4, 0, h_z1_x4);
      double h_z1_y4 = ((cse_p12_t155_g372375 * cse_p12_t163_g372383 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t178_g372398 * cse_p12_t213_g372433) + (cse_p12_t169_g372389 * cse_p14_t1_g372472) + (cse_p12_t178_g372398 * cse_p14_t1_g372470) + (cse_p12_t208_g372428 * cse_p14_t1_g372467) + (cse_p13_t1_g372458 * v2y * (-(v2z))) + (cse_p13_t1_g372459 * (-(v2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x4, 1, h_z1_y4);
      double h_z1_z4 = ((cse_p12_t157_g372377 * cse_p12_t163_g372383 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t178_g372398 * cse_p12_t214_g372434) + (cse_p12_t171_g372391 * cse_p14_t1_g372472) + (cse_p12_t178_g372398 * cse_p14_t1_g372468) + (cse_p12_t209_g372429 * cse_p14_t1_g372467) + (cse_p13_t1_g372458 * (((-(v2x)) * (-(v2x))) + (v2y * v2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x4, 2, h_z1_z4);
      double cse_p14_t1_g372475 = (cse_p11_t2_invsqrt372217 * (-(v2x)));
      double h_x2_x2 = ((cse_p12_t161_g372381 * cse_p12_t161_g372381 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p12_t147_g372367 * cse_p12_t202_g372422 * t2) + (cse_p12_t175_g372395 * ((cse_p11_t4_invr372219 * cse_p12_t147_g372367) + (cse_p11_t4_invr372219 * cse_p12_t216_g372436) + (cse_p12_t202_g372422 * (-(t1))))) + (cse_p12_t1_g372447 * (cse_p12_t1_g372444 + (cse_p12_t192_g372412 * cse_p14_t1_g372475) + (cse_p12_t192_g372412 * cse_p14_t1_g372475) + (-(cse_p12_t66_g372286)))) + (cse_p13_t1_g372458 * ((cse_p12_t188_g372408 * v3z) + (cse_p12_t188_g372408 * v3z) + (cse_p12_t196_g372416 * (-(v3y))) + (cse_p12_t196_g372416 * (-(v3y))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
      double cse_p14_t1_g372476 = (cse_p11_t2_invsqrt372217 * (-(v2y)));
      double h_x2_y2 = ((cse_p12_t156_g372376 * cse_p12_t161_g372381 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t147_g372367 * cse_p12_t170_g372390) + (cse_p11_t4_invr372219 * cse_p12_t175_g372395 * cse_p12_t217_g372437) + (cse_p12_t147_g372367 * cse_p12_t203_g372423 * t2) + (cse_p12_t175_g372395 * cse_p12_t203_g372423 * (-(t1))) + (cse_p12_t1_g372447 * ((cse_p12_t192_g372412 * cse_p14_t1_g372476) + (cse_p12_t193_g372413 * cse_p14_t1_g372475) + (-(cse_p12_t67_g372287)))) + (cse_p13_t1_g372458 * ((cse_p12_t186_g372406 * (-(v3y))) + (cse_p12_t196_g372416 * v3x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
      double cse_p15_t1_g372477 = (cse_p11_t2_invsqrt372217 * (-(v2z)));
      double h_x2_z2 = ((cse_p12_t158_g372378 * cse_p12_t161_g372381 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t147_g372367 * cse_p12_t172_g372392) + (cse_p11_t4_invr372219 * cse_p12_t175_g372395 * cse_p12_t218_g372438) + (cse_p12_t147_g372367 * cse_p12_t204_g372424 * t2) + (cse_p12_t175_g372395 * cse_p12_t204_g372424 * (-(t1))) + (cse_p12_t1_g372447 * ((cse_p12_t192_g372412 * cse_p15_t1_g372477) + (cse_p12_t194_g372414 * cse_p14_t1_g372475) + (-(cse_p12_t68_g372288)))) + (cse_p13_t1_g372458 * ((cse_p12_t188_g372408 * (-(v3x))) + (cse_p12_t195_g372415 * v3z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
      double h_x2_x3 = ((cse_p12_t159_g372379 * cse_p12_t161_g372381 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t147_g372367 * cse_p12_t173_g372393) + (cse_p11_t4_invr372219 * cse_p12_t175_g372395 * cse_p12_t219_g372439) + (cse_p12_t147_g372367 * cse_p12_t205_g372425 * t2) + (cse_p12_t175_g372395 * cse_p12_t205_g372425 * (-(t1))) + (cse_p12_t1_g372447 * (cse_p12_t201_g372421 + cse_p12_t66_g372286 + (cse_p11_t2_invsqrt372217 * cse_p12_t192_g372412 * v2x) + (cse_p12_t184_g372404 * cse_p14_t1_g372475))) + (cse_p13_t1_g372458 * (((-(v1y)) * (-(v3y))) + (cse_p12_t188_g372408 * cse_p12_t191_g372411) + (cse_p12_t196_g372416 * cse_p12_t199_g372419) + (v1z * v3z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 0, h_x2_x3);
      double h_x2_y3 = ((cse_p12_t154_g372374 * cse_p12_t161_g372381 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t147_g372367 * cse_p12_t168_g372388) + (cse_p11_t4_invr372219 * cse_p12_t175_g372395 * cse_p12_t220_g372440) + (cse_p12_t147_g372367 * cse_p12_t206_g372426 * t2) + (cse_p12_t175_g372395 * cse_p12_t206_g372426 * (-(t1))) + (cse_p12_t1_g372447 * (cse_p12_t67_g372287 + (cse_p11_t2_invsqrt372217 * cse_p12_t192_g372412 * v2y) + (cse_p12_t181_g372401 * cse_p14_t1_g372475) + (v2_len * (v1z + v2z + v3z)))) + (cse_p13_t1_g372458 * (c1z + c2z + (cse_p12_t189_g372409 * cse_p12_t196_g372416) + (v1x * (-(v3y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 1, h_x2_y3);
      double h_x2_z3 = ((cse_p12_t153_g372373 * cse_p12_t161_g372381 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t147_g372367 * cse_p12_t167_g372387) + (cse_p11_t4_invr372219 * cse_p12_t175_g372395 * cse_p12_t221_g372441) + (cse_p12_t147_g372367 * cse_p12_t207_g372427 * t2) + (cse_p12_t175_g372395 * cse_p14_t1_g372474) + (cse_p12_t1_g372447 * (cse_p12_t68_g372288 + (cse_p11_t2_invsqrt372217 * cse_p12_t192_g372412 * v2z) + (cse_p12_t180_g372400 * cse_p14_t1_g372475) + (v2_len * ((-(v1y)) + (-(v2y)) + (-(v3y)))))) + (cse_p13_t1_g372458 * ((cse_p12_t188_g372408 * cse_p12_t198_g372418) + (v3z * (-(v1x))) + (-(c1y)) + (-(c2y)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 2, h_x2_z3);
      double h_x2_x4 = ((cse_p12_t160_g372380 * cse_p12_t161_g372381 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t147_g372367 * cse_p12_t174_g372394) + (cse_p11_t4_invr372219 * cse_p12_t175_g372395 * cse_p12_t215_g372435) + (cse_p12_t147_g372367 * cse_p12_t210_g372430 * t2) + (cse_p12_t175_g372395 * cse_p14_t1_g372473) + (cse_p13_t1_g372458 * ((cse_p12_t188_g372408 * v2z) + (cse_p12_t196_g372416 * (-(v2y))))) + (-((cse_p12_t185_g372405 * cse_p13_t1_g372462 * v2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x4, 0, h_x2_x4);
      double h_x2_y4 = ((cse_p12_t155_g372375 * cse_p12_t161_g372381 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t147_g372367 * cse_p12_t169_g372389) + (cse_p11_t4_invr372219 * cse_p12_t175_g372395 * cse_p12_t213_g372433) + (cse_p12_t147_g372367 * cse_p12_t208_g372428 * t2) + (cse_p12_t175_g372395 * cse_p14_t1_g372470) + (cse_p12_t1_g372447 * ((cse_p12_t182_g372402 * cse_p14_t1_g372475) + (cse_p12_t188_g372408 * v2_len))) + (cse_p13_t1_g372458 * ((cse_p12_t196_g372416 * v2x) + (-(c1z)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x4, 1, h_x2_y4);
      double h_x2_z4 = ((cse_p12_t157_g372377 * cse_p12_t161_g372381 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t147_g372367 * cse_p12_t171_g372391) + (cse_p11_t4_invr372219 * cse_p12_t175_g372395 * cse_p12_t214_g372434) + (cse_p12_t147_g372367 * cse_p12_t209_g372429 * t2) + (cse_p12_t175_g372395 * cse_p14_t1_g372468) + (cse_p12_t1_g372447 * ((cse_p12_t183_g372403 * cse_p14_t1_g372475) + (cse_p12_t196_g372416 * v2_len))) + (cse_p13_t1_g372458 * (c1y + (cse_p12_t188_g372408 * (-(v2x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x4, 2, h_x2_z4);
      double h_y2_y2 = ((cse_p12_t156_g372376 * cse_p12_t156_g372376 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p12_t148_g372368 * cse_p12_t203_g372423 * t2) + (cse_p12_t170_g372390 * ((cse_p11_t4_invr372219 * cse_p12_t148_g372368) + (cse_p11_t4_invr372219 * cse_p12_t217_g372437) + (cse_p12_t203_g372423 * (-(t1))))) + (cse_p12_t1_g372447 * (cse_p12_t1_g372444 + (cse_p12_t193_g372413 * cse_p14_t1_g372476) + (cse_p12_t193_g372413 * cse_p14_t1_g372476) + (-(cse_p12_t69_g372289)))) + (cse_p13_t1_g372458 * ((cse_p12_t186_g372406 * v3x) + (cse_p12_t186_g372406 * v3x) + (cse_p12_t197_g372417 * (-(v3z))) + (cse_p12_t197_g372417 * (-(v3z))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
      double h_y2_z2 = ((cse_p12_t156_g372376 * cse_p12_t158_g372378 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t148_g372368 * cse_p12_t172_g372392) + (cse_p11_t4_invr372219 * cse_p12_t170_g372390 * cse_p12_t218_g372438) + (cse_p12_t148_g372368 * cse_p12_t204_g372424 * t2) + (cse_p12_t170_g372390 * cse_p12_t204_g372424 * (-(t1))) + (cse_p12_t1_g372447 * ((cse_p12_t193_g372413 * cse_p15_t1_g372477) + (cse_p12_t194_g372414 * cse_p14_t1_g372476) + (-(cse_p12_t70_g372290)))) + (cse_p13_t1_g372458 * ((cse_p12_t187_g372407 * (-(v3z))) + (cse_p12_t197_g372417 * v3y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
      double h_y2_x3 = ((cse_p12_t156_g372376 * cse_p12_t159_g372379 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t148_g372368 * cse_p12_t173_g372393) + (cse_p11_t4_invr372219 * cse_p12_t170_g372390 * cse_p12_t219_g372439) + (cse_p12_t148_g372368 * cse_p12_t205_g372425 * t2) + (cse_p12_t170_g372390 * cse_p12_t205_g372425 * (-(t1))) + (cse_p12_t1_g372447 * (cse_p12_t67_g372287 + (cse_p11_t2_invsqrt372217 * cse_p12_t193_g372413 * v2x) + (cse_p12_t184_g372404 * cse_p14_t1_g372476) + (v2_len * ((-(v1z)) + (-(v2z)) + (-(v3z)))))) + (cse_p13_t1_g372458 * ((cse_p12_t186_g372406 * cse_p12_t199_g372419) + (v3x * (-(v1y))) + (-(c1z)) + (-(c2z)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 0, h_y2_x3);
      double h_y2_y3 = ((cse_p12_t154_g372374 * cse_p12_t156_g372376 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t148_g372368 * cse_p12_t168_g372388) + (cse_p11_t4_invr372219 * cse_p12_t170_g372390 * cse_p12_t220_g372440) + (cse_p12_t148_g372368 * cse_p12_t206_g372426 * t2) + (cse_p12_t170_g372390 * cse_p12_t206_g372426 * (-(t1))) + (cse_p12_t1_g372447 * (cse_p12_t201_g372421 + cse_p12_t69_g372289 + (cse_p11_t2_invsqrt372217 * cse_p12_t193_g372413 * v2y) + (cse_p12_t181_g372401 * cse_p14_t1_g372476))) + (cse_p13_t1_g372458 * (((-(v1z)) * (-(v3z))) + (cse_p12_t186_g372406 * cse_p12_t189_g372409) + (cse_p12_t197_g372417 * cse_p12_t200_g372420) + (v1x * v3x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 1, h_y2_y3);
      double h_y2_z3 = ((cse_p12_t153_g372373 * cse_p12_t156_g372376 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t148_g372368 * cse_p12_t167_g372387) + (cse_p11_t4_invr372219 * cse_p12_t170_g372390 * cse_p12_t221_g372441) + (cse_p12_t148_g372368 * cse_p12_t207_g372427 * t2) + (cse_p12_t170_g372390 * cse_p14_t1_g372474) + (cse_p12_t1_g372447 * (cse_p12_t70_g372290 + (cse_p11_t2_invsqrt372217 * cse_p12_t193_g372413 * v2z) + (cse_p12_t180_g372400 * cse_p14_t1_g372476) + (v2_len * (v1x + v2x + v3x)))) + (cse_p13_t1_g372458 * (c1x + c2x + (cse_p12_t190_g372410 * cse_p12_t197_g372417) + (v1y * (-(v3z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 2, h_y2_z3);
      double h_y2_x4 = ((cse_p12_t156_g372376 * cse_p12_t160_g372380 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t148_g372368 * cse_p12_t174_g372394) + (cse_p11_t4_invr372219 * cse_p12_t170_g372390 * cse_p12_t215_g372435) + (cse_p12_t148_g372368 * cse_p12_t210_g372430 * t2) + (cse_p12_t170_g372390 * cse_p14_t1_g372473) + (cse_p12_t1_g372447 * ((cse_p12_t185_g372405 * cse_p14_t1_g372476) + (cse_p12_t197_g372417 * v2_len))) + (cse_p13_t1_g372458 * (c1z + (cse_p12_t186_g372406 * (-(v2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x4, 0, h_y2_x4);
      double h_y2_y4 = ((cse_p12_t155_g372375 * cse_p12_t156_g372376 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t148_g372368 * cse_p12_t169_g372389) + (cse_p11_t4_invr372219 * cse_p12_t170_g372390 * cse_p12_t213_g372433) + (cse_p12_t148_g372368 * cse_p12_t208_g372428 * t2) + (cse_p12_t170_g372390 * cse_p14_t1_g372470) + (cse_p13_t1_g372458 * ((cse_p12_t186_g372406 * v2x) + (cse_p12_t197_g372417 * (-(v2z))))) + (-((cse_p12_t182_g372402 * cse_p13_t1_g372462 * v2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x4, 1, h_y2_y4);
      double h_y2_z4 = ((cse_p12_t156_g372376 * cse_p12_t157_g372377 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t148_g372368 * cse_p12_t171_g372391) + (cse_p11_t4_invr372219 * cse_p12_t170_g372390 * cse_p12_t214_g372434) + (cse_p12_t148_g372368 * cse_p12_t209_g372429 * t2) + (cse_p12_t170_g372390 * cse_p14_t1_g372468) + (cse_p12_t1_g372447 * ((cse_p12_t183_g372403 * cse_p14_t1_g372476) + (cse_p12_t186_g372406 * v2_len))) + (cse_p13_t1_g372458 * ((cse_p12_t197_g372417 * v2y) + (-(c1x)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x4, 2, h_y2_z4);
      double h_z2_z2 = ((cse_p12_t158_g372378 * cse_p12_t158_g372378 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p12_t149_g372369 * cse_p12_t204_g372424 * t2) + (cse_p12_t172_g372392 * ((cse_p11_t4_invr372219 * cse_p12_t149_g372369) + (cse_p11_t4_invr372219 * cse_p12_t218_g372438) + (cse_p12_t204_g372424 * (-(t1))))) + (cse_p12_t1_g372447 * (cse_p12_t1_g372444 + (cse_p12_t194_g372414 * cse_p15_t1_g372477) + (cse_p12_t194_g372414 * cse_p15_t1_g372477) + (-(cse_p12_t71_g372291)))) + (cse_p13_t1_g372458 * ((cse_p12_t187_g372407 * v3y) + (cse_p12_t187_g372407 * v3y) + (cse_p12_t195_g372415 * (-(v3x))) + (cse_p12_t195_g372415 * (-(v3x))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      double h_z2_x3 = ((cse_p12_t158_g372378 * cse_p12_t159_g372379 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t149_g372369 * cse_p12_t173_g372393) + (cse_p11_t4_invr372219 * cse_p12_t172_g372392 * cse_p12_t219_g372439) + (cse_p12_t149_g372369 * cse_p12_t205_g372425 * t2) + (cse_p12_t172_g372392 * cse_p12_t205_g372425 * (-(t1))) + (cse_p12_t1_g372447 * (cse_p12_t68_g372288 + (cse_p11_t2_invsqrt372217 * cse_p12_t194_g372414 * v2x) + (cse_p12_t184_g372404 * cse_p15_t1_g372477) + (v2_len * (v1y + v2y + v3y)))) + (cse_p13_t1_g372458 * (c1y + c2y + (cse_p12_t191_g372411 * cse_p12_t195_g372415) + (v1z * (-(v3x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 0, h_z2_x3);
      double h_z2_y3 = ((cse_p12_t154_g372374 * cse_p12_t158_g372378 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t149_g372369 * cse_p12_t168_g372388) + (cse_p11_t4_invr372219 * cse_p12_t172_g372392 * cse_p12_t220_g372440) + (cse_p12_t149_g372369 * cse_p12_t206_g372426 * t2) + (cse_p12_t172_g372392 * cse_p12_t206_g372426 * (-(t1))) + (cse_p12_t1_g372447 * (cse_p12_t70_g372290 + (cse_p11_t2_invsqrt372217 * cse_p12_t194_g372414 * v2y) + (cse_p12_t181_g372401 * cse_p15_t1_g372477) + (v2_len * ((-(v1x)) + (-(v2x)) + (-(v3x)))))) + (cse_p13_t1_g372458 * ((cse_p12_t187_g372407 * cse_p12_t200_g372420) + (v3y * (-(v1z))) + (-(c1x)) + (-(c2x)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 1, h_z2_y3);
      double h_z2_z3 = ((cse_p12_t153_g372373 * cse_p12_t158_g372378 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t149_g372369 * cse_p12_t167_g372387) + (cse_p11_t4_invr372219 * cse_p12_t172_g372392 * cse_p12_t221_g372441) + (cse_p12_t149_g372369 * cse_p12_t207_g372427 * t2) + (cse_p12_t172_g372392 * cse_p14_t1_g372474) + (cse_p12_t1_g372447 * (cse_p12_t201_g372421 + cse_p12_t71_g372291 + (cse_p11_t2_invsqrt372217 * cse_p12_t194_g372414 * v2z) + (cse_p12_t180_g372400 * cse_p15_t1_g372477))) + (cse_p13_t1_g372458 * (((-(v1x)) * (-(v3x))) + (cse_p12_t187_g372407 * cse_p12_t190_g372410) + (cse_p12_t195_g372415 * cse_p12_t198_g372418) + (v1y * v3y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 2, h_z2_z3);
      double h_z2_x4 = ((cse_p12_t158_g372378 * cse_p12_t160_g372380 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t149_g372369 * cse_p12_t174_g372394) + (cse_p11_t4_invr372219 * cse_p12_t172_g372392 * cse_p12_t215_g372435) + (cse_p12_t149_g372369 * cse_p12_t210_g372430 * t2) + (cse_p12_t172_g372392 * cse_p14_t1_g372473) + (cse_p12_t1_g372447 * ((cse_p12_t185_g372405 * cse_p15_t1_g372477) + (cse_p12_t187_g372407 * v2_len))) + (cse_p13_t1_g372458 * ((cse_p12_t195_g372415 * v2z) + (-(c1y)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x4, 0, h_z2_x4);
      double h_z2_y4 = ((cse_p12_t155_g372375 * cse_p12_t158_g372378 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t149_g372369 * cse_p12_t169_g372389) + (cse_p11_t4_invr372219 * cse_p12_t172_g372392 * cse_p12_t213_g372433) + (cse_p12_t149_g372369 * cse_p12_t208_g372428 * t2) + (cse_p12_t172_g372392 * cse_p14_t1_g372470) + (cse_p12_t1_g372447 * ((cse_p12_t182_g372402 * cse_p15_t1_g372477) + (cse_p12_t195_g372415 * v2_len))) + (cse_p13_t1_g372458 * (c1x + (cse_p12_t187_g372407 * (-(v2z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x4, 1, h_z2_y4);
      double h_z2_z4 = ((cse_p12_t157_g372377 * cse_p12_t158_g372378 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t149_g372369 * cse_p12_t171_g372391) + (cse_p11_t4_invr372219 * cse_p12_t172_g372392 * cse_p12_t214_g372434) + (cse_p12_t149_g372369 * cse_p12_t209_g372429 * t2) + (cse_p12_t172_g372392 * cse_p14_t1_g372468) + (cse_p13_t1_g372458 * ((cse_p12_t187_g372407 * v2y) + (cse_p12_t195_g372415 * (-(v2x))))) + (-((cse_p12_t183_g372403 * cse_p13_t1_g372462 * v2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x4, 2, h_z2_z4);
      double h_x3_x3 = ((cse_p12_t159_g372379 * cse_p12_t159_g372379 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p12_t150_g372370 * cse_p12_t205_g372425 * t2) + (cse_p12_t173_g372393 * ((cse_p11_t4_invr372219 * cse_p12_t150_g372370) + (cse_p11_t4_invr372219 * cse_p12_t219_g372439) + (cse_p12_t205_g372425 * (-(t1))))) + (cse_p12_t1_g372447 * (cse_p12_t1_g372444 + (cse_p11_t2_invsqrt372217 * cse_p12_t184_g372404 * v2x) + (cse_p11_t2_invsqrt372217 * cse_p12_t184_g372404 * v2x) + (-(cse_p12_t66_g372286)))) + (cse_p13_t1_g372458 * ((cse_p12_t191_g372411 * v1z) + (cse_p12_t191_g372411 * v1z) + (cse_p12_t199_g372419 * (-(v1y))) + (cse_p12_t199_g372419 * (-(v1y))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 0, h_x3_x3);
      double h_x3_y3 = ((cse_p12_t154_g372374 * cse_p12_t159_g372379 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t150_g372370 * cse_p12_t168_g372388) + (cse_p11_t4_invr372219 * cse_p12_t173_g372393 * cse_p12_t220_g372440) + (cse_p12_t150_g372370 * cse_p12_t206_g372426 * t2) + (cse_p12_t173_g372393 * cse_p12_t206_g372426 * (-(t1))) + (cse_p12_t1_g372447 * ((cse_p11_t2_invsqrt372217 * cse_p12_t181_g372401 * v2x) + (cse_p11_t2_invsqrt372217 * cse_p12_t184_g372404 * v2y) + (-(cse_p12_t67_g372287)))) + (cse_p13_t1_g372458 * ((cse_p12_t189_g372409 * (-(v1y))) + (cse_p12_t199_g372419 * v1x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 1, h_x3_y3);
      double h_x3_z3 = ((cse_p12_t153_g372373 * cse_p12_t159_g372379 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t150_g372370 * cse_p12_t167_g372387) + (cse_p11_t4_invr372219 * cse_p12_t173_g372393 * cse_p12_t221_g372441) + (cse_p12_t150_g372370 * cse_p12_t207_g372427 * t2) + (cse_p12_t173_g372393 * cse_p14_t1_g372474) + (cse_p12_t1_g372447 * ((cse_p11_t2_invsqrt372217 * cse_p12_t180_g372400 * v2x) + (cse_p11_t2_invsqrt372217 * cse_p12_t184_g372404 * v2z) + (-(cse_p12_t68_g372288)))) + (cse_p13_t1_g372458 * ((cse_p12_t191_g372411 * (-(v1x))) + (cse_p12_t198_g372418 * v1z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 2, h_x3_z3);
      double h_x3_x4 = ((cse_p12_t159_g372379 * cse_p12_t160_g372380 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t150_g372370 * cse_p12_t174_g372394) + (cse_p11_t4_invr372219 * cse_p12_t173_g372393 * cse_p12_t215_g372435) + (cse_p12_t150_g372370 * cse_p12_t210_g372430 * t2) + (cse_p12_t173_g372393 * cse_p14_t1_g372473) + (cse_p12_t185_g372405 * cse_p13_t1_g372462 * v2x) + (cse_p13_t1_g372458 * (((-(v1y)) * (-(v2y))) + (v1z * v2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x4, 0, h_x3_x4);
      double h_x3_y4 = ((cse_p12_t155_g372375 * cse_p12_t159_g372379 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t150_g372370 * cse_p12_t169_g372389) + (cse_p11_t4_invr372219 * cse_p12_t173_g372393 * cse_p12_t213_g372433) + (cse_p12_t150_g372370 * cse_p12_t208_g372428 * t2) + (cse_p12_t173_g372393 * cse_p14_t1_g372470) + (cse_p12_t1_g372447 * ((cse_p11_t2_invsqrt372217 * cse_p12_t182_g372402 * v2x) + (v1z * v2_len))) + (cse_p13_t1_g372458 * (c1z + (v2x * (-(v1y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x4, 1, h_x3_y4);
      double h_x3_z4 = ((cse_p12_t157_g372377 * cse_p12_t159_g372379 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t150_g372370 * cse_p12_t171_g372391) + (cse_p11_t4_invr372219 * cse_p12_t173_g372393 * cse_p12_t214_g372434) + (cse_p12_t150_g372370 * cse_p12_t209_g372429 * t2) + (cse_p12_t173_g372393 * cse_p14_t1_g372468) + (cse_p12_t1_g372447 * ((cse_p11_t2_invsqrt372217 * cse_p12_t183_g372403 * v2x) + (v2_len * (-(v1y))))) + (cse_p13_t1_g372458 * ((v1z * (-(v2x))) + (-(c1y)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x4, 2, h_x3_z4);
      double h_y3_y3 = ((cse_p12_t154_g372374 * cse_p12_t154_g372374 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p12_t151_g372371 * cse_p12_t206_g372426 * t2) + (cse_p12_t168_g372388 * ((cse_p11_t4_invr372219 * cse_p12_t151_g372371) + (cse_p11_t4_invr372219 * cse_p12_t220_g372440) + (cse_p12_t206_g372426 * (-(t1))))) + (cse_p12_t1_g372447 * (cse_p12_t1_g372444 + (cse_p11_t2_invsqrt372217 * cse_p12_t181_g372401 * v2y) + (cse_p11_t2_invsqrt372217 * cse_p12_t181_g372401 * v2y) + (-(cse_p12_t69_g372289)))) + (cse_p13_t1_g372458 * ((cse_p12_t189_g372409 * v1x) + (cse_p12_t189_g372409 * v1x) + (cse_p12_t200_g372420 * (-(v1z))) + (cse_p12_t200_g372420 * (-(v1z))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 1, h_y3_y3);
      double h_y3_z3 = ((cse_p12_t153_g372373 * cse_p12_t154_g372374 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t151_g372371 * cse_p12_t167_g372387) + (cse_p11_t4_invr372219 * cse_p12_t168_g372388 * cse_p12_t221_g372441) + (cse_p12_t151_g372371 * cse_p12_t207_g372427 * t2) + (cse_p12_t168_g372388 * cse_p14_t1_g372474) + (cse_p12_t1_g372447 * ((cse_p11_t2_invsqrt372217 * cse_p12_t180_g372400 * v2y) + (cse_p11_t2_invsqrt372217 * cse_p12_t181_g372401 * v2z) + (-(cse_p12_t70_g372290)))) + (cse_p13_t1_g372458 * ((cse_p12_t190_g372410 * (-(v1z))) + (cse_p12_t200_g372420 * v1y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 2, h_y3_z3);
      double h_y3_x4 = ((cse_p12_t154_g372374 * cse_p12_t160_g372380 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t151_g372371 * cse_p12_t174_g372394) + (cse_p11_t4_invr372219 * cse_p12_t168_g372388 * cse_p12_t215_g372435) + (cse_p12_t151_g372371 * cse_p12_t210_g372430 * t2) + (cse_p12_t168_g372388 * cse_p14_t1_g372473) + (cse_p12_t1_g372447 * ((cse_p11_t2_invsqrt372217 * cse_p12_t185_g372405 * v2y) + (v2_len * (-(v1z))))) + (cse_p13_t1_g372458 * ((v1x * (-(v2y))) + (-(c1z)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x4, 0, h_y3_x4);
      double h_y3_y4 = ((cse_p12_t154_g372374 * cse_p12_t155_g372375 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t151_g372371 * cse_p12_t169_g372389) + (cse_p11_t4_invr372219 * cse_p12_t168_g372388 * cse_p12_t213_g372433) + (cse_p12_t151_g372371 * cse_p12_t208_g372428 * t2) + (cse_p12_t168_g372388 * cse_p14_t1_g372470) + (cse_p12_t182_g372402 * cse_p13_t1_g372462 * v2y) + (cse_p13_t1_g372458 * (((-(v1z)) * (-(v2z))) + (v1x * v2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x4, 1, h_y3_y4);
      double h_y3_z4 = ((cse_p12_t154_g372374 * cse_p12_t157_g372377 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t151_g372371 * cse_p12_t171_g372391) + (cse_p11_t4_invr372219 * cse_p12_t168_g372388 * cse_p12_t214_g372434) + (cse_p12_t151_g372371 * cse_p12_t209_g372429 * t2) + (cse_p12_t168_g372388 * cse_p14_t1_g372468) + (cse_p12_t1_g372447 * ((cse_p11_t2_invsqrt372217 * cse_p12_t183_g372403 * v2y) + (v1x * v2_len))) + (cse_p13_t1_g372458 * (c1x + (v2y * (-(v1z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x4, 2, h_y3_z4);
      double h_z3_z3 = ((cse_p12_t153_g372373 * cse_p12_t153_g372373 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p12_t152_g372372 * cse_p12_t207_g372427 * t2) + (cse_p12_t167_g372387 * (cse_p14_t1_g372474 + (cse_p11_t4_invr372219 * cse_p12_t152_g372372) + (cse_p11_t4_invr372219 * cse_p12_t221_g372441))) + (cse_p12_t1_g372447 * (cse_p12_t1_g372444 + (cse_p11_t2_invsqrt372217 * cse_p12_t180_g372400 * v2z) + (cse_p11_t2_invsqrt372217 * cse_p12_t180_g372400 * v2z) + (-(cse_p12_t71_g372291)))) + (cse_p13_t1_g372458 * ((cse_p12_t190_g372410 * v1y) + (cse_p12_t190_g372410 * v1y) + (cse_p12_t198_g372418 * (-(v1x))) + (cse_p12_t198_g372418 * (-(v1x))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x3, 2, h_z3_z3);
      double h_z3_x4 = ((cse_p12_t153_g372373 * cse_p12_t160_g372380 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t152_g372372 * cse_p12_t174_g372394) + (cse_p11_t4_invr372219 * cse_p12_t167_g372387 * cse_p12_t215_g372435) + (cse_p12_t152_g372372 * cse_p12_t210_g372430 * t2) + (cse_p12_t167_g372387 * cse_p14_t1_g372473) + (cse_p12_t1_g372447 * ((cse_p11_t2_invsqrt372217 * cse_p12_t185_g372405 * v2z) + (v1y * v2_len))) + (cse_p13_t1_g372458 * (c1y + (v2z * (-(v1x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x4, 0, h_z3_x4);
      double h_z3_y4 = ((cse_p12_t153_g372373 * cse_p12_t155_g372375 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t152_g372372 * cse_p12_t169_g372389) + (cse_p11_t4_invr372219 * cse_p12_t167_g372387 * cse_p12_t213_g372433) + (cse_p12_t152_g372372 * cse_p12_t208_g372428 * t2) + (cse_p12_t167_g372387 * cse_p14_t1_g372470) + (cse_p12_t1_g372447 * ((cse_p11_t2_invsqrt372217 * cse_p12_t182_g372402 * v2z) + (v2_len * (-(v1x))))) + (cse_p13_t1_g372458 * ((v1y * (-(v2z))) + (-(c1x)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x4, 1, h_z3_y4);
      double h_z3_z4 = ((cse_p12_t153_g372373 * cse_p12_t157_g372377 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t152_g372372 * cse_p12_t171_g372391) + (cse_p11_t4_invr372219 * cse_p12_t167_g372387 * cse_p12_t214_g372434) + (cse_p12_t152_g372372 * cse_p12_t209_g372429 * t2) + (cse_p12_t167_g372387 * cse_p14_t1_g372468) + (cse_p12_t183_g372403 * cse_p13_t1_g372462 * v2z) + (cse_p13_t1_g372458 * (((-(v1x)) * (-(v2x))) + (v1y * v2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x4, 2, h_z3_z4);
      double h_x4_x4 = ((cse_p12_t160_g372380 * cse_p12_t160_g372380 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p12_t174_g372394 * (cse_p14_t1_g372473 + (cse_p11_t4_invr372219 * cse_p12_t215_g372435) + (cse_p12_t185_g372405 * cse_p13_t1_g372457))) + (cse_p12_t185_g372405 * cse_p12_t1_g372450 * cse_p12_t210_g372430))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 0, term.i3x4, 0, h_x4_x4);
      double h_x4_y4 = ((cse_p12_t155_g372375 * cse_p12_t160_g372380 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t174_g372394 * cse_p12_t213_g372433) + (cse_p12_t169_g372389 * cse_p12_t185_g372405 * cse_p13_t1_g372457) + (cse_p12_t174_g372394 * cse_p14_t1_g372470) + (cse_p12_t185_g372405 * cse_p12_t1_g372450 * cse_p12_t208_g372428))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 0, term.i3x4, 1, h_x4_y4);
      double h_x4_z4 = ((cse_p12_t157_g372377 * cse_p12_t160_g372380 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t174_g372394 * cse_p12_t214_g372434) + (cse_p12_t171_g372391 * cse_p12_t185_g372405 * cse_p13_t1_g372457) + (cse_p12_t174_g372394 * cse_p14_t1_g372468) + (cse_p12_t185_g372405 * cse_p12_t1_g372450 * cse_p12_t209_g372429))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 0, term.i3x4, 2, h_x4_z4);
      double h_y4_y4 = ((cse_p12_t155_g372375 * cse_p12_t155_g372375 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p12_t169_g372389 * (cse_p14_t1_g372470 + (cse_p11_t4_invr372219 * cse_p12_t213_g372433) + (cse_p12_t182_g372402 * cse_p13_t1_g372457))) + (cse_p12_t182_g372402 * cse_p12_t1_g372450 * cse_p12_t208_g372428))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 1, term.i3x4, 1, h_y4_y4);
      double h_y4_z4 = ((cse_p12_t155_g372375 * cse_p12_t157_g372377 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p11_t4_invr372219 * cse_p12_t169_g372389 * cse_p12_t214_g372434) + (cse_p12_t169_g372389 * cse_p14_t1_g372468) + (cse_p12_t171_g372391 * cse_p12_t182_g372402 * cse_p13_t1_g372457) + (cse_p12_t182_g372402 * cse_p12_t1_g372450 * cse_p12_t209_g372429))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 1, term.i3x4, 2, h_y4_z4);
      double h_z4_z4 = ((cse_p12_t157_g372377 * cse_p12_t157_g372377 * cse_p12_t1_g372448) + (cse_p12_t1_g372449 * ((cse_p12_t171_g372391 * (cse_p14_t1_g372468 + (cse_p11_t4_invr372219 * cse_p12_t214_g372434) + (cse_p12_t183_g372403 * cse_p13_t1_g372457))) + (cse_p12_t183_g372403 * cse_p12_t1_g372450 * cse_p12_t209_g372429))));
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
