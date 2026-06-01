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
      double v2y = (y3 + (-(y2)));
      double v2z = (z3 + (-(z2)));
      double v3x = (x4 + (-(x3)));
      double v3y = (y4 + (-(y3)));
      double v3z = (z4 + (-(z3)));
      double c1x = ((v1y * v2z) + (-((v1z * v2y))));
      double c1y = ((v1z * v2x) + (-((v1x * v2z))));
      double c1z = ((v1x * v2y) + (-((v1y * v2x))));
      double c2x = ((v2y * v3z) + (-((v2z * v3y))));
      double cse_p12_t1_g4 = (c2x * v1x);
      double c2y = ((v2z * v3x) + (-((v2x * v3z))));
      double cse_p12_t2_g5 = (c2y * v1y);
      double c2z = ((v2x * v3y) + (-((v2y * v3x))));
      double cse_p12_t3_g6 = (c2z * v1z);
      double cse_p12_t4_g7 = (cse_p12_t1_g4 + cse_p12_t2_g5 + cse_p12_t3_g6);
      double v2_sq = ((v2x * v2x) + (v2y * v2y) + (v2z * v2z));
      double cse_p1_t1_sqrt1 = sqrt(v2_sq);
      double cse_p1_t2_invsqrt2 = (1.0 / (cse_p1_t1_sqrt1));
      double v2_len = cse_p1_t1_sqrt1;
      double t1 = (cse_p12_t4_g7 * v2_len);
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
      double cse_p1_t3_invr3 = (1.0 / ((((t1) * (t1)) + ((t2) * (t2)))));
      double cse_p12_t1_g8 = (term.n * sin_angle * (-(term.v)));
      double cse_p12_t1_g9 = (cse_p1_t3_invr3 * (-(t1)));
      double cse_p12_t1_g10 = (cse_p1_t3_invr3 * t2);
      double cse_p12_t1_g11 = (cse_p12_t1_g10 * v2_len);
      double g_x1 = (cse_p12_t1_g8 * ((cse_p12_t1_g11 * (-(c2x))) + (cse_p12_t1_g9 * ((c2y * v2z) + (c2z * (-(v2y)))))));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = (cse_p12_t1_g8 * ((cse_p12_t1_g11 * (-(c2y))) + (cse_p12_t1_g9 * ((c2x * (-(v2z))) + (c2z * v2x)))));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = (cse_p12_t1_g8 * ((cse_p12_t1_g11 * (-(c2z))) + (cse_p12_t1_g9 * ((c2x * v2y) + (c2y * (-(v2x)))))));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p12_t1_g12 = (cse_p12_t4_g7 * cse_p1_t2_invsqrt2);
      double g_x2 = (cse_p12_t1_g8 * ((cse_p12_t1_g10 * ((cse_p12_t1_g12 * (-(v2x))) + (v2_len * (c2x + (v1y * v3z) + (v1z * (-(v3y))))))) + (cse_p12_t1_g9 * ((c1y * v3z) + (c1z * (-(v3y))) + (c2y * ((-(v1z)) + (-(v2z)))) + (c2z * (v1y + v2y))))));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (cse_p12_t1_g8 * ((cse_p12_t1_g10 * ((cse_p12_t1_g12 * (-(v2y))) + (v2_len * (c2y + (v1x * (-(v3z))) + (v1z * v3x))))) + (cse_p12_t1_g9 * ((c1x * (-(v3z))) + (c1z * v3x) + (c2x * (v1z + v2z)) + (c2z * ((-(v1x)) + (-(v2x))))))));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (cse_p12_t1_g8 * ((cse_p12_t1_g10 * ((cse_p12_t1_g12 * (-(v2z))) + (v2_len * (c2z + (v1x * v3y) + (v1y * (-(v3x))))))) + (cse_p12_t1_g9 * ((c1x * v3y) + (c1y * (-(v3x))) + (c2x * ((-(v1y)) + (-(v2y)))) + (c2y * (v1x + v2x))))));
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double g_x3 = (cse_p12_t1_g8 * ((cse_p12_t1_g10 * ((cse_p12_t1_g12 * v2x) + (v2_len * ((v1y * ((-(v2z)) + (-(v3z)))) + (v1z * (v2y + v3y)))))) + (cse_p12_t1_g9 * ((c1y * ((-(v2z)) + (-(v3z)))) + (c1z * (v2y + v3y)) + (c2y * v1z) + (c2z * (-(v1y)))))));
      KernelGradientAcc(term.i3x3, 0, g_x3);
      double g_y3 = (cse_p12_t1_g8 * ((cse_p12_t1_g10 * ((cse_p12_t1_g12 * v2y) + (v2_len * ((v1x * (v2z + v3z)) + (v1z * ((-(v2x)) + (-(v3x)))))))) + (cse_p12_t1_g9 * ((c1x * (v2z + v3z)) + (c1z * ((-(v2x)) + (-(v3x)))) + (c2x * (-(v1z))) + (c2z * v1x)))));
      KernelGradientAcc(term.i3x3, 1, g_y3);
      double g_z3 = (cse_p12_t1_g8 * ((cse_p12_t1_g10 * ((cse_p12_t1_g12 * v2z) + (v2_len * ((v1x * ((-(v2y)) + (-(v3y)))) + (v1y * (v2x + v3x)))))) + (cse_p12_t1_g9 * ((c1x * ((-(v2y)) + (-(v3y)))) + (c1y * (v2x + v3x)) + (c2x * v1y) + (c2y * (-(v1x)))))));
      KernelGradientAcc(term.i3x3, 2, g_z3);
      double g_x4 = (cse_p12_t1_g8 * ((cse_p12_t1_g11 * ((v1y * v2z) + (v1z * (-(v2y))))) + (cse_p12_t1_g9 * ((c1y * v2z) + (c1z * (-(v2y)))))));
      KernelGradientAcc(term.i3x4, 0, g_x4);
      double g_y4 = (cse_p12_t1_g8 * ((cse_p12_t1_g11 * ((v1x * (-(v2z))) + (v1z * v2x))) + (cse_p12_t1_g9 * ((c1x * (-(v2z))) + (c1z * v2x)))));
      KernelGradientAcc(term.i3x4, 1, g_y4);
      double g_z4 = (cse_p12_t1_g8 * ((cse_p12_t1_g11 * ((v1x * v2y) + (v1y * (-(v2x))))) + (cse_p12_t1_g9 * ((c1x * v2y) + (c1y * (-(v2x)))))));
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
      double cse_p12_t118_g126 = (v1y * (-(v2x)));
      double cse_p12_t126_g134 = (v1z * v2x);
      double cse_p12_t186_g194 = ((-(v1x)) + (-(v2x)));
      double cse_p12_t195_g203 = (v1x + v2x);
      double v2y = (y3 + (-(y2)));
      double cse_p12_t114_g122 = (v1x * v2y);
      double cse_p12_t124_g132 = (v1z * (-(v2y)));
      double cse_p12_t183_g191 = (cse_p12_t114_g122 + cse_p12_t118_g126);
      double cse_p12_t187_g195 = ((-(v1y)) + (-(v2y)));
      double cse_p12_t196_g204 = (v1y + v2y);
      double v2z = (z3 + (-(z2)));
      double cse_p12_t112_g120 = (v1x * (-(v2z)));
      double cse_p12_t120_g128 = (v1y * v2z);
      double cse_p12_t182_g190 = (cse_p12_t112_g120 + cse_p12_t126_g134);
      double cse_p12_t185_g193 = (cse_p12_t120_g128 + cse_p12_t124_g132);
      double cse_p12_t188_g196 = ((-(v1z)) + (-(v2z)));
      double cse_p12_t197_g205 = (v1z + v2z);
      double v3x = (x4 + (-(x3)));
      double cse_p12_t117_g125 = (v1y * (v2x + v3x));
      double cse_p12_t119_g127 = (v1y * (-(v3x)));
      double cse_p12_t122_g130 = (v1z * ((-(v2x)) + (-(v3x))));
      double cse_p12_t127_g135 = (v1z * v3x);
      double cse_p12_t189_g197 = ((-(v2x)) + (-(v3x)));
      double cse_p12_t198_g206 = (v2x + v3x);
      double v3y = (y4 + (-(y3)));
      double cse_p12_t110_g118 = (v1x * ((-(v2y)) + (-(v3y))));
      double cse_p12_t115_g123 = (v1x * v3y);
      double cse_p12_t123_g131 = (v1z * (v2y + v3y));
      double cse_p12_t125_g133 = (v1z * (-(v3y)));
      double cse_p12_t180_g188 = (cse_p12_t110_g118 + cse_p12_t117_g125);
      double cse_p12_t190_g198 = ((-(v2y)) + (-(v3y)));
      double cse_p12_t199_g207 = (v2y + v3y);
      double v3z = (z4 + (-(z3)));
      double cse_p12_t111_g119 = (v1x * (v2z + v3z));
      double cse_p12_t113_g121 = (v1x * (-(v3z)));
      double cse_p12_t116_g124 = (v1y * ((-(v2z)) + (-(v3z))));
      double cse_p12_t121_g129 = (v1y * v3z);
      double cse_p12_t181_g189 = (cse_p12_t111_g119 + cse_p12_t122_g130);
      double cse_p12_t184_g192 = (cse_p12_t116_g124 + cse_p12_t123_g131);
      double cse_p12_t191_g199 = ((-(v2z)) + (-(v3z)));
      double cse_p12_t200_g208 = (v2z + v3z);
      double c1x = (cse_p12_t120_g128 + (-((v1z * v2y))));
      double cse_p12_t23_g31 = (c1x * cse_p12_t190_g198);
      double cse_p12_t24_g32 = (c1x * cse_p12_t200_g208);
      double cse_p12_t25_g33 = (c1x * (-(v2z)));
      double cse_p12_t26_g34 = (c1x * (-(v3z)));
      double cse_p12_t27_g35 = (c1x * v2y);
      double cse_p12_t28_g36 = (c1x * v3y);
      double c1y = (cse_p12_t126_g134 + (-((v1x * v2z))));
      double cse_p12_t29_g37 = (c1y * cse_p12_t191_g199);
      double cse_p12_t30_g38 = (c1y * cse_p12_t198_g206);
      double cse_p12_t31_g39 = (c1y * (-(v2x)));
      double cse_p12_t32_g40 = (c1y * (-(v3x)));
      double cse_p12_t33_g41 = (c1y * v2z);
      double cse_p12_t34_g42 = (c1y * v3z);
      double cse_p12_t171_g179 = (cse_p12_t27_g35 + cse_p12_t31_g39);
      double c1z = (cse_p12_t114_g122 + (-((v1y * v2x))));
      double cse_p12_t35_g43 = (c1z * cse_p12_t189_g197);
      double cse_p12_t36_g44 = (c1z * cse_p12_t199_g207);
      double cse_p12_t37_g45 = (c1z * (-(v2y)));
      double cse_p12_t38_g46 = (c1z * (-(v3y)));
      double cse_p12_t39_g47 = (c1z * v2x);
      double cse_p12_t40_g48 = (c1z * v3x);
      double cse_p12_t169_g177 = (cse_p12_t25_g33 + cse_p12_t39_g47);
      double cse_p12_t174_g182 = (cse_p12_t33_g41 + cse_p12_t37_g45);
      double c2x = ((v2y * v3z) + (-((v2z * v3y))));
      double cse_p12_t41_g49 = (c2x * cse_p12_t187_g195);
      double cse_p12_t42_g50 = (c2x * cse_p12_t197_g205);
      double cse_p12_t43_g51 = (c2x * (-(v1z)));
      double cse_p12_t44_g52 = (c2x * (-(v2z)));
      double cse_p12_t45_g53 = (c2x * v1x);
      double cse_p12_t46_g54 = (c2x * v1y);
      double cse_p12_t47_g55 = (c2x * v2y);
      double cse_p12_t192_g200 = (c2x + cse_p12_t121_g129 + cse_p12_t125_g133);
      double c2y = ((v2z * v3x) + (-((v2x * v3z))));
      double cse_p12_t48_g56 = (c2y * cse_p12_t188_g196);
      double cse_p12_t49_g57 = (c2y * cse_p12_t195_g203);
      double cse_p12_t50_g58 = (c2y * (-(v1x)));
      double cse_p12_t51_g59 = (c2y * (-(v2x)));
      double cse_p12_t52_g60 = (c2y * v1y);
      double cse_p12_t53_g61 = (c2y * v1z);
      double cse_p12_t54_g62 = (c2y * v2z);
      double cse_p12_t167_g175 = (cse_p12_t23_g31 + cse_p12_t30_g38 + cse_p12_t46_g54 + cse_p12_t50_g58);
      double cse_p12_t172_g180 = (cse_p12_t28_g36 + cse_p12_t32_g40 + cse_p12_t41_g49 + cse_p12_t49_g57);
      double cse_p12_t178_g186 = (cse_p12_t47_g55 + cse_p12_t51_g59);
      double cse_p12_t193_g201 = (c2y + cse_p12_t113_g121 + cse_p12_t127_g135);
      double c2z = ((v2x * v3y) + (-((v2y * v3x))));
      double cse_p12_t55_g63 = (c2z * cse_p12_t186_g194);
      double cse_p12_t56_g64 = (c2z * cse_p12_t196_g204);
      double cse_p12_t57_g65 = (c2z * (-(v1y)));
      double cse_p12_t58_g66 = (c2z * (-(v2y)));
      double cse_p12_t59_g67 = (c2z * v1x);
      double cse_p12_t60_g68 = (c2z * v1z);
      double cse_p12_t61_g69 = (c2z * v2x);
      double cse_p12_t168_g176 = (cse_p12_t24_g32 + cse_p12_t35_g43 + cse_p12_t43_g51 + cse_p12_t59_g67);
      double cse_p12_t170_g178 = (cse_p12_t26_g34 + cse_p12_t40_g48 + cse_p12_t42_g50 + cse_p12_t55_g63);
      double cse_p12_t173_g181 = (cse_p12_t29_g37 + cse_p12_t36_g44 + cse_p12_t53_g61 + cse_p12_t57_g65);
      double cse_p12_t175_g183 = (cse_p12_t34_g42 + cse_p12_t38_g46 + cse_p12_t48_g56 + cse_p12_t56_g64);
      double cse_p12_t176_g184 = (cse_p12_t44_g52 + cse_p12_t61_g69);
      double cse_p12_t177_g185 = (cse_p12_t45_g53 + cse_p12_t52_g60 + cse_p12_t60_g68);
      double cse_p12_t179_g187 = (cse_p12_t54_g62 + cse_p12_t58_g66);
      double cse_p12_t194_g202 = (c2z + cse_p12_t115_g123 + cse_p12_t119_g127);
      double v2_sq = ((v2x * v2x) + (v2y * v2y) + (v2z * v2z));
      double cse_p11_t1_sqrt4 = sqrt(v2_sq);
      double cse_p11_t2_invsqrt5 = (1.0 / (cse_p11_t1_sqrt4));
      double cse_p12_t1_g231 = (cse_p11_t2_invsqrt5 * cse_p12_t177_g185);
      double cse_p12_t62_g70 = (cse_p12_t1_g231 * (-(v2x)));
      double cse_p12_t63_g71 = (cse_p12_t1_g231 * (-(v2y)));
      double cse_p12_t64_g72 = (cse_p12_t1_g231 * (-(v2z)));
      double cse_p12_t72_g80 = (cse_p12_t1_g231 * v2x);
      double cse_p12_t73_g81 = (cse_p12_t1_g231 * v2y);
      double cse_p12_t74_g82 = (cse_p12_t1_g231 * v2z);
      double cse_p12_t201_g209 = (-(cse_p12_t1_g231));
      double cse_p11_t3_invr6 = (cse_p11_t2_invsqrt5 * cse_p11_t2_invsqrt5);
      double cse_p12_t66_g74 = (cse_p11_t3_invr6 * cse_p12_t1_g231 * v2x * v2x);
      double cse_p12_t67_g75 = (cse_p11_t3_invr6 * cse_p12_t1_g231 * v2x * v2y);
      double cse_p12_t68_g76 = (cse_p11_t3_invr6 * cse_p12_t1_g231 * v2x * v2z);
      double cse_p12_t69_g77 = (cse_p11_t3_invr6 * cse_p12_t1_g231 * v2y * v2y);
      double cse_p12_t70_g78 = (cse_p11_t3_invr6 * cse_p12_t1_g231 * v2y * v2z);
      double cse_p12_t71_g79 = (cse_p11_t3_invr6 * cse_p12_t1_g231 * v2z * v2z);
      double v2_len = cse_p11_t1_sqrt4;
      double cse_p12_t1_g237 = (cse_p12_t180_g188 * v2_len);
      double cse_p12_t1_g238 = (cse_p12_t181_g189 * v2_len);
      double cse_p12_t130_g138 = (cse_p12_t182_g190 * v2_len);
      double cse_p12_t131_g139 = (cse_p12_t183_g191 * v2_len);
      double cse_p12_t1_g239 = (cse_p12_t184_g192 * v2_len);
      double cse_p12_t133_g141 = (cse_p12_t185_g193 * v2_len);
      double cse_p12_t134_g142 = (cse_p12_t192_g200 * v2_len);
      double cse_p12_t135_g143 = (cse_p12_t193_g201 * v2_len);
      double cse_p12_t136_g144 = (cse_p12_t194_g202 * v2_len);
      double cse_p13_t1_g242 = (v2_len * (-(c2z)));
      double cse_p12_t147_g155 = (cse_p12_t134_g142 + cse_p12_t62_g70);
      double cse_p12_t148_g156 = (cse_p12_t135_g143 + cse_p12_t63_g71);
      double cse_p12_t149_g157 = (cse_p12_t136_g144 + cse_p12_t64_g72);
      double cse_p12_t150_g158 = (cse_p12_t1_g239 + cse_p12_t72_g80);
      double cse_p12_t151_g159 = (cse_p12_t1_g238 + cse_p12_t73_g81);
      double cse_p12_t152_g160 = (cse_p12_t1_g237 + cse_p12_t74_g82);
      double cse_p12_t213_g221 = (-(cse_p12_t130_g138));
      double cse_p12_t214_g222 = (-(cse_p12_t131_g139));
      double cse_p12_t215_g223 = (-(cse_p12_t133_g141));
      double cse_p12_t216_g224 = (-(cse_p12_t147_g155));
      double cse_p12_t217_g225 = (-(cse_p12_t148_g156));
      double cse_p12_t218_g226 = (-(cse_p12_t149_g157));
      double cse_p12_t219_g227 = (-(cse_p12_t150_g158));
      double cse_p12_t220_g228 = (-(cse_p12_t151_g159));
      double cse_p12_t221_g229 = (-(cse_p12_t152_g160));
      double cse_p12_t222_g230 = (-(cse_p13_t1_g242));
      double t1 = (cse_p12_t177_g185 * v2_len);
      double cse_p13_t1_g243 = (2.00000000000000000e+0 * t1);
      double cse_p12_t1_g9 = (cse_p12_t147_g155 * cse_p13_t1_g243);
      double cse_p12_t2_g10 = (cse_p12_t148_g156 * cse_p13_t1_g243);
      double cse_p12_t3_g11 = (cse_p12_t149_g157 * cse_p13_t1_g243);
      double cse_p12_t4_g12 = (cse_p12_t150_g158 * cse_p13_t1_g243);
      double cse_p12_t5_g13 = (cse_p12_t151_g159 * cse_p13_t1_g243);
      double cse_p12_t6_g14 = (cse_p12_t152_g160 * cse_p13_t1_g243);
      double cse_p12_t7_g15 = (cse_p12_t182_g190 * cse_p13_t1_g243 * v2_len);
      double cse_p12_t8_g16 = (cse_p12_t183_g191 * cse_p13_t1_g243 * v2_len);
      double cse_p12_t9_g17 = (cse_p12_t185_g193 * cse_p13_t1_g243 * v2_len);
      double cse_p12_t1_g240 = (v2_len * (-(c2y)));
      double cse_p12_t10_g18 = (cse_p12_t1_g240 * cse_p13_t1_g243);
      double cse_p12_t11_g19 = (cse_p13_t1_g242 * cse_p13_t1_g243);
      double t2 = ((c1x * c2x) + (c1y * c2y) + (c1z * c2z));
      double cse_p13_t1_g244 = (2.00000000000000000e+0 * t2);
      double cse_p12_t12_g20 = (cse_p12_t167_g175 * cse_p13_t1_g244);
      double cse_p12_t13_g21 = (cse_p12_t168_g176 * cse_p13_t1_g244);
      double cse_p12_t14_g22 = (cse_p12_t169_g177 * cse_p13_t1_g244);
      double cse_p12_t15_g23 = (cse_p12_t170_g178 * cse_p13_t1_g244);
      double cse_p12_t16_g24 = (cse_p12_t171_g179 * cse_p13_t1_g244);
      double cse_p12_t17_g25 = (cse_p12_t172_g180 * cse_p13_t1_g244);
      double cse_p12_t18_g26 = (cse_p12_t173_g181 * cse_p13_t1_g244);
      double cse_p12_t19_g27 = (cse_p12_t174_g182 * cse_p13_t1_g244);
      double cse_p12_t20_g28 = (cse_p12_t175_g183 * cse_p13_t1_g244);
      double cse_p12_t21_g29 = (cse_p12_t176_g184 * cse_p13_t1_g244);
      double cse_p12_t22_g30 = (cse_p12_t178_g186 * cse_p13_t1_g244);
      double cse_p12_t138_g146 = (cse_p12_t1_g9 + cse_p12_t20_g28);
      double cse_p12_t139_g147 = (cse_p12_t15_g23 + cse_p12_t2_g10);
      double cse_p12_t140_g148 = (cse_p12_t17_g25 + cse_p12_t3_g11);
      double cse_p12_t141_g149 = (cse_p12_t18_g26 + cse_p12_t4_g12);
      double cse_p12_t142_g150 = (cse_p12_t13_g21 + cse_p12_t5_g13);
      double cse_p12_t143_g151 = (cse_p12_t12_g20 + cse_p12_t6_g14);
      double cse_p12_t144_g152 = (cse_p12_t14_g22 + cse_p12_t7_g15);
      double cse_p12_t145_g153 = (cse_p12_t16_g24 + cse_p12_t8_g16);
      double cse_p12_t146_g154 = (cse_p12_t19_g27 + cse_p12_t9_g17);
      double cse_p12_t165_g173 = (cse_p12_t10_g18 + cse_p12_t21_g29);
      double cse_p12_t166_g174 = (cse_p12_t11_g19 + cse_p12_t22_g30);
      double phi = atan2(t1, t2);
      double nphi = (term.n * phi);
      double sin_nphi = sin(nphi);
      double cos_nphi = cos(nphi);
      double cos_angle = ((term.cosphase * cos_nphi) + (term.sinphase * sin_nphi));
      double sin_angle = ((term.cosphase * sin_nphi) + (-((cos_nphi * term.sinphase))));
      double energy = (term.v * (1.00000000000000000e+0 + cos_angle));
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t4_invr7 = (1.0 / ((((t1) * (t1)) + ((t2) * (t2)))));
      double cse_p13_t1_g247 = (cse_p11_t4_invr7 * (-(t1)));
      double cse_p12_t75_g83 = (cse_p12_t167_g175 * cse_p13_t1_g247);
      double cse_p12_t76_g84 = (cse_p12_t168_g176 * cse_p13_t1_g247);
      double cse_p12_t77_g85 = (cse_p12_t169_g177 * cse_p13_t1_g247);
      double cse_p12_t78_g86 = (cse_p12_t170_g178 * cse_p13_t1_g247);
      double cse_p12_t79_g87 = (cse_p12_t171_g179 * cse_p13_t1_g247);
      double cse_p12_t80_g88 = (cse_p12_t172_g180 * cse_p13_t1_g247);
      double cse_p12_t81_g89 = (cse_p12_t173_g181 * cse_p13_t1_g247);
      double cse_p12_t82_g90 = (cse_p12_t174_g182 * cse_p13_t1_g247);
      double cse_p12_t83_g91 = (cse_p12_t175_g183 * cse_p13_t1_g247);
      double cse_p12_t84_g92 = (cse_p12_t176_g184 * cse_p13_t1_g247);
      double cse_p12_t85_g93 = (cse_p12_t178_g186 * cse_p13_t1_g247);
      double cse_p12_t86_g94 = (cse_p12_t179_g187 * cse_p13_t1_g247);
      double cse_p12_t1_g234 = (cse_p11_t4_invr7 * t2);
      double cse_p12_t87_g95 = (cse_p12_t147_g155 * cse_p12_t1_g234);
      double cse_p12_t88_g96 = (cse_p12_t148_g156 * cse_p12_t1_g234);
      double cse_p12_t89_g97 = (cse_p12_t149_g157 * cse_p12_t1_g234);
      double cse_p12_t90_g98 = (cse_p12_t150_g158 * cse_p12_t1_g234);
      double cse_p12_t91_g99 = (cse_p12_t151_g159 * cse_p12_t1_g234);
      double cse_p12_t92_g100 = (cse_p12_t152_g160 * cse_p12_t1_g234);
      double cse_p13_t1_g246 = (cse_p12_t1_g234 * v2_len);
      double cse_p12_t93_g101 = (cse_p12_t182_g190 * cse_p13_t1_g246);
      double cse_p12_t94_g102 = (cse_p12_t183_g191 * cse_p13_t1_g246);
      double cse_p12_t95_g103 = (cse_p12_t185_g193 * cse_p13_t1_g246);
      double cse_p13_t1_g241 = (v2_len * (-(c2x)));
      double cse_p12_t96_g104 = (cse_p12_t1_g234 * cse_p13_t1_g241);
      double cse_p12_t97_g105 = (cse_p12_t1_g234 * cse_p12_t1_g240);
      double cse_p12_t98_g106 = (cse_p12_t1_g234 * cse_p13_t1_g242);
      double cse_p12_t153_g161 = (cse_p12_t75_g83 + cse_p12_t92_g100);
      double cse_p12_t154_g162 = (cse_p12_t76_g84 + cse_p12_t91_g99);
      double cse_p12_t155_g163 = (cse_p12_t77_g85 + cse_p12_t93_g101);
      double cse_p12_t156_g164 = (cse_p12_t78_g86 + cse_p12_t88_g96);
      double cse_p12_t157_g165 = (cse_p12_t79_g87 + cse_p12_t94_g102);
      double cse_p12_t158_g166 = (cse_p12_t80_g88 + cse_p12_t89_g97);
      double cse_p12_t159_g167 = (cse_p12_t81_g89 + cse_p12_t90_g98);
      double cse_p12_t160_g168 = (cse_p12_t82_g90 + cse_p12_t95_g103);
      double cse_p12_t161_g169 = (cse_p12_t83_g91 + cse_p12_t87_g95);
      double cse_p12_t162_g170 = (cse_p12_t84_g92 + cse_p12_t97_g105);
      double cse_p12_t163_g171 = (cse_p12_t85_g93 + cse_p12_t98_g106);
      double cse_p12_t164_g172 = (cse_p12_t86_g94 + cse_p12_t96_g104);
      double cse_p11_t5_invr28 = (cse_p11_t4_invr7 * cse_p11_t4_invr7);
      double cse_p12_t99_g107 = (cse_p11_t5_invr28 * cse_p12_t138_g146);
      double cse_p12_t100_g108 = (cse_p11_t5_invr28 * cse_p12_t139_g147);
      double cse_p12_t101_g109 = (cse_p11_t5_invr28 * cse_p12_t140_g148);
      double cse_p12_t102_g110 = (cse_p11_t5_invr28 * cse_p12_t141_g149);
      double cse_p12_t103_g111 = (cse_p11_t5_invr28 * cse_p12_t142_g150);
      double cse_p12_t104_g112 = (cse_p11_t5_invr28 * cse_p12_t143_g151);
      double cse_p12_t105_g113 = (cse_p11_t5_invr28 * cse_p12_t144_g152);
      double cse_p12_t106_g114 = (cse_p11_t5_invr28 * cse_p12_t145_g153);
      double cse_p12_t107_g115 = (cse_p11_t5_invr28 * cse_p12_t146_g154);
      double cse_p12_t108_g116 = (cse_p11_t5_invr28 * cse_p12_t165_g173);
      double cse_p12_t109_g117 = (cse_p11_t5_invr28 * cse_p12_t166_g174);
      double cse_p12_t202_g210 = (-(cse_p12_t99_g107));
      double cse_p12_t203_g211 = (-(cse_p12_t100_g108));
      double cse_p12_t204_g212 = (-(cse_p12_t101_g109));
      double cse_p12_t205_g213 = (-(cse_p12_t102_g110));
      double cse_p12_t206_g214 = (-(cse_p12_t103_g111));
      double cse_p12_t207_g215 = (-(cse_p12_t104_g112));
      double cse_p12_t208_g216 = (-(cse_p12_t105_g113));
      double cse_p12_t209_g217 = (-(cse_p12_t106_g114));
      double cse_p12_t210_g218 = (-(cse_p12_t107_g115));
      double cse_p12_t211_g219 = (-(cse_p12_t108_g116));
      double cse_p12_t212_g220 = (-(cse_p12_t109_g117));
      double cse_p12_t1_g232 = (term.n * (-(term.v)));
      double cse_p12_t1_g236 = (cse_p12_t1_g232 * sin_angle);
      double cse_p14_t1_g255 = (cse_p11_t4_invr7 * cse_p13_t1_g241);
      double g_x1 = (cse_p12_t1_g236 * (cse_p12_t86_g94 + (cse_p14_t1_g255 * t2)));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double cse_p14_t1_g256 = (cse_p11_t4_invr7 * cse_p12_t1_g240);
      double g_y1 = (cse_p12_t1_g236 * (cse_p12_t84_g92 + (cse_p14_t1_g256 * t2)));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double cse_p15_t1_g262 = (cse_p11_t4_invr7 * cse_p13_t1_g242);
      double g_z1 = (cse_p12_t1_g236 * (cse_p12_t85_g93 + (cse_p15_t1_g262 * t2)));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double g_x2 = (cse_p12_t1_g236 * (cse_p12_t83_g91 + (cse_p11_t4_invr7 * t2 * (cse_p12_t134_g142 + (cse_p11_t2_invsqrt5 * cse_p12_t177_g185 * (-(v2x)))))));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (cse_p12_t1_g236 * (cse_p12_t78_g86 + (cse_p11_t4_invr7 * t2 * (cse_p12_t135_g143 + (cse_p11_t2_invsqrt5 * cse_p12_t177_g185 * (-(v2y)))))));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (cse_p12_t1_g236 * (cse_p12_t80_g88 + (cse_p11_t4_invr7 * t2 * (cse_p12_t136_g144 + (cse_p11_t2_invsqrt5 * cse_p12_t177_g185 * (-(v2z)))))));
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double cse_p22_t1_g267 = (cse_p11_t2_invsqrt5 * v2x);
      double g_x3 = (cse_p12_t1_g236 * (cse_p12_t81_g89 + (cse_p11_t4_invr7 * t2 * (cse_p12_t1_g239 + (cse_p12_t177_g185 * cse_p22_t1_g267)))));
      KernelGradientAcc(term.i3x3, 0, g_x3);
      double cse_p22_t1_g268 = (cse_p11_t2_invsqrt5 * v2y);
      double g_y3 = (cse_p12_t1_g236 * (cse_p12_t76_g84 + (cse_p11_t4_invr7 * t2 * (cse_p12_t1_g238 + (cse_p12_t177_g185 * cse_p22_t1_g268)))));
      KernelGradientAcc(term.i3x3, 1, g_y3);
      double cse_p22_t1_g269 = (cse_p11_t2_invsqrt5 * v2z);
      double g_z3 = (cse_p12_t1_g236 * (cse_p12_t75_g83 + (cse_p11_t4_invr7 * t2 * (cse_p12_t1_g237 + (cse_p12_t177_g185 * cse_p22_t1_g269)))));
      KernelGradientAcc(term.i3x3, 2, g_z3);
      double g_x4 = (cse_p12_t1_g236 * (cse_p12_t82_g90 + (cse_p11_t4_invr7 * cse_p12_t185_g193 * t2 * v2_len)));
      KernelGradientAcc(term.i3x4, 0, g_x4);
      double g_y4 = (cse_p12_t1_g236 * (cse_p12_t77_g85 + (cse_p11_t4_invr7 * cse_p12_t182_g190 * t2 * v2_len)));
      KernelGradientAcc(term.i3x4, 1, g_y4);
      double g_z4 = (cse_p12_t1_g236 * (cse_p12_t79_g87 + (cse_p11_t4_invr7 * cse_p12_t183_g191 * t2 * v2_len)));
      KernelGradientAcc(term.i3x4, 2, g_z4);
      double cse_p12_t1_g235 = (cos_angle * cse_p12_t1_g232 * term.n);
      double h_x1_x1 = ((cse_p12_t164_g172 * cse_p12_t164_g172 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t179_g187 * (cse_p14_t1_g255 + (cse_p11_t4_invr7 * (-(cse_p13_t1_g241))) + (-((cse_p11_t5_invr28 * ((cse_p12_t179_g187 * cse_p13_t1_g244) + (cse_p13_t1_g241 * cse_p13_t1_g243)) * (-(t1))))))) + (-((cse_p11_t5_invr28 * cse_p13_t1_g241 * t2 * ((cse_p12_t179_g187 * cse_p13_t1_g244) + (cse_p13_t1_g241 * cse_p13_t1_g243))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
      double h_x1_y1 = ((cse_p12_t162_g170 * cse_p12_t164_g172 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t179_g187 * (-(cse_p12_t1_g240))) + (cse_p12_t176_g184 * cse_p14_t1_g255) + (cse_p12_t179_g187 * cse_p12_t211_g219 * (-(t1))) + (cse_p12_t211_g219 * cse_p13_t1_g241 * t2))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
      double h_x1_z1 = ((cse_p12_t163_g171 * cse_p12_t164_g172 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t179_g187 * cse_p12_t222_g230) + (cse_p12_t178_g186 * cse_p14_t1_g255) + (cse_p12_t179_g187 * cse_p12_t212_g220 * (-(t1))) + (cse_p12_t212_g220 * cse_p13_t1_g241 * t2))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
      double cse_p13_t1_g245 = (cse_p11_t2_invsqrt5 * cse_p12_t1_g234);
      double h_x1_x2 = ((cse_p12_t161_g169 * cse_p12_t164_g172 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t179_g187 * cse_p12_t216_g224) + (cse_p12_t175_g183 * cse_p14_t1_g255) + (cse_p12_t179_g187 * cse_p12_t202_g210 * (-(t1))) + (cse_p12_t202_g210 * cse_p13_t1_g241 * t2) + (cse_p13_t1_g247 * (((-(v2y)) * (-(v3y))) + (v2z * v3z))) + (-((cse_p13_t1_g245 * v2x * (-(c2x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
      double h_x1_y2 = ((cse_p12_t156_g164 * cse_p12_t164_g172 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t179_g187 * cse_p12_t217_g225) + (cse_p12_t170_g178 * cse_p14_t1_g255) + (cse_p12_t179_g187 * cse_p12_t203_g211 * (-(t1))) + (cse_p12_t203_g211 * cse_p13_t1_g241 * t2) + (cse_p13_t1_g246 * v3z) + (cse_p13_t1_g247 * (c2z + (v3x * (-(v2y))))) + (-((cse_p13_t1_g245 * v2y * (-(c2x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
      double h_x1_z2 = ((cse_p12_t158_g166 * cse_p12_t164_g172 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t179_g187 * cse_p12_t218_g226) + (cse_p12_t172_g180 * cse_p14_t1_g255) + (cse_p12_t179_g187 * cse_p12_t204_g212 * (-(t1))) + (cse_p12_t204_g212 * cse_p13_t1_g241 * t2) + (cse_p13_t1_g246 * (-(v3y))) + (cse_p13_t1_g247 * ((v2z * (-(v3x))) + (-(c2y)))) + (-((cse_p13_t1_g245 * v2z * (-(c2x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
      double cse_p21_t1_g265 = (cse_p12_t205_g213 * (-(t1)));
      double h_x1_x3 = ((cse_p12_t159_g167 * cse_p12_t164_g172 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t179_g187 * cse_p12_t219_g227) + (cse_p12_t173_g181 * cse_p14_t1_g255) + (cse_p12_t179_g187 * cse_p21_t1_g265) + (cse_p12_t205_g213 * cse_p13_t1_g241 * t2) + (cse_p13_t1_g245 * v2x * (-(c2x))) + (cse_p13_t1_g247 * ((cse_p12_t191_g199 * v2z) + (cse_p12_t199_g207 * (-(v2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 0, h_x1_x3);
      double cse_p14_t1_g260 = (cse_p12_t206_g214 * (-(t1)));
      double cse_p15_t1_g263 = (cse_p12_t154_g162 * cse_p12_t1_g235);
      double h_x1_y3 = ((cse_p12_t164_g172 * cse_p15_t1_g263) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t179_g187 * cse_p12_t220_g228) + (cse_p12_t168_g176 * cse_p14_t1_g255) + (cse_p12_t179_g187 * cse_p14_t1_g260) + (cse_p12_t191_g199 * cse_p13_t1_g246) + (cse_p12_t206_g214 * cse_p13_t1_g241 * t2) + (cse_p13_t1_g245 * v2y * (-(c2x))) + (cse_p13_t1_g247 * ((cse_p12_t189_g197 * (-(v2y))) + (-(c2z)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 1, h_x1_y3);
      double cse_p14_t1_g251 = (cse_p12_t207_g215 * (-(t1)));
      double cse_p14_t1_g258 = (cse_p12_t153_g161 * cse_p12_t1_g235);
      double h_x1_z3 = ((cse_p12_t164_g172 * cse_p14_t1_g258) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t179_g187 * cse_p12_t221_g229) + (cse_p12_t167_g175 * cse_p14_t1_g255) + (cse_p12_t179_g187 * cse_p14_t1_g251) + (cse_p12_t199_g207 * cse_p13_t1_g246) + (cse_p12_t207_g215 * cse_p13_t1_g241 * t2) + (cse_p13_t1_g245 * v2z * (-(c2x))) + (cse_p13_t1_g247 * (c2y + (cse_p12_t198_g206 * v2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 2, h_x1_z3);
      double cse_p13_t1_g250 = (cse_p12_t210_g218 * (-(t1)));
      double cse_p21_t1_g266 = (cse_p11_t4_invr7 * cse_p12_t215_g223);
      double h_x1_x4 = ((cse_p12_t160_g168 * cse_p12_t164_g172 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t174_g182 * cse_p14_t1_g255) + (cse_p12_t179_g187 * cse_p13_t1_g250) + (cse_p12_t179_g187 * cse_p21_t1_g266) + (cse_p12_t210_g218 * cse_p13_t1_g241 * t2) + (cse_p13_t1_g247 * (((-(v2y)) * (-(v2y))) + (v2z * v2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x4, 0, h_x1_x4);
      double cse_p13_t1_g249 = (cse_p12_t208_g216 * (-(t1)));
      double cse_p15_t1_g261 = (cse_p11_t4_invr7 * cse_p12_t213_g221);
      double cse_p21_t1_g264 = (cse_p12_t208_g216 * t2);
      double h_x1_y4 = ((cse_p12_t155_g163 * cse_p12_t164_g172 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t169_g177 * cse_p14_t1_g255) + (cse_p12_t179_g187 * cse_p13_t1_g249) + (cse_p12_t179_g187 * cse_p15_t1_g261) + (cse_p13_t1_g241 * cse_p21_t1_g264) + (cse_p13_t1_g246 * v2z) + (cse_p13_t1_g247 * v2x * (-(v2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x4, 1, h_x1_y4);
      double cse_p13_t1_g248 = (cse_p12_t209_g217 * (-(t1)));
      double cse_p14_t1_g257 = (cse_p11_t4_invr7 * cse_p12_t214_g222);
      double cse_p14_t1_g259 = (cse_p12_t209_g217 * t2);
      double h_x1_z4 = ((cse_p12_t157_g165 * cse_p12_t164_g172 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t171_g179 * cse_p14_t1_g255) + (cse_p12_t179_g187 * cse_p13_t1_g248) + (cse_p12_t179_g187 * cse_p14_t1_g257) + (cse_p13_t1_g241 * cse_p14_t1_g259) + (cse_p13_t1_g246 * (-(v2y))) + (cse_p13_t1_g247 * v2z * (-(v2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x4, 2, h_x1_z4);
      double h_y1_y1 = ((cse_p12_t162_g170 * cse_p12_t162_g170 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t176_g184 * (cse_p14_t1_g256 + (cse_p11_t4_invr7 * (-(cse_p12_t1_g240))) + (cse_p12_t211_g219 * (-(t1))))) + (cse_p12_t1_g240 * cse_p12_t211_g219 * t2))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
      double h_y1_z1 = ((cse_p12_t162_g170 * cse_p12_t163_g171 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t176_g184 * cse_p12_t222_g230) + (cse_p12_t176_g184 * cse_p12_t212_g220 * (-(t1))) + (cse_p12_t178_g186 * cse_p14_t1_g256) + (cse_p12_t1_g240 * cse_p12_t212_g220 * t2))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
      double h_y1_x2 = ((cse_p12_t161_g169 * cse_p12_t162_g170 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t176_g184 * cse_p12_t216_g224) + (cse_p12_t175_g183 * cse_p14_t1_g256) + (cse_p12_t176_g184 * cse_p12_t202_g210 * (-(t1))) + (cse_p12_t1_g240 * cse_p12_t202_g210 * t2) + (cse_p13_t1_g246 * (-(v3z))) + (cse_p13_t1_g247 * ((v2x * (-(v3y))) + (-(c2z)))) + (-((cse_p13_t1_g245 * v2x * (-(c2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
      double h_y1_y2 = ((cse_p12_t156_g164 * cse_p12_t162_g170 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t176_g184 * cse_p12_t217_g225) + (cse_p12_t170_g178 * cse_p14_t1_g256) + (cse_p12_t176_g184 * cse_p12_t203_g211 * (-(t1))) + (cse_p12_t1_g240 * cse_p12_t203_g211 * t2) + (cse_p13_t1_g247 * (((-(v2z)) * (-(v3z))) + (v2x * v3x))) + (-((cse_p13_t1_g245 * v2y * (-(c2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
      double h_y1_z2 = ((cse_p12_t158_g166 * cse_p12_t162_g170 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t176_g184 * cse_p12_t218_g226) + (cse_p12_t172_g180 * cse_p14_t1_g256) + (cse_p12_t176_g184 * cse_p12_t204_g212 * (-(t1))) + (cse_p12_t1_g240 * cse_p12_t204_g212 * t2) + (cse_p13_t1_g246 * v3x) + (cse_p13_t1_g247 * (c2x + (v3y * (-(v2z))))) + (-((cse_p13_t1_g245 * v2z * (-(c2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
      double h_y1_x3 = ((cse_p12_t159_g167 * cse_p12_t162_g170 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t176_g184 * cse_p12_t219_g227) + (cse_p12_t173_g181 * cse_p14_t1_g256) + (cse_p12_t176_g184 * cse_p21_t1_g265) + (cse_p12_t1_g240 * cse_p12_t205_g213 * t2) + (cse_p12_t200_g208 * cse_p13_t1_g246) + (cse_p13_t1_g245 * v2x * (-(c2y))) + (cse_p13_t1_g247 * (c2z + (cse_p12_t199_g207 * v2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 0, h_y1_x3);
      double h_y1_y3 = ((cse_p12_t162_g170 * cse_p15_t1_g263) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t176_g184 * cse_p12_t220_g228) + (cse_p12_t168_g176 * cse_p14_t1_g256) + (cse_p12_t176_g184 * cse_p14_t1_g260) + (cse_p12_t1_g240 * cse_p12_t206_g214 * t2) + (cse_p13_t1_g245 * v2y * (-(c2y))) + (cse_p13_t1_g247 * ((cse_p12_t189_g197 * v2x) + (cse_p12_t200_g208 * (-(v2z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 1, h_y1_y3);
      double h_y1_z3 = ((cse_p12_t162_g170 * cse_p14_t1_g258) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t176_g184 * cse_p12_t221_g229) + (cse_p12_t167_g175 * cse_p14_t1_g256) + (cse_p12_t176_g184 * cse_p14_t1_g251) + (cse_p12_t189_g197 * cse_p13_t1_g246) + (cse_p12_t1_g240 * cse_p12_t207_g215 * t2) + (cse_p13_t1_g245 * v2z * (-(c2y))) + (cse_p13_t1_g247 * ((cse_p12_t190_g198 * (-(v2z))) + (-(c2x)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 2, h_y1_z3);
      double h_y1_x4 = ((cse_p12_t160_g168 * cse_p12_t162_g170 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t174_g182 * cse_p14_t1_g256) + (cse_p12_t176_g184 * cse_p13_t1_g250) + (cse_p12_t176_g184 * cse_p21_t1_g266) + (cse_p12_t1_g240 * cse_p12_t210_g218 * t2) + (cse_p13_t1_g246 * (-(v2z))) + (cse_p13_t1_g247 * v2x * (-(v2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x4, 0, h_y1_x4);
      double h_y1_y4 = ((cse_p12_t155_g163 * cse_p12_t162_g170 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t169_g177 * cse_p14_t1_g256) + (cse_p12_t176_g184 * cse_p13_t1_g249) + (cse_p12_t176_g184 * cse_p15_t1_g261) + (cse_p12_t1_g240 * cse_p21_t1_g264) + (cse_p13_t1_g247 * (((-(v2z)) * (-(v2z))) + (v2x * v2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x4, 1, h_y1_y4);
      double h_y1_z4 = ((cse_p12_t157_g165 * cse_p12_t162_g170 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t171_g179 * cse_p14_t1_g256) + (cse_p12_t176_g184 * cse_p13_t1_g248) + (cse_p12_t176_g184 * cse_p14_t1_g257) + (cse_p12_t1_g240 * cse_p14_t1_g259) + (cse_p13_t1_g246 * v2x) + (cse_p13_t1_g247 * v2y * (-(v2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x4, 2, h_y1_z4);
      double h_z1_z1 = ((cse_p12_t163_g171 * cse_p12_t163_g171 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t178_g186 * (cse_p15_t1_g262 + (cse_p11_t4_invr7 * cse_p12_t222_g230) + (cse_p12_t212_g220 * (-(t1))))) + (cse_p12_t212_g220 * cse_p13_t1_g242 * t2))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
      double h_z1_x2 = ((cse_p12_t161_g169 * cse_p12_t163_g171 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t178_g186 * cse_p12_t216_g224) + (cse_p12_t175_g183 * cse_p15_t1_g262) + (cse_p12_t178_g186 * cse_p12_t202_g210 * (-(t1))) + (cse_p12_t202_g210 * cse_p13_t1_g242 * t2) + (cse_p13_t1_g246 * v3y) + (cse_p13_t1_g247 * (c2y + (v3z * (-(v2x))))) + (-((cse_p13_t1_g245 * v2x * (-(c2z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
      double h_z1_y2 = ((cse_p12_t156_g164 * cse_p12_t163_g171 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t178_g186 * cse_p12_t217_g225) + (cse_p12_t170_g178 * cse_p15_t1_g262) + (cse_p12_t178_g186 * cse_p12_t203_g211 * (-(t1))) + (cse_p12_t203_g211 * cse_p13_t1_g242 * t2) + (cse_p13_t1_g246 * (-(v3x))) + (cse_p13_t1_g247 * ((v2y * (-(v3z))) + (-(c2x)))) + (-((cse_p13_t1_g245 * v2y * (-(c2z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
      double h_z1_z2 = ((cse_p12_t158_g166 * cse_p12_t163_g171 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t178_g186 * cse_p12_t218_g226) + (cse_p12_t172_g180 * cse_p15_t1_g262) + (cse_p12_t178_g186 * cse_p12_t204_g212 * (-(t1))) + (cse_p12_t204_g212 * cse_p13_t1_g242 * t2) + (cse_p13_t1_g247 * (((-(v2x)) * (-(v3x))) + (v2y * v3y))) + (-((cse_p13_t1_g245 * v2z * (-(c2z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
      double h_z1_x3 = ((cse_p12_t159_g167 * cse_p12_t163_g171 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t178_g186 * cse_p12_t219_g227) + (cse_p12_t173_g181 * cse_p15_t1_g262) + (cse_p12_t178_g186 * cse_p21_t1_g265) + (cse_p12_t190_g198 * cse_p13_t1_g246) + (cse_p12_t205_g213 * cse_p13_t1_g242 * t2) + (cse_p13_t1_g245 * v2x * (-(c2z))) + (cse_p13_t1_g247 * ((cse_p12_t191_g199 * (-(v2x))) + (-(c2y)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 0, h_z1_x3);
      double h_z1_y3 = ((cse_p12_t163_g171 * cse_p15_t1_g263) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t178_g186 * cse_p12_t220_g228) + (cse_p12_t168_g176 * cse_p15_t1_g262) + (cse_p12_t178_g186 * cse_p14_t1_g260) + (cse_p12_t198_g206 * cse_p13_t1_g246) + (cse_p12_t206_g214 * cse_p13_t1_g242 * t2) + (cse_p13_t1_g245 * v2y * (-(c2z))) + (cse_p13_t1_g247 * (c2x + (cse_p12_t200_g208 * v2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 1, h_z1_y3);
      double h_z1_z3 = ((cse_p12_t163_g171 * cse_p14_t1_g258) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t178_g186 * cse_p12_t221_g229) + (cse_p12_t167_g175 * cse_p15_t1_g262) + (cse_p12_t178_g186 * cse_p14_t1_g251) + (cse_p12_t207_g215 * cse_p13_t1_g242 * t2) + (cse_p13_t1_g245 * v2z * (-(c2z))) + (cse_p13_t1_g247 * ((cse_p12_t190_g198 * v2y) + (cse_p12_t198_g206 * (-(v2x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 2, h_z1_z3);
      double h_z1_x4 = ((cse_p12_t160_g168 * cse_p12_t163_g171 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t174_g182 * cse_p15_t1_g262) + (cse_p12_t178_g186 * cse_p13_t1_g250) + (cse_p12_t178_g186 * cse_p21_t1_g266) + (cse_p12_t210_g218 * cse_p13_t1_g242 * t2) + (cse_p13_t1_g246 * v2y) + (cse_p13_t1_g247 * v2z * (-(v2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x4, 0, h_z1_x4);
      double h_z1_y4 = ((cse_p12_t155_g163 * cse_p12_t163_g171 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t169_g177 * cse_p15_t1_g262) + (cse_p12_t178_g186 * cse_p13_t1_g249) + (cse_p12_t178_g186 * cse_p15_t1_g261) + (cse_p13_t1_g242 * cse_p21_t1_g264) + (cse_p13_t1_g246 * (-(v2x))) + (cse_p13_t1_g247 * v2y * (-(v2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x4, 1, h_z1_y4);
      double h_z1_z4 = ((cse_p12_t157_g165 * cse_p12_t163_g171 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t171_g179 * cse_p15_t1_g262) + (cse_p12_t178_g186 * cse_p13_t1_g248) + (cse_p12_t178_g186 * cse_p14_t1_g257) + (cse_p13_t1_g242 * cse_p14_t1_g259) + (cse_p13_t1_g247 * (((-(v2x)) * (-(v2x))) + (v2y * v2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x4, 2, h_z1_z4);
      double cse_p14_t1_g252 = (cse_p11_t2_invsqrt5 * (-(v2x)));
      double h_x2_x2 = ((cse_p12_t161_g169 * cse_p12_t161_g169 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t147_g155 * cse_p12_t202_g210 * t2) + (cse_p12_t175_g183 * ((cse_p11_t4_invr7 * cse_p12_t147_g155) + (cse_p11_t4_invr7 * cse_p12_t216_g224) + (cse_p12_t202_g210 * (-(t1))))) + (cse_p12_t1_g234 * (cse_p12_t1_g231 + (cse_p12_t192_g200 * cse_p14_t1_g252) + (cse_p12_t192_g200 * cse_p14_t1_g252) + (-(cse_p12_t66_g74)))) + (cse_p13_t1_g247 * ((cse_p12_t188_g196 * v3z) + (cse_p12_t188_g196 * v3z) + (cse_p12_t196_g204 * (-(v3y))) + (cse_p12_t196_g204 * (-(v3y))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
      double cse_p14_t1_g253 = (cse_p11_t2_invsqrt5 * (-(v2y)));
      double h_x2_y2 = ((cse_p12_t156_g164 * cse_p12_t161_g169 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t147_g155 * cse_p12_t170_g178) + (cse_p11_t4_invr7 * cse_p12_t175_g183 * cse_p12_t217_g225) + (cse_p12_t147_g155 * cse_p12_t203_g211 * t2) + (cse_p12_t175_g183 * cse_p12_t203_g211 * (-(t1))) + (cse_p12_t1_g234 * ((cse_p12_t192_g200 * cse_p14_t1_g253) + (cse_p12_t193_g201 * cse_p14_t1_g252) + (-(cse_p12_t67_g75)))) + (cse_p13_t1_g247 * ((cse_p12_t186_g194 * (-(v3y))) + (cse_p12_t196_g204 * v3x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
      double cse_p14_t1_g254 = (cse_p11_t2_invsqrt5 * (-(v2z)));
      double h_x2_z2 = ((cse_p12_t158_g166 * cse_p12_t161_g169 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t147_g155 * cse_p12_t172_g180) + (cse_p11_t4_invr7 * cse_p12_t175_g183 * cse_p12_t218_g226) + (cse_p12_t147_g155 * cse_p12_t204_g212 * t2) + (cse_p12_t175_g183 * cse_p12_t204_g212 * (-(t1))) + (cse_p12_t1_g234 * ((cse_p12_t192_g200 * cse_p14_t1_g254) + (cse_p12_t194_g202 * cse_p14_t1_g252) + (-(cse_p12_t68_g76)))) + (cse_p13_t1_g247 * ((cse_p12_t188_g196 * (-(v3x))) + (cse_p12_t195_g203 * v3z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
      double h_x2_x3 = ((cse_p12_t159_g167 * cse_p12_t161_g169 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t147_g155 * cse_p12_t173_g181) + (cse_p11_t4_invr7 * cse_p12_t175_g183 * cse_p12_t219_g227) + (cse_p12_t147_g155 * cse_p12_t205_g213 * t2) + (cse_p12_t175_g183 * cse_p21_t1_g265) + (cse_p12_t1_g234 * (cse_p12_t201_g209 + cse_p12_t66_g74 + (cse_p12_t184_g192 * cse_p14_t1_g252) + (cse_p12_t192_g200 * cse_p22_t1_g267))) + (cse_p13_t1_g247 * (((-(v1y)) * (-(v3y))) + (cse_p12_t188_g196 * cse_p12_t191_g199) + (cse_p12_t196_g204 * cse_p12_t199_g207) + (v1z * v3z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 0, h_x2_x3);
      double h_x2_y3 = ((cse_p12_t161_g169 * cse_p15_t1_g263) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t147_g155 * cse_p12_t168_g176) + (cse_p11_t4_invr7 * cse_p12_t175_g183 * cse_p12_t220_g228) + (cse_p12_t147_g155 * cse_p12_t206_g214 * t2) + (cse_p12_t175_g183 * cse_p14_t1_g260) + (cse_p12_t1_g234 * (cse_p12_t67_g75 + (cse_p12_t181_g189 * cse_p14_t1_g252) + (cse_p12_t192_g200 * cse_p22_t1_g268) + (v2_len * (v1z + v2z + v3z)))) + (cse_p13_t1_g247 * (c1z + c2z + (cse_p12_t189_g197 * cse_p12_t196_g204) + (v1x * (-(v3y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 1, h_x2_y3);
      double h_x2_z3 = ((cse_p12_t161_g169 * cse_p14_t1_g258) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t147_g155 * cse_p12_t167_g175) + (cse_p11_t4_invr7 * cse_p12_t175_g183 * cse_p12_t221_g229) + (cse_p12_t147_g155 * cse_p12_t207_g215 * t2) + (cse_p12_t175_g183 * cse_p14_t1_g251) + (cse_p12_t1_g234 * (cse_p12_t68_g76 + (cse_p12_t180_g188 * cse_p14_t1_g252) + (cse_p12_t192_g200 * cse_p22_t1_g269) + (v2_len * ((-(v1y)) + (-(v2y)) + (-(v3y)))))) + (cse_p13_t1_g247 * ((cse_p12_t188_g196 * cse_p12_t198_g206) + (v3z * (-(v1x))) + (-(c1y)) + (-(c2y)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 2, h_x2_z3);
      double h_x2_x4 = ((cse_p12_t160_g168 * cse_p12_t161_g169 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t147_g155 * cse_p12_t174_g182) + (cse_p12_t147_g155 * cse_p12_t210_g218 * t2) + (cse_p12_t175_g183 * cse_p13_t1_g250) + (cse_p12_t175_g183 * cse_p21_t1_g266) + (cse_p13_t1_g247 * ((cse_p12_t188_g196 * v2z) + (cse_p12_t196_g204 * (-(v2y))))) + (-((cse_p12_t185_g193 * cse_p13_t1_g245 * v2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x4, 0, h_x2_x4);
      double h_x2_y4 = ((cse_p12_t155_g163 * cse_p12_t161_g169 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t147_g155 * cse_p12_t169_g177) + (cse_p12_t147_g155 * cse_p21_t1_g264) + (cse_p12_t175_g183 * cse_p13_t1_g249) + (cse_p12_t175_g183 * cse_p15_t1_g261) + (cse_p12_t1_g234 * ((cse_p12_t182_g190 * cse_p14_t1_g252) + (cse_p12_t188_g196 * v2_len))) + (cse_p13_t1_g247 * ((cse_p12_t196_g204 * v2x) + (-(c1z)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x4, 1, h_x2_y4);
      double h_x2_z4 = ((cse_p12_t157_g165 * cse_p12_t161_g169 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t147_g155 * cse_p12_t171_g179) + (cse_p12_t147_g155 * cse_p14_t1_g259) + (cse_p12_t175_g183 * cse_p13_t1_g248) + (cse_p12_t175_g183 * cse_p14_t1_g257) + (cse_p12_t1_g234 * ((cse_p12_t183_g191 * cse_p14_t1_g252) + (cse_p12_t196_g204 * v2_len))) + (cse_p13_t1_g247 * (c1y + (cse_p12_t188_g196 * (-(v2x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x4, 2, h_x2_z4);
      double h_y2_y2 = ((cse_p12_t156_g164 * cse_p12_t156_g164 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t148_g156 * cse_p12_t203_g211 * t2) + (cse_p12_t170_g178 * ((cse_p11_t4_invr7 * cse_p12_t148_g156) + (cse_p11_t4_invr7 * cse_p12_t217_g225) + (cse_p12_t203_g211 * (-(t1))))) + (cse_p12_t1_g234 * (cse_p12_t1_g231 + (cse_p12_t193_g201 * cse_p14_t1_g253) + (cse_p12_t193_g201 * cse_p14_t1_g253) + (-(cse_p12_t69_g77)))) + (cse_p13_t1_g247 * ((cse_p12_t186_g194 * v3x) + (cse_p12_t186_g194 * v3x) + (cse_p12_t197_g205 * (-(v3z))) + (cse_p12_t197_g205 * (-(v3z))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
      double h_y2_z2 = ((cse_p12_t156_g164 * cse_p12_t158_g166 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t148_g156 * cse_p12_t172_g180) + (cse_p11_t4_invr7 * cse_p12_t170_g178 * cse_p12_t218_g226) + (cse_p12_t148_g156 * cse_p12_t204_g212 * t2) + (cse_p12_t170_g178 * cse_p12_t204_g212 * (-(t1))) + (cse_p12_t1_g234 * ((cse_p12_t193_g201 * cse_p14_t1_g254) + (cse_p12_t194_g202 * cse_p14_t1_g253) + (-(cse_p12_t70_g78)))) + (cse_p13_t1_g247 * ((cse_p12_t187_g195 * (-(v3z))) + (cse_p12_t197_g205 * v3y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
      double h_y2_x3 = ((cse_p12_t156_g164 * cse_p12_t159_g167 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t148_g156 * cse_p12_t173_g181) + (cse_p11_t4_invr7 * cse_p12_t170_g178 * cse_p12_t219_g227) + (cse_p12_t148_g156 * cse_p12_t205_g213 * t2) + (cse_p12_t170_g178 * cse_p21_t1_g265) + (cse_p12_t1_g234 * (cse_p12_t67_g75 + (cse_p12_t184_g192 * cse_p14_t1_g253) + (cse_p12_t193_g201 * cse_p22_t1_g267) + (v2_len * ((-(v1z)) + (-(v2z)) + (-(v3z)))))) + (cse_p13_t1_g247 * ((cse_p12_t186_g194 * cse_p12_t199_g207) + (v3x * (-(v1y))) + (-(c1z)) + (-(c2z)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 0, h_y2_x3);
      double h_y2_y3 = ((cse_p12_t156_g164 * cse_p15_t1_g263) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t148_g156 * cse_p12_t168_g176) + (cse_p11_t4_invr7 * cse_p12_t170_g178 * cse_p12_t220_g228) + (cse_p12_t148_g156 * cse_p12_t206_g214 * t2) + (cse_p12_t170_g178 * cse_p14_t1_g260) + (cse_p12_t1_g234 * (cse_p12_t201_g209 + cse_p12_t69_g77 + (cse_p12_t181_g189 * cse_p14_t1_g253) + (cse_p12_t193_g201 * cse_p22_t1_g268))) + (cse_p13_t1_g247 * (((-(v1z)) * (-(v3z))) + (cse_p12_t186_g194 * cse_p12_t189_g197) + (cse_p12_t197_g205 * cse_p12_t200_g208) + (v1x * v3x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 1, h_y2_y3);
      double h_y2_z3 = ((cse_p12_t156_g164 * cse_p14_t1_g258) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t148_g156 * cse_p12_t167_g175) + (cse_p11_t4_invr7 * cse_p12_t170_g178 * cse_p12_t221_g229) + (cse_p12_t148_g156 * cse_p12_t207_g215 * t2) + (cse_p12_t170_g178 * cse_p14_t1_g251) + (cse_p12_t1_g234 * (cse_p12_t70_g78 + (cse_p12_t180_g188 * cse_p14_t1_g253) + (cse_p12_t193_g201 * cse_p22_t1_g269) + (v2_len * (v1x + v2x + v3x)))) + (cse_p13_t1_g247 * (c1x + c2x + (cse_p12_t190_g198 * cse_p12_t197_g205) + (v1y * (-(v3z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 2, h_y2_z3);
      double h_y2_x4 = ((cse_p12_t156_g164 * cse_p12_t160_g168 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t148_g156 * cse_p12_t174_g182) + (cse_p12_t148_g156 * cse_p12_t210_g218 * t2) + (cse_p12_t170_g178 * cse_p13_t1_g250) + (cse_p12_t170_g178 * cse_p21_t1_g266) + (cse_p12_t1_g234 * ((cse_p12_t185_g193 * cse_p14_t1_g253) + (cse_p12_t197_g205 * v2_len))) + (cse_p13_t1_g247 * (c1z + (cse_p12_t186_g194 * (-(v2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x4, 0, h_y2_x4);
      double h_y2_y4 = ((cse_p12_t155_g163 * cse_p12_t156_g164 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t148_g156 * cse_p12_t169_g177) + (cse_p12_t148_g156 * cse_p21_t1_g264) + (cse_p12_t170_g178 * cse_p13_t1_g249) + (cse_p12_t170_g178 * cse_p15_t1_g261) + (cse_p13_t1_g247 * ((cse_p12_t186_g194 * v2x) + (cse_p12_t197_g205 * (-(v2z))))) + (-((cse_p12_t182_g190 * cse_p13_t1_g245 * v2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x4, 1, h_y2_y4);
      double h_y2_z4 = ((cse_p12_t156_g164 * cse_p12_t157_g165 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t148_g156 * cse_p12_t171_g179) + (cse_p12_t148_g156 * cse_p14_t1_g259) + (cse_p12_t170_g178 * cse_p13_t1_g248) + (cse_p12_t170_g178 * cse_p14_t1_g257) + (cse_p12_t1_g234 * ((cse_p12_t183_g191 * cse_p14_t1_g253) + (cse_p12_t186_g194 * v2_len))) + (cse_p13_t1_g247 * ((cse_p12_t197_g205 * v2y) + (-(c1x)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x4, 2, h_y2_z4);
      double h_z2_z2 = ((cse_p12_t158_g166 * cse_p12_t158_g166 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t149_g157 * cse_p12_t204_g212 * t2) + (cse_p12_t172_g180 * ((cse_p11_t4_invr7 * cse_p12_t149_g157) + (cse_p11_t4_invr7 * cse_p12_t218_g226) + (cse_p12_t204_g212 * (-(t1))))) + (cse_p12_t1_g234 * (cse_p12_t1_g231 + (cse_p12_t194_g202 * cse_p14_t1_g254) + (cse_p12_t194_g202 * cse_p14_t1_g254) + (-(cse_p12_t71_g79)))) + (cse_p13_t1_g247 * ((cse_p12_t187_g195 * v3y) + (cse_p12_t187_g195 * v3y) + (cse_p12_t195_g203 * (-(v3x))) + (cse_p12_t195_g203 * (-(v3x))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      double h_z2_x3 = ((cse_p12_t158_g166 * cse_p12_t159_g167 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t149_g157 * cse_p12_t173_g181) + (cse_p11_t4_invr7 * cse_p12_t172_g180 * cse_p12_t219_g227) + (cse_p12_t149_g157 * cse_p12_t205_g213 * t2) + (cse_p12_t172_g180 * cse_p21_t1_g265) + (cse_p12_t1_g234 * (cse_p12_t68_g76 + (cse_p12_t184_g192 * cse_p14_t1_g254) + (cse_p12_t194_g202 * cse_p22_t1_g267) + (v2_len * (v1y + v2y + v3y)))) + (cse_p13_t1_g247 * (c1y + c2y + (cse_p12_t191_g199 * cse_p12_t195_g203) + (v1z * (-(v3x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 0, h_z2_x3);
      double h_z2_y3 = ((cse_p12_t158_g166 * cse_p15_t1_g263) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t149_g157 * cse_p12_t168_g176) + (cse_p11_t4_invr7 * cse_p12_t172_g180 * cse_p12_t220_g228) + (cse_p12_t149_g157 * cse_p12_t206_g214 * t2) + (cse_p12_t172_g180 * cse_p14_t1_g260) + (cse_p12_t1_g234 * (cse_p12_t70_g78 + (cse_p12_t181_g189 * cse_p14_t1_g254) + (cse_p12_t194_g202 * cse_p22_t1_g268) + (v2_len * ((-(v1x)) + (-(v2x)) + (-(v3x)))))) + (cse_p13_t1_g247 * ((cse_p12_t187_g195 * cse_p12_t200_g208) + (v3y * (-(v1z))) + (-(c1x)) + (-(c2x)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 1, h_z2_y3);
      double h_z2_z3 = ((cse_p12_t158_g166 * cse_p14_t1_g258) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t149_g157 * cse_p12_t167_g175) + (cse_p11_t4_invr7 * cse_p12_t172_g180 * cse_p12_t221_g229) + (cse_p12_t149_g157 * cse_p12_t207_g215 * t2) + (cse_p12_t172_g180 * cse_p14_t1_g251) + (cse_p12_t1_g234 * (cse_p12_t201_g209 + cse_p12_t71_g79 + (cse_p12_t180_g188 * cse_p14_t1_g254) + (cse_p12_t194_g202 * cse_p22_t1_g269))) + (cse_p13_t1_g247 * (((-(v1x)) * (-(v3x))) + (cse_p12_t187_g195 * cse_p12_t190_g198) + (cse_p12_t195_g203 * cse_p12_t198_g206) + (v1y * v3y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 2, h_z2_z3);
      double h_z2_x4 = ((cse_p12_t158_g166 * cse_p12_t160_g168 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t149_g157 * cse_p12_t174_g182) + (cse_p12_t149_g157 * cse_p12_t210_g218 * t2) + (cse_p12_t172_g180 * cse_p13_t1_g250) + (cse_p12_t172_g180 * cse_p21_t1_g266) + (cse_p12_t1_g234 * ((cse_p12_t185_g193 * cse_p14_t1_g254) + (cse_p12_t187_g195 * v2_len))) + (cse_p13_t1_g247 * ((cse_p12_t195_g203 * v2z) + (-(c1y)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x4, 0, h_z2_x4);
      double h_z2_y4 = ((cse_p12_t155_g163 * cse_p12_t158_g166 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t149_g157 * cse_p12_t169_g177) + (cse_p12_t149_g157 * cse_p21_t1_g264) + (cse_p12_t172_g180 * cse_p13_t1_g249) + (cse_p12_t172_g180 * cse_p15_t1_g261) + (cse_p12_t1_g234 * ((cse_p12_t182_g190 * cse_p14_t1_g254) + (cse_p12_t195_g203 * v2_len))) + (cse_p13_t1_g247 * (c1x + (cse_p12_t187_g195 * (-(v2z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x4, 1, h_z2_y4);
      double h_z2_z4 = ((cse_p12_t157_g165 * cse_p12_t158_g166 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t149_g157 * cse_p12_t171_g179) + (cse_p12_t149_g157 * cse_p14_t1_g259) + (cse_p12_t172_g180 * cse_p13_t1_g248) + (cse_p12_t172_g180 * cse_p14_t1_g257) + (cse_p13_t1_g247 * ((cse_p12_t187_g195 * v2y) + (cse_p12_t195_g203 * (-(v2x))))) + (-((cse_p12_t183_g191 * cse_p13_t1_g245 * v2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x4, 2, h_z2_z4);
      double h_x3_x3 = ((cse_p12_t159_g167 * cse_p12_t159_g167 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t150_g158 * cse_p12_t205_g213 * t2) + (cse_p12_t173_g181 * (cse_p21_t1_g265 + (cse_p11_t4_invr7 * cse_p12_t150_g158) + (cse_p11_t4_invr7 * cse_p12_t219_g227))) + (cse_p12_t1_g234 * (cse_p12_t1_g231 + (cse_p12_t184_g192 * cse_p22_t1_g267) + (cse_p12_t184_g192 * cse_p22_t1_g267) + (-(cse_p12_t66_g74)))) + (cse_p13_t1_g247 * ((cse_p12_t191_g199 * v1z) + (cse_p12_t191_g199 * v1z) + (cse_p12_t199_g207 * (-(v1y))) + (cse_p12_t199_g207 * (-(v1y))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 0, h_x3_x3);
      double h_x3_y3 = ((cse_p12_t159_g167 * cse_p15_t1_g263) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t150_g158 * cse_p12_t168_g176) + (cse_p11_t4_invr7 * cse_p12_t173_g181 * cse_p12_t220_g228) + (cse_p12_t150_g158 * cse_p12_t206_g214 * t2) + (cse_p12_t173_g181 * cse_p14_t1_g260) + (cse_p12_t1_g234 * ((cse_p12_t181_g189 * cse_p22_t1_g267) + (cse_p12_t184_g192 * cse_p22_t1_g268) + (-(cse_p12_t67_g75)))) + (cse_p13_t1_g247 * ((cse_p12_t189_g197 * (-(v1y))) + (cse_p12_t199_g207 * v1x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 1, h_x3_y3);
      double h_x3_z3 = ((cse_p12_t159_g167 * cse_p14_t1_g258) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t150_g158 * cse_p12_t167_g175) + (cse_p11_t4_invr7 * cse_p12_t173_g181 * cse_p12_t221_g229) + (cse_p12_t150_g158 * cse_p12_t207_g215 * t2) + (cse_p12_t173_g181 * cse_p14_t1_g251) + (cse_p12_t1_g234 * ((cse_p12_t180_g188 * cse_p22_t1_g267) + (cse_p12_t184_g192 * cse_p22_t1_g269) + (-(cse_p12_t68_g76)))) + (cse_p13_t1_g247 * ((cse_p12_t191_g199 * (-(v1x))) + (cse_p12_t198_g206 * v1z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 2, h_x3_z3);
      double h_x3_x4 = ((cse_p12_t159_g167 * cse_p12_t160_g168 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t150_g158 * cse_p12_t174_g182) + (cse_p12_t150_g158 * cse_p12_t210_g218 * t2) + (cse_p12_t173_g181 * cse_p13_t1_g250) + (cse_p12_t173_g181 * cse_p21_t1_g266) + (cse_p12_t185_g193 * cse_p13_t1_g245 * v2x) + (cse_p13_t1_g247 * (((-(v1y)) * (-(v2y))) + (v1z * v2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x4, 0, h_x3_x4);
      double h_x3_y4 = ((cse_p12_t155_g163 * cse_p12_t159_g167 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t150_g158 * cse_p12_t169_g177) + (cse_p12_t150_g158 * cse_p21_t1_g264) + (cse_p12_t173_g181 * cse_p13_t1_g249) + (cse_p12_t173_g181 * cse_p15_t1_g261) + (cse_p12_t1_g234 * ((cse_p12_t182_g190 * cse_p22_t1_g267) + (v1z * v2_len))) + (cse_p13_t1_g247 * (c1z + (v2x * (-(v1y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x4, 1, h_x3_y4);
      double h_x3_z4 = ((cse_p12_t157_g165 * cse_p12_t159_g167 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t150_g158 * cse_p12_t171_g179) + (cse_p12_t150_g158 * cse_p14_t1_g259) + (cse_p12_t173_g181 * cse_p13_t1_g248) + (cse_p12_t173_g181 * cse_p14_t1_g257) + (cse_p12_t1_g234 * ((cse_p12_t183_g191 * cse_p22_t1_g267) + (v2_len * (-(v1y))))) + (cse_p13_t1_g247 * ((v1z * (-(v2x))) + (-(c1y)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x4, 2, h_x3_z4);
      double h_y3_y3 = ((cse_p12_t154_g162 * cse_p15_t1_g263) + (cse_p12_t1_g236 * ((cse_p12_t151_g159 * cse_p12_t206_g214 * t2) + (cse_p12_t168_g176 * (cse_p14_t1_g260 + (cse_p11_t4_invr7 * cse_p12_t151_g159) + (cse_p11_t4_invr7 * cse_p12_t220_g228))) + (cse_p12_t1_g234 * (cse_p12_t1_g231 + (cse_p12_t181_g189 * cse_p22_t1_g268) + (cse_p12_t181_g189 * cse_p22_t1_g268) + (-(cse_p12_t69_g77)))) + (cse_p13_t1_g247 * ((cse_p12_t189_g197 * v1x) + (cse_p12_t189_g197 * v1x) + (cse_p12_t200_g208 * (-(v1z))) + (cse_p12_t200_g208 * (-(v1z))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 1, h_y3_y3);
      double h_y3_z3 = ((cse_p12_t154_g162 * cse_p14_t1_g258) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t151_g159 * cse_p12_t167_g175) + (cse_p11_t4_invr7 * cse_p12_t168_g176 * cse_p12_t221_g229) + (cse_p12_t151_g159 * cse_p12_t207_g215 * t2) + (cse_p12_t168_g176 * cse_p14_t1_g251) + (cse_p12_t1_g234 * ((cse_p12_t180_g188 * cse_p22_t1_g268) + (cse_p12_t181_g189 * cse_p22_t1_g269) + (-(cse_p12_t70_g78)))) + (cse_p13_t1_g247 * ((cse_p12_t190_g198 * (-(v1z))) + (cse_p12_t200_g208 * v1y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 2, h_y3_z3);
      double h_y3_x4 = ((cse_p12_t160_g168 * cse_p15_t1_g263) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t151_g159 * cse_p12_t174_g182) + (cse_p12_t151_g159 * cse_p12_t210_g218 * t2) + (cse_p12_t168_g176 * cse_p13_t1_g250) + (cse_p12_t168_g176 * cse_p21_t1_g266) + (cse_p12_t1_g234 * ((cse_p12_t185_g193 * cse_p22_t1_g268) + (v2_len * (-(v1z))))) + (cse_p13_t1_g247 * ((v1x * (-(v2y))) + (-(c1z)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x4, 0, h_y3_x4);
      double h_y3_y4 = ((cse_p12_t155_g163 * cse_p15_t1_g263) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t151_g159 * cse_p12_t169_g177) + (cse_p12_t151_g159 * cse_p21_t1_g264) + (cse_p12_t168_g176 * cse_p13_t1_g249) + (cse_p12_t168_g176 * cse_p15_t1_g261) + (cse_p12_t182_g190 * cse_p13_t1_g245 * v2y) + (cse_p13_t1_g247 * (((-(v1z)) * (-(v2z))) + (v1x * v2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x4, 1, h_y3_y4);
      double h_y3_z4 = ((cse_p12_t157_g165 * cse_p15_t1_g263) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t151_g159 * cse_p12_t171_g179) + (cse_p12_t151_g159 * cse_p14_t1_g259) + (cse_p12_t168_g176 * cse_p13_t1_g248) + (cse_p12_t168_g176 * cse_p14_t1_g257) + (cse_p12_t1_g234 * ((cse_p12_t183_g191 * cse_p22_t1_g268) + (v1x * v2_len))) + (cse_p13_t1_g247 * (c1x + (v2y * (-(v1z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x4, 2, h_y3_z4);
      double h_z3_z3 = ((cse_p12_t153_g161 * cse_p14_t1_g258) + (cse_p12_t1_g236 * ((cse_p12_t152_g160 * cse_p12_t207_g215 * t2) + (cse_p12_t167_g175 * (cse_p14_t1_g251 + (cse_p11_t4_invr7 * cse_p12_t152_g160) + (cse_p11_t4_invr7 * cse_p12_t221_g229))) + (cse_p12_t1_g234 * (cse_p12_t1_g231 + (cse_p12_t180_g188 * cse_p22_t1_g269) + (cse_p12_t180_g188 * cse_p22_t1_g269) + (-(cse_p12_t71_g79)))) + (cse_p13_t1_g247 * ((cse_p12_t190_g198 * v1y) + (cse_p12_t190_g198 * v1y) + (cse_p12_t198_g206 * (-(v1x))) + (cse_p12_t198_g206 * (-(v1x))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x3, 2, h_z3_z3);
      double h_z3_x4 = ((cse_p12_t160_g168 * cse_p14_t1_g258) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t152_g160 * cse_p12_t174_g182) + (cse_p12_t152_g160 * cse_p12_t210_g218 * t2) + (cse_p12_t167_g175 * cse_p13_t1_g250) + (cse_p12_t167_g175 * cse_p21_t1_g266) + (cse_p12_t1_g234 * ((cse_p12_t185_g193 * cse_p22_t1_g269) + (v1y * v2_len))) + (cse_p13_t1_g247 * (c1y + (v2z * (-(v1x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x4, 0, h_z3_x4);
      double h_z3_y4 = ((cse_p12_t155_g163 * cse_p14_t1_g258) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t152_g160 * cse_p12_t169_g177) + (cse_p12_t152_g160 * cse_p21_t1_g264) + (cse_p12_t167_g175 * cse_p13_t1_g249) + (cse_p12_t167_g175 * cse_p15_t1_g261) + (cse_p12_t1_g234 * ((cse_p12_t182_g190 * cse_p22_t1_g269) + (v2_len * (-(v1x))))) + (cse_p13_t1_g247 * ((v1y * (-(v2z))) + (-(c1x)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x4, 1, h_z3_y4);
      double h_z3_z4 = ((cse_p12_t157_g165 * cse_p14_t1_g258) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t152_g160 * cse_p12_t171_g179) + (cse_p12_t152_g160 * cse_p14_t1_g259) + (cse_p12_t167_g175 * cse_p13_t1_g248) + (cse_p12_t167_g175 * cse_p14_t1_g257) + (cse_p12_t183_g191 * cse_p13_t1_g245 * v2z) + (cse_p13_t1_g247 * (((-(v1x)) * (-(v2x))) + (v1y * v2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x4, 2, h_z3_z4);
      double h_x4_x4 = ((cse_p12_t160_g168 * cse_p12_t160_g168 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t174_g182 * (cse_p13_t1_g250 + cse_p21_t1_g266 + (cse_p11_t4_invr7 * cse_p12_t185_g193 * v2_len))) + (cse_p12_t185_g193 * cse_p12_t210_g218 * t2 * v2_len))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 0, term.i3x4, 0, h_x4_x4);
      double h_x4_y4 = ((cse_p12_t155_g163 * cse_p12_t160_g168 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t169_g177 * cse_p12_t185_g193 * v2_len) + (cse_p12_t174_g182 * cse_p13_t1_g249) + (cse_p12_t174_g182 * cse_p15_t1_g261) + (cse_p12_t185_g193 * cse_p21_t1_g264 * v2_len))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 0, term.i3x4, 1, h_x4_y4);
      double h_x4_z4 = ((cse_p12_t157_g165 * cse_p12_t160_g168 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t171_g179 * cse_p12_t185_g193 * v2_len) + (cse_p12_t174_g182 * cse_p13_t1_g248) + (cse_p12_t174_g182 * cse_p14_t1_g257) + (cse_p12_t185_g193 * cse_p14_t1_g259 * v2_len))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 0, term.i3x4, 2, h_x4_z4);
      double h_y4_y4 = ((cse_p12_t155_g163 * cse_p12_t155_g163 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t169_g177 * (cse_p13_t1_g249 + cse_p15_t1_g261 + (cse_p11_t4_invr7 * cse_p12_t182_g190 * v2_len))) + (cse_p12_t182_g190 * cse_p21_t1_g264 * v2_len))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 1, term.i3x4, 1, h_y4_y4);
      double h_y4_z4 = ((cse_p12_t155_g163 * cse_p12_t157_g165 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p11_t4_invr7 * cse_p12_t171_g179 * cse_p12_t182_g190 * v2_len) + (cse_p12_t169_g177 * cse_p13_t1_g248) + (cse_p12_t169_g177 * cse_p14_t1_g257) + (cse_p12_t182_g190 * cse_p14_t1_g259 * v2_len))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x4, 1, term.i3x4, 2, h_y4_z4);
      double h_z4_z4 = ((cse_p12_t157_g165 * cse_p12_t157_g165 * cse_p12_t1_g235) + (cse_p12_t1_g236 * ((cse_p12_t171_g179 * (cse_p13_t1_g248 + cse_p14_t1_g257 + (cse_p11_t4_invr7 * cse_p12_t183_g191 * v2_len))) + (cse_p12_t183_g191 * cse_p14_t1_g259 * v2_len))));
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
