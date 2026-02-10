template <typename HESSIAN>
struct Dihedral_Fast {
  static constexpr size_t PositionSize = 12;
  static std::string description() { return "mathkernel-dihedral_fast"; };
double energy(double v, double n, double sinphase, double cosphase, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
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
      double nphi = (n * phi);
      double sin_nphi = sin(nphi);
      double cos_nphi = cos(nphi);
      double cos_angle = ((cosphase * cos_nphi) + (sinphase * sin_nphi));
      double sin_angle = ((cosphase * sin_nphi) + (-((cos_nphi * sinphase))));
      double energy = (v * (1.0000000000000000     + cos_angle));
      energy_added += energy;
      *energy_accumulate += energy;
    }
  }
  return energy_added;
}
void energy_fd(double v, double n, double sinphase, double cosphase, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, energy_accumulate, force, hessian, dvec, hdvec);
}

double gradient(double v, double n, double sinphase, double cosphase, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
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
      double v1x = (x2 + (-(x1)));
      double v1y = (y2 + (-(y1)));
      double v1z = (z2 + (-(z1)));
      double v2x = (x3 + (-(x2)));
      double cse_p12_t6_g27976 = (v1z * v2x);
      double v2y = (y3 + (-(y2)));
      double cse_p12_t4_g27974 = (v1x * v2y);
      double v2z = (z3 + (-(z2)));
      double cse_p12_t5_g27975 = (v1y * v2z);
      double v3x = (x4 + (-(x3)));
      double cse_p12_t8_g27978 = ((-(v2x)) + (-(v3x)));
      double cse_p12_t11_g27981 = (v2x + v3x);
      double v3y = (y4 + (-(y3)));
      double cse_p12_t9_g27979 = ((-(v2y)) + (-(v3y)));
      double cse_p12_t12_g27982 = (v2y + v3y);
      double v3z = (z4 + (-(z3)));
      double cse_p12_t10_g27980 = ((-(v2z)) + (-(v3z)));
      double cse_p12_t13_g27983 = (v2z + v3z);
      double c1x = (cse_p12_t5_g27975 + (-((v1z * v2y))));
      double c1y = (cse_p12_t6_g27976 + (-((v1x * v2z))));
      double c1z = (cse_p12_t4_g27974 + (-((v1y * v2x))));
      double c2x = ((v2y * v3z) + (-((v2z * v3y))));
      double cse_p12_t1_g27971 = (c2x * v1x);
      double c2y = ((v2z * v3x) + (-((v2x * v3z))));
      double cse_p12_t2_g27972 = (c2y * v1y);
      double c2z = ((v2x * v3y) + (-((v2y * v3x))));
      double cse_p12_t3_g27973 = (c2z * v1z);
      double cse_p12_t7_g27977 = (cse_p12_t1_g27971 + cse_p12_t2_g27972 + cse_p12_t3_g27973);
      double v2_sq = ((v2x * v2x) + (v2y * v2y) + (v2z * v2z));
      double cse_p1_t1_sqrt27968 = sqrt(v2_sq);
      double cse_p1_t2_invsqrt27969 = (1.0 / (cse_p1_t1_sqrt27968));
      double t1 = (cse_p12_t7_g27977 * cse_p1_t1_sqrt27968);
      double t2 = ((c1x * c2x) + (c1y * c2y) + (c1z * c2z));
      double phi = atan2(t1, t2);
      double nphi = (n * phi);
      double cse_p12_t14_g27984 = (-(sin(nphi)));
      double sin_nphi = sin(nphi);
      double cos_nphi = cos(nphi);
      double cos_angle = ((cosphase * cos_nphi) + (sinphase * sin_nphi));
      double sin_angle = ((cosphase * sin_nphi) + (-((cos_nphi * sinphase))));
      double energy = (v * (1.0000000000000000     + cos_angle));
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p1_t3_invr27970 = (1.0 / ((((t1) * (t1)) + ((t2) * (t2)))));
      double cse_p12_t1_g27985 = (cse_p1_t3_invr27970 * (-(t1)));
      double cse_p12_t1_g27986 = (cse_p1_t3_invr27970 * t2);
      double cse_p12_t1_g27987 = (cos_nphi * n * sinphase);
      double cse_p12_t1_g27988 = (cosphase * cse_p12_t14_g27984 * n);
      double cse_p12_t1_g27989 = (cse_p12_t1_g27986 * cse_p1_t1_sqrt27968);
      double g_x1 = (v * ((cse_p12_t1_g27987 * ((cse_p12_t1_g27985 * ((c2y * v2z) + (c2z * (-(v2y))))) + (cse_p12_t1_g27989 * (-(c2x))))) + (cse_p12_t1_g27988 * ((cse_p12_t1_g27985 * ((c2y * v2z) + (c2z * (-(v2y))))) + (cse_p12_t1_g27989 * (-(c2x)))))));
      KernelGradientAcc(i3x1, 0, g_x1);
      double g_y1 = (v * ((cse_p12_t1_g27987 * ((cse_p12_t1_g27985 * ((c2x * (-(v2z))) + (c2z * v2x))) + (cse_p12_t1_g27989 * (-(c2y))))) + (cse_p12_t1_g27988 * ((cse_p12_t1_g27985 * ((c2x * (-(v2z))) + (c2z * v2x))) + (cse_p12_t1_g27989 * (-(c2y)))))));
      KernelGradientAcc(i3x1, 1, g_y1);
      double g_z1 = (v * ((cse_p12_t1_g27987 * ((cse_p12_t1_g27985 * ((c2x * v2y) + (c2y * (-(v2x))))) + (cse_p12_t1_g27989 * (-(c2z))))) + (cse_p12_t1_g27988 * ((cse_p12_t1_g27985 * ((c2x * v2y) + (c2y * (-(v2x))))) + (cse_p12_t1_g27989 * (-(c2z)))))));
      KernelGradientAcc(i3x1, 2, g_z1);
      double cse_p12_t1_g27990 = (cse_p12_t7_g27977 * cse_p1_t2_invsqrt27969);
      double g_x2 = (v * ((cse_p12_t1_g27987 * ((cse_p12_t1_g27985 * ((c1y * v3z) + (c1z * (-(v3y))) + (c2y * ((-(v1z)) + (-(v2z)))) + (c2z * (v1y + v2y)))) + (cse_p12_t1_g27986 * ((cse_p12_t1_g27990 * (-(v2x))) + (cse_p1_t1_sqrt27968 * (c2x + (v1y * v3z) + (v1z * (-(v3y))))))))) + (cse_p12_t1_g27988 * ((cse_p12_t1_g27985 * ((c1y * v3z) + (c1z * (-(v3y))) + (c2y * ((-(v1z)) + (-(v2z)))) + (c2z * (v1y + v2y)))) + (cse_p12_t1_g27986 * ((cse_p12_t1_g27990 * (-(v2x))) + (cse_p1_t1_sqrt27968 * (c2x + (v1y * v3z) + (v1z * (-(v3y)))))))))));
      KernelGradientAcc(i3x2, 0, g_x2);
      double g_y2 = (v * ((cse_p12_t1_g27987 * ((cse_p12_t1_g27985 * ((c1x * (-(v3z))) + (c1z * v3x) + (c2x * (v1z + v2z)) + (c2z * ((-(v1x)) + (-(v2x)))))) + (cse_p12_t1_g27986 * ((cse_p12_t1_g27990 * (-(v2y))) + (cse_p1_t1_sqrt27968 * (c2y + (v1x * (-(v3z))) + (v1z * v3x))))))) + (cse_p12_t1_g27988 * ((cse_p12_t1_g27985 * ((c1x * (-(v3z))) + (c1z * v3x) + (c2x * (v1z + v2z)) + (c2z * ((-(v1x)) + (-(v2x)))))) + (cse_p12_t1_g27986 * ((cse_p12_t1_g27990 * (-(v2y))) + (cse_p1_t1_sqrt27968 * (c2y + (v1x * (-(v3z))) + (v1z * v3x)))))))));
      KernelGradientAcc(i3x2, 1, g_y2);
      double g_z2 = (v * ((cse_p12_t1_g27987 * ((cse_p12_t1_g27985 * ((c1x * v3y) + (c1y * (-(v3x))) + (c2x * ((-(v1y)) + (-(v2y)))) + (c2y * (v1x + v2x)))) + (cse_p12_t1_g27986 * ((cse_p12_t1_g27990 * (-(v2z))) + (cse_p1_t1_sqrt27968 * (c2z + (v1x * v3y) + (v1y * (-(v3x))))))))) + (cse_p12_t1_g27988 * ((cse_p12_t1_g27985 * ((c1x * v3y) + (c1y * (-(v3x))) + (c2x * ((-(v1y)) + (-(v2y)))) + (c2y * (v1x + v2x)))) + (cse_p12_t1_g27986 * ((cse_p12_t1_g27990 * (-(v2z))) + (cse_p1_t1_sqrt27968 * (c2z + (v1x * v3y) + (v1y * (-(v3x)))))))))));
      KernelGradientAcc(i3x2, 2, g_z2);
      double g_x3 = (v * ((cse_p12_t1_g27987 * ((cse_p12_t1_g27985 * ((c1y * cse_p12_t10_g27980) + (c1z * cse_p12_t12_g27982) + (c2y * v1z) + (c2z * (-(v1y))))) + (cse_p12_t1_g27986 * ((cse_p12_t1_g27990 * v2x) + (cse_p1_t1_sqrt27968 * ((cse_p12_t10_g27980 * v1y) + (cse_p12_t12_g27982 * v1z))))))) + (cse_p12_t1_g27988 * ((cse_p12_t1_g27985 * ((c1y * cse_p12_t10_g27980) + (c1z * cse_p12_t12_g27982) + (c2y * v1z) + (c2z * (-(v1y))))) + (cse_p12_t1_g27986 * ((cse_p12_t1_g27990 * v2x) + (cse_p1_t1_sqrt27968 * ((cse_p12_t10_g27980 * v1y) + (cse_p12_t12_g27982 * v1z)))))))));
      KernelGradientAcc(i3x3, 0, g_x3);
      double g_y3 = (v * ((cse_p12_t1_g27987 * ((cse_p12_t1_g27985 * ((c1x * cse_p12_t13_g27983) + (c1z * cse_p12_t8_g27978) + (c2x * (-(v1z))) + (c2z * v1x))) + (cse_p12_t1_g27986 * ((cse_p12_t1_g27990 * v2y) + (cse_p1_t1_sqrt27968 * ((cse_p12_t13_g27983 * v1x) + (cse_p12_t8_g27978 * v1z))))))) + (cse_p12_t1_g27988 * ((cse_p12_t1_g27985 * ((c1x * cse_p12_t13_g27983) + (c1z * cse_p12_t8_g27978) + (c2x * (-(v1z))) + (c2z * v1x))) + (cse_p12_t1_g27986 * ((cse_p12_t1_g27990 * v2y) + (cse_p1_t1_sqrt27968 * ((cse_p12_t13_g27983 * v1x) + (cse_p12_t8_g27978 * v1z)))))))));
      KernelGradientAcc(i3x3, 1, g_y3);
      double g_z3 = (v * ((cse_p12_t1_g27987 * ((cse_p12_t1_g27985 * ((c1x * cse_p12_t9_g27979) + (c1y * cse_p12_t11_g27981) + (c2x * v1y) + (c2y * (-(v1x))))) + (cse_p12_t1_g27986 * ((cse_p12_t1_g27990 * v2z) + (cse_p1_t1_sqrt27968 * ((cse_p12_t11_g27981 * v1y) + (cse_p12_t9_g27979 * v1x))))))) + (cse_p12_t1_g27988 * ((cse_p12_t1_g27985 * ((c1x * cse_p12_t9_g27979) + (c1y * cse_p12_t11_g27981) + (c2x * v1y) + (c2y * (-(v1x))))) + (cse_p12_t1_g27986 * ((cse_p12_t1_g27990 * v2z) + (cse_p1_t1_sqrt27968 * ((cse_p12_t11_g27981 * v1y) + (cse_p12_t9_g27979 * v1x)))))))));
      KernelGradientAcc(i3x3, 2, g_z3);
      double g_x4 = (v * ((cse_p12_t1_g27987 * ((cse_p12_t1_g27985 * ((c1y * v2z) + (c1z * (-(v2y))))) + (cse_p12_t1_g27989 * (cse_p12_t5_g27975 + (v1z * (-(v2y))))))) + (cse_p12_t1_g27988 * ((cse_p12_t1_g27985 * ((c1y * v2z) + (c1z * (-(v2y))))) + (cse_p12_t1_g27989 * (cse_p12_t5_g27975 + (v1z * (-(v2y)))))))));
      KernelGradientAcc(i3x4, 0, g_x4);
      double g_y4 = (v * ((cse_p12_t1_g27987 * ((cse_p12_t1_g27985 * ((c1x * (-(v2z))) + (c1z * v2x))) + (cse_p12_t1_g27989 * (cse_p12_t6_g27976 + (v1x * (-(v2z))))))) + (cse_p12_t1_g27988 * ((cse_p12_t1_g27985 * ((c1x * (-(v2z))) + (c1z * v2x))) + (cse_p12_t1_g27989 * (cse_p12_t6_g27976 + (v1x * (-(v2z)))))))));
      KernelGradientAcc(i3x4, 1, g_y4);
      double g_z4 = (v * ((cse_p12_t1_g27987 * ((cse_p12_t1_g27985 * ((c1x * v2y) + (c1y * (-(v2x))))) + (cse_p12_t1_g27989 * (cse_p12_t4_g27974 + (v1y * (-(v2x))))))) + (cse_p12_t1_g27988 * ((cse_p12_t1_g27985 * ((c1x * v2y) + (c1y * (-(v2x))))) + (cse_p12_t1_g27989 * (cse_p12_t4_g27974 + (v1y * (-(v2x)))))))));
      KernelGradientAcc(i3x4, 2, g_z4);
    }
  }
  return energy_added;
}
void gradient_fd(double v, double n, double sinphase, double cosphase, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 12;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  double e0 = 0.0;
  energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
  {
    double saved = position[i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 0] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 0, d);
  }
  {
    double saved = position[i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 1] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 1, d);
  }
  {
    double saved = position[i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 2] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 2, d);
  }
  {
    double saved = position[i3x4 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 0] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 0, d);
  }
  {
    double saved = position[i3x4 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 1] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 1, d);
  }
  {
    double saved = position[i3x4 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 2] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 2, d);
  }
}

double hessian(double v, double n, double sinphase, double cosphase, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
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
      double v1x = (x2 + (-(x1)));
      double v1y = (y2 + (-(y1)));
      double v1z = (z2 + (-(z1)));
      double v2x = (x3 + (-(x2)));
      double cse_p12_t118_g28116 = (v1y * (-(v2x)));
      double cse_p12_t126_g28124 = (v1z * v2x);
      double cse_p12_t186_g28184 = ((-(v1x)) + (-(v2x)));
      double cse_p12_t195_g28193 = (v1x + v2x);
      double v2y = (y3 + (-(y2)));
      double cse_p12_t114_g28112 = (v1x * v2y);
      double cse_p12_t124_g28122 = (v1z * (-(v2y)));
      double cse_p12_t183_g28181 = (cse_p12_t114_g28112 + cse_p12_t118_g28116);
      double cse_p12_t187_g28185 = ((-(v1y)) + (-(v2y)));
      double cse_p12_t196_g28194 = (v1y + v2y);
      double v2z = (z3 + (-(z2)));
      double cse_p12_t112_g28110 = (v1x * (-(v2z)));
      double cse_p12_t120_g28118 = (v1y * v2z);
      double cse_p12_t182_g28180 = (cse_p12_t112_g28110 + cse_p12_t126_g28124);
      double cse_p12_t185_g28183 = (cse_p12_t120_g28118 + cse_p12_t124_g28122);
      double cse_p12_t188_g28186 = ((-(v1z)) + (-(v2z)));
      double cse_p12_t197_g28195 = (v1z + v2z);
      double v3x = (x4 + (-(x3)));
      double cse_p12_t117_g28115 = (v1y * (v2x + v3x));
      double cse_p12_t119_g28117 = (v1y * (-(v3x)));
      double cse_p12_t122_g28120 = (v1z * ((-(v2x)) + (-(v3x))));
      double cse_p12_t127_g28125 = (v1z * v3x);
      double cse_p12_t189_g28187 = ((-(v2x)) + (-(v3x)));
      double cse_p12_t198_g28196 = (v2x + v3x);
      double v3y = (y4 + (-(y3)));
      double cse_p12_t110_g28108 = (v1x * ((-(v2y)) + (-(v3y))));
      double cse_p12_t115_g28113 = (v1x * v3y);
      double cse_p12_t123_g28121 = (v1z * (v2y + v3y));
      double cse_p12_t125_g28123 = (v1z * (-(v3y)));
      double cse_p12_t180_g28178 = (cse_p12_t110_g28108 + cse_p12_t117_g28115);
      double cse_p12_t190_g28188 = ((-(v2y)) + (-(v3y)));
      double cse_p12_t199_g28197 = (v2y + v3y);
      double v3z = (z4 + (-(z3)));
      double cse_p12_t111_g28109 = (v1x * (v2z + v3z));
      double cse_p12_t113_g28111 = (v1x * (-(v3z)));
      double cse_p12_t116_g28114 = (v1y * ((-(v2z)) + (-(v3z))));
      double cse_p12_t121_g28119 = (v1y * v3z);
      double cse_p12_t181_g28179 = (cse_p12_t111_g28109 + cse_p12_t122_g28120);
      double cse_p12_t184_g28182 = (cse_p12_t116_g28114 + cse_p12_t123_g28121);
      double cse_p12_t191_g28189 = ((-(v2z)) + (-(v3z)));
      double cse_p12_t200_g28198 = (v2z + v3z);
      double c1x = (cse_p12_t120_g28118 + (-((v1z * v2y))));
      double cse_p12_t23_g28021 = (c1x * cse_p12_t190_g28188);
      double cse_p12_t24_g28022 = (c1x * cse_p12_t200_g28198);
      double cse_p12_t25_g28023 = (c1x * (-(v2z)));
      double cse_p12_t26_g28024 = (c1x * (-(v3z)));
      double cse_p12_t27_g28025 = (c1x * v2y);
      double cse_p12_t28_g28026 = (c1x * v3y);
      double c1y = (cse_p12_t126_g28124 + (-((v1x * v2z))));
      double cse_p12_t29_g28027 = (c1y * cse_p12_t191_g28189);
      double cse_p12_t30_g28028 = (c1y * cse_p12_t198_g28196);
      double cse_p12_t31_g28029 = (c1y * (-(v2x)));
      double cse_p12_t32_g28030 = (c1y * (-(v3x)));
      double cse_p12_t33_g28031 = (c1y * v2z);
      double cse_p12_t34_g28032 = (c1y * v3z);
      double cse_p12_t171_g28169 = (cse_p12_t27_g28025 + cse_p12_t31_g28029);
      double c1z = (cse_p12_t114_g28112 + (-((v1y * v2x))));
      double cse_p12_t35_g28033 = (c1z * cse_p12_t189_g28187);
      double cse_p12_t36_g28034 = (c1z * cse_p12_t199_g28197);
      double cse_p12_t37_g28035 = (c1z * (-(v2y)));
      double cse_p12_t38_g28036 = (c1z * (-(v3y)));
      double cse_p12_t39_g28037 = (c1z * v2x);
      double cse_p12_t40_g28038 = (c1z * v3x);
      double cse_p12_t169_g28167 = (cse_p12_t25_g28023 + cse_p12_t39_g28037);
      double cse_p12_t174_g28172 = (cse_p12_t33_g28031 + cse_p12_t37_g28035);
      double c2x = ((v2y * v3z) + (-((v2z * v3y))));
      double cse_p12_t41_g28039 = (c2x * cse_p12_t187_g28185);
      double cse_p12_t42_g28040 = (c2x * cse_p12_t197_g28195);
      double cse_p12_t43_g28041 = (c2x * (-(v1z)));
      double cse_p12_t44_g28042 = (c2x * (-(v2z)));
      double cse_p12_t45_g28043 = (c2x * v1x);
      double cse_p12_t46_g28044 = (c2x * v1y);
      double cse_p12_t47_g28045 = (c2x * v2y);
      double cse_p12_t192_g28190 = (c2x + cse_p12_t121_g28119 + cse_p12_t125_g28123);
      double c2y = ((v2z * v3x) + (-((v2x * v3z))));
      double cse_p12_t48_g28046 = (c2y * cse_p12_t188_g28186);
      double cse_p12_t49_g28047 = (c2y * cse_p12_t195_g28193);
      double cse_p12_t50_g28048 = (c2y * (-(v1x)));
      double cse_p12_t51_g28049 = (c2y * (-(v2x)));
      double cse_p12_t52_g28050 = (c2y * v1y);
      double cse_p12_t53_g28051 = (c2y * v1z);
      double cse_p12_t54_g28052 = (c2y * v2z);
      double cse_p12_t167_g28165 = (cse_p12_t23_g28021 + cse_p12_t30_g28028 + cse_p12_t46_g28044 + cse_p12_t50_g28048);
      double cse_p12_t172_g28170 = (cse_p12_t28_g28026 + cse_p12_t32_g28030 + cse_p12_t41_g28039 + cse_p12_t49_g28047);
      double cse_p12_t178_g28176 = (cse_p12_t47_g28045 + cse_p12_t51_g28049);
      double cse_p12_t193_g28191 = (c2y + cse_p12_t113_g28111 + cse_p12_t127_g28125);
      double c2z = ((v2x * v3y) + (-((v2y * v3x))));
      double cse_p12_t55_g28053 = (c2z * cse_p12_t186_g28184);
      double cse_p12_t56_g28054 = (c2z * cse_p12_t196_g28194);
      double cse_p12_t57_g28055 = (c2z * (-(v1y)));
      double cse_p12_t58_g28056 = (c2z * (-(v2y)));
      double cse_p12_t59_g28057 = (c2z * v1x);
      double cse_p12_t60_g28058 = (c2z * v1z);
      double cse_p12_t61_g28059 = (c2z * v2x);
      double cse_p12_t168_g28166 = (cse_p12_t24_g28022 + cse_p12_t35_g28033 + cse_p12_t43_g28041 + cse_p12_t59_g28057);
      double cse_p12_t170_g28168 = (cse_p12_t26_g28024 + cse_p12_t40_g28038 + cse_p12_t42_g28040 + cse_p12_t55_g28053);
      double cse_p12_t173_g28171 = (cse_p12_t29_g28027 + cse_p12_t36_g28034 + cse_p12_t53_g28051 + cse_p12_t57_g28055);
      double cse_p12_t175_g28173 = (cse_p12_t34_g28032 + cse_p12_t38_g28036 + cse_p12_t48_g28046 + cse_p12_t56_g28054);
      double cse_p12_t176_g28174 = (cse_p12_t44_g28042 + cse_p12_t61_g28059);
      double cse_p12_t177_g28175 = (cse_p12_t45_g28043 + cse_p12_t52_g28050 + cse_p12_t60_g28058);
      double cse_p12_t179_g28177 = (cse_p12_t54_g28052 + cse_p12_t58_g28056);
      double cse_p12_t194_g28192 = (c2z + cse_p12_t115_g28113 + cse_p12_t119_g28117);
      double v2_sq = ((v2x * v2x) + (v2y * v2y) + (v2z * v2z));
      double cse_p11_t1_sqrt27994 = sqrt(v2_sq);
      double cse_p11_t2_invsqrt27995 = (1.0 / (cse_p11_t1_sqrt27994));
      double cse_p12_t1_g28222 = (cse_p11_t2_invsqrt27995 * cse_p12_t177_g28175);
      double cse_p12_t62_g28060 = (cse_p12_t1_g28222 * (-(v2x)));
      double cse_p12_t63_g28061 = (cse_p12_t1_g28222 * (-(v2y)));
      double cse_p12_t64_g28062 = (cse_p12_t1_g28222 * (-(v2z)));
      double cse_p12_t72_g28070 = (cse_p12_t1_g28222 * v2x);
      double cse_p12_t73_g28071 = (cse_p12_t1_g28222 * v2y);
      double cse_p12_t74_g28072 = (cse_p12_t1_g28222 * v2z);
      double cse_p12_t201_g28199 = (-(cse_p12_t1_g28222));
      double cse_p11_t3_invr27996 = (cse_p11_t2_invsqrt27995 * cse_p11_t2_invsqrt27995);
      double cse_p12_t66_g28064 = (cse_p11_t3_invr27996 * cse_p12_t1_g28222 * v2x * v2x);
      double cse_p12_t67_g28065 = (cse_p11_t3_invr27996 * cse_p12_t1_g28222 * v2x * v2y);
      double cse_p12_t68_g28066 = (cse_p11_t3_invr27996 * cse_p12_t1_g28222 * v2x * v2z);
      double cse_p12_t69_g28067 = (cse_p11_t3_invr27996 * cse_p12_t1_g28222 * v2y * v2y);
      double cse_p12_t70_g28068 = (cse_p11_t3_invr27996 * cse_p12_t1_g28222 * v2y * v2z);
      double cse_p12_t71_g28069 = (cse_p11_t3_invr27996 * cse_p12_t1_g28222 * v2z * v2z);
      double cse_p12_t1_g28229 = (cse_p11_t1_sqrt27994 * cse_p12_t180_g28178);
      double cse_p12_t1_g28230 = (cse_p11_t1_sqrt27994 * cse_p12_t181_g28179);
      double cse_p12_t130_g28128 = (cse_p11_t1_sqrt27994 * cse_p12_t182_g28180);
      double cse_p12_t131_g28129 = (cse_p11_t1_sqrt27994 * cse_p12_t183_g28181);
      double cse_p12_t1_g28231 = (cse_p11_t1_sqrt27994 * cse_p12_t184_g28182);
      double cse_p12_t133_g28131 = (cse_p11_t1_sqrt27994 * cse_p12_t185_g28183);
      double cse_p12_t134_g28132 = (cse_p11_t1_sqrt27994 * cse_p12_t192_g28190);
      double cse_p12_t135_g28133 = (cse_p11_t1_sqrt27994 * cse_p12_t193_g28191);
      double cse_p12_t136_g28134 = (cse_p11_t1_sqrt27994 * cse_p12_t194_g28192);
      double cse_p12_t137_g28135 = (cse_p11_t1_sqrt27994 * (-(c2z)));
      double cse_p12_t147_g28145 = (cse_p12_t134_g28132 + cse_p12_t62_g28060);
      double cse_p12_t148_g28146 = (cse_p12_t135_g28133 + cse_p12_t63_g28061);
      double cse_p12_t149_g28147 = (cse_p12_t136_g28134 + cse_p12_t64_g28062);
      double cse_p12_t150_g28148 = (cse_p12_t1_g28231 + cse_p12_t72_g28070);
      double cse_p12_t151_g28149 = (cse_p12_t1_g28230 + cse_p12_t73_g28071);
      double cse_p12_t152_g28150 = (cse_p12_t1_g28229 + cse_p12_t74_g28072);
      double cse_p12_t213_g28211 = (-(cse_p12_t130_g28128));
      double cse_p12_t214_g28212 = (-(cse_p12_t131_g28129));
      double cse_p12_t215_g28213 = (-(cse_p12_t133_g28131));
      double cse_p12_t216_g28214 = (-(cse_p12_t147_g28145));
      double cse_p12_t217_g28215 = (-(cse_p12_t148_g28146));
      double cse_p12_t218_g28216 = (-(cse_p12_t149_g28147));
      double cse_p12_t219_g28217 = (-(cse_p12_t150_g28148));
      double cse_p12_t220_g28218 = (-(cse_p12_t151_g28149));
      double cse_p12_t221_g28219 = (-(cse_p12_t152_g28150));
      double cse_p12_t222_g28220 = (-(cse_p12_t137_g28135));
      double t1 = (cse_p11_t1_sqrt27994 * cse_p12_t177_g28175);
      double cse_p13_t1_g28238 = (2.0000000000000000     * t1);
      double cse_p12_t1_g27999 = (cse_p12_t147_g28145 * cse_p13_t1_g28238);
      double cse_p12_t2_g28000 = (cse_p12_t148_g28146 * cse_p13_t1_g28238);
      double cse_p12_t3_g28001 = (cse_p12_t149_g28147 * cse_p13_t1_g28238);
      double cse_p12_t4_g28002 = (cse_p12_t150_g28148 * cse_p13_t1_g28238);
      double cse_p12_t5_g28003 = (cse_p12_t151_g28149 * cse_p13_t1_g28238);
      double cse_p12_t6_g28004 = (cse_p12_t152_g28150 * cse_p13_t1_g28238);
      double cse_p12_t7_g28005 = (cse_p11_t1_sqrt27994 * cse_p12_t182_g28180 * cse_p13_t1_g28238);
      double cse_p12_t8_g28006 = (cse_p11_t1_sqrt27994 * cse_p12_t183_g28181 * cse_p13_t1_g28238);
      double cse_p12_t9_g28007 = (cse_p11_t1_sqrt27994 * cse_p12_t185_g28183 * cse_p13_t1_g28238);
      double cse_p12_t10_g28008 = (cse_p11_t1_sqrt27994 * cse_p13_t1_g28238 * (-(c2y)));
      double cse_p12_t11_g28009 = (cse_p11_t1_sqrt27994 * cse_p13_t1_g28238 * (-(c2z)));
      double t2 = ((c1x * c2x) + (c1y * c2y) + (c1z * c2z));
      double cse_p13_t1_g28239 = (2.0000000000000000     * t2);
      double cse_p12_t12_g28010 = (cse_p12_t167_g28165 * cse_p13_t1_g28239);
      double cse_p12_t13_g28011 = (cse_p12_t168_g28166 * cse_p13_t1_g28239);
      double cse_p12_t14_g28012 = (cse_p12_t169_g28167 * cse_p13_t1_g28239);
      double cse_p12_t15_g28013 = (cse_p12_t170_g28168 * cse_p13_t1_g28239);
      double cse_p12_t16_g28014 = (cse_p12_t171_g28169 * cse_p13_t1_g28239);
      double cse_p12_t17_g28015 = (cse_p12_t172_g28170 * cse_p13_t1_g28239);
      double cse_p12_t18_g28016 = (cse_p12_t173_g28171 * cse_p13_t1_g28239);
      double cse_p12_t19_g28017 = (cse_p12_t174_g28172 * cse_p13_t1_g28239);
      double cse_p12_t20_g28018 = (cse_p12_t175_g28173 * cse_p13_t1_g28239);
      double cse_p12_t21_g28019 = (cse_p12_t176_g28174 * cse_p13_t1_g28239);
      double cse_p12_t22_g28020 = (cse_p12_t178_g28176 * cse_p13_t1_g28239);
      double cse_p12_t138_g28136 = (cse_p12_t1_g27999 + cse_p12_t20_g28018);
      double cse_p12_t139_g28137 = (cse_p12_t15_g28013 + cse_p12_t2_g28000);
      double cse_p12_t140_g28138 = (cse_p12_t17_g28015 + cse_p12_t3_g28001);
      double cse_p12_t141_g28139 = (cse_p12_t18_g28016 + cse_p12_t4_g28002);
      double cse_p12_t142_g28140 = (cse_p12_t13_g28011 + cse_p12_t5_g28003);
      double cse_p12_t143_g28141 = (cse_p12_t12_g28010 + cse_p12_t6_g28004);
      double cse_p12_t144_g28142 = (cse_p12_t14_g28012 + cse_p12_t7_g28005);
      double cse_p12_t145_g28143 = (cse_p12_t16_g28014 + cse_p12_t8_g28006);
      double cse_p12_t146_g28144 = (cse_p12_t19_g28017 + cse_p12_t9_g28007);
      double cse_p12_t165_g28163 = (cse_p12_t10_g28008 + cse_p12_t21_g28019);
      double cse_p12_t166_g28164 = (cse_p12_t11_g28009 + cse_p12_t22_g28020);
      double phi = atan2(t1, t2);
      double nphi = (n * phi);
      double cse_p12_t223_g28221 = (-(sin(nphi)));
      double sin_nphi = sin(nphi);
      double cos_nphi = cos(nphi);
      double cos_angle = ((cosphase * cos_nphi) + (sinphase * sin_nphi));
      double sin_angle = ((cosphase * sin_nphi) + (-((cos_nphi * sinphase))));
      double energy = (v * (1.0000000000000000     + cos_angle));
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t4_invr27997 = (1.0 / ((((t1) * (t1)) + ((t2) * (t2)))));
      double cse_p13_t1_g28236 = (cse_p11_t4_invr27997 * (-(t1)));
      double cse_p12_t75_g28073 = (cse_p12_t167_g28165 * cse_p13_t1_g28236);
      double cse_p12_t76_g28074 = (cse_p12_t168_g28166 * cse_p13_t1_g28236);
      double cse_p12_t77_g28075 = (cse_p12_t169_g28167 * cse_p13_t1_g28236);
      double cse_p12_t78_g28076 = (cse_p12_t170_g28168 * cse_p13_t1_g28236);
      double cse_p12_t79_g28077 = (cse_p12_t171_g28169 * cse_p13_t1_g28236);
      double cse_p12_t80_g28078 = (cse_p12_t172_g28170 * cse_p13_t1_g28236);
      double cse_p12_t81_g28079 = (cse_p12_t173_g28171 * cse_p13_t1_g28236);
      double cse_p12_t82_g28080 = (cse_p12_t174_g28172 * cse_p13_t1_g28236);
      double cse_p12_t83_g28081 = (cse_p12_t175_g28173 * cse_p13_t1_g28236);
      double cse_p12_t84_g28082 = (cse_p12_t176_g28174 * cse_p13_t1_g28236);
      double cse_p12_t85_g28083 = (cse_p12_t178_g28176 * cse_p13_t1_g28236);
      double cse_p12_t86_g28084 = (cse_p12_t179_g28177 * cse_p13_t1_g28236);
      double cse_p12_t1_g28225 = (cse_p11_t4_invr27997 * t2);
      double cse_p12_t87_g28085 = (cse_p12_t147_g28145 * cse_p12_t1_g28225);
      double cse_p12_t88_g28086 = (cse_p12_t148_g28146 * cse_p12_t1_g28225);
      double cse_p12_t89_g28087 = (cse_p12_t149_g28147 * cse_p12_t1_g28225);
      double cse_p12_t90_g28088 = (cse_p12_t150_g28148 * cse_p12_t1_g28225);
      double cse_p12_t91_g28089 = (cse_p12_t151_g28149 * cse_p12_t1_g28225);
      double cse_p12_t92_g28090 = (cse_p12_t152_g28150 * cse_p12_t1_g28225);
      double cse_p13_t1_g28237 = (cse_p11_t1_sqrt27994 * cse_p12_t1_g28225);
      double cse_p12_t93_g28091 = (cse_p12_t182_g28180 * cse_p13_t1_g28237);
      double cse_p12_t94_g28092 = (cse_p12_t183_g28181 * cse_p13_t1_g28237);
      double cse_p12_t95_g28093 = (cse_p12_t185_g28183 * cse_p13_t1_g28237);
      double cse_p12_t96_g28094 = (cse_p13_t1_g28237 * (-(c2x)));
      double cse_p12_t97_g28095 = (cse_p13_t1_g28237 * (-(c2y)));
      double cse_p12_t98_g28096 = (cse_p13_t1_g28237 * (-(c2z)));
      double cse_p12_t153_g28151 = (cse_p12_t75_g28073 + cse_p12_t92_g28090);
      double cse_p12_t154_g28152 = (cse_p12_t76_g28074 + cse_p12_t91_g28089);
      double cse_p12_t155_g28153 = (cse_p12_t77_g28075 + cse_p12_t93_g28091);
      double cse_p12_t156_g28154 = (cse_p12_t78_g28076 + cse_p12_t88_g28086);
      double cse_p12_t157_g28155 = (cse_p12_t79_g28077 + cse_p12_t94_g28092);
      double cse_p12_t158_g28156 = (cse_p12_t80_g28078 + cse_p12_t89_g28087);
      double cse_p12_t159_g28157 = (cse_p12_t81_g28079 + cse_p12_t90_g28088);
      double cse_p12_t160_g28158 = (cse_p12_t82_g28080 + cse_p12_t95_g28093);
      double cse_p12_t161_g28159 = (cse_p12_t83_g28081 + cse_p12_t87_g28085);
      double cse_p12_t162_g28160 = (cse_p12_t84_g28082 + cse_p12_t97_g28095);
      double cse_p12_t163_g28161 = (cse_p12_t85_g28083 + cse_p12_t98_g28096);
      double cse_p12_t164_g28162 = (cse_p12_t86_g28084 + cse_p12_t96_g28094);
      double cse_p11_t5_invr227998 = (cse_p11_t4_invr27997 * cse_p11_t4_invr27997);
      double cse_p12_t99_g28097 = (cse_p11_t5_invr227998 * cse_p12_t138_g28136);
      double cse_p12_t100_g28098 = (cse_p11_t5_invr227998 * cse_p12_t139_g28137);
      double cse_p12_t101_g28099 = (cse_p11_t5_invr227998 * cse_p12_t140_g28138);
      double cse_p12_t102_g28100 = (cse_p11_t5_invr227998 * cse_p12_t141_g28139);
      double cse_p12_t103_g28101 = (cse_p11_t5_invr227998 * cse_p12_t142_g28140);
      double cse_p12_t104_g28102 = (cse_p11_t5_invr227998 * cse_p12_t143_g28141);
      double cse_p12_t105_g28103 = (cse_p11_t5_invr227998 * cse_p12_t144_g28142);
      double cse_p12_t106_g28104 = (cse_p11_t5_invr227998 * cse_p12_t145_g28143);
      double cse_p12_t107_g28105 = (cse_p11_t5_invr227998 * cse_p12_t146_g28144);
      double cse_p12_t108_g28106 = (cse_p11_t5_invr227998 * cse_p12_t165_g28163);
      double cse_p12_t109_g28107 = (cse_p11_t5_invr227998 * cse_p12_t166_g28164);
      double cse_p12_t202_g28200 = (-(cse_p12_t99_g28097));
      double cse_p12_t203_g28201 = (-(cse_p12_t100_g28098));
      double cse_p12_t204_g28202 = (-(cse_p12_t101_g28099));
      double cse_p12_t205_g28203 = (-(cse_p12_t102_g28100));
      double cse_p12_t206_g28204 = (-(cse_p12_t103_g28101));
      double cse_p12_t207_g28205 = (-(cse_p12_t104_g28102));
      double cse_p12_t208_g28206 = (-(cse_p12_t105_g28103));
      double cse_p12_t209_g28207 = (-(cse_p12_t106_g28104));
      double cse_p12_t210_g28208 = (-(cse_p12_t107_g28105));
      double cse_p12_t211_g28209 = (-(cse_p12_t108_g28106));
      double cse_p12_t212_g28210 = (-(cse_p12_t109_g28107));
      double cse_p12_t1_g28228 = (cse_p11_t1_sqrt27994 * t2);
      double cse_p13_t1_g28241 = (cse_p12_t1_g28228 * (-(c2x)));
      double cse_p13_t1_g28242 = (cos_nphi * n * sinphase);
      double cse_p14_t1_g28249 = (cosphase * cse_p12_t223_g28221 * n);
      double g_x1 = (v * ((cse_p13_t1_g28242 * (cse_p12_t86_g28084 + (cse_p11_t4_invr27997 * cse_p13_t1_g28241))) + (cse_p14_t1_g28249 * (cse_p12_t86_g28084 + (cse_p11_t4_invr27997 * cse_p13_t1_g28241)))));
      KernelGradientAcc(i3x1, 0, g_x1);
      double cse_p13_t1_g28243 = (cse_p12_t1_g28228 * (-(c2y)));
      double g_y1 = (v * ((cse_p13_t1_g28242 * (cse_p12_t84_g28082 + (cse_p11_t4_invr27997 * cse_p13_t1_g28243))) + (cse_p14_t1_g28249 * (cse_p12_t84_g28082 + (cse_p11_t4_invr27997 * cse_p13_t1_g28243)))));
      KernelGradientAcc(i3x1, 1, g_y1);
      double cse_p14_t1_g28245 = (cse_p12_t1_g28228 * (-(c2z)));
      double g_z1 = (v * ((cse_p13_t1_g28242 * (cse_p12_t85_g28083 + (cse_p11_t4_invr27997 * cse_p14_t1_g28245))) + (cse_p14_t1_g28249 * (cse_p12_t85_g28083 + (cse_p11_t4_invr27997 * cse_p14_t1_g28245)))));
      KernelGradientAcc(i3x1, 2, g_z1);
      double g_x2 = (v * ((cse_p13_t1_g28242 * (cse_p12_t83_g28081 + (cse_p11_t4_invr27997 * t2 * (cse_p12_t134_g28132 + (cse_p11_t2_invsqrt27995 * cse_p12_t177_g28175 * (-(v2x))))))) + (cse_p14_t1_g28249 * (cse_p12_t83_g28081 + (cse_p11_t4_invr27997 * t2 * (cse_p12_t134_g28132 + (cse_p11_t2_invsqrt27995 * cse_p12_t177_g28175 * (-(v2x)))))))));
      KernelGradientAcc(i3x2, 0, g_x2);
      double g_y2 = (v * ((cse_p13_t1_g28242 * (cse_p12_t78_g28076 + (cse_p11_t4_invr27997 * t2 * (cse_p12_t135_g28133 + (cse_p11_t2_invsqrt27995 * cse_p12_t177_g28175 * (-(v2y))))))) + (cse_p14_t1_g28249 * (cse_p12_t78_g28076 + (cse_p11_t4_invr27997 * t2 * (cse_p12_t135_g28133 + (cse_p11_t2_invsqrt27995 * cse_p12_t177_g28175 * (-(v2y)))))))));
      KernelGradientAcc(i3x2, 1, g_y2);
      double g_z2 = (v * ((cse_p13_t1_g28242 * (cse_p12_t80_g28078 + (cse_p11_t4_invr27997 * t2 * (cse_p12_t136_g28134 + (cse_p11_t2_invsqrt27995 * cse_p12_t177_g28175 * (-(v2z))))))) + (cse_p14_t1_g28249 * (cse_p12_t80_g28078 + (cse_p11_t4_invr27997 * t2 * (cse_p12_t136_g28134 + (cse_p11_t2_invsqrt27995 * cse_p12_t177_g28175 * (-(v2z)))))))));
      KernelGradientAcc(i3x2, 2, g_z2);
      double g_x3 = (v * ((cse_p13_t1_g28242 * (cse_p12_t81_g28079 + (cse_p11_t4_invr27997 * t2 * (cse_p12_t1_g28231 + (cse_p11_t2_invsqrt27995 * cse_p12_t177_g28175 * v2x))))) + (cse_p14_t1_g28249 * (cse_p12_t81_g28079 + (cse_p11_t4_invr27997 * t2 * (cse_p12_t1_g28231 + (cse_p11_t2_invsqrt27995 * cse_p12_t177_g28175 * v2x)))))));
      KernelGradientAcc(i3x3, 0, g_x3);
      double g_y3 = (v * ((cse_p13_t1_g28242 * (cse_p12_t76_g28074 + (cse_p11_t4_invr27997 * t2 * (cse_p12_t1_g28230 + (cse_p11_t2_invsqrt27995 * cse_p12_t177_g28175 * v2y))))) + (cse_p14_t1_g28249 * (cse_p12_t76_g28074 + (cse_p11_t4_invr27997 * t2 * (cse_p12_t1_g28230 + (cse_p11_t2_invsqrt27995 * cse_p12_t177_g28175 * v2y)))))));
      KernelGradientAcc(i3x3, 1, g_y3);
      double g_z3 = (v * ((cse_p13_t1_g28242 * (cse_p12_t75_g28073 + (cse_p11_t4_invr27997 * t2 * (cse_p12_t1_g28229 + (cse_p11_t2_invsqrt27995 * cse_p12_t177_g28175 * v2z))))) + (cse_p14_t1_g28249 * (cse_p12_t75_g28073 + (cse_p11_t4_invr27997 * t2 * (cse_p12_t1_g28229 + (cse_p11_t2_invsqrt27995 * cse_p12_t177_g28175 * v2z)))))));
      KernelGradientAcc(i3x3, 2, g_z3);
      double cse_p12_t1_g28234 = (cse_p11_t4_invr27997 * cse_p12_t185_g28183 * cse_p12_t1_g28228);
      double g_x4 = (v * ((cse_p13_t1_g28242 * (cse_p12_t1_g28234 + cse_p12_t82_g28080)) + (cse_p14_t1_g28249 * (cse_p12_t1_g28234 + cse_p12_t82_g28080))));
      KernelGradientAcc(i3x4, 0, g_x4);
      double cse_p12_t1_g28232 = (cse_p11_t4_invr27997 * cse_p12_t182_g28180 * cse_p12_t1_g28228);
      double g_y4 = (v * ((cse_p13_t1_g28242 * (cse_p12_t1_g28232 + cse_p12_t77_g28075)) + (cse_p14_t1_g28249 * (cse_p12_t1_g28232 + cse_p12_t77_g28075))));
      KernelGradientAcc(i3x4, 1, g_y4);
      double cse_p12_t1_g28233 = (cse_p11_t4_invr27997 * cse_p12_t183_g28181 * cse_p12_t1_g28228);
      double g_z4 = (v * ((cse_p13_t1_g28242 * (cse_p12_t1_g28233 + cse_p12_t79_g28077)) + (cse_p14_t1_g28249 * (cse_p12_t1_g28233 + cse_p12_t79_g28077))));
      KernelGradientAcc(i3x4, 2, g_z4);
      double cse_p12_t1_g28223 = (n * (-(v)));
      double cse_p12_t1_g28226 = (cos_angle * cse_p12_t1_g28223 * n);
      double cse_p12_t1_g28227 = (cse_p12_t1_g28223 * sin_angle);
      double cse_p13_t1_g28235 = (cse_p11_t1_sqrt27994 * cse_p11_t4_invr27997);
      double cse_p13_t1_g28244 = (cse_p13_t1_g28235 * (-(c2x)));
      double h_x1_x1 = ((cse_p12_t164_g28162 * cse_p12_t164_g28162 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p12_t179_g28177 * (cse_p13_t1_g28244 + (-((cse_p11_t1_sqrt27994 * cse_p11_t4_invr27997 * (-(c2x))))) + (-((cse_p11_t5_invr227998 * ((cse_p11_t1_sqrt27994 * cse_p13_t1_g28238 * (-(c2x))) + (cse_p12_t179_g28177 * cse_p13_t1_g28239)) * (-(t1))))))) + (-((cse_p11_t5_invr227998 * cse_p13_t1_g28241 * ((cse_p11_t1_sqrt27994 * cse_p13_t1_g28238 * (-(c2x))) + (cse_p12_t179_g28177 * cse_p13_t1_g28239))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 0, h_x1_x1);
      double h_x1_y1 = ((cse_p12_t162_g28160 * cse_p12_t164_g28162 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p12_t176_g28174 * cse_p13_t1_g28244) + (cse_p12_t179_g28177 * cse_p12_t211_g28209 * (-(t1))) + (cse_p12_t211_g28209 * cse_p13_t1_g28241) + (-((cse_p11_t1_sqrt27994 * cse_p11_t4_invr27997 * cse_p12_t179_g28177 * (-(c2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 1, h_x1_y1);
      double h_x1_z1 = ((cse_p12_t163_g28161 * cse_p12_t164_g28162 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t179_g28177 * cse_p12_t222_g28220) + (cse_p12_t178_g28176 * cse_p13_t1_g28244) + (cse_p12_t179_g28177 * cse_p12_t212_g28210 * (-(t1))) + (cse_p12_t212_g28210 * cse_p13_t1_g28241))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x1, 2, h_x1_z1);
      double cse_p13_t1_g28240 = (cse_p11_t2_invsqrt27995 * cse_p12_t1_g28225);
      double h_x1_x2 = ((cse_p12_t161_g28159 * cse_p12_t164_g28162 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t179_g28177 * cse_p12_t216_g28214) + (cse_p12_t175_g28173 * cse_p13_t1_g28244) + (cse_p12_t179_g28177 * cse_p12_t202_g28200 * (-(t1))) + (cse_p12_t202_g28200 * cse_p13_t1_g28241) + (cse_p13_t1_g28236 * (((-(v2y)) * (-(v3y))) + (v2z * v3z))) + (-((cse_p13_t1_g28240 * v2x * (-(c2x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 0, h_x1_x2);
      double h_x1_y2 = ((cse_p12_t156_g28154 * cse_p12_t164_g28162 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t179_g28177 * cse_p12_t217_g28215) + (cse_p12_t170_g28168 * cse_p13_t1_g28244) + (cse_p12_t179_g28177 * cse_p12_t203_g28201 * (-(t1))) + (cse_p12_t203_g28201 * cse_p13_t1_g28241) + (cse_p13_t1_g28236 * (c2z + (v3x * (-(v2y))))) + (cse_p13_t1_g28237 * v3z) + (-((cse_p13_t1_g28240 * v2y * (-(c2x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 1, h_x1_y2);
      double h_x1_z2 = ((cse_p12_t158_g28156 * cse_p12_t164_g28162 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t179_g28177 * cse_p12_t218_g28216) + (cse_p12_t172_g28170 * cse_p13_t1_g28244) + (cse_p12_t179_g28177 * cse_p12_t204_g28202 * (-(t1))) + (cse_p12_t204_g28202 * cse_p13_t1_g28241) + (cse_p13_t1_g28236 * ((v2z * (-(v3x))) + (-(c2y)))) + (cse_p13_t1_g28237 * (-(v3y))) + (-((cse_p13_t1_g28240 * v2z * (-(c2x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x2, 2, h_x1_z2);
      double h_x1_x3 = ((cse_p12_t159_g28157 * cse_p12_t164_g28162 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t179_g28177 * cse_p12_t219_g28217) + (cse_p12_t173_g28171 * cse_p13_t1_g28244) + (cse_p12_t179_g28177 * cse_p12_t205_g28203 * (-(t1))) + (cse_p12_t205_g28203 * cse_p13_t1_g28241) + (cse_p13_t1_g28236 * ((cse_p12_t191_g28189 * v2z) + (cse_p12_t199_g28197 * (-(v2y))))) + (cse_p13_t1_g28240 * v2x * (-(c2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x3, 0, h_x1_x3);
      double h_x1_y3 = ((cse_p12_t154_g28152 * cse_p12_t164_g28162 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t179_g28177 * cse_p12_t220_g28218) + (cse_p12_t168_g28166 * cse_p13_t1_g28244) + (cse_p12_t179_g28177 * cse_p12_t206_g28204 * (-(t1))) + (cse_p12_t191_g28189 * cse_p13_t1_g28237) + (cse_p12_t206_g28204 * cse_p13_t1_g28241) + (cse_p13_t1_g28236 * ((cse_p12_t189_g28187 * (-(v2y))) + (-(c2z)))) + (cse_p13_t1_g28240 * v2y * (-(c2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x3, 1, h_x1_y3);
      double cse_p14_t1_g28252 = (cse_p12_t207_g28205 * (-(t1)));
      double h_x1_z3 = ((cse_p12_t153_g28151 * cse_p12_t164_g28162 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t179_g28177 * cse_p12_t221_g28219) + (cse_p12_t167_g28165 * cse_p13_t1_g28244) + (cse_p12_t179_g28177 * cse_p14_t1_g28252) + (cse_p12_t199_g28197 * cse_p13_t1_g28237) + (cse_p12_t207_g28205 * cse_p13_t1_g28241) + (cse_p13_t1_g28236 * (c2y + (cse_p12_t198_g28196 * v2z))) + (cse_p13_t1_g28240 * v2z * (-(c2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x3, 2, h_x1_z3);
      double cse_p14_t1_g28251 = (cse_p12_t210_g28208 * (-(t1)));
      double h_x1_x4 = ((cse_p12_t160_g28158 * cse_p12_t164_g28162 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t179_g28177 * cse_p12_t215_g28213) + (cse_p12_t174_g28172 * cse_p13_t1_g28244) + (cse_p12_t179_g28177 * cse_p14_t1_g28251) + (cse_p12_t210_g28208 * cse_p13_t1_g28241) + (cse_p13_t1_g28236 * (((-(v2y)) * (-(v2y))) + (v2z * v2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x4, 0, h_x1_x4);
      double cse_p14_t1_g28248 = (cse_p12_t208_g28206 * (-(t1)));
      double h_x1_y4 = ((cse_p12_t155_g28153 * cse_p12_t164_g28162 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t179_g28177 * cse_p12_t213_g28211) + (cse_p12_t169_g28167 * cse_p13_t1_g28244) + (cse_p12_t179_g28177 * cse_p14_t1_g28248) + (cse_p12_t208_g28206 * cse_p13_t1_g28241) + (cse_p13_t1_g28236 * v2x * (-(v2y))) + (cse_p13_t1_g28237 * v2z))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x4, 1, h_x1_y4);
      double cse_p14_t1_g28246 = (cse_p12_t209_g28207 * (-(t1)));
      double h_x1_z4 = ((cse_p12_t157_g28155 * cse_p12_t164_g28162 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t179_g28177 * cse_p12_t214_g28212) + (cse_p12_t171_g28169 * cse_p13_t1_g28244) + (cse_p12_t179_g28177 * cse_p14_t1_g28246) + (cse_p12_t209_g28207 * cse_p13_t1_g28241) + (cse_p13_t1_g28236 * v2z * (-(v2x))) + (cse_p13_t1_g28237 * (-(v2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 0, i3x4, 2, h_x1_z4);
      double cse_p14_t1_g28247 = (cse_p13_t1_g28235 * (-(c2y)));
      double h_y1_y1 = ((cse_p12_t162_g28160 * cse_p12_t162_g28160 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p12_t176_g28174 * (cse_p14_t1_g28247 + (cse_p12_t211_g28209 * (-(t1))) + (-((cse_p11_t1_sqrt27994 * cse_p11_t4_invr27997 * (-(c2y))))))) + (cse_p12_t211_g28209 * cse_p13_t1_g28243))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 1, h_y1_y1);
      double h_y1_z1 = ((cse_p12_t162_g28160 * cse_p12_t163_g28161 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t176_g28174 * cse_p12_t222_g28220) + (cse_p12_t176_g28174 * cse_p12_t212_g28210 * (-(t1))) + (cse_p12_t178_g28176 * cse_p14_t1_g28247) + (cse_p12_t212_g28210 * cse_p13_t1_g28243))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x1, 2, h_y1_z1);
      double h_y1_x2 = ((cse_p12_t161_g28159 * cse_p12_t162_g28160 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t176_g28174 * cse_p12_t216_g28214) + (cse_p12_t175_g28173 * cse_p14_t1_g28247) + (cse_p12_t176_g28174 * cse_p12_t202_g28200 * (-(t1))) + (cse_p12_t202_g28200 * cse_p13_t1_g28243) + (cse_p13_t1_g28236 * ((v2x * (-(v3y))) + (-(c2z)))) + (cse_p13_t1_g28237 * (-(v3z))) + (-((cse_p13_t1_g28240 * v2x * (-(c2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 0, h_y1_x2);
      double h_y1_y2 = ((cse_p12_t156_g28154 * cse_p12_t162_g28160 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t176_g28174 * cse_p12_t217_g28215) + (cse_p12_t170_g28168 * cse_p14_t1_g28247) + (cse_p12_t176_g28174 * cse_p12_t203_g28201 * (-(t1))) + (cse_p12_t203_g28201 * cse_p13_t1_g28243) + (cse_p13_t1_g28236 * (((-(v2z)) * (-(v3z))) + (v2x * v3x))) + (-((cse_p13_t1_g28240 * v2y * (-(c2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 1, h_y1_y2);
      double h_y1_z2 = ((cse_p12_t158_g28156 * cse_p12_t162_g28160 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t176_g28174 * cse_p12_t218_g28216) + (cse_p12_t172_g28170 * cse_p14_t1_g28247) + (cse_p12_t176_g28174 * cse_p12_t204_g28202 * (-(t1))) + (cse_p12_t204_g28202 * cse_p13_t1_g28243) + (cse_p13_t1_g28236 * (c2x + (v3y * (-(v2z))))) + (cse_p13_t1_g28237 * v3x) + (-((cse_p13_t1_g28240 * v2z * (-(c2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x2, 2, h_y1_z2);
      double h_y1_x3 = ((cse_p12_t159_g28157 * cse_p12_t162_g28160 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t176_g28174 * cse_p12_t219_g28217) + (cse_p12_t173_g28171 * cse_p14_t1_g28247) + (cse_p12_t176_g28174 * cse_p12_t205_g28203 * (-(t1))) + (cse_p12_t200_g28198 * cse_p13_t1_g28237) + (cse_p12_t205_g28203 * cse_p13_t1_g28243) + (cse_p13_t1_g28236 * (c2z + (cse_p12_t199_g28197 * v2x))) + (cse_p13_t1_g28240 * v2x * (-(c2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x3, 0, h_y1_x3);
      double h_y1_y3 = ((cse_p12_t154_g28152 * cse_p12_t162_g28160 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t176_g28174 * cse_p12_t220_g28218) + (cse_p12_t168_g28166 * cse_p14_t1_g28247) + (cse_p12_t176_g28174 * cse_p12_t206_g28204 * (-(t1))) + (cse_p12_t206_g28204 * cse_p13_t1_g28243) + (cse_p13_t1_g28236 * ((cse_p12_t189_g28187 * v2x) + (cse_p12_t200_g28198 * (-(v2z))))) + (cse_p13_t1_g28240 * v2y * (-(c2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x3, 1, h_y1_y3);
      double h_y1_z3 = ((cse_p12_t153_g28151 * cse_p12_t162_g28160 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t176_g28174 * cse_p12_t221_g28219) + (cse_p12_t167_g28165 * cse_p14_t1_g28247) + (cse_p12_t176_g28174 * cse_p14_t1_g28252) + (cse_p12_t189_g28187 * cse_p13_t1_g28237) + (cse_p12_t207_g28205 * cse_p13_t1_g28243) + (cse_p13_t1_g28236 * ((cse_p12_t190_g28188 * (-(v2z))) + (-(c2x)))) + (cse_p13_t1_g28240 * v2z * (-(c2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x3, 2, h_y1_z3);
      double h_y1_x4 = ((cse_p12_t160_g28158 * cse_p12_t162_g28160 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t176_g28174 * cse_p12_t215_g28213) + (cse_p12_t174_g28172 * cse_p14_t1_g28247) + (cse_p12_t176_g28174 * cse_p14_t1_g28251) + (cse_p12_t210_g28208 * cse_p13_t1_g28243) + (cse_p13_t1_g28236 * v2x * (-(v2y))) + (cse_p13_t1_g28237 * (-(v2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x4, 0, h_y1_x4);
      double h_y1_y4 = ((cse_p12_t155_g28153 * cse_p12_t162_g28160 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t176_g28174 * cse_p12_t213_g28211) + (cse_p12_t169_g28167 * cse_p14_t1_g28247) + (cse_p12_t176_g28174 * cse_p14_t1_g28248) + (cse_p12_t208_g28206 * cse_p13_t1_g28243) + (cse_p13_t1_g28236 * (((-(v2z)) * (-(v2z))) + (v2x * v2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x4, 1, h_y1_y4);
      double h_y1_z4 = ((cse_p12_t157_g28155 * cse_p12_t162_g28160 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t176_g28174 * cse_p12_t214_g28212) + (cse_p12_t171_g28169 * cse_p14_t1_g28247) + (cse_p12_t176_g28174 * cse_p14_t1_g28246) + (cse_p12_t209_g28207 * cse_p13_t1_g28243) + (cse_p13_t1_g28236 * v2y * (-(v2z))) + (cse_p13_t1_g28237 * v2x))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 1, i3x4, 2, h_y1_z4);
      double cse_p14_t1_g28250 = (cse_p13_t1_g28235 * (-(c2z)));
      double h_z1_z1 = ((cse_p12_t163_g28161 * cse_p12_t163_g28161 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p12_t178_g28176 * (cse_p14_t1_g28250 + (cse_p11_t4_invr27997 * cse_p12_t222_g28220) + (cse_p12_t212_g28210 * (-(t1))))) + (cse_p12_t212_g28210 * cse_p14_t1_g28245))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x1, 2, h_z1_z1);
      double h_z1_x2 = ((cse_p12_t161_g28159 * cse_p12_t163_g28161 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t178_g28176 * cse_p12_t216_g28214) + (cse_p12_t175_g28173 * cse_p14_t1_g28250) + (cse_p12_t178_g28176 * cse_p12_t202_g28200 * (-(t1))) + (cse_p12_t202_g28200 * cse_p14_t1_g28245) + (cse_p13_t1_g28236 * (c2y + (v3z * (-(v2x))))) + (cse_p13_t1_g28237 * v3y) + (-((cse_p13_t1_g28240 * v2x * (-(c2z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 0, h_z1_x2);
      double h_z1_y2 = ((cse_p12_t156_g28154 * cse_p12_t163_g28161 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t178_g28176 * cse_p12_t217_g28215) + (cse_p12_t170_g28168 * cse_p14_t1_g28250) + (cse_p12_t178_g28176 * cse_p12_t203_g28201 * (-(t1))) + (cse_p12_t203_g28201 * cse_p14_t1_g28245) + (cse_p13_t1_g28236 * ((v2y * (-(v3z))) + (-(c2x)))) + (cse_p13_t1_g28237 * (-(v3x))) + (-((cse_p13_t1_g28240 * v2y * (-(c2z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 1, h_z1_y2);
      double h_z1_z2 = ((cse_p12_t158_g28156 * cse_p12_t163_g28161 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t178_g28176 * cse_p12_t218_g28216) + (cse_p12_t172_g28170 * cse_p14_t1_g28250) + (cse_p12_t178_g28176 * cse_p12_t204_g28202 * (-(t1))) + (cse_p12_t204_g28202 * cse_p14_t1_g28245) + (cse_p13_t1_g28236 * (((-(v2x)) * (-(v3x))) + (v2y * v3y))) + (-((cse_p13_t1_g28240 * v2z * (-(c2z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x2, 2, h_z1_z2);
      double h_z1_x3 = ((cse_p12_t159_g28157 * cse_p12_t163_g28161 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t178_g28176 * cse_p12_t219_g28217) + (cse_p12_t173_g28171 * cse_p14_t1_g28250) + (cse_p12_t178_g28176 * cse_p12_t205_g28203 * (-(t1))) + (cse_p12_t190_g28188 * cse_p13_t1_g28237) + (cse_p12_t205_g28203 * cse_p14_t1_g28245) + (cse_p13_t1_g28236 * ((cse_p12_t191_g28189 * (-(v2x))) + (-(c2y)))) + (cse_p13_t1_g28240 * v2x * (-(c2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x3, 0, h_z1_x3);
      double h_z1_y3 = ((cse_p12_t154_g28152 * cse_p12_t163_g28161 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t178_g28176 * cse_p12_t220_g28218) + (cse_p12_t168_g28166 * cse_p14_t1_g28250) + (cse_p12_t178_g28176 * cse_p12_t206_g28204 * (-(t1))) + (cse_p12_t198_g28196 * cse_p13_t1_g28237) + (cse_p12_t206_g28204 * cse_p14_t1_g28245) + (cse_p13_t1_g28236 * (c2x + (cse_p12_t200_g28198 * v2y))) + (cse_p13_t1_g28240 * v2y * (-(c2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x3, 1, h_z1_y3);
      double h_z1_z3 = ((cse_p12_t153_g28151 * cse_p12_t163_g28161 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t178_g28176 * cse_p12_t221_g28219) + (cse_p12_t167_g28165 * cse_p14_t1_g28250) + (cse_p12_t178_g28176 * cse_p14_t1_g28252) + (cse_p12_t207_g28205 * cse_p14_t1_g28245) + (cse_p13_t1_g28236 * ((cse_p12_t190_g28188 * v2y) + (cse_p12_t198_g28196 * (-(v2x))))) + (cse_p13_t1_g28240 * v2z * (-(c2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x3, 2, h_z1_z3);
      double h_z1_x4 = ((cse_p12_t160_g28158 * cse_p12_t163_g28161 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t178_g28176 * cse_p12_t215_g28213) + (cse_p12_t174_g28172 * cse_p14_t1_g28250) + (cse_p12_t178_g28176 * cse_p14_t1_g28251) + (cse_p12_t210_g28208 * cse_p14_t1_g28245) + (cse_p13_t1_g28236 * v2z * (-(v2x))) + (cse_p13_t1_g28237 * v2y))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x4, 0, h_z1_x4);
      double h_z1_y4 = ((cse_p12_t155_g28153 * cse_p12_t163_g28161 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t178_g28176 * cse_p12_t213_g28211) + (cse_p12_t169_g28167 * cse_p14_t1_g28250) + (cse_p12_t178_g28176 * cse_p14_t1_g28248) + (cse_p12_t208_g28206 * cse_p14_t1_g28245) + (cse_p13_t1_g28236 * v2y * (-(v2z))) + (cse_p13_t1_g28237 * (-(v2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x4, 1, h_z1_y4);
      double h_z1_z4 = ((cse_p12_t157_g28155 * cse_p12_t163_g28161 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t178_g28176 * cse_p12_t214_g28212) + (cse_p12_t171_g28169 * cse_p14_t1_g28250) + (cse_p12_t178_g28176 * cse_p14_t1_g28246) + (cse_p12_t209_g28207 * cse_p14_t1_g28245) + (cse_p13_t1_g28236 * (((-(v2x)) * (-(v2x))) + (v2y * v2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x1, 2, i3x4, 2, h_z1_z4);
      double cse_p14_t1_g28253 = (cse_p11_t2_invsqrt27995 * (-(v2x)));
      double h_x2_x2 = ((cse_p12_t161_g28159 * cse_p12_t161_g28159 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p12_t147_g28145 * cse_p12_t202_g28200 * t2) + (cse_p12_t175_g28173 * ((cse_p11_t4_invr27997 * cse_p12_t147_g28145) + (cse_p11_t4_invr27997 * cse_p12_t216_g28214) + (cse_p12_t202_g28200 * (-(t1))))) + (cse_p12_t1_g28225 * (cse_p12_t1_g28222 + (cse_p12_t192_g28190 * cse_p14_t1_g28253) + (cse_p12_t192_g28190 * cse_p14_t1_g28253) + (-(cse_p12_t66_g28064)))) + (cse_p13_t1_g28236 * ((cse_p12_t188_g28186 * v3z) + (cse_p12_t188_g28186 * v3z) + (cse_p12_t196_g28194 * (-(v3y))) + (cse_p12_t196_g28194 * (-(v3y))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 0, h_x2_x2);
      double cse_p14_t1_g28254 = (cse_p11_t2_invsqrt27995 * (-(v2y)));
      double h_x2_y2 = ((cse_p12_t156_g28154 * cse_p12_t161_g28159 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t147_g28145 * cse_p12_t170_g28168) + (cse_p11_t4_invr27997 * cse_p12_t175_g28173 * cse_p12_t217_g28215) + (cse_p12_t147_g28145 * cse_p12_t203_g28201 * t2) + (cse_p12_t175_g28173 * cse_p12_t203_g28201 * (-(t1))) + (cse_p12_t1_g28225 * ((cse_p12_t192_g28190 * cse_p14_t1_g28254) + (cse_p12_t193_g28191 * cse_p14_t1_g28253) + (-(cse_p12_t67_g28065)))) + (cse_p13_t1_g28236 * ((cse_p12_t186_g28184 * (-(v3y))) + (cse_p12_t196_g28194 * v3x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 1, h_x2_y2);
      double cse_p15_t1_g28255 = (cse_p11_t2_invsqrt27995 * (-(v2z)));
      double h_x2_z2 = ((cse_p12_t158_g28156 * cse_p12_t161_g28159 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t147_g28145 * cse_p12_t172_g28170) + (cse_p11_t4_invr27997 * cse_p12_t175_g28173 * cse_p12_t218_g28216) + (cse_p12_t147_g28145 * cse_p12_t204_g28202 * t2) + (cse_p12_t175_g28173 * cse_p12_t204_g28202 * (-(t1))) + (cse_p12_t1_g28225 * ((cse_p12_t192_g28190 * cse_p15_t1_g28255) + (cse_p12_t194_g28192 * cse_p14_t1_g28253) + (-(cse_p12_t68_g28066)))) + (cse_p13_t1_g28236 * ((cse_p12_t188_g28186 * (-(v3x))) + (cse_p12_t195_g28193 * v3z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x2, 2, h_x2_z2);
      double h_x2_x3 = ((cse_p12_t159_g28157 * cse_p12_t161_g28159 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t147_g28145 * cse_p12_t173_g28171) + (cse_p11_t4_invr27997 * cse_p12_t175_g28173 * cse_p12_t219_g28217) + (cse_p12_t147_g28145 * cse_p12_t205_g28203 * t2) + (cse_p12_t175_g28173 * cse_p12_t205_g28203 * (-(t1))) + (cse_p12_t1_g28225 * (cse_p12_t201_g28199 + cse_p12_t66_g28064 + (cse_p11_t2_invsqrt27995 * cse_p12_t192_g28190 * v2x) + (cse_p12_t184_g28182 * cse_p14_t1_g28253))) + (cse_p13_t1_g28236 * (((-(v1y)) * (-(v3y))) + (cse_p12_t188_g28186 * cse_p12_t191_g28189) + (cse_p12_t196_g28194 * cse_p12_t199_g28197) + (v1z * v3z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x3, 0, h_x2_x3);
      double h_x2_y3 = ((cse_p12_t154_g28152 * cse_p12_t161_g28159 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t147_g28145 * cse_p12_t168_g28166) + (cse_p11_t4_invr27997 * cse_p12_t175_g28173 * cse_p12_t220_g28218) + (cse_p12_t147_g28145 * cse_p12_t206_g28204 * t2) + (cse_p12_t175_g28173 * cse_p12_t206_g28204 * (-(t1))) + (cse_p12_t1_g28225 * (cse_p12_t67_g28065 + (cse_p11_t1_sqrt27994 * (v1z + v2z + v3z)) + (cse_p11_t2_invsqrt27995 * cse_p12_t192_g28190 * v2y) + (cse_p12_t181_g28179 * cse_p14_t1_g28253))) + (cse_p13_t1_g28236 * (c1z + c2z + (cse_p12_t189_g28187 * cse_p12_t196_g28194) + (v1x * (-(v3y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x3, 1, h_x2_y3);
      double h_x2_z3 = ((cse_p12_t153_g28151 * cse_p12_t161_g28159 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t147_g28145 * cse_p12_t167_g28165) + (cse_p11_t4_invr27997 * cse_p12_t175_g28173 * cse_p12_t221_g28219) + (cse_p12_t147_g28145 * cse_p12_t207_g28205 * t2) + (cse_p12_t175_g28173 * cse_p14_t1_g28252) + (cse_p12_t1_g28225 * (cse_p12_t68_g28066 + (cse_p11_t1_sqrt27994 * ((-(v1y)) + (-(v2y)) + (-(v3y)))) + (cse_p11_t2_invsqrt27995 * cse_p12_t192_g28190 * v2z) + (cse_p12_t180_g28178 * cse_p14_t1_g28253))) + (cse_p13_t1_g28236 * ((cse_p12_t188_g28186 * cse_p12_t198_g28196) + (v3z * (-(v1x))) + (-(c1y)) + (-(c2y)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x3, 2, h_x2_z3);
      double h_x2_x4 = ((cse_p12_t160_g28158 * cse_p12_t161_g28159 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t147_g28145 * cse_p12_t174_g28172) + (cse_p11_t4_invr27997 * cse_p12_t175_g28173 * cse_p12_t215_g28213) + (cse_p12_t147_g28145 * cse_p12_t210_g28208 * t2) + (cse_p12_t175_g28173 * cse_p14_t1_g28251) + (cse_p13_t1_g28236 * ((cse_p12_t188_g28186 * v2z) + (cse_p12_t196_g28194 * (-(v2y))))) + (-((cse_p12_t185_g28183 * cse_p13_t1_g28240 * v2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x4, 0, h_x2_x4);
      double h_x2_y4 = ((cse_p12_t155_g28153 * cse_p12_t161_g28159 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t147_g28145 * cse_p12_t169_g28167) + (cse_p11_t4_invr27997 * cse_p12_t175_g28173 * cse_p12_t213_g28211) + (cse_p12_t147_g28145 * cse_p12_t208_g28206 * t2) + (cse_p12_t175_g28173 * cse_p14_t1_g28248) + (cse_p12_t1_g28225 * ((cse_p11_t1_sqrt27994 * cse_p12_t188_g28186) + (cse_p12_t182_g28180 * cse_p14_t1_g28253))) + (cse_p13_t1_g28236 * ((cse_p12_t196_g28194 * v2x) + (-(c1z)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x4, 1, h_x2_y4);
      double h_x2_z4 = ((cse_p12_t157_g28155 * cse_p12_t161_g28159 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t147_g28145 * cse_p12_t171_g28169) + (cse_p11_t4_invr27997 * cse_p12_t175_g28173 * cse_p12_t214_g28212) + (cse_p12_t147_g28145 * cse_p12_t209_g28207 * t2) + (cse_p12_t175_g28173 * cse_p14_t1_g28246) + (cse_p12_t1_g28225 * ((cse_p11_t1_sqrt27994 * cse_p12_t196_g28194) + (cse_p12_t183_g28181 * cse_p14_t1_g28253))) + (cse_p13_t1_g28236 * (c1y + (cse_p12_t188_g28186 * (-(v2x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 0, i3x4, 2, h_x2_z4);
      double h_y2_y2 = ((cse_p12_t156_g28154 * cse_p12_t156_g28154 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p12_t148_g28146 * cse_p12_t203_g28201 * t2) + (cse_p12_t170_g28168 * ((cse_p11_t4_invr27997 * cse_p12_t148_g28146) + (cse_p11_t4_invr27997 * cse_p12_t217_g28215) + (cse_p12_t203_g28201 * (-(t1))))) + (cse_p12_t1_g28225 * (cse_p12_t1_g28222 + (cse_p12_t193_g28191 * cse_p14_t1_g28254) + (cse_p12_t193_g28191 * cse_p14_t1_g28254) + (-(cse_p12_t69_g28067)))) + (cse_p13_t1_g28236 * ((cse_p12_t186_g28184 * v3x) + (cse_p12_t186_g28184 * v3x) + (cse_p12_t197_g28195 * (-(v3z))) + (cse_p12_t197_g28195 * (-(v3z))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 1, h_y2_y2);
      double h_y2_z2 = ((cse_p12_t156_g28154 * cse_p12_t158_g28156 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t148_g28146 * cse_p12_t172_g28170) + (cse_p11_t4_invr27997 * cse_p12_t170_g28168 * cse_p12_t218_g28216) + (cse_p12_t148_g28146 * cse_p12_t204_g28202 * t2) + (cse_p12_t170_g28168 * cse_p12_t204_g28202 * (-(t1))) + (cse_p12_t1_g28225 * ((cse_p12_t193_g28191 * cse_p15_t1_g28255) + (cse_p12_t194_g28192 * cse_p14_t1_g28254) + (-(cse_p12_t70_g28068)))) + (cse_p13_t1_g28236 * ((cse_p12_t187_g28185 * (-(v3z))) + (cse_p12_t197_g28195 * v3y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x2, 2, h_y2_z2);
      double h_y2_x3 = ((cse_p12_t156_g28154 * cse_p12_t159_g28157 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t148_g28146 * cse_p12_t173_g28171) + (cse_p11_t4_invr27997 * cse_p12_t170_g28168 * cse_p12_t219_g28217) + (cse_p12_t148_g28146 * cse_p12_t205_g28203 * t2) + (cse_p12_t170_g28168 * cse_p12_t205_g28203 * (-(t1))) + (cse_p12_t1_g28225 * (cse_p12_t67_g28065 + (cse_p11_t1_sqrt27994 * ((-(v1z)) + (-(v2z)) + (-(v3z)))) + (cse_p11_t2_invsqrt27995 * cse_p12_t193_g28191 * v2x) + (cse_p12_t184_g28182 * cse_p14_t1_g28254))) + (cse_p13_t1_g28236 * ((cse_p12_t186_g28184 * cse_p12_t199_g28197) + (v3x * (-(v1y))) + (-(c1z)) + (-(c2z)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x3, 0, h_y2_x3);
      double h_y2_y3 = ((cse_p12_t154_g28152 * cse_p12_t156_g28154 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t148_g28146 * cse_p12_t168_g28166) + (cse_p11_t4_invr27997 * cse_p12_t170_g28168 * cse_p12_t220_g28218) + (cse_p12_t148_g28146 * cse_p12_t206_g28204 * t2) + (cse_p12_t170_g28168 * cse_p12_t206_g28204 * (-(t1))) + (cse_p12_t1_g28225 * (cse_p12_t201_g28199 + cse_p12_t69_g28067 + (cse_p11_t2_invsqrt27995 * cse_p12_t193_g28191 * v2y) + (cse_p12_t181_g28179 * cse_p14_t1_g28254))) + (cse_p13_t1_g28236 * (((-(v1z)) * (-(v3z))) + (cse_p12_t186_g28184 * cse_p12_t189_g28187) + (cse_p12_t197_g28195 * cse_p12_t200_g28198) + (v1x * v3x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x3, 1, h_y2_y3);
      double h_y2_z3 = ((cse_p12_t153_g28151 * cse_p12_t156_g28154 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t148_g28146 * cse_p12_t167_g28165) + (cse_p11_t4_invr27997 * cse_p12_t170_g28168 * cse_p12_t221_g28219) + (cse_p12_t148_g28146 * cse_p12_t207_g28205 * t2) + (cse_p12_t170_g28168 * cse_p14_t1_g28252) + (cse_p12_t1_g28225 * (cse_p12_t70_g28068 + (cse_p11_t1_sqrt27994 * (v1x + v2x + v3x)) + (cse_p11_t2_invsqrt27995 * cse_p12_t193_g28191 * v2z) + (cse_p12_t180_g28178 * cse_p14_t1_g28254))) + (cse_p13_t1_g28236 * (c1x + c2x + (cse_p12_t190_g28188 * cse_p12_t197_g28195) + (v1y * (-(v3z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x3, 2, h_y2_z3);
      double h_y2_x4 = ((cse_p12_t156_g28154 * cse_p12_t160_g28158 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t148_g28146 * cse_p12_t174_g28172) + (cse_p11_t4_invr27997 * cse_p12_t170_g28168 * cse_p12_t215_g28213) + (cse_p12_t148_g28146 * cse_p12_t210_g28208 * t2) + (cse_p12_t170_g28168 * cse_p14_t1_g28251) + (cse_p12_t1_g28225 * ((cse_p11_t1_sqrt27994 * cse_p12_t197_g28195) + (cse_p12_t185_g28183 * cse_p14_t1_g28254))) + (cse_p13_t1_g28236 * (c1z + (cse_p12_t186_g28184 * (-(v2y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x4, 0, h_y2_x4);
      double h_y2_y4 = ((cse_p12_t155_g28153 * cse_p12_t156_g28154 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t148_g28146 * cse_p12_t169_g28167) + (cse_p11_t4_invr27997 * cse_p12_t170_g28168 * cse_p12_t213_g28211) + (cse_p12_t148_g28146 * cse_p12_t208_g28206 * t2) + (cse_p12_t170_g28168 * cse_p14_t1_g28248) + (cse_p13_t1_g28236 * ((cse_p12_t186_g28184 * v2x) + (cse_p12_t197_g28195 * (-(v2z))))) + (-((cse_p12_t182_g28180 * cse_p13_t1_g28240 * v2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x4, 1, h_y2_y4);
      double h_y2_z4 = ((cse_p12_t156_g28154 * cse_p12_t157_g28155 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t148_g28146 * cse_p12_t171_g28169) + (cse_p11_t4_invr27997 * cse_p12_t170_g28168 * cse_p12_t214_g28212) + (cse_p12_t148_g28146 * cse_p12_t209_g28207 * t2) + (cse_p12_t170_g28168 * cse_p14_t1_g28246) + (cse_p12_t1_g28225 * ((cse_p11_t1_sqrt27994 * cse_p12_t186_g28184) + (cse_p12_t183_g28181 * cse_p14_t1_g28254))) + (cse_p13_t1_g28236 * ((cse_p12_t197_g28195 * v2y) + (-(c1x)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 1, i3x4, 2, h_y2_z4);
      double h_z2_z2 = ((cse_p12_t158_g28156 * cse_p12_t158_g28156 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p12_t149_g28147 * cse_p12_t204_g28202 * t2) + (cse_p12_t172_g28170 * ((cse_p11_t4_invr27997 * cse_p12_t149_g28147) + (cse_p11_t4_invr27997 * cse_p12_t218_g28216) + (cse_p12_t204_g28202 * (-(t1))))) + (cse_p12_t1_g28225 * (cse_p12_t1_g28222 + (cse_p12_t194_g28192 * cse_p15_t1_g28255) + (cse_p12_t194_g28192 * cse_p15_t1_g28255) + (-(cse_p12_t71_g28069)))) + (cse_p13_t1_g28236 * ((cse_p12_t187_g28185 * v3y) + (cse_p12_t187_g28185 * v3y) + (cse_p12_t195_g28193 * (-(v3x))) + (cse_p12_t195_g28193 * (-(v3x))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x2, 2, h_z2_z2);
      double h_z2_x3 = ((cse_p12_t158_g28156 * cse_p12_t159_g28157 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t149_g28147 * cse_p12_t173_g28171) + (cse_p11_t4_invr27997 * cse_p12_t172_g28170 * cse_p12_t219_g28217) + (cse_p12_t149_g28147 * cse_p12_t205_g28203 * t2) + (cse_p12_t172_g28170 * cse_p12_t205_g28203 * (-(t1))) + (cse_p12_t1_g28225 * (cse_p12_t68_g28066 + (cse_p11_t1_sqrt27994 * (v1y + v2y + v3y)) + (cse_p11_t2_invsqrt27995 * cse_p12_t194_g28192 * v2x) + (cse_p12_t184_g28182 * cse_p15_t1_g28255))) + (cse_p13_t1_g28236 * (c1y + c2y + (cse_p12_t191_g28189 * cse_p12_t195_g28193) + (v1z * (-(v3x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x3, 0, h_z2_x3);
      double h_z2_y3 = ((cse_p12_t154_g28152 * cse_p12_t158_g28156 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t149_g28147 * cse_p12_t168_g28166) + (cse_p11_t4_invr27997 * cse_p12_t172_g28170 * cse_p12_t220_g28218) + (cse_p12_t149_g28147 * cse_p12_t206_g28204 * t2) + (cse_p12_t172_g28170 * cse_p12_t206_g28204 * (-(t1))) + (cse_p12_t1_g28225 * (cse_p12_t70_g28068 + (cse_p11_t1_sqrt27994 * ((-(v1x)) + (-(v2x)) + (-(v3x)))) + (cse_p11_t2_invsqrt27995 * cse_p12_t194_g28192 * v2y) + (cse_p12_t181_g28179 * cse_p15_t1_g28255))) + (cse_p13_t1_g28236 * ((cse_p12_t187_g28185 * cse_p12_t200_g28198) + (v3y * (-(v1z))) + (-(c1x)) + (-(c2x)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x3, 1, h_z2_y3);
      double h_z2_z3 = ((cse_p12_t153_g28151 * cse_p12_t158_g28156 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t149_g28147 * cse_p12_t167_g28165) + (cse_p11_t4_invr27997 * cse_p12_t172_g28170 * cse_p12_t221_g28219) + (cse_p12_t149_g28147 * cse_p12_t207_g28205 * t2) + (cse_p12_t172_g28170 * cse_p14_t1_g28252) + (cse_p12_t1_g28225 * (cse_p12_t201_g28199 + cse_p12_t71_g28069 + (cse_p11_t2_invsqrt27995 * cse_p12_t194_g28192 * v2z) + (cse_p12_t180_g28178 * cse_p15_t1_g28255))) + (cse_p13_t1_g28236 * (((-(v1x)) * (-(v3x))) + (cse_p12_t187_g28185 * cse_p12_t190_g28188) + (cse_p12_t195_g28193 * cse_p12_t198_g28196) + (v1y * v3y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x3, 2, h_z2_z3);
      double h_z2_x4 = ((cse_p12_t158_g28156 * cse_p12_t160_g28158 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t149_g28147 * cse_p12_t174_g28172) + (cse_p11_t4_invr27997 * cse_p12_t172_g28170 * cse_p12_t215_g28213) + (cse_p12_t149_g28147 * cse_p12_t210_g28208 * t2) + (cse_p12_t172_g28170 * cse_p14_t1_g28251) + (cse_p12_t1_g28225 * ((cse_p11_t1_sqrt27994 * cse_p12_t187_g28185) + (cse_p12_t185_g28183 * cse_p15_t1_g28255))) + (cse_p13_t1_g28236 * ((cse_p12_t195_g28193 * v2z) + (-(c1y)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x4, 0, h_z2_x4);
      double h_z2_y4 = ((cse_p12_t155_g28153 * cse_p12_t158_g28156 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t149_g28147 * cse_p12_t169_g28167) + (cse_p11_t4_invr27997 * cse_p12_t172_g28170 * cse_p12_t213_g28211) + (cse_p12_t149_g28147 * cse_p12_t208_g28206 * t2) + (cse_p12_t172_g28170 * cse_p14_t1_g28248) + (cse_p12_t1_g28225 * ((cse_p11_t1_sqrt27994 * cse_p12_t195_g28193) + (cse_p12_t182_g28180 * cse_p15_t1_g28255))) + (cse_p13_t1_g28236 * (c1x + (cse_p12_t187_g28185 * (-(v2z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x4, 1, h_z2_y4);
      double h_z2_z4 = ((cse_p12_t157_g28155 * cse_p12_t158_g28156 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t149_g28147 * cse_p12_t171_g28169) + (cse_p11_t4_invr27997 * cse_p12_t172_g28170 * cse_p12_t214_g28212) + (cse_p12_t149_g28147 * cse_p12_t209_g28207 * t2) + (cse_p12_t172_g28170 * cse_p14_t1_g28246) + (cse_p13_t1_g28236 * ((cse_p12_t187_g28185 * v2y) + (cse_p12_t195_g28193 * (-(v2x))))) + (-((cse_p12_t183_g28181 * cse_p13_t1_g28240 * v2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x2, 2, i3x4, 2, h_z2_z4);
      double h_x3_x3 = ((cse_p12_t159_g28157 * cse_p12_t159_g28157 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p12_t150_g28148 * cse_p12_t205_g28203 * t2) + (cse_p12_t173_g28171 * ((cse_p11_t4_invr27997 * cse_p12_t150_g28148) + (cse_p11_t4_invr27997 * cse_p12_t219_g28217) + (cse_p12_t205_g28203 * (-(t1))))) + (cse_p12_t1_g28225 * (cse_p12_t1_g28222 + (cse_p11_t2_invsqrt27995 * cse_p12_t184_g28182 * v2x) + (cse_p11_t2_invsqrt27995 * cse_p12_t184_g28182 * v2x) + (-(cse_p12_t66_g28064)))) + (cse_p13_t1_g28236 * ((cse_p12_t191_g28189 * v1z) + (cse_p12_t191_g28189 * v1z) + (cse_p12_t199_g28197 * (-(v1y))) + (cse_p12_t199_g28197 * (-(v1y))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 0, i3x3, 0, h_x3_x3);
      double h_x3_y3 = ((cse_p12_t154_g28152 * cse_p12_t159_g28157 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t150_g28148 * cse_p12_t168_g28166) + (cse_p11_t4_invr27997 * cse_p12_t173_g28171 * cse_p12_t220_g28218) + (cse_p12_t150_g28148 * cse_p12_t206_g28204 * t2) + (cse_p12_t173_g28171 * cse_p12_t206_g28204 * (-(t1))) + (cse_p12_t1_g28225 * ((cse_p11_t2_invsqrt27995 * cse_p12_t181_g28179 * v2x) + (cse_p11_t2_invsqrt27995 * cse_p12_t184_g28182 * v2y) + (-(cse_p12_t67_g28065)))) + (cse_p13_t1_g28236 * ((cse_p12_t189_g28187 * (-(v1y))) + (cse_p12_t199_g28197 * v1x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 0, i3x3, 1, h_x3_y3);
      double h_x3_z3 = ((cse_p12_t153_g28151 * cse_p12_t159_g28157 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t150_g28148 * cse_p12_t167_g28165) + (cse_p11_t4_invr27997 * cse_p12_t173_g28171 * cse_p12_t221_g28219) + (cse_p12_t150_g28148 * cse_p12_t207_g28205 * t2) + (cse_p12_t173_g28171 * cse_p14_t1_g28252) + (cse_p12_t1_g28225 * ((cse_p11_t2_invsqrt27995 * cse_p12_t180_g28178 * v2x) + (cse_p11_t2_invsqrt27995 * cse_p12_t184_g28182 * v2z) + (-(cse_p12_t68_g28066)))) + (cse_p13_t1_g28236 * ((cse_p12_t191_g28189 * (-(v1x))) + (cse_p12_t198_g28196 * v1z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 0, i3x3, 2, h_x3_z3);
      double h_x3_x4 = ((cse_p12_t159_g28157 * cse_p12_t160_g28158 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t150_g28148 * cse_p12_t174_g28172) + (cse_p11_t4_invr27997 * cse_p12_t173_g28171 * cse_p12_t215_g28213) + (cse_p12_t150_g28148 * cse_p12_t210_g28208 * t2) + (cse_p12_t173_g28171 * cse_p14_t1_g28251) + (cse_p12_t185_g28183 * cse_p13_t1_g28240 * v2x) + (cse_p13_t1_g28236 * (((-(v1y)) * (-(v2y))) + (v1z * v2z))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 0, i3x4, 0, h_x3_x4);
      double h_x3_y4 = ((cse_p12_t155_g28153 * cse_p12_t159_g28157 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t150_g28148 * cse_p12_t169_g28167) + (cse_p11_t4_invr27997 * cse_p12_t173_g28171 * cse_p12_t213_g28211) + (cse_p12_t150_g28148 * cse_p12_t208_g28206 * t2) + (cse_p12_t173_g28171 * cse_p14_t1_g28248) + (cse_p12_t1_g28225 * ((cse_p11_t1_sqrt27994 * v1z) + (cse_p11_t2_invsqrt27995 * cse_p12_t182_g28180 * v2x))) + (cse_p13_t1_g28236 * (c1z + (v2x * (-(v1y))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 0, i3x4, 1, h_x3_y4);
      double h_x3_z4 = ((cse_p12_t157_g28155 * cse_p12_t159_g28157 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t150_g28148 * cse_p12_t171_g28169) + (cse_p11_t4_invr27997 * cse_p12_t173_g28171 * cse_p12_t214_g28212) + (cse_p12_t150_g28148 * cse_p12_t209_g28207 * t2) + (cse_p12_t173_g28171 * cse_p14_t1_g28246) + (cse_p12_t1_g28225 * ((cse_p11_t1_sqrt27994 * (-(v1y))) + (cse_p11_t2_invsqrt27995 * cse_p12_t183_g28181 * v2x))) + (cse_p13_t1_g28236 * ((v1z * (-(v2x))) + (-(c1y)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 0, i3x4, 2, h_x3_z4);
      double h_y3_y3 = ((cse_p12_t154_g28152 * cse_p12_t154_g28152 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p12_t151_g28149 * cse_p12_t206_g28204 * t2) + (cse_p12_t168_g28166 * ((cse_p11_t4_invr27997 * cse_p12_t151_g28149) + (cse_p11_t4_invr27997 * cse_p12_t220_g28218) + (cse_p12_t206_g28204 * (-(t1))))) + (cse_p12_t1_g28225 * (cse_p12_t1_g28222 + (cse_p11_t2_invsqrt27995 * cse_p12_t181_g28179 * v2y) + (cse_p11_t2_invsqrt27995 * cse_p12_t181_g28179 * v2y) + (-(cse_p12_t69_g28067)))) + (cse_p13_t1_g28236 * ((cse_p12_t189_g28187 * v1x) + (cse_p12_t189_g28187 * v1x) + (cse_p12_t200_g28198 * (-(v1z))) + (cse_p12_t200_g28198 * (-(v1z))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 1, i3x3, 1, h_y3_y3);
      double h_y3_z3 = ((cse_p12_t153_g28151 * cse_p12_t154_g28152 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t151_g28149 * cse_p12_t167_g28165) + (cse_p11_t4_invr27997 * cse_p12_t168_g28166 * cse_p12_t221_g28219) + (cse_p12_t151_g28149 * cse_p12_t207_g28205 * t2) + (cse_p12_t168_g28166 * cse_p14_t1_g28252) + (cse_p12_t1_g28225 * ((cse_p11_t2_invsqrt27995 * cse_p12_t180_g28178 * v2y) + (cse_p11_t2_invsqrt27995 * cse_p12_t181_g28179 * v2z) + (-(cse_p12_t70_g28068)))) + (cse_p13_t1_g28236 * ((cse_p12_t190_g28188 * (-(v1z))) + (cse_p12_t200_g28198 * v1y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 1, i3x3, 2, h_y3_z3);
      double h_y3_x4 = ((cse_p12_t154_g28152 * cse_p12_t160_g28158 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t151_g28149 * cse_p12_t174_g28172) + (cse_p11_t4_invr27997 * cse_p12_t168_g28166 * cse_p12_t215_g28213) + (cse_p12_t151_g28149 * cse_p12_t210_g28208 * t2) + (cse_p12_t168_g28166 * cse_p14_t1_g28251) + (cse_p12_t1_g28225 * ((cse_p11_t1_sqrt27994 * (-(v1z))) + (cse_p11_t2_invsqrt27995 * cse_p12_t185_g28183 * v2y))) + (cse_p13_t1_g28236 * ((v1x * (-(v2y))) + (-(c1z)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 1, i3x4, 0, h_y3_x4);
      double h_y3_y4 = ((cse_p12_t154_g28152 * cse_p12_t155_g28153 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t151_g28149 * cse_p12_t169_g28167) + (cse_p11_t4_invr27997 * cse_p12_t168_g28166 * cse_p12_t213_g28211) + (cse_p12_t151_g28149 * cse_p12_t208_g28206 * t2) + (cse_p12_t168_g28166 * cse_p14_t1_g28248) + (cse_p12_t182_g28180 * cse_p13_t1_g28240 * v2y) + (cse_p13_t1_g28236 * (((-(v1z)) * (-(v2z))) + (v1x * v2x))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 1, i3x4, 1, h_y3_y4);
      double h_y3_z4 = ((cse_p12_t154_g28152 * cse_p12_t157_g28155 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t151_g28149 * cse_p12_t171_g28169) + (cse_p11_t4_invr27997 * cse_p12_t168_g28166 * cse_p12_t214_g28212) + (cse_p12_t151_g28149 * cse_p12_t209_g28207 * t2) + (cse_p12_t168_g28166 * cse_p14_t1_g28246) + (cse_p12_t1_g28225 * ((cse_p11_t1_sqrt27994 * v1x) + (cse_p11_t2_invsqrt27995 * cse_p12_t183_g28181 * v2y))) + (cse_p13_t1_g28236 * (c1x + (v2y * (-(v1z))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 1, i3x4, 2, h_y3_z4);
      double h_z3_z3 = ((cse_p12_t153_g28151 * cse_p12_t153_g28151 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p12_t152_g28150 * cse_p12_t207_g28205 * t2) + (cse_p12_t167_g28165 * (cse_p14_t1_g28252 + (cse_p11_t4_invr27997 * cse_p12_t152_g28150) + (cse_p11_t4_invr27997 * cse_p12_t221_g28219))) + (cse_p12_t1_g28225 * (cse_p12_t1_g28222 + (cse_p11_t2_invsqrt27995 * cse_p12_t180_g28178 * v2z) + (cse_p11_t2_invsqrt27995 * cse_p12_t180_g28178 * v2z) + (-(cse_p12_t71_g28069)))) + (cse_p13_t1_g28236 * ((cse_p12_t190_g28188 * v1y) + (cse_p12_t190_g28188 * v1y) + (cse_p12_t198_g28196 * (-(v1x))) + (cse_p12_t198_g28196 * (-(v1x))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 2, i3x3, 2, h_z3_z3);
      double h_z3_x4 = ((cse_p12_t153_g28151 * cse_p12_t160_g28158 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t152_g28150 * cse_p12_t174_g28172) + (cse_p11_t4_invr27997 * cse_p12_t167_g28165 * cse_p12_t215_g28213) + (cse_p12_t152_g28150 * cse_p12_t210_g28208 * t2) + (cse_p12_t167_g28165 * cse_p14_t1_g28251) + (cse_p12_t1_g28225 * ((cse_p11_t1_sqrt27994 * v1y) + (cse_p11_t2_invsqrt27995 * cse_p12_t185_g28183 * v2z))) + (cse_p13_t1_g28236 * (c1y + (v2z * (-(v1x))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 2, i3x4, 0, h_z3_x4);
      double h_z3_y4 = ((cse_p12_t153_g28151 * cse_p12_t155_g28153 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t152_g28150 * cse_p12_t169_g28167) + (cse_p11_t4_invr27997 * cse_p12_t167_g28165 * cse_p12_t213_g28211) + (cse_p12_t152_g28150 * cse_p12_t208_g28206 * t2) + (cse_p12_t167_g28165 * cse_p14_t1_g28248) + (cse_p12_t1_g28225 * ((cse_p11_t1_sqrt27994 * (-(v1x))) + (cse_p11_t2_invsqrt27995 * cse_p12_t182_g28180 * v2z))) + (cse_p13_t1_g28236 * ((v1y * (-(v2z))) + (-(c1x)))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 2, i3x4, 1, h_z3_y4);
      double h_z3_z4 = ((cse_p12_t153_g28151 * cse_p12_t157_g28155 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t152_g28150 * cse_p12_t171_g28169) + (cse_p11_t4_invr27997 * cse_p12_t167_g28165 * cse_p12_t214_g28212) + (cse_p12_t152_g28150 * cse_p12_t209_g28207 * t2) + (cse_p12_t167_g28165 * cse_p14_t1_g28246) + (cse_p12_t183_g28181 * cse_p13_t1_g28240 * v2z) + (cse_p13_t1_g28236 * (((-(v1x)) * (-(v2x))) + (v1y * v2y))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x3, 2, i3x4, 2, h_z3_z4);
      double h_x4_x4 = ((cse_p12_t160_g28158 * cse_p12_t160_g28158 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p12_t174_g28172 * (cse_p14_t1_g28251 + (cse_p11_t4_invr27997 * cse_p12_t215_g28213) + (cse_p12_t185_g28183 * cse_p13_t1_g28235))) + (cse_p12_t185_g28183 * cse_p12_t1_g28228 * cse_p12_t210_g28208))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x4, 0, i3x4, 0, h_x4_x4);
      double h_x4_y4 = ((cse_p12_t155_g28153 * cse_p12_t160_g28158 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t174_g28172 * cse_p12_t213_g28211) + (cse_p12_t169_g28167 * cse_p12_t185_g28183 * cse_p13_t1_g28235) + (cse_p12_t174_g28172 * cse_p14_t1_g28248) + (cse_p12_t185_g28183 * cse_p12_t1_g28228 * cse_p12_t208_g28206))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x4, 0, i3x4, 1, h_x4_y4);
      double h_x4_z4 = ((cse_p12_t157_g28155 * cse_p12_t160_g28158 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t174_g28172 * cse_p12_t214_g28212) + (cse_p12_t171_g28169 * cse_p12_t185_g28183 * cse_p13_t1_g28235) + (cse_p12_t174_g28172 * cse_p14_t1_g28246) + (cse_p12_t185_g28183 * cse_p12_t1_g28228 * cse_p12_t209_g28207))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x4, 0, i3x4, 2, h_x4_z4);
      double h_y4_y4 = ((cse_p12_t155_g28153 * cse_p12_t155_g28153 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p12_t169_g28167 * (cse_p14_t1_g28248 + (cse_p11_t4_invr27997 * cse_p12_t213_g28211) + (cse_p12_t182_g28180 * cse_p13_t1_g28235))) + (cse_p12_t182_g28180 * cse_p12_t1_g28228 * cse_p12_t208_g28206))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x4, 1, i3x4, 1, h_y4_y4);
      double h_y4_z4 = ((cse_p12_t155_g28153 * cse_p12_t157_g28155 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p11_t4_invr27997 * cse_p12_t169_g28167 * cse_p12_t214_g28212) + (cse_p12_t169_g28167 * cse_p14_t1_g28246) + (cse_p12_t171_g28169 * cse_p12_t182_g28180 * cse_p13_t1_g28235) + (cse_p12_t182_g28180 * cse_p12_t1_g28228 * cse_p12_t209_g28207))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x4, 1, i3x4, 2, h_y4_z4);
      double h_z4_z4 = ((cse_p12_t157_g28155 * cse_p12_t157_g28155 * cse_p12_t1_g28226) + (cse_p12_t1_g28227 * ((cse_p12_t171_g28169 * (cse_p14_t1_g28246 + (cse_p11_t4_invr27997 * cse_p12_t214_g28212) + (cse_p12_t183_g28181 * cse_p13_t1_g28235))) + (cse_p12_t183_g28181 * cse_p12_t1_g28228 * cse_p12_t209_g28207))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  i3x4, 2, i3x4, 2, h_z4_z4);
    }
  }
  return energy_added;
}
void hessian_fd(double v, double n, double sinphase, double cosphase, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 12;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  const double invh2 = 1.0/((h*h));
  double e0 = 0.0;
  energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e0, 0, 0, 0, 0);
  if (energy_accumulate) { *energy_accumulate += e0; }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 0, d);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 1, d);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x1, 2, d);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 0, d);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 1, d);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x2, 2, d);
  }
  {
    double saved = position[i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 0] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 0, d);
  }
  {
    double saved = position[i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 1] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 1, d);
  }
  {
    double saved = position[i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 2] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x3, 2, d);
  }
  {
    double saved = position[i3x4 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 0] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 0, d);
  }
  {
    double saved = position[i3x4 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 1] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 1, d);
  }
  {
    double saved = position[i3x4 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 2] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved;
    double d = (e_plus - e_minus) * inv2h;
    KernelGradientAcc(i3x4, 2, d);
  }
  {
    double saved = position[i3x1 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 0] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 0, i3x1, 0, hval);
  }
  {
    double saved = position[i3x1 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 1] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 1, i3x1, 1, hval);
  }
  {
    double saved = position[i3x1 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x1 + 2] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x1 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x1, 2, i3x1, 2, hval);
  }
  {
    double saved = position[i3x2 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 0] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 0, i3x2, 0, hval);
  }
  {
    double saved = position[i3x2 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 1] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 1, i3x2, 1, hval);
  }
  {
    double saved = position[i3x2 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x2 + 2] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x2 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x2, 2, hval);
  }
  {
    double saved = position[i3x3 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 0] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 0, i3x3, 0, hval);
  }
  {
    double saved = position[i3x3 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 1] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 1, i3x3, 1, hval);
  }
  {
    double saved = position[i3x3 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x3 + 2] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x3 + 2] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x3, 2, i3x3, 2, hval);
  }
  {
    double saved = position[i3x4 + 0];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 0] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 0] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 0, i3x4, 0, hval);
  }
  {
    double saved = position[i3x4 + 1];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 1] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
    position[i3x4 + 1] = saved;
    double hval = (e_plus + e_minus - (2.0*e0)) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 1, i3x4, 1, hval);
  }
  {
    double saved = position[i3x4 + 2];
    double e_plus = 0.0;
    double e_minus = 0.0;
    position[i3x4 + 2] = saved + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_plus, 0, 0, 0, 0);
    position[i3x4 + 2] = saved - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_minus, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_i - h; position[i3x3 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x3 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_i - h; position[i3x3 + 2] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x3 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x3 + 2] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_i - h; position[i3x4 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x1 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x1 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x1 + 2] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x1 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x2 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x2 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x2 + 2] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x2 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x3 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x3 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x3 + 2] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x3 + 2] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x4 + 0] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x4 + 0] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
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
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pp, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_pm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i - h; position[i3x4 + 1] = saved_j + h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mp, 0, 0, 0, 0);
    position[i3x4 + 1] = saved_j - h;
    energy(v, n, sinphase, cosphase, i3x1, i3x2, i3x3, i3x4, position, &e_mm, 0, 0, 0, 0);
    position[i3x4 + 2] = saved_i; position[i3x4 + 1] = saved_j;
    double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, i3x4, 2, i3x4, 1, hval);
  }
}

};
