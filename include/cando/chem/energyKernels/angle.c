#include "angle.h"

template <typename HESSIAN>
struct Angle {
  static constexpr size_t PositionSize = 9;
  static constexpr size_t EnergyAccumulatorSize = 1;
  static std::string description() { return "mathkernel-angle"; };
double energy(const angle_term& term, double* position, double* energy_accumulate) {
  double energy_added = 0.0;
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
    {
      /* !BASE */
      double cse_p11_t1_g1 = (-(x2));
      double cse_p11_t2_g2 = (-(y2));
      double cse_p11_t3_g3 = (-(z2));
      double vx1 = (cse_p11_t1_g1 + x1);
      double vy1 = (cse_p11_t2_g2 + y1);
      double vz1 = (cse_p11_t3_g3 + z1);
      double vx2 = (cse_p11_t1_g1 + x3);
      double vy2 = (cse_p11_t2_g2 + y3);
      double vz2 = (cse_p11_t3_g3 + z3);
      double dot = ((vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2));
      double n1_sq = ((vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1));
      double n2_sq = ((vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2));
      double n1 = sqrt(n1_sq);
      double n2 = sqrt(n2_sq);
      double cos_theta = (dot * (1.0 / ((n1 * n2))));
      if (fabs(cos_theta)>(1.0-VERYSMALL)) cos_theta = linear_angle_clamp(cos_theta);
      double sin_theta = sqrt((1.0 + (-(((cos_theta) * (cos_theta))))));
      double theta = acos(cos_theta);
      double dtheta = (theta + (-(term.t0)));
      double energy = (dtheta * dtheta * term.kt);
      energy_added += energy;
      *energy_accumulate += energy;
    }
  }
  return energy_added;
}
void energy_fd(const angle_term& term, double* position, double* energy_accumulate)
{
  energy(term, position, energy_accumulate);
}

double gradient(const angle_term& term, double* position, double* energy_accumulate, double* force) {
  double energy_added = 0.0;
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
    {
      /* !BASE */
      double cse_p11_t1_g9 = (-(x2));
      double cse_p11_t2_g10 = (-(y2));
      double cse_p11_t3_g11 = (-(z2));
      double vx1 = (cse_p11_t1_g9 + x1);
      double vy1 = (cse_p11_t2_g10 + y1);
      double vz1 = (cse_p11_t3_g11 + z1);
      double vx2 = (cse_p11_t1_g9 + x3);
      double vy2 = (cse_p11_t2_g10 + y3);
      double vz2 = (cse_p11_t3_g11 + z3);
      double dot = ((vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2));
      double n1_sq = ((vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1));
      double n2_sq = ((vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2));
      double cse_p1_t1_sqrt1 = sqrt(n1_sq);
      double cse_p1_t2_invsqrt2 = (1.0 / (cse_p1_t1_sqrt1));
      double n1 = cse_p1_t1_sqrt1;
      double cse_p1_t3_sqrt3 = sqrt(n2_sq);
      double cse_p1_t4_invsqrt4 = (1.0 / (cse_p1_t3_sqrt3));
      double n2 = cse_p1_t3_sqrt3;
      double cse_p1_t5_invr5 = (1.0 / ((n1 * n2)));
      double cse_p1_t6_invr26 = (cse_p1_t5_invr5 * cse_p1_t5_invr5);
      double cse_p113_t1_g20 = (cse_p1_t6_invr26 * dot);
      double cos_theta = (cse_p1_t5_invr5 * dot);
      if (fabs(cos_theta)>(1.0-VERYSMALL)) cos_theta = linear_angle_clamp(cos_theta);
      double cse_p1_t7_sqrt7 = sqrt((1.0 + (-(((cos_theta) * (cos_theta))))));
      double cse_p1_t8_invsqrt8 = (1.0 / (cse_p1_t7_sqrt7));
      double sin_theta = cse_p1_t7_sqrt7;
      double theta = acos(cos_theta);
      double dtheta = (theta + (-(term.t0)));
      double cse_p11_t1_g18 = (dtheta * term.kt);
      double energy = (cse_p11_t1_g18 * dtheta);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p11_t1_g12 = (2.0 * cse_p11_t1_g18);
      double cse_p11_t1_g13 = (cse_p1_t2_invsqrt2 * n2);
      double cse_p11_t1_g16 = (cse_p113_t1_g20 * cse_p11_t1_g13);
      double cse_p62_t1_g19 = (cse_p11_t1_g12 * cse_p1_t8_invsqrt8);
      double g_x1 = (-((cse_p62_t1_g19 * ((cse_p1_t5_invr5 * vx2) + (-((cse_p11_t1_g16 * vx1)))))));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = (-((cse_p62_t1_g19 * ((cse_p1_t5_invr5 * vy2) + (-((cse_p11_t1_g16 * vy1)))))));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = (-((cse_p62_t1_g19 * ((cse_p1_t5_invr5 * vz2) + (-((cse_p11_t1_g16 * vz1)))))));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p11_t1_g14 = (cse_p1_t4_invsqrt4 * n1);
      double g_x2 = (-((cse_p62_t1_g19 * ((cse_p1_t5_invr5 * ((-(vx1)) + (-(vx2)))) + (-((cse_p113_t1_g20 * ((-((cse_p11_t1_g13 * vx1))) + (-((cse_p11_t1_g14 * vx2)))))))))));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (-((cse_p62_t1_g19 * ((cse_p1_t5_invr5 * ((-(vy1)) + (-(vy2)))) + (-((cse_p113_t1_g20 * ((-((cse_p11_t1_g13 * vy1))) + (-((cse_p11_t1_g14 * vy2)))))))))));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (-((cse_p62_t1_g19 * ((cse_p1_t5_invr5 * ((-(vz1)) + (-(vz2)))) + (-((cse_p113_t1_g20 * ((-((cse_p11_t1_g13 * vz1))) + (-((cse_p11_t1_g14 * vz2)))))))))));
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double cse_p11_t1_g17 = (cse_p113_t1_g20 * cse_p11_t1_g14);
      double g_x3 = (-((cse_p62_t1_g19 * ((cse_p1_t5_invr5 * vx1) + (-((cse_p11_t1_g17 * vx2)))))));
      KernelGradientAcc(term.i3x3, 0, g_x3);
      double g_y3 = (-((cse_p62_t1_g19 * ((cse_p1_t5_invr5 * vy1) + (-((cse_p11_t1_g17 * vy2)))))));
      KernelGradientAcc(term.i3x3, 1, g_y3);
      double g_z3 = (-((cse_p62_t1_g19 * ((cse_p1_t5_invr5 * vz1) + (-((cse_p11_t1_g17 * vz2)))))));
      KernelGradientAcc(term.i3x3, 2, g_z3);
    }
  }
  return energy_added;
}
void gradient_fd(const angle_term& term, double* position, double* energy_accumulate, double* force)
{
  constexpr size_t PositionSize = 9;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  double e0[EnergyAccumulatorSize] = {0.0};
  energy(term, position, e0);
  if (energy_accumulate) { *energy_accumulate += e0[0]; }
  {
    double saved = position[term.i3x1 + 0];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x1 + 0] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x1 + 0] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x1 + 0] = saved;
    double d = (e_plus[0] - e_minus[0]) * inv2h;
    KernelGradientAcc(term.i3x1, 0, d);
  }
  {
    double saved = position[term.i3x1 + 1];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x1 + 1] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x1 + 1] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x1 + 1] = saved;
    double d = (e_plus[0] - e_minus[0]) * inv2h;
    KernelGradientAcc(term.i3x1, 1, d);
  }
  {
    double saved = position[term.i3x1 + 2];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x1 + 2] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x1 + 2] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x1 + 2] = saved;
    double d = (e_plus[0] - e_minus[0]) * inv2h;
    KernelGradientAcc(term.i3x1, 2, d);
  }
  {
    double saved = position[term.i3x2 + 0];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 0] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x2 + 0] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x2 + 0] = saved;
    double d = (e_plus[0] - e_minus[0]) * inv2h;
    KernelGradientAcc(term.i3x2, 0, d);
  }
  {
    double saved = position[term.i3x2 + 1];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 1] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x2 + 1] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x2 + 1] = saved;
    double d = (e_plus[0] - e_minus[0]) * inv2h;
    KernelGradientAcc(term.i3x2, 1, d);
  }
  {
    double saved = position[term.i3x2 + 2];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 2] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x2 + 2] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x2 + 2] = saved;
    double d = (e_plus[0] - e_minus[0]) * inv2h;
    KernelGradientAcc(term.i3x2, 2, d);
  }
  {
    double saved = position[term.i3x3 + 0];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 0] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x3 + 0] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x3 + 0] = saved;
    double d = (e_plus[0] - e_minus[0]) * inv2h;
    KernelGradientAcc(term.i3x3, 0, d);
  }
  {
    double saved = position[term.i3x3 + 1];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 1] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x3 + 1] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x3 + 1] = saved;
    double d = (e_plus[0] - e_minus[0]) * inv2h;
    KernelGradientAcc(term.i3x3, 1, d);
  }
  {
    double saved = position[term.i3x3 + 2];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 2] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x3 + 2] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x3 + 2] = saved;
    double d = (e_plus[0] - e_minus[0]) * inv2h;
    KernelGradientAcc(term.i3x3, 2, d);
  }
}

double hessian(const angle_term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
  double energy_added = 0.0;
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
    {
      /* !BASE */
      double cse_p11_t186_g194 = (-(x2));
      double cse_p11_t187_g195 = (-(y2));
      double cse_p11_t188_g196 = (-(z2));
      double vx1 = (cse_p11_t186_g194 + x1);
      double cse_p11_t180_g188 = (-(vx1));
      double vy1 = (cse_p11_t187_g195 + y1);
      double cse_p11_t182_g190 = (-(vy1));
      double vz1 = (cse_p11_t188_g196 + z1);
      double cse_p11_t184_g192 = (-(vz1));
      double vx2 = (cse_p11_t186_g194 + x3);
      double cse_p11_t138_g146 = (cse_p11_t180_g188 + (-(vx2)));
      double cse_p11_t181_g189 = (-(vx2));
      double vy2 = (cse_p11_t187_g195 + y3);
      double cse_p11_t139_g147 = (cse_p11_t182_g190 + (-(vy2)));
      double cse_p11_t183_g191 = (-(vy2));
      double vz2 = (cse_p11_t188_g196 + z3);
      double cse_p11_t140_g148 = (cse_p11_t184_g192 + (-(vz2)));
      double cse_p11_t185_g193 = (-(vz2));
      double cse_p167_t1_g304 = (vx1 * vx2);
      double cse_p167_t1_g308 = (vy1 * vy2);
      double cse_p167_t1_g311 = (vz1 * vz2);
      double dot = (cse_p167_t1_g304 + cse_p167_t1_g308 + cse_p167_t1_g311);
      double n1_sq = ((vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1));
      double cse_p579_t1_sqrt398 = sqrt(n1_sq);
      double cse_p579_t2_invr399 = (1.0 / (n1_sq));
      double cse_p322_t1_g335 = (cse_p579_t2_invr399 * cse_p579_t2_invr399);
      double cse_p114_t1_g263 = (-2.0 * cse_p322_t1_g335);
      double cse_p63_t13_g236 = (cse_p114_t1_g263 + cse_p322_t1_g335);
      double cse_p167_t1_g305 = (vx2 * vx2);
      double cse_p167_t1_g309 = (vy2 * vy2);
      double cse_p167_t1_g312 = (vz2 * vz2);
      double n2_sq = (cse_p167_t1_g305 + cse_p167_t1_g309 + cse_p167_t1_g312);
      double cse_p579_t3_sqrt400 = sqrt(n2_sq);
      double cse_p579_t4_invr401 = (1.0 / (n2_sq));
      double cse_p61_t6_invr2216 = (cse_p579_t4_invr401 * cse_p579_t4_invr401);
      double cse_p323_t1_g336 = (cse_p579_t3_sqrt400 * cse_p579_t4_invr401);
      double cse_p11_t115_g123 = (cse_p167_t1_g304 * cse_p323_t1_g336);
      double cse_p63_t1_g244 = (cse_p323_t1_g336 * cse_p579_t2_invr399);
      double cse_p11_t116_g124 = (cse_p63_t1_g244 * vx2);
      double cse_p11_t117_g125 = (cse_p167_t1_g308 * cse_p323_t1_g336);
      double cse_p11_t118_g126 = (cse_p63_t1_g244 * vy2);
      double cse_p11_t119_g127 = (cse_p167_t1_g311 * cse_p323_t1_g336);
      double cse_p11_t120_g128 = (cse_p63_t1_g244 * vz2);
      double cse_p1_t2_invsqrt2 = (1.0 / (cse_p579_t1_sqrt398));
      double n1 = cse_p579_t1_sqrt398;
      double cse_p63_t1_g237 = (cse_p579_t3_sqrt400 * n1);
      double cse_p478_t1_g379 = (cse_p579_t4_invr401 * cse_p63_t1_g237);
      double cse_p64_t1_g252 = (cse_p61_t6_invr2216 * cse_p63_t1_g237);
      double cse_p116_t1_g282 = (2.0 * cse_p64_t1_g252);
      double cse_p167_t1_g306 = (vx2 * vy2);
      double cse_p11_t20_g28 = (cse_p116_t1_g282 * cse_p167_t1_g306);
      double cse_p167_t1_g307 = (vx2 * vz2);
      double cse_p11_t21_g29 = (cse_p116_t1_g282 * cse_p167_t1_g307);
      double cse_p167_t1_g310 = (vy2 * vz2);
      double cse_p11_t22_g30 = (cse_p116_t1_g282 * cse_p167_t1_g310);
      double cse_p11_t91_g99 = (cse_p478_t1_g379 * vx2);
      double cse_p116_t1_g276 = (cse_p478_t1_g379 * cse_p579_t4_invr401);
      double cse_p11_t92_g100 = (cse_p116_t1_g276 * cse_p167_t1_g305);
      double cse_p11_t94_g102 = (cse_p116_t1_g276 * cse_p167_t1_g306);
      double cse_p11_t96_g104 = (cse_p116_t1_g276 * cse_p167_t1_g307);
      double cse_p11_t97_g105 = (cse_p478_t1_g379 * vy2);
      double cse_p11_t98_g106 = (cse_p116_t1_g276 * cse_p167_t1_g309);
      double cse_p11_t100_g108 = (cse_p116_t1_g276 * cse_p167_t1_g310);
      double cse_p11_t101_g109 = (cse_p478_t1_g379 * vz2);
      double cse_p11_t102_g110 = (cse_p116_t1_g276 * cse_p167_t1_g312);
      double cse_p11_t162_g170 = (-(cse_p11_t91_g99));
      double cse_p11_t164_g172 = (-(cse_p11_t94_g102));
      double cse_p11_t166_g174 = (-(cse_p11_t96_g104));
      double cse_p11_t167_g175 = (-(cse_p11_t97_g105));
      double cse_p11_t169_g177 = (-(cse_p11_t100_g108));
      double cse_p11_t170_g178 = (-(cse_p11_t101_g109));
      double cse_p11_t177_g185 = (-(cse_p478_t1_g379));
      double cse_p1_t4_invsqrt4 = (1.0 / (cse_p579_t3_sqrt400));
      double n2 = cse_p579_t3_sqrt400;
      double cse_p116_t1_g283 = (cse_p579_t1_sqrt398 * cse_p579_t2_invr399);
      double cse_p11_t1_g207 = (cse_p116_t1_g283 * n2);
      double cse_p426_t1_g371 = (cse_p11_t1_g207 * vz1);
      double cse_p270_t1_g328 = (cse_p11_t1_g207 * vy1);
      double cse_p166_t1_g303 = (cse_p11_t1_g207 * vx1);
      double cse_p61_t7_invr217 = (1.0 / ((n1 * n2)));
      double cse_p61_t8_invr2218 = (cse_p61_t7_invr217 * cse_p61_t7_invr217);
      double cse_p374_t1_g343 = (cse_p61_t8_invr2218 * dot);
      double cse_p61_t9_invr3219 = (cse_p61_t7_invr217 * cse_p61_t8_invr2218);
      double cse_p63_t1_g241 = (2.0 * cse_p61_t9_invr3219);
      double cse_p63_t1_g243 = (cse_p478_t1_g379 * cse_p63_t1_g241);
      double cse_p63_t5_g228 = (cse_p579_t4_invr401 * cse_p63_t1_g243);
      double cse_p116_t1_g281 = (2.0 * cse_p61_t8_invr2218);
      double cse_p63_t12_g235 = (cse_p63_t5_g228 + (cse_p116_t1_g281 * cse_p61_t6_invr2216));
      double cse_p11_t1_g9 = (cse_p61_t8_invr2218 * (cse_p11_t162_g170 + (-(cse_p166_t1_g303))));
      double cse_p11_t2_g10 = (cse_p61_t8_invr2218 * (cse_p11_t167_g175 + (-(cse_p270_t1_g328))));
      double cse_p11_t3_g11 = (cse_p61_t8_invr2218 * (cse_p11_t170_g178 + (-(cse_p426_t1_g371))));
      double cse_p115_t1_g273 = (n2 * vx1);
      double cse_p11_t13_g21 = (cse_p114_t1_g263 * cse_p115_t1_g273);
      double cse_p116_t1_g279 = (n2 * vy1);
      double cse_p11_t14_g22 = (cse_p114_t1_g263 * cse_p116_t1_g279);
      double cse_p117_t1_g293 = (n2 * vz1);
      double cse_p11_t15_g23 = (cse_p114_t1_g263 * cse_p117_t1_g293);
      double cse_p11_t16_g24 = (cse_p116_t1_g281 * cse_p322_t1_g335);
      double cse_p11_t17_g25 = (cse_p116_t1_g281 * cse_p61_t6_invr2216);
      double cse_p11_t19_g27 = (cse_p116_t1_g276 * cse_p63_t1_g241);
      double cse_p117_t1_g294 = (cse_p579_t2_invr399 * cse_p63_t1_g241);
      double cse_p11_t23_g31 = (cse_p117_t1_g294 * cse_p11_t1_g207);
      double cse_p11_t33_g41 = (dot * (-(cse_p11_t1_g9)));
      double cse_p11_t34_g42 = (dot * (-(cse_p11_t2_g10)));
      double cse_p11_t35_g43 = (dot * (-(cse_p11_t3_g11)));
      double cse_p218_t1_g321 = (cse_p374_t1_g343 * cse_p478_t1_g379);
      double cse_p11_t37_g45 = (cse_p218_t1_g321 * vx2);
      double cse_p11_t38_g46 = (cse_p218_t1_g321 * vy2);
      double cse_p11_t39_g47 = (cse_p218_t1_g321 * vz2);
      double cse_p11_t40_g48 = (cse_p11_t1_g207 * cse_p374_t1_g343);
      double cse_p11_t41_g49 = (cse_p166_t1_g303 * cse_p374_t1_g343);
      double cse_p11_t42_g50 = (cse_p270_t1_g328 * cse_p374_t1_g343);
      double cse_p11_t43_g51 = (cse_p374_t1_g343 * cse_p426_t1_g371);
      double cse_p11_t90_g98 = (n1 * n2);
      double cse_p377_t1_g364 = (cse_p478_t1_g379 * cse_p61_t8_invr2218);
      double cse_p11_t93_g101 = (cse_p167_t1_g306 * cse_p377_t1_g364);
      double cse_p11_t95_g103 = (cse_p167_t1_g307 * cse_p377_t1_g364);
      double cse_p11_t99_g107 = (cse_p167_t1_g310 * cse_p377_t1_g364);
      double cse_p11_t103_g111 = (cse_p115_t1_g273 * cse_p322_t1_g335);
      double cse_p11_t105_g113 = (cse_p115_t1_g273 * cse_p579_t2_invr399 * vx1);
      double cse_p11_t106_g114 = (cse_p166_t1_g303 * cse_p61_t8_invr2218 * vy1);
      double cse_p11_t107_g115 = (cse_p166_t1_g303 * cse_p61_t8_invr2218 * vz1);
      double cse_p11_t108_g116 = (cse_p116_t1_g279 * cse_p322_t1_g335);
      double cse_p11_t110_g118 = (cse_p116_t1_g279 * cse_p579_t2_invr399 * vy1);
      double cse_p11_t111_g119 = (cse_p270_t1_g328 * cse_p61_t8_invr2218 * vz1);
      double cse_p11_t112_g120 = (cse_p117_t1_g293 * cse_p322_t1_g335);
      double cse_p11_t114_g122 = (cse_p117_t1_g293 * cse_p579_t2_invr399 * vz1);
      double cse_p11_t121_g129 = (cse_p11_t162_g170 + (-(cse_p166_t1_g303)));
      double cse_p63_t1_g245 = (cse_p117_t1_g294 * cse_p11_t121_g129);
      double cse_p63_t2_g225 = (cse_p63_t1_g245 * n2);
      double cse_p11_t122_g130 = (cse_p11_t167_g175 + (-(cse_p270_t1_g328)));
      double cse_p63_t1_g246 = (cse_p117_t1_g294 * cse_p11_t122_g130);
      double cse_p63_t3_g226 = (cse_p63_t1_g246 * n2);
      double cse_p11_t123_g131 = (cse_p11_t170_g178 + (-(cse_p426_t1_g371)));
      double cse_p63_t1_g247 = (cse_p117_t1_g294 * cse_p11_t123_g131);
      double cse_p63_t4_g227 = (cse_p63_t1_g247 * n2);
      double cse_p11_t124_g132 = (cse_p11_t103_g111 + cse_p11_t116_g124 + cse_p11_t13_g21);
      double cse_p63_t6_g229 = (cse_p11_t124_g132 * cse_p61_t8_invr2218);
      double cse_p63_t9_g232 = (cse_p63_t2_g225 + cse_p63_t6_g229);
      double cse_p11_t125_g133 = (cse_p11_t108_g116 + cse_p11_t118_g126 + cse_p11_t14_g22);
      double cse_p63_t7_g230 = (cse_p11_t125_g133 * cse_p61_t8_invr2218);
      double cse_p63_t10_g233 = (cse_p63_t3_g226 + cse_p63_t7_g230);
      double cse_p11_t126_g134 = (cse_p11_t112_g120 + cse_p11_t120_g128 + cse_p11_t15_g23);
      double cse_p63_t8_g231 = (cse_p11_t126_g134 * cse_p61_t8_invr2218);
      double cse_p63_t11_g234 = (cse_p63_t4_g227 + cse_p63_t8_g231);
      double cse_p11_t127_g135 = (cse_p11_t16_g24 + cse_p11_t23_g31);
      double cse_p11_t128_g136 = (cse_p11_t17_g25 + cse_p11_t19_g27);
      double cse_p11_t142_g150 = (-((cse_p11_t121_g129 * cse_p61_t8_invr2218)));
      double cse_p11_t143_g151 = (-((cse_p11_t122_g130 * cse_p61_t8_invr2218)));
      double cse_p11_t144_g152 = (-((cse_p11_t123_g131 * cse_p61_t8_invr2218)));
      double cse_p11_t154_g162 = (-(cse_p218_t1_g321));
      double cse_p11_t155_g163 = (-(cse_p11_t37_g45));
      double cse_p11_t156_g164 = (-(cse_p11_t38_g46));
      double cse_p11_t157_g165 = (-(cse_p11_t39_g47));
      double cse_p11_t158_g166 = (-(cse_p11_t40_g48));
      double cse_p11_t159_g167 = (-(cse_p11_t41_g49));
      double cse_p11_t160_g168 = (-(cse_p11_t42_g50));
      double cse_p11_t161_g169 = (-(cse_p11_t43_g51));
      double cse_p11_t163_g171 = (-(cse_p11_t93_g101));
      double cse_p11_t165_g173 = (-(cse_p11_t95_g103));
      double cse_p11_t168_g176 = (-(cse_p11_t99_g107));
      double cse_p11_t171_g179 = (-(cse_p166_t1_g303));
      double cse_p11_t172_g180 = (-(cse_p11_t106_g114));
      double cse_p11_t173_g181 = (-(cse_p11_t107_g115));
      double cse_p11_t174_g182 = (-(cse_p270_t1_g328));
      double cse_p11_t175_g183 = (-(cse_p11_t111_g119));
      double cse_p11_t176_g184 = (-(cse_p426_t1_g371));
      double cse_p1_t5_invr5 = (1.0 / (cse_p11_t90_g98));
      double cse_p11_t18_g26 = (2.0 * cse_p1_t5_invr5);
      double cse_p11_t24_g32 = (cse_p1_t5_invr5 * (cse_p11_t180_g188 + cse_p11_t181_g189));
      double cse_p11_t25_g33 = (cse_p1_t5_invr5 * (cse_p11_t182_g190 + cse_p11_t183_g191));
      double cse_p11_t26_g34 = (cse_p1_t5_invr5 * (cse_p11_t184_g192 + cse_p11_t185_g193));
      double cse_p11_t27_g35 = (cse_p1_t5_invr5 * vx1);
      double cse_p11_t28_g36 = (cse_p1_t5_invr5 * vx2);
      double cse_p11_t29_g37 = (cse_p1_t5_invr5 * vy1);
      double cse_p11_t30_g38 = (cse_p1_t5_invr5 * vy2);
      double cse_p11_t31_g39 = (cse_p1_t5_invr5 * vz1);
      double cse_p11_t32_g40 = (cse_p1_t5_invr5 * vz2);
      double cse_p11_t129_g137 = (cse_p11_t24_g32 + (-((cse_p374_t1_g343 * (cse_p11_t162_g170 + cse_p11_t171_g179)))));
      double cse_p11_t130_g138 = (cse_p11_t25_g33 + (-((cse_p374_t1_g343 * (cse_p11_t167_g175 + cse_p11_t174_g182)))));
      double cse_p11_t131_g139 = (cse_p11_t26_g34 + (-((cse_p374_t1_g343 * (cse_p11_t170_g178 + cse_p11_t176_g184)))));
      double cse_p11_t132_g140 = (cse_p11_t155_g163 + cse_p11_t27_g35);
      double cse_p11_t133_g141 = (cse_p11_t159_g167 + cse_p11_t28_g36);
      double cse_p11_t134_g142 = (cse_p11_t156_g164 + cse_p11_t29_g37);
      double cse_p11_t135_g143 = (cse_p11_t160_g168 + cse_p11_t30_g38);
      double cse_p11_t136_g144 = (cse_p11_t157_g165 + cse_p11_t31_g39);
      double cse_p11_t137_g145 = (cse_p11_t161_g169 + cse_p11_t32_g40);
      double cse_p11_t179_g187 = (-(cse_p1_t5_invr5));
      double cse_p1_t6_invr26 = (cse_p1_t5_invr5 * cse_p1_t5_invr5);
      double cos_theta = (cse_p1_t5_invr5 * dot);
      double cse_p61_t10_invsqrt220 = (1.0 / sqrt((1.0 + (-(((cos_theta) * (cos_theta)))))));
      double cse_p11_t4_g12 = (cse_p11_t129_g137 * cse_p61_t10_invsqrt220);
      double cse_p11_t5_g13 = (cse_p11_t130_g138 * cse_p61_t10_invsqrt220);
      double cse_p11_t6_g14 = (cse_p11_t131_g139 * cse_p61_t10_invsqrt220);
      double cse_p11_t7_g15 = (cse_p11_t132_g140 * cse_p61_t10_invsqrt220);
      double cse_p11_t8_g16 = (cse_p11_t133_g141 * cse_p61_t10_invsqrt220);
      double cse_p11_t9_g17 = (cse_p11_t134_g142 * cse_p61_t10_invsqrt220);
      double cse_p11_t10_g18 = (cse_p11_t135_g143 * cse_p61_t10_invsqrt220);
      double cse_p11_t11_g19 = (cse_p11_t136_g144 * cse_p61_t10_invsqrt220);
      double cse_p11_t12_g20 = (cse_p11_t137_g145 * cse_p61_t10_invsqrt220);
      double cse_p63_t1_g238 = (cse_p61_t10_invsqrt220 * cse_p61_t10_invsqrt220 * term.kt);
      double cse_p115_t1_g266 = (cse_p11_t129_g137 * cse_p63_t1_g238);
      double cse_p11_t44_g52 = (cse_p115_t1_g266 * cse_p11_t129_g137);
      double cse_p11_t45_g53 = (cse_p115_t1_g266 * cse_p11_t130_g138);
      double cse_p11_t46_g54 = (cse_p115_t1_g266 * cse_p11_t131_g139);
      double cse_p11_t47_g55 = (cse_p115_t1_g266 * cse_p11_t132_g140);
      double cse_p11_t48_g56 = (cse_p115_t1_g266 * cse_p11_t133_g141);
      double cse_p11_t49_g57 = (cse_p115_t1_g266 * cse_p11_t134_g142);
      double cse_p11_t50_g58 = (cse_p115_t1_g266 * cse_p11_t135_g143);
      double cse_p11_t51_g59 = (cse_p115_t1_g266 * cse_p11_t136_g144);
      double cse_p11_t52_g60 = (cse_p115_t1_g266 * cse_p11_t137_g145);
      double cse_p115_t1_g270 = (cse_p11_t130_g138 * cse_p63_t1_g238);
      double cse_p11_t53_g61 = (cse_p115_t1_g270 * cse_p11_t130_g138);
      double cse_p11_t54_g62 = (cse_p115_t1_g270 * cse_p11_t131_g139);
      double cse_p11_t55_g63 = (cse_p115_t1_g270 * cse_p11_t132_g140);
      double cse_p11_t56_g64 = (cse_p115_t1_g270 * cse_p11_t133_g141);
      double cse_p11_t57_g65 = (cse_p115_t1_g270 * cse_p11_t134_g142);
      double cse_p11_t58_g66 = (cse_p115_t1_g270 * cse_p11_t135_g143);
      double cse_p11_t59_g67 = (cse_p115_t1_g270 * cse_p11_t136_g144);
      double cse_p11_t60_g68 = (cse_p115_t1_g270 * cse_p11_t137_g145);
      double cse_p116_t1_g277 = (cse_p11_t131_g139 * cse_p63_t1_g238);
      double cse_p11_t61_g69 = (cse_p116_t1_g277 * cse_p11_t131_g139);
      double cse_p11_t62_g70 = (cse_p116_t1_g277 * cse_p11_t132_g140);
      double cse_p11_t63_g71 = (cse_p116_t1_g277 * cse_p11_t133_g141);
      double cse_p11_t64_g72 = (cse_p116_t1_g277 * cse_p11_t134_g142);
      double cse_p11_t65_g73 = (cse_p116_t1_g277 * cse_p11_t135_g143);
      double cse_p11_t66_g74 = (cse_p116_t1_g277 * cse_p11_t136_g144);
      double cse_p11_t67_g75 = (cse_p116_t1_g277 * cse_p11_t137_g145);
      double cse_p116_t1_g284 = (cse_p11_t132_g140 * cse_p63_t1_g238);
      double cse_p11_t68_g76 = (cse_p116_t1_g284 * cse_p11_t132_g140);
      double cse_p11_t69_g77 = (cse_p116_t1_g284 * cse_p11_t133_g141);
      double cse_p11_t70_g78 = (cse_p116_t1_g284 * cse_p11_t134_g142);
      double cse_p11_t71_g79 = (cse_p116_t1_g284 * cse_p11_t135_g143);
      double cse_p11_t72_g80 = (cse_p116_t1_g284 * cse_p11_t136_g144);
      double cse_p11_t73_g81 = (cse_p116_t1_g284 * cse_p11_t137_g145);
      double cse_p117_t1_g295 = (cse_p11_t133_g141 * cse_p63_t1_g238);
      double cse_p11_t74_g82 = (cse_p117_t1_g295 * cse_p11_t133_g141);
      double cse_p11_t75_g83 = (cse_p117_t1_g295 * cse_p11_t134_g142);
      double cse_p11_t76_g84 = (cse_p117_t1_g295 * cse_p11_t135_g143);
      double cse_p11_t77_g85 = (cse_p117_t1_g295 * cse_p11_t136_g144);
      double cse_p11_t78_g86 = (cse_p117_t1_g295 * cse_p11_t137_g145);
      double cse_p167_t1_g313 = (cse_p11_t134_g142 * cse_p63_t1_g238);
      double cse_p11_t79_g87 = (cse_p11_t134_g142 * cse_p167_t1_g313);
      double cse_p11_t80_g88 = (cse_p11_t135_g143 * cse_p167_t1_g313);
      double cse_p11_t81_g89 = (cse_p11_t136_g144 * cse_p167_t1_g313);
      double cse_p11_t82_g90 = (cse_p11_t137_g145 * cse_p167_t1_g313);
      double cse_p376_t1_g362 = (cse_p11_t135_g143 * cse_p63_t1_g238);
      double cse_p11_t83_g91 = (cse_p11_t135_g143 * cse_p376_t1_g362);
      double cse_p11_t84_g92 = (cse_p11_t136_g144 * cse_p376_t1_g362);
      double cse_p11_t85_g93 = (cse_p11_t137_g145 * cse_p376_t1_g362);
      double cse_p11_t86_g94 = (cse_p11_t136_g144 * cse_p11_t136_g144 * cse_p63_t1_g238);
      double cse_p11_t87_g95 = (cse_p11_t136_g144 * cse_p11_t137_g145 * cse_p63_t1_g238);
      double cse_p11_t88_g96 = (cse_p11_t137_g145 * cse_p11_t137_g145 * cse_p63_t1_g238);
      double cse_p11_t141_g149 = (1.0 + (-(((cos_theta) * (cos_theta)))));
      double cse_p11_t145_g153 = (-(cse_p11_t4_g12));
      double cse_p11_t146_g154 = (-(cse_p11_t5_g13));
      double cse_p11_t147_g155 = (-(cse_p11_t6_g14));
      double cse_p11_t148_g156 = (-(cse_p11_t7_g15));
      double cse_p11_t149_g157 = (-(cse_p11_t8_g16));
      double cse_p11_t150_g158 = (-(cse_p11_t9_g17));
      double cse_p11_t151_g159 = (-(cse_p11_t10_g18));
      double cse_p11_t152_g160 = (-(cse_p11_t11_g19));
      double cse_p11_t153_g161 = (-(cse_p11_t12_g20));
      double cse_p11_t178_g186 = (-(((cos_theta) * (cos_theta))));
      double cse_p579_t5_sqrt402 = sqrt(cse_p11_t141_g149);
      double cse_p579_t6_invr403 = (1.0 / (cse_p11_t141_g149));
      double cse_p11_t193_g201 = ((cos_theta) * (cos_theta));
      if (fabs(cos_theta)>(1.0-VERYSMALL)) cos_theta = linear_angle_clamp(cos_theta);
      double cse_p1_t8_invsqrt8 = (1.0 / (cse_p579_t5_sqrt402));
      double sin_theta = cse_p579_t5_sqrt402;
      double theta = acos(cos_theta);
      double dtheta = (theta + (-(term.t0)));
      double cse_p62_t1_g223 = (2.0 * dtheta * term.kt);
      double energy = (dtheta * dtheta * term.kt);
      energy_added += energy;
      *energy_accumulate += energy;
      double cse_p64_t1_g248 = (dot * n2);
      double cse_p64_t1_g256 = (cse_p1_t2_invsqrt2 * cse_p1_t6_invr26 * cse_p64_t1_g248);
      double cse_p115_t1_g268 = (cse_p1_t8_invsqrt8 * cse_p62_t1_g223);
      double g_x1 = (-((cse_p115_t1_g268 * (cse_p11_t28_g36 + (-((cse_p64_t1_g256 * vx1)))))));
      KernelGradientAcc(term.i3x1, 0, g_x1);
      double g_y1 = (-((cse_p115_t1_g268 * (cse_p11_t30_g38 + (-((cse_p64_t1_g256 * vy1)))))));
      KernelGradientAcc(term.i3x1, 1, g_y1);
      double g_z1 = (-((cse_p115_t1_g268 * (cse_p11_t32_g40 + (-((cse_p64_t1_g256 * vz1)))))));
      KernelGradientAcc(term.i3x1, 2, g_z1);
      double cse_p117_t1_g287 = (cse_p1_t4_invsqrt4 * n1);
      double cse_p117_t1_g288 = (cse_p1_t6_invr26 * dot);
      double g_x2 = (-((cse_p115_t1_g268 * (cse_p11_t24_g32 + (-((cse_p117_t1_g288 * ((-((cse_p115_t1_g273 * cse_p1_t2_invsqrt2))) + (-((cse_p117_t1_g287 * vx2)))))))))));
      KernelGradientAcc(term.i3x2, 0, g_x2);
      double g_y2 = (-((cse_p115_t1_g268 * (cse_p11_t25_g33 + (-((cse_p117_t1_g288 * ((-((cse_p116_t1_g279 * cse_p1_t2_invsqrt2))) + (-((cse_p117_t1_g287 * vy2)))))))))));
      KernelGradientAcc(term.i3x2, 1, g_y2);
      double g_z2 = (-((cse_p115_t1_g268 * (cse_p11_t26_g34 + (-((cse_p117_t1_g288 * ((-((cse_p117_t1_g287 * vz2))) + (-((cse_p117_t1_g293 * cse_p1_t2_invsqrt2)))))))))));
      KernelGradientAcc(term.i3x2, 2, g_z2);
      double cse_p376_t1_g355 = (cse_p117_t1_g287 * cse_p117_t1_g288);
      double g_x3 = (-((cse_p115_t1_g268 * (cse_p11_t27_g35 + (-((cse_p376_t1_g355 * vx2)))))));
      KernelGradientAcc(term.i3x3, 0, g_x3);
      double g_y3 = (-((cse_p115_t1_g268 * (cse_p11_t29_g37 + (-((cse_p376_t1_g355 * vy2)))))));
      KernelGradientAcc(term.i3x3, 1, g_y3);
      double g_z3 = (-((cse_p115_t1_g268 * (cse_p11_t31_g39 + (-((cse_p376_t1_g355 * vz2)))))));
      KernelGradientAcc(term.i3x3, 2, g_z3);
      double cse_p11_t1_g209 = (cse_p116_t1_g283 * cse_p61_t8_invr2218);
      double cse_p63_t1_g239 = (cos_theta * cse_p579_t6_invr403);
      double cse_p63_t1_g240 = (cse_p61_t10_invsqrt220 * cse_p62_t1_g223);
      double cse_p64_t1_g254 = (cse_p11_t1_g209 * cse_p579_t2_invr399 * cse_p64_t1_g248);
      double cse_p580_t1_g404 = (cse_p64_t1_g254 * vx1);
      double cse_p116_t1_g275 = (cse_p579_t1_sqrt398 * cse_p64_t1_g248);
      double cse_p64_t1_g255 = (cse_p116_t1_g275 * cse_p116_t1_g283 * cse_p117_t1_g294);
      double cse_p115_t1_g265 = (cse_p579_t1_sqrt398 * dot);
      double cse_p116_t1_g285 = (cse_p115_t1_g265 * cse_p116_t1_g281 * cse_p322_t1_g335);
      double cse_p118_t1_g296 = (cse_p11_t133_g141 * cse_p63_t1_g239);
      double cse_p375_t1_g350 = (-2.0 * cse_p11_t1_g209);
      double h_x1_x1 = (cse_p11_t74_g82 + cse_p11_t74_g82 + (-((cse_p63_t1_g240 * (cse_p11_t158_g166 + (cse_p115_t1_g273 * ((cse_p116_t1_g285 * vx1) + (cse_p375_t1_g350 * vx2) + (cse_p64_t1_g255 * vx1))) + (cse_p118_t1_g296 * cse_p11_t133_g141) + (-((cse_p580_t1_g404 * vx1))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 0, h_x1_x1);
      double cse_p64_t1_g249 = (cse_p11_t1_g209 * n2);
      double cse_p117_t1_g290 = (cse_p64_t1_g249 * vx1);
      double cse_p117_t1_g291 = (cse_p64_t1_g249 * vy1);
      double cse_p376_t1_g354 = (cse_p116_t1_g275 * vx1);
      double h_x1_y1 = (cse_p11_t76_g84 + cse_p11_t76_g84 + (-((cse_p63_t1_g240 * ((cse_p118_t1_g296 * cse_p11_t135_g143) + (cse_p11_t127_g135 * cse_p376_t1_g354 * vy1) + (-((cse_p117_t1_g290 * vy2))) + (-((cse_p117_t1_g291 * vx2))) + (-((cse_p580_t1_g404 * vy1))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 1, h_x1_y1);
      double cse_p117_t1_g292 = (cse_p64_t1_g249 * vz1);
      double h_x1_z1 = (cse_p11_t78_g86 + cse_p11_t78_g86 + (-((cse_p63_t1_g240 * ((cse_p118_t1_g296 * cse_p11_t137_g145) + (cse_p11_t127_g135 * cse_p376_t1_g354 * vz1) + (-((cse_p117_t1_g290 * vz2))) + (-((cse_p117_t1_g292 * vx2))) + (-((cse_p580_t1_g404 * vz1))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x1, 2, h_x1_z1);
      double cse_p115_t1_g267 = (cse_p11_t129_g137 * cse_p63_t1_g239);
      double cse_p375_t1_g352 = (cse_p114_t1_g263 * cse_p61_t8_invr2218);
      double cse_p427_t1_g372 = (cse_p11_t1_g209 * dot);
      double h_x1_x2 = (cse_p11_t48_g56 + cse_p11_t48_g56 + (-((cse_p63_t1_g240 * (cse_p11_t179_g187 + (cse_p115_t1_g267 * cse_p11_t133_g141) + (cse_p11_t142_g150 * vx2) + (cse_p376_t1_g354 * (cse_p63_t1_g245 + (cse_p375_t1_g352 * vx1))) + (cse_p427_t1_g372 * (cse_p11_t105_g113 + cse_p11_t115_g123 + n2)) + (-((cse_p117_t1_g290 * cse_p11_t138_g146))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 0, h_x1_x2);
      double cse_p115_t1_g271 = (cse_p11_t130_g138 * cse_p63_t1_g239);
      double h_x1_y2 = (cse_p11_t56_g64 + cse_p11_t56_g64 + (-((cse_p63_t1_g240 * ((cse_p115_t1_g265 * cse_p63_t10_g233 * vx1) + (cse_p115_t1_g271 * cse_p11_t133_g141) + (cse_p11_t143_g151 * vx2) + (-((cse_p117_t1_g290 * cse_p11_t139_g147))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 1, h_x1_y2);
      double cse_p116_t1_g278 = (cse_p11_t131_g139 * cse_p63_t1_g239);
      double h_x1_z2 = (cse_p11_t63_g71 + cse_p11_t63_g71 + (-((cse_p63_t1_g240 * ((cse_p115_t1_g265 * cse_p63_t11_g234 * vx1) + (cse_p116_t1_g278 * cse_p11_t133_g141) + (cse_p11_t144_g152 * vx2) + (-((cse_p117_t1_g290 * cse_p11_t140_g148))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x2, 2, h_x1_z2);
      double cse_p11_t1_g210 = (cse_p579_t1_sqrt398 * cse_p63_t1_g244);
      double cse_p64_t1_g250 = (cse_p579_t4_invr401 * dot);
      double cse_p64_t1_g251 = (cse_p11_t1_g210 * cse_p63_t1_g241 * cse_p64_t1_g248 * n1);
      double cse_p64_t1_g253 = (cse_p11_t1_g209 * cse_p579_t3_sqrt400 * cse_p64_t1_g250);
      double cse_p115_t1_g264 = (cse_p377_t1_g364 * vx2);
      double cse_p117_t1_g286 = (cse_p11_t132_g140 * cse_p63_t1_g239);
      double h_x1_x3 = (cse_p11_t69_g77 + cse_p11_t69_g77 + (-((cse_p63_t1_g240 * (cse_p1_t5_invr5 + (cse_p117_t1_g286 * cse_p11_t133_g141) + (cse_p167_t1_g304 * cse_p64_t1_g251) + (-((cse_p115_t1_g264 * vx2))) + (-((cse_p117_t1_g290 * vx1))) + (-((cse_p167_t1_g304 * cse_p64_t1_g253))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 0, h_x1_x3);
      double cse_p375_t1_g344 = (vx1 * vy2);
      double h_x1_y3 = (cse_p11_t75_g83 + cse_p11_t75_g83 + (-((cse_p63_t1_g240 * (cse_p11_t163_g171 + cse_p11_t172_g180 + (cse_p118_t1_g296 * cse_p11_t134_g142) + (cse_p375_t1_g344 * cse_p64_t1_g251) + (-((cse_p375_t1_g344 * cse_p64_t1_g253))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 1, h_x1_y3);
      double cse_p375_t1_g345 = (vx1 * vz2);
      double h_x1_z3 = (cse_p11_t77_g85 + cse_p11_t77_g85 + (-((cse_p63_t1_g240 * (cse_p11_t165_g173 + cse_p11_t173_g181 + (cse_p118_t1_g296 * cse_p11_t136_g144) + (cse_p375_t1_g345 * cse_p64_t1_g251) + (-((cse_p375_t1_g345 * cse_p64_t1_g253))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 0, term.i3x3, 2, h_x1_z3);
      double cse_p376_t1_g363 = (cse_p11_t135_g143 * cse_p63_t1_g239);
      double h_y1_y1 = (cse_p11_t83_g91 + cse_p11_t83_g91 + (-((cse_p63_t1_g240 * (cse_p11_t158_g166 + (cse_p116_t1_g279 * ((cse_p116_t1_g285 * vy1) + (cse_p375_t1_g350 * vy2) + (cse_p64_t1_g255 * vy1))) + (cse_p11_t135_g143 * cse_p376_t1_g363) + (-((cse_p64_t1_g254 * vy1 * vy1))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 1, h_y1_y1);
      double h_y1_z1 = (cse_p11_t85_g93 + cse_p11_t85_g93 + (-((cse_p63_t1_g240 * ((cse_p116_t1_g275 * cse_p11_t127_g135 * vy1 * vz1) + (cse_p11_t137_g145 * cse_p376_t1_g363) + (-((cse_p117_t1_g291 * vz2))) + (-((cse_p117_t1_g292 * vy2))) + (-((cse_p64_t1_g254 * vy1 * vz1))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x1, 2, h_y1_z1);
      double h_y1_x2 = (cse_p11_t50_g58 + cse_p11_t50_g58 + (-((cse_p63_t1_g240 * ((cse_p115_t1_g265 * cse_p63_t9_g232 * vy1) + (cse_p115_t1_g267 * cse_p11_t135_g143) + (cse_p11_t142_g150 * vy2) + (-((cse_p117_t1_g291 * cse_p11_t138_g146))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 0, h_y1_x2);
      double h_y1_y2 = (cse_p11_t58_g66 + cse_p11_t58_g66 + (-((cse_p63_t1_g240 * (cse_p11_t179_g187 + (cse_p115_t1_g271 * cse_p11_t135_g143) + (cse_p116_t1_g275 * vy1 * (cse_p63_t1_g246 + (cse_p375_t1_g352 * vy1))) + (cse_p11_t143_g151 * vy2) + (cse_p427_t1_g372 * (cse_p11_t110_g118 + cse_p11_t117_g125 + n2)) + (-((cse_p117_t1_g291 * cse_p11_t139_g147))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 1, h_y1_y2);
      double h_y1_z2 = (cse_p11_t65_g73 + cse_p11_t65_g73 + (-((cse_p63_t1_g240 * ((cse_p115_t1_g265 * cse_p63_t11_g234 * vy1) + (cse_p116_t1_g278 * cse_p11_t135_g143) + (cse_p11_t144_g152 * vy2) + (-((cse_p117_t1_g291 * cse_p11_t140_g148))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x2, 2, h_y1_z2);
      double cse_p375_t1_g346 = (vx2 * vy1);
      double h_y1_x3 = (cse_p11_t71_g79 + cse_p11_t71_g79 + (-((cse_p63_t1_g240 * (cse_p11_t163_g171 + cse_p11_t172_g180 + (cse_p117_t1_g286 * cse_p11_t135_g143) + (cse_p375_t1_g346 * cse_p64_t1_g251) + (-((cse_p375_t1_g346 * cse_p64_t1_g253))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 0, h_y1_x3);
      double cse_p115_t1_g269 = (cse_p377_t1_g364 * vy2);
      double cse_p168_t1_g314 = (cse_p11_t134_g142 * cse_p63_t1_g239);
      double h_y1_y3 = (cse_p11_t80_g88 + cse_p11_t80_g88 + (-((cse_p63_t1_g240 * (cse_p1_t5_invr5 + (cse_p11_t135_g143 * cse_p168_t1_g314) + (cse_p167_t1_g308 * cse_p64_t1_g251) + (-((cse_p115_t1_g269 * vy2))) + (-((cse_p117_t1_g291 * vy1))) + (-((cse_p167_t1_g308 * cse_p64_t1_g253))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 1, h_y1_y3);
      double cse_p375_t1_g348 = (vy1 * vz2);
      double h_y1_z3 = (cse_p11_t84_g92 + cse_p11_t84_g92 + (-((cse_p63_t1_g240 * (cse_p11_t168_g176 + cse_p11_t175_g183 + (cse_p11_t136_g144 * cse_p376_t1_g363) + (cse_p375_t1_g348 * cse_p64_t1_g251) + (-((cse_p375_t1_g348 * cse_p64_t1_g253))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 1, term.i3x3, 2, h_y1_z3);
      double h_z1_z1 = (cse_p11_t88_g96 + cse_p11_t88_g96 + (-((cse_p63_t1_g240 * (cse_p11_t158_g166 + (cse_p117_t1_g293 * ((cse_p116_t1_g285 * vz1) + (cse_p375_t1_g350 * vz2) + (cse_p64_t1_g255 * vz1))) + (cse_p11_t137_g145 * cse_p11_t137_g145 * cse_p63_t1_g239) + (-((cse_p64_t1_g254 * vz1 * vz1))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x1, 2, h_z1_z1);
      double h_z1_x2 = (cse_p11_t52_g60 + cse_p11_t52_g60 + (-((cse_p63_t1_g240 * ((cse_p115_t1_g265 * cse_p63_t9_g232 * vz1) + (cse_p115_t1_g267 * cse_p11_t137_g145) + (cse_p11_t142_g150 * vz2) + (-((cse_p117_t1_g292 * cse_p11_t138_g146))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 0, h_z1_x2);
      double h_z1_y2 = (cse_p11_t60_g68 + cse_p11_t60_g68 + (-((cse_p63_t1_g240 * ((cse_p115_t1_g265 * cse_p63_t10_g233 * vz1) + (cse_p115_t1_g271 * cse_p11_t137_g145) + (cse_p11_t143_g151 * vz2) + (-((cse_p117_t1_g292 * cse_p11_t139_g147))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 1, h_z1_y2);
      double h_z1_z2 = (cse_p11_t67_g75 + cse_p11_t67_g75 + (-((cse_p63_t1_g240 * (cse_p11_t179_g187 + (cse_p116_t1_g275 * vz1 * (cse_p63_t1_g247 + (cse_p375_t1_g352 * vz1))) + (cse_p116_t1_g278 * cse_p11_t137_g145) + (cse_p11_t144_g152 * vz2) + (cse_p427_t1_g372 * (cse_p11_t114_g122 + cse_p11_t119_g127 + n2)) + (-((cse_p117_t1_g292 * cse_p11_t140_g148))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x2, 2, h_z1_z2);
      double cse_p375_t1_g347 = (vx2 * vz1);
      double h_z1_x3 = (cse_p11_t73_g81 + cse_p11_t73_g81 + (-((cse_p63_t1_g240 * (cse_p11_t165_g173 + cse_p11_t173_g181 + (cse_p117_t1_g286 * cse_p11_t137_g145) + (cse_p375_t1_g347 * cse_p64_t1_g251) + (-((cse_p375_t1_g347 * cse_p64_t1_g253))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 0, h_z1_x3);
      double cse_p375_t1_g349 = (vy2 * vz1);
      double h_z1_y3 = (cse_p11_t82_g90 + cse_p11_t82_g90 + (-((cse_p63_t1_g240 * (cse_p11_t168_g176 + cse_p11_t175_g183 + (cse_p11_t137_g145 * cse_p168_t1_g314) + (cse_p375_t1_g349 * cse_p64_t1_g251) + (-((cse_p375_t1_g349 * cse_p64_t1_g253))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 1, h_z1_y3);
      double cse_p115_t1_g272 = (cse_p377_t1_g364 * vz2);
      double h_z1_z3 = (cse_p11_t87_g95 + cse_p11_t87_g95 + (-((cse_p63_t1_g240 * (cse_p1_t5_invr5 + (cse_p11_t136_g144 * cse_p11_t137_g145 * cse_p63_t1_g239) + (cse_p167_t1_g311 * cse_p64_t1_g251) + (-((cse_p115_t1_g272 * vz2))) + (-((cse_p117_t1_g292 * vz1))) + (-((cse_p167_t1_g311 * cse_p64_t1_g253))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x1, 2, term.i3x3, 2, h_z1_z3);
      double cse_p116_t1_g280 = (-2.0 * cse_p64_t1_g252);
      double cse_p375_t1_g351 = (cse_p114_t1_g263 * cse_p579_t1_sqrt398);
      double cse_p376_t1_g358 = (cse_p11_t121_g129 * cse_p63_t1_g241);
      double h_x2_x2 = (cse_p11_t44_g52 + cse_p11_t44_g52 + (-((cse_p63_t1_g240 * (cse_p11_t18_g26 + (2.0 * cse_p11_t138_g146 * cse_p11_t142_g150) + (cse_p115_t1_g267 * cse_p11_t129_g137) + (dot * ((cse_p11_t121_g129 * cse_p376_t1_g358) + (-((cse_p61_t8_invr2218 * (cse_p11_t92_g100 + cse_p478_t1_g379 + (cse_p115_t1_g273 * cse_p375_t1_g351 * vx1) + (cse_p116_t1_g280 * cse_p167_t1_g305) + (cse_p116_t1_g283 * (cse_p11_t105_g113 + cse_p11_t115_g123 + cse_p11_t115_g123 + n2)))))))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 0, h_x2_x2);
      double cse_p115_t1_g274 = (cse_p115_t1_g273 * cse_p579_t1_sqrt398 * cse_p63_t13_g236);
      double h_x2_y2 = (cse_p11_t45_g53 + cse_p11_t45_g53 + (-((cse_p63_t1_g240 * ((cse_p115_t1_g267 * cse_p11_t130_g138) + (cse_p11_t138_g146 * cse_p11_t143_g151) + (cse_p11_t139_g147 * cse_p11_t142_g150) + (dot * ((cse_p11_t122_g130 * cse_p376_t1_g358) + (-((cse_p61_t8_invr2218 * (cse_p11_t94_g102 + (cse_p115_t1_g274 * vy1) + (cse_p116_t1_g280 * cse_p167_t1_g306) + (cse_p11_t1_g210 * (cse_p375_t1_g344 + cse_p375_t1_g346)))))))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 1, h_x2_y2);
      double h_x2_z2 = (cse_p11_t46_g54 + cse_p11_t46_g54 + (-((cse_p63_t1_g240 * ((cse_p115_t1_g267 * cse_p11_t131_g139) + (cse_p11_t138_g146 * cse_p11_t144_g152) + (cse_p11_t140_g148 * cse_p11_t142_g150) + (dot * ((cse_p11_t123_g131 * cse_p376_t1_g358) + (-((cse_p61_t8_invr2218 * (cse_p11_t96_g104 + (cse_p115_t1_g274 * vz1) + (cse_p116_t1_g280 * cse_p167_t1_g307) + (cse_p11_t1_g210 * (cse_p375_t1_g345 + cse_p375_t1_g347)))))))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x2, 2, h_x2_z2);
      double cse_p376_t1_g359 = (cse_p11_t121_g129 * cse_p63_t1_g243);
      double h_x2_x3 = (cse_p11_t47_g55 + cse_p11_t47_g55 + (-((cse_p63_t1_g240 * (cse_p11_t179_g187 + (cse_p115_t1_g267 * cse_p11_t132_g140) + (cse_p11_t142_g150 * vx1) + (dot * ((cse_p376_t1_g359 * vx2) + (-((cse_p61_t8_invr2218 * (cse_p11_t177_g185 + (cse_p116_t1_g282 * cse_p167_t1_g305) + (-((cse_p11_t1_g210 * cse_p167_t1_g304))) + (-(cse_p11_t92_g100)))))))) + (-((cse_p115_t1_g264 * cse_p11_t138_g146))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 0, h_x2_x3);
      double h_x2_y3 = (cse_p11_t49_g57 + cse_p11_t49_g57 + (-((cse_p63_t1_g240 * ((cse_p115_t1_g267 * cse_p11_t134_g142) + (cse_p11_t142_g150 * vy1) + (dot * ((cse_p376_t1_g359 * vy2) + (-((cse_p61_t8_invr2218 * (cse_p11_t164_g172 + cse_p11_t20_g28 + (-((cse_p11_t1_g210 * cse_p375_t1_g344))))))))) + (-((cse_p115_t1_g269 * cse_p11_t138_g146))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 1, h_x2_y3);
      double h_x2_z3 = (cse_p11_t51_g59 + cse_p11_t51_g59 + (-((cse_p63_t1_g240 * ((cse_p115_t1_g267 * cse_p11_t136_g144) + (cse_p11_t142_g150 * vz1) + (dot * ((cse_p376_t1_g359 * vz2) + (-((cse_p61_t8_invr2218 * (cse_p11_t166_g174 + cse_p11_t21_g29 + (-((cse_p11_t1_g210 * cse_p375_t1_g345))))))))) + (-((cse_p115_t1_g272 * cse_p11_t138_g146))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 0, term.i3x3, 2, h_x2_z3);
      double h_y2_y2 = (cse_p11_t53_g61 + cse_p11_t53_g61 + (-((cse_p63_t1_g240 * (cse_p11_t18_g26 + (2.0 * cse_p11_t139_g147 * cse_p11_t143_g151) + (cse_p115_t1_g271 * cse_p11_t130_g138) + (dot * ((cse_p11_t122_g130 * cse_p11_t122_g130 * cse_p63_t1_g241) + (-((cse_p61_t8_invr2218 * (cse_p11_t98_g106 + cse_p478_t1_g379 + (cse_p116_t1_g279 * cse_p375_t1_g351 * vy1) + (cse_p116_t1_g280 * cse_p167_t1_g309) + (cse_p116_t1_g283 * (cse_p11_t110_g118 + cse_p11_t117_g125 + cse_p11_t117_g125 + n2)))))))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 1, h_y2_y2);
      double h_y2_z2 = (cse_p11_t54_g62 + cse_p11_t54_g62 + (-((cse_p63_t1_g240 * ((cse_p115_t1_g271 * cse_p11_t131_g139) + (cse_p11_t139_g147 * cse_p11_t144_g152) + (cse_p11_t140_g148 * cse_p11_t143_g151) + (dot * ((cse_p11_t122_g130 * cse_p11_t123_g131 * cse_p63_t1_g241) + (-((cse_p61_t8_invr2218 * (cse_p11_t100_g108 + (cse_p116_t1_g279 * cse_p579_t1_sqrt398 * cse_p63_t13_g236 * vz1) + (cse_p116_t1_g280 * cse_p167_t1_g310) + (cse_p11_t1_g210 * (cse_p375_t1_g348 + cse_p375_t1_g349)))))))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x2, 2, h_y2_z2);
      double cse_p376_t1_g360 = (cse_p11_t122_g130 * cse_p63_t1_g243);
      double h_y2_x3 = (cse_p11_t55_g63 + cse_p11_t55_g63 + (-((cse_p63_t1_g240 * ((cse_p115_t1_g271 * cse_p11_t132_g140) + (cse_p11_t143_g151 * vx1) + (dot * ((cse_p376_t1_g360 * vx2) + (-((cse_p61_t8_invr2218 * (cse_p11_t164_g172 + cse_p11_t20_g28 + (-((cse_p11_t1_g210 * cse_p375_t1_g346))))))))) + (-((cse_p115_t1_g264 * cse_p11_t139_g147))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 0, h_y2_x3);
      double h_y2_y3 = (cse_p11_t57_g65 + cse_p11_t57_g65 + (-((cse_p63_t1_g240 * (cse_p11_t179_g187 + (cse_p115_t1_g271 * cse_p11_t134_g142) + (cse_p11_t143_g151 * vy1) + (dot * ((cse_p376_t1_g360 * vy2) + (-((cse_p61_t8_invr2218 * (cse_p11_t177_g185 + (cse_p116_t1_g282 * cse_p167_t1_g309) + (-((cse_p11_t1_g210 * cse_p167_t1_g308))) + (-(cse_p11_t98_g106)))))))) + (-((cse_p115_t1_g269 * cse_p11_t139_g147))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 1, h_y2_y3);
      double h_y2_z3 = (cse_p11_t59_g67 + cse_p11_t59_g67 + (-((cse_p63_t1_g240 * ((cse_p115_t1_g271 * cse_p11_t136_g144) + (cse_p11_t143_g151 * vz1) + (dot * ((cse_p376_t1_g360 * vz2) + (-((cse_p61_t8_invr2218 * (cse_p11_t169_g177 + cse_p11_t22_g30 + (-((cse_p11_t1_g210 * cse_p375_t1_g348))))))))) + (-((cse_p115_t1_g272 * cse_p11_t139_g147))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 1, term.i3x3, 2, h_y2_z3);
      double h_z2_z2 = (cse_p11_t61_g69 + cse_p11_t61_g69 + (-((cse_p63_t1_g240 * (cse_p11_t18_g26 + (2.0 * cse_p11_t140_g148 * cse_p11_t144_g152) + (cse_p116_t1_g278 * cse_p11_t131_g139) + (dot * ((cse_p11_t123_g131 * cse_p11_t123_g131 * cse_p63_t1_g241) + (-((cse_p61_t8_invr2218 * (cse_p11_t102_g110 + cse_p478_t1_g379 + (cse_p116_t1_g280 * cse_p167_t1_g312) + (cse_p116_t1_g283 * (cse_p11_t114_g122 + cse_p11_t119_g127 + cse_p11_t119_g127 + n2)) + (cse_p117_t1_g293 * cse_p375_t1_g351 * vz1))))))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x2, 2, h_z2_z2);
      double cse_p376_t1_g361 = (cse_p11_t123_g131 * cse_p63_t1_g243);
      double h_z2_x3 = (cse_p11_t62_g70 + cse_p11_t62_g70 + (-((cse_p63_t1_g240 * ((cse_p116_t1_g278 * cse_p11_t132_g140) + (cse_p11_t144_g152 * vx1) + (dot * ((cse_p376_t1_g361 * vx2) + (-((cse_p61_t8_invr2218 * (cse_p11_t166_g174 + cse_p11_t21_g29 + (-((cse_p11_t1_g210 * cse_p375_t1_g347))))))))) + (-((cse_p115_t1_g264 * cse_p11_t140_g148))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 0, h_z2_x3);
      double h_z2_y3 = (cse_p11_t64_g72 + cse_p11_t64_g72 + (-((cse_p63_t1_g240 * ((cse_p116_t1_g278 * cse_p11_t134_g142) + (cse_p11_t144_g152 * vy1) + (dot * ((cse_p376_t1_g361 * vy2) + (-((cse_p61_t8_invr2218 * (cse_p11_t169_g177 + cse_p11_t22_g30 + (-((cse_p11_t1_g210 * cse_p375_t1_g349))))))))) + (-((cse_p115_t1_g269 * cse_p11_t140_g148))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 1, h_z2_y3);
      double h_z2_z3 = (cse_p11_t66_g74 + cse_p11_t66_g74 + (-((cse_p63_t1_g240 * (cse_p11_t179_g187 + (cse_p116_t1_g278 * cse_p11_t136_g144) + (cse_p11_t144_g152 * vz1) + (dot * ((cse_p376_t1_g361 * vz2) + (-((cse_p61_t8_invr2218 * (cse_p11_t177_g185 + (cse_p116_t1_g282 * cse_p167_t1_g312) + (-((cse_p11_t1_g210 * cse_p167_t1_g311))) + (-(cse_p11_t102_g110)))))))) + (-((cse_p115_t1_g272 * cse_p11_t140_g148))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x2, 2, term.i3x3, 2, h_z2_z3);
      double cse_p117_t1_g289 = (cse_p63_t1_g237 * dot);
      double cse_p375_t1_g353 = (cse_p115_t1_g264 * cse_p64_t1_g250);
      double cse_p376_t1_g357 = (cse_p117_t1_g289 * cse_p63_t12_g235);
      double h_x3_x3 = (cse_p11_t68_g76 + cse_p11_t68_g76 + (-((cse_p63_t1_g240 * (cse_p11_t154_g162 + (-2.0 * cse_p115_t1_g264 * vx1) + (cse_p117_t1_g286 * cse_p11_t132_g140) + (cse_p167_t1_g305 * cse_p376_t1_g357) + (-((cse_p375_t1_g353 * vx2))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 0, h_x3_x3);
      double cse_p376_t1_g356 = (cse_p117_t1_g289 * cse_p11_t128_g136);
      double h_x3_y3 = (cse_p11_t70_g78 + cse_p11_t70_g78 + (-((cse_p63_t1_g240 * ((cse_p117_t1_g286 * cse_p11_t134_g142) + (cse_p167_t1_g306 * cse_p376_t1_g356) + (-((cse_p115_t1_g264 * vy1))) + (-((cse_p115_t1_g269 * vx1))) + (-((cse_p375_t1_g353 * vy2))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 1, h_x3_y3);
      double h_x3_z3 = (cse_p11_t72_g80 + cse_p11_t72_g80 + (-((cse_p63_t1_g240 * ((cse_p117_t1_g286 * cse_p11_t136_g144) + (cse_p167_t1_g307 * cse_p376_t1_g356) + (-((cse_p115_t1_g264 * vz1))) + (-((cse_p115_t1_g272 * vx1))) + (-((cse_p375_t1_g353 * vz2))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 0, term.i3x3, 2, h_x3_z3);
      double h_y3_y3 = (cse_p11_t79_g87 + cse_p11_t79_g87 + (-((cse_p63_t1_g240 * (cse_p11_t154_g162 + (-2.0 * cse_p115_t1_g269 * vy1) + (cse_p11_t134_g142 * cse_p168_t1_g314) + (cse_p167_t1_g309 * cse_p376_t1_g357) + (-((cse_p115_t1_g269 * cse_p64_t1_g250 * vy2))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 1, h_y3_y3);
      double h_y3_z3 = (cse_p11_t81_g89 + cse_p11_t81_g89 + (-((cse_p63_t1_g240 * ((cse_p11_t136_g144 * cse_p168_t1_g314) + (cse_p167_t1_g310 * cse_p376_t1_g356) + (-((cse_p115_t1_g269 * cse_p64_t1_g250 * vz2))) + (-((cse_p115_t1_g269 * vz1))) + (-((cse_p115_t1_g272 * vy1))))))));
      KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 1, term.i3x3, 2, h_y3_z3);
      double h_z3_z3 = (cse_p11_t86_g94 + cse_p11_t86_g94 + (-((cse_p63_t1_g240 * (cse_p11_t154_g162 + (-2.0 * cse_p115_t1_g272 * vz1) + (cse_p11_t136_g144 * cse_p11_t136_g144 * cse_p63_t1_g239) + (cse_p167_t1_g312 * cse_p376_t1_g357) + (-((cse_p115_t1_g272 * cse_p64_t1_g250 * vz2))))))));
      KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec,  term.i3x3, 2, term.i3x3, 2, h_z3_z3);
    }
  }
  return energy_added;
}
void hessian_fd(const angle_term& term, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
{
  constexpr size_t PositionSize = 9;
  const double h = 1.0e-5;
  const double inv2h = 1.0/(2.0*h);
  const double invh2 = 1.0/((h*h));
  double e0[EnergyAccumulatorSize] = {0.0};
  energy(term, position, e0);
  if (energy_accumulate) { *energy_accumulate += e0[0]; }
  {
    double saved = position[term.i3x1 + 0];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x1 + 0] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x1 + 0] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x1 + 0] = saved;
    double d = (e_plus[0] - e_minus[0]) * inv2h;
    KernelGradientAcc(term.i3x1, 0, d);
  }
  {
    double saved = position[term.i3x1 + 1];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x1 + 1] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x1 + 1] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x1 + 1] = saved;
    double d = (e_plus[0] - e_minus[0]) * inv2h;
    KernelGradientAcc(term.i3x1, 1, d);
  }
  {
    double saved = position[term.i3x1 + 2];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x1 + 2] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x1 + 2] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x1 + 2] = saved;
    double d = (e_plus[0] - e_minus[0]) * inv2h;
    KernelGradientAcc(term.i3x1, 2, d);
  }
  {
    double saved = position[term.i3x2 + 0];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 0] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x2 + 0] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x2 + 0] = saved;
    double d = (e_plus[0] - e_minus[0]) * inv2h;
    KernelGradientAcc(term.i3x2, 0, d);
  }
  {
    double saved = position[term.i3x2 + 1];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 1] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x2 + 1] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x2 + 1] = saved;
    double d = (e_plus[0] - e_minus[0]) * inv2h;
    KernelGradientAcc(term.i3x2, 1, d);
  }
  {
    double saved = position[term.i3x2 + 2];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 2] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x2 + 2] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x2 + 2] = saved;
    double d = (e_plus[0] - e_minus[0]) * inv2h;
    KernelGradientAcc(term.i3x2, 2, d);
  }
  {
    double saved = position[term.i3x3 + 0];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 0] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x3 + 0] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x3 + 0] = saved;
    double d = (e_plus[0] - e_minus[0]) * inv2h;
    KernelGradientAcc(term.i3x3, 0, d);
  }
  {
    double saved = position[term.i3x3 + 1];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 1] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x3 + 1] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x3 + 1] = saved;
    double d = (e_plus[0] - e_minus[0]) * inv2h;
    KernelGradientAcc(term.i3x3, 1, d);
  }
  {
    double saved = position[term.i3x3 + 2];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 2] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x3 + 2] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x3 + 2] = saved;
    double d = (e_plus[0] - e_minus[0]) * inv2h;
    KernelGradientAcc(term.i3x3, 2, d);
  }
  {
    double saved = position[term.i3x1 + 0];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x1 + 0] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x1 + 0] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x1 + 0] = saved;
    double hval = (e_plus[0] + e_minus[0] - (2.0*e0[0])) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 0, term.i3x1, 0, hval);
  }
  {
    double saved = position[term.i3x1 + 1];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x1 + 1] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x1 + 1] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x1 + 1] = saved;
    double hval = (e_plus[0] + e_minus[0] - (2.0*e0[0])) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 1, term.i3x1, 1, hval);
  }
  {
    double saved = position[term.i3x1 + 2];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x1 + 2] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x1 + 2] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x1 + 2] = saved;
    double hval = (e_plus[0] + e_minus[0] - (2.0*e0[0])) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 2, term.i3x1, 2, hval);
  }
  {
    double saved = position[term.i3x2 + 0];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 0] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x2 + 0] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x2 + 0] = saved;
    double hval = (e_plus[0] + e_minus[0] - (2.0*e0[0])) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 0, term.i3x2, 0, hval);
  }
  {
    double saved = position[term.i3x2 + 1];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 1] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x2 + 1] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x2 + 1] = saved;
    double hval = (e_plus[0] + e_minus[0] - (2.0*e0[0])) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 1, term.i3x2, 1, hval);
  }
  {
    double saved = position[term.i3x2 + 2];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 2] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x2 + 2] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x2 + 2] = saved;
    double hval = (e_plus[0] + e_minus[0] - (2.0*e0[0])) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 2, term.i3x2, 2, hval);
  }
  {
    double saved = position[term.i3x3 + 0];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 0] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x3 + 0] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x3 + 0] = saved;
    double hval = (e_plus[0] + e_minus[0] - (2.0*e0[0])) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x3, 0, hval);
  }
  {
    double saved = position[term.i3x3 + 1];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 1] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x3 + 1] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x3 + 1] = saved;
    double hval = (e_plus[0] + e_minus[0] - (2.0*e0[0])) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x3, 1, hval);
  }
  {
    double saved = position[term.i3x3 + 2];
    double e_plus[EnergyAccumulatorSize] = {0.0};
    double e_minus[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 2] = saved + h;
    energy(term, position, e_plus);
    position[term.i3x3 + 2] = saved - h;
    energy(term, position, e_minus);
    position[term.i3x3 + 2] = saved;
    double hval = (e_plus[0] + e_minus[0] - (2.0*e0[0])) * invh2;
    KernelHessDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x3, 2, hval);
  }
  {
    double saved_i = position[term.i3x1 + 1];
    double saved_j = position[term.i3x1 + 0];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x1 + 1] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x1 + 1] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x1 + 1] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 1, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x1 + 2];
    double saved_j = position[term.i3x1 + 0];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x1 + 2] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x1 + 2] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x1 + 2] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 2, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x1 + 2];
    double saved_j = position[term.i3x1 + 1];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x1 + 2] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x1 + 2] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x1 + 2] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x1, 2, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x2 + 0];
    double saved_j = position[term.i3x1 + 0];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 0] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x2 + 0] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x2 + 0] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 0, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x2 + 0];
    double saved_j = position[term.i3x1 + 1];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 0] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x2 + 0] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x2 + 0] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 0, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x2 + 0];
    double saved_j = position[term.i3x1 + 2];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 0] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x2 + 0] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x2 + 0] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 0, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x2 + 1];
    double saved_j = position[term.i3x1 + 0];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 1] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x2 + 1] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x2 + 1] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 1, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x2 + 1];
    double saved_j = position[term.i3x1 + 1];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 1] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x2 + 1] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x2 + 1] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 1, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x2 + 1];
    double saved_j = position[term.i3x1 + 2];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 1] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x2 + 1] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x2 + 1] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 1, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x2 + 1];
    double saved_j = position[term.i3x2 + 0];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 1] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x2 + 1] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x2 + 1] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 1, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x2 + 2];
    double saved_j = position[term.i3x1 + 0];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 2] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x2 + 2] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 2, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x2 + 2];
    double saved_j = position[term.i3x1 + 1];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 2] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x2 + 2] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 2, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x2 + 2];
    double saved_j = position[term.i3x1 + 2];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 2] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x2 + 2] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 2, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x2 + 2];
    double saved_j = position[term.i3x2 + 0];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 2] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x2 + 2] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 2, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x2 + 2];
    double saved_j = position[term.i3x2 + 1];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x2 + 2] = saved_i + h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x2 + 2] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x2 + 2] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x2, 2, term.i3x2, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 0];
    double saved_j = position[term.i3x1 + 0];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 0] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 0] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 0] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 0];
    double saved_j = position[term.i3x1 + 1];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 0] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 0] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 0] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 0];
    double saved_j = position[term.i3x1 + 2];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 0] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 0] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 0] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x3 + 0];
    double saved_j = position[term.i3x2 + 0];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 0] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 0] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 0] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 0];
    double saved_j = position[term.i3x2 + 1];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 0] = saved_i + h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 0] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 0] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x2, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 0];
    double saved_j = position[term.i3x2 + 2];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 0] = saved_i + h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 0] = saved_i - h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 0] = saved_i; position[term.i3x2 + 2] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 0, term.i3x2, 2, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x1 + 0];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x1 + 1];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x1 + 2];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x2 + 0];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x2 + 1];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x2, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x2 + 2];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x2 + 2] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x2, 2, hval);
  }
  {
    double saved_i = position[term.i3x3 + 1];
    double saved_j = position[term.i3x3 + 0];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 1] = saved_i + h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 1] = saved_i - h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 1] = saved_i; position[term.i3x3 + 0] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 1, term.i3x3, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x1 + 0];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x1 + 0] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 0] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x1 + 0] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x1, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x1 + 1];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x1 + 1] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 1] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x1 + 1] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x1, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x1 + 2];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x1 + 2] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x1 + 2] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x1 + 2] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x1, 2, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x2 + 0];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x2 + 0] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x2 + 0] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x2 + 0] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x2, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x2 + 1];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x2 + 1] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x2 + 1] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x2 + 1] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x2, 1, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x2 + 2];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x2 + 2] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x2 + 2] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x2 + 2] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x2, 2, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x3 + 0];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x3 + 0] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x3 + 0] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x3 + 0] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x3, 0, hval);
  }
  {
    double saved_i = position[term.i3x3 + 2];
    double saved_j = position[term.i3x3 + 1];
    double e_pp[EnergyAccumulatorSize] = {0.0};
    double e_pm[EnergyAccumulatorSize] = {0.0};
    double e_mp[EnergyAccumulatorSize] = {0.0};
    double e_mm[EnergyAccumulatorSize] = {0.0};
    position[term.i3x3 + 2] = saved_i + h; position[term.i3x3 + 1] = saved_j + h;
    energy(term, position, e_pp);
    position[term.i3x3 + 1] = saved_j - h;
    energy(term, position, e_pm);
    position[term.i3x3 + 2] = saved_i - h; position[term.i3x3 + 1] = saved_j + h;
    energy(term, position, e_mp);
    position[term.i3x3 + 1] = saved_j - h;
    energy(term, position, e_mm);
    position[term.i3x3 + 2] = saved_i; position[term.i3x3 + 1] = saved_j;
    double hval = (e_pp[0] - e_pm[0] - e_mp[0] + e_mm[0]) * (0.25*invh2);
    KernelHessOffDiagAcc( PositionSize, hessian, dvec, hdvec, term.i3x3, 2, term.i3x3, 1, hval);
  }
}

};
