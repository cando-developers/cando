#define NDIM 9

typedef struct {
  double val;
  double g[NDIM];
  double h[NDIM][NDIM];
} Dual2;

/* Small helpers for second-order dual numbers */
static inline Dual2 dual_const(double v) {
  Dual2 d;
  d.val = v;
  for (int i = 0; i < NDIM; ++i) {
    d.g[i] = 0.0;
    for (int j = 0; j < NDIM; ++j) {
      d.h[i][j] = 0.0;
    }
  }
  return d;
}

static inline Dual2 dual_var(double v, int idx) {
  Dual2 d = dual_const(v);
  d.g[idx] = 1.0;
  return d;
}

static inline Dual2 dual_add(Dual2 a, Dual2 b) {
  Dual2 r;
  r.val = a.val + b.val;
  for (int i = 0; i < NDIM; ++i) {
    r.g[i] = a.g[i] + b.g[i];
    for (int j = 0; j < NDIM; ++j) {
      r.h[i][j] = a.h[i][j] + b.h[i][j];
    }
  }
  return r;
}

static inline Dual2 dual_sub(Dual2 a, Dual2 b) {
  Dual2 r;
  r.val = a.val - b.val;
  for (int i = 0; i < NDIM; ++i) {
    r.g[i] = a.g[i] - b.g[i];
    for (int j = 0; j < NDIM; ++j) {
      r.h[i][j] = a.h[i][j] - b.h[i][j];
    }
  }
  return r;
}

static inline Dual2 dual_mul(Dual2 a, Dual2 b) {
  Dual2 r;
  r.val = a.val * b.val;
  for (int i = 0; i < NDIM; ++i) {
    r.g[i] = a.g[i] * b.val + b.g[i] * a.val;
  }
  for (int i = 0; i < NDIM; ++i) {
    for (int j = 0; j < NDIM; ++j) {
      r.h[i][j] = a.h[i][j] * b.val + b.h[i][j] * a.val + a.g[i] * b.g[j] + b.g[i] * a.g[j];
    }
  }
  return r;
}

static inline Dual2 dual_mul_const(Dual2 a, double c) {
  Dual2 r;
  r.val = a.val * c;
  for (int i = 0; i < NDIM; ++i) {
    r.g[i] = a.g[i] * c;
    for (int j = 0; j < NDIM; ++j) {
      r.h[i][j] = a.h[i][j] * c;
    }
  }
  return r;
}

static inline Dual2 dual_neg(Dual2 a) {
  return dual_mul_const(a, -1.0);
}

static inline Dual2 dual_dot(const Dual2 a[3], const Dual2 b[3]) {
  Dual2 r = dual_const(0.0);
  for (int k = 0; k < 3; ++k) {
    r = dual_add(r, dual_mul(a[k], b[k]));
  }
  return r;
}

static inline void dual_cross(const Dual2 a[3], const Dual2 b[3], Dual2 out[3]) {
  out[0] = dual_sub(dual_mul(a[1], b[2]), dual_mul(a[2], b[1]));
  out[1] = dual_sub(dual_mul(a[2], b[0]), dual_mul(a[0], b[2]));
  out[2] = dual_sub(dual_mul(a[0], b[1]), dual_mul(a[1], b[0]));
}

static inline Dual2 dual_sqrt(Dual2 a) {
  double root = sqrt(a.val);
  double inv_root = 1.0 / root;
  double g_scale = 0.5 * inv_root;
  double h_scale = -0.25 * inv_root * inv_root * inv_root;

  Dual2 r = dual_const(0.0);
  r.val = root;
  for (int i = 0; i < NDIM; ++i) {
    r.g[i] = a.g[i] * g_scale;
  }
  for (int i = 0; i < NDIM; ++i) {
    for (int j = 0; j < NDIM; ++j) {
      r.h[i][j] = a.h[i][j] * g_scale + h_scale * a.g[i] * a.g[j];
    }
  }
  return r;
}

static inline Dual2 dual_atan2(Dual2 y, Dual2 x) {
  Dual2 r = dual_const(0.0);
  r.val = atan2(y.val, x.val);
  double den = x.val * x.val + y.val * y.val;

  double N[NDIM];
  for (int i = 0; i < NDIM; ++i) {
    N[i] = x.val * y.g[i] - y.val * x.g[i];
    r.g[i] = N[i] / den;
  }

  for (int i = 0; i < NDIM; ++i) {
    for (int j = 0; j < NDIM; ++j) {
      double Nij = x.g[j] * y.g[i] + x.val * y.h[i][j] - y.g[j] * x.g[i] - y.val * x.h[i][j];
      double den_j = 2.0 * x.val * x.g[j] + 2.0 * y.val * y.g[j];
      r.h[i][j] = (Nij * den - N[i] * den_j) / (den * den);
    }
  }
  return r;
}

void dihedral_hessian_ai(double v, double n, double sinphase, double cosphase, size_t i3x1, size_t i3x2, size_t i3x3, size_t i3x4, double* position, double* energy_accumulate, double* force, double* hessian, double* dvec, double* hdvec) {
  /* Compact AD-based dihedral energy/gradient/Hessian */
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

  Dual2 v1[3];
  Dual2 v2[3];
  Dual2 v3[3];
  v1[0] = dual_var(x2 - x1, 0);
  v1[1] = dual_var(y2 - y1, 1);
  v1[2] = dual_var(z2 - z1, 2);
  v2[0] = dual_var(x3 - x2, 3);
  v2[1] = dual_var(y3 - y2, 4);
  v2[2] = dual_var(z3 - z2, 5);
  v3[0] = dual_var(x4 - x3, 6);
  v3[1] = dual_var(y4 - y3, 7);
  v3[2] = dual_var(z4 - z3, 8);

  Dual2 c1[3];
  Dual2 c2[3];
  dual_cross(v1, v2, c1);
  dual_cross(v2, v3, c2);

  Dual2 v2_sq = dual_dot(v2, v2);
  Dual2 v2_len = dual_sqrt(v2_sq);
  Dual2 dot_v1_c2 = dual_dot(v1, c2);
  Dual2 t1 = dual_mul(v2_len, dot_v1_c2);
  Dual2 t2 = dual_dot(c1, c2);
  Dual2 phi_dual = dual_atan2(t1, t2);
  double phi = phi_dual.val;

  double nphi = n * phi;
  double sin_nphi = sin(nphi);
  double cos_nphi = cos(nphi);
  double cos_angle = (cosphase * cos_nphi) + (sinphase * sin_nphi);
  double sin_angle = (cosphase * sin_nphi) - (sinphase * cos_nphi);
  double energy = v * (1.0 + cos_angle);
  *energy_accumulate += energy;

  double dE_dphi = v * n * (-sin_angle);
  double d2E_dphi2 = v * n * n * (-cos_angle);

  double g_v[NDIM];
  for (int i = 0; i < NDIM; ++i) {
    g_v[i] = dE_dphi * phi_dual.g[i];
  }

  double h_v[NDIM][NDIM];
  for (int i = 0; i < NDIM; ++i) {
    for (int j = 0; j < NDIM; ++j) {
      h_v[i][j] = d2E_dphi2 * phi_dual.g[i] * phi_dual.g[j] + dE_dphi * phi_dual.h[i][j];
    }
  }

  /* Map derivatives from v1/v2/v3 space to atom coordinates */
  static const double C[12][9] = {
    { -1,  0,  0,  0,  0,  0,  0,  0,  0 },
    {  0, -1,  0,  0,  0,  0,  0,  0,  0 },
    {  0,  0, -1,  0,  0,  0,  0,  0,  0 },
    {  1,  0,  0, -1,  0,  0,  0,  0,  0 },
    {  0,  1,  0,  0, -1,  0,  0,  0,  0 },
    {  0,  0,  1,  0,  0, -1,  0,  0,  0 },
    {  0,  0,  0,  1,  0,  0, -1,  0,  0 },
    {  0,  0,  0,  0,  1,  0,  0, -1,  0 },
    {  0,  0,  0,  0,  0,  1,  0,  0, -1 },
    {  0,  0,  0,  0,  0,  0,  1,  0,  0 },
    {  0,  0,  0,  0,  0,  0,  0,  1,  0 },
    {  0,  0,  0,  0,  0,  0,  0,  0,  1 }
  };

  double g_atoms[12] = {0};
  for (int i = 0; i < 12; ++i) {
    for (int k = 0; k < NDIM; ++k) {
      g_atoms[i] += C[i][k] * g_v[k];
    }
  }

  double h_atoms[12][12];
  for (int i = 0; i < 12; ++i) {
    for (int j = 0; j < 12; ++j) {
      h_atoms[i][j] = 0.0;
    }
  }

  for (int i = 0; i < 12; ++i) {
    for (int j = i; j < 12; ++j) {
      double val = 0.0;
      for (int a = 0; a < NDIM; ++a) {
        double cia = C[i][a];
        if (cia == 0.0) continue;
        for (int b = 0; b < NDIM; ++b) {
          double cjb = C[j][b];
          if (cjb == 0.0) continue;
          val += cia * h_v[a][b] * cjb;
        }
      }
      h_atoms[i][j] = val;
      h_atoms[j][i] = val;
    }
  }

  /* Accumulate gradient without branches */
  size_t base_idx[12] = { i3x1, i3x1, i3x1, i3x2, i3x2, i3x2, i3x3, i3x3, i3x3, i3x4, i3x4, i3x4 };
  int axis_idx[12]  = { 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2 };
  for (int i = 0; i < 12; ++i) {
    KernelGradientAcc(base_idx[i], axis_idx[i], g_atoms[i]);
  }

  /* Accumulate Hessian with pre-baked pairs (no ifs in the inner loop) */
  static const unsigned char off_i[66] = {
    0,0,0,0,0,0,0,0,0,0,0,
    1,1,1,1,1,1,1,1,1,1,
    2,2,2,2,2,2,2,2,2,
    3,3,3,3,3,3,3,3,
    4,4,4,4,4,4,4,
    5,5,5,5,5,5,
    6,6,6,6,6,
    7,7,7,7,
    8,8,8,
    9,9,
    10
  };
  static const unsigned char off_j[66] = {
    1,2,3,4,5,6,7,8,9,10,11,
    2,3,4,5,6,7,8,9,10,11,
    3,4,5,6,7,8,9,10,11,
    4,5,6,7,8,9,10,11,
    5,6,7,8,9,10,11,
    6,7,8,9,10,11,
    7,8,9,10,11,
    8,9,10,11,
    9,10,11,
    10,11,
    11
  };

  for (int i = 0; i < 12; ++i) {
    KernelDiagHessAcc(base_idx[i], axis_idx[i], base_idx[i], axis_idx[i], h_atoms[i][i]);
  }
  for (int k = 0; k < 66; ++k) {
    int i = off_i[k];
    int j = off_j[k];
    KernelOffDiagHessAcc(base_idx[i], axis_idx[i], base_idx[j], axis_idx[j], h_atoms[i][j]);
  }
}
