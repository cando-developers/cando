
template <typename HESSIAN>
struct Nonbond_Dd_Cutoff {
  static constexpr size_t PositionSize = 6;

  static std::string description() { return "gpt-nonbond_dd_cutoff"; };
  static inline void switch_terms(double r, double r_switch, double inv_range, double* s, double* ds_dr, double* d2s_dr2) {
    double t = (r - r_switch) * inv_range;
    double t2 = t * t;
    double t3 = t2 * t;
    *s = 1.0 + t3 * (-10.0 + t * (15.0 - 6.0 * t));
    double ds_dt = t2 * (-30.0 + t * (60.0 - 30.0 * t));
    *ds_dr = ds_dt * inv_range;
    if (d2s_dr2) {
      double d2s_dt2 = t * (-60.0 + t * (180.0 - 120.0 * t));
      *d2s_dr2 = d2s_dt2 * inv_range * inv_range;
    }
  }

  static inline double energy(double a, double b, double qq, double invdd, double r_switch, double r_switch2, double r_cut, double r_cut2, double inv_range, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
    DOUBLE x1 = position[i3x1 + 0];
    DOUBLE y1 = position[i3x1 + 1];
    DOUBLE z1 = position[i3x1 + 2];
    DOUBLE x2 = position[i3x2 + 0];
    DOUBLE y2 = position[i3x2 + 1];
    DOUBLE z2 = position[i3x2 + 2];

    double dx = x1 - x2;
    double dy = y1 - y2;
    double dz = z1 - z2;
    double r2 = (dx * dx) + (dy * dy) + (dz * dz);
    if (r2 >= r_cut2) {
      return 0.0;
    }

    double invr = 1.0 / sqrt(r2);
    double invr2 = invr * invr;
    double invr6 = invr2 * invr2 * invr2;
    double invr12 = invr6 * invr6;
    double k = invdd * qq;
    double e_base = (a * invr12 - b * invr6) + (k * invr2);

    if (r2 < r_switch2) {
      *energy_accumulate += e_base;
      return e_base;
    }

    double r = r2 * invr;
    double t = (r - r_switch) * inv_range;
    double t2 = t * t;
    double t3 = t2 * t;
    double s = 1.0 + t3 * (-10.0 + t * (15.0 - 6.0 * t));
    *energy_accumulate += e_base * s;
    return e_base * s;
  }

















  
  static inline void energy_fd(double a, double b, double qq, double invdd, double r_switch, double r_switch2, double r_cut, double r_cut2, double inv_range, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
  {
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, energy_accumulate, force, hessian, dvec, hdvec);
  }

  static inline double gradient(double a, double b, double qq, double invdd, double r_switch, double r_switch2, double r_cut, double r_cut2, double inv_range, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
    DOUBLE x1 = position[i3x1 + 0];
    DOUBLE y1 = position[i3x1 + 1];
    DOUBLE z1 = position[i3x1 + 2];
    DOUBLE x2 = position[i3x2 + 0];
    DOUBLE y2 = position[i3x2 + 1];
    DOUBLE z2 = position[i3x2 + 2];

    double dx = x1 - x2;
    double dy = y1 - y2;
    double dz = z1 - z2;
    double r2 = (dx * dx) + (dy * dy) + (dz * dz);
    if (r2 >= r_cut2) {
      return 0.0;
    }

    double invr = 1.0 / sqrt(r2);
    double invr2 = invr * invr;
    double invr3 = invr2 * invr;
    double invr6 = invr2 * invr2 * invr2;
    double invr12 = invr6 * invr6;

    double k = invdd * qq;
    double e_lj = (a * invr12) - (b * invr6);
    double e_coul = k * invr2;
    double e_base = e_lj + e_coul;
    double de_base_dr = ((-12.0 * a * invr12) + (6.0 * b * invr6)) * invr - (2.0 * k * invr3);

    double s = 1.0;
    double ds_dr = 0.0;
    if (r2 >= r_switch2) {
      double r = r2 * invr;
      switch_terms(r, r_switch, inv_range, &s, &ds_dr, nullptr);
    }

    double energy = e_base * s;
    *energy_accumulate += energy;

    double de_dr = de_base_dr * s + e_base * ds_dr;
    double g = de_dr * invr;

    double g_x1 = g * dx;
    double g_y1 = g * dy;
    double g_z1 = g * dz;
    double g_x2 = -g_x1;
    double g_y2 = -g_y1;
    double g_z2 = -g_z1;

    KernelGradientAcc(i3x1, 0, g_x1);
    KernelGradientAcc(i3x1, 1, g_y1);
    KernelGradientAcc(i3x1, 2, g_z1);
    KernelGradientAcc(i3x2, 0, g_x2);
    KernelGradientAcc(i3x2, 1, g_y2);
    KernelGradientAcc(i3x2, 2, g_z2);
    return energy;
  }

  static inline void  gradient_fd(double a, double b, double qq, double invdd, double r_switch, double r_switch2, double r_cut, double r_cut2, double inv_range, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
  {
    const double h = 1.0e-5;
    const double inv2h = 1.0 / (2.0 * h);
    double e0 = 0.0;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e0, 0, 0, 0, 0);
    if (energy_accumulate) {
      *energy_accumulate += e0;
    }

    const size_t coord_base[6] = {i3x1, i3x1, i3x1, i3x2, i3x2, i3x2};
    const size_t coord_off[6] = {0, 1, 2, 0, 1, 2};

    for (size_t idx = 0; idx < 6; ++idx) {
      size_t base = coord_base[idx];
      size_t off = coord_off[idx];
      double saved = position[base + off];
      double e_plus = 0.0;
      double e_minus = 0.0;
      position[base + off] = saved + h;
      energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
      position[base + off] = saved - h;
      energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
      position[base + off] = saved;
      double d = (e_plus - e_minus) * inv2h;
      KernelGradientAcc(base, off, d);
    }
  }

  static inline double hessian(double a, double b, double qq, double invdd, double r_switch, double r_switch2, double r_cut, double r_cut2, double inv_range, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec) {
    DOUBLE x1 = position[i3x1 + 0];
    DOUBLE y1 = position[i3x1 + 1];
    DOUBLE z1 = position[i3x1 + 2];
    DOUBLE x2 = position[i3x2 + 0];
    DOUBLE y2 = position[i3x2 + 1];
    DOUBLE z2 = position[i3x2 + 2];

    double dx = x1 - x2;
    double dy = y1 - y2;
    double dz = z1 - z2;
    double r2 = (dx * dx) + (dy * dy) + (dz * dz);
    if (r2 >= r_cut2) {
      return 0.0;
    }

    double invr = 1.0 / sqrt(r2);
    double invr2 = invr * invr;
    double invr3 = invr2 * invr;
    double invr4 = invr2 * invr2;
    double invr6 = invr4 * invr2;
    double invr8 = invr6 * invr2;
    double invr12 = invr6 * invr6;
    double invr14 = invr12 * invr2;

    double k = invdd * qq;
    double e_lj = (a * invr12) - (b * invr6);
    double e_coul = k * invr2;
    double e_base = e_lj + e_coul;
    double de_base_dr = ((-12.0 * a * invr12) + (6.0 * b * invr6)) * invr - (2.0 * k * invr3);
    double d2e_base_dr2 = (156.0 * a * invr14) - (42.0 * b * invr8) + (6.0 * k * invr4);

    double s = 1.0;
    double ds_dr = 0.0;
    double d2s_dr2 = 0.0;
    if (r2 >= r_switch2) {
      double r = r2 * invr;
      switch_terms(r, r_switch, inv_range, &s, &ds_dr, &d2s_dr2);
    }

    double energy = e_base * s;
    *energy_accumulate += energy;

    double de_dr = de_base_dr * s + e_base * ds_dr;
    double d2e_dr2 = d2e_base_dr2 * s + (2.0 * de_base_dr * ds_dr) + (e_base * d2s_dr2);

    double g = de_dr * invr;
    double h = (d2e_dr2 - (de_dr * invr)) * invr2;

    double g_x1 = g * dx;
    double g_y1 = g * dy;
    double g_z1 = g * dz;
    double g_x2 = -g_x1;
    double g_y2 = -g_y1;
    double g_z2 = -g_z1;

    KernelGradientAcc(i3x1, 0, g_x1);
    KernelGradientAcc(i3x1, 1, g_y1);
    KernelGradientAcc(i3x1, 2, g_z1);
    KernelGradientAcc(i3x2, 0, g_x2);
    KernelGradientAcc(i3x2, 1, g_y2);
    KernelGradientAcc(i3x2, 2, g_z2);

    double hxx = h * dx * dx + g;
    double hxy = h * dx * dy;
    double hxz = h * dx * dz;
    double hyy = h * dy * dy + g;
    double hyz = h * dy * dz;
    double hzz = h * dz * dz + g;

    double h_x1_x1 = hxx;
    double h_x1_y1 = hxy;
    double h_x1_z1 = hxz;
    double h_y1_y1 = hyy;
    double h_y1_z1 = hyz;
    double h_z1_z1 = hzz;

    double h_x1_x2 = -hxx;
    double h_x1_y2 = -hxy;
    double h_x1_z2 = -hxz;
    double h_y1_x2 = -hxy;
    double h_y1_y2 = -hyy;
    double h_y1_z2 = -hyz;
    double h_z1_x2 = -hxz;
    double h_z1_y2 = -hyz;
    double h_z1_z2 = -hzz;

    double h_x2_x2 = hxx;
    double h_x2_y2 = hxy;
    double h_x2_z2 = hxz;
    double h_y2_y2 = hyy;
    double h_y2_z2 = hyz;
    double h_z2_z2 = hzz;

    KernelHessDiagAcc(PositionSize, hessian, dvec, hdvec, i3x1, 0, i3x1, 0, h_x1_x1);
    KernelHessOffDiagAcc(PositionSize, hessian, dvec, hdvec, i3x1, 0, i3x1, 1, h_x1_y1);
    KernelHessOffDiagAcc(PositionSize, hessian, dvec, hdvec, i3x1, 0, i3x1, 2, h_x1_z1);
    KernelHessOffDiagAcc(PositionSize, hessian, dvec, hdvec, i3x1, 0, i3x2, 0, h_x1_x2);
    KernelHessOffDiagAcc(PositionSize, hessian, dvec, hdvec, i3x1, 0, i3x2, 1, h_x1_y2);
    KernelHessOffDiagAcc(PositionSize, hessian, dvec, hdvec, i3x1, 0, i3x2, 2, h_x1_z2);
    KernelHessDiagAcc(PositionSize, hessian, dvec, hdvec, i3x1, 1, i3x1, 1, h_y1_y1);
    KernelHessOffDiagAcc(PositionSize, hessian, dvec, hdvec, i3x1, 1, i3x1, 2, h_y1_z1);
    KernelHessOffDiagAcc(PositionSize, hessian, dvec, hdvec, i3x1, 1, i3x2, 0, h_y1_x2);
    KernelHessOffDiagAcc(PositionSize, hessian, dvec, hdvec, i3x1, 1, i3x2, 1, h_y1_y2);
    KernelHessOffDiagAcc(PositionSize, hessian, dvec, hdvec, i3x1, 1, i3x2, 2, h_y1_z2);
    KernelHessDiagAcc(PositionSize, hessian, dvec, hdvec, i3x1, 2, i3x1, 2, h_z1_z1);
    KernelHessOffDiagAcc(PositionSize, hessian, dvec, hdvec, i3x1, 2, i3x2, 0, h_z1_x2);
    KernelHessOffDiagAcc(PositionSize, hessian, dvec, hdvec, i3x1, 2, i3x2, 1, h_z1_y2);
    KernelHessOffDiagAcc(PositionSize, hessian, dvec, hdvec, i3x1, 2, i3x2, 2, h_z1_z2);
    KernelHessDiagAcc(PositionSize, hessian, dvec, hdvec, i3x2, 0, i3x2, 0, h_x2_x2);
    KernelHessOffDiagAcc(PositionSize, hessian, dvec, hdvec, i3x2, 0, i3x2, 1, h_x2_y2);
    KernelHessOffDiagAcc(PositionSize, hessian, dvec, hdvec, i3x2, 0, i3x2, 2, h_x2_z2);
    KernelHessDiagAcc(PositionSize, hessian, dvec, hdvec, i3x2, 1, i3x2, 1, h_y2_y2);
    KernelHessOffDiagAcc(PositionSize, hessian, dvec, hdvec, i3x2, 1, i3x2, 2, h_y2_z2);
    KernelHessDiagAcc(PositionSize, hessian, dvec, hdvec, i3x2, 2, i3x2, 2, h_z2_z2);
    return energy;
  }

  static inline void hessian_fd(double a, double b, double qq, double invdd, double r_switch, double r_switch2, double r_cut, double r_cut2, double inv_range, size_t i3x1, size_t i3x2, double* position, double* energy_accumulate, double* force, HESSIAN hessian, double* dvec, double* hdvec)
  {
    const double h = 1.0e-5;
    const double inv2h = 1.0 / (2.0 * h);
    const double invh2 = 1.0 / (h * h);
    double e0 = 0.0;
    energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e0, 0, 0, 0, 0);
    if (energy_accumulate) {
      *energy_accumulate += e0;
    }

    const size_t coord_base[6] = {i3x1, i3x1, i3x1, i3x2, i3x2, i3x2};
    const size_t coord_off[6] = {0, 1, 2, 0, 1, 2};

    for (size_t idx = 0; idx < 6; ++idx) {
      size_t base = coord_base[idx];
      size_t off = coord_off[idx];
      double saved = position[base + off];
      double e_plus = 0.0;
      double e_minus = 0.0;
      position[base + off] = saved + h;
      energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_plus, 0, 0, 0, 0);
      position[base + off] = saved - h;
      energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_minus, 0, 0, 0, 0);
      position[base + off] = saved;
      double d = (e_plus - e_minus) * inv2h;
      KernelGradientAcc(base, off, d);
      double hval = (e_plus + e_minus - (2.0 * e0)) * invh2;
      KernelHessDiagAcc(PositionSize, hessian, dvec, hdvec, base, off, base, off, hval);
    }

    for (size_t i = 1; i < 6; ++i) {
      size_t base_i = coord_base[i];
      size_t off_i = coord_off[i];
      double saved_i = position[base_i + off_i];
      for (size_t j = 0; j < i; ++j) {
        size_t base_j = coord_base[j];
        size_t off_j = coord_off[j];
        double saved_j = position[base_j + off_j];
        double e_pp = 0.0;
        double e_pm = 0.0;
        double e_mp = 0.0;
        double e_mm = 0.0;
        position[base_i + off_i] = saved_i + h;
        position[base_j + off_j] = saved_j + h;
        energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pp, 0, 0, 0, 0);
        position[base_j + off_j] = saved_j - h;
        energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_pm, 0, 0, 0, 0);
        position[base_i + off_i] = saved_i - h;
        position[base_j + off_j] = saved_j + h;
        energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mp, 0, 0, 0, 0);
        position[base_j + off_j] = saved_j - h;
        energy(a, b, qq, invdd, r_switch, r_switch2, r_cut, r_cut2, inv_range, i3x1, i3x2, position, &e_mm, 0, 0, 0, 0);
        position[base_i + off_i] = saved_i;
        position[base_j + off_j] = saved_j;
        double hval = (e_pp - e_pm - e_mp + e_mm) * (0.25 * invh2);
        KernelHessOffDiagAcc(PositionSize, hessian, dvec, hdvec, base_i, off_i, base_j, off_j, hval);
      }
    }
  }
};
