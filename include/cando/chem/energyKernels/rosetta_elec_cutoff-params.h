#pragma once
#include <math.h>

struct rosetta_elec_parameters {
  double elec_weight;
  double eps_core;
  double eps_solvent;
  double rmin;
  double rlow;
  double rhi;
  double rcut;
  double rpairlist;
  rosetta_elec_parameters()
    : elec_weight(1.00000000000000000e+0), eps_core(6.00000000000000000e+0), eps_solvent(8.00000000000000000e+1), rmin(1.45000004768371550e+0), rlow(1.85000002384185800e+0), rhi(4.50000000000000000e+0), rcut(5.50000000000000000e+0), rpairlist(7.50000000000000000e+0)
  {}
};

