#pragma once
#include <math.h>

struct rosetta_nonbond_parameters {
  double rep_weight;
  double r_switch;
  double r_cut;
  double r_pairlist;
  rosetta_nonbond_parameters()
    : rep_weight(1.00000000000000000e+0), r_switch(4.50000000000000000e+0), r_cut(6.00000000000000000e+0), r_pairlist(9.00000000000000000e+0)
  {}
};

