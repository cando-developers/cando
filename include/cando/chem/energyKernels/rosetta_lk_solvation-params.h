#pragma once
#include <math.h>

struct rosetta_lk_solvation_parameters {
  double c0;
  double c1;
  double r_solv_low;
  double r_solv_high;
  double rpairlist;
  rosetta_lk_solvation_parameters()
    : c0(3.00000011920928960e-1), c1(2.00000002980232240e-1), r_solv_low(4.50000000000000000e+0), r_solv_high(6.00000000000000000e+0), rpairlist(8.00000000000000000e+0)
  {}
  void do_apply(SetupAccumulator& setupAcc) {
    setupAcc.maybe_apply(INTERN_(kw,c0),this->c0);
    setupAcc.maybe_apply(INTERN_(kw,c1),this->c1);
    setupAcc.maybe_apply(INTERN_(kw,r_solv_low),this->r_solv_low);
    setupAcc.maybe_apply(INTERN_(kw,r_solv_high),this->r_solv_high);
    setupAcc.maybe_apply(INTERN_(kw,rpairlist),this->rpairlist);
  }
  void fields(core::Record_sp node) {
    node->field(INTERN_(kw, c0), this->c0);
    node->field(INTERN_(kw, c1), this->c1);
    node->field(INTERN_(kw, r_solv_low), this->r_solv_low);
    node->field(INTERN_(kw, r_solv_high), this->r_solv_high);
    node->field(INTERN_(kw, rpairlist), this->rpairlist);
  }
};

