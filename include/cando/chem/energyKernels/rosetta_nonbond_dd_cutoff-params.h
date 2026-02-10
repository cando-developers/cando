#pragma once
#include <math.h>

struct rosetta_nonbond_parameters {
  double rep_weight;
  double rswitch;
  double rcut;
  double rpairlist;
  rosetta_nonbond_parameters()
    : rep_weight(1.00000000000000000e+0), rswitch(4.50000000000000000e+0), rcut(6.00000000000000000e+0), rpairlist(9.00000000000000000e+0)
  {}
  void do_apply(SetupAccumulator& setupAcc) {
    setupAcc.maybe_apply(INTERN_(kw,rep_weight),this->rep_weight);
    setupAcc.maybe_apply(INTERN_(kw,rswitch),this->rswitch);
    setupAcc.maybe_apply(INTERN_(kw,rcut),this->rcut);
    setupAcc.maybe_apply(INTERN_(kw,rpairlist),this->rpairlist);
  }
  void fields(core::Record_sp node) {
    node->field(INTERN_(kw, rep_weight), this->rep_weight);
    node->field(INTERN_(kw, rswitch), this->rswitch);
    node->field(INTERN_(kw, rcut), this->rcut);
    node->field(INTERN_(kw, rpairlist), this->rpairlist);
  }
};

