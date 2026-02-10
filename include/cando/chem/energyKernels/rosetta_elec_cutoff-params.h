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
  void do_apply(SetupAccumulator& setupAcc) {
    setupAcc.maybe_apply(INTERN_(kw,elec_weight),this->elec_weight);
    setupAcc.maybe_apply(INTERN_(kw,eps_core),this->eps_core);
    setupAcc.maybe_apply(INTERN_(kw,eps_solvent),this->eps_solvent);
    setupAcc.maybe_apply(INTERN_(kw,rmin),this->rmin);
    setupAcc.maybe_apply(INTERN_(kw,rlow),this->rlow);
    setupAcc.maybe_apply(INTERN_(kw,rhi),this->rhi);
    setupAcc.maybe_apply(INTERN_(kw,rcut),this->rcut);
    setupAcc.maybe_apply(INTERN_(kw,rpairlist),this->rpairlist);
  }
  void fields(core::Record_sp node) {
    node->field(INTERN_(kw, elec_weight), this->elec_weight);
    node->field(INTERN_(kw, eps_core), this->eps_core);
    node->field(INTERN_(kw, eps_solvent), this->eps_solvent);
    node->field(INTERN_(kw, rmin), this->rmin);
    node->field(INTERN_(kw, rlow), this->rlow);
    node->field(INTERN_(kw, rhi), this->rhi);
    node->field(INTERN_(kw, rcut), this->rcut);
    node->field(INTERN_(kw, rpairlist), this->rpairlist);
  }
};

