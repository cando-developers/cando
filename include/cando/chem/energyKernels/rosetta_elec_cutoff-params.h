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
    : elec_weight(1.0), eps_core(6.0), eps_solvent(80.0), rmin(1.45), rlow(1.85), rhi(4.5), rcut(5.5), rpairlist(7.5)
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

