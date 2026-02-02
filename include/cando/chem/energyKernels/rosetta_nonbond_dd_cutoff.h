#pragma once
#include <math.h>

struct rosetta_nonbond_term {
  double sigma;
  double epsilon;
  double mij;
  double bij;
  double aa;
  double bb;
  double cc;
  double dd;
  int i3x1;
  int i3x2;
  rosetta_nonbond_term() {}
  rosetta_nonbond_term(const rosetta_nonbond_parameters& params, double a, double c, int i3x1, int i3x2) {
    this->i3x1 = i3x1;
    this->i3x2 = i3x2;
    double sigma6 = (2.00000000000000000e+0 * a * (1.0 / (c)));
    double sigma = pow(sigma6, (1.00000000000000000e+0 * (1.0 / (6.00000000000000000e+0))));
    this->sigma = sigma;
    double epsilon = (c * c * (1.0 / ((4.00000000000000000e+0 * a))));
    this->epsilon = epsilon;
    double rrep = (6.00000023841857900e-1 * sigma);
    double sr_rep = (sigma * (1.0 / (rrep)));
    double sr2_rep = (sr_rep * sr_rep);
    double sr6_rep = (sr2_rep * sr2_rep * sr2_rep);
    double sr12_rep = (sr6_rep * sr6_rep);
    double rep_rrep = (epsilon * (1.00000000000000000e+0 + sr12_rep + (-2.00000000000000000e+0 * sr6_rep)));
    double drep_rrep = (epsilon * ((-1.20000000000000000e+1 * sr12_rep) + (1.20000000000000000e+1 * sr6_rep)) * (1.0 / (rrep)));
    double mij = drep_rrep;
    this->mij = mij;
    double bij = (rep_rrep + (-((mij * rrep))));
    this->bij = bij;
    double sr_sw = (sigma * (1.0 / (params.r_switch)));
    double sr2_sw = (sr_sw * sr_sw);
    double sr6_sw = (sr2_sw * sr2_sw * sr2_sw);
    double sr12_sw = (sr6_sw * sr6_sw);
    double atr_sw = (epsilon * (sr12_sw + (-2.00000000000000000e+0 * sr6_sw)));
    double d_atr_sw = (epsilon * ((-1.20000000000000000e+1 * sr12_sw) + (1.20000000000000000e+1 * sr6_sw)) * (1.0 / (params.r_switch)));
    double dr = (params.r_cut + (-(params.r_switch)));
    double aa_shift = (((2.00000000000000000e+0 * atr_sw) + (dr * d_atr_sw)) * (1.0 / ((dr * dr * dr))));
    double bb_shift = (((-2.00000000000000000e+0 * dr * d_atr_sw) + (-3.00000000000000000e+0 * atr_sw)) * (1.0 / ((dr * dr))));
    double cc_shift = d_atr_sw;
    double dd_shift = atr_sw;
    double aa = aa_shift;
    this->aa = aa;
    double bb = (bb_shift + (-3.00000000000000000e+0 * aa_shift * params.r_switch));
    this->bb = bb;
    double cc = (cc_shift + (-2.00000000000000000e+0 * bb_shift * params.r_switch) + (3.00000000000000000e+0 * aa_shift * params.r_switch * params.r_switch));
    this->cc = cc;
    double dd = (dd_shift + (bb_shift * params.r_switch * params.r_switch) + (-((aa_shift * params.r_switch * params.r_switch * params.r_switch))) + (-((cc_shift * params.r_switch))));
    this->dd = dd;
  }
};
