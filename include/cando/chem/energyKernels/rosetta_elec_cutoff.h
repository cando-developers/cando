#pragma once
#include <math.h>

struct rosetta_elec_term {
  double kqq;
  double e_rmin;
  double aa_low;
  double bb_low;
  double cc_low;
  double dd_low;
  double aa_high;
  double bb_high;
  double cc_high;
  double dd_high;
  int i3x1;
  int i3x2;
  rosetta_elec_term() {}
  rosetta_elec_term(const rosetta_elec_parameters& params, double kqq, int i3x1, int i3x2) {
    this->kqq = kqq;
    this->i3x1 = i3x1;
    this->i3x2 = i3x2;
    double rx_rmin = (2.50000000000000000e-1 * params.rmin);
    double ex_rmin = exp((-(rx_rmin)));
    double gx_rmin = (ex_rmin * (1.00000000000000000e+0 + rx_rmin + (5.00000000000000000e-1 * rx_rmin * rx_rmin)));
    double ep_rmin = (params.eps_solvent + (gx_rmin * (params.eps_core + (-(params.eps_solvent)))));
    double e_rmin = (kqq * (1.0 / (ep_rmin)));
    this->e_rmin = e_rmin;
    double rx_rlow = (2.50000000000000000e-1 * params.rlow);
    double ex_rlow = exp((-(rx_rlow)));
    double gx_rlow = (ex_rlow * (1.00000000000000000e+0 + rx_rlow + (5.00000000000000000e-1 * rx_rlow * rx_rlow)));
    double ep_rlow = (params.eps_solvent + (gx_rlow * (params.eps_core + (-(params.eps_solvent)))));
    double e_rlow = (kqq * (1.0 / (ep_rlow)));
    double gx_p_rlow = (-5.00000000000000000e-1 * ex_rlow * rx_rlow * rx_rlow);
    double ep_p_rlow = (2.50000000000000000e-1 * gx_p_rlow * (params.eps_core + (-(params.eps_solvent))));
    double de_rlow = (ep_p_rlow * (1.0 / ((ep_rlow * ep_rlow))) * (-(kqq)));
    double rx_rhi = (2.50000000000000000e-1 * params.rhi);
    double ex_rhi = exp((-(rx_rhi)));
    double gx_rhi = (ex_rhi * (1.00000000000000000e+0 + rx_rhi + (5.00000000000000000e-1 * rx_rhi * rx_rhi)));
    double ep_rhi = (params.eps_solvent + (gx_rhi * (params.eps_core + (-(params.eps_solvent)))));
    double e_rhi = (kqq * (1.0 / (ep_rhi)));
    double gx_p_rhi = (-5.00000000000000000e-1 * ex_rhi * rx_rhi * rx_rhi);
    double ep_p_rhi = (2.50000000000000000e-1 * gx_p_rhi * (params.eps_core + (-(params.eps_solvent))));
    double de_rhi = (ep_p_rhi * (1.0 / ((ep_rhi * ep_rhi))) * (-(kqq)));
    double dr_low = (params.rlow + (-(params.rmin)));
    double dr_high = (params.rcut + (-(params.rhi)));
    double aa_low_shifted = (((-2.00000000000000000e+0 * e_rlow) + (2.00000000000000000e+0 * e_rmin) + (de_rlow * dr_low)) * (1.0 / ((dr_low * dr_low * dr_low))));
    double bb_low_shifted = (((-3.00000000000000000e+0 * e_rmin) + (3.00000000000000000e+0 * e_rlow) + (-((de_rlow * dr_low)))) * (1.0 / ((dr_low * dr_low))));
    double cc_low_shifted = 0.00000000000000000e+0;
    double dd_low_shifted = e_rmin;
    double aa_low = aa_low_shifted;
    this->aa_low = aa_low;
    double bb_low = (bb_low_shifted + (-3.00000000000000000e+0 * aa_low_shifted * params.rmin));
    this->bb_low = bb_low;
    double cc_low = (cc_low_shifted + (-2.00000000000000000e+0 * bb_low_shifted * params.rmin) + (3.00000000000000000e+0 * aa_low_shifted * params.rmin * params.rmin));
    this->cc_low = cc_low;
    double dd_low = (dd_low_shifted + (bb_low_shifted * params.rmin * params.rmin) + (-((aa_low_shifted * params.rmin * params.rmin * params.rmin))) + (-((cc_low_shifted * params.rmin))));
    this->dd_low = dd_low;
    double aa_high_shifted = (((2.00000000000000000e+0 * e_rhi) + (de_rhi * dr_high)) * (1.0 / ((dr_high * dr_high * dr_high))));
    double bb_high_shifted = (((-2.00000000000000000e+0 * de_rhi * dr_high) + (-3.00000000000000000e+0 * e_rhi)) * (1.0 / ((dr_high * dr_high))));
    double cc_high_shifted = de_rhi;
    double dd_high_shifted = e_rhi;
    double aa_high = aa_high_shifted;
    this->aa_high = aa_high;
    double bb_high = (bb_high_shifted + (-3.00000000000000000e+0 * aa_high_shifted * params.rhi));
    this->bb_high = bb_high;
    double cc_high = (cc_high_shifted + (-2.00000000000000000e+0 * bb_high_shifted * params.rhi) + (3.00000000000000000e+0 * aa_high_shifted * params.rhi * params.rhi));
    this->cc_high = cc_high;
    double dd_high = (dd_high_shifted + (bb_high_shifted * params.rhi * params.rhi) + (-((aa_high_shifted * params.rhi * params.rhi * params.rhi))) + (-((cc_high_shifted * params.rhi))));
    this->dd_high = dd_high;
  }
};
