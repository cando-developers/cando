#pragma once
#include <math.h>

struct rosetta_lk_solvation_term {
  double lj_radius_i;
  double lj_radius_j;
  double sigmaij;
  double r0_low;
  double r1_low;
  double r2_high;
  double r3_high;
  double dr_low;
  double dr_high;
  double inv2pi32;
  double lambda2_i;
  double lambda2_j;
  double pref_i_j;
  double pref_j_i;
  double f_contact_i;
  double f_contact_j;
  double f_r1_i;
  double f_r1_j;
  double df_r1_i;
  double df_r1_j;
  double f_r2_i;
  double f_r2_j;
  double df_r2_i;
  double df_r2_j;
  double aa_low_i;
  double bb_low_i;
  double cc_low_i;
  double dd_low_i;
  double aa_low_j;
  double bb_low_j;
  double cc_low_j;
  double dd_low_j;
  double aa_high_i;
  double bb_high_i;
  double cc_high_i;
  double dd_high_i;
  double aa_high_j;
  double bb_high_j;
  double cc_high_j;
  double dd_high_j;
  int i3x1;
  int i3x2;
  rosetta_lk_solvation_term() {}
  rosetta_lk_solvation_term(const rosetta_lk_solvation_parameters& params, double lk_dgfree_i, double lk_lambda_i, double lj_radius_i, double lk_volume_i, double lk_dgfree_j, double lk_lambda_j, double lj_radius_j, double lk_volume_j, int i3x1, int i3x2) {
    this->lj_radius_i = lj_radius_i;
    this->lj_radius_j = lj_radius_j;
    this->i3x1 = i3x1;
    this->i3x2 = i3x2;
    double sigmaij = (lj_radius_i + lj_radius_j);
    this->sigmaij = sigmaij;
    double r0_low = (sigmaij + (-(params.c0)));
    this->r0_low = r0_low;
    double r1_low = (params.c1 + sigmaij);
    this->r1_low = r1_low;
    double r2_high = params.r_solv_low;
    this->r2_high = r2_high;
    double r3_high = params.r_solv_high;
    this->r3_high = r3_high;
    double dr_low = (r1_low + (-(r0_low)));
    this->dr_low = dr_low;
    double dr_high = (r3_high + (-(r2_high)));
    this->dr_high = dr_high;
    double inv2pi32 = (1.00000000000000000e+0 * (1.0 / ((2.00000000000000000e+0 * pi * sqrt(pi)))));
    this->inv2pi32 = inv2pi32;
    double lambda2_i = (lk_lambda_i * lk_lambda_i);
    this->lambda2_i = lambda2_i;
    double lambda2_j = (lk_lambda_j * lk_lambda_j);
    this->lambda2_j = lambda2_j;
    double pref_i_j = (1.00000000000000000e+0 * lk_dgfree_i * lk_volume_j * (1.0 / (lk_lambda_i)) * (-(inv2pi32)));
    this->pref_i_j = pref_i_j;
    double pref_j_i = (1.00000000000000000e+0 * lk_dgfree_j * lk_volume_i * (1.0 / (lk_lambda_j)) * (-(inv2pi32)));
    this->pref_j_i = pref_j_i;
    double delta_contact_i = (sigmaij + (-(lj_radius_i)));
    double exp_contact_i = exp(((1.0 / (lambda2_i)) * (-((delta_contact_i * delta_contact_i)))));
    double f_contact_i = (exp_contact_i * pref_i_j * (1.0 / ((sigmaij * sigmaij))));
    this->f_contact_i = f_contact_i;
    double delta_contact_j = (sigmaij + (-(lj_radius_j)));
    double exp_contact_j = exp(((1.0 / (lambda2_j)) * (-((delta_contact_j * delta_contact_j)))));
    double f_contact_j = (exp_contact_j * pref_j_i * (1.0 / ((sigmaij * sigmaij))));
    this->f_contact_j = f_contact_j;
    double delta_r1_i = (r1_low + (-(lj_radius_i)));
    double exp_r1_i = exp(((1.0 / (lambda2_i)) * (-((delta_r1_i * delta_r1_i)))));
    double f_r1_i = (exp_r1_i * pref_i_j * (1.0 / ((r1_low * r1_low))));
    this->f_r1_i = f_r1_i;
    double df_r1_i = (-2.00000000000000000e+0 * f_r1_i * ((1.00000000000000000e+0 * (1.0 / (r1_low))) + (delta_r1_i * (1.0 / (lambda2_i)))));
    this->df_r1_i = df_r1_i;
    double delta_r1_j = (r1_low + (-(lj_radius_j)));
    double exp_r1_j = exp(((1.0 / (lambda2_j)) * (-((delta_r1_j * delta_r1_j)))));
    double f_r1_j = (exp_r1_j * pref_j_i * (1.0 / ((r1_low * r1_low))));
    this->f_r1_j = f_r1_j;
    double df_r1_j = (-2.00000000000000000e+0 * f_r1_j * ((1.00000000000000000e+0 * (1.0 / (r1_low))) + (delta_r1_j * (1.0 / (lambda2_j)))));
    this->df_r1_j = df_r1_j;
    double delta_r2_i = (r2_high + (-(lj_radius_i)));
    double exp_r2_i = exp(((1.0 / (lambda2_i)) * (-((delta_r2_i * delta_r2_i)))));
    double f_r2_i = (exp_r2_i * pref_i_j * (1.0 / ((r2_high * r2_high))));
    this->f_r2_i = f_r2_i;
    double df_r2_i = (-2.00000000000000000e+0 * f_r2_i * ((1.00000000000000000e+0 * (1.0 / (r2_high))) + (delta_r2_i * (1.0 / (lambda2_i)))));
    this->df_r2_i = df_r2_i;
    double delta_r2_j = (r2_high + (-(lj_radius_j)));
    double exp_r2_j = exp(((1.0 / (lambda2_j)) * (-((delta_r2_j * delta_r2_j)))));
    double f_r2_j = (exp_r2_j * pref_j_i * (1.0 / ((r2_high * r2_high))));
    this->f_r2_j = f_r2_j;
    double df_r2_j = (-2.00000000000000000e+0 * f_r2_j * ((1.00000000000000000e+0 * (1.0 / (r2_high))) + (delta_r2_j * (1.0 / (lambda2_j)))));
    this->df_r2_j = df_r2_j;
    double aa_low_i_shift = (((-2.00000000000000000e+0 * f_r1_i) + (2.00000000000000000e+0 * f_contact_i) + (df_r1_i * dr_low)) * (1.0 / ((dr_low * dr_low * dr_low))));
    double bb_low_i_shift = (((-3.00000000000000000e+0 * f_contact_i) + (3.00000000000000000e+0 * f_r1_i) + (-((df_r1_i * dr_low)))) * (1.0 / ((dr_low * dr_low))));
    double cc_low_i_shift = 0.00000000000000000e+0;
    double dd_low_i_shift = f_contact_i;
    double aa_low_j_shift = (((-2.00000000000000000e+0 * f_r1_j) + (2.00000000000000000e+0 * f_contact_j) + (df_r1_j * dr_low)) * (1.0 / ((dr_low * dr_low * dr_low))));
    double bb_low_j_shift = (((-3.00000000000000000e+0 * f_contact_j) + (3.00000000000000000e+0 * f_r1_j) + (-((df_r1_j * dr_low)))) * (1.0 / ((dr_low * dr_low))));
    double cc_low_j_shift = 0.00000000000000000e+0;
    double dd_low_j_shift = f_contact_j;
    double aa_low_i = aa_low_i_shift;
    this->aa_low_i = aa_low_i;
    double bb_low_i = (bb_low_i_shift + (-3.00000000000000000e+0 * aa_low_i_shift * r0_low));
    this->bb_low_i = bb_low_i;
    double cc_low_i = (cc_low_i_shift + (-2.00000000000000000e+0 * bb_low_i_shift * r0_low) + (3.00000000000000000e+0 * aa_low_i_shift * r0_low * r0_low));
    this->cc_low_i = cc_low_i;
    double dd_low_i = (dd_low_i_shift + (bb_low_i_shift * r0_low * r0_low) + (-((aa_low_i_shift * r0_low * r0_low * r0_low))) + (-((cc_low_i_shift * r0_low))));
    this->dd_low_i = dd_low_i;
    double aa_low_j = aa_low_j_shift;
    this->aa_low_j = aa_low_j;
    double bb_low_j = (bb_low_j_shift + (-3.00000000000000000e+0 * aa_low_j_shift * r0_low));
    this->bb_low_j = bb_low_j;
    double cc_low_j = (cc_low_j_shift + (-2.00000000000000000e+0 * bb_low_j_shift * r0_low) + (3.00000000000000000e+0 * aa_low_j_shift * r0_low * r0_low));
    this->cc_low_j = cc_low_j;
    double dd_low_j = (dd_low_j_shift + (bb_low_j_shift * r0_low * r0_low) + (-((aa_low_j_shift * r0_low * r0_low * r0_low))) + (-((cc_low_j_shift * r0_low))));
    this->dd_low_j = dd_low_j;
    double aa_high_i_shift = (((2.00000000000000000e+0 * f_r2_i) + (df_r2_i * dr_high)) * (1.0 / ((dr_high * dr_high * dr_high))));
    double bb_high_i_shift = (((-2.00000000000000000e+0 * df_r2_i * dr_high) + (-3.00000000000000000e+0 * f_r2_i)) * (1.0 / ((dr_high * dr_high))));
    double cc_high_i_shift = df_r2_i;
    double dd_high_i_shift = f_r2_i;
    double aa_high_j_shift = (((2.00000000000000000e+0 * f_r2_j) + (df_r2_j * dr_high)) * (1.0 / ((dr_high * dr_high * dr_high))));
    double bb_high_j_shift = (((-2.00000000000000000e+0 * df_r2_j * dr_high) + (-3.00000000000000000e+0 * f_r2_j)) * (1.0 / ((dr_high * dr_high))));
    double cc_high_j_shift = df_r2_j;
    double dd_high_j_shift = f_r2_j;
    double aa_high_i = aa_high_i_shift;
    this->aa_high_i = aa_high_i;
    double bb_high_i = (bb_high_i_shift + (-3.00000000000000000e+0 * aa_high_i_shift * r2_high));
    this->bb_high_i = bb_high_i;
    double cc_high_i = (cc_high_i_shift + (-2.00000000000000000e+0 * bb_high_i_shift * r2_high) + (3.00000000000000000e+0 * aa_high_i_shift * r2_high * r2_high));
    this->cc_high_i = cc_high_i;
    double dd_high_i = (dd_high_i_shift + (bb_high_i_shift * r2_high * r2_high) + (-((aa_high_i_shift * r2_high * r2_high * r2_high))) + (-((cc_high_i_shift * r2_high))));
    this->dd_high_i = dd_high_i;
    double aa_high_j = aa_high_j_shift;
    this->aa_high_j = aa_high_j;
    double bb_high_j = (bb_high_j_shift + (-3.00000000000000000e+0 * aa_high_j_shift * r2_high));
    this->bb_high_j = bb_high_j;
    double cc_high_j = (cc_high_j_shift + (-2.00000000000000000e+0 * bb_high_j_shift * r2_high) + (3.00000000000000000e+0 * aa_high_j_shift * r2_high * r2_high));
    this->cc_high_j = cc_high_j;
    double dd_high_j = (dd_high_j_shift + (bb_high_j_shift * r2_high * r2_high) + (-((aa_high_j_shift * r2_high * r2_high * r2_high))) + (-((cc_high_j_shift * r2_high))));
    this->dd_high_j = dd_high_j;
  }

  void encode(ql::list& ll) const {
    ll
        << INTERN_(kw, lj_radius_i) << core::clasp_make_double_float(this->lj_radius_i)
        << INTERN_(kw, lj_radius_j) << core::clasp_make_double_float(this->lj_radius_j)
        << INTERN_(kw, sigmaij) << core::clasp_make_double_float(this->sigmaij)
        << INTERN_(kw, r0_low) << core::clasp_make_double_float(this->r0_low)
        << INTERN_(kw, r1_low) << core::clasp_make_double_float(this->r1_low)
        << INTERN_(kw, r2_high) << core::clasp_make_double_float(this->r2_high)
        << INTERN_(kw, r3_high) << core::clasp_make_double_float(this->r3_high)
        << INTERN_(kw, dr_low) << core::clasp_make_double_float(this->dr_low)
        << INTERN_(kw, dr_high) << core::clasp_make_double_float(this->dr_high)
        << INTERN_(kw, inv2pi32) << core::clasp_make_double_float(this->inv2pi32)
        << INTERN_(kw, lambda2_i) << core::clasp_make_double_float(this->lambda2_i)
        << INTERN_(kw, lambda2_j) << core::clasp_make_double_float(this->lambda2_j)
        << INTERN_(kw, pref_i_j) << core::clasp_make_double_float(this->pref_i_j)
        << INTERN_(kw, pref_j_i) << core::clasp_make_double_float(this->pref_j_i)
        << INTERN_(kw, f_contact_i) << core::clasp_make_double_float(this->f_contact_i)
        << INTERN_(kw, f_contact_j) << core::clasp_make_double_float(this->f_contact_j)
        << INTERN_(kw, f_r1_i) << core::clasp_make_double_float(this->f_r1_i)
        << INTERN_(kw, f_r1_j) << core::clasp_make_double_float(this->f_r1_j)
        << INTERN_(kw, df_r1_i) << core::clasp_make_double_float(this->df_r1_i)
        << INTERN_(kw, df_r1_j) << core::clasp_make_double_float(this->df_r1_j)
        << INTERN_(kw, f_r2_i) << core::clasp_make_double_float(this->f_r2_i)
        << INTERN_(kw, f_r2_j) << core::clasp_make_double_float(this->f_r2_j)
        << INTERN_(kw, df_r2_i) << core::clasp_make_double_float(this->df_r2_i)
        << INTERN_(kw, df_r2_j) << core::clasp_make_double_float(this->df_r2_j)
        << INTERN_(kw, aa_low_i) << core::clasp_make_double_float(this->aa_low_i)
        << INTERN_(kw, bb_low_i) << core::clasp_make_double_float(this->bb_low_i)
        << INTERN_(kw, cc_low_i) << core::clasp_make_double_float(this->cc_low_i)
        << INTERN_(kw, dd_low_i) << core::clasp_make_double_float(this->dd_low_i)
        << INTERN_(kw, aa_low_j) << core::clasp_make_double_float(this->aa_low_j)
        << INTERN_(kw, bb_low_j) << core::clasp_make_double_float(this->bb_low_j)
        << INTERN_(kw, cc_low_j) << core::clasp_make_double_float(this->cc_low_j)
        << INTERN_(kw, dd_low_j) << core::clasp_make_double_float(this->dd_low_j)
        << INTERN_(kw, aa_high_i) << core::clasp_make_double_float(this->aa_high_i)
        << INTERN_(kw, bb_high_i) << core::clasp_make_double_float(this->bb_high_i)
        << INTERN_(kw, cc_high_i) << core::clasp_make_double_float(this->cc_high_i)
        << INTERN_(kw, dd_high_i) << core::clasp_make_double_float(this->dd_high_i)
        << INTERN_(kw, aa_high_j) << core::clasp_make_double_float(this->aa_high_j)
        << INTERN_(kw, bb_high_j) << core::clasp_make_double_float(this->bb_high_j)
        << INTERN_(kw, cc_high_j) << core::clasp_make_double_float(this->cc_high_j)
        << INTERN_(kw, dd_high_j) << core::clasp_make_double_float(this->dd_high_j)
        << INTERN_(kw, i1) << core::make_fixnum(this->i3x1)
        << INTERN_(kw, i2) << core::make_fixnum(this->i3x2)
        ;
  }
};
