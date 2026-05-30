#pragma once
#include <math.h>

struct angle_term {
  double kt;
  double t0;
  int i3x1;
  int i3x2;
  int i3x3;
  angle_term() {}
  angle_term(double kt, double t0, int i3x1, int i3x2, int i3x3) {
    this->kt = kt;
    this->t0 = t0;
    this->i3x1 = i3x1;
    this->i3x2 = i3x2;
    this->i3x3 = i3x3;
  }

  void encode(ql::list& ll) const {
    ll
        << INTERN_(kw, kt) << core::clasp_make_double_float(this->kt)
        << INTERN_(kw, t0) << core::clasp_make_double_float(this->t0)
        << INTERN_(kw, i1) << core::make_fixnum(this->i3x1)
        << INTERN_(kw, i2) << core::make_fixnum(this->i3x2)
        << INTERN_(kw, i3) << core::make_fixnum(this->i3x3)
        ;
  }
};
