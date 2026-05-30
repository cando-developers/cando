#pragma once
#include <math.h>

struct linear_angle_term {
  double kt;
  int i3x1;
  int i3x2;
  int i3x3;
  linear_angle_term() {}
  linear_angle_term(double kt, int i3x1, int i3x2, int i3x3) {
    this->kt = kt;
    this->i3x1 = i3x1;
    this->i3x2 = i3x2;
    this->i3x3 = i3x3;
  }

  void encode(ql::list& ll) const {
    ll
        << INTERN_(kw, kt) << core::clasp_make_double_float(this->kt)
        << INTERN_(kw, i1) << core::make_fixnum(this->i3x1)
        << INTERN_(kw, i2) << core::make_fixnum(this->i3x2)
        << INTERN_(kw, i3) << core::make_fixnum(this->i3x3)
        ;
  }
};
