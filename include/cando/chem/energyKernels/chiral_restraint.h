#pragma once
#include <math.h>

struct chiral_term {
  double k;
  double co;
  int i3x1;
  int i3x2;
  int i3x3;
  int i3x4;
  chiral_term() {}
  chiral_term(double k, double co, int i3x1, int i3x2, int i3x3, int i3x4) {
    this->k = k;
    this->co = co;
    this->i3x1 = i3x1;
    this->i3x2 = i3x2;
    this->i3x3 = i3x3;
    this->i3x4 = i3x4;
  }

  void encode(ql::list& ll) const {
    ll
        << INTERN_(kw, k) << core::clasp_make_double_float(this->k)
        << INTERN_(kw, co) << core::clasp_make_double_float(this->co)
        << INTERN_(kw, i1) << core::make_fixnum(this->i3x1)
        << INTERN_(kw, i2) << core::make_fixnum(this->i3x2)
        << INTERN_(kw, i3) << core::make_fixnum(this->i3x3)
        << INTERN_(kw, i4) << core::make_fixnum(this->i3x4)
        ;
  }
};
