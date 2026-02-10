#pragma once
#include <math.h>

struct AmberNonbond14Term {
  double a;
  double c;
  double kqq;
  int i3x1;
  int i3x2;
  AmberNonbond14Term() {}
  AmberNonbond14Term(double a, double c, double kqq, int i3x1, int i3x2) {
    this->a = a;
    this->c = c;
    this->kqq = kqq;
    this->i3x1 = i3x1;
    this->i3x2 = i3x2;
  }

  void encode(ql::list& ll) const {
    ll
        << INTERN_(kw, a) << core::clasp_make_double_float(this->a)
        << INTERN_(kw, c) << core::clasp_make_double_float(this->c)
        << INTERN_(kw, kqq) << core::clasp_make_double_float(this->kqq)
        << INTERN_(kw, i1) << core::make_fixnum(this->i3x1)
        << INTERN_(kw, i2) << core::make_fixnum(this->i3x2)
        ;
  }
};
