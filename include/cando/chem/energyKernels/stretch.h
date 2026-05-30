#pragma once
#include <math.h>

struct TermStretch {
  double kb;
  double r0;
  int i3x1;
  int i3x2;
  TermStretch() {}
  TermStretch(double kb, double r0, int i3x1, int i3x2) {
    this->kb = kb;
    this->r0 = r0;
    this->i3x1 = i3x1;
    this->i3x2 = i3x2;
  }

  void encode(ql::list& ll) const {
    ll
        << INTERN_(kw, kb) << core::clasp_make_double_float(this->kb)
        << INTERN_(kw, r0) << core::clasp_make_double_float(this->r0)
        << INTERN_(kw, i1) << core::make_fixnum(this->i3x1)
        << INTERN_(kw, i2) << core::make_fixnum(this->i3x2)
        ;
  }
};
