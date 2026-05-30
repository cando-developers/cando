#pragma once
#include <math.h>

struct anchor_term {
  double ka;
  double xa;
  double ya;
  double za;
  int i3x1;
  anchor_term() {}
  anchor_term(double ka, double xa, double ya, double za, int i3x1) {
    this->ka = ka;
    this->xa = xa;
    this->ya = ya;
    this->za = za;
    this->i3x1 = i3x1;
  }

  void encode(ql::list& ll) const {
    ll
        << INTERN_(kw, ka) << core::clasp_make_double_float(this->ka)
        << INTERN_(kw, xa) << core::clasp_make_double_float(this->xa)
        << INTERN_(kw, ya) << core::clasp_make_double_float(this->ya)
        << INTERN_(kw, za) << core::clasp_make_double_float(this->za)
        << INTERN_(kw, i1) << core::make_fixnum(this->i3x1)
        ;
  }
};
