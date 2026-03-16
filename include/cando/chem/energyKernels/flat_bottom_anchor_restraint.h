#pragma once
#include <math.h>

struct flat_bottom_anchor_restraint_term {
  double ka;
  double r0;
  double r02;
  double xa;
  double ya;
  double za;
  int i3x1;
  flat_bottom_anchor_restraint_term() {}
  flat_bottom_anchor_restraint_term(double ka, double r0, double r02, double xa, double ya, double za, int i3x1) {
    this->ka = ka;
    this->r0 = r0;
    this->r02 = r02;
    this->xa = xa;
    this->ya = ya;
    this->za = za;
    this->i3x1 = i3x1;
  }

  void encode(ql::list& ll) const {
    ll
        << INTERN_(kw, ka) << core::clasp_make_double_float(this->ka)
        << INTERN_(kw, r0) << core::clasp_make_double_float(this->r0)
        << INTERN_(kw, r02) << core::clasp_make_double_float(this->r02)
        << INTERN_(kw, xa) << core::clasp_make_double_float(this->xa)
        << INTERN_(kw, ya) << core::clasp_make_double_float(this->ya)
        << INTERN_(kw, za) << core::clasp_make_double_float(this->za)
        << INTERN_(kw, i1) << core::make_fixnum(this->i3x1)
        ;
  }
};
