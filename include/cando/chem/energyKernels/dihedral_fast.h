#pragma once
#include <math.h>

struct dihedral_term {
  double v;
  double n;
  double sinphase;
  double cosphase;
  int i3x1;
  int i3x2;
  int i3x3;
  int i3x4;
  dihedral_term() {}
  dihedral_term(double v, double n, double sinphase, double cosphase, int i3x1, int i3x2, int i3x3, int i3x4) {
    this->v = v;
    this->n = n;
    this->sinphase = sinphase;
    this->cosphase = cosphase;
    this->i3x1 = i3x1;
    this->i3x2 = i3x2;
    this->i3x3 = i3x3;
    this->i3x4 = i3x4;
  }

  void encode(ql::list& ll) const {
    ll
        << INTERN_(kw, v) << core::clasp_make_double_float(this->v)
        << INTERN_(kw, n) << core::clasp_make_double_float(this->n)
        << INTERN_(kw, sinphase) << core::clasp_make_double_float(this->sinphase)
        << INTERN_(kw, cosphase) << core::clasp_make_double_float(this->cosphase)
        << INTERN_(kw, i1) << core::make_fixnum(this->i3x1)
        << INTERN_(kw, i2) << core::make_fixnum(this->i3x2)
        << INTERN_(kw, i3) << core::make_fixnum(this->i3x3)
        << INTERN_(kw, i4) << core::make_fixnum(this->i3x4)
        ;
  }
};
