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
};
