/*
    File: quantity.fwd.h
*/
/*
Open Source License
Copyright (c) 2016, Christian E. Schafmeister
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University
at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
/* @(#)quantity.fwd.h
 */

#ifndef _units_Quantity_fwd_H
#define _units_Quantity_fwd_H

namespace units {

inline int gcd(int a, int b) {
  if (a == 0)
    return b;
  return gcd(b % a, a);
}

struct SimpleRational {
  int _Num;
  int _Denom;

  SimpleRational() : _Num(0), _Denom(1){};
  SimpleRational(int num) : _Num(num), _Denom(1){};
  SimpleRational(int num, int denom) : _Num(num), _Denom(denom){};

  SimpleRational operator*(const SimpleRational &other) const {
    return SimpleRational(this->_Num * other._Num, this->_Denom * other._Denom).maybe_simplify();
  }

  SimpleRational operator*(int other) const {
    return SimpleRational(this->_Num * other, this->_Denom).maybe_simplify();
  }

  SimpleRational operator/(const SimpleRational &other) const {
    return SimpleRational(this->_Num * other._Denom, this->_Denom * other._Num).maybe_simplify();
  }

  SimpleRational operator+(const SimpleRational &other) const {
    if (this->_Denom == other._Denom)
      return SimpleRational(this->_Num + other._Num, this->_Denom);
    auto result = SimpleRational(this->_Num * other._Denom + other._Num * this->_Denom, this->_Denom * other._Denom);
    return result.maybe_simplify();
  }

  SimpleRational operator-(const SimpleRational &other) const {
    if (this->_Denom == other._Denom)
      return SimpleRational(this->_Num - other._Num, this->_Denom);
    auto result = SimpleRational(this->_Num * other._Denom - other._Num * this->_Denom, this->_Denom * other._Denom);
    return result.maybe_simplify();
  }
  bool operator==(const SimpleRational &other) const { return (this->_Num == other._Num && this->_Denom == other._Denom); }
  bool operator!=(const SimpleRational &other) const {
    // Assume everything is normalized
    return (this->_Num != other._Num) || (this->_Denom != other._Denom);
  }

  bool operator!=(int other) const {
    return (this->_Denom * other) == this->_Num;
  }

  double as_double() const { return (double)this->_Num / (double)this->_Denom; }

  SimpleRational maybe_simplify() const {
    if (this->_Num % this->_Denom == 0) {
      return SimpleRational(this->_Num / this->_Denom, 1);
    }
    int gcd_ = gcd(this->_Num, this->_Denom);
    if (gcd_ > 1) {
      return SimpleRational(this->_Num / gcd_, this->_Denom / gcd_);
    }
    return *this;
  }
};

}; // namespace units

namespace units {
FORWARD(Quantity);

}; // namespace units

#endif /* _units_Quantity_fwd_H */
