/*
    File: boundingBox.h
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
#ifndef BoundingCuboid_H //[
#define BoundingCuboid_H

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <cando/geom/geomPackage.fwd.h>
#include <cando/geom/matrix.h>
#include <cando/geom/vector3.h>

namespace geom {

FORWARD(BoundingCuboid);
class BoundingCuboid_O : public core::CxxObject_O {
  LISP_CLASS(geom, GeomPkg, BoundingCuboid_O, "BoundingCuboid", core::CxxObject_O);

public:
  void initialize();

private:
  // instance variables
  bool _Defined;
  Vector3 _MinCorner;
  Vector3 _MaxCorner;

public:
  bool fieldsp() const { return true; }
  void fields(core::Record_sp node);

public:
  bool isDefined() { return this->_Defined; };
  Vector3 &minCorner() { return this->_MinCorner; };
  Vector3 &maxCorner() { return this->_MaxCorner; };

  double getMinX() const;
  double getMaxX() const;
  double getMinY() const;
  double getMaxY() const;
  double getMinZ() const;
  double getMaxZ() const;

  Vector3 getCenter() const;
  double getExtentX() const;
  double getExtentY() const;
  double getExtentZ() const;
  Vector3 getHalfWidth() const;

  string __repr__() const;

  bool encompasses(const Vector3 &point);
  void expandToEncompassPoint(const Vector3 &point);
  void pad(double add);

  BoundingCuboid_O(const BoundingCuboid_O &ss); //!< Copy constructor

  CL_LISPIFY_NAME(make-bounding-cuboid);
  CL_DEF_CLASS_METHOD static BoundingCuboid_sp make()
  {
    return create();
  }

  DEFAULT_CTOR_DTOR(BoundingCuboid_O);
};

}; // namespace geom
#endif //]
