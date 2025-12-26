/*
    File: atomGrid.h
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
 
This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
       
       

#ifndef	AtomGrid_H //[
#define AtomGrid_H



#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/smallMap.h>
#include <clasp/core/array.h>
#include <clasp/core/array_double.h>
#include <cando/chem/nVector.h>

#include <cando/chem/chemPackage.h>

namespace chem {


SMART(AtomGrid);
class AtomGrid_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,AtomGrid_O,"AtomGrid",core::CxxObject_O);

public:
private:
  double                        _OriginX;
  double                        _OriginY;
  double                        _OriginZ;
  double                        _CellSize;
  size_t                        _Nx;
  size_t                        _Ny;
  size_t                        _Nz;
  gc::Vec0<core::T_sp>          _Cells;
  gc::Vec0<double>              _CellMaxRadius;
  NVector_sp                    _Positions;
  core::SimpleVector_double_sp  _EffRadii;
  core::ComplexVector_double_sp _HeapBounds;
  core::ComplexVector_size_t_sp _HeapCells;
  uint8_t                       _VisitedFlag;
  gc::Vec0<uint8_t>             _Visited;
public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);

  static AtomGrid_sp makeAtomGrid(NVector_sp positions,
                                  core::SimpleVector_double_sp effRadii,
                                  double cellSize,
                                  double padding );

  size_t cellIndex(size_t ix, size_t iy, size_t iz) const;
  void decodeIndex(size_t idx, size_t& cx, size_t& cy, size_t& cz) const;
  double cellMinDist(double x, double y, double z, size_t ix, size_t iy, size_t iz) const;

  double sdfGridDouble_(double x, double y, double z, size_t& bestAtom );
  core::T_mv sdfGridDouble(double x, double y, double z );
  core::T_mv sdfGridSingle(float x, float y, float z );

  void buildAtomGrid(NVector_sp positions,
                     core::SimpleVector_double_sp effRadii,
                     double cellSize,
                     double padding = 0.0);
};




};
#endif //]
