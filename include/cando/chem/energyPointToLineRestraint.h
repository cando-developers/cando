/*
    File: energyPointToLineRestraint.h
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
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	energyPointToLineRestraint.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyPointToLineRestraint_H  //[
#define	EnergyPointToLineRestraint_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>
#include <cando/chem/energySketchStretch.h>



namespace chem
{


/*! Store a PointToLineRestraint energy term
 */

double	_evaluateEnergyOnly_PointToLineRestraint(
		double x1, double y1, double z1,
		double xa, double ya, double za,
		double xb, double yb, double zb,
		double ka, double ra );
FORWARD(EnergyPointToLineRestraint);
class EnergyPointToLineRestraint_O : public EnergyComponent_O
{
    LISP_CLASS(chem,ChemPkg,EnergyPointToLineRestraint_O,"EnergyPointToLineRestraint",EnergyComponent_O);
public: // virtual functions inherited from Object
    void	initialize();
public: // instance variables
  double                  _Bond_div_2;  // cutoff - if _ForceConstant is negative then less than this force is off
  EnergySketchStretch_sp  _Stretch;
  double                  _ForceConstant; // positive repels, negative attracts

public:
  static EnergyPointToLineRestraint_sp create(EnergySketchStretch_sp stretch);
public:
  virtual double evaluateAll( ScoringFunction_sp scorer,
                              NVector_sp 	pos,
                              bool 		calcForce,
                              gc::Nilable<NVector_sp> 	force,
                              bool		calcDiagonalHessian,
                              bool		calcOffDiagonalHessian,
                              gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                              gc::Nilable<NVector_sp>	hdvec,
                              gc::Nilable<NVector_sp> dvec);

public:
  EnergyPointToLineRestraint_O(EnergySketchStretch_sp stretch) : _Stretch(stretch), _ForceConstant(0.5) {};
  void reset();



};


};

#endif //]
