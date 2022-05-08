/*
    File: energyNonbond.h
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
 *	energyNonbond.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergySketchNonbond_H  //[
#define	EnergySketchNonbond_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>
//#include "geom/render.fwd.h"// energyNonbond.h wants DisplayList needs render.fwd.h

namespace chem {



struct EnergySketchNonbond
{
  size_t        _FreezeFlags;
  double       _Constant;
  int		I1; //!< i*3 index into coordinate vector, must match Mathematica code!
  int		I2; //!< i*3 index into coordinate vector, must match Mathematica code!
  EnergySketchNonbond(size_t flags, double constant, int i1, int i2) : _FreezeFlags(flags), _Constant(constant), I1(i1), I2(i2) {};
  EnergySketchNonbond() {};
    core::List_sp encode() const;
  void decode(core::List_sp alist);

};


};

namespace translate {

template <>
struct	to_object<chem::EnergySketchNonbond >
{
  typedef	core::Cons_sp ExpectedType;
  typedef	core::Cons_sp DeclareType;
  static core::T_sp convert(const chem::EnergySketchNonbond& sketchnb)
  {
    return sketchnb.encode();
  }
};

template <>
struct	from_object<chem::EnergySketchNonbond>
{
  typedef	chem::EnergySketchNonbond	ExpectedType;
  typedef	ExpectedType 		DeclareType;
	DeclareType _v;
  from_object(core::T_sp o)
	{
          SIMPLE_ERROR(("Implement me"));
        }
};
};

namespace chem {

double	_evaluateEnergyOnly_Nonbond(
		double x1, double y1, double z1,
		double x2, double y2, double z2,
		double dC );

FORWARD(EnergySketchNonbond);
class EnergySketchNonbond_O : public EnergyComponent_O
{
  LISP_CLASS(chem,ChemPkg,EnergySketchNonbond_O,"EnergySketchNonbond",EnergyComponent_O);
  
 public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
  double _ScaleSketchNonbond;
  double _LongDistanceCutoff;
  bool   _IgnoreHydrogensAndLps;
  double _Energy;
  size_t _FreezeFlags;
  gctools::Vec0<EnergySketchNonbond> _Terms;
 public: // virtual functions inherited from Object
  void	initialize();
 public:
    virtual	void	zeroEnergy();
  virtual	double	getEnergy();

  void setScaleSketchNonbond(double d);
  double	getScaleSketchNonbond();

  void setFreezeFlags(size_t freezeFlags);
  
  void addSketchNonbondTerm(size_t coordinate1IndexTimes3, size_t coordinate2IndexTimes3, size_t freezeFlags, double constant);
  
  virtual void setupHessianPreconditioner(NVector_sp nvPosition,
                                          AbstractLargeSquareMatrix_sp m );
    
  virtual double evaluateAllComponent( ScoringFunction_sp scorer,
                              NVector_sp 	pos,
                              bool 		calcForce,
                              gc::Nilable<NVector_sp> 	force,
                              bool		calcDiagonalHessian,
                              bool		calcOffDiagonalHessian,
                              gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                              gc::Nilable<NVector_sp>	hdvec,
                              gc::Nilable<NVector_sp> dvec);

  void evaluateTerms( NVector_sp 	pos,
                      bool 		calcForce,
                      gc::Nilable<NVector_sp> 	force,
                      bool		calcDiagonalHessian,
                      bool		calcOffDiagonalHessian,
                      gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                      gc::Nilable<NVector_sp>	hdvec,
                      gc::Nilable<NVector_sp> dvec);

  void setLongDistanceCutoff(float cutoff);
  void walkSketchNonbondTerms(core::T_sp callback);
  void modifySketchNonbondTermConstant(size_t index, float constant);

 public:
  EnergySketchNonbond_O() : _LongDistanceCutoff(80.0), _ScaleSketchNonbond(1.0), _IgnoreHydrogensAndLps(false), _FreezeFlags() {};
  void reset();
};

};

TRANSLATE(chem::EnergySketchNonbond_O);
#endif //]
