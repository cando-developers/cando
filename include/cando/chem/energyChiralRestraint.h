/*
    File: energyChiralRestraint.h
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
 *	energyChiralRestraint.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyChiralRestraint_H  //[
#define	EnergyChiralRestraint_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>


namespace chem
{


struct TermChiralRestraint
{
	REAL		K;
	REAL		CO;
	INT		I1;
	INT		I2;
	INT		I3;
	INT		I4;
};

inline	string	XmlTag_ChiralRestraint() { return "ChiralRestraint"; };



class	EnergyChiralRestraint : public EnergyTerm {
public:
	string	className()	{ return "EnergyChiralRestraint"; };
public:
		// Variables
        Atom_sp      _Atom1;
        Atom_sp      _Atom2;
        Atom_sp      _Atom3;
        Atom_sp      _Atom4;
		// Threshold
	TermChiralRestraint	term;

#if TURN_ENERGY_FUNCTION_DEBUG_ON
	bool		_calcForce;
	bool		_calcDiagonalHessian;
	bool		_calcOffDiagonalHessian;
#include <cando/chem/energy_functions/_ChiralRestraint_debugEvalDeclares.cc>
#endif
        Atom_sp	getAtom1() { return this->_Atom1; };
        Atom_sp	getAtom2() { return this->_Atom2; };
        Atom_sp	getAtom3() { return this->_Atom3; };
        Atom_sp	getAtom4() { return this->_Atom4; };
public:
public:
  core::List_sp encode() const;
  adapt::QDomNode_sp	asXml();
  void	parseFromXmlUsingAtomTable(adapt::QDomNode_sp xml, AtomTable_sp atomTable );
public:
  string	description();
};
};
namespace translate {

template <>
struct	to_object<chem::EnergyChiralRestraint >
{
  typedef	core::Cons_sp ExpectedType;
  typedef	core::Cons_sp DeclareType;
  static core::T_sp convert(const chem::EnergyChiralRestraint& tt)
  {
    return tt.encode();
  }
};

template <>
struct	from_object<chem::EnergyChiralRestraint>
{
  typedef	chem::EnergyChiralRestraint	ExpectedType;
  typedef	ExpectedType 		DeclareType;
	DeclareType _v;
	from_object(core::T_sp o)
	{
          SIMPLE_ERROR("Implement me");
        }
};
};

namespace chem {


num_real	_evaluateEnergyOnly_ChiralRestraint(
		num_real x1, num_real y1, num_real z1,
		num_real x2, num_real y2, num_real z2,
		num_real x3, num_real y3, num_real z3,
		num_real x4, num_real y4, num_real z4,
		num_real K, num_real CO );




FORWARD(EnergyChiralRestraint);
class EnergyChiralRestraint_O : public EnergyComponent_O
{
    LISP_CLASS(chem,ChemPkg,EnergyChiralRestraint_O,"EnergyChiralRestraint",EnergyComponent_O);
public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
public: // virtual functions inherited from Object
    void	initialize();
//	string	__repr__() const;

public:
    typedef EnergyChiralRestraint	TermType;
public: // instance variables
    gctools::Vec0<TermType>	_Terms;
    gctools::Vec0<TermType>	_BeyondThresholdTerms;

public:	// Creation class functions
public:	
    typedef gctools::Vec0<TermType>::iterator iterator;
    iterator begin() { return this->_Terms.begin(); };
    iterator end() { return this->_Terms.end(); };
//added by G 7.19.2011
public:
    virtual size_t numberOfTerms() { return this->_Terms.size();};
    

public:
    void addTerm(const TermType& term);
    virtual void dumpTerms(core::HashTable_sp atomTypes);
    virtual bool is_restraint() const { return true; };

    virtual void setupHessianPreconditioner(NVector_sp nvPosition,
					    AbstractLargeSquareMatrix_sp m,
                                            core::T_sp activeAtomMask );
  virtual num_real evaluateAllComponent( ScoringFunction_sp scorer,
                                         NVector_sp 	pos,
                                         core::T_sp componentEnergy,
                                         bool 		calcForce,
                                         gc::Nilable<NVector_sp> 	force,
                                         bool		calcDiagonalHessian,
                                         bool		calcOffDiagonalHessian,
                                         gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                         gc::Nilable<NVector_sp>	hdvec,
                                         gc::Nilable<NVector_sp> dvec,
                                         core::T_sp activeAtomMask,
                                         core::T_sp debugInteractions );

    virtual	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm(
	NVector_sp pos );

    // virtual	int	checkForBeyondThresholdInteractions( stringstream& info, NVector_sp pos );

    virtual string	beyondThresholdInteractionsAsString();

  EnergyChiralRestraint_sp copyFilter(core::T_sp keepInteraction);

public:
    EnergyChiralRestraint_O( const EnergyChiralRestraint_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(EnergyChiralRestraint_O);
};

};

TRANSLATE(chem::EnergyChiralRestraint_O);
#endif //]
