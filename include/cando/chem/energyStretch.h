/*
    File: energyStretch.h
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
 *	energyStretch.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyStretch1_H  //[
#define	EnergyStretch1_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>

namespace       chem {

class EnergyAtom;
SMART(FFStretch);

struct TermStretch {
  REAL	kb;	//!< Stretch force constant, this must match Mathematica code!
  REAL	r0;	//!< Stretch equilibrium distance, this must match Mathematica code!
  INT	I1;	//!< i*3 index into coordinate vector for atom1, this must match Mathematica code!
  INT	I2;	//!< i*3 index into coordinate vector for atom2, this must match Mathematica code!
  TermStretch(REAL k,REAL r, INT i1, INT i2) : kb(k), r0(r),I1(i1),I2(i2) {};
  TermStretch() {};
};


/*! Store a Stretch energy term
 */
class EnergyStretch : public EnergyTerm {
public:
  string	className()	{ return "EnergyStretch"; };
public:
  TermStretch	term;
  Atom_sp		_Atom1;
  Atom_sp		_Atom2;
#if TURN_ENERGY_FUNCTION_DEBUG_ON
  bool		_calcForce;
  bool		_calcDiagonalHessian;
  bool		_calcOffDiagonalHessian;
#include <cando/chem/energy_functions/_Stretch_debugEvalDeclares.cc>
#endif


  Atom_sp	getAtom1() {return this->_Atom1;};
  Atom_sp	getAtom2() {return this->_Atom2;};
  double	getR0()	{return this->term.r0;};
  double	getR();

public:
  EnergyStretch(Atom_sp a1, Atom_sp a2, size_t i1, size_t i2, double kb, double r0) : term(kb,r0,i1,i2),_Atom1(a1),_Atom2(a2) {};
  EnergyStretch() {};
  
public:
  adapt::QDomNode_sp	asXml();
  void	parseFromXmlUsingAtomTable(adapt::QDomNode_sp xml, AtomTable_sp atomTable );


		//! Initialize self from a FFStretch parameter and two EnergyAtom(s)
  void defineFrom( FFStretch_sp stretch, EnergyAtom *ea1, EnergyAtom *ea2, double scale);
		//! Handle situation when FFStretch term could not be found for these two atoms
  void defineMissing( EnergyAtom *ea1, EnergyAtom *ea2);
//        void defineFrom( EnergyAtom *ea1, EnergyAtom *ea2);
  core::List_sp encode() const;
  void decode(core::List_sp alist);
};
};


namespace translate {

template <>
struct	to_object<chem::EnergyStretch >
{
  typedef	core::Cons_sp ExpectedType;
  typedef	core::Cons_sp DeclareType;
  static core::T_sp convert(const chem::EnergyStretch& stretch)
  {
    return stretch.encode();
  }
};

template <>
struct	from_object<chem::EnergyStretch>
{
  typedef	chem::EnergyStretch	ExpectedType;
  typedef	ExpectedType 		DeclareType;
	DeclareType _v;
	from_object(core::T_sp o)
	{
          SIMPLE_ERROR("Implement me");
        }
};
};



namespace chem {

num_real	_evaluateEnergyOnly_Stretch (
		num_real x1,
		num_real y1,
		num_real z1,
		num_real x2,
		num_real y2,
		num_real z2,
		num_real r0,
		num_real kb );

FORWARD(EnergyStretch);

class EnergyStretch_O : public EnergyComponent_O
{
    LISP_CLASS(chem,ChemPkg,EnergyStretch_O,"EnergyStretch",EnergyComponent_O);
public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
public: // virtual functions inherited from Object
    void	initialize();
public:
    typedef EnergyStretch	TermType;
public: // instance variables
    gctools::Vec0<TermType>	_Terms;
    gctools::Vec0<TermType>	_BeyondThresholdTerms;
public:	
    typedef gctools::Vec0<TermType>::iterator iterator;
    iterator begin() { return this->_Terms.begin(); };
    iterator end() { return this->_Terms.end(); };
	
public:
    virtual size_t numberOfTerms() { return this->_Terms.size();};
    void addTerm(const TermType& term);
    virtual void dumpTerms(core::HashTable_sp atomTypes);

    CL_DEFMETHOD core::T_mv safe_amber_energy_stretch_term(size_t index) {
      if (index >= this->numberOfTerms() ) {
        SIMPLE_ERROR("Illegal term index {} must be less than {}" , index , this->_Terms.size() );
      }
      return Values(core::DoubleFloat_O::create(this->_Terms[index].term.kb),
                    core::DoubleFloat_O::create(this->_Terms[index].term.r0),
                    core::make_fixnum(this->_Terms[index].term.I1),
                    core::make_fixnum(this->_Terms[index].term.I2));
    }

    /*! Extract the components of the energy term into a collection of parallel vectors.
        The vectors are returned as an alist with keyword keys.
        The keys are :kb, :r0, :i1, :i2, :atom1 :atom2 */
    virtual core::List_sp extract_vectors_as_alist() const;
    virtual void fill_from_vectors_in_alist(core::List_sp vectors);
    
    core::T_sp stretchTermBetweenAtoms(Atom_sp a1, Atom_sp a2);
    
    virtual void setupHessianPreconditioner(NVector_sp nvPosition,
					    AbstractLargeSquareMatrix_sp m );
  virtual num_real evaluateAllComponent( ScoringFunction_sp scorer,
                                         NVector_sp 	pos,
                                         core::T_sp     componentEnergy,
                                         bool 		calcForce,
                                         gc::Nilable<NVector_sp> 	force,
                                         bool		calcDiagonalHessian,
                                         bool		calcOffDiagonalHessian,
                                         gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                         gc::Nilable<NVector_sp>	hdvec,
                                         gc::Nilable<NVector_sp> dvec,
                                         core::T_sp     activeAtomMask
                                         );

    virtual	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm(
	NVector_sp pos );

  virtual	core::List_sp checkForBeyondThresholdInteractionsWithPosition(NVector_sp pos, double threshold );

    virtual string	beyondThresholdInteractionsAsString();

  void addStretchTerm(AtomTable_sp at, Atom_sp a1, Atom_sp a2, double kb, double r0);

  void walkStretchTerms(core::T_sp callback);
  void modifyStretchTermKb(size_t index, float kb);
  void modifyStretchTermR0(size_t index, float r0);

  core::List_sp lookupStretchTerms(AtomTable_sp at, Atom_sp a1, Atom_sp a2, core::HashTable_sp atomTypes );

  EnergyStretch_sp copyFilter(core::T_sp keepInteraction);

  void reset();
public:
    EnergyStretch_O( const EnergyStretch_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(EnergyStretch_O);
};

};

#endif //]
