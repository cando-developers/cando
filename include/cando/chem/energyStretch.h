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
//	void	archive(core::ArchiveP node);
public:
	adapt::QDomNode_sp	asXml();
	void	parseFromXmlUsingAtomTable(adapt::QDomNode_sp xml, AtomTable_sp atomTable );


		//! Initialize self from a FFStretch parameter and two EnergyAtom(s)
	void defineFrom( FFStretch_sp stretch, EnergyAtom *ea1, EnergyAtom *ea2, double scale);
		//! Handle situation when FFStretch term could not be found for these two atoms
	void defineMissing( EnergyAtom *ea1, EnergyAtom *ea2);
//        void defineFrom( EnergyAtom *ea1, EnergyAtom *ea2);
};



double	_evaluateEnergyOnly_Stretch (
		double x1,
		double y1,
		double z1,
		double x2,
		double y2,
		double z2,
		double r0,
		double kb );






class EnergyStretch_O : public EnergyComponent_O
{
    LISP_CLASS(chem,ChemPkg,EnergyStretch_O,"EnergyStretch",EnergyComponent_O);
public:
public: // virtual functions inherited from Object
    void	initialize();
//    void	archiveBase(core::ArchiveP node);
//	string	__repr__() const;

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
    virtual void dumpTerms();

    CL_DEFMETHOD core::T_mv safe_amber_energy_stretch_term(size_t index) {
      if (index >= this->numberOfTerms() ) {
        SIMPLE_ERROR(BF("Illegal term index %zu must be less than %zu") % index % this->_Terms.size() );
      }
      return Values(core::DoubleFloat_O::create(this->_Terms[index].term.kb),
                    core::DoubleFloat_O::create(this->_Terms[index].term.r0),
                    core::make_fixnum(this->_Terms[index].term.I1),
                    core::make_fixnum(this->_Terms[index].term.I2));
    }
    
    core::T_sp stretchTermBetweenAtoms(Atom_sp a1, Atom_sp a2);
    
    virtual void setupHessianPreconditioner(NVector_sp nvPosition,
					    AbstractLargeSquareMatrix_sp m );
    virtual double evaluateAll( NVector_sp 	pos,
                                bool 		calcForce,
                                gc::Nilable<NVector_sp> 	force,
                                bool		calcDiagonalHessian,
                                bool		calcOffDiagonalHessian,
                                gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                gc::Nilable<NVector_sp>	hdvec,
                                gc::Nilable<NVector_sp> dvec);

    virtual	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm(
	NVector_sp pos );

    virtual	int	checkForBeyondThresholdInteractions( stringstream& info, NVector_sp pos );

    virtual string	beyondThresholdInteractionsAsString();


public:
    EnergyStretch_O( const EnergyStretch_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(EnergyStretch_O);
};







};

TRANSLATE(chem::EnergyStretch_O);
#endif //]
