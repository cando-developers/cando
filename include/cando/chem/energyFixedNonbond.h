/*
    File: energyFixedNonbond.h
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
 *	energyFixedNonbond.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyFixedNonbond_H  //[
#define	EnergyFixedNonbond_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>

//#include "geom/render.fwd.h"// energyFixedNonbond.h wants DisplayList needs render.fwd.h

namespace chem
{

SMART(ForceField);
SMART(FFNonbondCrossTermTable);
SMART(FFNonbondDb);
SMART(EnergyFixedNonbondRestraint);
//SMART(DisplayList);


/*! A FixedNonbond term, this is subclassed in AmberFunction to deal with Amber specific stuff.
 * I should combine energyFunction.cc and amberFunction.cc together.
 */
class FixedNonbondRestraint : public EnergyTerm {
public:
	string	className()	{ return "FixedNonbondRestraint"; };
public:
                // Parameters
                // Variables
        Atom_sp		_FixedAtom;
	double		_FixedCharge;
	uint		_FixedType;
	Vector3		_FixedPosition;
//	int		I1; //!< i*3 index into coordinate vector, must match Mathematica code!
#if TURN_ENERGY_FUNCTION_DEBUG_ON
	bool		_calcForce;
	bool		_calcDiagonalHessian;
	bool		_calcOffDiagonalHessian;
#include <cando/chem/energy_functions/_FixedNonbond_debugEvalDeclares.cc>
#endif
        Atom_sp	getFixedAtom() { return this->_FixedAtom; };
	double	getDistance();
	bool	defineFrom(Atom_sp a);

public:
//	void	archive(core::ArchiveP node);
public:
	adapt::QDomNode_sp	asXml();
	void	parseFromXmlRelativeToContainer( adapt::QDomNode_sp xml, Matter_sp parent );
	void	parseFromXmlUsingAtomTable(adapt::QDomNode_sp xml, AtomTable_sp atomTable );

};

class BeyondThresholdFixedNonbondRestraint : public BeyondThresholdEnergyTerm 
{
 public:
  Atom_sp	_MobileAtom;
  Atom_sp	_FixedAtom;
  double	_NonbondDistance;
};






double	_evaluateEnergyOnly_FixedNonbond(
		double x1, double y1, double z1,
		double x2, double y2, double z2,
		double dA, double dC, double dQ1Q2 );


class EnergyFixedNonbondRestraint_O : public EnergyComponent_O
{
    LISP_CLASS(chem,ChemPkg,EnergyFixedNonbondRestraint_O,"EnergyFixedNonbondRestraint",EnergyComponent_O);
public:
public: // virtual functions inherited from Object
    void	initialize();
//    void	archiveBase(core::ArchiveP node);

public:
    typedef FixedNonbondRestraint TermType;
    typedef BeyondThresholdFixedNonbondRestraint BeyondThresholdTermType;
public: // instance variables
    AtomTable_sp			_MobileAtomTable;
    FFNonbondCrossTermTable_sp	_NonbondCrossTermTable;
    double				_DielectricConstant;
    double				_ScaleVdw;
    double				_ScaleElectrostatic;
    double				_EnergyVdw;
    double				_EnergyElectrostatic;
    gctools::Vec0<TermType>		_Terms;
    gctools::Vec0<BeyondThresholdTermType>	_BeyondThresholdTerms;

public:	// Creation class functions
public:	
    typedef gctools::Vec0<TermType>::iterator iterator;
    iterator begin() { return this->_Terms.begin(); };
    iterator end() { return this->_Terms.end(); };
//added by G 7.19.2011
public:
    virtual size_t numberOfTerms() { return this->_Terms.size();};
    


public:
    /*! Add a fixed atom to this FixedNonbond component
     */
    void	addFixedAtom(FFNonbondDb_sp nonbondDb, Atom_sp a);
    void	setupForEvaluation(AtomTable_sp m, FFNonbondCrossTermTable_sp crossTerms);

public:
    void	setDielectricConstant(double d) { this->_DielectricConstant = d; };
    double	getDielectricConstant() { return this->_DielectricConstant; };

    void	setVdwScale(double d) { this->_ScaleVdw = d; };
    double	getVdwScale()	{return this->_ScaleVdw; };
    void	setElectrostaticScale(double d) { this->_ScaleElectrostatic = d; };
    double	getElectrostaticScale()	{return this->_ScaleElectrostatic; };

    double	getVdwEnergy() { return this->_EnergyVdw; };
    double	getElectrostaticEnergy() { return this->_EnergyElectrostatic; };


public:

    void addTerm(const TermType& term);
    virtual void dumpTerms();

    virtual	void	zeroEnergy();
//virtual void setupHessianPreconditioner(NVector_sp nvPosition,
//				AbstractLargeSquareMatrix_sp m );
    virtual double evaluateAll( NVector_sp 	pos,
                                bool 		calcForce,
                                gc::Nilable<NVector_sp> 	force,
                                bool		calcDiagonalHessian,
                                bool		calcOffDiagonalHessian,
                                gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                gc::Nilable<NVector_sp>	hdvec,
                                gc::Nilable<NVector_sp> dvec) override;

    virtual	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm(
	NVector_sp pos );

    virtual	int	checkForBeyondThresholdInteractions( stringstream& info, NVector_sp pos );

    virtual string	beyondThresholdInteractionsAsString();

//    int countBadVdwOverlaps(double scaleSumOfVdwRadii, NVector_sp pos, geom::DisplayList_sp displayIn );

    virtual	double	getEnergy();


public:

    EnergyFixedNonbondRestraint_O( const EnergyFixedNonbondRestraint_O& ss ); //!< Copy constructor


    DEFAULT_CTOR_DTOR(EnergyFixedNonbondRestraint_O);
};

};

TRANSLATE(chem::EnergyFixedNonbondRestraint_O);
#endif //]
