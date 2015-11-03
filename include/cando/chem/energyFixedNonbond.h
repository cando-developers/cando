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
        gc::Nilable<Atom_sp>		_FixedAtom;
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
        gc::Nilable<Atom_sp>	getFixedAtom() { return this->_FixedAtom; };
	double	getDistance();
	bool	defineFrom(Atom_sp a);

public:
//	void	archive(core::ArchiveP node);
public:
	adapt::QDomNode_sp	asXml();
	void	parseFromXmlRelativeToContainer( adapt::QDomNode_sp xml, Matter_sp parent );
	void	parseFromXmlUsingAtomTable(adapt::QDomNode_sp xml, AtomTable_sp atomTable );

    FixedNonbondRestraint();
	virtual ~FixedNonbondRestraint();

};

class BeyondThresholdFixedNonbondRestraint : public BeyondThresholdEnergyTerm 
{
    public:
  gc::Nilable<Atom_sp>	_MobileAtom;
  gc::Nilable<Atom_sp>	_FixedAtom;
	double	_NonbondDistance;
};






double	_evaluateEnergyOnly_FixedNonbond(
		double x1, double y1, double z1,
		double x2, double y2, double z2,
		double dA, double dC, double dQ1Q2 );


class EnergyFixedNonbondRestraint_O : public EnergyComponent_O
{
    LISP_BASE1(EnergyComponent_O);
    LISP_CLASS(chem,ChemPkg,EnergyFixedNonbondRestraint_O,"EnergyFixedNonbondRestraint");
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
    virtual int numberOfTerms() { return this->_Terms.size();};
    


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
//    virtual void dumpTerms();

    virtual	void	zeroEnergy();
//virtual void setupHessianPreconditioner(NVector_sp nvPosition,
//				AbstractLargeSquareMatrix_sp m );
    virtual	void evaluateAll ( NVector_sp 	pos,
				   bool 		calcForce,
				   gc::Nilable<NVector_sp> 	force,
				   bool		calcDiagonalHessian,
				   bool		calcOffDiagonalHessian,
				   gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
				   NVector_sp	hdvec,
				   NVector_sp dvec);

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


