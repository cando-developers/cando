//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	energyAnchorRestraint.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyAnchorRestraint_H  //[
#define	EnergyAnchorRestraint_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>



namespace chem
{


    typedef struct {
	double		ka;
	int		I1;
	double		xa;
	double		ya;
	double		za;
    } TermAnchorRestraint;

/*! Store a AnchorRestraint energy term
 */
class	EnergyAnchorRestraint : public EnergyTerm {
public:
	string	className()	{ return "EnergyAnchorRestraint"; };
public:
	TermAnchorRestraint	term;
		// Variables
        Atom_sp      _Atom1;
		// Threshold
#if TURN_ENERGY_FUNCTION_DEBUG_ON
	bool		_calcForce;
	bool		_calcDiagonalHessian;
	bool		_calcOffDiagonalHessian;
#include <cando/chem/_AnchorRestraint_debugEvalDeclares.cc>
#endif
	Atom_sp	getAtom() { return this->_Atom1; };
	double		getXa() { return this->term.xa; };
	double		getYa() { return this->term.ya; };
	double		getZa() { return this->term.za; };
public:
//	void	archive(core::ArchiveP node);
public:
	geom::QDomNode_sp	asXml(core::Lisp_sp);
	void	parseFromXmlUsingAtomTable(geom::QDomNode_sp xml, AtomTable_sp atomTable );

    EnergyAnchorRestraint();
	virtual ~EnergyAnchorRestraint();

};



double	_evaluateEnergyOnly_AnchorRestraint(
		double x1, double y1, double z1,
		double xa, double ya, double za,
		double kb );

class EnergyAnchorRestraint_O : public EnergyComponent_O
{
    LISP_BASE1(EnergyComponent_O);
    LISP_CLASS(chem,ChemPkg,EnergyAnchorRestraint_O,"EnergyAnchorRestraint");
public:
public: // virtual functions inherited from Object
    void	initialize();
//    void	archiveBase(core::ArchiveP node);
public:
    typedef EnergyAnchorRestraint	TermType;
public: // instance variables
    gctools::Vec0<TermType>		_Terms;
    gctools::Vec0<TermType>		_BeyondThresholdTerms;

public:	// Creation class functions
public:	
    typedef gctools::Vec0<TermType>::iterator iterator;
    iterator begin() { return this->_Terms.begin(); };
    iterator end() { return this->_Terms.end(); };
//added by G 7.19.2011
public:
    virtual int numberOfTerms() { return this->_Terms.size();};


public:
    void addTerm(const TermType& term);
    virtual void dumpTerms();

    virtual void setupHessianPreconditioner(NVector_sp nvPosition,
					    AbstractLargeSquareMatrix_sp m );
    virtual	void evaluateAll ( NVector_sp 	pos,
				   bool 		calcForce,
				   NVector_sp 	force,
				   bool		calcDiagonalHessian,
				   bool		calcOffDiagonalHessian,
				   AbstractLargeSquareMatrix_sp	hessian,
				   NVector_sp	hdvec,
				   NVector_sp dvec);

    virtual	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm(
	NVector_sp pos );

    virtual	int	checkForBeyondThresholdInteractions( stringstream& info, NVector_sp pos );

    virtual string	beyondThresholdInteractionsAsString();


public:
    EnergyAnchorRestraint_O( const EnergyAnchorRestraint_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(EnergyAnchorRestraint_O);
};




};

TRANSLATE(chem::EnergyAnchorRestraint_O);
#endif //]


