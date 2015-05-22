       
       
#define	DEBUG_LEVEL_NONE


//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	energyFunction.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyFunction_H
#define	EnergyFunction_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/common.h"
#include "addon/vector3.h"
#include "atom.h"




#include "addon/quickDom.fwd.h"// energyFunction.h wants QDomNode needs quickDom.fwd.h
//#include "addon/render.fwd.h"// energyFunction.h wants DisplayList needs render.fwd.h
#include "core/iterator.fwd.h"// energyFunction.h wants Iterator needs iterator.fwd.h


#include "chemPackage.h"

namespace       chem
{

class EnergyAtom;

SMART(AbstractLargeSquareMatrix);
SMART(FFNonbondCrossTermTable);
SMART(QDomNode);
SMART(Atom);
SMART(NVector);
SMART(Matter);
SMART(ForceField);
SMART(AtomTable);
SMART(EnergyStretch);
SMART(EnergyAngle);
SMART(EnergyDihedral);
SMART(EnergyNonbond);
SMART(EnergyImproperRestraint);
SMART(EnergyChiralRestraint);
SMART(EnergyAnchorRestraint);
SMART(EnergyFixedNonbondRestraint);
//SMART(DisplayList);
SMART(Iterator);
SMART(ForceMatchReport);



#define	DefaultChiralRestraintOffset	0.2
#define	DefaultChiralRestraintWeight	100000.0
#define	DefaultAnchorRestraintWeight	10.0


//
// Remove InteractionProblem
// once we have split up the energy components
//
//
typedef	enum { linearAngle, linearDihedral, linearImproperRestraint, overlappingNonbond } InteractionProblemType;
class	InteractionProblem 
{
public:
    string _Message;
	InteractionProblemType	_Type;
	Atom_sp		_Atom1;
	Atom_sp		_Atom2;
	Atom_sp		_Atom3;
	Atom_sp		_Atom4;

	string	message() { return this->_Message;};
	InteractionProblem() {};
	virtual ~InteractionProblem() throw() {};
};


















SMART(EnergyFunction);
class EnergyFunction_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,EnergyFunction_O,"EnergyFunction");
#if INIT_TO_FACTORIES
 public:
    static EnergyFunction_sp make(Matter_sp matter, ForceField_sp forceField);
#else
    DECLARE_INIT();
#endif
public:
    static void lisp_initGlobals(core::Lisp_sp lisp);
    void initialize();
public:
//    void	archive(core::ArchiveP node);

private:
    string					_Name;
    Matter_sp				_Matter;	// Aggregate or Molecule
    /*! Stores cross terms for evaluating nonbond interactions
     */
    FFNonbondCrossTermTable_sp		_NonbondCrossTermTable;
    AtomTable_sp			_AtomTable;
    EnergyStretch_sp			_Stretch;
#if USE_ALL_ENERGY_COMPONENTS
    EnergyAngle_sp			_Angle;
    EnergyDihedral_sp			_Dihedral;
    EnergyNonbond_sp			_Nonbond;
    EnergyImproperRestraint_sp		_ImproperRestraint;
    EnergyChiralRestraint_sp		_ChiralRestraint;
    EnergyAnchorRestraint_sp		_AnchorRestraint;
    EnergyFixedNonbondRestraint_sp	_FixedNonbondRestraint;
#endif
    /*! If true then secondary amides are
     * automatically restrainted to be trans
     */
    bool					_RestrainSecondaryAmides;


    double					_ChiralRestraintWeight;
    double					_ChiralRestraintOffset;
    double					_AnchorRestraintWeight;
    double					_FixedNonbondRestraintWeight;
    double					_TotalEnergy;
    string					_Message;
    double					_DielectricConstant;
    std::set<string>				_MissingParameters;
public:
public:
    void	_eraseMissingParameters() { this->_MissingParameters.clear();};
    void	_addMissingParameter(const string& s) { this->_MissingParameters.insert(s);};
    void __createSecondaryAmideRestraints(VectorAtom& nitrogens);





    void	flagImproperRestraintsAboveThreshold(NVector_sp nvPosition);

private:
    void _applyRestraints(ForceField_sp forceField, core::Iterator_sp restraintIterator);
    void _applyDihedralRestraint(Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4, double min, double max, double weight);

public:

    AtomTable_sp atomTable() const { return this->_AtomTable;};

    string	energyTermsEnabled() ;
    void	extractCoordinatesFromAtoms(NVector_sp pos);
    void	writeCoordinatesToAtoms(NVector_sp pos);
    void	writeCoordinatesAndForceToAtoms(NVector_sp pos, NVector_sp force);
    uint	getNVectorSize();
    double	evaluateRaw( NVector_sp pos, NVector_sp force );
    double	evaluate( NVector_sp pos, NVector_sp force, bool calculateForce );
    addon::QDomNode_sp	identifyTermsBeyondThreshold();
//    uint	countBadVdwInteractions(double scaleSumOfVdwRadii, addon::DisplayList_sp displayIn);

    ForceMatchReport_sp checkIfAnalyticalForceMatchesNumericalForce( NVector_sp pos, NVector_sp force );

    Matter_sp	getMatter() { return this->_Matter;};

    void	useDefaultSettings();


    /*! Set a single options */
    void	setOption( core::Symbol_sp option, core::T_sp val);


    /*! Set the energy function options. List the options as a flat list of keyword/value pairs */
    void	setOptions( core::Cons_sp options );

    EnergyStretch_sp	getStretchComponent() { return this->_Stretch; };
#if USE_ALL_ENERGY_COMPONENTS
    EnergyAngle_sp	getAngleComponent() { return this->_Angle; };
    EnergyDihedral_sp	getDihedralComponent() { return this->_Dihedral; };
    EnergyNonbond_sp	getNonbondComponent() { return this->_Nonbond; };
    EnergyChiralRestraint_sp	getChiralRestraintComponent() { return this->_ChiralRestraint; };
    EnergyAnchorRestraint_sp	getAnchorRestraintComponent() { return this->_AnchorRestraint; };
    EnergyImproperRestraint_sp	getImproperRestraintComponent() { return this->_ImproperRestraint; };
    EnergyFixedNonbondRestraint_sp	getFixedNonbondRestraintComponent() { return this->_FixedNonbondRestraint; };
#endif

    double getDihedralComponentEnergy(); // { return this->_Dihedral->getEnergy(); };
    double getImproperRestraintComponentEnergy(); // { return this->_ImproperRestraint->getEnergy(); };

    bool	hasMissingParameters();
    string	getMissingParameters();

    double	getTotalEnergy() { return this->_TotalEnergy; };
    void	setupHessianPreconditioner( NVector_sp pos, AbstractLargeSquareMatrix_sp hessian);

    /*! Enable debugging on all energy components
     */
    void	enableDebug();
    /*! Disable debugging on all energy components
     */
    void	disableDebug();

    void	dumpTerms();
    string	getMessage() { return this->_Message;};

    /*! Print the energy components as a single, multi-line string
     */
    string	energyComponentsAsString();


    int	compareAnalyticalAndNumericalForceAndHessianTermByTerm(NVector_sp pos);
    int	compareAnalyticalAndNumericalForceAndHessianTermByTermAtCurrentPosition();



    void	setName(const string& nm) { this->_Name = nm; };
    string	getName() { return this->_Name; };
    void		writeForceToAtoms(NVector_sp f);
    EnergyAtom*     getEnergyAtomPointer(Atom_sp a);

    double		calculateEnergy();
    double		calculateEnergyAndForce();

    void	defineForMatter(Matter_sp agg, ForceField_sp forceField );


    /*! Add the restraints to the energy function.
     * This allows restraints to be applied to the system
     * without having to add them to the molecule/aggregate.
     */
    void	addTermsForListOfRestraints( ForceField_sp forceField,  core::Cons_sp restraintList );


    double	calculateNumericalDerivative(NVector_sp pos, double delta, uint i );
    double	calculateNumericalSecondDerivative(NVector_sp pos, double delta, uint i, uint j );
    double	evaluateAll( 	NVector_sp pos,
				bool calcForce,
				NVector_sp force,
       				bool calcDiagonalHessian,
				bool calcOffDiagonalHessian,
				AbstractLargeSquareMatrix_sp	hessian,
				NVector_sp hdvec, NVector_sp dvec	);
    double	evaluateEnergy( NVector_sp pos );
    double	evaluateEnergyForce( NVector_sp pos, bool calcForce, NVector_sp force );
    double	evaluateEnergyForceFullHessian(
	NVector_sp pos,
	bool calcForce, NVector_sp force,
	bool calcDiagonalHessian,
	bool calcOffDiagonalHessian,
	AbstractLargeSquareMatrix_sp hessian );
    double	evaluateEnergyForceFullHessianForDebugging();

    std::set<string>::iterator	begin_MissingParameters() {
	return this->_MissingParameters.begin(); };
    std::set<string>::iterator	end_MissingParameters() {
	return this->_MissingParameters.end(); };




    string	summarizeBeyondThresholdInteractionsAsString();
    string	summarizeEnergyAsString();


//		addon::QDomNode_sp	rawAccumulateTermsBeyondThresholdAsXml(uint& count);
//		addon::QDomNode_sp	accumulateTermsBeyondThresholdAsXml();
    uint		countTermsBeyondThreshold();

    void	evaluateNumericalForce(NVector_sp pos, NVector_sp numForce, double delta );
    void	evaluateNumericalHessian(NVector_sp pos, AbstractLargeSquareMatrix_sp numHessian, bool calcOffDiagonalElements, double delta);

    uint	checkForBeyondThresholdInteractions();

    string	debugLogAsString();

    void	dealWithProblem(const InteractionProblem& problem);
    DEFAULT_CTOR_DTOR(EnergyFunction_O);
};



#define	FINITE_DIFFERENCE_TOLERANCE 0.10
int	_areValuesClose(double numVal, double analVal, const char* funcName, const char* termName, int index );

};

TRANSLATE(chem::EnergyFunction_O);
#endif
