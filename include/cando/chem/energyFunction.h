/*
    File: energyFunction.h
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
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/chem/scoringFunction.h>
#include <cando/chem/atom.h>




#include <cando/adapt/quickDom.fwd.h>// energyFunction.h wants QDomNode needs quickDom.fwd.h
//#include "geom/render.fwd.h"// energyFunction.h wants DisplayList needs render.fwd.h
#include <clasp/core/iterator.fwd.h>// energyFunction.h wants Iterator needs iterator.fwd.h


#include <cando/chem/chemPackage.h>

namespace       chem
{
  SMART(FFNonbondDb);
  
  class EnergyAtom;

  SMART(FFParameter);
  SMART(AbstractLargeSquareMatrix);
  SMART(FFNonbondCrossTermTable);
  SMART(QDomNode);
  SMART(Atom);
  SMART(Matter);
  SMART(ForceField);
  SMART(CombinedForceField);
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


 /*! Throw this exception from Cando if you want the minimizer to deal with a bad interaction */
  class	InteractionCondition 
  {
  public:
    core::Symbol_sp  condition_name;
    core::List_sp    condition_arguments;
  };

  FORWARD(EnergyFunction);
};

template <>
struct gctools::GCInfo<chem::EnergyFunction_O> {
  static bool constexpr NeedsInitialization = true;
  static bool constexpr NeedsFinalization = false;
  static GCInfo_policy constexpr Policy = normal;
};

namespace chem {
  SMART(EnergyFunction);
  class EnergyFunction_O : public ScoringFunction_O
  {
    LISP_CLASS(chem,ChemPkg,EnergyFunction_O,"EnergyFunction",ScoringFunction_O);
  public:
    static EnergyFunction_sp make(Matter_sp matter, core::T_sp system, bool useExcludedAtoms=false, core::T_sp activeAtoms=_Nil<core::T_O>(), bool progress=false, bool assign_types=false);
  public:
    void initialize();
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
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
    core::List_sp			_MissingParameters;
  public:
  public:
    void	_eraseMissingParameters() { this->_MissingParameters = _Nil<core::T_O>();};
    void	_addMissingParameter(FFParameter_sp p) { this->_MissingParameters = core::Cons_O::create(p,this->_MissingParameters);};
    void __createSecondaryAmideRestraints(VectorAtom& nitrogens, core::T_sp activeAtoms );

    void	flagImproperRestraintsAboveThreshold(NVector_sp nvPosition);

  private:
    int _applyRestraints(core::T_sp forceField, core::Iterator_sp restraintIterator, core::T_sp activeAtoms );
    void _applyDihedralRestraint(Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4, double min, double max, double weight, core::T_sp activeAtoms );

  public:

    CL_LISPIFY_NAME("atomTable");
    CL_DEFMETHOD     AtomTable_sp atomTable() const { return this->_AtomTable;};

    string	energyTermsEnabled() ;
    void	loadCoordinatesIntoVector(NVector_sp pos);
    void	saveCoordinatesFromVector(NVector_sp pos);
    void	saveCoordinatesAndForcesFromVectors(NVector_sp pos, NVector_sp force);
    size_t	getNVectorSize();
    double	evaluateRaw( NVector_sp pos, NVector_sp force );
//    double	evaluate( NVector_sp pos, NVector_sp force, bool calculateForce );
    adapt::QDomNode_sp	identifyTermsBeyondThreshold();
//    uint	countBadVdwInteractions(double scaleSumOfVdwRadii, geom::DisplayList_sp displayIn);

    ForceMatchReport_sp checkIfAnalyticalForceMatchesNumericalForce( NVector_sp pos, NVector_sp force );

    CL_LISPIFY_NAME("getMatter");
    CL_DEFMETHOD     Matter_sp	getMatter() { return this->_Matter;};

    void	useDefaultSettings();


    /*! Set a single options */
    void	setOption( core::Symbol_sp option, core::T_sp val);


    /*! Set the energy function options. List the options as a flat list of keyword/value pairs */
    void	setOptions( core::List_sp options );

    CL_LISPIFY_NAME("getStretchComponent");
    CL_DEFMETHOD     EnergyStretch_sp	getStretchComponent() { return this->_Stretch; };
#if USE_ALL_ENERGY_COMPONENTS
    CL_LISPIFY_NAME("getAngleComponent");
    CL_DEFMETHOD     EnergyAngle_sp	getAngleComponent() { return this->_Angle; };
    CL_LISPIFY_NAME("getDihedralComponent");
    CL_DEFMETHOD     EnergyDihedral_sp	getDihedralComponent() { return this->_Dihedral; };
    CL_LISPIFY_NAME("getNonbondComponent");
    CL_DEFMETHOD     EnergyNonbond_sp	getNonbondComponent() { return this->_Nonbond; };
    CL_LISPIFY_NAME("getChiralRestraintComponent");
    CL_DEFMETHOD     EnergyChiralRestraint_sp	getChiralRestraintComponent() { return this->_ChiralRestraint; };
    CL_LISPIFY_NAME("getAnchorRestraintComponent");
    CL_DEFMETHOD     EnergyAnchorRestraint_sp	getAnchorRestraintComponent() { return this->_AnchorRestraint; };
    CL_LISPIFY_NAME("getImproperRestraintComponent");
    CL_DEFMETHOD     EnergyImproperRestraint_sp	getImproperRestraintComponent() { return this->_ImproperRestraint; };
    CL_LISPIFY_NAME("getFixedNonbondRestraintComponent");
    CL_DEFMETHOD     EnergyFixedNonbondRestraint_sp	getFixedNonbondRestraintComponent() { return this->_FixedNonbondRestraint; };
#endif

    double getDihedralComponentEnergy(); // { return this->_Dihedral->getEnergy(); };
    double getImproperRestraintComponentEnergy(); // { return this->_ImproperRestraint->getEnergy(); };

    CL_DEFMETHOD bool hasMissingParameters();
    CL_DEFMETHOD core::List_sp getMissingParameters();

    CL_LISPIFY_NAME("getTotalEnergy");
    CL_DEFMETHOD     double	getTotalEnergy() { return this->_TotalEnergy; };
    void	setupHessianPreconditioner( NVector_sp pos, AbstractLargeSquareMatrix_sp hessian);

    /*! Enable debugging on all energy components
     */
    void	enableDebug();
    /*! Disable debugging on all energy components
     */
    void	disableDebug();

    void	summarizeTerms();
    void	dumpTerms();
    CL_DEFMETHOD     string	getMessage() { return this->_Message;};

    /*! Print the energy components as a single, multi-line string
     */
    string	energyComponentsAsString();

    int	compareAnalyticalAndNumericalForceAndHessianTermByTerm(NVector_sp pos);
    int	compareAnalyticalAndNumericalForceAndHessianTermByTermAtCurrentPosition();


    void		writeForceToAtoms(NVector_sp f);
    EnergyAtom*     getEnergyAtomPointer(Atom_sp a);

    void assignAtomTypes(Matter_sp matter, bool show_progress);
    void defineForMatter(Matter_sp agg, core::T_sp forceFieldSystem, bool useExcludedAtoms, core::T_sp activeAtomSet=_Nil<core::T_O>(), bool show_progress=false, bool assign_types=true );
    void defineForMatterWithAtomTypes(Matter_sp matter, core::T_sp system, bool useExcludedAtoms, core::T_sp activeAtoms, bool show_progress);
    void generateStandardEnergyFunctionTables(Matter_sp mol, CombinedForceField_sp forceField, core::T_sp activeAtomSet=_Nil<core::T_O>(), bool show_progress=false );
    void generateNonbondEnergyFunctionTables(bool useExcludedAtoms, Matter_sp agg, core::T_sp forceField, core::T_sp activeAtomSet=_Nil<core::T_O>(), bool show_progress=false );
    void generateRestraintEnergyFunctionTables(Matter_sp agg, core::T_sp nonbonds, core::T_sp activeAtomSet=_Nil<core::T_O>(), bool show_progress=false );


    /*! Add the restraints to the energy function.
     * This allows restraints to be applied to the system
     * without having to add them to the molecule/aggregate.
     */
    void	addTermsForListOfRestraints( ForceField_sp forceField,  core::List_sp restraintList, core::T_sp activeAtoms );


    double	calculateNumericalDerivative(NVector_sp pos, double delta, uint i );
    double	calculateNumericalSecondDerivative(NVector_sp pos, double delta, uint i, uint j );
    double	evaluateAll( 	NVector_sp pos,
				bool calcForce,
				gc::Nilable<NVector_sp> force,
       				bool calcDiagonalHessian,
				bool calcOffDiagonalHessian,
				gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
				gc::Nilable<NVector_sp> hdvec,
                                gc::Nilable<NVector_sp> dvec	);

    string	summarizeBeyondThresholdInteractionsAsString();
    string	summarizeEnergyAsString();


//		adapt::QDomNode_sp	rawAccumulateTermsBeyondThresholdAsXml(uint& count);
//		adapt::QDomNode_sp	accumulateTermsBeyondThresholdAsXml();
    uint		countTermsBeyondThreshold();

    void	evaluateNumericalForce(NVector_sp pos, NVector_sp numForce, double delta );
    void	evaluateNumericalHessian(NVector_sp pos, AbstractLargeSquareMatrix_sp numHessian, bool calcOffDiagonalElements, double delta);

    uint	checkForBeyondThresholdInteractions();

    string	debugLogAsString();

    void	dealWithProblem(core::Symbol_sp error_symbol, core::T_sp arguments);
    
  EnergyFunction_O() :
    _Matter(_Unbound<Matter_O>())
      , _NonbondCrossTermTable(_Unbound<FFNonbondCrossTermTable_O>())
      , _AtomTable(_Unbound<AtomTable_O>())
      , _Stretch(_Unbound<EnergyStretch_O>())
#if USE_ALL_ENERGY_COMPONENTS
      , _Angle(_Unbound<EnergyAngle_O>())
      , _Dihedral(_Unbound<EnergyDihedral_O>())
      , _Nonbond(_Unbound<EnergyNonbond_O>())
      , _ImproperRestraint(_Unbound<EnergyImproperRestraint_O>())
      , _ChiralRestraint(_Unbound<EnergyChiralRestraint_O>())
      , _AnchorRestraint(_Unbound<EnergyAnchorRestraint_O>())
      , _FixedNonbondRestraint(_Unbound<EnergyFixedNonbondRestraint_O>())
#endif
//      , _MissingParameters(_Unbound<core::List_O>())
    {};
  };


  bool inAtomSet(core::T_sp atomSet, Atom_sp atom);

#define	FINITE_DIFFERENCE_TOLERANCE 0.10
  int	_areValuesClose(double numVal, double analVal, const char* funcName, const char* termName, int index );

  void energyFunction_initializeSmarts();
};

#endif
