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
#include <clasp/core/hashTableEq.h>
#include <cando/geom/vector3.h>
#include <cando/chem/scoringFunction.h>
#include <cando/chem/energyComponent.fwd.h>
#include <cando/chem/atom.h>




#include <cando/adapt/quickDom.fwd.h>// energyFunction.h wants QDomNode needs quickDom.fwd.h
//#include "geom/render.fwd.h"// energyFunction.h wants DisplayList needs render.fwd.h
#include <clasp/core/iterator.fwd.h>// energyFunction.h wants Iterator needs iterator.fwd.h


#include <cando/chem/chemPackage.h>

namespace       chem
{
  FORWARD(FFStretchDb);
  FORWARD(FFAngleDb);
  FORWARD(FFPtorDb);
  FORWARD(FFItorDb);
  FORWARD(FFNonbondDb);
  
  class EnergyAtom;

  FORWARD(FFParameter);
  FORWARD(AbstractLargeSquareMatrix);
  FORWARD(FFNonbondCrossTermTable);
  FORWARD(QDomNode);
  FORWARD(Atom);
  FORWARD(Matter);
  FORWARD(ForceField);
  FORWARD(CombinedForceField);
  FORWARD(AtomTable);
  FORWARD(EnergyStretch);
  FORWARD(EnergyAngle);
  FORWARD(EnergyDihedral);
  FORWARD(EnergyNonbond);
  FORWARD(EnergyPeriodicBoundaryConditionNonbond);
  FORWARD(EnergyDihedralRestraint);
  FORWARD(EnergyChiralRestraint);
  FORWARD(EnergyAnchorRestraint);
  FORWARD(EnergyFixedNonbondRestraint);
  FORWARD(Iterator);
  FORWARD(ForceMatchReport);
FORWARD(BoundingBox);


#define	DefaultChiralRestraintOffset	0.2
#define	DefaultChiralRestraintWeight	100000.0
#define	DefaultAnchorRestraintWeight	10.0


bool skipInteraction( core::T_sp keepInteractionFunction,
                      core::T_sp interactionClass,
                      core::T_sp atomA = nil<core::T_O>(),
                      core::T_sp atomB = nil<core::T_O>(),
                      core::T_sp atomC = nil<core::T_O>(),
                      core::T_sp atomD = nil<core::T_O>(),
                      core::T_sp option = nil<core::T_O>() );
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
  FORWARD(EnergyFunction);
  class EnergyFunction_O : public ScoringFunction_O
  {
    LISP_CLASS(chem,ChemPkg,EnergyFunction_O,"EnergyFunction",ScoringFunction_O);
  public:
    static EnergyFunction_sp make(core::T_sp matter,
                                  core::T_sp disableComponents,
                                  core::List_sp enableComponents,
                                  bool useExcludedAtoms=false,
                                  core::T_sp keepInteraction=nil<core::T_O>(),
                                  bool assign_types=false );
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
    EnergyAngle_sp			_Angle;
    EnergyDihedral_sp			_Dihedral;
    EnergyNonbond_sp			_Nonbond;
    EnergyDihedralRestraint_sp		_DihedralRestraint;
    EnergyChiralRestraint_sp		_ChiralRestraint;
    EnergyAnchorRestraint_sp		_AnchorRestraint;
    EnergyFixedNonbondRestraint_sp	_FixedNonbondRestraint;
    core::List_sp                       _OtherEnergyComponents; // alist of additional (name . energy-function) pairs
    BoundingBox_sp                      _BoundingBox;
    double		                _DielectricConstant;
    double		                _NonbondCutoff;
    /*! If true then secondary amides are
     * automatically restrainted to be trans
     */
    bool					_RestrainSecondaryAmides;
    double					_ChiralRestraintWeight;
    double					_ChiralRestraintOffset;
    double					_AnchorRestraintWeight;
    double					_FixedNonbondRestraintWeight;
    core::T_sp                   _Message;
    core::List_sp			_MissingParameters;
  public:
  public:
    void	_eraseMissingParameters() { this->_MissingParameters = nil<core::T_O>();};
    void	_addMissingParameter(FFParameter_sp p) { this->_MissingParameters = core::Cons_O::create(p,this->_MissingParameters);};
    void __createSecondaryAmideRestraints(VectorAtom& nitrogens, core::T_sp keepInteraction );

    void	flagDihedralRestraintsAboveThreshold(NVector_sp nvPosition);
    core::T_mv enabledDisabled() const;
  private:
    int _applyRestraints(core::T_sp forceField, core::Iterator_sp restraintIterator, core::T_sp keepInteraction, core::HashTable_sp atomTypes );
    void _addDihedralRestraint(Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4, double min, double max, double weight, core::T_sp keepInteraction );

  public:

    CL_LISPIFY_NAME("atomTable");
    CL_DEFMETHOD     AtomTable_sp atomTable() const { return this->_AtomTable;};

    core::List_sp allComponents() const;

    string	energyTermsEnabled() ;
    void	loadCoordinatesIntoVector(NVector_sp pos);
    void	saveCoordinatesFromVector(NVector_sp pos);
    void	saveCoordinatesAndForcesFromVectors(NVector_sp pos, NVector_sp force);
    size_t	getNVectorSize() const;
    double	evaluateRaw( NVector_sp pos, NVector_sp force );
//    double	evaluate( NVector_sp pos, NVector_sp force, bool calculateForce );
    adapt::QDomNode_sp	identifyTermsBeyondThreshold();
//    uint	countBadVdwInteractions(double scaleSumOfVdwRadii, geom::DisplayList_sp displayIn);

    BoundingBox_sp boundingBox() const;
    bool boundingBoxBoundP() const;
    void setBoundingBox(BoundingBox_sp bounding_box);
    void makUnboundBoundingBox();
    
    ForceMatchReport_sp checkIfAnalyticalForceMatchesNumericalForce( NVector_sp pos, NVector_sp force, core::T_sp activeAtomMask );

    CL_LISPIFY_NAME("getMatter");
    CL_DEFMETHOD     Matter_sp	getMatter() { return this->_Matter;};

    void	useDefaultSettings();


    /*! Set a single options */
    void	setOption( core::Symbol_sp option, core::T_sp val);


    /*! Set the energy function options. List the options as a flat list of keyword/value pairs */
    void	setOptions( core::List_sp options );

    CL_LISPIFY_NAME("getStretchComponent");
    CL_DEFMETHOD     EnergyStretch_sp	getStretchComponent() { return this->_Stretch; };
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
    CL_LISPIFY_NAME("getDihedralRestraintComponent");
    CL_DEFMETHOD     EnergyDihedralRestraint_sp	getDihedralRestraintComponent() { return this->_DihedralRestraint; };
    CL_LISPIFY_NAME("getFixedNonbondRestraintComponent");
    CL_DEFMETHOD     EnergyFixedNonbondRestraint_sp	getFixedNonbondRestraintComponent() { return this->_FixedNonbondRestraint; };

    core::List_sp allEnergyComponents() const;
    
    core::List_sp otherEnergyComponents() const;
    void          pushOtherEnergyComponent(EnergyComponent_sp component);


    CL_DEFMETHOD bool hasMissingParameters();
    CL_DEFMETHOD core::List_sp getMissingParameters();

    void	setupHessianPreconditioner( NVector_sp pos, AbstractLargeSquareMatrix_sp hessian,
                                            core::T_sp activeAtomMask );

    /*! Enable debugging on all energy components
     */
    void	enableDebug();
    /*! Disable debugging on all energy components
     */
    void	disableDebug();

    void	summarizeTerms();
    void	dumpTerms();
    CL_DEFMETHOD     core::T_sp	getMessage() { return this->_Message;};

    int	compareAnalyticalAndNumericalForceAndHessianTermByTerm(NVector_sp pos);
    int	compareAnalyticalAndNumericalForceAndHessianTermByTermAtCurrentPosition();


    void		writeForceToAtoms(NVector_sp f);
    EnergyAtom*     getEnergyAtomPointer(Atom_sp a);

    void assignAtomTypes(Matter_sp matter, bool show_progress);
    void defineForMatter(Matter_sp agg, bool useExcludedAtoms, core::T_sp keepInteraction=nil<core::T_O>(), bool assign_types=true );
    void defineForMatterWithAtomTypes(Matter_sp matter, bool useExcludedAtoms, core::T_sp keepInteraction, core::T_sp cip_priorities, core::HashTable_sp atomTypes );
    void generateStandardEnergyFunctionTables(Matter_sp mol,
                                              FFStretchDb_sp stretchDb,
                                              FFAngleDb_sp angleDb,
                                              FFPtorDb_sp ptorDb,
                                              FFItorDb_sp itorDb,
                                              core::T_sp keepInteraction,
                                              core::HashTable_sp atomTypes);
    void generateNonbondEnergyFunctionTables(bool useExcludedAtoms, Matter_sp agg, core::T_sp forceField, core::T_sp keepInteraction, core::HashTable_sp atomTypes );
    void generateRestraintEnergyFunctionTables(Matter_sp agg, core::T_sp nonbonds, core::T_sp keepInteraction, core::T_sp cip_priorities, core::HashTable_sp atomTypes );

    CL_DEFMETHOD void	setDielectricConstant(double d) { this->_DielectricConstant = d; };
    CL_DEFMETHOD double	getDielectricConstant() { return this->_DielectricConstant; };
    CL_DEFMETHOD void	setNonbondCutoff(double d) { this->_NonbondCutoff = d; };
    CL_DEFMETHOD double	getNonbondCutoff() { return this->_NonbondCutoff; };

    /*! Add the restraints to the energy function.
     * This allows restraints to be applied to the system
     * without having to add them to the molecule/aggregate.
     */
    void	addTermsForListOfRestraints( ForceField_sp forceField,  core::List_sp restraintList, core::T_sp keepInteraction, core::HashTable_sp atomTypes );


    double	calculateNumericalDerivative(NVector_sp pos, double delta, uint i, core::T_sp activeAtomMask );
    double	calculateNumericalSecondDerivative(NVector_sp pos, double delta, uint i, uint j, core::T_sp activeAtomMask );
    double	evaluateAll(NVector_sp pos,
                            core::T_sp componentEnergy,
                            bool calcForce,
                            gc::Nilable<NVector_sp> force,
                            bool calcDiagonalHessian,
                            bool calcOffDiagonalHessian,
                            gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                            gc::Nilable<NVector_sp> hdvec,
                            gc::Nilable<NVector_sp> dvec,
                            core::T_sp activeAtomMask,
                            core::T_sp debugInteractions );

    string	summarizeBeyondThresholdInteractionsAsString();
    string	summarizeEnergyAsString();


//		adapt::QDomNode_sp	rawAccumulateTermsBeyondThresholdAsXml(uint& count);
//		adapt::QDomNode_sp	accumulateTermsBeyondThresholdAsXml();
    uint		countTermsBeyondThreshold();

    void	evaluateNumericalForce(NVector_sp pos, NVector_sp numForce, double delta, core::T_sp activeAtomMask );
    void	evaluateNumericalHessian(NVector_sp pos, AbstractLargeSquareMatrix_sp numHessian, bool calcOffDiagonalElements, double delta, core::T_sp activeAtomMask);

    core::List_sp checkForBeyondThresholdInteractions(double threshold);

    string	debugLogAsString();

    void	dealWithProblem(core::Symbol_sp error_symbol, core::T_sp arguments);

    EnergyFunction_sp copyFilter(core::T_sp keepInteraction);

    EnergyFunction_O(BoundingBox_sp bounding_box) :
        _Matter(unbound<Matter_O>())
        , _NonbondCrossTermTable(unbound<FFNonbondCrossTermTable_O>())
        , _AtomTable(unbound<AtomTable_O>())
        , _Stretch(unbound<EnergyStretch_O>())
        , _Angle(unbound<EnergyAngle_O>())
        , _Dihedral(unbound<EnergyDihedral_O>())
        , _Nonbond(unbound<EnergyNonbond_O>())
        , _DihedralRestraint(unbound<EnergyDihedralRestraint_O>())
        , _ChiralRestraint(unbound<EnergyChiralRestraint_O>())
        , _AnchorRestraint(unbound<EnergyAnchorRestraint_O>())
        , _FixedNonbondRestraint(unbound<EnergyFixedNonbondRestraint_O>())
        ,_OtherEnergyComponents(nil<core::T_O>())
        ,_BoundingBox(bounding_box)
        ,_DielectricConstant(1.0)
        ,_NonbondCutoff(16.0)
//      , _MissingParameters(unbound<core::List_O>())
    {};

    EnergyFunction_O() :
        _Matter(unbound<Matter_O>())
        , _NonbondCrossTermTable(unbound<FFNonbondCrossTermTable_O>())
        , _AtomTable(unbound<AtomTable_O>())
        , _Stretch(unbound<EnergyStretch_O>())
        , _Angle(unbound<EnergyAngle_O>())
        , _Dihedral(unbound<EnergyDihedral_O>())
        , _Nonbond(unbound<EnergyNonbond_O>())
        , _DihedralRestraint(unbound<EnergyDihedralRestraint_O>())
        , _ChiralRestraint(unbound<EnergyChiralRestraint_O>())
        , _AnchorRestraint(unbound<EnergyAnchorRestraint_O>())
        , _FixedNonbondRestraint(unbound<EnergyFixedNonbondRestraint_O>())
        ,_OtherEnergyComponents(nil<core::T_O>())
        ,_BoundingBox(unbound<BoundingBox_O>())
        ,_DielectricConstant(1.0)
        ,_NonbondCutoff(16.0)
//      , _MissingParameters(unbound<core::List_O>())
    {};
    EnergyFunction_O( const EnergyFunction_O& ef ) :
        _Matter(unbound<Matter_O>())
        , _NonbondCrossTermTable(unbound<FFNonbondCrossTermTable_O>())
        , _AtomTable(unbound<AtomTable_O>())
        , _Stretch(unbound<EnergyStretch_O>())
        , _Angle(unbound<EnergyAngle_O>())
        , _Dihedral(unbound<EnergyDihedral_O>())
        , _Nonbond(unbound<EnergyNonbond_O>())
        , _DihedralRestraint(unbound<EnergyDihedralRestraint_O>())
        , _ChiralRestraint(unbound<EnergyChiralRestraint_O>())
        , _AnchorRestraint(unbound<EnergyAnchorRestraint_O>())
        , _FixedNonbondRestraint(unbound<EnergyFixedNonbondRestraint_O>())
        ,_OtherEnergyComponents(nil<core::T_O>())
        ,_BoundingBox(unbound<BoundingBox_O>())
//      , _MissingParameters(unbound<core::List_O>())
    {};
  };

    bool inAtomSet(core::T_sp atomSet, Atom_sp atom);

#define	FINITE_DIFFERENCE_TOLERANCE 0.10
  int	_areValuesClose(double numVal, double analVal, const char* funcName, const char* termName, int index );

  void energyFunction_initializeSmarts();


};

#endif
