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


core::T_sp specializeKeepInteractionFactory( core::T_sp keepInteractionFactory, core::T_sp aclass );

bool skipInteraction_EnergyNonbond( core::T_sp keepInteractionFunction,
                                   Atom_sp atomA,
                                   Atom_sp atomB,
                                   core::T_sp ti3A,
                                   core::T_sp ti3B );

bool skipInteraction( core::T_sp keepInteractionFunction,
                      Atom_sp atomA = unbound<Atom_O>(),
                      Atom_sp atomB = unbound<Atom_O>(),
                      Atom_sp atomC = unbound<Atom_O>(),
                      Atom_sp atomD = unbound<Atom_O>() );
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
                                  core::T_sp disableComponents=nil<core::T_O>(),
                                  core::List_sp enableComponents=nil<core::T_O>(),
                                  bool useExcludedAtoms=false,
                                  core::T_sp keepInteractionFactory=nil<core::T_O>(),
                                  bool assign_types=false,
                                  core::T_sp setup = nil<core::T_O>());
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
    double                              _MonomerCorrectionEnergy;
    core::List_sp                       _EnergyComponents; // list of EnergyComponent_sp objects
    BoundingBox_sp                      _BoundingBox;
    /*! If true then secondary amides are
     * automatically restrainted to be trans
     */
    bool					_RestrainSecondaryAmides;
    core::T_sp                   _Message;
    core::List_sp			_MissingParameters;
  public:
  public:
    void	_eraseMissingParameters() { this->_MissingParameters = nil<core::T_O>();};
    void	_addMissingParameter(FFParameter_sp p) { this->_MissingParameters = core::Cons_O::create(p,this->_MissingParameters);};
    void __createSecondaryAmideRestraints(VectorAtom& nitrogens, core::T_sp keepInteractionFactory );

    void	flagDihedralRestraintsAboveThreshold(NVector_sp nvPosition);
    core::T_mv enabledDisabled() const;
  private:
    int _applyRestraints(core::T_sp forceField, core::Iterator_sp restraintIterator, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes );
    void _addDihedralRestraintDegrees(EnergyDihedralRestraint_sp rest, Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4, double phi0_degrees, double kdh, core::T_sp keepInteractionFactory );

  public:

    CL_LISPIFY_NAME("atomTable");
    CL_DEFMETHOD     AtomTable_sp atomTable() const { return this->_AtomTable;};

    CL_DEFMETHOD void setMonomerCorrectionEnergy(double en) { this->_MonomerCorrectionEnergy = en; };

    core::List_sp allComponents() const;

    string	energyTermsEnabled() ;
    void	loadCoordinatesIntoVector(NVector_sp pos);
    void	saveCoordinatesFromVector(NVector_sp pos);
    void	saveCoordinatesAndForcesFromVectors(NVector_sp pos, NVector_sp force);
    size_t	getNVectorSize() const;
    double	evaluateRaw( NVector_sp pos, NVector_sp force );

    BoundingBox_sp boundingBox() const;
    bool boundingBoxBoundP() const;
    void setBoundingBox(BoundingBox_sp bounding_box);
    void makUnboundBoundingBox();

    ForceMatchReport_sp checkIfAnalyticalForceMatchesNumericalForce( NVector_sp pos, core::T_sp energyScale, NVector_sp force, core::T_sp activeAtomMask );

    CL_LISPIFY_NAME("getMatter");
    CL_DEFMETHOD     Matter_sp	getMatter() { return this->_Matter;};

    void	useDefaultSettings();


    /*! Set a single options */
    void	setOption( core::Symbol_sp option, core::T_sp val);


    /*! Set the energy function options. List the options as a flat list of keyword/value pairs */
    void	setOptions( core::List_sp options );

    CL_DEFMETHOD core::T_sp findComponentOrNil(core::Symbol_sp className);
    CL_DEFMETHOD EnergyComponent_sp findComponent(core::Symbol_sp className);

    void addComponentIfMissing(EnergyComponent_sp comp);
    void ensureBaseComponents();

    CL_LISPIFY_NAME("getStretchComponent");
    CL_DEFMETHOD     EnergyStretch_sp	getStretchComponent();
    CL_LISPIFY_NAME("getAngleComponent");
    CL_DEFMETHOD     EnergyAngle_sp	getAngleComponent();
    CL_LISPIFY_NAME("getDihedralComponent");
    CL_DEFMETHOD EnergyDihedral_sp getDihedralComponent();
    CL_LISPIFY_NAME("getNonbondComponent");
    CL_DEFMETHOD EnergyNonbond_sp getNonbondComponent();
    CL_LISPIFY_NAME("getChiralRestraintComponent");
    CL_DEFMETHOD EnergyChiralRestraint_sp getChiralRestraintComponent();
    CL_LISPIFY_NAME("getAnchorRestraintComponent");
    CL_DEFMETHOD EnergyAnchorRestraint_sp getAnchorRestraintComponent();
    CL_LISPIFY_NAME("getDihedralRestraintComponent");
    CL_DEFMETHOD EnergyDihedralRestraint_sp getDihedralRestraintComponent();
    CL_LISPIFY_NAME("getFixedNonbondRestraintComponent");
    CL_DEFMETHOD EnergyFixedNonbondRestraint_sp getFixedNonbondRestraintComponent();

    core::List_sp allEnergyComponents() const;
    void          pushEnergyComponent(EnergyComponent_sp component) {
      this->_EnergyComponents = core::Cons_O::create(component,this->_EnergyComponents);
    }


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
    void defineForMatter(Matter_sp agg, bool useExcludedAtoms, core::T_sp keepInteractionFactory=nil<core::T_O>(), bool assign_types=true, core::T_sp coordinates=nil<core::T_O>(), core::T_sp setup=nil<core::T_O>() );
    void defineForMatterWithAtomTypes(Matter_sp matter, bool useExcludedAtoms, core::T_sp keepInteractionFactory, core::T_sp cip_priorities, core::HashTable_sp atomTypes, core::T_sp coordinates, core::T_sp setup );
    void generateStandardEnergyFunctionTables(Matter_sp mol,
                                              FFStretchDb_sp stretchDb,
                                              FFAngleDb_sp angleDb,
                                              FFPtorDb_sp ptorDb,
                                              FFItorDb_sp itorDb,
                                              core::T_sp keepInteractionFactory,
                                              core::HashTable_sp atomTypes);
    void generateNonbondEnergyFunctionTables(bool useExcludedAtoms, Matter_sp agg, core::T_sp forceField, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes, core::T_sp coordinates, core::T_sp setup );
    void generateRestraintEnergyFunctionTables(Matter_sp agg, core::T_sp nonbonds, core::T_sp keepInteractionFactory, core::T_sp cip_priorities, core::HashTable_sp atomTypes );

    /*! Add the restraints to the energy function.
     * This allows restraints to be applied to the system
     * without having to add them to the molecule/aggregate.
     */
    void	addTermsForListOfRestraints( ForceField_sp forceField,  core::List_sp restraintList, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes );


    double	calculateNumericalDerivative(NVector_sp pos, core::T_sp energyScale, double delta, uint i, core::T_sp activeAtomMask );
    double	calculateNumericalSecondDerivative(NVector_sp pos, core::T_sp energyScale, double delta, uint i, uint j, core::T_sp activeAtomMask );
    double	evaluateAll(NVector_sp pos,
                            core::T_sp energyScale,
                            core::T_sp energyComponents,
                            bool calcForce,
                            gc::Nilable<NVector_sp> force,
                            bool calcDiagonalHessian,
                            bool calcOffDiagonalHessian,
                            gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                            gc::Nilable<NVector_sp> hdvec,
                            gc::Nilable<NVector_sp> dvec,
                            core::T_sp activeAtomMask,
                            core::T_sp debugInteractions,
                            bool disableRestraints );

    string	summarizeEnergyAsString();


    void	evaluateNumericalForce(NVector_sp pos, core::T_sp energyScale, NVector_sp numForce, double delta, core::T_sp activeAtomMask );
    void	evaluateNumericalHessian(NVector_sp pos, core::T_sp energyScale, AbstractLargeSquareMatrix_sp numHessian, bool calcOffDiagonalElements, double delta, core::T_sp activeAtomMask);

    string	debugLogAsString();

    void	dealWithProblem(core::Symbol_sp error_symbol, core::T_sp arguments);

    EnergyFunction_sp copyFilter(core::T_sp keepInteractionFactory, core::List_sp setup );

    EnergyFunction_O(BoundingBox_sp bounding_box) :
        _Matter(unbound<Matter_O>())
        , _NonbondCrossTermTable(unbound<FFNonbondCrossTermTable_O>())
        , _AtomTable(unbound<AtomTable_O>())
        , _MonomerCorrectionEnergy(0.0)
        ,_EnergyComponents(nil<core::T_O>())
        ,_BoundingBox(bounding_box)
//      , _MissingParameters(unbound<core::List_O>())
    {};

    EnergyFunction_O() :
        _Matter(unbound<Matter_O>())
        , _NonbondCrossTermTable(unbound<FFNonbondCrossTermTable_O>())
        , _AtomTable(unbound<AtomTable_O>())
        , _MonomerCorrectionEnergy(0.0)
        ,_EnergyComponents(nil<core::T_O>())
        ,_BoundingBox(unbound<BoundingBox_O>())
//      , _MissingParameters(unbound<core::List_O>())
    {};
    EnergyFunction_O( const EnergyFunction_O& ef ) :
        _Matter(unbound<Matter_O>())
        , _NonbondCrossTermTable(unbound<FFNonbondCrossTermTable_O>())
        , _AtomTable(unbound<AtomTable_O>())
        , _MonomerCorrectionEnergy(0.0)
        ,_EnergyComponents(nil<core::T_O>())
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
