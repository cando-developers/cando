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

#ifndef SketchFunction_H
#define	SketchFunction_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/chem/scoringFunction.h>
#include <cando/chem/energySketchNonbond.h>
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
  SMART(EnergySketchStretch);
  SMART(EnergyPointToLineRestraint);
  SMART(EnergyOutOfZPlane);
  SMART(ForceMatchReport);



#define	DefaultChiralRestraintOffset	0.2
#define	DefaultChiralRestraintWeight	100000.0
#define	DefaultAnchorRestraintWeight	10.0
#define	DefaultPointToLineRestraintWeight	10.0
#define	DefaultOutOfZPlaneWeight	10.0


 /*! Throw this exception from Cando if you want the minimizer to deal with a bad interaction */
  FORWARD(SketchFunction);
};

template <>
struct gctools::GCInfo<chem::SketchFunction_O> {
  static bool constexpr NeedsInitialization = true;
  static bool constexpr NeedsFinalization = false;
  static GCInfo_policy constexpr Policy = normal;
};

namespace chem {
  SMART(SketchFunction);
  class SketchFunction_O : public ScoringFunction_O
  {
    LISP_CLASS(chem,ChemPkg,SketchFunction_O,"SketchFunction",ScoringFunction_O);
  public:
    static SketchFunction_sp make(core::T_sp graph, core::T_sp sketchNonbondForceField, core::HashTable_sp atomTypes);
  public:
    void initialize();
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    core::T_sp				_Graph;	// Graph object
    /*! Stores cross terms for evaluating nonbond interactions
     */
    core::T_sp			        _NodeTable;
    EnergySketchStretch_sp		_Stretch;
    EnergyOutOfZPlane_sp                _OutOfZPlane;
    EnergyPointToLineRestraint_sp       _PointToLineRestraint;
    EnergySketchNonbond_sp              _Nonbond;
    /*! If true then secondary amides are
     * automatically restrainted to be trans
     */
    double					_PointToLineRestraintWeight;
    double					_OutOfZPlaneWeight;
    double					_TotalEnergy;
    core::SimpleBaseString_sp			_Message;
    core::T_sp                                  _Frozen;
  public:
    SketchFunction_O(core::T_sp graph) : _Graph(graph), _Frozen(nil<core::T_O>()) {};
  public:
    CL_LISPIFY_NAME("nodeTable");
    CL_DEFMETHOD     core::T_sp nodeTable() const { return this->_NodeTable;};
    AtomTable_sp atomTable() const;
    void setf_velocity_scale(double xscale, double yscale, double zscale);
    
    string	energyTermsEnabled() ;
    void	loadCoordinatesIntoVector(NVector_sp pos);
    void	saveCoordinatesFromVector(NVector_sp pos);
    void	saveCoordinatesAndForcesFromVectors(NVector_sp pos, NVector_sp force);
    size_t	getNVectorSize() const;
    double	evaluateRaw( NVector_sp pos, NVector_sp force );
//    double	evaluate( NVector_sp pos, NVector_sp force, bool calculateForce );
    adapt::QDomNode_sp	identifyTermsBeyondThreshold();
//    uint	countBadVdwInteractions(double scaleSumOfVdwRadii, geom::DisplayList_sp displayIn);

    ForceMatchReport_sp checkIfAnalyticalForceMatchesNumericalForce( NVector_sp pos, NVector_sp force );

    CL_LISPIFY_NAME("getGraph");
    CL_DEFMETHOD     core::T_sp	getGraph() { return this->_Graph;};
    virtual Matter_sp getMatter();
    
    void	useDefaultSettings();

    core::List_sp allComponents() const;


    /*! Set a single options */
    void	setOption( core::Symbol_sp option, core::T_sp val);


    /*! Set the energy function options. List the options as a flat list of keyword/value pairs */
    void	setOptions( core::List_sp options );

    CL_LISPIFY_NAME("getSketchNonbondComponent");
    CL_DEFMETHOD     EnergySketchNonbond_sp	getSketchNonbondComponent() { return this->_Nonbond; };
    CL_LISPIFY_NAME("getStretchComponent");
    CL_DEFMETHOD     EnergySketchStretch_sp	getSketchStretchComponent() { return this->_Stretch; };
    CL_LISPIFY_NAME("getPointToLineRestraintComponent");
    CL_DEFMETHOD     EnergyPointToLineRestraint_sp	getPointToLineRestraintComponent() { return this->_PointToLineRestraint; };
    CL_LISPIFY_NAME("getOutOfZPlaneComponent");
    CL_DEFMETHOD     EnergyOutOfZPlane_sp	getOutOfZPlaneComponent() { return this->_OutOfZPlane; };

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
    void	dumpTerms(core::HashTable_sp atomTypes);
    CL_DEFMETHOD     core::T_sp	getMessage() { return this->_Message;};

    /*! Print the energy components as a single, multi-line string
     */
    string	energyComponentsAsString();

//    void		writeForceToAtoms(NVector_sp f);
//    EnergyAtom*     getEnergyAtomPointer(Atom_sp a);

    double	calculateNumericalDerivative(NVector_sp pos, double delta, uint i );
    double	calculateNumericalSecondDerivative(NVector_sp pos, double delta, uint i, uint j );
    double	evaluateAll(NVector_sp pos,
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

    string	debugLogAsString();

    virtual void	dealWithProblem(core::Symbol_sp error_symbol, core::T_sp arguments);
    void resetSketchFunction();


  };


  bool inAtomSet(core::T_sp atomSet, Atom_sp atom);

  void energyFunction_initializeSmarts();

};

#endif
