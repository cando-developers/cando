/*
    File: scoringFunction.h
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
 *	scoringFunction.h
 *
 *	Maintains a database of stretch types
 */

#ifndef ScoringFunction_H
#define	ScoringFunction_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/chem/atom.h>
#include <cando/chem/nVector.h>




#include <cando/adapt/quickDom.fwd.h>// scoringFunction.h wants QDomNode needs quickDom.fwd.h
//#include "geom/render.fwd.h"// scoringFunction.h wants DisplayList needs render.fwd.h
#include <clasp/core/iterator.fwd.h>// scoringFunction.h wants Iterator needs iterator.fwd.h


#include <cando/chem/chemPackage.h>

namespace chem {
  SMART(AbstractLargeSquareMatrix);
  FORWARD(ScoringFunction);
};

template <>
struct gctools::GCInfo<chem::ScoringFunction_O> {
  static bool constexpr NeedsInitialization = false;
  static bool constexpr NeedsFinalization = false;
  static GCInfo_policy constexpr Policy = normal;
};

namespace       chem
{

  SMART(FFParameter);
  SMART(AbstractLargeSquareMatrix);
  SMART(Atom);
  SMART(Matter);
//SMART(DisplayList);
  SMART(Iterator);
  SMART(ForceMatchReport);


//
// Remove InteractionProblem
// once we have split up the scoring components
//
//
  SMART(ScoringFunction);
  class ScoringFunction_O : public core::CxxObject_O
  {
    LISP_CLASS(chem,ChemPkg,ScoringFunction_O,"ScoringFunction",core::CxxObject_O);
  public:
    string	    _Name;
    bool        _Debug;
  public:
    CL_DEFMETHOD     void	setScoringFunctionName(const string& nm) { this->_Name = nm; };
    CL_DEFMETHOD     string	scoringFunctionName() { return this->_Name; };
//    void		writeForceToAtoms(NVector_sp f);

  public:
    CL_LISPIFY_NAME("energyTermsEnabled");
    CL_DEFMETHOD virtual string	energyTermsEnabled() = 0;
 
    CL_DEFMETHOD virtual size_t	getNVectorSize() = 0;
//    virtual double	evaluate( NVector_sp pos, NVector_sp force, bool calculateForce ) = 0;

    virtual ForceMatchReport_sp checkIfAnalyticalForceMatchesNumericalForce( NVector_sp pos, NVector_sp force ) = 0;

    CL_LISPIFY_NAME("useDefaultSettings");
    CL_DEFMETHOD    virtual void	useDefaultSettings() = 0;
    /*! Set a single options */
    CL_LISPIFY_NAME("setOption");
    CL_DEFMETHOD     virtual void	setOption( core::Symbol_sp option, core::T_sp val) = 0;


    /*! Set the scoring function options. List the options as a flat list of keyword/value pairs */
    CL_LISPIFY_NAME("setOptions");
    CL_DEFMETHOD     virtual void	setOptions( core::List_sp options ) = 0;

    CL_LISPIFY_NAME("setupHessianPreconditioner");
    CL_DEFMETHOD virtual void	setupHessianPreconditioner( NVector_sp pos, AbstractLargeSquareMatrix_sp hessian) = 0;

            /*! Enable debugging on all scoring components
     */
    CL_LISPIFY_NAME("enableDebug");
    CL_DEFMETHOD virtual void	enableDebug() = 0;
    /*! Disable debugging on all scoring components
     */ 
    CL_LISPIFY_NAME("disableDebug");
    CL_DEFMETHOD virtual void	disableDebug() = 0;


//    virtual void	summarizeTerms() = 0;
    CL_LISPIFY_NAME("dumpTerms");
    CL_DEFMETHOD virtual void	dumpTerms() = 0;

    /*! Print the scoring components as a single, multi-line string
     */
//    virtual string	scoringComponentsAsString() = 0;

#if 0
    int	compareAnalyticalAndNumericalForceAndHessianTermByTerm(NVector_sp pos);
    CL_LISPIFY_NAME("compareAnalyticalAndNumericalForceAndHessianTermByTermAtCurrentPosition");
    CL_DEFMETHOD int	compareAnalyticalAndNumericalForceAndHessianTermByTermAtCurrentPosition();
#endif

    CL_LISPIFY_NAME("energyComponentsAsString");
    CL_DEFMETHOD virtual string energyComponentsAsString() = 0;
 
    CL_DEFMETHOD virtual void	loadCoordinatesIntoVector(NVector_sp pos) = 0;
    CL_LISPIFY_NAME("saveCoordinatesFromVector");
    CL_DEFMETHOD virtual void	saveCoordinatesFromVector(NVector_sp pos) = 0;
    CL_LISPIFY_NAME("saveCoordinatesAndForcesFromVectors");
    CL_DEFMETHOD virtual void	saveCoordinatesAndForcesFromVectors(NVector_sp pos, NVector_sp force) = 0;

    /*! Add the restraints to the scoring function.
     * This allows restraints to be applied to the system
     * without having to add them to the molecule/aggregate.
     */
//    void	addTermsForListOfRestraints( ForceField_sp forceField,  core::List_sp restraintList, core::T_sp activeAtoms );


#if 0
    CL_LISPIFY_NAME("calculateNumericalDerivative");
    CL_DEFMETHOD double	calculateNumericalDerivative(NVector_sp pos, double delta, uint i );
    CL_LISPIFY_NAME("calculateNumericalSecondDerivative");
    CL_DEFMETHOD double	calculateNumericalSecondDerivative(NVector_sp pos, double delta, uint i, uint j );
#endif

    double calculateEnergy();

    double calculateEnergyAndForce();

 
    CL_LISPIFY_NAME("evaluateAll");
    CL_LAMBDA(pos &optional calc-force force calc-diagonal-hessian calc-off-diagonal-hessian hessian hdvec dvec);
    CL_DEFMETHOD virtual double	evaluateAll( 	NVector_sp pos,
                                                bool calcForce,
                                                gc::Nilable<NVector_sp> force,
                                                bool calcDiagonalHessian,
                                                bool calcOffDiagonalHessian,
                                                gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                                gc::Nilable<NVector_sp> hdvec,
                                                gc::Nilable<NVector_sp> dvec	) = 0;
    virtual double	evaluateEnergy( NVector_sp pos );
    virtual double	evaluateEnergyForce( NVector_sp pos, bool calcForce, NVector_sp force );
    virtual double	evaluateEnergyForceFullHessian(NVector_sp pos,
                                                       bool calcForce, NVector_sp force,
                                                       bool calcDiagonalHessian,
                                                       bool calcOffDiagonalHessian,
                                                       AbstractLargeSquareMatrix_sp hessian );
    virtual double	evaluateEnergyForceFullHessianForDebugging();
#if 0

    string	summarizeBeyondThresholdInteractionsAsString();
    string	summarizeEnergyAsString();


//		adapt::QDomNode_sp	rawAccumulateTermsBeyondThresholdAsXml(uint& count);
//		adapt::QDomNode_sp	accumulateTermsBeyondThresholdAsXml();
    CL_LISPIFY_NAME("countTermsBeyondThreshold");
    CL_DEFMETHOD uint		countTermsBeyondThreshold();

    CL_LISPIFY_NAME("evaluateNumericalForce");
    CL_DEFMETHOD void	evaluateNumericalForce(NVector_sp pos, NVector_sp numForce, double delta );
    CL_LISPIFY_NAME("evaluateNumericalHessian");
    CL_DEFMETHOD void	evaluateNumericalHessian(NVector_sp pos, AbstractLargeSquareMatrix_sp numHessian, bool calcOffDiagonalElements, double delta);

    string	debugLogAsString();
#endif

    CL_LISPIFY_NAME("checkForBeyondThresholdInteractions");
    CL_DEFMETHOD virtual uint	checkForBeyondThresholdInteractions() = 0;

    virtual void	dealWithProblem(core::Symbol_sp error_symbol, core::T_sp arguments) = 0;

    DEFAULT_CTOR_DTOR(ScoringFunction_O);
  };

};

#endif
