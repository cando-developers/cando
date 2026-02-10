/*
    File: rigidBodyEnergyFunction.h
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

#ifndef RigidBodyEnergyFunction_H
#define	RigidBodyEnergyFunction_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/chem/scoringFunction.h>
#include <cando/chem/energyRigidBodyComponent.h>




#include <cando/adapt/quickDom.fwd.h>// energyFunction.h wants QDomNode needs quickDom.fwd.h
//#include "geom/render.fwd.h"// energyFunction.h wants DisplayList needs render.fwd.h
#include <clasp/core/iterator.fwd.h>// energyFunction.h wants Iterator needs iterator.fwd.h


#include <cando/chem/chemPackage.h>

namespace       chem
{


  SMART(FFParameter);
  SMART(AbstractLargeSquareMatrix);
  SMART(FFNonbondCrossTermTable);
  SMART(QDomNode);
  SMART(Atom);
  SMART(Matter);
  SMART(ForceField);
  SMART(AtomTable);

  FORWARD(RigidBodyEnergyFunction);
  FORWARD(EnergyRigidBodyNonbond);
  FORWARD(EnergyRigidBodyStaple);
};


template <>
struct gctools::GCInfo<chem::RigidBodyEnergyFunction_O> {
  static bool constexpr NeedsInitialization = false;
  static bool constexpr NeedsFinalization = false;
  static GCInfo_policy constexpr Policy = normal;
};

namespace chem {
FORWARD(BoundingBox);
FORWARD(RigidBodyEnergyFunction);
class RigidBodyEnergyFunction_O : public ScoringFunction_O
{
  LISP_CLASS(chem,ChemPkg,RigidBodyEnergyFunction_O,"RigidBodyEnergyFunction",ScoringFunction_O);
public:
  static RigidBodyEnergyFunction_sp make(size_t number_of_rigid_bodies, BoundingBox_sp boundingBox );
public:
  size_t              _RigidBodies;
  NVector_sp          _SavedCoordinates;
  core::List_sp       _Terms;
  BoundingBox_sp      _BoundingBox;
public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
public:
  BoundingBox_sp boundingBox() const;
  bool boundingBoxBoundP() const;
  void setBoundingBox(BoundingBox_sp bounding_box);
  void makUnboundBoundingBox();

    
  CL_LISPIFY_NAME("rigid-body-energy-function-add-term");
  CL_DEFMETHOD void addTerm(EnergyRigidBodyComponent_sp comp) { this->_Terms = core::Cons_O::create(comp,this->_Terms);};
    
  CL_LISPIFY_NAME("rigid-body-energy-function-terms")
  CL_DEFMETHOD core::List_sp allComponents() const { return this->_Terms;};
public:
    /*! 4 quaternion and 3 cartesian coordinates for each rigid body */
  size_t numberOfRigidBodies() const;
  virtual size_t getNVectorSize() const override { return this->_RigidBodies*7; };
  CL_DEFMETHOD NVector_sp get_coordinates() const { return this->_SavedCoordinates; };
  CL_DEFMETHOD void set_coordinates(NVector_sp coords);

  virtual string	energyTermsEnabled() override;

  virtual void	setupHessianPreconditioner( NVector_sp pos, AbstractLargeSquareMatrix_sp hessian,
                                            core::T_sp activeAtomMask ) override;

    /*! Load the coordinates in the ScoringFunction into the position vector */
  virtual void	loadCoordinatesIntoVector(NVector_sp pos);
    /*! Save the coordinates in the pos vector into the ScoringFunction */
  virtual void	saveCoordinatesFromVector(NVector_sp pos);
    /*! Save the coordinates in the pos vector and forces in force into the ScoringFunction */
  virtual void	saveCoordinatesAndForcesFromVectors(NVector_sp pos, NVector_sp force);
  virtual double	evaluateRaw( NVector_sp pos, NVector_sp force ) ;

  ForceMatchReport_sp checkIfAnalyticalForceMatchesNumericalForce( NVector_sp pos, core::T_sp energyScale, NVector_sp force, core::T_sp activeAtomMask );

  void	useDefaultSettings();

    /*! Set a single options */
  void	setOption( core::Symbol_sp option, core::T_sp val);


    /*! Set the energy function options. List the options as a flat list of keyword/value pairs */
  void	setOptions( core::List_sp options );

  double	evaluateAll(NVector_sp pos,
                            core::T_sp energyScale,
                            core::T_sp componentEnergy,
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

  void	dealWithProblem(core::Symbol_sp error_symbol, core::T_sp arguments);
  CL_LISPIFY_NAME("rigid-body-energy-function-set-position");
  CL_DEFMETHOD void setPosition(size_t index, double a, double b, double c, double d, double x, double y, double z);

  string	energyComponentsAsString();

  CL_LISPIFY_NAME("rigid-body-energy-function-get-position");
  CL_DEFMETHOD core::T_mv getPosition(size_t index);
  CL_LISPIFY_NAME("rigid-body-energy-function-normalize-position");
  CL_DEFMETHOD void normalizePosition(NVector_sp pos);

  Aggregate_sp buildPseudoAggregate(NVector_sp pos) const;
  core::SimpleVector_float_sp buildPseudoAggregateCoordinates(NVector_sp pos) const;
  size_t numberOfPoints() const;
    
  void dumpTerms();
 
  RigidBodyEnergyFunction_O(size_t number_of_rigid_bodies, BoundingBox_sp boundingBox)
      : _RigidBodies(number_of_rigid_bodies),
        _Terms(nil<core::T_O>()),
        _BoundingBox(boundingBox){};
};

};

#endif
