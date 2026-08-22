/*
    File: energyComponentGroup.h
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
 *	energyComponentGroup.h
 *
 *	An EnergyComponent_O that CONTAINS other EnergyComponent_O objects and forwards to them.
 *
 *	WHY
 *	  The energy function has historically held exactly one EnergyStretch, one EnergyAngle and
 *	  so on - see ensureComponent in energyFunction.h.  A blueprint fan-out materializes every
 *	  rotamer of every sidechain simultaneously and needs to evaluate ONE SLOT at a time, which
 *	  means one set of components per slot rather than one per energy function.
 *
 *	  A group is created per generation pass and handed to ensureComponent, which then scopes
 *	  its find-or-create to that group.  So a pass driven with a fresh group produces a fresh
 *	  stretch/angle/dihedral holding only that pass's terms.  The group is both the lookup
 *	  scope during generation and the handle afterwards - it is what a ROTAMER-SCAN stores.
 *
 *	EVALUATION
 *	  evaluateAllComponent sums its enabled children, so evaluating one slot is evaluating one
 *	  group, and ENABLE/DISABLE - already on EnergyComponent_O - turns a whole slot on or off.
 *
 *	  Groups ARE pushed onto the energy function's component list, so they are visible to
 *	  serialization, setRepWeight walks and everything else that iterates components.  That is
 *	  safe here only because nothing evaluates the blueprint's energy function as a whole: a
 *	  total over all groups would sum every mutually exclusive rotamer at once.  A Monte Carlo
 *	  step reads cached per-slot energies selected by MCSTATE instead.
 *
 *	FINDING COMPONENTS
 *	  Once groups exist, findComponentOrNil is no longer meaningful for the grouped classes -
 *	  it returns the FIRST component of a class, which is the backbone's.  Consumers of slot
 *	  components must hold direct references, and any walk that must reach children needs to
 *	  descend through COMPONENTS below rather than only iterating the top level.
 */

#ifndef EnergyComponentGroup_H  //[
#define	EnergyComponentGroup_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/energyComponent.h>

namespace       chem {

  FORWARD(EnergyFunction);

SMART(EnergyComponentGroup);
class EnergyComponentGroup_O : public EnergyComponent_O
{
  LISP_CLASS(chem,ChemPkg,EnergyComponentGroup_O,"EnergyComponentGroup",EnergyComponent_O);

public: // instance variables
  gctools::Vec0<EnergyComponent_sp>	_Components;
  /*! A label for this group - anything: a symbol, a string, a list.
   *
   * Groups are otherwise indistinguishable in a backtrace or an inspector: a blueprint holds one
   * backbone-intra, one backbone-inter, one pair-scan and 940 per-slot groups, and every one of
   * them prints as "#<EnergyComponentGroup N components>".  Naming them is what makes an error or
   * a profile line say WHICH group it came from. */
  core::T_sp                            _GroupName;

public:
  typedef gctools::Vec0<EnergyComponent_sp>::iterator iterator;
  iterator begin() { return this->_Components.begin(); };
  iterator end() { return this->_Components.end(); };

public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);

public:
  static EnergyComponentGroup_sp make();

  CL_DEFMETHOD void addComponent(EnergyComponent_sp c);
  /*! The children, as a list - for walkers that must descend into groups. */
  CL_DEFMETHOD core::List_sp components() const;
  CL_LISPIFY_NAME("energy-component-group-name");
  CL_DEFMETHOD core::T_sp groupName() const { return this->_GroupName; };
  CL_LISPIFY_NAME("setf-energy-component-group-name");
  CL_DEFMETHOD void setGroupName(core::T_sp name) { this->_GroupName = name; };
  CL_DEFMETHOD size_t numberOfComponents() const { return this->_Components.size(); };

  /*! TRUE if ANY child is a restraint.  A group is as much a restraint as its children. */
  virtual bool restraintp() const;
  virtual bool is_restraint() const;

  virtual std::string descriptionOfContents() const;
  virtual std::string implementation_details() const;
  /*! Sum over children - NOT the number of children. */
  virtual size_t numberOfTerms();

  virtual void dumpTerms(core::HashTable_sp atomTypes);
  virtual void atomsForEachTerm(core::Function_sp callback);
  virtual core::List_sp extract_vectors_as_alist() const;

  virtual size_t runTestCalls(core::T_sp stream, chem::NVector_sp pos) const;
  virtual void emitTestCalls(core::T_sp stream, chem::NVector_sp pos) const;

  virtual void setupHessianPreconditioner(NVector_sp nvPosition,
                                          AbstractLargeSquareMatrix_sp m,
                                          core::T_sp activeAtomMask );

  /*! Sum the ENABLED children.  A disabled child contributes nothing, which is how one slot is
      evaluated without touching the others. */
  virtual double evaluateAllComponent( ScoringFunction_sp scorer,
                                       NVector_sp 	pos,
                                       core::T_sp energyScale,
                                       core::T_sp componentEnergy,
                                       bool 		calcForce,
                                       gc::Nilable<NVector_sp> 	force,
                                       bool		calcDiagonalHessian,
                                       bool		calcOffDiagonalHessian,
                                       gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                       gc::Nilable<NVector_sp>	hdvec,
                                       gc::Nilable<NVector_sp> dvec,
                                       core::T_sp activeAtomMask,
                                       core::T_sp debugInteractions );

  virtual void compareAnalyticalAndNumericalForceAndHessianTermByTerm(NVector_sp pos );

  EnergyComponent_sp copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc);

  virtual string __repr__() const;

public:
  EnergyComponentGroup_O( const EnergyComponentGroup_O& ss ); //!< Copy constructor

  EnergyComponentGroup_O() : _GroupName(nil<core::T_O>()) {};
  virtual ~EnergyComponentGroup_O() {};
};

};

#endif //]