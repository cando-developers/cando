/*
    File: monomerNode.h
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
#ifndef	_kinematicsMonomerNode_H
#define _kinematicsMonomerNode_H

#include <clasp/core/common.h>
#include <cando/adapt/symbolMap.h>
#include <cando/kinematics/monomerBaseNode.h>
#include <cando/kinematics/chainNode.fwd.h>
#include <cando/kinematics/monomerId.h>
#include <cando/chem/coupling.fwd.h>
#include <cando/chem/candoDatabase.fwd.h>
#include <cando/chem/monomer.fwd.h>
#include <cando/chem/constitution.fwd.h>
#include <cando/chem/topology.fwd.h>
#include <cando/kinematics/monomerNode.fwd.h>
#include <cando/kinematics/kinematicsPackage.h>

namespace kinematics
{
FORWARD(Joint);
FORWARD(MonomerNode);
class MonomerNode_O : public MonomerBaseNode_O
{
  LISP_CLASS(kinematics,KinPkg,MonomerNode_O,"MonomerNode",MonomerBaseNode_O);
public:
  static MonomerNode_sp make();
  MonomerNode_O() : _Parent(unbound<MonomerNode_O>()), _ParentPlugName(unbound<core::Symbol_O>()), _StereoisomerName(unbound<core::Symbol_O>()), _Topology(unbound<chem::Topology_O>()), _ConformationIndex(-1) {};
public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);	
  void initialize();
public:	// instance variables here
	/*! Store the parent ChainNode to this one */
  MonomerNode_sp		_Parent;
	/*! Store the name of the plug to the parent */
  core::Symbol_sp		_ParentPlugName;
	/*! Store the monomerId/residueId of this MonomerNode */
  MonomerId			_Id;
	/*! Store the children of this MonomerNode indexed by PlugNames */
  gc::SmallMap<core::Symbol_sp,MonomerNode_sp>	_Children;
  /*! Store the name of the stereoisomer */
  core::Symbol_sp		_StereoisomerName;
	/*! Store the Topology of the monomer */
  chem::Topology_sp	_Topology;
	/*! Store the current conformation index of this monomer */
  size_t			_ConformationIndex;
  /*! Store the JointTree joints for this MonomerNode_sp */
  gc::Vec0<Joint_sp>            _Joints;
private:
	/*! Return the MonomerName */
  core::Symbol_sp	stereoisomerName() const;
  
	/*! Return the Topology */
CL_DEFMETHOD  chem::Topology_sp topology() const { return this->_Topology;};

	/*! Describe recursively into a stringstream */
  void describeRecursivelyIntoStringStream(const string& prefix, stringstream& output) const;

  string __repr__() const;
	

public:
  static MonomerNode_sp create(MonomerId monomerId);
public: // Functions here
	/*! Recursively build children Nodes for the OutPlugs of this children and hook them in as my children */
  void recursivelyBuildChildren(ChainNode_sp chainNode, RingClosingMonomerMap monomerMap,
                                gctools::Nilable<MonomerNode_sp> parent,
                                gctools::Nilable<chem::DirectionalCoupling_sp> coupling,
                                chem::Monomer_sp monomer);

	/*! Return the Constitution and Topology for this ChainNode
	 as a two element Cons */
  chem::Constitution_mv identifyConstitutionAndTopology();

	/*! Return a description of this MonomerNode and its children */
  string asString() const;
  void addJoint(size_t index, Joint_sp joint);
  Joint_sp jointAt(size_t index) const;

  size_t conformationIndex() const;
  void setConformationIndex(size_t i);
  core::T_sp parent() const;
  core::Symbol_sp parentPlugName() const;
  
  core::List_sp children() const;

  MonomerId monomerId() const;

  void walkJoints(core::Function_sp callback);
  
};

}; /* kinematics */

TRANSLATE(kinematics::MonomerNode_O);

#endif /* _kinematicsMonomerNode_H */
