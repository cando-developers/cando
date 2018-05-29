/*
    File: jointTemplate.h
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

#ifndef	kinematics_jointTemplate_H
#define kinematics_jointTemplate_H

#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/adapt/symbolMap.h>
#include <cando/chem/plug.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/kinematics/kinematicsPackage.h>
#include <cando/kinematics/atomTree.fwd.h>
#include <cando/kinematics/jointTemplate.fwd.h>
#include <cando/kinematics/atom.fwd.h>
#include <cando/chem/atomId.h>
#include <cando/kinematics/bondId.h>
#include <cando/kinematics/pool.h>


// Currently a maximum of 5 children hanging off of any atom
// means we can handle octahedral coordination but no more

namespace kinematics
{
  FORWARD(DelayedBondedJoint);
  FORWARD(JointTemplate);
  FORWARD(Checkpoint);
  class Checkpoint_O : public core::CxxObject_O
  {
    LISP_CLASS(kinematics,KinPkg,Checkpoint_O,"Checkpoint",core::CxxObject_O);
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
    DEFAULT_CTOR_DTOR(Checkpoint_O);
  protected:
    core::Symbol_sp	_ConstitutionName;
    core::Symbol_sp	_TopologyName;
  public:
    static Checkpoint_sp make(const core::Symbol_sp& constitutionName,
                              const core::Symbol_sp& topologyName);
  public:
	/*! Set up the DelayedBondedAtom */
    virtual void setupDelayedBondedAtom(DelayedBondedJoint_sp atom) const;
    CL_DEFMETHOD virtual core::Symbol_sp atomName() const {_OF(); SUBCLASS_MUST_IMPLEMENT();};
  };


  class CheckpointJoint_O : public Checkpoint_O
  {
    LISP_CLASS(kinematics,KinPkg,CheckpointJoint_O,"CheckpointJoint",Checkpoint_O);
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    static CheckpointJoint_sp make(core::Symbol_sp atomName);
    DEFAULT_CTOR_DTOR(CheckpointJoint_O);
  private:
	/*! The name of the atom that must be built (the checkpoint) before a
	  DelayedBondedAtom is built */
    core::Symbol_sp	_AtomName;
  public:
    virtual void setupDelayedBondedAtom(DelayedBondedJoint_sp atom) const;
    virtual core::Symbol_sp atomName() const {return this->_AtomName;};
  };


  class CheckpointOutPlugJoint_O : public Checkpoint_O
  {
    LISP_CLASS(kinematics,KinPkg,CheckpointOutPlugJoint_O,"CheckpointOutPlugJoint",Checkpoint_O);
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    static CheckpointOutPlugJoint_sp make(chem::OutPlug_sp outPlug);
  CheckpointOutPlugJoint_O() : _Plug(_Nil<chem::OutPlug_O>()) {};
    virtual ~CheckpointOutPlugJoint_O() {};
  public:
	/*! The name of the plug on the other side of which the Bond1 atom
	  must be built (the checkpoint) before a
	  DelayedBondedAtom is built */
    chem::OutPlug_sp	_Plug;
  public:
    virtual void setupDelayedBondedAtom(DelayedBondedJoint_sp atom) const;
    virtual core::Symbol_sp atomName() const {return this->_Plug->getB0();};
  };




/*! Builds Atoms within an AtomTree
 */
  class JointTemplate_O : public core::CxxObject_O
  {
    LISP_CLASS(kinematics,KinPkg,JointTemplate_O,"JointTemplate",core::CxxObject_O);
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    static JointTemplate_sp make(const int id, const string& comment, const JointTemplate_sp parent);
  protected:
	//! Point to the parent atom (also used to contruct linked list of unused PoolMembers)
    gc::Nilable<JointTemplate_sp>			_Parent;
    chem::ConstitutionAtomIndex0N	_Id;
    string				_Comment;
  public:
	/*! Return the parent */
    gc::Nilable<JointTemplate_sp> parent() const;

	/*! Return the atomName for this Template by looking up its index
	  in the ConstitutionAtoms object */
    core::Symbol_sp atomName(chem::ConstitutionAtoms_sp ca) const;

	/*! Return the ConstitutionAtomIndex0N id for the JointTemplate */
    chem::ConstitutionAtomIndex0N id() const { return this->_Id;};

	/*! Return the comment associated with this object */
    string comment() const { return this->_Comment;};


    typedef adapt::SymbolMap<BondId_O> PlugNamesToBondIdMap;
    virtual Joint_sp writeIntoAtomTree(const AtomTree_sp& atomTree,
                                       uint moleculeId,
                                       uint residueId,
                                       const BondId_sp& incoming,
                                       const PlugNamesToBondIdMap& outgoing,
                                       bool rootNode = false )
    {_OF();
      SUBCLASS_MUST_IMPLEMENT();
    }


	/*! Extract the internal coordinates from the atom */
    virtual void extractInternalCoords(Joint_sp const& atom) {THROW_HARD_ERROR(BF("Subclass must implement"));};

  JointTemplate_O() : _Parent(_Nil<JointTemplate_O>()), _Id(-1), _Comment("") {};
    virtual ~JointTemplate_O() {};
  };


    /*!  A template that builds a BondedAtom within an AtomTree
     */
  class BondedJointTemplate_O : public JointTemplate_O
  {
    LISP_CLASS(kinematics,KinPkg,BondedJointTemplate_O,"BondedJointTemplate",JointTemplate_O);
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    static BondedJointTemplate_sp make(chem::OutPlug_sp outPlug);
  public:
    typedef gctools::Vec0<JointTemplate_sp>	ChildList;
  protected:
    ChildList		_Children;
    double			_Distance;
    double			_Theta;
    double			_Phi;
    chem::OutPlug_sp	_OutPlug;
  public:
    chem::OutPlug_sp outPlug() const { return this->_OutPlug;};

    CL_DEFMETHOD void addChild(JointTemplate_sp child) {this->_Children.push_back(child);};
	/*! Return a Cons of children */
    core::List_sp children();

    void addChildren(Joint_sp me,
                     uint moleculeId,
                     uint residueId,
                     const AtomTree_sp& atomTree,
                     const BondId_sp& incoming,
                     const PlugNamesToBondIdMap& outgoing);


    virtual Joint_sp writeIntoAtomTree(const AtomTree_sp& atomTree,
                                       uint moleculeId,
                                       uint residueId,
                                       const BondId_sp& incoming,
                                       const PlugNamesToBondIdMap& outgoing,
                                       bool rootNode = false);
    void setupOutPlugAtomTree(Joint_sp owned,
                              const AtomTree_sp& atomTree,
                              uint moleculeId,
                              uint residueId,
                              const BondId_sp& incoming,
                              const PlugNamesToBondIdMap& outgoing );


	/*! Extract the internal coordinates from the atom */
    virtual void extractInternalCoords(Joint_sp const& atom);



//	DEFAULT_CTOR_DTOR(BondedJointTemplate_O);
  BondedJointTemplate_O() : _Distance(0.0), _Theta(0.0), _Phi(0.0), _OutPlug(_Unbound<chem::OutPlug_O>()) {};
  };


    //! Builds DelayedBondedAtom within an AtomTree.
    /*! DelayedBondedAtoms don't update their children until the second pass of
    the build process */
  class DelayedBondedJointTemplate_O : public BondedJointTemplate_O
  {
    LISP_CLASS(kinematics,KinPkg,DelayedBondedJointTemplate_O,"DelayedBondedJointTemplate",BondedJointTemplate_O);
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    static DelayedBondedJointTemplate_sp make(const Checkpoint_sp& checkpoint);
  protected:
	/*! Defines the atom, either in the current residue or in the following
	  one that must be built before this one can create a stub that it gives
	  to its children */
    Checkpoint_sp	_Checkpoint;
  public:

    virtual Joint_sp writeIntoAtomTree(const AtomTree_sp& atomTree,
                                       uint moleculeId,
                                       uint residueId,
                                       const BondId_sp& incoming,
                                       const PlugNamesToBondIdMap& outgoing,
                                       bool rootNode = false);



//	DEFAULT_CTOR_DTOR(DelayedBondedJointTemplate_O);	
  DelayedBondedJointTemplate_O() : _Checkpoint(_Nil<CheckpointJoint_O>()) {};
    virtual ~DelayedBondedJointTemplate_O() {};
  };



    /*! A template that builds an atom that will be the root of a residue
      within an AtomTree.
      The type of atom that is built depends on the Parent of the new atom.
      If the parent is a JumpAtom then a JumpAtom is created.
      If the parent is a BondedAtom then a RootBondedAtom is created.
    In both cases, information about the residue is written into the Atom
    that is created for use by the AtomTree as Atoms are built. */
  class RootBondedJointTemplate_O : public BondedJointTemplate_O
  {
    LISP_CLASS(kinematics,KinPkg,RootBondedJointTemplate_O,"RootBondedJointTemplate",BondedJointTemplate_O);
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    static RootBondedJointTemplate_sp make(core::Symbol_sp constitutionName, core::Symbol_sp topologyName, chem::Plug_sp inPlug);
  protected:
    core::Symbol_sp		_ConstitutionName;
    core::Symbol_sp		_TopologyName;
	/*! Can be InPlug or JumpPlug */
    chem::Plug_sp		_InPlug;
  public:
    virtual Joint_sp writeIntoAtomTree(const AtomTree_sp& atomTree,
                                       uint moleculeId,
                                       uint residueId,
                                       const BondId_sp& incoming,
                                       const PlugNamesToBondIdMap& outgoing,
                                       bool rootNode = false);

    DEFAULT_CTOR_DTOR(RootBondedJointTemplate_O);
  };








};

#endif
