/*
    File: atomTemplate.h
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

#ifndef	kinematics_atomTemplate_H
#define kinematics_atomTemplate_H

#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/adapt/symbolMap.h>
#include <clasp/core/holder.h>
#include <cando/chem/plug.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/kinematics/kinematicsPackage.h>
#include <cando/kinematics/atomTree.fwd.h>
#include <cando/kinematics/atomTemplate.fwd.h>
#include <cando/kinematics/atom.fwd.h>
#include <cando/chem/atomId.h>
#include <cando/kinematics/bondId.h>
#include <cando/kinematics/pool.h>


// Currently a maximum of 5 children hanging off of any atom
// means we can handle octahedral coordination but no more

namespace kinematics
{

    class DelayedBondedAtom;



    FORWARD(AtomTemplate);
    FORWARD(Checkpoint);
    class Checkpoint_O : public core::General_O
    {
	LISP_CLASS(kinematics,KinPkg,Checkpoint_O,"Checkpoint",core::General_O);
//	DECLARE_STANDARD_LISP_FUNCTIONS();
//	DECLARE_ARCHIVE();
	DEFAULT_CTOR_DTOR(Checkpoint_O);
    protected:
	core::Symbol_sp	_ConstitutionName;
	core::Symbol_sp	_TopologyName;
    public:
	static Checkpoint_sp make(const core::Symbol_sp& constitutionName,
				  const core::Symbol_sp& topologyName);
    public:
	/*! Set up the DelayedBondedAtom */
	virtual void setupDelayedBondedAtom(DelayedBondedAtom* atom) const;
	virtual core::Symbol_sp atomName() const {_OF(); SUBCLASS_MUST_IMPLEMENT();};
    };


    class CheckpointAtom_O : public Checkpoint_O
    {
	LISP_CLASS(kinematics,KinPkg,CheckpointAtom_O,"CheckpointAtom",Checkpoint_O);
#if INIT_TO_FACTORIES
    public:
	static CheckpointAtom_sp make(core::Symbol_sp atomName);
#endif
    public:
//	DECLARE_STANDARD_LISP_FUNCTIONS();
//	DECLARE_ARCHIVE();
	DEFAULT_CTOR_DTOR(CheckpointAtom_O);
    private:
	/*! The name of the atom that must be built (the checkpoint) before a
	  DelayedBondedAtom is built */
        core::Symbol_sp	_AtomName;
    public:
	virtual void setupDelayedBondedAtom(DelayedBondedAtom* atom) const;
	virtual core::Symbol_sp atomName() const {return this->_AtomName;};
    };


    class CheckpointOutPlugAtom_O : public Checkpoint_O
    {
	LISP_CLASS(kinematics,KinPkg,CheckpointOutPlugAtom_O,"CheckpointOutPlugAtom",Checkpoint_O);
#if INIT_TO_FACTORIES
    public:
	static CheckpointOutPlugAtom_sp make(const chem::OutPlug_sp& outPlug);
#endif
    public:
//	DECLARE_STANDARD_LISP_FUNCTIONS();
//	DECLARE_ARCHIVE();
	CheckpointOutPlugAtom_O() : _Plug(_Nil<chem::OutPlug_O>()) {};
	virtual ~CheckpointOutPlugAtom_O() {};
    private:
	/*! The name of the plug on the other side of which the Bond1 atom
	  must be built (the checkpoint) before a
	  DelayedBondedAtom is built */
	chem::OutPlug_sp	_Plug;
    public:
	virtual void setupDelayedBondedAtom(DelayedBondedAtom* atom) const;
	virtual core::Symbol_sp atomName() const {return this->_Plug->getB0();};
};




/*! Builds Atoms within an AtomTree
 */
    class AtomTemplate_O : public core::General_O
    {
	LISP_CLASS(kinematics,KinPkg,AtomTemplate_O,"AtomTemplate",core::General_O);
#if INIT_TO_FACTORIES
    public:
	static AtomTemplate_sp make(const int id, const string& comment, const AtomTemplate_sp parent);
#endif
    public:
//	DECLARE_STANDARD_LISP_FUNCTIONS();
//	DECLARE_ARCHIVE();
    protected:
	//! Point to the parent atom (also used to contruct linked list of unused PoolMembers)
        gc::Nilable<AtomTemplate_sp>			_Parent;
	chem::ConstitutionAtomIndex0N	_Id;
	string				_Comment;
    private:

    public:

	/*! Return the parent */
        gc::Nilable<AtomTemplate_sp> parent() const;

	/*! Return the atomName for this Template by looking up its index
	  in the ConstitutionAtoms object */
        core::Symbol_sp atomName(chem::ConstitutionAtoms_sp ca) const;

	/*! Return the ConstitutionAtomIndex0N id for the AtomTemplate */
	chem::ConstitutionAtomIndex0N id() const { return this->_Id;};

	/*! Return the comment associated with this object */
	string comment() const { return this->_Comment;};


	typedef adapt::SymbolMap<BondId_O> PlugNamesToBondIdMap;
	virtual Atom_sp writeIntoAtomTree(const AtomTree_sp& atomTree,
					  uint moleculeId,
					  uint residueId,
					  const BondId_sp& incoming,
					  const PlugNamesToBondIdMap& outgoing,
					  bool rootNode = false )
	{_OF();
	    SUBCLASS_MUST_IMPLEMENT();
	}


	/*! Extract the internal coordinates from the atom */
	virtual void extractInternalCoords(Atom_sp const& atom) {THROW_HARD_ERROR(BF("Subclass must implement"));};

	AtomTemplate_O() : _Parent(_Nil<AtomTemplate_O>()), _Id(-1), _Comment("") {};
	virtual ~AtomTemplate_O() {};
    };


    /*!  A template that builds a BondedAtom within an AtomTree
     */
    class BondedAtomTemplate_O : public AtomTemplate_O
    {
	LISP_CLASS(kinematics,KinPkg,BondedAtomTemplate_O,"BondedAtomTemplate",AtomTemplate_O);
#if INIT_TO_FACTORIES
    public:
	static BondedAtomTemplate_sp make(const chem::OutPlug_sp& outPlug);
#endif
    public:
//	DECLARE_STANDARD_LISP_FUNCTIONS();
//	DECLARE_ARCHIVE();
    public:
	typedef gctools::Vec0<AtomTemplate_sp>	ChildList;
    protected:
	ChildList		_Children;
	double			_Distance;
	double			_Theta;
	double			_Phi;
	chem::OutPlug_sp	_OutPlug;
    public:
	chem::OutPlug_sp outPlug() const { return this->_OutPlug;};

	void addChild(AtomTemplate_sp child) {this->_Children.push_back(child);};
	/*! Return a Cons of children */
	core::List_sp children();

	void addChildren(Atom_sp& me,
			 uint moleculeId,
			 uint residueId,
			 const AtomTree_sp& atomTree,
			 const BondId_sp& incoming,
			 const PlugNamesToBondIdMap& outgoing);


	virtual Atom_sp writeIntoAtomTree(const AtomTree_sp& atomTree,
					  uint moleculeId,
					  uint residueId,
					  const BondId_sp& incoming,
					  const PlugNamesToBondIdMap& outgoing,
					  bool rootNode = false);
	void setupOutPlugAtomTree(Atom_sp owned,
				  const AtomTree_sp& atomTree,
				  uint moleculeId,
				  uint residueId,
				  const BondId_sp& incoming,
				  const PlugNamesToBondIdMap& outgoing );


	/*! Extract the internal coordinates from the atom */
	virtual void extractInternalCoords(Atom_sp const& atom);



//	DEFAULT_CTOR_DTOR(BondedAtomTemplate_O);
	BondedAtomTemplate_O() : _Distance(0.0), _Theta(0.0), _Phi(0.0), _OutPlug(_Nil<chem::OutPlug_O>()) {};
    };


    //! Builds DelayedBondedAtom within an AtomTree.
    /*! DelayedBondedAtoms don't update their children until the second pass of
    the build process */
    class DelayedBondedAtomTemplate_O : public BondedAtomTemplate_O
    {
	LISP_CLASS(kinematics,KinPkg,DelayedBondedAtomTemplate_O,"DelayedBondedAtomTemplate",BondedAtomTemplate_O);
#if INIT_TO_FACTORIES
    public:
	static DelayedBondedAtomTemplate_sp make(const Checkpoint_sp& checkpoint);
#endif
    public:
//	DECLARE_STANDARD_LISP_FUNCTIONS();
//	DECLARE_ARCHIVE();
    protected:
	/*! Defines the atom, either in the current residue or in the following
	  one that must be built before this one can create a stub that it gives
	  to its children */
	Checkpoint_sp	_Checkpoint;
    public:

	virtual Atom_sp writeIntoAtomTree(const AtomTree_sp& atomTree,
					  uint moleculeId,
					  uint residueId,
					  const BondId_sp& incoming,
					  const PlugNamesToBondIdMap& outgoing,
					  bool rootNode = false);



//	DEFAULT_CTOR_DTOR(DelayedBondedAtomTemplate_O);	
	DelayedBondedAtomTemplate_O() : _Checkpoint(_Nil<CheckpointAtom_O>()) {};
	virtual ~DelayedBondedAtomTemplate_O() {};
    };



    /*! A template that builds an atom that will be the root of a residue
      within an AtomTree.
      The type of atom that is built depends on the Parent of the new atom.
      If the parent is a JumpAtom then a JumpAtom is created.
      If the parent is a BondedAtom then a RootBondedAtom is created.
    In both cases, information about the residue is written into the Atom
    that is created for use by the AtomTree as Atoms are built. */
    class RootBondedAtomTemplate_O : public BondedAtomTemplate_O
    {
	LISP_CLASS(kinematics,KinPkg,RootBondedAtomTemplate_O,"RootBondedAtomTemplate",BondedAtomTemplate_O);
#if INIT_TO_FACTORIES
    public:
	static RootBondedAtomTemplate_sp make(const core::Symbol_sp& constitutionName, const core::Symbol_sp topologyName, const chem::Plug_sp& inPlug);
#endif
    public:
//	DECLARE_STANDARD_LISP_FUNCTIONS();
//	DECLARE_ARCHIVE();
    protected:
	core::Symbol_sp		_ConstitutionName;
	core::Symbol_sp		_TopologyName;
	/*! Can be InPlug or JumpPlug */
	chem::Plug_sp		_InPlug;
    public:
	virtual Atom_sp writeIntoAtomTree(const AtomTree_sp& atomTree,
					  uint moleculeId,
					  uint residueId,
					  const BondId_sp& incoming,
					  const PlugNamesToBondIdMap& outgoing,
					  bool rootNode = false);

	DEFAULT_CTOR_DTOR(RootBondedAtomTemplate_O);
    };








};

TRANSLATE(kinematics::Checkpoint_O);
TRANSLATE(kinematics::CheckpointAtom_O);
TRANSLATE(kinematics::CheckpointOutPlugAtom_O);
TRANSLATE(kinematics::AtomTemplate_O);
TRANSLATE(kinematics::BondedAtomTemplate_O);
TRANSLATE(kinematics::DelayedBondedAtomTemplate_O);
TRANSLATE(kinematics::RootBondedAtomTemplate_O);


#endif
