
#ifndef	kinematics_atomTemplate_H
#define kinematics_atomTemplate_H

#include "core/common.h"
#include "addon/vector3.h"
#include "core/holder.h"
#include "chem/plug.h"
#include "chem/constitutionAtoms.h"
#include "kinematicsPackage.h"
#include "atomTree.fwd.h"
#include "atomTemplate.fwd.h"
#include "atom.fwd.h"
#include "chem/atomId.h"
#include "bondId.h"
#include "pool.h"


// Currently a maximum of 5 children hanging off of any atom
// means we can handle octahedral coordination but no more

namespace kinematics
{

    class DelayedBondedAtom;



    FORWARD(AtomTemplate);
    FORWARD(Checkpoint);
    class Checkpoint_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(kinematics,KinPkg,Checkpoint_O,"Checkpoint");
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
	virtual string atomName() const {_OF(); SUBCLASS_MUST_IMPLEMENT();};
    };


    class CheckpointAtom_O : public Checkpoint_O
    {
	LISP_BASE1(Checkpoint_O);
	LISP_CLASS(kinematics,KinPkg,CheckpointAtom_O,"CheckpointAtom");
#if INIT_TO_FACTORIES
    public:
	static CheckpointAtom_sp make(const string& atomName);
#endif
    public:
//	DECLARE_STANDARD_LISP_FUNCTIONS();
//	DECLARE_ARCHIVE();
	DEFAULT_CTOR_DTOR(CheckpointAtom_O);
    private:
	/*! The name of the atom that must be built (the checkpoint) before a
	  DelayedBondedAtom is built */
	string	_AtomName;
    public:
	virtual void setupDelayedBondedAtom(DelayedBondedAtom* atom) const;
	virtual string atomName() const {return this->_AtomName;};
    };


    class CheckpointOutPlugAtom_O : public Checkpoint_O
    {
	LISP_BASE1(Checkpoint_O);
	LISP_CLASS(kinematics,KinPkg,CheckpointOutPlugAtom_O,"CheckpointOutPlugAtom");
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
	virtual string atomName() const {return this->_Plug->getB0();};
};




/*! Builds Atoms within an AtomTree
 */
    class AtomTemplate_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(kinematics,KinPkg,AtomTemplate_O,"AtomTemplate");
#if INIT_TO_FACTORIES
    public:
	static AtomTemplate_sp make(const int id, const string& comment, const AtomTemplate_sp& parent);
#endif
    public:
//	DECLARE_STANDARD_LISP_FUNCTIONS();
//	DECLARE_ARCHIVE();
    protected:
	//! Point to the parent atom (also used to contruct linked list of unused PoolMembers)
	AtomTemplate_wp			_Parent;
	chem::ConstitutionAtomIndex0N	_Id;
	string				_Comment;
    private:

    public:

	/*! Return the parent */
	AtomTemplate_sp parent() const;

	/*! Return the atomName for this Template by looking up its index
	  in the ConstitutionAtoms object */
	string atomName(chem::ConstitutionAtoms_sp ca) const;

	/*! Return the ConstitutionAtomIndex0N id for the AtomTemplate */
	chem::ConstitutionAtomIndex0N id() const { return this->_Id;};

	/*! Return the comment associated with this object */
	string comment() const { return this->_Comment;};


	typedef map<core::Symbol_sp,BondId_sp> PlugNamesToBondIdMap;
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
	LISP_BASE1(AtomTemplate_O);
	LISP_CLASS(kinematics,KinPkg,BondedAtomTemplate_O,"BondedAtomTemplate");
#if INIT_TO_FACTORIES
    public:
	static BondedAtomTemplate_sp make(const chem::OutPlug_sp& outPlug);
#endif
    public:
//	DECLARE_STANDARD_LISP_FUNCTIONS();
//	DECLARE_ARCHIVE();
    public:
	typedef core::Vector0<AtomTemplate_O>	ChildList;
    protected:
	ChildList		_Children;
	double			_Distance;
	double			_Theta;
	double			_Phi;
	chem::OutPlug_sp	_OutPlug;
    public:
	chem::OutPlug_sp outPlug() const { return this->_OutPlug;};

	void addChild(AtomTemplate_sp child) {this->_Children.append(child);};
	/*! Return a Cons of children */
	core::Cons_sp children() { return this->_Children.asCons();};

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
	LISP_BASE1(BondedAtomTemplate_O);
	LISP_CLASS(kinematics,KinPkg,DelayedBondedAtomTemplate_O,"DelayedBondedAtomTemplate");
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
	LISP_BASE1(BondedAtomTemplate_O);
	LISP_CLASS(kinematics,KinPkg,RootBondedAtomTemplate_O,"RootBondedAtomTemplate");
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
