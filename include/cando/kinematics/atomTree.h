/*
    File: atomTree.h
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
#ifndef kinematics_atomTree_H
#define kinematics_atomTree_H

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <cando/chem/oligomer.fwd.h>
#include <cando/kinematics/kinematicsPackage.h>
#include <cando/kinematics/pool.h>
#include <cando/chem/atomIdMap.h>
#include <cando/kinematics/atomTree.fwd.h>
#include <cando/kinematics/bondId.fwd.h>
#include <cando/chem/plug.fwd.h>
#include <cando/kinematics/atom.fwd.h>
#include <cando/kinematics/bondedAtom.h>
#include <cando/kinematics/rootBondedAtom.h>
#include <cando/kinematics/chainNode.fwd.h>
#include <cando/kinematics/delayedBondedAtom.h>
#include <cando/kinematics/jumpAtom.h>
#include <cando/kinematics/originJumpAtom.h>
#include <cando/chem/atomId.h>
#include <cando/kinematics/bondId.h>


namespace kinematics
{

    FORWARD(AtomTemplate);

    class AtomHolder : public PoolMember
    {
	friend class AtomTree_O;
	friend class Atom_O;
	friend class Atom;
	friend class WeakAtomHandle;
	friend class RefCountedAtomHandle;
    private:
	//! Number of references to this resource
	uint		_RefCountOrNextUnused;
	//! The type of the resource
	NodeType	_Type;
	//! The index into the Pool
	uint		_NodeIndex;
    public:
	virtual void setNextUnusedMember(uint idx) { this->_RefCountOrNextUnused = idx;};
	virtual uint nextUnusedMember() const { return this->_RefCountOrNextUnused;};
    public:

    AtomHolder() : _Type(unused),_NodeIndex(UndefinedUnsignedInt) {};

	


	void setRefCount(uint c) { this->_RefCountOrNextUnused = c;};

        uint nodeIndex() const { return this->_NodeIndex;};
	int refCount() const { return this->_RefCountOrNextUnused;};
	void refAdd()
	{
	    this->_RefCountOrNextUnused++;
	}

	/*! Dec the reference count, if the RefCount got to zero then return true */
	bool refDec()
	{
	    this->_RefCountOrNextUnused--;
	    return this->_RefCountOrNextUnused==0;
	}

	string typeAsString() const;

    };





    class ExecutableAtomTreeWalkFunctor : public AtomTreeWalkFunctor
    {
    private:
	core::Function_sp	_Callback;
    public:
	ExecutableAtomTreeWalkFunctor(core::Function_sp callback) :
	AtomTreeWalkFunctor(lisp), _Callback(callback) {};
	void operator()(kinematics::Atom_sp atom) const;
    };



    class AtomTree_O : public core::General_O
    {
	friend class Atom_O;
	friend class Atom;
	friend class RefCountedAtomHandle;
	friend class WeakAtomHandle;
	LISP_CLASS(kinematics,KinPkg,AtomTree_O,"AtomTree",core::General_O);
    public:
	void initialize();
    private:
	//! The root of the tree
	Atom_sp				_Root;
	/*! chem::AtomIdMap maps AtomId's to Atoms in the tree */
	chem::AtomIdMap<WeakAtomHandle>	_AtomMap;

	//------- The following store atoms in pools for better cache performance -------

	//! All Atoms are indirectly referenced from this table that is indexed with WeakAtomHandle
	Pool<AtomHolder>	_AtomHolders;
	//! All JumpAtoms are stored in this vector
	Pool<JumpAtom>		_JumpAtoms;
	//! All OriginJumpAtoms are stored in this vector
	Pool<OriginJumpAtom>	_OriginJumpAtoms;
	//! All BondedAtoms are stored in this vector
	Pool<BondedAtom>	_BondedAtoms;
	//! All RootBondedAtoms are stored in this vector
	Pool<RootBondedAtom>	_RootBondedAtoms;
	//! All BondedAtoms are stored in this vector
	Pool<DelayedBondedAtom>	_DelayedBondedAtoms;
    public:
    AtomTree_O() :
	_AtomHolders(_lisp,"AtomHolders"),
	_JumpAtoms(_lisp,"JumpAtoms"),
	_OriginJumpAtoms(_lisp,"OriginJumpAtoms"),
	_BondedAtoms(_lisp,"BondedAtoms"),
	_RootBondedAtoms(_lisp,"RootBondedAtoms"),
	_DelayedBondedAtoms(_lisp,"DelayedBondedAtoms") {};

    protected:
	/*! Release the AtomHandle pointed to by handle */
	void releaseHandle(RefCountedAtomHandle& handle);
    public:


	/*! Lookup the atom by atomId */
	RefCountedAtomHandle lookup(const chem::AtomId& atomId) const;

	/*! Lookup the atom by atomId and return it as a shared Atom_O */
	Atom_sp atomTreeLookupAtomid(const chem::AtomId& atomId) const;

	int numberOfMolecules() const { return this->_AtomMap.numberOfMolecules();};

	/*! Given the molecule Id return the number of residues in it */
	int numberOfResiduesInMolecule(int mid) const { return this->_AtomMap.numberOfResidues(mid);};

	/*! Given the molecule Id and residue Id return the number of atoms in it */
	int numberOfAtomsInResidueInMolecule(int mid,int rid) const { return this->_AtomMap.numberOfAtoms(mid,rid);};


	/*! Return the total number of entries in the atom tree */
	int numberOfEntries() const { return this->_AtomHolders.size();};

	/*! Return the total number of entries in the atom tree */
	int numberOfAtoms() const {return this->_AtomHolders.size();};


	/*! Return the root of the AtomTree */
	Atom_sp root() const;


	/*! Print a description of the kin:Atom */
	string asString() const;


	/*! Create a shared Atom_sp for the atom */
	Atom_sp _initializeNewAtom(uint holderIndex, const chem::AtomId& atomId );

	template <class T>
	    Atom_sp _newAtom(Pool<T>& pool, const chem::AtomId& atomId, const string& comment)
	{_OF();
	    LOG(BF("BEFORE_ALLOC AtomTree:%s") % this->__repr__());
	    T entry(atomId,comment);
	    AtomHolder holder;
	    holder._RefCountOrNextUnused = 0;
	    holder._Type = T::nodeType;
	    holder._NodeIndex = pool.allocate(entry);
	    uint holderIndex = this->_AtomHolders.allocate(holder);
	    LOG(BF("Allocated AtomHolder[%d] for %s") % holderIndex % comment );
	    return this->_initializeNewAtom(holderIndex,atomId);
	};


	/*! Return a new JumpAtom */
	Atom_sp newJumpAtom(const chem::AtomId& atomId, const string& comment);


	/*! Return a new OriginJumpAtom */
	Atom_sp newOriginJumpAtom(const chem::AtomId& atomId, const string& comment);


	/*! allocate a new BondedAtom */
	Atom_sp newBondedAtom(const chem::AtomId& atomId, const string& comment);

	/*! allocate a new RootBondedAtom */
	Atom_sp newRootBondedAtom(const chem::AtomId& atomId, const string& comment,
				  core::Symbol_sp constitutionName,
				  core::Symbol_sp topologyName,
				  chem::Plug_sp inPlug );

	/*! allocate a new DelayedBondedAtom */
	Atom_sp newDelayedBondedAtom(const chem::AtomId& atomId, const string& comment);


	/*! Resize the number of molecules */
	void resizeMolecules(int numMolecules);

	/*! Resize the number of residues in a molecule */
	void resizeResidues(int molecule, int numResidues);

	/*! Resize the number of atoms in a residue in a molecule */
	void resizeAtoms(int molecule, int residue, int numAtoms );


	/*! Update the mapping between the AtomId and the Atom */
	void updateAtomIdMap(const chem::AtomId& atomId, const WeakAtomHandle& atomHandle );



	/*! Replace an existing sub-tree with a sub-tree derived from an AtomTemplate that
	  is the root of a new sub-tree.  (incoming) represents a bond within the AtomTree
	  whose child atom  will be replaced by the new sub-tree.
	  This function will identify all of the connections going out of the monomer sub-tree
	  being replaced and hook them into the new sub-tree.
	  There are a couple of special instances.
	  1) (incoming) is nil - in that case the new sub-tree
	  is being inserted as the root of the AtomTree.
	  2) (incoming) parent is nil - in that case we are replacing
	  the root monomer sub-tree of the AtomTree
	  3) (incoming) parent is defined but child is nil 
	      - in that case we are appending a sub-tree to an existing sub-tree.
	  4) (incoming) parent is defined and child is defined
	      - in that case we are replacing a sub-tree with the new monomer sub-tree
	*/
	void recursivelyBuildMolecule(int moleculeId,
				      int residueId,
				      chem::CandoDatabase_sp db,
				      MonomerNode_sp chainNode,
				      Atom_sp parent,
				      bool rootNode = false);

	/*! Build a Molecule according to the plan in ChainNode */
	void buildMoleculeUsingChainNode( int moleculeId,
					     ChainNode_sp chainNode,
					     chem::Oligomer_sp oligomer);

	/*! Dump the AtomHolders */
	string __repr__() const;


	/*! Walk the AtomTree and for each atom evaluate the Executable
	  with the atom as an argument */
	void walk(core::Function_sp exec);

	/*! Walk the AtomTree and for each atom evaluate the Functor
	  with the atom as an argument */
	void walkTree(AtomTreeWalkFunctor& functor);


	/*! Update internal coordinates recursively */
	void updateInternalCoords();

	/*! Walk the atoms in a Residue and evaluate the executable for each Atom
	 The executable should evaluate a single argument which is the atom */
	void walkResidue(int residueId, Atom_sp const& rootAtom, core::Function_sp exec);


    }; // AtomTree





};


TRANSLATE(kinematics::AtomTree_O);


#endif
