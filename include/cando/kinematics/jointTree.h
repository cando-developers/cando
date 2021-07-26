/*
    File: JointTree.h
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
#ifndef kinematics_JointTree_H
#define kinematics_JointTree_H

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <cando/chem/oligomer.fwd.h>
#include <cando/kinematics/kinematicsPackage.h>
#include <cando/kinematics/pool.h>
#include <cando/chem/atomIdMap.h>
#include <cando/kinematics/jointTree.fwd.h>
#include <cando/kinematics/bondId.fwd.h>
#include <cando/chem/plug.fwd.h>
#include <cando/kinematics/joint.fwd.h>
#include <cando/kinematics/monomerId.h>
#include <cando/kinematics/bondedJoint.h>
#include <cando/kinematics/rootBondedJoint.h>
#include <cando/kinematics/chainNode.fwd.h>
#include <cando/kinematics/delayedBondedJoint.h>
#include <cando/kinematics/jumpJoint.h>
#include <cando/kinematics/originJumpJoint.h>
#include <cando/chem/atomId.h>
#include <cando/kinematics/bondId.h>


namespace kinematics
{

  FORWARD(JointTemplate);
  FORWARD(JumpJoint);

  class AtomHolder : public PoolMember
  {
  public:
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






  class JointTree_O : public core::CxxObject_O
  {
    LISP_CLASS(kinematics,KinPkg,JointTree_O,"JointTree",core::CxxObject_O);
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
    void initialize();
  public:
	//! The root of the tree
    Joint_sp				_Root;
	/*! chem::AtomIdMap maps AtomId's to Atoms in the tree */
    chem::AtomIdMap<Joint_sp>	_AtomMap;
        
  public:
  JointTree_O() : _Root(unbound<Joint_O>()) {};

	/*! Lookup the atom by atomId */
    Joint_sp lookup(const chem::AtomId& atomId) const;

	/*! Lookup the atom by atomId and return it as a shared Joint_O */
    Joint_sp JointTreeLookupAtomid(const chem::AtomId& atomId) const;

    size_t numberOfMolecules() const { return this->_AtomMap.numberOfMolecules();};

	/*! Given the molecule Id return the number of residues in it */
    int numberOfResiduesInMolecule(int mid) const { return this->_AtomMap.numberOfResidues(mid);};

	/*! Given the molecule Id and residue Id return the number of atoms in it */
    int numberOfAtomsInResidueInMolecule(int mid,int rid) const { return this->_AtomMap.numberOfAtoms(mid,rid);};


#if 0
	/*! Return the total number of entries in the atom tree */
    int numberOfEntries() const { return this->_AtomHolders.size();};

	/*! Return the total number of entries in the atom tree */
    int numberOfAtoms() const {return this->_AtomHolders.size();};
#endif

	/*! Return the root of the JointTree */
    CL_DEFMETHOD
    Joint_sp jointTreeRoot() const { return this->_Root; };


	/*! Print a description of the kin:Atom */
    string asString() const;


	/*! Create a shared Joint_sp for the atom */
    Joint_sp _initializeNewAtom(Joint_sp atom, const chem::AtomId& atomId );

    template <class T>
      Joint_sp _newAtom(const chem::AtomId& atomId, core::T_sp name, const string& comment)
    {_OF();
      auto atom = gctools::GC<T>::allocate(atomId,name,comment);
      return this->_initializeNewAtom(atom,atomId);
    };


	/*! Return a new JumpJoint */
    Joint_sp newJumpJoint(const chem::AtomId& atomId, core::T_sp name, const string& comment);


	/*! Return a new OriginJumpJoint */
    Joint_sp newOriginJumpJoint(const chem::AtomId& atomId, core::T_sp name, const string& comment);


	/*! allocate a new BondedJoint */
    Joint_sp newBondedJoint(const chem::AtomId& atomId, core::T_sp name, const string& comment);

	/*! allocate a new RootBondedJoint */
    Joint_sp newRootBondedJoint(const chem::AtomId& atomId, core::T_sp name, const string& comment,
                               core::Symbol_sp constitutionName,
                               core::Symbol_sp topologyName,
                               chem::Plug_sp inPlug );

	/*! allocate a new DelayedBondedJoint */
    Joint_sp newDelayedBondedJoint(const chem::AtomId& atomId, core::T_sp name, const string& comment);


	/*! Resize the number of molecules */
    void resizeMolecules(int numMolecules);

	/*! Resize the number of residues in a molecule */
    void resizeResidues(int molecule, int numResidues);

	/*! Resize the number of atoms in a residue in a molecule */
    void resizeAtoms(int molecule, int residue, int numAtoms );


	/*! Update the mapping between the AtomId and the Atom */
    void updateAtomIdMap(const chem::AtomId& atomId, Joint_sp atomHandle );



	/*! Replace an existing sub-tree with a sub-tree derived from an JointTemplate that
	  is the root of a new sub-tree.  (incoming) represents a bond within the JointTree
	  whose child atom  will be replaced by the new sub-tree.
	  This function will identify all of the connections going out of the monomer sub-tree
	  being replaced and hook them into the new sub-tree.
	  There are a couple of special instances.
	  1) (incoming) is nil - in that case the new sub-tree
	  is being inserted as the root of the JointTree.
	  2) (incoming) parent is nil - in that case we are replacing
	  the root monomer sub-tree of the JointTree
	  3) (incoming) parent is defined but child is nil 
	      - in that case we are appending a sub-tree to an existing sub-tree.
	  4) (incoming) parent is defined and child is defined
	      - in that case we are replacing a sub-tree with the new monomer sub-tree
	*/
    void recursivelyBuildMolecule(MonomerId monomerId,
                                  MonomerNode_sp chainNode,
                                  Joint_sp parent,
                                  bool rootNode = false);

	/*! Build a Molecule according to the plan in ChainNode */
    void buildMoleculeUsingChainNode( int moleculeId,
                                      ChainNode_sp chainNode,
                                      chem::Oligomer_sp oligomer);

	/*! Walk the JointTree and for each atom evaluate the Executable
	  with the atom as an argument */
    void walk(core::Function_sp exec);

	/*! Walk the JointTree and for each atom evaluate the Functor
	  with the atom as an argument */
    void walkTree(core::Function_sp callback_one_arg);


	/*! Update internal coordinates recursively */
    void updateInternalCoords();

	/*! Walk the atoms in a Residue and evaluate the executable for each Atom
	 The executable should evaluate a single argument which is the atom */
    void walkResidue(int residueId, Joint_sp const& rootJoint, core::Function_sp exec);


  }; // JointTree





};


TRANSLATE(kinematics::JointTree_O);


#endif
