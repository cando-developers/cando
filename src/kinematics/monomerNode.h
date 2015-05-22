#ifndef	_kinematicsMonomerNode_H
#define _kinematicsMonomerNode_H

#include "core/common.h"
#include "monomerBaseNode.h"
#include "chainNode.fwd.h"
#include "chem/coupling.fwd.h"
#include "chem/candoDatabase.fwd.h"
#include "chem/monomer.fwd.h"
#include "chem/constitution.fwd.h"
#include "chem/topology.fwd.h"
#include "monomerNode.fwd.h"
#include "kinematicsPackage.h"

namespace kinematics
{

    FORWARD(MonomerNode);
class MonomerNode_O : public MonomerBaseNode_O
{
    friend class AtomTree_O;
    LISP_BASE1(MonomerBaseNode_O);
    LISP_CLASS(kinematics,KinPkg,MonomerNode_O,"MonomerNode");
#if INIT_TO_FACTORIES
 public:
    static MonomerNode_sp make();
#else
    DECLARE_INIT();
#endif
//    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(MonomerNode_O);
public:
	void initialize();
protected:	// instance variables here
	/*! Store the parent ChainNode to this one */
	MonomerNode_wp		_Parent;
	/*! Store the name of the plug to the parent */
	core::Symbol_sp		_ParentPlugName;
	/*! Store the monomerId/residueId of this MonomerNode */
	int			_MonomerId;
	/*! Store the children of this MonomerNode indexed by PlugNames */
	core::SymbolMap<MonomerNode_O>	_Children;
	/*! Store the index of the monomer */
	core::Symbol_sp		_MonomerName;
	/*! Store the Topology of the monomer */
	chem::Topology_sp	_Topology;
	/*! Store the current conformation index of this monomer */
	int			_ConformationIndex;
private:
	/*! Return the MonomerName */
	core::Symbol_sp	monomerName() const { return this->_MonomerName;};
	/*! Return the Topology */
	chem::Topology_sp topology() const { return this->_Topology;};

	/*! Describe recursively into a stringstream */
	void describeRecursivelyIntoStringStream(const string& prefix, stringstream& output) const;

	

public:
	static MonomerNode_sp create(int monomerId);
public: // Functions here
	/*! Recursively build children Nodes for the OutPlugs of this children and hook them in as my children */
	void recursivelyBuildChildren(ChainNode_sp chainNode, RingClosingMonomerMap& monomerMap, MonomerNode_sp parent,
				      chem::DirectionalCoupling_sp coupling, chem::Monomer_sp monomer);

	/*! Return the Constitution and Topology for this ChainNode
	 as a two element Cons */
	core::Cons_sp identifyConstitutionAndTopology(chem::CandoDatabase_sp db);

	/*! Return a description of this MonomerNode and its children */
	string asString() const;

};

}; /* kinematics */

TRANSLATE(kinematics::MonomerNode_O);

#endif /* _kinematicsMonomerNode_H */


