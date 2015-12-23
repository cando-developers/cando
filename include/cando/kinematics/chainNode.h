#ifndef	_kinematicsChainNode_H
#define _kinematicsChainNode_H

#include <clasp/core/common.h>
#include <cando/kinematics/monomerBaseNode.h>
#include <cando/chem/oligomer.fwd.h>
#include <cando/kinematics/chainNode.fwd.h>
#include <cando/chem/monomer.fwd.h>
#include <cando/kinematics/monomerNode.fwd.h>
#include <cando/kinematics/kinematicsPackage.h>

namespace kinematics
{


    class ChainNode_O : public MonomerBaseNode_O
    {
	friend class AtomTree_O;
	friend class Conformation_O;
	LISP_CLASS(kinematics,KinPkg,ChainNode_O,"ChainNode",MonomerBaseNode_O);
#if INIT_TO_FACTORIES
    public:
	static ChainNode_sp make();
#else
	DECLARE_INIT();
#endif
//	DECLARE_ARCHIVE();
	DEFAULT_CTOR_DTOR(ChainNode_O);
    public:
	void initialize();
    private:	// instance variables here
	/*! Store the root MonomerNode */
	MonomerNode_sp		_RootMonomerNode;
	/*! monomerId index into the MonomerNodes */
	gctools::Vec0<MonomerNode_sp>	_IndexedMonomerNodes;
    public: // Functions here

	void resizeMonomers(int numMonomers);

	/*! Build a graph of monomers according to the oligomer */
	void buildUsingOligomer(chem::Oligomer_sp oligomer, int chainId);


	/*! Build a MonomerNode or RingClosingMonomerNode depending on whether (monomer)
	  has a RingClosingOutPlug or not.  If it does have a RingClosingOutPlug then 
	  update the RingClosingMonomerMap to map the Monomer_sp to the MonomerNode_sp 
	  so that we can make RingClosing connections after everything is built */
	static MonomerNode_sp monomerNodeFactory(ChainNode_sp chainNode, RingClosingMonomerMap ringClosingMonomerMap, chem::Monomer_sp monomer);

	/*! Connect all of the ring closing connections */
	void makeRingClosingConnections(RingClosingMonomerMap ringClosings);

	/*! Return the monomer with the MonomerId */
	MonomerNode_sp lookupMonomerId(int monomerId) const;

 
    };

}; /* kinematics */

TRANSLATE(kinematics::ChainNode_O);

#endif /* _kinematicsChainNode_H */


