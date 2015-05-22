#ifndef	_kinematicsChainNode_H
#define _kinematicsChainNode_H

#include "core/common.h"
#include "monomerBaseNode.h"
#include "chem/oligomer.fwd.h"
#include "chainNode.fwd.h"
#include "chem/monomer.fwd.h"
#include "monomerNode.fwd.h"
#include "kinematicsPackage.h"

namespace kinematics
{


    class ChainNode_O : public MonomerBaseNode_O
    {
	friend class AtomTree_O;
	friend class Conformation_O;
	LISP_BASE1(MonomerBaseNode_O);
	LISP_CLASS(kinematics,KinPkg,ChainNode_O,"ChainNode");
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
	core::Vector0<MonomerNode_O>	_IndexedMonomerNodes;
    public: // Functions here

	void resizeMonomers(int numMonomers);

	/*! Build a graph of monomers according to the oligomer */
	void buildUsingOligomer(chem::Oligomer_sp oligomer, int chainId);


	/*! Build a MonomerNode or RingClosingMonomerNode depending on whether (monomer)
	  has a RingClosingOutPlug or not.  If it does have a RingClosingOutPlug then 
	  update the RingClosingMonomerMap to map the Monomer_sp to the MonomerNode_sp 
	  so that we can make RingClosing connections after everything is built */
	static MonomerNode_sp monomerNodeFactory(ChainNode_O* chainNode, RingClosingMonomerMap& ringClosingMonomerMap, chem::Monomer_sp monomer, core::Lisp_sp lisp);

	/*! Connect all of the ring closing connections */
	void makeRingClosingConnections(RingClosingMonomerMap& ringClosings);

	/*! Return the monomer with the MonomerId */
	MonomerNode_sp lookupMonomerId(int monomerId) const;

 
    };

}; /* kinematics */

TRANSLATE(kinematics::ChainNode_O);

#endif /* _kinematicsChainNode_H */


