#ifndef	_kinematicsFoldTree_H
#define _kinematicsFoldTree_H

#include <clasp/core/common.h>
#include <cando/kinematics/chainNode.fwd.h>
#include <cando/kinematics/monomerIdMap.h>
#include <cando/kinematics/kinematicsPackage.h>

namespace kinematics
{

FORWARD(AggregateNode);


FORWARD(FoldTree);
class FoldTree_O : public core::T_O
{
    friend class ChainNode_O;
    LISP_BASE1(core::T_O);
    LISP_CLASS(kinematics,KinPkg,FoldTree_O,"FoldTree");
#if INIT_TO_FACTORIES
 public:
    static FoldTree_sp make();
#else
    DECLARE_INIT();
#endif
//    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(FoldTree_O);
public:
	void initialize();
public:

	/*! allocate a number of ChainNodes */
	void resizeChains(int numChains);

	/*! Return the number of Chains */
	int numberOfChains() const;

	/*! Resize the number of Monomers in a Chain */
	void resizeMonomers(int chainId, int numMonomers );

	/*! Build a Chain graph using the oligomer as a plan */
    ChainNode_sp buildChainUsingOligomer(int chainId, chem::Oligomer_sp oligomer);

    /*! Return the ChainNode for chainId */
    ChainNode_sp getChainNode(int chainId) const;

    /*! Return the ChainNode for the MonomerId */
    MonomerNode_sp lookupMonomerId(MonomerId const& monomerId) const;


private:
	/*! Point to the AggregateNode that contains the ChainNodes */
	AggregateNode_sp		_AggregateNode;

};

}; /* kinematics */

TRANSLATE(kinematics::FoldTree_O);

#endif /* _kinematicsFoldTree_H */


