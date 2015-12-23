#ifndef	_kinematicsAggregateNode_H
#define _kinematicsAggregateNode_H

#include <clasp/core/common.h>
#include <cando/kinematics/monomerBaseNode.h>
#include <cando/kinematics/kinematicsPackage.h>

namespace kinematics
{

    FORWARD(ChainNode);

class AggregateNode_O : public MonomerBaseNode_O
{
    friend class FoldTree_O;
    LISP_CLASS(kinematics,KinPkg,AggregateNode_O,"AggregateNode",MonomerBaseNode_O);
//    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(AggregateNode_O);
public:
	void initialize();
public:	// Functions here
	/*! Resize the number of molecules */
	void resizeChains(int numChains);

private:
	// instance variables here
	/*! Point to the Chains */
	gctools::Vec0<ChainNode_sp>	_Chains;
};

}; /* kinematics */

TRANSLATE(kinematics::AggregateNode_O);

#endif /* _kinematicsAggregateNode_H */


