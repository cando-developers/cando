#ifndef	_kinematicsAggregateNode_H
#define _kinematicsAggregateNode_H

#include "core/common.h"
#include "monomerBaseNode.h"
#include "kinematicsPackage.h"

namespace kinematics
{

    FORWARD(ChainNode);

class AggregateNode_O : public MonomerBaseNode_O
{
    friend class FoldTree_O;
    LISP_BASE1(MonomerBaseNode_O);
    LISP_CLASS(kinematics,KinPkg,AggregateNode_O,"AggregateNode");
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
	core::Vector0<ChainNode_O>	_Chains;
};

}; /* kinematics */

TRANSLATE(kinematics::AggregateNode_O);

#endif /* _kinematicsAggregateNode_H */


