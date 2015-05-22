#ifndef	_kinematicsSequenceNode_H
#define _kinematicsSequenceNode_H

#include "core/common.h"
#include "sequenceBaseNode.h"
#include "kinematicsPackage.h"

namespace kinematics
{

    FORWARD(SequenceNode);

    /*! @class This is a base class for MonomerNode, MoleculeNode, AggregateNode
      It stores the parent */
class SequenceNode_O : public MonomerBaseNode_O
{
    LISP_CLASS(KinPkg,SequenceNode_O,SequenceNode,"MonomerBaseNode_O");
    DECLARE_INIT();
    DECLARE_ARCHIVE();
    DECLARE_EXPOSE_CANDO();
    DEFAULT_CTOR_DTOR(SequenceNode_O);
public:
	void initialize();
public:
	// Functions here
protected:
	// instance variables here
};

}; /* kinematics */

TRANSLATE(kinematics::SequenceNode_O);

#endif /* _kinematicsSEQUENCENODE_H */


