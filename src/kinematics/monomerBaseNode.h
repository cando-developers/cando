#ifndef	_kinematicsMonomerBaseNode_H
#define _kinematicsMonomerBaseNode_H

#include <clasp/core/common.h>
#include <cando/kinematics/kinematicsPackage.h>
#include <cando/kinematics/monomerBaseNode.fwd.h>

namespace kinematics
{


    class MonomerBaseNode_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(kinematics,KinPkg,MonomerBaseNode_O,"MonomerBaseNode");
#if INIT_TO_FACTORIES
 public:
    static MonomerBaseNode_sp make();
#else
    DECLARE_INIT();
#endif
//    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(MonomerBaseNode_O);
public:
	void initialize();
public:
	// Functions here
private:
	// instance variables here
};

}; /* kinematics */

TRANSLATE(kinematics::MonomerBaseNode_O);

#endif /* _kinematicsMonomerBaseNode_H */


