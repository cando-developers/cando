#ifndef	_kinematicsRingClosingMonomerNode_H
#define _kinematicsRingClosingMonomerNode_H

#include <clasp/core/common.h>
#include <cando/kinematics/monomerNode.h>
#include <cando/kinematics/kinematicsPackage.h>

namespace kinematics
{

    FORWARD(RingClosingMonomerNode);
class RingClosingMonomerNode_O : public MonomerNode_O
{
    LISP_BASE1(MonomerNode_O);
    LISP_CLASS(kinematics,KinPkg,RingClosingMonomerNode_O,"RingClosingMonomerNode");
#if INIT_TO_FACTORIES
 public:
    static RingClosingMonomerNode_sp make();
#else
    DECLARE_INIT();
#endif
//    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(RingClosingMonomerNode_O);
public:
	void initialize();
public:
	static RingClosingMonomerNode_sp create(int monomerId);
private:
	/*! Store weak pointers to the other ring closing monomers indexed
	  by the name of the ring closing plug */
	core::WeakSymbolMap<RingClosingMonomerNode_O>	_RingClosingChildren;
public:
	// Functions here

};

}; /* kinematics */

TRANSLATE(kinematics::RingClosingMonomerNode_O);

#endif /* _kinematicsRingClosingMonomerNode_H */


