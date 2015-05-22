/* @(#)chainNode.fwd.h
 */

#ifndef _kinematics_ChainNode_fwd_H
#define _kinematics_ChainNode_fwd_H

#include "chem/monomer.fwd.h"
#include "monomerNode.fwd.h"

namespace kinematics
{
    FORWARD(ChainNode);

    typedef map<chem::Monomer_sp,MonomerNode_sp>	RingClosingMonomerMap;


}; /* kinematics */

#endif /* _kinematics_ChainNode_fwd_H */
