/*
    File: chainNode.cc
*/
/*
Open Source License
Copyright (c) 2016, Christian E. Schafmeister
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 
This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
#define	DEBUG_LEVEL_NONE

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
#include <cando/chem/monomer.h>
#include <cando/chem/coupling.h>
#include <cando/chem/oligomer.h>
#include <cando/kinematics/monomerNode.h>
#include <cando/kinematics/ringClosingMonomerNode.h>
#include <cando/kinematics/chainNode.h>
#include <clasp/core/wrappers.h>
namespace kinematics
{

// ----------------------------------------------------------------------
//





#if INIT_TO_FACTORIES

  ChainNode_sp ChainNode_O::make()
    {
      IMPLEMENT_ME();
    };

#else
    core::T_sp ChainNode_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {
//      this->Base::__init__(exec,args,env,lisp);
//      arg = translate::from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	return _Nil<core::T_O>();
    }
#endif

#ifdef XML_ARCHIVE
    void ChainNode_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
	IMPLEMENT_ME();
    }
#endif


    void ChainNode_O::initialize()
    {_OF();
        this->Base::initialize();
    }


    MonomerNode_sp ChainNode_O::monomerNodeFactory(ChainNode_sp chainNode,
						   RingClosingMonomerMap ringClosingMonomerMap,
						   chem::Monomer_sp monomer )
    {
	int seqId = monomer->getSequenceNumber();
	if ( monomer->hasRingClosingOutPlug() )
	{
	    RingClosingMonomerNode_sp monomerNode = RingClosingMonomerNode_O::create(seqId);
	    ringClosingMonomerMap->setf_gethash(monomer,monomerNode);
	    chainNode->_IndexedMonomerNodes[seqId] = monomerNode;
	    return monomerNode;
	}
	MonomerNode_sp monomerNode = MonomerNode_O::create(seqId);
	chainNode->_IndexedMonomerNodes[seqId] = monomerNode;
	return monomerNode;
    }

    void ChainNode_O::makeRingClosingConnections(RingClosingMonomerMap ringClosings)
    {_OF();
	if ( ringClosings->size() == 0 ) return;
	IMPLEMENT_ME(); // Make connections between RingClosingMonomerNodes
    }




    void ChainNode_O::buildUsingOligomer(chem::Oligomer_sp oligomer, int chainId)
    {_OF();
	chem::Monomer_sp rootMonomer = oligomer->rootMonomer();
	RingClosingMonomerMap	ringClosingMonomerMap = core::HashTableEq_O::create_default();
	MonomerNode_sp rootMonomerNode = this->monomerNodeFactory(this->asSmartPtr(),ringClosingMonomerMap,rootMonomer);
	rootMonomerNode->recursivelyBuildChildren(this->sharedThis<ChainNode_O>(),
						  ringClosingMonomerMap,_Nil<core::T_O>(),
						  _Nil<core::T_O>(),rootMonomer);
	this->makeRingClosingConnections(ringClosingMonomerMap);
	this->_RootMonomerNode = rootMonomerNode;
    }

    void ChainNode_O::resizeMonomers(int numMonomers)
    {_OF();
	this->_IndexedMonomerNodes.resize(numMonomers);
    }

    MonomerNode_sp ChainNode_O::lookupMonomerId(int monomerId) const
    {_OF();
	ASSERTF(monomerId < (int)this->_IndexedMonomerNodes.size(),
		BF("Illegal monomerId[%d] - there are only %d monomerNodes")
		% monomerId
		% this->_IndexedMonomerNodes.size());
	return this->_IndexedMonomerNodes[monomerId];
    }
    

}; /* kinematics */
