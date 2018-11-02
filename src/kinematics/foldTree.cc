/*
    File: foldTree.cc
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
#include <clasp/core/evaluator.h>
#include <cando/kinematics/aggregateNode.h>
#include <cando/chem/oligomer.h>
#include <cando/kinematics/chainNode.h>
#include <cando/kinematics/monomerNode.h>
#include <cando/kinematics/foldTree.h>
#include <clasp/core/wrappers.h>
namespace kinematics
{

// ----------------------------------------------------------------------
//


void FoldTree_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,aggregate_node),this->_AggregateNode);
}


FoldTree_sp FoldTree_O::make()
{
  IMPLEMENT_ME();
};


void FoldTree_O::initialize()
{_OF();
  this->Base::initialize();
  this->_AggregateNode = AggregateNode_O::create();
}


int FoldTree_O::numberOfChains() const
{
  return this->_AggregateNode->_Chains.size();
}

void FoldTree_O::resizeChains(int numChains)
{
  this->_AggregateNode->resizeChains(numChains);
}

void FoldTree_O::resizeMonomers(int chainId, int numMonomers)
{_OF();
  ASSERTF(chainId<(int)this->_AggregateNode->_Chains.size(),
          BF("Illegal chainId[%d] - there are only %d chains")
          % chainId
          % this->_AggregateNode->_Chains.size());
  this->_AggregateNode->_Chains[chainId]->resizeMonomers(numMonomers);
}


ChainNode_sp FoldTree_O::buildChainUsingOligomer(int chainId, chem::Oligomer_sp oligomer)
{_OF();
  ASSERTF(chainId<this->numberOfChains(),BF("Illegal chainId[%d]") % chainId );
  this->resizeMonomers(chainId,oligomer->numberOfMonomers());
  this->_AggregateNode->_Chains[chainId]->buildUsingOligomer(oligomer,chainId);
  return this->_AggregateNode->_Chains[chainId];
}
	
ChainNode_sp FoldTree_O::getChainNode(int chainId) const
{_OF();
  ASSERTF(chainId>0 && chainId<this->numberOfChains(),
          BF("Illegal chainId[%d]") % chainId );
  return this->_AggregateNode->_Chains[chainId];
}

CL_DEFMETHOD AggregateNode_sp FoldTree_O::root() const { return this->_AggregateNode; };


CL_DOCSTRING("Return the monomer-node that corresponds to the monomer-id");
CL_DEFMETHOD MonomerNode_sp FoldTree_O::lookupMonomerId(MonomerId const& monomerId) const
{_OF();
  ASSERTF(monomerId.chainId()<(int)this->_AggregateNode->_Chains.size(),
          BF("Illegal chainId[%d] - there are only %d chains") % monomerId.chainId()
          % this->_AggregateNode->_Chains.size() );
  if (monomerId.chainId()<this->_AggregateNode->_Chains.size()) {
    ChainNode_sp chainNode = this->_AggregateNode->_Chains[monomerId.chainId()];
    MonomerNode_sp monomer = chainNode->lookupMonomerId(monomerId.monomerId());
    return monomer;
  }
  SIMPLE_ERROR(BF("Out of bound chain node id %lu must be less than %lu") % monomerId.chainId() % this->_AggregateNode->_Chains.size());
}

CL_DEFMETHOD void FoldTree_O::walkMonomerNodes(core::Function_sp callback)
{
  this->_AggregateNode->walkMonomerNodes(callback);
}

}; /* kinematics */
