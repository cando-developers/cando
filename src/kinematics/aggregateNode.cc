/*
    File: aggregateNode.cc
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

#include <clasp/core/record.h>
#include <clasp/core/environment.h>
#include <cando/kinematics/kinematicsPackage.h>
#include <cando/kinematics/chainNode.h>
#include <cando/kinematics/aggregateNode.h>
#include <clasp/core/wrappers.h>
namespace kinematics
{

// ----------------------------------------------------------------------
//


void AggregateNode_O::fields(core::Record_sp node)
{
  node->field_if_not_empty(INTERN_(kw,chains),this->_Chains);
  this->Base::fields(node);
}


void AggregateNode_O::initialize()
{_OF();
  this->Base::initialize();
}

void AggregateNode_O::resizeChains(int numChains)
{_OF();
  int numMolsOld = this->_Chains.size();
  this->_Chains.resize(numChains);
  if ( numChains > numMolsOld )
  {
    for ( int i=numMolsOld; i<numChains; i++ )
    {
      this->_Chains[i] = ChainNode_O::create();
    }
  }
}

    

}; /* kinematics */
