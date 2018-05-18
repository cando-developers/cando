/*
    File: foldTree.h
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
class FoldTree_O : public core::General_O
{
    LISP_CLASS(kinematics,KinPkg,FoldTree_O,"FoldTree",core::General_O);
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
