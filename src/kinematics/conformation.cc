/*
    File: conformation.cc
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
#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <cando/kinematics/kinematicsPackage.h>
#include <clasp/core/environment.h>
#include <cando/kinematics/conformation.h>
#include <cando/kinematics/foldTree.h>
#include <cando/kinematics/chainNode.h>
#include <cando/kinematics/monomerNode.h>
#include <cando/kinematics/jointTree.h>
#include <cando/chem/oligomer.h>
#include <clasp/core/wrappers.h>


namespace kinematics
{


SYMBOL_EXPORT_SC_(KinPkg,KinPkg_resizeMoleculesEvent);
SYMBOL_EXPORT_SC_(KinPkg,KinPkg_buildMoleculeUsingOligomerEvent);
SYMBOL_EXPORT_SC_(ChemKwPkg,newValue);
SYMBOL_EXPORT_SC_(ChemKwPkg,oldValue);
SYMBOL_EXPORT_SC_(ChemKwPkg,moleculeId);
SYMBOL_EXPORT_SC_(ChemKwPkg,residueId);


void Conformation_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw,fold_tree),this->_FoldTree);
  node->field(INTERN_(kw,joint_tree),this->_JointTree);
}

// ----------------------------------------------------------------------
//
CL_LISPIFY_NAME(make-conformation);
CL_DEF_CLASS_METHOD Conformation_sp Conformation_O::make(core::List_sp oligomers)
{_OF();
  auto  me  = gctools::GC<Conformation_O>::allocate_with_default_constructor();
  me->resizeMolecules(core::cl__length(oligomers));
  int moleculeId = 0;
  for ( auto cur : oligomers ) {
    chem::Oligomer_sp oligomer = cur->car<chem::Oligomer_O>();
    me->buildMoleculeUsingOligomer(moleculeId,oligomer);
    ++moleculeId;
  }
  return me;
};

Conformation_sp Conformation_O::create(chem::Oligomer_sp oligomer)
{
  auto  conf  = gctools::GC<Conformation_O>::allocate_with_default_constructor();
  conf->resizeMolecules(1);
  conf->buildMoleculeUsingOligomer(0,oligomer);
  return conf;
}


void Conformation_O::initialize()
{_OF();
  this->Base::initialize();
  this->_FoldTree = FoldTree_O::create();
  this->_JointTree = JointTree_O::create();
}


int Conformation_O::numberOfMolecules() const
{_OF();
  return this->_FoldTree->numberOfChains();
}

void Conformation_O::resizeMolecules(int numMolecules)
{_OF();
  ASSERTF(numMolecules > 0, BF("You must allocate more than 0 molecules"));
  int oldNumberOfMolecules = this->_FoldTree->numberOfChains();
  this->_FoldTree->resizeChains(numMolecules);
  this->_JointTree->resizeMolecules(numMolecules);
#if 0
#if 1
  SIMPLE_ERROR(BF("What are resizeMoleculesEvents for?"));
#else
  this->notify(_sym_resizeMoleculesEvent,
               core::Cons_O::createList(kw::_sym_newValue,
                                        core::clasp_make_fixnum(numMolecules),
                                        kw::_sym_oldValue,
                                        core::clasp_make_fixnum(oldNumberOfMolecules)));
#endif
#endif
}


void Conformation_O::buildMoleculeUsingOligomer(int moleculeId, chem::Oligomer_sp oligomer)
{_OF();
  ASSERTF(moleculeId<this->numberOfMolecules(),
          BF("Illegal moleculeId[%d]") % moleculeId );
  ChainNode_sp chainNode =
    this->_FoldTree->buildChainUsingOligomer(moleculeId,oligomer);
  LOG(BF("Built FoldTree--->\n%s")%chainNode->_RootMonomerNode->asString());
  this->_JointTree->buildMoleculeUsingChainNode(moleculeId,chainNode,oligomer);
}


}; /* kinematics */
