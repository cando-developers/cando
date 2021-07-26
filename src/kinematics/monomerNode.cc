/*
    File: monomerNode.cc
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
#include <clasp/core/lispStream.h>
#include <cando/adapt/symbolSet.h>
#include <cando/kinematics/joint.h>
#include <cando/chem/monomer.h>
#include <cando/chem/coupling.h>
#include <cando/chem/constitution.h>
#include <cando/chem/topology.h>
#include <cando/chem/candoDatabase.h>
#include <cando/kinematics/chainNode.h>
#include <cando/kinematics/monomerNode.h>
#include <clasp/core/wrappers.h>
namespace kinematics
{

// ----------------------------------------------------------------------
//



void MonomerNode_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw,parent),this->_Parent);
  node->field(INTERN_(kw,parent_plug_name),this->_ParentPlugName);
  node->field(INTERN_(kw,Id),this->_Id);
  node->field(INTERN_(kw,Children),this->_Children);
  node->field(INTERN_(kw,StereoisomerName),this->_StereoisomerName);
  node->field(INTERN_(kw,Topology),this->_Topology);
  node->field(INTERN_(kw,ConformationIndex),this->_ConformationIndex);
  node->field(INTERN_(kw,Joints),this->_Joints);
  this->Base::fields(node);
}


MonomerNode_sp MonomerNode_O::create(MonomerId monomerId)
{
  auto  monomerNode  = gctools::GC<MonomerNode_O>::allocate_with_default_constructor();
  monomerNode->_Id = monomerId;
  return monomerNode;
}

CL_LISPIFY_NAME(make-MonomerNode);
MonomerNode_sp MonomerNode_O::make()
{
  IMPLEMENT_ME();
};

void MonomerNode_O::initialize()
{_OF();
  this->Base::initialize();
}


CL_LISPIFY_NAME(stereoisomer-name);
CL_DEFMETHOD core::Symbol_sp	MonomerNode_O::stereoisomerName() const { return this->_StereoisomerName;};

void MonomerNode_O::recursivelyBuildChildren(ChainNode_sp chainNode,
                                             RingClosingMonomerMap ringClosingMonomerMap,
                                             gctools::Nilable<MonomerNode_sp> parent,
                                             gctools::Nilable<chem::DirectionalCoupling_sp> coupling,
                                             chem::Monomer_sp monomer)
{
  LOG(BF("recursivelyBuildChildren on monomer[%s]") % _rep_(monomer->getName()) );
  this->_Parent = parent;
  if ( coupling.notnilp() )
  {
    this->_ParentPlugName = coupling->getTargetMonomerPlugName();
    ASSERTF(chem::DirectionalCoupling_O::isInPlugName(this->_ParentPlugName),
            BF("This has to be an inPlugName[%s] - it isnt")
            % _rep_(this->_ParentPlugName) );
  } else
  {
    this->_ParentPlugName = nil<core::Symbol_O>();
  }
  this->_StereoisomerName = monomer->currentStereoisomerName();
  this->_Topology = monomer->currentTopology();
  this->_ConformationIndex = 0;
  for ( chem::Monomer_O::Couplings::const_iterator it=monomer->_Couplings.begin();
        it!=monomer->_Couplings.end(); it++ )
  {
    chem::Coupling_sp coupling = it->second;
    if ( coupling->isRingClosing() ) continue;
	    /*! Skip InPlugs */
    chem::DirectionalCoupling_sp dirCoupling = coupling.as<chem::DirectionalCoupling_O>();
    if ( dirCoupling->getSourceMonomer() != monomer ) continue;
    chem::Monomer_sp otherMonomer = dirCoupling->getTargetMonomer();
    MonomerNode_sp otherMonomerNode = ChainNode_O::monomerNodeFactory(chainNode,ringClosingMonomerMap,otherMonomer);
    core::Symbol_sp outPlugName = dirCoupling->getSourceMonomerPlugName();
    ASSERTF(chem::DirectionalCoupling_O::isOutPlugName(outPlugName),
            BF("Problem - this[%s] should be an OutPlug name but it isnt")
            %_rep_(outPlugName) );
    this->_Children.set(outPlugName,otherMonomerNode);
    otherMonomerNode->recursivelyBuildChildren(chainNode,
                                               ringClosingMonomerMap,
                                               this->sharedThis<MonomerNode_O>(),
                                               dirCoupling,otherMonomer);
  }
}


chem::Constitution_mv MonomerNode_O::identifyConstitutionAndTopology()
{
  chem::Constitution_sp constitution = this->_Topology->getConstitution();
  return Values(constitution,this->_Topology);
#if 0  
  core::T_sp tconstitution = core::eval::funcall(chem::_sym_constitutionForNameOrPdb,_Nil<core::T_O>(),this->stereoisomerName());
  if (tconstitution.nilp()) {
    SIMPLE_ERROR(BF("Could not find constitution for monomer %s") % _rep_(this->_StereoisomerName));
  }
  chem::Constitution_sp constitution  = gc::As<chem::Constitution_sp>(tconstitution);
  adapt::SymbolSet_sp myPlugNameSet = adapt::SymbolSet_O::create();
  if ( this->_ParentPlugName.notnilp() )
  {
    myPlugNameSet->insert(this->_ParentPlugName);
  }
  for ( adapt::SymbolMap<MonomerNode_O>::iterator it=this->_Children.begin();
        it!=this->_Children.end(); it++ )
  {
    myPlugNameSet->insert(it->first);
  }
  chem::Topology_sp topology; // default is 0x0
  chem::Constitution_O::TopologyMap::iterator it;
  {_BLOCK_TRACEF(BF("Looking for Topology with plugs: %s") % myPlugNameSet->asString() );
    for ( it= constitution->begin_Topologies();
          it!=constitution->end_Topologies(); it++ ) {
      if ( (it->second)->hasMatchingPlugsWithMates(myPlugNameSet))
      {
        topology = it->second;
        break;
      }
    }
  }
  if (!topology) {
    SIMPLE_ERROR(BF("No topology could be found for monomer %s with plugs %s")
                 % _rep_(this->_StereoisomerName) % myPlugNameSet->asString() );
  }
  return Values(constitution, topology);
#endif
}


void MonomerNode_O::describeRecursivelyIntoStringStream(const string& prefix, stringstream& output) const
{_OF();
  LOG(BF("Describing %s[%s]") % this->className() % _rep_(this->_StereoisomerName) );
  output << prefix;
  if ( this->_ParentPlugName.notnilp() )
  {
    output << _rep_(this->_ParentPlugName) << " ";
  }
  output << this->className()
         << "["
         << _rep_(this->_StereoisomerName)
         << "]" << std::endl;
  for ( adapt::SymbolMap<MonomerNode_O>::const_iterator it=this->_Children.begin();
        it!=this->_Children.end(); it++ )
  {
    it->second->describeRecursivelyIntoStringStream(prefix+"  ",output);
  }
}

string MonomerNode_O::asString() const
{_OF();
  stringstream out;
  this->describeRecursivelyIntoStringStream("",out);
  return out.str();
}

string MonomerNode_O::__repr__() const {
  stringstream ss;
  ss << "#<";
  ss << this->className();
  ss << " " << _rep_(this->_StereoisomerName) << ">";
  return ss.str();
}

void MonomerNode_O::addJoint(size_t index, Joint_sp joint)
{
//  core::write_bf_stream(BF("%s:%d:%s  joint.id = %s  index -> %lu  joint = %s\n") % __FILE__ % __LINE__ % __FUNCTION__ % joint->atomId().asString() % index % core::_rep_(joint));
  if (this->_Joints.size() <= index ) {
    this->_Joints.resize(index+1);
  }
  this->_Joints[index] = joint;
}

core::List_sp MonomerNode_O::children() const
{
  ql::list l;
  for ( auto i : this->_Children ) {
    l << i.second;
  }
  return l.cons();
}

CL_DEFMETHOD MonomerId MonomerNode_O::monomerId() const
{
  return this->_Id;
}

CL_DOCSTRING("Walk the joints of this monomer-node. Provide a callback that takes two arguments, the joint atom-id index and the joint.");
CL_DEFMETHOD void MonomerNode_O::walkJoints(core::Function_sp callback) {
  for ( int i=0, iEnd(this->_Joints.size()); i<iEnd; ++i ) {
    core::eval::funcall(callback,core::make_fixnum(i),this->_Joints[i]);
  }
}

CL_DEFMETHOD size_t MonomerNode_O::conformationIndex() const
{
  return this->_ConformationIndex;
};
CL_DEFMETHOD void MonomerNode_O::setConformationIndex(size_t i) {
  this->_ConformationIndex = i;
};

CL_DEFMETHOD core::T_sp MonomerNode_O::parent() const {
  if (this->_Parent.boundp()) return this->_Parent;
  return nil<core::T_O>();
}

CL_DEFMETHOD core::Symbol_sp MonomerNode_O::parentPlugName() const {
  return this->_ParentPlugName;
}

CL_DEFMETHOD Joint_sp MonomerNode_O::jointAt(size_t i) const {
  if (i< this->_Joints.size()) {
    return this->_Joints[i];
  }
  SIMPLE_ERROR(BF("Joint index %d must be less than %d") % i % this->_Joints.size());
}

    
}; /* kinematics */
