#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "core/environment.h"
#include "core/symbolSet.h"
#include "chem/monomer.h"
#include "chem/coupling.h"
#include "chem/constitution.h"
#include "chem/topology.h"
#include "chem/candoDatabase.h"
#include "chainNode.h"
#include "monomerNode.h"
#include "core/wrappers.h"
namespace kinematics
{

// ----------------------------------------------------------------------
//

    EXPOSE_CLASS(kinematics,MonomerNode_O);

    void MonomerNode_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<MonomerNode_O>()
//	.def_raw("__init__",&MonomerNode_O::__init__,"(self)")
	    .def("monomerName",&MonomerNode_O::monomerName)
	    .def("monomer-node-topology",&MonomerNode_O::topology)
	;
    }

    void MonomerNode_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(KinPkg,MonomerNode,"","",_lisp)
//	.def_raw("__init__",&"(self)")
	    .def("monomerName",&MonomerNode_O::monomerName)
	    .def("monomer-node-topology",&MonomerNode_O::topology)
	;
#endif
    }


    MonomerNode_sp MonomerNode_O::create(int monomerId)
    {_G();
	GC_ALLOCATE(MonomerNode_O, monomerNode );
	monomerNode->_MonomerId = monomerId;
	return monomerNode;
    }

#if INIT_TO_FACTORIES

#define ARGS_MonomerNode_O_make "()"
#define DECL_MonomerNode_O_make ""
#define DOCS_MonomerNode_O_make "make MonomerNode"
  MonomerNode_sp MonomerNode_O::make()
    {_G();
      IMPLEMENT_ME();
    };

#else
    core::T_sp MonomerNode_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_G();
//      this->Base::__init__(exec,args,env,lisp);
//      arg = translate::from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	return _Nil<core::T_O>();
    }
#endif

#ifdef XML_ARCHIVE
    void MonomerNode_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
    }
#endif


    void MonomerNode_O::initialize()
    {_OF();
        this->Base::initialize();
    }



    void MonomerNode_O::recursivelyBuildChildren(ChainNode_sp chainNode,
						 RingClosingMonomerMap& ringClosingMonomerMap,
						 MonomerNode_sp parent,
						 chem::DirectionalCoupling_sp coupling,
						 chem::Monomer_sp monomer)
    {_G();
	LOG(BF("recursivelyBuildChildren on monomer[%s]") % monomer->getName()->__repr__() );
	this->_Parent = parent;
	if ( coupling.notnilp() )
	{
	    this->_ParentPlugName = coupling->getOutMonomerPlugName();
	    ASSERTF(chem::DirectionalCoupling_O::isInPlugName(this->_ParentPlugName),
		    BF("This has to be an inPlugName[%s] - it isnt")
		    % this->_ParentPlugName->__repr__() );
	} else
	{
	    this->_ParentPlugName = core::Symbol_O::_nil;
	}
	this->_MonomerName = monomer->getName();
	this->_Topology = monomer->getTopology();
	this->_ConformationIndex = 0;
	for ( chem::Monomer_O::Couplings::const_iterator it=monomer->begin_WeakCouplings();
	      it!=monomer->end_WeakCouplings(); it++ )
	{
	    chem::Coupling_sp coupling = it->second.lock();
	    if ( coupling->isRingClosing() ) continue;
	    /*! Skip InPlugs */
	    chem::DirectionalCoupling_sp dirCoupling = coupling.as<chem::DirectionalCoupling_O>();
	    if ( dirCoupling->getInMonomer() != monomer ) continue;
	    chem::Monomer_sp otherMonomer = dirCoupling->getOutMonomer();
	    MonomerNode_sp otherMonomerNode = ChainNode_O::monomerNodeFactory(chainNode.get(),ringClosingMonomerMap,otherMonomer,_lisp);
	    core::Symbol_sp outPlugName = dirCoupling->getInMonomerPlugName();
	    ASSERTF(chem::DirectionalCoupling_O::isOutPlugName(outPlugName),
		    BF("Problem - this[%s] should be an OutPlug name but it isnt")
		    %outPlugName->__repr__() );
	    this->_Children.set(outPlugName,otherMonomerNode);
	    otherMonomerNode->recursivelyBuildChildren(chainNode,
						       ringClosingMonomerMap,
						       this->sharedThis<MonomerNode_O>(),
						       dirCoupling,otherMonomer);
	}
    }


    core::Cons_sp MonomerNode_O::identifyConstitutionAndTopology(chem::CandoDatabase_sp db)
    {_OF();
	chem::Constitution_sp constitution = db->constitutionForNameOrPdb(this->_MonomerName);
	if ( constitution.nilp() )
	{
	    SIMPLE_ERROR(BF("Could not find Constitution for monomer[%s]")
			       % this->_MonomerName->__repr__() );
	}
	core::SymbolSet_sp myPlugNameSet = core::SymbolSet_O::create();
	if ( this->_ParentPlugName.notnilp() )
	{
	    myPlugNameSet->insert(this->_ParentPlugName);
	}
	for ( core::SymbolMap<MonomerNode_O>::iterator it=this->_Children.begin();
	      it!=this->_Children.end(); it++ )
	{
	    myPlugNameSet->insert(it->first);
	}
	chem::Topology_sp topology = _Nil<chem::Topology_O>();
	chem::Constitution_O::TopologyMap::iterator it;
	{_BLOCK_TRACEF(BF("Looking for Topology with plugs: %s") % myPlugNameSet->asString() );
	    for ( it= constitution->begin_Topologies();
		  it!=constitution->end_Topologies(); it++ )
	    {
		if ( (it->second)->hasMatchingPlugsWithMates(myPlugNameSet))
		{
		    topology = it->second;
		    break;
		}
	    }
	}
	if ( topology.nilp() )
	{
	    SIMPLE_ERROR(BF("No topology could be found for monomer[%s] with plugs[%s]")
			       % this->_MonomerName->__repr__() % myPlugNameSet->asString() );
	}
	return core::Cons_O::createList(constitution,topology);
    }


    void MonomerNode_O::describeRecursivelyIntoStringStream(const string& prefix, stringstream& output) const
    {_OF();
	LOG(BF("Describing %s[%s]") % this->className() % this->_MonomerName->__repr__() );
	output << prefix;
	if ( this->_ParentPlugName.notnilp() )
	{
	    output << this->_ParentPlugName->__repr__() << " ";
	}
	output << this->className()
	       << "["
	       << this->_MonomerName->__repr__()
	       << "]" << endl;
	for ( core::SymbolMap<MonomerNode_O>::const_iterator it=this->_Children.begin();
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


}; /* kinematics */
