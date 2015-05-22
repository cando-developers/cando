#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "entity.h"
#include "candoDatabase.h"
#include "core/hashTableEq.h"
#include "core/binder.h"
#include "representedEntityNameSet.h"
#include "constitution.h"
#include "core/wrappers.h"


namespace chem {


//
// Constructor
//

void	Entity_O::initialize()
{
    this->Base::initialize();
}


//
// Copy Constructor
//
    Entity_O::Entity_O(const Entity_O& emr ) : core::T_O(emr)
{
}

//
// Destructor
//

#ifdef XML_ARCHIVE
    void	Entity_O::archiveBase(core::ArchiveP node)
{_G();
    this->Base::archiveBase(node);
}
#endif


Constitution_sp Entity_O::constitution()
{_G();
    Constitution_sp x = _Nil<Constitution_O>();
    return x;
}


    string Entity_O::__repr__() const
    {
	stringstream ss;
	ss << "[" << this->className() << "/name=" << this->getName()->__repr__() << "]";
	return ss.str();
    }



    RepresentativeList_sp Entity_O::minimalRepresentativeList() const
    {_OF();
	SHOUT(BF("Generating minimalRepresentativeList of: %s") % this->__repr__() );
	RepresentativeList_sp expandedList = this->expandedRepresentativeList();
	SHOUT(BF("   after expandedRepresentativeList: %s") % expandedList->__repr__() );
	core::HashTableEq_sp binder = core::HashTableEq_O::create_default();
	for ( int i(0), iEnd(expandedList->length()); i<iEnd; ++i ) {
	    RepresentedEntityNameSet_sp curNameSet = expandedList->elt(i).as<RepresentedEntityNameSet_O>();
	    core::Symbol_sp representor = curNameSet->getRepresentative();
	    RepresentedEntityNameSet_sp minimalRepresentedEntityNameSet;
            core::T_mv bi = binder->gethash(representor);
	    if ( bi.second().notnilp() )
	    {
		minimalRepresentedEntityNameSet = bi.as<RepresentedEntityNameSet_O>(); //binder->indexed_value(bi).as<RepresentedEntityNameSet_O>();
		minimalRepresentedEntityNameSet->mergeEntityNames(curNameSet);
	    } else
	    {
		GC_COPY(RepresentedEntityNameSet_O, temp , *curNameSet); // = RP_Copy<RepresentedEntityNameSet_O>(curNameSet);
		minimalRepresentedEntityNameSet = temp;
		minimalRepresentedEntityNameSet->expandToTerminalEntityNames();
		binder->setf_gethash(representor,minimalRepresentedEntityNameSet);
	    }
	}
	RepresentativeList_sp minimalList = RepresentativeList_O::create();
        binder->maphash( [&minimalList] (core::T_sp key, core::T_sp val) {
                minimalList->vectorPushExtend(val); //binder->indexed_value(val));
            } );
	SHOUT(BF("    resulting minimalRepresentativeList of: %s") % minimalList->__repr__() );
	return minimalList;
    }


void Entity_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<Entity_O>()
	.def("hasConstitution",&Entity_O::hasConstitution)
	.def("constitution",&Entity_O::constitution)
	.def("getName",&Entity_O::getName,"","","",false)
	.def("expandedRepresentativeList",&Entity_O::expandedRepresentativeList)
	.def("minimalRepresentativeList",&Entity_O::minimalRepresentativeList)
    ;
}

    void Entity_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,Entity,"","",_lisp)
	.def("hasConstitution",&Entity_O::hasConstitution)
	.def("constitution",&Entity_O::constitution)
//	.def("getName",&Entity_O::getName,"","","",false)
	.def("expandedRepresentativeList",&Entity_O::expandedRepresentativeList)
	.def("minimalRepresentativeList",&Entity_O::minimalRepresentativeList)
    ;
#endif
}




    EXPOSE_CLASS(chem,Entity_O);
};


