

#include "core/common.h"
#include "chem/plug.h"
#include "chem/candoDatabase.h"
#include "chem/constitution.h"
#include "chem/constitutionAtoms.h"
#include "rootAtomInfo.h"

namespace kinematics
{

    void RootAtomInfo::setup(core::Symbol_sp constitutionName,
			     core::Symbol_sp topologyName,
			     chem::Plug_sp plug)
    {_G();
	LOG(BF("Setting up RootAtomInfo for Constitution[%s] Topology[%s]")
	    % constitutionName->__repr__()
	    % topologyName->__repr__() );
	this->_ConstitutionName = constitutionName;
	this->_TopologyName = topologyName;
	this->_Bond1Id = UndefinedUnsignedInt;
	if ( plug->isAssignableTo<chem::InPlug_O>() )
	{
	    chem::InPlug_sp inPlug = plug.as<chem::InPlug_O>();
	    if ( inPlug->hasB1() )
	    {
		string bond1AtomName = inPlug->getB1();
		LOG(BF("InPlug has a Bond1 atom[%s]") % bond1AtomName); 
		chem::CandoDatabase_sp db = chem::getCandoDatabase(_lisp);
		chem::Constitution_sp constitution = db->getEntity(constitutionName).as<chem::Constitution_O>();
		chem::ConstitutionAtoms_sp constitutionAtoms = constitution->getConstitutionAtoms();
		chem::ConstitutionAtom_sp bond1ConstitutionAtom = constitutionAtoms->atomWithName(bond1AtomName);
		this->_Bond1Id = bond1ConstitutionAtom->atomIndex();
		LOG(BF("The Bond1Id is %d") % this->_Bond1Id );
	    } else
	    {
		LOG(BF("There is no Bond1 atom"));
	    }
	}
    }
};

