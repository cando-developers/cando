

#include <clasp/core/common.h>
#include <cando/chem/plug.h>
#include <cando/chem/candoDatabase.h>
#include <cando/chem/constitution.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/kinematics/rootAtomInfo.h>

namespace kinematics
{

    void RootAtomInfo::setup(core::Symbol_sp constitutionName,
			     core::Symbol_sp topologyName,
			     chem::Plug_sp plug)
    {
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
            core::Symbol_sp bond1AtomName = inPlug->getB1();
            LOG(BF("InPlug has a Bond1 atom[%s]") % bond1AtomName); 
            chem::CandoDatabase_sp db = chem::getCandoDatabase();
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

