/*
    File: rootAtomInfo.cc
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
	    % _rep_(constitutionName)
	    % _rep_(topologyName) );
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
