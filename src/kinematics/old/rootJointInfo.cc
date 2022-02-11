/*
    File: rootJointInfo.cc
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
#include <clasp/core/evaluator.h>
#include <cando/chem/plug.h>
#include <cando/chem/candoDatabase.h>
#include <cando/chem/constitution.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/kinematics/rootJointInfo.h>

namespace kinematics
{

    void RootJointInfo::setup(core::Symbol_sp constitutionName,
                              core::Symbol_sp topologyName,
                              chem::Plug_sp plug)
    {
	LOG(BF("Setting up RootAtomInfo for Constitution[%s] Topology[%s]")
	    % _rep_(constitutionName)
	    % _rep_(topologyName) );
	this->_ConstitutionName = constitutionName;
	this->_TopologyName = topologyName;
	this->_Bond1Id = UndefinedUnsignedInt;
	if ( gc::IsA<chem::InPlug_sp>(plug) ) {
          chem::InPlug_sp inPlug = gc::As_unsafe<chem::InPlug_sp>(plug);
          if ( inPlug->getB1().notnilp() ) {
            core::Symbol_sp bond1AtomName = inPlug->getB1();
            LOG(BF("InPlug has a Bond1 atom[%s]") % bond1AtomName); 
            chem::CandoDatabase_sp db = chem::getCandoDatabase();
            chem::Constitution_sp constitution = gc::As<chem::Constitution_sp>(core::eval::funcall(chem::_sym_constitutionForNameOrPdb,db,constitutionName));
            chem::ConstitutionAtoms_sp constitutionAtoms = constitution->getConstitutionAtoms();
            this->_Bond1Id = constitutionAtoms->index(bond1AtomName);
          } else
          {
            LOG(BF("There is no Bond1 atom"));
          }
	}
    }
};
