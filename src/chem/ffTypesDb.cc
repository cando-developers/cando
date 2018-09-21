/*
    File: ffTypesDb.cc
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
//
// (C) 2004 Christian E. Schafmeister
//

/*
 *	ffTypes.cc
 *
 *	Maintains databases and structures to store types and type assignement
 *	rules.
 */
#include <clasp/core/foundation.h>
#include <clasp/core/hashTable.h>
#include <clasp/core/hashTableEq.h>
#include <clasp/core/record.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/loop.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/wrappers.h>
#include <clasp/core/bformat.h>





namespace chem {






void	FFTypesDb_O::initialize()
{
    this->Base::initialize();
    this->_TypeAssignmentRules.clear();
}


SYMBOL_EXPORT_SC_(ChemPkg,assignType);

CL_LISPIFY_NAME("assignType");
CL_LAMBDA((types-db !) atom &key verbose)
CL_DEFMETHOD core::Symbol_sp FFTypesDb_O::assignType(chem::Atom_sp atom, bool verbose) {
  LOG(BF("Got atom") );
  LOG(BF("atom name: %s") % atom->getName().c_str() );
  LOG(BF("Assigning type for atom: %s") % atom->description().c_str()  );
  {_BLOCK_TRACE("Testing every type rule");
    for ( auto it=this->_TypeAssignmentRules.begin();
          it!=this->_TypeAssignmentRules.end(); it++ ) 
    {_BLOCK_TRACEF(BF("Testing rule code(%s)") % (*it)->getCode().c_str() );
//		LOG(BF("as xml: %s") % ((*it)->asXmlString().c_str() ) );
      ChemInfo_sp cheminfo = gc::As<ChemInfo_sp>(*it);
      if ( cheminfo->matches(atom)) {
        AntechamberRoot_sp antechamberRoot = gc::As<AntechamberRoot_sp>(cheminfo->_Root);
        LOG(BF("Rule MATCH!!!") );
        if (verbose) write_bf_stream(BF("Matched %s\n") % _rep_(cheminfo));
        return antechamberRoot->_AssignType;
      } else {
        if (verbose) write_bf_stream(BF("Did not match %s\n") % _rep_(cheminfo));
      }
      LOG(BF("Rule does not match, keep going") );
    }
  }
  return _Nil<core::Symbol_O>();
}


CL_LISPIFY_NAME("assignTypes");
CL_DEFMETHOD void    FFTypesDb_O::assignTypes(chem::Matter_sp matter)
{ 
  chem::Loop    				lAtoms;
  chem::Atom_sp  				atom;
  chem::Matter_sp				c;
  gctools::Vec0<chem::ChemInfo_sp>::iterator  it;
  chem::Molecule_sp                     	mol;
  chem::Residue_sp				res;
  string 		                        name;
  if (this->_TypeAssignmentRules.size()==0) return;
  lAtoms.loopTopGoal(matter,ATOMS);
  LOG(BF("defined loop") );
  while ( lAtoms.advanceLoopAndProcess() ) {
    LOG(BF("Getting container") );
    c = lAtoms.getMatter();
//        rawGet = c.get();
//        castGet = dynamic_cast<chem::Atom_O*>(c.get());
//        LOG(BF("rawGet = %X") % rawGet );
//        LOG(BF("castGet = %X") % castGet );
//        LOG(BF("getting first atom in loop") );
    atom = (c).as<chem::Atom_O>();
    core::Symbol_sp type = this->assignType(atom,false);
    atom->setType(type);
  }
}


CL_DOCSTRING("Return the atom types as a hash-table with the atoms as keys and types as values");
CL_DEFMETHOD core::HashTableEq_sp    FFTypesDb_O::atomTypes(chem::Matter_sp matter)
{ 
  chem::Loop    				lAtoms;
  chem::Atom_sp  				atom;
  core::HashTableEq_sp atomTypes = core::HashTableEq_O::create_default();
  if (this->_TypeAssignmentRules.size()==0)  return atomTypes;
  lAtoms.loopTopGoal(matter,ATOMS);
  LOG(BF("defined loop") );
  while ( lAtoms.advanceLoopAndProcess() ) {
    LOG(BF("Getting container") );
    atom = lAtoms.getAtom();
    core::Symbol_sp type = this->assignType(atom,false);
    atomTypes->setf_gethash(atom,type);
  }
  return atomTypes;
}


void FFTypesDb_O::fields(core::Record_sp node)
{	
  node->field(INTERN_(kw,type_rules), this->_TypeAssignmentRules);
  this->FFBaseDb_O::fields(node);
}



CL_LISPIFY_NAME("FFTypes-getRule");
CL_DEFMETHOD chem::ChemInfo_sp	FFTypesDb_O::getRule(uint index)
{_OF();
  if ( index < this->_TypeAssignmentRules.size() ) {
    return this->_TypeAssignmentRules[index];
  }
  SIMPLE_ERROR(BF("Rule index is out of bounds"));
}


CL_DEFUN FFTypesDb_sp chem__make_fftypes_db()
{
  return FFTypesDb_O::create();
}










};
