#define	DEBUG_LEVEL_FULL       
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
#include <clasp/core/record.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/loop.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/wrappers.h>





namespace chem {






void	FFTypesDb_O::initialize()
{
    this->Base::initialize();
    this->_TypeAssignmentRules.clear();
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
        atom = downcast<chem::Atom_O>(c);
        LOG(BF("Got atom") );
        LOG(BF("atom name: %s") % atom->getName().c_str() );
        LOG(BF("Assigning type for atom: %s") % atom->description().c_str()  );
	{_BLOCK_TRACE("Testing every type rule");
	    for ( it=this->_TypeAssignmentRules.begin();
		    it!=this->_TypeAssignmentRules.end(); it++ ) 
	    {_BLOCK_TRACEF(BF("Testing rule code(%s)") % (*it)->getCode().c_str() );
//		LOG(BF("as xml: %s") % ((*it)->asXmlString().c_str() ) );
		if ( (*it)->matches(atom) ) 
		{
		    LOG(BF("Rule MATCH!!!") );
		    break;
		}
		LOG(BF("Rule does not match, keep going") );
	    }
	} 
    }
}


void FFTypesDb_O::fields(core::Record_sp node)
    {	
      this->FFBaseDb_O::fields(node);
      node->field(INTERN_(kw,type_rules), this->_TypeAssignmentRules);
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




