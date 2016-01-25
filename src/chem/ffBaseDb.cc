#define	DEBUG_LEVEL_NONE


//
// (C) 2004 Christian E. Schafmeister
//


/*!
 *	ffBase.cc
 *
 *	Maintain the base class for parameters
 */
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/forceField.h>
//#include "core/archiveNode.h"
#include <clasp/core/symbolTable.h>
#include <clasp/core/wrappers.h>

SYMBOL_EXPORT_SC_(ChemKwPkg,parameterized);
SYMBOL_EXPORT_SC_(ChemKwPkg,estimated);
SYMBOL_EXPORT_SC_(ChemKwPkg,rough);

namespace chem {

CL_DEFMETHOD string	FFParameter_O::levelDescription()
{
  switch ( this->_Level ) {
  case parameterized:
      return "parameterized ";
  case estimated:
      return "estimated ";
  case rough:
      return "guessed ";
  case unknown:
      return "unknown ";
  }
  return "unknownParameterLevel ";
}


    core::NullTerminatedEnumAssociation paramLevelEnum[] = {
	{ "parameterized", parameterized },
	{ "estimated", estimated },
	{ "rough", rough },
	{ "", -1 }
};

void FFParameter_O::fields(core::Record_sp node)
{
  //this->Base::fields(node); // T_O
  node->field( INTERN_(kw,level), this->_Level);
}


void	FFBaseDb_O::fields(core::Record_sp node)
{
  //this->Base::fields(node); // T_O
  node->field( INTERN_(kw,forceField), this->_ForceField );
}


void	FFBaseDb_O::setForceField(ForceField_sp ff)
{
    ASSERTNOTNULL(ff);
    this->_ForceField = ff;
};















};

