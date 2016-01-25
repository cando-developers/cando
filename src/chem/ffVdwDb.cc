       

#define	DEBUG_LEVEL_NONE

#include <cando/chem/ffVdwDb.h>
//#include "core/archiveNode.h"
#include <clasp/core/wrappers.h>



namespace chem {


void	FFVdw_O::fields(core::Record_sp node) {
  this->Base::fields(node);
  node->field(INTERN_(kw,type1),this->_T1);
  node->field(INTERN_(kw,type2),this->_T2);
  node->field(INTERN_(kw,R),this->_R);
  node->field(INTERN_(kw,Eps),this->_Eps);
  node->field(INTERN_(kw,M),this->_M);
  node->field(INTERN_(kw,N),this->_N);
  node->field(INTERN_(kw,Edep),this->_Edep);
};

void FFVdwDb_O::fields(core::Record_sp node)
{
  this->Base::fields(node);
  node->field(INTERN_(kw,vdw),this->_Terms);
};


};
