



#ifndef gaff_TypeParser_H
#define gaff_TypeParser_H

#include <clasp/core/foundation.h>
#include <clasp/core/lisp.h>
#include <clasp/core/object.h>
#include <clasp/core/cons.h>


typedef union
{
  char                            carr[255];
  int 				ival;
  chem::BondEnum              benum;
  chem::Hold<chem::AntechamberRoot_O>	*antechamberRoot;
  chem::Hold<chem::BondListMatchNode_O>		*bondListMatchNode;
  chem::Hold<chem::ResidueList_O>	*residueList;
  chem::Hold<chem::Logical_O>		*logical;
  chem::Hold<chem::AtomTest_O>		*atomTest;
  chem::Hold<chem::AfterMatchBondTest_O>	*afterMatchBondTest;
  chem::Hold<chem::Chain_O>		*chain;
  chem::Hold<chem::AntechamberBondTest_O> *antechamberBondTest;
} gaff_STypeParser;

#define YYSTYPE gaff_STypeParser


#endif
