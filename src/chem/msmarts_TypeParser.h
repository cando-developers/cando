



#ifndef msmarts_TypeParser_H
#define msmarts_TypeParser_H

#include "core/common.h"
#include "core/cons.h"
#include "hold.h"

typedef union
{
    const char*				eval;
    char					cval[1024];
    int 					ival;
    chem::BondEnum             		benum;
    chem::Hold<chem::Root_O>*			root;
    chem::Hold<chem::SmartsRoot_O>*		smartsRoot;
    chem::Hold<chem::BondListMatchNode_O>*	bondListMatchNode;
    chem::Hold<chem::AtomOrBondMatchNode_O>*	atomOrBondMatchNode;
//		chem::Hold<chem::AtomTest_O>*		atomTest;
    chem::Hold<chem::BondTest_O>*		bondTest;
    chem::Hold<chem::Logical_O>*			logical;
} msmarts_STypeParser;

#define YYSTYPE msmarts_STypeParser


#endif
