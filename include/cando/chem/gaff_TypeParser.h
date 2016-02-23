/*
    File: gaff_TypeParser.h
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
