/*
    File: msmarts_TypeParser.h
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




#ifndef msmarts_TypeParser_H
#define msmarts_TypeParser_H

#include <clasp/core/common.h>
#include <clasp/core/cons.h>
#include <cando/chem/hold.h>

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
