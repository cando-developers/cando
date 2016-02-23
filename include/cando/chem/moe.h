/*
    File: moe.h
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
       
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	moe.h
 *
 *	Functions that read and write MOE files.
 *
 */

#ifndef	MOE_H
#define MOE_H

#include <iostream>
#include <fstream>
#include <string>
#include <cando/chem/matter.h>
#include <cando/chem/molecule.h>
#include <cando/chem/residue.h>
#include <cando/chem/atom.h>

#include <cando/chem/chemPackage.h>


namespace chem {

class	Aggregate_O;
    typedef	gctools::smart_ptr<Aggregate_O>	Aggregate_sp;


void	moeReadAggregateFromFileName(  Aggregate_sp	a,
				const char* 	sFileName );

inline void	moeReadAggregateFromStringFileName(  Aggregate_sp	a,
						const string&	sFileName ) 
{
	moeReadAggregateFromFileName(a,sFileName.c_str());
};


Aggregate_sp	moeReadAggregate( const string& fileName);


    void	moeReadAggregateFromStream( Aggregate_sp 	agg, core::T_sp stream);

void	moeReadAggregateWithAtomTypesFromFileName(  Aggregate_sp	a,
                                                    const char* 	sFileName );

    Aggregate_sp	moeReadAggregateWithAtomTypes(  const string& nm );





    void	moeWriteAggregateStream( Aggregate_sp 	agg, core::T_sp stream);
string	moeWriteAggregateReturnAsString( Aggregate_sp 	agg );
void	moeWriteAggregateFileName( Aggregate_sp a, const char* sFileName );
void	moeWrite( Aggregate_sp a, const string& sFileName );
void	moeRead( Aggregate_sp a, const string& sFileName );

};



#endif
