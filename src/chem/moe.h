       
       
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
#include "matter.h"
#include "molecule.h"
#include "residue.h"
#include "atom.h"

#include "chemPackage.h"


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
