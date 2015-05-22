       
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	moe.h
 *
 *	Functions that read and write MOE files.
 *
 */

#ifndef	MOL2_H
#define MOL2_H

#include "core/common.h"


#include <iostream>
#include <fstream>
#include <string>
#include "matter.h"
#include "aggregate.h"
#include "molecule.h"
#include "residue.h"
#include "atom.h"


#include "chemPackage.h"


namespace chem {

void	_calculateElementAndHybridization(Atom_sp a);

void	mol2ReadAggregateFromFileName( Aggregate_sp	a,
                                       const string&	sFileName);

//Aggregate_sp	mol2ReadAggregate( string fileName );


    void	mol2ReadAggregateFromStream( Aggregate_sp 	agg, core::T_sp in );

string	mol2AtomType(Atom_sp a);


    void	mol2WriteAggregateStream( Aggregate_sp 	agg, core::T_sp out);
string	mol2WriteAggregateReturnAsString( Aggregate_sp agg );
void	mol2WriteAggregateToFileName( Aggregate_sp a, const string& sFileName );

void	mol2WriteDumbAggregateToFileName( Aggregate_O& a, const string& sFileName );

void	mol2WriteMatterToFileName( Matter_sp matter, const string& fileName );


};


#endif
