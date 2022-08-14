/*
    File: mol2.h
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

#ifndef	MOL2_H
#define MOL2_H

#include <clasp/core/common.h>


#include <iostream>
#include <fstream>
#include <string>
#include <cando/chem/matter.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/molecule.h>
#include <cando/chem/residue.h>
#include <cando/chem/atom.h>

#include <cando/chem/typeAssignmentRules.fwd.h>
#include <cando/chem/chemPackage.h>

namespace chem {

extern chem::TypeAssignmentRules_sp	sybylRules;

class	Mol2File
{
public:
  core::T_sp	fIn;
  bool _eof;
  std::stringstream	mLine;
public:
  bool eof() { return this->_eof;};
  void	advanceLine();
  stringstream&	line() { return this->mLine; };
  std::queue<string>	splitLine();
  bool	hasDataLine();
  void	openFileName(core::T_sp pn);
  Mol2File() : _eof(false) {};
  Mol2File(core::T_sp stream) : fIn(stream), _eof(false) {};
  ~Mol2File();
};

void	_calculateElementAndHybridization(Atom_sp a);

core::T_sp mol2Read( Mol2File& fin);

//Aggregate_sp	mol2ReadAggregate( string fileName );
string	mol2AtomType(Atom_sp a);


 void	mol2WriteAggregateStream( Aggregate_sp agg, std::ostream &out, bool useSybylTypes );
 void	mol2WriteAggregateToFileName( Aggregate_sp a, core::T_sp sFileName, bool useSybylTypes );

 void	mol2WriteDumbAggregateToFileName( Aggregate_O& a, core::T_sp sFileName );

 void	mol2WriteMatterToFileName( Matter_sp matter, core::T_sp fileName, bool useSybylTypes );


 void initialize_mol2_type_rules();

};


#endif
