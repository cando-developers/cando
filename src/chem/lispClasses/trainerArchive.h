/*
    File: trainerArchive.h
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

#ifndef TrainerArchive_H
#define	TrainerArchive_H

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/common.h"
#include "core/holder.h"
#include "candoDatabase.h"
#include "monomerContext.h"
#include "conformationExplorer.h"


//#include "jobHistory.h"


namespace       chem
{



    SMART(Trainer);
SMART(TrainerArchive);
class TrainerArchive_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,TrainerArchive_O,"TrainerArchive");

public:
	void	initialize();
public:
//	void	archive(core::ArchiveP node);
public:
	CandoDatabase_sp	_CandoDatabase;
	core::StringMap<Trainer_O>	_Trainers;

public:
	typedef core::StringMap<Trainer_O>::iterator	jobIterator;
public:

	void	addTrainer(Trainer_sp job);

	bool	hasTrainerWithKey(const string& key);
	Trainer_sp	getTrainerWithKey(const string& key);


	jobIterator	begin_Trainers() { return this->_Trainers.begin(); }
	jobIterator	end_Trainers() { return this->_Trainers.end(); }

	DEFAULT_CTOR_DTOR(TrainerArchive_O);
};




};


TRANSLATE(chem::TrainerArchive_O);
#endif
