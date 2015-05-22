       
       
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
