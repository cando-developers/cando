/*
    File: trainerArchive.cc
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
       
#define	DEBUG_LEVEL_FULL

#include "trainerArchive.h"
//#include "core/archive.h"
//#include "core/xmlSaveArchive.h"
//#include "core/xmlLoadArchive.h"
#include "trainer.h"
#include "core/wrappers.h"



namespace chem {





void	TrainerArchive_O::initialize()
{
    this->Base::initialize();
    this->_CandoDatabase = _Nil<CandoDatabase_O>();
    this->_Trainers.clear();
}

#ifdef XML_ARCHIVE
    void	TrainerArchive_O::archive(core::ArchiveP node)
{
    node->attribute("CandoDatabase",this->_CandoDatabase);
    node->attributeStringMap("Trainers",this->_Trainers);
}
#endif


#if 0
void	TrainerArchive_O::setCandoDatabase(CandoDatabase_sp hold)
{
    ASSERTNOTNULL(hold);
    this->_CandoDatabase = hold;
}
#endif

void	TrainerArchive_O::addTrainer(Trainer_sp job)
{
string		key;
    key = job->getHeader()->getContextKey();
    this->_Trainers.set(key,job);
}


bool	TrainerArchive_O::hasTrainerWithKey(const string& key)
{
    return this->_Trainers.contains(key);
}

Trainer_sp	TrainerArchive_O::getTrainerWithKey(const string& key)
{
    return this->_Trainers.get(key);
}



};
