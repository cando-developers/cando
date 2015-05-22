       
       
//
// (C) 2004 Christian E. Schafmeister
//

#ifndef TrainerHistory_H
#define	TrainerHistory_H

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/common.h"
#include "core/holder.h"
#include "chemPackage.h"


namespace       chem {



SMART(TrainerHistoryEntry);
class TrainerHistoryEntry_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,TrainerHistoryEntry_O,"TrainerHistoryEntry");

public:
void initialize();
//void	archive(core::ArchiveP node);

public:
	string		_Command;
	string		_State;
	int		_StartSeconds;
	int		_StopSeconds;
	string		_Outcome;
	string		_Comment;


public:
	void setCommand(const string& cmd) {this->_Command=cmd;};
	string getCommand() {return this->_Command;};

	void setState(const string& cmd) {this->_State=cmd;};
	string getState() {return this->_State;};

	void start(const string& commandName);
	void finish(const string& stateName );

	void setStartSeconds(int i) {this->_StartSeconds=i;};
	int getStartSeconds() {return this->_StartSeconds;};

	void setStopSeconds(int i) {this->_StopSeconds=i;};
	int getStopSeconds() {return this->_StopSeconds;};

	void setComment(const string& txt) {this->_Comment=txt;};
	string getComment() {return this->_Comment;};

	void setOutcome(const string& txt) {this->_Outcome=txt;};
	string getOutcome() {return this->_Outcome;};
public:
	int	getElapsedSeconds() { return this->_StopSeconds-this->_StartSeconds; };

	DEFAULT_CTOR_DTOR(TrainerHistoryEntry_O);
};



class TrainerHistory_O;
    typedef	gctools::smart_ptr<TrainerHistory_O>	TrainerHistory_sp;
SMART(TrainerHistory);
class TrainerHistory_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,TrainerHistory_O,"TrainerHistory");

public:
//    void	archive(core::ArchiveP node);

private:
    gctools::Vec0<TrainerHistoryEntry_sp>	_Entries;

public:

public:
	void add(TrainerHistoryEntry_sp ent) { this->_Entries.append(ent); };
	void clear() { this->_Entries.clear(); };

    gctools::Vec0<TrainerHistoryEntry_sp>::iterator	begin_Entries()
		{ return this->_Entries.begin(); };

    gctools::Vec0<TrainerHistoryEntry_sp>::iterator	end_Entries()
		{ return this->_Entries.end(); };

	DEFAULT_CTOR_DTOR(TrainerHistory_O);
};




};


TRANSLATE(chem::TrainerHistoryEntry_O);
TRANSLATE(chem::TrainerHistory_O);
#endif
