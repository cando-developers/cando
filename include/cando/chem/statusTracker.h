       
       

#ifndef	StatusTracker_H //[
#define StatusTracker_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>


#include <cando/chem/chemPackage.h>


namespace chem {


SMART(StatusTracker);

SMART(StatusTracker);
SMART(StatusTracker );
class StatusTracker_O : public core::CxxObject_O
{
    LISP_BASE1(core::CxxObject_O);
    LISP_CLASS(chem,ChemPkg,StatusTracker_O,"StatusTracker");

public:
	void initialize();
public:
//	void	archive(core::ArchiveP node);
private:
	stringstream	_MessageStream;
	bool		_HasError;
	stringstream	_ErrorStream;
public:

	void	addMessage(const string& msg);
	void	addError(const string& err);
	bool	hasError() { return this->_HasError; };
	void	reset();

	string	getStatus();

	StatusTracker_sp copy();


	StatusTracker_O( const StatusTracker_O& ss ); //!< Copy constructor


	DEFAULT_CTOR_DTOR(StatusTracker_O);
};



};
TRANSLATE(chem::StatusTracker_O);
#endif //]
