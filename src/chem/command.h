       
       

#ifndef	Command_H //[
#define Command_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>


#include <cando/chem/chemPackage.h>


namespace chem {


/*! Commands encapsulate a bunch of primitive operations
 * and provide the capability of undoing them.
 */


SMART(Command);
SMART(Command);
SMART(Command );
class Command_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,Command_O,"Command");

public:
//    void	archiveBase(core::ArchiveP node);
public:
	virtual	void doIt() {};
	virtual void undoIt() {};

	Command_O( const Command_O& ss ); //!< Copy constructor


	DEFAULT_CTOR_DTOR(Command_O);
};




};
TRANSLATE(chem::Command_O);
#endif //]
