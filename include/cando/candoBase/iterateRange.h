#ifndef	IterateRange_H //[
#define IterateRange_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/numbers.fwd.h>
#include <clasp/core/iterator.h>
#include <cando/candoBase/candoBasePackage.h>

namespace candoBase {

SMART(IterateRange );
class IterateRange_O : public core::Iterator_O
{
    LISP_BASE1(core::Iterator_O);
    LISP_CLASS(candoBase,CandoBasePkg,IterateRange_O,"IterateRange");
public:
	void	initialize();
private:
	// instance variables
	int	_Begin;
	int	_Step;
	int	_End;
	core::Fixnum_sp	_Cur;
public:
	static core::T_sp prim_range(core::Function_sp e,core::Cons_sp args, core::Environment_sp environ, core::Lisp_sp lisp);


	static IterateRange_sp create(core::Lisp_sp e,int ibegin, int iend);
	static IterateRange_sp create(core::Lisp_sp e,int ibegin, int iend, int istep);
public:
	void initRange(int b,int e);
	void first();
	void next();
	bool isDone();
	core::T_sp currentObject();

	IterateRange_O( const IterateRange_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(IterateRange_O);
};


};
TRANSLATE(candoBase::IterateRange_O);
#endif //]
