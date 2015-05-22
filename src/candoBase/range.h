#ifndef	Range_H //[
#define Range_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/foundation.h"
#include "core/object.h"
#include "core/numerics.h"
#include "candoBasePackage.fwd.h"

namespace candoBase {


    SMART(Range);
    class Range_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(candoBase,CandoBasePkg,Range_O,"Range");
	DECLARE_INIT();
    public:
//	void	archive(core::ArchiveP node);
	void	initialize();
    private:
	core::LongLongInt	_Begin;
	core::LongLongInt	_End;
	core::LongLongInt	_Step;
    public:
	static Range_sp create(core::Lisp_sp e,core::LongLongInt begin, core::LongLongInt end);
	static Range_sp create3(core::Lisp_sp e,core::LongLongInt begin, core::LongLongInt end, core::LongLongInt step);
    public:

	string __repr__() const;

	core::LongLongInt	begin() { return this->_Begin; };
	core::LongLongInt	end() { return this->_End; };
	core::LongLongInt	step() { return this->_Step; };
	core::Cons_sp asCons();

	Range_O( const Range_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(Range_O);
    };


};
TRANSLATE(candoBase::Range_O);
#endif //]
