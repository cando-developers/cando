#ifndef	Range_H //[
#define Range_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/numerics.h>
#include <cando/geom/geomPackage.fwd.h>

namespace geom {


    SMART(Range);
    class Range_O : public core::CxxObject_O
    {
	LISP_CLASS(geom,GeomPkg,Range_O,"Range",core::CxxObject_O);
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
TRANSLATE(geom::Range_O);
#endif //]
