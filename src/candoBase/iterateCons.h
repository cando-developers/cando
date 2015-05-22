#ifndef	IterateCons_H //[
#define IterateCons_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/foundation.h"
#include "core/object.h"
#include "core/iterator.h"
#include "candoBasePackage.fwd.h"

namespace candoBase {


SMART(Cons);


SMART(IterateCons );
    class IterateCons_O : public core::Iterator_O
{
    LISP_BASE1(core::Iterator_O);
    LISP_CLASS(candoBase,CandoBasePkg,IterateCons_O,"IterateCons");
public:
	void	initialize();
private:
	// instance variables
    core::Cons_sp	_Begin;
    core::Cons_sp	_Cur;
public:
    static IterateCons_sp create(core::Cons_sp beg)
    {
	GC_ALLOCATE(IterateCons_O, r );
	r->init(beg);
	return r;
    };
public:

    void init(core::Cons_sp c) { this->_Begin = c; };
	void first();
	void next();
	bool isDone();
    core::T_sp currentObject();



	DEFAULT_CTOR_DTOR(IterateCons_O);
};


};


TRANSLATE(candoBase::IterateCons_O);

#endif //]
