#ifndef	IterateCons_H //[
#define IterateCons_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/iterator.h>
#include <cando/adapt/adaptPackage.fwd.h>

namespace adapt {


SMART(Cons);


SMART(IterateCons );
    class IterateCons_O : public core::Iterator_O
{
    LISP_CLASS(adapt,AdaptPkg,IterateCons_O,"IterateCons",core::Iterator_O);
public:
	void	initialize();
private:
	// instance variables
    core::List_sp	_Begin;
    core::List_sp	_Cur;
public:
    static IterateCons_sp create(core::List_sp beg)
    {
	GC_ALLOCATE(IterateCons_O, r );
	r->init(beg);
	return r;
    };
public:

CL_NAME("init");
CL_DEFMETHOD     void init(core::List_sp c) { this->_Begin = c; };
	void first();
	void next();
	bool isDone();
    core::T_sp currentObject();



	DEFAULT_CTOR_DTOR(IterateCons_O);
};


};


TRANSLATE(adapt::IterateCons_O);

#endif //]
