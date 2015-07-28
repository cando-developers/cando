#ifndef	Container_H //[
#define Container_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <cando/geom/geomPackage.fwd.h>

namespace geom {

/*! Virtual class for containers
 */
    SMART(Container );
    class Container_O : public core::CxxObject_O
    {
	LISP_BASE1(core::CxxObject_O);
	LISP_CLASS(geom,GeomPkg,Container_O,"Container");
public:
public:

    virtual bool isEmpty() { return this->size()==0; };
    virtual bool notEmpty() { return this->size()!=0;};

    virtual int size() {_OF(); SUBCLASS_MUST_IMPLEMENT();};

//	Container_O( const Container_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(Container_O);
    };



};
TRANSLATE(geom::Container_O);
#endif //]
