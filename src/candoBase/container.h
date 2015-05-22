#ifndef	Container_H //[
#define Container_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/foundation.h"
#include "core/object.h"
#include "candoBasePackage.fwd.h"

namespace candoBase {

/*! Virtual class for containers
 */
    SMART(Container );
    class Container_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(candoBase,CandoBasePkg,Container_O,"Container");
public:
public:

    virtual bool isEmpty() { return this->size()==0; };
    virtual bool notEmpty() { return this->size()!=0;};

    virtual int size() {_OF(); SUBCLASS_MUST_IMPLEMENT();};

//	Container_O( const Container_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(Container_O);
    };



};
TRANSLATE(candoBase::Container_O);
#endif //]
