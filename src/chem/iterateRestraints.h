
#ifndef	IterateRestraints_H //[
#define IterateRestraints_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/common.h"
#include "core/iterator.h"
#include "iterateMatter.h"



#include "chemPackage.h"


namespace chem {


    SMART(IterateRestraints );
    class IterateRestraints_O : public core::Iterator_O
    {
	LISP_BASE1(core::Iterator_O);
	LISP_CLASS(chem,ChemPkg,IterateRestraints_O,"IterateRestraints");
public:
	void	initialize();
private:
	Matter_sp			_Top;
	core::Cons_sp				_Restraints;
	core::Cons_sp				_CurRestraint;
private:
    void _accumulateRestraints(Matter_sp m);

public:
    static IterateRestraints_sp create(core::Lisp_sp e,Matter_sp matter);
public:
	void initTop(Matter_sp matter);
	void first();
	void next();
	bool isDone();
	core::T_sp currentObject();
	
//	IterateRestraints_O( const IterateRestraints_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(IterateRestraints_O);
    };



};
TRANSLATE(chem::IterateRestraints_O);
#endif //]
