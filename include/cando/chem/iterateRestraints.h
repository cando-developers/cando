
#ifndef	IterateRestraints_H //[
#define IterateRestraints_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/iterator.h>
#include <cando/chem/iterateMatter.h>



#include <cando/chem/chemPackage.h>


namespace chem {


    SMART(IterateRestraints );
    class IterateRestraints_O : public core::Iterator_O
    {
	LISP_CLASS(chem,ChemPkg,IterateRestraints_O,"IterateRestraints",core::Iterator_O);
public:
	void	initialize();
private:
	Matter_sp			_Top;
	core::List_sp				_Restraints;
	core::List_sp				_CurRestraint;
private:
    void _accumulateRestraints(Matter_sp m);

public:
    static IterateRestraints_sp create(Matter_sp matter);
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
