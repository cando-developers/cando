#ifndef	O_Loop_H //[
#define O_Loop_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/loop.h>
#include <clasp/core/iterator.h>


#include <cando/chem/chemPackage.h>


namespace chem 
{




SMART(IterateMatter);
class IterateMatter_O : public core::Iterator_O
{
    LISP_BASE1(core::Iterator_O);
    LISP_CLASS(chem,ChemPkg,IterateMatter_O,"IterateMatter");
#if INIT_TO_FACTORIES
 public:
    static IterateMatter_sp make();
#else
    DECLARE_INIT();
#endif
    public:
static void lisp_initGlobals(core::Lisp_sp lisp);
public:
//void	archive(core::ArchiveP node);
	void	initialize();
protected:
	Matter_sp _Top;
	int	_Goal;
	Loop	_Loop;
	bool	_IsDone;
public:
    static IterateMatter_sp create(Matter_sp top, int goal );

public:

		//! Initialize the top of the loop and the goal
	void initTopAndGoal(Matter_sp top, int goal );
	void first();
	void next();
	core::T_sp currentObject();
	bool isDone();
	bool advance();
	Atom_sp getAtom();
	Atom_sp getBondA1();
	Atom_sp getBondA2();
	int	getBondOrder();

	IterateMatter_O( const IterateMatter_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(IterateMatter_O);
};




SMART(IterateAtoms);
class IterateAtoms_O : public IterateMatter_O
{
    LISP_BASE1(IterateMatter_O);
    LISP_CLASS(chem,ChemPkg,IterateAtoms_O,"IterateAtoms");
#if INIT_TO_FACTORIES
 public:
    static IterateAtoms_sp make(Matter_sp matter);
#else
    DECLARE_INIT();
#endif
public:
public:
    static IterateAtoms_sp create(Matter_sp top);
    static core::T_sp createForMatter(Matter_sp matter);
public:

    core::T_sp currentObject();

	IterateAtoms_O( const IterateAtoms_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(IterateAtoms_O);
};




SMART(IterateResidues);
class IterateResidues_O : public IterateMatter_O
{
    LISP_BASE1(IterateMatter_O);
    LISP_CLASS(chem,ChemPkg,IterateResidues_O,"IterateResidues");
public:
private:
public:
    static IterateResidues_sp create(Matter_sp top);
    static core::T_sp createForMatter(Matter_sp matter);

public:

    core::T_sp currentObject();


	DEFAULT_CTOR_DTOR(IterateResidues_O);
};



SMART(IterateBonds);
class IterateBonds_O : public IterateMatter_O
{
    LISP_BASE1(IterateMatter_O);
    LISP_CLASS(chem,ChemPkg,IterateBonds_O,"IterateBonds");
public:
private:
public:
    static IterateBonds_sp create(Matter_sp top);
    static core::T_sp createForMatter(Matter_sp matter);

public:

    core::T_sp currentObject();
	Atom_sp	getAtom1();
	Atom_sp getAtom2();
	int getBondOrder();



	DEFAULT_CTOR_DTOR(IterateBonds_O);
};




};
TRANSLATE(chem::IterateMatter_O);
TRANSLATE(chem::IterateAtoms_O);
TRANSLATE(chem::IterateResidues_O);
TRANSLATE(chem::IterateBonds_O);
#endif //]
