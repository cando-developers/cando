       
       
#ifndef	candoScript_H //[
#define candoScript_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/common.h"


#include "chemPackage.h"


namespace chem {

extern	void	setupCandoPrimitives(core::Lisp_sp);
extern	void	setupPythonPrimitives(core::Lisp_sp);

SMART(MultiMonomer);
SMART(CandoDatabase);
SMART(OligomerPart_Base);
/*! A simple object to store and create a monomer within an oligomer
 */
SMART(OligomerPart_Base);
class OligomerPart_Base_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,OligomerPart_Base_O,"OligomerPart-Base");

public:
//    void archiveBase(core::ArchiveP node);
public:
	virtual MultiMonomer_sp createMonomer(CandoDatabase_sp bdb) {_OF();SUBCLASS_MUST_IMPLEMENT();};
	DEFAULT_CTOR_DTOR(OligomerPart_Base_O);
};


SMART(OligomerPart_Monomer);
/*! A simple object to store and create a monomer within an oligomer
 */
SMART(OligomerPart_Monomer);
class OligomerPart_Monomer_O : public OligomerPart_Base_O
{
    LISP_BASE1(OligomerPart_Base_O);
    LISP_CLASS(chem,ChemPkg,OligomerPart_Monomer_O,"OligomerPart-Monomer");

public:
//    void archive(core::ArchiveP node);
public:	// Leave everything exposed to simplify getting values in and out
	core::Symbol_sp _MonomerId;
	core::Symbol_sp	_GroupName;
	core::Cons_sp	_MonomerAliases;
	string	_Comment;
public:
public:
    virtual MultiMonomer_sp createMonomer(CandoDatabase_sp bdb);
    DEFAULT_CTOR_DTOR(OligomerPart_Monomer_O);
};


SMART(OligomerPart_Link);
/*! A simple object to store and create a monomer within an oligomer
 */
SMART(OligomerPart_Link);
class OligomerPart_Link_O : public OligomerPart_Base_O
{
    LISP_BASE1(OligomerPart_Base_O);
    LISP_CLASS(chem,ChemPkg,OligomerPart_Link_O,"OligomerPart-Link");

public:
//    void archive(core::ArchiveP node);
public:	// Leave everything exposed to simplify getting values in and out
	core::Symbol_sp		_Monomer1Id;
	core::Symbol_sp		_Coupling;
	OligomerPart_Monomer_sp	_Monomer2;
public:
public:
    virtual MultiMonomer_sp createMonomer(CandoDatabase_sp bdb);
    DEFAULT_CTOR_DTOR(OligomerPart_Link_O);
};


extern void setCandoDatabase(CandoDatabase_sp bdb, core::Lisp_sp lisp);







};
TRANSLATE(chem::OligomerPart_Base_O);
TRANSLATE(chem::OligomerPart_Monomer_O);
TRANSLATE(chem::OligomerPart_Link_O);
#endif //]
