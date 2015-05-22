       
       
#ifndef	TypeAssignmentRules_H //[
#define TypeAssignmentRules_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/chemInfo.h>

#include <cando/chem/chemPackage.h>


namespace chem {

    SMART(OneTypeRule );
    class OneTypeRule_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,OneTypeRule_O,"OneTypeRule");
public:
//	void	archive(core::ArchiveP node);
	void	initialize();
private:
	core::Symbol_sp _AssignTypeName;
	ChemInfo_sp	_Match;
public:
public:

	void setAssignTypeName(core::Symbol_sp nm) { this->_AssignTypeName = nm;};
	core::Symbol_sp  getAssignTypeName()	{return this->_AssignTypeName;};
	void setCode(ChemInfo_sp code) { this->_Match = code; };

	bool testMatch(Atom_sp a);


	OneTypeRule_O( const OneTypeRule_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(OneTypeRule_O);
    };

    SMART(TypeAssignmentRules );
    class TypeAssignmentRules_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,TypeAssignmentRules_O,"TypeAssignmentRules");
public:
//	void	archive(core::ArchiveP node);
	void	initialize();
private:
	WildElementDict_sp	_WildElementDict;
        gctools::Vec0<OneTypeRule_sp>	_Rules;
public:
public:

	void setWildElementDict(WildElementDict_sp dict) { this->_WildElementDict = dict;};
	WildElementDict_sp getWildElementDict(WildElementDict_sp dict) { return this->_WildElementDict;};

	void appendRule(OneTypeRule_sp rule) { this->_Rules.push_back(rule);};

	core::Symbol_sp calculateType(Atom_sp a);

	TypeAssignmentRules_O( const TypeAssignmentRules_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(TypeAssignmentRules_O);
    };



};
TRANSLATE(chem::OneTypeRule_O);
TRANSLATE(chem::TypeAssignmentRules_O);
#endif //]
