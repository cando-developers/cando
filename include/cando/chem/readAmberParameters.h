#ifndef	ReadAmberParameters_H //[
#define ReadAmberParameters_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>


#include <cando/chem/chemPackage.h>


namespace chem {

SMART(ForceField);
SMART(FFTypesDb);
SMART(FFNonbondDb);
SMART(FFStretchDb);
SMART(FFAngleDb);
SMART(FFPtorDb);
SMART(FFItorDb);

SMART(ReadAmberParameters );
class ReadAmberParameters_O : public core::CxxObject_O
{
    LISP_BASE1(core::CxxObject_O);
    LISP_CLASS(chem,ChemPkg,ReadAmberParameters_O,"ReadAmberParameters");
private:
	FFTypesDb_sp	_Types;
	ForceField_sp	_ForceField;
public:
	void	initialize();
private:
	// instance variables
public:
private:
        FFTypesDb_sp parseTypeRules(core::T_sp fin);
    FFNonbondDb_sp parseFirstNonbondDb(core::T_sp fin);
    FFStretchDb_sp parseStretchDb(core::T_sp fin);
    FFAngleDb_sp parseAngleDb(core::T_sp fin);
    FFPtorDb_sp parsePtorDb(core::T_sp fin);
    FFItorDb_sp parseItorDb(core::T_sp fin);
    void parseNonbondDb(core::T_sp fin, FFNonbondDb_sp );
    ForceField_sp parseAmberFormattedForceField(core::T_sp strm);
public:

    void readTypes(core::T_sp fin);
    void readParameters(core::T_sp fileName);

    FFTypesDb_sp getTypeRules() { return this->_Types; };
    ForceField_sp getForceField();

    ReadAmberParameters_O( const ReadAmberParameters_O& ss ); //!< Copy constructor

    DEFAULT_CTOR_DTOR(ReadAmberParameters_O);
};


};
TRANSLATE(chem::ReadAmberParameters_O);
#endif //]
