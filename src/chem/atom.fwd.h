#ifndef atom_FWD_H
#define atom_FWD_H

#include <cando/chem/chemPackage.h>

/*! Defines whether each atom has a single built in Anchor restraint or
if it maintains a list of restraints which can contain mulitple
AnchorRestraints */
#define ATOMIC_ANCHOR 0


namespace chem
{
    typedef unsigned int AtomFlags;

    FORWARD(Atom);

    typedef gctools::Vec0<Atom_sp> VectorAtom;

    typedef core::Symbol_sp AtomAliasName;
    typedef core::Symbol_sp AtomType;
    typedef core::Symbol_sp MoeType;

    typedef	enum {
	undefinedConfiguration,
	S_Configuration,
	R_Configuration
    } ConfigurationEnum;

    
    typedef	enum { undefinedCenter, chiralCenter, prochiralCenter} StereochemistryType;



    extern core::Symbol_sp _sym__PLUS_configurationEnumConverter_PLUS_;
    extern core::Symbol_sp _sym__PLUS_stereochemistryTypeConverter_PLUS_;

};


    DECLARE_ENUM_SYMBOL_TRANSLATOR(chem::ConfigurationEnum,chem::_sym__PLUS_configurationEnumConverter_PLUS_);
    DECLARE_ENUM_SYMBOL_TRANSLATOR(chem::StereochemistryType,chem::_sym__PLUS_stereochemistryTypeConverter_PLUS_);


#endif // atom_FWD_H
