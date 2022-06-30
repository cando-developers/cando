/*
    File: atom.fwd.h
*/
/*
Open Source License
Copyright (c) 2016, Christian E. Schafmeister
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 
This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
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
    typedef core::T_sp AtomType;
    typedef core::Symbol_sp MoeType;

    /*!
      RightHanded_Configuration and
      LeftHanded_Configuration
      means bonds to atoms (a1 a2 a3 a4)
      point a4 away then RightHanded_Configuration a1->a2->a3 right handed
      point a4 away then LeftHanded_Configuration a1->a2->a3 left handed
    */

    typedef	enum {
	undefinedConfiguration,
	S_Configuration,
	R_Configuration,
        RightHanded_Configuration,
        LeftHanded_Configuration
    } ConfigurationEnum;

    typedef	enum { undefinedCenter, chiralCenter, prochiralCenter} StereochemistryType;

    extern core::Symbol_sp& _sym__PLUS_configurationEnumConverter_PLUS_;
    extern core::Symbol_sp& _sym__PLUS_stereochemistryTypeConverter_PLUS_;

};


    DECLARE_ENUM_SYMBOL_TRANSLATOR(chem::ConfigurationEnum,chem::_sym__PLUS_configurationEnumConverter_PLUS_);
    DECLARE_ENUM_SYMBOL_TRANSLATOR(chem::StereochemistryType,chem::_sym__PLUS_stereochemistryTypeConverter_PLUS_);


#endif // atom_FWD_H
