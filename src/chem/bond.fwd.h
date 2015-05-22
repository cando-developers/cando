#ifndef	bond_FWD_H
#define bond_FWD_H

#include "core/common.h"

#include "chemPackage.h"


namespace chem
{

FORWARD(Bond);

typedef	enum 	BondOrderEnum 
{
    noBond,
    singleBond,
    doubleBond,
    tripleBond,
    aromaticBond,
    hydrogenBond,
    virtualBond,
    dashedSingleBond,
    dashedDoubleBond
} BondOrder;


    extern core::Symbol_sp _sym__PLUS_bondOrderToSymbolConverter_PLUS_;
};



DECLARE_ENUM_SYMBOL_TRANSLATOR(chem::BondOrder,chem::_sym__PLUS_bondOrderToSymbolConverter_PLUS_);

#endif // bond_FWD_H
