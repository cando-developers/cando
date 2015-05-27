#ifndef	_chem_symbolTable_H
#define _chem_symbolTable_H

#include <clasp/core/foundation.h>

namespace chem
{

#define	ChemPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkg,lispname,export) extern core::Symbol_sp cname;
#include <cando/chem/generated/symbols_scraped_inc.h>
#undef DO_SYMBOL
#undef ChemPkg_SYMBOLS

}; /* chem */


namespace chemkw
{
#define	ChemKwPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkg,lispname,export) extern core::Symbol_sp cname;
#include <cando/chem/generated/symbols_scraped_inc.h>
#undef DO_SYMBOL
#undef ChemKwPkg_SYMBOLS
}; /* chemkw */

namespace kw {
#define KeywordPkg_SYMBOLS
#define DO_SYMBOL(cname, idx, pkg, lispname, export) extern core::Symbol_sp cname;
#include <cando/chem/generated/symbols_scraped_inc.h>
#undef DO_SYMBOL
#undef KeywordPkg_SYMBOLS
}; /* kw package */

#endif /* _chem_symbolTable_H */
