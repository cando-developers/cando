#ifndef	_candoBase_symbolTable_H
#define _candoBase_symbolTable_H

#include <clasp/core/foundation.h>

namespace candoBase
{

#define	CandoBasePkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkg,lispname,export) extern core::Symbol_sp cname;
#include <symbols_scraped_inc.h>
#undef DO_SYMBOL
#undef CandoBasePkg_SYMBOLS

}; /* candoBase */


namespace kw
{
#define	KeywordPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkg,lispname,export) extern core::Symbol_sp cname;
#include <symbols_scraped_inc.h>
#undef DO_SYMBOL
#undef KeywordPkg_SYMBOLS
}; /* kw */

#endif /* _candoBase_symbolTable_H */

