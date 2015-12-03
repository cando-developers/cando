#ifndef	_geom_symbolTable_H
#define _geom_symbolTable_H

#include <clasp/core/foundation.h>

namespace geom
{

#define	GeomPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkg,lispname,export) extern core::Symbol_sp cname;
  #ifndef SCRAPING
    #include SYMBOLS_SCRAPED_INC_H
  #endif
#undef DO_SYMBOL
#undef GeomPkg_SYMBOLS

}; /* geom */


namespace kw
{
#define	KeywordPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkg,lispname,export) extern core::Symbol_sp cname;
  #ifndef SCRAPING
    #include SYMBOLS_SCRAPED_INC_H
  #endif
#undef DO_SYMBOL
#undef KeywordPkg_SYMBOLS
}; /* kw */

#endif /* _geom_symbolTable_H */

