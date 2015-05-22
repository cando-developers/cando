#ifndef gaff_LocationType_H
#define gaff_LocationType_H

#include "core/foundation.h"
#include "core/object.h"
#include "core/lisp.h"
#include "core/cons.h"



// define the type for flex and bison                                                      
// Location aware type
//
// Defined in cons.h
//
#define YYLTYPE ::core::LispParserPos



#define YYLLOC_DEFAULT(Current, Rhs, N)                                \
         do                                                                  \
	     if (N)                                                            \
             {                                                               \
		 (Current).first_line   = YYRHSLOC(Rhs, 1).first_line;         \
		 (Current).first_column = YYRHSLOC(Rhs, 1).first_column;       \
		 (Current).last_line    = YYRHSLOC(Rhs, N).last_line;          \
		 (Current).last_column  = YYRHSLOC(Rhs, N).last_column;        \
             }                                                               \
	     else                                                              \
             {                                                               \
		 (Current).first_line   = (Current).last_line   =              \
		     YYRHSLOC(Rhs, 0).last_line;                                 \
		 (Current).first_column = (Current).last_column =              \
		     YYRHSLOC(Rhs, 0).last_column;                               \
             }                                                               \
         while (0)


#endif // gaff_LocationType_H
