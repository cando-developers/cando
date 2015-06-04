#ifndef gaff_ParserParam_H
#define gaff_ParserParam_H

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <clasp/core/cons.h>


typedef struct tagSGaffParserParam
{
    /* The scanner extra pointer will point back to this lisp_SParserParam structure
     */
    void* scanner; // was yyscan_t

    chem::AntechamberRoot_sp	expression;
    chem::WildElementDict_sp wildElementDictionary;
    stringstream*		gaffErrorStream;


  std::istream*	acLexStream;		
// This defines the stream that we are reading from



} gaff_SParserParam;



#endif
