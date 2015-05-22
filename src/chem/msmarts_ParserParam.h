#ifndef msmarts_ParserParam_H
#define msmarts_ParserParam_H

#include "core/foundation.h"
#include "core/object.h"
#include "core/lisp.h"
#include "core/cons.h"


typedef struct tagSMsmartsParserParam
{
    /* The scanner extra pointer will point back to this lisp_SParserParam structure
     */
    void* scanner; // was yyscan_t

    core::Lisp_sp lisp;

    chem::SmartsRoot_sp	expression;
    map<string,int> atomTagCounter;
    stringstream*		msmartsErrorStream;


    istream*	lexStream;		
// This defines the stream that we are reading from



} msmarts_SParserParam;



#endif
