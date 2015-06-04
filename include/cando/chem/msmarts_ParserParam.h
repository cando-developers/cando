#ifndef msmarts_ParserParam_H
#define msmarts_ParserParam_H

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <clasp/core/cons.h>


typedef struct tagSMsmartsParserParam
{
    /* The scanner extra pointer will point back to this lisp_SParserParam structure
     */
    void* scanner; // was yyscan_t

    core::Lisp_sp lisp;

    chem::SmartsRoot_sp	expression;
  std::map<string,int> atomTagCounter;
    stringstream*		msmartsErrorStream;


  std::istream*	lexStream;		
// This defines the stream that we are reading from



} msmarts_SParserParam;



#endif
