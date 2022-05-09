/*
    File: mySaxInterface.h
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
       
       
//
// (C) 2004 Christian E. Schafmeister
//


#ifndef	_MYSAXINTERFACE
#define	_MYSAXINTERFACE

#include <expat.h>
#include <string>

#include <clasp/core/foundation.h>
#include <clasp/core/lisp.h>
#include <clasp/core/lispStream.fwd.h>


namespace adapt {


/////////////////////////////////////////////////////////////////////////////
//
// MySaxHandlerDefault
//
// The class that my code implements to handle the SAX messages in a system
// independant way.
//

class	MySaxHandlerDefault {
public:


	void	setLineNumber(int l) {this->iLineNumber = l; };
	int	getLineNumber() {return this->iLineNumber; };

	void	setLocalName(const char* s) { this->cPLocalName = s; };
	const char*	getLocalName() {return this->cPLocalName;}

	void	setData(const char* s) { this->cPData = s;};
	const char*	getData() {return this->cPData;}

	void	setDataLength(uint l) {this->iDataLength = l;}
	int	getDataLength() {return this->iDataLength;}

	string	getAttribute(char* str);
	int	getAttributeCount();
	void	getAttributeAtIndex(int ind, string& attrName, string& attrVal );

	void	setAttributes(const XML_Char ** att) { this->attributes = att; };

	// Sub-class these
    virtual void	startElement() {THROW_HARD_ERROR("Subclass must implement");};
  virtual void	endElement() {return;};
  virtual void	characters() {return;};
  virtual void	ignorableWhiteSpace() {return;};

	virtual ~MySaxHandlerDefault() {};

private:
// Common stuff
    const XML_Char		**attributes;
    uint			iLineNumber;
    const char*			cPLocalName;
    const char*			cPData;
    uint			iDataLength;
};



////////////////////////////////////////////////////////////////////////////////
//
// MySaxParser
//
//  System independent SAX parser
//
//

class	MySaxParser {
private:
	XML_Parser		parser;
//	ProgressMonitor			_Progress;
public:
        MySaxParser();
	virtual ~MySaxParser() {};
        MySaxHandlerDefault	*handler;
	bool	initialize();
	void	setContentHandler(MySaxHandlerDefault* ahand );
	XML_Parser	getXMLParser() { return this->parser; };
        bool	parse(core::T_sp fIn);
	void	shutDown();

};


};

#endif
