       
       
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


namespace candoBase {


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
    virtual void	startElement() {THROW_HARD_ERROR(BF("Subclass must implement"));};
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
    bool	parse(core::T_sp fIn, const string& fileName );
	bool	parseString(const char* xmlContent);
	void	shutDown();


};


};

#endif
