/*
    File: quickDom.cc
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
#define		DEBUG_LEVEL_FULL
//
// (C) 2004 Christian E. Schafmeister
//


#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <clasp/core/foundation.h>
#include <clasp/core/lisp.h>
#include <clasp/core/multipleValues.h>
#include <clasp/core/str.h>
#include <clasp/core/fileSystem.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/pathname.h>
#include <cando/adapt/mySaxInterface.h>
#include <cando/adapt/quickDom.h>
//#i n c l u d e "stringSet.h"
#include <clasp/core/wrappers.h>


namespace adapt {


class MySaxDomHandler : public MySaxHandlerDefault // , private XMLFormatTarget
{
public:
    // -----------------------------------------------------------------------
    //  Constructors and Destructor
    // -----------------------------------------------------------------------
    virtual ~MySaxDomHandler() {};



    // -----------------------------------------------------------------------
    //  Handlers for the SAX ContentHandler interface
    // -----------------------------------------------------------------------
    virtual void startElement();
    virtual void endElement();
    virtual void characters();
    virtual void ignorableWhitespace();
    QDomNode_O	*getTopNode() {return this->topNode;}
    void		setTopNode(QDomNode_O* n);
private:
    QDomNode_O	*_currentNode;
    QDomNode_O	*topNode;
};



void MySaxDomHandler::setTopNode(QDomNode_O* n)
{_OF();
    this->topNode = n;
    LOG(BF("Setting _currentNode to NULL"));
    this->_currentNode=NULL;
}

void	MySaxDomHandler::startElement() 
{
QDomNode_sp	newNode;
QDomNode_O	*nd;
string		attrName, attrVal;
int		i;

    if ( this->_currentNode == NULL ) 
    {
	nd = this->topNode;
//	LOG(BF("In top node address=%lx") % (nd ) );
    } else {
		    // Create a new node, give it its name and parent
	newNode = QDomNode_O::create();
	nd = newNode.get();
//	LOG(BF("Created new node address=%lx") % (nd ) );
        this->_currentNode->addChild(newNode);
    }

//    LOG(BF("Setting data for node: %lx  localName=%s") % (nd) % this->getLocalName() );
    nd->setLocalName(this->getLocalName());
    nd->setFileName(this->topNode->getFileName());
    nd->setLineNumber(this->getLineNumber());
    nd->setCharacters("");
    nd->setParent(this->_currentNode);
		// Set up the attributes
    for ( i = 0; i<this->getAttributeCount(); i++ ) 
    {
	this->getAttributeAtIndex( i, attrName, attrVal );
//	LOG(BF("Adding attribute: %s = %s") % (attrName.c_str()) % (attrVal.c_str() ) );
	nd->addAttribute( attrName, attrVal );
    }
//    LOG(BF("Setting _currentNode to %X") % (void*)(nd) );
    this->_currentNode = nd;
}


void	MySaxDomHandler::endElement() 
{
    if ( this->_currentNode != this->topNode ) 
    {
        this->_currentNode = this->_currentNode->getParent();
//	LOG(BF("Restored _currentNode to %X") % (void*)(this->_currentNode));
    }
}



void	MySaxDomHandler::characters() 
{
//    LOG(BF("About to appendCharacters(%s)") % this->getData() );
    ASSERTP(this->_currentNode!=NULL,"currentNode is NULL and you are trying to write characters to it");
    this->_currentNode->appendCharacters(this->getData());
}


void	MySaxDomHandler::ignorableWhitespace() 
{
}



    QDomNode_sp	QDomNode_O::create(const string& name)
{
    QDomNode_sp node = QDomNode_O::create();
    node->setLocalName(name);
    return node;
}


void	QDomNode_O::throwErrorForChildren()
{_OF();
VectorQDomNodes		children;
QDomNode_sp		child;

    LOG(BF("Entering throwErrorForChildren") );
    children = this->getChildren();
    LOG(BF("got children") );
    if ( children.size() != 0 ) {
        LOG(BF("children.size() != 0") );
        child = *(children.begin());
        LOG(BF("Defining error for unknown xml command: %s") % (child->getLocalName().c_str() ) );
	SIMPLE_ERROR(BF("There is an unknown command: %s in file(%s) line(%d)") % child->getLocalName() % this->getFileName() % this->getLineNumber() );
    }
}




void	QDomNode_O::throwErrorForChildrenWithoutName(string nm)
{_OF();
    VectorQDomNodes		children;
    QDomNode_sp		child;

    LOG(BF("Entering throwErrorForChildren") );
    children = this->getChildren();
    LOG(BF("got children") );
    for ( VectorQDomNodes::iterator x=children.begin(); x!= children.end(); x++ )
    {
	child = *x;
	if ( child->getLocalName() != nm ) {
            LOG(BF("Defining error for unknown xml command: %s") % (child->getLocalName().c_str() ) );
	    SIMPLE_ERROR(BF("There is an unknown command: %s in file(%s) line(%d)") % child->getLocalName() % this->getFileName() % this->getLineNumber());
	}
    }
}

//
//      dataIsAllWhiteSpace
//
//      Return true if the data is all white space
bool    QDomNode_O::dataIsAllWhiteSpace()
{_OF();
    bool                        sawChar;
    string                      val;
    string::iterator    it;
    sawChar = false;
    val = this->characters.str();
    for (it=val.begin();it!=val.end();it++)
    {
        if (!isspace(*it))
        {
            LOG(BF( "isspace was false char=|%c|\n")% *it );
            return false;
        }
    }
    return true;
}



int	QDomNode_O::dataCountNewLines()
{
string			val;
string::iterator	it;
int			newLines;
    newLines = 0;
    val = this->characters.str();
    for (it=val.begin();it!=val.end();it++){
	if (*it == '\n' ) {
	    newLines++;
	}
    }
    return newLines;
}

CL_LISPIFY_NAME("isLeaf");
CL_DEFMETHOD     bool QDomNode_O::isLeaf()
    {
	return this->_children.size() == 0;
    }

CL_LISPIFY_NAME("children");
CL_DEFMETHOD core::List_sp QDomNode_O::children()
{
    core::Cons_sp first, cur;
    first = core::Cons_O::create();
    cur = first;
    for ( QDomNode_O::iterator ci=this->begin_Children(); ci!=this->end_Children(); ci++ )
    {
	core::Cons_sp one = core::Cons_O::create(*ci);
	cur->setCdr(one);
	cur = one;
    }
    return oCdr(first);
}



CL_LISPIFY_NAME("writeToFileName");
CL_DEFMETHOD void	QDomNode_O::writeToFileName( string fileName )
{
    std::ofstream	out;
    out.open(fileName.c_str(), std::ios::out );
    this->writeXml( "", out );
    out.close();
}

CL_LISPIFY_NAME("appendToFileName");
CL_DEFMETHOD void	QDomNode_O::appendToFileName( string fileName )
{
    std::ofstream	out;
    out.open(fileName.c_str(), std::ios::app );
    this->writeXml( "", out );
    out.close();
}


void	QDomNode_O::writeXml( string prefix, std::ostream& out )
{
    gctools::Vec0<QDomNode_sp>::iterator	it;
map<string,string>::iterator	attr;
map<string,string>::iterator	abegin;
map<string,string>::iterator	aend;
bool				dataWhiteSpace;
int				dataNewLines;
bool				usePrefix;
    out << prefix << "<" << this->getLocalName();
	out.flush();
    abegin = this->attributes.begin();
    aend = this->attributes.end();
    LOG(BF("A") );
    for ( attr=abegin; attr!=aend; attr++ ) {
	out << " " << (*attr).first << "=\"" << (*attr).second << "\"";
	out.flush();
    }
    LOG(BF("B") );
    dataWhiteSpace = this->dataIsAllWhiteSpace();
    dataNewLines = 0;
    LOG(BF("C") );
    if ( !dataWhiteSpace ) {
        LOG(BF("CA") );
	dataNewLines = this->dataCountNewLines();
    }
    if ( this->_children.begin() == this->_children.end() && dataWhiteSpace ) {
        LOG(BF("CB") );
	out << "/>\n";
	out.flush();
    } else {
        LOG(BF("CC") );
	out << ">";
		//
		// If there are children then start them
		// on the next line
		//
        LOG(BF("CD") );
	usePrefix = false;
	if ( this->_children.begin() != this->_children.end() ) {
	    out << std::endl;
	    usePrefix = true;
	}
        LOG(BF("Ce") );
	out.flush();
	for ( it = this->_children.begin(); it!=this->_children.end(); it++ ) {
            LOG(BF("Cea") );
	    if ( !(*it) ) {
		SIMPLE_ERROR(BF("Bad child"));
	    }
	    (*it)->writeXml(prefix+" ",out);
	}
        LOG(BF("Cf") );
	out.flush();
	if ( !dataWhiteSpace ) {
	    if ( dataNewLines > 0 ) {
		out << std::endl;
		usePrefix = true;
	    }
	    out << this->getCharacters();
	    out.flush();
	}
        LOG(BF("Cg") );
	if ( usePrefix ) {
	    out << prefix;
	}
        LOG(BF("Ch") );
	out << "</" << this->getLocalName() << ">\n";
        LOG(BF("Ci") );
	out.flush();
    }
    LOG(BF("D") );
}

#ifdef	_DEBUG
_CrtMemState	memoryState;
#endif


#if 0
QDomNode_O::QDomNode_O(string ln) 
{
    REF_CREATE("QDomNode_O");
    this->localName = ln;
#ifdef	SAVE_XML_FILE_NAME
    this->_FileName = "";
#endif
}

QDomNode_O::QDomNode_O(string ln,string data) {
    REF_CREATE("QDomNode_O");
    this->localName = ln;
    this->characters.str(data);
#ifdef	SAVE_XML_FILE_NAME
    this->_FileName = "";
#endif
}
#endif




void	QDomNode_O::dumpChildToStream( std::ostream& o, const string& prefix )
{
    this->writeXml(prefix,o);
}

vector<string>	QDomNode_O::getDataAsVectorOfStrings()
{
string		s;
vector<string>	vs;
unsigned int	cur;
int		start, stop;
    s = this->getCharacters();
    cur = 0;
    while (cur <s.size()) {
	if ( isspace(s[cur]) || s[cur]<' ' ) {
	    cur++;
	    continue;
	}
	start = cur;
	while (cur<s.size() && !isspace(s[cur])) {
	    cur++;
	}
	stop = cur-1;
	vs.push_back(s.substr(start,stop-start+1));
    }
    return vs;
}

#if 0 //[
StringSet_sp	QDomNode_O::getDataAsStringSet()
{
string		s, str;
StringSet_sp	ss;
unsigned int	cur;
int		start, stop;

    ss = StringSet_O::create(p);
    s = this->getCharacters();
    cur = 0;
    while (cur <s.size()) {
	if ( isspace(s[cur]) || s[cur]<' ' ) {
	    cur++;
	    continue;
	}
	start = cur;
	while (cur<s.size() && !isspace(s[cur])) {
	    cur++;
	}
	stop = cur-1;
	str = s.substr(start,stop-start+1);
	LOG(BF("adding data (%s) to set of strings") % (str.c_str() ) );
	ss->insert(str);
    }
    return ss;
}
#endif //]


vector<int>	QDomNode_O::getDataAsVectorOfInts()
{
string		s;
vector<int>	vs;
unsigned int	cur;
int		start, stop;
    s = this->getCharacters();
    cur = 0;
    while (cur <s.size()) {
	if ( isspace(s[cur]) || s[cur] <' ' ) {
	    cur++;
	    continue;
	}
	start = cur;
	while (cur<s.size() && !isspace(s[cur])) {
	    cur++;
	}
	stop = cur-1;
	vs.push_back(atoi((s.substr(start,stop-start+1)).c_str()));
    }
    return vs;
}

vector<double>	QDomNode_O::getDataAsVectorOfDoubles()
{_OF();
string		s;
vector<double>	vs;
uint	cur;
int		start, stop;
double		dVal;
    s = this->getCharacters();
    cur = 0;
    while (cur <s.size()) {
	if ( isspace(s[cur]) || s[cur] < ' ' ) {
	    cur++;
	    continue;
	}
	start = cur;
	while (cur<s.size() && !isspace(s[cur])) {
	    cur++;
	}
	stop = cur-1;
	LOG(BF("convertToVectorOfDoubles converting: |%s|") % (s.substr(start,stop-start+1).c_str() ) );
	dVal = atof((s.substr(start,stop-start+1)).c_str());
	LOG(BF("convertToVectorOfDoubles value = %lf") % (dVal ) );
	vs.push_back(dVal);
	LOG(BF("convertToVectorOfDoubles pushed value") );
    }
    return vs;
}












std::ostream& QDomNode_O::dumpToStream(std::ostream& o) {
    this->dumpChildToStream(o, "");
    return o;
}




CL_LISPIFY_NAME("dump");
CL_DEFMETHOD void	QDomNode_O::dump() {
    this->dumpChildToStream(std::cout, "");
}


/////////////////////////////////////////////////////////////////////////
//
// QDomNode methods
//
#if 0
iterator	QDomNode_O::privateNextChildWithName( const string& name, const iterator& cur ) {
iterator	newCur;
    newCur = cur;
    while ( 1 ) {
        if ( newCur == this->_children.end() ) goto FAILED;
	if ( (*newCur)->getLocalName() == name ) {
	    return newCur;
	}
	newCur ++;
    }
    return newCur;
FAILED:
    return NULL;
}
#endif

VectorQDomNodes QDomNode_O::getChildrenWithName( const string& name ) {
VectorQDomNodes childs;
VectorQDomNodes::iterator	it;

    for ( it=this->_children.begin(); it!=this->_children.end(); it++ ) {
	if ( (*it)->getLocalName() == name ) {
	    childs.push_back(*it);
	}
    }
    return childs;
}

CL_LISPIFY_NAME("childWithName");
CL_DEFMETHOD QDomNode_sp	QDomNode_O::childWithName( const string& name ) 
{_OF();
QDomNode_sp			child;
iterator	it;
int				children;

    children = 0;
    for ( it=this->_children.begin(); it!=this->_children.end(); it++ ) {
	if ( (*it)->getLocalName() == name ) {
	    child = *it;
	    children++;
	}
    }
    if ( children > 1 )
    {
	SIMPLE_ERROR(BF("There can only be one child with the name(%s) in XML file(%s) at line(%d)") % name % this->getFileName() % this->getLineNumber() );
    }
    if ( children < 1 ) {
	SIMPLE_ERROR(BF("There must be one child with the name(%s) in XML file(%s) at line(%d)") % name % this->getFileName() % this->getLineNumber() );
    }
    return child;
}


CL_LISPIFY_NAME("onlyChild");
CL_DEFMETHOD QDomNode_sp	QDomNode_O::onlyChild()
{_OF();
QDomNode_sp			child;
    if ( this->_children.size() != 1 ) {
        SIMPLE_ERROR(BF("This node must have one and only one child file(%s) line(%d) ") % this->getFileName() % this->getLineNumber() );
    }
    child = *(this->_children.begin());
    return child;
}



VectorQDomNodes QDomNode_O::getChildrenWithNameAndRemoveThem( const string& name )
{_OF();
VectorQDomNodes 			childs;
iterator		it;
    LOG(BF("QDomNode_O::getChildrenWithNameAndRemoveThem") );
    LOG(BF("Coming in there are %d children") % (this->_children.size() ) );
    it = this->_children.begin();
    while ( it != this->_children.end() ) {
	if ( (*it)->getLocalName() == name ) {
	    childs.push_back(*it);
	    it = this->_children.erase(it);
	    continue;
	}
	it++;
    }
    LOG(BF("I found %d children with the name: %s") % (childs.size()) % (name.c_str() ) );
    LOG(BF("Going out there are %d children") % (this->_children.size() ) );

    return childs;
}





CL_LISPIFY_NAME("hasChildrenWithName");
CL_DEFMETHOD bool	QDomNode_O::hasChildrenWithName( const string& name )
{
    return this->countChildrenWithName(name) > 0;
}



int	QDomNode_O::countChildrenWithName( const string& name ) {
VectorQDomNodes childs;
iterator	it;
int				num;
    num = 0;
    for ( it=this->_children.begin(); it!=this->_children.end(); it++ ) {
	if ( (*it)->getLocalName() == name ) num++;
    }
    return num;
}



#ifdef	USEBOOSTPYTHON
#if 0
boost::python::list	QDomNode_O::python_gatherSubNodesWithName( const string& nm ) {
boost::python::list			vals;
VectorQDomNodes			childs;
VectorQDomNodes::iterator	it;

    LOG(BF("Entered python_childrenWithName, looking for name: %s") % (nm.c_str() ) );
    childs = this->gatherSubNodesWithName(nm);
    for ( it = childs.begin(); it!=childs.end(); it++ ) {
	vals.append(*it);
    }
    return vals;
}

boost::python::list	QDomNode_O::python_childrenWithName( const string& nm ) {
boost::python::list			vals;
VectorQDomNodes			childs;
VectorQDomNodes::iterator	it;

    LOG(BF("Entered python_childrenWithName, looking for name: %s") % (nm.c_str() ) );
    childs = this->getChildrenWithName(nm);
    for ( it = childs.begin(); it!=childs.end(); it++ ) {
	vals.append(*it);
    }
    return vals;
}

boost::python::long_ QDomNode_O::python_hasChildrenWithName( const string& nm ) {
boost::python::list			vals;
VectorQDomNodes			childs;
VectorQDomNodes::iterator	it;

    childs = this->getChildrenWithName(nm);
    for ( it = childs.begin(); it!=childs.end(); it++ ) {
	return boost::python::long_(1);
    }
    return boost::python::long_(0);
}

boost::python::list	QDomNode_O::python_getChildren() {
boost::python::list			vals;
VectorQDomNodes			childs;
VectorQDomNodes::iterator	it;

    childs = this->getChildren();
    for ( it = childs.begin(); it!=childs.end(); it++ ) {
	vals.append(*it);
    }
    return vals;
}
#endif
#endif	// USEBOOSTPYTHON




void	QDomNode_O::eraseAll() {
    this->parent = NULL;
    this->localName = "";
    this->characters.str("");
    this->attributes.erase(this->attributes.begin(), this->attributes.end());
    this->_children.clear();
}





//////////////////////////////////////////////////////////////////////////
//
//
//
// Parsing routines.
//
//
//


#if 0
void	QDomNode_O::parseStream(core::Stream_sp fIn)
{_OF();
    MySaxDomHandler	handler;
    bool		res;

		//
		// First erase everything from the Node
		//
    this->eraseAll();

#ifdef	SAVE_XML_FILE_NAME
    this->_FileName = fileName;
#endif
    handler.setTopNode(this);

    LOG(BF("myDomParseFileName about to parse file") );

    LOG(BF("Initializing parser"));
    MySaxParser parser;
    res = parser.initialize();
    if ( !res ) {
	printf ("error during XML initialization\n" );
	LOG(BF("Error during XML initialization") );
	parser.shutDown();
    }
    LOG(BF("myDomParseFileName About to parse XML file") );
    parser.setContentHandler(&handler);
    res = parser.parse(fIn,fileName);
    if ( !res ) {
	printf( "Could not parse file\n" );
	LOG(BF("myDomParseFileName failed to parse file") );
    }
    parser.shutDown();

}
#endif

gc::Nilable<QDomNode_sp>  QDomNode_O::parse(core::T_sp stream)
{
    MySaxDomHandler	handler;
    bool		res;
    QDomNode_sp	node;

    node = QDomNode_O::create();
    handler.setTopNode(node.get());
    LOG(BF("QDomParseString about to parse string") );
    MySaxParser parser;
    res = parser.initialize();
    if ( !res ) {
	printf ("error durint XML initialization\n" );
	LOG(BF("Error during XML initialization") );
	parser.shutDown();
	return _Nil<core::T_O>();
    }
//    handler.setExecutionContext(this);
    LOG(BF("QDomParseString About to parse XML file") );
    parser.setContentHandler(&handler);
    res = parser.parse(stream);
    if ( !res ) {
	printf( "Could not parse file\n" );
	LOG(BF("QDomParseString failed to parse file") );
    }
    LOG(BF("QDomParseString shutting down") );
    parser.shutDown();
    LOG(BF("About to return") );
    return node;
}



CL_LISPIFY_NAME("asString");
CL_DEFMETHOD string	QDomNode_O::asString()
{_OF();
stringstream	ss;
    this->dumpToStream(ss);
    return ss.str();
}


CL_LISPIFY_NAME("addChild");
CL_DEFMETHOD void	QDomNode_O::addChild(QDomNode_sp child)
{
    if ( !child ) {
	SIMPLE_ERROR(BF("You tried to add an empty child"));
    }
    this->_children.push_back(child);
}


void	QDomNode_O::fillVectorQDomNodesIfNameIs(int depth,QDomNode_sp me, VectorQDomNodes& vnodes, const string& name)
{//
    VectorQDomNodes	children;
    QDomNode_sp		child;
    if ( me->getLocalName() == name ) {
//	LOG(BF("Found name: %s") % name.c_str()  );
	vnodes.push_back(me);
	return;
    }
    children = this->getChildren();
    for ( VectorQDomNodes::iterator x=children.begin(); x!= children.end(); x++ ) {
	child = *x;
	child->fillVectorQDomNodesIfNameIs(depth+1,child,vnodes,name);
    }
}


VectorQDomNodes QDomNode_O::gatherSubNodesWithName(const string& name)
{//
    VectorQDomNodes	vnodes;
    VectorQDomNodes	children;
    QDomNode_sp		child;
    children = this->getChildren();
    for ( VectorQDomNodes::iterator x=children.begin(); x!= children.end(); x++ ) {
	child = *x;
	child->fillVectorQDomNodesIfNameIs(1,child,vnodes,name);
    }
    return vnodes;
}






};
