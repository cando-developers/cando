/*
    File: quickDom.h
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


#ifndef	adapt_QUICKDOM_H
#define	adapt_QUICKDOM_H

#define	SAVE_XML_FILE_NAME
//#define	DONT_SAVE_XML_FILE_NAME

#include <stdio.h>
#include <math.h>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <map>
#include <set>
#include <exception>
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <cando/adapt/adaptPackage.fwd.h>

namespace adapt {

SMART(QDomNode);


    typedef	gctools::Vec0<QDomNode_sp>	VectorQDomNodes;


SMART(QDomNode );
    class QDomNode_O : public core::CxxObject_O
{
    LISP_CLASS(adapt,AdaptPkg,QDomNode_O,"QDomNode",core::CxxObject_O);

private:
    void	dumpChild( const string& prefix);
    void	dumpChildToStream( std::ostream& o, const string& prefix);

    QDomNode_O			*parent;
    string			localName;
#ifdef	SAVE_XML_FILE_NAME
    string			_FileName;
#endif
    int				lineNumber;
    stringstream		characters;
    VectorQDomNodes		_children;
    map<string,string>		attributes;

    void	fillVectorQDomNodesIfNameIs(int depth, QDomNode_sp me, VectorQDomNodes& vnodes, const string& name);

public:
    static QDomNode_sp create(const string& nodeName);
    typedef VectorQDomNodes::iterator	iterator;
    static gc::Nilable<QDomNode_sp> parse(core::T_sp stream);
public:
    
    QDomNode_O	*getParent()		{return this->parent;};
    void	setParent(QDomNode_O *par)	{this->parent = par;};

    /*! call proc with each child as an argument, return a tree of objects */
//    T_sp walk(core::Function_sp exec, core::Cons_sp args, Environment_sp environ, Lisp_sp lisp);

    iterator begin_Children() { return this->_children.begin(); };
    iterator end_Children() { return this->_children.end(); };

CL_LISPIFY_NAME("isNamed");
CL_DEFMETHOD     bool	isNamed(const string& nm)	{return this->localName==nm;};
CL_LISPIFY_NAME("getLocalName");
CL_DEFMETHOD     string	getLocalName()			{return this->localName;};
CL_LISPIFY_NAME("setLocalName");
CL_DEFMETHOD     void	setLocalName(const string& nm)		{this->localName = nm;};
#ifdef	SAVE_XML_FILE_NAME
    void	setFileName(const string& fn)	{ this->_FileName = fn; };
    string	getFileName()			{ return this->_FileName; };
#else
    void	setFileName(const string& fn)	{};
    string	getFileName()			{ return "-no file name, turn debugging on-";};
#endif
    void	setLineNumber(int line)		{this->lineNumber = line;};
CL_LISPIFY_NAME("getLineNumber");
CL_DEFMETHOD     int		getLineNumber()			{return this->lineNumber;};
CL_LISPIFY_NAME("getCharacters");
CL_DEFMETHOD     string	getCharacters()			{return this->characters.str();};
    stringstream&	charactersStream()	{return this->characters;};
		//
		// getData same as getCharacters
		//
    bool	dataIsAllWhiteSpace();
    int		dataCountNewLines();
    string	getData()			{return this->characters.str();};
    vector<string>	getDataAsVectorOfStrings();
//    StringSet_sp		getDataAsStringSet();
    vector<int>		getDataAsVectorOfInts();
    vector<double>	getDataAsVectorOfDoubles();


CL_LISPIFY_NAME("setCharacters");
CL_DEFMETHOD     void	setCharacters(const string& chrs)	{this->characters.str(""); this->characters << chrs;};
    void	appendCharacters(const string& chrs)	{this->characters << chrs; };

    void	writeXml( string prefix, std::ostream& out );
    void	writeToFileName( string fileName );
    void	appendToFileName( string fileName );
CL_LISPIFY_NAME("saveAs");
CL_DEFMETHOD     void	saveAs( string fileName ) { this->writeToFileName(fileName);};
    std::ostream&	dumpToStream(std::ostream& o);
    void	dump();

    bool	isLeaf();
    bool		hasChildrenWithName(const string& nm);
    int			countChildrenWithName(const string& nm);
    QDomNode_sp		childWithName(const string& nm);
    QDomNode_sp		onlyChild();


    core::List_sp	childrenWithName(const string& nm);
    VectorQDomNodes	getChildrenWithName(const string& nm);
    VectorQDomNodes	getChildrenWithNameAndRemoveThem(const string& nm);
    VectorQDomNodes	getChildren()	{return this->_children;};
#if 0
#ifdef	USEBOOSTPYTHON
    boost::python::list		python_childrenWithName(const string& nm);
    boost::python::list		python_getChildren();
    boost::python::list		python_gatherSubNodesWithName(const string& nm);
    boost::python::long_	python_hasChildrenWithName(const string& nm);
#endif
#endif
    void	addChild(QDomNode_sp child );
CL_LISPIFY_NAME("numberOfChildren");
CL_DEFMETHOD     int		numberOfChildren() { return this->_children.size();};

    core::List_sp children();

CL_LISPIFY_NAME("hasAttribute");
CL_DEFMETHOD     bool	hasAttribute(const string& at ) {
	return this->attributes.count(at)>0;
    };

CL_LISPIFY_NAME("removeAttribute");
CL_DEFMETHOD     void	removeAttribute(const string& at)
    {
        this->attributes.erase(at);
    }

CL_LISPIFY_NAME("getAttribute");
CL_DEFMETHOD     string	getAttributeValueNoError(const string& at) {
	if ( this->attributes.count(at)>0 ) {
	    return this->attributes[at];
	} else {
	    return "";	// Don't throw an error, just return a blank string
	}
    };
    string	getAttributeValueDefault(const string& at,const string& df) {
	if ( this->attributes.count(at)>0 ) {
	    return this->attributes[at];
	} else {
	    return df;	// return the default value
	}
    };
    string	getAttributeValue(const string& at)
    {_OF();
	if ( this->attributes.count(at)>0 ) {
	    return this->attributes[at];
	} else {
	    SIMPLE_ERROR(BF("Missing attribute: %s file(%s) line(%d)") % at % this->getFileName() % this->getLineNumber());
	}
    };

	//
	// Overloaded versions of getAttribute
	//
    void getAttribute(const string& name, string& v) {
	v = this->getAttributeString(name);
    }
    void getAttribute(const string& name, string& v, string defval) {
	v = this->getAttributeStringDefault(name,defval);
    }


    void getAttribute(const string& name, int& v) {
	v = this->getAttributeInt(name);
    }
    void getAttribute(const string& name, int& v, int defval) {
	v = this->getAttributeIntDefault(name,defval);
    }

    void getAttribute(const string& name, double& v) {
	v = this->getAttributeDouble(name);
    }
    void getAttribute(const string& name, double& v, double defval) {
	v = this->getAttributeDoubleDefault(name,defval);
    }

    void getAttribute(const string& name, float& v) {
	v = this->getAttributeFloat(name);
    }
    void getAttribute(const string& name, float& v, float defval) {
	v = this->getAttributeFloatDefault(name,defval);
    }

    void getAttribute(const string& name, bool& v ) {
	v = this->getAttributeBool(name);
    }
    void getAttribute(const string& name, bool& v, bool defval) {
	v = this->getAttributeBoolDefault(name,defval);
    }




		//
		// Attribute type conversion routines.
		// These routines could throw exceptions if they
		// found invalid data
		//
CL_LISPIFY_NAME("getAttributeDouble");
CL_DEFMETHOD     double	getAttributeDouble(const string& at) {
		return atof(this->getAttributeValue(at).c_str());
    };
    double	getAttributeFloat(const string& at) {
		return atof(this->getAttributeValue(at).c_str());
    };
    double	getAttributeDoubleDefault(const string& at,double df) {
			if ( this->attributes.count(at) ) {
			    return atof(this->getAttributeValue(at).c_str());
			} else {
			    return df;
			}
    };
    float	getAttributeFloatDefault(const string& at,float df) {
			if ( this->attributes.count(at) ) {
			    return atof(this->getAttributeValue(at).c_str());
			} else {
			    return df;
			}
    };
    string	getAttributeString(const string& at) {
			return this->getAttributeValue(at);
    };
    string	getAttributeStringDefault(const string& at, const string& df) {
			if ( this->attributes.count(at) ) {
			    return this->getAttributeValue(at);
			}
			return df;
    };
    int		getAttributeInt(const string& at) {
			return atoi(this->getAttributeValue(at).c_str());
    };
CL_LISPIFY_NAME("getAttributeIntDefault");
CL_DEFMETHOD     int		getAttributeIntDefault(const string& at,int df) {
			if ( this->attributes.count(at) ) {
			    return atoi(this->getAttributeValue(at).c_str());
			} else {
			    return df;
			}
    };
    bool	getAttributeBool(const string& at)
    {_OF();
	string val;
	val = this->getAttributeValue(at);
	if ( val== "true") return true;
	if ( val== "True") return true;
	if ( val== "TRUE") return true;
	if ( val== "yes") return true;
	if ( val== "Yes") return true;
	if ( val== "YES") return true;
	if ( val== "1" ) return true;
	if ( val== "false") return false;
	if ( val== "False") return false;
	if ( val== "FALSE") return false;
	if ( val== "no") return false;
	if ( val== "No") return false;
	if ( val== "NO") return false;
	if ( val== "0" ) return false;
	SIMPLE_ERROR(BF("Illegal boolean value: %s file(%s) line(%d) ") % val % this->getFileName() % this->getLineNumber());
    };

    bool	getAttributeBoolDefault(const string& at,bool d) {
		string val;
		if ( this->hasAttribute(at) ) {
		    val = this->getAttributeValue(at);
		    if ( val== "true") return true;
		    if ( val== "True") return true;
		    if ( val== "TRUE") return true;
		    if ( val== "yes") return true;
		    if ( val== "Yes") return true;
		    if ( val== "YES") return true;
		    if ( val== "1" ) return true;
		    if ( val== "false") return false;
		    if ( val== "False") return false;
		    if ( val== "FALSE") return false;
		    if ( val== "no") return false;
		    if ( val== "No") return false;
		    if ( val== "NO") return false;
		    if ( val== "0" ) return false;
		}
		return d;
    };



    void addAttribute(const string& at, const string& val) {
		this->attributes[at] = val;
    };

    void addAttribute(const string& at, double val ) {
	this->addAttributeDoubleScientific(at,val);
    };

    void addAttribute(const string& at, float val ) {
	this->addAttributeFloatScientific(at,val);
    };

    void addAttribute(const string& at, int val ) {
	this->addAttributeInt(at,val);
    };

    void addAttribute(const string& at, bool val ) {
	this->addAttributeBool(at,val);
    };

CL_LISPIFY_NAME("addAttributeString");
CL_DEFMETHOD     void addAttributeString(const string& at, const string& val) {
		this->addAttribute(at,val);
    };
CL_LISPIFY_NAME("addAttributeDouble");
CL_DEFMETHOD     void addAttributeDouble(const string& at, double val,int width,int prec ) {
        std::ostringstream temp;
	temp << std::setiosflags(std::ios::fixed)
             << std::setw(width)
             << std::setprecision(prec)
		<< val;
	this->addAttribute(at,temp.str());
    };
    void addAttributeDoubleDefault(const string& at, double val,int width,int prec ) {
	if ( val != 0.0 ) {
	    this->addAttributeDouble(at,val,width,prec);
	}
    };
    void addAttributeDoubleScientificWidth(const string& at, double val, int width ) {
        std::ostringstream temp;
	temp << std::setiosflags(std::ios::scientific)
	    	<< std::setprecision(width)
		<< val;
	this->addAttribute(at,temp.str());
    };
    void addAttributeFloatScientificWidth(const string& at, float val, int width ) {
	std::ostringstream temp;
	temp << std::setiosflags(std::ios::scientific)
	    	<< std::setprecision(width)
		<< val;
	this->addAttribute(at,temp.str());
    };
CL_LISPIFY_NAME("addAttributeDoubleScientific");
CL_DEFMETHOD     void addAttributeDoubleScientific(const string& at, double val ) {
	std::ostringstream temp;
	temp << std::setiosflags(std::ios::scientific) << val;
	this->addAttribute(at,temp.str());
    };
    void addAttributeFloatScientific(const string& at, float val ) {
	std::ostringstream temp;
	temp << std::setiosflags(std::ios::scientific) << val;
	this->addAttribute(at,temp.str());
    };
    void addAttributeDoubleScientificDefault(const string& at, double val ) {
	if ( val != 0.0 ) {
	    this->addAttributeDoubleScientific(at,val);
	}
    };
    void addAttributeFloatScientificDefault(const string& at, double val ) {
	if ( val != 0.0 ) {
	    this->addAttributeFloatScientific(at,val);
	}
    };

CL_LISPIFY_NAME("addAttributeInt");
CL_DEFMETHOD     void addAttributeInt(const string& at, int val) {
	std::ostringstream temp;
	temp << std::dec << val;
	this->addAttribute(at,temp.str());
    };
    void addAttributeIntDefault(const string& at, int val ) {
	if ( val != 0 ) {
	    this->addAttributeInt(at,val);
	}
    };

    void addAttributeHex(const string& at, int val) {
	std::ostringstream temp;
	temp << std::hex << val << std::dec;
	this->addAttribute(at,temp.str());
    };

    void addAttributeBool(const string& at, bool b ) {
	std::ostringstream temp;
	if ( b ) {
		temp << "yes";
	} else {
		temp << "no";
	}
	this->addAttribute(at,temp.str());
    };

    string asString();

    void throwErrorForChildren();
    void throwErrorForChildrenWithoutName(string name);

    void eraseAll();

		/*! Lock the file and parse the fileName
			if lock fails return false otherwise true
		*/
    bool lockAndParseFileName(const string fn);
		/*! Lock the file and write to the fileName
			if lock fails return false otherwise true
		*/
    bool lockAndWriteToFileName( string fileName );

    VectorQDomNodes	gatherSubNodesWithName(const string& name);


       QDomNode_O(string localName);
    QDomNode_O(string localName, string data);




    DEFAULT_CTOR_DTOR(QDomNode_O);
};




};

TRANSLATE(adapt::QDomNode_O);
#endif
