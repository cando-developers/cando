/*
    File: integerKeyObjectDictionary.cc
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
#define	DEBUG_LEVEL_NONE

#include <clasp/core/common.h>
#include <cando/geom/integerKeyObjectDictionary.h>
#include <clasp/core/numbers.h>
//#i n c l u d e "render.h"
#include <clasp/core/wrappers.h>

namespace geom {





    core::T_sp IntegerKeyObjectDictionary_O::get(int ikey)
{_OF();
    core::T_sp obj;
    if ( this->_Map.count(ikey) == 0 ) {
        stringstream ss;
	ss << "Could not find key(" << ikey << ")";
        SIMPLE_ERROR((ss.str()));
    }
    return this->_Map[ikey];
}

core::T_sp IntegerKeyObjectDictionary_O::getOrDefault(int ikey, core::T_sp defVal)
{
    if ( this->_Map.count(ikey) == 0 ) 
    {
	return defVal;
    }
    return this->get(ikey);
}


string	IntegerKeyObjectDictionary_O::__repr__() const
{
    stringstream ss;
    map<int,core::T_sp>::const_iterator	vi;
    for ( vi=this->_Map.begin(); vi!=this->_Map.end(); vi++ )
    {
	ss << "Entry{ " << vi->first << ": ";
	ss << vi->second->_instanceClass()->classNameAsString();
	if ( vi->second.isA<core::Number_O>() )
	{
	    ss <<":" <<_rep_(vi->second);
	}
	ss << std::endl;
    }
    return ss.str();
}



core::T_sp IntegerKeyObjectDictionary_O::getAndRemove(int ikey)
{_OF();
    core::T_sp	obj;
    if ( this->_Map.count(ikey) == 0 ) 
    {
        stringstream ss;
	ss << "Could not find key(" << ikey << ")";
        SIMPLE_ERROR((ss.str()));
    }
    obj = this->_Map[ikey];
    this->_Map.erase(ikey);
    return obj;
}

core::T_sp IntegerKeyObjectDictionary_O::getAndRemoveOrDefault(int ikey, core::T_sp dobj)
{
    core::T_sp	obj;
    if ( this->_Map.count(ikey) == 0 ) {
	return dobj;
    }
    obj = this->_Map[ikey];
    this->_Map.erase(ikey);
    return obj;
}

#if 0
Render_sp	IntegerKeyObjectDictionary_O::rendered(core::Cons_sp opts)
{
DisplayList_sp	dl;
Render_sp		dlOne;
map<int,core::T_sp>::iterator	vi;
    dl = DisplayList_O::create();
    for ( vi=this->_Map.begin(); vi!=this->_Map.end(); vi++ )
    {
	LOG("Looking to render %s" , vi->second->description().c_str()  );
	if ( vi->second->canRender() )
	{
	    LOG("YES,Rendering object" );
	    dlOne = vi->second->rendered(opts);
	    Symbol_sp isymbol = _lisp->internKeyword((BF("i%d")%vi->first).str());
	    dlOne->setName(isymbol);
	    dl->append(dlOne);
	} else
	{
	    LOG("NO can't render" );
	}
    }
    return dl;
}
#endif


string	IntegerKeyObjectDictionary_O::allKeys()
{
map<int,core::T_sp>::iterator	vi;
stringstream		ss;

    for ( vi=this->_Map.begin(); vi!=this->_Map.end(); vi++ ) {
        ss << vi->first << " ";
    }
    return ss.str();
}

core::Cons_sp IntegerKeyObjectDictionary_O::keysAsCons()
{
    map<int,core::T_sp>::iterator	vi;
    stringstream		ss;
    core::Cons_sp first = core::Cons_O::create(_Nil<core::T_O>(),_Nil<core::Cons_O>());
    core::Cons_sp cur = first;
    for ( vi=this->_Map.begin(); vi!=this->_Map.end(); vi++ ) {
	core::Cons_sp one = core::Cons_O::create(core::clasp_make_fixnum(vi->first),_Nil<core::Cons_O>());
	cur->setCdr(one);
	cur=one;
    }
    return cCdr(first);
}

core::Cons_sp IntegerKeyObjectDictionary_O::valuesAsCons()
{
    map<int,core::T_sp>::iterator	vi;
    stringstream		ss;
    core::Cons_sp first = core::Cons_O::create(_Nil<core::T_O>(), _Nil<core::Cons_O>());
    core::Cons_sp cur = first;
    for ( vi=this->_Map.begin(); vi!=this->_Map.end(); vi++ ) {
	core::Cons_sp one = core::Cons_O::create(vi->second,_Nil<core::Cons_O>());
	cur->setCdr(one);
	cur=one;
    }
    return cCdr(first);
}

string	IntegerKeyObjectDictionary_O::summary()
{
map<int,core::T_sp>::iterator	vi;
stringstream		ss;

    for ( vi=this->_Map.begin(); vi!=this->_Map.end(); vi++ ) {
        ss << vi->first << ":";
	ss << vi->second->description() << std::endl;
    }
    return ss.str();
}





#if defined(XML_ARCHIVE)
void	IntegerKeyObjectDictionary_O::archiveBase(ArchiveP node)
{
    if ( node->saving() )
    { 
	if ( this->_Map.size() != 0 )
	{
	    int i = 0;
	    stringstream suid;
	    IntegerKeyObjectDictionary_O::iterator oi;
	    for ( oi=this->_Map.begin(); oi!=this->_Map.end(); i++,oi++ )
	    {
	    	LOG( "Archiving entry key(%d)"% oi->first );
		suid.str("");
		suid << i;
		boost::format fmt = boost::format("%d");
		fmt % oi->first;
		TRY_BOOST_FORMAT_STRING(fmt,fmt_str);
		string key = fmt_str; // (boost::format("%d") % oi->first).str();
		node->archiveObject(key,oi->second);
	    }
	}
    } else
    { 
	VectorNodes::iterator	ci;
	core::T_sp object;
	this->_Map.clear();
	for ( ci=node->begin_Children(); ci!=node->end_Children(); ci++ )
	{
	    (*ci)->setRecognized(true);
	    object = node->getArchive()->loadObjectDirectly((*ci));
	    ASSERTNOTNULL(object);
	    LOG( "Adding to the IntegerKeyObjectDictionary key(%d)"% (*ci)->getUniqueIdNumeric() );
	    uint ikey = (*ci)->getUniqueIdNumeric();
	    this->_Map[ikey] = object;
	}
    }
}
#endif // defined(XML_ARCHIVE)


#if defined(XML_ARCHIVE)
void	IntegerKeyObjectDictionary_O::archive(ArchiveP node)
{
    this->archiveBase(node);
}
#endif // defined(XML_ARCHIVE)


core::Cons_sp IntegerKeyObjectDictionary_O::keysAsList() 
{
map<int,core::T_sp>::iterator	vi;
core::Cons_sp list = _Nil<core::Cons_O>();
    for ( vi=this->_Map.begin(); vi!=this->_Map.end(); vi++ ) {
        list = core::Cons_O::create(core::clasp_make_fixnum(vi->first),list);
    }
    return list;
}



#ifdef	USEBOOSTPYTHON //[
boost::python::list IntegerKeyObjectDictionary_O::boost_keys() {
boost::python::list	res;
map<int,core::T_sp>::iterator	vi;
    for ( vi=this->_Map.begin(); vi!=this->_Map.end(); vi++ ) {
	res.append(boost::python::str(vi->first));
    }
    return res;
}
#endif //]








};
