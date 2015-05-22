#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <cando/candoBase/integerKeyObjectDictionary.h>
#include <clasp/core/numbers.h>
//#i n c l u d e "render.h"
#include <clasp/core/wrappers.h>

namespace candoBase {





    core::T_sp IntegerKeyObjectDictionary_O::get(int ikey)
{_OF();
    core::T_sp obj;
    if ( this->_Map.count(ikey) == 0 ) {
        stringstream ss;
	ss << "Could not find key(" << ikey << ")";
        SIMPLE_ERROR(BF(ss.str()));
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
        SIMPLE_ERROR(BF(ss.str()));
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
{_G();
DisplayList_sp	dl;
Render_sp		dlOne;
map<int,core::T_sp>::iterator	vi;
    dl = DisplayList_O::create();
    for ( vi=this->_Map.begin(); vi!=this->_Map.end(); vi++ )
    {
	LOG(BF("Looking to render %s") % vi->second->description().c_str()  );
	if ( vi->second->canRender() )
	{
	    LOG(BF("YES,Rendering object") );
	    dlOne = vi->second->rendered(opts);
	    Symbol_sp isymbol = _lisp->internKeyword((BF("i%d")%vi->first).str());
	    dlOne->setName(isymbol);
	    dl->append(dlOne);
	} else
	{
	    LOG(BF("NO can't render") );
	}
    }
    return dl;
}
#endif


string	IntegerKeyObjectDictionary_O::allKeys()
{_G();
map<int,core::T_sp>::iterator	vi;
stringstream		ss;

    for ( vi=this->_Map.begin(); vi!=this->_Map.end(); vi++ ) {
        ss << vi->first << " ";
    }
    return ss.str();
}

core::Cons_sp IntegerKeyObjectDictionary_O::keysAsCons()
{_G();
    map<int,core::T_sp>::iterator	vi;
    stringstream		ss;
    core::Cons_sp first = core::Cons_O::create(_Nil<core::T_O>(),_Nil<core::Cons_O>());
    core::Cons_sp cur = first;
    for ( vi=this->_Map.begin(); vi!=this->_Map.end(); vi++ ) {
	core::Cons_sp one = core::Cons_O::create(core::Fixnum_O::create(vi->first),_Nil<core::Cons_O>());
	cur->setCdr(one);
	cur=one;
    }
    return cCdr(first);
}

core::Cons_sp IntegerKeyObjectDictionary_O::valuesAsCons()
{_G();
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
{_G();
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
{_G();
    if ( node->saving() )
    { _BLOCK_TRACEF(BF("Saving IntegerKeyObjectDictionary - it contains %d objects") % this->_Map.size() );
	if ( this->_Map.size() != 0 )
	{
	    int i = 0;
	    stringstream suid;
	    IntegerKeyObjectDictionary_O::iterator oi;
	    for ( oi=this->_Map.begin(); oi!=this->_Map.end(); i++,oi++ )
	    {
	    	LOG(BF( "Archiving entry key(%d)")% oi->first );
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
    { _BLOCK_TRACEF(BF("Loading IntegerKeyObjectDictionary with uid(%d)") % node->getUniqueIdNumeric() );
	VectorNodes::iterator	ci;
	core::T_sp object;
	this->_Map.clear();
	for ( ci=node->begin_Children(); ci!=node->end_Children(); ci++ )
	{
	    (*ci)->setRecognized(true);
	    object = node->getArchive()->loadObjectDirectly((*ci));
	    ASSERTNOTNULL(object);
	    LOG(BF( "Adding to the IntegerKeyObjectDictionary key(%d)")% (*ci)->getUniqueIdNumeric() );
	    uint ikey = (*ci)->getUniqueIdNumeric();
	    this->_Map[ikey] = object;
	}
    }
}
#endif // defined(XML_ARCHIVE)


#if defined(XML_ARCHIVE)
void	IntegerKeyObjectDictionary_O::archive(ArchiveP node)
{_G();
    this->archiveBase(node);
}
#endif // defined(XML_ARCHIVE)


core::Cons_sp IntegerKeyObjectDictionary_O::keysAsList() 
{
map<int,core::T_sp>::iterator	vi;
core::Cons_sp list = _Nil<core::Cons_O>();
    for ( vi=this->_Map.begin(); vi!=this->_Map.end(); vi++ ) {
        list = core::Cons_O::create(core::Fixnum_O::create(vi->first),list);
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



    void IntegerKeyObjectDictionary_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<IntegerKeyObjectDictionary_O>()
	.def("put",&IntegerKeyObjectDictionary_O::put)
	.def("keys", &IntegerKeyObjectDictionary_O::keysAsList)
	.def("valuesAsCons", &IntegerKeyObjectDictionary_O::valuesAsCons)
	.def("allKeysAsString", &IntegerKeyObjectDictionary_O::allKeys)
	.def("contains",&IntegerKeyObjectDictionary_O::contains)
	.def("get",&IntegerKeyObjectDictionary_O::get)
	.def("getOrDefault",&IntegerKeyObjectDictionary_O::getOrDefault)
	.def("size",&IntegerKeyObjectDictionary_O::size)
    ;
}



    void IntegerKeyObjectDictionary_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(CandoBasePkg,IntegerKeyObjectDictionary,"","",_lisp)
	.def("__len__", &IntegerKeyObjectDictionary_O::size)
	.def("__getitem__", &IntegerKeyObjectDictionary_O::get)
	.def("keys", &IntegerKeyObjectDictionary_O::boost_keys)
	.def("set",&IntegerKeyObjectDictionary_O::set)
	.def("contains",&IntegerKeyObjectDictionary_O::contains)
	.def("get",&IntegerKeyObjectDictionary_O::get)
	.def("getOrDefault",&IntegerKeyObjectDictionary_O::getOrDefault)
	.def("size",&IntegerKeyObjectDictionary_O::size)
	;

#endif //]
}


    EXPOSE_CLASS(candoBase, IntegerKeyObjectDictionary_O);
};
