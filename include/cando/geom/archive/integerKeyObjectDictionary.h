       
       


#ifndef	IntegerKeyObjectDictionary_H //[
#define	IntegerKeyObjectDictionary_H


#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <cando/geom/geomPackage.fwd.h>
#include <clasp/core/holder.h>

namespace geom {


//SMART(IterateIntegerKeyObjectDictionary);

SMART(IntegerKeyObjectDictionary);
    class IntegerKeyObjectDictionary_O : public core::CxxObject_O
{
    LISP_CLASS(geom,GeomPkg,IntegerKeyObjectDictionary_O,"IntegerKeyObjectDictionary",core::CxxObject_O);
public:
#if defined(XML_ARCHIVE)
	void	archiveBase(ArchiveP node);
	void	archive(ArchiveP node);
#endif // defined(XML_ARCHIVE)
public:
    typedef gctools::smart_ptr<core::T_O>	valueType;
    typedef map<int,core::T_sp>::iterator	iterator;

protected:
    map<int,core::T_sp>	_Map;

public:
public:
#if 0
	KeyValueIterator begin_keyValue() {_G();
	    KeyValueIterator end_keyValue() {_G();
#endif

	string	__repr__() const;

	iterator begin() { return this->_Map.begin(); };
	iterator end() { return this->_Map.end(); };
public:

    virtual bool	contains( int c) { return this->_Map.count(c)>0; };
	void	clear() { this->_Map.clear(); };
    virtual valueType get(int name);
    valueType getOrDefault(int name, core::T_sp obj);

    		/*! Get and return the object and
		 * then remove it from the Dictionary
		 */
    valueType getAndRemove(int name);

    		/*! Get and return the object and
		 * then remove it from the Dictionary
		 * otherwise return the defaultObj
		 */
    valueType getAndRemoveOrDefault(int name, core::T_sp defaultObj );
		
    virtual void set(int s, valueType obj) { this->_Map[s] = obj; };


	void	put(int s, valueType obj) { this->set(s,obj);};


    virtual int	size() { return this->_Map.size(); };

		/*!Return a multiline string that summarizes the contents
		 * of this dictionary in human readable format
		 */
	string	summary();

    virtual string	allKeys();

    core::Cons_sp  valuesAsCons();
    core::Cons_sp  keysAsCons();
        core::Cons_sp keysAsList();

#ifdef	USEBOOSTPYTHON //[
	boost::python::list boost_keys();
#endif //]
public:

//	bool	canRender() { return true; };
//	Render_sp rendered(core::Cons_sp options);

	DEFAULT_CTOR_DTOR(IntegerKeyObjectDictionary_O);
	    };


};

TRANSLATE(geom::IntegerKeyObjectDictionary_O);

#endif //]
