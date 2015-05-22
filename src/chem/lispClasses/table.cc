
#define	DEBUG_LEVEL_FULL

#include "table.h"
//#include "core/archiveNode.h"
#include "core/metaClass.h"
#include "core/builtInClass.h"
#include "core/hashTableEq.h"
//#include "core/archive.h"
#include "core/environment.h"

// last include is wrappers.h
#include "core/wrappers.h"


namespace chem
{


    void TableEntry_O::exposeCando(core::Lisp_sp e)
    {
	core::class_<TableEntry_O>()
	    .def("read",&TableEntry_O::read,"","","",false)
	    .def("write",&TableEntry_O::write,"","","",false)
	    .def("getField",&TableEntry_O::getField)
	    .def("putField",&TableEntry_O::putField)
	    .def("appendField",&TableEntry_O::appendField)
	    .def("insertField",&TableEntry_O::insertField)
	    .def("getTable",&TableEntry_O::getTableNotConst)
	    ;
    }





    void TableEntry_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(ChemPkg,TableEntry,"","",_lisp)
	    ;
#endif //]
    }




    TableEntry_sp TableEntry_O::createWithFields(uint numFields, Table_sp owner, core::Lisp_sp lisp)
    {_G();
	TableEntry_sp entry = TableEntry_O::create();
	IMPLEMENT_MEF(BF("handle setOwner"));
//	entry->setOwner(owner);
	entry->_Fields.resize(numFields,_Nil<core::T_O>());
	return entry;
    }

#if INIT_TO_FACTORIES

#define ARGS_TableEntry_O_make "()"
#define DECL_TableEntry_O_make ""
#define DOCS_TableEntry_O_make "make TableEntry"
  TableEntry_sp TableEntry_O::make()
  {_G();
    IMPLEMENT_ME();
  };

#else

    core::T_sp 	TableEntry_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {
	IMPLEMENT_ME();
    	// your stuff here
    }

#endif

    void	TableEntry_O::initialize()
    {
	this->Base::initialize();
    }

    Table_sp TableEntry_O::getTable() const
    {_G();
	IMPLEMENT_MEF(BF("Handle owner"));
//	return this->owner().as<Table_O>();
    }



#ifdef XML_ARCHIVE
    void	TableEntry_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	node->archiveVector0("fields",this->_Fields);
    }
#endif


    core::T_sp TableEntry_O::read(core::Symbol_sp fieldSymbol)
    {_G();
	Table_sp table = this->getTable();
	uint fieldIdx = table->indexOfField(fieldSymbol);
	return this->_Fields[fieldIdx];
    }

    core::T_sp TableEntry_O::getField(uint fieldIndex )
    {_G();
	ASSERT_lt(fieldIndex,this->_Fields.size());
	return this->_Fields[fieldIndex];
    }

    void TableEntry_O::write(core::Symbol_sp fieldSymbol, core::T_sp obj)
    {_G();
	Table_sp table = this->getTable();
	uint fieldIdx = table->indexOfField(fieldSymbol);
	this->_Fields[fieldIdx] = obj;
    }

    void TableEntry_O::putField(uint fieldIndex, core::T_sp obj)
    {_G();
	ASSERT_lt(fieldIndex,this->_Fields.size());
	this->_Fields[fieldIndex] = obj;
    }




    uint TableEntry_O::appendField(core::T_sp val)
    {_G();
	this->_Fields.push_back(val);
	return this->_Fields.size()-1;
    }



    uint TableEntry_O::insertField(uint position, core::T_sp val)
    {_G();
	ASSERT_lt(position,this->_Fields.size());
        gctools::Vec0<core::T_sp>::iterator it;
	it = this->_Fields.begin()+position;
	this->_Fields.insert(it,val);
	return position;
    }

#if 0
#ifdef RENDER
    candoBase::Render_sp TableEntry_O::rendered(core::Cons_sp opts)
    {_G();
	candoBase::DisplayList_sp dl;
	dl = candoBase::DisplayList_O::create();
        gctools::Vec0<core::T_sp>::iterator si;
	uint idx;
	for ( si=this->_Fields.begin(), idx=0;
	      si!=this->_Fields.end(); si++,idx++ )
	{
	    if ( (*si)->canRender() )
	    {
		candoBase::Render_sp rend = (*si)->rendered(opts);
		rend->setName(this->getTable()->fieldSymbol(idx));
		dl->append(rend);
	    } else 
	    {
		stringstream ss; 
		ss << this->getTable()->fieldSymbol(idx)->fullName() << " -> ";
		ss << (*si)->__repr__() << std::endl;
		candoBase::RenderInformation_sp info = candoBase::RenderInformation_O::create(ss.str(),_lisp);
		dl->append(info);
	    }
	}
	return dl;
    }
#endif
#endif

    string TableEntry_O::__repr__() const
    {_G();
	stringstream ss;
	uint idx=0;
	for ( gctools::Vec0<core::T_sp>::const_iterator fi=this->_Fields.begin(); fi!=this->_Fields.end(); fi++, idx++ )
	{
	    ss << this->getTable()->fieldSymbol(idx)->fullName() << " -> " << (*fi)->__repr__() << std::endl;
	}
	return ss.str();
    }




    void Table_O::exposeCando(core::Lisp_sp e)
    {
	core::class_<Table_O>()
	    .def("appendField",&Table_O::appendField)
	    .def("insertField",&Table_O::insertField)
	    .def("writeField",&Table_O::writeField)
	    .def("write",&Table_O::write,"","","",false)
	    .def("readField",&Table_O::readField)
	    .def("read",&Table_O::read,"","","",false)
	    .def("hasField",&Table_O::hasField)
	    .def("fieldClass",&Table_O::fieldClass)
	    .def("appendEntry",&Table_O::appendEntry)
	    .def("appendWrite",&Table_O::appendWrite)
	    .def("insertEntryBefore",&Table_O::insertEntryBefore)
	    .def("numberOfEntries",&Table_O::numberOfEntries)
	    .def("asCons",&Table_O::asCons,"","","",false)
	    .def("merge",&Table_O::merge)
	    .def("fieldSymbol",&Table_O::fieldSymbol)
	    .def("contentsAsString",&Table_O::contentsAsString)
	    ;
    }

    void Table_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(ChemPkg,Table,"","",_lisp)
	    ;
#endif //]
    }


#if INIT_TO_FACTORIES

#define ARGS_Table_O_make "()"
#define DECL_Table_O_make ""
#define DOCS_Table_O_make "make Table"
  Table_sp Table_O::make()
  {_G();
    IMPLEMENT_ME();
  };

#else

    core::T_sp 	Table_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {
	IMPLEMENT_ME();
    	// your stuff here
    }

#endif

    void	Table_O::initialize()
    {
	this->Base::initialize();
        this->_FieldIndices = core::HashTableEq_O::create_default();
    }

#ifdef XML_ARCHIVE
    void	Table_O::archiveBase(core::ArchiveP node)
    {_OF();
	this->Base::archiveBase(node);
	if ( node->saving() )
	{
	    core::VectorStrings symbolNames;
	    for ( gctools::Vec0<core::Symbol_sp>::iterator si=this->_FieldSymbols.begin(); si!=this->_FieldSymbols.end(); si++ )
	    {
		symbolNames.push_back((*si)->fullName());
	    }
	    node->archiveVectorStrings("symbolNames",symbolNames);
	    core::VectorStrings classNames;
	    for ( gctools::Vec0<core::Class_sp>::iterator ci=this->_FieldClasses.begin(); ci!=this->_FieldClasses.end(); ci++ )
	    {
		classNames.push_back((*ci)->instanceClassName());
	    }
	    node->archiveVectorStrings("classNames",classNames);
	} else // loading
	{
	    core::VectorStrings symbolNames;
	    node->archiveVectorStrings("symbolNames",symbolNames);
	    core::VectorStrings classNames;
	    node->archiveVectorStrings("classNames",classNames);
	    core::VectorStrings::iterator	si,ci;
	    this->_FieldSymbols.clear();
	    this->_FieldClasses.clear();
	    this->_FieldIndices.clear();
	    uint idx = 0;
	    for ( si=symbolNames.begin(),ci=classNames.begin(); si!=symbolNames.end(); si++, ci++, idx++ )
	    {
		core::Symbol_sp fieldSymbol = _lisp->intern(*si);
		this->_FieldSymbols.push_back(fieldSymbol);
		core::Class_sp mc = _lisp->classFromClassName(*ci);
		if ( mc.nilp() )
		{
		    core::Symbol_sp classSymbol = _lisp->intern(*ci);
		    IMPLEMENT_ME(); // Handle new way of creating classes
#if 0
		    mc = _lisp->defineStandardClass(classSymbol,_Nil<core::T_O>(),_Nil<core::Cons_O>());
#endif
		}
		this->_FieldClasses.push_back(mc);
		this->_FieldIndices[fieldSymbol] = idx;
	    }
	}
	node->archiveVector0("entries",this->_Entries);
	if ( node->loading() )
	{
	    this->_Entries.setOwnerOfAll(this->sharedThis<Table_O>());
	}
    }
#endif

    core::Symbol_sp Table_O::fieldSymbol(uint idx)
    {_G();
	ASSERT_lt(idx,this->_FieldSymbols.size());
	return this->_FieldSymbols[idx];
    }


    uint Table_O::indexOfField(core::Symbol_sp positionSymbol)
    {_OF();
        core::T_mv pi = this->_FieldIndices->gethash(positionSymbol);
	if ( pi.second().nilp() ) 
	{
	    SIMPLE_ERROR(BF("This table doesn't contain the field(%s)") % _rep_(positionSymbol));
	}
	return pi.as<core::Fixnum_O>()->get();
    }


    bool Table_O::hasField(core::Symbol_sp fieldSymbol)
    {_G();
	return (this->_FieldIndices->gethash(fieldSymbol).second().notnilp());
    }

    core::Class_sp Table_O::fieldClass(core::Symbol_sp fieldSymbol)
    {_G();
	if ( this->hasField(fieldSymbol) )
	{
	    return this->_FieldClasses[this->_FieldIndices[fieldSymbol]];
	}
	return _Nil<core::Class_O>();
    }


    void Table_O::appendField(core::Symbol_sp fieldSymbol, core::Class_sp fieldClass )
    {_G();
	if ( this->_FieldIndices.count(fieldSymbol) != 0 )
	{
	    SIMPLE_ERROR(BF("This table already contains a field named(%s)") % fieldSymbol->fullName() );
	}
	uint idx = this->_FieldSymbols.size();
	this->_FieldSymbols.push_back(fieldSymbol);
	this->_FieldClasses.push_back(fieldClass);
	for ( uint idx=0; idx<this->numberOfEntries(); idx++ )
	{
	    this->_Entries.get(idx)->appendField(_Nil<core::T_O>());
	}
	this->_FieldIndices[fieldSymbol] = idx;
    }


    void Table_O::insertField(core::Symbol_sp positionSymbol, core::Symbol_sp fieldSymbol, core::Class_sp fieldClass )
    {_G();
	uint pos = this->indexOfField(positionSymbol);
	this->_FieldSymbols.insert(this->_FieldSymbols.begin()+pos,fieldSymbol);
	this->_FieldClasses.insert(this->_FieldClasses.begin()+pos,fieldClass);
    	//
	// Slide all of the field indices at or above pos up one
	//
	map<core::Symbol_sp,uint>::iterator pi;
	for ( pi=this->_FieldIndices.begin(); pi!=this->_FieldIndices.end(); pi++ )
	{
	    if ( pos <= pi->second )
	    {
		pi->second++;
	    }
	}
	this->_FieldIndices[fieldSymbol] = pos;
    	//
	// Open up a column in the table for this
	// new field
	//
	for ( uint idx=0; idx<this->numberOfEntries(); idx++ )
	{
	    this->_Entries.get(idx)->insertField(pos,_Nil<core::T_O>());
	}
    }



    void Table_O::writeField(uint entryIdx, core::Symbol_sp field, core::T_sp val )
    {_OF();
	ASSERT_lt(entryIdx,this->_Entries.size());
	uint fieldIdx =  this->indexOfField(field);
	TableEntry_sp entry = this->_Entries[entryIdx];
	if ( !val->isAssignableToClass(this->_FieldClasses[fieldIdx]) )
	{
	    SIMPLE_ERROR(BF("You tried to put a value of class(%s) into a field of a table that expected a %s") % val->className() % this->_FieldClasses[fieldIdx]->instanceClassName() );
	}
	entry->putField(fieldIdx,val);
    }

    void Table_O::write(uint entryIdx, core::Cons_sp values )
    {_OF();
	ASSERT_lt(entryIdx,this->_Entries.size());
	while ( values.notnilp() )
	{
	    core::Symbol_sp sym = core::oCar(values).as<core::Symbol_O>();
	    values = values->cdr();
	    core::T_sp obj = core::oCar(values);
	    values = values->cdr();
	    this->writeField(entryIdx,sym,obj);
	}
    }


    TableEntry_sp Table_O::appendEntry()
    {_G();
	TableEntry_sp entry = TableEntry_O::createWithFields(this->_FieldSymbols.size(),this->sharedThis<Table_O>(),_lisp);
	this->_Entries.push_back(entry);
	return entry;
    }

    void Table_O::appendWrite(core::Cons_sp values)
    {_G();
	TableEntry_sp entry = this->appendEntry();
	while ( values.notnilp() )
	{
	    core::Symbol_sp sym = core::oCar(values).as<core::Symbol_O>();
	    values = values->cdr();
	    core::T_sp obj = core::oCar(values);
	    values = values->cdr();
	    entry->write(sym,obj);
	}
    }


    core::Cons_sp Table_O::asCons()
    {_G();
	return core::Cons_O::createFromVec0(this->_Entries);
    }


    TableEntry_sp Table_O::insertEntryBefore(uint position)
    {_G();
	ASSERT_lt(position,this->_Entries.size());
        gctools::Vec0<TableEntry_sp>::iterator it = this->_Entries.begin()+position;
	TableEntry_sp entry = TableEntry_O::createWithFields(this->_FieldSymbols.size(),this->sharedThis<Table_O>(),_lisp);
	this->_Entries.insert(it,entry);
	return entry;
    }

    uint Table_O::numberOfEntries()
    {_G();
	return this->_Entries.size();
    }


    TableEntry_sp Table_O::read(uint idx)
    {_G();
	ASSERT_lt(idx,this->_Entries.size());
	return this->_Entries.get(idx);
    }

    core::T_sp Table_O::readField(uint entryIdx, core::Symbol_sp field)
    {_G();
	ASSERT_lt(entryIdx,this->_Entries.size());
	uint fieldIdx =  this->indexOfField(field);
	TableEntry_sp entry = this->_Entries[entryIdx];
	return entry->getField(fieldIdx);
    }


    bool Table_O::fieldsMatch(Table_sp other)
    {_G();
	if ( this->_FieldSymbols.size() != other->_FieldSymbols.size() ) return false;
        gctools::Vec0<core::Symbol_sp>::iterator	ms,os;
        gctools::Vec0<core::Class_sp>::iterator mc,oc;
	for ( ms=this->_FieldSymbols.begin(), os = other->_FieldSymbols.begin(), 
		  mc = this->_FieldClasses.begin(), oc = other->_FieldClasses.begin();
	      ms < this->_FieldSymbols.end(); ms++, os++, mc++, oc++ )
	{
	    if ( (*ms) != (*os) ) return false;
	    if ( (*mc) != (*oc) ) return false;
	}
	return true;
    }

    void Table_O::merge(Table_sp other)
    {_G();
	if ( !this->fieldsMatch(other) )
	{
	    stringstream ss;
	    ss << "Tables to merge don't have matching fields" << std::endl;
	    ss << "this: " << this->__repr__() << std::endl;
	    ss << "other: " << other->__repr__() << std::endl;
	    SIMPLE_ERROR(BF("%s")%BF("%s")%ss.str());
	}
        gctools::Vec0<TableEntry_sp>::iterator it;
	for ( it=other->_Entries.begin(); it!=other->_Entries.end(); it++ )
	{
	    TableEntry_sp one = this->appendEntry();
	    for ( uint i=0; i<this->_FieldSymbols.size(); i++ )
	    {
		one->putField(i,(*it)->getField(i));
	    }
	}
    }




    string Table_O::__repr__() const
    {_G();
	stringstream sout;
	sout << (BF("Table number of fields: %d") % this->_FieldSymbols.size() ) << std::endl;
	sout << "Symbols: ";
	uint currentIdx=0;
	for ( gctools::Vec0<core::Symbol_sp>::const_iterator it = this->_FieldSymbols.begin(); it!=this->_FieldSymbols.end(); it++ )
	{
	    core::Symbol_sp sym = *it;
	    sout << " " << (*it)->fullName();
	    map<core::Symbol_sp,uint>::const_iterator fi = this->_FieldIndices.find(sym);
	    ASSERT(fi!=this->_FieldIndices.end());
	    uint fieldIdx = fi->second;
	    sout << "("<<fieldIdx<<")";
	    if ( fieldIdx != currentIdx )
	    {
		sout << "<-INDEX_WRONG!";
	    }
	    currentIdx++;
	}
	sout << std::endl;
	sout << "Classes: ";
	for ( gctools::Vec0<core::Class_sp>::const_iterator mi = this->_FieldClasses.begin(); mi!=this->_FieldClasses.end(); mi++ )
	{
	    sout << " " << (*mi)->instanceClassName();
	}
	sout << "Number of entries: " << this->_Entries.size() << std::endl;
	return sout.str();
    }

#if 0
    candoBase::Render_sp	Table_O::rendered(core::Cons_sp opts)
    {_G();
	candoBase::FrameList_sp	frames;
	frames = candoBase::FrameList_O::create();
	frames->setName(_lisp->internKeyword("entry"));
        gctools::Vec0<TableEntry_sp>::iterator si;
	for ( si=this->_Entries.begin(); si!=this->_Entries.end(); si++ )
	{
	    if ( (*si)->canRender() )
	    {
		candoBase::Render_sp rend = (*si)->rendered(opts);
		frames->append(rend);
	    }
	}
	return frames;
    }
#endif

    string	Table_O::contentsAsString()
    {_G();
	stringstream ss;
	uint idx = 0;
	for (gctools::Vec0<TableEntry_sp>::iterator ei=this->_Entries.begin(); ei!=this->_Entries.end(); ei++,idx++ )
	{
	    ss << "entry("<<idx<<"):: " << (*ei)->__repr__();
	}
	return ss.str();
    }


    EXPOSE_CLASS(chem,TableEntry_O);
    EXPOSE_CLASS(chem,Table_O);

};




