       
#define	DEBUG_LEVEL_NONE
#include <clasp/core/common.h>
#include <clasp/core/str.h>
#include <cando/adapt/adapters.h>
#include <cando/chem/atomIndexer.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/adapt/stringList.h>
#include <clasp/core/wrappers.h>


namespace chem {

    REGISTER_CLASS(chem,AtomIndexer_O);
    REGISTER_CLASS(chem,MapOfMonomerNamesToAtomIndexers_O);




//
// Constructor
//

#ifdef XML_ARCHIVE
    void	AtomIndexer_O::archive(core::ArchiveP node)
{_G();
//    node->attribute("indexAll",this->_IndexAll);
    node->archiveVectorStrings("AtomNames",this->_Names);
}
#endif


int	AtomIndexer_O::addAtomName(core::Symbol_sp name)
{_G();
//    this->_IndexAll = false;
    if ( !this->containsAtomName(name) )
    {
        int id = this->_Names.size();
        this->_Names.push_back(name);
	LOG(BF("Added atom with name(%s) gave it id(%d)") % name.c_str() % id );
	return id;
    }
    SIMPLE_ERROR(BF("The atom name (%s) is already in the AtomIndexer") % _rep_(name));
}

void	AtomIndexer_O::setFromAtomNames(const string& names)
{_G();
vector<string>	ss;
vector<string>::iterator	it;
    LOG(BF("Setting from atom names: %s") % names.c_str()  );
    core::tokenize(names,ss," \n\t");
    this->_Names.clear();
    for ( it = ss.begin(); it!=ss.end(); it++ )
    {
	if ( (*it).size()!= 0 )
	{
          this->_Names.push_back(chemkw_intern(*it));
	}
    }
}


void	AtomIndexer_O::setFromStringList(adapt::StringList_sp sl)
{_G();
adapt::StringList_O::iterator		it;
    LOG(BF("Setting from atom names: %s") % sl->asString().c_str()  );
    for ( it = sl->begin(); it!=sl->end(); it++ )
    {
	if ( (*it).size()!= 0 )
	{
          this->_Names.push_back(chemkw_intern(*it));
	}
    }
}

void	AtomIndexer_O::appendConsOfTexts(core::List_sp sl)
{_G();
    LOG(BF("Appending atom names: %s") % sl->__repr__().c_str()  );
    for ( auto cur : sl ) {
      this->_Names.push_back(oCar(cur).as<core::Symbol_O>());
    }
}





bool	AtomIndexer_O::containsAtomName(core::Symbol_sp name)
{_G();
//    if ( this->_IndexAll ) return true;
    for ( auto mi=this->_Names.begin(); mi!=this->_Names.end(); mi++ ) {
        if ( (*mi)==name ) {
	    LOG(BF("Yes") );
	    return true;
	}
    }
    LOG(BF("No") );
    return false;
}

int AtomIndexer_O::indexForAtomName(core::Symbol_sp name)
{_G();
int				idx;
    idx = 0;
    for ( auto mi=this->_Names.begin(); mi!=this->_Names.end(); mi++ ) {
      if ( (*mi)==name ) return idx;
      idx++;
    }
    SIMPLE_ERROR(BF("Could not find name(%s) in AtomIndexer") % _rep_(name));
}


string AtomIndexer_O::asString() const
{_G();
  stringstream ss;
    ss.str("");
    if ( this->_Names.size() > 0 )
    {
        ss << *(this->_Names.begin());
	for ( auto mi=this->_Names.begin()+1; mi!=this->_Names.end(); mi++ )
	{
          ss << " " << _rep_(*mi);
	}
    }
    return ss.str();
}

string	AtomIndexer_O::description() const
{
    stringstream ss;
    ss << "( " << this->className() << " " << this->asString() << " )";
    return ss.str();
}


#if 0 //[
void	AtomIndexer_O::setIndexAll(bool b)
{_G();
    this->_IndexAll = b;
    if ( this->_IndexAll )
    {
        this->_Names.clear();
    }
}
#endif //]






//
// Constructor
//


#ifdef XML_ARCHIVE
    void	MapOfMonomerNamesToAtomIndexers_O::archive(core::ArchiveP node)
{_G();
    node->archiveSymbolMap("atomIndexers",this->_AtomIndexers);
}
#endif

void	MapOfMonomerNamesToAtomIndexers_O::set(core::Symbol_sp nm, AtomIndexer_O::smart_ptr ai)
{_G();
    this->_AtomIndexers.set(nm,ai);
}

    adapt::SymbolSet_sp MapOfMonomerNamesToAtomIndexers_O::getMonomerNamesAsSymbolSet() const
    {_OF();
	return adapt::SymbolSet_O::createFromKeysOfSymbolMap(this->_AtomIndexers);
    }

    bool MapOfMonomerNamesToAtomIndexers_O::recognizesMonomerName(core::Symbol_sp nm) const
    {_OF();
	return this->_AtomIndexers.contains(nm);
    }


    AtomIndexer_sp MapOfMonomerNamesToAtomIndexers_O::getAtomIndexerForMonomerName(core::Symbol_sp name) const
    {_OF();
	return this->_AtomIndexers.get(name);
    }

#if 0
    bool MapOfMonomerNamesToAtomIndexers_O::getAtomIndexerForMonomerName(core::Symbol_sp name) const
    {_OF();
	return this->_AtomIndexers.contains(name);
    }
#endif



void	MapOfMonomerNamesToAtomIndexers_O::eraseEntryForMonomer(core::Symbol_sp nm)
{_G();
    this->_AtomIndexers.remove(nm);
}

};


