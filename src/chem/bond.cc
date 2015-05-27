#define	DEBUG_LEVEL_FULL
//
// (C) 2004 Christian E. Schafmeister
//



// Comment
#include <clasp/core/common.h>
#include <clasp/core/exceptions.h>
#include <cando/adapt/stringSet.h>
#include <cando/chem/bond.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <clasp/core/hashTableEq.h>
//#include "core/serialize.h"
//#include "core/archive.h"
//#include "core/archiveNode.h"
#include <cando/chem/symbolTable.h>
#include <cando/chem/constitutionAtoms.h>
#include <clasp/core/wrappers.h>



namespace chem 
{
    SYMBOL_EXPORT_SC_(ChemPkg,_PLUS_bondOrderToSymbolConverter_PLUS_);
    SYMBOL_EXPORT_SC_(ChemKwPkg,_singleBond);
    SYMBOL_EXPORT_SC_(ChemKwPkg,_doubleBond);
    SYMBOL_EXPORT_SC_(ChemKwPkg,_tripleBond);
    SYMBOL_EXPORT_SC_(ChemKwPkg,_aromaticBond);
    SYMBOL_EXPORT_SC_(ChemKwPkg,_hydrogenBond);
    SYMBOL_EXPORT_SC_(ChemKwPkg,_virtualBond);
    SYMBOL_EXPORT_SC_(ChemKwPkg,_dashedSingleBond);
    SYMBOL_EXPORT_SC_(ChemKwPkg,_dashedDoubleBond);


    string	XmlName_Bond = "bond";


    Bond_sp	Bond_O::create(Atom_sp from, Atom_sp to, BondOrder o)
    {_G();
	GC_ALLOCATE(Bond_O, bond );
	bond->_Atom1 = from;
	bond->_Atom2 = to;
	bond->setOrder(o);
	LOG(BF("created bond from=%s to=%s")
	    % bond->_Atom1->description()
	    % bond->_Atom2->description() );
	return bond;
    }


    void Bond_O::initialize()
    {_OF();
	this->Base::initialize();
	this->order = singleBond;
	LOG(BF("Setting fromAtom and toAtom to nil") ); //
	this->_Atom1 = _Nil<Atom_O>();
	this->_Atom2 = _Nil<Atom_O>();
    }

    Bond_O::Bond_O(const Bond_O& bb)  : core::T_O(bb)
    {_G();
	LOG(BF("Copying Bond_O") ); //
	this->order = bb.order;
	this->_Atom1 = bb._Atom1;
	this->_Atom2 = bb._Atom2;
	LOG(BF("copy _Atom1=%s _Atom2=%s")
	    % this->_Atom1->description()
	    % this->_Atom2->description() );
    }


    bool Bond_O::isAtom1(Atom_sp a) const
    {_OF();
	ASSERTNOTNULL(this->_Atom1);
	if (this->_Atom1.get()==a.get()) return true;
	ASSERTNOTNULL(this->_Atom2);
	if ( this->_Atom2.get()==a.get()) return false;
	SIMPLE_ERROR(BF("Atom[%s] is neither atom1[%s] or atom2[%s] of bond")
			   % a->__repr__() % this->_Atom1->__repr__() % this->_Atom2->__repr__() );
    }




    Bond_sp Bond_O::copyDontRedirectAtoms()
    {_OF();
	GC_COPY(Bond_O, bondNew , *this); // = RP_Copy<Bond_O>(this);
	return bondNew;
    }

    void Bond_O::redirectAtoms()
    {_OF();
	ASSERTNOTNULL(this->_Atom1);
	this->_Atom1 = this->_Atom1->getCopyAtom();
	ASSERTNOTNULL(this->_Atom2);
	this->_Atom2 = this->_Atom2->getCopyAtom();
	this->_Atom2->addBond(this->sharedThis<Bond_O>());
    }




//
//	isInterResidueBond
//
//	Return true if this is an inter-residue bond
//
    bool	Bond_O::isInterResidueBond()
    {_G();
	Matter_wp	wfromCont, wtoCont;
	Matter_sp	fromCont, toCont;
	Atom_sp a1 = this->_Atom1;
	ASSERTP(a1.objectp(), "Bond_O::isInterResidueBond from atom is undefined!");
	LOG(BF("from atom = %s") % a1->description() );
	wfromCont = a1->getResidueContainedBy();
	ASSERTNOTNULL(wfromCont);
	if ( wfromCont.nilp() )
	{
	    SIMPLE_ERROR(BF("From atom isnt in a residue"));
	}
	Atom_sp a2 = this->_Atom2;
	ASSERTNOTNULLP(a2, "Bond_O::isInterResidueBond to atom is undefined!");
	wtoCont = a2->getResidueContainedBy();
	ASSERTNOTNULL(wtoCont);
	if ( wtoCont.nilp() )
	{
	    SIMPLE_ERROR(BF("To atom isnt in a residue"));
	}
	fromCont = wfromCont;
	toCont = wtoCont;
	return ( fromCont != toCont );
    }




    core::NullTerminatedEnumAssociation bondOrderKeys[] = {
	{ "noBond", noBond },
	{ "single", singleBond },
	{ "double", doubleBond },
	{ "triple", tripleBond },
	{ "aromatic", aromaticBond },
	{ "hbond", hydrogenBond },
	{ "virtual", virtualBond },
	{ "dashedSingle", dashedSingleBond },
	{ "dashedDouble", dashedDoubleBond },
	{ "", -1 }
    };


/*!
 * Saving: The Archive getCurrentNode() will be a Node that
 * 	we will fill the node with attributes and children.
 * 	If it is then we just create a node that links to the existing one.
 * Loading: The Archive _CurrentNode will be pointing at the deep node for this object.
 * 	This function will extract the data.
 *
 */
#ifdef OLD_SERIALIZE
    void Bond_O::serialize(serialize::SNode node)
    {_G();
	LOG(BF("archive direction = %s") % (node->loading()?"loading":"saving") ); //
	node->attributeSymbolEnumHiddenConverter("order",this->order,_sym__PLUS_bondOrderToSymbolConverter_PLUS_);
	node->archiveSymbolMapIfDefined("properties",this->_Properties);
#ifdef	DEBUG_ON
	if ( !node->loading() )
	{
	    LOG(BF("before save this->_Atom1 = %s") % this->_Atom1->description() );
	    LOG(BF("before save this->_Atom2 = %s") % this->_Atom2->description() );
	}
#endif
	node->archiveWeakPointer<Atom_O>( "a1", this->_Atom1  );
	node->archiveWeakPointer<Atom_O>( "a2", this->_Atom2 );
#ifdef	DEBUG_ON
	if ( node->loading() )
	{
	    LOG(BF("after load this->_Atom1 = %s") % this->_Atom1->description() );
	    LOG(BF("after load this->_Atom2 = %s") % this->_Atom2->description() );
	}
#endif
    };
#endif



/*!
 * Saving: The Archive getCurrentNode() will be a Node that
 * 	we will fill the node with attributes and children.
 * 	If it is then we just create a node that links to the existing one.
 * Loading: The Archive _CurrentNode will be pointing at the deep node for this object.
 * 	This function will extract the data.
 *
 */
#ifdef XML_ARCHIVE
    void	Bond_O::archiveBase(core::ArchiveP node)
    {_G();
	LOG(BF("archive direction = %s") % (node->loading()?"loading":"saving") ); //
	node->attributeSymbolEnumHiddenConverter("order",this->order,_sym__PLUS_bondOrderToSymbolConverter_PLUS_);
	node->archiveSymbolMapIfDefined("properties",this->_Properties);
#ifdef	DEBUG_ON
	if ( !node->loading() )
	{
	    LOG(BF("before save this->_Atom1 = %s") % this->_Atom1->description() );
	    LOG(BF("before save this->_Atom2 = %s") % this->_Atom2->description() );
	}
#endif
	node->archiveWeakPointer<Atom_O>( "a1", this->_Atom1  );
	node->archiveWeakPointer<Atom_O>( "a2", this->_Atom2 );
#ifdef	DEBUG_ON
	if ( node->loading() )
	{
	    LOG(BF("after load this->_Atom1 = %s") % this->_Atom1->description() );
	    LOG(BF("after load this->_Atom2 = %s") % this->_Atom2->description() );
	}
#endif
    };
#endif






    string Bond_O::propertiesAsString() const
    {_OF();
        return core::alist_asString(this->_Properties);
    }
	


    void Bond_O::clearProperty(core::Symbol_sp prop)
    {_G();
        this->_Properties = core::alist_erase(this->_Properties,prop);
    }

    void Bond_O::pushProperty(core::Symbol_sp prop, core::T_sp val)
    {_G();
        this->_Properties = core::alist_push(this->_Properties,prop,val);
    }

    void Bond_O::pushPropertyTrue(core::Symbol_sp prop)
    {_G();
        this->_Properties = core::alist_push(this->_Properties,prop,_lisp->_true());
    }

    core::T_sp Bond_O::getProperty(core::Symbol_sp prop, core::T_sp defval)
    {_G();
        core::List_sp match = core::alist_get(this->_Properties,prop);
        if ( match.nilp() ) {
            return defval;
        }
        return oCdr(match);
    }

    bool Bond_O::hasProperty(core::Symbol_sp prop)
    {_G();
        core::List_sp match = core::alist_get(this->_Properties,prop);
        return match.notnilp();
    }



    Atom_sp Bond_O::getOtherAtom(Atom_sp atom) const
    {_OF();
	ASSERTNOTNULL(this->_Atom1);
	ASSERTNOTNULL(this->_Atom2);
	if ( atom.get() == this->_Atom1.get() ) return this->_Atom2;
	if ( atom.get() == this->_Atom2.get() ) return this->_Atom1;
	SIMPLE_ERROR(BF("Illegal getOtherAtom because from atom[%s] is not part of bond: %s") % atom->description() % this->description() );
    }


    bool Bond_O::isBondBetween(Atom_sp a, Atom_sp b) const
    {_OF();
	ASSERTNOTNULL(this->_Atom1);
	ASSERTNOTNULL(this->_Atom2);
	if ( this->_Atom1.get() == a.get() )
	{
	    if ( this->_Atom2.get() == b.get() ) return true;
	}
	if ( this->_Atom2.get() == a.get() )
	{
	    if ( this->_Atom1.get() == b.get() ) return true;
	}
	return false;
    }


    string  Bond_O::description() const
    {_OF();
	stringstream    ss;
	ss << "Bond(";
	Atom_wp wa1 = this->_Atom1;
	ANN(wa1);
	ss << wa1->description();
	ss << bondOrderToChar(this->order);
	Atom_wp wa2 = this->_Atom2;
	ANN(wa2);
	ss << wa2->description();
	ss << ")@"<<std::hex<<this<<std::dec;
	return ss.str();
    }

    string  Bond_O::describeOther(Atom_sp from) const
    {_OF();
	ASSERTF(from.get()==this->_Atom1.get()||from.get()==this->_Atom2.get(),BF("describeFrom failed because from atom[%s] is not part of bond: %s") % from->description() % this->description() );
	stringstream    ss;
	ss << "Bond(";
	ss << from->description();
	ss << bondOrderToChar(this->order);
	Atom_sp a2 = this->getOtherAtom(from);
	ss << a2->description();
	ss << ")";
	return ss.str();
    }



    bool	Bond_O::invalid(Atom_sp a)
    {_OF();
	ASSERTNOTNULL(this->_Atom1);
	ASSERTNOTNULL(this->_Atom2);
	if ( this->_Atom1.nilp() ) return true;
	if ( this->_Atom2.nilp() ) return true;
	return false;
    }

    void	Bond_O::failIfInvalid(Atom_sp a)
    {_G();
	if ( this->invalid(a) ) {
	    SIMPLE_ERROR(BF("INVALID %s")%this->description());
	}
    }



    void Bond_O::redirectToAtomCopies()
    {_G();
	LOG(BF("Redirecting bond@%p") % this ); //
	Atom_sp fa =  this->_Atom1;
	Atom_sp ta = this->_Atom2;
#ifdef	DEBUG_ON
	if ( !fa.objectp() ) {
	    SIMPLE_ERROR(BF("redirectToAtomCopies _Atom1 is NULL"));
	}
	LOG(BF("  original from atom@%p") % fa.get() ); //
	if ( !ta.objectp() ) {
	    SIMPLE_ERROR(BF("redirectToAtomCopies _Atom2 is NULL"));
	}
	LOG(BF("  original   to atom@%p") % ta.get() ); //
#endif
	Atom_wp fc = fa->getCopyAtom();
	ASSERTNOTNULL(fc);
	LOG(BF("    new from %s") % fc->description() );
	Atom_wp tc = ta->getCopyAtom();
	ASSERTNOTNULL(tc);
	LOG(BF("    new   to %s") % tc->description() );
	this->_Atom1 = fc;
	this->_Atom2 = tc;
    }

    void	Bond_O::imposeYourself()
    {_G();
	Atom_sp a1 = this->_Atom1;
	Atom_sp a2 = this->_Atom2;
	a1->_addExistingBond(this->sharedThis<Bond_O>());
	a2->_addExistingBond(this->sharedThis<Bond_O>());
    }



    ConstitutionBond_sp Bond_O::asConstitutionBond(Atom_sp from, MapAtomsToConstitutionAtomIndex0N atomMap)
    {_OF();
	Atom_sp to = this->getOtherAtom(from);
        MapAtomsToConstitutionAtomIndex0N::iterator it = atomMap.find(to);
	if ( it == atomMap.end() )
	{
	    SIMPLE_ERROR(BF("Could not find atom[%s] in atomMap") % _rep_(to));
	}
        ConstitutionAtomIndex0N index = it->second;
	ConstitutionBond_sp cb = ConstitutionBond_O::create(this->getOrder(),index);
	return cb;
    }




    string	Bond_O::getOrderAsString()
    {_G();
	return bondOrderToString(this->order);
    }

    void	Bond_O::joinYourAtoms()
    {_OF();
	IMPLEMENT_ME(); // handle new way of storing bonds
#if 0
	ASSERTNOTNULLP(this->fromAtom,"Bond_O::joinYourAtoms from atom is undefined!");
	ASSERTNOTNULLP(this->toAtom,"Bond_O::joinYourAtoms to atom is undefined!");
	this->fromAtom->bondTo(this->toAtom,this->order);
#endif
    }




//typedef	enum 	BondOrderEnum {
//		noBond=NoBond,
//		singleBond=SingleBond,
//    		doubleBond = DoubleBond,
//		tripleBond = TripleBond,
//		aromaticBond = AromaticBond } BondOrder;

    string	bondOrderToString(BondOrder bo) {

	switch (bo) {
	case noBond:
	    return "none";
	    break;
	case singleBond:
	    return "single";
	    break;
	case doubleBond:
	    return "double";
	    break;
	case tripleBond:
	    return "triple";
	    break;
	case aromaticBond:
	    return "aromatic";
	    break;
	case hydrogenBond:
	    return "HBond";
	    break;
	default:
	    return "-unknown-";
	}
	printf( "Illegal bond order: %d\n", (int)(bo));
	exit(1);
    }



    BondOrder	stringToBondOrder(string bos) {
	if ( bos=="none" ) return noBond;
	if ( bos=="single" ) return singleBond;
	if ( bos=="double" ) return doubleBond;
	if ( bos=="triple" ) return tripleBond;
	if ( bos=="aromatic" ) return aromaticBond;
	if ( bos=="HBond" ) return hydrogenBond;
	printf( "Illegal bond order: %s\n", bos.c_str() );
	exit(1);
    }


    char bondOrderToChar(BondOrder bo) 
    {
	switch (bo) {
	case noBond:
	    return '?';
	    break;
	case singleBond:
	    return '-';
	    break;
	case doubleBond:
	    return '=';
	    break;
	case tripleBond:
	    return '#';
	    break;
	case aromaticBond:
	    return '@';
	    break;
	case dashedSingleBond:
	    return '!';
	case dashedDoubleBond:
	    return '%';
	    break;
	default:
	    THROW_HARD_ERROR(BF("No char available for bond"));
	}
	printf( "Illegal bond order: %d\n", (int)(bo));
	exit(1);
    }


    BondList_O::BondList_O(const BondList_O& orig) : Base(orig), _Bonds(orig._Bonds)
    {_OF();
    }



    void	BondList_O::initialize()
    {
	this->Base::initialize();
	this->_Bonds.clear();
    };




    void	BondList_O::imposeYourself()
    { _G();
        gctools::Vec0<Bond_sp>::iterator	bi;
	for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ ) {
	    (*bi)->imposeYourself();
	}
    }



    void	BondList_O::removeBond(Bond_sp b)
    { _G();
	BondList_O::iterator	bi;
	for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ ) {
	    if ( (*bi) == b ) {
		this->_Bonds.erase(bi);
		return;
	    }
	}
	// If the bond wasn't found then just return
    }

    void	BondList_O::removeBondBetween(Atom_sp a,Atom_sp b)
    {_OF();
	BondList_O::iterator	bi;
	for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ )
	{
	    if ( (*bi)->isBondBetween(a,b) )
	    {
		this->_Bonds.erase(bi);
		return;
	    }
	}
	SIMPLE_ERROR(BF("Could not find bond between %s and %s in bondList: %s") % a->description() % b->description() % this->description() );
	// If the bond wasn't found then just return
    }


    void	BondList_O::addBond(Bond_sp b)
    {_G();
	this->_Bonds.push_back(b);
    };

#if 0
/*! Shallow copy of a bond list */
    BondList_sp	BondList_O::copy()
    {
	BondList_O::iterator	bi;
	GC_ALLOCATE(BondList_O, bl );
	for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ )
	{
	    bl->addBond(*bi);
	}
	return bl;
    }
#endif

    string	BondList_O::description() const
    {
	BondList_O::const_iterator	bi;
	stringstream ss;
	ss << "BondList[[[" << std::endl;
	for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ ) {
	    ss << (*bi)->description() << "," << std::endl;
	}
	ss << "]]]" << std::endl;
	return ss.str();
    }

    string	BondList_O::describeOthers(Atom_sp from) const
    {
	BondList_O::const_iterator	bi;
	stringstream ss;
	ss << "BondList[[[" << std::endl;
	for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ ) {
	    ss << (*bi)->describeOther(from) << "," << std::endl;
	}
	ss << "]]]" << std::endl;
	return ss.str();
    }



#ifdef OLD_SERIALIZE
    void BondList_O::serialize(serialize::SNode node)
    {_G();
	node->orderedVector0(this->_Bonds);
    }
#endif





#ifdef XML_ARCHIVE
    void	BondList_O::archiveBase(core::ArchiveP node)
    {_G();
	BondList_O::iterator	ib;
	if ( node->saving() ) {
	    LOG(BF("About to save BondList with %d members") % this->size()  ); //
	    for ( ib=this->_Bonds.begin(); ib!=this->_Bonds.end(); ib++ ) {
		_BLOCK_TRACE("Archiving one bond of bondList");
		node->archiveSaveObjectAsChildAssignAutoUniqueId<Bond_O>(*ib);
		LOG(BF("After archiveObject just print") ); //
	    }
	} else {
	    LOG(BF("About to load BondList") ); //
	    this->_Bonds.clear();
	    core::VectorNodes::iterator	ic;
	    Bond_sp			b,bEmpty;
	    for ( ic=node->begin_Children(); ic!=node->end_Children(); ic++ ) {
		(*ic)->archiveLoadObjectDirectly<Bond_O>(b);
		ASSERTNOTNULL(b);
		this->_Bonds.push_back(b);
	    }
	}
    }
#endif


    void Bond_O::lisp_initGlobals(core::Lisp_sp lisp)
    {_G();
	SYMBOL_EXPORT_SC_(ChemKwPkg,noBond);
	SYMBOL_EXPORT_SC_(ChemKwPkg,singleBond);
	SYMBOL_EXPORT_SC_(ChemKwPkg,doubleBond);
	SYMBOL_EXPORT_SC_(ChemKwPkg,tripleBond);
	SYMBOL_EXPORT_SC_(ChemKwPkg,aromaticBond);
	SYMBOL_EXPORT_SC_(ChemKwPkg,hydrogenBond);
	SYMBOL_EXPORT_SC_(ChemKwPkg,virtualBond);
	SYMBOL_EXPORT_SC_(ChemKwPkg,dashedSingleBond);
	SYMBOL_EXPORT_SC_(ChemKwPkg,dashedDoubleBond);
        FIX_ME();
#if 0
	core::enum_<BondOrder> en(_sym__PLUS_bondOrderToSymbolConverter_PLUS_,"BondOrder");
	en.value(chemkw::_sym_noBond, _lisp->internKeyword("none"), singleBond );
	en.value(chemkw::_sym_singleBond, _lisp->internKeyword("s"), singleBond );
	en.value(chemkw::_sym_doubleBond, _lisp->internKeyword("d"), doubleBond );
	en.value(chemkw::_sym_tripleBond, _lisp->internKeyword("t"), tripleBond );
	en.value(chemkw::_sym_aromaticBond, _lisp->internKeyword("a"), aromaticBond );
	en.value(_sym_hydrogenBond, _lisp->internKeyword("h"), hydrogenBond );
	en.value(_sym_virtualBond,_lisp->internKeyword("v"),virtualBond);
	en.value(_sym_dashedSingleBond,_lisp->internKeyword("ds"),dashedSingleBond);
	en.value(_sym_dashedDoubleBond,_lisp->internKeyword("dd"),dashedDoubleBond);
#endif
    }	


    void Bond_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<Bond_O>()
	    .def("getAtom1",&Bond_O::getAtom1)
	    .def("getAtom2",&Bond_O::getAtom2)
	    .def("getOrder",&Bond_O::getOrder)
	    .def("getOrderAsString",&Bond_O::getOrderAsString)
	    .def("setOrder",&Bond_O::setOrderFromInt)
	    .def("getOtherAtom",&Bond_O::getOtherAtom)
	    .def("hasProperty",&Bond_O::hasProperty)
	    .def("pushProperty",&Bond_O::pushProperty)
	    .def("getProperty",&Bond_O::getProperty)
	    ;
    }

    void Bond_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,Bond,"","",_lisp)
	    .def("getAtom1",&Bond_O::getAtom1)
	    .def("getAtom2",&Bond_O::getAtom2)
	    .def("getOrder",&Bond_O::getOrder)
	    .def("getOrderAsString",&Bond_O::getOrderAsString)
	    .def("setOrder",&Bond_O::setOrderFromInt)
	    .def("getOtherAtom",&Bond_O::getOtherAtom)
	    .def("hasProperty",&Bond_O::hasProperty)
	    .def("setProperty",&Bond_O::setProperty)
	    .def("getProperty",&Bond_O::getProperty)
	    .def("getPropertyOrDefault",&Bond_O::getPropertyOrDefault)
	    ;
#endif
    }



    void BondList_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<BondList_O>();
    }

    void BondList_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,BondList,"","",_lisp)
	    ;
#endif
    }


    EXPOSE_CLASS_AND_GLOBALS(chem,Bond_O);
    EXPOSE_CLASS(chem,BondList_O);


}; // namespace chem


