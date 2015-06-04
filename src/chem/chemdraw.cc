#define	DEBUG_LEVEL_FULL

#include <string.h>
#include <clasp/core/common.h>
#include <clasp/core/str.h>
#include <cando/adapt/adapters.h>
#include <cando/adapt/stringSet.h>
#include <clasp/core/binder.h>
#include <cando/chem/chemdraw.h>
#include <cando/adapt/quickDom.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/bond.h>
#include <cando/chem/residue.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/aggregate.h>
#include <clasp/core/environment.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/lispReader.h>
#include <cando/chem/symbolTable.h>
#include <cando/chem/spanningLoop.h>
#include <clasp/core/reader.h>
#include <cando/chem/cipPrioritizer.h>
#include <cando/chem/constitutionAtoms.h>
#include <clasp/core/wrappers.h>


namespace chem
{

    SYMBOL_EXPORT_SC_(ChemPkg,name);



    void	CDNode_O::initialize()
    {
	this->Base::initialize();
	this->_Atom = _Nil<Atom_O>();
    }

    string	CDNode_O::_extractLabel(adapt::QDomNode_sp node)
    {
	string name;
	if ( node->hasChildrenWithName("t") )
	{
	    adapt::QDomNode_sp text = node->childWithName("t");
	    if ( text->hasChildrenWithName("s") )
	    {
		adapt::QDomNode_sp xmls = text->childWithName("s");
		string name = xmls->getData();
		return name;
	    }
	}
	return "C_"+node->getAttributeString("id");
    }


    void CDNode_O::getParsedLabel(string& name, int& ionization) const
    {
      vector<string> parts = core::split(this->getLabel(),"/");
      name = parts[0];
      ionization = 0;
      string ionstr = "";
      if ( parts.size() == 2 )
      {
        ionstr = parts[1];
        for ( size_t i=0; i<ionstr.size(); i++ )
        {
          if ( ionstr[i] == '+' ) ionization++;
          if ( ionstr[i] == '-' ) ionization--;
        }
      }
      LOG(BF("parsed[%s] into name[%s] ionization[%d]") % this->getLabel() % name % ionization );
    }
    void	CDNode_O::parseFromXml(adapt::QDomNode_sp xml)
    {_OF();
	this->_Id = xml->getAttributeInt("id");
	this->_Label = this->_extractLabel(xml);
	LOG(BF("Parsing CDNode with label: %s") % this->_Label);
	this->_StereochemistryType = undefinedCenter;
	this->_Configuration = undefinedConfiguration;
	if ( xml->hasAttribute("Geometry") )
	{
	    this->_StereochemistryType = chiralCenter;
	    string geometry = xml->getAttributeString("Geometry");
	    if ( xml->hasAttribute("AS") )
	    {
		string as = xml->getAttributeString("AS");
		if ( as == "S" )
		{
		    LOG(BF("Atom has geometry AS=%s") % as );
		    this->_Configuration = S_Configuration;
		} else if ( as == "R" )
		{
		    LOG(BF("Atom has geometry AS=%s") % as );
		    this->_Configuration = R_Configuration;
		} else
		{
		    LOG(BF("Could not interpret geometry AS[%s]") % as );
		}
	    }
	}
    }


    void	CDNode_O::unidirectionalBondTo(CDNode_sp target, CDBondOrder o)
    {
	pair< CDNode_wp, CDBondOrder> bd;
	bd.first = target;
	bd.second = o;
	this->_Neighbors.push_back(bd);
    }

    void	CDNode_O::bondTo(CDNode_sp other, CDBondOrder o )
    {
	this->bondTo(other,o);
	other->bondTo(this->sharedThis<CDNode_O>(),o);
    }


    void	CDBond_O::initialize()
    {
	this->Base::initialize();
	this->_BeginNode = _Nil<CDNode_O>();
	this->_EndNode = _Nil<CDNode_O>();
    }

    string	CDBond_O::getOrderAsString()
    {
	switch ( this->_Order )
	{
	case singleCDBond: return "single";
	case doubleCDBond: return "double";
	case tripleCDBond: return "triple";
	case dativeCDBond: return "dative";
	case singleDashCDBond: return "singleDash";
	case doubleDashCDBond: return "doubleDash";
	case tripleDashCDBond: return "tripleDash";
	case hashCDBond: return "hash";
	case unknownCDBond: return "unknown";
	}
	return "unknown";
    }

    BondOrder CDBond_O::getOrderAsBondOrder()
    {
	switch ( this->_Order )
	{
	case singleDashCDBond: return dashedSingleBond;
	case singleCDBond: return singleBond;
	case doubleDashCDBond: return dashedDoubleBond;
	case doubleCDBond: return doubleBond;
	case tripleCDBond: return tripleBond;
	case dativeCDBond: return noBond;
	case hashCDBond: return noBond;
	case unknownCDBond: return noBond;
	case tripleDashCDBond: return noBond;
	}
	return noBond;
    }


    void	CDBond_O::parseFromXml(adapt::QDomNode_sp xml)
    {
	this->_IdBegin = xml->getAttributeInt("B");
	this->_IdEnd = xml->getAttributeInt("E");
	string order = xml->getAttributeStringDefault("Order","1");
	string display = xml->getAttributeStringDefault("Display","");
	if ( order == "1" )
	{
	    if ( display == "Dash" )
	    {
		this->_Order = singleDashCDBond;
	    } else if ( display == "Hash" )
	    {
		this->_Order = hashCDBond;
	    } else
	    {
		this->_Order = singleCDBond;
	    }
	} else if ( order == "1.5") // Chemdraw saves a double dashed bond as Order=1.5
	{
	    if ( display == "Dash" ) this->_Order = doubleDashCDBond;
	} else if ( order == "2" )
	{
	    this->_Order = doubleCDBond;
	} else if ( order == "3" ) 
	{
	    this->_Order = tripleCDBond;
	} else if ( order == "dative" )
	{
	    this->_Order = dativeCDBond;
	} else 
	{
	    this->_Order = unknownCDBond;
	}
    }


    void	CDFragment_O::initialize()
    {_OF();
	this->Base::initialize();
	this->_ConstitutionName = _Nil<core::Symbol_O>();
	this->_Nodes.clear();
	this->_AtomsToNodes.clear();
	this->_Bonds.clear();
	this->_LargestId = -1;
	this->_RootNode = _Nil<CDNode_O>();
    }


void CDFragment_O::createImplicitHydrogen(CDNode_sp fromNode, const string& name)
    {_G();
	Atom_sp fromAtom = fromNode->getAtom();
	LOG(BF("From %s") % fromAtom->description() );
	Atom_sp a = Atom_O::create();
	a->setContainedBy(_Nil<core::T_O>());
	a->setName(_lisp->intern(name,ChemKwPkg));
	a->setElement(element_H);
	LOG(BF("Created implicit hydrogen (%s)") % name  );
	fromAtom->bondTo(a,singleBond);
	GC_ALLOCATE(CDNode_O, toNode );
	toNode->setAtom(a);
	toNode->setLabel(name);
	toNode->setId(this->_LargestId+1);
	this->_LargestId++;
	this->_Nodes[toNode->getId()] = toNode;
	this->_AtomsToNodes[a] = toNode;
	GC_ALLOCATE(CDBond_O, bond );
	bond->setIdBegin(fromNode->getId());
	bond->setIdEnd(toNode->getId());
	bond->setBeginNode(fromNode);
	bond->setEndNode(toNode);
	bond->setOrder(singleCDBond);
	this->_Bonds.push_back(bond);
    }



    void	CDFragment_O::parseFromXml(adapt::QDomNode_sp fragment)
    {_G();
	adapt::QDomNode_O::iterator	it;
	this->_Nodes.clear();
	this->_AtomsToNodes.clear();
	for ( it=fragment->begin_Children(); it!=fragment->end_Children(); it++ )
	{
	    adapt::QDomNode_sp child = (*it);
	    if ( child->getLocalName() == "n" )
	    {
		GC_ALLOCATE(CDNode_O, node );
		node->parseFromXml(child);
		int id = node->getId();
		if ( id > this->_LargestId )
		    this->_LargestId = id;
		this->_Nodes[id] = node;
		LOG(BF("Processed node label(%s)") % node->getLabel()  );
	    }
	}
	for ( it=fragment->begin_Children(); it!=fragment->end_Children(); it++ )
	{
	    adapt::QDomNode_sp child = (*it);
	    if ( child->getLocalName() == "b" )
	    {
		GC_ALLOCATE(CDBond_O, bond );
		bond->parseFromXml(child);
		uint idBegin = bond->getIdBegin();
		uint idEnd = bond->getIdEnd();
		ASSERT(this->_Nodes.count(idBegin)>0);
		ASSERT(this->_Nodes.count(idEnd)>0);
		CDNode_sp nodeBegin = this->_Nodes[idBegin];
		CDNode_sp nodeEnd = this->_Nodes[idEnd];
		bond->setBeginNode( nodeBegin );
		bond->setEndNode( nodeEnd);
		this->_Bonds.push_back(bond);
		LOG(BF("Processed bond order(%s) from(%s) to (%s)")
		    % bond->getOrderAsString() % nodeBegin->getLabel() % nodeEnd->getLabel()  );
	    }
	}
    }



    void CDFragment_O::addProperties(core::HashTableEq_sp d)
    {_G();
      d->mapHash([this] (core::T_sp key, core::T_sp value) {
          this->_Properties->setf_gethash(key,value);
        } );
    }


/*!
 * Take a label of the form "xxxx:yyyy" and convert the xxxx part into a keyword symbol and 
 * the yyyy part into a string object
 * return false if there it wasn't in xxxx:yyyy form
 */
    bool CDFragment_O::_asKeyedObject(const string& label, core::Symbol_sp& keyword, core::T_sp& obj)
    {_G();
	core::VectorStrings keyValue;
	LOG(BF("tokenizing string [%s]") % label  );
	core::tokenize( label, keyValue, " ");
	if ( keyValue.size() != 2 )
	{
	    return false;
	}
	if ( keyValue.size() == 1 )
	{
	    SIMPLE_ERROR(BF("There is an error in a ChemDraw file, the string["+label+"] should be a key:value pair"));
	}

	// now we have a key value pair
	// trim off white-space and put it in the _Properties
	string trimKey = core::trimWhiteSpace(keyValue[0]);
	string trimValue = core::trimWhiteSpace(keyValue[1]);
	FIX_ME();
        keyword = _lisp->internKeyword(trimKey);
	obj = core::Str_O::create(trimValue);
	return true;
    }




    bool	CDFragment_O::hasProperty(core::Symbol_sp key)
    {_G();
	return this->getProperties()->contains(key);
    }


    core::T_sp CDFragment_O::setProperty(core::Symbol_sp key,core::T_sp obj)
    {_G();
	if ( this->hasProperty(key) )
	{
	    SIMPLE_ERROR(BF("You have already set the property(%s) it has the value(%s) and you are trying to set it with(%s)") % key->__repr__() % this->_Properties->gethash(key)->__repr__() % obj->__repr__()  );
	}
	LOG(BF("Setting variable(%s) class(%s) to %s") % key % obj->className() % obj->__repr__() );
	this->_Properties->setf_gethash(key,obj);
	return obj;
    }



    string	CDFragment_O::describeProperties()
    {_G();
	stringstream ss;
        this->_Properties->mapHash( [&ss] (core::T_sp key, core::T_sp value) {
	    ss << _rep_(key) << "  value=" << _rep_(value);
          } );
	LOG(BF("%s")%ss.str());
	return ss.str();
    }











/*!
 * Look for properties (key:value) pairs in the nodes and
 * create properties for them.
 */
    bool CDFragment_O::interpret()
    {_G();
	if ( this->_Bonds.size() == 0 )
	{
	    LOG(BF("There are no bonds in this fragment"));
	    return false;
	}
	LOG(BF("Starting a new fragment, number of bonds = %d") % this->_Bonds.size() );
	CDBonds::iterator bi;
	bool foundHashedBond = false;
	adapt::StringSet_sp allNames = adapt::StringSet_O::create();
	this->_Properties = core::HashTableEq_O::create_default();
	{_BLOCK_TRACEF(BF("Processing bonds"));
	    for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ )
	    { _BLOCK_TRACEF(BF("Processing bond with order: %s") % (*bi)->getOrderAsString() );
		LOG(BF("label Begin = %s") % ((*bi)->getBeginNode()->getLabel())  );
		LOG(BF("label End   = %s") % ((*bi)->getEndNode()->getLabel())  );
		allNames->insert((*bi)->getBeginNode()->getLabel());
		allNames->insert((*bi)->getEndNode()->getLabel());
		//
		// Dative bonds must have form begin(XXX:)->end(YYY)
		// A property with key XXX is created  with the value YYY
		//
		if ( (*bi)->getOrder() == dativeCDBond )
		{_BLOCK_TRACE("dative bond");
                  string beginLabel = (*bi)->getBeginNode()->getLabel();
                  string endLabel = (*bi)->getEndNode()->getLabel();
                  string propertyStr = core::trimWhiteSpace(beginLabel);
                  string valueStr = core::trimWhiteSpace(endLabel);
                  this->_Properties->setf_gethash(_lisp->intern(propertyStr,ChemKwPkg),core::Str_O::create(valueStr));
		    LOG(BF("Created property(%s)") % propertyStr  );
		} else if ( (*bi)->getOrder() == hashCDBond )
		{ _BLOCK_TRACE("hash bond");
		    //
		    // There should only be one hashed bond, one side has the
		    // value ":name XXX" where XXX is a constitution name
		    // and the other side is any atom of the constitution.
		    core::Symbol_sp kwBegin; core::T_sp objBegin;
		    core::Symbol_sp kwEnd; core::T_sp objEnd;
		    bool koBegin = this->_asKeyedObject((*bi)->getBeginNode()->getLabel(),kwBegin,objBegin);
		    bool koEnd = this->_asKeyedObject((*bi)->getEndNode()->getLabel(),kwEnd,objEnd);
		    core::Symbol_sp kwPart = _Nil<core::Symbol_O>(); 
		    core::T_sp objPart = _Nil<core::T_O>();
		    CDNode_sp rootNode = _Nil<CDNode_O>();
		    if ( koBegin )
		    {
			// If the bond Begin side is the name:: XXX node
			// then the End side is the root node
			//
			kwPart = kwBegin;
			objPart = objBegin;
			rootNode = (*bi)->getEndNode();
		    }
		    if ( koEnd )
		    {
			if ( kwPart.notnilp() )
			{
			    SIMPLE_ERROR(BF("Invalid hash bond ("
								   + (*bi)->getBeginNode()->getLabel()
								   + " - "
								   + (*bi)->getEndNode()->getLabel()
					    + ")" ));
			}
			kwPart = kwEnd;
			objPart = objEnd;
			rootNode = (*bi)->getBeginNode();
		    }
		    if ( kwPart.nilp() )
		    {
			stringstream serr;
			serr << "A ChemDraw file has a hashed bond with no :name xxx value on one side and a fragment on the other"<<std::endl;
			serr << "On one side the label is ["<< (*bi)->getBeginNode()->getLabel() << "]"<<std::endl;
			serr << "On the other side the label is ["<< (*bi)->getEndNode()->getLabel() << "]"<<std::endl;
			SIMPLE_ERROR(BF(serr.str()));
		    }
		    FIX_ME();
                    core::Symbol_sp sym = _lisp->intern(objPart.as<core::Str_O>()->get());
		    this->_Properties->setf_gethash(kwPart,sym);
		    this->_RootNode = rootNode;
		    LOG(BF("Assigned rootNode: %s") % sym->__repr__() );
		    foundHashedBond = true;
		} else
		{
		    LOG(BF("Doing nothing with bond type(%s)") % (*bi)->getOrderAsString()  );
		}
	    }
	}
	ASSERTF(foundHashedBond,BF("There was no hashed bond for fragment containing the following atoms: %s") % allNames->asString() );
	ASSERTF(this->_RootNode.notnilp(),BF("Read a fragment that did not have a _RootNode defined"));
	this->createAtoms();
	// Now build the residue and fill in the implicit
	// hydrogens on carbon
	Residue_sp res = this->getEntireResidue();
	core::List_sp carbons = res->allAtomsOfElementAsList(element_C);
	{_BLOCK_TRACEF(BF("Creating implicit hydrogens"));
          for ( auto cur : carbons ) {
		Atom_sp c = cur->car<Atom_O>();
		ASSERT(this->_AtomsToNodes.count(c)>0);
		CDNode_sp fromNode = this->_AtomsToNodes[c];
		core::List_sp hydrogens = c->createImplicitHydrogenNamesOnCarbon();
		LOG(BF("Creating implicit hydrogens for %s named %s") % c->description() % hydrogens->__repr__()  );
		for ( auto ih : hydrogens ) {
		    string name = ih->car<core::Str_O>()->get();
		    this->createImplicitHydrogen(fromNode,name);
		}
	    }
	}
	Residue_sp everything = this->getEntireResidue();

    	//
	// Here we want to assign cipPriorities and
	// determine whether atoms will be chiral or prochiral
	//
	{_BLOCK_TRACEF(BF("Assigning cipPriorities"));
	    CipPrioritizer_sp cip = CipPrioritizer_O::create();
	    cip->defineStereochemicalConfigurationsForAllAtoms(everything);
	}

	//
	// Now build the :constitutionAtoms property
	// and the :builtResidue property
	//
	ConstitutionAtoms_sp ca = this->asConstitutionAtoms();
	FIX_ME();
        this->setProperty(_lisp->internKeyword("constitutionAtoms"),ca);
	Residue_sp builtResidue = ca->makeResidue();
	FIX_ME();
        this->setProperty(_lisp->internKeyword("builtResidue"),builtResidue);
	return true;
    }


    Atom_sp CDFragment_O::createOneAtom(CDNode_sp n)
    {_OF();
	Atom_sp a = Atom_O::create();
	a->setContainedBy(_Nil<core::T_O>());
	n->setAtom(a);
	string name;
	int ionization;
	n->getParsedLabel(name,ionization);
	a->setName(_lisp->intern(name,ChemKwPkg));
	a->setIonization(ionization);
	a->setElementFromAtomName();
	a->setStereochemistryType(n->_StereochemistryType);
	a->setConfiguration(n->_Configuration);
	LOG(BF("Just set configuration of atom[%s] to config[%s]")
		     % a->__repr__() % a->getConfigurationAsString() );
	ASSERT(elementIsRealElement(a->getElement()));
	return a;
    }

    void CDFragment_O::createAtoms()
    {_G();
	// First create atoms that are on the ends
	// of solid and dashed bonds
	// 
	CDBonds::iterator bi;
	for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ )
	{
	    CDBondOrder o = ( (*bi)->getOrder() );
	    if ( o == singleCDBond || o == doubleCDBond || o == tripleCDBond
		 || o == singleDashCDBond || o == doubleDashCDBond || o == tripleDashCDBond )
	    {
		CDNode_sp n = (*bi)->getBeginNode();
		if ( n->getAtom().nilp() ) 
		{
		    Atom_sp a = this->createOneAtom(n);
		    this->_AtomsToNodes[a] = n;
		}
		n = (*bi)->getEndNode();
		if ( n->getAtom().nilp() ) 
		{
		    Atom_sp a = this->createOneAtom(n);
		    this->_AtomsToNodes[a] = n;
		}
	    }
	}
    }


    core::T_sp	CDFragment_O::getProperty(core::Symbol_sp s)
    {_G();
	core::HashTableEq_sp dict = this->getProperties();
	if ( !dict->contains(s) )
	{
	    stringstream serr;
	    serr << "Fragment(" << this->getConstitutionName()->__repr__() << ") is missing property: " << s->__repr__() << std::endl;
	    serr << " available property names are: " << std::endl;
            dict->mapHash([&serr](core::T_sp key, core::T_sp val) {
                serr << _rep_(key) << " : " << _rep_(val) << std::endl; } );
	    SIMPLE_ERROR(BF(serr.str()));
	}
	return dict->gethash(s);
    }

    core::T_sp	CDFragment_O::getPropertyOrDefault(core::Symbol_sp s, core::T_sp df)
    {_G();
	core::HashTableEq_sp dict = this->getProperties();
	if ( dict->contains(s) )
	{
	    return dict->gethash(s);
	}
	return df;
    }


    void	CDFragment_O::createBonds(bool selectedAtomsOnly)
    {_G();
	for ( CDBonds::iterator bi=this->_Bonds.begin();
	      bi!=this->_Bonds.end(); bi++ )
	{
//	    CDBondOrder o = ( (*bi)->getOrder() );
	    BondOrder bo = (*bi)->getOrderAsBondOrder();
	    Atom_sp a1 = (*bi)->getBeginNode()->getAtom();
	    ASSERTNOTNULL(a1);
	    if ( a1.nilp() ) continue;
	    Atom_sp a2 = (*bi)->getEndNode()->getAtom();
	    ASSERTNOTNULL(a2);
	    if ( a2.nilp() ) continue;
	    bool createBond = true;
	    if ( selectedAtomsOnly && (!a1->testAllFlags(SELECTED) || !a2->testAllFlags(SELECTED)) )
	    {
		createBond = false;
	    }
	    if ( createBond )
	    {
		LOG(BF("Creating a bond between %s and %s order: %s") % a1->description() % a2->description() % bondOrderToString(bo)  );
		a1->bondTo(a2,bo);
	    }
	}
    }

    void CDFragment_O::removeAllBonds()
    {_G();
        gctools::SmallMap<int,CDNode_sp>::iterator	ni;
	for ( ni = this->_Nodes.begin(); ni!=this->_Nodes.end(); ni++ )
	{
	    Atom_sp a = (*ni).second->getAtom();
	    if ( a.notnilp() )
	    {
		LOG(BF("Removing bonds for atom: %s") % a->description()  );
		a->removeAllBonds();
	    }
	}
    }


    void CDFragment_O::clearAtomSelected()
    {_G();
        gctools::SmallMap<int,CDNode_sp>::iterator	ni;
	for ( ni = this->_Nodes.begin(); ni!=this->_Nodes.end(); ni++ )
	{
	    Atom_sp a = (*ni).second->getAtom();
	    if ( a.notnilp() )
	    {
		a->resetFlags();
	    }
	}
    }


    Residue_sp CDFragment_O::createResidueOfSelectedAtoms()
    {_G();
	Residue_sp res = Residue_O::create();
        gctools::SmallMap<int,CDNode_sp>::iterator	ni;
	for ( ni = this->_Nodes.begin(); ni!=this->_Nodes.end(); ni++ )
	{
	    Atom_sp a = (*ni).second->getAtom();
	    if ( a.notnilp() )
	    {
		if ( a->testAllFlags(SELECTED) ) res->addAtom(a);
	    }
	}
	return res;
    }


/*! Build the residue.  If constitutionOnly is false then 
 * all atoms connected by solid or dashed bonds are included
 * If it is true then only atoms connected by solid bonds are included
 */
    Residue_sp	CDFragment_O::_buildResidue( bool constitutionOnly )
    {_G();
	uint atomCount = 0;
  	// First all of the bonds between atoms
	//
	this->clearAtomSelected();
	this->removeAllBonds();
	this->createBonds(false);

	// Select the atoms that we want to put into
	// the residue
	{_BLOCK_TRACE("Selecting residue atoms");
	    SpanningLoop_sp span = SpanningLoop_O::create();
	    span->setOnlyFollowRealBonds(constitutionOnly);
	    ASSERTNOTNULL(this->_RootNode);
	    ASSERTP(this->_RootNode.notnilp(),"Every fragment must have a root node");
	    span->setTop(this->_RootNode->getAtom());
	    while ( span->advance() )
	    {
		Atom_sp a = span->getAtom();
		a->turnOnFlags(SELECTED);
		LOG(BF("Selected %s") % a->description()  );
		atomCount++;
	    }
	}
	ASSERTP(atomCount>0, "There are no atoms in the residue!");

	this->removeAllBonds();
    	//
    	// Now create bonds just between selected atoms
	//
	this->createBonds(true);
	Residue_sp res = this->createResidueOfSelectedAtoms();
	ASSERTP(res->contentSize()>0,"The residue is empty");
	return res;
    }










    Residue_sp	CDFragment_O::getEntireResidue()
    {_G();
	return this->_buildResidue(false);
    }

    ConstitutionAtoms_sp	CDFragment_O::asConstitutionAtoms()
    {_G();
	Residue_sp residue = this->_buildResidue(true);
	return ConstitutionAtoms_O::create(residue);
    }



    string CDFragment_O::__repr__() const
    {
	stringstream ss;
	ss << "#S(" << this->className() << " ";
	ss << this->_Properties->__repr__();
	ss << ")" << std::endl;
	return ss.str();
    }

#if INIT_TO_FACTORIES

#define ARGS_CDText_O_make "(kprops)"
#define DECL_CDText_O_make ""
#define DOCS_CDText_O_make "make CDText"
  CDText_sp CDText_O::make(core::HashTableEq_sp kprops)
  {_G();
    IMPLEMENT_ME();
#if 0
    GC_ALLOCATE(CDText_O, me );
    me->_CDText._Properties = kprops;
    me->_CDText._Text = "n/a";
    return me;
#endif
  };

#else

    core::T_sp  CDText_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_OF();
	IMPLEMENT_ME();
#if 0
	this->Base::__init__(exec,args,env,lisp);
	core::Cons_sp props = env->lookup(lisp->internWithPackageName(ChemPkg,"props")).as<core::Cons_O>();
	core::HashTableEq_sp kprops = env->lookup(lisp->internWithPackageName(ChemPkg,"kprops")).as<core::HashTableEq_O>();
	this->_Properties = kprops;
	this->_Text = "n/a";
	return _Nil<core::T_O>();
#endif
    }

#endif

    void	CDText_O::initialize()
    {
	this->Base::initialize();
	this->_Properties = core::HashTableEq_O::create_default();
    }


    bool	CDText_O::hasProperties() 
    { 
	return this->_Properties->size()>0; 
    };



/*!
 * Text blocks should be list of key: value pairs separated by line feeds
 */
void	CDText_O::parseFromXml(adapt::QDomNode_sp text)
    {
      adapt::QDomNode_sp sub = text->childWithName("s");
      this->_Text = core::trimWhiteSpace(sub->getData());
      if ( this->_Text[0] != '(' )
      {
        LOG(BF("Text block is not code") );
        return;
      }
      core::StringInputStream_sp sin = core::cl_make_string_input_stream(core::Str_O::create(this->_Text)
                                                                         ,core::clasp_make_fixnum(0)
                                                                         ,_Nil<core::T_O>());
#if 0
      core::Reader_sp reader = core::Reader_O::create(sin,lisp);
      core::Cons_sp block = reader->read(true,_Nil<core::T_O>()).as<core::Cons_O>();
#endif
      core::List_sp block = read_lisp_object(sin,true,_Nil<core::T_O>(),false);
      core::cl_close(sin);

      LOG(BF("Parsed text block: %s\n") % this->_Text);
      if ( block.nilp() )
      {
        SIMPLE_ERROR(BF("Error compiling code:\n"+this->_Text));
      }
      LOG(BF("About to evaluate CDText: %s") % block->__repr__() );
      core::List_sp result = core::eval::evaluate(block,_Nil<core::Environment_O>());
      this->_Properties = core::HashTableEq_O::create_default();
      while (result.notnilp()) {
        core::T_sp key = oCar(result);
        core::T_sp value = oCadr(result);
        result = oCddr(result);
        if ( chem::_sym__PLUS_validChemdrawKeywords_PLUS_->symbolValue().as<core::HashTable_O>()->gethash(key).notnilp() ) {
          this->_Properties->setf_gethash(key,value);
        } else SIMPLE_ERROR(BF("Illegal chemdraw keyword: %s value: %s") % _rep_(key) % _rep_(value) );
      }
    }



    void ChemDraw_O::exposeCando(core::Lisp_sp lisp)
    {
      core::class_<ChemDraw_O>()
//	    .def_raw("core:__init__",&ChemDraw_O::__init__,"(self fileName)")
        .def("getFragments",&ChemDraw_O::getFragments)
        .def("allFragmentsAsCons",&ChemDraw_O::allFragmentsAsCons)
        .def("getSubSetOfFragments",&ChemDraw_O::getSubSetOfFragments)
        .def("setFragmentProperties",&ChemDraw_O::setFragmentProperties)
        .def("asAggregate",&ChemDraw_O::asAggregate)
        ;
      SYMBOL_EXPORT_SC_(KeywordPkg,comment);
      SYMBOL_EXPORT_SC_(KeywordPkg,chiral_centers);
      SYMBOL_EXPORT_SC_(KeywordPkg,group);
      SYMBOL_EXPORT_SC_(KeywordPkg,name_template);
      SYMBOL_EXPORT_SC_(KeywordPkg,pdb_template);
      SYMBOL_EXPORT_SC_(KeywordPkg,restraints);
      SYMBOL_EXPORT_SC_(KeywordPkg,residue_charge);
      SYMBOL_EXPORT_SC_(KeywordPkg,restrained_pi_bonds);
      SYMBOL_EXPORT_SC_(KeywordPkg,caps);
      SYMBOL_EXPORT_SC_(ChemPkg,_PLUS_validChemdrawKeywords_PLUS_);
      adapt::SymbolSet_sp vck = adapt::SymbolSet_O::create();
      vck->insert(kw::_sym_name);
      vck->insert(kw::_sym_comment);
      vck->insert(kw::_sym_chiral_centers);
      vck->insert(kw::_sym_group);
      vck->insert(kw::_sym_name_template);
      vck->insert(kw::_sym_pdb_template);
      vck->insert(kw::_sym_restraints);
      vck->insert(kw::_sym_residue_charge);
      vck->insert(kw::_sym_restrained_pi_bonds);
      vck->insert(kw::_sym_caps);
      chem::_sym__PLUS_validChemdrawKeywords_PLUS_->defparameter(vck);
    }

    void ChemDraw_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,ChemDraw,"","",_lisp)
	    .def("getFragments",&ChemDraw_O::getFragments)
	    .def("allFragmentsAsCons",&ChemDraw_O::allFragmentsAsCons)
	    .def("getSubSetOfFragments",&ChemDraw_O::getSubSetOfFragments)
	    .def("setFragmentProperties",&ChemDraw_O::setFragmentProperties)
	    .def("asAggregate",&ChemDraw_O::asAggregate)
	    ;
#endif
    }




    void ChemDraw_O::lisp_initGlobals(core::Lisp_sp lisp)
    {_G();
    }




    void	ChemDraw_O::initialize()
    {
	this->Base::initialize();
    }


/*
  __BEGIN_DOC(classes.ChemDraw.!class.ChemDraw)
  \requiredKeyed{fileName:}{Text::name}

  Define a ChemDraw object.  Load a cdxml file from \sa{name} and return the ChemDraw object define by it.
  __END_DOC
*/
#if INIT_TO_FACTORIES

#define ARGS_ChemDraw_O_make "(file_name)"
#define DECL_ChemDraw_O_make ""
#define DOCS_ChemDraw_O_make "make ChemDraw"
ChemDraw_sp ChemDraw_O::make(core::T_sp stream)
  {_G();
      GC_ALLOCATE(ChemDraw_O, me );
      me->parse(stream); // me->parse(stream);
      return me;
  };

#else

    core::T_sp 	ChemDraw_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_G();
	string fileName = env->lookup(lisp->internWithPackageName(ChemPkg,"fileName")).as<core::Str_O>()->get();
	LOG(BF("Loading Chemdraw file from: %s") % fileName );
	this->parseFromFileName(fileName);
	return _Nil<core::T_O>();
    }

#endif

#ifdef XML_ARCHIVE
    void	ChemDraw_O::archive(core::ArchiveP node)
    {
	IMPLEMENT_ME();
    }
#endif



#define ARGS_core_ "(name &key comment chiral-centers group name-template pdb-template restraints residue-charge restrained-pi-bonds caps)"
#define DECL_core_ ""
#define DOCS_core_ ""
void ChemDraw_O::setFragmentProperties(core::Symbol_sp name
                                       , core::T_sp comment
                                       , core::T_sp chiral_centers
                                       , core::T_sp group
                                       , core::T_sp name_template
                                       , core::T_sp pdb_template
                                       , core::T_sp restraints
                                       , core::T_sp residue_charge
                                       , core::T_sp restrained_pi_bonds
                                       , core::T_sp caps
                                       )
    {_OF();
	CDFragment_sp frag = this->_NamedFragments.get(name);
	core::HashTableEq_sp properties = core::HashTableEq_O::create_default();
        if ( comment.notnilp() ) properties->setf_gethash(kw::_sym_comment, comment );
        if ( chiral_centers.notnilp() ) properties->setf_gethash(kw::_sym_chiral_centers, chiral_centers );
        if ( group.notnilp() ) properties->setf_gethash(kw::_sym_group, group );
        if ( name_template.notnilp() ) properties->setf_gethash(kw::_sym_name_template, name_template );
        if ( pdb_template.notnilp() ) properties->setf_gethash(kw::_sym_pdb_template, pdb_template );
        if ( restraints.notnilp() ) properties->setf_gethash(kw::_sym_restraints, restraints );
        if ( residue_charge.notnilp() ) properties->setf_gethash(kw::_sym_residue_charge, residue_charge );
        if ( restrained_pi_bonds.notnilp() ) properties->setf_gethash(kw::_sym_restrained_pi_bonds, restrained_pi_bonds );
        if ( caps.notnilp() ) properties->setf_gethash(kw::_sym_caps, caps );
        frag->addProperties(properties);
    }

#if 0
    void	ChemDraw_O::setFragmentProperties(core::List_sp props)
    {_G();
	DEPRECIATED();
	core::HashTableEq_sp kargs = core::HashTableEq_O::createFromKeywordCons(props,validChemdrawKeywords);
	if ( !kargs->contains(_kw_name) )
	{
	    stringstream ss;
	    ss << "Properties for ChemDraw Fragments must have a (name) property, what was passed has: " << kargs->description();
	    SIMPLE_ERROR(BF("%s")%ss.str());
	}
	string name = kargs->lookup(_kw_name).as<core::Str_O>()->get();
	if ( !this->_NamedFragments.contains(name) )
	{
	    stringstream ss;
	    ss << "There are properties for a fragment named("<<name<<") but there is no ChemDraw fragment with that name";
	    SIMPLE_ERROR(BF("%s")%ss.str());
	}
	CDFragment_sp frag = this->_NamedFragments.get(name);
	frag->addProperties(kargs);
    }
#endif

void	ChemDraw_O::parse( core::T_sp strm )
{_G();
  adapt::QDomNode_sp xml = adapt::QDomNode_O::parse(strm);
  if ( !xml->hasChildrenWithName("page") )
    SIMPLE_ERROR(BF("Not a cdxml file" ));
  adapt::QDomNode_sp page = xml->childWithName("page");
  adapt::QDomNode_O::iterator	it;
  this->_NamedFragments.clear();
  for ( it=page->begin_Children(); it!=page->end_Children(); it++ ) {
    adapt::QDomNode_sp child= (*it);
    if ( child->getLocalName() == "fragment" )
    {_BLOCK_TRACE("Processing fragment node");
      GC_ALLOCATE(CDFragment_O, fragment );
      fragment->parseFromXml(child);
      if ( fragment->interpret() )
      {
        core::HashTableEq_sp properties = fragment->getProperties();
        FIX_ME();
        if ( !properties->contains(_lisp->internKeyword("name") ))
        {
          SIMPLE_ERROR(BF("Every fragment must have a property(name:)"));
        }
        FIX_ME();
        core::Symbol_sp constitutionName = properties->gethash(_lisp->internKeyword("name")).as<core::Symbol_O>();
        fragment->setConstitutionName(constitutionName);
        this->_NamedFragments.set(constitutionName,fragment);
        this->_AllFragments.push_back(fragment);
      }
    }
  }
  for ( it=page->begin_Children(); it!=page->end_Children(); it++ )
  {
    adapt::QDomNode_sp child= (*it);
    if ( child->getLocalName() == "t" )
    {_BLOCK_TRACEF(BF("Processing text block"));
      GC_ALLOCATE(CDText_O, text );
      text->parseFromXml(child);
      if ( text->hasProperties() )
      {
        LOG(BF("Found properties: %s") % text->__repr__() );
        core::HashTableEq_sp properties = text->getProperties();
        FIX_ME(); if (!properties->contains(_lisp->internKeyword("name")) )
                  {
                    SIMPLE_ERROR(BF("Every properties block must have a property(name:)"));
                  }
        FIX_ME(); core::Symbol_sp constitutionName = properties->gethash(_lisp->internKeyword("name")).as<core::Symbol_O>();
        if ( !this->_NamedFragments.contains(constitutionName) )
        {
          SIMPLE_ERROR(BF("Could not find fragment with name("+constitutionName->__repr__()+")"));
        }
        CDFragment_sp fragment = this->_NamedFragments.get(constitutionName);
        fragment->addProperties(properties);
      }
    }
  }
}


    Aggregate_sp ChemDraw_O::asAggregate()
    {_OF();
	core::List_sp fragments = this->allFragmentsAsCons();
	Aggregate_sp agg = Aggregate_O::create();
	for ( auto cur : fragments ) {
	    Molecule_sp mol = Molecule_O::create();
	    CDFragment_sp frag = core::oCar(cur).as<CDFragment_O>();
	    Residue_sp res = frag->getEntireResidue();
            core::Symbol_sp name = frag->getProperty(kw::_sym_name).as<core::Symbol_O>();
	    res->setName(name);
	    mol->setName(name);
	    mol->addMatter(res);
	    agg->addMatter(mol);
	}
	return agg;
    }

    core::List_sp	ChemDraw_O::getFragments()
    {_G();
	core::List_sp	frags = _Nil<core::T_O>();
	Fragments::iterator	fi;
	for ( fi=this->_AllFragments.begin(); fi!=this->_AllFragments.end(); fi++ ) {
	    core::Cons_sp n = core::Cons_O::create(*fi,frags);
	    frags = n;
	}
	return frags;
    }


    core::List_sp	ChemDraw_O::getSubSetOfFragments(adapt::SymbolSet_sp namesOfSubSet)
    {_G();
	core::List_sp	frags = _Nil<core::T_O>();
	adapt::SymbolSet_sp namesChosen = adapt::SymbolSet_O::create();
	NamedFragments::iterator	fi;
	for ( fi=this->_NamedFragments.begin(); fi!=this->_NamedFragments.end(); fi++ )
	{
	    if ( namesOfSubSet->contains(fi->first) )
	    {
		namesChosen->insert(fi->first);
		core::Cons_sp n = core::Cons_O::create(fi->second,frags);
		frags = n;
	    }
	}
	if ( namesChosen->size() != namesOfSubSet->size() )
	{
	    adapt::SymbolSet_sp diff = namesOfSubSet->relativeComplement(namesChosen);
	    stringstream ss;
	    ss << "The following names were not found in the ChemDraw object: ";
	    ss << diff->asString();
	    SIMPLE_ERROR(BF("%s")%ss.str());
	}
	return frags;
    }




    
    void CDFragment_O::exposeCando(core::Lisp_sp lisp)
	{
	    core::class_<CDFragment_O>()
		.def("getConstitutionName",&CDFragment_O::getConstitutionName)
		.def("CDFragment-getProperties",&CDFragment_O::getProperties)
		.def("CDFragment-getProperty",&CDFragment_O::getProperty)
		.def("CDFragment-hasProperty",&CDFragment_O::hasProperty)
		.def("CDFragment-setProperty",&CDFragment_O::setProperty)
		.def("CDFragment-getPropertyOrDefault",&CDFragment_O::getPropertyOrDefault)
		.def("getEntireResidue",&CDFragment_O::getEntireResidue)
		.def("asConstitutionAtoms",&CDFragment_O::asConstitutionAtoms)
		.def("describeProperties",&CDFragment_O::describeProperties)
		;
	}

    void CDFragment_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(ChemPkg,CDFragment,"","",_lisp)
		.def("getPropertyOrDefault",&CDFragment_O::getPropertyOrDefault)
		;
#endif //]
	}





    REGISTER_CLASS(chem,CDNode_O);
    REGISTER_CLASS(chem,CDBond_O);
    EXPOSE_CLASS(chem,CDFragment_O);
    REGISTER_CLASS(chem,CDText_O);
    EXPOSE_CLASS_AND_GLOBALS(chem,ChemDraw_O);
};



