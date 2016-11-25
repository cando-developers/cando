/*
    File: chemdraw.cc
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
#define	DEBUG_LEVEL_FULL

#include <string.h>
#include <clasp/core/common.h>
#include <clasp/core/str.h>
#include <clasp/core/primitives.h>
#include <clasp/core/bformat.h>
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
#include <clasp/core/symbolTable.h>
#include <cando/chem/spanningLoop.h>
#include <clasp/core/reader.h>
#include <cando/chem/cipPrioritizer.h>
#include <cando/chem/constitutionAtoms.h>
#include <clasp/core/wrappers.h>


namespace chem
{

SYMBOL_EXPORT_SC_(ChemPkg,name);

// Ensure that property is a Keyword symbol
void set_property(core::HashTableEq_sp properties, core::T_sp property, core::T_sp val)
{
  core::Symbol_sp psym = gc::As<core::Symbol_sp>(property);
  if (!psym->isKeywordSymbol()) {
    SIMPLE_ERROR(BF("All properties must be keyword symbols property: %s") % _rep_(psym));
  }
  properties->setf_gethash(psym,val);
}

void	CDNode_O::initialize()
{
  this->Base::initialize();
  this->_Atom = _Nil<Atom_O>();
}

#if 0
void CDNode_O::fieldsp(core::Record_sp node)
{
  node->field(INTERN_(kw,atom_properties),this->_AtomProperties);
  node->field(INTERN_(kw,residue_properties),this->_ResidueProperties);
  node->field(INTERN_(kw,molecule_properties),this->_MoleculeProperties);
  this->Base(node);
}
#endif


std::string CDNode_O::_extractLabel(adapt::QDomNode_sp node)
{
  string name;
  if ( node->hasChildrenWithName("t") )
  {
    adapt::QDomNode_sp text = node->childWithName("t");
    if ( text->hasChildrenWithName("s") )
    {
      core::List_sp xmls = text->childrenWithName("s");
      stringstream ss;
      for ( auto c : xmls ) {
        adapt::QDomNode_sp child = gctools::As<adapt::QDomNode_sp>(oCar(c));
        ss << child->getData();
      }
      string name = ss.str();
      return name;
    }
  }
  return "C_"+node->getAttributeString("id");
}


void CDNode_O::getParsedLabel(string& name, int& ionization) const
{
  string label = this->_Label;
  vector<string> parts = core::split(label,"/");
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
void	CDNode_O::parseFromXml(adapt::QDomNode_sp xml, bool print)
{_OF();
  this->_Id = xml->getAttributeInt("id");
  this->_Color = xml->getAttributeIntDefault("color",3);
  this->_Label = this->_extractLabel(xml);
  if (print) BFORMAT_T(BF("CDNode id(%s) color(%s) label(%s)\n") % this->_Id % this->_Color % this->_Label );
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
        if (print) BFORMAT_T(BF("    configuration(S)\n"));
        LOG(BF("Atom has geometry AS=%s") % as );
        this->_Configuration = S_Configuration;
      } else if ( as == "R" )
      {
        LOG(BF("Atom has geometry AS=%s") % as );
        if (print) BFORMAT_T(BF("    configuration(R)\n"));
        this->_Configuration = R_Configuration;
      } else
      {
        if (print) BFORMAT_T(BF("    could not determin configuration from AS\n"));
        LOG(BF("Could not interpret geometry AS[%s]") % as );
      }
    }
  }
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
  case hollowWedgeCDBond: return "hollowWedge";
  case wedgeHashCDBond: return "wedgeHash";
  case wedgeCDBond: return "wedge";
  case wavyCDBond: return "wavy";
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
  default:
      return noBond;
  }
  return noBond;
}


void	CDBond_O::parseFromXml(adapt::QDomNode_sp xml, bool print)
{
  this->_IdBegin = xml->getAttributeInt("B");
  this->_IdEnd = xml->getAttributeInt("E");
  string order = xml->getAttributeStringDefault("Order","1");
  string display = xml->getAttributeStringDefault("Display","");
  if ( order == "1" ) {
    if ( display == "Dash" ) {
      this->_Order = singleDashCDBond;
    } else if ( display == "Hash" ) {
      this->_Order = hashCDBond;
    } else if ( display == "HollowWedgeBegin" ) {
      this->_Order = hollowWedgeCDBond;
    } else if ( display == "WedgedHashBegin" ) {
      this->_Order = wedgeHashCDBond;
    } else if ( display == "WedgedBegin" ) {
      this->_Order = wedgeCDBond;
    } else if ( display == "Wavy" ) {
      this->_Order = wavyCDBond;
    } else {
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
    Warn(core::Str_O::create((BF("Unknown bond order %s") % this->_Order)),
         _Nil<core::T_O>());
    this->_Order = unknownCDBond;
  }
  if (print) BFORMAT_T(BF("CDBond _IdBegin(%s) _IdEnd(%s) order(%s) display(%s)\n") % this->_IdBegin % this->_IdEnd % order % display );
}


void	CDFragment_O::initialize()
{_OF();
  this->Base::initialize();
  this->_ConstitutionName = _Nil<core::Symbol_O>();
  this->_Nodes.clear();
  this->_AtomsToNodes.clear();
  this->_Bonds.clear();
  this->_LargestId = -1;
}

#if 0
void CDFragment_O::createImplicitHydrogen(CDNode_sp fromNode, const std::string& name)
{
  Atom_sp fromAtom = fromNode->getAtom();
  LOG(BF("From %s") % fromAtom->description() );
  Atom_sp a = Atom_O::create();
  a->setContainedBy(_Nil<core::T_O>());
  a->setName(chemkw_intern(name));
  a->setElement(element_H);
  LOG(BF("Created implicit hydrogen (%s)") % name  );
  fromAtom->bondTo(a,singleBond);
  GC_ALLOCATE(CDNode_O, toNode );
  toNode->setAtom(a);
  toNode->_Label = name;
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
#endif


void	CDFragment_O::parseFromXml(adapt::QDomNode_sp fragment, bool print)
{
  adapt::QDomNode_O::iterator	it;
  this->_Nodes.clear();
  this->_AtomsToNodes.clear();
  if (print) BFORMAT_T(BF("CDFragment - starting\n"));
  for ( it=fragment->begin_Children(); it!=fragment->end_Children(); it++ ) {
    adapt::QDomNode_sp child = (*it);
    if ( child->getLocalName() == "n" ) {
      GC_ALLOCATE(CDNode_O, node );
      node->parseFromXml(child, print);
      int id = node->getId();
      if ( id > this->_LargestId )
        this->_LargestId = id;
      this->_Nodes[id] = node;
      LOG(BF("Processed node label(%s)") % node->getLabel()  );
    }
  }
  for ( it=fragment->begin_Children(); it!=fragment->end_Children(); it++ ) {
    adapt::QDomNode_sp child = (*it);
    if ( child->getLocalName() == "b" ) {
      GC_ALLOCATE(CDBond_O, bond );
      bond->parseFromXml(child,print);
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
  if (print) BFORMAT_T(BF("CDFragment - done.\n"));
}

int CDFragment_O::countNeighbors(CDNode_sp node)
{
  int num = 0;
  for ( auto cur : this->_Bonds ) {
    if ( cur->getBeginNode() == node ) num++;
    if ( cur->getEndNode() == node ) num++;
  }
  return num;
}

#if 0
void CDFragment_O::addProperties(core::HashTableEq_sp d)
{
  d->mapHash([this] (core::T_sp key, core::T_sp value) {
      set_property(this->_Properties,key,value);
    } );
}
#endif

/*!
 * Take a label of the form "xxxx:yyyy" and convert the xxxx part into a keyword symbol and 
 * the yyyy part into a string object
 * return false if there it wasn't in xxxx:yyyy form
 */
bool CDFragment_O::_asKeyedObject(core::Symbol_sp labelSym, core::Symbol_sp& keyword, core::T_sp& obj)
{
  string label = labelSym->symbolNameAsString();
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
  keyword = _lisp->intern(trimKey,KeywordPkg);
  obj = core::Str_O::create(trimValue);
  return true;
}




#if 0
bool	CDFragment_O::hasProperty(core::Symbol_sp key)
{
  return this->getProperties()->contains(key);
}


core::T_sp CDFragment_O::setProperty(core::Symbol_sp key,core::T_sp obj)
{
  if ( this->hasProperty(key) ) {
    SIMPLE_ERROR(BF("You have already set the property(%s) it has the value(%s) and you are trying to set it with(%s)") % _rep_(key) % _rep_(this->_Properties->gethash(key)) % _rep_(obj)  );
  }
  LOG(BF("Setting variable(%s) class(%s) to %s") % key % obj->className() % _rep_(obj) );
  set_property(this->_Properties,key,obj);
  return obj;
}



string	CDFragment_O::describeProperties()
{
  stringstream ss;
  this->_Properties->mapHash( [&ss] (core::T_sp key, core::T_sp value) {
      ss << _rep_(key) << "  value=" << _rep_(value) << " ";
    } );
  LOG(BF("%s")%ss.str());
  return ss.str();
}
#endif





core::Symbol_mv parse_property(const string& propertyValue, CDBond_sp bond, const string& otherSideValue)
{
  core::T_sp stream = core::cl__make_string_input_stream(core::Str_O::create(propertyValue),core::clasp_make_fixnum(0),_Nil<core::T_O>());
//  printf("%s:%d Parsing property string: %s\n", __FILE__, __LINE__, propertyValue.c_str());
  core::T_sp eof = core::Cons_O::create(_Nil<core::T_O>());
  core::DynamicScopeManager scope(cl::_sym_STARpackageSTAR,_lisp->findPackage(ChemKwPkg));
  core::T_sp property = core::cl__read(stream,_Nil<core::T_O>(),eof);
  if ( property == eof ) {
    SIMPLE_ERROR(BF("Could not parse first part of \"%s\" as a (symbol value) pair - in property bond of order %s other side of bond is \"%s\"") % propertyValue % bond->getOrderAsString() % otherSideValue );
  }
//  printf("%s:%d Parsed property: %s\n", __FILE__, __LINE__, _rep_(property).c_str());
  core::T_sp value = core::cl__read(stream,_Nil<core::T_O>(),eof);
  if ( value == eof ) {
    SIMPLE_ERROR(BF("Could not parse second part of \"%s\" as a (symbol value) pair - in property bond of order %s other side of bond is \"%s\"") % propertyValue % bond->getOrderAsString() % otherSideValue );
  }
//  printf("%s:%d Parsed value: %s\n", __FILE__, __LINE__, _rep_(value).c_str());
  core::cl__close(stream);
  if ( core::Symbol_sp key = property.asOrNull<core::Symbol_O>() ) {
    return Values(property,value);
  }
  SIMPLE_ERROR(BF("Could not parse \"%s\" as a (symbol value) pair - in property bond of order %s other side of bond is \"%s\"") % propertyValue % bond->getOrderAsString() % otherSideValue );
}




/*!
 * Look for edges that specify properties and move them into the CDNodes that
* they target
 */
bool CDFragment_O::interpret()
{
  printf("%s:%d  Interpreting a fragment\n", __FILE__, __LINE__ );
  int nextFragmentNameIndex = 1;
  if ( this->_Bonds.size() == 0 ) {return false;}
  CDBonds::iterator bi;
  bool foundHashedBond = false;
//  adapt::SymbolSet_sp allNames = adapt::SymbolSet_O::create();
  {
    for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ ) {
      auto cdorder = (*bi)->getOrder();
      //
      // Dative bonds have a property on one end of the bond
      // and CDNode on the other end
      // The property is added to the _AtomProperties of the CDNode
      // and when a molecule is created later the atom generated
      // from the node will get those properties.
      //
      if ( cdorder == wedgeHashCDBond
           || cdorder == singleCDBond
           || cdorder == doubleCDBond
           || cdorder == tripleCDBond
           || cdorder == singleDashCDBond
           || cdorder == doubleDashCDBond
           || cdorder == tripleDashCDBond ) {
        // Do nothing
      } else if ( cdorder == dativeCDBond || // Atom
                  cdorder == hollowWedgeCDBond || // Residue
                  cdorder == wavyCDBond ) { // Molecule property
        // define properties
        // On one side it will have no other nodes connected to it and the
        // other side has a CDNode that will represent an atom in a residue in a molecule.
        // Add the property _AtomProperties, _ResidueProperties or _MoleculeProperties
        // depending on the type of the bond
        int beginNeighbors = this->countNeighbors((*bi)->getBeginNode());
        int endNeighbors = this->countNeighbors((*bi)->getEndNode());
        CDNode_sp targetNode, propertyNode;
//        printf("%s:%d %s bond has %d and %d neighbors\n", __FILE__, __LINE__, (*bi)->getOrderAsString().c_str(), beginNeighbors, endNeighbors);
        if ( beginNeighbors == 1 ) {
          targetNode = (*bi)->getEndNode();
          propertyNode = (*bi)->getBeginNode();
        } else if ( endNeighbors == 1 ) {
          targetNode = (*bi)->getBeginNode();
          propertyNode = (*bi)->getEndNode();
        } else {
          SIMPLE_ERROR(BF("The %s bond must have one end with only one neighbor - the one in question has %d and %d neighbors\n") % (*bi)->getOrderAsString() % beginNeighbors % endNeighbors );
        }
        string propertyCode = propertyNode->_Label;
        core::Symbol_mv parsedProperty = parse_property(propertyCode, *bi, targetNode->_Label);
        core::T_sp value = parsedProperty.second();
        if ( parsedProperty.number_of_values() == 2 ) {
          if ( cdorder == dativeCDBond ) {
            targetNode->_AtomProperties = core__put_f(targetNode->_AtomProperties,value,parsedProperty);
          } else if ( cdorder == hollowWedgeCDBond ) {
            targetNode->_ResidueProperties = core__put_f(targetNode->_ResidueProperties,value,parsedProperty);
          } else if ( cdorder == wavyCDBond ) {
            targetNode->_MoleculeProperties = core__put_f(targetNode->_MoleculeProperties,value,parsedProperty);
          } else {
            SIMPLE_ERROR(BF("Cannot interpret bond %s in terms of where to put the property") % (*bi)->getOrderAsString().c_str());
          }
        } else {
          SIMPLE_ERROR(BF("What should I do when parse_property returns only one value from: %s") % propertyCode );
        }
      } else {
        Warn(core::Str_O::create((BF("Doing nothing with bond type %d")%(*bi)->getOrderAsString()).str()),
             _Nil<core::T_O>());
      }
    }
  }

  // Create an atom for every node that is on the end of a
  // single/double/triple solid or dashed bond
  this->createAtomsAndBonds();
  // Now build the residues and fill in the implicit
  // hydrogens on carbon
  Molecule_sp mol = this->createMolecule();
// IF YOU COMMENT OUT THE NEXT LINE THEN A TRULY FRIGHTENING BUG WILL HAPPEN
  this->_Molecule = mol;
  core::List_sp carbons = mol->allAtomsOfElementAsList(element_C);
  for ( core::List_sp cur = carbons; cur.consp(); cur = oCdr(cur) ) {
    core::T_sp ta = oCar(cur);
    Atom_sp a = gc::As<Atom_sp>(ta);
    a->fillInImplicitHydrogens();
  }

  core::List_sp nitrogens = mol->allAtomsOfElementAsList(element_N);
  for ( core::List_sp cur = nitrogens; cur.consp(); cur = oCdr(cur) ) {
    core::T_sp ta = oCar(cur);
    Atom_sp a = gc::As<Atom_sp>(ta);
    a->fillInImplicitHydrogens();
  }

  core::List_sp oxygens = mol->allAtomsOfElementAsList(element_O);
  for ( core::List_sp cur = oxygens; cur.consp(); cur = oCdr(cur) ) {
    core::T_sp ta = oCar(cur);
    Atom_sp a = gc::As<Atom_sp>(ta);
    a->fillInImplicitHydrogens();
  }

  {_BLOCK_TRACEF(BF("Assigning cipPriorities"));
    CipPrioritizer_sp cip = CipPrioritizer_O::create();
    cip->defineStereochemicalConfigurationsForAllAtoms(mol);
  }
#if 0    

  Residue_sp res = this->getEntireResidue();
  core::List_sp carbons = res->allAtomsOfElementAsList(element_C);
  {_BLOCK_TRACEF(BF("Creating implicit hydrogens"));
    for ( auto cur : carbons ) {
      Atom_sp c = cur->car<Atom_O>();
      ASSERT(this->_AtomsToNodes.count(c)>0);
      CDNode_sp fromNode = this->_AtomsToNodes[c];
      core::List_sp hydrogens = c->createImplicitHydrogenNamesOnCarbon();
      LOG(BF("Creating implicit hydrogens for %s named %s") % c->description() % _rep_(hydrogens)  );
      for ( auto ih : hydrogens ) {
        core::Symbol_sp name = gc::As<core::Symbol_sp>(oCar(ih));
        this->createImplicitHydrogen(fromNode,name);
      }
    }
  }
#endif

#if 0
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
  this->setProperty(INTERN_(kw,constitutionAtoms),ca);
  Residue_sp builtResidue = ca->makeResidue();
  this->setProperty(INTERN_(kw,builtResidue),builtResidue);
  return true;
#endif
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
  a->setName(chemkw_intern(name));
  a->setIonization(ionization);
  a->setElementFromAtomName();
  a->setProperty(INTERN_(kw,chemdraw_color),core::clasp_make_fixnum(n->_Color));
  a->setStereochemistryType(n->_StereochemistryType);
  a->setConfiguration(n->_Configuration);
  LOG(BF("Just set configuration of atom[%s] to config[%s]")
      % _rep_(a) % a->getConfigurationAsString() );
  ASSERT(elementIsRealElement(a->getElement()));
  return a;
}

void CDFragment_O::createAtomsAndBonds()
{
	// First create atoms that are on the ends
	// of solid and dashed bonds
	// 
  CDBonds::iterator bi;
  for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ )
  {
    CDBondOrder o = ( (*bi)->getOrder() );
    gc::Nilable<Atom_sp> beginAtom;
    gc::Nilable<Atom_sp> endAtom;
    if ( o == wedgeHashCDBond
         || o == singleCDBond || o == doubleCDBond || o == tripleCDBond
         || o == singleDashCDBond || o == doubleDashCDBond || o == tripleDashCDBond )
    {
      CDNode_sp n = (*bi)->getBeginNode();
      if ( n->getAtom().nilp() ) {
        beginAtom = this->createOneAtom(n);
        this->_AtomsToNodes[beginAtom] = n;
      } else beginAtom = n->getAtom();
      n = (*bi)->getEndNode();
      if ( n->getAtom().nilp() ) {
        endAtom = this->createOneAtom(n);
        this->_AtomsToNodes[endAtom] = n;
      } else endAtom = n->getAtom();
      Bond_sp bond = Bond_O::create(beginAtom,endAtom);
      switch (o) {
      case singleDashCDBond:
          bond->setProperty(INTERN_(kw,chemdraw_dashed_bond),_lisp->_true());
      case wedgeHashCDBond:
      case singleCDBond:
          bond->setOrder(singleBond);
          break;
      case doubleDashCDBond:
          bond->setProperty(INTERN_(kw,chemdraw_dashed_bond),_lisp->_true());
      case doubleCDBond:
          bond->setOrder(doubleBond);
          break;
      case tripleDashCDBond:
          bond->setProperty(INTERN_(kw,chemdraw_dashed_bond),_lisp->_true());
      case tripleCDBond:
          bond->setOrder(tripleBond);
          break;
      default:
          Warn(core::Str_O::create((BF("Add support for ChemDraw bond: %s") % (*bi)->getOrderAsString()).str()), _Nil<core::T_O>() );
      }
      beginAtom->addBond(bond);
      endAtom->addBond(bond);
    }
  }
}


Molecule_sp CDFragment_O::createMolecule()
{
  core::HashTable_sp atomsToNodes = core::HashTableEq_O::create_default();
  // Gather all of the nodes that have atoms
  for ( auto cur : this->_Nodes ) {
    if ( cur.second->_Atom.notnilp() ) {
      atomsToNodes->setf_gethash(cur.second->_Atom,cur.second);
      cur.second->_Atom->applyPropertyList(cur.second->_AtomProperties);
    }
  }
  // Loop, pick the first atom remaining in the hash table and build a
  // spanning tree, avoiding inter-residue bonds
  Molecule_sp mol = Molecule_O::create();
  do {
    Atom_sp firstAtom;
    atomsToNodes->map_while_true( [&firstAtom] (core::T_sp atom, core::T_sp value ) -> bool {
        firstAtom = gc::As<Atom_sp>(atom);
        return false;
      } );
    if ( firstAtom ) {
    // Build a spanning loop from this atom and a Residue to hold the atoms
      Residue_sp res = Residue_O::create();
      mol->addMatter(res);
      SpanningLoop_sp span = SpanningLoop_O::create();
//      printf("%s:%d starting spanning tree from atom: %s\n", __FILE__, __LINE__, _rep_(firstAtom).c_str());
      span->setTop(firstAtom);
      while (span->advanceLoopAndProcessWhenTestTrue( [] (Atom_sp a, Bond_sp b) -> bool {
            bool advance = !b->hasProperty(INTERN_(kw,chemdraw_dashed_bond));
//            printf("%s:%d looking at atom: %s  bond: %s\n", __FILE__, __LINE__, _rep_(a).c_str(), _rep_(b).c_str());
            return advance;
          } ) ) {
        Atom_sp a = span->getAtom();
//        printf("%s:%d adding to residue atom: %s\n", __FILE__, __LINE__, _rep_(a).c_str());
        res->addAtom(a);
        CDNode_sp atomNode = gc::As<CDNode_sp>(atomsToNodes->gethash(a));
//        printf("%s:%d Applying residue and molecule properties\n", __FILE__, __LINE__ );
//        printf("%s:%d    residue properties: %s\n", __FILE__, __LINE__, _rep_(atomNode->_ResidueProperties).c_str());
//        printf("%s:%d    molecule properties: %s\n", __FILE__, __LINE__, _rep_(atomNode->_MoleculeProperties).c_str());
        res->applyPropertyList(atomNode->_ResidueProperties);
        mol->applyPropertyList(atomNode->_MoleculeProperties);
        atomsToNodes->remhash(a);
      }
    }
  } while (atomsToNodes->hashTableCount()>0);
  return mol;
};

#if 0
core::T_sp	CDFragment_O::getProperty(core::Symbol_sp s)
{
  core::HashTableEq_sp dict = this->getProperties();
  if ( !dict->contains(s) )
  {
    stringstream serr;
    serr << "Fragment(" << _rep_(this->getConstitutionName()) << ") is missing property: " << s->__repr__() << std::endl;
    serr << " available property names are: " << std::endl;
    dict->mapHash([&serr](core::T_sp key, core::T_sp val) {
        serr << _rep_(key) << " : " << _rep_(val) << std::endl; } );
    SIMPLE_ERROR(BF(serr.str()));
  }
  return dict->gethash(s);
}

core::T_sp	CDFragment_O::getPropertyOrDefault(core::Symbol_sp s, core::T_sp df)
{
  core::HashTableEq_sp dict = this->getProperties();
  if ( dict->contains(s) )
  {
    return dict->gethash(s);
  }
  return df;
}
#endif

void	CDFragment_O::createBonds(bool selectedAtomsOnly)
{
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
{
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
{
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
{
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
#if 0
Residue_sp	CDFragment_O::_buildResidue( bool constitutionOnly )
{
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
#endif







#if 0
Residue_sp	CDFragment_O::getEntireResidue()
{
  return this->_buildResidue(false);
}
#endif
CL_LISPIFY_NAME("asConstitutionAtoms");
CL_DEFMETHOD     ConstitutionAtoms_sp	CDFragment_O::asConstitutionAtoms()
{
  IMPLEMENT_ME();
#if 0
  Residue_sp residue = this->_buildResidue(true);
  return ConstitutionAtoms_O::create(residue);
#endif
}


string CDFragment_O::__repr__() const
{
  stringstream ss;
  ss << "#S(" << this->className() << " ";
//	ss << this->_Properties->__repr__();
  ss << ")" << std::endl;
  return ss.str();
}

CDText_sp CDText_O::make(core::HashTableEq_sp kprops)
{
  IMPLEMENT_ME();
};

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
void	CDText_O::parseFromXml(adapt::QDomNode_sp text, bool print)
{
  core::List_sp xmls = text->childrenWithName("s");
  stringstream ss;
  for ( auto c : xmls ) {
    adapt::QDomNode_sp child = gctools::As<adapt::QDomNode_sp>(oCar(c));
    ss << child->getData();
  }
  string name = ss.str();
  this->_Text = core::trimWhiteSpace(name);
  if ( this->_Text[0] != '(' )
  {
    LOG(BF("Text block is not code") );
    return;
  }
  if (print) BFORMAT_T(BF("CDText parsed: %s") % this->_Text);
  core::StringInputStream_sp sin = core::cl__make_string_input_stream(core::Str_O::create(this->_Text)
                                                                      ,core::clasp_make_fixnum(0)
                                                                      ,_Nil<core::T_O>());
#if 0
  core::Reader_sp reader = core::Reader_O::create(sin,lisp);
  core::Cons_sp block = reader->read(true,_Nil<core::T_O>()).as<core::Cons_O>();
#endif
  core::List_sp block = read_lisp_object(sin,true,_Nil<core::T_O>(),false);
  core::cl__close(sin);
  LOG(BF("Parsed text block: %s\n") % this->_Text);
  if ( block.nilp() ) {
    SIMPLE_ERROR(BF("Error compiling code:\n"+this->_Text));
  }
  LOG(BF("About to evaluate CDText: %s") % _rep_(block) );
  core::List_sp result = core::eval::evaluate(block,_Nil<core::Environment_O>());
  this->_Properties = core::HashTableEq_O::create_default();
  while (result.notnilp()) {
    core::T_sp key = oCar(result);
    core::T_sp value = oCadr(result);
    result = oCddr(result);
    if ( chem::_sym__PLUS_validChemdrawKeywords_PLUS_->symbolValue().as<core::HashTable_O>()->gethash(key).notnilp() ) {
      set_property(this->_Properties,key,value);
      if (print) BFORMAT_T(BF("   setting property %s to %s\n") % core::_rep_(key) % core::_rep_(value));
    } else SIMPLE_ERROR(BF("Illegal chemdraw keyword: %s value: %s") % _rep_(key) % _rep_(value) );
  }
}

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
CL_INITIALIZE void initialize_valid_chemdraw_keywords()
{
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
};

  








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

CL_LISPIFY_NAME(make-chem-draw);
CL_LAMBDA(file_name &optional print);
CL_DOCSTRING("Make a chem:chem-draw object from a string.  If print is T then print info to *standard-output*.");
CL_DEFUN ChemDraw_sp ChemDraw_O::make(core::T_sp stream, bool print)
{
  GC_ALLOCATE(ChemDraw_O, me );
  me->parse(stream,print); // me->parse(stream);
  return me;
};

#ifdef XML_ARCHIVE
void	ChemDraw_O::archive(core::ArchiveP node)
{
  IMPLEMENT_ME();
}
#endif


#if 0
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
#endif


#if 0
void	ChemDraw_O::setFragmentProperties(core::List_sp props)
{
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



void ChemDraw_O::parseChild( adapt::QDomNode_sp child, bool print )
{
  if (print) BFORMAT_T(BF("ChemDraw_O::parse child of page with name(%s)\n") % child->getLocalName());
  if ( child->getLocalName() == "fragment" ) {
    GC_ALLOCATE(CDFragment_O, fragment );
    fragment->parseFromXml(child,print);
    if ( fragment->interpret() ) {
#if 0
      core::HashTableEq_sp properties = fragment->getProperties();
      if ( !properties->contains(INTERN_(kw,name)))
      {
        SIMPLE_ERROR(BF("Every fragment must have a property(%s) available properties: %s") % _rep_(INTERN_(kw,name)) % properties->keysAsString());
      }
      core::Symbol_sp constitutionName = properties->gethash(INTERN_(kw,name)).as<core::Symbol_O>();
      fragment->setConstitutionName(constitutionName);
      this->_NamedFragments.set(constitutionName,fragment);
#endif
      this->_AllFragments.push_back(fragment);
    } else {
      SIMPLE_ERROR(BF("Could not interpret a ChemDraw CDFragment"));
    }
  } else if ( child->getLocalName() == "t" ) {
    GC_ALLOCATE(CDText_O, text );
    text->parseFromXml(child,print);
#if 0
    if ( text->hasProperties() )
    {
      LOG(BF("Found properties: %s") % _rep_(text) );
      core::HashTableEq_sp properties = text->getProperties();
      if (!properties->contains(INTERN_(kw,name)))
      {
        SIMPLE_ERROR(BF("Every properties block must have a property(name:)"));
      }
      core::Symbol_sp constitutionName = properties->gethash(INTERN_(kw,name)).as<core::Symbol_O>();
      if ( !this->_NamedFragments.contains(constitutionName) )
      {
        SIMPLE_ERROR(BF("Could not find fragment with name("+_rep_(constitutionName)+")"));
      }
      CDFragment_sp fragment = this->_NamedFragments.get(constitutionName);
      fragment->addProperties(properties);
    }
#endif
  } else if ( child->getLocalName() == "group" ) {
    if (print) BFORMAT_T(BF("ChemDraw_O::parsing group start...\n"));
    for ( adapt::QDomNode_O::iterator it=child->begin_Children(); it!=child->end_Children(); it++ ) {
      this->parseChild(*it,print);
    }
    if (print) BFORMAT_T(BF("ChemDraw_O::parsing group done.\n"));
  }
}


void	ChemDraw_O::parse( core::T_sp strm, bool print )
{
  if (print) BFORMAT_T(BF("ChemDraw_O::parse starting\n"));
  adapt::QDomNode_sp xml = adapt::QDomNode_O::parse(strm);
  if ( !xml->hasChildrenWithName("page") )
    SIMPLE_ERROR(BF("Not a cdxml file" ));
  adapt::QDomNode_sp page = xml->childWithName("page");
  adapt::QDomNode_O::iterator	it;
  this->_NamedFragments.clear();
  for ( it=page->begin_Children(); it!=page->end_Children(); it++ ) {
    this->parseChild(*it,print);
  }
  if (print) BFORMAT_T(BF("ChemDraw_O::parse done.\n"));
}


CL_LISPIFY_NAME("asAggregate");
CL_DEFMETHOD     Aggregate_sp ChemDraw_O::asAggregate()
{_OF();
  core::List_sp fragments = this->allFragmentsAsList();
  Aggregate_sp agg = Aggregate_O::create();
  for ( auto cur : fragments ) {
    CDFragment_sp frag = core::oCar(cur).as<CDFragment_O>();
    Molecule_sp mol = gc::As<Molecule_sp>(frag->getMolecule());
    agg->addMatter(mol);
  }
  return agg;
}

CL_LISPIFY_NAME("getFragments");
CL_DEFMETHOD     core::List_sp	ChemDraw_O::getFragments()
{
  core::List_sp	frags = _Nil<core::T_O>();
  Fragments::iterator	fi;
  for ( fi=this->_AllFragments.begin(); fi!=this->_AllFragments.end(); fi++ ) {
    core::Cons_sp n = core::Cons_O::create(*fi,frags);
    frags = n;
  }
  return frags;
}


CL_LISPIFY_NAME("getSubSetOfFragments");
CL_DEFMETHOD     core::List_sp	ChemDraw_O::getSubSetOfFragments(adapt::SymbolSet_sp namesOfSubSet)
{
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




    







};
