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
#define	DEBUG_LEVEL_NONE

#include <string.h>
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <clasp/core/primitives.h>
#include <clasp/core/bformat.h>
#include <cando/adapt/adapters.h>
#include <cando/adapt/stringSet.h>
#include <cando/chem/chemdraw.h>
#include <cando/adapt/quickDom.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/bond.h>
#include <cando/chem/residue.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/aggregate.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/lispReader.h>
#include <clasp/core/symbolTable.h>
#include <cando/geom/vector3.h>
#include <cando/chem/loop.h>
#include <cando/chem/spanningLoop.h>
#include <cando/chem/cipPrioritizer.h>
#include <cando/chem/constitutionAtoms.h>
#include <clasp/core/wrappers.h>
#include <cando/geom/ovector2.h>

namespace chem
{

SYMBOL_EXPORT_SC_(ChemPkg,name);

// Ensure that property is a Keyword symbol
void set_property(core::HashTable_sp properties, core::T_sp property, core::T_sp val)
{
  core::Symbol_sp psym = gc::As<core::Symbol_sp>(property);
  if (!psym->isKeywordSymbol()) {
    SIMPLE_ERROR("All properties must be keyword symbols property: {}" , _rep_(psym));
  }
  properties->setf_gethash(psym,val);
}

void	CDNode_O::initialize()
{
  this->Base::initialize();
  this->_Atom = nil<Atom_O>();
}
void CDNode_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,id),this->_Id);
  node->field(INTERN_(kw,color),this->_Color);
  node->field(INTERN_(kw,pos),this->_Pos);
  node->field(INTERN_(kw,label),this->_Label);
  node->/*pod_*/field_if_not_default( INTERN_(kw,stereochemistryType), this->_StereochemistryType, undefinedCenter );
  node->/*pod_*/field_if_not_default( INTERN_(kw,configuration), this->_Configuration, undefinedConfiguration  );
  node->field(INTERN_(kw,atom_properties),this->_AtomProperties);
  node->field(INTERN_(kw,residue_properties),this->_ResidueProperties);
  node->field(INTERN_(kw,molecule_properties),this->_MoleculeProperties);
  this->Base::fields(node);
}

std::string CDNode_O::__repr__() const {
  stringstream ss;
  ss << "#<" << this->className() << " ";
  ss << this->_Label << " " << this->_Pos.asString() << ">";
  return ss.str();
}

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


void CDNode_O::getParsedLabel(string& name, bool& saw_ionization, int& ionization, bool& saw_id, int& id) const
{
  string label = this->_Label;
  vector<string> parts = core::split(label,"/");
  name = parts[0];
  ionization = 0;
  string ionstr = "";
  saw_ionization = false;
  saw_id = false;
  if ( parts.size() == 2 ) {
    ionstr = parts[1];
    for ( size_t i=0; i<ionstr.size(); i++ ) {
      if ( ionstr[i] == '+' ) {
        ionization++;
        saw_ionization = true;
      }
      if ( ionstr[i] == '-' ) {
        ionization--;
        saw_ionization = true;
      }
      if ( ionstr[i]>='0' && ionstr[i]<='9') {
        id = id*10+(ionstr[i]-'0');  // goofy way to parse number +1+2+3 --> +3 123
        saw_id = true;
      }
    }
  }
  LOG("parsed[{}] into name[{}] ionization[{}]" , this->getLabel() , name , ionization );
}

void	CDNode_O::parseFromXml(adapt::QDomNode_sp xml, bool verbose)
{
  this->_Id = xml->getAttributeInt("id");
  this->_Color = xml->getAttributeIntDefault("color",3);
  std::string pos = xml->getAttributeString("p");
  std::vector<std::string> posparts = core::split(pos," ");
  float x = std::stof(posparts[0]) * 0.1;
  float y = std::stof(posparts[1]) * 0.1;
  this->_Pos = geom::Vector2(x,y);
  this->_Label = this->_extractLabel(xml);
  if (verbose) core::clasp_write_string(fmt::format("CDNode id({}) color({}) label({})\n" , this->_Id , this->_Color , this->_Label ));
  LOG("Parsing CDNode with label: {}" , this->_Label);
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
        if (verbose) core::clasp_write_string("    configuration(S)\n");
        LOG("Atom has geometry AS={}" , as );
        this->_Configuration = S_Configuration;
      } else if ( as == "R" )
      {
        LOG("Atom has geometry AS={}" , as );
        if (verbose) core::clasp_write_string("    configuration(R)\n");
        this->_Configuration = R_Configuration;
      } else
      {
        if (verbose) core::clasp_write_string("    could not determin configuration from AS\n");
        LOG("Could not interpret geometry AS[{}]" , as );
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
  this->_BeginNode = nil<CDNode_O>();
  this->_EndNode = nil<CDNode_O>();
}

void CDBond_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,idbegin),this->_IdBegin);
  node->field(INTERN_(kw,idend),this->_IdEnd);
  node->field(INTERN_(kw,nodebegin),this->_BeginNode);
  node->field(INTERN_(kw,nodeend),this->_EndNode);
  node->field(INTERN_(kw,order),this->_Order);
  this->Base::fields(node);
}

std::string CDBond_O::__repr__() const {
  stringstream ss;
  ss << "#<" << this->className() << " ";
  ss << this->_BeginNode->_Label << " " << this->getOrderAsString() << " " << this->_EndNode->_Label << ">";
  return ss.str();
}


string	CDBond_O::getOrderAsString() const
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
  case WedgedHashBeginCDBond: return "WedgedHashBegin";
  case WedgedHashEndCDBond: return "WedgedHashEnd";
  case WedgeBeginCDBond: return "WedgeBegin";
  case WedgeEndCDBond: return "WedgeEnd";
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


void	CDBond_O::parseFromXml(adapt::QDomNode_sp xml, bool verbose)
{
  this->_IdBegin = xml->getAttributeInt("B");
  this->_IdEnd = xml->getAttributeInt("E");
  string order = xml->getAttributeValueDefault("Order","1");
  string display = xml->getAttributeValueDefault("Display","");
  if ( order == "1" ) {
    if ( display == "Dash" ) {
      this->_Order = singleDashCDBond;
    } else if ( display == "Hash" ) {
      this->_Order = hashCDBond;
    } else if ( display == "HollowWedgeBegin" ) {
      this->_Order = hollowWedgeCDBond;
    } else if ( display == "WedgedHashBegin" ) {
      this->_Order = WedgedHashBeginCDBond;
    } else if ( display == "WedgedHashEnd" ) {
      this->_Order = WedgedHashEndCDBond;
    } else if ( display == "WedgeBegin" ) {
      this->_Order = WedgeBeginCDBond;
    } else if ( display == "WedgeEnd" ) {
      this->_Order = WedgeEndCDBond;
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
    Warn(core::SimpleBaseString_O::make(fmt::format("Unknown bond order {}" , (int)this->_Order)),
         nil<core::T_O>());
    this->_Order = unknownCDBond;
  }
  if (verbose) core::clasp_write_string(fmt::format("CDBond _IdBegin({}) _IdEnd({}) order({}) display({})\n" , this->_IdBegin , this->_IdEnd , order , display ));
}



void CDFragment_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,constitutionName),this->_ConstitutionName);
  node->field(INTERN_(kw,nodes),this->_Nodes);
  node->field(INTERN_(kw,atomsToBonds),this->_AtomsToNodes);
  node->field(INTERN_(kw,bonds),this->_Bonds);
  node->field(INTERN_(kw,largestId),this->_LargestId);
  node->field_if_not_nil(INTERN_(kw,molecule),this->_Molecule);
  this->Base::fields(node);
}

void	CDFragment_O::initialize()
{
  this->Base::initialize();
  this->_ConstitutionName = nil<core::Symbol_O>();
  this->_Nodes.clear();
  this->_AtomsToNodes.clear();
  this->_Bonds.clear();
  this->_LargestId = -1;
}

#if 0
void CDFragment_O::createImplicitHydrogen(CDNode_sp fromNode, const std::string& name)
{
  Atom_sp fromAtom = fromNode->getAtom();
  LOG("From {}" , fromAtom->description() );
  Atom_sp a = Atom_O::create();
  a->setContainedBy(_Nil<core::T_O>());
  a->setName(chemkw_intern(name));
  a->setElement(element_H);
  LOG("Created implicit hydrogen ({})" , name  );
  fromAtom->bondTo(a,singleBond);
  auto  toNode  = gctools::GC<CDNode_O>::allocate_with_default_constructor();
  toNode->setAtom(a);
  toNode->_Label = name;
  toNode->setId(this->_LargestId+1);
  this->_LargestId++;
  this->_Nodes[toNode->getId()] = toNode;
  this->_AtomsToNodes[a] = toNode;
  auto  bond  = gctools::GC<CDBond_O>::allocate_with_default_constructor();
  bond->setIdBegin(fromNode->getId());
  bond->setIdEnd(toNode->getId());
  bond->setBeginNode(fromNode);
  bond->setEndNode(toNode);
  bond->setOrder(singleCDBond);
  this->_Bonds.push_back(bond);
}
#endif


void	CDFragment_O::parseFromXml(adapt::QDomNode_sp fragment, bool verbose)
{
  adapt::QDomNode_O::iterator	it;
  this->_Nodes.clear();
  this->_AtomsToNodes.clear();
  if (verbose) core::clasp_write_string("CDFragment - starting\n");
  for ( it=fragment->begin_Children(); it!=fragment->end_Children(); it++ ) {
    adapt::QDomNode_sp child = (*it);
    if ( child->getLocalName() == "n" ) {
      auto  node  = gctools::GC<CDNode_O>::allocate_with_default_constructor();
      node->parseFromXml(child, verbose);
      int id = node->getId();
      if ( id > this->_LargestId )
        this->_LargestId = id;
      this->_Nodes[id] = node;
      LOG("Processed node label({})" , node->getLabel()  );
    }
  }
  for ( it=fragment->begin_Children(); it!=fragment->end_Children(); it++ ) {
    adapt::QDomNode_sp child = (*it);
    if ( child->getLocalName() == "b" ) {
      auto  bond  = gctools::GC<CDBond_O>::allocate_with_default_constructor();
      bond->parseFromXml(child,verbose);
      uint idBegin = bond->getIdBegin();
      uint idEnd = bond->getIdEnd();
      ASSERT(this->_Nodes.count(idBegin)>0);
      ASSERT(this->_Nodes.count(idEnd)>0);
      CDNode_sp nodeBegin = this->_Nodes[idBegin];
      CDNode_sp nodeEnd = this->_Nodes[idEnd];
      bond->setBeginNode( nodeBegin );
      bond->setEndNode( nodeEnd);
      this->_Bonds.push_back(bond);
      LOG("Processed bond order({}) from({}) to ({})"
          , bond->getOrderAsString() , nodeBegin->getLabel() , nodeEnd->getLabel()  );
    }
  }
  if (verbose) core::clasp_write_string("CDFragment - done.\n");
}

CL_DEFMETHOD core::List_sp CDFragment_O::getBonds() const {
  ql::list result;
  for ( size_t idx=0; idx<this->_Bonds.size(); ++idx) {
    result << this->_Bonds[idx];
  }
  return result.cons();
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
void CDFragment_O::addProperties(core::HashTable_sp d)
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
  LOG("tokenizing string [{}]" , label  );
  core::tokenize( label, keyValue, " ");
  if ( keyValue.size() != 2 )
  {
    return false;
  }
  if ( keyValue.size() == 1 )
  {
    SIMPLE_ERROR("There is an error in a ChemDraw file, the string[{}] should be a key:value pair", label);
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
    SIMPLE_ERROR("You have already set the property({}) it has the value({}) and you are trying to set it with({})" , _rep_(key) , _rep_(this->_Properties->gethash(key)) , _rep_(obj)  );
  }
  LOG("Setting variable({}) class({}) to {}" , key , obj->className() , _rep_(obj) );
  set_property(this->_Properties,key,obj);
  return obj;
}



string	CDFragment_O::describeProperties()
{
  stringstream ss;
  this->_Properties->mapHash( [&ss] (core::T_sp key, core::T_sp value) {
      ss << _rep_(key) << "  value=" << _rep_(value) << " ";
    } );
  LOG("{}"%ss.str());
  return ss.str();
}
#endif





core::Symbol_mv parse_property(core::T_sp stream, const string& propertyValue, CDBond_sp bond, const string& otherSideValue)
{
  core::T_sp eof = core::Cons_O::create(nil<core::T_O>(),nil<core::T_O>());
  core::DynamicScopeManager scope(cl::_sym_STARpackageSTAR,_lisp->findPackage(ChemKwPkg));
  core::T_sp property = core::cl__read(stream,nil<core::T_O>(),eof);
  if ( property == eof ) {
    return Values(nil<core::T_O>(),nil<core::T_O>());
  }
//  printf("%s:%d Parsed property: %s\n", __FILE__, __LINE__, _rep_(property).c_str());
  core::T_sp value = core::cl__read(stream,nil<core::T_O>(),eof);
  if ( value == eof ) {
    // If no value is provided then it will default to NIL
    value = nil<core::T_O>();
    // SIMPLE_ERROR(("Could not parse second part of \"{}\" as a (symbol value) pair - in property bond of order {} other side of bond is \"{}\"") , propertyValue , bond->getOrderAsString() , otherSideValue );
  }
//  printf("%s:%d Parsed value: %s\n", __FILE__, __LINE__, _rep_(value).c_str());
  if ( core::Symbol_sp key = property.asOrNull<core::Symbol_O>() ) {
    return Values(property,value);
  }
  SIMPLE_ERROR(("Could not parse \"{}\" as a (symbol value) pair - in property bond of order {} other side of bond is \"{}\"") , propertyValue , bond->getOrderAsString() , otherSideValue );
}






/*!
 * Look for edges that specify properties and move them into the CDNodes that
* they target
 */
bool CDFragment_O::interpret(bool verbose, bool addHydrogens)
{
//  printf("%s:%d  Interpreting a fragment\n", __FILE__, __LINE__ );
  if ( this->_Bonds.size() == 0 ) {return false;}
  CDBonds::iterator bi;
  core::MultipleValues &values = core::lisp_multipleValues();
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
      if ( cdorder == WedgedHashBeginCDBond
           || cdorder == WedgedHashEndCDBond
           || cdorder == WedgeBeginCDBond
           || cdorder == WedgeEndCDBond
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
          SIMPLE_ERROR("The {} bond must have one end with only one neighbor - the one in question has {} and {} neighbors\n" , (*bi)->getOrderAsString() , beginNeighbors , endNeighbors );
        }
        string propertyCode = propertyNode->_Label;
        if (cdorder == dativeCDBond) {
          targetNode->_AtomProperties = core__put_f(targetNode->_AtomProperties,geom::OVector2_O::createFromVector2(propertyNode->_Pos),INTERN_(kw,property_position));
        }
        core::T_sp stream = core::cl__make_string_input_stream(core::Str_O::create(propertyCode),0,nil<core::T_O>());
        core::Symbol_mv parsedProperty;
        do {
          parsedProperty = parse_property(stream, propertyCode, *bi, targetNode->_Label);
          core::T_sp value = values.second(parsedProperty.number_of_values());
          if ( parsedProperty.notnilp() && parsedProperty.number_of_values() == 2 ) {
            if ( cdorder == dativeCDBond ) {
              if (verbose) core::clasp_write_string(fmt::format("Adding atom property {} value: {}\n" , _rep_(parsedProperty) , _rep_(value)));
              targetNode->_AtomProperties = core__put_f(targetNode->_AtomProperties,value,parsedProperty);
            } else if ( cdorder == hollowWedgeCDBond ) {
              if (verbose) core::clasp_write_string(fmt::format("Adding residue property {} value: {}\n" , _rep_(parsedProperty) , _rep_(value)));
              targetNode->_ResidueProperties = core__put_f(targetNode->_ResidueProperties,value,parsedProperty);
            } else if ( cdorder == wavyCDBond ) {
              if (verbose) core::clasp_write_string(fmt::format("Adding molecule property {} value: {}\n" , _rep_(parsedProperty) , _rep_(value)));
              targetNode->_MoleculeProperties = core__put_f(targetNode->_MoleculeProperties,value,parsedProperty);
            } else {
              SIMPLE_ERROR("Cannot interpret bond {} in terms of where to put the property" , (*bi)->getOrderAsString().c_str());
            }
          } else if (parsedProperty.number_of_values() != 2) {
            SIMPLE_ERROR("The property on atom {} must be a list of keyword value pairs - got: {}" , targetNode->_Label , core::_rep_(targetNode->_AtomProperties));
          }
        } while (parsedProperty.notnilp());
        core::cl__close(stream);
      } else {
        Warn(core::Str_O::create(fmt::format("Doing nothing with bond type {}" , (*bi)->getOrderAsString())),
             nil<core::T_O>());
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
#if 1
  if (addHydrogens) {
    mol->fillInImplicitHydrogens();
  }
#else
  core::List_sp carbons = mol->allAtomsOfElementAsList(element_C);
  if (addHydrogens) {
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
  }
#endif

  {
    CipPrioritizer_sp cip = CipPrioritizer_O::create();
    core::HashTable_mv cips = cip->calculateStereochemistryTypeForAllAtoms(mol);
    core::MultipleValues &values = core::lisp_multipleValues();
    core::T_sp ret2 = values.second( cips.number_of_values());
    if (ret2.nilp()) {
      SIMPLE_ERROR("NIL was returned rather than a hash-table for stereochemistryType");
    }
    core::HashTable_sp stereochemistryType = gc::As<core::HashTable_sp>(ret2);
    stereochemistryType->mapHash([](core::T_sp key, core::T_sp val) {
      StereochemistryType st = translate::from_object<StereochemistryType>(val)._v;
        gc::As<Atom_sp>(key)->_StereochemistryType = st;
      }
      );
  }

  this->_Molecule->makeAllAtomNamesInEachResidueUnique();

#if 0

  Residue_sp res = this->getEntireResidue();
  core::List_sp carbons = res->allAtomsOfElementAsList(element_C);
  {
    for ( auto cur : carbons ) {
      Atom_sp c = cur->car<Atom_O>();
      ASSERT(this->_AtomsToNodes.count(c)>0);
      CDNode_sp fromNode = this->_AtomsToNodes[c];
      core::List_sp hydrogens = c->createImplicitHydrogenNamesOnCarbon();
      LOG("Creating implicit hydrogens for {} named {}" , c->description() , _rep_(hydrogens)  );
      for ( auto ih : hydrogens ) {
        core::Symbol_sp name = gc::As<core::Symbol_sp>(oCar(ih));
        this->createImplicitHydrogen(fromNode,name);
      }
    }
  }
#endif

  return true;
}


Atom_sp CDFragment_O::createOneAtom(CDNode_sp n)
{
  Atom_sp a = Atom_O::create();
  Vector3 pos3;
  pos3.set(n->_Pos.getX(),n->_Pos.getY(), 0.0 );
  a->setf_needs_build(false);
  a->setPosition(pos3);
  n->setAtom(a);
  string name;
  bool saw_ionization = false;
  int ionization;
  bool saw_id = false;
  int id = 0;
  n->getParsedLabel(name,saw_ionization,ionization,saw_id,id);
  a->setName(chemkw_intern(name));
  a->setIonization(ionization);
  a->setElementFromAtomName();
  a->setProperty(INTERN_(kw,chemdraw_color),core::clasp_make_fixnum(n->_Color));
  if (saw_id) {
    a->setProperty(INTERN_(kw,id),core::clasp_make_fixnum(id));
  }
  a->setStereochemistryType(n->_StereochemistryType);
  a->setConfiguration(n->_Configuration);
  LOG("Just set configuration of atom[{}] to config[{}]"
      , _rep_(a) , a->getConfigurationAsString() );
#if 0
  if (!elementIsRealElement(a->getElement())) {
    SIMPLE_ERROR("The element for atom {} must be a real element", _rep_(a) );
  }
#endif
  return a;
}

void CDFragment_O::createAtomsAndBonds()
{
	// First create atoms that are on the ends
	// of solid and dashed bonds or the arrow-head of a dative bond
	//
  CDBonds::iterator bi;
  for ( bi=this->_Bonds.begin(); bi!=this->_Bonds.end(); bi++ )
  {
    CDBondOrder o = ( (*bi)->getOrder() );
    gc::Nilable<Atom_sp> beginAtom;
    gc::Nilable<Atom_sp> endAtom;
    if ( o == WedgedHashBeginCDBond || o == WedgedHashEndCDBond
         || o == WedgeBeginCDBond || o == WedgeEndCDBond
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
          bond->setOrder(dashedSingleBond);
          break;
      case WedgedHashBeginCDBond:
//          printf("%s:%d wedgeHashCDBond %s\n", __FILE__, __LINE__, _rep_(bond).c_str());
          bond->setOrder(singleHashBegin);
          break;
      case WedgedHashEndCDBond:
//          printf("%s:%d wedgeHashCDBond %s\n", __FILE__, __LINE__, _rep_(bond).c_str());
          bond->setOrder(singleHashEnd);
          break;
      case WedgeBeginCDBond:
//          printf("%s:%d WedgeBeginCDBond %s\n", __FILE__, __LINE__, _rep_(bond).c_str());
          bond->setOrder(singleWedgeBegin);
          break;
      case WedgeEndCDBond:
//          printf("%s:%d WedgeEndCDBond %s\n", __FILE__, __LINE__, _rep_(bond).c_str());
          bond->setOrder(singleWedgeEnd);
          break;
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
          Warn(core::Str_O::create(fmt::format("Add support for ChemDraw bond: {}" , (*bi)->getOrderAsString())), nil<core::T_O>() );
      }
      beginAtom->addBond(bond);
      endAtom->addBond(bond);
    } else if (o == dativeCDBond) {
      printf("%s:%d:%s  Handle dativeCDBond with beginNode->%s  endNode->%s\n", __FILE__, __LINE__, __FUNCTION__, _rep_((*bi)->getBeginNode()).c_str(), _rep_((*bi)->getEndNode()).c_str());
    }
  }
}


Molecule_sp CDFragment_O::createMolecule()
{
  core::HashTable_sp atomsToNodes = core::HashTable_O::createEq();
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
  core::HashTable_sp dict = this->getProperties();
  if ( !dict->contains(s) )
  {
    stringstream serr;
    serr << "Fragment(" << _rep_(this->getConstitutionName()) << ") is missing property: " << s->__repr__() << std::endl;
    serr << " available property names are: " << std::endl;
    dict->mapHash([&serr](core::T_sp key, core::T_sp val) {
        serr << _rep_(key) << " : " << _rep_(val) << std::endl; } );
    SIMPLE_ERROR("{}", serr.str());
  }
  return dict->gethash(s);
}

core::T_sp	CDFragment_O::getPropertyOrDefault(core::Symbol_sp s, core::T_sp df)
{
  core::HashTable_sp dict = this->getProperties();
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
      LOG("Creating a bond between {} and {} order: {}" , a1->description() , a2->description() , bondOrderToString(bo)  );
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
      LOG("Removing bonds for atom: {}" , a->description()  );
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
  {
    SpanningLoop_sp span = SpanningLoop_O::create();
    span->setOnlyFollowRealBonds(constitutionOnly);
    ASSERTNOTNULL(this->_RootNode);
    ASSERTP(this->_RootNode.notnilp(),"Every fragment must have a root node");
    span->setTop(this->_RootNode->getAtom());
    while ( span->advance() )
    {
      Atom_sp a = span->getAtom();
      a->turnOnFlags(SELECTED);
      LOG("Selected {}" , a->description()  );
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
CL_DEFMETHOD     core::T_sp	CDFragment_O::asConstitutionAtoms()
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


void CDText_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,code),this->_Code);
}

/*!
 * Text blocks should be list of key: value pairs separated by line feeds
 */
bool CDText_O::parseFromXml(adapt::QDomNode_sp text, bool verbose)
{
  core::List_sp xmls = text->childrenWithName("s");
  stringstream ss;
  for ( auto c : xmls ) {
    adapt::QDomNode_sp child = gctools::As<adapt::QDomNode_sp>(oCar(c));
    ss << child->getData();
  }
  string name = ss.str();
  string stext = core::trimWhiteSpace(name);
  if ( stext[0] != '(' )
  {
    LOG("Text block is not code" );
    return false;
  }
  if (verbose) core::clasp_write_string(fmt::format("CDText parsed: {}" , stext));
  core::StringInputStream_sp sin =
    gc::As<core::StringInputStream_sp>(core::cl__make_string_input_stream(core::Str_O::create(stext)
                                                                          ,0
                                                                          ,nil<core::T_O>()));

  core::DynamicScopeManager scope(cl::_sym_STARpackageSTAR,_lisp->findPackage(ChemPkg));
  core::List_sp block = read_lisp_object(sin,true,nil<core::T_O>(),false);
  core::cl__close(sin);
  LOG("Parsed text block: {}\n" , stext);
  if ( block.nilp() ) {
    SIMPLE_ERROR("Error compiling code:\n{}", stext);
  }
  this->_Code = block;
  return true;
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







void ChemDraw_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,allFragments),this->_AllFragments);
  node->field(INTERN_(kw,namedFragments),this->_NamedFragments);
  node->field(INTERN_(kw,code),this->_Code);
}



/*
  __BEGIN_DOC(classes.ChemDraw.!class.ChemDraw)
  \requiredKeyed{fileName:}{Text::name}

  Define a ChemDraw object.  Load a cdxml file from \sa{name} and return the ChemDraw object define by it.
  __END_DOC
*/

CL_LISPIFY_NAME(make-chem-draw);
CL_LAMBDA(file-name &optional verbose (add-hydrogens t));
CL_DOCSTRING(R"dx(Make a chem:chem-draw object from a string.  If verbose is T then print info to *standard-output*.)dx");
DOCGROUP(cando);
CL_DEFUN ChemDraw_sp ChemDraw_O::make(core::T_sp stream, bool verbose, bool addHydrogens)
{
  auto  me  = gctools::GC<ChemDraw_O>::allocate_with_default_constructor();
  me->parse(stream,verbose,addHydrogens); // me->parse(stream);
  return me;
};

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
{
  CDFragment_sp frag = this->_NamedFragments.get(name);
  core::HashTable_sp properties = core::HashTable_O::createEq();
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
  DEPRECATED();
  core::HashTable_sp kargs = core::HashTable_O::createFromKeywordCons(props,validChemdrawKeywords);
  if ( !kargs->contains(_kw_name) )
  {
    stringstream ss;
    ss << "Properties for ChemDraw Fragments must have a (name) property, what was passed has: " << kargs->description();
    SIMPLE_ERROR("{}" , ss.str());
  }
  string name = kargs->lookup(_kw_name).as<core::Str_O>()->get();
  if ( !this->_NamedFragments.contains(name) )
  {
    stringstream ss;
    ss << "There are properties for a fragment named("<<name<<") but there is no ChemDraw fragment with that name";
    SIMPLE_ERROR("{}" , ss.str());
  }
  CDFragment_sp frag = this->_NamedFragments.get(name);
  frag->addProperties(kargs);
}
#endif



void ChemDraw_O::parseChild( adapt::QDomNode_sp child, bool verbose, bool addHydrogens )
{
  if (verbose) core::clasp_write_string(fmt::format("ChemDraw_O::parse child of page with name({})\n" , child->getLocalName()));
  if ( child->getLocalName() == "fragment" ) {
    auto  fragment  = gctools::GC<CDFragment_O>::allocate_with_default_constructor();
    fragment->parseFromXml(child,verbose);
    if ( fragment->interpret(verbose, addHydrogens) ) {
#if 0
      core::HashTable_sp properties = fragment->getProperties();
      if ( !properties->contains(INTERN_(kw,name)))
      {
        SIMPLE_ERROR("Every fragment must have a property({}) available properties: {}" , _rep_(INTERN_(kw,name)) , properties->keysAsString());
      }
      core::Symbol_sp constitutionName = properties->gethash(INTERN_(kw,name)).as<core::Symbol_O>();
      fragment->setConstitutionName(constitutionName);
      this->_NamedFragments.set(constitutionName,fragment);
#endif
      this->_AllFragments.push_back(fragment);
    } else {
      SIMPLE_ERROR("Could not interpret a ChemDraw CDFragment");
    }
  } else if ( child->getLocalName() == "t" ) {
    auto  text  = gctools::GC<CDText_O>::allocate_with_default_constructor();
    if (text->parseFromXml(child,verbose)) {
      this->_Code = core::Cons_O::create(text->_Code,this->_Code);
    }
  } else if ( child->getLocalName() == "group" ) {
    if (verbose) core::clasp_write_string("ChemDraw_O::parsing group start...\n");
    for ( adapt::QDomNode_O::iterator it=child->begin_Children(); it!=child->end_Children(); it++ ) {
      this->parseChild(*it,verbose,addHydrogens);
    }
    if (verbose) core::clasp_write_string("ChemDraw_O::parsing group done.\n");
  }
}


void	ChemDraw_O::parse( core::T_sp strm, bool verbose, bool addHydrogens )
{
  adapt::QDomNode_sp xml = adapt::QDomNode_O::parse(strm);
  if ( !xml->hasChildrenWithName("page") )
    SIMPLE_ERROR("Not a cdxml file");
  adapt::QDomNode_sp page = xml->childWithName("page");
  adapt::QDomNode_O::iterator	it;
  this->_NamedFragments.clear();
  for ( it=page->begin_Children(); it!=page->end_Children(); it++ ) {
    this->parseChild(*it,verbose,addHydrogens);
  }
  if (verbose) core::clasp_write_string("ChemDraw_O::parse done.\n");
}


CL_LISPIFY_NAME("asAggregate");
CL_DEFMETHOD     Aggregate_sp ChemDraw_O::asAggregate()
{
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
  core::List_sp	frags = nil<core::T_O>();
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
  core::List_sp	frags = nil<core::T_O>();
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
    SIMPLE_ERROR("{}" , ss.str());
  }
  return frags;
}


SYMBOL_EXPORT_SC_(ChemPkg,singleCDBond);
SYMBOL_EXPORT_SC_(ChemPkg,doubleCDBond);
SYMBOL_EXPORT_SC_(ChemPkg,tripleCDBond);
SYMBOL_EXPORT_SC_(ChemPkg,dativeCDBond);
SYMBOL_EXPORT_SC_(ChemPkg,singleDashCDBond);
SYMBOL_EXPORT_SC_(ChemPkg,doubleDashCDBond);
SYMBOL_EXPORT_SC_(ChemPkg,tripleDashCDBond);
SYMBOL_EXPORT_SC_(ChemPkg,hashCDBond);
SYMBOL_EXPORT_SC_(ChemPkg,hollowWedgeCDBond);
SYMBOL_EXPORT_SC_(ChemPkg,wedgedHashBeginCDBond);
SYMBOL_EXPORT_SC_(ChemPkg,wedgedHashEndCDBond);
SYMBOL_EXPORT_SC_(ChemPkg,wedgeBeginCDBond);
SYMBOL_EXPORT_SC_(ChemPkg,wedgeEndCDBond);
SYMBOL_EXPORT_SC_(ChemPkg,wavyCDBond);
SYMBOL_EXPORT_SC_(ChemPkg,unknownCDBond);

SYMBOL_EXPORT_SC_(ChemPkg,_PLUS_CDBondOrderConverter_PLUS_);

CL_BEGIN_ENUM(CDBondOrder,_sym__PLUS_CDBondOrderConverter_PLUS_,"CDBondOrder");
CL_VALUE_ENUM(_sym_singleCDBond,singleCDBond);
CL_VALUE_ENUM(_sym_doubleCDBond,doubleCDBond);
CL_VALUE_ENUM(_sym_tripleCDBond,tripleCDBond);
CL_VALUE_ENUM(_sym_dativeCDBond,dativeCDBond);
CL_VALUE_ENUM(_sym_singleDashCDBond,singleDashCDBond);
CL_VALUE_ENUM(_sym_doubleDashCDBond,doubleDashCDBond);
CL_VALUE_ENUM(_sym_tripleDashCDBond,tripleDashCDBond);
CL_VALUE_ENUM(_sym_hashCDBond,hashCDBond);
CL_VALUE_ENUM(_sym_hollowWedgeCDBond,hollowWedgeCDBond);
CL_VALUE_ENUM(_sym_wedgedHashBeginCDBond,WedgedHashBeginCDBond);
CL_VALUE_ENUM(_sym_wedgedHashEndCDBond,WedgedHashEndCDBond);
CL_VALUE_ENUM(_sym_wedgeBeginCDBond,WedgeBeginCDBond);
CL_VALUE_ENUM(_sym_wedgeEndCDBond,WedgeEndCDBond);
CL_VALUE_ENUM(_sym_wavyCDBond,wavyCDBond);
CL_VALUE_ENUM(_sym_unknownCDBond,unknownCDBond);
CL_END_ENUM(_sym__PLUS_CDBondOrderConverter_PLUS_);









};
